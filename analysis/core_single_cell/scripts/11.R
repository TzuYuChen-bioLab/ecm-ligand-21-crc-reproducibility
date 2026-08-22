# 11.R — external recurrence validation with explicit overlap exclusion
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
source(file.path(.script_dir, "helpers_bulk.R"))
require_packages(c("dplyr", "survival", "broom", "metafor", "ggplot2"))

cohorts <- c("GSE14333", "GSE17536", "GSE17537", "GSE33113")
force <- identical(tolower(Sys.getenv("CRC_FORCE_REPROCESS", unset = "false")), "true")
genes <- read_signature_definition()
table_dir <- ensure_dir(project_path("08_tables", "bulk", "corrected_external_validation"))
fig_dir <- ensure_dir(project_path("07_figures", "bulk", "corrected_external_validation"))
diag_dir <- ensure_dir(file.path(table_dir, "diagnostics"))

normalise_sex_external <- function(x) {
  z <- toupper(trimws(as.character(x)))
  factor(ifelse(z %in% c("F", "FEMALE"), "Female", ifelse(z %in% c("M", "MALE"), "Male", NA_character_)),
         levels = c("Female", "Male"))
}

endpoint_spec <- list(
  GSE14333 = list(endpoint = "DFS", time = "parsed_dfs_time", event = "parsed_dfs_cens", censor_inverted = TRUE, days = FALSE,
                  age = c("parsed_age_diag"), sex = c("parsed_gender"), stage = c("parsed_dukesstage")),
  GSE17536 = list(endpoint = "DFS", time = "parsed_dfs_time", event = "parsed_dfs_event_disease_free_survival_cancer_recurrence", censor_inverted = FALSE, days = FALSE,
                  age = c("parsed_age"), sex = c("parsed_gender"), stage = c("parsed_ajcc_stage")),
  GSE17537 = list(endpoint = "DFS", time = "parsed_dfs_time", event = "parsed_dfs_event_disease_free_survival_cancer_recurrence", censor_inverted = FALSE, days = FALSE,
                  age = c("parsed_age"), sex = c("parsed_gender"), stage = c("parsed_ajcc_stage")),
  GSE33113 = list(endpoint = "RFS", time = "parsed_time_to_meta_or_recurrence", event = "parsed_meta_or_recurrence_within_3_years", censor_inverted = FALSE, days = TRUE,
                  age = c("parsed_age_at_diagnosis"), sex = c("parsed_sex"), stage = character(0))
)

make_external_data <- function(gse, score_data) {
  spec <- endpoint_spec[[gse]]
  time_col <- first_existing(score_data, spec$time, required = TRUE, label = paste(gse, "recurrence time"))
  event_col <- first_existing(score_data, spec$event, required = TRUE, label = paste(gse, "recurrence event"))
  age_col <- first_existing(score_data, spec$age)
  sex_col <- first_existing(score_data, spec$sex)
  stage_col <- first_existing(score_data, spec$stage)
  time <- extract_numeric(score_data[[time_col]])
  if (isTRUE(spec$days)) time <- time / 30.4375
  if (isTRUE(spec$censor_inverted)) {
    censor <- extract_numeric(score_data[[event_col]])
    event <- ifelse(censor %in% c(0, 1), 1 - censor, NA_real_)
  } else {
    event <- convert_event(score_data[[event_col]])
  }
  out <- data.frame(
    sample = score_data$sample,
    cohort = gse,
    endpoint = spec$endpoint,
    time = time,
    event = event,
    score_raw = as.numeric(score_data$ECM_LIGAND_21),
    age = if (is.na(age_col)) NA_real_ else extract_numeric(score_data[[age_col]]),
    sex = if (is.na(sex_col)) factor(rep(NA_character_, nrow(score_data))) else normalise_sex_external(score_data[[sex_col]]),
    stage = if (is.na(stage_col)) factor(rep(NA_character_, nrow(score_data))) else simplify_stage(score_data[[stage_col]]),
    stringsAsFactors = FALSE
  )
  out$score_z <- safe_z(out$score_raw)$values
  out <- out[is.finite(out$time) & out$time > 0 & out$event %in% c(0, 1) & is.finite(out$score_z), , drop = FALSE]
  rownames(out) <- out$sample
  out
}

analysis_data <- list()
univ_results <- list()
adjusted_results <- list()
for (gse in cohorts) {
  processed <- process_geo_bulk_corrected(gse, force = force, requantile = FALSE)
  scored <- score_fixed_signature(gse, processed$expr, processed$pheno, genes)
  dat <- make_external_data(gse, scored)
  analysis_data[[gse]] <- dat
  write_csv_atomic(dat, file.path(table_dir, paste0(gse, "_recurrence_analysis_ready.csv")))

  model_name <- paste0(gse, "_", unique(dat$endpoint), "_univariable_per_SD")
  uni <- fit_cox_with_diagnostics(dat, survival::Surv(time, event) ~ score_z, model_name, file.path(diag_dir, gse, "univariable"))
  if (!uni$estimable) stop(gse, " failed the prespecified univariable Cox estimability gate: ", uni$status$reason)
  row <- uni$tidy[uni$tidy$term == "score_z", , drop = FALSE]
  row$cohort <- gse
  row$endpoint <- unique(dat$endpoint)
  row$logHR <- log(row$estimate)
  row$SE_logHR <- row$std.error
  row$model_role <- "primary external per-SD association"
  univ_results[[gse]] <- row

  terms <- character(0)
  if (sum(is.finite(dat$age)) >= 30) terms <- c(terms, "age")
  if (length(unique(stats::na.omit(dat$sex))) >= 2) terms <- c(terms, "sex")
  if (length(unique(stats::na.omit(dat$stage))) >= 2) terms <- c(terms, "stage")
  if (length(terms) > 0) {
    formula <- stats::as.formula(paste("survival::Surv(time, event) ~ score_z +", paste(terms, collapse = " + ")))
    model_vars <- all.vars(formula)
    complete <- droplevels(dat[stats::complete.cases(dat[, model_vars, drop = FALSE]), , drop = FALSE])
    design <- if (nrow(complete) > 0) {
      tryCatch(stats::model.matrix(formula, complete), error = function(e) NULL)
    } else NULL
    parameter_count <- if (is.null(design)) Inf else ncol(design) - 1L
    if (sum(complete$event) >= max(10, 5 * parameter_count)) {
      adj <- fit_cox_with_diagnostics(dat, formula, paste0(gse, "_adjusted_", paste(terms, collapse = "_")), file.path(diag_dir, gse, "adjusted"))
      if (adj$estimable) {
        adj_row <- adj$tidy[adj$tidy$term == "score_z", , drop = FALSE]
        adj_row$cohort <- gse
        adj_row$endpoint <- unique(dat$endpoint)
        adj_row$model_role <- paste("secondary adjusted for", paste(terms, collapse = ", "))
        adjusted_results[[gse]] <- adj_row
      }
    }
  }
}

effects <- dplyr::bind_rows(univ_results)
adjusted <- dplyr::bind_rows(adjusted_results)
write_csv_atomic(effects, file.path(table_dir, "external_ECM_LIGAND_21_recurrence_effects.csv"))
write_csv_atomic(adjusted, file.path(table_dir, "external_ECM_LIGAND_21_adjusted_sensitivity.csv"))

# Detect the known GSE14333/GSE17536 duplication using clinical/endpoint
# fingerprints. The cohorts remain useful, but never enter one meta-analysis.
make_fingerprint <- function(dat, suffix) {
  complete <- dat[is.finite(dat$age) & !is.na(dat$sex) & !is.na(dat$stage), , drop = FALSE]
  fingerprint <- paste(round(complete$age, 3), as.character(complete$sex), as.character(complete$stage),
                       round(complete$time, 2), complete$event, sep = "|")
  out <- data.frame(fingerprint = fingerprint, sample = complete$sample, score = complete$score_raw, stringsAsFactors = FALSE)
  colnames(out)[2:3] <- paste0(c("sample_", "score_"), suffix)
  out
}
fp14333 <- make_fingerprint(analysis_data$GSE14333, "GSE14333")
fp17536 <- make_fingerprint(analysis_data$GSE17536, "GSE17536")
overlap <- merge(fp14333, fp17536, by = "fingerprint", all = FALSE, sort = FALSE)
if (nrow(overlap) == 0) stop("No GSE14333/GSE17536 overlap was detected; inspect endpoint parsing before meta-analysis.")
overlap$score_pair_correlation <- stats::cor(overlap$score_GSE14333, overlap$score_GSE17536, use = "complete.obs")
write_csv_atomic(overlap, file.path(table_dir, "GSE14333_GSE17536_exact_clinical_fingerprint_overlap.csv"))
overlap_summary <- data.frame(
  exact_match_rows = nrow(overlap),
  unique_GSE14333_samples = length(unique(overlap$sample_GSE14333)),
  unique_GSE17536_samples = length(unique(overlap$sample_GSE17536)),
  score_Pearson_r = unique(overlap$score_pair_correlation)[1],
  decision = "mutually exclusive meta-analysis alternatives",
  stringsAsFactors = FALSE
)
write_csv_atomic(overlap_summary, file.path(table_dir, "GSE14333_GSE17536_overlap_summary.csv"))
if (overlap_summary$unique_GSE14333_samples < 100) warning("Fewer than 100 exact overlapping patients detected; verify GEO version/metadata.")

run_meta <- function(effect_data, analysis_name) {
  if (all(c("GSE14333", "GSE17536") %in% effect_data$cohort)) stop("Overlapping GSE14333 and GSE17536 cannot enter the same meta-analysis.")
  if (nrow(effect_data) < 3) stop(analysis_name, " needs at least three independent cohorts.")
  fe <- metafor::rma.uni(yi = effect_data$logHR, sei = effect_data$SE_logHR, method = "FE")
  re <- metafor::rma.uni(yi = effect_data$logHR, sei = effect_data$SE_logHR, method = "REML")
  prediction <- predict(re)
  summary <- data.frame(
    analysis = analysis_name, k = nrow(effect_data), cohorts = paste(effect_data$cohort, collapse = ";"),
    FE_HR = exp(as.numeric(fe$b)), FE_CI_low = exp(fe$ci.lb), FE_CI_high = exp(fe$ci.ub), FE_p = fe$pval,
    RE_HR = exp(as.numeric(re$b)), RE_CI_low = exp(re$ci.lb), RE_CI_high = exp(re$ci.ub), RE_p = re$pval,
    prediction_low = exp(prediction$pi.lb), prediction_high = exp(prediction$pi.ub),
    Q = re$QE, Q_df = re$k - 1, Q_p = re$QEp, I2_percent = re$I2, tau2_REML = re$tau2,
    stringsAsFactors = FALSE
  )
  loo <- as.data.frame(metafor::leave1out(re))
  loo$omitted_cohort <- effect_data$cohort
  loo$analysis <- analysis_name
  list(summary = summary, leave_one_out = loo)
}

primary <- effects[effects$cohort %in% c("GSE17536", "GSE17537", "GSE33113"), , drop = FALSE]
sensitivity <- effects[effects$cohort %in% c("GSE14333", "GSE17537", "GSE33113"), , drop = FALSE]
meta_primary <- run_meta(primary, "Primary: GSE17536 retained; overlapping GSE14333 excluded")
meta_sensitivity <- run_meta(sensitivity, "Sensitivity: GSE14333 replaces overlapping GSE17536")
meta_summaries <- dplyr::bind_rows(meta_primary$summary, meta_sensitivity$summary)
meta_loo <- dplyr::bind_rows(meta_primary$leave_one_out, meta_sensitivity$leave_one_out)

# GSE38832 is a separate sensitivity extension when script 9 produced it.
gse38832_file <- project_path("08_tables", "bulk", "corrected_survival", "GSE38832", "GSE38832_ECM_LIGAND_21_survival_effects.csv")
if (file.exists(gse38832_file)) {
  extra_all <- utils::read.csv(gse38832_file, stringsAsFactors = FALSE)
  extra <- extra_all[extra_all$endpoint == "DFS" & extra_all$term == "score_z", , drop = FALSE]
  if (nrow(extra) == 1) {
    extra$cohort <- "GSE38832"
    extra$logHR <- log(extra$estimate)
    extra$SE_logHR <- extra$std.error
    extended_primary <- run_meta(dplyr::bind_rows(primary, extra), "Extended sensitivity: primary independent cohorts plus GSE38832")
    extended_swap <- run_meta(dplyr::bind_rows(sensitivity, extra), "Extended sensitivity: GSE14333 swap plus GSE38832")
    meta_summaries <- dplyr::bind_rows(meta_summaries, extended_primary$summary, extended_swap$summary)
    meta_loo <- dplyr::bind_rows(meta_loo, extended_primary$leave_one_out, extended_swap$leave_one_out)
  }
}
write_csv_atomic(meta_summaries, file.path(table_dir, "external_recurrence_meta_summary_independent_cohorts.csv"))
write_csv_atomic(meta_loo, file.path(table_dir, "external_recurrence_meta_leave_one_out.csv"))

forest_data <- primary |>
  dplyr::transmute(label = paste(cohort, endpoint), HR = estimate, low = conf.low, high = conf.high, type = "Cohort")
forest_data <- dplyr::bind_rows(
  forest_data,
  data.frame(label = "Random-effects meta", HR = meta_primary$summary$RE_HR,
             low = meta_primary$summary$RE_CI_low, high = meta_primary$summary$RE_CI_high, type = "Meta")
)
forest_data$label <- factor(forest_data$label, levels = rev(forest_data$label))
p <- ggplot2::ggplot(forest_data, ggplot2::aes(HR, label, shape = type)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = low, xmax = high), width = 0.15, orientation = "y") +
  ggplot2::geom_point(size = 2.8) + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
  ggplot2::labs(x = "Hazard ratio per 1 SD ECM_LIGAND_21", y = NULL, shape = NULL,
                title = "Independent external recurrence cohorts")
ggplot2::ggsave(file.path(fig_dir, "external_recurrence_primary_forest.pdf"), p, width = 8, height = 5)
write_session_info(project_path("11_logs", "bulk_corrected", "11_external_validation_sessionInfo.txt"))
message("External validation completed with explicit GSE14333/GSE17536 mutual exclusion and REML meta-analysis.")
