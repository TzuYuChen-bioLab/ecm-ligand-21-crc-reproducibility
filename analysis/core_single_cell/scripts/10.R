# 10.R — estimability-aware GSE39582 recurrence models with full diagnostics
# The primary adjustment set and frozen score are never changed after looking at outcomes.
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
source(file.path(.script_dir, "helpers_bulk.R"))
require_packages(c("dplyr", "survival", "broom", "ggplot2"))

# Keep the special term unqualified in formulas so it is recognized by coxph
# across survival package versions.
strata <- survival::strata

score_file <- project_path("04_bulk_model", "corrected_scores", "GSE39582_ECM_LIGAND_21_scores.rds")
if (!file.exists(score_file)) stop("Run 9.R first: ", score_file)
raw <- readRDS(score_file)
out_dir <- ensure_dir(project_path("08_tables", "bulk", "corrected_survival", "GSE39582"))
fig_dir <- ensure_dir(project_path("07_figures", "bulk", "corrected_survival", "GSE39582"))
diag_path <- file.path(out_dir, "diagnostics")
if (dir.exists(diag_path)) unlink(diag_path, recursive = TRUE, force = TRUE)
diag_dir <- ensure_dir(diag_path)

pick <- function(candidates, label) first_existing(raw, candidates, required = TRUE, label = label)
time_col <- pick(c("parsed_rfs_delay", "rfs.delay:ch1"), "RFS time")
event_col <- pick(c("parsed_rfs_event", "rfs.event:ch1"), "RFS event")
split_col <- pick(c("parsed_dataset", "dataset:ch1"), "discovery/validation split")
age_col <- pick(c("parsed_age_at_diagnosis_year", "age.at.diagnosis (year):ch1"), "age")
sex_col <- pick(c("parsed_sex", "parsed_Sex", "Sex:ch1"), "sex")
stage_col <- pick(c("parsed_tnm_stage", "tnm.stage:ch1"), "TNM stage")
mmr_col <- pick(c("parsed_mmr_status", "mmr.status:ch1"), "MMR status")
site_col <- first_existing(raw, c("parsed_tumor_location", "tumor.location:ch1"))
cms_col <- first_existing(raw, c("parsed_cit_molecularsubtype", "cit.molecularsubtype:ch1"))
braf_col <- first_existing(raw, c("parsed_braf_mutation", "braf.mutation:ch1"))
kras_col <- first_existing(raw, c("parsed_kras_mutation", "kras.mutation:ch1"))
tp53_col <- first_existing(raw, c("parsed_tp53_mutation", "tp53.mutation:ch1"))

normalise_sex <- function(x) {
  z <- toupper(trimws(as.character(x)))
  out <- ifelse(z %in% c("F", "FEMALE", "WOMAN"), "Female",
                ifelse(z %in% c("M", "MALE", "MAN"), "Male", NA_character_))
  factor(out, levels = c("Female", "Male"))
}
normalise_mmr <- function(x) {
  z <- toupper(trimws(as.character(x)))
  out <- ifelse(grepl("DMMR|DEFICIENT|MSI", z), "dMMR",
                ifelse(grepl("PMMR|PROFICIENT|MSS", z), "pMMR", NA_character_))
  factor(out, levels = c("pMMR", "dMMR"))
}
optional_value <- function(column, transform = identity) {
  if (is.na(column)) return(rep(NA, nrow(raw)))
  transform(raw[[column]])
}
clean_optional_category <- function(x) {
  z <- trimws(as.character(x))
  z[toupper(z) %in% c("", "NA", "N/A", "UNKNOWN", "NOT AVAILABLE")] <- NA_character_
  z
}

dat <- data.frame(
  sample = raw$sample,
  cohort_split = tools::toTitleCase(tolower(trimws(as.character(raw[[split_col]])))),
  time = extract_numeric(raw[[time_col]]),
  event = convert_event(raw[[event_col]]),
  score_raw = as.numeric(raw$ECM_LIGAND_21),
  age = extract_numeric(raw[[age_col]]),
  sex = normalise_sex(raw[[sex_col]]),
  stage = simplify_stage(raw[[stage_col]]),
  mmr = normalise_mmr(raw[[mmr_col]]),
  site = factor(optional_value(site_col, clean_optional_category)),
  cms = factor(optional_value(cms_col, clean_optional_category)),
  braf = optional_value(braf_col, simplify_binary_status),
  kras = optional_value(kras_col, simplify_binary_status),
  tp53 = optional_value(tp53_col, simplify_binary_status),
  stringsAsFactors = FALSE
)
if (!all(c("Discovery", "Validation") %in% unique(dat$cohort_split))) {
  stop("Both GSE39582 discovery and validation labels are required.")
}
discovery_score <- dat$score_raw[dat$cohort_split == "Discovery"]
z <- safe_z(dat$score_raw, centre = mean(discovery_score, na.rm = TRUE),
            scale_value = stats::sd(discovery_score, na.rm = TRUE))
dat$score_z <- z$values
dat <- dat[is.finite(dat$time) & dat$time > 0 & dat$event %in% c(0, 1) & is.finite(dat$score_z), , drop = FALSE]
assert_no_duplicate_samples(dat$sample, "GSE39582 RFS analysis data")
dat$cohort_split <- factor(dat$cohort_split, levels = c("Discovery", "Validation"))
rownames(dat) <- make.unique(as.character(dat$sample))
write_csv_atomic(dat, file.path(out_dir, "GSE39582_tumour_only_RFS_analysis_ready.csv"))

registry <- data.frame(
  model = c(
    "univariable", "primary_prespecified", "stratified_sensitivity",
    "secondary_site", "secondary_CMS", "secondary_mutations",
    "pooled_interaction", "pooled_interaction_stratified"
  ),
  formula = c(
    "score", "score + age + sex + TNM stage + MMR",
    "score + age + MMR + strata(sex, TNM stage)",
    "primary + tumour site", "primary + CMS", "primary + BRAF + KRAS + TP53",
    "primary + score*cohort split",
    "score*cohort split + age + MMR + strata(sex, TNM stage)"
  ),
  role = c(
    "supportive", "PRIMARY in discovery; attempted with estimability gate in validation",
    "separation/PH-robust sensitivity", "secondary discovery only", "secondary discovery only",
    "secondary discovery only", "transportability test", "transportability sensitivity"
  ),
  selection_rule = "specified before refitting; never chosen by minimum P value",
  estimability_rule = c(
    "fit if minimum n/events are met",
    "reject output if convergence warning, zero-event coefficient level, non-finite or extreme coefficient/SE is detected",
    "stage and sex define separate baseline hazards; score remains the target coefficient",
    rep("validation not attempted because sparse events cannot support the multi-parameter secondary model", 3),
    "reject output if estimability gate fails",
    "stage and sex define separate baseline hazards; score-by-cohort remains the target interaction"
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(registry, file.path(out_dir, "GSE39582_model_registry_prespecified.csv"))

primary_formula <- survival::Surv(time, event) ~ score_z + age + sex + stage + mmr
univ_formula <- survival::Surv(time, event) ~ score_z
stratified_formula <- survival::Surv(time, event) ~ score_z + age + mmr + strata(sex, stage)
clinical_primary_formula <- survival::Surv(time, event) ~ age + sex + stage + mmr
clinical_stratified_formula <- survival::Surv(time, event) ~ age + mmr + strata(sex, stage)

all_results <- list()
incremental <- list()
model_statuses <- list()
event_audits <- list()
register_fit <- function(key, result) {
  model_statuses[[key]] <<- result$status
  event_audits[[key]] <<- result$event_audit
  invisible(result)
}

for (split in c("Discovery", "Validation")) {
  cohort <- droplevels(dat[dat$cohort_split == split, , drop = FALSE])
  base_name <- paste0("GSE39582_", split)

  uni <- fit_cox_with_diagnostics(
    cohort, univ_formula, paste0(base_name, "_univariable"),
    file.path(diag_dir, split, "univariable")
  )
  register_fit(paste0(split, "_univariable"), uni)
  if (!uni$estimable) stop(split, " cohort did not support the univariable score model.")
  all_results[[paste0(split, "_uni")]] <- transform(
    uni$tidy, cohort_split = split, model_role = "supportive"
  )

  primary <- fit_cox_with_diagnostics(
    cohort, primary_formula, paste0(base_name, "_primary_prespecified"),
    file.path(diag_dir, split, "primary")
  )
  register_fit(paste0(split, "_primary"), primary)
  if (split == "Discovery" && !primary$estimable) {
    stop("Discovery prespecified primary model failed its estimability gate: ", primary$status$reason)
  }
  if (primary$estimable) {
    all_results[[paste0(split, "_primary")]] <- transform(
      primary$tidy, cohort_split = split, model_role = "PRIMARY"
    )
  }

  stratified <- fit_cox_with_diagnostics(
    cohort, stratified_formula, paste0(base_name, "_stratified_sensitivity"),
    file.path(diag_dir, split, "stratified_sensitivity")
  )
  register_fit(paste0(split, "_stratified_sensitivity"), stratified)
  if (!stratified$estimable) {
    stop(split, " stratified sensitivity model failed its estimability gate: ", stratified$status$reason)
  }
  all_results[[paste0(split, "_stratified_sensitivity")]] <- transform(
    stratified$tidy, cohort_split = split,
    model_role = "sensitivity: sex/stage-stratified"
  )

  if (split == "Discovery" && primary$estimable) {
    test_score_nonlinearity(
      cohort, c("age", "sex", "stage", "mmr"), paste0(base_name, "_primary"),
      file.path(diag_dir, split, "primary")
    )
    incremental[[paste0(split, "_primary")]] <- compare_incremental_cox(
      cohort, clinical_primary_formula, primary_formula, split, "prespecified_primary"
    )
  }
  test_score_nonlinearity(
    cohort, c("age", "mmr", "strata(sex, stage)"),
    paste0(base_name, "_stratified_sensitivity"),
    file.path(diag_dir, split, "stratified_sensitivity")
  )
  incremental[[paste0(split, "_stratified")]] <- compare_incremental_cox(
    cohort, clinical_stratified_formula, stratified_formula, split, "stratified_sensitivity"
  )
}

# Secondary models remain exploratory. Validation has only 30 recurrence events
# and contains zero-event levels, so its multi-parameter secondary models are
# recorded as not attempted instead of being coerced to convergence.
secondary_specs <- list(
  secondary_site = survival::Surv(time, event) ~ score_z + age + sex + stage + mmr + site,
  secondary_CMS = survival::Surv(time, event) ~ score_z + age + sex + stage + mmr + cms,
  secondary_mutations = survival::Surv(time, event) ~ score_z + age + sex + stage + mmr + braf + kras + tp53
)
for (split in c("Discovery", "Validation")) {
  cohort <- droplevels(dat[dat$cohort_split == split, , drop = FALSE])
  for (name in names(secondary_specs)) {
    model_name <- paste0("GSE39582_", split, "_", name)
    output_path <- file.path(diag_dir, split, name)
    if (split == "Validation") {
      result <- record_unfitted_cox_model(
        cohort, secondary_specs[[name]], model_name, output_path,
        paste0(
          "not attempted by design in validation: ", sum(cohort$event == 1),
          " recurrence events overall plus sparse/zero-event covariate levels cannot support this multi-parameter secondary model"
        )
      )
    } else {
      vars <- all.vars(secondary_specs[[name]])
      factor_vars <- vars[!vars %in% c("time", "event", "score_z", "age")]
      available <- vapply(factor_vars, function(v) {
        length(unique(stats::na.omit(cohort[[v]]))) >= 2
      }, logical(1))
      if (all(available)) {
        result <- fit_cox_with_diagnostics(cohort, secondary_specs[[name]], model_name, output_path)
      } else {
        result <- record_unfitted_cox_model(
          cohort, secondary_specs[[name]], model_name, output_path,
          "not attempted: at least one required categorical covariate has fewer than two observed levels"
        )
      }
      if (result$estimable) {
        all_results[[paste(split, name, sep = "_")]] <- transform(
          result$tidy, cohort_split = split, model_role = "secondary discovery-only"
        )
      }
    }
    register_fit(paste(split, name, sep = "_"), result)
  }
}

interaction_formula <- survival::Surv(time, event) ~ score_z * cohort_split + age + sex + stage + mmr
interaction <- fit_cox_with_diagnostics(
  dat, interaction_formula, "GSE39582_pooled_score_by_cohort_interaction",
  file.path(diag_dir, "pooled_interaction")
)
register_fit("pooled_interaction", interaction)
if (interaction$estimable) {
  all_results$interaction <- transform(
    interaction$tidy, cohort_split = "Pooled", model_role = "transportability"
  )
}

interaction_stratified_formula <- survival::Surv(time, event) ~
  score_z * cohort_split + age + mmr + strata(sex, stage)
interaction_stratified <- fit_cox_with_diagnostics(
  dat, interaction_stratified_formula, "GSE39582_pooled_score_by_cohort_interaction_stratified",
  file.path(diag_dir, "pooled_interaction_stratified")
)
register_fit("pooled_interaction_stratified", interaction_stratified)
if (!interaction_stratified$estimable) {
  stop("Pooled stratified interaction model failed its estimability gate: ", interaction_stratified$status$reason)
}
all_results$interaction_stratified <- transform(
  interaction_stratified$tidy, cohort_split = "Pooled",
  model_role = "transportability sensitivity: sex/stage-stratified"
)

results <- dplyr::bind_rows(all_results)
status_table <- dplyr::bind_rows(model_statuses)
event_audit_table <- dplyr::bind_rows(event_audits)
incremental_table <- dplyr::bind_rows(incremental)
write_csv_atomic(results, file.path(out_dir, "GSE39582_all_estimable_Cox_results.csv"))
# Retain the historical filename for downstream compatibility, but it now
# contains estimable models only and includes model_status/model_role columns.
write_csv_atomic(results, file.path(out_dir, "GSE39582_all_prespecified_Cox_results.csv"))
write_csv_atomic(status_table, file.path(out_dir, "GSE39582_model_estimability_status.csv"))
write_csv_atomic(event_audit_table, file.path(out_dir, "GSE39582_model_event_level_audit.csv"))
write_csv_atomic(incremental_table, file.path(out_dir, "GSE39582_incremental_value_clinical_vs_score.csv"))

# A discovery median is carried unchanged into validation for descriptive KM only.
cutpoint <- stats::median(dat$score_raw[dat$cohort_split == "Discovery"], na.rm = TRUE)
dat$risk_group_frozen <- factor(ifelse(dat$score_raw >= cutpoint, "High", "Low"), levels = c("Low", "High"))
km_source <- list()
for (split in c("Discovery", "Validation")) {
  cohort_data <- dat[dat$cohort_split == split, , drop = FALSE]
  fit <- survival::survfit(survival::Surv(time, event) ~ risk_group_frozen, data = cohort_data)
  s <- summary(fit)
  km_source[[split]] <- data.frame(
    cohort_split = split, time = s$time, n_risk = s$n.risk, n_event = s$n.event,
    survival = s$surv, lower = s$lower, upper = s$upper,
    risk_group = base::sub("risk_group_frozen=", "", s$strata), stringsAsFactors = FALSE
  )
}
km_source <- dplyr::bind_rows(km_source)
write_csv_atomic(km_source, file.path(out_dir, "GSE39582_frozen_discovery_median_KM_source.csv"))
write_csv_atomic(
  data.frame(cutpoint_source = "Discovery", raw_score_cutpoint = cutpoint, purpose = "descriptive only"),
  file.path(out_dir, "GSE39582_frozen_KM_cutpoint.csv")
)
p_km <- ggplot2::ggplot(km_source, ggplot2::aes(time, survival, colour = risk_group)) +
  ggplot2::geom_step(linewidth = 0.8) + ggplot2::facet_wrap(~ cohort_split) + ggplot2::theme_bw() +
  ggplot2::labs(
    x = "RFS time (months)", y = "Recurrence-free survival", colour = "ECM_LIGAND_21",
    title = "Descriptive KM using the discovery-set median in both cohorts"
  )
ggplot2::ggsave(file.path(fig_dir, "GSE39582_frozen_cutpoint_KM.pdf"), p_km, width = 9, height = 5)
write_session_info(project_path("11_logs", "bulk_corrected", "10_GSE39582_diagnostics_sessionInfo.txt"))
message(
  "GSE39582 estimability-aware Cox models completed. Review GSE39582_model_estimability_status.csv ",
  "before running step 11."
)
