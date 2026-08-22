# 9.R — tumour-only bulk processing, locked scoring, and GSE38832 diagnostics
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
source(file.path(.script_dir, "helpers_bulk.R"))
require_packages(c("dplyr", "survival", "broom", "ggplot2"))

genes <- read_signature_definition()
if (length(genes) != 21) stop("The locked ECM_LIGAND_21 signature must contain exactly 21 unique genes.")
force <- identical(tolower(Sys.getenv("CRC_FORCE_REPROCESS", unset = "false")), "true")

processed <- lapply(c("GSE39582", "GSE38832"), function(gse) process_geo_bulk_corrected(gse, force = force, requantile = FALSE))
names(processed) <- c("GSE39582", "GSE38832")

gse39582 <- processed$GSE39582
split_col <- first_existing(gse39582$pheno, c("parsed_dataset", "dataset:ch1"), required = TRUE, label = "GSE39582 cohort split")
discovery_samples <- rownames(gse39582$pheno)[tolower(trimws(as.character(gse39582$pheno[[split_col]]))) == "discovery"]
if (length(discovery_samples) != 443) stop("Expected 443 GSE39582 discovery tumours; found ", length(discovery_samples), ".")
score39582 <- score_fixed_signature("GSE39582", gse39582$expr, gse39582$pheno, genes, reference_samples = discovery_samples)
score38832 <- score_fixed_signature("GSE38832", processed$GSE38832$expr, processed$GSE38832$pheno, genes)

endpoint_map <- list(
  DFS = c("parsed_dfs_time_disease_free_survival_time_months", "parsed_dfs_event_disease_free_survival"),
  DSS = c("parsed_dss_time_disease_specific_survival_time_months", "parsed_dss_event_disease_specific_survival")
)
diagnostic_dir <- ensure_dir(project_path("08_tables", "bulk", "corrected_survival", "GSE38832"))
effects <- list()
for (endpoint in names(endpoint_map)) {
  time_col <- endpoint_map[[endpoint]][1]
  event_col <- endpoint_map[[endpoint]][2]
  if (!all(c(time_col, event_col) %in% colnames(score38832))) stop("GSE38832 ", endpoint, " columns are missing.")
  dat <- data.frame(
    sample = score38832$sample,
    time = extract_numeric(score38832[[time_col]]),
    event = convert_event(score38832[[event_col]]),
    score_raw = as.numeric(score38832$ECM_LIGAND_21),
    stringsAsFactors = FALSE
  )
  z <- safe_z(dat$score_raw)
  dat$score_z <- z$values
  dat <- dat[is.finite(dat$time) & dat$time > 0 & dat$event %in% c(0, 1) & is.finite(dat$score_z), , drop = FALSE]
  model_name <- paste0("GSE38832_", endpoint, "_univariable_per_SD")
  fit_result <- fit_cox_with_diagnostics(
    dat, survival::Surv(time, event) ~ score_z, model_name,
    ensure_dir(file.path(diagnostic_dir, endpoint))
  )
  if (!fit_result$estimable) {
    stop(model_name, " failed the Cox estimability gate: ", fit_result$status$reason)
  }
  test_score_nonlinearity(dat, character(0), model_name, file.path(diagnostic_dir, endpoint))
  row <- fit_result$tidy[fit_result$tidy$term == "score_z", , drop = FALSE]
  row$cohort <- "GSE38832"
  row$endpoint <- endpoint
  row$logHR <- log(row$estimate)
  row$SE_logHR <- row$std.error
  effects[[endpoint]] <- row
}
effects <- dplyr::bind_rows(effects)
write_csv_atomic(effects, file.path(diagnostic_dir, "GSE38832_ECM_LIGAND_21_survival_effects.csv"))

audit <- data.frame(
  cohort = c("GSE39582", "GSE38832"),
  samples_scored = c(nrow(score39582), nrow(score38832)),
  population = c("443 discovery + 123 validation primary colon cancers; 19 non-tumoral mucosa excluded before mapping/scoring",
                 "primary colorectal tumour samples after phenotype audit"),
  expression_standardisation = c("gene means/SD and probe selection learned in the 443-sample discovery set",
                                 "within-cohort gene means/SD"),
  repeated_quantile_normalisation = FALSE,
  stringsAsFactors = FALSE
)
write_csv_atomic(audit, project_path("08_tables", "bulk", "corrected_scores", "bulk_processing_population_audit.csv"))
write_session_info(project_path("11_logs", "bulk_corrected", "9_bulk_processing_scoring_sessionInfo.txt"))
message("Tumour-only processing and locked ECM_LIGAND_21 scoring completed.")
