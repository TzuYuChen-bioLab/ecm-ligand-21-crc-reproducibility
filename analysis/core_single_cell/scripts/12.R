# 12.R — publication figures from corrected bulk analyses only
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
require_packages(c("dplyr", "ggplot2", "patchwork", "scales"))

gse_dir <- project_path("08_tables", "bulk", "corrected_survival", "GSE39582")
external_dir <- project_path("08_tables", "bulk", "corrected_external_validation")
fig_dir <- ensure_dir(project_path("07_figures", "publication_bulk_validation_corrected"))
source_dir <- ensure_dir(project_path("08_tables", "publication_bulk_validation_corrected_source_data"))

read_required <- function(path) {
  if (!file.exists(path)) stop("Required corrected result is missing: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
save_both <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggplot2::ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 350)
}

cox <- read_required(file.path(gse_dir, "GSE39582_all_prespecified_Cox_results.csv"))
incremental <- read_required(file.path(gse_dir, "GSE39582_incremental_value_clinical_vs_score.csv"))
km <- read_required(file.path(gse_dir, "GSE39582_frozen_discovery_median_KM_source.csv"))
status <- read_required(file.path(gse_dir, "GSE39582_model_estimability_status.csv"))
discovery_primary <- cox |>
  dplyr::filter(term == "score_z", model_role == "PRIMARY", cohort_split == "Discovery") |>
  dplyr::mutate(analysis_label = "Discovery primary")
stratified_hr <- cox |>
  dplyr::filter(term == "score_z", model_role == "sensitivity: sex/stage-stratified") |>
  dplyr::mutate(analysis_label = paste(cohort_split, "stratified"))
primary_hr <- dplyr::bind_rows(discovery_primary, stratified_hr) |>
  dplyr::transmute(analysis_label, cohort = cohort_split, HR = estimate, low = conf.low, high = conf.high,
                   p_value = p.value, n = n, events = n_event, C_index_apparent = concordance,
                   model_role = model_role)
expected_hr_labels <- c("Discovery primary", "Discovery stratified", "Validation stratified")
if (!setequal(primary_hr$analysis_label, expected_hr_labels) || nrow(primary_hr) != 3) {
  stop("Expected Discovery primary plus Discovery/Validation stratified score effects.")
}
if (any(status$status == "not_estimable" & grepl("Validation_primary_prespecified", status$model))) {
  message("Validation ordinary primary model is correctly excluded from Figure 5 because it is not estimable.")
}

ph_files <- list.files(file.path(gse_dir, "diagnostics"), pattern = "_PH_test\\.csv$", recursive = TRUE, full.names = TRUE)
ph <- dplyr::bind_rows(lapply(ph_files, function(path) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  x$file <- basename(path)
  x
}))
primary_ph <- ph[
  grepl("Discovery_primary_prespecified|Discovery_stratified_sensitivity|Validation_stratified_sensitivity", ph$model) &
    ph$term %in% c("score_z", "sex", "GLOBAL"), , drop = FALSE
]
p_col <- grep("^p$|p.value|Pr", colnames(primary_ph), value = TRUE, ignore.case = TRUE)[1]
if (is.na(p_col)) primary_ph$PH_p <- NA_real_ else primary_ph$PH_p <- as.numeric(primary_ph[[p_col]])
primary_ph$analysis_label <- ifelse(
  grepl("Discovery_primary_prespecified", primary_ph$model), "Discovery primary",
  ifelse(grepl("Discovery_stratified_sensitivity", primary_ph$model), "Discovery stratified", "Validation stratified")
)
incremental_plot <- incremental |>
  dplyr::filter(status == "estimated") |>
  dplyr::mutate(
    analysis_label = ifelse(
      analysis_type == "prespecified_primary", "Discovery primary",
      paste(cohort_split, "stratified")
    )
  )
if (!setequal(incremental_plot$analysis_label, expected_hr_labels) || nrow(incremental_plot) != 3) {
  stop("Expected three estimable incremental-value rows for Figure 5.")
}

write_csv_atomic(primary_hr, file.path(source_dir, "Figure5_primary_prespecified_HR.csv"))
write_csv_atomic(incremental, file.path(source_dir, "Figure5_incremental_value.csv"))
write_csv_atomic(km, file.path(source_dir, "Figure5_frozen_cutpoint_KM.csv"))
write_csv_atomic(primary_ph, file.path(source_dir, "Figure5_PH_diagnostics.csv"))
write_csv_atomic(status, file.path(source_dir, "Figure5_model_estimability_status.csv"))

p5a <- ggplot2::ggplot(km, ggplot2::aes(time, survival, colour = risk_group)) +
  ggplot2::geom_step(linewidth = 0.8) + ggplot2::facet_wrap(~ cohort_split) + ggplot2::theme_bw() +
  ggplot2::labs(x = "RFS time (months)", y = "Recurrence-free survival", colour = "Frozen group",
                title = "A  Discovery median carried into validation")
primary_hr$analysis_label <- factor(primary_hr$analysis_label, levels = rev(expected_hr_labels))
p5b <- ggplot2::ggplot(primary_hr, ggplot2::aes(HR, analysis_label)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = low, xmax = high), width = 0.14, orientation = "y") +
  ggplot2::geom_point(size = 3, colour = "#B2182B") + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
  ggplot2::labs(x = "Adjusted HR per 1 SD", y = NULL, title = "B  Primary and separation-robust sensitivity")
incremental_plot$analysis_label <- factor(incremental_plot$analysis_label, levels = expected_hr_labels)
p5c <- ggplot2::ggplot(incremental_plot, ggplot2::aes(analysis_label, delta_C_index_apparent, fill = cohort_split)) +
  ggplot2::geom_col(width = 0.65, show.legend = FALSE) + ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
  ggplot2::labs(x = NULL, y = "Apparent Δ C-index", title = "C  Incremental discrimination")
primary_ph$analysis_label <- factor(primary_ph$analysis_label, levels = expected_hr_labels)
p5d <- ggplot2::ggplot(primary_ph, ggplot2::aes(analysis_label, -log10(pmax(PH_p, .Machine$double.xmin)), shape = term)) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey55") +
  ggplot2::geom_point(size = 3) + ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
  ggplot2::labs(x = NULL, y = expression(-log[10](italic(P)[PH])), shape = NULL,
                title = "D  PH checks (sex stratified in sensitivity models)")
figure5 <- (p5a / (p5b + p5c + p5d)) + patchwork::plot_annotation(title = "Figure 5. GSE39582 tumour-only recurrence analysis")
save_both(figure5, "Figure5_GSE39582_corrected", 14, 10)

external <- read_required(file.path(external_dir, "external_ECM_LIGAND_21_recurrence_effects.csv"))
meta <- read_required(file.path(external_dir, "external_recurrence_meta_summary_independent_cohorts.csv"))
overlap <- read_required(file.path(external_dir, "GSE14333_GSE17536_overlap_summary.csv"))
primary_external <- external[external$cohort %in% c("GSE17536", "GSE17537", "GSE33113"), , drop = FALSE]
if (any(c("GSE14333", "GSE17536") %in% primary_external$cohort) && all(c("GSE14333", "GSE17536") %in% primary_external$cohort)) {
  stop("Overlapping GSE14333 and GSE17536 entered the Figure 6 primary set.")
}
primary_meta <- meta[grepl("^Primary:", meta$analysis), , drop = FALSE]
forest <- primary_external |>
  dplyr::transmute(label = paste(cohort, endpoint), HR = estimate, low = conf.low, high = conf.high, type = "Independent cohort")
forest <- dplyr::bind_rows(forest, data.frame(label = "REML random effects", HR = primary_meta$RE_HR,
                                              low = primary_meta$RE_CI_low, high = primary_meta$RE_CI_high, type = "Meta-analysis"))
forest$label <- factor(forest$label, levels = rev(forest$label))
write_csv_atomic(forest, file.path(source_dir, "Figure6_independent_external_forest.csv"))
write_csv_atomic(meta, file.path(source_dir, "Figure6_meta_primary_and_sensitivity.csv"))
write_csv_atomic(overlap, file.path(source_dir, "Figure6_overlap_exclusion_audit.csv"))

p6a <- ggplot2::ggplot(forest, ggplot2::aes(HR, label, shape = type)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = low, xmax = high), width = 0.14, orientation = "y") +
  ggplot2::geom_point(size = 3) + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
  ggplot2::labs(x = "HR per 1 SD ECM_LIGAND_21", y = NULL, shape = NULL,
                title = "A  Independent external recurrence cohorts")
meta_plot <- meta[grepl("^Primary:|^Sensitivity:", meta$analysis), , drop = FALSE]
meta_plot$short <- ifelse(grepl("^Primary", meta_plot$analysis), "Primary: GSE17536", "Sensitivity: GSE14333 swap")
meta_plot$short <- factor(meta_plot$short, levels = rev(meta_plot$short))
p6b <- ggplot2::ggplot(meta_plot, ggplot2::aes(RE_HR, short)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = RE_CI_low, xmax = RE_CI_high), width = 0.14, orientation = "y") +
  ggplot2::geom_point(size = 3, colour = "#2166AC") + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
  ggplot2::labs(x = "REML pooled HR", y = NULL, title = "B  Mutually exclusive overlap sensitivity")
figure6 <- (p6a | p6b) + patchwork::plot_annotation(title = "Figure 6. External validation without duplicate patients")
save_both(figure6, "Figure6_external_validation_corrected", 13, 6)

inventory <- data.frame(
  figure = c("Figure 5", "Figure 6"),
  population = c("GSE39582 primary tumours only", "independent external recurrence cohorts"),
  key_guardrail = c("Discovery primary plus estimability-gated sex/stage-stratified sensitivity and full diagnostics",
                    "GSE14333 and GSE17536 mutually exclusive; one recurrence endpoint per cohort"),
  stringsAsFactors = FALSE
)
write_csv_atomic(inventory, file.path(source_dir, "bulk_figure_method_inventory.csv"))
write_session_info(project_path("11_logs", "bulk_corrected", "12_publication_figures_sessionInfo.txt"))
message("Corrected bulk publication figures and source data completed.")
