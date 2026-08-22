# 8.R — freeze a prespecified, bulk-compatible 21-gene ECM ligand score
# The signature membership and equal weights are fixed before survival testing.
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
require_packages(c("dplyr", "ggplot2"))

signature_dir <- ensure_dir(project_path("03_cellchat_nichenet", "GSE132465_corrected", "signature"))
table_dir <- ensure_dir(project_path("08_tables", "GSE132465_corrected", "signature"))
fig_dir <- ensure_dir(project_path("07_figures", "GSE132465_corrected", "signature"))
nichenet_dir <- project_path("08_tables", "GSE132465_corrected", "nichenet")
cellchat_dir <- project_path("08_tables", "GSE132465_corrected", "cellchat")

genes <- c(
  "COL1A1", "COL1A2", "COL4A1", "COL4A2", "COL4A5", "COL6A1", "COL6A2", "COL6A3",
  "FN1", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMC1", "THBS1", "THBS2", "HSPG2",
  "TNC", "TNXB", "PTN", "MDK"
)
signature <- data.frame(
  signature = "ECM_LIGAND_21", gene = genes, weight = 1 / length(genes),
  role = "extracellular_matrix_or_growth_factor_ligand", locked_before_survival = TRUE,
  stringsAsFactors = FALSE
)
write_csv_atomic(signature, file.path(signature_dir, "CRC_ECM_LIGAND_21.csv"))
write_csv_atomic(signature, file.path(table_dir, "CRC_ECM_LIGAND_21.csv"))
writeLines(paste(c("ECM_LIGAND_21", "fixed_equal_weight_z_mean", genes), collapse = "\t"),
           file.path(signature_dir, "CRC_ECM_LIGAND_21.gmt"))

activities <- utils::read.csv(file.path(nichenet_dir, "GSE132465_NicheNet_ligand_activities_paired_pseudobulk.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
support <- utils::read.csv(file.path(cellchat_dir, "GSE132465_CellChat_donor_supported_comparison.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
act_cols <- intersect(c("test_ligand", "pearson", "aupr_corrected", "candidate_tier"), colnames(activities))
activity_gene <- activities[, act_cols, drop = FALSE]
colnames(activity_gene)[colnames(activity_gene) == "test_ligand"] <- "gene"
support_gene <- support |>
  dplyr::filter(ligand %in% genes) |>
  dplyr::group_by(gene = ligand) |>
  dplyr::summarise(
    n_CellChat_edges = dplyr::n(),
    n_donor_supported_edges = sum(donor_supported %in% TRUE, na.rm = TRUE),
    max_paired_donor_median_difference = suppressWarnings(max(median_difference, na.rm = TRUE)),
    max_pooled_CellChat_probability_difference = suppressWarnings(max(prob_difference, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with("max_"), ~ ifelse(is.infinite(.x), NA_real_, .x)))
evidence <- signature |>
  dplyr::left_join(activity_gene, by = "gene") |>
  dplyr::left_join(support_gene, by = "gene") |>
  dplyr::mutate(
    evidence_is_descriptive = TRUE,
    membership_changed_by_current_data = FALSE,
    interpretation = "Evidence columns describe support; they did not select or reweight signature genes."
  )
write_csv_atomic(evidence, file.path(table_dir, "CRC_ECM_LIGAND_21_single_cell_evidence.csv"))

metadata <- data.frame(
  field = c("signature_name", "version", "freeze_date", "n_genes", "weights", "bulk_score",
            "missing_gene_rule", "training_on_survival", "primary_bulk_endpoint", "scope"),
  value = c("ECM_LIGAND_21", "2.0.0", "2026-07-12", length(genes), "equal",
            "mean of per-gene z-scores", "require at least 80% gene coverage and use available locked genes",
            "none", "recurrence/DFS per-SD continuous score", "prognostic association; not a diagnostic classifier"),
  stringsAsFactors = FALSE
)
write_csv_atomic(metadata, file.path(signature_dir, "CRC_ECM_LIGAND_21_metadata.csv"))
write_csv_atomic(metadata, file.path(table_dir, "CRC_ECM_LIGAND_21_metadata.csv"))

plot_data <- evidence[!is.na(evidence$pearson), , drop = FALSE]
if (nrow(plot_data) > 0) {
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(stats::reorder(gene, pearson), pearson)) +
    ggplot2::geom_col(fill = "#2166AC") + ggplot2::coord_flip() + ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = "NicheNet Pearson activity", title = "Descriptive support for fixed ECM_LIGAND_21 genes")
  ggplot2::ggsave(file.path(fig_dir, "ECM_LIGAND_21_NicheNet_support.pdf"), p, width = 7, height = 6)
}
write_session_info(project_path("11_logs", "GSE132465_corrected", "8_signature_freeze_sessionInfo.txt"))
message("Frozen ECM_LIGAND_21 definition written with 21 equal-weight genes.")
