# 6.R — condition-specific CellChat plus paired donor expression support
# CellChat is run separately in tumour and normal tissue, following the official
# comparison workflow. No significance-filtered missing edge is replaced by zero.
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
source(file.path(.script_dir, "helpers_cellchat.R"))
require_packages(c("Seurat", "CellChat", "dplyr", "tidyr", "ggplot2", "patchwork", "future", "Matrix"))
future::plan("sequential")
options(future.globals.maxSize = 20 * 1024^3)

in_dir <- project_path("03_cellchat_nichenet", "GSE132465_corrected")
out_dir <- ensure_dir(file.path(in_dir, "cellchat"))
table_dir <- ensure_dir(project_path("08_tables", "GSE132465_corrected", "cellchat"))
fig_dir <- ensure_dir(project_path("07_figures", "GSE132465_corrected", "cellchat"))

seu_all <- readRDS(file.path(in_dir, "GSE132465_paired_condition_communication_seurat.rds"))
seu_tumour <- readRDS(file.path(in_dir, "GSE132465_Tumor_communication_seurat.rds"))
seu_normal <- readRDS(file.path(in_dir, "GSE132465_Normal_communication_seurat.rds"))
paired_patients <- sort(intersect(unique(as.character(seu_tumour$Patient)), unique(as.character(seu_normal$Patient))))
common_groups <- intersect(levels(droplevels(factor(seu_tumour$comm_group))), levels(droplevels(factor(seu_normal$comm_group))))
source_groups <- setdiff(common_groups, "Epithelial")

cellchat_tumour <- run_cellchat_condition(seu_tumour, "Tumor", min_cells = 20)
cellchat_normal <- run_cellchat_condition(seu_normal, "Normal", min_cells = 20)
object_list <- list(Normal = cellchat_normal, Tumor = cellchat_tumour)
cellchat_merged <- CellChat::mergeCellChat(object_list, add.names = names(object_list))
saveRDS(cellchat_tumour, file.path(out_dir, "CellChat_Tumor_paired.rds"))
saveRDS(cellchat_normal, file.path(out_dir, "CellChat_Normal_paired.rds"))
saveRDS(cellchat_merged, file.path(out_dir, "CellChat_Normal_Tumor_merged.rds"))

raw_tumour <- extract_cellchat_raw_probabilities(cellchat_tumour) |>
  dplyr::filter(source %in% source_groups, target == "Epithelial") |>
  dplyr::rename(prob_tumor = prob, pval_tumor = pval)
raw_normal <- extract_cellchat_raw_probabilities(cellchat_normal) |>
  dplyr::filter(source %in% source_groups, target == "Epithelial") |>
  dplyr::rename(prob_normal = prob, pval_normal = pval)
join_keys <- intersect(c("source", "target", "interaction_name", "ligand", "receptor", "pathway_name", "annotation"),
                       intersect(colnames(raw_tumour), colnames(raw_normal)))
comparison <- dplyr::full_join(raw_tumour, raw_normal, by = join_keys) |>
  dplyr::mutate(
    prob_difference = dplyr::if_else(!is.na(prob_tumor) & !is.na(prob_normal), prob_tumor - prob_normal, NA_real_),
    significant_tumor = !is.na(pval_tumor) & pval_tumor < 0.05,
    significant_normal = !is.na(pval_normal) & pval_normal < 0.05,
    comparison_status = dplyr::case_when(
      !is.na(prob_tumor) & !is.na(prob_normal) ~ "estimated_in_both_conditions",
      !is.na(prob_tumor) ~ "tumor_only_estimate_missing_normal",
      !is.na(prob_normal) ~ "normal_only_estimate_missing_tumor",
      TRUE ~ "missing_both"
    ),
    evidence_scope = "paired_patients_but_pooled_condition_specific_CellChat"
  ) |>
  dplyr::arrange(dplyr::desc(prob_difference))
write_csv_atomic(comparison, file.path(table_dir, "GSE132465_CellChat_condition_comparison_raw.csv"))

lr_union <- unique(dplyr::bind_rows(
  as.data.frame(cellchat_tumour@LR$LRsig, stringsAsFactors = FALSE) |> tibble::rownames_to_column("interaction_name_row"),
  as.data.frame(cellchat_normal@LR$LRsig, stringsAsFactors = FALSE) |> tibble::rownames_to_column("interaction_name_row")
))
if (!"interaction_name" %in% colnames(lr_union)) lr_union$interaction_name <- lr_union$interaction_name_row
donor_tests <- paired_lr_expression_tests(seu_all, lr_union, source_groups, paired_patients)
write_csv_atomic(donor_tests, file.path(table_dir, "GSE132465_donor_paired_LR_expression_tests.csv"))

support_keys <- intersect(c("source", "target", "interaction_name", "ligand", "receptor", "pathway_name"),
                          intersect(colnames(comparison), colnames(donor_tests)))
supported <- dplyr::left_join(comparison, donor_tests, by = support_keys) |>
  dplyr::mutate(
    donor_supported = n_pairs >= 5 & !is.na(FDR) & FDR < 0.10 & median_difference > 0,
    interpretation = "CellChat probabilities are pooled exploratory estimates; paired donor statistics refer to an LR expression score, not a CellChat probability test."
  )
write_csv_atomic(supported, file.path(table_dir, "GSE132465_CellChat_donor_supported_comparison.csv"))

pdf(file.path(fig_dir, "CellChat_condition_network_comparison.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2), xpd = TRUE)
CellChat::netVisual_diffInteraction(cellchat_merged, comparison = c(1, 2), weight.scale = TRUE)
CellChat::netVisual_diffInteraction(cellchat_merged, comparison = c(1, 2), weight.scale = TRUE, measure = "weight")
dev.off()

p <- CellChat::netVisual_bubble(cellchat_merged, sources.use = source_groups, targets.use = "Epithelial",
                                comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE)
ggplot2::ggsave(file.path(fig_dir, "CellChat_LR_comparison_bubble.pdf"), p, width = 12, height = 7)
write_session_info(project_path("11_logs", "GSE132465_corrected", "6_CellChat_sessionInfo.txt"))
message("Condition-specific CellChat and paired donor LR-expression support completed.")
