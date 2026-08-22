# 13.R — publication figures from corrected single-cell analyses only
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
require_packages(c("Seurat", "dplyr", "tidyr", "ggplot2", "patchwork", "scales"))

fig_dir <- ensure_dir(project_path("07_figures", "publication_single_cell_corrected"))
source_dir <- ensure_dir(project_path("08_tables", "publication_single_cell_corrected_source_data"))
annotation_file <- project_path("02_scRNA_analysis", "GSE132465_corrected", "GSE132465_author_annotated_corrected.rds")
communication_file <- project_path("03_cellchat_nichenet", "GSE132465_corrected", "GSE132465_paired_condition_communication_seurat.rds")
cellchat_file <- project_path("08_tables", "GSE132465_corrected", "cellchat", "GSE132465_CellChat_donor_supported_comparison.csv")
deg_file <- project_path("08_tables", "GSE132465_corrected", "nichenet", "GSE132465_epithelial_paired_pseudobulk_DEG.csv")
activity_file <- project_path("08_tables", "GSE132465_corrected", "nichenet", "GSE132465_NicheNet_ligand_activities_paired_pseudobulk.csv")
evidence_file <- project_path("08_tables", "GSE132465_corrected", "signature", "CRC_ECM_LIGAND_21_single_cell_evidence.csv")

required <- c(annotation_file, communication_file, cellchat_file, deg_file, activity_file, evidence_file)
if (any(!file.exists(required))) stop("Missing corrected input(s): ", paste(required[!file.exists(required)], collapse = ", "))
save_both <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggplot2::ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 350)
}
seu <- readRDS(annotation_file)
comm <- readRDS(communication_file)
if (!"umap" %in% names(seu@reductions)) stop("Corrected annotated object has no UMAP reduction.")
umap <- as.data.frame(Seurat::Embeddings(seu, "umap"))
umap$cell <- rownames(umap)
meta <- seu@meta.data
meta$cell <- rownames(meta)
umap <- dplyr::left_join(umap, meta[, c("cell", "Patient", "condition", "major_celltype", "epithelial_provenance")], by = "cell")
colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
counts <- as.data.frame(table(patient = comm$Patient, condition = comm$condition, group = comm$comm_group), stringsAsFactors = FALSE)
colnames(counts)[4] <- "n_cells"
write_csv_atomic(umap, file.path(source_dir, "Figure2_UMAP_source.csv"))
write_csv_atomic(counts, file.path(source_dir, "Figure2_paired_group_counts.csv"))

p2a <- ggplot2::ggplot(umap, ggplot2::aes(UMAP_1, UMAP_2, colour = major_celltype)) +
  ggplot2::geom_point(size = 0.12, alpha = 0.65) + ggplot2::theme_void() +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
  ggplot2::labs(colour = NULL, title = "A  Author-guided major cell types")
p2b <- ggplot2::ggplot(umap, ggplot2::aes(UMAP_1, UMAP_2, colour = condition)) +
  ggplot2::geom_point(size = 0.12, alpha = 0.6) + ggplot2::facet_wrap(~ condition) + ggplot2::theme_void() +
  ggplot2::scale_colour_manual(values = c("Normal" = "#2166AC", "Tumor" = "#B2182B")) +
  ggplot2::labs(colour = NULL, title = "B  Specimen condition (not malignant-cell status)")
p2c <- counts |>
  dplyr::filter(n_cells > 0) |>
  ggplot2::ggplot(ggplot2::aes(condition, n_cells, colour = condition)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.12, alpha = 0.55) +
  ggplot2::facet_wrap(~ group, scales = "free_y") + ggplot2::theme_bw() +
  ggplot2::scale_colour_manual(values = c("Normal" = "#2166AC", "Tumor" = "#B2182B"), guide = "none") +
  ggplot2::labs(x = NULL, y = "Cells per paired patient", title = "C  Paired communication-cell coverage")
figure2 <- (p2a | p2b) / p2c + patchwork::plot_annotation(title = "Figure 2. Single-cell landscape with explicit tissue provenance")
save_both(figure2, "Figure2_single_cell_landscape_corrected", 14, 10)

cellchat <- utils::read.csv(cellchat_file, stringsAsFactors = FALSE, check.names = FALSE)
cellchat_plot <- cellchat |>
  dplyr::filter(target == "Epithelial", is.finite(median_difference), is.finite(prob_difference)) |>
  dplyr::mutate(
    donor_support = ifelse(donor_supported %in% TRUE, "Paired-donor FDR<0.10", "Exploratory only"),
    axis = paste(source, interaction_name, sep = " → ")
  ) |>
  dplyr::arrange(dplyr::desc(donor_supported), dplyr::desc(median_difference)) |>
  dplyr::slice_head(n = 40)
write_csv_atomic(cellchat_plot, file.path(source_dir, "Figure3_condition_specific_CellChat_paired_donor_support.csv"))

p3a <- ggplot2::ggplot(cellchat_plot, ggplot2::aes(prob_difference, median_difference, colour = donor_support)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey65") +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey65") +
  ggplot2::geom_point(ggplot2::aes(size = n_pairs), alpha = 0.8) + ggplot2::theme_bw() +
  ggplot2::scale_colour_manual(values = c("Paired-donor FDR<0.10" = "#B2182B", "Exploratory only" = "grey65")) +
  ggplot2::labs(x = "Tumour − normal pooled CellChat probability",
                y = "Paired-donor median LR-expression-score difference", colour = NULL, size = "Patient pairs",
                title = "A  Separate-condition CellChat with donor-level support")
top_supported <- cellchat_plot |>
  dplyr::filter(donor_support == "Paired-donor FDR<0.10") |>
  dplyr::slice_max(median_difference, n = 20, with_ties = FALSE)
if (nrow(top_supported) == 0) top_supported <- cellchat_plot |> dplyr::slice_max(median_difference, n = 20, with_ties = FALSE)
top_supported$axis <- factor(top_supported$axis, levels = rev(top_supported$axis))
p3b <- ggplot2::ggplot(top_supported, ggplot2::aes(median_difference, axis, colour = donor_support)) +
  ggplot2::geom_point(size = 2.5) + ggplot2::theme_bw() +
  ggplot2::scale_colour_manual(values = c("Paired-donor FDR<0.10" = "#B2182B", "Exploratory only" = "grey65"), guide = "none") +
  ggplot2::labs(x = "Paired median difference", y = NULL,
                title = "B  Ligand–receptor axes (receiver: epithelial)")
figure3 <- (p3a | p3b) + patchwork::plot_annotation(title = "Figure 3. Tumour-versus-normal communication analysis")
save_both(figure3, "Figure3_CellChat_corrected", 14, 7)

deg <- utils::read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
activity <- utils::read.csv(activity_file, stringsAsFactors = FALSE, check.names = FALSE)
evidence <- utils::read.csv(evidence_file, stringsAsFactors = FALSE, check.names = FALSE)
cell_cycle_pattern <- paste(c("^MKI67$", "^TOP2A$", "^UBE2C$", "^AURK[AB]$", "^BIRC5$", "^CDK1$", "^CDC",
                              "^CCNA", "^CCNB", "^CCND", "^MCM", "^PCNA$", "^NUSAP1$", "^CENP", "^PLK1$",
                              "^TK1$", "^TYMS$", "^HMMR$", "^KIF", "^TPX2$"), collapse = "|")
deg$is_cell_cycle_like <- grepl(cell_cycle_pattern, deg$gene)
activity_top <- activity |>
  dplyr::slice_max(pearson, n = 20, with_ties = FALSE)
evidence_plot <- evidence[!is.na(evidence$pearson), , drop = FALSE]
write_csv_atomic(deg, file.path(source_dir, "Figure4_paired_pseudobulk_DEG.csv"))
write_csv_atomic(activity_top, file.path(source_dir, "Figure4_NicheNet_top_ligands.csv"))
write_csv_atomic(evidence, file.path(source_dir, "Figure4_fixed_signature_evidence.csv"))

p4a <- ggplot2::ggplot(deg, ggplot2::aes(logFC, -log10(pmax(FDR, .Machine$double.xmin)))) +
  ggplot2::geom_point(ggplot2::aes(colour = FDR < 0.05 & abs(logFC) > 0.25), alpha = 0.55, size = 1) +
  ggplot2::scale_colour_manual(values = c("grey72", "#B2182B"), guide = "none") + ggplot2::theme_bw() +
  ggplot2::labs(x = "Tumour vs normal log2FC", y = expression(-log[10](FDR)),
                title = "A  Paired patient-level pseudobulk")
activity_top$test_ligand <- factor(activity_top$test_ligand, levels = rev(activity_top$test_ligand[order(activity_top$pearson)]))
p4b <- ggplot2::ggplot(activity_top, ggplot2::aes(pearson, test_ligand, fill = candidate_tier)) +
  ggplot2::geom_col() + ggplot2::theme_bw() + ggplot2::labs(x = "NicheNet Pearson activity", y = NULL, fill = NULL,
                                                            title = "B  Ligand activity")
p4c <- ggplot2::ggplot(evidence_plot, ggplot2::aes(stats::reorder(gene, pearson), pearson,
                                                   fill = n_donor_supported_edges > 0)) +
  ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(values = c("grey70", "#2166AC"), na.value = "grey85") +
  ggplot2::labs(x = NULL, y = "NicheNet Pearson activity", fill = "Donor-supported\nCellChat edge",
                title = "C  Descriptive evidence for locked 21-gene score")
figure4 <- (p4a | p4b | p4c) + patchwork::plot_annotation(title = "Figure 4. Pseudobulk-to-NicheNet mechanism and fixed score")
save_both(figure4, "Figure4_pseudobulk_NicheNet_signature_corrected", 17, 7)

inventory <- data.frame(
  figure = c("Figure 2", "Figure 3", "Figure 4"),
  experimental_unit = c("cell for visualization; patient coverage shown", "paired patient for statistical support", "paired patient pseudobulk"),
  prohibited_claim = c("tumour-tissue epithelial is not proof of malignancy", "LR-expression test is not a CellChat probability test",
                       "signature evidence is descriptive; survival did not select genes"),
  stringsAsFactors = FALSE
)
write_csv_atomic(inventory, file.path(source_dir, "single_cell_figure_method_inventory.csv"))
write_session_info(project_path("11_logs", "GSE132465_corrected", "13_publication_figures_sessionInfo.txt"))
message("Corrected single-cell publication figures and source data completed.")
