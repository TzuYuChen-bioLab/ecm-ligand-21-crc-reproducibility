# 7.R — patient-level pseudobulk differential expression and NicheNet
# The experimental unit is the patient, not the cell. Tumour and normal
# epithelial pseudobulks are analysed with a paired edgeR design.
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
require_packages(c("Seurat", "Matrix", "edgeR", "dplyr", "tidyr", "tibble", "ggplot2", "nichenetr"), bioconductor = TRUE)

in_dir <- project_path("03_cellchat_nichenet", "GSE132465_corrected")
table_dir <- ensure_dir(project_path("08_tables", "GSE132465_corrected", "nichenet"))
fig_dir <- ensure_dir(project_path("07_figures", "GSE132465_corrected", "nichenet"))
resource_dir <- project_path("03_cellchat_nichenet", "nichenet_resources")
seu <- readRDS(file.path(in_dir, "GSE132465_paired_condition_communication_seurat.rds"))
Seurat::DefaultAssay(seu) <- "RNA"

required_meta <- c("Patient", "condition", "comm_group")
if (length(setdiff(required_meta, colnames(seu@meta.data))) > 0) stop("Patient/condition/comm_group metadata is missing.")
epi_cells <- colnames(seu)[seu$comm_group == "Epithelial"]
epi <- subset(seu, cells = epi_cells)
counts <- get_assay_data_compat(epi, "RNA", "counts")
meta <- epi@meta.data[colnames(counts), , drop = FALSE]

min_cells <- as.integer(Sys.getenv("CRC_MIN_EPITHELIAL_CELLS_PER_PSEUDOBULK", unset = "20"))
sample_key <- paste(meta$Patient, meta$condition, sep = "__")
cell_counts <- as.data.frame(table(patient = meta$Patient, condition = meta$condition), stringsAsFactors = FALSE)
colnames(cell_counts)[3] <- "n_cells"
eligible <- cell_counts |>
  dplyr::filter(n_cells >= min_cells) |>
  dplyr::group_by(patient) |>
  dplyr::filter(all(c("Normal", "Tumor") %in% condition)) |>
  dplyr::ungroup()
paired_patients <- sort(unique(as.character(eligible$patient)))
if (length(paired_patients) < 5) {
  stop("Only ", length(paired_patients), " paired patients have at least ", min_cells,
       " epithelial cells in each condition; five are required.")
}
keep_cells <- rownames(meta)[as.character(meta$Patient) %in% paired_patients]
meta <- meta[keep_cells, , drop = FALSE]
counts <- counts[, keep_cells, drop = FALSE]
sample_key <- paste(meta$Patient, meta$condition, sep = "__")
split_cells <- split(seq_len(ncol(counts)), sample_key)
pseudobulk <- vapply(split_cells, function(index) as.numeric(Matrix::rowSums(counts[, index, drop = FALSE])), numeric(nrow(counts)))
rownames(pseudobulk) <- rownames(counts)

sample_meta <- do.call(rbind, strsplit(colnames(pseudobulk), "__", fixed = TRUE))
sample_meta <- data.frame(sample = colnames(pseudobulk), patient = sample_meta[, 1], condition = sample_meta[, 2], stringsAsFactors = FALSE)
sample_meta$condition <- factor(sample_meta$condition, levels = c("Normal", "Tumor"))
sample_meta$patient <- factor(sample_meta$patient)
sample_meta <- sample_meta[order(sample_meta$patient, sample_meta$condition), , drop = FALSE]
pseudobulk <- pseudobulk[, sample_meta$sample, drop = FALSE]
assert_no_duplicate_samples(sample_meta$sample, "Epithelial pseudobulk")

design <- stats::model.matrix(~ patient + condition, data = sample_meta)
y <- edgeR::DGEList(counts = pseudobulk, samples = sample_meta)
keep_gene <- edgeR::filterByExpr(y, design = design)
if (sum(keep_gene) < 1000) stop("Fewer than 1,000 genes passed pseudobulk expression filtering.")
y <- edgeR::calcNormFactors(y[keep_gene, , keep.lib.sizes = FALSE], method = "TMM")
y <- edgeR::estimateDisp(y, design = design, robust = TRUE)
fit <- edgeR::glmQLFit(y, design = design, robust = TRUE)
coef_name <- "conditionTumor"
if (!coef_name %in% colnames(design)) stop("Paired design did not create conditionTumor coefficient.")
test <- edgeR::glmQLFTest(fit, coef = which(colnames(design) == coef_name))
deg <- edgeR::topTags(test, n = Inf, sort.by = "PValue")$table |>
  tibble::rownames_to_column("gene")

# Cell-level detection rates are descriptive annotations only and are not used
# as replicates or for inferential P values.
norm <- get_assay_data_compat(epi, "RNA", "data")
pct_tumor <- Matrix::rowMeans(norm[, epi$condition == "Tumor", drop = FALSE] > 0)
pct_normal <- Matrix::rowMeans(norm[, epi$condition == "Normal", drop = FALSE] > 0)
deg$pct_tumor_cells_descriptive <- as.numeric(pct_tumor[deg$gene])
deg$pct_normal_cells_descriptive <- as.numeric(pct_normal[deg$gene])
deg$experimental_unit <- "paired_patient_pseudobulk"
write_csv_atomic(deg, file.path(table_dir, "GSE132465_epithelial_paired_pseudobulk_DEG.csv"))
write_csv_atomic(cell_counts, file.path(table_dir, "GSE132465_epithelial_cells_per_patient_condition.csv"))
write_csv_atomic(sample_meta, file.path(table_dir, "GSE132465_epithelial_pseudobulk_sample_design.csv"))

geneset_oi <- deg$gene[deg$FDR < 0.05 & deg$logFC > 0.25]
if (length(geneset_oi) < 10) {
  stop("Only ", length(geneset_oi), " tumour-upregulated genes met prespecified FDR/logFC thresholds; NicheNet was not run.")
}

support_file <- project_path("08_tables", "GSE132465_corrected", "cellchat", "GSE132465_CellChat_donor_supported_comparison.csv")
support <- utils::read.csv(support_file, stringsAsFactors = FALSE, check.names = FALSE)
support_ligand <- support |>
  dplyr::filter(!is.na(ligand), nzchar(ligand)) |>
  dplyr::group_by(ligand) |>
  dplyr::summarise(
    donor_supported = any(donor_supported %in% TRUE, na.rm = TRUE),
    max_pairs = suppressWarnings(max(n_pairs, na.rm = TRUE)),
    max_median_difference = suppressWarnings(max(median_difference, na.rm = TRUE)),
    max_pooled_probability_difference = suppressWarnings(max(prob_difference, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with("max_"), ~ ifelse(is.infinite(.x), NA_real_, .x)))

ligand_target_file <- file.path(resource_dir, "ligand_target_matrix.rds")
lr_network_file <- file.path(resource_dir, "lr_network.rds")
missing_resources <- c(ligand_target_file, lr_network_file)[!file.exists(c(ligand_target_file, lr_network_file))]
if (length(missing_resources) > 0) {
  stop("Missing local NicheNet resource(s): ", paste(missing_resources, collapse = ", "),
       ". Download/version them explicitly; this pipeline does not silently fetch changing resources.")
}
ligand_target_matrix <- readRDS(ligand_target_file)
lr_network <- as.data.frame(readRDS(lr_network_file), stringsAsFactors = FALSE)
if (!is.matrix(ligand_target_matrix)) ligand_target_matrix <- as.matrix(ligand_target_matrix)

sender_cells <- colnames(seu)[seu$condition == "Tumor" & seu$comm_group != "Epithelial"]
receiver_cells <- colnames(seu)[seu$condition == "Tumor" & seu$comm_group == "Epithelial"]
norm_all <- get_assay_data_compat(seu, "RNA", "data")
sender_expressed <- names(which(Matrix::rowMeans(norm_all[, sender_cells, drop = FALSE] > 0) >= 0.10))
receiver_expressed <- names(which(Matrix::rowMeans(norm_all[, receiver_cells, drop = FALSE] > 0) >= 0.10))
background <- intersect(receiver_expressed, rownames(ligand_target_matrix))
geneset_oi <- intersect(geneset_oi, background)
if (length(geneset_oi) < 10) stop("Fewer than 10 prespecified DE genes overlap the NicheNet target background.")

tier1 <- support_ligand$ligand[support_ligand$donor_supported]
tier2 <- support_ligand$ligand[
  support_ligand$max_pairs >= 5 & support_ligand$max_median_difference > 0 &
    support_ligand$max_pooled_probability_difference > 0
]
ligand_col <- c("from", "ligand", "source")[c("from", "ligand", "source") %in% colnames(lr_network)][1]
receptor_col <- c("to", "receptor", "target")[c("to", "receptor", "target") %in% colnames(lr_network)][1]
if (is.na(ligand_col) || is.na(receptor_col)) stop("NicheNet lr_network lacks recognizable ligand/receptor columns.")
feasible_lr_ligands <- unique(as.character(lr_network[[ligand_col]][
  as.character(lr_network[[ligand_col]]) %in% sender_expressed &
    as.character(lr_network[[receptor_col]]) %in% receiver_expressed
]))
feasible <- intersect(feasible_lr_ligands, colnames(ligand_target_matrix))
tier1_feasible <- intersect(tier1, feasible)
tier2_feasible <- intersect(tier2, feasible)
candidate_tier <- if (length(tier1_feasible) >= 3) "tier1_donor_supported" else "tier2_directionally_concordant_exploratory"
potential_ligands <- if (candidate_tier == "tier1_donor_supported") tier1_feasible else tier2_feasible
if (length(potential_ligands) < 3) stop("Fewer than three CellChat-supported, sender-expressed NicheNet ligands remain.")

activities <- nichenetr::predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
) |>
  dplyr::left_join(support_ligand, by = c("test_ligand" = "ligand")) |>
  dplyr::mutate(candidate_tier = candidate_tier, DEG_definition = "paired pseudobulk FDR<0.05 and logFC>0.25") |>
  dplyr::arrange(dplyr::desc(pearson))
write_csv_atomic(activities, file.path(table_dir, "GSE132465_NicheNet_ligand_activities_paired_pseudobulk.csv"))

top_ligands <- head(activities$test_ligand, 20)
links <- expand.grid(target = geneset_oi, ligand = top_ligands, stringsAsFactors = FALSE)
links$regulatory_potential <- mapply(function(target, ligand) ligand_target_matrix[target, ligand], links$target, links$ligand)
links <- links[is.finite(links$regulatory_potential) & links$regulatory_potential > 0, , drop = FALSE]
links <- links |>
  dplyr::group_by(ligand) |>
  dplyr::slice_max(regulatory_potential, n = 100, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::left_join(deg[, c("gene", "logFC", "FDR")], by = c("target" = "gene"))
write_csv_atomic(links, file.path(table_dir, "GSE132465_NicheNet_ligand_target_links_paired_pseudobulk.csv"))
write_csv_atomic(data.frame(gene = geneset_oi), file.path(table_dir, "GSE132465_NicheNet_tumor_up_geneset_paired_pseudobulk.csv"))

p_volcano <- ggplot2::ggplot(deg, ggplot2::aes(logFC, -log10(pmax(FDR, .Machine$double.xmin)))) +
  ggplot2::geom_point(ggplot2::aes(colour = FDR < 0.05 & abs(logFC) > 0.25), alpha = 0.55, size = 1) +
  ggplot2::scale_colour_manual(values = c("grey70", "#B2182B"), guide = "none") +
  ggplot2::theme_bw() + ggplot2::labs(x = "Tumour vs normal log2 fold change", y = expression(-log[10](FDR)),
                                      title = "Paired epithelial pseudobulk differential expression")
ggplot2::ggsave(file.path(fig_dir, "paired_pseudobulk_DEG_volcano.pdf"), p_volcano, width = 7, height = 6)
saveRDS(list(dge = y, design = design, fit = fit, test = test, sample_meta = sample_meta),
        file.path(in_dir, "GSE132465_epithelial_paired_pseudobulk_edgeR.rds"))
write_session_info(project_path("11_logs", "GSE132465_corrected", "7_pseudobulk_NicheNet_sessionInfo.txt"))
message("Paired patient-level pseudobulk and NicheNet analysis completed for ", length(paired_patients), " patients.")
