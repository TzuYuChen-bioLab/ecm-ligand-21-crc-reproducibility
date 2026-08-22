# 5.R — construct matched, condition-specific communication cohorts
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
require_packages(c("Seurat", "dplyr", "ggplot2", "tidyr"))

in_file <- project_path("02_scRNA_analysis", "GSE132465_corrected", "GSE132465_author_annotated_corrected.rds")
out_dir <- ensure_dir(project_path("03_cellchat_nichenet", "GSE132465_corrected"))
table_dir <- ensure_dir(project_path("08_tables", "GSE132465_corrected", "communication_prepare"))
seu <- readRDS(in_file)
Seurat::DefaultAssay(seu) <- "RNA"
rna_data <- tryCatch(get_assay_data_compat(seu, "RNA", "data"), error = function(e) NULL)
if (is.null(rna_data) || nrow(rna_data) == 0 || ncol(rna_data) == 0) {
  seu <- Seurat::NormalizeData(seu, assay = "RNA", verbose = FALSE)
}

seu$comm_group <- NA_character_
stromal_map <- c("Myofibroblasts" = "Stromal_Myofibroblasts", "Stromal 1" = "Stromal_1",
                 "Stromal 2" = "Stromal_2", "Stromal 3" = "Stromal_3")
myeloid_map <- c("SPP1+" = "Myeloid_SPP1pos", "Pro-inflammatory" = "Myeloid_Proinflammatory",
                 "cDC" = "Myeloid_cDC", "Proliferating" = "Myeloid_Proliferating")
for (subtype in names(stromal_map)) seu$comm_group[seu$major_celltype == "Fibroblast_stromal" & seu$Cell_subtype == subtype] <- stromal_map[[subtype]]
for (subtype in names(myeloid_map)) seu$comm_group[seu$major_celltype == "Myeloid" & seu$Cell_subtype == subtype] <- myeloid_map[[subtype]]
seu$comm_group[seu$major_celltype == "Epithelial"] <- "Epithelial"
seu_comm <- subset(seu, cells = colnames(seu)[!is.na(seu$comm_group)])

# Define pairing from the receiver population itself so both CellChat conditions
# contain the same patients with adequate epithelial coverage.
epi_counts <- as.data.frame(table(patient = seu_comm$Patient, condition = seu_comm$condition,
                                  is_epithelial = seu_comm$comm_group == "Epithelial"), stringsAsFactors = FALSE)
colnames(epi_counts)[4] <- "n_cells"
epi_counts$is_epithelial <- as.character(epi_counts$is_epithelial) == "TRUE"
paired_patients <- epi_counts |>
  dplyr::filter(is_epithelial, n_cells >= 20) |>
  dplyr::group_by(patient) |>
  dplyr::filter(all(c("Tumor", "Normal") %in% condition)) |>
  dplyr::ungroup() |>
  dplyr::pull(patient) |>
  unique() |>
  as.character()
if (length(paired_patients) < 8) stop("Fewer than 8 paired tumour-normal patients; found ", length(paired_patients), ".")
seu_paired <- subset(seu_comm, cells = colnames(seu_comm)[seu_comm$Patient %in% paired_patients])

counts_condition <- as.data.frame(table(seu_paired$comm_group, seu_paired$condition))
colnames(counts_condition) <- c("comm_group", "condition", "n_cells")
wide <- tidyr::pivot_wider(counts_condition, names_from = condition, values_from = n_cells, values_fill = 0)
common_groups <- wide$comm_group[wide$Tumor >= 20 & wide$Normal >= 20]
if (!"Epithelial" %in% common_groups) stop("Epithelial cells do not meet the minimum in both paired conditions.")
seu_paired <- subset(seu_paired, cells = colnames(seu_paired)[seu_paired$comm_group %in% common_groups])
seu_paired$comm_group <- factor(seu_paired$comm_group, levels = sort(common_groups))

write_csv_atomic(data.frame(patient = sort(paired_patients)), file.path(table_dir, "paired_patients.csv"))
write_csv_atomic(epi_counts, file.path(table_dir, "epithelial_pairing_audit.csv"))
write_csv_atomic(wide, file.path(table_dir, "group_counts_by_condition_before_common_filter.csv"))
write_csv_atomic(as.data.frame(table(seu_paired$comm_group, seu_paired$condition, seu_paired$Patient)),
                 file.path(table_dir, "group_counts_by_condition_patient.csv"))
saveRDS(seu_paired, file.path(out_dir, "GSE132465_paired_condition_communication_seurat.rds"))
saveRDS(subset(seu_paired, cells = colnames(seu_paired)[seu_paired$condition == "Tumor"]),
        file.path(out_dir, "GSE132465_Tumor_communication_seurat.rds"))
saveRDS(subset(seu_paired, cells = colnames(seu_paired)[seu_paired$condition == "Normal"]),
        file.path(out_dir, "GSE132465_Normal_communication_seurat.rds"))
message("Prepared separate tumour and normal communication datasets for ", length(paired_patients), " paired patients.")
