# 2.R — metadata alignment and QC
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
require_packages(c("Seurat", "data.table", "ggplot2"))

in_dir <- project_path("01_processed_data", "GSE132465_corrected")
out_dir <- ensure_dir(project_path("02_scRNA_analysis", "GSE132465_corrected"))
fig_dir <- ensure_dir(project_path("07_figures", "GSE132465_corrected", "QC"))
seu <- readRDS(file.path(in_dir, "GSE132465_raw_seurat.rds"))
annotation <- as.data.frame(readRDS(file.path(in_dir, "GSE132465_annotation.rds")), stringsAsFactors = FALSE)

cell_candidates <- c("cell", "Cell", "cell_id", "Cell_ID", "barcode", "Barcode", "CellName", "name")
cell_col <- cell_candidates[cell_candidates %in% colnames(annotation)][1]
if (is.na(cell_col)) cell_col <- colnames(annotation)[1]
rownames(annotation) <- as.character(annotation[[cell_col]])
common <- intersect(colnames(seu), rownames(annotation))
if (length(common) < 0.9 * ncol(seu)) stop("Fewer than 90% of count-matrix cells matched annotation.")
seu <- subset(seu, cells = common)
seu <- Seurat::AddMetaData(seu, metadata = annotation[common, , drop = FALSE])

required <- c("Class", "Sample", "Patient", "Cell_type", "Cell_subtype")
missing <- setdiff(required, colnames(seu@meta.data))
if (length(missing) > 0) stop("Author annotation is missing: ", paste(missing, collapse = ", "))
seu$condition <- normalise_condition(seu$Class)
if (anyNA(seu$condition)) stop("Some Class values could not be mapped to Tumor/Normal.")
seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")

qc_before <- data.frame(cell = colnames(seu), patient = seu$Patient, condition = seu$condition,
                        nCount_RNA = seu$nCount_RNA, nFeature_RNA = seu$nFeature_RNA, percent_mt = seu$percent.mt)
keep <- seu$nFeature_RNA > 200 & seu$nFeature_RNA < 6000 & seu$nCount_RNA > 500 & seu$percent.mt < 20
seu <- subset(seu, cells = colnames(seu)[keep])
qc_after <- data.frame(cell = colnames(seu), patient = seu$Patient, condition = seu$condition,
                       nCount_RNA = seu$nCount_RNA, nFeature_RNA = seu$nFeature_RNA, percent_mt = seu$percent.mt)
write_csv_atomic(qc_before, file.path(out_dir, "QC_before.csv"))
write_csv_atomic(qc_after, file.path(out_dir, "QC_after.csv"))
write_csv_atomic(as.data.frame.matrix(table(seu$Patient, seu$condition)), file.path(out_dir, "cells_by_patient_condition.csv"), row.names = TRUE)

p <- Seurat::VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "condition", ncol = 3, pt.size = 0)
ggplot2::ggsave(file.path(fig_dir, "QC_after_by_condition.pdf"), p, width = 11, height = 4)
saveRDS(seu, file.path(out_dir, "GSE132465_QC_filtered_seurat.rds"))
