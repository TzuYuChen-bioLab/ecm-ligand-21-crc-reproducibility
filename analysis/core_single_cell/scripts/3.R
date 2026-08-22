# 3.R — normalisation, PCA, UMAP and clustering
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
require_packages(c("Seurat", "ggplot2", "future"))
future::plan("sequential")
options(future.globals.maxSize = 20 * 1024^3)

out_dir <- ensure_dir(project_path("02_scRNA_analysis", "GSE132465_corrected"))
fig_dir <- ensure_dir(project_path("07_figures", "GSE132465_corrected", "cluster"))
seu <- readRDS(file.path(out_dir, "GSE132465_QC_filtered_seurat.rds"))
seu <- Seurat::SCTransform(seu, vars.to.regress = "percent.mt", conserve.memory = TRUE,
                          return.only.var.genes = TRUE, verbose = TRUE)
seu <- Seurat::RunPCA(seu, verbose = FALSE)
seu <- Seurat::RunUMAP(seu, dims = 1:30, seed.use = 1234)
seu <- Seurat::FindNeighbors(seu, dims = 1:30)
seu <- Seurat::FindClusters(seu, resolution = 0.5, random.seed = 1234)

p <- Seurat::DimPlot(seu, reduction = "umap", group.by = "Cell_type", label = TRUE, repel = TRUE)
ggplot2::ggsave(file.path(fig_dir, "UMAP_author_Cell_type.pdf"), p, width = 9, height = 7)
saveRDS(seu, file.path(out_dir, "GSE132465_SCT_clustered_seurat.rds"))
write_session_info(project_path("11_logs", "GSE132465_corrected", "3_sessionInfo.txt"))
