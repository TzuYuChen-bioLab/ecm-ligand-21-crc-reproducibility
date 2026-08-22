# 1.R — create Seurat object without inventing duplicate gene identifiers
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))
require_packages(c("data.table", "Matrix", "Seurat"))

raw_dir <- project_path("00_raw_data", "scRNA", "GSE132465")
out_dir <- ensure_dir(project_path("01_processed_data", "GSE132465_corrected"))
annotation_file <- file.path(raw_dir, "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")
count_file <- file.path(raw_dir, "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")
stopifnot(file.exists(annotation_file), file.exists(count_file))

annotation <- data.table::fread(annotation_file)
count_dt <- data.table::fread(count_file)
gene_column <- names(count_dt)[1]
genes <- as.character(count_dt[[gene_column]])
if (anyNA(genes) || any(!nzchar(genes))) stop("Count matrix contains missing/blank gene identifiers.")
cell_names <- names(count_dt)[-1]
if (anyDuplicated(cell_names)) stop("Count matrix contains duplicate cell-barcode columns.")

# Build sparse triplets directly from the data.table. This avoids materialising
# the entire genes × cells matrix as a dense R matrix. Mapping duplicate source
# rows to one index makes sparseMatrix sum their counts automatically.
gene_levels <- unique(genes)
gene_index <- match(genes, gene_levels)
i_list <- vector("list", length(cell_names))
x_list <- vector("list", length(cell_names))
for (j in seq_along(cell_names)) {
  values <- as.numeric(count_dt[[j + 1L]])
  if (any(!is.finite(values)) || any(values < 0)) stop("Non-finite or negative UMI count in cell column: ", cell_names[j])
  nonzero <- which(values != 0)
  i_list[[j]] <- gene_index[nonzero]
  x_list[[j]] <- values[nonzero]
}
nonzero_per_cell <- lengths(i_list)
count_matrix <- Matrix::sparseMatrix(
  i = unlist(i_list, use.names = FALSE),
  j = rep.int(seq_along(cell_names), nonzero_per_cell),
  x = unlist(x_list, use.names = FALSE),
  dims = c(length(gene_levels), length(cell_names)),
  dimnames = list(gene_levels, cell_names),
  giveCsparse = TRUE
)
rm(count_dt, i_list, x_list)
invisible(gc())

seu <- Seurat::CreateSeuratObject(
  counts = count_matrix, project = "GSE132465_CRC", min.cells = 3, min.features = 200
)
saveRDS(annotation, file.path(out_dir, "GSE132465_annotation.rds"))
saveRDS(seu, file.path(out_dir, "GSE132465_raw_seurat.rds"))
write_csv_atomic(
  data.frame(n_genes = nrow(seu), n_cells = ncol(seu), duplicate_gene_rows_aggregated = anyDuplicated(genes) > 0),
  file.path(out_dir, "create_seurat_summary.csv")
)
message("Saved corrected raw Seurat object: ", file.path(out_dir, "GSE132465_raw_seurat.rds"))
