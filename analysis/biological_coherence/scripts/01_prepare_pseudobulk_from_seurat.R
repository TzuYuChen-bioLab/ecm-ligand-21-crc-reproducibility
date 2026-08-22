#!/usr/bin/env Rscript

# Prepare patient-level pseudobulk raw counts from an annotated Seurat object.
# The script does not use survival outcomes and does not select panel genes.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) y else x

parse_cli <- function(x) {
  out <- list()
  for (item in x) {
    if (!grepl("^--[^=]+=", item)) next
    key <- sub("^--([^=]+)=.*$", "\\1", item)
    value <- sub("^--[^=]+=", "", item)
    out[[key]] <- value
  }
  out
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
required_args <- c("rds", "dataset", "patient_col", "condition_col", "celltype_col", "outdir")
missing_args <- required_args[!vapply(required_args, function(z) nzchar(args[[z]] %||% ""), logical(1))]

if (length(missing_args)) {
  stop(
    "Missing argument(s): ", paste(missing_args, collapse = ", "), "\n",
    "Example:\n",
    "Rscript 01_prepare_pseudobulk_from_seurat.R \\\n",
    "  --rds=C:/CRC/GSE132465_author_annotated_corrected.rds \\\n",
    "  --dataset=GSE132465 --patient_col=Patient --condition_col=condition \\\n",
    "  --celltype_col=major_celltype --outdir=C:/CRC/provenance_input"
  )
}

needed <- c("Matrix")
missing_packages <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install missing package(s): ", paste(missing_packages, collapse = ", "))

if (!file.exists(args$rds)) stop("RDS file not found: ", args$rds)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

dataset <- gsub("[^A-Za-z0-9.-]+", "-", args$dataset)
assay <- args$assay %||% "RNA"
minimum_cells <- as.integer(args$minimum_cells %||% "20")
if (!is.finite(minimum_cells) || minimum_cells < 1) stop("minimum_cells must be a positive integer.")

clean_token <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x <- gsub("[^A-Za-z0-9.-]+", "-", x)
  x
}

standardize_condition <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(raw))
  out[grepl("tumou?r|cancer|carcinoma|(^|[-_ ])t($|[-_ ])", raw)] <- "Tumor"
  out[grepl("border|margin|interface|(^|[-_ ])b($|[-_ ])", raw)] <- "Border"
  out[grepl("normal|adjacent|non[-_ ]?tumou?r|(^|[-_ ])n($|[-_ ])", raw)] <- "Normal"
  out
}

standardize_cell_class <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  out <- rep("Other", length(raw))
  out[grepl("pericyte", raw)] <- "Pericyte"
  out[grepl("smooth.?muscle", raw)] <- "Smooth_muscle"
  out[grepl("endothel|vascular", raw)] <- "Endothelial"
  out[grepl("epithel|malignan|tumou?r.?cell|cancer.?cell|colonocyte|goblet|enterocyte", raw)] <- "Epithelial"
  out[grepl("myelo|macro|monocyte|dendritic|neutroph", raw)] <- "Myeloid"
  out[grepl("t.?cell|nk.?cell|lymphoid", raw)] <- "T_NK"
  out[grepl("b.?cell|plasma", raw)] <- "B_Plasma"
  out[grepl("mast", raw)] <- "Mast"
  out[grepl("fibro|myofibro|stromal|(^|[^a-z])caf([^a-z]|$)", raw)] <- "Fibroblast"
  out
}

get_counts <- function(object, assay_name) {
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    ans <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay_name, layer = "counts"),
      error = function(e) tryCatch(
        SeuratObject::GetAssayData(object, assay = assay_name, slot = "counts"),
        error = function(e2) NULL
      )
    )
    if (!is.null(ans)) return(ans)
  }
  if (requireNamespace("Seurat", quietly = TRUE)) {
    ans <- tryCatch(
      Seurat::GetAssayData(object, assay = assay_name, layer = "counts"),
      error = function(e) tryCatch(
        Seurat::GetAssayData(object, assay = assay_name, slot = "counts"),
        error = function(e2) NULL
      )
    )
    if (!is.null(ans)) return(ans)
  }
  stop("Cannot obtain the raw-count layer. Install SeuratObject/Seurat and confirm that assay '", assay_name, "' contains counts.")
}

object <- readRDS(args$rds)
metadata <- tryCatch(as.data.frame(object[[]]), error = function(e) NULL)
if (is.null(metadata) || !nrow(metadata)) {
  metadata <- tryCatch(as.data.frame(object@meta.data), error = function(e) NULL)
}
if (is.null(metadata) || !nrow(metadata)) stop("The RDS does not expose Seurat cell metadata.")

required_columns <- c(args$patient_col, args$condition_col, args$celltype_col)
if (length(setdiff(required_columns, names(metadata)))) {
  stop(
    "Metadata column(s) not found: ", paste(setdiff(required_columns, names(metadata)), collapse = ", "),
    "\nAvailable columns: ", paste(names(metadata), collapse = ", ")
  )
}

counts <- get_counts(object, assay)
if (!inherits(counts, "Matrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)
if (is.null(rownames(counts)) || is.null(colnames(counts))) stop("Count matrix needs gene and cell names.")
source_cell_count <- ncol(counts)
count_values <- if ("x" %in% methods::slotNames(counts)) counts@x else as.numeric(counts)
if (length(count_values) && (any(!is.finite(count_values)) || any(count_values < 0))) stop("Counts contain non-finite or negative values.")
if (length(count_values) && max(abs(count_values - round(count_values))) > 1e-6) {
  stop("The selected assay layer is not raw integer counts. Use the RNA counts layer, not normalized/SCT data.")
}

cell_index <- match(colnames(counts), rownames(metadata))
if (anyNA(cell_index)) stop("Count-matrix cells do not all match metadata row names.")
metadata <- metadata[cell_index, , drop = FALSE]

cell_meta <- data.frame(
  cell_id = colnames(counts),
  patient_id = clean_token(metadata[[args$patient_col]]),
  condition = standardize_condition(metadata[[args$condition_col]]),
  cell_class = standardize_cell_class(metadata[[args$celltype_col]]),
  stringsAsFactors = FALSE
)

valid <- !is.na(cell_meta$patient_id) & !is.na(cell_meta$condition) & cell_meta$cell_class != "Other"
if (!any(valid)) stop("No cells retained after patient/condition/cell-class standardization.")

counts <- counts[, valid, drop = FALSE]
cell_meta <- cell_meta[valid, , drop = FALSE]
group_id <- paste(dataset, cell_meta$patient_id, cell_meta$condition, cell_meta$cell_class, sep = "__")
group_sizes <- table(group_id)
keep_groups <- names(group_sizes)[group_sizes >= minimum_cells]
keep_cells <- group_id %in% keep_groups
if (!any(keep_cells)) stop("No pseudobulk group has at least ", minimum_cells, " cells.")

counts <- counts[, keep_cells, drop = FALSE]
cell_meta <- cell_meta[keep_cells, , drop = FALSE]
group_id <- group_id[keep_cells]
group_factor <- factor(group_id, levels = unique(group_id))
aggregation_matrix <- Matrix::sparseMatrix(
  i = seq_along(group_factor),
  j = as.integer(group_factor),
  x = 1,
  dims = c(length(group_factor), nlevels(group_factor)),
  dimnames = list(NULL, levels(group_factor))
)
pseudobulk <- counts %*% aggregation_matrix
colnames(pseudobulk) <- levels(group_factor)

split_id <- strsplit(colnames(pseudobulk), "__", fixed = TRUE)
pb_meta <- do.call(rbind, lapply(split_id, function(z) {
  if (length(z) != 4) stop("Internal sample identifier could not be parsed: ", paste(z, collapse = "__"))
  z
}))
pb_meta <- data.frame(
  sample_id = colnames(pseudobulk),
  dataset = pb_meta[, 1],
  patient_id = pb_meta[, 2],
  condition = pb_meta[, 3],
  cell_class = pb_meta[, 4],
  n_cells = as.integer(table(group_id)[colnames(pseudobulk)]),
  paired_id = paste(pb_meta[, 1], pb_meta[, 2], sep = "__"),
  stringsAsFactors = FALSE
)

counts_path <- file.path(args$outdir, paste0(dataset, "_pseudobulk_counts.csv.gz"))
metadata_path <- file.path(args$outdir, paste0(dataset, "_pseudobulk_metadata.csv"))
audit_path <- file.path(args$outdir, paste0(dataset, "_pseudobulk_preparation_audit.csv"))

counts_df <- data.frame(gene = rownames(pseudobulk), as.matrix(pseudobulk), check.names = FALSE)
connection <- gzfile(counts_path, open = "wt")
utils::write.table(counts_df, connection, sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
close(connection)
utils::write.csv(pb_meta, metadata_path, row.names = FALSE, na = "")

audit <- data.frame(
  dataset = dataset,
  source_rds = normalizePath(args$rds, winslash = "/", mustWork = TRUE),
  assay = assay,
  minimum_cells = minimum_cells,
  source_cells = source_cell_count,
  retained_cells = sum(pb_meta$n_cells),
  pseudobulk_samples = ncol(pseudobulk),
  genes = nrow(pseudobulk),
  fibroblast_samples = sum(pb_meta$cell_class == "Fibroblast"),
  paired_fibroblast_patients = length(intersect(
    pb_meta$paired_id[pb_meta$cell_class == "Fibroblast" & pb_meta$condition == "Tumor"],
    pb_meta$paired_id[pb_meta$cell_class == "Fibroblast" & pb_meta$condition == "Normal"]
  )),
  stringsAsFactors = FALSE
)
utils::write.csv(audit, audit_path, row.names = FALSE, na = "")

message("Created:\n", counts_path, "\n", metadata_path, "\n", audit_path)
