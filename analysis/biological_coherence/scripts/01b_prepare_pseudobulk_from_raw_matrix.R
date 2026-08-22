#!/usr/bin/env Rscript

# Prepare patient-level pseudobulk raw counts from a gene-by-cell raw UMI
# matrix plus a cell-annotation table. Intended for deposits such as GSE144735.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) y else x

parse_cli <- function(x) {
  out <- list()
  for (item in x) {
    if (!grepl("^--[^=]+=", item)) next
    out[[sub("^--([^=]+)=.*$", "\\1", item)]] <- sub("^--[^=]+=", "", item)
  }
  out
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
required_args <- c("counts", "metadata", "dataset", "patient_col", "condition_col", "celltype_col", "outdir")
missing_args <- required_args[!vapply(required_args, function(z) nzchar(args[[z]] %||% ""), logical(1))]
if (length(missing_args)) {
  stop(
    "Missing argument(s): ", paste(missing_args, collapse = ", "), "\n",
    "Use explicit annotation column names. cell_id_col is optional and will be inferred by barcode overlap."
  )
}

needed <- c("data.table", "Matrix")
missing_packages <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install missing package(s): ", paste(missing_packages, collapse = ", "))
if (!file.exists(args$counts)) stop("Count matrix not found: ", args$counts)
if (!file.exists(args$metadata)) stop("Metadata table not found: ", args$metadata)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

minimum_cells <- as.integer(args$minimum_cells %||% "20")
if (!is.finite(minimum_cells) || minimum_cells < 1) stop("minimum_cells must be a positive integer.")
dataset <- gsub("[^A-Za-z0-9.-]+", "-", args$dataset)

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\.[0-9]+$", "", x)
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

clean_token <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  gsub("[^A-Za-z0-9.-]+", "-", x)
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

message("Reading raw gene-by-cell matrix; a 64-bit R session and substantial RAM may be required.")
count_table <- data.table::fread(args$counts, data.table = FALSE, check.names = FALSE)
if (ncol(count_table) < 2) stop("Count matrix needs one gene column plus cell columns.")
gene_col <- args$gene_col %||% names(count_table)[1]
if (!gene_col %in% names(count_table)) stop("gene_col not found: ", gene_col)
genes <- clean_gene(count_table[[gene_col]])
count_table[[gene_col]] <- NULL
cell_ids <- names(count_table)
if (anyDuplicated(cell_ids)) stop("Count matrix has duplicated cell columns.")
count_table[] <- lapply(count_table, function(z) suppressWarnings(as.numeric(z)))
counts <- as.matrix(count_table)
rm(count_table)
invisible(gc())
if (length(genes) != nrow(counts)) stop("Gene vector does not match count-matrix rows.")
if (anyNA(genes)) stop("Count matrix contains blank/missing gene identifiers.")
if (any(!is.finite(counts)) || any(counts < 0) || max(abs(counts - round(counts))) > 1e-6) {
  stop("Count matrix must contain finite non-negative integer UMI counts.")
}
rownames(counts) <- genes
colnames(counts) <- cell_ids
if (anyDuplicated(genes)) counts <- rowsum(counts, group = genes, reorder = FALSE)
counts <- Matrix::Matrix(counts, sparse = TRUE)

annotation <- data.table::fread(args$metadata, data.table = FALSE, check.names = FALSE)
required_columns <- c(args$patient_col, args$condition_col, args$celltype_col)
if (length(setdiff(required_columns, names(annotation)))) {
  stop(
    "Metadata column(s) not found: ", paste(setdiff(required_columns, names(annotation)), collapse = ", "),
    "\nAvailable columns: ", paste(names(annotation), collapse = ", ")
  )
}

if (!is.null(args$cell_id_col) && nzchar(args$cell_id_col)) {
  if (!args$cell_id_col %in% names(annotation)) stop("cell_id_col not found: ", args$cell_id_col)
  cell_id_col <- args$cell_id_col
} else {
  overlap <- vapply(annotation, function(z) sum(as.character(z) %in% cell_ids), integer(1))
  if (!length(overlap) || max(overlap) < 0.90 * length(cell_ids)) {
    stop("Could not infer a metadata cell-ID column with at least 90% barcode overlap. Supply --cell_id_col=...")
  }
  cell_id_col <- names(which.max(overlap))
}

annotation_cell_ids <- as.character(annotation[[cell_id_col]])
index <- match(cell_ids, annotation_cell_ids)
if (sum(!is.na(index)) < 0.90 * length(cell_ids)) stop("Fewer than 90% of count-matrix cells matched metadata.")
keep_matched <- !is.na(index)
counts <- counts[, keep_matched, drop = FALSE]
cell_ids <- cell_ids[keep_matched]
annotation <- annotation[index[keep_matched], , drop = FALSE]

cell_meta <- data.frame(
  cell_id = cell_ids,
  patient_id = clean_token(annotation[[args$patient_col]]),
  condition = standardize_condition(annotation[[args$condition_col]]),
  cell_class = standardize_cell_class(annotation[[args$celltype_col]]),
  stringsAsFactors = FALSE
)
valid <- !is.na(cell_meta$patient_id) & !is.na(cell_meta$condition) & cell_meta$cell_class != "Other"
counts <- counts[, valid, drop = FALSE]
cell_meta <- cell_meta[valid, , drop = FALSE]
if (!ncol(counts)) stop("No cells retained after metadata standardization.")

group_id <- paste(dataset, cell_meta$patient_id, cell_meta$condition, cell_meta$cell_class, sep = "__")
group_sizes <- table(group_id)
keep_groups <- names(group_sizes)[group_sizes >= minimum_cells]
keep_cells <- group_id %in% keep_groups
counts <- counts[, keep_cells, drop = FALSE]
cell_meta <- cell_meta[keep_cells, , drop = FALSE]
group_id <- group_id[keep_cells]
if (!length(group_id)) stop("No pseudobulk group has at least ", minimum_cells, " cells.")

group_factor <- factor(group_id, levels = unique(group_id))
aggregation_matrix <- Matrix::sparseMatrix(
  i = seq_along(group_factor), j = as.integer(group_factor), x = 1,
  dims = c(length(group_factor), nlevels(group_factor)),
  dimnames = list(NULL, levels(group_factor))
)
pseudobulk <- counts %*% aggregation_matrix
colnames(pseudobulk) <- levels(group_factor)

split_id <- strsplit(colnames(pseudobulk), "__", fixed = TRUE)
pb_fields <- do.call(rbind, lapply(split_id, function(z) {
  if (length(z) != 4) stop("Internal sample identifier could not be parsed.")
  z
}))
pb_meta <- data.frame(
  sample_id = colnames(pseudobulk),
  dataset = pb_fields[, 1],
  patient_id = pb_fields[, 2],
  condition = pb_fields[, 3],
  cell_class = pb_fields[, 4],
  n_cells = as.integer(table(group_id)[colnames(pseudobulk)]),
  paired_id = paste(pb_fields[, 1], pb_fields[, 2], sep = "__"),
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
  count_source = normalizePath(args$counts, winslash = "/", mustWork = TRUE),
  metadata_source = normalizePath(args$metadata, winslash = "/", mustWork = TRUE),
  cell_id_column = cell_id_col,
  patient_column = args$patient_col,
  condition_column = args$condition_col,
  celltype_column = args$celltype_col,
  minimum_cells = minimum_cells,
  retained_cells = sum(pb_meta$n_cells),
  pseudobulk_samples = ncol(pseudobulk),
  genes = nrow(pseudobulk),
  paired_fibroblast_patients = length(intersect(
    pb_meta$paired_id[pb_meta$cell_class == "Fibroblast" & pb_meta$condition == "Tumor"],
    pb_meta$paired_id[pb_meta$cell_class == "Fibroblast" & pb_meta$condition == "Normal"]
  )),
  stringsAsFactors = FALSE
)
utils::write.csv(audit, audit_path, row.names = FALSE, na = "")
message("Created:\n", counts_path, "\n", metadata_path, "\n", audit_path)
