# 0.R — preflight and input checksum inventory
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "config.R"))

raw_files <- c(
  GSE132465_annotation = project_path("00_raw_data", "scRNA", "GSE132465", "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"),
  GSE132465_counts = project_path("00_raw_data", "scRNA", "GSE132465", "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"),
  GSE39582 = project_path("00_raw_data", "bulk", "GSE39582", "GSE39582_series_matrix.txt.gz"),
  GSE38832 = project_path("00_raw_data", "bulk", "GSE38832", "GSE38832_series_matrix.txt.gz"),
  GSE14333 = project_path("00_raw_data", "bulk", "GSE14333", "GSE14333_series_matrix.txt.gz"),
  GSE17536 = project_path("00_raw_data", "bulk", "GSE17536", "GSE17536_series_matrix.txt.gz"),
  GSE17537 = project_path("00_raw_data", "bulk", "GSE17537", "GSE17537_series_matrix.txt.gz"),
  GSE33113 = project_path("00_raw_data", "bulk", "GSE33113", "GSE33113_series_matrix.txt.gz")
)

inventory <- data.frame(
  item = names(raw_files), path = unname(raw_files), exists = file.exists(raw_files),
  size_bytes = ifelse(file.exists(raw_files), file.info(raw_files)$size, NA_real_),
  md5 = ifelse(file.exists(raw_files), unname(tools::md5sum(raw_files)), NA_character_),
  stringsAsFactors = FALSE
)
out_dir <- ensure_dir(project_path("11_logs", "preflight"))
write_csv_atomic(inventory, file.path(out_dir, "input_file_inventory.csv"))

required_now <- c("GSE132465_annotation", "GSE132465_counts")
if (any(!inventory$exists[inventory$item %in% required_now])) {
  stop("Required GSE132465 input files are missing. See ", file.path(out_dir, "input_file_inventory.csv"))
}

message("Preflight complete. Missing bulk files are allowed until scripts 9.R/11.R are run.")
