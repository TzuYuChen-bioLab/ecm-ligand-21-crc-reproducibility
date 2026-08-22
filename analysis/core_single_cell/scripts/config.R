# CRC extracellular-matrix ligand programme: shared configuration
# Version: 2.1.0 (2026-07-16)

options(stringsAsFactors = FALSE)
options(timeout = max(7200, getOption("timeout")))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active)) return(dirname(normalizePath(active, mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- get_script_dir()

resolve_project_root <- function() {
  from_env <- Sys.getenv("CRC_PROJECT_ROOT", unset = "")
  if (nzchar(from_env)) return(normalizePath(from_env, mustWork = FALSE))

  candidates <- unique(c(
    normalizePath(getwd(), mustWork = FALSE),
    normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE),
    normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
  ))
  has_project_layout <- vapply(
    candidates,
    function(x) dir.exists(file.path(x, "00_raw_data")) || dir.exists(file.path(x, "01_processed_data")),
    logical(1)
  )
  if (any(has_project_layout)) return(candidates[which(has_project_layout)[1]])

  stop(
    "Cannot determine the project root. Set environment variable CRC_PROJECT_ROOT ",
    "to the directory containing 00_raw_data, then rerun."
  )
}

PROJECT_ROOT <- resolve_project_root()
set.seed(1234)

project_path <- function(...) file.path(PROJECT_ROOT, ...)

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}

require_packages <- function(packages, bioconductor = FALSE) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  installer <- if (bioconductor) {
    paste0("BiocManager::install(c(", paste(sprintf("'%s'", missing), collapse = ", "), "))")
  } else {
    paste0("install.packages(c(", paste(sprintf("'%s'", missing), collapse = ", "), "))")
  }
  stop(
    "Missing required package(s): ", paste(missing, collapse = ", "),
    ".\nInstall explicitly in the locked project environment, for example:\n", installer,
    "\nThe pipeline never installs or upgrades packages automatically."
  )
}

get_assay_data_compat <- function(object, assay = "RNA", layer = "data") {
  tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = layer)
  )
}

normalise_condition <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("tumou?r|cancer|crc", x)] <- "Tumor"
  out[grepl("normal|non[- ]?tumou?r|adjacent", x)] <- "Normal"
  out
}

clean_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

extract_numeric <- function(x, decimal_comma = TRUE) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "n/a", "null", "unknown")] <- NA_character_
  if (decimal_comma) x <- gsub("(?<=\\d),(?=\\d)", ".", x, perl = TRUE)
  suppressWarnings(as.numeric(sub(".*?(-?[0-9]+(?:\\.[0-9]+)?).*", "\\1", x)))
}

convert_event <- function(x) {
  raw <- trimws(tolower(as.character(x)))
  raw[raw %in% c("", "na", "n/a", "unknown", "null")] <- NA_character_
  out <- rep(NA_real_, length(raw))
  out[grepl("^(0)(\\b|\\s)|^no$|^false$|no recurrence|no relapse|no death|alive|censor|disease free", raw)] <- 0
  out[grepl("^(1)(\\b|\\s)|^yes$|^true$|recurrence|relapse|death|dead|deceased|progression", raw) & is.na(out)] <- 1
  num <- suppressWarnings(as.numeric(raw))
  out[is.na(out) & num %in% c(0, 1)] <- num[is.na(out) & num %in% c(0, 1)]
  out
}

safe_z <- function(x, centre = NULL, scale_value = NULL) {
  x <- as.numeric(x)
  if (is.null(centre)) centre <- mean(x, na.rm = TRUE)
  if (is.null(scale_value)) scale_value <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(scale_value) || scale_value <= 0) {
    return(list(values = rep(NA_real_, length(x)), centre = centre, scale = scale_value))
  }
  list(values = (x - centre) / scale_value, centre = centre, scale = scale_value)
}

write_csv_atomic <- function(x, path, row.names = FALSE) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  utils::write.csv(x, tmp, row.names = row.names, na = "")
  if (!file.rename(tmp, path)) stop("Failed to publish output: ", path)
  invisible(path)
}

write_session_info <- function(path) {
  ensure_dir(dirname(path))
  writeLines(capture.output(sessionInfo()), path)
}

assert_no_duplicate_samples <- function(sample_ids, label) {
  dup <- unique(sample_ids[duplicated(sample_ids)])
  if (length(dup) > 0) stop(label, " contains duplicated sample identifiers: ", paste(head(dup, 10), collapse = ", "))
  invisible(TRUE)
}

message("CRC corrected pipeline project root: ", PROJECT_ROOT)
