# =============================================================================
# GSE39582 ECM21 module Cox analysis and paired bootstrap delta C-index
# Standalone repair/continuation script for the ECM21 module analysis
#
# Author-facing purpose
#   1. Reuse the official GSE39582 processed series matrix (RMA + ComBat).
#   2. Learn probe choice and gene standardisation in Discovery tumours only.
#   3. Score M1-M4, the gene-equal 21-gene score, the module-balanced score,
#      the 19-gene structural core, and leave-one-module-out scores.
#   4. Fit continuous-score RFS Cox models per Discovery-reference SD.
#   5. Compare paired Harrell C-indices by patient-level bootstrap.
#   6. Apply a conservative claim gate: a positive gain, bootstrap CI support,
#      and stable gain direction in Discovery and Validation are all required
#      before the word "enhanced" is permitted.
#
# Recommended location
#   PROJECT_ROOT
#
# Run
#   Open this file in RStudio and click Source, or run:
#   Rscript RUN_GSE39582_MODULE_COX_AND_BOOTSTRAP_DELTA_CINDEX.R
#
# Required input
#   This run is fixed to the author's intact local GSE39582 series matrix:
#   PROJECT_ROOT
#   GSE39582_series_matrix.txt.gz
#   The file is validated before parsing. No GEO download is attempted.
#
# Optional CAF input (needed for CAF-adjusted rows)
#   A) a sample-level table containing a sample ID and CAF_PROXY_30 score, or
#   B) a text/CSV file containing the exact frozen CAF_PROXY_30 gene symbols.
#   Set CAF_SCORE_FILE or CAF_GENE_FILE explicitly below. The script never
#   invents or substitutes a CAF signature. Clinical-only analyses still run
#   if neither CAF resource is available, and the omission is recorded.
# =============================================================================

# ------------------------------- USER SETTINGS -------------------------------
USER_ROOT <- Sys.getenv("CBC_REPRO_ROOT", unset = "")

# Fixed local series-matrix path supplied by the author.
SERIES_MATRIX_FILE <- ""
CAF_SCORE_FILE <- ""
CAF_GENE_FILE <- ""

AUTO_INSTALL_PACKAGES <- FALSE
ALLOW_OFFICIAL_GEO_DOWNLOAD <- FALSE
REQUIRE_CAF_INPUT <- FALSE

# The original master script used 500. Use 1000 for the final manuscript run;
# 2000 is reasonable if runtime permits.
CINDEX_BOOTSTRAP_REPS <- 1000L
BOOTSTRAP_SEED <- 20260812L
MIN_BOOTSTRAP_SUCCESS_FRACTION <- 0.80

# Ordinary clinical model and the prespecified sparse-level sensitivity model.
# CAF-adjusted variants are added automatically when a valid CAF resource exists.
BOOTSTRAP_MODEL_KEYS <- c("clinical_primary", "clinical_stratified", "caf_primary", "caf_stratified")

# An intact 585-sample series matrix is much larger than the previous truncated
# 0.96 MB cache. This threshold is only an early corruption screen; GEOquery
# parsing and sample-count checks provide the definitive validation.
MIN_SERIES_MATRIX_GZIP_BYTES <- 10 * 1024^2

# ----------------------------- FIXED GENE MODULES -----------------------------
ECM_MODULES <- list(
  M1 = c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "COL6A3"),
  M2 = c("COL4A1", "COL4A2", "COL4A5", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMC1", "HSPG2"),
  M3 = c("FN1", "THBS1", "THBS2", "TNC", "TNXB"),
  M4 = c("PTN", "MDK")
)
ECM21_GENES <- unique(unlist(ECM_MODULES, use.names = FALSE))
ECM19_GENES <- unique(unlist(ECM_MODULES[c("M1", "M2", "M3")], use.names = FALSE))
stopifnot(length(ECM21_GENES) == 21L, length(ECM19_GENES) == 19L)

# ------------------------------- BASIC HELPERS --------------------------------
options(stringsAsFactors = FALSE)
options(timeout = max(7200, getOption("timeout", 60)))

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) {
    return(dirname(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE)))
  }
  frames <- sys.frames()
  ofiles <- vapply(frames, function(x) if (!is.null(x$ofile)) as.character(x$ofile)[1] else "", character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles)) return(dirname(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = FALSE)))
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

SCRIPT_DIR <- script_dir()
PACKAGE_ROOT <- normalizePath(file.path(SCRIPT_DIR, '..', '..', '..'), winslash = '/', mustWork = TRUE)
if (!nzchar(USER_ROOT)) USER_ROOT <- file.path(PACKAGE_ROOT, 'results', 'generated')
SERIES_MATRIX_FILE <- file.path(PACKAGE_ROOT, 'data', 'available', 'GSE39582', 'GSE39582_series_matrix.txt.gz')
dir.create(USER_ROOT, recursive = TRUE, showWarnings = FALSE)
if (FALSE) {
  warning("Configured USER_ROOT does not exist; using the script directory: ", SCRIPT_DIR)
  USER_ROOT <- SCRIPT_DIR
}
USER_ROOT <- normalizePath(USER_ROOT, winslash = "/", mustWork = FALSE)
RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(USER_ROOT, paste0("GSE39582_ECM21_MODULE_COX_BOOTSTRAP_RESULTS_", RUN_ID))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
CACHE_DIR <- file.path(PACKAGE_ROOT, "data", "available", "GSE39582")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUTPUT_DIR, "MASTER_RUN_LOG.txt")

log_line <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOG_FILE, append = TRUE)
}

write_csv <- function(x, filename) {
  utils::write.csv(x, filename, row.names = FALSE, na = "")
  invisible(filename)
}

normalise_token <- function(x) {
  x <- trimws(as.character(x))
  x[tolower(x) %in% c("", "na", "n/a", "nan", "null", "not available", "not determined", "nd", "-")] <- NA_character_
  x
}

safe_numeric <- function(x) {
  x <- normalise_token(x)
  x <- gsub(",", ".", x, fixed = TRUE)
  out <- suppressWarnings(as.numeric(sub(".*?(-?[0-9]+(?:\\.[0-9]+)?).*", "\\1", x, perl = TRUE)))
  out[is.na(x)] <- NA_real_
  out
}

normalise_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

install_and_load <- function() {
  cran <- c("survival", "curl")
  bioc <- c("GEOquery", "Biobase", "AnnotationDbi", "hgu133plus2.db")
  missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if ((length(missing_cran) || length(missing_bioc)) && !isTRUE(AUTO_INSTALL_PACKAGES)) {
    stop("Missing packages: ", paste(c(missing_cran, missing_bioc), collapse = ", "))
  }
  if (length(missing_cran)) {
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  }
  if (length(missing_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
  still_missing <- c(cran, bioc)[!vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing)) stop("Package installation failed: ", paste(still_missing, collapse = ", "))
}

cat(c(
  "GSE39582 ECM21 module Cox and bootstrap delta C-index run",
  paste0("Run ID: ", RUN_ID),
  paste0("User root: ", USER_ROOT),
  paste0("Output directory: ", OUTPUT_DIR),
  paste0("R version: ", R.version.string),
  ""
), sep = "\n", file = LOG_FILE)

install_and_load()

# -------------------------- SERIES MATRIX ACQUISITION -------------------------
looks_like_intact_gzip <- function(path, minimum_bytes = MIN_SERIES_MATRIX_GZIP_BYTES) {
  if (!nzchar(path) || !file.exists(path)) return(FALSE)
  info <- file.info(path)
  if (!is.finite(info$size) || info$size < minimum_bytes) return(FALSE)
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  first <- tryCatch(readLines(con, n = 8L, warn = FALSE), error = function(e) character())
  length(first) > 0L && any(grepl("Series_", first, fixed = TRUE))
}

find_local_series_matrix <- function() {
  explicit <- normalise_token(SERIES_MATRIX_FILE)
  if (!is.na(explicit) && nzchar(explicit)) {
    explicit <- normalizePath(explicit, winslash = "/", mustWork = FALSE)
    if (!looks_like_intact_gzip(explicit)) {
      stop("SERIES_MATRIX_FILE is missing, too small, or not a readable GEO series matrix: ", explicit)
    }
    return(explicit)
  }
  candidates <- list.files(
    USER_ROOT,
    pattern = "GSE39582.*series[_-]?matrix.*\\.txt\\.gz$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(candidates)) return(NA_character_)
  valid <- candidates[vapply(candidates, looks_like_intact_gzip, logical(1))]
  if (!length(valid)) return(NA_character_)
  valid[which.max(file.info(valid)$size)]
}

download_official_series_matrix <- function() {
  if (!isTRUE(ALLOW_OFFICIAL_GEO_DOWNLOAD)) return(NA_character_)
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39582/matrix/GSE39582_series_matrix.txt.gz"
  target <- file.path(CACHE_DIR, "GSE39582_series_matrix.txt.gz")
  if (looks_like_intact_gzip(target)) return(target)
  if (file.exists(target)) unlink(target)
  part <- paste0(target, ".part")
  for (attempt in seq_len(4L)) {
    if (file.exists(part)) unlink(part)
    log_line("Official GEO download attempt ", attempt, "/4")
    ok <- tryCatch({
      handle <- curl::new_handle(
        connecttimeout = 120,
        timeout = 7200,
        low_speed_limit = 1024,
        low_speed_time = 180,
        followlocation = TRUE
      )
      curl::curl_download(url, part, quiet = FALSE, mode = "wb", handle = handle)
      file.rename(part, target)
      looks_like_intact_gzip(target)
    }, error = function(e) {
      log_line("Download attempt failed: ", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) return(target)
    if (file.exists(target) && !looks_like_intact_gzip(target)) unlink(target)
    if (attempt < 4L) Sys.sleep(5 * attempt)
  }
  NA_character_
}

matrix_file <- find_local_series_matrix()
if (is.na(matrix_file)) matrix_file <- download_official_series_matrix()
if (is.na(matrix_file) || !looks_like_intact_gzip(matrix_file)) {
  stop(
    "No intact GSE39582 series matrix was available. Manually download the official file from:\n",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39582/matrix/GSE39582_series_matrix.txt.gz\n",
    "Place it under ", USER_ROOT, " and rerun. Do not reuse the previous 0.96 MB truncated cache."
  )
}
matrix_file <- normalizePath(matrix_file, winslash = "/", mustWork = TRUE)
log_line("Using series matrix: ", matrix_file, " (", round(file.info(matrix_file)$size / 1024^2, 1), " MB)")

# --------------------------- IMPORT AND PHENOTYPE -----------------------------
log_line("Parsing GSE39582 series matrix with GEOquery")
gse_raw <- GEOquery::getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE)
if (inherits(gse_raw, "ExpressionSet")) {
  eset <- gse_raw
} else if (is.list(gse_raw) && length(gse_raw)) {
  sizes <- vapply(gse_raw, function(x) if (inherits(x, "ExpressionSet")) ncol(Biobase::exprs(x)) else 0L, integer(1))
  if (!any(sizes > 0L)) stop("GEOquery returned no ExpressionSet.")
  eset <- gse_raw[[which.max(sizes)]]
} else {
  stop("Unexpected GEOquery return type: ", paste(class(gse_raw), collapse = "/"))
}
expr_probe <- Biobase::exprs(eset)
pd <- Biobase::pData(eset)
if (ncol(expr_probe) != nrow(pd)) stop("Expression and phenotype sample counts do not align.")
if (is.null(colnames(expr_probe)) || is.null(rownames(pd))) stop("Sample identifiers are missing.")
pd <- pd[colnames(expr_probe), , drop = FALSE]

extract_geo_field <- function(pd, aliases) {
  aliases <- tolower(aliases)
  norm_names <- tolower(gsub("[^a-z0-9]", "", names(pd)))
  norm_alias <- gsub("[^a-z0-9]", "", aliases)
  candidate_cols <- unique(unlist(lapply(norm_alias, function(a) which(grepl(a, norm_names, fixed = TRUE)))))
  strip_prefix <- function(x) sub("^[^:]+:\\s*", "", as.character(x), perl = TRUE)
  if (length(candidate_cols)) {
    for (j in candidate_cols) {
      value <- normalise_token(strip_prefix(pd[[j]]))
      if (sum(!is.na(value)) >= max(3L, floor(0.25 * nrow(pd)))) return(value)
    }
  }
  out <- rep(NA_character_, nrow(pd))
  alias_pattern <- paste(gsub("([.\\^$*+?(){}|\\[\\]\\\\])", "\\\\\\1", aliases), collapse = "|")
  rx <- paste0("^\\s*(?:", alias_pattern, ")\\s*:")
  for (i in seq_len(nrow(pd))) {
    vals <- as.character(pd[i, , drop = TRUE])
    hit <- which(grepl(rx, vals, ignore.case = TRUE, perl = TRUE))[1]
    if (length(hit) && !is.na(hit)) out[i] <- sub(rx, "", vals[hit], ignore.case = TRUE, perl = TRUE)
  }
  normalise_token(out)
}

parse_event <- function(x) {
  x0 <- tolower(normalise_token(x))
  out <- rep(NA_integer_, length(x0))
  out[x0 %in% c("1", "yes", "y", "event", "recurred", "recurrence")] <- 1L
  out[x0 %in% c("0", "no", "n", "censored", "no event", "non-recurrence")] <- 0L
  unresolved <- is.na(out) & !is.na(x0)
  out[unresolved] <- suppressWarnings(as.integer(x0[unresolved]))
  out[!out %in% c(0L, 1L)] <- NA_integer_
  out
}

parse_stage <- function(x) {
  x <- toupper(normalise_token(x))
  x <- gsub("TNM|STAGE|AJCC|[^IVX0-9]", "", x)
  roman <- c(I = "1", II = "2", III = "3", IV = "4")
  x[x %in% names(roman)] <- roman[x[x %in% names(roman)]]
  x <- sub("^([1-4]).*$", "\\1", x)
  x[!x %in% as.character(1:4)] <- NA_character_
  factor(x, levels = as.character(1:4))
}

dataset_raw <- extract_geo_field(pd, c("dataset"))
sex_raw <- extract_geo_field(pd, c("sex"))
age_raw <- extract_geo_field(pd, c("age.at.diagnosis", "age at diagnosis"))
stage_raw <- extract_geo_field(pd, c("tnm.stage", "tnm stage"))
rfs_event_raw <- extract_geo_field(pd, c("rfs.event", "rfs event"))
rfs_time_raw <- extract_geo_field(pd, c("rfs.delay", "rfs delay"))
mmr_raw <- extract_geo_field(pd, c("mmr.status", "mmr status"))
source_raw <- extract_geo_field(pd, c("source_name", "source name"))
title_raw <- if ("title" %in% names(pd)) as.character(pd$title) else rownames(pd)

dataset <- tolower(dataset_raw)
dataset[grepl("discov", dataset)] <- "Discovery"
dataset[grepl("valid", dataset)] <- "Validation"
dataset[!dataset %in% c("Discovery", "Validation")] <- NA_character_
sex <- toupper(substr(normalise_token(sex_raw), 1, 1))
sex[!sex %in% c("F", "M")] <- NA_character_
sex <- factor(sex, levels = c("F", "M"))
mmr <- tolower(normalise_token(mmr_raw))
mmr[grepl("dmmr|deficient|msi", mmr)] <- "dMMR"
mmr[grepl("pmmr|proficient|mss", mmr)] <- "pMMR"
mmr[!mmr %in% c("dMMR", "pMMR")] <- NA_character_
mmr <- factor(mmr, levels = c("pMMR", "dMMR"))
source_text <- tolower(paste(source_raw, title_raw))
is_tumour <- !grepl("non[- ]?tum|normal|mucosa", source_text, perl = TRUE)

clinical <- data.frame(
  geo_accession = rownames(pd),
  sample_title = title_raw,
  dataset = factor(dataset, levels = c("Discovery", "Validation")),
  is_tumour = is_tumour,
  age = safe_numeric(age_raw),
  sex = sex,
  stage = parse_stage(stage_raw),
  mmr = mmr,
  rfs_time_months = safe_numeric(rfs_time_raw),
  rfs_event = parse_event(rfs_event_raw),
  stringsAsFactors = FALSE,
  row.names = rownames(pd)
)

audit <- data.frame(
  metric = c(
    "series_matrix_bytes", "expression_probes", "all_samples", "tumour_samples",
    "non_tumour_samples", "discovery_tumours", "validation_tumours",
    "RFS_evaluable_tumours", "RFS_events", "missing_age", "missing_sex",
    "missing_stage", "missing_MMR"
  ),
  value = c(
    file.info(matrix_file)$size, nrow(expr_probe), ncol(expr_probe), sum(clinical$is_tumour),
    sum(!clinical$is_tumour), sum(clinical$is_tumour & clinical$dataset == "Discovery", na.rm = TRUE),
    sum(clinical$is_tumour & clinical$dataset == "Validation", na.rm = TRUE),
    sum(clinical$is_tumour & is.finite(clinical$rfs_time_months) & clinical$rfs_time_months > 0 & !is.na(clinical$rfs_event)),
    sum(clinical$rfs_event[clinical$is_tumour & is.finite(clinical$rfs_time_months) & clinical$rfs_time_months > 0], na.rm = TRUE),
    sum(clinical$is_tumour & is.na(clinical$age)), sum(clinical$is_tumour & is.na(clinical$sex)),
    sum(clinical$is_tumour & is.na(clinical$stage)), sum(clinical$is_tumour & is.na(clinical$mmr))
  ),
  stringsAsFactors = FALSE
)
write_csv(audit, file.path(OUTPUT_DIR, "tables", "GSE39582_input_audit.csv"))

if (ncol(expr_probe) < 580L || nrow(expr_probe) < 50000L) {
  stop("The parsed object is unexpectedly small (", nrow(expr_probe), " probes x ", ncol(expr_probe), " samples).")
}
if (sum(clinical$is_tumour) < 560L) stop("Tumour detection failed; fewer than 560 tumours were identified.")
if (sum(clinical$is_tumour & clinical$dataset == "Discovery", na.rm = TRUE) < 440L) stop("Discovery labels were not parsed correctly.")
if (sum(clinical$is_tumour & clinical$dataset == "Validation", na.rm = TRUE) < 120L) stop("Validation labels were not parsed correctly.")
log_line("Input audit passed: ", sum(clinical$is_tumour), " tumours; ", sum(!clinical$is_tumour), " non-tumour samples")

# ------------------------------ OPTIONAL CAF INPUT ----------------------------
read_delimited_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(readRDS(path))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  x <- tryCatch(utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || ncol(x) <= 1L) {
    x <- tryCatch(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  }
  x
}

find_optional_file <- function(explicit, pattern) {
  explicit <- normalise_token(explicit)
  if (!is.na(explicit) && nzchar(explicit)) {
    if (!file.exists(explicit)) stop("Configured optional input does not exist: ", explicit)
    return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  }
  hits <- list.files(USER_ROOT, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  hits <- hits[!grepl("GSE39582_ECM21_MODULE_COX_BOOTSTRAP_RESULTS_", hits, fixed = TRUE)]
  if (!length(hits)) return(NA_character_)
  normalizePath(hits[which.max(file.info(hits)$mtime)], winslash = "/", mustWork = TRUE)
}

caf_score_path <- find_optional_file(CAF_SCORE_FILE, "CAF.*(PROXY.*SCORE|SCORE.*PROXY|SAMPLE.*SCORE|SCORE.*SAMPLE).*(csv|tsv|txt|rds)$")
caf_gene_path <- find_optional_file(CAF_GENE_FILE, "CAF.*(GENE|SIGNATURE).*(csv|tsv|txt)$")

read_caf_genes <- function(path) {
  if (is.na(path)) return(character())
  obj <- read_delimited_auto(path)
  if (is.data.frame(obj)) {
    norm <- tolower(gsub("[^a-z0-9]", "", names(obj)))
    gene_col <- which(grepl("gene|symbol", norm))[1]
    vals <- if (!is.na(gene_col)) obj[[gene_col]] else unlist(obj, use.names = FALSE)
  } else {
    vals <- obj
  }
  vals <- toupper(trimws(as.character(vals)))
  vals <- vals[grepl("^[A-Z][A-Z0-9.-]*$", vals)]
  vals <- setdiff(unique(vals), c("GENE", "GENES", "SYMBOL", "GENE_SYMBOL", "CAF_PROXY_30"))
  vals
}

caf_genes <- read_caf_genes(caf_gene_path)
if (length(caf_genes) && length(intersect(caf_genes, ECM21_GENES))) {
  stop("CAF gene input overlaps ECM21: ", paste(intersect(caf_genes, ECM21_GENES), collapse = ", "))
}
if (length(caf_genes) && length(caf_genes) != 30L) {
  stop("CAF_GENE_FILE must contain the exact frozen 30-gene set; detected ", length(caf_genes), " unique symbols.")
}

# ------------------------- PROBE MAPPING AND SCORING --------------------------
genes_needed <- unique(c(ECM21_GENES, caf_genes))
log_line("Mapping GPL570 probes and selecting one probe per gene in Discovery only")
probe_map <- AnnotationDbi::select(
  hgu133plus2.db::hgu133plus2.db,
  keys = rownames(expr_probe),
  keytype = "PROBEID",
  columns = "SYMBOL"
)
probe_map <- probe_map[!is.na(probe_map$SYMBOL) & nzchar(probe_map$SYMBOL), c("PROBEID", "SYMBOL")]
probe_map$SYMBOL <- toupper(probe_map$SYMBOL)
symbol_count <- aggregate(SYMBOL ~ PROBEID, probe_map, function(x) length(unique(x)))
unambiguous_ids <- symbol_count$PROBEID[symbol_count$SYMBOL == 1L]
probe_map <- probe_map[probe_map$PROBEID %in% unambiguous_ids, , drop = FALSE]
probe_map <- probe_map[!duplicated(probe_map$PROBEID), , drop = FALSE]
probe_map <- probe_map[probe_map$SYMBOL %in% genes_needed, , drop = FALSE]

disc_ids_all <- rownames(clinical)[clinical$is_tumour & clinical$dataset == "Discovery"]
disc_ids_all <- intersect(disc_ids_all, colnames(expr_probe))
if (!length(disc_ids_all)) stop("No Discovery tumour samples aligned to the expression matrix.")

probe_iqr <- apply(expr_probe[probe_map$PROBEID, disc_ids_all, drop = FALSE], 1, stats::IQR, na.rm = TRUE)
probe_map$discovery_IQR <- probe_iqr[probe_map$PROBEID]
probe_map <- probe_map[order(probe_map$SYMBOL, -probe_map$discovery_IQR, probe_map$PROBEID), , drop = FALSE]
selected_map <- probe_map[!duplicated(probe_map$SYMBOL), , drop = FALSE]
missing_ecm <- setdiff(ECM21_GENES, selected_map$SYMBOL)
if (length(missing_ecm)) stop("Missing unambiguous GPL570 mapping for ECM21 genes: ", paste(missing_ecm, collapse = ", "))
if (length(caf_genes)) {
  missing_caf <- setdiff(caf_genes, selected_map$SYMBOL)
  if (length(missing_caf)) stop("Missing unambiguous GPL570 mapping for CAF genes: ", paste(missing_caf, collapse = ", "))
}
write_csv(selected_map, file.path(OUTPUT_DIR, "tables", "selected_probe_per_gene_discovery_reference.csv"))

gene_expr <- expr_probe[selected_map$PROBEID, , drop = FALSE]
rownames(gene_expr) <- selected_map$SYMBOL
gene_expr <- gene_expr[genes_needed, , drop = FALSE]
if (stats::quantile(gene_expr, 0.99, na.rm = TRUE) > 100) gene_expr <- log2(gene_expr + 1)

disc_mean <- rowMeans(gene_expr[, disc_ids_all, drop = FALSE], na.rm = TRUE)
disc_sd <- apply(gene_expr[, disc_ids_all, drop = FALSE], 1, stats::sd, na.rm = TRUE)
if (any(!is.finite(disc_sd) | disc_sd <= 0)) {
  stop("Non-positive Discovery SD for: ", paste(names(disc_sd)[!is.finite(disc_sd) | disc_sd <= 0], collapse = ", "))
}
gene_z <- sweep(sweep(gene_expr, 1, disc_mean, "-"), 1, disc_sd, "/")
standardisation <- data.frame(gene = names(disc_mean), discovery_mean = disc_mean, discovery_sd = disc_sd, stringsAsFactors = FALSE)
write_csv(standardisation, file.path(OUTPUT_DIR, "tables", "gene_standardisation_discovery_reference.csv"))

score_mean <- function(genes) colMeans(gene_z[genes, , drop = FALSE], na.rm = FALSE)
score_raw <- data.frame(
  geo_accession = colnames(gene_z),
  M1 = score_mean(ECM_MODULES$M1),
  M2 = score_mean(ECM_MODULES$M2),
  M3 = score_mean(ECM_MODULES$M3),
  M4 = score_mean(ECM_MODULES$M4),
  ECM21_gene_equal = score_mean(ECM21_GENES),
  ECM19_structural = score_mean(ECM19_GENES),
  LOMO_minus_M1 = score_mean(setdiff(ECM21_GENES, ECM_MODULES$M1)),
  LOMO_minus_M2 = score_mean(setdiff(ECM21_GENES, ECM_MODULES$M2)),
  LOMO_minus_M3 = score_mean(setdiff(ECM21_GENES, ECM_MODULES$M3)),
  LOMO_minus_M4 = score_mean(setdiff(ECM21_GENES, ECM_MODULES$M4)),
  stringsAsFactors = FALSE,
  row.names = colnames(gene_z)
)
score_raw$ECM21_module_balanced <- rowMeans(score_raw[, c("M1", "M2", "M3", "M4")])

score_columns <- c(
  "M1", "M2", "M3", "M4", "ECM21_gene_equal", "ECM21_module_balanced",
  "ECM19_structural", "LOMO_minus_M1", "LOMO_minus_M2", "LOMO_minus_M3", "LOMO_minus_M4"
)
score_sd <- score_raw
score_scaling <- list()
for (nm in score_columns) {
  mu <- mean(score_raw[disc_ids_all, nm], na.rm = TRUE)
  sig <- stats::sd(score_raw[disc_ids_all, nm], na.rm = TRUE)
  if (!is.finite(sig) || sig <= 0) stop("Invalid Discovery-reference SD for score ", nm)
  score_sd[[nm]] <- (score_raw[[nm]] - mu) / sig
  score_scaling[[nm]] <- data.frame(score = nm, discovery_mean = mu, discovery_sd = sig, stringsAsFactors = FALSE)
}
write_csv(do.call(rbind, score_scaling), file.path(OUTPUT_DIR, "tables", "score_standardisation_discovery_reference.csv"))

# CAF score: prefer explicit/sample-level score; otherwise compute from exact genes.
extract_caf_score_table <- function(path, clinical) {
  if (is.na(path)) return(rep(NA_real_, nrow(clinical)))
  obj <- read_delimited_auto(path)
  if (is.numeric(obj) && !is.null(names(obj))) {
    obj <- data.frame(sample_id = names(obj), CAF_PROXY_30_score = as.numeric(obj), stringsAsFactors = FALSE)
  }
  if (is.list(obj) && !is.data.frame(obj)) {
    data_frames <- Filter(is.data.frame, obj)
    if (length(data_frames)) obj <- data_frames[[1]]
  }
  if (!is.data.frame(obj)) stop("CAF score input must resolve to a data frame.")
  norm <- tolower(gsub("[^a-z0-9]", "", names(obj)))
  score_hit <- which(grepl("caf", norm) & grepl("score|proxy30|cafproxy", norm))[1]
  id_hit <- which(norm %in% c("geoaccession", "gsm", "sampleid", "sample", "title", "sampletitle", "patientid", "citid"))[1]
  if (is.na(score_hit) || is.na(id_hit)) stop("CAF score table needs a sample-ID column and a CAF score column.")
  ids <- normalise_id(obj[[id_hit]])
  vals <- suppressWarnings(as.numeric(obj[[score_hit]]))
  gsm_match <- match(normalise_id(clinical$geo_accession), ids)
  title_match <- match(normalise_id(clinical$sample_title), ids)
  idx <- gsm_match
  idx[is.na(idx)] <- title_match[is.na(idx)]
  vals[idx]
}

caf_source <- "not_available"
caf_raw <- rep(NA_real_, nrow(clinical))
names(caf_raw) <- rownames(clinical)
if (!is.na(caf_score_path)) {
  caf_raw <- extract_caf_score_table(caf_score_path, clinical)
  names(caf_raw) <- rownames(clinical)
  caf_source <- paste0("sample_level_score:", caf_score_path)
} else if (length(caf_genes) == 30L) {
  caf_raw <- score_mean(caf_genes)
  caf_source <- paste0("computed_from_frozen_30_genes:", caf_gene_path)
}
if (sum(is.finite(caf_raw[disc_ids_all])) >= 0.80 * length(disc_ids_all)) {
  caf_mu <- mean(caf_raw[disc_ids_all], na.rm = TRUE)
  caf_sig <- stats::sd(caf_raw[disc_ids_all], na.rm = TRUE)
  if (is.finite(caf_sig) && caf_sig > 0) caf_sd <- (caf_raw - caf_mu) / caf_sig else caf_sd <- rep(NA_real_, length(caf_raw))
} else {
  caf_sd <- rep(NA_real_, length(caf_raw))
}
names(caf_sd) <- rownames(clinical)
caf_available <- sum(is.finite(caf_sd[rownames(clinical)[clinical$is_tumour]])) >= 0.75 * sum(clinical$is_tumour)
if (!caf_available) {
  caf_source <- "not_available"
  if (isTRUE(REQUIRE_CAF_INPUT)) stop("A valid frozen CAF_PROXY_30 score or exact 30-gene file is required but was not found.")
  warning("CAF_PROXY_30 resource was not found or had inadequate coverage. CAF-adjusted models will be explicitly skipped.")
}
write_csv(data.frame(caf_available = caf_available, caf_source = caf_source, stringsAsFactors = FALSE), file.path(OUTPUT_DIR, "tables", "CAF_input_status.csv"))

# Assemble the tumour-only analysis table.
analysis_data <- clinical[clinical$is_tumour, , drop = FALSE]
analysis_data <- cbind(analysis_data, score_sd[rownames(analysis_data), score_columns, drop = FALSE])
analysis_data$caf_value <- caf_sd[rownames(analysis_data)]
analysis_data$dataset <- droplevels(analysis_data$dataset)
write_csv(cbind(analysis_data[, c("geo_accession", "sample_title", "dataset", "age", "sex", "stage", "mmr", "rfs_time_months", "rfs_event")], analysis_data[, score_columns], CAF_PROXY_30 = analysis_data$caf_value), file.path(OUTPUT_DIR, "tables", "GSE39582_sample_level_module_scores_and_clinical.csv"))

# ------------------------------- COX FUNCTIONS --------------------------------
model_terms <- function(model_key, include_score = TRUE) {
  score <- if (include_score) "score_value" else NULL
  switch(
    model_key,
    clinical_primary = c(score, "age", "sex", "stage", "mmr"),
    clinical_stratified = c(score, "age", "mmr", "strata(sex, stage)"),
    caf_primary = c(score, "age", "sex", "stage", "mmr", "caf_value"),
    caf_stratified = c(score, "age", "mmr", "caf_value", "strata(sex, stage)"),
    stop("Unknown model key: ", model_key)
  )
}

model_required_columns <- function(model_key, include_score = TRUE) {
  terms <- model_terms(model_key, include_score)
  terms <- gsub("strata\\(([^,]+),\\s*([^\\)]+)\\)", "\\1+\\2", terms)
  unique(c("rfs_time_months", "rfs_event", unlist(strsplit(terms, "\\+", fixed = FALSE))))
}

safe_cox_fit <- function(dat, score_name = NULL, model_key, include_score = TRUE) {
  d <- dat
  if (include_score) {
    if (is.null(score_name) || !score_name %in% names(d)) return(list(ok = FALSE, reason = "score_missing"))
    d$score_value <- d[[score_name]]
  }
  req <- trimws(model_required_columns(model_key, include_score))
  req <- req[nzchar(req)]
  cc <- stats::complete.cases(d[, req, drop = FALSE])
  d <- d[cc & d$rfs_time_months > 0 & d$rfs_event %in% c(0, 1), , drop = FALSE]
  if (nrow(d) < 30L || sum(d$rfs_event) < 10L) return(list(ok = FALSE, reason = "insufficient_rows_or_events", n = nrow(d), events = sum(d$rfs_event)))
  rhs <- paste(model_terms(model_key, include_score), collapse = " + ")
  form <- stats::as.formula(paste("survival::Surv(rfs_time_months, rfs_event) ~", rhs))
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      survival::coxph(form, data = d, x = TRUE, y = TRUE, model = TRUE, ties = "efron"),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(ok = FALSE, reason = conditionMessage(fit), n = nrow(d), events = sum(d$rfs_event)))
  bad_warning <- any(grepl("converg|infinite|did not converge|loglik converged before", warnings, ignore.case = TRUE))
  cf <- stats::coef(fit)
  vv <- diag(stats::vcov(fit))
  if (bad_warning || any(!is.finite(cf)) || any(!is.finite(vv)) || any(vv <= 0) || any(abs(cf) > 20) || any(sqrt(vv) > 20)) {
    return(list(ok = FALSE, reason = paste(c("estimability_gate_failed", warnings), collapse = " | "), n = nrow(d), events = sum(d$rfs_event)))
  }
  lp <- stats::predict(fit, newdata = d, type = "lp")
  conc <- tryCatch(survival::concordance(survival::Surv(d$rfs_time_months, d$rfs_event) ~ lp, reverse = TRUE)$concordance, error = function(e) NA_real_)
  list(ok = is.finite(conc), fit = fit, data = d, n = nrow(d), events = sum(d$rfs_event), cindex = unname(conc), warnings = warnings)
}

cohort_subset <- function(name) {
  d <- analysis_data[is.finite(analysis_data$rfs_time_months) & analysis_data$rfs_time_months > 0 & analysis_data$rfs_event %in% c(0, 1), , drop = FALSE]
  if (name == "Discovery") d <- d[d$dataset == "Discovery", , drop = FALSE]
  if (name == "Validation") d <- d[d$dataset == "Validation", , drop = FALSE]
  droplevels(d)
}

model_keys <- c("clinical_primary", "clinical_stratified")
if (caf_available) model_keys <- c(model_keys, "caf_primary", "caf_stratified")

cox_rows <- list()
cindex_rows <- list()
for (cohort_name in c("Discovery", "Validation", "Pooled")) {
  dd <- cohort_subset(cohort_name)
  for (mk in model_keys) {
    base_fit <- safe_cox_fit(dd, model_key = mk, include_score = FALSE)
    for (score_name in score_columns) {
      sf <- safe_cox_fit(dd, score_name = score_name, model_key = mk, include_score = TRUE)
      if (!isTRUE(sf$ok)) {
        cox_rows[[length(cox_rows) + 1L]] <- data.frame(
          cohort = cohort_name, model = mk, score = score_name, status = "failed",
          reason = sf$reason, n = ifelse(is.null(sf$n), NA, sf$n), events = ifelse(is.null(sf$events), NA, sf$events),
          HR = NA, CI_lower = NA, CI_upper = NA, p_value = NA, stringsAsFactors = FALSE
        )
        next
      }
      beta <- stats::coef(sf$fit)["score_value"]
      se <- sqrt(stats::vcov(sf$fit)["score_value", "score_value"])
      p <- summary(sf$fit)$coefficients["score_value", "Pr(>|z|)"]
      cox_rows[[length(cox_rows) + 1L]] <- data.frame(
        cohort = cohort_name, model = mk, score = score_name, status = "completed", reason = "",
        n = sf$n, events = sf$events, HR = exp(beta), CI_lower = exp(beta - 1.96 * se),
        CI_upper = exp(beta + 1.96 * se), p_value = p, stringsAsFactors = FALSE
      )
      cindex_rows[[length(cindex_rows) + 1L]] <- data.frame(
        cohort = cohort_name, model = mk, score = score_name,
        n = sf$n, events = sf$events,
        base_Cindex = if (isTRUE(base_fit$ok)) base_fit$cindex else NA_real_,
        score_Cindex = sf$cindex,
        incremental_Cindex_vs_base = if (isTRUE(base_fit$ok)) sf$cindex - base_fit$cindex else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
}
cox_results <- do.call(rbind, cox_rows)
cindex_results <- do.call(rbind, cindex_rows)

# BH correction is restricted to M1-M4 within each cohort/model endpoint.
cox_results$BH_FDR_M1_to_M4 <- NA_real_
groups <- interaction(cox_results$cohort, cox_results$model, drop = TRUE)
for (g in levels(groups)) {
  idx <- which(groups == g & cox_results$score %in% names(ECM_MODULES) & is.finite(cox_results$p_value))
  if (length(idx)) cox_results$BH_FDR_M1_to_M4[idx] <- p.adjust(cox_results$p_value[idx], method = "BH")
}
write_csv(cox_results, file.path(OUTPUT_DIR, "tables", "GSE39582_module_and_aggregate_Cox_results.csv"))
write_csv(cindex_results, file.path(OUTPUT_DIR, "tables", "GSE39582_apparent_incremental_Cindex.csv"))

# Select the best single module once in Discovery using the prespecified primary
# clinical model. Freeze that identity for Validation and all bootstrap contrasts.
selection_pool <- cindex_results[
  cindex_results$cohort == "Discovery" & cindex_results$model == "clinical_primary" &
    cindex_results$score %in% names(ECM_MODULES) & is.finite(cindex_results$score_Cindex), , drop = FALSE
]
if (!nrow(selection_pool)) stop("No estimable Discovery primary model was available for best-module freezing.")
selection_pool <- selection_pool[order(-selection_pool$score_Cindex, selection_pool$score), , drop = FALSE]
BEST_MODULE <- selection_pool$score[1]
selection_table <- transform(selection_pool, selected_and_frozen = score == BEST_MODULE)
write_csv(selection_table, file.path(OUTPUT_DIR, "tables", "best_single_module_selection_discovery_only.csv"))
log_line("Best single module frozen from Discovery apparent C-index: ", BEST_MODULE)

# ------------------------ PAIRED BOOTSTRAP C-INDEX ----------------------------
bootstrap_cindices <- function(dat, model_key, scores, reps, seed) {
  required <- unique(c(model_required_columns(model_key, include_score = FALSE), scores))
  required <- trimws(required[nzchar(required)])
  d <- dat[stats::complete.cases(dat[, required, drop = FALSE]) & dat$rfs_time_months > 0 & dat$rfs_event %in% c(0, 1), , drop = FALSE]
  apparent <- setNames(rep(NA_real_, length(scores)), scores)
  for (s in scores) {
    fit <- safe_cox_fit(d, score_name = s, model_key = model_key, include_score = TRUE)
    if (isTRUE(fit$ok)) apparent[s] <- fit$cindex
  }
  set.seed(seed)
  boot <- matrix(NA_real_, nrow = reps, ncol = length(scores), dimnames = list(NULL, scores))
  n <- nrow(d)
  for (b in seq_len(reps)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    db <- d[idx, , drop = FALSE]
    if (sum(db$rfs_event) < 10L) next
    for (s in scores) {
      fit <- safe_cox_fit(db, score_name = s, model_key = model_key, include_score = TRUE)
      if (isTRUE(fit$ok)) boot[b, s] <- fit$cindex
    }
    if (b %% 100L == 0L) log_line("Bootstrap ", model_key, ": ", b, "/", reps)
  }
  list(data = d, apparent = apparent, bootstrap = boot)
}

comparison_specs <- data.frame(
  score_A = c(
    rep("ECM21_gene_equal", 7L),
    "ECM21_module_balanced"
  ),
  score_B = c(
    BEST_MODULE, "M1", "M2", "M3", "M4", "ECM19_structural", "ECM21_module_balanced",
    BEST_MODULE
  ),
  comparison = c(
    "ECM21_gene_equal_minus_best_module_frozen",
    "ECM21_gene_equal_minus_M1", "ECM21_gene_equal_minus_M2",
    "ECM21_gene_equal_minus_M3", "ECM21_gene_equal_minus_M4",
    "ECM21_gene_equal_minus_ECM19_structural",
    "ECM21_gene_equal_minus_ECM21_module_balanced",
    "ECM21_module_balanced_minus_best_module_frozen"
  ),
  stringsAsFactors = FALSE
)
bootstrap_scores <- unique(c(comparison_specs$score_A, comparison_specs$score_B))
bootstrap_rows <- list()
bootstrap_model_keys <- intersect(BOOTSTRAP_MODEL_KEYS, model_keys)
seed_counter <- 0L
for (cohort_name in c("Discovery", "Validation", "Pooled")) {
  dd <- cohort_subset(cohort_name)
  for (mk in bootstrap_model_keys) {
    seed_counter <- seed_counter + 1L
    log_line("Starting paired bootstrap: ", cohort_name, " / ", mk)
    bobj <- bootstrap_cindices(dd, mk, bootstrap_scores, CINDEX_BOOTSTRAP_REPS, BOOTSTRAP_SEED + seed_counter)
    for (i in seq_len(nrow(comparison_specs))) {
      a <- comparison_specs$score_A[i]
      b <- comparison_specs$score_B[i]
      delta <- bobj$bootstrap[, a] - bobj$bootstrap[, b]
      good <- is.finite(delta)
      success_fraction <- mean(good)
      q <- if (sum(good) >= 100L) stats::quantile(delta[good], c(0.025, 0.5, 0.975), names = FALSE, type = 6) else rep(NA_real_, 3)
      bootstrap_rows[[length(bootstrap_rows) + 1L]] <- data.frame(
        cohort = cohort_name, model = mk, comparison = comparison_specs$comparison[i],
        score_A = a, score_B = b, n = nrow(bobj$data), events = sum(bobj$data$rfs_event),
        bootstrap_reps_requested = CINDEX_BOOTSTRAP_REPS,
        bootstrap_reps_successful = sum(good), bootstrap_success_fraction = success_fraction,
        Cindex_A_apparent = bobj$apparent[a], Cindex_B_apparent = bobj$apparent[b],
        delta_Cindex_apparent = bobj$apparent[a] - bobj$apparent[b],
        delta_Cindex_boot_median = q[2], CI_lower = q[1], CI_upper = q[3],
        bootstrap_probability_delta_gt_0 = if (sum(good)) mean(delta[good] > 0) else NA_real_,
        reliable_bootstrap = success_fraction >= MIN_BOOTSTRAP_SUCCESS_FRACTION,
        stringsAsFactors = FALSE
      )
    }
  }
}
bootstrap_results <- do.call(rbind, bootstrap_rows)
write_csv(bootstrap_results, file.path(OUTPUT_DIR, "tables", "paired_bootstrap_delta_Cindex_results.csv"))

# -------------------------- CROSS-COHORT CLAIM GATES --------------------------
make_claim_gate <- function(model_key, comparison_name) {
  x <- bootstrap_results[bootstrap_results$model == model_key & bootstrap_results$comparison == comparison_name & bootstrap_results$cohort %in% c("Discovery", "Validation", "Pooled"), , drop = FALSE]
  get_row <- function(cohort) x[x$cohort == cohort, , drop = FALSE]
  d <- get_row("Discovery")
  v <- get_row("Validation")
  p <- get_row("Pooled")
  enough <- nrow(d) == 1L && nrow(v) == 1L && nrow(p) == 1L
  direction_stable <- enough && isTRUE(d$delta_Cindex_apparent > 0) && isTRUE(v$delta_Cindex_apparent > 0)
  positive_gain <- enough && direction_stable && isTRUE(p$delta_Cindex_apparent > 0)
  ci_support_all <- enough && isTRUE(d$CI_lower > 0) && isTRUE(v$CI_lower > 0) && isTRUE(p$CI_lower > 0)
  reliable <- enough && isTRUE(d$reliable_bootstrap) && isTRUE(v$reliable_bootstrap) && isTRUE(p$reliable_bootstrap)
  allow_enhancement <- positive_gain && ci_support_all && direction_stable && reliable
  data.frame(
    model = model_key,
    comparison = comparison_name,
    discovery_delta = if (nrow(d)) d$delta_Cindex_apparent else NA,
    discovery_CI_lower = if (nrow(d)) d$CI_lower else NA,
    discovery_CI_upper = if (nrow(d)) d$CI_upper else NA,
    validation_delta = if (nrow(v)) v$delta_Cindex_apparent else NA,
    validation_CI_lower = if (nrow(v)) v$CI_lower else NA,
    validation_CI_upper = if (nrow(v)) v$CI_upper else NA,
    pooled_delta = if (nrow(p)) p$delta_Cindex_apparent else NA,
    pooled_CI_lower = if (nrow(p)) p$CI_lower else NA,
    pooled_CI_upper = if (nrow(p)) p$CI_upper else NA,
    gain_positive = positive_gain,
    bootstrap_CI_support_all = ci_support_all,
    cross_cohort_direction_stable = direction_stable,
    bootstrap_reliability_gate = reliable,
    enhancement_word_allowed = allow_enhancement,
    decision = if (allow_enhancement) {
      "All prespecified gates passed; an enhancement statement is permitted for this exact comparison and model."
    } else {
      "Do not write enhanced: at least one required gain, CI, cross-cohort-direction, or bootstrap-reliability gate failed."
    },
    stringsAsFactors = FALSE
  )
}

gate_rows <- list()
for (mk in bootstrap_model_keys) {
  gate_rows[[length(gate_rows) + 1L]] <- make_claim_gate(mk, "ECM21_gene_equal_minus_best_module_frozen")
  gate_rows[[length(gate_rows) + 1L]] <- make_claim_gate(mk, "ECM21_gene_equal_minus_ECM19_structural")
}
claim_gates <- do.call(rbind, gate_rows)
write_csv(claim_gates, file.path(OUTPUT_DIR, "tables", "CLAIM_GATE_enhancement_decision.csv"))

# Cox-direction concordance is reported separately and is not substituted for
# a positive, CI-supported delta C-index.
direction_rows <- list()
for (mk in model_keys) {
  for (s in score_columns) {
    d <- cox_results[cox_results$cohort == "Discovery" & cox_results$model == mk & cox_results$score == s, , drop = FALSE]
    v <- cox_results[cox_results$cohort == "Validation" & cox_results$model == mk & cox_results$score == s, , drop = FALSE]
    direction_rows[[length(direction_rows) + 1L]] <- data.frame(
      model = mk, score = s,
      discovery_HR = if (nrow(d)) d$HR else NA,
      validation_HR = if (nrow(v)) v$HR else NA,
      same_logHR_direction = nrow(d) == 1L && nrow(v) == 1L && is.finite(d$HR) && is.finite(v$HR) && sign(log(d$HR)) == sign(log(v$HR)),
      stringsAsFactors = FALSE
    )
  }
}
write_csv(do.call(rbind, direction_rows), file.path(OUTPUT_DIR, "tables", "cross_cohort_Cox_direction_concordance.csv"))

# --------------------------------- FIGURES ------------------------------------
forest_plot <- function(df, estimate, lower, upper, labels, null, xlab, title, file_stub, log_x = FALSE) {
  df <- df[is.finite(df[[estimate]]) & is.finite(df[[lower]]) & is.finite(df[[upper]]), , drop = FALSE]
  if (!nrow(df)) return(invisible(NULL))
  draw <- function(device_fun, path) {
    device_fun(path, width = 1800, height = max(1200, 85 * nrow(df) + 500), res = 180)
    on.exit(grDevices::dev.off(), add = TRUE)
    op <- par(mar = c(5, 15, 4, 2))
    on.exit(par(op), add = TRUE)
    y <- rev(seq_len(nrow(df)))
    xr <- range(c(df[[lower]], df[[upper]], null), finite = TRUE)
    if (log_x) xr <- range(log(c(df[[lower]], df[[upper]], null)), finite = TRUE)
    est <- if (log_x) log(df[[estimate]]) else df[[estimate]]
    lo <- if (log_x) log(df[[lower]]) else df[[lower]]
    hi <- if (log_x) log(df[[upper]]) else df[[upper]]
    null_draw <- if (log_x) log(null) else null
    plot(est, y, xlim = xr, ylim = c(0.5, nrow(df) + 0.5), yaxt = "n", xaxt = if (log_x) "n" else "s", ylab = "", xlab = xlab, pch = 19, main = title)
    segments(lo, y, hi, y, lwd = 2)
    abline(v = null_draw, lty = 2, col = "grey40")
    axis(2, at = y, labels = df[[labels]], las = 1, cex.axis = 0.75)
    if (log_x) {
      ticks <- pretty(c(df[[lower]], df[[upper]], null))
      ticks <- ticks[ticks > 0]
      axis(1, at = log(ticks), labels = format(ticks, digits = 3))
    }
  }
  png_path <- file.path(OUTPUT_DIR, "figures", paste0(file_stub, ".png"))
  draw(function(path, width, height, res) grDevices::png(path, width = width, height = height, res = res), png_path)
  pdf_path <- file.path(OUTPUT_DIR, "figures", paste0(file_stub, ".pdf"))
  grDevices::pdf(pdf_path, width = 11, height = max(7, 0.42 * nrow(df) + 2))
  op <- par(mar = c(5, 15, 4, 2))
  y <- rev(seq_len(nrow(df)))
  est <- if (log_x) log(df[[estimate]]) else df[[estimate]]
  lo <- if (log_x) log(df[[lower]]) else df[[lower]]
  hi <- if (log_x) log(df[[upper]]) else df[[upper]]
  null_draw <- if (log_x) log(null) else null
  xr <- range(c(lo, hi, null_draw), finite = TRUE)
  plot(est, y, xlim = xr, ylim = c(0.5, nrow(df) + 0.5), yaxt = "n", xaxt = if (log_x) "n" else "s", ylab = "", xlab = xlab, pch = 19, main = title)
  segments(lo, y, hi, y, lwd = 2)
  abline(v = null_draw, lty = 2, col = "grey40")
  axis(2, at = y, labels = df[[labels]], las = 1, cex.axis = 0.75)
  if (log_x) {
    ticks <- pretty(c(df[[lower]], df[[upper]], null))
    ticks <- ticks[ticks > 0]
    axis(1, at = log(ticks), labels = format(ticks, digits = 3))
  }
  par(op)
  grDevices::dev.off()
}

cox_plot_data <- cox_results[
  cox_results$model == "clinical_stratified" & cox_results$score %in% c(names(ECM_MODULES), "ECM21_gene_equal", "ECM21_module_balanced", "ECM19_structural"), , drop = FALSE
]
cox_plot_data$plot_label <- paste(cox_plot_data$cohort, cox_plot_data$score, sep = " | ")
forest_plot(
  cox_plot_data, "HR", "CI_lower", "CI_upper", "plot_label", 1,
  "Hazard ratio per Discovery-reference SD", "GSE39582 module and aggregate RFS Cox models",
  "Figure_GSE39582_module_and_aggregate_Cox_forest", log_x = TRUE
)

delta_plot_data <- bootstrap_results[
  bootstrap_results$model == "clinical_stratified" & bootstrap_results$comparison %in% c(
    "ECM21_gene_equal_minus_best_module_frozen",
    "ECM21_gene_equal_minus_ECM19_structural",
    "ECM21_gene_equal_minus_ECM21_module_balanced"
  ), , drop = FALSE
]
delta_plot_data$plot_label <- paste(delta_plot_data$cohort, delta_plot_data$comparison, sep = " | ")
forest_plot(
  delta_plot_data, "delta_Cindex_apparent", "CI_lower", "CI_upper", "plot_label", 0,
  "Paired delta C-index", "Paired bootstrap comparison of ECM21 aggregation",
  "Figure_GSE39582_bootstrap_delta_Cindex_forest", log_x = FALSE
)

# ----------------------------- TEXT HANDOFF FILES -----------------------------
primary_cox <- cox_results[cox_results$model == "clinical_stratified" & cox_results$score %in% c(names(ECM_MODULES), "ECM21_gene_equal", "ECM21_module_balanced", "ECM19_structural"), , drop = FALSE]
fmt <- function(x, digits = 3) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
methods_text <- c(
  "GSE39582 MODULE SURVIVAL METHODS - READY TO EDIT",
  "",
  paste0("The official GSE39582 GPL570 processed series matrix was imported with GEOquery. Non-tumour mucosa samples were excluded before scoring. Unambiguous probes were mapped with hgu133plus2.db; for genes represented by multiple probes, the probe with the largest interquartile range in the Discovery tumours was selected, with probe ID as a deterministic tie-breaker."),
  "",
  "Gene expression was standardised using Discovery-tumour means and standard deviations and transferred unchanged to Validation. M1-M4 scores were equal-weight means of their member-gene z values. ECM21_gene_equal was the equal-weight mean of all 21 genes; ECM21_module_balanced was mean(M1,M2,M3,M4); ECM19_structural omitted PTN and MDK. All Cox effects were reported per one Discovery-reference score SD.",
  "",
  "The primary clinical Cox model included score, age, sex, TNM stage and mismatch-repair status. A prespecified sensitivity model retained score, age and mismatch-repair status and stratified the baseline hazard jointly by sex and stage. When the exact frozen CAF_PROXY_30 resource was available, otherwise identical CAF-adjusted models were added. Benjamini-Hochberg correction was applied only across M1-M4 within each cohort/model endpoint.",
  "",
  paste0("The single-module comparator was selected once in Discovery as the module with the largest apparent Harrell C-index under the primary clinical model (", BEST_MODULE, ") and was then frozen. Patients were resampled with replacement within each cohort. In every bootstrap replicate, comparator models were fitted to the same resampled patients and delta C-index was calculated as C-index(score A) minus C-index(score B). Percentile 95% confidence intervals used ", CINDEX_BOOTSTRAP_REPS, " requested replicates. An enhancement statement required positive gains in Discovery, Validation and Pooled analyses, bootstrap CIs above zero in all three, stable Discovery/Validation direction, and at least ", round(100 * MIN_BOOTSTRAP_SUCCESS_FRACTION), "% successful replicates."),
  ""
)
writeLines(methods_text, file.path(OUTPUT_DIR, "METHODS_READY_TO_EDIT.txt"), useBytes = TRUE)

results_lines <- c(
  "GSE39582 MODULE SURVIVAL NUMERIC HANDOFF",
  "",
  paste0("Best single module frozen from Discovery: ", BEST_MODULE),
  paste0("CAF-adjusted analysis available: ", caf_available, "; source: ", caf_source),
  "",
  "Clinical-stratified Cox results:",
  apply(primary_cox, 1, function(r) paste0(
    r[["cohort"]], " | ", r[["score"]], " | HR ", fmt(as.numeric(r[["HR"]])),
    " (95% CI ", fmt(as.numeric(r[["CI_lower"]])), "-", fmt(as.numeric(r[["CI_upper"]])),
    "); P=", fmt(as.numeric(r[["p_value"]]), 4)
  )),
  "",
  "Enhancement decision gates:",
  apply(claim_gates, 1, function(r) paste0(
    r[["model"]], " | ", r[["comparison"]], " | allowed=", r[["enhancement_word_allowed"]],
    " | ", r[["decision"]]
  ))
)
writeLines(results_lines, file.path(OUTPUT_DIR, "AUTOMATED_RESULTS_NUMERIC_SUMMARY.txt"), useBytes = TRUE)

status <- data.frame(
  component = c("GSE39582 input", "module Cox", "paired bootstrap delta C-index", "CAF-adjusted analysis", "claim gate"),
  status = c("completed", "completed", "completed", if (caf_available) "completed" else "not_run_missing_frozen_input", "completed"),
  detail = c(
    matrix_file,
    paste0(nrow(cox_results), " rows"),
    paste0(nrow(bootstrap_results), " rows; ", CINDEX_BOOTSTRAP_REPS, " requested replicates per cohort/model"),
    caf_source,
    paste0(sum(claim_gates$enhancement_word_allowed, na.rm = TRUE), " comparisons passed all strict gates")
  ),
  stringsAsFactors = FALSE
)
write_csv(status, file.path(OUTPUT_DIR, "RUN_STATUS.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "R_sessionInfo.txt"), useBytes = TRUE)

# Archive creation is handled by the repository release workflow.
zip_ok <- FALSE
log_line("Analysis completed")
message("\n============================================================")
message("GSE39582 module Cox and bootstrap delta C-index completed.")
message("Results folder: ", OUTPUT_DIR)
if (zip_ok) message("Return ZIP: ", zip_path)
message("Upload RUN_STATUS.csv, AUTOMATED_RESULTS_NUMERIC_SUMMARY.txt, and the tables folder for manuscript updating.")
message("============================================================\n")
