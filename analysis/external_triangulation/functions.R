options(stringsAsFactors = FALSE)

require_packages <- function(packages) {
  ok <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    stop(
      "Missing R packages: ", paste(packages[!ok], collapse = ", "),
      "\nInstall them manually using the commands in README_FIRST_CN.md, then rerun."
    )
  }
  invisible(TRUE)
}

result_dir <- function(name) {
  out <- file.path(RESULT_ROOT, name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.frame(x), path, na = "NA")
  invisible(path)
}

write_status <- function(step, status, detail = "") {
  path <- file.path(LOG_ROOT, paste0(step, "_status.csv"))
  write_csv(
    data.frame(
      step = step,
      status = status,
      detail = detail,
      timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    path
  )
}

assert_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop("Missing ", label, ":\n", path, "\nDownload it manually from DOWNLOAD_LINKS.md.")
  }
  invisible(path)
}

read_table_auto <- function(path, ...) {
  assert_file(path)
  data.table::fread(path, data.table = FALSE, check.names = FALSE, ...)
}

clean_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\.[0-9]+$", "", x)
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

read_gene_set <- function(path, gene_column = "gene") {
  x <- read.csv(path, check.names = FALSE)
  if (!gene_column %in% names(x)) stop("Gene column not found in ", path)
  unique(stats::na.omit(clean_symbol(x[[gene_column]])))
}

guess_gene_column <- function(x) {
  nms <- names(x)
  hit <- grep("^(gene|symbol|gene_symbol|hgnc|external_gene_name|protein)$", nms, ignore.case = TRUE)
  if (length(hit)) return(nms[hit[1]])
  character_cols <- nms[vapply(x, function(z) is.character(z) || is.factor(z), logical(1))]
  if (length(character_cols)) return(character_cols[1])
  nms[1]
}

numeric_expression_matrix <- function(path, duplicate_action = c("sum", "mean")) {
  duplicate_action <- match.arg(duplicate_action)
  x <- read_table_auto(path)
  gene_col <- guess_gene_column(x)
  genes <- clean_symbol(x[[gene_col]])
  x[[gene_col]] <- NULL
  keep_numeric <- vapply(x, is.numeric, logical(1))
  if (!any(keep_numeric)) {
    x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
    keep_numeric <- vapply(x, function(z) sum(is.finite(z)) > 0, logical(1))
  }
  mat <- as.matrix(x[, keep_numeric, drop = FALSE])
  storage.mode(mat) <- "numeric"
  keep <- !is.na(genes) & rowSums(is.finite(mat)) > 0
  genes <- genes[keep]
  mat <- mat[keep, , drop = FALSE]
  split_rows <- split(seq_along(genes), genes)
  collapsed <- vapply(
    split_rows,
    function(ii) {
      if (duplicate_action == "sum") colSums(mat[ii, , drop = FALSE], na.rm = TRUE)
      else colMeans(mat[ii, , drop = FALSE], na.rm = TRUE)
    },
    numeric(ncol(mat))
  )
  collapsed <- t(collapsed)
  colnames(collapsed) <- colnames(mat)
  collapsed
}

log_cpm <- function(counts, prior.count = 1) {
  edgeR::cpm(counts, log = TRUE, prior.count = prior.count)
}

zscore_rows <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sig <- apply(mat, 1, stats::sd, na.rm = TRUE)
  sig[!is.finite(sig) | sig == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sig, "/")
}

score_gene_set <- function(expression, genes, minimum_coverage = 0.60, weights = NULL) {
  genes <- unique(clean_symbol(genes))
  rownames(expression) <- clean_symbol(rownames(expression))
  present <- intersect(genes, rownames(expression))
  coverage <- length(present) / length(genes)
  if (coverage < minimum_coverage) {
    stop(
      sprintf(
        "Gene-set coverage %.1f%% is below the pre-specified %.1f%% threshold.",
        100 * coverage, 100 * minimum_coverage
      )
    )
  }
  z <- zscore_rows(expression[present, , drop = FALSE])
  if (is.null(weights)) {
    score <- colMeans(z, na.rm = TRUE)
  } else {
    w <- weights[match(present, genes)]
    w[!is.finite(w)] <- 0
    w <- w / sum(w)
    score <- as.numeric(crossprod(w, z))
    names(score) <- colnames(z)
  }
  list(
    score = score,
    present = present,
    missing = setdiff(genes, present),
    coverage = coverage
  )
}

save_gene_set_audit <- function(score_object, path, set_name) {
  all_genes <- c(score_object$present, score_object$missing)
  write_csv(
    data.frame(
      gene_set = set_name,
      gene = all_genes,
      detected = all_genes %in% score_object$present,
      coverage = score_object$coverage
    ),
    path
  )
}

fit_edger <- function(counts, design, contrast, robust = TRUE) {
  y <- edgeR::DGEList(counts = round(counts))
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design, robust = robust)
  fit <- edgeR::glmQLFit(y, design, robust = robust)
  test <- edgeR::glmQLFTest(fit, contrast = contrast)
  tab <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  tab[, c("gene", setdiff(names(tab), "gene")), drop = FALSE]
}

paired_effect <- function(scores, patient, condition, reference, comparison) {
  x <- data.frame(score = as.numeric(scores), patient = patient, condition = condition)
  x <- aggregate(score ~ patient + condition, x, mean, na.rm = TRUE)
  wide <- reshape(x, idvar = "patient", timevar = "condition", direction = "wide")
  ref <- paste0("score.", reference)
  cmp <- paste0("score.", comparison)
  if (!all(c(ref, cmp) %in% names(wide))) stop("Both paired conditions are not available.")
  wide$delta <- wide[[cmp]] - wide[[ref]]
  wide <- wide[is.finite(wide$delta), , drop = FALSE]
  if (nrow(wide) < 2) stop("Fewer than two complete patient pairs.")
  test <- stats::t.test(wide$delta, mu = 0)
  list(
    pairs = wide,
    summary = data.frame(
      n_pairs = nrow(wide),
      comparison = paste(comparison, "minus", reference),
      mean_difference = mean(wide$delta),
      ci_low = unname(test$conf.int[1]),
      ci_high = unname(test$conf.int[2]),
      p_value = test$p.value,
      sign_concordance = mean(wide$delta > 0)
    )
  )
}

best_matching_column <- function(annotation, identifiers) {
  scores <- vapply(
    annotation,
    function(z) sum(as.character(z) %in% identifiers),
    integer(1)
  )
  if (!length(scores) || max(scores) == 0) return(NA_character_)
  names(which.max(scores))
}

first_name_matching <- function(x, patterns, required = TRUE, label = "column") {
  nms <- names(x)
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  if (required) stop("Could not infer ", label, ". Available columns: ", paste(nms, collapse = ", "))
  NA_character_
}

standardize_condition <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x0))
  out[grepl("normal|adjacent|(^|[-_])n($|[-_])", x0)] <- "Normal"
  out[grepl("border|margin|(^|[-_])b($|[-_])", x0)] <- "Border"
  out[grepl("tumou?r|core|cancer|(^|[-_])t($|[-_])", x0)] <- "Tumor"
  out
}

standardize_major_type <- function(x) {
  x0 <- tolower(as.character(x))
  out <- rep("Other", length(x0))
  out[grepl("fibro|stromal|myofibro|caf", x0)] <- "Fibroblast"
  out[grepl("epithelial|tumou?r|malignant|colonocyte|goblet|enterocyte|stem", x0)] <- "Epithelial"
  out
}

pseudobulk_counts <- function(counts, sample_id, minimum_cells = 20L) {
  stopifnot(length(sample_id) == ncol(counts))
  groups <- split(seq_len(ncol(counts)), sample_id)
  n_cells <- lengths(groups)
  keep <- n_cells >= minimum_cells
  groups <- groups[keep]
  if (!length(groups)) stop("No pseudobulk group met the cell-count threshold.")
  pb <- vapply(groups, function(ii) Matrix::rowSums(counts[, ii, drop = FALSE]), numeric(nrow(counts)))
  rownames(pb) <- rownames(counts)
  list(counts = pb, n_cells = data.frame(sample_id = names(groups), n_cells = n_cells[keep]))
}

open_parquet_maybe_gz <- function(path) {
  assert_file(path)
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) return(arrow::read_parquet(path))
  require_packages("R.utils")
  cache <- file.path(dirname(path), sub("\\.gz$", "", basename(path), ignore.case = TRUE))
  if (!file.exists(cache)) R.utils::gunzip(path, destname = cache, overwrite = FALSE, remove = FALSE)
  arrow::read_parquet(cache)
}

save_plot <- function(plot, stem, width = 8, height = 6) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 400, bg = "white")
}

theme_publication <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot",
      strip.background = ggplot2::element_rect(fill = "grey92"),
      legend.position = "right"
    )
}

bh <- function(p) stats::p.adjust(p, method = "BH")

cohen_d_paired <- function(delta) {
  if (length(delta) < 2 || stats::sd(delta) == 0) return(NA_real_)
  mean(delta) / stats::sd(delta)
}

manual_mapping_stop <- function(x, out_file, message) {
  write_csv(x, out_file)
  stop(message, "\nAn audit/template was written to: ", out_file)
}

