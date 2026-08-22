PIPELINE_ROOT <- normalizePath(Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))
require_packages(c("Seurat", "Matrix", "arrow", "R.utils", "RANN", "data.table", "ggplot2"))
set.seed(ANALYSIS_OPTIONS$random_seed)
out <- result_dir("04_GSE280315_spatial")
signature <- read_gene_set(SIGNATURE_FILE)
target <- read_gene_set(TARGET_FILE)

normalize_subset <- function(counts, genes, columns) {
  present <- intersect(clean_symbol(genes), clean_symbol(rownames(counts)))
  index <- match(present, clean_symbol(rownames(counts)))
  sub <- counts[index, columns, drop = FALSE]
  rownames(sub) <- present
  lib <- Matrix::colSums(counts[, columns, drop = FALSE])
  lib[lib <= 0] <- 1
  norm <- Matrix::t(Matrix::t(sub) / lib * 10000)
  norm@x <- log1p(norm@x)
  norm
}

section_summaries <- list()
neighbour_rows <- list()
label_audits <- list()
section_plots <- list()

for (i in seq_len(nrow(SPATIAL_SAMPLES))) {
  section <- SPATIAL_SAMPLES$section[i]
  h5 <- SPATIAL_SAMPLES$h5[i]
  parquet <- SPATIAL_SAMPLES$metadata[i]
  assert_file(h5, paste(section, "H5"))
  assert_file(parquet, paste(section, "metadata"))
  message("Spatial section: ", section)

  counts <- Seurat::Read10X_h5(h5, use.names = TRUE, unique.features = TRUE)
  if (is.list(counts)) {
    assay_name <- if ("Gene Expression" %in% names(counts)) "Gene Expression" else names(counts)[1]
    counts <- counts[[assay_name]]
  }
  meta <- as.data.frame(open_parquet_maybe_gz(parquet), check.names = FALSE)
  barcode_col <- first_name_matching(meta, c("^barcode$", "barcodes", "cell"), TRUE, "spatial barcode")
  barcode <- as.character(meta[[barcode_col]])
  idx <- match(colnames(counts), barcode)
  keep_match <- !is.na(idx)
  if (mean(keep_match) < 0.50) {
    # Accommodate deposits that differ only by the terminal -1 suffix.
    idx <- match(sub("-1$", "", colnames(counts)), sub("-1$", "", barcode))
    keep_match <- !is.na(idx)
  }
  if (mean(keep_match) < 0.50) stop(section, ": fewer than 50% of H5 barcodes match Metadata.parquet.")
  counts <- counts[, keep_match, drop = FALSE]
  meta <- meta[idx[keep_match], , drop = FALSE]

  x_col <- first_name_matching(meta, c("^X$", "x_coord", "pxl_col", "array_col", "column"), TRUE, "X coordinate")
  y_col <- first_name_matching(meta, c("^Y$", "y_coord", "pxl_row", "array_row", "row"), TRUE, "Y coordinate")
  label_cols <- grep("Deconvolution|Label|Unsupervised|cluster|class", names(meta), ignore.case = TRUE, value = TRUE)
  if (!length(label_cols)) {
    manual_mapping_stop(
      data.frame(column = names(meta), example = vapply(meta, function(z) paste(head(z, 3), collapse = " | "), character(1))),
      file.path(out, paste0(section, "_metadata_columns_NEEDS_REVIEW.csv")),
      paste(section, "has no recognizable cell-label columns.")
    )
  }
  primary_candidates <- grep("^DeconvolutionLabel1$|^Label1$|primary.*label|cell.?type", names(meta), ignore.case = TRUE, value = TRUE)
  primary_label_col <- if (length(primary_candidates)) primary_candidates[1] else label_cols[1]
  primary_label <- as.character(meta[[primary_label_col]])
  cell_class <- rep("Other", length(primary_label))
  cell_class[grepl(SPATIAL_LABEL_REGEX$fibroblast, primary_label, ignore.case = TRUE)] <- "Fibroblast"
  cell_class[grepl(SPATIAL_LABEL_REGEX$epithelial, primary_label, ignore.case = TRUE)] <- "Epithelial"
  class_candidates <- grep("^DeconvolutionClass$|purity|singlet", names(meta), ignore.case = TRUE, value = TRUE)
  if (length(class_candidates)) {
    class_text <- as.character(meta[[class_candidates[1]]])
    if (sum(grepl("singlet", class_text, ignore.case = TRUE), na.rm = TRUE) >= 100) {
      cell_class[!grepl("singlet", class_text, ignore.case = TRUE)] <- "Other"
    }
  }

  label_audits[[section]] <- cbind(
    section = section,
    primary_label_column = primary_label_col,
    as.data.frame(table(primary_label = primary_label, assigned_class = cell_class))
  )
  write_csv(
    label_audits[[section]],
    file.path(out, paste0(section, "_label_frequency_audit.csv"))
  )
  fib <- which(cell_class == "Fibroblast")
  epi <- which(cell_class == "Epithelial")
  if (length(fib) < 50 || length(epi) < 50) {
    stop(section, ": fewer than 50 fibroblast or epithelial bins after pre-specified label mapping. Review label audit.")
  }
  # Fixed-size computational subsample is drawn before scores/outcomes are calculated.
  if (length(epi) > ANALYSIS_OPTIONS$spatial_max_epithelial_bins_per_section) {
    epi <- sort(sample(epi, ANALYSIS_OPTIONS$spatial_max_epithelial_bins_per_section))
  }

  fib_norm <- normalize_subset(counts, signature, fib)
  epi_norm <- normalize_subset(counts, target, epi)
  fib_score_obj <- score_gene_set(as.matrix(fib_norm), signature, ANALYSIS_OPTIONS$minimum_signature_coverage)
  epi_score_obj <- score_gene_set(as.matrix(epi_norm), target, ANALYSIS_OPTIONS$minimum_target_coverage)
  save_gene_set_audit(fib_score_obj, file.path(out, paste0(section, "_ECM_coverage.csv")), "ECM_LIGAND_21")
  save_gene_set_audit(epi_score_obj, file.path(out, paste0(section, "_target_coverage.csv")), "EPITHELIAL_TUMOR_UP_50")

  fib_xy <- cbind(as.numeric(meta[[x_col]][fib]), as.numeric(meta[[y_col]][fib]))
  epi_xy <- cbind(as.numeric(meta[[x_col]][epi]), as.numeric(meta[[y_col]][epi]))
  finite_fib <- rowSums(!is.finite(fib_xy)) == 0
  finite_epi <- rowSums(!is.finite(epi_xy)) == 0
  fib_xy <- fib_xy[finite_fib, , drop = FALSE]
  fib_scores <- fib_score_obj$score[finite_fib]
  epi_xy <- epi_xy[finite_epi, , drop = FALSE]
  epi_scores <- epi_score_obj$score[finite_epi]

  k <- min(ANALYSIS_OPTIONS$spatial_nearest_fibroblast_bins, nrow(fib_xy))
  nn <- RANN::nn2(data = fib_xy, query = epi_xy, k = k)
  nn_index <- if (k == 1) matrix(nn$nn.idx, ncol = 1) else nn$nn.idx
  neighbour_ecm <- rowMeans(matrix(fib_scores[nn_index], nrow = nrow(nn_index)))
  observed <- suppressWarnings(stats::cor(neighbour_ecm, epi_scores, method = "spearman", use = "complete.obs"))
  permuted <- replicate(
    ANALYSIS_OPTIONS$spatial_permutations,
    {
      shuffled <- sample(fib_scores, replace = FALSE)
      local <- rowMeans(matrix(shuffled[nn_index], nrow = nrow(nn_index)))
      suppressWarnings(stats::cor(local, epi_scores, method = "spearman", use = "complete.obs"))
    }
  )
  p_emp <- (1 + sum(permuted >= observed, na.rm = TRUE)) / (1 + sum(is.finite(permuted)))
  section_summaries[[section]] <- data.frame(
    section = section,
    n_fibroblast_bins = nrow(fib_xy),
    n_epithelial_bins_analyzed = nrow(epi_xy),
    nearest_k = k,
    spearman_rho = observed,
    one_sided_permutation_p = p_emp,
    ecm_coverage = fib_score_obj$coverage,
    target_coverage = epi_score_obj$coverage,
    inference_note = "section-level exploratory spatial association; bins are not patients"
  )
  neighbour_rows[[section]] <- data.frame(
    section = section,
    x = epi_xy[, 1], y = epi_xy[, 2],
    epithelial_target_score = epi_scores,
    nearest_fibroblast_ECM_score = neighbour_ecm
  )
  plot_data <- neighbour_rows[[section]]
  if (nrow(plot_data) > 5000) plot_data <- plot_data[sample(seq_len(nrow(plot_data)), 5000), ]
  section_plots[[section]] <- ggplot2::ggplot(plot_data, ggplot2::aes(nearest_fibroblast_ECM_score, epithelial_target_score)) +
    ggplot2::geom_point(alpha = 0.18, size = 0.7, colour = "#365E96") +
    ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "#C9223A") +
    ggplot2::labs(title = section, subtitle = sprintf("Spearman rho = %.3f; permutation P = %.4f", observed, p_emp), x = "Nearest fibroblast ECM score", y = "Epithelial target score") +
    theme_publication()
  rm(counts, meta, fib_norm, epi_norm, nn, permuted)
  gc()
}

summary_df <- do.call(rbind, section_summaries)
neighbour_df <- do.call(rbind, neighbour_rows)
label_df <- do.call(rbind, label_audits)
write_csv(summary_df, file.path(out, "spatial_section_summary.csv"))
write_csv(neighbour_df, file.path(out, "spatial_neighbour_source_data.csv"))
write_csv(label_df, file.path(out, "spatial_label_audit.csv"))

# Fisher combination is shown only as an exploratory across-section summary.
combined_p <- stats::pchisq(-2 * sum(log(pmax(summary_df$one_sided_permutation_p, 1e-12))), df = 2 * nrow(summary_df), lower.tail = FALSE)
combined <- data.frame(
  n_sections = nrow(summary_df),
  positive_sections = sum(summary_df$spearman_rho > 0),
  median_rho = median(summary_df$spearman_rho),
  fisher_combined_permutation_p = combined_p,
  interpretation = "exploratory; sections are the reporting units and may not represent independent patients"
)
write_csv(combined, file.path(out, "spatial_combined_exploratory_summary.csv"))

p <- patchwork::wrap_plots(section_plots, nrow = 1)
save_plot(p, file.path(out, "GSE280315_spatial_neighbour_association"), 13, 4.2)
write_status("04_GSE280315_spatial", "PASS", paste("Sections:", paste(summary_df$section, collapse = ", ")))
