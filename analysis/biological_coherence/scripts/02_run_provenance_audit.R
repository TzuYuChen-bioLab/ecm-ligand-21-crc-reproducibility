#!/usr/bin/env Rscript

# Outcome-blind retrospective biological-coherence audit of ECM_CORE_19 and
# ECM_LIGAND_21. This script evaluates an unchanged panel; it does not recreate
# or claim the historical feature-selection process.

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

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
script_dir <- get_script_dir()
input_dir <- normalizePath(args$input_dir %||% file.path(script_dir, "input"), winslash = "/", mustWork = FALSE)
output_dir <- normalizePath(args$output_dir %||% file.path(script_dir, "results"), winslash = "/", mustWork = FALSE)
panel_file <- normalizePath(args$panel_file %||% file.path(script_dir, "panel_21_annotation.csv"), winslash = "/", mustWork = FALSE)
background_file <- args$background_file %||% file.path(input_dir, "naba_background_genes.csv")
n_random <- as.integer(args$n_random %||% "10000")
random_seed <- as.integer(args$seed %||% "20260808")

if (!dir.exists(input_dir)) stop("Input directory not found: ", input_dir)
if (!file.exists(panel_file)) stop("Panel annotation not found: ", panel_file)
if (!is.finite(n_random) || n_random < 1000) stop("n_random must be at least 1000; 10000 is recommended.")
if (!is.finite(random_seed)) stop("seed must be an integer.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(random_seed)

needed <- c("data.table", "edgeR", "ggplot2")
missing_packages <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Install missing package(s) before running: ", paste(missing_packages, collapse = ", "),
    "\nFor edgeR: if (!requireNamespace('BiocManager')) install.packages('BiocManager'); BiocManager::install('edgeR')"
  )
}

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\.[0-9]+$", "", x)
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

safe_cor <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3 || stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0) return(NA_real_)
  unname(stats::cor(x[keep], y[keep], method = method))
}

safe_t_test <- function(delta) {
  delta <- delta[is.finite(delta)]
  if (length(delta) < 2 || stats::sd(delta) == 0) return(c(p = NA_real_, low = NA_real_, high = NA_real_))
  ans <- tryCatch(stats::t.test(delta, mu = 0), error = function(e) NULL)
  if (is.null(ans)) return(c(p = NA_real_, low = NA_real_, high = NA_real_))
  c(p = ans$p.value, low = unname(ans$conf.int[1]), high = unname(ans$conf.int[2]))
}

safe_wilcox <- function(delta) {
  delta <- delta[is.finite(delta)]
  if (length(delta) < 2 || all(delta == 0)) return(NA_real_)
  tryCatch(
    suppressWarnings(stats::wilcox.test(delta, mu = 0, exact = length(delta) <= 50)$p.value),
    error = function(e) NA_real_
  )
}

summarise_delta <- function(delta, dataset, analysis, panel) {
  delta <- as.numeric(delta[is.finite(delta)])
  tt <- safe_t_test(delta)
  data.frame(
    dataset = dataset,
    analysis = analysis,
    panel = panel,
    n_pairs = length(delta),
    mean_difference = if (length(delta)) mean(delta) else NA_real_,
    ci_low_t = unname(tt["low"]),
    ci_high_t = unname(tt["high"]),
    median_difference = if (length(delta)) stats::median(delta) else NA_real_,
    paired_t_p = unname(tt["p"]),
    wilcoxon_p = safe_wilcox(delta),
    positive_pair_fraction = if (length(delta)) mean(delta > 0) else NA_real_,
    stringsAsFactors = FALSE
  )
}

row_or_na <- function(mat, gene) {
  if (gene %in% rownames(mat)) return(as.numeric(mat[gene, ]))
  rep(NA_real_, ncol(mat))
}

read_count_matrix <- function(path) {
  x <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(x) < 2) stop("Count file needs a gene column and sample columns: ", path)
  gene_col <- names(x)[1]
  genes <- clean_gene(x[[gene_col]])
  x[[gene_col]] <- NULL
  x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
  mat <- as.matrix(x)
  storage.mode(mat) <- "numeric"
  keep <- !is.na(genes)
  genes <- genes[keep]
  mat <- mat[keep, , drop = FALSE]
  if (any(!is.finite(mat)) || any(mat < 0)) stop("Counts must be finite and non-negative: ", path)
  if (max(abs(mat - round(mat))) > 1e-6) stop("Input must contain raw integer pseudobulk counts, not log-normalized expression: ", path)
  if (anyDuplicated(genes)) mat <- rowsum(mat, group = genes, reorder = FALSE)
  else rownames(mat) <- genes
  round(mat)
}

zscore_rows <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sigma <- apply(mat, 1, stats::sd, na.rm = TRUE)
  sigma[!is.finite(sigma) | sigma == 0] <- NA_real_
  sweep(sweep(mat, 1, mu, "-"), 1, sigma, "/")
}

paired_gene_delta <- function(expression, metadata, target_class = "Fibroblast") {
  use <- metadata$cell_class == target_class & metadata$condition %in% c("Tumor", "Normal")
  m <- metadata[use, , drop = FALSE]
  e <- expression[, use, drop = FALSE]
  patients <- intersect(m$patient_id[m$condition == "Tumor"], m$patient_id[m$condition == "Normal"])
  if (!length(patients)) return(matrix(numeric(0), nrow = nrow(expression), dimnames = list(rownames(expression), NULL)))
  delta_list <- lapply(patients, function(patient) {
    as.numeric(
      rowMeans(e[, m$patient_id == patient & m$condition == "Tumor", drop = FALSE]) -
        rowMeans(e[, m$patient_id == patient & m$condition == "Normal", drop = FALSE])
    )
  })
  delta <- matrix(
    unlist(delta_list, use.names = FALSE),
    nrow = nrow(expression),
    ncol = length(patients),
    dimnames = list(rownames(expression), patients)
  )
  delta
}

fibroblast_specificity_delta <- function(expression, metadata) {
  use <- metadata$condition == "Tumor"
  m <- metadata[use, , drop = FALSE]
  e <- expression[, use, drop = FALSE]
  patients <- intersect(m$patient_id[m$cell_class == "Fibroblast"], m$patient_id[m$cell_class != "Fibroblast"])
  if (!length(patients)) return(matrix(numeric(0), nrow = nrow(expression), dimnames = list(rownames(expression), NULL)))
  delta_list <- lapply(patients, function(patient) {
    as.numeric(
      rowMeans(e[, m$patient_id == patient & m$cell_class == "Fibroblast", drop = FALSE]) -
        rowMeans(e[, m$patient_id == patient & m$cell_class != "Fibroblast", drop = FALSE])
    )
  })
  delta <- matrix(
    unlist(delta_list, use.names = FALSE),
    nrow = nrow(expression),
    ncol = length(patients),
    dimnames = list(rownames(expression), patients)
  )
  delta
}

score_panel <- function(z_expression, genes, minimum_coverage = 0.80) {
  present <- intersect(genes, rownames(z_expression))
  coverage <- length(present) / length(genes)
  if (coverage < minimum_coverage) {
    stop(sprintf("Panel coverage %.1f%% is below the pre-specified %.1f%% threshold.", 100 * coverage, 100 * minimum_coverage))
  }
  list(score = colMeans(z_expression[present, , drop = FALSE], na.rm = TRUE), present = present, coverage = coverage)
}

panel <- data.table::fread(panel_file, data.table = FALSE)
required_panel_columns <- c("gene", "panel", "module", "naba_gene_set")
if (length(setdiff(required_panel_columns, names(panel)))) stop("Panel annotation is missing required columns.")
panel$gene <- clean_gene(panel$gene)
if (anyNA(panel$gene) || anyDuplicated(panel$gene)) stop("Panel genes must be unique and non-missing.")
core19 <- panel$gene[panel$panel == "ECM_CORE_19"]
ligand21 <- panel$gene
if (length(core19) != 19 || length(ligand21) != 21 || !all(c("PTN", "MDK") %in% setdiff(ligand21, core19))) {
  stop("Panel annotation must define 19 core genes and a 21-gene extension containing PTN and MDK.")
}

background_available <- file.exists(background_file)
if (background_available) {
  background <- data.table::fread(background_file, data.table = FALSE)
  if (length(setdiff(c("gene", "naba_gene_set"), names(background)))) {
    stop("Background file needs columns 'gene' and 'naba_gene_set': ", background_file)
  }
  background$gene <- clean_gene(background$gene)
  background <- unique(background[!is.na(background$gene), c("gene", "naba_gene_set")])
}

count_files <- list.files(input_dir, pattern = "_pseudobulk_counts\\.csv(\\.gz)?$", full.names = TRUE)
if (!length(count_files)) {
  stop(
    "No *_pseudobulk_counts.csv.gz files found in ", input_dir,
    ". Run 01_prepare_pseudobulk_from_seurat.R first or follow the templates."
  )
}

all_gene_panel <- list()
all_sample_scores <- list()
all_panel_effects <- list()
all_similarity <- list()
all_benchmark_summary <- list()
all_benchmark_null <- list()
all_matching_audit <- list()
all_input_audit <- list()

run_random_benchmark <- function(dataset, panel_name, panel_genes, gene_stats, background, background_type, B) {
  metric_columns <- c("fibroblast_specificity", "fibroblast_tumor_normal_effect")
  eligible <- gene_stats[is.finite(gene_stats[[metric_columns[1]]]) & is.finite(gene_stats[[metric_columns[2]]]), , drop = FALSE]
  if (background_type == "NABA_category_matched") {
    eligible <- merge(eligible, background, by = "gene", all = FALSE)
    eligible <- eligible[!duplicated(eligible$gene), , drop = FALSE]
  } else {
    eligible$naba_gene_set <- "TRANSCRIPTOME"
  }

  target <- merge(
    data.frame(gene = panel_genes, stringsAsFactors = FALSE),
    panel[, c("gene", "naba_gene_set")],
    by = "gene", all.x = TRUE, sort = FALSE
  )
  target <- target[match(panel_genes, target$gene), , drop = FALSE]
  if (background_type != "NABA_category_matched") target$naba_gene_set <- "TRANSCRIPTOME"
  target <- merge(target, gene_stats[, c("gene", "mean_logCPM", "detection_rate", metric_columns)], by = "gene", all.x = TRUE, sort = FALSE)
  target <- target[match(panel_genes, target$gene), , drop = FALSE]
  target <- target[is.finite(target[[metric_columns[1]]]) & is.finite(target[[metric_columns[2]]]), , drop = FALSE]
  if (nrow(target) < ceiling(0.80 * length(panel_genes))) stop("Fewer than 80% of ", panel_name, " genes have both audit metrics in ", dataset, ".")

  eligible <- eligible[!eligible$gene %in% ligand21, , drop = FALSE]
  if (nrow(eligible) < 10 * nrow(target)) stop("Too few eligible background genes for matched random panels in ", dataset, ".")

  match_frame <- rbind(
    data.frame(gene = eligible$gene, mean_logCPM = eligible$mean_logCPM, detection_rate = eligible$detection_rate),
    data.frame(gene = target$gene, mean_logCPM = target$mean_logCPM, detection_rate = target$detection_rate)
  )
  scaled <- scale(match_frame[, c("mean_logCPM", "detection_rate")])
  scaled[!is.finite(scaled)] <- 0
  rownames(scaled) <- match_frame$gene

  pools <- vector("list", nrow(target))
  pool_audit <- vector("list", nrow(target))
  for (i in seq_len(nrow(target))) {
    same_category <- eligible$naba_gene_set == target$naba_gene_set[i]
    candidates <- eligible$gene[same_category]
    if (!length(candidates)) stop("No measured background genes for category ", target$naba_gene_set[i], " in ", dataset, ".")
    distance <- sqrt(rowSums((scaled[candidates, , drop = FALSE] - matrix(scaled[target$gene[i], ], nrow = length(candidates), ncol = 2, byrow = TRUE))^2))
    candidates <- candidates[order(distance)]
    pools[[i]] <- head(candidates, min(200, length(candidates)))
    pool_audit[[i]] <- data.frame(
      dataset = dataset,
      panel = panel_name,
      target_gene = target$gene[i],
      naba_gene_set = target$naba_gene_set[i],
      eligible_same_category = length(candidates),
      nearest_pool_size = length(pools[[i]]),
      stringsAsFactors = FALSE
    )
  }

  metric_index <- match(eligible$gene, gene_stats$gene)
  names(metric_index) <- eligible$gene
  null_spec <- numeric(B)
  null_tn <- numeric(B)
  for (b in seq_len(B)) {
    chosen <- character(0)
    for (i in sample(seq_len(nrow(target)))) {
      candidates <- setdiff(pools[[i]], chosen)
      if (!length(candidates)) {
        candidates <- setdiff(eligible$gene[eligible$naba_gene_set == target$naba_gene_set[i]], chosen)
      }
      if (!length(candidates)) stop("Unable to sample a unique category-matched panel.")
      chosen <- c(chosen, sample(candidates, 1))
    }
    idx <- unname(metric_index[chosen])
    null_spec[b] <- mean(gene_stats$fibroblast_specificity[idx])
    null_tn[b] <- mean(gene_stats$fibroblast_tumor_normal_effect[idx])
  }

  observed_spec <- mean(target$fibroblast_specificity)
  observed_tn <- mean(target$fibroblast_tumor_normal_effect)
  spec_sd <- stats::sd(null_spec)
  tn_sd <- stats::sd(null_tn)
  if (!is.finite(spec_sd) || spec_sd == 0 || !is.finite(tn_sd) || tn_sd == 0) {
    stop("Matched-random null distribution has zero/undefined variance in ", dataset, " for ", panel_name, ".")
  }
  combined_null <- (null_spec - mean(null_spec)) / spec_sd + (null_tn - mean(null_tn)) / tn_sd
  combined_observed <- (observed_spec - mean(null_spec)) / spec_sd + (observed_tn - mean(null_tn)) / tn_sd
  add_one_p <- function(null, observed) (1 + sum(null >= observed)) / (1 + length(null))

  summary <- data.frame(
    dataset = dataset,
    panel = panel_name,
    background_type = background_type,
    panel_genes_requested = length(panel_genes),
    panel_genes_tested = nrow(target),
    n_random = B,
    observed_fibroblast_specificity = observed_spec,
    null_mean_fibroblast_specificity = mean(null_spec),
    empirical_p_fibroblast_specificity = add_one_p(null_spec, observed_spec),
    observed_fibroblast_tumor_normal_effect = observed_tn,
    null_mean_fibroblast_tumor_normal_effect = mean(null_tn),
    empirical_p_fibroblast_tumor_normal_effect = add_one_p(null_tn, observed_tn),
    observed_combined_index = combined_observed,
    empirical_p_combined_index = add_one_p(combined_null, combined_observed),
    stringsAsFactors = FALSE
  )
  null <- data.frame(
    dataset = dataset,
    panel = panel_name,
    iteration = seq_len(B),
    fibroblast_specificity = null_spec,
    fibroblast_tumor_normal_effect = null_tn,
    combined_index = combined_null,
    observed_combined_index = combined_observed,
    stringsAsFactors = FALSE
  )
  list(summary = summary, null = null, matching = do.call(rbind, pool_audit))
}

for (count_path in count_files) {
  dataset_from_file <- sub("_pseudobulk_counts\\.csv(\\.gz)?$", "", basename(count_path))
  metadata_path <- file.path(input_dir, paste0(dataset_from_file, "_pseudobulk_metadata.csv"))
  if (!file.exists(metadata_path)) stop("Matching metadata file not found: ", metadata_path)

  counts <- read_count_matrix(count_path)
  metadata <- data.table::fread(metadata_path, data.table = FALSE)
  required_meta <- c("sample_id", "dataset", "patient_id", "condition", "cell_class", "n_cells", "paired_id")
  if (length(setdiff(required_meta, names(metadata)))) stop("Metadata missing required columns: ", metadata_path)
  if (anyDuplicated(metadata$sample_id)) stop("Duplicate sample_id in ", metadata_path)
  if (!all(colnames(counts) %in% metadata$sample_id)) stop("Some count columns are absent from metadata: ", count_path)
  metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
  if (anyNA(metadata$sample_id)) stop("Failed to align counts and metadata for ", dataset_from_file)
  dataset <- unique(metadata$dataset)
  if (length(dataset) != 1) stop("Each input file pair must contain exactly one dataset.")
  dataset <- as.character(dataset)

  if (!all(metadata$condition %in% c("Tumor", "Normal", "Border"))) stop("condition must be Tumor, Normal, or Border.")
  if (sum(metadata$cell_class == "Fibroblast" & metadata$condition == "Tumor") < 3 ||
      sum(metadata$cell_class == "Fibroblast" & metadata$condition == "Normal") < 3) {
    stop(dataset, " needs at least three Tumor and three Normal fibroblast pseudobulks.")
  }

  y <- edgeR::DGEList(counts = counts)
  raw_cpm <- edgeR::cpm(y, log = FALSE)
  minimum_samples <- min(ncol(counts), max(3L, ceiling(0.20 * ncol(counts))))
  keep_gene <- rowSums(raw_cpm >= 1) >= minimum_samples
  if (sum(keep_gene) < 1000) stop("Fewer than 1,000 genes pass CPM >= 1 in at least ", minimum_samples, " pseudobulks for ", dataset, ".")
  y <- y[keep_gene, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  normalized_cpm <- edgeR::cpm(y, log = FALSE)
  log_cpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  z_expression <- zscore_rows(log_cpm)

  core_score <- score_panel(z_expression, core19)
  ligand_score <- score_panel(z_expression, ligand21)
  sample_scores <- cbind(
    metadata,
    ECM_CORE_19 = as.numeric(core_score$score),
    ECM_LIGAND_21 = as.numeric(ligand_score$score),
    extension_contribution = as.numeric(ligand_score$score - core_score$score)
  )
  all_sample_scores[[dataset]] <- sample_scores

  tn_delta <- paired_gene_delta(log_cpm, metadata, "Fibroblast")
  specificity_delta <- fibroblast_specificity_delta(log_cpm, metadata)
  if (ncol(tn_delta) < 3) stop(dataset, " has fewer than three complete Tumor-Normal fibroblast pairs.")
  if (ncol(specificity_delta) < 3) stop(dataset, " has fewer than three patients with both Tumor fibroblast and Tumor non-fibroblast pseudobulks.")

  gene_stats <- data.frame(
    gene = rownames(log_cpm),
    mean_logCPM = rowMeans(log_cpm),
    detection_rate = rowMeans(normalized_cpm >= 1),
    fibroblast_specificity = rowMeans(specificity_delta),
    n_specificity_patients = ncol(specificity_delta),
    fibroblast_tumor_normal_effect = rowMeans(tn_delta),
    n_tumor_normal_pairs = ncol(tn_delta),
    stringsAsFactors = FALSE
  )

  panel_gene_stats <- merge(panel, gene_stats, by = "gene", all.x = TRUE, sort = FALSE)
  panel_gene_stats <- panel_gene_stats[match(panel$gene, panel_gene_stats$gene), , drop = FALSE]
  panel_gene_stats$dataset <- dataset
  panel_gene_stats$specificity_t_p <- vapply(panel_gene_stats$gene, function(g) {
    if (!g %in% rownames(specificity_delta)) return(NA_real_)
    unname(safe_t_test(specificity_delta[g, ])["p"])
  }, numeric(1))
  panel_gene_stats$tumor_normal_t_p <- vapply(panel_gene_stats$gene, function(g) {
    if (!g %in% rownames(tn_delta)) return(NA_real_)
    unname(safe_t_test(tn_delta[g, ])["p"])
  }, numeric(1))
  all_gene_panel[[dataset]] <- panel_gene_stats

  score_tn_core <- paired_gene_delta(matrix(core_score$score, nrow = 1, dimnames = list("ECM_CORE_19", names(core_score$score))), metadata)
  score_tn_ligand <- paired_gene_delta(matrix(ligand_score$score, nrow = 1, dimnames = list("ECM_LIGAND_21", names(ligand_score$score))), metadata)
  all_panel_effects[[dataset]] <- rbind(
    summarise_delta(score_tn_core[1, ], dataset, "Tumor_minus_Normal_fibroblast_score", "ECM_CORE_19"),
    summarise_delta(score_tn_ligand[1, ], dataset, "Tumor_minus_Normal_fibroblast_score", "ECM_LIGAND_21"),
    summarise_delta(row_or_na(tn_delta, "PTN"), dataset, "Tumor_minus_Normal_fibroblast_logCPM", "PTN"),
    summarise_delta(row_or_na(tn_delta, "MDK"), dataset, "Tumor_minus_Normal_fibroblast_logCPM", "MDK")
  )

  fibroblast_samples <- metadata$cell_class == "Fibroblast"
  all_similarity[[dataset]] <- data.frame(
    dataset = dataset,
    n_fibroblast_pseudobulks = sum(fibroblast_samples),
    pearson_core19_vs_ligand21 = safe_cor(core_score$score[fibroblast_samples], ligand_score$score[fibroblast_samples], "pearson"),
    spearman_core19_vs_ligand21 = safe_cor(core_score$score[fibroblast_samples], ligand_score$score[fibroblast_samples], "spearman"),
    mean_extension_contribution = mean(sample_scores$extension_contribution[fibroblast_samples]),
    sd_extension_contribution = stats::sd(sample_scores$extension_contribution[fibroblast_samples]),
    maximum_absolute_extension_contribution = max(abs(sample_scores$extension_contribution[fibroblast_samples])),
    stringsAsFactors = FALSE
  )

  if (background_available) {
    background_for_dataset <- background
    background_type <- "NABA_category_matched"
  } else {
    background_for_dataset <- data.frame(gene = gene_stats$gene, naba_gene_set = "TRANSCRIPTOME")
    background_type <- "expression_detection_matched_transcriptome"
  }

  bench_core <- run_random_benchmark(dataset, "ECM_CORE_19", core19, gene_stats, background_for_dataset, background_type, n_random)
  bench_ligand <- run_random_benchmark(dataset, "ECM_LIGAND_21", ligand21, gene_stats, background_for_dataset, background_type, n_random)
  all_benchmark_summary[[dataset]] <- rbind(bench_core$summary, bench_ligand$summary)
  all_benchmark_null[[dataset]] <- rbind(bench_core$null, bench_ligand$null)
  all_matching_audit[[dataset]] <- rbind(bench_core$matching, bench_ligand$matching)

  all_input_audit[[dataset]] <- data.frame(
    dataset = dataset,
    count_file = basename(count_path),
    metadata_file = basename(metadata_path),
    raw_genes = nrow(counts),
    retained_genes = nrow(log_cpm),
    pseudobulk_samples = ncol(counts),
    tumor_normal_fibroblast_pairs = ncol(tn_delta),
    tumor_fibroblast_specificity_patients = ncol(specificity_delta),
    core19_coverage = core_score$coverage,
    ligand21_coverage = ligand_score$coverage,
    background_type = background_type,
    stringsAsFactors = FALSE
  )
}

gene_panel_results <- do.call(rbind, all_gene_panel)
sample_score_results <- do.call(rbind, all_sample_scores)
panel_effect_results <- do.call(rbind, all_panel_effects)
similarity_results <- do.call(rbind, all_similarity)
benchmark_summary <- do.call(rbind, all_benchmark_summary)
benchmark_null <- do.call(rbind, all_benchmark_null)
matching_audit <- do.call(rbind, all_matching_audit)
input_audit <- do.call(rbind, all_input_audit)

data.table::fwrite(gene_panel_results, file.path(output_dir, "Table_gene_level_panel_audit.csv"), na = "NA")
data.table::fwrite(sample_score_results, file.path(output_dir, "Table_sample_level_core19_ligand21_scores.csv"), na = "NA")
data.table::fwrite(panel_effect_results, file.path(output_dir, "Table_panel_and_extension_effects.csv"), na = "NA")
data.table::fwrite(similarity_results, file.path(output_dir, "Table_core19_vs_ligand21_similarity.csv"), na = "NA")
data.table::fwrite(benchmark_summary, file.path(output_dir, "Table_biological_coherence_random_benchmark.csv"), na = "NA")
data.table::fwrite(benchmark_null, file.path(output_dir, "Random_benchmark_null_distributions.csv.gz"), na = "NA")
data.table::fwrite(matching_audit, file.path(output_dir, "Random_matching_pool_audit.csv"), na = "NA")
data.table::fwrite(input_audit, file.path(output_dir, "Input_and_coverage_audit.csv"), na = "NA")

module_coverage <- aggregate(gene ~ module, panel, length)
names(module_coverage)[2] <- "n_genes_ECM_LIGAND_21"
module_coverage$n_genes_ECM_CORE_19 <- vapply(module_coverage$module, function(z) sum(panel$module == z & panel$gene %in% core19), integer(1))
module_coverage$interpretation <- ifelse(
  module_coverage$module == "ECM_linked_secreted_factor",
  "Added only by PTN/MDK; functional-scope extension, not a demonstrated performance gain",
  "Represented in the structural 19-gene core"
)
data.table::fwrite(module_coverage, file.path(output_dir, "Table_functional_module_coverage.csv"), na = "NA")

paired_plot_data <- sample_score_results[
  sample_score_results$cell_class == "Fibroblast" & sample_score_results$condition %in% c("Normal", "Tumor"),
  c("dataset", "patient_id", "condition", "ECM_CORE_19", "ECM_LIGAND_21")
]
paired_plot_data <- reshape(
  paired_plot_data,
  varying = c("ECM_CORE_19", "ECM_LIGAND_21"),
  v.names = "score", timevar = "panel", times = c("ECM_CORE_19", "ECM_LIGAND_21"),
  direction = "long"
)
paired_plot_data$condition <- factor(paired_plot_data$condition, levels = c("Normal", "Tumor"))
p1 <- ggplot2::ggplot(
  paired_plot_data,
  ggplot2::aes(condition, score, group = interaction(dataset, patient_id), colour = condition)
) +
  ggplot2::geom_line(colour = "grey65", linewidth = 0.45) +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::facet_grid(panel ~ dataset, scales = "free_y") +
  ggplot2::scale_colour_manual(values = c(Normal = "#2C73B9", Tumor = "#C9223A")) +
  ggplot2::labs(x = NULL, y = "Mean gene-wise z score", title = "Outcome-blind 19-vs-21 paired fibroblast audit") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(file.path(output_dir, "Figure_core19_vs_ligand21_paired.pdf"), p1, width = 8.2, height = 5.8)
ggplot2::ggsave(file.path(output_dir, "Figure_core19_vs_ligand21_paired.png"), p1, width = 8.2, height = 5.8, dpi = 300)

p2 <- ggplot2::ggplot(
  gene_panel_results,
  ggplot2::aes(fibroblast_specificity, fibroblast_tumor_normal_effect, colour = module, label = gene)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
  ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_text(size = 2.3, vjust = -0.55, check_overlap = TRUE, show.legend = FALSE) +
  ggplot2::facet_wrap(~ dataset, scales = "free") +
  ggplot2::labs(
    x = "Tumor fibroblast minus non-fibroblast mean logCPM",
    y = "Fibroblast tumor minus normal mean logCPM",
    title = "Gene-level biological-coherence audit"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
ggplot2::ggsave(file.path(output_dir, "Figure_gene_level_coherence_map.pdf"), p2, width = 9.0, height = 5.2)
ggplot2::ggsave(file.path(output_dir, "Figure_gene_level_coherence_map.png"), p2, width = 9.0, height = 5.2, dpi = 300)

p3 <- ggplot2::ggplot(benchmark_null, ggplot2::aes(combined_index)) +
  ggplot2::geom_histogram(bins = 45, fill = "grey75", colour = "white") +
  ggplot2::geom_vline(ggplot2::aes(xintercept = observed_combined_index), colour = "#C9223A", linewidth = 0.8) +
  ggplot2::facet_grid(panel ~ dataset, scales = "free_y") +
  ggplot2::labs(x = "Empirical combined coherence index", y = "Random panels", title = "Expression/detection-matched random-panel benchmark") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(file.path(output_dir, "Figure_random_panel_benchmark.pdf"), p3, width = 8.2, height = 5.8)
ggplot2::ggsave(file.path(output_dir, "Figure_random_panel_benchmark.png"), p3, width = 8.2, height = 5.8, dpi = 300)

fmt <- function(x, digits = 3) ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
methods_text <- paste(
  "Retrospective biological-coherence audit (manuscript template)",
  "",
  "ECM_LIGAND_21 was evaluated unchanged and without access to recurrence outcomes. Raw UMI counts were summed within patient, condition, and major cell class. Pseudobulk libraries were TMM-normalized, and log2 counts per million were used for gene-level comparisons. Two outcome-independent quantities were calculated within each cohort: (i) fibroblast specificity, defined as tumor-fibroblast expression minus the mean expression of non-fibroblast tumor pseudobulks from the same patient; and (ii) the paired tumor-minus-normal effect within fibroblast pseudobulks. ECM_CORE_19 and ECM_LIGAND_21 scores were calculated as equal-weight means of gene-wise standardized expression. The unchanged panels were compared with random panels matched on expression and detection rate; when the NABA background file was supplied, matching was additionally constrained by Matrisome category. This retrospective audit evaluated biological coherence and the effect of adding PTN/MDK; it was not used to select, remove, or reweight genes and does not reconstruct the original historical inclusion rule.",
  sep = "\n"
)

result_lines <- c("Outcome-blind audit results (insert values only after review)", "")
for (dataset in unique(similarity_results$dataset)) {
  sim <- similarity_results[similarity_results$dataset == dataset, ][1, ]
  eff19 <- panel_effect_results[panel_effect_results$dataset == dataset & panel_effect_results$panel == "ECM_CORE_19", ][1, ]
  eff21 <- panel_effect_results[panel_effect_results$dataset == dataset & panel_effect_results$panel == "ECM_LIGAND_21", ][1, ]
  b21 <- benchmark_summary[benchmark_summary$dataset == dataset & benchmark_summary$panel == "ECM_LIGAND_21", ][1, ]
  result_lines <- c(
    result_lines,
    paste0(
      dataset, ": ECM_CORE_19 and ECM_LIGAND_21 scores were highly similar across fibroblast pseudobulks (Pearson r = ",
      fmt(sim$pearson_core19_vs_ligand21), "; Spearman rho = ", fmt(sim$spearman_core19_vs_ligand21), "). "
    ),
    paste0(
      "The paired tumor-minus-normal score difference was ", fmt(eff19$mean_difference), " for ECM_CORE_19 and ",
      fmt(eff21$mean_difference), " for ECM_LIGAND_21 (", eff21$n_pairs, " patient pairs). "
    ),
    paste0(
      "For the 21-gene panel, empirical one-sided P values were ", fmt(b21$empirical_p_fibroblast_specificity, 4),
      " for fibroblast specificity, ", fmt(b21$empirical_p_fibroblast_tumor_normal_effect, 4),
      " for the fibroblast tumor-normal effect, and ", fmt(b21$empirical_p_combined_index, 4),
      " for the combined outcome-independent coherence index using ", b21$background_type, "."
    ),
    "These results may support biological coherence or scope coverage, but they do not establish that 21 is an optimal size, that PTN/MDK are necessary, or that the historical panel was generated by this audit.",
    ""
  )
}

writeLines(methods_text, file.path(output_dir, "Manuscript_methods_template.txt"), useBytes = TRUE)
writeLines(result_lines, file.path(output_dir, "Manuscript_results_template.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"), useBytes = TRUE)

manifest_paths <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  file = basename(manifest_paths),
  bytes = file.info(manifest_paths)$size,
  md5 = unname(tools::md5sum(manifest_paths)),
  stringsAsFactors = FALSE
)
data.table::fwrite(manifest, file.path(output_dir, "OUTPUT_MANIFEST.csv"), na = "NA")

message("Audit complete. Results written to: ", output_dir)
