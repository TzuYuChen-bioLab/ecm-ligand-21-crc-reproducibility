PIPELINE_ROOT <- normalizePath(
  Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))

require_packages(c("data.table", "ggplot2"))

out <- result_dir("07_GSE155343")

if (
  !ANALYSIS_OPTIONS$run_optional_GSE155343 ||
    !file.exists(INPUT_FILES$GSE155343_expression)
) {
  write_status(
    "07_GSE155343",
    "SKIPPED",
    "Optional file absent or analysis disabled."
  )
} else {
  # --------------------------------------------------------------------
  # Read the deposited processed matrix without discarding its gene_name
  # annotation. The GEO file contains four annotation columns followed by
  # Sample 1--Sample 25. The generic numeric_expression_matrix() helper is
  # intentionally not used here because its first-column-ID convention is
  # incompatible with this multi-annotation-column file.
  # --------------------------------------------------------------------
  deposited <- data.table::fread(
    INPUT_FILES$GSE155343_expression,
    check.names = FALSE,
    data.table = FALSE,
    showProgress = FALSE
  )

  normalise_header <- function(x) {
    tolower(gsub("[^[:alnum:]]+", "", x))
  }

  header_key <- normalise_header(names(deposited))

  symbol_column <- which(
    header_key %in% c(
      "genename",
      "genesymbol",
      "symbol",
      "hgncsymbol"
    )
  )

  ensembl_column <- which(
    header_key %in% c(
      "ensemblegeneid",
      "ensemblgeneid",
      "ensemblid"
    )
  )

  chromosome_column <- which(header_key == "chromosome")
  gene_type_column <- which(header_key == "genetype")
  sample_columns <- grep(
    "^sample[[:space:]]*[0-9]+$",
    names(deposited),
    ignore.case = TRUE
  )

  if (length(symbol_column) != 1L) {
    stop(
      "Expected exactly one gene_name/gene-symbol column in the deposited ",
      "GSE155343 matrix; detected ", length(symbol_column), "."
    )
  }

  if (length(sample_columns) != 25L) {
    stop(
      "Expected Sample 1--Sample 25 in the deposited GSE155343 matrix; ",
      "detected ", length(sample_columns), " sample columns."
    )
  }

  sample_number <- suppressWarnings(
    as.integer(
      sub(
        "^sample[[:space:]]*([0-9]+)$",
        "\\1",
        names(deposited)[sample_columns],
        ignore.case = TRUE
      )
    )
  )

  if (
    anyNA(sample_number) ||
      !identical(sort(sample_number), seq_len(25L))
  ) {
    stop(
      "The GSE155343 sample columns do not form an unambiguous Sample 1--25 ",
      "sequence. No order-only fallback was used."
    )
  }

  expression_raw <- as.matrix(
    deposited[, sample_columns, drop = FALSE]
  )
  suppressWarnings(storage.mode(expression_raw) <- "double")

  if (any(!is.finite(expression_raw))) {
    stop(
      "Non-finite values were introduced while reading the 25 rLog sample ",
      "columns. Check the deposited file integrity."
    )
  }

  gene_symbol_raw <- trimws(
    as.character(deposited[[symbol_column]])
  )
  gene_symbol <- toupper(gene_symbol_raw)
  valid_symbol <- (
    !is.na(gene_symbol) &
      nzchar(gene_symbol) &
      gene_symbol != "NA"
  )

  if (!any(valid_symbol)) {
    stop("No usable gene symbols were found in the gene_name column.")
  }

  expression_raw <- expression_raw[valid_symbol, , drop = FALSE]
  gene_symbol <- gene_symbol[valid_symbol]

  # Average rLog values when more than one deposited Ensembl row has the
  # same gene symbol. This is a descriptive transformed-expression analysis,
  # not count-level differential expression.
  symbol_levels <- unique(gene_symbol)
  symbol_index <- match(gene_symbol, symbol_levels)
  expression_sum <- rowsum(
    expression_raw,
    group = symbol_index,
    reorder = FALSE
  )
  expression <- sweep(
    expression_sum,
    MARGIN = 1L,
    STATS = tabulate(symbol_index, nbins = length(symbol_levels)),
    FUN = "/"
  )
  rownames(expression) <- symbol_levels
  colnames(expression) <- names(deposited)[sample_columns]

  signature <- unique(
    toupper(trimws(as.character(read_gene_set(SIGNATURE_FILE))))
  )
  target <- unique(
    toupper(trimws(as.character(read_gene_set(TARGET_FILE))))
  )
  signature <- signature[!is.na(signature) & nzchar(signature)]
  target <- target[!is.na(target) & nzchar(target)]

  if (length(signature) != 21L) {
    stop(
      "The locked ECM_LIGAND_21 set must contain exactly 21 unique genes; ",
      "detected ", length(signature), "."
    )
  }
  if (length(target) != 50L) {
    stop(
      "The locked epithelial target set must contain exactly 50 unique ",
      "genes; detected ", length(target), "."
    )
  }

  make_coverage_audit <- function(genes, gene_set, minimum_required) {
    detected <- genes %in% rownames(expression)
    coverage <- mean(detected)
    data.frame(
      gene_set = gene_set,
      gene = genes,
      detected = detected,
      overall_coverage = coverage,
      minimum_required = minimum_required,
      passes_prespecified_gate = coverage >= minimum_required,
      stringsAsFactors = FALSE
    )
  }

  signature_coverage <- make_coverage_audit(
    signature,
    "ECM_LIGAND_21",
    ANALYSIS_OPTIONS$minimum_signature_coverage
  )
  target_coverage <- make_coverage_audit(
    target,
    "EPITHELIAL_TUMOR_UP_50",
    ANALYSIS_OPTIONS$minimum_target_coverage
  )

  write_csv(
    signature_coverage,
    file.path(out, "ECM_LIGAND_21_coverage.csv")
  )
  write_csv(
    target_coverage,
    file.path(out, "EPITHELIAL_TUMOR_UP_50_coverage.csv")
  )

  signature_fraction <- unique(signature_coverage$overall_coverage)
  target_fraction <- unique(target_coverage$overall_coverage)

  if (
    signature_fraction <
      ANALYSIS_OPTIONS$minimum_signature_coverage
  ) {
    stop(
      sprintf(
        paste0(
          "Gene-symbol-resolved ECM_LIGAND_21 coverage %.1f%% is below ",
          "the pre-specified %.1f%% threshold. See ",
          "ECM_LIGAND_21_coverage.csv."
        ),
        100 * signature_fraction,
        100 * ANALYSIS_OPTIONS$minimum_signature_coverage
      )
    )
  }
  if (
    target_fraction <
      ANALYSIS_OPTIONS$minimum_target_coverage
  ) {
    stop(
      sprintf(
        paste0(
          "Gene-symbol-resolved epithelial-target coverage %.1f%% is below ",
          "the pre-specified %.1f%% threshold. See ",
          "EPITHELIAL_TUMOR_UP_50_coverage.csv."
        ),
        100 * target_fraction,
        100 * ANALYSIS_OPTIONS$minimum_target_coverage
      )
    )
  }

  ensembl_raw <- if (length(ensembl_column) == 1L) {
    as.character(deposited[[ensembl_column]])
  } else {
    rep(NA_character_, nrow(deposited))
  }
  chromosome_raw <- if (length(chromosome_column) == 1L) {
    as.character(deposited[[chromosome_column]])
  } else {
    rep(NA_character_, nrow(deposited))
  }
  gene_type_raw <- if (length(gene_type_column) == 1L) {
    as.character(deposited[[gene_type_column]])
  } else {
    rep(NA_character_, nrow(deposited))
  }

  feature_audit <- data.frame(
    matrix_row = seq_len(nrow(deposited)),
    chromosome = chromosome_raw,
    gene_type = gene_type_raw,
    gene_name_deposited = gene_symbol_raw,
    gene_symbol_used = ifelse(valid_symbol, toupper(gene_symbol_raw), NA),
    ensembl_gene_id = ensembl_raw,
    ensembl_without_version = sub(
      "[.][0-9]+$",
      "",
      ensembl_raw
    ),
    retained_after_symbol_check = valid_symbol,
    n_deposited_rows_for_symbol = ifelse(
      valid_symbol,
      as.integer(table(toupper(gene_symbol_raw[valid_symbol]))[
        toupper(gene_symbol_raw)
      ]),
      NA_integer_
    ),
    in_ECM_LIGAND_21 = valid_symbol &
      toupper(gene_symbol_raw) %in% signature,
    in_EPITHELIAL_TUMOR_UP_50 = valid_symbol &
      toupper(gene_symbol_raw) %in% target,
    stringsAsFactors = FALSE
  )

  write_csv(
    feature_audit,
    file.path(out, "matrix_feature_identifier_audit.csv")
  )

  provenance <- data.frame(
    field = c(
      "deposited_file",
      "deposited_value_type",
      "deposited_annotation_columns",
      "identifier_used_for_scoring",
      "duplicate_symbol_rule",
      "technical_replicate_rule",
      "inferential_scope",
      "ECM_LIGAND_21_detected",
      "ECM_LIGAND_21_total",
      "ECM_LIGAND_21_coverage",
      "EPITHELIAL_TUMOR_UP_50_detected",
      "EPITHELIAL_TUMOR_UP_50_total",
      "EPITHELIAL_TUMOR_UP_50_coverage"
    ),
    value = c(
      basename(INPUT_FILES$GSE155343_expression),
      "DESeq2 rLog transformed expression; not raw counts",
      paste(names(deposited)[-sample_columns], collapse = ";"),
      names(deposited)[symbol_column],
      "mean rLog value across deposited Ensembl rows sharing one gene symbol",
      "mean within each biological group after retaining library-level audit",
      "descriptive only; no inferential P values from technical libraries",
      sum(signature_coverage$detected),
      nrow(signature_coverage),
      signature_fraction,
      sum(target_coverage$detected),
      nrow(target_coverage),
      target_fraction
    ),
    stringsAsFactors = FALSE
  )
  write_csv(
    provenance,
    file.path(out, "identifier_mapping_provenance.csv")
  )

  # --------------------------------------------------------------------
  # Explicit Sample N -> official GEO series-order mapping.
  # All deposited replicates are technical and therefore are collapsed
  # before any biological interpretation.
  # --------------------------------------------------------------------
  official <- data.frame(
    sample_number = seq_len(25L),
    gsm = paste0("GSM", 4699516:4699540),
    biological_group = c(
      rep("HT29_mono", 3),
      rep("SW480_mono", 3),
      rep("NF_mono", 3),
      rep("HT29_after_coculture", 3),
      rep("NF_after_HT29", 3),
      rep("SW480_after_coculture", 3),
      rep("NF_after_SW480", 3),
      rep("CAF1", 2),
      rep("CAF2", 2)
    ),
    compartment = c(
      rep("epithelial", 6),
      rep("fibroblast", 3),
      rep("epithelial", 3),
      rep("fibroblast", 3),
      rep("epithelial", 3),
      rep("fibroblast", 7)
    ),
    replicate_type = "technical",
    stringsAsFactors = FALSE
  )

  idx <- match(sample_number, official$sample_number)
  mapping_is_valid <- (
    length(idx) == 25L &&
      all(!is.na(idx)) &&
      length(unique(idx)) == 25L
  )
  if (!mapping_is_valid) {
    manual_mapping_stop(
      data.frame(
        matrix_column = colnames(expression),
        parsed_sample_number = sample_number,
        stringsAsFactors = FALSE
      ),
      file.path(out, "sample_mapping_NEEDS_REVIEW.csv"),
      "Could not map Sample 1--25 to the 25 official GSE155343 libraries."
    )
  }

  meta <- official[idx, , drop = FALSE]
  meta$matrix_column <- colnames(expression)
  meta$mapping_rule <- paste(
    "deposited_Sample_number_to_official_GEO_series_order"
  )
  meta <- meta[
    ,
    c(
      "matrix_column",
      "sample_number",
      "gsm",
      "biological_group",
      "compartment",
      "replicate_type",
      "mapping_rule"
    )
  ]
  write_csv(meta, file.path(out, "sample_mapping_audit.csv"))

  ecm <- score_gene_set(
    expression,
    signature,
    ANALYSIS_OPTIONS$minimum_signature_coverage
  )
  epithelial <- score_gene_set(
    expression,
    target,
    ANALYSIS_OPTIONS$minimum_target_coverage
  )

  libraries <- cbind(
    meta,
    ECM_LIGAND_21 = as.numeric(ecm$score),
    ECM_LIGAND_21_coverage = as.numeric(ecm$coverage),
    epithelial_target = as.numeric(epithelial$score),
    epithelial_target_coverage = as.numeric(epithelial$coverage)
  )
  libraries$inference_note <- paste(
    "deposited libraries are technical replicates and are not independent",
    "biological units"
  )
  write_csv(
    libraries,
    file.path(out, "technical_library_scores.csv")
  )

  summarise_group <- function(x) {
    data.frame(
      biological_group = x$biological_group[1],
      compartment = x$compartment[1],
      n_technical_libraries = nrow(x),
      ECM_LIGAND_21_mean = mean(x$ECM_LIGAND_21),
      ECM_LIGAND_21_sd = stats::sd(x$ECM_LIGAND_21),
      epithelial_target_mean = mean(x$epithelial_target),
      epithelial_target_sd = stats::sd(x$epithelial_target),
      stringsAsFactors = FALSE
    )
  }

  collapsed <- do.call(
    rbind,
    lapply(
      split(libraries, libraries$biological_group),
      summarise_group
    )
  )
  rownames(collapsed) <- NULL
  # Preserve the original Step 07 column names for downstream compatibility.
  collapsed$ECM_LIGAND_21 <- collapsed$ECM_LIGAND_21_mean
  collapsed$epithelial_target <- collapsed$epithelial_target_mean
  collapsed$inference_note <- paste(
    "descriptive means after collapsing technical replicates;",
    "no inferential P values"
  )
  write_csv(
    collapsed,
    file.path(out, "collapsed_descriptive_scores.csv")
  )

  contrast_spec <- data.frame(
    contrast = c(
      "NF_after_HT29 minus NF_mono",
      "NF_after_SW480 minus NF_mono",
      "CAF1 minus NF_mono",
      "CAF2 minus NF_mono",
      "HT29_after_coculture minus HT29_mono",
      "SW480_after_coculture minus SW480_mono"
    ),
    program = c(
      rep("ECM_LIGAND_21", 4),
      rep("epithelial_target", 2)
    ),
    comparison_group = c(
      "NF_after_HT29",
      "NF_after_SW480",
      "CAF1",
      "CAF2",
      "HT29_after_coculture",
      "SW480_after_coculture"
    ),
    reference_group = c(
      rep("NF_mono", 4),
      "HT29_mono",
      "SW480_mono"
    ),
    biological_scope = c(
      rep("fibroblast ECM-program response", 4),
      rep("CRC-cell epithelial tumour-up-program response", 2)
    ),
    stringsAsFactors = FALSE
  )

  group_mean <- function(group, program) {
    column <- paste0(program, "_mean")
    hit <- collapsed[collapsed$biological_group == group, column]
    if (length(hit) != 1L || !is.finite(hit)) {
      stop("Could not obtain one collapsed value for ", group, " / ", program)
    }
    as.numeric(hit)
  }

  contrast_spec$comparison_mean <- mapply(
    group_mean,
    contrast_spec$comparison_group,
    contrast_spec$program
  )
  contrast_spec$reference_mean <- mapply(
    group_mean,
    contrast_spec$reference_group,
    contrast_spec$program
  )
  contrast_spec$descriptive_difference <- (
    contrast_spec$comparison_mean - contrast_spec$reference_mean
  )
  contrast_spec$inference_note <- paste(
    "descriptive difference between technical-replicate-collapsed means;",
    "not an inferential effect estimate"
  )
  write_csv(
    contrast_spec,
    file.path(out, "prespecified_descriptive_contrasts.csv")
  )

  fibroblast_groups <- c(
    "NF_mono",
    "NF_after_HT29",
    "NF_after_SW480",
    "CAF1",
    "CAF2"
  )
  epithelial_groups <- c(
    "HT29_mono",
    "HT29_after_coculture",
    "SW480_mono",
    "SW480_after_coculture"
  )

  plot_data <- rbind(
    data.frame(
      biological_group = collapsed$biological_group[
        match(fibroblast_groups, collapsed$biological_group)
      ],
      program = "Fibroblast: ECM_LIGAND_21",
      score = collapsed$ECM_LIGAND_21_mean[
        match(fibroblast_groups, collapsed$biological_group)
      ],
      stringsAsFactors = FALSE
    ),
    data.frame(
      biological_group = collapsed$biological_group[
        match(epithelial_groups, collapsed$biological_group)
      ],
      program = "CRC cells: epithelial tumour-up program",
      score = collapsed$epithelial_target_mean[
        match(epithelial_groups, collapsed$biological_group)
      ],
      stringsAsFactors = FALSE
    )
  )

  if (anyNA(plot_data$biological_group) || any(!is.finite(plot_data$score))) {
    stop("Plot data are incomplete after biological-group collapsing.")
  }

  plot_levels <- c(fibroblast_groups, epithelial_groups)
  plot_data$biological_group <- factor(
    plot_data$biological_group,
    levels = rev(plot_levels)
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(score, biological_group)
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      colour = "#BDBDBD",
      linewidth = 0.4
    ) +
    ggplot2::geom_col(fill = "#777777", width = 0.72) +
    ggplot2::facet_wrap(
      ~ program,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::labs(
      x = "Locked program score",
      y = NULL,
      title = "CRC-fibroblast 3D coculture (GSE155343)",
      subtitle = sprintf(
        paste0(
          "Descriptive only: technical replicates collapsed; ",
          "ECM coverage %.0f%%, epithelial-target coverage %.0f%%"
        ),
        100 * signature_fraction,
        100 * target_fraction
      )
    ) +
    theme_publication()

  save_plot(
    p,
    file.path(out, "GSE155343_descriptive_program_scores"),
    10,
    5.5
  )
  # Compatibility copy for any downstream script that expects the original
  # Step 07 plot basename.
  save_plot(
    p,
    file.path(out, "GSE155343_descriptive_ECM_score"),
    10,
    5.5
  )

  write_status(
    "07_GSE155343",
    "PASS",
    sprintf(
      paste0(
        "Technical replicates collapsed; descriptive analysis only. ",
        "ECM_LIGAND_21 coverage %d/%d (%.1f%%); epithelial-target ",
        "coverage %d/%d (%.1f%%)."
      ),
      sum(signature_coverage$detected),
      nrow(signature_coverage),
      100 * signature_fraction,
      sum(target_coverage$detected),
      nrow(target_coverage),
      100 * target_fraction
    )
  )
}
