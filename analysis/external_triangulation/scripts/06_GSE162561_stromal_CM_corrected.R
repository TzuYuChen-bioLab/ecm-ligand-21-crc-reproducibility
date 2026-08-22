PIPELINE_ROOT <- normalizePath(
  Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))

require_packages(
  c(
    "data.table",
    "edgeR",
    "limma",
    "ggplot2",
    "AnnotationDbi",
    "org.Hs.eg.db"
  )
)

out <- result_dir("06_GSE162561")

if (
  !ANALYSIS_OPTIONS$run_optional_GSE162561 ||
    !file.exists(INPUT_FILES$GSE162561_counts)
) {
  write_status(
    "06_GSE162561",
    "SKIPPED",
    "Optional file absent or analysis disabled."
  )
} else {
  counts <- numeric_expression_matrix(
    INPUT_FILES$GSE162561_counts,
    "sum"
  )
  target <- read_gene_set(TARGET_FILE)

  target <- unique(trimws(as.character(target)))
  target <- target[!is.na(target) & nzchar(target)]

  if (length(target) != 50L) {
    stop(
      "The locked epithelial target set must contain exactly 50 unique genes; ",
      "detected ", length(target), "."
    )
  }

  # The GEO series contains 12 libraries: two patient-derived CRC models,
  # three conditions and two deposited replicates per model-condition pair.
  # The supplementary count matrix uses G605/ME59 labels rather than GSM IDs.
  # This explicit map prevents order-only or silently guessed sample assignment.
  official <- data.frame(
    gsm = paste0("GSM", 4954445:4954456),
    matrix_column_expected = c(
      "G605-control-1",
      "G605-S-ASC-1",
      "G605-V-ASC-1",
      "G605-control-2",
      "G605-S-ASC-2",
      "G605-V-ASC-2",
      "ME59-control-1",
      "ME59-S-ASC-1",
      "ME59-V-ASC-1",
      "ME59-control-2",
      "ME59-S-ASC-2",
      "ME59-V-ASC-2"
    ),
    model = rep(c("CRC8", "CRC9"), each = 6),
    condition = rep(
      c("Control", "S_ASC", "V_ASC", "Control", "S_ASC", "V_ASC"),
      2
    ),
    replicate = rep(c(1, 1, 1, 2, 2, 2), 2),
    stringsAsFactors = FALSE
  )

  idx <- match(
    tolower(colnames(counts)),
    tolower(official$matrix_column_expected)
  )

  # Also accept matrices whose columns contain the official GSM accessions.
  if (anyNA(idx)) {
    gsm_idx <- vapply(
      colnames(counts),
      function(nm) {
        hit <- which(
          vapply(
            official$gsm,
            function(g) grepl(g, nm, ignore.case = TRUE),
            logical(1)
          )
        )
        if (length(hit) == 1L) hit else NA_integer_
      },
      integer(1)
    )
    idx[is.na(idx)] <- gsm_idx[is.na(idx)]
  }

  mapping_is_valid <- (
    ncol(counts) == 12L &&
      length(idx) == 12L &&
      all(!is.na(idx)) &&
      length(unique(idx)) == 12L
  )

  if (!mapping_is_valid) {
    manual_mapping_stop(
      data.frame(
        matrix_column = colnames(counts),
        matched_official_row = idx,
        stringsAsFactors = FALSE
      ),
      file.path(out, "sample_mapping_NEEDS_REVIEW.csv"),
      paste(
        "Expected the 12 official GSE162561 libraries and would not use an",
        "order-only fallback."
      )
    )
  }

  meta <- official[idx, , drop = FALSE]
  meta$matrix_column <- colnames(counts)
  meta$mapping_rule <- ifelse(
    tolower(meta$matrix_column) == tolower(meta$matrix_column_expected),
    "explicit_deposited_matrix_label",
    "official_GSM_accession"
  )
  meta <- meta[
    ,
    c(
      "matrix_column",
      "matrix_column_expected",
      "gsm",
      "model",
      "condition",
      "replicate",
      "mapping_rule"
    )
  ]
  write_csv(meta, file.path(out, "sample_mapping_audit.csv"))

  meta$model <- factor(meta$model)
  meta$condition <- factor(
    meta$condition,
    levels = c("Control", "S_ASC", "V_ASC")
  )

  design <- stats::model.matrix(~ model + condition, meta)
  contrasts <- list(
    S_ASC_vs_Control = "conditionS_ASC",
    V_ASC_vs_Control = "conditionV_ASC"
  )

  de_all <- lapply(
    names(contrasts),
    function(label) {
      con <- rep(0, ncol(design))
      contrast_column <- match(
        contrasts[[label]],
        colnames(design)
      )
      if (is.na(contrast_column)) {
        stop("Required design column not found for ", label, ".")
      }
      con[contrast_column] <- 1
      x <- fit_edger(counts, design, con)
      x$contrast <- label
      x$experimental_unit <- paste(
        "two_patient_derived_models_with_within_model_replicates;",
        "exploratory_model_adjusted_library_level_DE"
      )
      x
    }
  )
  write_csv(
    do.call(rbind, de_all),
    file.path(out, "conditioned_medium_DEG.csv")
  )

  # ----------------------------------------------------------------------
  # Locked target-set ID reconciliation
  # ----------------------------------------------------------------------
  # The deposited matrix is indexed entirely by Ensembl gene IDs, whereas
  # EPITHELIAL_TUMOR_UP_50.csv contains HGNC symbols. Mapping is performed
  # locally with org.Hs.eg.db; no study data or annotations are downloaded.
  # Direct SYMBOL mappings are preferred. ALIAS mappings are used only for
  # locked symbols that have no direct current-symbol mapping. Every mapping
  # is written to disk before the pre-specified 60% coverage gate is applied.
  orgdb <- org.Hs.eg.db::org.Hs.eg.db

  valid_symbols <- intersect(
    target,
    AnnotationDbi::keys(orgdb, keytype = "SYMBOL")
  )

  direct_raw <- if (length(valid_symbols) > 0L) {
    AnnotationDbi::select(
      orgdb,
      keys = valid_symbols,
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ENSEMBL")
    )
  } else {
    data.frame(
      SYMBOL = character(0),
      ENSEMBL = character(0),
      stringsAsFactors = FALSE
    )
  }

  direct_map <- data.frame(
    target_gene = direct_raw$SYMBOL,
    current_symbol = direct_raw$SYMBOL,
    ensembl = direct_raw$ENSEMBL,
    mapping_route = "direct_SYMBOL",
    stringsAsFactors = FALSE
  )
  direct_map <- direct_map[
    !is.na(direct_map$ensembl) & nzchar(direct_map$ensembl),
    ,
    drop = FALSE
  ]

  directly_mapped <- unique(direct_map$target_gene)
  alias_candidates <- setdiff(target, directly_mapped)
  valid_aliases <- intersect(
    alias_candidates,
    AnnotationDbi::keys(orgdb, keytype = "ALIAS")
  )

  alias_raw <- if (length(valid_aliases) > 0L) {
    AnnotationDbi::select(
      orgdb,
      keys = valid_aliases,
      keytype = "ALIAS",
      columns = c("SYMBOL", "ENSEMBL")
    )
  } else {
    data.frame(
      ALIAS = character(0),
      SYMBOL = character(0),
      ENSEMBL = character(0),
      stringsAsFactors = FALSE
    )
  }

  alias_map <- data.frame(
    target_gene = alias_raw$ALIAS,
    current_symbol = alias_raw$SYMBOL,
    ensembl = alias_raw$ENSEMBL,
    mapping_route = "fallback_ALIAS",
    stringsAsFactors = FALSE
  )
  alias_map <- alias_map[
    !is.na(alias_map$ensembl) & nzchar(alias_map$ensembl),
    ,
    drop = FALSE
  ]

  target_map <- unique(rbind(direct_map, alias_map))
  target_map$ensembl <- sub(
    "[.][0-9]+$",
    "",
    as.character(target_map$ensembl)
  )

  matrix_ensembl <- sub(
    "[.][0-9]+$",
    "",
    rownames(counts)
  )

  if (anyDuplicated(matrix_ensembl)) {
    stop(
      "Version-stripped Ensembl identifiers are duplicated in the count ",
      "matrix. Resolve them at the raw-count level before scoring."
    )
  }

  rownames(counts) <- matrix_ensembl
  target_map$present_in_matrix <- (
    !is.na(target_map$ensembl) &
      target_map$ensembl %in% rownames(counts)
  )
  target_map$used_for_scoring <- target_map$present_in_matrix

  completely_unmapped <- setdiff(
    target,
    unique(target_map$target_gene)
  )
  if (length(completely_unmapped) > 0L) {
    target_map <- rbind(
      target_map,
      data.frame(
        target_gene = completely_unmapped,
        current_symbol = NA_character_,
        ensembl = NA_character_,
        mapping_route = "unmapped",
        present_in_matrix = FALSE,
        used_for_scoring = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }

  target_map <- target_map[
    order(
      match(target_map$target_gene, target),
      target_map$mapping_route,
      target_map$ensembl
    ),
    ,
    drop = FALSE
  ]

  write_csv(
    target_map,
    file.path(out, "EPITHELIAL_TUMOR_UP_50_symbol_Ensembl_mapping_audit.csv")
  )

  coverage_audit <- do.call(
    rbind,
    lapply(
      target,
      function(gene) {
        rows <- target_map[
          target_map$target_gene == gene,
          ,
          drop = FALSE
        ]
        mapped_ids <- unique(
          rows$ensembl[
            !is.na(rows$ensembl) & nzchar(rows$ensembl)
          ]
        )
        detected_ids <- unique(
          rows$ensembl[rows$present_in_matrix %in% TRUE]
        )
        current_symbols <- unique(
          rows$current_symbol[
            !is.na(rows$current_symbol) & nzchar(rows$current_symbol)
          ]
        )
        data.frame(
          target_gene = gene,
          current_symbols = paste(current_symbols, collapse = ";"),
          n_Ensembl_mapped = length(mapped_ids),
          mapped_Ensembl = paste(mapped_ids, collapse = ";"),
          n_Ensembl_detected = length(detected_ids),
          detected_Ensembl = paste(detected_ids, collapse = ";"),
          detected = length(detected_ids) > 0L,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  target_coverage <- mean(coverage_audit$detected)
  coverage_audit$overall_coverage <- target_coverage
  coverage_audit$minimum_required <- (
    ANALYSIS_OPTIONS$minimum_target_coverage
  )
  coverage_audit$passes_prespecified_gate <- (
    target_coverage >= ANALYSIS_OPTIONS$minimum_target_coverage
  )

  write_csv(
    coverage_audit,
    file.path(out, "EPITHELIAL_TUMOR_UP_50_coverage.csv")
  )

  annotation_provenance <- data.frame(
    field = c(
      "matrix_identifier_type",
      "locked_target_identifier_type",
      "mapping_package",
      "mapping_package_version",
      "mapping_precedence",
      "multi_Ensembl_rule",
      "target_genes_detected",
      "target_genes_total",
      "target_coverage",
      "minimum_required_coverage"
    ),
    value = c(
      "Ensembl gene ID",
      "HGNC symbol",
      "org.Hs.eg.db",
      as.character(utils::packageVersion("org.Hs.eg.db")),
      "direct SYMBOL first; ALIAS only when direct mapping absent",
      "mean log-CPM across detected Ensembl IDs so each locked symbol has equal weight",
      sum(coverage_audit$detected),
      nrow(coverage_audit),
      target_coverage,
      ANALYSIS_OPTIONS$minimum_target_coverage
    ),
    stringsAsFactors = FALSE
  )
  write_csv(
    annotation_provenance,
    file.path(out, "identifier_mapping_provenance.csv")
  )

  if (target_coverage < ANALYSIS_OPTIONS$minimum_target_coverage) {
    stop(
      sprintf(
        paste0(
          "Mapped target coverage %.1f%% is below the pre-specified %.1f%% ",
          "threshold. See EPITHELIAL_TUMOR_UP_50_coverage.csv."
        ),
        100 * target_coverage,
        100 * ANALYSIS_OPTIONS$minimum_target_coverage
      )
    )
  }

  lcpm <- log_cpm(counts)
  detected_targets <- coverage_audit$target_gene[
    coverage_audit$detected
  ]

  target_lcpm <- matrix(
    NA_real_,
    nrow = length(detected_targets),
    ncol = ncol(lcpm),
    dimnames = list(detected_targets, colnames(lcpm))
  )

  for (gene in detected_targets) {
    detected_ids <- unique(
      target_map$ensembl[
        target_map$target_gene == gene &
          target_map$present_in_matrix %in% TRUE
      ]
    )
    detected_ids <- detected_ids[
      !is.na(detected_ids) & detected_ids %in% rownames(lcpm)
    ]

    if (length(detected_ids) == 1L) {
      target_lcpm[gene, ] <- lcpm[detected_ids, ]
    } else {
      target_lcpm[gene, ] <- colMeans(
        lcpm[detected_ids, , drop = FALSE]
      )
    }
  }

  sc <- score_gene_set(
    target_lcpm,
    target,
    ANALYSIS_OPTIONS$minimum_target_coverage
  )

  scores <- cbind(
    meta,
    epithelial_target_score = as.numeric(sc$score),
    coverage = sc$coverage
  )

  model_means <- aggregate(
    epithelial_target_score ~ model + condition,
    scores,
    mean
  )
  wide <- reshape(
    model_means,
    idvar = "model",
    timevar = "condition",
    direction = "wide"
  )
  wide$delta_S_ASC <- (
    wide$epithelial_target_score.S_ASC -
      wide$epithelial_target_score.Control
  )
  wide$delta_V_ASC <- (
    wide$epithelial_target_score.V_ASC -
      wide$epithelial_target_score.Control
  )
  wide$inference_note <- paste(
    "two patient-derived models; replicate libraries are averaged within",
    "model-condition and are not treated as independent patients"
  )

  write_csv(
    scores,
    file.path(out, "epithelial_target_library_scores.csv")
  )
  write_csv(
    wide,
    file.path(out, "epithelial_target_model_level_effects.csv")
  )

  p <- ggplot2::ggplot(
    model_means,
    ggplot2::aes(
      condition,
      epithelial_target_score,
      group = model,
      colour = model
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2.7) +
    ggplot2::labs(
      x = NULL,
      y = "Epithelial tumour-up program",
      title = "Stromal conditioned-medium response (GSE162561)",
      subtitle = paste(
        "Replicates averaged within each patient-derived CRC model;",
        sprintf("target coverage %.0f%%", 100 * target_coverage)
      )
    ) +
    theme_publication()

  save_plot(
    p,
    file.path(out, "GSE162561_model_level_response"),
    6.5,
    4.5
  )

  write_status(
    "06_GSE162561",
    "PASS",
    sprintf(
      paste0(
        "Optional two-model conditioned-medium analysis completed. ",
        "Locked target coverage: %d/%d (%.1f%%)."
      ),
      sum(coverage_audit$detected),
      nrow(coverage_audit),
      100 * target_coverage
    )
  )
}
