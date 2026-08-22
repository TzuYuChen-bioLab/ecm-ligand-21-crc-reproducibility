PIPELINE_ROOT <- normalizePath(
  Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))

require_packages(c("readxl", "ggplot2"))

out <- result_dir("08_ECM_proteomics")
assert_file(INPUT_FILES$ECM_proteomics)

signature <- unique(
  toupper(trimws(as.character(read_gene_set(SIGNATURE_FILE))))
)
signature <- signature[!is.na(signature) & nzchar(signature)]

if (length(signature) != 21L) {
  stop(
    "The locked ECM_LIGAND_21 set must contain exactly 21 unique genes; ",
    "detected ", length(signature), "."
  )
}

# ----------------------------------------------------------------------
# 1. Workbook and sheet audit
# ----------------------------------------------------------------------
sheets <- readxl::excel_sheets(INPUT_FILES$ECM_proteomics)
expected_sheet <- "C. Normalized intensity"

sheet_audit <- do.call(
  rbind,
  lapply(
    sheets,
    function(sheet_name) {
      preview <- tryCatch(
        readxl::read_excel(
          INPUT_FILES$ECM_proteomics,
          sheet = sheet_name,
          n_max = 8,
          na = c("", "NA", "#N/A"),
          .name_repair = "minimal"
        ),
        error = function(e) NULL
      )

      if (is.null(preview)) {
        return(
          data.frame(
            sheet = sheet_name,
            n_columns = NA_integer_,
            tissue_channel_columns = NA_integer_,
            reference_channel_columns = NA_integer_,
            selected = FALSE,
            columns = NA_character_,
            stringsAsFactors = FALSE
          )
        )
      }

      header <- names(preview)
      data.frame(
        sheet = sheet_name,
        n_columns = length(header),
        tissue_channel_columns = sum(
          grepl("^SEV[0-9]+_[NT]$", header, ignore.case = TRUE)
        ),
        reference_channel_columns = sum(
          grepl("^R[1-4]-[1-3]$", header, ignore.case = TRUE)
        ),
        selected = identical(sheet_name, expected_sheet),
        columns = paste(header, collapse = " | "),
        stringsAsFactors = FALSE
      )
    }
  )
)

write_csv(
  sheet_audit,
  file.path(out, "workbook_sheet_audit.csv")
)

if (!expected_sheet %in% sheets) {
  manual_mapping_stop(
    sheet_audit,
    file.path(out, "workbook_sheet_mapping_NEEDS_REVIEW.csv"),
    paste0(
      "The source-defined normalized-intensity sheet was not found: ",
      expected_sheet
    )
  )
}

proteomics <- as.data.frame(
  readxl::read_excel(
    INPUT_FILES$ECM_proteomics,
    sheet = expected_sheet,
    na = c("", "NA", "#N/A"),
    .name_repair = "minimal"
  ),
  check.names = FALSE
)

header <- names(proteomics)
gene_column <- which(toupper(trimws(header)) == "GENE")
locus_column <- which(toupper(trimws(header)) == "LOCUS")
tissue_columns <- grep(
  "^SEV[0-9]+_[NT]$",
  header,
  ignore.case = TRUE
)
reference_columns <- grep(
  "^R[1-4]-[1-3]$",
  header,
  ignore.case = TRUE
)

structure_is_valid <- (
  length(gene_column) == 1L &&
    length(locus_column) == 1L &&
    length(tissue_columns) == 32L &&
    length(reference_columns) == 12L &&
    ncol(proteomics) == 46L
)

if (!structure_is_valid) {
  manual_mapping_stop(
    data.frame(
      column_index = seq_along(header),
      column = header,
      is_tissue_channel = seq_along(header) %in% tissue_columns,
      is_reference_channel = seq_along(header) %in% reference_columns,
      stringsAsFactors = FALSE
    ),
    file.path(out, "proteomics_workbook_structure_NEEDS_REVIEW.csv"),
    paste(
      "Expected LOCUS, GENE, 32 tissue channels and 12 internal-reference",
      "channels in 'C. Normalized intensity'."
    )
  )
}

channel_set <- function(column_index) {
  findInterval(column_index, c(3L, 14L, 25L, 36L))
}

channel_audit <- data.frame(
  column_index = seq_along(header),
  original_column = header,
  channel_type = ifelse(
    seq_along(header) %in% tissue_columns,
    "tissue",
    ifelse(
      seq_along(header) %in% reference_columns,
      "pooled_internal_reference",
      "annotation"
    )
  ),
  TMT_set = ifelse(
    seq_along(header) >= 3L,
    channel_set(seq_along(header)),
    NA_integer_
  ),
  included_in_tissue_analysis = seq_along(header) %in% tissue_columns,
  stringsAsFactors = FALSE
)

write_csv(
  channel_audit,
  file.path(out, "proteomics_channel_role_audit.csv")
)

# ----------------------------------------------------------------------
# 2. Explicit tissue-channel mapping and internal-reference exclusion
# ----------------------------------------------------------------------
tissue_labels <- toupper(trimws(header[tissue_columns]))
tissue_patient <- sub("_.*$", "", tissue_labels)
tissue_condition <- ifelse(
  grepl("_N$", tissue_labels),
  "Normal",
  "Tumor"
)
tissue_set <- channel_set(tissue_columns)
tissue_channel_id <- paste0(
  tissue_labels,
  "__TMT",
  tissue_set,
  "__COL",
  tissue_columns
)

sample_meta <- data.frame(
  matrix_column = tissue_channel_id,
  original_column = tissue_labels,
  column_index = tissue_columns,
  TMT_set = tissue_set,
  patient = tissue_patient,
  condition = tissue_condition,
  channel_role = "biological_tissue_measurement",
  stringsAsFactors = FALSE
)

if (
  any(!grepl("^SEV[0-9]+$", sample_meta$patient)) ||
    any(!sample_meta$condition %in% c("Normal", "Tumor"))
) {
  manual_mapping_stop(
    sample_meta,
    file.path(out, "proteomics_sample_mapping_NEEDS_REVIEW.csv"),
    "Some tissue channels lack an unambiguous patient or condition."
  )
}

write_csv(
  sample_meta,
  file.path(out, "proteomics_sample_mapping_audit.csv")
)

reference_audit <- channel_audit[
  channel_audit$channel_type == "pooled_internal_reference",
  ,
  drop = FALSE
]
reference_audit$exclusion_reason <- paste(
  "pooled internal TMT reference; technical normalization/QC channel,",
  "not a patient tissue sample"
)
write_csv(
  reference_audit,
  file.path(out, "proteomics_internal_reference_exclusion_audit.csv")
)

# ----------------------------------------------------------------------
# 3. Gene-symbol matrix on the source-normalized log2 abundance scale
# ----------------------------------------------------------------------
normalise_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("[;|, ].*$", "", x)
  x <- toupper(x)
  x[is.na(x) | !nzchar(x) | x == "NA"] <- NA_character_
  x
}

genes_raw <- as.character(proteomics[[gene_column]])
genes <- normalise_symbol(genes_raw)
locus <- as.character(proteomics[[locus_column]])

linear_abundance <- do.call(
  cbind,
  lapply(
    tissue_columns,
    function(column_index) {
      suppressWarnings(as.numeric(proteomics[[column_index]]))
    }
  )
)
colnames(linear_abundance) <- tissue_channel_id

linear_abundance[
  !is.finite(linear_abundance) | linear_abundance <= 0
] <- NA_real_
log2_abundance <- log2(linear_abundance)

valid_gene <- !is.na(genes) & nzchar(genes)
genes_valid <- genes[valid_gene]
log2_valid <- log2_abundance[valid_gene, , drop = FALSE]

gene_rows <- split(seq_along(genes_valid), genes_valid)
gene_expression <- t(
  vapply(
    gene_rows,
    function(row_index) {
      colMeans(
        log2_valid[row_index, , drop = FALSE],
        na.rm = TRUE
      )
    },
    numeric(ncol(log2_valid))
  )
)
gene_expression[!is.finite(gene_expression)] <- NA_real_
colnames(gene_expression) <- tissue_channel_id

gene_counts <- table(genes_valid)
n_rows_for_symbol <- rep(NA_integer_, length(genes))
n_rows_for_symbol[valid_gene] <- as.integer(
  gene_counts[genes[valid_gene]]
)

feature_audit <- data.frame(
  workbook_row = seq_len(nrow(proteomics)) + 1L,
  LOCUS = locus,
  GENE_deposited = genes_raw,
  gene_symbol_used = genes,
  valid_gene_symbol = valid_gene,
  n_rows_for_gene_symbol = n_rows_for_symbol,
  duplicate_gene_rule = ifelse(
    valid_gene & n_rows_for_symbol > 1L,
    "mean log2 abundance across deposited rows sharing the symbol",
    ifelse(valid_gene, "single deposited row", "invalid/missing symbol")
  ),
  in_ECM_LIGAND_21 = valid_gene & genes %in% signature,
  stringsAsFactors = FALSE
)

write_csv(
  feature_audit,
  file.path(out, "proteomics_feature_identifier_audit.csv")
)

# ----------------------------------------------------------------------
# 4. Locked signature coverage with the original 70% gate retained
# ----------------------------------------------------------------------
coverage_audit <- do.call(
  rbind,
  lapply(
    signature,
    function(gene) {
      present <- gene %in% rownames(gene_expression)
      n_finite <- if (present) {
        sum(is.finite(gene_expression[gene, ]))
      } else {
        0L
      }
      data.frame(
        gene = gene,
        present_in_workbook = present,
        finite_tissue_channels = n_finite,
        total_tissue_channels = ncol(gene_expression),
        detected_for_scoring = present && n_finite == ncol(gene_expression),
        stringsAsFactors = FALSE
      )
    }
  )
)

coverage <- mean(coverage_audit$detected_for_scoring)
coverage_audit$overall_coverage <- coverage
coverage_audit$minimum_required <- (
  ANALYSIS_OPTIONS$minimum_signature_coverage
)
coverage_audit$passes_prespecified_gate <- (
  coverage >= ANALYSIS_OPTIONS$minimum_signature_coverage
)

write_csv(
  coverage_audit,
  file.path(out, "ECM_LIGAND_21_protein_coverage.csv")
)

if (coverage < ANALYSIS_OPTIONS$minimum_signature_coverage) {
  stop(
    sprintf(
      paste0(
        "ECM_LIGAND_21 protein coverage %.1f%% is below the locked %.1f%% ",
        "threshold. See ECM_LIGAND_21_protein_coverage.csv."
      ),
      100 * coverage,
      100 * ANALYSIS_OPTIONS$minimum_signature_coverage
    )
  )
}

detected_signature <- coverage_audit$gene[
  coverage_audit$detected_for_scoring
]

# ----------------------------------------------------------------------
# 5. Collapse repeated TMT measurements of the same tissue sample
# ----------------------------------------------------------------------
sample_meta$sample <- paste(
  sample_meta$patient,
  sample_meta$condition,
  sep = "|"
)
sample_levels <- unique(sample_meta$sample)

collapsed <- vapply(
  sample_levels,
  function(sample_id) {
    column_index <- which(sample_meta$sample == sample_id)
    rowMeans(
      gene_expression[, column_index, drop = FALSE],
      na.rm = TRUE
    )
  },
  numeric(nrow(gene_expression))
)
collapsed[!is.finite(collapsed)] <- NA_real_
rownames(collapsed) <- rownames(gene_expression)
colnames(collapsed) <- sample_levels

collapsed_meta <- do.call(
  rbind,
  lapply(
    sample_levels,
    function(sample_id) {
      rows <- sample_meta[sample_meta$sample == sample_id, , drop = FALSE]
      data.frame(
        sample = sample_id,
        patient = rows$patient[1],
        condition = rows$condition[1],
        n_TMT_channels_collapsed = nrow(rows),
        TMT_sets = paste(sort(unique(rows$TMT_set)), collapse = ";"),
        original_columns = paste(rows$original_column, collapse = ";"),
        matrix_columns = paste(rows$matrix_column, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  )
)
rownames(collapsed_meta) <- NULL

if (nrow(collapsed_meta) != 24L) {
  manual_mapping_stop(
    collapsed_meta,
    file.path(out, "collapsed_proteomics_design_NEEDS_REVIEW.csv"),
    paste(
      "Expected 24 unique tissue samples after collapsing repeated TMT",
      "measurements."
    )
  )
}

write_csv(
  collapsed_meta,
  file.path(out, "collapsed_proteomics_design.csv")
)

technical_duplicate_audit <- collapsed_meta[
  collapsed_meta$n_TMT_channels_collapsed > 1L,
  ,
  drop = FALSE
]
technical_duplicate_audit$collapse_rule <- paste(
  "mean of source-normalized log2 abundance across repeated TMT",
  "measurements; not treated as independent tissues"
)
write_csv(
  technical_duplicate_audit,
  file.path(out, "proteomics_repeated_TMT_measurement_audit.csv")
)

# ----------------------------------------------------------------------
# 6. Source-defined pathological exclusions
# ----------------------------------------------------------------------
collapsed_meta$source_excluded <- FALSE
collapsed_meta$source_exclusion_reason <- NA_character_

exclusion_map <- data.frame(
  sample = c(
    "SEV01|Tumor",
    "SEV04|Tumor",
    "SEV09|Normal"
  ),
  reason = c(
    "perforation",
    "chemotherapy/pre-treatment",
    "stent insertion"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(exclusion_map))) {
  hit <- collapsed_meta$sample == exclusion_map$sample[i]
  collapsed_meta$source_excluded[hit] <- TRUE
  collapsed_meta$source_exclusion_reason[hit] <- exclusion_map$reason[i]
}
collapsed_meta$source_eligible <- !collapsed_meta$source_excluded

if (
  sum(collapsed_meta$source_excluded) != 3L ||
    sum(collapsed_meta$source_eligible) != 21L
) {
  manual_mapping_stop(
    collapsed_meta,
    file.path(out, "source_exclusion_mapping_NEEDS_REVIEW.csv"),
    "The three source-defined pathological exclusions were not recovered."
  )
}

write_csv(
  collapsed_meta,
  file.path(out, "proteomics_source_exclusion_audit.csv")
)

# Rewrite the legacy design file with eligibility fields included.
write_csv(
  collapsed_meta,
  file.path(out, "collapsed_proteomics_design.csv")
)

# ----------------------------------------------------------------------
# 7. Locked sample-level protein score
# ----------------------------------------------------------------------
score_object <- score_gene_set(
  collapsed,
  signature,
  minimum_coverage = ANALYSIS_OPTIONS$minimum_signature_coverage
)

scores <- cbind(
  collapsed_meta,
  score = as.numeric(score_object$score),
  detected_signature_proteins = length(detected_signature),
  total_signature_proteins = length(signature),
  coverage = coverage
)
scores$inference_note <- paste(
  "sample-level locked protein score after repeated-TMT collapse;",
  "source exclusions and partial pairing must be respected"
)

write_csv(
  scores,
  file.path(out, "ECM_LIGAND_21_protein_sample_scores.csv")
)

eligible_scores <- scores[scores$source_eligible, , drop = FALSE]

clean_group_summary <- do.call(
  rbind,
  lapply(
    c("Normal", "Tumor"),
    function(condition_label) {
      value <- eligible_scores$score[
        eligible_scores$condition == condition_label
      ]
      data.frame(
        condition = condition_label,
        n_tissues = length(value),
        mean = mean(value),
        sd = stats::sd(value),
        median = stats::median(value),
        IQR = stats::IQR(value),
        minimum = min(value),
        maximum = max(value),
        stringsAsFactors = FALSE
      )
    }
  )
)

clean_group_difference <- data.frame(
  analysis_scope = "source_eligible_all_tissues",
  n_normal = sum(eligible_scores$condition == "Normal"),
  n_tumor = sum(eligible_scores$condition == "Tumor"),
  mean_difference_tumor_minus_normal = (
    mean(eligible_scores$score[eligible_scores$condition == "Tumor"]) -
      mean(eligible_scores$score[eligible_scores$condition == "Normal"])
  ),
  median_difference_tumor_minus_normal = (
    stats::median(
      eligible_scores$score[eligible_scores$condition == "Tumor"]
    ) -
      stats::median(
        eligible_scores$score[eligible_scores$condition == "Normal"]
      )
  ),
  inference_note = paste(
    "descriptive full-cohort contrast after the three source exclusions;",
    "the cohort is partially paired, so no independence-based P value is",
    "reported"
  ),
  stringsAsFactors = FALSE
)

write_csv(
  clean_group_summary,
  file.path(out, "ECM_LIGAND_21_clean_cohort_group_summary.csv")
)
write_csv(
  clean_group_difference,
  file.path(out, "ECM_LIGAND_21_clean_cohort_group_difference.csv")
)

# ----------------------------------------------------------------------
# 8. Strict two-pair analysis and flagged three-pair sensitivity
# ----------------------------------------------------------------------
make_score_pairs <- function(score_data, patients, analysis_scope, warning) {
  pairs <- do.call(
    rbind,
    lapply(
      patients,
      function(patient_id) {
        normal <- score_data$score[
          score_data$patient == patient_id &
            score_data$condition == "Normal"
        ]
        tumor <- score_data$score[
          score_data$patient == patient_id &
            score_data$condition == "Tumor"
        ]

        if (length(normal) != 1L || length(tumor) != 1L) {
          stop(
            "Expected one Normal and one Tumor score for paired patient ",
            patient_id, "."
          )
        }

        data.frame(
          analysis_scope = analysis_scope,
          patient = patient_id,
          score_Normal = normal,
          score_Tumor = tumor,
          delta_Tumor_minus_Normal = tumor - normal,
          source_warning = warning,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  summary <- data.frame(
    analysis_scope = analysis_scope,
    n_pairs = nrow(pairs),
    patients = paste(pairs$patient, collapse = ";"),
    mean_difference = mean(pairs$delta_Tumor_minus_Normal),
    median_difference = stats::median(pairs$delta_Tumor_minus_Normal),
    minimum_difference = min(pairs$delta_Tumor_minus_Normal),
    maximum_difference = max(pairs$delta_Tumor_minus_Normal),
    patients_with_increase = sum(pairs$delta_Tumor_minus_Normal > 0),
    patients_with_decrease = sum(pairs$delta_Tumor_minus_Normal < 0),
    P_value = NA_real_,
    inference_note = paste(
      "descriptive paired differences only; sample size is insufficient for",
      "a confirmatory P value;",
      warning
    ),
    stringsAsFactors = FALSE
  )

  list(pairs = pairs, summary = summary)
}

strict_pair_patients <- c("SEV13", "SEV14")
flagged_pair_patients <- c("SEV09", "SEV13", "SEV14")

strict_pairs <- make_score_pairs(
  scores,
  strict_pair_patients,
  "strict_source_aligned_pairs",
  "all included tissues are source-eligible"
)
flagged_pairs <- make_score_pairs(
  scores,
  flagged_pair_patients,
  "all_three_pairs_flagged_sensitivity",
  "SEV09 Normal tissue was excluded by the source study because of stent insertion"
)

write_csv(
  rbind(strict_pairs$pairs, flagged_pairs$pairs),
  file.path(out, "ECM_LIGAND_21_protein_paired_differences.csv")
)
write_csv(
  rbind(strict_pairs$summary, flagged_pairs$summary),
  file.path(out, "ECM_LIGAND_21_protein_effect_summary.csv")
)

make_protein_pair_statistics <- function(
  expression,
  patients,
  analysis_scope,
  source_warning
) {
  delta_list <- lapply(
    patients,
    function(patient_id) {
      tumor_column <- paste(patient_id, "Tumor", sep = "|")
      normal_column <- paste(patient_id, "Normal", sep = "|")

      if (!all(c(tumor_column, normal_column) %in% colnames(expression))) {
        stop("Missing paired protein columns for ", patient_id, ".")
      }

      expression[, tumor_column] - expression[, normal_column]
    }
  )
  delta <- do.call(cbind, delta_list)
  colnames(delta) <- paste0("delta_", patients)

  involved_columns <- unlist(
    lapply(
      patients,
      function(patient_id) {
        c(
          paste(patient_id, "Normal", sep = "|"),
          paste(patient_id, "Tumor", sep = "|")
        )
      }
    )
  )

  n_observed <- rowSums(is.finite(delta))
  mean_delta <- rowMeans(delta, na.rm = TRUE)
  mean_delta[n_observed == 0L] <- NA_real_
  median_delta <- apply(
    delta,
    1,
    function(value) {
      value <- value[is.finite(value)]
      if (length(value)) stats::median(value) else NA_real_
    }
  )

  result <- data.frame(
    gene = rownames(expression),
    logFC = mean_delta,
    median_logFC = median_delta,
    AveExpr = rowMeans(
      expression[, involved_columns, drop = FALSE],
      na.rm = TRUE
    ),
    n_pairs_observed = n_observed,
    n_tumor_higher = rowSums(delta > 0, na.rm = TRUE),
    n_tumor_lower = rowSums(delta < 0, na.rm = TRUE),
    t = NA_real_,
    P.Value = NA_real_,
    adj.P.Val = NA_real_,
    B = NA_real_,
    contrast = paste0(analysis_scope, "_Tumor_vs_Normal"),
    analysis_scope = analysis_scope,
    source_warning = source_warning,
    inference_note = paste(
      "mean paired log2-abundance difference; descriptive only;",
      "no limma P value fitted at this sample size"
    ),
    stringsAsFactors = FALSE
  )
  cbind(result, as.data.frame(delta, check.names = FALSE))
}

strict_statistics <- make_protein_pair_statistics(
  collapsed,
  strict_pair_patients,
  "strict_source_aligned_two_pairs",
  "none"
)
flagged_statistics <- make_protein_pair_statistics(
  collapsed,
  flagged_pair_patients,
  "all_three_pairs_flagged_sensitivity",
  "SEV09 Normal was source-excluded because of stent insertion"
)

# Legacy filename retained, but its contents are explicitly descriptive.
write_csv(
  strict_statistics,
  file.path(out, "paired_proteomics_statistics.csv")
)
write_csv(
  flagged_statistics,
  file.path(out, "all_three_pair_proteomics_sensitivity.csv")
)

# ----------------------------------------------------------------------
# 9. Author-reported DEP overlap from Supplementary Table 1E
# ----------------------------------------------------------------------
dep_sheet <- "E. DEP list"
if (!dep_sheet %in% sheets) {
  stop("The author-reported DEP sheet was not found: ", dep_sheet)
}

dep_table <- as.data.frame(
  readxl::read_excel(
    INPUT_FILES$ECM_proteomics,
    sheet = dep_sheet,
    na = c("", "NA", "#N/A"),
    .name_repair = "minimal"
  ),
  check.names = FALSE
)

if (ncol(dep_table) < 2L) {
  stop("The author-reported DEP sheet has fewer than two columns.")
}

normal_dep <- unique(normalise_symbol(dep_table[[1]]))
tumor_dep <- unique(normalise_symbol(dep_table[[2]]))
normal_dep <- normal_dep[!is.na(normal_dep)]
tumor_dep <- tumor_dep[!is.na(tumor_dep)]

dep_overlap <- data.frame(
  gene = signature,
  detected_in_normalized_intensity = (
    signature %in% detected_signature
  ),
  author_DEP_classification = ifelse(
    signature %in% normal_dep,
    "Normal_enriched_DEP",
    ifelse(
      signature %in% tumor_dep,
      "Tumor_enriched_DEP",
      "Not_author_reported_DEP"
    )
  ),
  strict_two_pair_logFC = strict_statistics$logFC[
    match(signature, strict_statistics$gene)
  ],
  flagged_three_pair_logFC = flagged_statistics$logFC[
    match(signature, flagged_statistics$gene)
  ],
  stringsAsFactors = FALSE
)
dep_overlap$author_DEP_direction <- ifelse(
  dep_overlap$author_DEP_classification == "Tumor_enriched_DEP",
  1L,
  ifelse(
    dep_overlap$author_DEP_classification == "Normal_enriched_DEP",
    -1L,
    0L
  )
)
dep_overlap$strict_pair_direction <- sign(
  dep_overlap$strict_two_pair_logFC
)
dep_overlap$strict_pair_concordant_with_author_DEP <- ifelse(
  dep_overlap$author_DEP_direction == 0L ||
    is.na(dep_overlap$strict_pair_direction),
  NA,
  dep_overlap$author_DEP_direction == dep_overlap$strict_pair_direction
)

write_csv(
  dep_overlap,
  file.path(out, "ECM_LIGAND_21_author_DEP_overlap.csv")
)

sig_stats <- strict_statistics[
  match(detected_signature, strict_statistics$gene),
  ,
  drop = FALSE
]
sig_stats$detected_signature_count <- length(detected_signature)
sig_stats$total_signature_count <- length(signature)
sig_stats$coverage <- coverage
sig_stats$author_DEP_classification <- dep_overlap$author_DEP_classification[
  match(sig_stats$gene, dep_overlap$gene)
]

write_csv(
  sig_stats,
  file.path(out, "ECM_LIGAND_21_individual_protein_statistics.csv")
)

# ----------------------------------------------------------------------
# 10. Figures and provenance
# ----------------------------------------------------------------------
plot_stats <- sig_stats[is.finite(sig_stats$logFC), , drop = FALSE]
plot_stats$direction <- ifelse(
  plot_stats$logFC > 0,
  "Tumor higher",
  "Normal higher"
)

p_direction <- ggplot2::ggplot(
  plot_stats,
  ggplot2::aes(
    reorder(gene, logFC),
    logFC,
    fill = direction
  )
) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(
    values = c(
      "Tumor higher" = "#D55E00",
      "Normal higher" = "#0072B2"
    )
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey35",
    linewidth = 0.4
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Mean paired tumour-minus-normal log2 abundance",
    title = "Independent decellularized-ECM proteomics",
    subtitle = paste(
      length(detected_signature),
      "of 21 locked proteins detected; strict source-aligned pairs:",
      paste(strict_pair_patients, collapse = ", "),
      "(descriptive)"
    ),
    fill = NULL
  ) +
  theme_publication()

save_plot(
  p_direction,
  file.path(out, "ECM_proteomics_signature_directions"),
  8,
  max(5.5, 0.30 * nrow(plot_stats))
)

eligible_scores$condition <- factor(
  eligible_scores$condition,
  levels = c("Normal", "Tumor")
)
p_score <- ggplot2::ggplot(
  eligible_scores,
  ggplot2::aes(condition, score, fill = condition)
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.45
  ) +
  ggplot2::geom_jitter(
    width = 0.08,
    height = 0,
    size = 2.2,
    alpha = 0.85
  ) +
  ggplot2::scale_fill_manual(
    values = c("Normal" = "#0072B2", "Tumor" = "#D55E00"),
    guide = "none"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "ECM_LIGAND_21 protein score",
    title = "Source-eligible decellularized-ECM tissues",
    subtitle = paste(
      "Descriptive, partially paired cohort;",
      sum(eligible_scores$condition == "Normal"),
      "normal and",
      sum(eligible_scores$condition == "Tumor"),
      "tumour tissues"
    )
  ) +
  theme_publication()

save_plot(
  p_score,
  file.path(out, "ECM_proteomics_sample_scores"),
  6.5,
  5
)

provenance <- data.frame(
  field = c(
    "source_article",
    "selected_sheet",
    "abundance_scale_input",
    "analysis_scale",
    "tissue_channels",
    "internal_reference_channels_excluded",
    "unique_tissues_after_TMT_collapse",
    "source_excluded_tissues",
    "source_eligible_tissues",
    "strict_paired_patients",
    "flagged_sensitivity_patients",
    "locked_signature_detected",
    "locked_signature_total",
    "locked_signature_coverage",
    "locked_minimum_coverage",
    "primary_inference_scope"
  ),
  value = c(
    "Lee et al., British Journal of Cancer 2025; doi:10.1038/s41416-025-02964-z",
    expected_sheet,
    "source-normalized linear TMT intensity",
    "log2 abundance; repeated TMT measurements averaged within tissue",
    length(tissue_columns),
    length(reference_columns),
    nrow(collapsed_meta),
    sum(collapsed_meta$source_excluded),
    sum(collapsed_meta$source_eligible),
    paste(strict_pair_patients, collapse = ";"),
    paste(flagged_pair_patients, collapse = ";"),
    length(detected_signature),
    length(signature),
    coverage,
    ANALYSIS_OPTIONS$minimum_signature_coverage,
    paste(
      "protein-level and sample-score evidence is descriptive; author DEP",
      "classification is retained as source-reported evidence"
    )
  ),
  stringsAsFactors = FALSE
)

write_csv(
  provenance,
  file.path(out, "ECM_proteomics_analysis_provenance.csv")
)

write_status(
  "08_ECM_proteomics",
  "PASS",
  sprintf(
    paste0(
      "%d/%d locked proteins detected (%.1f%%); 32 tissue channels ",
      "collapsed to 24 tissues; 21 source-eligible tissues; strict paired ",
      "analysis descriptive for SEV13/SEV14; SEV09 retained only as a ",
      "flagged sensitivity."
    ),
    length(detected_signature),
    length(signature),
    100 * coverage
  )
)
