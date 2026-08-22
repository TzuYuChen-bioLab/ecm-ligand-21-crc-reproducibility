PIPELINE_ROOT <- normalizePath(
  Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))
require_packages(c("data.table", "edgeR", "limma", "ggplot2"))

out <- result_dir("05_GSE160686_TGFb")
counts <- numeric_expression_matrix(INPUT_FILES$GSE160686_counts, "sum")
deposit_meta <- read_table_auto(INPUT_FILES$GSE160686_metadata)
signature <- read_gene_set(SIGNATURE_FILE)

# Explicit GEO sample-to-library map. The deposited count-matrix columns are
# formatted as "library_id.cell_barcode.numeric_suffix"; they do not contain
# GSM accessions. Keeping this map explicit prevents silent sample relabelling.
official <- data.frame(
  gsm = paste0("GSM", 4877056:4877073),
  library_id = c(
    "AS_IL_LIB_19_H996_TB_AAC252_07",
    "AS_IL_LIB_19_I004_TB_AAC252_15",
    "AS_IL_LIB_19_H994_TB_AAC252_05",
    "AS_IL_LIB_19_I002_TB_AAC252_13",
    "AS_IL_LIB_19_H995_TB_AAC252_06",
    "AS_IL_LIB_19_I003_TB_AAC252_14",
    "AS_IL_LIB_19_I012_TB_AAC252_23",
    "AS_IL_LIB_19_I010_TB_AAC252_21",
    "AS_IL_LIB_19_I018_TB_AAC252_29",
    "AS_IL_LIB_19_I011_TB_AAC252_22",
    "AS_IL_LIB_19_I019_TB_AAC252_30",
    "AS_IL_LIB_19_I020_TB_AAC252_31",
    "GX_IL_LIB_20_A332_TB_AAC513_06",
    "GX_IL_LIB_20_A340_TB_AAC513_14",
    "GX_IL_LIB_20_A333_TB_AAC513_07",
    "GX_IL_LIB_20_A341_TB_AAC513_15",
    "GX_IL_LIB_20_A334_TB_AAC513_08",
    "GX_IL_LIB_20_A342_TB_AAC513_16"
  ),
  patient = c(rep("P1", 6), rep("P2", 6), rep("P3", 6)),
  treatment = c(
    "DSE", "DSE", "Isotype", "Isotype", "NIS", "NIS",
    "DSE", "Isotype", "Isotype", "NIS", "NIS", "DSE",
    "IgG2", "IgG2", "NIS793", "NIS793", "DSE", "DSE"
  ),
  replicate = c(1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2),
  stringsAsFactors = FALSE
)
official$condition <- ifelse(
  official$treatment %in% c("Isotype", "IgG2"),
  "Control",
  ifelse(
    official$treatment %in% c("NIS", "NIS793"),
    "TGFb_blockade",
    "DSE_reference"
  )
)

required_deposit_columns <- c("library_id", "patient", "treatment")
if (!all(required_deposit_columns %in% colnames(deposit_meta))) {
  manual_mapping_stop(
    data.frame(
      required_column = required_deposit_columns,
      detected = required_deposit_columns %in% colnames(deposit_meta)
    ),
    file.path(out, "deposited_metadata_NEEDS_REVIEW.csv"),
    "The deposited GSE160686 metadata lacks a required mapping column."
  )
}

normalise_patient <- function(x) {
  number <- gsub("[^0-9]", "", as.character(x))
  ifelse(nzchar(number), paste0("P", number), NA_character_)
}

classify_deposited_treatment <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(
    grepl("isotype|igg2", x),
    "Control",
    ifelse(
      grepl("nis", x),
      "TGFb_blockade",
      ifelse(grepl("dse", x), "DSE_reference", NA_character_)
    )
  )
}

# Verify the explicit map against the author-provided metadata before reading
# biological labels from it. This turns any future repository change into an
# auditable stop rather than a silent remapping.
deposit_check <- merge(
  official,
  deposit_meta,
  by = "library_id",
  all = TRUE,
  suffixes = c("_official", "_deposited"),
  sort = FALSE
)
deposit_check$patient_deposited_normalised <- normalise_patient(
  deposit_check$patient_deposited
)
deposit_check$condition_deposited_normalised <- classify_deposited_treatment(
  deposit_check$treatment_deposited
)
deposit_check$patient_concordant <- with(
  deposit_check,
  patient_official == patient_deposited_normalised
)
deposit_check$condition_concordant <- with(
  deposit_check,
  condition_official == condition_deposited_normalised
)
write_csv(deposit_check, file.path(out, "deposited_metadata_mapping_audit.csv"))

metadata_mapping_pass <- (
  nrow(deposit_check) == 18 &&
    all(!is.na(deposit_check$gsm)) &&
    all(deposit_check$patient_concordant %in% TRUE) &&
    all(deposit_check$condition_concordant %in% TRUE)
)
if (!metadata_mapping_pass) {
  stop(
    "The explicit 18-library GEO map did not agree with the deposited ",
    "GSE160686 metadata. See deposited_metadata_mapping_audit.csv."
  )
}

extract_library_id <- function(nm) {
  nm <- as.character(nm)

  # Actual deposited cell-column format:
  # AS_IL_LIB_..._05.AAACCCAAGTCCGTCG.20
  prefix <- sub("\\..*$", "", nm)
  if (prefix %in% official$library_id) return(prefix)

  # Permit already aggregated matrices named directly by library ID.
  if (nm %in% official$library_id) return(nm)

  # Permit matrices named by GSM accession.
  gsm_hit <- which(vapply(
    official$gsm,
    function(g) grepl(g, nm, fixed = TRUE),
    logical(1)
  ))
  if (length(gsm_hit) == 1) return(official$library_id[gsm_hit])

  # Final safe fallback: one and only one embedded official library ID.
  library_hit <- which(vapply(
    official$library_id,
    function(id) grepl(id, nm, fixed = TRUE),
    logical(1)
  ))
  if (length(library_hit) == 1) return(official$library_id[library_hit])

  NA_character_
}

matrix_library_id <- vapply(
  colnames(counts),
  extract_library_id,
  character(1)
)
idx <- match(matrix_library_id, official$library_id)

library_count_audit <- as.data.frame(
  table(matrix_library_id, useNA = "ifany"),
  stringsAsFactors = FALSE
)
colnames(library_count_audit) <- c("matrix_library_id", "n_matrix_columns")
library_count_audit$matched_gsm <- official$gsm[
  match(library_count_audit$matrix_library_id, official$library_id)
]
write_csv(library_count_audit, file.path(out, "matrix_library_prefix_audit.csv"))

matrix_mapping_pass <- (
  all(!is.na(idx)) &&
    setequal(unique(matrix_library_id), official$library_id) &&
    length(unique(matrix_library_id)) == 18
)
if (!matrix_mapping_pass) {
  manual_mapping_stop(
    library_count_audit,
    file.path(out, "sample_mapping_NEEDS_REVIEW.csv"),
    paste0(
      "Expected every count-matrix column to map to exactly one of the 18 ",
      "official GSE160686 libraries."
    )
  )
}

original_columns <- colnames(counts)
gene_names <- rownames(counts)

# First aggregation level: cell columns -> 18 GEO libraries.
official_index <- seq_len(nrow(official))
library_groups <- lapply(
  official_index,
  function(j) which(idx == j)
)
names(library_groups) <- official$gsm

library_counts <- vapply(
  library_groups,
  function(ii) rowSums(counts[, ii, drop = FALSE]),
  numeric(nrow(counts))
)
rownames(library_counts) <- gene_names
colnames(library_counts) <- official$gsm

library_meta <- official
library_meta$n_cells_aggregated <- lengths(library_groups)
library_meta$matrix_column_examples <- vapply(
  library_groups,
  function(ii) paste(utils::head(original_columns[ii], 3), collapse = ";"),
  character(1)
)
write_csv(library_meta, file.path(out, "sample_mapping_audit.csv"))
write_csv(deposit_meta, file.path(out, "deposited_metadata_copy.csv"))

if (any(library_meta$n_cells_aggregated < 1)) {
  stop("At least one official GSE160686 library contains no mapped cells.")
}

aggregate_libraries <- function(count_matrix, metadata, grouping_columns) {
  grouping_key <- do.call(
    paste,
    c(metadata[grouping_columns], list(sep = "__"))
  )
  grouping_levels <- unique(grouping_key)
  grouping_indices <- split(
    seq_along(grouping_key),
    factor(grouping_key, levels = grouping_levels)
  )

  aggregated_counts <- vapply(
    grouping_indices,
    function(ii) rowSums(count_matrix[, ii, drop = FALSE]),
    numeric(nrow(count_matrix))
  )
  rownames(aggregated_counts) <- rownames(count_matrix)

  aggregated_meta <- metadata[
    vapply(grouping_indices, function(ii) ii[1], integer(1)),
    grouping_columns,
    drop = FALSE
  ]
  aggregated_meta$analysis_sample <- names(grouping_indices)
  aggregated_meta$n_libraries_aggregated <- lengths(grouping_indices)
  aggregated_meta$library_ids <- vapply(
    grouping_indices,
    function(ii) paste(metadata$library_id[ii], collapse = ";"),
    character(1)
  )
  aggregated_meta$gsms <- vapply(
    grouping_indices,
    function(ii) paste(metadata$gsm[ii], collapse = ";"),
    character(1)
  )
  rownames(aggregated_meta) <- NULL
  colnames(aggregated_counts) <- aggregated_meta$analysis_sample

  list(counts = aggregated_counts, meta = aggregated_meta)
}

# Library-level signature scores are retained only for replication/QC display.
# They are not used as independent biological observations.
library_lcpm <- log_cpm(library_counts)
library_sc <- score_gene_set(
  library_lcpm,
  signature,
  ANALYSIS_OPTIONS$minimum_signature_coverage
)
library_score_df <- cbind(
  library_meta,
  score = as.numeric(library_sc$score),
  signature_coverage = library_sc$coverage
)
write_csv(
  library_score_df,
  file.path(out, "ECM_LIGAND_21_library_scores_descriptive.csv")
)
# Backward-compatible filename retained for downstream scripts.
write_csv(
  library_score_df,
  file.path(out, "ECM_LIGAND_21_library_scores.csv")
)
write_csv(
  library_score_df,
  file.path(out, "all_conditions_score_sensitivity.csv")
)

# Primary prespecified contrast: NIS/NIS793 versus matched isotype/IgG2.
# The two libraries within each patient-condition are summed before modelling,
# leaving six pseudobulk samples and three independent paired patients.
primary_library <- library_meta$condition %in% c(
  "Control",
  "TGFb_blockade"
)
primary_pb <- aggregate_libraries(
  library_counts[, primary_library, drop = FALSE],
  library_meta[primary_library, , drop = FALSE],
  c("patient", "condition")
)
cts <- primary_pb$counts
m <- primary_pb$meta

if (
  ncol(cts) != 6 ||
    nrow(m) != 6 ||
    !all(m$n_libraries_aggregated == 2) ||
    !identical(sort(unique(m$patient)), c("P1", "P2", "P3"))
) {
  stop(
    "Primary pseudobulk design must contain 3 patients x 2 conditions, ",
    "with 2 libraries summed in every patient-condition."
  )
}
write_csv(m, file.path(out, "primary_patient_condition_design.csv"))

m$condition <- factor(
  m$condition,
  levels = c("Control", "TGFb_blockade")
)
m$patient <- factor(m$patient, levels = c("P1", "P2", "P3"))
design <- stats::model.matrix(~ patient + condition, m)
if (qr(design)$rank != ncol(design)) {
  stop("The patient-paired primary design matrix is not full rank.")
}
contrast_column <- grep(
  "^conditionTGFb_blockade$",
  colnames(design)
)
if (length(contrast_column) != 1) {
  stop("Could not identify the prespecified TGFb-blockade contrast column.")
}
contrast <- rep(0, ncol(design))
contrast[contrast_column] <- 1

de <- fit_edger(cts, design, contrast)
de$contrast <- "TGFb_blockade_vs_isotype_or_IgG2_control"
de$experimental_unit <- paste0(
  "patient; two within-patient treatment libraries were summed before ",
  "modelling"
)
de$pseudobulk_scope <- paste0(
  "all deposited cells; no cell-type annotation was provided in the ",
  "combined metadata"
)
write_csv(de, file.path(out, "TGFb_blockade_patient_pseudobulk_DEG.csv"))
# Backward-compatible filename; contents now use patient-condition pseudobulk.
write_csv(de, file.path(out, "TGFb_blockade_paired_blocked_DEG.csv"))

lcpm <- log_cpm(cts)
sc <- score_gene_set(
  lcpm,
  signature,
  ANALYSIS_OPTIONS$minimum_signature_coverage
)
save_gene_set_audit(
  sc,
  file.path(out, "ECM_LIGAND_21_coverage.csv"),
  "ECM_LIGAND_21"
)
score_df <- cbind(
  m,
  score = as.numeric(sc$score),
  signature_coverage = sc$coverage
)
write_csv(
  score_df,
  file.path(out, "ECM_LIGAND_21_patient_condition_scores.csv")
)

pe <- paired_effect(
  score_df$score,
  score_df$patient,
  score_df$condition,
  "Control",
  "TGFb_blockade"
)
if (nrow(pe$pairs) != 3 || any(!is.finite(pe$pairs$delta))) {
  stop("The primary score analysis did not yield three finite paired-patient effects.")
}
pe$summary$expected_direction <- "negative"
pe$summary$patients_with_decrease <- sum(pe$pairs$delta < 0)
pe$summary$patient_level_inference_note <- paste0(
  "n=3 independent patients; libraries and cells are not treated as ",
  "independent biological replicates"
)
pe$summary$pseudobulk_scope <- paste0(
  "all deposited cells; interpretation is a whole-sample perturbation ",
  "response rather than a fibroblast-specific estimate"
)
write_csv(
  pe$pairs,
  file.path(out, "ECM_LIGAND_21_patient_paired_differences.csv")
)
write_csv(
  pe$summary,
  file.path(out, "ECM_LIGAND_21_patient_effect_summary.csv")
)

# DSE is a separately labelled concordance/sensitivity analysis and is never
# pooled with NIS/NIS793 in the prespecified primary contrast.
all_pb <- aggregate_libraries(
  library_counts,
  library_meta,
  c("patient", "condition")
)
all_lcpm <- log_cpm(all_pb$counts)
all_sc <- score_gene_set(
  all_lcpm,
  signature,
  ANALYSIS_OPTIONS$minimum_signature_coverage
)
all_scores <- cbind(
  all_pb$meta,
  score = as.numeric(all_sc$score),
  signature_coverage = all_sc$coverage
)
write_csv(
  all_scores,
  file.path(out, "all_conditions_patient_scores_sensitivity.csv")
)

dse_subset <- all_scores$condition %in% c("Control", "DSE_reference")
dse_effect <- paired_effect(
  all_scores$score[dse_subset],
  all_scores$patient[dse_subset],
  all_scores$condition[dse_subset],
  "Control",
  "DSE_reference"
)
if (nrow(dse_effect$pairs) != 3 || any(!is.finite(dse_effect$pairs$delta))) {
  stop("The DSE sensitivity did not yield three finite paired-patient effects.")
}
dse_effect$summary$expected_direction <- "negative"
dse_effect$summary$patients_with_decrease <- sum(
  dse_effect$pairs$delta < 0
)
dse_effect$summary$analysis_role <- paste0(
  "separately labelled DSE concordance sensitivity; not pooled into primary"
)
write_csv(
  dse_effect$pairs,
  file.path(out, "DSE_reference_patient_paired_differences.csv")
)
write_csv(
  dse_effect$summary,
  file.path(out, "DSE_reference_patient_effect_summary.csv")
)

primary_plot_df <- score_df
primary_plot_df$condition <- factor(
  primary_plot_df$condition,
  levels = c("Control", "TGFb_blockade")
)
p <- ggplot2::ggplot(
  primary_plot_df,
  ggplot2::aes(condition, score, group = patient, colour = patient)
) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::labs(
    x = NULL,
    y = "ECM_LIGAND_21 score",
    title = "Human CRC TGF-beta blockade (GSE160686)",
    subtitle = paste0(
      "Two libraries summed within patient-condition; ",
      "patient is the inference unit"
    )
  ) +
  theme_publication()
save_plot(
  p,
  file.path(out, "GSE160686_patient_level_TGFb_response"),
  6.5,
  4.5
)

all_plot_df <- all_scores
all_plot_df$condition <- factor(
  all_plot_df$condition,
  levels = c("Control", "TGFb_blockade", "DSE_reference")
)
p_all <- ggplot2::ggplot(
  all_plot_df,
  ggplot2::aes(condition, score, group = patient, colour = patient)
) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::labs(
    x = NULL,
    y = "ECM_LIGAND_21 score",
    title = "GSE160686 primary and DSE sensitivity conditions",
    subtitle = "DSE is displayed separately and is not pooled into the primary contrast"
  ) +
  theme_publication()
save_plot(
  p_all,
  file.path(out, "GSE160686_all_conditions_patient_level_sensitivity"),
  7.2,
  4.5
)

write_status(
  "05_GSE160686_TGFb",
  "PASS",
  paste0(
    "Primary patient-level decreases: ",
    sum(pe$pairs$delta < 0),
    " of ",
    nrow(pe$pairs),
    "; DSE sensitivity decreases: ",
    sum(dse_effect$pairs$delta < 0),
    " of ",
    nrow(dse_effect$pairs),
    "."
  )
)
