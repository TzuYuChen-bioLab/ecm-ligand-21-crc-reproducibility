# CRC ECM mechanistic-triangulation extension
# Edit only this file before running the analysis.

PROJECT_ROOT <- Sys.getenv(
  "CRC_PROJECT_ROOT",
  unset = getwd()
)

PIPELINE_ROOT <- Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd())
PIPELINE_ROOT <- normalizePath(PIPELINE_ROOT, winslash = "/", mustWork = TRUE)
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)

DATA_ROOT <- file.path(PIPELINE_ROOT, "data")
RESULT_ROOT <- file.path(PIPELINE_ROOT, "results")
RESOURCE_ROOT <- file.path(PIPELINE_ROOT, "resources")
LOG_ROOT <- file.path(PIPELINE_ROOT, "logs")

dir.create(RESULT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_ROOT, recursive = TRUE, showWarnings = FALSE)

ANALYSIS_OPTIONS <- list(
  minimum_cells_per_pseudobulk = 20L,
  minimum_signature_coverage = 0.70,
  minimum_target_coverage = 0.60,
  spatial_nearest_fibroblast_bins = 5L,
  spatial_max_epithelial_bins_per_section = 20000L,
  spatial_permutations = 1000L,
  random_gene_sets = 10000L,
  random_seed = 20260721L,
  run_optional_GSE162561 = TRUE,
  run_optional_GSE155343 = TRUE
)

INPUT_FILES <- list(
  GSE144735_counts = file.path(DATA_ROOT, "GSE144735", "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"),
  GSE144735_annotation = file.path(DATA_ROOT, "GSE144735", "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"),
  GSE92945_counts = file.path(DATA_ROOT, "GSE92945", "GSE92945_Fibroblast_RNAseq_counts.csv.gz"),
  GSE160686_counts = file.path(DATA_ROOT, "GSE160686", "GSE160686_combined_count_matrix_TGFb.csv.gz"),
  GSE160686_metadata = file.path(DATA_ROOT, "GSE160686", "GSE160686_combined_metadata_TGFb.csv.gz"),
  GSE162561_counts = file.path(DATA_ROOT, "GSE162561", "GSE162561_raw_counts.txt.gz"),
  GSE155343_expression = file.path(DATA_ROOT, "GSE155343", "GSE155343_DESeq2_rLog_Normalized_counts.txt.gz"),
  ECM_proteomics = file.path(DATA_ROOT, "ECM_proteomics", "41416_2025_2964_MOESM3_ESM.xlsx")
)

SPATIAL_SAMPLES <- data.frame(
  section = c("P1CRC", "P2CRC", "P5CRC"),
  h5 = file.path(
    DATA_ROOT,
    "GSE280315",
    c(
      "GSM8594567_P1CRC_filtered_feature_bc_matrix.h5",
      "GSM8594568_P2CRC_filtered_feature_bc_matrix.h5",
      "GSM8594569_P5CRC_filtered_feature_bc_matrix.h5"
    )
  ),
  metadata = file.path(
    DATA_ROOT,
    "GSE280315",
    c(
      "GSM8594567_P1CRC_Metadata.parquet.gz",
      "GSM8594568_P2CRC_Metadata.parquet.gz",
      "GSM8594569_P5CRC_Metadata.parquet.gz"
    )
  ),
  stringsAsFactors = FALSE
)

# Pre-specified label patterns. Inspect the exported label audit before accepting
# spatial results; do not change these after viewing outcome associations.
SPATIAL_LABEL_REGEX <- list(
  fibroblast = "fibro|stromal|myofibro|caf|pericyte",
  epithelial = "epithelial|tumou?r|malignant|colonocyte|goblet|stem|enterocyte"
)

# Local raw GSE39582 series matrix. The script searches these folders if blank.
GSE39582_LOCAL_MATRIX <- ""

SIGNATURE_FILE <- file.path(RESOURCE_ROOT, "CRC_ECM_LIGAND_21.csv")
TARGET_FILE <- file.path(RESOURCE_ROOT, "EPITHELIAL_TUMOR_UP_50.csv")
