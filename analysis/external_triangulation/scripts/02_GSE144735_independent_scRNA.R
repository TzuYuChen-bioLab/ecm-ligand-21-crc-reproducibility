PIPELINE_ROOT <- normalizePath(Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))
require_packages(c("data.table", "Matrix", "edgeR", "limma", "ggplot2"))
set.seed(ANALYSIS_OPTIONS$random_seed)
out <- result_dir("02_GSE144735")

assert_file(INPUT_FILES$GSE144735_counts)
assert_file(INPUT_FILES$GSE144735_annotation)
signature <- read_gene_set(SIGNATURE_FILE)

message("Reading GSE144735 raw UMI matrix. A 64-bit R session and substantial RAM are recommended.")
raw <- read_table_auto(INPUT_FILES$GSE144735_counts)
gene_col <- guess_gene_column(raw)
genes <- clean_symbol(raw[[gene_col]])
raw[[gene_col]] <- NULL
raw[] <- lapply(raw, function(z) suppressWarnings(as.numeric(z)))
counts <- as.matrix(raw)
rm(raw)
gc()

# Detect orientation before collapsing duplicate genes.
annotation <- read_table_auto(INPUT_FILES$GSE144735_annotation)
if (nrow(counts) == nrow(annotation) && ncol(counts) != nrow(annotation)) {
  counts <- t(counts)
  genes <- clean_symbol(rownames(counts))
}
if (length(genes) != nrow(counts)) stop("Could not establish gene-by-cell orientation for GSE144735.")
rownames(counts) <- make.unique(genes)
counts <- Matrix::Matrix(counts, sparse = TRUE)

cell_col <- best_matching_column(annotation, colnames(counts))
if (is.na(cell_col)) {
  # Some deposits store cell identifiers as annotation row names or in the first column.
  candidate <- names(annotation)[1]
  if (nrow(annotation) == ncol(counts)) cell_col <- candidate
  else manual_mapping_stop(
    data.frame(column = names(annotation), example = vapply(annotation, function(z) paste(head(z, 3), collapse = " | "), character(1))),
    file.path(out, "annotation_column_audit_NEEDS_REVIEW.csv"),
    "Cell identifiers in annotation could not be matched to count-matrix columns."
  )
}
cell_ids <- as.character(annotation[[cell_col]])
idx <- match(colnames(counts), cell_ids)
if (sum(!is.na(idx)) < 0.90 * ncol(counts)) {
  if (nrow(annotation) == ncol(counts)) idx <- seq_len(nrow(annotation))
  else stop("Fewer than 90% of count-matrix cells matched the annotation.")
}
annotation <- annotation[idx, , drop = FALSE]

patient_col <- first_name_matching(annotation, c("patient", "donor", "subject", "sample"), FALSE)
condition_col <- first_name_matching(annotation, c("condition", "tissue", "region", "sample", "patient"), FALSE)
type_col <- first_name_matching(annotation, c("cell.?type", "celltype", "annotation", "cluster", "subtype"), FALSE)

if (is.na(patient_col)) {
  token <- sub("[-_].*$", "", cell_ids[idx])
  annotation$inferred_patient <- token
  patient_col <- "inferred_patient"
}
if (is.na(condition_col)) {
  annotation$inferred_condition <- standardize_condition(cell_ids[idx])
  condition_col <- "inferred_condition"
}
if (is.na(type_col)) {
  manual_mapping_stop(
    data.frame(column = names(annotation), example = vapply(annotation, function(z) paste(head(z, 3), collapse = " | "), character(1))),
    file.path(out, "annotation_column_audit_NEEDS_REVIEW.csv"),
    "Cell-type column could not be inferred."
  )
}

meta <- data.frame(
  cell = colnames(counts),
  patient = as.character(annotation[[patient_col]]),
  condition_raw = as.character(annotation[[condition_col]]),
  cell_type_raw = as.character(annotation[[type_col]])
)
meta$condition <- standardize_condition(meta$condition_raw)
meta$cell_class <- standardize_major_type(meta$cell_type_raw)

# Recover KUL patient and T/B/N region from any informative string if needed.
combined_text <- apply(annotation, 1, paste, collapse = "|")
has_kul <- grepl("KUL[0-9]+", combined_text, ignore.case = TRUE)
kul <- rep(NA_character_, length(combined_text))
kul[has_kul] <- toupper(sub(".*(KUL[0-9]+).*", "\\1", combined_text[has_kul], ignore.case = TRUE))
meta$patient[has_kul] <- kul[has_kul]
meta$patient <- sub("[-_](T|B|N)$", "", meta$patient, ignore.case = TRUE)
fallback_condition <- standardize_condition(combined_text)
meta$condition[is.na(meta$condition)] <- fallback_condition[is.na(meta$condition)]

write_csv(
  data.frame(
    inferred_role = c("cell_id", "patient", "condition", "cell_type"),
    source_column = c(cell_col, patient_col, condition_col, type_col)
  ),
  file.path(out, "annotation_mapping_audit.csv")
)
write_csv(as.data.frame(table(meta$patient, meta$condition, meta$cell_class)), file.path(out, "cell_count_audit.csv"))

keep <- !is.na(meta$condition) & !is.na(meta$patient) & meta$cell_class %in% c("Fibroblast", "Epithelial")
counts <- counts[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
sample_id <- paste(meta$patient, meta$condition, meta$cell_class, sep = "|")
pb <- pseudobulk_counts(counts, sample_id, ANALYSIS_OPTIONS$minimum_cells_per_pseudobulk)
pb_meta <- do.call(rbind, strsplit(colnames(pb$counts), "\\|", fixed = FALSE))
pb_meta <- data.frame(
  sample_id = colnames(pb$counts),
  patient = pb_meta[, 1], condition = pb_meta[, 2], cell_class = pb_meta[, 3]
)
pb_meta <- merge(pb_meta, pb$n_cells, by = "sample_id", sort = FALSE)
pb_meta <- pb_meta[match(colnames(pb$counts), pb_meta$sample_id), ]
write_csv(pb_meta, file.path(out, "pseudobulk_design_and_cell_counts.csv"))

all_de <- list()
score_rows <- list()
effect_rows <- list()
for (cell_class in c("Fibroblast", "Epithelial")) {
  for (comparison in c("Tumor", "Border")) {
    use <- pb_meta$cell_class == cell_class & pb_meta$condition %in% c("Normal", comparison)
    m <- pb_meta[use, , drop = FALSE]
    tab <- table(m$patient, m$condition)
    if (!all(c("Normal", comparison) %in% colnames(tab))) next
    complete <- rownames(tab)[tab[, "Normal"] > 0 & tab[, comparison] > 0]
    use <- use & pb_meta$patient %in% complete
    m <- pb_meta[use, , drop = FALSE]
    cts <- pb$counts[, use, drop = FALSE]
    if (length(unique(m$patient)) < 3) next
    m$condition <- factor(m$condition, levels = c("Normal", comparison))
    m$patient <- factor(m$patient)
    design <- stats::model.matrix(~ patient + condition, m)
    contrast <- rep(0, ncol(design)); contrast[grep("^condition", colnames(design))] <- 1
    de <- fit_edger(cts, design, contrast)
    de$cell_class <- cell_class
    de$comparison <- paste(comparison, "vs Normal")
    de$experimental_unit <- "paired_patient_pseudobulk"
    all_de[[paste(cell_class, comparison)]] <- de

    lcpm <- log_cpm(cts)
    sc <- score_gene_set(lcpm, signature, ANALYSIS_OPTIONS$minimum_signature_coverage)
    sr <- data.frame(
      patient = m$patient,
      condition = m$condition,
      cell_class = cell_class,
      comparison = comparison,
      score = as.numeric(sc$score),
      signature_coverage = sc$coverage
    )
    score_rows[[paste(cell_class, comparison)]] <- sr
    pe <- paired_effect(sr$score, sr$patient, sr$condition, "Normal", comparison)
    eff <- pe$summary
    eff$cell_class <- cell_class
    eff$effect_size_paired_d <- cohen_d_paired(pe$pairs$delta)
    effect_rows[[paste(cell_class, comparison)]] <- eff
  }
}
de_all <- do.call(rbind, all_de)
scores <- do.call(rbind, score_rows)
effects <- do.call(rbind, effect_rows)
write_csv(de_all, file.path(out, "paired_pseudobulk_differential_expression.csv"))
write_csv(scores, file.path(out, "ECM_LIGAND_21_patient_scores.csv"))
write_csv(effects, file.path(out, "ECM_LIGAND_21_paired_effects.csv"))

# Derive a non-overlapping fibroblast abundance proxy without any survival data.
use <- pb_meta$condition %in% c("Normal", "Tumor")
m <- pb_meta[use, , drop = FALSE]
cts <- pb$counts[, use, drop = FALSE]
m$cell_class <- factor(m$cell_class, levels = c("Epithelial", "Fibroblast"))
m$patient <- factor(m$patient)
m$condition <- factor(m$condition)
design <- stats::model.matrix(~ patient + condition + cell_class, m)
contrast <- rep(0, ncol(design)); contrast[grep("cell_classFibroblast", colnames(design))] <- 1
proxy_de <- fit_edger(cts, design, contrast)
proxy_de <- proxy_de[order(proxy_de$FDR, -proxy_de$logFC), ]
proxy <- proxy_de[proxy_de$FDR < 0.05 & proxy_de$logFC > 1 & !proxy_de$gene %in% signature, ]
proxy <- head(proxy, 30)
proxy$selection_rule <- "fibroblast_vs_epithelial_FDR_lt_0.05_logFC_gt_1_excluding_ECM_LIGAND_21"
write_csv(proxy_de, file.path(out, "fibroblast_vs_epithelial_marker_statistics.csv"))
write_csv(proxy, file.path(out, "CAF_PROXY_30_nonoverlapping.csv"))

p <- ggplot2::ggplot(scores[scores$comparison == "Tumor", ], ggplot2::aes(condition, score, group = patient)) +
  ggplot2::geom_line(alpha = 0.55, colour = "grey45") +
  ggplot2::geom_point(ggplot2::aes(colour = condition), size = 2.4) +
  ggplot2::facet_wrap(~ cell_class, scales = "free_y") +
  ggplot2::scale_colour_manual(values = c(Normal = "#2C73B9", Tumor = "#C9223A")) +
  ggplot2::labs(x = NULL, y = "ECM_LIGAND_21 score", title = "Independent paired GSE144735 pseudobulk validation") +
  theme_publication()
save_plot(p, file.path(out, "GSE144735_paired_ECM_score"), 8, 4.5)

write_status("02_GSE144735", "PASS", paste("Pseudobulk groups:", nrow(pb_meta)))
message("GSE144735 completed. Review annotation_mapping_audit.csv and cell_count_audit.csv before interpreting effects.")
