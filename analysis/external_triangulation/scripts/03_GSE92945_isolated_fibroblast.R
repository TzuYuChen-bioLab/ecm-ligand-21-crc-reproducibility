PIPELINE_ROOT <- normalizePath(Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))
require_packages(c("data.table", "edgeR", "limma", "ggplot2"))
out <- result_dir("03_GSE92945")
counts <- numeric_expression_matrix(INPUT_FILES$GSE92945_counts, "sum")
signature <- read_gene_set(SIGNATURE_FILE)
signature_map_file <- file.path(RESOURCE_ROOT, "CRC_ECM_LIGAND_21_ENSEMBL.csv")
assert_file(signature_map_file, "locked ECM_LIGAND_21 Ensembl crosswalk")
signature_map <- read.csv(signature_map_file, check.names = FALSE, stringsAsFactors = FALSE)
signature_map$gene <- clean_symbol(signature_map$gene)
signature_map$ensembl_gene_id <- clean_symbol(signature_map$ensembl_gene_id)
if (
  nrow(signature_map) != 21L ||
  length(unique(signature_map$gene)) != 21L ||
  length(unique(signature_map$ensembl_gene_id)) != 21L ||
  !setequal(signature_map$gene, signature)
) {
  stop("The bundled ECM_LIGAND_21 Ensembl crosswalk does not exactly match the locked 21-gene signature.")
}

# The official supplementary CSV contains two columns per RNA sample:
# an integer raw-count column (e.g. CA1) and a normalized FPKM column
# (e.g. CA1_fpkm). edgeR must receive raw counts only.  Keep a complete
# audit and exclude FPKM columns before matching the ten RNA samples.
all_matrix_columns <- colnames(counts)
is_fpkm <- grepl("(^|[._-])fpkm($|[._-])", all_matrix_columns, ignore.case = TRUE)
column_role_audit <- data.frame(
  matrix_column = all_matrix_columns,
  analysis_role = ifelse(is_fpkm, "excluded_normalized_FPKM", "candidate_raw_count"),
  stringsAsFactors = FALSE
)
write_csv(column_role_audit, file.path(out, "matrix_column_role_audit.csv"))

counts <- counts[, !is_fpkm, drop = FALSE]
if (ncol(counts) != 10L) {
  manual_mapping_stop(
    column_role_audit,
    file.path(out, "matrix_column_role_NEEDS_REVIEW.csv"),
    paste0(
      "Expected 10 raw-count columns after excluding FPKM columns, but found ",
      ncol(counts), "."
    )
  )
}
if (any(!is.finite(counts)) || any(abs(counts - round(counts)) > 1e-8)) {
  stop("GSE92945 candidate raw-count columns are not finite integer-like counts; do not run edgeR.")
}
counts <- round(counts)

# GSE92945 rows are Ensembl stable gene IDs, whereas the locked signature uses
# HGNC symbols. Use the bundled, version-audited 21-gene crosswalk locally;
# do not query an online annotation service during analysis and do not change
# the locked signature or its coverage threshold.
signature_map$detected_in_GSE92945 <- signature_map$ensembl_gene_id %in% rownames(counts)
signature_map$total_raw_count <- NA_real_
present_map <- which(signature_map$detected_in_GSE92945)
signature_map$total_raw_count[present_map] <- base::rowSums(
  counts[signature_map$ensembl_gene_id[present_map], , drop = FALSE]
)
write_csv(
  signature_map,
  file.path(out, "ECM_LIGAND_21_Ensembl_mapping_audit.csv")
)
mapping_coverage <- mean(signature_map$detected_in_GSE92945)
if (mapping_coverage < ANALYSIS_OPTIONS$minimum_signature_coverage) {
  stop(
    sprintf(
      "Mapped ECM_LIGAND_21 coverage %.1f%% is below the pre-specified %.1f%% threshold.",
      100 * mapping_coverage,
      100 * ANALYSIS_OPTIONS$minimum_signature_coverage
    )
  )
}

sample_map <- data.frame(
  gsm = paste0("GSM", 2440822:2440831),
  sample = c("CA1", "CA2", "CA3", "CA4", "CO2", "CO4", "CO5", "N1", "N2", "N3"),
  group = c(rep("Cancer", 4), rep("Colitis", 3), rep("Normal", 3))
)
match_sample <- function(nm) {
  hit <- which(vapply(seq_len(nrow(sample_map)), function(i) grepl(sample_map$gsm[i], nm, ignore.case = TRUE) || grepl(paste0("(^|[^A-Z0-9])", sample_map$sample[i], "([^A-Z0-9]|$)"), nm, ignore.case = TRUE), logical(1)))
  if (length(hit) == 1) hit else NA_integer_
}
map_idx <- vapply(colnames(counts), match_sample, integer(1))
keep <- !is.na(map_idx)
counts <- counts[, keep, drop = FALSE]
meta <- sample_map[map_idx[keep], ]
meta$matrix_column <- colnames(counts)
if (ncol(counts) != 10 || any(table(meta$group) != c(Cancer = 4, Colitis = 3, Normal = 3))) {
  manual_mapping_stop(
    data.frame(matrix_column = colnames(counts), inferred_group = meta$group, inferred_sample = meta$sample),
    file.path(out, "sample_mapping_NEEDS_REVIEW.csv"),
    "Expected 10 RNA samples (4 Cancer, 3 Colitis, 3 Normal)."
  )
}
write_csv(meta, file.path(out, "sample_mapping_audit.csv"))

meta$group <- factor(meta$group, levels = c("Normal", "Colitis", "Cancer"))
design <- stats::model.matrix(~ 0 + group, meta)
colnames(design) <- levels(meta$group)
de_cancer <- fit_edger(counts, design, limma::makeContrasts(Cancer - Normal, levels = design))
de_cancer$contrast <- "Cancer_vs_Normal"
de_colitis <- fit_edger(counts, design, limma::makeContrasts(Colitis - Normal, levels = design))
de_colitis$contrast <- "Colitis_vs_Normal"
annotate_locked_signature <- function(x) {
  x$gene_ensembl <- x$gene
  x$gene_symbol <- signature_map$gene[match(x$gene, signature_map$ensembl_gene_id)]
  x$ECM_LIGAND_21_member <- !is.na(x$gene_symbol)
  x
}
de_cancer <- annotate_locked_signature(de_cancer)
de_colitis <- annotate_locked_signature(de_colitis)
write_csv(rbind(de_cancer, de_colitis), file.path(out, "fibroblast_differential_expression.csv"))

lcpm <- log_cpm(counts)
signature_lcpm <- lcpm[
  signature_map$ensembl_gene_id[signature_map$detected_in_GSE92945],
  ,
  drop = FALSE
]
rownames(signature_lcpm) <- signature_map$gene[signature_map$detected_in_GSE92945]
sc <- score_gene_set(signature_lcpm, signature, ANALYSIS_OPTIONS$minimum_signature_coverage)
save_gene_set_audit(sc, file.path(out, "ECM_LIGAND_21_coverage.csv"), "ECM_LIGAND_21")
score_df <- cbind(meta, score = as.numeric(sc$score), signature_coverage = sc$coverage)
write_csv(score_df, file.path(out, "ECM_LIGAND_21_sample_scores.csv"))

contrasts <- lapply(c("Cancer", "Colitis"), function(g) {
  tt <- stats::t.test(score_df$score[score_df$group == g], score_df$score[score_df$group == "Normal"])
  data.frame(
    contrast = paste(g, "vs Normal"),
    n_group = sum(score_df$group == g), n_normal = sum(score_df$group == "Normal"),
    mean_difference = mean(score_df$score[score_df$group == g]) - mean(score_df$score[score_df$group == "Normal"]),
    ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], p_value = tt$p.value,
    test = "Welch_t_test_on_prespecified_sample_level_score"
  )
})
contrasts <- do.call(rbind, contrasts)
contrasts$FDR <- bh(contrasts$p_value)
write_csv(contrasts, file.path(out, "ECM_LIGAND_21_group_contrasts.csv"))

p <- ggplot2::ggplot(score_df, ggplot2::aes(group, score, colour = group)) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.08, size = 2.5) +
  ggplot2::scale_colour_manual(values = c(Normal = "#2C73B9", Colitis = "#E69F00", Cancer = "#C9223A")) +
  ggplot2::labs(x = NULL, y = "ECM_LIGAND_21 score", title = "Isolated colon fibroblasts (GSE92945)") +
  theme_publication() + ggplot2::theme(legend.position = "none")
save_plot(p, file.path(out, "GSE92945_isolated_fibroblast_score"), 5.5, 4.5)
write_status("03_GSE92945", "PASS", "Cancer=4, Colitis=3, Normal=3 independent fibroblast samples.")
