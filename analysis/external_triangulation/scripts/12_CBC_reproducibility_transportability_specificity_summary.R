#!/usr/bin/env Rscript
# CBC framing audit: reproducibility, transportability, and specificity
# This script uses only supplied, precomputed source-data files. It does not
# download data and does not replace the primary analysis pipelines.

options(stringsAsFactors = FALSE)

script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit)) return(dirname(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE)))
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

PACKAGE_ROOT <- normalizePath(file.path(script_dir(), "..", "..", ".."), winslash = "/", mustWork = TRUE)
default_source <- file.path(dirname(PACKAGE_ROOT), "ESM_4_Figure_Source_Data_and_Redraw_CBC_FINAL", "source_data")
SOURCE_ROOT <- Sys.getenv("CBC_SOURCE_ROOT", unset = default_source)
if (!dir.exists(SOURCE_ROOT)) {
  stop("Could not locate the ESM 4 source_data folder. Set CBC_SOURCE_ROOT to that folder.")
}
OUT_DIR <- file.path(PACKAGE_ROOT, "results", "generated", "CBC_summary_outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_required <- function(rel) {
  path <- file.path(SOURCE_ROOT, rel)
  if (!file.exists(path)) stop("Missing required source file: ", path)
  read.csv(path, check.names = FALSE)
}

paired <- read_required("Figure_06/ECM_LIGAND_21_paired_effects.csv")
isolated <- read_required("Figure_06/ECM_LIGAND_21_group_contrasts.csv")
spatial <- read_required("Shared/spatial/primary_vs_strict_no_pericytes_comparison.csv")
perturb <- read_required("Figure_06/ECM_LIGAND_21_patient_effect_summary.csv")
protein <- read_required("Figure_06/ECM_LIGAND_21_protein_effect_summary.csv")
random_summary <- read_required("Figure_07/matched_random_gene_set_benchmark_summary.csv")
ledger_path <- file.path(SOURCE_ROOT, "Figure_01", "Figure1_guardrails_source.csv")
ledger <- if (file.exists(ledger_path)) read.csv(ledger_path, check.names = FALSE) else data.frame(
  note = "Claim-control logic is specified in the manuscript, ESM 1, and pipeline scripts.",
  stringsAsFactors = FALSE
)

write.csv(paired, file.path(OUT_DIR, "01_reproducibility_paired_effects.csv"), row.names = FALSE)
write.csv(isolated, file.path(OUT_DIR, "02a_transportability_isolated_fibroblast.csv"), row.names = FALSE)
write.csv(spatial, file.path(OUT_DIR, "02b_transportability_spatial.csv"), row.names = FALSE)
write.csv(perturb, file.path(OUT_DIR, "02c_transportability_TGFb_perturbation.csv"), row.names = FALSE)
write.csv(protein, file.path(OUT_DIR, "02d_transportability_ECM_proteomics.csv"), row.names = FALSE)
write.csv(random_summary, file.path(OUT_DIR, "03_specificity_random_benchmark.csv"), row.names = FALSE)
write.csv(ledger, file.path(OUT_DIR, "04_claim_control_ledger.csv"), row.names = FALSE)

# Machine-readable manifest for manuscript cross-checking.
manifest <- data.frame(
  domain = c("reproducibility", "transportability", "specificity", "claim_control"),
  output = c("01_reproducibility_paired_effects.csv", "02a-02d_transportability_outputs",
             "03_specificity_random_benchmark.csv", "04_claim_control_ledger.csv"),
  interpretation_boundary = c(
    "Patient-level directional replication; not secretion or causality",
    "Cross-context evidence; effect magnitudes are not pooled across modalities",
    "GSE39582 matched-random benchmark; not a universal clinical validation",
    "Categorical claim restriction; not a weighted statistical meta-analysis"
  )
)
write.csv(manifest, file.path(OUT_DIR, "CBC_summary_manifest.csv"), row.names = FALSE)
message("CBC summary outputs written to: ", normalizePath(OUT_DIR, winslash = "/"))
