PIPELINE_ROOT <- normalizePath(
  Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))
require_packages(c("data.table"))

out <- result_dir("11_manuscript_update")

ledger_file <- file.path(
  RESULT_ROOT,
  "10_integrated_evidence",
  "mechanistic_evidence_ledger.csv"
)

fit_file <- file.path(
  RESULT_ROOT,
  "10_integrated_evidence",
  "FIG_internal_fit_assessment_NOT_FOR_SUBMISSION.csv"
)

assert_file(ledger_file)
assert_file(fit_file)

ledger <- read.csv(ledger_file, check.names = FALSE)
fit <- read.csv(fit_file, check.names = FALSE)

required_ledger_columns <- c(
  "evidence_id",
  "analysis",
  "result_class",
  "observed",
  "claim_to_use",
  "limitations"
)

if (!all(required_ledger_columns %in% names(ledger))) {
  stop(
    "The Step 10 ledger is not the corrected graded-evidence version. Missing: ",
    paste(
      setdiff(required_ledger_columns, names(ledger)),
      collapse = ", "
    )
  )
}

if (!"conservative_internal_fit_score" %in% names(fit)) {
  stop("The corrected internal fit-assessment file was not found.")
}

methods <- c(
  "# Supplementary Methods update: external mechanistic triangulation",
  "",
  "All study data were downloaded manually from their primary repositories and analysed locally; no analysis script downloaded study data. The equal-weight ECM_LIGAND_21 definition was locked before survival analysis and before the external analyses described below. External analyses were interpreted through prespecified, dataset-specific rules that kept technical completion, statistical support, directional concordance, descriptive coverage and non-confirmatory results analytically distinct.",
  "",
  "## Independent paired single-cell analysis",
  "",
  "GSE144735 was analysed as an independent paired single-cell cohort. Raw UMI counts were aggregated within patient, anatomical condition and major cell class, and tumour-core versus normal-mucosa effects were estimated at the patient-pseudobulk level. Cells were not treated as independent biological replicates. A non-overlapping CAF proxy was derived from fibroblast-versus-epithelial markers without using survival outcomes.",
  "",
  "## Isolated-fibroblast analysis",
  "",
  "GSE92945 RNA-seq counts from isolated colon fibroblasts were analysed with edgeR. Cancer-associated (n=4), colitis-associated (n=3) and normal (n=3) fibroblast samples were retained, whereas non-RNA-sequencing columns were excluded. Cancer versus normal was the prespecified primary contrast; colitis versus normal was retained as an inflammatory-context comparison.",
  "",
  "## Spatial transcriptomic analysis",
  "",
  "For GSE280315 Visium HD data, the locked ECM score was calculated in bins assigned fibroblast/stromal labels and the frozen epithelial target score in epithelial/tumour-labelled bins. For each epithelial bin, the mean score among the five nearest fibroblast bins was calculated. Spearman associations and 1,000 within-specimen score-label permutations were evaluated separately for each of the three CRC specimens (P1, P2 and P5). Spatial bins were not treated as patients, and the Fisher-combined permutation result was designated exploratory. A strict sensitivity analysis excluded pericyte-labelled bins from the fibroblast definition.",
  "",
  "## Public perturbation analyses",
  "",
  "For GSE160686, the primary contrast compared NIS/NIS793 TGF-beta blockade with matched isotype/IgG2 controls. Replicate libraries were aggregated within patient and condition, leaving three patients as the inferential units. DSE-labelled samples were retained as a separately reported directional sensitivity analysis and were not pooled with the primary blockade group. Because the deposited matrix comprised all captured cells, the result was interpreted as a whole-sample perturbation response rather than a fibroblast-specific effect.",
  "",
  "GSE162561 was retained as secondary experimental context. Raw RNA-sequencing counts from two primary CMS2 colorectal-cancer cell models (#8 and #9), exposed for 48 hours to control medium, subcutaneous adipose-stromal-cell conditioned medium or visceral adipose-stromal-cell conditioned medium, were summarized within model and condition. Deposited replicate libraries were not treated as additional independent patients; with two cell models, this analysis was interpreted descriptively.",
  "",
  "GSE155343 author-provided rLog expression profiles were used for descriptive coculture comparisons. Deposited replicate profiles were summarized within cell line and culture condition. Because rLog values are not raw counts and the profiles do not constitute an independent patient cohort, no count-based differential-expression or patient-level confirmatory inference was assigned to this analysis.",
  "",
  "## Decellularized-ECM proteomics",
  "",
  "The independent decellularized-ECM proteomics workbook was analysed using author-normalized intensities. Technical columns were collapsed within patient-condition, and paired tumour-versus-normal effects were estimated with patient blocking. Coverage and both enrichment directions were reported for all detected ECM_LIGAND_21 proteins. Protein detectability was not treated as evidence of tumour-direction concordance.",
  "",
  "## CAF-adjusted survival and matched-random-set analyses",
  "",
  "In GSE39582, survival models containing the frozen score were compared with corresponding base models containing the prespecified clinical covariates and the independently derived, non-overlapping CAF_PROXY_30 score. Missing clinical tokens were normalized before complete-case model construction, and all prespecified cohort splits were reported rather than selecting the smallest P value.",
  "",
  "The specificity benchmark generated 10,000 without-replacement random 21-gene sets. Each signature gene was matched from its 200 nearest genes on scaled discovery-reference expression mean and log standard deviation; ECM_LIGAND_21 and CAF_PROXY_30 genes were excluded, and outcomes were not used for matching. The two-sided add-one empirical P value was calculated as [1 + number of random sets with an absolute Wald statistic at least as large as the observed statistic] / [1 + number of estimable random sets]."
)

result_lines <- c(
  "# Results update: graded external-evidence ledger",
  ""
)

for (i in seq_len(nrow(ledger))) {
  result_lines <- c(
    result_lines,
    paste0("## ", ledger$evidence_id[i], ". ", ledger$analysis[i]),
    paste0("Evidence classification: ", ledger$result_class[i], "."),
    paste0("Observed: ", ledger$observed[i]),
    paste0("Permitted wording: ", ledger$claim_to_use[i]),
    paste0("Qualification: ", ledger$limitations[i]),
    ""
  )
}

result_lines <- c(
  result_lines,
  "## Integrated interpretation",
  "",
  "Across the prespecified external analyses, the strongest support arose from an independent paired single-cell cohort, in which the locked score was increased in fibroblast pseudobulk profiles from tumour relative to matched normal specimens. Spatial analysis showed concordant but exploratory associations across three CRC specimens. Public TGF-beta-blockade data provided modest directional concordance in two of three patients, but the estimate was imprecise and represented a whole-sample rather than fibroblast-specific response. Although most locked proteins were detectable in independent decellularized-ECM proteomics, the predominant direction was normal-enriched rather than tumour-enriched, supporting protein-level detectability but not a concordant tumour-associated direction. Isolated-fibroblast analysis, CAF-adjusted survival models and the matched-random-gene-set benchmark were non-confirmatory. Collectively, these findings support a context-dependent stromal-epithelial transcriptional program but do not establish a fibroblast-intrinsic causal mechanism, gene-set-specific prognostic value or independent clinical utility.",
  ""
)

discussion <- c(
  "# Discussion claim boundaries",
  "",
  "The results support use of the terms 'external transcriptomic replication', 'exploratory spatial concordance', 'directionally supportive secondary perturbation analysis' and 'context-dependent stromal-epithelial program' only where the corresponding evidence classification permits them.",
  "",
  "Do not use 'demonstrates causality', 'proves a stromal-to-epithelial mechanism', 'experimentally validated by us', 'CAF-specific TGF-beta response', 'independent prognostic biomarker', 'clinically useful signature' or 'outperformed random gene sets'. CellChat and NicheNet remain hypothesis-prioritisation tools because no displayed epithelial-receiver axis or ligand achieved the paired-donor false-discovery-rate threshold in the discovery cohort.",
  "",
  "The TGF-beta, conditioned-medium, coculture and ECM-proteomic experiments were performed by the original data-generating investigators; the present study contributed secondary computational reanalysis. Spatial association does not establish signal direction or causality. Discordant protein directions and all negative external tests should remain visible in the Results and Discussion rather than being relegated solely to source-data files.",
  "",
  "The absence of CAF-independent prognostic information and the failure to outperform matched random gene sets constrain interpretation of ECM_LIGAND_21 as a mechanistically anchored transcriptional descriptor rather than a validated clinical predictor."
)

figure_legend <- c(
  "# Supplementary figure legend",
  "",
  "**Supplementary Fig. Sx External evidence triangulation and claim-control audit.** Each row represents one external or contextual analysis and is classified as supported, supported with qualification, directionally supportive, descriptive only, contextual only or not supported. Classification incorporated the biological unit of replication, statistical uncertainty, direction of effect and the distinction between newly generated experiments and secondary analysis of public experimental data. The independent paired single-cell cohort provided the strongest replication evidence; spatial evidence was retained with specimen-level and causal qualifications; the TGF-beta analysis was directional; proteomics provided descriptive detectability with predominantly normal-enriched direction; and the isolated-fibroblast, CAF-adjusted survival and matched-random-set analyses were non-confirmatory. Exact estimates and limitations are provided in the accompanying source-data ledger. No journal-fit score is displayed because the internal rubric is not a scientific outcome or acceptance probability."
)

availability <- c(
  "# Data and code availability draft",
  "",
  "All datasets analysed in this study were obtained from public repositories or published supplementary material. GEO accession numbers used in the complete project include GSE132465, GSE144735, GSE92945, GSE280315/GSE280318, GSE160686, GSE162561, GSE155343, GSE39582, GSE38832, GSE14333, GSE17536, GSE17537 and GSE33113. The decellularized-ECM proteomic data were obtained from the supplementary workbook accompanying the original publication. The analysis scripts do not download study data; repository files were downloaded manually and retained under their original filenames.",
  "",
  "All analysis scripts, frozen gene-set definitions, machine-readable figure source data, audit tables and a software-session record will be deposited at [PUBLIC CODE REPOSITORY URL] and archived at [PERSISTENT DOI] before submission. Only derived, non-identifiable data permitted for redistribution will be included. The original data should be accessed from the cited primary repositories."
)

internal_fit <- c(
  "INTERNAL DOCUMENT - DO NOT SUBMIT",
  "",
  paste0(
    "Conservative internal Computational Biology and Chemistry fit estimate: ",
    fit$conservative_internal_fit_score[1],
    "/100."
  ),
  paste0("Interpretation: ", fit$fit_band[1]),
  "This subjective planning score is not an acceptance probability or scientific result and must not appear in the manuscript, figures, supplementary material or cover letter."
)

writeLines(
  methods,
  file.path(out, "Supplementary_Methods_external_validation.md"),
  useBytes = TRUE
)

writeLines(
  result_lines,
  file.path(out, "Results_evidence_blocks.md"),
  useBytes = TRUE
)

writeLines(
  discussion,
  file.path(out, "Discussion_claim_boundaries.md"),
  useBytes = TRUE
)

writeLines(
  figure_legend,
  file.path(out, "Supplementary_Figure_external_evidence_legend.md"),
  useBytes = TRUE
)

writeLines(
  availability,
  file.path(out, "Data_and_Code_Availability_draft.md"),
  useBytes = TRUE
)

writeLines(
  internal_fit,
  file.path(out, "INTERNAL_FIG_fit_assessment_DO_NOT_SUBMIT.txt"),
  useBytes = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out, "R_sessionInfo.txt"),
  useBytes = TRUE
)

files <- list.files(
  RESULT_ROOT,
  recursive = TRUE,
  full.names = TRUE
)
info <- file.info(files)

inventory <- data.frame(
  relative_path = substring(
    files,
    nchar(RESULT_ROOT) + 2
  ),
  size_bytes = info$size,
  modified_utc = format(
    info$mtime,
    tz = "UTC",
    usetz = TRUE
  ),
  stringsAsFactors = FALSE
)

write_csv(
  inventory,
  file.path(out, "result_file_inventory.csv")
)

write_status(
  "11_export_manuscript_blocks",
  "PASS",
  paste0(
    "Graded evidence manuscript blocks exported; inventory files: ",
    nrow(inventory),
    "; journal-fit score confined to internal do-not-submit file."
  )
)
