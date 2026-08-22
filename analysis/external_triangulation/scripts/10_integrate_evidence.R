PIPELINE_ROOT <- normalizePath(
  Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))
require_packages(c("data.table", "ggplot2"))

out <- result_dir("10_integrated_evidence")

# This ledger intentionally separates:
#   1) technical completion,
#   2) statistically supported evidence,
#   3) directionally supportive evidence,
#   4) descriptive/contextual evidence, and
#   5) non-confirmatory evidence.
# Journal-fit scoring is retained for internal planning only and is never plotted
# or exported as manuscript text.

ledger <- data.frame(
  evidence_id = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
  analysis = c(
    "GSE144735 independent paired single-cell cohort",
    "GSE92945 isolated fibroblasts",
    "GSE280315 spatial transcriptomics",
    "GSE160686 TGF-beta perturbation",
    "Decellularized-ECM proteomics",
    "GSE39582 CAF-adjusted survival model",
    "Matched random 21-gene sets",
    "GSE162561 stromal conditioned-medium response",
    "GSE155343 coculture expression profiles"
  ),
  evidence_role = c(
    "external paired-patient replication",
    "fibroblast-intrinsic specificity",
    "exploratory spatial concordance",
    "directional perturbational concordance",
    "protein-level detectability and direction",
    "CAF-independent prognostic information",
    "gene-set specificity benchmark",
    "secondary experimental context",
    "descriptive coculture context"
  ),
  prespecified_rule = c(
    "Paired fibroblast tumour-versus-normal score difference is positive with P < 0.10",
    "Cancer-versus-normal isolated-fibroblast score difference is positive with P < 0.05",
    "At least two of three CRC specimens have positive rho and the exploratory combined permutation P < 0.05",
    "The score decreases in at least two of three paired patients after TGF-beta blockade",
    "Report coverage and both protein directions without treating detectability as directional validation",
    "The prespecified score adds information beyond clinical covariates and a non-overlapping CAF proxy",
    "The observed prognostic statistic exceeds expression- and variability-matched random 21-gene sets",
    "Report the completed secondary analysis without using it as an inferential gate",
    "Report averaged deposited replicate profiles descriptively without patient-level inference"
  ),
  maximum_internal_points = c(4, 3, 4, 4, 3, 3, 3, 0, 0),
  stringsAsFactors = FALSE
)

ledger$evaluated <- FALSE
ledger$criterion_met <- FALSE
ledger$result_class <- "Not evaluated"
ledger$observed <- NA_character_
ledger$source_file <- NA_character_
ledger$earned_internal_points <- 0
ledger$claim_to_use <- "not evaluated"
ledger$limitations <- NA_character_

set_evidence <- function(
  id,
  result_class,
  criterion_met,
  observed,
  source,
  earned_points,
  claim,
  limitations
) {
  row <- which(ledger$evidence_id == id)

  if (length(row) != 1) {
    stop("Evidence ID must identify exactly one row: ", id)
  }

  if (!earned_points %in% 0:ledger$maximum_internal_points[row]) {
    stop("Invalid internal point assignment for evidence ", id)
  }

  ledger$evaluated[row] <<- TRUE
  ledger$criterion_met[row] <<- isTRUE(criterion_met)
  ledger$result_class[row] <<- as.character(result_class)
  ledger$observed[row] <<- as.character(observed)
  ledger$source_file[row] <<- as.character(source)
  ledger$earned_internal_points[row] <<- as.numeric(earned_points)
  ledger$claim_to_use[row] <<- as.character(claim)
  ledger$limitations[row] <<- as.character(limitations)
}

fmt_p <- function(x) {
  if (!length(x) || !is.finite(x[1])) return("not available")
  format.pval(x[1], digits = 4, eps = 1e-04)
}

fmt_n <- function(x, digits = 3) {
  if (!length(x) || !is.finite(x[1])) return("not available")
  formatC(x[1], format = "f", digits = digits)
}

get_first <- function(x, column, default = NA_real_) {
  if (!column %in% names(x) || !nrow(x)) return(default)
  x[[column]][1]
}

# A: independent paired-patient single-cell replication.
f <- file.path(
  RESULT_ROOT,
  "02_GSE144735",
  "ECM_LIGAND_21_paired_effects.csv"
)

if (file.exists(f)) {
  x <- read.csv(f, check.names = FALSE)
  x <- x[
    x$cell_class == "Fibroblast" &
      grepl("Tumor", x$comparison, ignore.case = TRUE),
    ,
    drop = FALSE
  ]

  if (nrow(x)) {
    passed <- (
      x$mean_difference[1] > 0 &&
        x$p_value[1] < 0.10
    )
    n_pairs <- get_first(x, "n_pairs", NA_real_)

    set_evidence(
      "A",
      if (passed) "Supported" else "Not supported",
      passed,
      sprintf(
        "Fibroblast mean paired difference %s (95%% CI %s to %s), P=%s%s",
        fmt_n(x$mean_difference),
        fmt_n(x$ci_low),
        fmt_n(x$ci_high),
        fmt_p(x$p_value),
        if (is.finite(n_pairs)) paste0("; ", n_pairs, " paired patients") else ""
      ),
      f,
      if (passed) 4 else 0,
      if (passed) {
        "The locked score was increased in fibroblast pseudobulk profiles in an independent paired single-cell cohort."
      } else {
        "The locked score was not replicated in fibroblast pseudobulk profiles in the independent paired cohort."
      },
      "This is transcriptomic replication of a score difference, not direct validation of ligand secretion, causality or clinical utility."
    )
  }
}

# B: isolated fibroblasts.
f <- file.path(
  RESULT_ROOT,
  "03_GSE92945",
  "ECM_LIGAND_21_group_contrasts.csv"
)

if (file.exists(f)) {
  x <- read.csv(f, check.names = FALSE)
  x <- x[x$contrast == "Cancer vs Normal", , drop = FALSE]

  if (nrow(x)) {
    passed <- (
      x$mean_difference[1] > 0 &&
        x$p_value[1] < 0.05
    )

    set_evidence(
      "B",
      if (passed) "Supported" else "Not supported",
      passed,
      sprintf(
        "Cancer-minus-normal difference %s (95%% CI %s to %s), P=%s",
        fmt_n(x$mean_difference),
        fmt_n(x$ci_low),
        fmt_n(x$ci_high),
        fmt_p(x$p_value)
      ),
      f,
      if (passed) 3 else 0,
      if (passed) {
        "The isolated-fibroblast comparison supported a fibroblast-intrinsic association."
      } else {
        "The isolated-fibroblast comparison was directionally evaluated but statistically non-confirmatory."
      },
      "The dataset contains 4 cancer-associated, 3 colitis-associated and 3 normal independent fibroblast samples and has limited precision."
    )
  }
}

# C: exploratory specimen-level spatial association.
f <- file.path(
  RESULT_ROOT,
  "04_GSE280315_spatial",
  "spatial_combined_exploratory_summary.csv"
)

if (file.exists(f)) {
  x <- read.csv(f, check.names = FALSE)
  passed <- (
    x$positive_sections[1] >= 2 &&
      x$fisher_combined_permutation_p[1] < 0.05
  )

  set_evidence(
    "C",
    if (passed) "Supported with qualification" else "Not supported",
    passed,
    sprintf(
      "%d/%d CRC specimens had positive section-level rho; median rho %s; exploratory combined permutation P=%s",
      x$positive_sections[1],
      x$n_sections[1],
      fmt_n(x$median_rho),
      fmt_p(x$fisher_combined_permutation_p)
    ),
    f,
    if (passed) 3 else 0,
    if (passed) {
      "The scores showed a concordant exploratory spatial association across three CRC specimens."
    } else {
      "The exploratory spatial evidence was not concordant across the analysed CRC specimens."
    },
    "The reporting units are three CRC specimens (P1, P2 and P5); spatial bins are not independent patients, and association does not establish signalling direction or causality."
  )
}

# D: patient-level directional response to public TGF-beta perturbation data.
f <- file.path(
  RESULT_ROOT,
  "05_GSE160686_TGFb",
  "ECM_LIGAND_21_patient_effect_summary.csv"
)

if (file.exists(f)) {
  x <- read.csv(f, check.names = FALSE)
  directional <- (
    x$patients_with_decrease[1] >= 2 &&
      x$mean_difference[1] < 0
  )

  set_evidence(
    "D",
    if (directional) "Directionally supportive" else "Not supported",
    directional,
    sprintf(
      "Mean blockade-minus-control difference %s (95%% CI %s to %s), P=%s; decreased in %d/%d patients",
      fmt_n(x$mean_difference),
      fmt_n(x$ci_low),
      fmt_n(x$ci_high),
      fmt_p(x$p_value),
      x$patients_with_decrease[1],
      x$n_pairs[1]
    ),
    f,
    if (directional) 1 else 0,
    if (directional) {
      "Public perturbation data provided modest directional concordance with TGF-beta sensitivity."
    } else {
      "Public perturbation data did not provide a consistent directional response."
    },
    "Only three patients were available; the confidence interval crossed zero, and the all-cell pseudobulk estimate is not fibroblast-specific. The present study did not perform the experiment."
  )
}

# E: protein detectability is descriptive and is not scored as directional support.
f <- file.path(
  RESULT_ROOT,
  "08_ECM_proteomics",
  "ECM_LIGAND_21_individual_protein_statistics.csv"
)

if (file.exists(f)) {
  x <- read.csv(f, check.names = FALSE)
  n_detected <- length(unique(x$gene[!is.na(x$gene) & nzchar(x$gene)]))
  n_tumour <- sum(x$logFC > 0, na.rm = TRUE)
  n_normal <- sum(x$logFC < 0, na.rm = TRUE)
  coverage_met <- n_detected >= 10

  set_evidence(
    "E",
    if (coverage_met) "Descriptive only" else "Not supported",
    coverage_met,
    sprintf(
      "%d/21 locked proteins detected; %d tumour-enriched and %d normal-enriched",
      n_detected,
      n_tumour,
      n_normal
    ),
    f,
    0,
    if (coverage_met) {
      "Most locked proteins were detectable in independent decellularized-ECM proteomics, but the predominant direction was normal-enriched rather than tumour-enriched."
    } else {
      "Protein-level coverage was insufficient for a descriptive comparison."
    },
    "Detectability is not directional validation. Published proteomic measurements were reanalysed; no proteomic experiment was performed in the present study."
  )
}

# F: prespecified CAF-adjusted models. Report every split rather than selecting
# the smallest P value.
f <- file.path(
  RESULT_ROOT,
  "09_bulk_CAF_random",
  "CAF_adjusted_incremental_Cox_results.csv"
)

if (file.exists(f)) {
  x <- read.csv(f, check.names = FALSE)
  x <- x[x$estimable %in% TRUE, , drop = FALSE]

  if (nrow(x)) {
    split_order <- match(x$split, c("All", "Discovery", "Validation"))
    x <- x[order(split_order, na.last = TRUE), , drop = FALSE]

    observed_parts <- vapply(
      seq_len(nrow(x)),
      function(i) {
        sprintf(
          "%s HR/SD %s (95%% CI %s to %s), incremental LRT P=%s",
          x$split[i],
          fmt_n(x$HR_per_SD[i]),
          fmt_n(x$CI_low[i]),
          fmt_n(x$CI_high[i]),
          fmt_p(x$incremental_LRT_p[i])
        )
      },
      character(1)
    )

    discovery <- x[x$split == "Discovery", , drop = FALSE]
    validation <- x[x$split == "Validation", , drop = FALSE]

    discovery_supported <- (
      nrow(discovery) == 1 &&
        discovery$HR_per_SD[1] > 1 &&
        discovery$incremental_LRT_p[1] < 0.05
    )
    validation_direction <- (
      nrow(validation) == 1 &&
        validation$HR_per_SD[1] > 1
    )
    supported <- discovery_supported && validation_direction

    set_evidence(
      "F",
      if (supported) "Supported with qualification" else "Not supported",
      supported,
      paste(observed_parts, collapse = "; "),
      f,
      if (supported) 2 else 0,
      if (supported) {
        "The score added prognostic information beyond clinical covariates and the non-overlapping CAF proxy in the prespecified discovery model, with concordant validation direction."
      } else {
        "The score did not add prognostic information beyond clinical covariates and the non-overlapping CAF proxy in the prespecified analyses."
      },
      "Complete-case analysis reduces the available sample size. Incremental association does not by itself establish clinical utility."
    )
  }
}

# G: matched random-gene-set specificity benchmark.
f <- file.path(
  RESULT_ROOT,
  "09_bulk_CAF_random",
  "matched_random_gene_set_benchmark_summary.csv"
)

if (file.exists(f)) {
  x <- read.csv(f, check.names = FALSE)
  passed <- x$two_sided_empirical_p[1] < 0.05

  set_evidence(
    "G",
    if (passed) "Supported" else "Not supported",
    passed,
    sprintf(
      "Observed Wald z %s; two-sided add-one empirical P=%s versus %d estimable matched sets",
      fmt_n(x$observed_score_wald_z),
      fmt_p(x$two_sided_empirical_p),
      x$random_sets_estimable[1]
    ),
    f,
    if (passed) 3 else 0,
    if (passed) {
      "The observed prognostic statistic exceeded that expected from matched random 21-gene sets."
    } else {
      "The observed prognostic statistic did not outperform expression- and variability-matched random 21-gene sets."
    },
    "The benchmark evaluates gene-set specificity for the selected clinical statistic; it does not test every possible biological property of the score."
  )
}

set_contextual_status <- function(id, status_name, claim, limitations) {
  status_file <- file.path(PIPELINE_ROOT, "logs", status_name)

  if (!file.exists(status_file)) return(invisible(NULL))

  x <- read.csv(status_file, check.names = FALSE)
  passed <- (
    "status" %in% names(x) &&
      nrow(x) >= 1 &&
      all(toupper(trimws(as.character(x$status))) == "PASS")
  )
  detail <- if ("detail" %in% names(x)) {
    paste(as.character(x$detail), collapse = " | ")
  } else {
    paste("Technical status:", paste(unique(x$status), collapse = "; "))
  }

  set_evidence(
    id,
    if (passed) "Contextual only" else "Not evaluated",
    FALSE,
    detail,
    status_file,
    0,
    if (passed) claim else "The contextual analysis was not technically complete.",
    limitations
  )
}

# H and I are retained in the ledger so completed optional analyses are not
# silently omitted. They do not contribute to the internal journal-fit score.
set_contextual_status(
  "H",
  "06_GSE162561_status.csv",
  "The public stromal conditioned-medium dataset was retained as secondary experimental context and not as confirmatory patient-level evidence.",
  "Only two primary CRC models were available; deposited replicate libraries do not create additional independent patients."
)

set_contextual_status(
  "I",
  "07_GSE155343_status.csv",
  "The coculture dataset was retained as descriptive context and not as inferential validation.",
  "The deposited expression profiles are not an independent patient cohort, and rLog values are unsuitable for count-based differential-expression inference."
)

allowed_classes <- c(
  "Supported",
  "Supported with qualification",
  "Directionally supportive",
  "Descriptive only",
  "Contextual only",
  "Not supported",
  "Not evaluated"
)

if (!all(ledger$result_class %in% allowed_classes)) {
  stop("Unexpected evidence classification detected.")
}

if (any(ledger$earned_internal_points > ledger$maximum_internal_points)) {
  stop("Earned internal points exceeded their prespecified maxima.")
}

# Internal planning rubric only. It must not appear in manuscript text, figure
# legends, cover letters or supplementary scientific results.
base_score <- 78
earned_points <- sum(ledger$earned_internal_points)
fit_score <- base_score + earned_points

external_replication_gate <- (
  ledger$result_class[ledger$evidence_id == "A"] == "Supported"
)
secondary_experimental_context <- any(
  ledger$evidence_id %in% c("D", "E", "H", "I") & ledger$evaluated
)
de_novo_experimental_validation <- FALSE

fit_band <- if (fit_score >= 86) {
  "Conservative conditional fit; editorial risk remains because no de novo experimental validation was performed"
} else if (fit_score >= 82) {
  "Promising conditional fit, with material evidence limitations"
} else {
  "Current evidence is insufficient for a strong journal-fit claim"
}

write_csv(
  ledger,
  file.path(out, "mechanistic_evidence_ledger.csv")
)

write_csv(
  data.frame(
    journal = "Computational Biology and Chemistry",
    baseline_internal_score = base_score,
    earned_internal_points = earned_points,
    conservative_internal_fit_score = fit_score,
    external_replication_gate = external_replication_gate,
    secondary_experimental_context = secondary_experimental_context,
    de_novo_experimental_validation = de_novo_experimental_validation,
    fit_band = fit_band,
    publication_use = "Internal planning only; exclude this score from the manuscript, figures, supplementary material and cover letter.",
    scoring_note = "Subjective manuscript-journal fit rubric, not an acceptance probability or scientific outcome.",
    stringsAsFactors = FALSE
  ),
  file.path(out, "FIG_internal_fit_assessment_NOT_FOR_SUBMISSION.csv")
)

plot_ledger <- ledger[ledger$evaluated, , drop = FALSE]
plot_ledger$label <- paste(plot_ledger$evidence_id, plot_ledger$analysis)
plot_ledger$label <- factor(
  plot_ledger$label,
  levels = rev(plot_ledger$label)
)
plot_ledger$result_class <- factor(
  plot_ledger$result_class,
  levels = c(
    "Not supported",
    "Contextual only",
    "Descriptive only",
    "Directionally supportive",
    "Supported with qualification",
    "Supported"
  )
)

class_colours <- c(
  "Not supported" = "#666666",
  "Contextual only" = "#CC79A7",
  "Descriptive only" = "#56B4E9",
  "Directionally supportive" = "#E69F00",
  "Supported with qualification" = "#0072B2",
  "Supported" = "#004488"
)

class_shapes <- c(
  "Not supported" = 4,
  "Contextual only" = 1,
  "Descriptive only" = 0,
  "Directionally supportive" = 2,
  "Supported with qualification" = 17,
  "Supported" = 16
)

p <- ggplot2::ggplot(
  plot_ledger,
  ggplot2::aes(
    x = result_class,
    y = label,
    colour = result_class,
    shape = result_class
  )
) +
  ggplot2::geom_point(size = 4.2, stroke = 1.1) +
  ggplot2::scale_colour_manual(
    values = class_colours,
    drop = FALSE
  ) +
  ggplot2::scale_shape_manual(
    values = class_shapes,
    drop = FALSE
  ) +
  ggplot2::labs(
    x = "Evidence classification",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(
      angle = 25,
      hjust = 1
    ),
    panel.grid.major.y = ggplot2::element_line(
      colour = "grey90",
      linewidth = 0.3
    )
  )

save_plot(
  p,
  file.path(
    out,
    "Supplementary_Figure_external_evidence_claim_control"
  ),
  10.5,
  5.8
)

write_status(
  "10_integrate_evidence",
  "PASS",
  paste0(
    "Graded evidence ledger generated; internal planning score ",
    fit_score,
    "/100; de novo experimental validation: no; journal-fit score excluded from manuscript figure."
  )
)
