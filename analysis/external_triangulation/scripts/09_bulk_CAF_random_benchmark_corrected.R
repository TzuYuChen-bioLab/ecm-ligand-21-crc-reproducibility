PIPELINE_ROOT <- normalizePath(
  Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))

require_packages(c("data.table", "survival", "ggplot2"))

out <- result_dir("09_bulk_CAF_random")

# This final benchmark is deliberately more precise than the 1,000-set pilot.
# It must be chosen before outcome inspection and must not be increased after
# seeing the empirical P value.
N_RANDOM <- suppressWarnings(as.integer(
  Sys.getenv("CRC_RANDOM_GENE_SETS", unset = "10000")
))
if (!is.finite(N_RANDOM) || N_RANDOM < 1000L) {
  stop("CRC_RANDOM_GENE_SETS must be an integer >= 1000; use 10000 for the final run.")
}

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(ANALYSIS_OPTIONS$random_seed)

signature <- read_gene_set(SIGNATURE_FILE)
if (length(signature) != 21L) {
  stop("The locked ECM signature must contain exactly 21 unique genes; found ", length(signature), ".")
}

proxy_file <- file.path(
  RESULT_ROOT,
  "02_GSE144735",
  "CAF_PROXY_30_nonoverlapping.csv"
)

# Reuse the exact corrected GSE39582 objects created by the parent pipeline.
# Do not reselect probes or renormalize the GEO matrix here.
input_paths <- c(
  expression_rds = file.path(
    PROJECT_ROOT,
    "01_processed_data", "bulk_corrected", "GSE39582",
    "GSE39582_primary_tumour_expr_gene.rds"
  ),
  full_frozen_scores = file.path(
    PROJECT_ROOT,
    "08_tables", "bulk", "corrected_scores",
    "GSE39582_ECM_LIGAND_21_scores.csv"
  ),
  frozen_gene_parameters = file.path(
    PROJECT_ROOT,
    "08_tables", "bulk", "corrected_scores",
    "GSE39582_ECM_LIGAND_21_standardisation_parameters.csv"
  ),
  survival_ready = file.path(
    PROJECT_ROOT,
    "08_tables", "bulk", "corrected_survival", "GSE39582",
    "GSE39582_tumour_only_RFS_analysis_ready.csv"
  ),
  caf_proxy = proxy_file
)

input_audit <- data.frame(
  input = names(input_paths),
  path = unname(input_paths),
  exists = file.exists(input_paths),
  size_bytes = ifelse(file.exists(input_paths), file.info(input_paths)$size, NA_real_),
  stringsAsFactors = FALSE
)
write_csv(input_audit, file.path(out, "local_input_audit.csv"))
if (!all(input_audit$exists)) {
  stop(
    "Step 09 canonical input files are incomplete. Inspect: ",
    file.path(out, "local_input_audit.csv")
  )
}

clean_sample <- function(x) toupper(trimws(as.character(x)))

assert_unique <- function(x, label) {
  if (anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) {
    stop(label, " contains missing, blank, or duplicate identifiers.")
  }
  invisible(TRUE)
}

assert_named_counts <- function(observed, expected, label) {
  observed <- observed[names(expected)]
  observed[is.na(observed)] <- 0L
  if (!identical(unname(as.integer(observed)), unname(as.integer(expected)))) {
    stop(
      label, " counts differ from the frozen corrected cohort. Expected ",
      paste(names(expected), expected, sep = "=", collapse = ", "),
      "; observed ",
      paste(names(observed), observed, sep = "=", collapse = ", "), "."
    )
  }
  invisible(TRUE)
}

first_existing_name <- function(x, candidates, label) {
  hit <- candidates[candidates %in% names(x)]
  if (!length(hit)) {
    stop("Required ", label, " column was not found. Tried: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

expression <- readRDS(input_paths[["expression_rds"]])
if (!is.matrix(expression) || !nrow(expression) || !ncol(expression)) {
  stop("The corrected GSE39582 gene-expression RDS is not a non-empty matrix.")
}
storage.mode(expression) <- "numeric"
rownames(expression) <- clean_symbol(rownames(expression))
colnames(expression) <- clean_sample(colnames(expression))
assert_unique(rownames(expression), "GSE39582 gene symbols")
assert_unique(colnames(expression), "GSE39582 expression samples")

full_scores <- read.csv(
  input_paths[["full_frozen_scores"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!all(c("sample", "ECM_LIGAND_21") %in% names(full_scores))) {
  stop("The full frozen-score table lacks sample or ECM_LIGAND_21.")
}
split_col <- first_existing_name(
  full_scores,
  c("parsed_dataset", "dataset:ch1"),
  "GSE39582 discovery/validation split"
)
full_scores$sample <- clean_sample(full_scores$sample)
full_scores$cohort_split_fixed <- tools::toTitleCase(
  tolower(trimws(as.character(full_scores[[split_col]])))
)
assert_unique(full_scores$sample, "Full frozen-score samples")
assert_named_counts(
  table(full_scores$cohort_split_fixed),
  c(Discovery = 443L, Validation = 123L),
  "Full tumour-only GSE39582"
)
if (nrow(full_scores) != 566L) {
  stop("The full tumour-only GSE39582 score table must contain 566 samples.")
}

full_expression_index <- match(full_scores$sample, colnames(expression))
if (anyNA(full_expression_index) || length(unique(full_expression_index)) != 566L) {
  write_csv(
    data.frame(
      sample = full_scores$sample,
      expression_column = full_expression_index,
      matched = !is.na(full_expression_index)
    ),
    file.path(out, "full_sample_alignment_NEEDS_REVIEW.csv")
  )
  stop("All 566 corrected GSE39582 tumour samples must match expression exactly once.")
}
expression <- expression[, full_expression_index, drop = FALSE]
colnames(expression) <- full_scores$sample

surv_data <- read.csv(
  input_paths[["survival_ready"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_survival <- c(
  "sample", "cohort_split", "time", "event", "score_raw", "score_z",
  "age", "sex", "stage", "mmr"
)
if (!all(required_survival %in% names(surv_data))) {
  stop(
    "The GSE39582 survival-ready table lacks: ",
    paste(setdiff(required_survival, names(surv_data)), collapse = ", ")
  )
}
surv_data$sample <- clean_sample(surv_data$sample)
surv_data$cohort_split <- tools::toTitleCase(
  tolower(trimws(as.character(surv_data$cohort_split)))
)
assert_unique(surv_data$sample, "GSE39582 survival-ready samples")
assert_named_counts(
  table(surv_data$cohort_split),
  c(Discovery = 410L, Validation = 109L),
  "RFS analysis-ready GSE39582"
)
assert_named_counts(
  tapply(surv_data$event, surv_data$cohort_split, sum),
  c(Discovery = 115L, Validation = 30L),
  "RFS event"
)
if (
  nrow(surv_data) != 519L ||
  any(!is.finite(surv_data$time)) ||
  any(surv_data$time <= 0) ||
  any(!surv_data$event %in% c(0, 1)) ||
  any(!is.finite(surv_data$score_raw)) ||
  any(!is.finite(surv_data$score_z))
) {
  stop("The frozen RFS analysis-ready population failed its integrity checks.")
}

survival_full_index <- match(surv_data$sample, full_scores$sample)
if (anyNA(survival_full_index) || length(unique(survival_full_index)) != 519L) {
  stop("All 519 RFS analysis samples must map uniquely into the 566 tumour samples.")
}

# Reproduce the frozen raw ECM score exactly from the saved reference parameters.
frozen_parameters <- read.csv(
  input_paths[["frozen_gene_parameters"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!all(c("gene", "centre", "scale") %in% names(frozen_parameters))) {
  stop("Frozen ECM score parameters must contain gene, centre, and scale.")
}
frozen_parameters$gene <- clean_symbol(frozen_parameters$gene)
if (
  nrow(frozen_parameters) != 21L ||
  anyDuplicated(frozen_parameters$gene) ||
  !setequal(frozen_parameters$gene, signature) ||
  any(!is.finite(frozen_parameters$centre)) ||
  any(!is.finite(frozen_parameters$scale)) ||
  any(frozen_parameters$scale <= 0) ||
  !all(frozen_parameters$gene %in% rownames(expression))
) {
  stop("The saved 21-gene frozen standardisation parameters failed validation.")
}
frozen_parameters <- frozen_parameters[match(signature, frozen_parameters$gene), ]
frozen_z <- sweep(
  expression[signature, , drop = FALSE],
  1,
  frozen_parameters$centre,
  "-"
)
frozen_z <- sweep(frozen_z, 1, frozen_parameters$scale, "/")
ecm_recomputed_raw <- colMeans(frozen_z)
frozen_saved_raw <- as.numeric(full_scores$ECM_LIGAND_21)

reproduction_audit <- data.frame(
  n = length(frozen_saved_raw),
  pearson_frozen_vs_recomputed = suppressWarnings(stats::cor(
    frozen_saved_raw, ecm_recomputed_raw, method = "pearson"
  )),
  spearman_frozen_vs_recomputed = suppressWarnings(stats::cor(
    frozen_saved_raw, ecm_recomputed_raw, method = "spearman"
  )),
  maximum_absolute_difference = max(abs(frozen_saved_raw - ecm_recomputed_raw)),
  tolerance = 1e-6,
  pass = max(abs(frozen_saved_raw - ecm_recomputed_raw)) <= 1e-6
)
write_csv(
  reproduction_audit,
  file.path(out, "frozen_score_reproduction_audit.csv")
)
if (!isTRUE(reproduction_audit$pass[1])) {
  stop("Frozen ECM score could not be reproduced to tolerance 1e-6.")
}

# The CAF proxy is derived independently in Step 02 and is fixed before this
# survival analysis. All centring/scaling uses all 443 discovery tumours, not
# the smaller outcome-complete subset.
proxy_table <- read.csv(proxy_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"gene" %in% names(proxy_table)) {
  stop("CAF proxy file must contain a gene column.")
}
proxy <- unique(stats::na.omit(clean_symbol(proxy_table$gene)))
if (length(proxy) != 30L) {
  stop("The locked non-overlapping CAF proxy must contain exactly 30 unique genes; found ", length(proxy), ".")
}
if (length(intersect(proxy, signature))) {
  stop("CAF proxy overlaps the locked ECM_LIGAND_21 signature.")
}

discovery_samples <- full_scores$sample[
  full_scores$cohort_split_fixed == "Discovery"
]
if (length(discovery_samples) != 443L) {
  stop("Discovery reference must contain exactly 443 tumour samples.")
}

score_with_discovery_reference <- function(expr, genes, reference_samples, minimum_coverage) {
  genes <- unique(clean_symbol(genes))
  present <- intersect(genes, rownames(expr))
  initial_coverage <- length(present) / length(genes)
  if (initial_coverage < minimum_coverage) {
    stop(sprintf(
      "Gene-set coverage %.1f%% is below the locked %.1f%% threshold.",
      100 * initial_coverage,
      100 * minimum_coverage
    ))
  }
  ref <- expr[present, reference_samples, drop = FALSE]
  centres <- rowMeans(ref, na.rm = TRUE)
  scales <- apply(ref, 1, stats::sd, na.rm = TRUE)
  valid <- is.finite(centres) & is.finite(scales) & scales > 0
  present <- present[valid]
  coverage <- length(present) / length(genes)
  if (coverage < minimum_coverage) {
    stop("Too many detected genes had zero or invalid Discovery-reference variance.")
  }
  centres <- centres[present]
  scales <- scales[present]
  z <- sweep(expr[present, , drop = FALSE], 1, centres, "-")
  z <- sweep(z, 1, scales, "/")
  raw_score <- colMeans(z, na.rm = TRUE)
  reference_score <- raw_score[match(reference_samples, colnames(expr))]
  composite_centre <- mean(reference_score)
  composite_scale <- stats::sd(reference_score)
  if (!is.finite(composite_scale) || composite_scale <= 0) {
    stop("Discovery-reference composite score has invalid variance.")
  }
  list(
    raw = raw_score,
    z = (raw_score - composite_centre) / composite_scale,
    present = present,
    missing = setdiff(genes, present),
    coverage = coverage,
    parameters = data.frame(
      gene = genes,
      detected_and_variable = genes %in% present,
      centre = unname(centres[match(genes, names(centres))]),
      scale = unname(scales[match(genes, names(scales))]),
      reference_n = length(reference_samples),
      stringsAsFactors = FALSE
    ),
    composite_centre = composite_centre,
    composite_scale = composite_scale
  )
}

caf_object <- score_with_discovery_reference(
  expression,
  proxy,
  discovery_samples,
  minimum_coverage = ANALYSIS_OPTIONS$minimum_signature_coverage
)
write_csv(
  transform(
    caf_object$parameters,
    gene_set = "CAF_PROXY_30_nonoverlapping",
    coverage = caf_object$coverage
  ),
  file.path(out, "CAF_PROXY_coverage.csv")
)
write_csv(
  data.frame(
    gene_set = "CAF_PROXY_30_nonoverlapping",
    reference = "all 443 GSE39582 discovery tumours",
    reference_n = 443L,
    composite_centre = caf_object$composite_centre,
    composite_scale = caf_object$composite_scale,
    coverage = caf_object$coverage
  ),
  file.path(out, "CAF_PROXY_standardisation_summary.csv")
)

dat <- surv_data
dat$frozen_score_recomputed_raw <- ecm_recomputed_raw[survival_full_index]
dat$caf_proxy_raw <- caf_object$raw[survival_full_index]
dat$caf_proxy_z <- caf_object$z[survival_full_index]
if (max(abs(dat$score_raw - dat$frozen_score_recomputed_raw)) > 1e-6) {
  stop("The survival-ready score_raw values do not match the reproduced frozen scores.")
}

dat$age <- suppressWarnings(as.numeric(dat$age))
dat$sex <- factor(dat$sex)
dat$stage <- factor(dat$stage, levels = c("I", "II", "III", "IV"))
dat$mmr <- factor(dat$mmr)
dat$cohort_split <- factor(
  dat$cohort_split,
  levels = c("Discovery", "Validation")
)

write_csv(dat, file.path(out, "GSE39582_with_nonoverlapping_CAF_proxy.csv"))
write_csv(
  data.frame(
    n = nrow(dat),
    pearson = stats::cor(dat$score_z, dat$caf_proxy_z, use = "complete.obs"),
    spearman = stats::cor(
      dat$score_z, dat$caf_proxy_z,
      method = "spearman", use = "complete.obs"
    ),
    signature_proxy_gene_overlap = length(intersect(signature, proxy))
  ),
  file.path(out, "ECM_score_CAF_proxy_correlation_audit.csv")
)

strata <- survival::strata

capture_cox <- function(formula, data, full_output = TRUE) {
  warnings <- character(0)
  fit <- withCallingHandlers(
    tryCatch(
      survival::coxph(
        formula,
        data = data,
        ties = "efron",
        x = full_output,
        y = full_output,
        model = full_output
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = unique(warnings))
}

coefficient_is_stable <- function(fit, term) {
  if (inherits(fit, "error") || is.null(fit)) return(FALSE)
  co <- stats::coef(fit)
  if (!term %in% names(co)) return(FALSE)
  vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(vc) || !term %in% rownames(vc)) return(FALSE)
  se <- sqrt(diag(vc))[term]
  is.finite(co[term]) && is.finite(se) && se > 0 && abs(co[term]) < 10 && se < 10
}

fit_incremental_pair <- function(d, split, analysis_type, base_formula, full_formula) {
  variables <- unique(c(all.vars(base_formula), all.vars(full_formula)))
  keep <- stats::complete.cases(d[, variables, drop = FALSE]) &
    is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1)
  model_data <- droplevels(d[keep, , drop = FALSE])
  n <- nrow(model_data)
  events <- sum(model_data$event)
  if (n < 40L || events < 15L) {
    return(data.frame(
      split = split, analysis_type = analysis_type, estimable = FALSE,
      n = n, events = events, reason = "insufficient_n_or_events",
      stringsAsFactors = FALSE
    ))
  }

  base_capture <- capture_cox(base_formula, model_data)
  full_capture <- capture_cox(full_formula, model_data)
  bad_warning <- grepl(
    "did not converge|infinite|singular",
    paste(c(base_capture$warnings, full_capture$warnings), collapse = " | "),
    ignore.case = TRUE
  )
  estimable <- !inherits(base_capture$fit, "error") &&
    !inherits(full_capture$fit, "error") &&
    !bad_warning &&
    coefficient_is_stable(full_capture$fit, "score_z")

  if (!estimable) {
    reason <- paste(
      c(
        if (inherits(base_capture$fit, "error")) conditionMessage(base_capture$fit),
        if (inherits(full_capture$fit, "error")) conditionMessage(full_capture$fit),
        base_capture$warnings,
        full_capture$warnings,
        if (!coefficient_is_stable(full_capture$fit, "score_z")) "unstable_score_coefficient"
      ),
      collapse = " | "
    )
    return(data.frame(
      split = split, analysis_type = analysis_type, estimable = FALSE,
      n = n, events = events, reason = reason,
      stringsAsFactors = FALSE
    ))
  }

  full <- full_capture$fit
  base <- base_capture$fit
  co <- summary(full)$coefficients
  ci <- summary(full)$conf.int
  row <- match("score_z", rownames(co))
  lrt <- stats::anova(base, full, test = "LRT")
  p_column <- grep("^Pr", colnames(lrt), value = TRUE)[1]
  ph <- tryCatch(survival::cox.zph(full), error = function(e) NULL)
  score_ph <- if (!is.null(ph) && "score_z" %in% rownames(ph$table)) {
    ph$table["score_z", "p"]
  } else {
    NA_real_
  }
  global_ph <- if (!is.null(ph) && "GLOBAL" %in% rownames(ph$table)) {
    ph$table["GLOBAL", "p"]
  } else {
    NA_real_
  }

  data.frame(
    split = split,
    analysis_type = analysis_type,
    estimable = TRUE,
    n = n,
    events = events,
    covariates = paste(
      setdiff(all.vars(full_formula), c("time", "event", "score_z")),
      collapse = ";"
    ),
    HR_per_SD = ci[row, "exp(coef)"],
    CI_low = ci[row, "lower .95"],
    CI_high = ci[row, "upper .95"],
    wald_p = co[row, "Pr(>|z|)"],
    wald_z = co[row, "z"],
    incremental_LRT_p = lrt[2, p_column],
    c_index_base = unname(summary(base)$concordance[1]),
    c_index_full = unname(summary(full)$concordance[1]),
    delta_c_index = unname(
      summary(full)$concordance[1] - summary(base)$concordance[1]
    ),
    score_PH_p = score_ph,
    global_PH_p = global_ph,
    base_formula = paste(deparse(base_formula), collapse = " "),
    full_formula = paste(deparse(full_formula), collapse = " "),
    warnings = paste(unique(c(base_capture$warnings, full_capture$warnings)), collapse = " | "),
    reason = NA_character_,
    stringsAsFactors = FALSE
  )
}

# One prespecified estimable model per reporting population. The validation
# cohort uses the same stage/sex-stratified sensitivity structure as the parent
# corrected pipeline because sparse/zero-event factor levels make a fully
# coefficient-adjusted validation model unstable.
model_registry <- data.frame(
  split = c("All", "Discovery", "Validation"),
  analysis_type = c(
    "pooled_cohort_sex_stage_stratified",
    "primary_prespecified_plus_CAF",
    "validation_sex_stage_stratified_plus_CAF"
  ),
  base_formula = c(
    "Surv(time,event) ~ age + mmr + caf_proxy_z + strata(cohort_split,sex,stage)",
    "Surv(time,event) ~ age + sex + stage + mmr + caf_proxy_z",
    "Surv(time,event) ~ age + mmr + caf_proxy_z + strata(sex,stage)"
  ),
  full_formula = c(
    "base + score_z",
    "base + score_z",
    "base + score_z"
  ),
  role = c("supportive pooled", "primary", "prespecified estimability sensitivity"),
  stringsAsFactors = FALSE
)
write_csv(model_registry, file.path(out, "CAF_adjusted_model_registry.csv"))

all_base <- survival::Surv(time, event) ~
  age + mmr + caf_proxy_z + strata(cohort_split, sex, stage)
all_full <- survival::Surv(time, event) ~
  age + mmr + caf_proxy_z + score_z + strata(cohort_split, sex, stage)
discovery_base <- survival::Surv(time, event) ~
  age + sex + stage + mmr + caf_proxy_z
discovery_full <- survival::Surv(time, event) ~
  age + sex + stage + mmr + caf_proxy_z + score_z
validation_base <- survival::Surv(time, event) ~
  age + mmr + caf_proxy_z + strata(sex, stage)
validation_full <- survival::Surv(time, event) ~
  age + mmr + caf_proxy_z + score_z + strata(sex, stage)

models <- list(
  fit_incremental_pair(
    dat, "All", "pooled_cohort_sex_stage_stratified",
    all_base, all_full
  ),
  fit_incremental_pair(
    dat[dat$cohort_split == "Discovery", , drop = FALSE],
    "Discovery", "primary_prespecified_plus_CAF",
    discovery_base, discovery_full
  ),
  fit_incremental_pair(
    dat[dat$cohort_split == "Validation", , drop = FALSE],
    "Validation", "validation_sex_stage_stratified_plus_CAF",
    validation_base, validation_full
  )
)
models <- data.table::rbindlist(models, fill = TRUE)
write_csv(models, file.path(out, "CAF_adjusted_incremental_Cox_results.csv"))
if (nrow(models) != 3L || any(!models$estimable)) {
  stop("At least one required CAF-adjusted Cox model failed its estimability gate.")
}

# Expression/variance-matched random 21-gene benchmark. Gene matching and
# standardisation use all 443 Discovery samples and never use outcomes.
benchmark_dat <- dat[dat$cohort_split == "Discovery", , drop = FALSE]
benchmark_variables <- c(
  "time", "event", "score_z", "caf_proxy_z", "age", "sex", "stage", "mmr", "sample"
)
benchmark_dat <- droplevels(benchmark_dat[
  stats::complete.cases(benchmark_dat[, benchmark_variables, drop = FALSE]) &
    benchmark_dat$time > 0 & benchmark_dat$event %in% c(0, 1),
  , drop = FALSE
])
if (nrow(benchmark_dat) < 300L || sum(benchmark_dat$event) < 80L) {
  stop("Discovery primary benchmark population is unexpectedly small.")
}

reference_expression <- expression[, discovery_samples, drop = FALSE]
gene_mean <- rowMeans(reference_expression, na.rm = TRUE)
gene_sd <- apply(reference_expression, 1, stats::sd, na.rm = TRUE)
valid_gene <- is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0
eligible <- names(gene_mean)[valid_gene]
eligible <- setdiff(eligible, union(signature, proxy))
if (length(eligible) < 5000L) {
  stop("Fewer than 5,000 eligible non-signature genes are available for matching.")
}
if (!all(signature %in% names(gene_mean)[valid_gene])) {
  stop("All 21 locked signature genes must have finite Discovery-reference mean and SD.")
}

metric_genes <- union(eligible, signature)
mean_metric <- as.numeric(scale(gene_mean[metric_genes]))
names(mean_metric) <- metric_genes
log_sd_metric <- as.numeric(scale(log(gene_sd[metric_genes])))
names(log_sd_metric) <- metric_genes
if (any(!is.finite(mean_metric)) || any(!is.finite(log_sd_metric))) {
  stop("Expression-matching metrics contain non-finite values.")
}

MATCH_POOL_SIZE <- min(200L, length(eligible))
ordered_candidates <- vector("list", length(signature))
names(ordered_candidates) <- signature
candidate_distances <- vector("list", length(signature))
names(candidate_distances) <- signature
pool_audit <- vector("list", length(signature))

for (g in signature) {
  distance <- sqrt(
    (mean_metric[eligible] - mean_metric[g])^2 +
      (log_sd_metric[eligible] - log_sd_metric[g])^2
  )
  names(distance) <- eligible
  ord <- order(distance, names(distance))
  ordered_candidates[[g]] <- names(distance)[ord]
  candidate_distances[[g]] <- distance
  pool <- ordered_candidates[[g]][seq_len(MATCH_POOL_SIZE)]
  pool_audit[[g]] <- data.frame(
    signature_gene = g,
    signature_expression_mean = gene_mean[g],
    signature_expression_sd = gene_sd[g],
    candidate_pool_size = length(pool),
    minimum_scaled_distance = min(distance[pool]),
    median_scaled_distance = stats::median(distance[pool]),
    maximum_scaled_distance = max(distance[pool]),
    stringsAsFactors = FALSE
  )
}
pool_audit <- do.call(rbind, pool_audit)
write_csv(pool_audit, file.path(out, "random_gene_matching_pool_audit.csv"))

sample_matched_set <- function() {
  selected <- setNames(rep(NA_character_, length(signature)), signature)
  for (g in sample(signature, length(signature), replace = FALSE)) {
    used <- unname(selected[!is.na(selected)])
    pool <- setdiff(
      ordered_candidates[[g]][seq_len(MATCH_POOL_SIZE)],
      used
    )
    if (!length(pool)) {
      pool <- setdiff(ordered_candidates[[g]], used)
    }
    if (!length(pool)) {
      stop("Unable to construct a unique 21-gene matched random set.")
    }
    selected[g] <- sample(pool, 1L)
  }
  unname(selected[signature])
}

# Standardise every eligible gene on all 443 reference samples once.
reference_z <- zscore_rows(reference_expression)
benchmark_reference_index <- match(benchmark_dat$sample, colnames(reference_z))
if (anyNA(benchmark_reference_index)) {
  stop("Benchmark survival samples did not align to the Discovery expression reference.")
}

observed_capture <- capture_cox(discovery_full, benchmark_dat, full_output = FALSE)
if (
  inherits(observed_capture$fit, "error") ||
  !coefficient_is_stable(observed_capture$fit, "score_z")
) {
  stop("Observed Discovery CAF-adjusted score model is not estimable.")
}
observed_co <- summary(observed_capture$fit)$coefficients
observed_z <- observed_co["score_z", "z"]
observed_hr <- exp(stats::coef(observed_capture$fit)["score_z"])

random_formula <- survival::Surv(time, event) ~
  age + sex + stage + mmr + caf_proxy_z + random_score_z
random_results <- vector("list", N_RANDOM)

for (i in seq_len(N_RANDOM)) {
  genes <- sample_matched_set()
  if (length(genes) != 21L || anyDuplicated(genes)) {
    stop("Internal error: a random set was not exactly 21 unique genes.")
  }
  random_raw_reference <- colMeans(reference_z[genes, , drop = FALSE])
  random_scale <- stats::sd(random_raw_reference)
  if (!is.finite(random_scale) || random_scale <= 0) {
    random_results[[i]] <- data.frame(
      iteration = i,
      wald_z = NA_real_,
      HR_per_SD = NA_real_,
      mean_matching_distance = NA_real_,
      maximum_matching_distance = NA_real_,
      n_genes = 21L,
      unique_genes = 21L,
      genes = paste(genes, collapse = ";"),
      warning = "invalid_random_composite_variance",
      stringsAsFactors = FALSE
    )
    next
  }
  random_z_reference <- (
    random_raw_reference - mean(random_raw_reference)
  ) / random_scale
  tmp <- benchmark_dat
  tmp$random_score_z <- random_z_reference[benchmark_reference_index]
  fit_capture <- capture_cox(random_formula, tmp, full_output = FALSE)
  stable <- coefficient_is_stable(fit_capture$fit, "random_score_z")
  if (stable) {
    co <- stats::coef(fit_capture$fit)["random_score_z"]
    se <- sqrt(diag(stats::vcov(fit_capture$fit)))["random_score_z"]
    z <- co / se
    hr <- exp(co)
  } else {
    z <- NA_real_
    hr <- NA_real_
  }
  distances <- vapply(
    seq_along(signature),
    function(j) candidate_distances[[signature[j]]][genes[j]],
    numeric(1)
  )
  random_results[[i]] <- data.frame(
    iteration = i,
    wald_z = unname(z),
    HR_per_SD = unname(hr),
    mean_matching_distance = mean(distances),
    maximum_matching_distance = max(distances),
    n_genes = length(genes),
    unique_genes = length(unique(genes)),
    genes = paste(genes, collapse = ";"),
    warning = paste(fit_capture$warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

random_results <- data.table::rbindlist(random_results, fill = TRUE)
write_csv(
  random_results,
  file.path(out, "matched_random_gene_set_results.csv")
)

estimable <- is.finite(random_results$wald_z)
random_sets_estimable <- sum(estimable)
if (random_sets_estimable < ceiling(0.99 * N_RANDOM)) {
  stop("Fewer than 99% of the matched random gene-set Cox models were estimable.")
}
exceedances <- sum(
  abs(random_results$wald_z[estimable]) >= abs(observed_z)
)
empirical_p <- (exceedances + 1) / (random_sets_estimable + 1)
binomial_ci <- stats::binom.test(exceedances, random_sets_estimable)$conf.int

set_keys <- vapply(
  strsplit(random_results$genes, ";", fixed = TRUE),
  function(x) paste(sort(x), collapse = ";"),
  character(1)
)

benchmark <- data.frame(
  cohort_split = "Discovery",
  n = nrow(benchmark_dat),
  events = sum(benchmark_dat$event),
  observed_score_HR_per_SD = unname(observed_hr),
  observed_score_wald_z = unname(observed_z),
  random_sets_requested = N_RANDOM,
  random_sets_estimable = random_sets_estimable,
  exceedances_two_sided = exceedances,
  two_sided_empirical_p = empirical_p,
  binomial_tail_probability_CI_low = unname(binomial_ci[1]),
  binomial_tail_probability_CI_high = unname(binomial_ci[2]),
  minimum_attainable_empirical_p = 1 / (random_sets_estimable + 1),
  unique_random_sets = length(unique(set_keys)),
  matching = paste0(
    "without-replacement 21-gene sets; each signature gene matched from its ",
    MATCH_POOL_SIZE,
    " nearest genes on scaled Discovery-reference expression mean and log(SD); ",
    "ECM_LIGAND_21 and CAF_PROXY_30 genes excluded; outcomes unused"
  ),
  standardisation_reference = "all 443 GSE39582 discovery tumours",
  empirical_p_definition = "(1 + count(abs(random z) >= abs(observed z))) / (1 + estimable random sets)",
  rng_kind = paste(RNGkind(), collapse = ";"),
  seed = ANALYSIS_OPTIONS$random_seed,
  stringsAsFactors = FALSE
)
write_csv(
  benchmark,
  file.path(out, "matched_random_gene_set_benchmark_summary.csv")
)

p_random <- ggplot2::ggplot(
  random_results[estimable, ],
  ggplot2::aes(x = wald_z)
) +
  ggplot2::geom_histogram(
    bins = 60,
    fill = "#8FB9D5",
    colour = "white",
    linewidth = 0.2
  ) +
  ggplot2::geom_vline(
    xintercept = c(-abs(observed_z), abs(observed_z)),
    colour = "#C51B33",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    x = "CAF-adjusted Cox Wald z",
    y = "Matched random gene sets",
    title = "Expression/variance-matched 21-gene null benchmark",
    subtitle = paste0(
      "Discovery cohort; B = ", random_sets_estimable,
      "; observed |z| = ", sprintf("%.2f", abs(observed_z)),
      "; empirical P = ", format.pval(empirical_p, digits = 3)
    )
  )
save_plot(
  p_random,
  file.path(out, "GSE39582_matched_random_gene_set_benchmark"),
  width = 7.2,
  height = 5.2
)

final_gates <- data.frame(
  gate = c(
    "canonical_inputs_present",
    "full_tumour_population_443_123",
    "RFS_population_410_109",
    "RFS_events_115_30",
    "frozen_score_exactly_reproduced",
    "CAF_proxy_30_genes_nonoverlapping",
    "three_CAF_adjusted_models_estimable",
    "random_sets_at_least_99_percent_estimable",
    "every_random_set_21_unique_genes",
    "at_least_95_percent_random_sets_unique"
  ),
  pass = c(
    all(input_audit$exists),
    TRUE,
    TRUE,
    TRUE,
    isTRUE(reproduction_audit$pass[1]),
    length(proxy) == 30L && !length(intersect(proxy, signature)),
    nrow(models) == 3L && all(models$estimable),
    random_sets_estimable >= ceiling(0.99 * N_RANDOM),
    all(random_results$n_genes == 21L & random_results$unique_genes == 21L),
    length(unique(set_keys)) >= ceiling(0.95 * N_RANDOM)
  ),
  stringsAsFactors = FALSE
)
write_csv(final_gates, file.path(out, "Step09_final_gate_audit.csv"))
if (!all(final_gates$pass)) {
  stop("Step 09 completed computations but failed at least one final audit gate.")
}

write_status(
  "09_bulk_CAF_random",
  "PASS",
  paste0(
    "Frozen-score reproduction passed; three CAF-adjusted models estimable; ",
    random_sets_estimable, "/", N_RANDOM,
    " matched random sets estimable; two-sided add-one empirical P = ",
    signif(empirical_p, 5), "."
  )
)

