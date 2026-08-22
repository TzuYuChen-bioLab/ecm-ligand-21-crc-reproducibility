#!/usr/bin/env Rscript

# =============================================================================
# GSE39582 ECM21 module CAF adjustment and Validation robustness continuation
#
# This script deliberately DOES NOT redefine the 21-gene panel, the four
# modules, the Discovery/Validation split, or the primary enhancement gate.
# It patches and runs the supplied base script, adds the exact frozen
# CAF_PROXY_30, handles the historical SDPR/current CAVIN2 gene-name alias,
# and then adds prespecified-looking, clearly labelled sensitivity analyses.
#
# Put these three files in the same folder:
#   1) RUN_GSE39582_MODULE_CAF_AND_VALIDATION_SENSITIVITY.R  (this file)
#   2) RUN_GSE39582_MODULE_COX_AND_BOOTSTRAP_DELTA_CINDEX.R (base script)
#   3) CAF_PROXY_30_nonoverlapping.csv                      (frozen resource)
#
# Required expression input (already present on the author's computer):
#   PROJECT_ROOT
#   GSE39582_series_matrix.txt.gz
#
# Run in RStudio with Source, or:
#   Rscript RUN_GSE39582_MODULE_CAF_AND_VALIDATION_SENSITIVITY.R
# =============================================================================

options(stringsAsFactors = FALSE)

# ------------------------------- USER SETTINGS -------------------------------
USER_ROOT <- Sys.getenv("CBC_REPRO_ROOT", unset = "")
BASE_SCRIPT_NAME <- "RUN_GSE39582_MODULE_COX_AND_BOOTSTRAP_DELTA_CINDEX.R"
CAF_GENE_FILE_NAME <- "CAF_PROXY_30_nonoverlapping.csv"

# The base analysis uses this number for its ordinary patient bootstrap.
BASE_BOOTSTRAP_REPS <- 1000L

# This additional Validation-only bootstrap preserves event/covariate blocks.
# It is a sensitivity analysis and must not replace the ordinary bootstrap.
VALIDATION_STRATIFIED_BOOTSTRAP_REPS <- 2000L
VALIDATION_STRATIFIED_BOOTSTRAP_SEED <- 20260813L
MIN_BOOTSTRAP_SUCCESS_FRACTION <- 0.80

# ------------------------------- PATH HELPERS --------------------------------
script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) {
    return(dirname(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE)))
  }
  frames <- sys.frames()
  ofiles <- vapply(frames, function(x) if (!is.null(x$ofile)) as.character(x$ofile)[1] else "", character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles)) return(dirname(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = FALSE)))
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

SCRIPT_DIR <- script_dir()
PACKAGE_ROOT <- normalizePath(file.path(SCRIPT_DIR, '..', '..', '..'), winslash = '/', mustWork = TRUE)
if (!nzchar(USER_ROOT)) USER_ROOT <- file.path(PACKAGE_ROOT, 'results', 'generated')
dir.create(USER_ROOT, recursive = TRUE, showWarnings = FALSE)
BASE_SCRIPT <- file.path(SCRIPT_DIR, BASE_SCRIPT_NAME)
CAF_GENE_FILE <- file.path(PACKAGE_ROOT, "analysis", "external_triangulation", "resources", CAF_GENE_FILE_NAME)

if (!file.exists(BASE_SCRIPT)) stop("Base script not found beside this file: ", BASE_SCRIPT)
if (!file.exists(CAF_GENE_FILE)) stop("Frozen CAF gene file not found beside this file: ", CAF_GENE_FILE)

series_matrix <- file.path(PACKAGE_ROOT, "data", "available", "GSE39582", "GSE39582_series_matrix.txt.gz")
if (!file.exists(series_matrix)) stop("Official GSE39582 series matrix not found: ", series_matrix)

path_literal <- function(x) {
  x <- normalizePath(x, winslash = "/", mustWork = TRUE)
  paste0('"', gsub('"', '\\\"', x, fixed = TRUE), '"')
}

# ------------------------- PATCH THE SUPPLIED BASE RUN ------------------------
base_lines <- readLines(BASE_SCRIPT, warn = FALSE)

replace_setting <- function(lines, setting, value_expression) {
  # Only search the initial USER SETTINGS block. For example, the base script
  # later contains `USER_ROOT <- normalizePath(...)`, which is executable
  # setup code rather than the user-supplied path that should be patched.
  settings_end <- grep("^# ----------------------------- FIXED GENE MODULES", lines)
  if (!length(settings_end)) stop("Could not locate the end of the USER SETTINGS block in the base script.")
  search_index <- seq_len(settings_end[1] - 1L)
  idx <- search_index[grep(paste0("^", setting, "[[:space:]]*<-"), lines[search_index])]
  if (length(idx) != 1L) {
    stop("Expected exactly one USER SETTINGS assignment named ", setting, "; found ", length(idx))
  }
  lines[idx] <- paste0(setting, " <- ", value_expression)
  lines
}

base_lines <- replace_setting(base_lines, "USER_ROOT", paste0('"', gsub("\\\\", "/", USER_ROOT), '"'))
base_lines <- replace_setting(base_lines, "CAF_SCORE_FILE", '""')
base_lines <- replace_setting(base_lines, "CAF_GENE_FILE", path_literal(CAF_GENE_FILE))
base_lines <- replace_setting(base_lines, "REQUIRE_CAF_INPUT", "TRUE")
base_lines <- replace_setting(base_lines, "CINDEX_BOOTSTRAP_REPS", paste0(as.integer(BASE_BOOTSTRAP_REPS), "L"))

# Preserve the frozen historical gene identity while accepting the current
# HGNC/annotation-package symbol. SDPR and CAVIN2 refer to the same gene; this
# is an annotation alias repair, not a replacement CAF signature.
symbol_line <- grep("^probe_map\\$SYMBOL <- toupper\\(probe_map\\$SYMBOL\\)$", base_lines)
if (length(symbol_line) != 1L) stop("Could not locate the probe-symbol normalization line in the base script.")
alias_patch <- c(
  "probe_map$annotation_symbol_current <- probe_map$SYMBOL",
  "FROZEN_TO_CURRENT_GENE_ALIAS <- c(SDPR = \"CAVIN2\")",
  "for (frozen_name in names(FROZEN_TO_CURRENT_GENE_ALIAS)) {",
  "  current_name <- FROZEN_TO_CURRENT_GENE_ALIAS[[frozen_name]]",
  "  probe_map$SYMBOL[probe_map$SYMBOL == current_name] <- frozen_name",
  "}"
)
base_lines <- append(base_lines, alias_patch, after = symbol_line)

patched_base <- tempfile(pattern = "RUN_GSE39582_BASE_WITH_CAF_", fileext = ".R")
writeLines(base_lines, patched_base, useBytes = TRUE)
on.exit(unlink(patched_base), add = TRUE)

message("Running the locked base analysis with the frozen CAF_PROXY_30 resource...")
sys.source(patched_base, envir = globalenv(), keep.source = TRUE)

# Objects below are created by the base script. Fail closed if its internals
# change, rather than silently producing a different analysis.
required_objects <- c(
  "OUTPUT_DIR", "analysis_data", "score_columns", "ECM_MODULES",
  "BEST_MODULE", "safe_cox_fit", "cohort_subset", "cox_results",
  "bootstrap_results", "caf_available"
)
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), envir = globalenv(), inherits = FALSE)]
if (length(missing_objects)) stop("Base script did not create required objects: ", paste(missing_objects, collapse = ", "))
if (!isTRUE(caf_available)) stop("Frozen CAF_PROXY_30 was supplied but CAF scoring was not available; inspect the base run log.")

TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

write_csv_local <- function(x, filename) {
  utils::write.csv(x, file.path(TABLE_DIR, filename), row.names = FALSE, na = "")
}

rbind_fill <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, rows)
}

# Record the one approved historical/current symbol mapping.
write_csv_local(
  data.frame(
    frozen_symbol = "SDPR",
    current_annotation_symbol = "CAVIN2",
    relationship = "same_gene_historical_alias",
    signature_membership_changed = FALSE,
    stringsAsFactors = FALSE
  ),
  "CAF_frozen_gene_alias_audit.csv"
)

# ------------------------ VALIDATION COX SENSITIVITIES ------------------------
# These models improve interpretability/estimability. They do not reselect a
# module, alter a score, optimize a cut point, or convert a null Validation
# result into confirmation.

validation_scores <- c(names(ECM_MODULES), "ECM21_gene_equal", "ECM21_module_balanced", "ECM19_structural")
strata <- survival::strata

validation_formulas <- list(
  univariable = survival::Surv(rfs_time_months, rfs_event) ~ score_value,
  age_mmr = survival::Surv(rfs_time_months, rfs_event) ~ score_value + age + mmr,
  clinical_primary = survival::Surv(rfs_time_months, rfs_event) ~ score_value + age + sex + stage + mmr,
  clinical_stratified = survival::Surv(rfs_time_months, rfs_event) ~ score_value + age + mmr + strata(sex, stage),
  caf_age_mmr = survival::Surv(rfs_time_months, rfs_event) ~ score_value + age + mmr + caf_value,
  caf_stratified = survival::Surv(rfs_time_months, rfs_event) ~ score_value + age + mmr + caf_value + strata(sex, stage)
)

fit_custom_cox <- function(dat, score_name, formula) {
  d <- dat
  d$score_value <- d[[score_name]]
  req <- unique(all.vars(formula))
  d <- d[stats::complete.cases(d[, req, drop = FALSE]) & d$rfs_time_months > 0 & d$rfs_event %in% c(0, 1), , drop = FALSE]
  d <- droplevels(d)
  if (nrow(d) < 30L || sum(d$rfs_event) < 10L) {
    return(list(ok = FALSE, reason = "insufficient_rows_or_events", n = nrow(d), events = sum(d$rfs_event)))
  }
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      survival::coxph(formula, data = d, ties = "efron", x = TRUE, y = TRUE, model = TRUE),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, reason = conditionMessage(fit), n = nrow(d), events = sum(d$rfs_event)))
  }
  if (!"score_value" %in% names(stats::coef(fit))) {
    return(list(ok = FALSE, reason = "score_coefficient_missing", n = nrow(d), events = sum(d$rfs_event)))
  }
  beta <- unname(stats::coef(fit)["score_value"])
  se <- sqrt(unname(stats::vcov(fit)["score_value", "score_value"]))
  if (!is.finite(beta) || !is.finite(se) || se <= 0 || abs(beta) > 20 || se > 20) {
    return(list(ok = FALSE, reason = paste(c("score_estimability_gate_failed", warnings), collapse = " | "), n = nrow(d), events = sum(d$rfs_event)))
  }
  lp <- stats::predict(fit, type = "lp")
  cindex <- tryCatch(
    survival::concordance(survival::Surv(d$rfs_time_months, d$rfs_event) ~ lp, reverse = TRUE)$concordance,
    error = function(e) NA_real_
  )
  ph <- tryCatch(survival::cox.zph(fit, transform = "km"), error = function(e) NULL)
  list(
    ok = TRUE, fit = fit, data = d, n = nrow(d), events = sum(d$rfs_event),
    beta = beta, se = se, cindex = unname(cindex), warnings = paste(unique(warnings), collapse = " | "),
    score_ph_p = if (!is.null(ph) && "score_value" %in% rownames(ph$table)) ph$table["score_value", "p"] else NA_real_,
    global_ph_p = if (!is.null(ph) && "GLOBAL" %in% rownames(ph$table)) ph$table["GLOBAL", "p"] else NA_real_
  )
}

validation_dat <- cohort_subset("Validation")
validation_rows <- list()
for (model_name in names(validation_formulas)) {
  for (score_name in validation_scores) {
    z <- fit_custom_cox(validation_dat, score_name, validation_formulas[[model_name]])
    if (!isTRUE(z$ok)) {
      validation_rows[[length(validation_rows) + 1L]] <- data.frame(
        model = model_name, score = score_name, status = "failed", reason = z$reason,
        n = z$n, events = z$events, HR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
        p_value = NA_real_, Cindex = NA_real_, score_PH_p = NA_real_, global_PH_p = NA_real_,
        warnings = "", stringsAsFactors = FALSE
      )
    } else {
      validation_rows[[length(validation_rows) + 1L]] <- data.frame(
        model = model_name, score = score_name, status = "completed", reason = "",
        n = z$n, events = z$events, HR = exp(z$beta), CI_lower = exp(z$beta - 1.96 * z$se),
        CI_upper = exp(z$beta + 1.96 * z$se), p_value = 2 * stats::pnorm(-abs(z$beta / z$se)),
        Cindex = z$cindex, score_PH_p = z$score_ph_p, global_PH_p = z$global_ph_p,
        warnings = z$warnings, stringsAsFactors = FALSE
      )
    }
  }
}
validation_results <- do.call(rbind, validation_rows)
validation_results$BH_FDR_M1_to_M4 <- NA_real_
for (model_name in unique(validation_results$model)) {
  idx <- which(
    validation_results$model == model_name &
      validation_results$score %in% names(ECM_MODULES) &
      is.finite(validation_results$p_value)
  )
  if (length(idx)) validation_results$BH_FDR_M1_to_M4[idx] <- stats::p.adjust(validation_results$p_value[idx], method = "BH")
}
write_csv_local(validation_results, "Validation_sensitivity_Cox_results.csv")

# ------------------ DISCOVERY-VALIDATION EFFECT INTERACTION ------------------
# A score-by-cohort interaction asks whether the transferred Validation effect
# differs from Discovery. It is an effect-heterogeneity diagnostic, not a way to
# pool away a failed Validation confirmation.

interaction_fit <- function(dat, score_name, add_caf = FALSE) {
  d <- dat
  d$score_value <- d[[score_name]]
  d$dataset <- stats::relevel(factor(d$dataset), ref = "Discovery")
  rhs <- if (add_caf) {
    "score_value * dataset + age + mmr + caf_value + strata(sex, stage)"
  } else {
    "score_value * dataset + age + mmr + strata(sex, stage)"
  }
  formula <- stats::as.formula(paste("survival::Surv(rfs_time_months, rfs_event) ~", rhs))
  req <- unique(all.vars(formula))
  d <- droplevels(d[stats::complete.cases(d[, req, drop = FALSE]) & d$rfs_time_months > 0 & d$rfs_event %in% c(0, 1), , drop = FALSE])
  fit <- tryCatch(survival::coxph(formula, data = d, ties = "efron", x = TRUE, y = TRUE), error = function(e) e)
  if (inherits(fit, "error")) return(data.frame(status = "failed", reason = conditionMessage(fit)))
  int_name <- grep("^score_value:datasetValidation$|^datasetValidation:score_value$", names(stats::coef(fit)), value = TRUE)
  if (length(int_name) != 1L) return(data.frame(status = "failed", reason = "interaction_term_missing"))
  b <- stats::coef(fit)
  v <- stats::vcov(fit)
  b_d <- unname(b["score_value"])
  b_i <- unname(b[int_name])
  var_d <- unname(v["score_value", "score_value"])
  var_v <- var_d + unname(v[int_name, int_name]) + 2 * unname(v["score_value", int_name])
  se_d <- sqrt(var_d)
  se_v <- sqrt(var_v)
  data.frame(
    status = "completed", reason = "", n = nrow(d), events = sum(d$rfs_event),
    discovery_HR = exp(b_d), discovery_CI_lower = exp(b_d - 1.96 * se_d), discovery_CI_upper = exp(b_d + 1.96 * se_d),
    validation_HR = exp(b_d + b_i), validation_CI_lower = exp(b_d + b_i - 1.96 * se_v), validation_CI_upper = exp(b_d + b_i + 1.96 * se_v),
    interaction_logHR = b_i,
    interaction_p = summary(fit)$coefficients[int_name, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

pooled_dat <- cohort_subset("Pooled")
interaction_rows <- list()
for (add_caf in c(FALSE, TRUE)) {
  model_name <- if (add_caf) "clinical_plus_CAF" else "clinical"
  for (score_name in validation_scores) {
    z <- interaction_fit(pooled_dat, score_name, add_caf)
    z$model <- model_name
    z$score <- score_name
    interaction_rows[[length(interaction_rows) + 1L]] <- z
  }
}
interaction_results <- rbind_fill(interaction_rows)
interaction_results$BH_FDR_interaction_M1_to_M4 <- NA_real_
for (model_name in unique(interaction_results$model)) {
  idx <- which(
    interaction_results$model == model_name & interaction_results$score %in% names(ECM_MODULES) &
      is.finite(interaction_results$interaction_p)
  )
  if (length(idx)) interaction_results$BH_FDR_interaction_M1_to_M4[idx] <- stats::p.adjust(interaction_results$interaction_p[idx], method = "BH")
}
write_csv_local(interaction_results, "Discovery_Validation_score_interaction_results.csv")

# --------------- VALIDATION EVENT/COVARIATE-STRATIFIED BOOTSTRAP -------------
# The ordinary bootstrap remains primary. This sensitivity preserves the exact
# event, sex, stage and MMR cell counts in each replicate to reduce failures
# caused solely by resampling sparse Validation covariate cells.

validation_comparisons <- data.frame(
  comparison = c(
    "ECM21_gene_equal_minus_best_module_frozen",
    "ECM21_gene_equal_minus_ECM19_structural",
    "ECM21_gene_equal_minus_ECM21_module_balanced"
  ),
  score_A = "ECM21_gene_equal",
  score_B = c(BEST_MODULE, "ECM19_structural", "ECM21_module_balanced"),
  stringsAsFactors = FALSE
)

validation_bootstrap_models <- c("clinical_stratified", "caf_stratified")

block_bootstrap_one_model <- function(dat, model_key, reps, seed) {
  scores <- unique(c(validation_comparisons$score_A, validation_comparisons$score_B))
  required <- c("rfs_time_months", "rfs_event", "age", "sex", "stage", "mmr", scores)
  if (model_key == "caf_stratified") required <- c(required, "caf_value")
  d <- droplevels(dat[stats::complete.cases(dat[, unique(required), drop = FALSE]) & dat$rfs_time_months > 0 & dat$rfs_event %in% c(0, 1), , drop = FALSE])
  block <- interaction(d$rfs_event, d$sex, d$stage, d$mmr, drop = TRUE, lex.order = TRUE)
  block_indices <- split(seq_len(nrow(d)), block)

  apparent <- setNames(rep(NA_real_, length(scores)), scores)
  for (s in scores) {
    fit <- safe_cox_fit(d, score_name = s, model_key = model_key, include_score = TRUE)
    if (isTRUE(fit$ok)) apparent[s] <- fit$cindex
  }

  set.seed(seed)
  boot <- matrix(NA_real_, nrow = reps, ncol = length(scores), dimnames = list(NULL, scores))
  for (b in seq_len(reps)) {
    idx <- unlist(lapply(block_indices, function(ii) sample(ii, length(ii), replace = TRUE)), use.names = FALSE)
    db <- droplevels(d[idx, , drop = FALSE])
    for (s in scores) {
      fit <- safe_cox_fit(db, score_name = s, model_key = model_key, include_score = TRUE)
      if (isTRUE(fit$ok)) boot[b, s] <- fit$cindex
    }
    if (b %% 100L == 0L) message("Validation block bootstrap ", model_key, ": ", b, "/", reps)
  }
  list(data = d, apparent = apparent, boot = boot, n_blocks = length(block_indices))
}

validation_boot_rows <- list()
for (m in seq_along(validation_bootstrap_models)) {
  model_key <- validation_bootstrap_models[m]
  bobj <- block_bootstrap_one_model(
    validation_dat, model_key,
    VALIDATION_STRATIFIED_BOOTSTRAP_REPS,
    VALIDATION_STRATIFIED_BOOTSTRAP_SEED + m
  )
  for (i in seq_len(nrow(validation_comparisons))) {
    a <- validation_comparisons$score_A[i]
    b <- validation_comparisons$score_B[i]
    delta <- bobj$boot[, a] - bobj$boot[, b]
    good <- is.finite(delta)
    q <- if (sum(good) >= 100L) stats::quantile(delta[good], c(0.025, 0.5, 0.975), names = FALSE, type = 6) else rep(NA_real_, 3L)
    apparent_delta <- bobj$apparent[a] - bobj$apparent[b]
    reliable <- mean(good) >= MIN_BOOTSTRAP_SUCCESS_FRACTION
    validation_boot_rows[[length(validation_boot_rows) + 1L]] <- data.frame(
      cohort = "Validation", model = model_key,
      bootstrap_type = "event_sex_stage_MMR_block_stratified_sensitivity",
      comparison = validation_comparisons$comparison[i], score_A = a, score_B = b,
      n = nrow(bobj$data), events = sum(bobj$data$rfs_event), blocks = bobj$n_blocks,
      reps_requested = VALIDATION_STRATIFIED_BOOTSTRAP_REPS,
      reps_successful = sum(good), success_fraction = mean(good),
      delta_Cindex_apparent = apparent_delta, delta_Cindex_boot_median = q[2],
      CI_lower = q[1], CI_upper = q[3], probability_delta_gt_0 = if (sum(good)) mean(delta[good] > 0) else NA_real_,
      reliable_bootstrap = reliable,
      positive_CI_supported_gain = isTRUE(reliable) && is.finite(apparent_delta) && apparent_delta > 0 && is.finite(q[1]) && q[1] > 0,
      primary_inference = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
validation_boot_results <- do.call(rbind, validation_boot_rows)
write_csv_local(validation_boot_results, "Validation_event_covariate_stratified_bootstrap_delta_Cindex.csv")

# -------------------------- MACHINE-READABLE DECISION ------------------------
ordinary_validation <- bootstrap_results[
  bootstrap_results$cohort == "Validation" &
    bootstrap_results$model %in% c("clinical_stratified", "caf_stratified") &
    bootstrap_results$comparison %in% c(
      "ECM21_gene_equal_minus_best_module_frozen",
      "ECM21_gene_equal_minus_ECM19_structural"
    ), , drop = FALSE
]

validation_decision <- data.frame(
  question = c(
    "Was Validation analysis technically completed?",
    "Did the ordinary Validation analysis confirm an association?",
    "Did the ordinary paired bootstrap support ECM21 enhancement?",
    "Can the block-stratified sensitivity replace the ordinary bootstrap?",
    "Was module-level CAF adjustment completed in this continuation?"
  ),
  answer = c(
    "YES",
    "NO unless the exported frozen results meet the prespecified CI/FDR rule",
    if (isTRUE(nrow(ordinary_validation) > 0L) && isTRUE(all(ordinary_validation$delta_Cindex_apparent > 0 & ordinary_validation$CI_lower > 0))) "YES" else "NO",
    "NO; sensitivity analysis only",
    "YES"
  ),
  interpretation = c(
    "Technical completion is distinct from statistical confirmation.",
    "A directionally positive HR with a CI crossing 1 is not confirmation.",
    "Enhancement requires positive deltas, CIs above zero, stable direction and reliable bootstrap in every required cohort.",
    "It diagnoses sparse-cell bootstrap instability and cannot be substituted post hoc for the locked primary analysis.",
    "Use caf_primary/caf_stratified rows to determine whether each module adds information beyond the frozen CAF proxy."
  ),
  stringsAsFactors = FALSE
)
write_csv_local(validation_decision, "VALIDATION_AND_CAF_DECISION_SUMMARY.csv")

readme_lines <- c(
  "GSE39582 MODULE CAF AND VALIDATION CONTINUATION",
  "",
  "Primary facts:",
  paste0("- Best module remains frozen from Discovery: ", BEST_MODULE, "."),
  "- Validation was already technically completed; non-confirmation must not be described as a missing analysis.",
  "- This run adds the exact frozen CAF_PROXY_30 to every module/aggregate model.",
  "- SDPR is mapped to the current annotation symbol CAVIN2 as the same gene; the 30-gene membership is unchanged.",
  "- The event/sex/stage/MMR block bootstrap is a labelled sensitivity analysis only.",
  "- Do not write 'enhanced' unless the original locked claim gate is TRUE.",
  "- Do not infer equivalence merely because a confidence interval crosses zero.",
  "",
  "Key outputs:",
  "- GSE39582_module_and_aggregate_Cox_results.csv: includes caf_primary and caf_stratified rows.",
  "- paired_bootstrap_delta_Cindex_results.csv: includes ordinary CAF-adjusted paired bootstraps.",
  "- CLAIM_GATE_enhancement_decision.csv: locked enhancement decisions.",
  "- Validation_sensitivity_Cox_results.csv: univariable, reduced, ordinary, stratified and CAF-adjusted Validation models.",
  "- Discovery_Validation_score_interaction_results.csv: score-by-cohort heterogeneity diagnostics.",
  "- Validation_event_covariate_stratified_bootstrap_delta_Cindex.csv: sparse-cell bootstrap sensitivity.",
  "- VALIDATION_AND_CAF_DECISION_SUMMARY.csv: machine-readable claim boundaries."
)
writeLines(readme_lines, file.path(OUTPUT_DIR, "README_CAF_AND_VALIDATION_CONTINUATION.txt"), useBytes = TRUE)

utils::write.csv(
  data.frame(
    component = c("base_clinical_and_CAF_run", "validation_sensitivity_cox", "cohort_interaction", "validation_block_bootstrap"),
    status = "completed",
    stringsAsFactors = FALSE
  ),
  file.path(OUTPUT_DIR, "RUN_STATUS_CAF_AND_VALIDATION_CONTINUATION.csv"),
  row.names = FALSE
)

capture.output(utils::sessionInfo(), file = file.path(OUTPUT_DIR, "R_sessionInfo_CAF_AND_VALIDATION_CONTINUATION.txt"))

message("Completed. Output directory: ", OUTPUT_DIR)
