# Executable synthetic test for the Step 10 Cox estimability gate.
# Run from the corrected_pipeline directory with:
# Rscript tests/test_step10_estimability.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
test_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_step10_estimability.R"
root <- normalizePath(file.path(dirname(test_file), ".."), mustWork = TRUE)

PROJECT_ROOT <- tempdir()
require_packages <- function(packages, bioconductor = FALSE) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Synthetic test requires: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}
write_csv_atomic <- function(x, path, row.names = FALSE) {
  ensure_dir(dirname(path))
  utils::write.csv(x, path, row.names = row.names, na = "")
  invisible(path)
}
source(file.path(root, "scripts", "helpers_bulk.R"))
require_packages(c("survival", "broom"))
strata <- survival::strata

set.seed(1234)
n <- 80L
stage <- factor(c(rep("I", 20), rep("II", 60)), levels = c("I", "II"))
event <- c(rep(0L, 20), stats::rbinom(60, 1, 0.40))
if (sum(event) < 10) stop("Synthetic event generation unexpectedly produced too few events.")
dat <- data.frame(
  sample = paste0("S", seq_len(n)),
  time = stats::rexp(n, rate = 0.05) + 0.1,
  event = event,
  score_z = stats::rnorm(n),
  stage = stage,
  stringsAsFactors = FALSE
)

bad_dir <- tempfile("cox_bad_")
bad <- fit_cox_with_diagnostics(
  dat, survival::Surv(time, event) ~ score_z + stage,
  "synthetic_zero_event_coefficient", bad_dir
)
stopifnot(!bad$estimable)
stopifnot(bad$status$status == "not_estimable")
stopifnot(grepl("zero-event coefficient level", bad$status$reason, fixed = TRUE))

stable_dir <- tempfile("cox_stable_")
stable <- fit_cox_with_diagnostics(
  dat, survival::Surv(time, event) ~ score_z + strata(stage),
  "synthetic_stage_stratified", stable_dir
)
stopifnot(stable$estimable)
dfbeta <- utils::read.csv(file.path(stable_dir, "synthetic_stage_stratified_dfbeta.csv"), check.names = FALSE)
stopifnot(identical(names(dfbeta), c("sample", "score_z")))
stopifnot(identical(dfbeta$sample, dat$sample))

message("PASS: zero-event coefficient model rejected; stratified sensitivity estimated; DFBETA labels preserved.")
