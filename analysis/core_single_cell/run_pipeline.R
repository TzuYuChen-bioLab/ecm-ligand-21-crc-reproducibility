# Run corrected scripts in isolated R sessions.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "run_pipeline.R"
pipeline_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
from_step <- as.integer(Sys.getenv("CRC_FROM_STEP", unset = "0"))
to_step <- as.integer(Sys.getenv("CRC_TO_STEP", unset = "13"))
if (!is.finite(from_step) || !is.finite(to_step) || from_step < 0 || to_step > 13 || from_step > to_step) {
  stop("Set CRC_FROM_STEP/CRC_TO_STEP to a valid inclusive range from 0 to 13.")
}
if (!nzchar(Sys.getenv("CRC_PROJECT_ROOT", unset = ""))) {
  stop("Set CRC_PROJECT_ROOT to the directory containing 00_raw_data before running the pipeline.")
}
rscript <- file.path(R.home("bin"), "Rscript")
for (step in seq.int(from_step, to_step)) {
  script <- file.path(pipeline_dir, "scripts", paste0(step, ".R"))
  message("\n===== Running corrected step ", step, ": ", script, " =====")
  status <- system2(rscript, args = shQuote(script), stdout = "", stderr = "")
  if (!identical(status, 0L)) stop("Corrected step ", step, " failed with status ", status, ".")
}
message("Corrected pipeline completed steps ", from_step, "–", to_step, ".")
