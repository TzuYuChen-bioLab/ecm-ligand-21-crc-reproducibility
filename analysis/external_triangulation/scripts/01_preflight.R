PIPELINE_ROOT <- normalizePath(Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))
require_packages(c("data.table"))

manifest <- read.csv(file.path(PIPELINE_ROOT, "DATA_DOWNLOAD_MANIFEST.csv"), check.names = FALSE)
manifest$local_path <- file.path(PIPELINE_ROOT, manifest$relative_path)
manifest$exists <- file.exists(manifest$local_path)
manifest$size_bytes <- ifelse(manifest$exists, file.info(manifest$local_path)$size, NA_real_)
manifest$required_now <- manifest$priority %in% c("CORE", "HIGH")

write_csv(manifest, file.path(LOG_ROOT, "preflight_file_audit.csv"))
print(manifest[, c("dataset", "file", "priority", "exists", "size_bytes")])

signature <- read.csv(SIGNATURE_FILE, check.names = FALSE)
target <- read.csv(TARGET_FILE, check.names = FALSE)
stopifnot(nrow(signature) == 21L, length(unique(signature$gene)) == 21L)
stopifnot(nrow(target) == 50L, length(unique(target$gene)) == 50L)
stopifnot(length(intersect(signature$gene, target$gene)) == 0L)

missing_core <- manifest$file[manifest$required_now & !manifest$exists]
if (length(missing_core)) {
  write_status("01_preflight", "BLOCKED", paste(missing_core, collapse = "; "))
  stop("Core/high-priority files are missing:\n", paste("-", missing_core, collapse = "\n"))
}
write_status("01_preflight", "PASS", "All core/high-priority files and frozen resources are present.")
cat("Preflight passed. You may run scripts 02-05 and 08. Optional scripts 06-07 depend on optional files.\n")

