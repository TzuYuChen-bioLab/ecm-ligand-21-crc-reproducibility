PIPELINE_ROOT <- normalizePath(Sys.getenv("CRC_TRIANGULATION_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
source(file.path(PIPELINE_ROOT, "config.R"))
source(file.path(PIPELINE_ROOT, "R", "functions.R"))

cran <- c("data.table", "Matrix", "ggplot2", "patchwork", "arrow", "readxl", "R.utils", "RANN", "survival", "Seurat")
bioc <- c("edgeR", "limma", "GEOquery", "AnnotationDbi", "hgu133plus2.db")
optional <- c("metafor")
all_packages <- c(cran, bioc, optional)
installed <- vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)

report <- data.frame(
  package = all_packages,
  source = c(rep("CRAN", length(cran)), rep("Bioconductor", length(bioc)), rep("CRAN", length(optional))),
  installed = installed,
  version = vapply(
    all_packages,
    function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_,
    character(1)
  )
)
write_csv(report, file.path(LOG_ROOT, "package_check.csv"))
print(report)
if (!all(installed)) stop("Install the missing packages manually; see README_FIRST_CN.md.")
write_status("00_package_check", "PASS", "All required packages are installed.")
