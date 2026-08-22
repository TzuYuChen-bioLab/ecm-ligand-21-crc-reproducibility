#!/usr/bin/env Rscript

# Optional helper: export the four Human Matrisome gene sets used as the
# category-matched random-panel background. No outcome data are accessed.

options(stringsAsFactors = FALSE)

parse_cli <- function(x) {
  out <- list()
  for (item in x) {
    if (!grepl("^--[^=]+=", item)) next
    out[[sub("^--([^=]+)=.*$", "\\1", item)]] <- sub("^--[^=]+=", "", item)
  }
  out
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
out_path <- if (!is.null(args$out)) args$out else "naba_background_genes.csv"

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop(
    "Package 'msigdbr' is required. Install it explicitly with install.packages('msigdbr'), ",
    "then rerun this helper."
  )
}

msig <- msigdbr::msigdbr(species = "Homo sapiens")
set_col <- intersect(c("gs_name", "gene_set_name"), names(msig))[1]
symbol_col <- intersect(c("gene_symbol", "human_gene_symbol"), names(msig))[1]
if (is.na(set_col) || is.na(symbol_col)) {
  stop("Could not identify gene-set and gene-symbol columns in the installed msigdbr version.")
}

target_sets <- c(
  "NABA_COLLAGENS",
  "NABA_ECM_GLYCOPROTEINS",
  "NABA_PROTEOGLYCANS",
  "NABA_SECRETED_FACTORS"
)

keep <- msig[[set_col]] %in% target_sets
background <- unique(data.frame(
  gene = toupper(trimws(as.character(msig[[symbol_col]][keep]))),
  naba_gene_set = as.character(msig[[set_col]][keep]),
  stringsAsFactors = FALSE
))
background <- background[nzchar(background$gene) & !is.na(background$gene), , drop = FALSE]

missing_sets <- setdiff(target_sets, unique(background$naba_gene_set))
if (length(missing_sets)) {
  stop(
    "The installed msigdbr release did not contain: ", paste(missing_sets, collapse = ", "),
    ". Do not silently substitute a different collection."
  )
}

dir.create(dirname(normalizePath(out_path, winslash = "/", mustWork = FALSE)), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(background[order(background$naba_gene_set, background$gene), ], out_path, row.names = FALSE, na = "")
message("Created ", out_path, " with ", nrow(background), " unique gene-to-set rows.")
