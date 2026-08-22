# Shared bulk-transcriptomic and survival helpers.
# Callers must source config.R first; failing explicitly is safer than guessing a
# working directory and silently reading a different project.
if (!exists("PROJECT_ROOT", inherits = TRUE)) {
  stop("Source scripts/config.R before scripts/helpers_bulk.R.")
}

parse_geo_characteristics <- function(pdata) {
  char_cols <- grep("^characteristics_ch1", colnames(pdata), value = TRUE, ignore.case = TRUE)
  parsed <- list()
  if (length(char_cols) == 0) return(pdata)

  for (i in seq_len(nrow(pdata))) {
    fields <- unlist(lapply(char_cols, function(col) {
      value <- as.character(pdata[i, col])
      if (is.na(value) || !nzchar(value)) return(character(0))
      trimws(unlist(strsplit(value, ";\\s*", perl = TRUE)))
    }), use.names = FALSE)
    for (field in fields) {
      if (!grepl(":", field, fixed = TRUE)) next
      key <- clean_key(sub(":.*$", "", field))
      value <- trimws(sub("^[^:]*:\\s*", "", field))
      if (!nzchar(key)) next
      column <- paste0("parsed_", key)
      if (is.null(parsed[[column]])) parsed[[column]] <- rep(NA_character_, nrow(pdata))
      if (is.na(parsed[[column]][i]) || !nzchar(parsed[[column]][i])) parsed[[column]][i] <- value
    }
  }

  if (length(parsed) > 0) {
    parsed_df <- as.data.frame(parsed, stringsAsFactors = FALSE, check.names = FALSE)
    rownames(parsed_df) <- rownames(pdata)
    pdata <- cbind(pdata, parsed_df)
  }
  pdata
}

as_expression_set_compat <- function(x) {
  if (inherits(x, "ExpressionSet")) return(x)
  if (is.list(x)) {
    hit <- which(vapply(x, inherits, logical(1), what = "ExpressionSet"))
    if (length(hit) > 0) return(x[[hit[1]]])
  }
  stop("GEOquery did not return an ExpressionSet.")
}

read_local_geo <- function(gse_id) {
  require_packages(c("GEOquery", "Biobase"), bioconductor = TRUE)
  raw_dir <- project_path("00_raw_data", "bulk", gse_id)
  matrix_file <- file.path(raw_dir, paste0(gse_id, "_series_matrix.txt.gz"))
  if (!file.exists(matrix_file)) {
    candidates <- list.files(raw_dir, pattern = "series_matrix.*\\.txt\\.gz$", full.names = TRUE, ignore.case = TRUE)
    if (length(candidates) == 1) matrix_file <- candidates[1]
  }
  if (!file.exists(matrix_file)) stop("Missing GEO series matrix for ", gse_id, ": ", matrix_file)
  obj <- GEOquery::getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE, parseCharacteristics = FALSE)
  list(eset = as_expression_set_compat(obj), matrix_file = matrix_file)
}

classify_primary_tumour <- function(pdata, gse_id) {
  all_text <- apply(pdata, 1, function(row) paste(tolower(as.character(row)), collapse = " | "))
  keep <- rep(TRUE, nrow(pdata))
  reason <- rep("assumed primary tumour", nrow(pdata))

  if (gse_id == "GSE39582") {
    role_col <- c("parsed_dataset", "dataset:ch1")
    role_col <- role_col[role_col %in% colnames(pdata)][1]
    if (is.na(role_col)) stop("GSE39582 dataset field not found; cannot exclude the 19 non-tumour mucosa samples safely.")
    role <- tolower(trimws(as.character(pdata[[role_col]])))
    keep <- role %in% c("discovery", "validation")
    reason <- ifelse(keep, role, "excluded non-tumoral mucosa")
  } else if (gse_id == "GSE33113") {
    tissue_col <- c("parsed_tissue", "tissue:ch1")
    tissue_col <- tissue_col[tissue_col %in% colnames(pdata)][1]
    if (is.na(tissue_col)) stop("GSE33113 tissue field not found; cannot exclude normal colon mucosa safely.")
    tissue <- tolower(trimws(as.character(pdata[[tissue_col]])))
    keep <- grepl("primary.*tumou?r|tumou?r.*resection", tissue) & !grepl("normal|mucosa", tissue)
    reason <- ifelse(keep, "primary tumour", "excluded normal mucosa")
  } else {
    non_tumour <- grepl("non[- ]?tumou?r|normal (colon|colorectal|mucosa|tissue)|adjacent normal", all_text)
    keep <- !non_tumour
    reason <- ifelse(keep, "not labelled non-tumour", "excluded non-tumour")
  }

  data.frame(sample = rownames(pdata), keep_primary_tumour = keep, classification_reason = reason, stringsAsFactors = FALSE)
}

map_probes_to_symbols <- function(probes, platform_id, fdata) {
  require_packages(c("AnnotationDbi", "hgu133plus2.db"), bioconductor = TRUE)
  is_gpl570 <- grepl("GPL570|hgu133plus2|HG-U133", paste(platform_id, collapse = " "), ignore.case = TRUE)
  if (!is_gpl570) stop("Only GPL570/HG-U133 Plus 2.0 is supported by this corrected helper. Platform: ", platform_id)
  mapped <- AnnotationDbi::select(
    hgu133plus2.db::hgu133plus2.db,
    keys = unique(as.character(probes)), columns = "SYMBOL", keytype = "PROBEID"
  )
  mapped <- mapped[!is.na(mapped$SYMBOL) & nzchar(mapped$SYMBOL), c("PROBEID", "SYMBOL")]
  # Ambiguous one-probe-to-many-symbol annotations are excluded rather than
  # resolved by an arbitrary first row.
  symbol_count <- vapply(split(mapped$SYMBOL, mapped$PROBEID), function(x) length(unique(x)), integer(1))
  unambiguous <- names(symbol_count)[symbol_count == 1]
  mapped <- mapped[mapped$PROBEID %in% unambiguous, , drop = FALSE]
  mapped <- mapped[!duplicated(mapped$PROBEID), , drop = FALSE]
  symbols <- setNames(rep(NA_character_, length(probes)), probes)
  symbols[mapped$PROBEID] <- mapped$SYMBOL
  symbols
}

collapse_probes_by_reference <- function(expr, symbols, reference_samples) {
  reference_samples <- intersect(reference_samples, colnames(expr))
  if (length(reference_samples) < 10) stop("Fewer than 10 reference samples for probe selection.")
  keep <- !is.na(symbols[rownames(expr)]) & nzchar(symbols[rownames(expr)])
  expr <- expr[keep, , drop = FALSE]
  symbols_use <- symbols[rownames(expr)]
  probe_iqr <- apply(expr[, reference_samples, drop = FALSE], 1, stats::IQR, na.rm = TRUE)
  mapping <- data.frame(probe = rownames(expr), symbol = unname(symbols_use), reference_iqr = probe_iqr, stringsAsFactors = FALSE)
  mapping <- mapping[order(mapping$symbol, -mapping$reference_iqr, mapping$probe), ]
  mapping <- mapping[!duplicated(mapping$symbol), ]
  out <- expr[mapping$probe, , drop = FALSE]
  rownames(out) <- mapping$symbol
  list(expr = out, mapping = mapping)
}

process_geo_bulk_corrected <- function(gse_id, force = FALSE, requantile = FALSE) {
  require_packages(c("limma"), bioconductor = TRUE)
  out_dir <- ensure_dir(project_path("01_processed_data", "bulk_corrected", gse_id))
  expr_path <- file.path(out_dir, paste0(gse_id, "_primary_tumour_expr_gene.rds"))
  pheno_path <- file.path(out_dir, paste0(gse_id, "_primary_tumour_pheno.rds"))
  if (!force && file.exists(expr_path) && file.exists(pheno_path)) {
    return(list(expr = readRDS(expr_path), pheno = readRDS(pheno_path), out_dir = out_dir))
  }

  geo <- read_local_geo(gse_id)
  expr <- Biobase::exprs(geo$eset)
  pdata <- parse_geo_characteristics(Biobase::pData(geo$eset))
  fdata <- Biobase::fData(geo$eset)
  platform_id <- Biobase::annotation(geo$eset)
  if (length(platform_id) == 0 || is.na(platform_id) || !nzchar(platform_id)) {
    platform_id <- unique(as.character(pdata$platform_id))
  }
  common <- intersect(colnames(expr), rownames(pdata))
  expr <- expr[, common, drop = FALSE]
  pdata <- pdata[common, , drop = FALSE]

  audit <- classify_primary_tumour(pdata, gse_id)
  write_csv_atomic(audit, file.path(out_dir, paste0(gse_id, "_tissue_filter_audit.csv")))
  keep_samples <- audit$sample[audit$keep_primary_tumour]
  if (gse_id == "GSE39582" && length(keep_samples) != 566) {
    stop("GSE39582 tumour filter retained ", length(keep_samples), " samples; expected 566 (443 discovery + 123 validation).")
  }
  expr <- expr[, keep_samples, drop = FALSE]
  pdata <- pdata[keep_samples, , drop = FALSE]

  q99 <- stats::quantile(expr, 0.99, na.rm = TRUE)
  if (is.finite(q99) && q99 > 100) expr <- log2(expr + 1)
  if (isTRUE(requantile)) expr <- limma::normalizeBetweenArrays(expr, method = "quantile")

  if (gse_id == "GSE39582") {
    split_col <- c("parsed_dataset", "dataset:ch1")
    split_col <- split_col[split_col %in% colnames(pdata)][1]
    reference_samples <- rownames(pdata)[tolower(trimws(as.character(pdata[[split_col]]))) == "discovery"]
    if (length(reference_samples) != 443) stop("Expected 443 GSE39582 discovery tumours; found ", length(reference_samples), ".")
  } else {
    reference_samples <- colnames(expr)
  }

  symbols <- map_probes_to_symbols(rownames(expr), platform_id, fdata)
  collapsed <- collapse_probes_by_reference(expr, symbols, reference_samples)
  expr_gene <- collapsed$expr
  write_csv_atomic(collapsed$mapping, file.path(out_dir, paste0(gse_id, "_probe_selection_reference_only.csv")))
  saveRDS(expr_gene, expr_path)
  saveRDS(pdata, pheno_path)
  write_csv_atomic(
    data.frame(gse_id = gse_id, n_primary_tumour = ncol(expr_gene), n_genes = nrow(expr_gene),
               reference_n = length(reference_samples), requantile = requantile,
               source_md5 = unname(tools::md5sum(geo$matrix_file))),
    file.path(out_dir, paste0(gse_id, "_processing_summary.csv"))
  )
  list(expr = expr_gene, pheno = pdata, out_dir = out_dir)
}

read_signature_definition <- function() {
  path <- project_path("03_cellchat_nichenet", "GSE132465_corrected", "signature", "CRC_ECM_LIGAND_21.csv")
  if (!file.exists(path)) stop("Run 8.R first; signature definition not found: ", path)
  sig <- utils::read.csv(path, stringsAsFactors = FALSE)
  unique(sig$gene[!is.na(sig$gene) & nzchar(sig$gene)])
}

score_fixed_signature <- function(gse_id, expr, pheno, genes, reference_samples = NULL) {
  genes_use <- intersect(genes, rownames(expr))
  if (length(genes_use) < ceiling(0.8 * length(genes))) {
    stop(gse_id, " matches only ", length(genes_use), "/", length(genes), " signature genes.")
  }
  if (is.null(reference_samples)) reference_samples <- colnames(expr)
  reference_samples <- intersect(reference_samples, colnames(expr))
  ref <- expr[genes_use, reference_samples, drop = FALSE]
  centres <- rowMeans(ref, na.rm = TRUE)
  scales <- apply(ref, 1, stats::sd, na.rm = TRUE)
  valid <- is.finite(scales) & scales > 0
  genes_use <- genes_use[valid]
  z <- sweep(expr[genes_use, , drop = FALSE], 1, centres[genes_use], "-")
  z <- sweep(z, 1, scales[genes_use], "/")
  score <- colMeans(z, na.rm = TRUE)
  out <- data.frame(sample = colnames(expr), ECM_LIGAND_21 = as.numeric(score), stringsAsFactors = FALSE)
  pheno_add <- pheno[out$sample, , drop = FALSE]
  pheno_add <- pheno_add[, setdiff(colnames(pheno_add), "sample"), drop = FALSE]
  out <- cbind(out, pheno_add)
  params <- data.frame(gene = genes_use, centre = centres[genes_use], scale = scales[genes_use], stringsAsFactors = FALSE)
  score_dir <- ensure_dir(project_path("04_bulk_model", "corrected_scores"))
  table_dir <- ensure_dir(project_path("08_tables", "bulk", "corrected_scores"))
  saveRDS(out, file.path(score_dir, paste0(gse_id, "_ECM_LIGAND_21_scores.rds")))
  write_csv_atomic(out, file.path(table_dir, paste0(gse_id, "_ECM_LIGAND_21_scores.csv")))
  write_csv_atomic(params, file.path(table_dir, paste0(gse_id, "_ECM_LIGAND_21_standardisation_parameters.csv")))
  write_csv_atomic(
    data.frame(gse_id = gse_id, signature_size = length(genes), matched = length(genes_use),
               reference_n = length(reference_samples), matched_genes = paste(genes_use, collapse = ";")),
    file.path(table_dir, paste0(gse_id, "_ECM_LIGAND_21_coverage.csv"))
  )
  out
}

first_existing <- function(data, candidates, required = FALSE, label = "column") {
  hit <- candidates[candidates %in% colnames(data)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("Required ", label, " not found. Tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

simplify_stage <- function(x) {
  raw <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(raw))
  out[grepl("DUKES.?A|^A$|^1$|^I$", raw)] <- "I"
  out[grepl("DUKES.?B|^B$|^2$|^II$", raw)] <- "II"
  out[grepl("DUKES.?C|^C$|^3$|^III$", raw)] <- "III"
  out[grepl("DUKES.?D|^D$|^4$|^IV$", raw)] <- "IV"
  factor(out, levels = c("I", "II", "III", "IV"), ordered = FALSE)
}

simplify_binary_status <- function(x) {
  raw <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(raw))
  out[raw %in% c("M", "MUT", "MUTATED", "MUTATION", "1", "POS", "POSITIVE", "DMMR", "MSI")] <- "Positive"
  out[raw %in% c("WT", "WILD TYPE", "WILDTYPE", "0", "NEG", "NEGATIVE", "PMMR", "MSS")] <- "Negative"
  factor(out, levels = c("Negative", "Positive"))
}

collapse_messages <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) "" else paste(x, collapse = " | ")
}

capture_coxph_fit <- function(formula, data) {
  warnings <- character(0)
  fit <- withCallingHandlers(
    tryCatch(
      survival::coxph(formula, data = data, ties = "efron", x = TRUE, y = TRUE, model = TRUE),
      error = function(e) e
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = unique(warnings))
}

cox_strata_variables <- function(formula) {
  labels <- attr(stats::terms(formula), "term.labels")
  strata_terms <- labels[grepl("strata\\s*\\(", labels)]
  if (length(strata_terms) == 0) return(character(0))
  unique(unlist(lapply(strata_terms, function(x) {
    all.vars(stats::as.formula(paste("~", x), env = environment(formula)))
  }), use.names = FALSE))
}

cox_event_level_audit <- function(data, formula, model_name) {
  vars <- setdiff(all.vars(formula), c("time", "event"))
  vars <- vars[vars %in% colnames(data)]
  strata_vars <- cox_strata_variables(formula)
  categorical <- vars[vapply(vars, function(v) {
    x <- data[[v]]
    is.factor(x) || is.character(x) || is.logical(x) || length(unique(stats::na.omit(x))) <= 10
  }, logical(1))]
  template <- data.frame(
    model = character(0), variable = character(0), level = character(0),
    model_role = character(0), n = integer(0), events = integer(0),
    non_events = integer(0), zero_event_level = logical(0),
    zero_non_event_level = logical(0), stringsAsFactors = FALSE
  )
  if (length(categorical) == 0 || nrow(data) == 0) return(template)
  out <- lapply(categorical, function(v) {
    value <- as.character(data[[v]])
    levels_present <- unique(value[!is.na(value)])
    if (length(levels_present) == 0) return(NULL)
    do.call(rbind, lapply(levels_present, function(level) {
      hit <- !is.na(value) & value == level
      n <- sum(hit)
      events <- sum(data$event[hit] == 1, na.rm = TRUE)
      data.frame(
        model = model_name, variable = v, level = level,
        model_role = if (v %in% strata_vars) "stratification" else "coefficient",
        n = n, events = events, non_events = n - events,
        zero_event_level = events == 0, zero_non_event_level = (n - events) == 0,
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) template else do.call(rbind, out)
}

assess_cox_fit <- function(captured, event_audit = NULL) {
  fit <- captured$fit
  warning_text <- collapse_messages(captured$warnings)
  if (inherits(fit, "error")) {
    return(list(
      estimable = FALSE, reason = paste0("coxph error: ", conditionMessage(fit)),
      warning_text = warning_text, parameters = NA_integer_,
      max_abs_log_coefficient = NA_real_, max_standard_error = NA_real_,
      zero_event_coefficient_levels = NA_integer_
    ))
  }
  coefficients <- stats::coef(fit)
  standard_errors <- suppressWarnings(sqrt(diag(stats::vcov(fit))))
  serious_warning <- grepl(
    "coefficient may be infinite|did not converge|ran out of iterations|loglik converged before",
    warning_text, ignore.case = TRUE
  )
  nonfinite <- any(!is.finite(coefficients)) || any(!is.finite(standard_errors))
  extreme <- any(abs(coefficients) > 10, na.rm = TRUE) || any(standard_errors > 10, na.rm = TRUE)
  finite_abs_coefficients <- abs(coefficients[is.finite(coefficients)])
  finite_standard_errors <- standard_errors[is.finite(standard_errors)]
  zero_rows <- if (is.null(event_audit) || nrow(event_audit) == 0) {
    data.frame(variable = character(0), level = character(0), stringsAsFactors = FALSE)
  } else {
    event_audit[event_audit$model_role == "coefficient" & event_audit$zero_event_level, , drop = FALSE]
  }
  zero_event_levels <- nrow(zero_rows)
  reasons <- character(0)
  if (serious_warning) reasons <- c(reasons, paste0("coxph convergence warning: ", warning_text))
  if (nonfinite) reasons <- c(reasons, "non-finite coefficient or standard error")
  if (extreme) reasons <- c(reasons, "absolute log coefficient or standard error exceeded 10")
  if (zero_event_levels > 0) {
    labels <- paste0(zero_rows$variable, "=", zero_rows$level)
    reasons <- c(reasons, paste0("zero-event coefficient level(s): ", paste(labels, collapse = "; ")))
  }
  list(
    estimable = length(reasons) == 0,
    reason = if (length(reasons) == 0) "estimated without convergence or separation flags" else collapse_messages(reasons),
    warning_text = warning_text,
    parameters = length(coefficients),
    max_abs_log_coefficient = if (length(finite_abs_coefficients)) max(finite_abs_coefficients) else NA_real_,
    max_standard_error = if (length(finite_standard_errors)) max(finite_standard_errors) else NA_real_,
    zero_event_coefficient_levels = zero_event_levels
  )
}

cox_status_table <- function(model_name, status, reason, n, events, parameters = NA_integer_,
                             warning_messages = "", max_abs_log_coefficient = NA_real_,
                             max_standard_error = NA_real_, zero_event_coefficient_levels = NA_integer_) {
  data.frame(
    model = model_name, status = status, reason = reason,
    n = as.integer(n), events = as.integer(events), parameters = as.integer(parameters),
    events_per_parameter = if (is.finite(parameters) && parameters > 0) events / parameters else NA_real_,
    warning_messages = warning_messages,
    max_abs_log_coefficient = max_abs_log_coefficient,
    max_standard_error = max_standard_error,
    zero_event_coefficient_levels = as.integer(zero_event_coefficient_levels),
    stringsAsFactors = FALSE
  )
}

fit_cox_with_diagnostics <- function(data, formula, model_name, output_dir, min_n = 30, min_events = 5) {
  require_packages(c("survival", "broom"))
  ensure_dir(output_dir)
  vars <- all.vars(formula)
  missing_vars <- setdiff(vars, colnames(data))
  if (length(missing_vars) > 0) stop(model_name, " is missing variables: ", paste(missing_vars, collapse = ", "))
  dat <- droplevels(data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE])
  event_audit <- cox_event_level_audit(dat, formula, model_name)
  write_csv_atomic(event_audit, file.path(output_dir, paste0(model_name, "_event_level_audit.csv")))
  n_events <- sum(dat$event == 1)
  if (nrow(dat) < min_n || n_events < min_events) {
    reason <- paste0("not fitted: requires at least ", min_n, " complete observations and ", min_events,
                     " events; observed n=", nrow(dat), ", events=", n_events)
    status <- cox_status_table(model_name, "not_estimable", reason, nrow(dat), n_events)
    write_csv_atomic(status, file.path(output_dir, paste0(model_name, "_model_status.csv")))
    return(list(fit = NULL, tidy = NULL, ph = NULL, data = dat, estimable = FALSE,
                status = status, event_audit = event_audit))
  }

  captured <- capture_coxph_fit(formula, dat)
  assessment <- assess_cox_fit(captured, event_audit)
  status <- cox_status_table(
    model_name,
    if (assessment$estimable) "estimated" else "not_estimable",
    assessment$reason, nrow(dat), n_events, assessment$parameters,
    assessment$warning_text, assessment$max_abs_log_coefficient,
    assessment$max_standard_error, assessment$zero_event_coefficient_levels
  )
  write_csv_atomic(status, file.path(output_dir, paste0(model_name, "_model_status.csv")))
  if (!assessment$estimable) {
    return(list(fit = captured$fit, tidy = NULL, ph = NULL, data = dat, estimable = FALSE,
                status = status, event_audit = event_audit))
  }

  fit <- captured$fit
  tidy <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tidy$model <- model_name
  tidy$model_status <- "estimated"
  tidy$n <- fit$n
  tidy$n_event <- fit$nevent
  tidy$concordance <- unname(summary(fit)$concordance[1])
  ph <- tryCatch(survival::cox.zph(fit, transform = "km"), error = function(e) e)
  if (inherits(ph, "error")) {
    ph_table <- data.frame(term = "NOT_ESTIMABLE", model = model_name, message = conditionMessage(ph), stringsAsFactors = FALSE)
  } else {
    ph_table <- as.data.frame(ph$table)
    ph_table$term <- rownames(ph_table)
    ph_table$model <- model_name
    rownames(ph_table) <- NULL
  }
  write_csv_atomic(tidy, file.path(output_dir, paste0(model_name, "_coefficients.csv")))
  write_csv_atomic(ph_table, file.path(output_dir, paste0(model_name, "_PH_test.csv")))
  if (!inherits(ph, "error")) {
    grDevices::pdf(file.path(output_dir, paste0(model_name, "_PH_plots.pdf")), width = 7, height = 6)
    try(plot(ph), silent = TRUE)
    grDevices::dev.off()
  }

  dfbeta_raw <- residuals(fit, type = "dfbeta")
  if (is.null(dim(dfbeta_raw))) dfbeta_raw <- matrix(dfbeta_raw, ncol = 1)
  coefficient_names <- names(stats::coef(fit))
  if (ncol(dfbeta_raw) == length(coefficient_names)) colnames(dfbeta_raw) <- coefficient_names
  sample_ids <- if ("sample" %in% colnames(dat)) as.character(dat$sample) else rownames(dat)
  if (length(sample_ids) != nrow(dfbeta_raw)) sample_ids <- rownames(dfbeta_raw)
  dfbeta <- cbind(data.frame(sample = sample_ids, stringsAsFactors = FALSE),
                  as.data.frame(dfbeta_raw, check.names = FALSE))
  write_csv_atomic(dfbeta, file.path(output_dir, paste0(model_name, "_dfbeta.csv")))
  threshold <- 2 / sqrt(nrow(dat))
  dfbeta_summary <- do.call(rbind, lapply(seq_len(ncol(dfbeta_raw)), function(j) {
    values <- abs(dfbeta_raw[, j])
    hit <- if (all(!is.finite(values))) NA_integer_ else which.max(replace(values, !is.finite(values), -Inf))
    data.frame(
      model = model_name, term = colnames(dfbeta_raw)[j],
      max_abs_dfbeta = if (is.na(hit)) NA_real_ else values[hit],
      max_abs_sample = if (is.na(hit)) NA_character_ else sample_ids[hit],
      heuristic_threshold_2_over_sqrt_n = threshold,
      exceeds_heuristic_threshold = if (is.na(hit)) NA else values[hit] > threshold,
      stringsAsFactors = FALSE
    )
  }))
  write_csv_atomic(dfbeta_summary, file.path(output_dir, paste0(model_name, "_dfbeta_summary.csv")))
  list(fit = fit, tidy = tidy, ph = ph_table, data = dat, estimable = TRUE,
       status = status, event_audit = event_audit, dfbeta_summary = dfbeta_summary)
}

record_unfitted_cox_model <- function(data, formula, model_name, output_dir, reason) {
  ensure_dir(output_dir)
  vars <- all.vars(formula)
  missing_vars <- setdiff(vars, colnames(data))
  if (length(missing_vars) > 0) stop(model_name, " is missing variables: ", paste(missing_vars, collapse = ", "))
  dat <- droplevels(data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE])
  audit <- cox_event_level_audit(dat, formula, model_name)
  status <- cox_status_table(model_name, "not_attempted", reason, nrow(dat), sum(dat$event == 1),
                             zero_event_coefficient_levels = sum(audit$model_role == "coefficient" & audit$zero_event_level))
  write_csv_atomic(audit, file.path(output_dir, paste0(model_name, "_event_level_audit.csv")))
  write_csv_atomic(status, file.path(output_dir, paste0(model_name, "_model_status.csv")))
  list(fit = NULL, tidy = NULL, ph = NULL, data = dat, estimable = FALSE,
       status = status, event_audit = audit)
}

test_score_nonlinearity <- function(data, clinical_terms, model_name, output_dir) {
  require_packages(c("survival"))
  ensure_dir(output_dir)
  linear_rhs <- paste(c("score_z", clinical_terms), collapse = " + ")
  spline_rhs <- paste(c("splines::ns(score_z, df = 3)", clinical_terms), collapse = " + ")
  f_linear <- stats::as.formula(paste("survival::Surv(time, event) ~", linear_rhs), env = parent.frame())
  f_spline <- stats::as.formula(paste("survival::Surv(time, event) ~", spline_rhs), env = parent.frame())
  vars <- unique(all.vars(f_linear))
  dat <- droplevels(data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE])
  n_events <- sum(dat$event == 1)
  if (nrow(dat) < 30 || n_events < 8) {
    out <- data.frame(model = model_name, n = nrow(dat), events = n_events,
                      linear_loglik = NA_real_, spline_loglik = NA_real_, LRT_p = NA_real_,
                      status = "not_estimable", reason = "fewer than 30 observations or 8 events",
                      warning_messages = "", stringsAsFactors = FALSE)
    write_csv_atomic(out, file.path(output_dir, paste0(model_name, "_score_nonlinearity.csv")))
    return(out)
  }
  audit <- cox_event_level_audit(dat, f_linear, model_name)
  linear_capture <- capture_coxph_fit(f_linear, dat)
  spline_capture <- capture_coxph_fit(f_spline, dat)
  linear_assessment <- assess_cox_fit(linear_capture, audit)
  spline_assessment <- assess_cox_fit(spline_capture, audit)
  estimable <- linear_assessment$estimable && spline_assessment$estimable
  warnings <- collapse_messages(c(linear_capture$warnings, spline_capture$warnings))
  if (!estimable) {
    reason <- collapse_messages(c(linear_assessment$reason, spline_assessment$reason))
    out <- data.frame(model = model_name, n = nrow(dat), events = n_events,
                      linear_loglik = NA_real_, spline_loglik = NA_real_, LRT_p = NA_real_,
                      status = "not_estimable", reason = reason,
                      warning_messages = warnings, stringsAsFactors = FALSE)
  } else {
    linear <- linear_capture$fit
    spline <- spline_capture$fit
    lrt <- stats::anova(linear, spline, test = "LRT")
    p_col <- grep("Pr|P\\(", colnames(lrt), value = TRUE)[1]
    lrt_p <- if (is.na(p_col)) NA_real_ else as.numeric(lrt[2, p_col])
    out <- data.frame(model = model_name, n = nrow(dat), events = n_events,
                      linear_loglik = as.numeric(stats::logLik(linear)),
                      spline_loglik = as.numeric(stats::logLik(spline)), LRT_p = lrt_p,
                      status = "estimated", reason = "linear and spline models passed estimability checks",
                      warning_messages = warnings, stringsAsFactors = FALSE)
  }
  write_csv_atomic(out, file.path(output_dir, paste0(model_name, "_score_nonlinearity.csv")))
  out
}

compare_incremental_cox <- function(data, clinical_formula, score_formula, cohort_split, analysis_type) {
  vars <- unique(c(all.vars(clinical_formula), all.vars(score_formula)))
  dat <- droplevels(data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE])
  n_events <- sum(dat$event == 1)
  clinical_capture <- capture_coxph_fit(clinical_formula, dat)
  score_capture <- capture_coxph_fit(score_formula, dat)
  clinical_audit <- cox_event_level_audit(dat, clinical_formula, paste0(cohort_split, "_", analysis_type, "_clinical"))
  score_audit <- cox_event_level_audit(dat, score_formula, paste0(cohort_split, "_", analysis_type, "_clinical_plus_score"))
  clinical_assessment <- assess_cox_fit(clinical_capture, clinical_audit)
  score_assessment <- assess_cox_fit(score_capture, score_audit)
  estimable <- clinical_assessment$estimable && score_assessment$estimable
  base <- data.frame(
    cohort_split = cohort_split, analysis_type = analysis_type,
    status = if (estimable) "estimated" else "not_estimable",
    reason = if (estimable) "both nested models passed estimability checks" else
      collapse_messages(c(clinical_assessment$reason, score_assessment$reason)),
    n = nrow(dat), events = n_events, stringsAsFactors = FALSE
  )
  if (!estimable) {
    return(cbind(base, data.frame(
      clinical_C_index = NA_real_, clinical_plus_score_C_index = NA_real_,
      delta_C_index_apparent = NA_real_, LRT_chisq = NA_real_, LRT_df = NA_real_, LRT_p = NA_real_
    )))
  }
  clinical_fit <- clinical_capture$fit
  score_fit <- score_capture$fit
  lrt <- stats::anova(clinical_fit, score_fit, test = "LRT")
  p_col <- grep("Pr|P\\(", colnames(lrt), value = TRUE)[1]
  clinical_c <- unname(summary(clinical_fit)$concordance[1])
  score_c <- unname(summary(score_fit)$concordance[1])
  cbind(base, data.frame(
    clinical_C_index = clinical_c,
    clinical_plus_score_C_index = score_c,
    delta_C_index_apparent = score_c - clinical_c,
    LRT_chisq = as.numeric(lrt[2, "Chisq"]), LRT_df = as.numeric(lrt[2, "Df"]),
    LRT_p = if (is.na(p_col)) NA_real_ else as.numeric(lrt[2, p_col])
  ))
}
