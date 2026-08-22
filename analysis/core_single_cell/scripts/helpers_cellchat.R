# Shared CellChat helpers. Source config.R before this file.

run_cellchat_condition <- function(seu, condition_label, min_cells = 20) {
  require_packages(c("Seurat", "CellChat", "Matrix", "future"))
  Seurat::DefaultAssay(seu) <- "RNA"
  expression <- get_assay_data_compat(seu, "RNA", "data")
  metadata <- seu@meta.data[colnames(expression), , drop = FALSE]
  metadata$comm_group <- droplevels(factor(metadata$comm_group))
  if (!"Patient" %in% colnames(metadata)) {
    stop("Patient metadata is required for CellChat samples.")
  }
  patient_id <- trimws(as.character(metadata$Patient))
  if (anyNA(patient_id) || any(!nzchar(patient_id))) {
    stop("Patient metadata contains missing or blank values.")
  }
  metadata$samples <- factor(patient_id)
  object <- CellChat::createCellChat(expression, meta = metadata, group.by = "comm_group")
  object@DB <- CellChat::CellChatDB.human
  object <- CellChat::subsetData(object)
  object <- CellChat::identifyOverExpressedGenes(object)
  object <- CellChat::identifyOverExpressedInteractions(object)
  object <- CellChat::computeCommunProb(object, type = "triMean", population.size = FALSE, raw.use = TRUE)
  object <- CellChat::filterCommunication(object, min.cells = min_cells)
  object <- CellChat::computeCommunProbPathway(object)
  object <- CellChat::aggregateNet(object)
  attr(object, "condition_label") <- condition_label
  object
}

extract_cellchat_raw_probabilities <- function(object) {
  prob <- object@net$prob
  pval <- object@net$pval
  if (length(dim(prob)) != 3) stop("CellChat net$prob is not a source × target × interaction array.")
  dims <- dimnames(prob)
  grid <- expand.grid(
    source = dims[[1]], target = dims[[2]], interaction_name = dims[[3]],
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid$prob <- as.numeric(prob)
  grid$pval <- as.numeric(pval)
  lr <- as.data.frame(object@LR$LRsig, stringsAsFactors = FALSE)
  if (!"interaction_name" %in% colnames(lr)) lr$interaction_name <- rownames(lr)
  keep <- intersect(c("interaction_name", "ligand", "receptor", "pathway_name", "annotation"), colnames(lr))
  dplyr::left_join(grid, unique(lr[, keep, drop = FALSE]), by = "interaction_name")
}

average_expression_by_patient_group <- function(seu) {
  require_packages(c("Matrix"))
  expression <- get_assay_data_compat(seu, "RNA", "data")
  metadata <- seu@meta.data[colnames(expression), , drop = FALSE]
  key <- paste(metadata$Patient, metadata$condition, metadata$comm_group, sep = "||")
  split_columns <- split(seq_len(ncol(expression)), key)
  average <- vapply(
    split_columns,
    function(index) as.numeric(Matrix::rowMeans(expression[, index, drop = FALSE])),
    numeric(nrow(expression))
  )
  rownames(average) <- rownames(expression)
  average
}

complex_expression <- function(symbol, average_matrix, column_name) {
  if (is.na(symbol) || !nzchar(symbol) || !column_name %in% colnames(average_matrix)) return(NA_real_)
  subunits <- unique(strsplit(as.character(symbol), "_", fixed = TRUE)[[1]])
  if (!all(subunits %in% rownames(average_matrix))) return(NA_real_)
  values <- average_matrix[subunits, column_name, drop = TRUE]
  if (any(!is.finite(values))) return(NA_real_)
  min(values)
}

paired_lr_expression_tests <- function(seu, lr_table, source_groups, paired_patients) {
  average <- average_expression_by_patient_group(seu)
  lr_table <- unique(lr_table[, intersect(c("interaction_name", "ligand", "receptor", "pathway_name"), colnames(lr_table)), drop = FALSE])
  lr_table <- lr_table[!is.na(lr_table$ligand) & !is.na(lr_table$receptor), , drop = FALSE]
  results <- vector("list", nrow(lr_table) * length(source_groups))
  counter <- 0L

  for (source in source_groups) {
    for (i in seq_len(nrow(lr_table))) {
      row <- lr_table[i, , drop = FALSE]
      tumour <- normal <- rep(NA_real_, length(paired_patients))
      for (j in seq_along(paired_patients)) {
        patient <- paired_patients[j]
        tumour_source <- paste(patient, "Tumor", source, sep = "||")
        tumour_target <- paste(patient, "Tumor", "Epithelial", sep = "||")
        normal_source <- paste(patient, "Normal", source, sep = "||")
        normal_target <- paste(patient, "Normal", "Epithelial", sep = "||")
        lt <- complex_expression(row$ligand, average, tumour_source)
        rt <- complex_expression(row$receptor, average, tumour_target)
        ln <- complex_expression(row$ligand, average, normal_source)
        rn <- complex_expression(row$receptor, average, normal_target)
        if (all(is.finite(c(lt, rt)))) tumour[j] <- sqrt(pmax(lt, 0) * pmax(rt, 0))
        if (all(is.finite(c(ln, rn)))) normal[j] <- sqrt(pmax(ln, 0) * pmax(rn, 0))
      }
      complete <- is.finite(tumour) & is.finite(normal)
      n_pairs <- sum(complete)
      p <- NA_real_
      if (n_pairs >= 5 && any(abs(tumour[complete] - normal[complete]) > 0)) {
        p <- tryCatch(
          stats::wilcox.test(tumour[complete], normal[complete], paired = TRUE, exact = FALSE)$p.value,
          error = function(e) NA_real_
        )
      }
      counter <- counter + 1L
      results[[counter]] <- data.frame(
        source = source, target = "Epithelial", interaction_name = row$interaction_name,
        ligand = row$ligand, receptor = row$receptor,
        pathway_name = if ("pathway_name" %in% colnames(row)) row$pathway_name else NA_character_,
        n_pairs = n_pairs, median_tumor = if (n_pairs) median(tumour[complete]) else NA_real_,
        median_normal = if (n_pairs) median(normal[complete]) else NA_real_,
        median_difference = if (n_pairs) median(tumour[complete] - normal[complete]) else NA_real_,
        p_value = p, stringsAsFactors = FALSE
      )
    }
  }
  out <- dplyr::bind_rows(results[seq_len(counter)])
  out$FDR <- p.adjust(out$p_value, method = "BH")
  out
}
