#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

task11_normalize_ensembl <- function(x) {
  out <- sub("[.][0-9]+$", "", as.character(x))
  if (any(!grepl("^ENSG[0-9]+$", out))) stop("Invalid normalized Ensembl ID", call. = FALSE)
  out
}

task11_prepare_rnaseq_background <- function(unique_ensembl_counts, design) {
  unique_ensembl_counts <- as.matrix(unique_ensembl_counts)
  if (is.null(rownames(unique_ensembl_counts)) || anyDuplicated(task11_normalize_ensembl(rownames(unique_ensembl_counts)))) {
    stop("Counts must already be aggregated to unique normalized Ensembl rows", call. = FALSE)
  }
  if (ncol(unique_ensembl_counts) != nrow(design)) stop("RNA-seq design does not match counts", call. = FALSE)
  dge <- edgeR::DGEList(counts = unique_ensembl_counts)
  keep <- edgeR::filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logcpm <- edgeR::cpm(
    dge, log = TRUE, prior.count = 0.25,
    normalized.lib.sizes = TRUE
  )
  list(dge = dge, mean_logcpm = rowMeans(logcpm), background_ensembl = rownames(dge))
}

task11_build_alias_map <- function(targets) {
  required <- c("entrez_id", "ensembl_gene_ids")
  if (!all(required %in% names(targets))) stop("Target alias input is incomplete", call. = FALSE)
  rows <- lapply(seq_len(nrow(targets)), function(i) {
    aliases <- strsplit(as.character(targets$ensembl_gene_ids[[i]]), ";", fixed = TRUE)[[1L]]
    aliases <- aliases[nzchar(aliases)]
    if (!length(aliases)) return(NULL)
    data.frame(
      entrez_id = rep(as.character(targets$entrez_id[[i]]), length(aliases)),
      ensembl_id = task11_normalize_ensembl(aliases),
      stringsAsFactors = FALSE
    )
  })
  alias_map <- unique(do.call(rbind, rows))
  n_entrez <- vapply(
    split(alias_map$entrez_id, alias_map$ensembl_id),
    function(x) length(unique(x)), integer(1)
  )
  ambiguous <- sort(names(n_entrez)[n_entrez > 1L], method = "radix")
  list(
    unambiguous = alias_map[!alias_map$ensembl_id %in% ambiguous, , drop = FALSE],
    ambiguous_ensembl = ambiguous
  )
}

task11_select_rnaseq_target_units <- function(filtered_ensembl, mean_logcpm, targets) {
  filtered_ensembl <- task11_normalize_ensembl(filtered_ensembl)
  if (anyDuplicated(filtered_ensembl)) stop("Filtered background must contain unique Ensembl rows", call. = FALSE)
  if (length(mean_logcpm) != length(filtered_ensembl) || any(!is.finite(mean_logcpm))) {
    stop("Blind mean logCPM does not match the filtered background", call. = FALSE)
  }
  aliases <- task11_build_alias_map(targets)
  candidates <- aliases$unambiguous[
    aliases$unambiguous$ensembl_id %in% filtered_ensembl, , drop = FALSE
  ]
  candidates$mean_logcpm <- mean_logcpm[match(candidates$ensembl_id, filtered_ensembl)]
  candidates <- candidates[order(
    candidates$entrez_id, -candidates$mean_logcpm, candidates$ensembl_id,
    method = "radix"
  ), , drop = FALSE]
  selected <- candidates[!duplicated(candidates$entrez_id), , drop = FALSE]
  rownames(selected) <- NULL
  if (anyDuplicated(selected$entrez_id) || anyDuplicated(selected$ensembl_id)) {
    stop("RNA-seq target selection is not one Entrez to one row", call. = FALSE)
  }
  list(selection = selected, ambiguous_ensembl = aliases$ambiguous_ensembl)
}

task11_target_index <- function(selection, target_entrez, background_ensembl) {
  target_entrez <- unique(as.character(target_entrez))
  selected <- selection[selection$entrez_id %in% target_entrez, , drop = FALSE]
  index <- match(selected$ensembl_id, background_ensembl)
  if (anyNA(index) || anyDuplicated(selected$entrez_id) || anyDuplicated(index)) {
    stop("Target coverage/index identity failed", call. = FALSE)
  }
  if (nrow(selected) != length(index)) stop("Coverage must equal camera index length", call. = FALSE)
  list(index = unname(index), coverage_unique_entrez = nrow(selected), selected = selected)
}

task11_camera <- function(y, index, design, contrast) {
  if (ncol(y) != nrow(design)) stop("Full design does not match the expression object", call. = FALSE)
  if (length(contrast) != ncol(design)) stop("Contrast does not match the full design", call. = FALSE)
  limma::camera(
    y = y, index = index, design = design, contrast = contrast,
    directional = TRUE, inter.gene.cor = NA_real_
  )
}

task11_signed_stouffer <- function(p_value, direction, n_independent) {
  if (length(p_value) < 2L || length(direction) != length(p_value) || length(n_independent) != length(p_value)) {
    stop("Stouffer inputs have incompatible lengths", call. = FALSE)
  }
  if (any(!is.finite(p_value)) || any(p_value <= 0 | p_value > 1)) stop("P values must be in (0,1]", call. = FALSE)
  if (any(!direction %in% c("Up", "Down"))) stop("Direction must be Up or Down", call. = FALSE)
  if (any(!is.finite(n_independent)) || any(n_independent <= 0)) stop("Independent sample counts must be positive", call. = FALSE)
  z_unsigned <- stats::qnorm(p_value / 2, lower.tail = FALSE)
  z_signed <- ifelse(direction == "Up", z_unsigned, -z_unsigned)
  weight <- sqrt(n_independent)
  z_combined <- sum(weight * z_signed) / sqrt(sum(weight^2))
  p_two_sided <- 2 * stats::pnorm(abs(z_combined), lower.tail = FALSE)
  list(z_i = z_signed, weight = weight, z_combined = z_combined, p_two_sided = p_two_sided)
}

task11_top_variable_gene_ids <- function(expression, n_top = 5000L) {
  if (is.null(rownames(expression)) || anyDuplicated(rownames(expression))) {
    stop("Gene IDs must be present and unique", call. = FALSE)
  }
  gene_variance <- apply(expression, 1L, stats::var)
  if (any(!is.finite(gene_variance))) stop("Gene variance is non-finite", call. = FALSE)
  ordered <- order(-gene_variance, rownames(expression), method = "radix")
  rownames(expression)[ordered[seq_len(min(as.integer(n_top), length(ordered)))]]
}

task11_load_task5_qc_gate <- function() {
  if (exists("apply_two_metric_qc_gate", mode = "function", inherits = TRUE)) return(invisible(TRUE))
  old <- Sys.getenv("TASK5_FUNCTIONS_ONLY", unset = NA_character_)
  Sys.setenv(TASK5_FUNCTIONS_ONLY = "true")
  on.exit(if (is.na(old)) Sys.unsetenv("TASK5_FUNCTIONS_ONLY") else Sys.setenv(TASK5_FUNCTIONS_ONLY = old), add = TRUE)
  sys.source(file.path("analysis", "R", "04_prepare_microarray.R"), envir = .GlobalEnv)
  invisible(TRUE)
}

task11_gse188427_blind_qc <- function(expression) {
  task11_load_task5_qc_gate()
  expression <- as.matrix(expression)
  storage.mode(expression) <- "double"
  if (is.null(colnames(expression)) || anyDuplicated(colnames(expression)) || any(!is.finite(expression))) {
    stop("Blind expression input is invalid", call. = FALSE)
  }
  top_gene_ids <- task11_top_variable_gene_ids(expression, 5000L)
  top_expression <- expression[top_gene_ids, , drop = FALSE]
  correlations <- stats::cor(top_expression, method = "spearman")
  diag(correlations) <- NA_real_
  pca <- stats::prcomp(t(top_expression), center = TRUE, scale. = FALSE)
  k <- min(5L, ncol(pca$x))
  if (any(!is.finite(pca$sdev[seq_len(k)]) | pca$sdev[seq_len(k)] <= 0)) stop("Blind PCA scale is invalid", call. = FALSE)
  standardized <- sweep(pca$x[, seq_len(k), drop = FALSE], 2L, pca$sdev[seq_len(k)], "/")
  pca_distance <- sqrt(rowSums(standardized^2))
  intercept_only <- matrix(1, nrow = ncol(expression), ncol = 1L)
  array_weight <- limma::arrayWeights(expression, design = intercept_only)
  metrics <- data.frame(
    sample_id = colnames(expression),
    missing_rate = colMeans(!is.finite(expression)),
    median_expression = apply(expression, 2L, stats::median),
    expression_iqr = apply(expression, 2L, stats::IQR),
    median_sample_correlation = apply(correlations, 2L, stats::median, na.rm = TRUE),
    pca_distance_5pc = unname(pca_distance[colnames(expression)]),
    log_array_weight = log(as.numeric(array_weight)),
    stringsAsFactors = FALSE
  )
  gated <- apply_two_metric_qc_gate(metrics, setdiff(names(metrics), "sample_id"))
  forbidden <- c("case_status", "status", "severity", "arm", "pathogen")
  if (any(names(gated) %in% forbidden)) stop("Phenotype leaked into blind QC", call. = FALSE)
  list(qc = gated, top_gene_ids = top_gene_ids, pca_scores = pca$x[, seq_len(k), drop = FALSE])
}
