#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

forbidden_qc_fields <- c(
  "case_status", "condition", "status", "severity", "pathogen", "coinfection",
  "treatment", "analysis_role", "inclusion_scope", "outcome_previously_observed"
)

make_blind_sample_data <- function(x) {
  required <- c("sample_id", "platform")
  if (!all(required %in% names(x))) stop("Blind sample input lacks sample_id/platform", call. = FALSE)
  technical_batch <- if ("technical_batch" %in% names(x)) x$technical_batch else if ("batch" %in% names(x)) x$batch else NA_character_
  out <- data.frame(
    sample_id = as.character(x$sample_id),
    platform = as.character(x$platform),
    technical_batch = as.character(technical_batch),
    stringsAsFactors = FALSE
  )
  if (any(names(out) %in% forbidden_qc_fields)) stop("Phenotype field leaked into blind QC object", call. = FALSE)
  if (anyDuplicated(out$sample_id)) stop("Blind QC sample IDs are not unique", call. = FALSE)
  out
}

transform_processed_expression <- function(expression) {
  expression <- as.matrix(expression)
  storage.mode(expression) <- "double"
  if (!length(expression) || any(!is.finite(expression))) {
    stop("Expression matrix contains non-finite values", call. = FALSE)
  }
  finite_values <- as.numeric(expression)
  intensity_scale <- unname(stats::quantile(finite_values, 0.99, na.rm = TRUE)) > 100
  if (intensity_scale) {
    transformed <- log2(pmax(expression, 1))
    action <- "log2_pmax1_only"
  } else {
    transformed <- expression
    action <- "already_log2_no_transform"
  }
  list(
    expression = transformed,
    action = action,
    between_array_normalization_applied = FALSE,
    raw_min = min(finite_values), raw_max = max(finite_values),
    raw_q99 = unname(stats::quantile(finite_values, 0.99, na.rm = TRUE))
  )
}

remove_expression_component <- function(x) {
  x$expression <- NULL
  x
}

strip_trailing_whitespace <- function(x) sub("[[:space:]]+$", "", x)

robust_flag <- function(x) {
  if (!length(x) || any(!is.finite(x))) stop("Blind QC metric contains non-finite values", call. = FALSE)
  center <- stats::median(x)
  spread <- stats::mad(x, center = center, constant = 1.4826)
  if (!is.finite(center) || !is.finite(spread)) stop("Blind QC center/scale is non-finite", call. = FALSE)
  if (spread == 0) return(abs(x - center) > 0)
  abs(x - center) > 5 * spread
}

apply_two_metric_qc_gate <- function(metrics, metric_names) {
  if (!all(c("sample_id", metric_names) %in% names(metrics))) stop("QC metrics are incomplete", call. = FALSE)
  if (any(!vapply(metrics[metric_names], function(x) all(is.finite(x)), logical(1)))) {
    stop("QC gate input contains non-finite values", call. = FALSE)
  }
  flags <- lapply(metric_names, function(name) robust_flag(metrics[[name]]))
  names(flags) <- paste0(metric_names, "_flag")
  flags <- as.data.frame(flags, check.names = FALSE)
  out <- cbind(metrics, flags)
  out$n_blind_flags <- rowSums(flags, na.rm = TRUE)
  out$exclude_qc <- out$n_blind_flags >= 2L
  out
}

collapse_probes_blind <- function(expression, probe_map) {
  expression <- as.matrix(expression)
  if (is.null(rownames(expression))) stop("Probe expression lacks row names", call. = FALSE)
  required <- c("probe_id", "entrez_id", "official_symbol", "mapping_status")
  if (!all(required %in% names(probe_map))) stop("Probe map is incomplete", call. = FALSE)
  map <- probe_map[match(rownames(expression), probe_map$probe_id), , drop = FALSE]
  keep <- !is.na(map$probe_id) & map$mapping_status == "unique_current" &
    grepl("^[0-9]+$", map$entrez_id)
  map <- map[keep, , drop = FALSE]
  expr <- expression[keep, , drop = FALSE]
  map$all_sample_mean <- rowMeans(expr, na.rm = TRUE)
  ord <- order(map$entrez_id, -map$all_sample_mean, map$probe_id, method = "radix")
  map <- map[ord, , drop = FALSE]
  expr <- expr[ord, , drop = FALSE]
  selected <- !duplicated(map$entrez_id)
  map <- map[selected, , drop = FALSE]
  expr <- expr[selected, , drop = FALSE]
  rownames(expr) <- map$entrez_id
  rownames(map) <- NULL
  list(expression = expr, selected_probes = map)
}

read_series_matrix_blind <- function(path) {
  if (!grepl("series_matrix[.]txt[.]gz$", path)) stop("Only processed series matrix is allowed", call. = FALSE)
  con <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  allowed <- list()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("series_matrix_table_begin not found", call. = FALSE)
    if (identical(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_geo_accession\t")) allowed$sample_id <- line
    if (startsWith(line, "!Sample_description\t")) allowed$array_id <- line
    if (startsWith(line, "!Sample_platform_id\t")) allowed$platform <- line
    if (startsWith(line, "!Sample_data_processing\t")) allowed$data_processing <- line
  }
  parse_values <- function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
    gsub('^"|"$', "", fields)
  }
  if (!all(c("sample_id", "array_id", "platform", "data_processing") %in% names(allowed))) {
    stop("Required non-phenotype series metadata is absent", call. = FALSE)
  }
  sample_id <- parse_values(allowed$sample_id)
  array_id <- parse_values(allowed$array_id)
  platform <- parse_values(allowed$platform)
  processing <- parse_values(allowed$data_processing)
  if (!all(lengths(list(array_id, platform, processing)) == length(sample_id))) {
    stop("Allowed series metadata widths differ", call. = FALSE)
  }
  table <- read.delim(
    con, sep = "\t", header = TRUE, quote = "\"", comment.char = "!",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!identical(names(table)[-1L], sample_id)) stop("Series expression columns do not match allowed sample IDs", call. = FALSE)
  probe_id <- as.character(table[[1L]])
  expression <- as.matrix(table[-1L])
  storage.mode(expression) <- "double"
  rownames(expression) <- probe_id
  list(
    expression = expression,
    sample_data = data.frame(sample_id = sample_id, array_id = array_id, platform = platform, stringsAsFactors = FALSE),
    data_processing = unique(processing)
  )
}

read_detection_pvalues <- function(path) {
  if (!grepl("non-normalized.*[.]txt[.]gz$|non-normalized_data[.]txt[.]gz$", basename(path), ignore.case = TRUE)) {
    stop("Detection P values must come from the processed non-normalized table", call. = FALSE)
  }
  tab <- read.delim(path, sep = "\t", header = TRUE, quote = "", check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(tab) < 3L || (ncol(tab) - 1L) %% 2L != 0L) stop("Unexpected non-normalized column structure", call. = FALSE)
  value_columns <- seq.int(2L, ncol(tab), by = 2L)
  p_columns <- value_columns + 1L
  if (!all(names(tab)[p_columns] == "Detection Pval")) stop("Detection Pval columns do not alternate with intensity columns", call. = FALSE)
  p <- as.matrix(tab[p_columns])
  storage.mode(p) <- "double"
  rownames(p) <- as.character(tab[[1L]])
  colnames(p) <- names(tab)[value_columns]
  if (any(!is.finite(p))) stop("Detection P table contains non-finite values", call. = FALSE)
  if (any(p < 0 | p > 1)) stop("Detection P values outside [0,1]", call. = FALSE)
  p
}

read_gpl10558_annotation <- function(path) {
  con <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("GPL10558 annotation table header not found", call. = FALSE)
    if (startsWith(line, "ID\t")) break
  }
  headers <- strsplit(line, "\t", fixed = TRUE)[[1L]]
  tab <- read.delim(
    con, sep = "\t", header = FALSE, col.names = make.names(headers, unique = TRUE),
    quote = "", comment.char = "!", fill = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(tab) <- headers
  tab <- tab[nzchar(tab$ID), c("ID", "Gene symbol", "Gene ID"), drop = FALSE]
  names(tab) <- c("probe_id", "source_gene_symbol", "source_entrez_id")
  tab
}

build_current_probe_map <- function(annotation) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Frozen AnnotationDbi and org.Hs.eg.db are required", call. = FALSE)
  }
  source_id <- trimws(annotation$source_entrez_id)
  syntactic_unique <- grepl("^[0-9]+$", source_id)
  keys <- unique(source_id[syntactic_unique])
  official <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys = keys, keytype = "ENTREZID",
    column = "SYMBOL", multiVals = "first"
  )
  official_symbol <- unname(official[source_id])
  status <- ifelse(
    !nzchar(source_id), "empty",
    ifelse(!syntactic_unique, "multi_entrez", ifelse(is.na(official_symbol) | !nzchar(official_symbol), "obsolete_or_unresolved", "unique_current"))
  )
  data.frame(
    probe_id = annotation$probe_id,
    source_gene_symbol = annotation$source_gene_symbol,
    source_entrez_id = source_id,
    entrez_id = ifelse(status == "unique_current", source_id, NA_character_),
    official_symbol = ifelse(status == "unique_current", official_symbol, NA_character_),
    mapping_status = status,
    symbol_conflict = status == "unique_current" & nzchar(annotation$source_gene_symbol) & annotation$source_gene_symbol != official_symbol,
    stringsAsFactors = FALSE
  )
}

calculate_blind_qc <- function(expression, detection_p) {
  if (any(!is.finite(expression)) || any(!is.finite(detection_p))) {
    stop("Non-finite base data cannot enter blind QC", call. = FALSE)
  }
  variable <- apply(expression, 1L, stats::var, na.rm = TRUE)
  variable[!is.finite(variable)] <- -Inf
  top <- order(variable, decreasing = TRUE)[seq_len(min(5000L, sum(is.finite(variable))))]
  top_expr <- expression[top, , drop = FALSE]
  correlations <- stats::cor(top_expr, use = "pairwise.complete.obs", method = "spearman")
  diag(correlations) <- NA_real_
  pca_fit <- stats::prcomp(t(top_expr), center = TRUE, scale. = FALSE)
  k <- min(5L, ncol(pca_fit$x))
  standardized_scores <- sweep(pca_fit$x[, seq_len(k), drop = FALSE], 2L, pca_fit$sdev[seq_len(k)], "/")
  pca_distance <- sqrt(rowSums(standardized_scores^2))
  if (!requireNamespace("limma", quietly = TRUE)) stop("Frozen limma is required", call. = FALSE)
  weights <- limma::arrayWeights(expression, design = matrix(1, ncol(expression), 1L))
  out <- data.frame(
    sample_id = colnames(expression),
    missing_rate = colMeans(!is.finite(expression)),
    median_expression = apply(expression, 2L, stats::median, na.rm = TRUE),
    expression_iqr = apply(expression, 2L, stats::IQR, na.rm = TRUE),
    median_sample_correlation = apply(correlations, 2L, stats::median, na.rm = TRUE),
    pca_distance_5pc = unname(pca_distance[colnames(expression)]),
    log_array_weight = log(as.numeric(weights)),
    detection_rate_p_lt_0_05 = colMeans(detection_p < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  if (any(!vapply(out[-1L], function(x) all(is.finite(x)), logical(1)))) {
    stop("Derived blind QC metric contains non-finite values", call. = FALSE)
  }
  out
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 failed for ", path, call. = FALSE)
  sub("[[:space:]].*$", "", out[[1L]])
}

write_qc_figure <- function(cohort, expression, qc, pca_scores, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = 10, height = 8, useDingbats = FALSE, timestamp = FALSE, bg = "white")
  old <- graphics::par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  on.exit({
    graphics::par(old)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::boxplot(as.data.frame(expression), outline = FALSE, xaxt = "n", main = paste(cohort, "pre-QC log2 expression"), ylab = "log2 intensity")
  graphics::hist(qc$median_sample_correlation, breaks = 20, main = "Median sample correlation", xlab = "Spearman correlation")
  graphics::hist(qc$detection_rate_p_lt_0_05, breaks = 20, main = "Detection rate", xlab = "fraction P < 0.05")
  color <- ifelse(qc$exclude_qc, "#D55E00", "#0072B2")
  graphics::plot(pca_scores[, 1], pca_scores[, 2], col = color, pch = 19, main = "Blind PCA", xlab = "PC1", ylab = "PC2")
  graphics::plot(qc$log_array_weight, qc$median_sample_correlation, col = color, pch = 19, main = "Blind array quality", xlab = "log array weight", ylab = "median correlation")
  graphics::barplot(table(qc$n_blind_flags), main = "Blind flags per sample", xlab = "number of flags", ylab = "samples")
}

run_cohort <- function(cohort, series_path, detection_path, annotation_map, blind_samples, targets) {
  series <- read_series_matrix_blind(series_path)
  if (length(series$data_processing) != 1L || !grepl("background.*(scale|normal)", series$data_processing, ignore.case = TRUE)) {
    stop(cohort, ": processed-data statement does not confirm background subtraction and scaling", call. = FALSE)
  }
  transformed <- transform_processed_expression(series$expression)
  expression <- transformed$expression
  series <- remove_expression_component(series)
  sample_index <- match(blind_samples$sample_id, colnames(expression))
  if (anyNA(sample_index)) stop(cohort, ": frozen samples absent from series matrix", call. = FALSE)
  expression <- expression[, sample_index, drop = FALSE]
  series_sample <- series$sample_data[match(colnames(expression), series$sample_data$sample_id), , drop = FALSE]
  if (!identical(series_sample$sample_id, blind_samples$sample_id)) stop(cohort, ": blind sample order mismatch", call. = FALSE)
  detection_by_array <- read_detection_pvalues(detection_path)
  array_index <- match(series_sample$array_id, colnames(detection_by_array))
  if (anyNA(array_index)) stop(cohort, ": array IDs cannot be matched to processed Detection P values", call. = FALSE)
  probe_index <- match(rownames(expression), rownames(detection_by_array))
  if (anyNA(probe_index)) stop(cohort, ": expression probes absent from Detection P table", call. = FALSE)
  detection <- detection_by_array[probe_index, array_index, drop = FALSE]
  rownames(detection) <- rownames(expression)
  colnames(detection) <- colnames(expression)
  rm(detection_by_array)
  if (any(!is.finite(expression)) || any(!is.finite(detection))) {
    stop(cohort, ": non-finite expression or Detection P value", call. = FALSE)
  }

  qc_base <- calculate_blind_qc(expression, detection)
  metric_names <- setdiff(names(qc_base), "sample_id")
  qc <- apply_two_metric_qc_gate(qc_base, metric_names)
  qc$mad_constant <- 1.4826
  qc$mad_rule <- "two_sided_abs_x_minus_median_gt_5_scaled_MAD"
  qc$mad_zero_rule <- "all_equal_no_flag_otherwise_any_nonzero_deviation_flags"
  qc$nonfinite_policy <- "stop_entire_cohort"
  qc$qc_decision <- ifelse(qc$exclude_qc, "exclude_two_or_more_blind_scaled_5MAD_flags", "retain")
  qc$qc_boundary <- "processed-only: no CEL image, RNA degradation, or raw background reprocessing claims"
  keep_samples <- !qc$exclude_qc

  variable <- apply(expression, 1L, stats::var, na.rm = TRUE)
  top <- order(variable, decreasing = TRUE)[seq_len(min(5000L, length(variable)))]
  pca_scores <- stats::prcomp(t(expression[top, , drop = FALSE]), center = TRUE, scale. = FALSE)$x[, 1:2, drop = FALSE]

  cohort_map <- annotation_map[match(rownames(expression), annotation_map$probe_id), , drop = FALSE]
  absent <- is.na(cohort_map$probe_id)
  if (any(absent)) {
    cohort_map[absent, ] <- data.frame(
      probe_id = rownames(expression)[absent], source_gene_symbol = NA_character_, source_entrez_id = NA_character_,
      entrez_id = NA_character_, official_symbol = NA_character_, mapping_status = "annotation_absent",
      symbol_conflict = FALSE, stringsAsFactors = FALSE
    )
  }
  cohort_map$probe_id <- rownames(expression)
  cohort_map$all_sample_mean_log2 <- rowMeans(expression, na.rm = TRUE)
  cohort_map$detection_rate_p_lt_0_05 <- rowMeans(detection < 0.05, na.rm = TRUE)
  collapsed <- collapse_probes_blind(expression[, keep_samples, drop = FALSE], cohort_map)
  selected <- collapsed$selected_probes
  selected_detection <- detection[match(selected$probe_id, rownames(detection)), keep_samples, drop = FALSE]
  detection_rate <- rowMeans(selected_detection < 0.05, na.rm = TRUE)
  feature_data <- data.frame(
    entrez_id = selected$entrez_id,
    official_symbol = selected$official_symbol,
    selected_probe_id = selected$probe_id,
    all_sample_mean_log2 = selected$all_sample_mean,
    detection_rate_p_lt_0_05 = detection_rate,
    detectable_blind = detection_rate >= 0.10,
    stringsAsFactors = FALSE
  )
  rownames(feature_data) <- feature_data$entrez_id
  detectable_expression <- collapsed$expression[feature_data$detectable_blind, , drop = FALSE]
  detectable_features <- feature_data[feature_data$detectable_blind, , drop = FALSE]
  target_detectable <- targets$entrez_id %in% detectable_features$entrez_id
  target_detail <- data.frame(
    entrez_id = targets$entrez_id,
    gene_symbol = targets$gene_symbol,
    detectable = target_detectable,
    selected_probe_id = detectable_features$selected_probe_id[match(targets$entrez_id, detectable_features$entrez_id)],
    detection_rate_p_lt_0_05 = detectable_features$detection_rate_p_lt_0_05[match(targets$entrez_id, detectable_features$entrez_id)],
    stringsAsFactors = FALSE
  )
  coverage_n <- sum(target_detectable)
  coverage_fraction <- coverage_n / nrow(targets)
  target_detail$cohort <- cohort
  target_detail$n_primary_targets <- nrow(targets)
  target_detail$n_detectable_targets <- coverage_n
  target_detail$coverage_fraction <- coverage_fraction
  target_detail$interpretability <- if (coverage_n < 10L || coverage_fraction < 0.50) "uninterpretable_by_preregistered_rule" else "interpretable"

  derived <- list(
    cohort = cohort,
    expression = detectable_expression,
    sample_data = blind_samples[keep_samples, , drop = FALSE],
    feature_data = detectable_features,
    qc_excluded_sample_ids = qc$sample_id[qc$exclude_qc],
    preprocessing = list(
      source = "GEO series matrix processed intensities",
      author_processing = series$data_processing,
      transform = transformed$action,
      between_array_normalization_applied = FALSE,
      detection_source = basename(detection_path),
      detection_rule = "Detection P < 0.05 in at least 10% of retained samples",
      qc_mad_constant = 1.4826,
      qc_mad_rule = "two_sided_abs_x_minus_median_gt_5_scaled_MAD",
      qc_nonfinite_policy = "stop_entire_cohort",
      phenotype_fields_available_in_qc_object = FALSE,
      processed_only_boundary = "No RAW/CEL/FASTQ/SRA read; no CEL image, RNA degradation, or raw background claims"
    )
  )
  if (any(names(derived$sample_data) %in% forbidden_qc_fields)) stop("Phenotype leakage in derived QC object", call. = FALSE)
  list(
    derived = derived, qc = qc, pca_scores = pca_scores, figure_expression = expression,
    probe_audit = cohort_map, target_detail = target_detail,
    counts = data.frame(
      cohort = cohort,
      samples_frozen = nrow(blind_samples), samples_qc_excluded = sum(qc$exclude_qc),
      samples_retained = sum(keep_samples), probes_input = nrow(expression),
      probes_unique_current_entrez = sum(cohort_map$mapping_status == "unique_current"),
      genes_after_blind_collapse = nrow(feature_data), detectable_genes = nrow(detectable_features),
      primary_targets_total = nrow(targets), primary_targets_detectable = coverage_n,
      primary_target_coverage = coverage_fraction,
      preprocessing_action = transformed$action,
      stringsAsFactors = FALSE
    )
  )
}

main <- function() {
  required_packages <- c("AnnotationDbi", "org.Hs.eg.db", "limma")
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Restore frozen renv; missing: ", paste(missing, collapse = ", "), call. = FALSE)
  manifest_path <- file.path("data", "clean", "sample_manifest_frozen.tsv")
  target_path <- file.path("data", "clean", "bfff_targets_primary.tsv")
  annotation_path <- file.path("data", "raw", "GPL10558", "GPL10558.annot.gz")
  required_files <- c(manifest_path, target_path, annotation_path)
  if (!all(file.exists(required_files))) stop("Task 5 immutable inputs are missing", call. = FALSE)

  # Extract only the previously frozen inclusion decision and technical fields,
  # then destroy the phenotype-bearing manifest before any expression is loaded.
  manifest <- read.delim(manifest_path, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
  cohort_samples <- lapply(c("GSE105450", "GSE103842"), function(cohort) {
    idx <- manifest$series_accession == cohort & manifest$meta_eligible_confirmatory & manifest$include
    make_blind_sample_data(manifest[idx, c("sample_id", "platform", "batch"), drop = FALSE])
  })
  names(cohort_samples) <- c("GSE105450", "GSE103842")
  expected <- c(GSE105450 = 122L, GSE103842 = 74L)
  observed <- vapply(cohort_samples, nrow, integer(1))
  if (!identical(observed, expected)) stop("Frozen confirmatory sample count mismatch", call. = FALSE)
  rm(manifest)
  invisible(gc(verbose = FALSE))

  targets <- read.delim(target_path, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
  targets$entrez_id <- as.character(targets$entrez_id)
  if (nrow(targets) != 520L || anyDuplicated(targets$entrez_id)) stop("Expected 520 unique frozen primary Entrez targets", call. = FALSE)
  annotation <- read_gpl10558_annotation(annotation_path)
  annotation_map <- build_current_probe_map(annotation)

  paths <- list(
    GSE105450 = list(
      series = file.path("data", "raw", "GSE105450", "GSE105450_series_matrix.txt.gz"),
      detection = file.path("data", "raw", "GSE105450", "GSE105450_non-normalized.txt.gz")
    ),
    GSE103842 = list(
      series = file.path("data", "raw", "GSE103842", "GSE103842_series_matrix.txt.gz"),
      detection = file.path("data", "raw", "GSE103842", "GSE103842_non-normalized_data.txt.gz")
    )
  )
  outputs <- list()
  written <- character()
  for (cohort in names(paths)) {
    if (!all(file.exists(unlist(paths[[cohort]])))) stop(cohort, ": processed input absent", call. = FALSE)
    result <- run_cohort(
      cohort, paths[[cohort]]$series, paths[[cohort]]$detection,
      annotation_map, cohort_samples[[cohort]], targets
    )
    derived_path <- file.path("data", "derived", paste0(cohort, "_expression.rds"))
    qc_path <- file.path("results", "tables", paste0(cohort, "_qc.tsv"))
    map_path <- file.path("results", "tables", paste0(cohort, "_probe_entrez_mapping_audit.tsv"))
    coverage_path <- file.path("results", "tables", paste0(cohort, "_target_coverage.tsv"))
    figure_path <- file.path("results", "figures", "qc", paste0(cohort, "_qc.pdf"))
    dir.create(dirname(derived_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result$derived, derived_path, compress = "xz")
    write_tsv(result$qc, qc_path)
    write_tsv(result$probe_audit, map_path)
    write_tsv(result$target_detail, coverage_path)
    write_qc_figure(cohort, result$figure_expression, result$qc, result$pca_scores, figure_path)
    outputs[[cohort]] <- result$counts
    written <- c(written, derived_path, qc_path, map_path, coverage_path, figure_path)
  }
  counts <- do.call(rbind, outputs)
  write_tsv(counts, file.path("results", "tables", "confirmatory_qc_mapping_summary.tsv"))
  written <- c(written, file.path("results", "tables", "confirmatory_qc_mapping_summary.tsv"))

  annotation_manifest <- data.frame(
    platform = "GPL10558",
    official_url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
    relative_path = annotation_path,
    bytes = as.numeric(file.info(annotation_path)$size),
    sha256 = sha256_file(annotation_path),
    downloaded_at = format(file.info(annotation_path)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    immutable_mode = "0444",
    stringsAsFactors = FALSE
  )
  annotation_manifest_path <- file.path("logs", "checksums", "gpl10558_annotation_sha256.tsv")
  write_tsv(annotation_manifest, annotation_manifest_path)
  written <- c(written, annotation_manifest_path)
  checksum_inputs <- c(
    manifest_path, target_path, annotation_path,
    paths$GSE105450$series, paths$GSE105450$detection,
    paths$GSE103842$series, paths$GSE103842$detection,
    file.path("analysis", "R", "04_prepare_microarray.R"),
    file.path("analysis", "tests", "test_prepare_microarray.R"),
    file.path("analysis", "tests", "test_prepare_microarray_outputs.R"),
    file.path("analysis", "tests", "integration_prepare_microarray_determinism.R"),
    file.path("docs", "science-superpowers", "preregistrations", "amendments", "2026-08-05-processed-qc-detection-filter.md"),
    file.path("docs", "science-superpowers", "preregistrations", "amendments", "2026-08-05-scaled-mad-blind-qc.md"),
    file.path("logs", "session_info", "task5_supersession.md")
  )
  checksum <- data.frame(
    artifact = c(checksum_inputs, written),
    role = c(rep("input", length(checksum_inputs)), rep("output", length(written))),
    bytes = as.numeric(file.info(c(checksum_inputs, written))$size),
    sha256 = vapply(c(checksum_inputs, written), sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_tsv(checksum, file.path("logs", "checksums", "task5_sha256.tsv"))
  dir.create(file.path("logs", "session_info"), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    strip_trailing_whitespace(capture.output(sessionInfo())),
    file.path("logs", "session_info", "task5_session_info.txt"), useBytes = TRUE
  )
  writeLines(
    c(
      "Task 5: processed-only blind QC and Entrez mapping",
      "Confirmatory cohorts only: GSE105450, GSE103842",
      "Expression source: author background-subtracted and average-scaled GEO series matrix",
      "Transformation: log2(pmax(intensity,1)); no second between-array normalization",
      "Detection source: processed non-normalized supplementary table Detection Pval columns only",
      "Blindness: phenotype-bearing manifest removed before expression load; derived sample_data contains sample_id/platform/technical_batch only",
      "QC metric scale: stats::mad constant=1.4826; two-sided abs(x-median)>5*MAD; MAD=0 uses any nonzero deviation",
      "QC exclusion: >=2 blind scaled-MAD flags; any non-finite base/metric stops the cohort",
      "Probe mapping: GPL10558 official Gene ID -> current org.Hs.eg.db Entrez; multi/empty/obsolete excluded",
      "Probe collapse: highest all-sample mean expression, no case label/effect/P value",
      "QC figure first panel: all pre-QC samples and probes; PDF timestamp=FALSE, bg=white",
      "Source archives read: 0; network downloads beyond GPL10558.annot.gz: 0",
      "Supersession: commit 5948947 raw-MAD QC artifacts invalidated and regenerated under bf603ad"
    ),
    file.path("logs", "session_info", "task5_commands.log"), useBytes = TRUE
  )
  print(counts)
}

if (!identical(Sys.getenv("TASK5_FUNCTIONS_ONLY", unset = "false"), "true")) main()
