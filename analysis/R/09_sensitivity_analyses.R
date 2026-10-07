#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

task11_contract_path <- file.path("analysis", "R", "task11_design_contract.R")
if (!file.exists(task11_contract_path)) stop("Task 11 design contract is absent", call. = FALSE)
source(task11_contract_path, local = FALSE)

task7_script <- file.path("analysis", "R", "05_fit_cohort_models.R")
old_task7_skip <- Sys.getenv("TASK7_SKIP_MAIN", unset = NA_character_)
Sys.setenv(TASK7_SKIP_MAIN = "1")
source(task7_script, local = FALSE)
if (is.na(old_task7_skip)) Sys.unsetenv("TASK7_SKIP_MAIN") else Sys.setenv(TASK7_SKIP_MAIN = old_task7_skip)

task11_direction <- "RSV_minus_healthy"
task11_primary_n <- 520L
task11_expanded_n <- 1270L

task11_read_tsv <- function(path) read.delim(
  path, sep = "\t", quote = "", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8"
)

task11_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    na = "NA", fileEncoding = "UTF-8"
  )
}

task11_sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 failed for ", path, call. = FALSE)
  sub("[[:space:]].*$", "", output[[1L]])
}

task11_read_series_expression <- function(path) {
  if (!grepl("series_matrix[.]txt[.]gz$", path)) stop("Only a processed series matrix is allowed", call. = FALSE)
  con <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sample_line <- NULL
  processing <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("series_matrix_table_begin not found", call. = FALSE)
    if (identical(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_geo_accession\t")) sample_line <- line
    if (startsWith(line, "!Sample_data_processing\t")) processing <- c(processing, line)
  }
  parse_values <- function(x) {
    values <- strsplit(x, "\t", fixed = TRUE)[[1L]][-1L]
    gsub('^"|"$', "", values)
  }
  if (is.null(sample_line)) stop("Series sample IDs are absent", call. = FALSE)
  sample_id <- parse_values(sample_line)
  table <- read.delim(
    con, sep = "\t", header = TRUE, quote = "\"", comment.char = "!",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!identical(names(table)[-1L], sample_id)) stop("Series expression/sample identity mismatch", call. = FALSE)
  expression <- as.matrix(table[-1L])
  storage.mode(expression) <- "double"
  rownames(expression) <- as.character(table[[1L]])
  if (any(!is.finite(expression))) stop("Processed expression contains non-finite values", call. = FALSE)
  list(
    expression = expression,
    sample_id = sample_id,
    processing = unique(unlist(lapply(processing, parse_values), use.names = FALSE))
  )
}

task11_fit_limma <- function(y, design, coefficient, cohort, direction) {
  y <- as.matrix(y)
  design <- as.matrix(design)
  if (any(!is.finite(y)) || any(!is.finite(design)) || ncol(y) != nrow(design)) {
    stop(cohort, ": invalid limma input", call. = FALSE)
  }
  if (qr(design)$rank != ncol(design)) stop(cohort, ": design is not full rank", call. = FALSE)
  coefficient_index <- if (is.character(coefficient)) match(coefficient, colnames(design)) else coefficient
  if (length(coefficient_index) != 1L || is.na(coefficient_index)) stop(cohort, ": coefficient is absent", call. = FALSE)
  fit <- limma::eBayes(limma::lmFit(y, design))
  se <- fit$stdev.unscaled[, coefficient_index] * sqrt(fit$s2.post)
  if (any(!is.finite(se) | se <= 0)) stop(cohort, ": invalid moderated SE", call. = FALSE)
  p <- unname(fit$p.value[, coefficient_index])
  effects <- data.frame(
    gene_id = rownames(y), cohort = cohort, direction = direction,
    log2FC = unname(fit$coefficients[, coefficient_index]),
    SE = unname(se), moderated_t = unname(fit$t[, coefficient_index]),
    p_value = p, fdr_bh = p.adjust(p, method = "BH"),
    average_expression = rowMeans(y), stringsAsFactors = FALSE
  )
  if (anyDuplicated(effects$gene_id) || any(!is.finite(effects$log2FC)) ||
      !isTRUE(all.equal(effects$moderated_t, effects$log2FC / effects$SE, tolerance = 1e-10))) {
    stop(cohort, ": gene-effect identity failed", call. = FALSE)
  }
  list(fit = fit, gene_effects = effects)
}

task11_make_contrast <- function(design, coefficient) {
  index <- match(coefficient, colnames(design))
  if (is.na(index)) stop("Coefficient is absent from design", call. = FALSE)
  out <- rep(0, ncol(design))
  out[[index]] <- 1
  out
}

task11_camera_row <- function(y, design, coefficient, targets, cohort, target_set,
                              contrast_label, n_model, role, family) {
  targets <- unique(as.character(targets))
  index <- which(rownames(y) %in% targets)
  n_detectable <- length(index)
  coverage <- n_detectable / length(targets)
  interpretable <- n_detectable >= 10L && coverage >= 0.50
  if (!interpretable) {
    return(list(
      result = data.frame(
        family = family, analysis_type = "camera", analysis_id = paste(cohort, target_set, sep = "__"),
        cohort_or_combination = cohort, target_set = target_set, contrast = contrast_label,
        n_independent = n_model, n_frozen = length(targets), n_detectable = n_detectable,
        coverage_fraction = coverage, direction = NA_character_, p_value = NA_real_,
        inter_gene_correlation = NA_real_, role = role,
        interpretability = "uninterpretable_by_frozen_coverage_rule",
        primary_H1_rescue_allowed = FALSE, stringsAsFactors = FALSE
      ),
      coverage = data.frame(
        cohort = cohort, target_set = target_set, n_frozen = length(targets),
        n_detectable = n_detectable, coverage_fraction = coverage,
        camera_index_length = n_detectable,
        interpretability = "uninterpretable_by_frozen_coverage_rule", stringsAsFactors = FALSE
      )
    ))
  }
  camera <- task11_camera(y, index, design, task11_make_contrast(design, coefficient))
  if (nrow(camera) != 1L || !camera$Direction[[1L]] %in% c("Up", "Down")) {
    stop(cohort, ": unexpected camera result", call. = FALSE)
  }
  list(
    result = data.frame(
      family = family, analysis_type = "camera", analysis_id = paste(cohort, target_set, sep = "__"),
      cohort_or_combination = cohort, target_set = target_set, contrast = contrast_label,
      n_independent = n_model, n_frozen = length(targets), n_detectable = n_detectable,
      coverage_fraction = coverage, direction = as.character(camera$Direction[[1L]]),
      p_value = as.numeric(camera$PValue[[1L]]),
      inter_gene_correlation = if ("Correlation" %in% names(camera)) as.numeric(camera$Correlation[[1L]]) else NA_real_,
      role = role,
      interpretability = "interpretable", primary_H1_rescue_allowed = FALSE,
      stringsAsFactors = FALSE
    ),
    coverage = data.frame(
      cohort = cohort, target_set = target_set, n_frozen = length(targets),
      n_detectable = n_detectable, coverage_fraction = coverage,
      camera_index_length = n_detectable, interpretability = "interpretable",
      stringsAsFactors = FALSE
    )
  )
}

task11_stouffer_row <- function(camera_rows, analysis_id, family, role) {
  if (nrow(camera_rows) != 2L || anyNA(camera_rows$p_value) ||
      any(camera_rows$interpretability != "interpretable")) {
    stop(analysis_id, ": Stouffer requires two interpretable camera rows", call. = FALSE)
  }
  combined <- task11_signed_stouffer(
    camera_rows$p_value, camera_rows$direction, camera_rows$n_independent
  )
  data.frame(
    family = family, analysis_type = "signed_stouffer", analysis_id = analysis_id,
    cohort_or_combination = paste(camera_rows$cohort_or_combination, collapse = "+"),
    target_set = unique(camera_rows$target_set), contrast = unique(camera_rows$contrast),
    n_independent = sum(camera_rows$n_independent), n_frozen = unique(camera_rows$n_frozen),
    n_detectable = NA_integer_, coverage_fraction = NA_real_,
    direction = ifelse(combined$z_combined >= 0, "Up", "Down"),
    p_value = combined$p_two_sided, inter_gene_correlation = NA_real_, role = role,
    interpretability = "interpretable", primary_H1_rescue_allowed = FALSE,
    z_combined = combined$z_combined,
    component_z = paste(format(combined$z_i, digits = 17), collapse = "|"),
    component_weights = paste(format(combined$weight, digits = 17), collapse = "|"),
    stringsAsFactors = FALSE
  )
}

task11_rbind_fill <- function(parts) {
  columns <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(x) {
    missing <- setdiff(columns, names(x))
    for (name in missing) x[[name]] <- NA
    x[columns]
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

task11_prepare_gse105450_hospital <- function(manifest) {
  object <- readRDS(file.path("data", "derived", "GSE105450_expression.rds"))
  eligible <- manifest[
    manifest$series_accession == "GSE105450" & manifest$include &
      manifest$severity %in% c("hospitalized", "healthy"), , drop = FALSE
  ]
  ids <- intersect(colnames(object$expression), eligible$sample_id)
  if (length(ids) != 89L) stop("GSE105450 hospital frozen sample count mismatch", call. = FALSE)
  sample <- eligible[match(ids, eligible$sample_id), , drop = FALSE]
  sample$case_status <- factor(ifelse(sample$status == "RSV_case", "RSV", "healthy"), levels = c("healthy", "RSV"))
  sample$sex <- factor(sample$sex, levels = c("female", "male"))
  sample$technical_batch <- factor(sample$batch)
  design <- model.matrix(~ case_status + age_months + sex + technical_batch, data = sample)
  rownames(design) <- ids
  if (qr(design)$rank != 6L || ncol(design) != 6L) stop("GSE105450 hospital design is not frozen 6/6", call. = FALSE)
  expression <- object$expression[, ids, drop = FALSE]
  fit <- task11_fit_limma(expression, design, "case_statusRSV", "GSE105450_hospital", task11_direction)
  list(expression = expression, sample = sample, design = design, coefficient = "case_statusRSV", fit = fit)
}

task11_prepare_gse103842_main <- function(manifest) {
  input <- load_cohort_inputs("GSE103842", manifest)
  prepared <- prepare_model_data(input$expression, input$sample_data, "GSE103842")
  fit <- task11_fit_limma(prepared$expression, prepared$design, prepared$coefficient, "GSE103842", task11_direction)
  list(
    expression = prepared$expression, sample = prepared$sample_data,
    design = prepared$design, coefficient = prepared$coefficient, fit = fit
  )
}

task11_prepare_gse105450_main <- function(manifest) {
  input <- load_cohort_inputs("GSE105450", manifest)
  prepared <- prepare_model_data(input$expression, input$sample_data, "GSE105450")
  fit <- task11_fit_limma(prepared$expression, prepared$design, prepared$coefficient, "GSE105450", task11_direction)
  list(
    expression = prepared$expression, sample = prepared$sample_data,
    design = prepared$design, coefficient = prepared$coefficient, fit = fit
  )
}

task11_prepare_gse188427 <- function(manifest) {
  path <- file.path("data", "raw", "GSE188427", "GSE188427_series_matrix.txt.gz")
  series <- task11_read_series_expression(path)
  ids <- manifest$sample_id[manifest$series_accession == "GSE188427" & manifest$include]
  expression <- series$expression[, match(ids, series$sample_id), drop = FALSE]
  rownames(expression) <- sub("_at$", "", rownames(expression))
  if (ncol(expression) != 198L || nrow(expression) != 18604L || anyDuplicated(rownames(expression)) ||
      any(!grepl("^[0-9]+$", rownames(expression)))) {
    stop("GSE188427 processed ENTREZG matrix identity failed", call. = FALSE)
  }
  qc_result <- task11_gse188427_blind_qc(expression)
  qc <- qc_result$qc
  keep <- !qc$exclude_qc
  expression <- expression[, keep, drop = FALSE]
  sample <- manifest[match(colnames(expression), manifest$sample_id), , drop = FALSE]
  sample$case_status <- factor(ifelse(sample$status == "RSV_case", "RSV", "healthy"), levels = c("healthy", "RSV"))
  design <- model.matrix(~ case_status, data = sample)
  rownames(design) <- sample$sample_id
  fit <- task11_fit_limma(expression, design, "case_statusRSV", "GSE188427", task11_direction)
  qc$cohort <- "GSE188427"
  qc$qc_decision <- ifelse(qc$exclude_qc, "exclude_two_or_more_blind_scaled_5MAD_flags", "retain")
  list(
    expression = expression, sample = sample, design = design,
    coefficient = "case_statusRSV", fit = fit, qc = qc,
    preprocessing = paste(series$processing, collapse = " | ")
  )
}

task11_prepare_gse103119 <- function(manifest) {
  series_path <- file.path("data", "raw", "GSE103119", "GSE103119_series_matrix.txt.gz")
  detection_path <- file.path("data", "raw", "GSE103119", "GSE103119_non-normalized.txt.gz")
  annotation_path <- file.path("data", "raw", "GPL10558", "GPL10558.annot.gz")
  series <- task11_read_series_expression(series_path)
  ids <- manifest$sample_id[manifest$series_accession == "GSE103119" & manifest$include]
  expression <- log2(pmax(series$expression[, match(ids, series$sample_id), drop = FALSE], 1))
  detection_all <- read_detection_pvalues(detection_path)
  sample_manifest <- manifest[match(ids, manifest$sample_id), , drop = FALSE]
  array_id <- sub("^.* ", "", sample_manifest$title)
  detection <- detection_all[match(rownames(expression), rownames(detection_all)), match(array_id, colnames(detection_all)), drop = FALSE]
  rownames(detection) <- rownames(expression)
  colnames(detection) <- colnames(expression)
  if (any(!is.finite(expression)) || any(!is.finite(detection))) stop("GSE103119 non-finite processed input", call. = FALSE)
  qc_base <- calculate_blind_qc(expression, detection)
  qc <- apply_two_metric_qc_gate(qc_base, setdiff(names(qc_base), "sample_id"))
  keep <- !qc$exclude_qc
  annotation_map <- build_current_probe_map(read_gpl10558_annotation(annotation_path))
  probe_map <- annotation_map[match(rownames(expression), annotation_map$probe_id), , drop = FALSE]
  absent <- is.na(probe_map$probe_id)
  if (any(absent)) {
    probe_map[absent, ] <- data.frame(
      probe_id = rownames(expression)[absent], source_gene_symbol = NA_character_,
      source_entrez_id = NA_character_, entrez_id = NA_character_, official_symbol = NA_character_,
      mapping_status = "annotation_absent", symbol_conflict = FALSE, stringsAsFactors = FALSE
    )
  }
  probe_map$probe_id <- rownames(expression)
  collapsed <- collapse_probes_blind(expression[, keep, drop = FALSE], probe_map)
  selected <- collapsed$selected_probes
  selected_detection <- detection[match(selected$probe_id, rownames(detection)), keep, drop = FALSE]
  detectable <- rowMeans(selected_detection < 0.05) >= 0.10
  expression_gene <- collapsed$expression[detectable, , drop = FALSE]
  sample <- manifest[match(colnames(expression_gene), manifest$sample_id), , drop = FALSE]
  sample$case_status <- factor(ifelse(sample$status == "RSV_case", "RSV", "healthy"), levels = c("healthy", "RSV"))
  sample$sex <- factor(sample$sex, levels = c("female", "male"))
  design <- model.matrix(~ case_status + age_months + sex, data = sample)
  rownames(design) <- sample$sample_id
  if (qr(design)$rank != ncol(design) || ncol(design) != 4L) stop("GSE103119 design is not frozen 4/4", call. = FALSE)
  fit <- task11_fit_limma(expression_gene, design, "case_statusRSV", "GSE103119", task11_direction)
  qc$cohort <- "GSE103119"
  qc$qc_decision <- ifelse(qc$exclude_qc, "exclude_two_or_more_blind_scaled_5MAD_flags", "retain")
  list(
    expression = expression_gene, sample = sample, design = design,
    coefficient = "case_statusRSV", fit = fit, qc = qc,
    probe_mapping = selected, preprocessing = "log2_pmax1_only; no second between-array normalization"
  )
}

task11_prepare_gse155925 <- function(manifest, targets_expanded) {
  path <- file.path("data", "raw", "GSE155925", "GSE155925_Raw_counts_matrix.txt.gz")
  table <- read.delim(path, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
  raw_id <- sub("^.*:", "", as.character(table[[1L]]))
  normalized <- task11_normalize_ensembl(raw_id)
  counts <- as.matrix(table[-1L])
  storage.mode(counts) <- "double"
  if (any(!is.finite(counts)) || any(counts < 0) || any(abs(counts - round(counts)) > 0)) {
    stop("GSE155925 count matrix is not finite non-negative integer data", call. = FALSE)
  }
  counts <- rowsum(counts, group = normalized, reorder = FALSE)
  all_manifest <- manifest[manifest$series_accession == "GSE155925", , drop = FALSE]
  if (!identical(colnames(counts), as.character(all_manifest$title))) {
    stop("GSE155925 Case-column/title/GSM ordering mismatch", call. = FALSE)
  }
  sample <- all_manifest[all_manifest$include, , drop = FALSE]
  counts <- counts[, match(sample$title, colnames(counts)), drop = FALSE]
  colnames(counts) <- sample$sample_id
  sample$virus_group <- factor(ifelse(sample$status == "RSV_case", "RSV", "other"), levels = c("other", "RSV"))
  sample$sex <- factor(sample$sex, levels = c("female", "male"))
  batch <- do.call(rbind, strsplit(sample$batch, ";", fixed = TRUE))
  sample$hospital_batch <- factor(batch[, 1L])
  sample$enrollment_batch <- factor(batch[, 2L])
  design <- model.matrix(~ virus_group + age_months + sex + hospital_batch + enrollment_batch, data = sample)
  rownames(design) <- sample$sample_id
  if (nrow(design) != 48L || qr(design)$rank != 6L || ncol(design) != 6L) {
    stop("GSE155925 design is not frozen 6/6", call. = FALSE)
  }
  background <- task11_prepare_rnaseq_background(counts, design)
  selection <- task11_select_rnaseq_target_units(
    background$background_ensembl, background$mean_logcpm, targets_expanded
  )
  voom <- limma::voom(background$dge, design = design, plot = FALSE)
  fit <- task11_fit_limma(voom$E, design, "virus_groupRSV", "GSE155925", "single_RSV_minus_other_single_virus")
  list(
    expression = voom, sample = sample, design = design, coefficient = "virus_groupRSV",
    fit = fit, selection = selection$selection,
    ambiguous_ensembl = selection$ambiguous_ensembl,
    counts_input_rows = nrow(table), unique_ensembl_rows = nrow(counts),
    filtered_ensembl_rows = nrow(voom$E)
  )
}

task11_assert_meta_pair <- function(cohorts) {
  cohorts <- as.character(cohorts)
  allowed <- list(
    c("GSE105450_hospital", "GSE103842"),
    c("GSE188427", "GSE103842")
  )
  if (length(cohorts) != 2L || !any(vapply(allowed, identical, logical(1), cohorts))) {
    stop("Only the frozen hospital or replacement pair is allowed; GSE105450 and GSE188427 cannot coexist", call. = FALSE)
  }
  cohorts
}

task11_fit_meta_pair <- function(effects, targets, cohorts, combo) {
  cohorts <- task11_assert_meta_pair(cohorts)
  ids <- Reduce(intersect, c(list(as.character(targets$entrez_id)), lapply(cohorts, function(x) {
    as.character(effects[[x]]$gene_id)
  })))
  ids <- as.character(targets$entrez_id[targets$entrez_id %in% ids])
  if (!length(ids)) stop(combo, ": no common primary targets", call. = FALSE)
  symbol <- setNames(as.character(targets$gene_symbol), as.character(targets$entrez_id))
  rows <- lapply(ids, function(id) {
    one <- do.call(rbind, lapply(cohorts, function(cohort) {
      effects[[cohort]][effects[[cohort]]$gene_id == id, , drop = FALSE]
    }))
    one <- one[match(cohorts, one$cohort), , drop = FALSE]
    models <- list(
      REML_HKSJ = metafor::rma.uni(yi = one$log2FC, sei = one$SE, method = "REML", test = "knha"),
      REML_Wald = metafor::rma.uni(yi = one$log2FC, sei = one$SE, method = "REML", test = "z"),
      fixed_inverse_variance = metafor::rma.uni(yi = one$log2FC, sei = one$SE, method = "FE", test = "z")
    )
    lapply(names(models), function(method) {
      model <- models[[method]]
      prediction <- if (method == "REML_HKSJ") predict(model, level = 95) else NULL
      data.frame(
        combination = combo, gene_id = id, gene_symbol = unname(symbol[[id]]),
        cohorts = paste(cohorts, collapse = "|"), k = 2L, method = method,
        estimate = as.numeric(model$b[[1L]]), meta_se = as.numeric(model$se),
        ci_lb = as.numeric(model$ci.lb), ci_ub = as.numeric(model$ci.ub),
        prediction_lb = if (is.null(prediction)) NA_real_ else as.numeric(prediction$pi.lb),
        prediction_ub = if (is.null(prediction)) NA_real_ else as.numeric(prediction$pi.ub),
        p_value = as.numeric(model$pval), tau2 = as.numeric(model$tau2),
        Q = as.numeric(model$QE), Q_p_value = as.numeric(model$QEp), I2 = as.numeric(model$I2),
        direction = task11_direction,
        heterogeneity_caveat = "k=2; tau2, Q, I2 and prediction interval are unstable and descriptive",
        leave_one_out = "not_applicable_k2", primary_H1_rescue_allowed = FALSE,
        stringsAsFactors = FALSE
      )
    })
  })
  all <- do.call(rbind, unlist(rows, recursive = FALSE))
  all$fdr_bh <- ave(all$p_value, all$method, FUN = function(x) p.adjust(x, method = "BH"))
  split(all, all$method)
}

task11_model_diagnostic <- function(cohort, data, formula, contrast, role) data.frame(
  cohort = cohort, role = role, model_formula = formula, contrast = contrast,
  n_model = nrow(data$design), n_reference = sum(data$sample[[if (cohort == "GSE155925") "virus_group" else "case_status"]] %in% c("healthy", "other")),
  n_case = sum(data$sample[[if (cohort == "GSE155925") "virus_group" else "case_status"]] == "RSV"),
  design_rank = qr(data$design)$rank, design_columns = ncol(data$design),
  genes_modeled = nrow(data$fit$gene_effects), primary_H1_rescue_allowed = FALSE,
  stringsAsFactors = FALSE
)

main_task11 <- function() {
  required_packages <- c("limma", "edgeR", "metafor", "AnnotationDbi", "org.Hs.eg.db")
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Restore frozen renv; missing: ", paste(missing, collapse = ", "), call. = FALSE)

  amendment_path <- file.path(
    "docs", "science-superpowers", "preregistrations", "amendments",
    "2026-08-05-task11-sensitivity-operationalization.md"
  )
  manifest_path <- file.path("data", "clean", "sample_manifest_frozen.tsv")
  primary_path <- file.path("data", "clean", "bfff_targets_primary.tsv")
  expanded_path <- file.path("data", "clean", "bfff_targets_sensitivity.tsv")
  decision_path <- file.path("results", "cohort", "target_set_primary_decision_partial.tsv")
  required_inputs <- c(amendment_path, manifest_path, primary_path, expanded_path, decision_path)
  if (!all(file.exists(required_inputs))) stop("Task 11 frozen inputs are incomplete", call. = FALSE)
  if (task11_sha256(manifest_path) != "5434ee45d31e80a38db81696118b823e45b5aab2cd2dc628269ab29ee417a859" ||
      task11_sha256(expanded_path) != "573c5f577dc071af45289ff68e051d762f7b4403288dc78267e62dab8755abcd") {
    stop("Task 11 frozen input hash mismatch", call. = FALSE)
  }
  decision <- task11_read_tsv(decision_path)
  if (nrow(decision) != 1L || !identical(decision$h1_supported, FALSE)) {
    stop("Failed primary H1 boundary is not intact", call. = FALSE)
  }
  manifest <- task11_read_tsv(manifest_path)
  primary <- task11_read_tsv(primary_path)
  expanded <- task11_read_tsv(expanded_path)
  primary$entrez_id <- as.character(primary$entrez_id)
  expanded$entrez_id <- as.character(expanded$entrez_id)
  # The frozen target table encodes an absent alias as an empty TSV field.
  # Restore that representation after the generic reader maps empty fields to NA.
  expanded$ensembl_gene_ids[is.na(expanded$ensembl_gene_ids)] <- ""
  if (nrow(primary) != task11_primary_n || nrow(expanded) != task11_expanded_n ||
      anyDuplicated(primary$entrez_id) || anyDuplicated(expanded$entrez_id) ||
      !all(primary$entrez_id %in% expanded$entrez_id)) {
    stop("Frozen target identity failed", call. = FALSE)
  }

  hospital <- task11_prepare_gse105450_hospital(manifest)
  gse105450 <- task11_prepare_gse105450_main(manifest)
  gse103842 <- task11_prepare_gse103842_main(manifest)
  replacement <- task11_prepare_gse188427(manifest)
  gse103119 <- task11_prepare_gse103119(manifest)
  gse155925 <- task11_prepare_gse155925(manifest, expanded)

  camera_parts <- list()
  coverage_parts <- list()
  add_camera <- function(data, targets, cohort, target_set, contrast, role, family, target_index = NULL) {
    if (is.null(target_index)) {
      result <- task11_camera_row(
        data$expression, data$design, data$coefficient, targets, cohort,
        target_set, contrast, nrow(data$design), role, family
      )
    } else {
      index_info <- task11_target_index(target_index, targets, rownames(data$expression$E))
      n_detectable <- index_info$coverage_unique_entrez
      coverage <- n_detectable / length(targets)
      if (n_detectable < 10L || coverage < 0.50) stop(cohort, ": uninterpretable RNA-seq target coverage", call. = FALSE)
      cam <- task11_camera(
        data$expression, index_info$index, data$design,
        task11_make_contrast(data$design, data$coefficient)
      )
      result <- list(
        result = data.frame(
          family = family, analysis_type = "camera", analysis_id = paste(cohort, target_set, sep = "__"),
          cohort_or_combination = cohort, target_set = target_set, contrast = contrast,
          n_independent = nrow(data$design), n_frozen = length(targets), n_detectable = n_detectable,
          coverage_fraction = coverage, direction = as.character(cam$Direction[[1L]]),
          p_value = as.numeric(cam$PValue[[1L]]),
          inter_gene_correlation = if ("Correlation" %in% names(cam)) as.numeric(cam$Correlation[[1L]]) else NA_real_,
          role = role,
          interpretability = "interpretable", primary_H1_rescue_allowed = FALSE,
          stringsAsFactors = FALSE
        ),
        coverage = data.frame(
          cohort = cohort, target_set = target_set, n_frozen = length(targets),
          n_detectable = n_detectable, coverage_fraction = coverage,
          camera_index_length = length(index_info$index), interpretability = "interpretable",
          stringsAsFactors = FALSE
        )
      )
    }
    camera_parts[[length(camera_parts) + 1L]] <<- result$result
    coverage_parts[[length(coverage_parts) + 1L]] <<- result$coverage
    result$result
  }

  s1_hospital <- add_camera(hospital, primary$entrez_id, "GSE105450_hospital", "primary_520", task11_direction, "sensitivity", "S1")
  s1_replacement <- add_camera(replacement, primary$entrez_id, "GSE188427", "primary_520", task11_direction, "exploratory_nonH1_unadjusted", "S1")
  primary_103842 <- task11_camera_row(
    gse103842$expression, gse103842$design, gse103842$coefficient, primary$entrez_id,
    "GSE103842", "primary_520", task11_direction, nrow(gse103842$design),
    "frozen_confirmatory_component_reused_for_sensitivity_combination", "component_only"
  )$result
  s1_hospital_meta <- task11_stouffer_row(
    rbind(s1_hospital, primary_103842), "hospital_plus_GSE103842__primary_520", "S1", "sensitivity"
  )
  s1_replacement_meta <- task11_stouffer_row(
    rbind(s1_replacement, primary_103842), "replacement_plus_GSE103842__primary_520", "S1", "exploratory_nonH1"
  )
  camera_parts <- c(camera_parts, list(s1_hospital_meta, s1_replacement_meta))

  s2_105450 <- add_camera(gse105450, expanded$entrez_id, "GSE105450", "expanded_1270", task11_direction, "sensitivity", "S2")
  s2_103842 <- add_camera(gse103842, expanded$entrez_id, "GSE103842", "expanded_1270", task11_direction, "sensitivity", "S2")
  camera_parts[[length(camera_parts) + 1L]] <- task11_stouffer_row(
    rbind(s2_105450, s2_103842), "GSE105450_plus_GSE103842__expanded_1270", "S2", "sensitivity"
  )

  add_camera(gse103119, primary$entrez_id, "GSE103119", "primary_520", task11_direction, "small_sample_exploratory", "E1")
  add_camera(gse103119, expanded$entrez_id, "GSE103119", "expanded_1270", task11_direction, "small_sample_exploratory", "E1")
  add_camera(gse155925, primary$entrez_id, "GSE155925", "primary_520", "single_RSV_minus_other_single_virus", "cross_virus_exploratory", "E1", gse155925$selection)
  add_camera(gse155925, expanded$entrez_id, "GSE155925", "expanded_1270", "single_RSV_minus_other_single_virus", "cross_virus_exploratory", "E1", gse155925$selection)

  set_results <- task11_rbind_fill(camera_parts)
  set_results$p_holm <- ave(set_results$p_value, set_results$family, FUN = function(x) p.adjust(x, method = "holm"))
  set_results <- set_results[order(match(set_results$family, c("S1", "S2", "E1")), set_results$analysis_id), , drop = FALSE]
  if (nrow(set_results) != 11L || anyNA(set_results$p_value) ||
      any(set_results$primary_H1_rescue_allowed)) stop("Task 11 set-family completeness failed", call. = FALSE)

  effects <- list(
    GSE105450_hospital = hospital$fit$gene_effects,
    GSE188427 = replacement$fit$gene_effects,
    GSE103842 = gse103842$fit$gene_effects
  )
  meta_hospital <- task11_fit_meta_pair(
    effects, primary, c("GSE105450_hospital", "GSE103842"), "hospital_plus_GSE103842"
  )
  meta_replacement <- task11_fit_meta_pair(
    effects, primary, c("GSE188427", "GSE103842"), "replacement_plus_GSE103842"
  )

  model_diagnostics <- rbind(
    task11_model_diagnostic("GSE105450_hospital", hospital, "~ case_status + age_months + sex + technical_batch", task11_direction, "sensitivity"),
    task11_model_diagnostic("GSE188427", replacement, "~ case_status", task11_direction, "exploratory_nonH1_unadjusted"),
    task11_model_diagnostic("GSE103119", gse103119, "~ case_status + age_months + sex", task11_direction, "small_sample_exploratory"),
    task11_model_diagnostic("GSE155925", gse155925, "~ virus_group + age_months + sex + hospital_batch + enrollment_batch", "single_RSV_minus_other_single_virus", "cross_virus_exploratory")
  )
  qc_audit <- task11_rbind_fill(list(replacement$qc, gse103119$qc))
  coverage <- do.call(rbind, coverage_parts)
  rownames(coverage) <- NULL

  alias_map <- task11_build_alias_map(expanded)
  alias_audit <- rbind(
    data.frame(
      ensembl_id = alias_map$ambiguous_ensembl,
      mapping_status = "globally_ambiguous_removed", note = "maps_to_multiple_Entrez_in_full_1270_alias_table",
      stringsAsFactors = FALSE
    ),
    data.frame(
      ensembl_id = gse155925$selection$ensembl_id,
      mapping_status = "unambiguous_selected", note = "highest_blind_mean_logCPM_per_Entrez_after_filterByExpr",
      stringsAsFactors = FALSE
    )
  )
  not_implemented <- data.frame(
    analysis = c("CBC_adjustment", "cell_proportion_adjustment", "GSE155925_severity"),
    status = "not_implemented_by_frozen_contract",
    reason = c(
      "No CBC/WBC/absolute differential fields in the five local series matrices",
      "No measured cell proportions or frozen pediatric whole-blood deconvolution reference",
      "No usable severity gradient in the fixed 48 single-virus hospitalized children"
    ),
    primary_H1_rescue_allowed = FALSE, stringsAsFactors = FALSE
  )
  input_paths <- c(
    amendment_path, manifest_path, primary_path, expanded_path, decision_path,
    file.path("data", "derived", "GSE105450_expression.rds"),
    file.path("data", "derived", "GSE103842_expression.rds"),
    file.path("data", "raw", "GSE188427", "GSE188427_series_matrix.txt.gz"),
    file.path("data", "raw", "GSE103119", "GSE103119_series_matrix.txt.gz"),
    file.path("data", "raw", "GSE103119", "GSE103119_non-normalized.txt.gz"),
    file.path("data", "raw", "GPL10558", "GPL10558.annot.gz"),
    file.path("data", "raw", "GSE155925", "GSE155925_Raw_counts_matrix.txt.gz")
  )
  if (any(grepl("RAW[.]tar|[.]CEL|FASTQ|SRA", input_paths, ignore.case = TRUE))) {
    stop("Forbidden source file entered the Task 11 input list", call. = FALSE)
  }
  access_audit <- data.frame(
    path = input_paths, input_type = c(
      rep("frozen_contract_or_table", 5L), rep("processed_derived_expression", 2L),
      "processed_series_matrix", "processed_series_matrix", "processed_detection_p_table",
      "small_official_platform_annotation", "author_gene_level_count_matrix"
    ),
    bytes = as.numeric(file.info(input_paths)$size), sha256 = vapply(input_paths, task11_sha256, character(1)),
    network_download = FALSE, source_archive_read = FALSE, stringsAsFactors = FALSE
  )

  output_dir <- file.path("results", "sensitivity")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  outputs <- character()
  emit <- function(x, name) {
    path <- file.path(output_dir, name)
    task11_write_tsv(x, path)
    outputs <<- c(outputs, path)
  }
  emit(set_results, "all_sensitivity_results.tsv")
  emit(set_results, "camera_and_stouffer_results.tsv")
  emit(coverage, "target_coverage.tsv")
  emit(model_diagnostics, "cohort_model_diagnostics.tsv")
  emit(qc_audit, "cohort_qc_audit.tsv")
  emit(hospital$fit$gene_effects, "gene_effects_GSE105450_hospital.tsv")
  emit(replacement$fit$gene_effects, "gene_effects_GSE188427.tsv")
  emit(gse103119$fit$gene_effects, "gene_effects_GSE103119.tsv")
  emit(gse155925$fit$gene_effects, "gene_effects_GSE155925.tsv")
  emit(meta_hospital$REML_HKSJ, "meta_hospital_reml_hksj.tsv")
  emit(meta_hospital$REML_Wald, "meta_hospital_reml_wald.tsv")
  emit(meta_hospital$fixed_inverse_variance, "meta_hospital_fixed.tsv")
  emit(meta_replacement$REML_HKSJ, "meta_replacement_reml_hksj.tsv")
  emit(meta_replacement$REML_Wald, "meta_replacement_reml_wald.tsv")
  emit(meta_replacement$fixed_inverse_variance, "meta_replacement_fixed.tsv")
  emit(gse155925$selection, "rnaseq_target_selection.tsv")
  emit(alias_audit, "rnaseq_alias_audit.tsv")
  emit(not_implemented, "not_implemented.tsv")
  emit(access_audit, "input_access_audit.tsv")

  session_path <- file.path("logs", "session_info", "task11_session_info.txt")
  commands_path <- file.path("logs", "session_info", "task11_commands.log")
  dir.create(dirname(session_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())), session_path, useBytes = TRUE)
  writeLines(c(
    "Task 11: frozen sensitivity and exploratory analyses",
    "Command: RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_CONFIG_NAMESPACES_CHECK=FALSE Rscript analysis/R/09_sensitivity_analyses.R",
    "Inputs: existing processed matrices, Task 5 derived RDS, GPL10558 annotation, and GSE155925 author gene-level count matrix only",
    "Forbidden source archives read: 0; RAW.tar/CEL/FASTQ/SRA read: 0; network downloads: 0",
    "GSE188427: blind processed-only QC then unadjusted ~case_status; exploratory non-H1",
    "GSE105450 hospital: fixed adjusted model; GSE103119: fixed small-sample adjusted model",
    "GSE155925: unique Ensembl aggregation, full-design filterByExpr, TMM, voom, single RSV minus other single virus",
    "Camera wrapper: directional=TRUE, inter.gene.cor=NA_real_; no cameraPR or fixed correlation",
    "Set multiplicity: Holm within S1, S2, and E1; gene-level BH within each meta combination/method",
    "CBC/cell proportions and GSE155925 severity: not implemented by frozen contract",
    "Primary H1 remains not supported and no Task 11 result may rescue it"
  ), commands_path, useBytes = TRUE)
  outputs <- c(outputs, session_path, commands_path)

  checksum_inputs <- c(
    file.path("analysis", "R", "09_sensitivity_analyses.R"),
    task11_contract_path, task7_script,
    file.path("analysis", "tests", "test_task11_sensitivity_analysis.R"),
    file.path("analysis", "tests", "test_task11_sensitivity_outputs.R"),
    input_paths
  )
  checksum_artifacts <- c(checksum_inputs, outputs)
  if (!all(file.exists(checksum_artifacts))) stop("Task 11 checksum artifact is absent", call. = FALSE)
  checksum <- data.frame(
    artifact = checksum_artifacts,
    role = c(rep("input", length(checksum_inputs)), rep("output", length(outputs))),
    bytes = as.numeric(file.info(checksum_artifacts)$size),
    sha256 = vapply(checksum_artifacts, task11_sha256, character(1)),
    stringsAsFactors = FALSE
  )
  task11_write_tsv(checksum, file.path("logs", "checksums", "task11_sha256.tsv"))

  print(set_results[, c(
    "family", "analysis_id", "direction", "p_value", "p_holm",
    "n_independent", "n_detectable", "coverage_fraction"
  )], row.names = FALSE)
  invisible(list(set_results = set_results, coverage = coverage, diagnostics = model_diagnostics))
}

if (!identical(Sys.getenv("TASK11_SKIP_MAIN"), "1")) main_task11()
