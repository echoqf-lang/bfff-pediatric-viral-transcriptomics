#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

revised_main_cohorts <- c("GSE77087", "GSE103842")
revised_analysis_label <- "post_outcome_investigator_selected_revised_main"
revised_direction <- "RSV_minus_healthy"
revised_seed <- 20260805L
revised_nrot <- 99999L
revised_root <- file.path("results", "revised_main_gse77087")

old_task5 <- Sys.getenv("TASK5_FUNCTIONS_ONLY", unset = NA_character_)
old_task7 <- Sys.getenv("TASK7_SKIP_MAIN", unset = NA_character_)
old_task8 <- Sys.getenv("TASK8_SKIP_MAIN", unset = NA_character_)
Sys.setenv(TASK5_FUNCTIONS_ONLY = "true", TASK7_SKIP_MAIN = "1", TASK8_SKIP_MAIN = "1")
source(file.path("analysis", "R", "04_prepare_microarray.R"), local = FALSE)
source(file.path("analysis", "R", "05_fit_cohort_models.R"), local = FALSE)
source(file.path("analysis", "R", "06_test_target_set.R"), local = FALSE)
if (is.na(old_task5)) Sys.unsetenv("TASK5_FUNCTIONS_ONLY") else Sys.setenv(TASK5_FUNCTIONS_ONLY = old_task5)
if (is.na(old_task7)) Sys.unsetenv("TASK7_SKIP_MAIN") else Sys.setenv(TASK7_SKIP_MAIN = old_task7)
if (is.na(old_task8)) Sys.unsetenv("TASK8_SKIP_MAIN") else Sys.setenv(TASK8_SKIP_MAIN = old_task8)
task8_seed <- revised_seed

read_header_values <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
  gsub('^"|"$', "", fields)
}

read_gse77087_metadata <- function(path) {
  con <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  fields <- list(characteristics = list())
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("GSE77087 series matrix metadata is incomplete", call. = FALSE)
    if (identical(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_geo_accession\t")) fields$sample_id <- read_header_values(line)
    if (startsWith(line, "!Sample_title\t")) fields$title <- read_header_values(line)
    if (startsWith(line, "!Sample_platform_id\t")) fields$platform <- read_header_values(line)
    if (startsWith(line, "!Sample_characteristics_ch1\t")) {
      values <- read_header_values(line)
      key <- tolower(trimws(sub(":.*$", "", values[[1L]])))
      fields$characteristics[[key]] <- trimws(sub("^[^:]+:", "", values))
    }
  }
  required <- c("sex", "age (mos)", "disease", "group", "batch", "tissue")
  if (!all(required %in% names(fields$characteristics))) stop("GSE77087 required phenotype fields are absent", call. = FALSE)
  n <- length(fields$sample_id)
  if (!all(c(length(fields$title), length(fields$platform), lengths(fields$characteristics)) == n)) {
    stop("GSE77087 metadata fields have unequal widths", call. = FALSE)
  }
  disease <- fields$characteristics$disease
  group <- tolower(fields$characteristics$group)
  out <- data.frame(
    sample_id = fields$sample_id,
    subject_id = fields$sample_id,
    cohort = "GSE77087",
    platform = fields$platform,
    case_status = ifelse(disease == "RSV", "RSV", ifelse(disease == "Healthy", "healthy", NA_character_)),
    age_months = suppressWarnings(as.numeric(fields$characteristics[["age (mos)"]])),
    sex = tolower(fields$characteristics$sex),
    technical_batch = fields$characteristics$batch,
    clinical_group = ifelse(group == "healthy control", "healthy", group),
    tissue = tolower(fields$characteristics$tissue),
    title = fields$title,
    stringsAsFactors = FALSE
  )
  if (nrow(out) != 104L || anyDuplicated(out$sample_id) || anyNA(out[c("case_status", "age_months", "sex")])) {
    stop("GSE77087 metadata fails the frozen 104-sample contract", call. = FALSE)
  }
  out
}

calculate_processed_blind_qc <- function(expression) {
  variable <- apply(expression, 1L, stats::var)
  top <- order(variable, decreasing = TRUE)[seq_len(min(5000L, length(variable)))]
  top_expression <- expression[top, , drop = FALSE]
  correlations <- stats::cor(top_expression, method = "spearman")
  diag(correlations) <- NA_real_
  pca <- stats::prcomp(t(top_expression), center = TRUE, scale. = FALSE)
  k <- min(5L, ncol(pca$x))
  standardized <- sweep(pca$x[, seq_len(k), drop = FALSE], 2L, pca$sdev[seq_len(k)], "/")
  weights <- limma::arrayWeights(expression, design = matrix(1, ncol(expression), 1L))
  data.frame(
    sample_id = colnames(expression),
    missing_rate = colMeans(!is.finite(expression)),
    median_expression = apply(expression, 2L, stats::median),
    expression_iqr = apply(expression, 2L, stats::IQR),
    median_sample_correlation = apply(correlations, 2L, stats::median, na.rm = TRUE),
    pca_distance_5pc = sqrt(rowSums(standardized^2)),
    log_array_weight = log(as.numeric(weights)),
    stringsAsFactors = FALSE
  )
}

filter_processed_features_blind <- function(expression) {
  expression <- as.matrix(expression)
  means <- rowMeans(expression)
  iqrs <- apply(expression, 1L, stats::IQR)
  mean_cutoff <- unname(stats::quantile(means, 0.10, names = FALSE, type = 7))
  keep <- is.finite(means) & is.finite(iqrs) & iqrs > 0 & means >= mean_cutoff
  list(
    expression = expression[keep, , drop = FALSE],
    feature_qc = data.frame(gene_id = rownames(expression), all_sample_mean = means,
                            all_sample_iqr = iqrs, mean_cutoff = mean_cutoff,
                            retained = keep, stringsAsFactors = FALSE)
  )
}

prepare_revised_model_data <- function(expression, sample_data, cohort_id) {
  prepare_model_data(expression, sample_data, cohort_id)
}

prepare_gse77087 <- function() {
  series_path <- file.path("data", "raw", "GSE77087", "GSE77087_series_matrix.txt.gz")
  metadata <- read_gse77087_metadata(series_path)
  series <- read_series_matrix_blind(series_path)
  transformed <- transform_processed_expression(series$expression)
  expression <- transformed$expression
  expression <- expression[, match(metadata$sample_id, colnames(expression)), drop = FALSE]
  if (anyNA(expression) || !identical(colnames(expression), metadata$sample_id)) stop("GSE77087 expression/metadata mismatch", call. = FALSE)

  qc_base <- calculate_processed_blind_qc(expression)
  qc <- apply_two_metric_qc_gate(qc_base, setdiff(names(qc_base), "sample_id"))
  qc$qc_decision <- ifelse(qc$exclude_qc, "exclude_two_or_more_blind_5_scaled_MAD_flags", "retain")
  retained <- !qc$exclude_qc

  annotation <- read_gpl10558_annotation(file.path("data", "raw", "GPL10558", "GPL10558.annot.gz"))
  probe_map <- build_current_probe_map(annotation)
  collapsed <- collapse_probes_blind(expression[, retained, drop = FALSE], probe_map)
  filtered <- filter_processed_features_blind(collapsed$expression)
  feature_data <- collapsed$selected_probes[match(rownames(filtered$expression), collapsed$selected_probes$entrez_id), , drop = FALSE]
  sample_data <- metadata[retained, c("sample_id", "subject_id", "cohort", "platform", "case_status",
                                     "age_months", "sex", "technical_batch", "clinical_group"), drop = FALSE]
  list(expression = filtered$expression, sample_data = sample_data, metadata = metadata,
       qc = qc, feature_data = feature_data, feature_qc = filtered$feature_qc,
       transform = transformed$action)
}

fit_revised_cohort <- function(expression, sample_data, cohort) {
  fit_task7_cohort(expression, sample_data, cohort, revised_direction)
}

load_gse103842_revised <- function() {
  manifest <- read_tsv(file.path("data", "clean", "sample_manifest_frozen.tsv"))
  input <- load_cohort_inputs("GSE103842", manifest)
  fit <- fit_revised_cohort(input$expression, input$sample_data, "GSE103842")
  prepared <- prepare_revised_model_data(input$expression, input$sample_data, "GSE103842")
  list(expression = prepared$expression, sample_data = prepared$sample_data, design = prepared$design,
       coefficient = prepared$coefficient, fit = fit)
}

make_model_bundle <- function(expression, sample_data, cohort) {
  prepared <- prepare_revised_model_data(expression, sample_data, cohort)
  fit <- fit_revised_cohort(expression, sample_data, cohort)
  list(expression = prepared$expression, sample_data = prepared$sample_data, design = prepared$design,
       coefficient = prepared$coefficient, formula = prepared$formula, fit = fit)
}

combine_revised_stouffer <- function(camera, sample_n) {
  if (!identical(as.character(camera$cohort), revised_main_cohorts)) stop("Revised Stouffer cohort order is invalid", call. = FALSE)
  weights <- sqrt(as.numeric(sample_n[revised_main_cohorts]))
  z <- sum(weights * camera$signed_z) / sqrt(sum(weights^2))
  data.frame(
    analysis_label = revised_analysis_label,
    cohorts = paste(revised_main_cohorts, collapse = "|"),
    direction_consistent = length(unique(camera$direction)) == 1L,
    weighted_stouffer_z = z,
    weighted_stouffer_two_sided_p = 2 * pnorm(-abs(z)),
    both_camera_p_lt_0_05 = all(camera$p_value < 0.05),
    original_preregistered_h1_changed = FALSE,
    stringsAsFactors = FALSE
  )
}

fit_revised_target_meta <- function(effects, targets) {
  common <- Reduce(intersect, lapply(effects, function(x) as.character(x$gene_id)))
  common <- intersect(as.character(targets$entrez_id), common)
  symbol <- setNames(as.character(targets$gene_symbol), as.character(targets$entrez_id))
  rows <- lapply(common, function(id) {
    one <- do.call(rbind, lapply(revised_main_cohorts, function(cohort) {
      x <- effects[[cohort]]
      x[x$gene_id == id, c("cohort", "log2FC", "SE"), drop = FALSE]
    }))
    fit <- metafor::rma.uni(yi = one$log2FC, sei = one$SE, method = "REML", test = "knha")
    data.frame(gene_id = id, gene_symbol = unname(symbol[id]), k = fit$k,
               estimate = as.numeric(fit$b), meta_se = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub,
               p_value = fit$pval, tau2 = fit$tau2, Q = fit$QE, Q_p_value = fit$QEp, I2 = fit$I2,
               method = "REML_HKSJ", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$fdr_bh <- p.adjust(out$p_value, method = "BH")
  out
}

write_revised_tsv <- function(x, name) {
  path <- file.path(revised_root, name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  path
}

write_revised_target_figure <- function(camera, stouffer) {
  path <- file.path(revised_root, "figures", "Fig2_revised_target_set.pdf")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  p_values <- c(camera$p_value, stouffer$weighted_stouffer_two_sided_p)
  labels <- c(as.character(camera$cohort), "Weighted Stouffer")
  colors <- c("#0072B2", "#56B4E9", "#444444")
  grDevices::pdf(path, width = 7.2, height = 5.4, useDingbats = FALSE,
                 timestamp = FALSE, bg = "white")
  old <- graphics::par(mar = c(4.5, 4.5, 1.5, 1))
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  heights <- -log10(p_values)
  mid <- graphics::barplot(heights, names.arg = labels, las = 1, col = colors,
                           border = NA, ylab = expression(-log[10](P)),
                           ylim = c(0, max(c(heights, -log10(0.05))) * 1.35))
  graphics::abline(h = -log10(0.05), lty = 2, col = "#D55E00")
  graphics::text(mid, heights, labels = formatC(p_values, format = "f", digits = 3),
                 pos = 3, cex = 0.9)
  path
}

main_gse77087_revised <- function() {
  if (!all(vapply(c("limma", "fgsea", "metafor"), requireNamespace, logical(1), quietly = TRUE))) {
    stop("Required frozen R packages are missing", call. = FALSE)
  }
  g77 <- prepare_gse77087()
  models <- list(
    GSE77087 = make_model_bundle(g77$expression, g77$sample_data, "GSE77087"),
    GSE103842 = load_gse103842_revised()
  )
  targets_table <- read_tsv(file.path("data", "clean", "bfff_targets_primary.tsv"))
  targets <- as.character(targets_table$entrez_id)

  effects <- lapply(models, function(x) x$fit$gene_effects)
  for (cohort in names(effects)) write_revised_tsv(effects[[cohort]], file.path("cohort", paste0(cohort, "_gene_effects.tsv")))
  diagnostics <- do.call(rbind, lapply(names(models), function(cohort) {
    x <- models[[cohort]]
    data.frame(cohort = cohort, analysis_label = revised_analysis_label, model_formula = x$fit$formula,
               n_complete = x$fit$n_complete, n_healthy = x$fit$n_healthy, n_RSV = x$fit$n_RSV,
               batch_included = x$fit$batch_included, genes_tested = nrow(x$fit$gene_effects), stringsAsFactors = FALSE)
  }))
  write_revised_tsv(diagnostics, file.path("cohort", "model_diagnostics.tsv"))
  write_revised_tsv(g77$metadata, file.path("cohort", "GSE77087_sample_manifest.tsv"))
  write_revised_tsv(g77$qc, file.path("cohort", "GSE77087_blind_qc.tsv"))
  write_revised_tsv(g77$feature_qc, file.path("cohort", "GSE77087_feature_filter.tsv"))

  camera <- do.call(rbind, lapply(revised_main_cohorts, function(cohort) {
    x <- models[[cohort]]
    run_camera_prepared(x$expression, x$design, x$coefficient, targets, cohort)
  }))
  rownames(camera) <- NULL
  sample_n <- setNames(vapply(models, function(x) nrow(x$sample_data), integer(1)), names(models))
  stouffer <- combine_revised_stouffer(camera, sample_n)
  write_revised_tsv(camera, file.path("target_set", "camera.tsv"))
  write_revised_tsv(stouffer, file.path("target_set", "weighted_stouffer.tsv"))
  table2 <- data.frame(
    analysis = c(as.character(camera$cohort), "Weighted Stouffer"),
    test = c(rep("CAMERA two-sided", nrow(camera)), "sample-size-weighted signed Stouffer"),
    direction = c(as.character(camera$direction), if (stouffer$direction_consistent) "concordant Up" else "discordant"),
    n_targets_evaluable = c(camera$n_detectable, NA_integer_),
    p_value = c(camera$p_value, stouffer$weighted_stouffer_two_sided_p),
    p_lt_0_05 = c(camera$p_value < 0.05, stouffer$weighted_stouffer_two_sided_p < 0.05),
    evidence_role = revised_analysis_label,
    stringsAsFactors = FALSE
  )
  write_revised_tsv(table2, file.path("tables", "Table2_revised_target_set.tsv"))
  write_revised_target_figure(camera, stouffer)

  set.seed(revised_seed)
  roast <- do.call(rbind, lapply(revised_main_cohorts, function(cohort) {
    x <- models[[cohort]]
    run_roast_prepared(x$expression, x$design, x$coefficient, targets, cohort, revised_nrot)
  }))
  set.seed(revised_seed)
  fgsea <- do.call(rbind, lapply(revised_main_cohorts, function(cohort) {
    run_fgsea_exploratory(models[[cohort]]$fit$gene_effects, targets, cohort)
  }))
  write_revised_tsv(roast, file.path("target_set", "roast.tsv"))
  write_revised_tsv(fgsea, file.path("target_set", "fgsea_exploratory.tsv"))

  meta <- fit_revised_target_meta(effects, targets_table)
  write_revised_tsv(meta, file.path("meta", "target_gene_reml_hksj.tsv"))

  hospital_data <- g77$sample_data$clinical_group %in% c("healthy", "inpatient")
  hospital <- make_model_bundle(g77$expression[, hospital_data, drop = FALSE],
                                g77$sample_data[hospital_data, , drop = FALSE], "GSE77087")
  hospital_camera <- run_camera_prepared(hospital$expression, hospital$design, hospital$coefficient, targets,
                                         "GSE77087_hospitalized_only")
  write_revised_tsv(hospital_camera, file.path("sensitivity", "GSE77087_hospitalized_only_camera.tsv"))

  summary <- data.frame(
    analysis_label = revised_analysis_label,
    gse77087_pre_qc_n = nrow(g77$metadata), gse77087_qc_excluded_n = sum(g77$qc$exclude_qc),
    gse77087_model_n = models$GSE77087$fit$n_complete,
    gse103842_model_n = models$GSE103842$fit$n_complete,
    gse77087_camera_p = camera$p_value[camera$cohort == "GSE77087"],
    gse103842_camera_p = camera$p_value[camera$cohort == "GSE103842"],
    direction_consistent = stouffer$direction_consistent,
    weighted_stouffer_p = stouffer$weighted_stouffer_two_sided_p,
    hksj_fdr_lt_0_05_n = sum(meta$fdr_bh < 0.05, na.rm = TRUE),
    original_preregistered_h1_changed = FALSE,
    stringsAsFactors = FALSE
  )
  write_revised_tsv(summary, "analysis_summary.tsv")
  dir.create(file.path("logs", "session_info"), recursive = TRUE, showWarnings = FALSE)
  writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())),
             file.path("logs", "session_info", "gse77087_revised_main_session_info.txt"), useBytes = TRUE)
  print(summary, row.names = FALSE)
  invisible(summary)
}

if (!identical(Sys.getenv("GSE77087_REVISED_SKIP_MAIN"), "1")) main_gse77087_revised()
