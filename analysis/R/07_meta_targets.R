#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

task9_cohorts <- c("GSE105450", "GSE103842")
task9_direction <- "RSV_minus_healthy"

read_tsv_task9 <- function(path) {
  read.delim(
    path, sep = "\t", quote = "", check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8"
  )
}

write_tsv_task9 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    na = "NA", fileEncoding = "UTF-8"
  )
}

sha256_file_task9 <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 failed for ", path, call. = FALSE)
  sub("[[:space:]].*$", "", output[[1L]])
}

assert_manifest_hash_task9 <- function(manifest, artifact) {
  row <- manifest[manifest$artifact == artifact, , drop = FALSE]
  if (nrow(row) != 1L) stop("Frozen checksum row absent or duplicated: ", artifact, call. = FALSE)
  if (!file.exists(artifact)) stop("Frozen input absent: ", artifact, call. = FALSE)
  observed <- sha256_file_task9(artifact)
  if (!identical(observed, as.character(row$sha256))) {
    stop("Frozen SHA-256 mismatch: ", artifact, call. = FALSE)
  }
  observed
}

validate_meta_input <- function(gene_effects, expected_cohorts = task9_cohorts) {
  if (!identical(as.character(expected_cohorts), task9_cohorts)) {
    stop("Task 9 expected_cohorts must be exactly GSE105450 and GSE103842", call. = FALSE)
  }
  required <- c(
    "gene_id", "gene_symbol", "cohort", "direction", "log2FC", "SE", "moderated_t"
  )
  absent <- setdiff(required, names(gene_effects))
  if (length(absent)) stop("Meta input missing: ", paste(absent, collapse = ", "), call. = FALSE)
  if (!nrow(gene_effects)) stop("Meta input is empty", call. = FALSE)
  if (anyNA(gene_effects$gene_id) || any(!nzchar(as.character(gene_effects$gene_id)))) {
    stop("Meta input requires non-empty gene_id", call. = FALSE)
  }
  if (anyNA(gene_effects$gene_symbol) || any(!nzchar(as.character(gene_effects$gene_symbol)))) {
    stop("Meta input requires non-empty gene_symbol", call. = FALSE)
  }
  if (anyDuplicated(gene_effects[c("gene_id", "cohort")])) {
    stop("Meta input contains duplicate gene-cohort rows", call. = FALSE)
  }
  if (any(!gene_effects$cohort %in% expected_cohorts)) {
    stop("Task 9 permits confirmatory cohorts only", call. = FALSE)
  }
  if (!setequal(unique(as.character(gene_effects$cohort)), expected_cohorts)) {
    stop("Meta input must contain exactly the two expected cohorts", call. = FALSE)
  }
  if (anyNA(gene_effects$direction) || any(gene_effects$direction != task9_direction)) {
    stop("All effects must use RSV_minus_healthy direction", call. = FALSE)
  }
  if (any(!is.finite(gene_effects$log2FC))) {
    stop("Meta input requires finite log2FC", call. = FALSE)
  }
  if (any(!is.finite(gene_effects$SE) | gene_effects$SE <= 0)) {
    stop("Meta input requires positive finite SE", call. = FALSE)
  }
  if (any(!is.finite(gene_effects$moderated_t))) {
    stop("Meta input requires finite moderated_t", call. = FALSE)
  }
  if (!isTRUE(all.equal(
    as.numeric(gene_effects$moderated_t),
    as.numeric(gene_effects$log2FC / gene_effects$SE),
    tolerance = 1e-10, check.attributes = FALSE
  ))) {
    stop("Task 7 posterior SE does not correspond to log2FC and moderated_t", call. = FALSE)
  }
  per_gene <- split(gene_effects$cohort, as.character(gene_effects$gene_id))
  complete <- vapply(
    per_gene,
    function(x) length(x) == 2L && setequal(as.character(x), expected_cohorts),
    logical(1)
  )
  if (!all(complete)) {
    stop("Every gene must contain exactly the two expected cohorts", call. = FALSE)
  }
  symbols <- split(as.character(gene_effects$gene_symbol), as.character(gene_effects$gene_id))
  if (any(vapply(symbols, function(x) length(unique(x)) != 1L, logical(1)))) {
    stop("A gene_id cannot map to conflicting official symbols", call. = FALSE)
  }
  gene_effects
}

extract_meta_row_task9 <- function(fit, gene_id, gene_symbol, method_label,
                                   expected_cohorts, include_prediction = FALSE) {
  prediction_lb <- NA_real_
  prediction_ub <- NA_real_
  if (include_prediction) {
    prediction <- stats::predict(fit, level = 95)
    prediction_lb <- as.numeric(prediction$pi.lb)
    prediction_ub <- as.numeric(prediction$pi.ub)
  }
  data.frame(
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    k = as.integer(fit$k),
    estimate = as.numeric(fit$b[[1L]]),
    meta_se = as.numeric(fit$se),
    ci_lb = as.numeric(fit$ci.lb),
    ci_ub = as.numeric(fit$ci.ub),
    prediction_lb = prediction_lb,
    prediction_ub = prediction_ub,
    p_value = as.numeric(fit$pval),
    tau2 = as.numeric(fit$tau2),
    Q = as.numeric(fit$QE),
    Q_p_value = as.numeric(fit$QEp),
    I2 = as.numeric(fit$I2),
    method = method_label,
    direction = task9_direction,
    expected_cohorts = paste(expected_cohorts, collapse = "|"),
    heterogeneity_caveat = "k=2; tau2, Q, I2 and prediction interval are unstable and descriptive",
    stringsAsFactors = FALSE
  )
}

fit_target_gene_meta <- function(gene_effects, expected_cohorts = task9_cohorts) {
  if (!requireNamespace("metafor", quietly = TRUE)) stop("Frozen metafor package is unavailable", call. = FALSE)
  gene_effects <- validate_meta_input(gene_effects, expected_cohorts)
  genes <- unique(as.character(gene_effects$gene_id))
  rows <- lapply(genes, function(gene_id) {
    one <- gene_effects[gene_effects$gene_id == gene_id, , drop = FALSE]
    one <- one[match(expected_cohorts, one$cohort), , drop = FALSE]
    symbol <- unique(as.character(one$gene_symbol))
    knha <- metafor::rma.uni(
      yi = one$log2FC, sei = one$SE, method = "REML", test = "knha"
    )
    wald <- metafor::rma.uni(
      yi = one$log2FC, sei = one$SE, method = "REML", test = "z"
    )
    fixed <- metafor::rma.uni(
      yi = one$log2FC, sei = one$SE, method = "FE", test = "z"
    )
    list(
      knha = extract_meta_row_task9(
        knha, gene_id, symbol, "REML_HKSJ", expected_cohorts, include_prediction = TRUE
      ),
      wald = extract_meta_row_task9(
        wald, gene_id, symbol, "REML_Wald", expected_cohorts, include_prediction = FALSE
      ),
      fixed = extract_meta_row_task9(
        fixed, gene_id, symbol, "fixed_inverse_variance", expected_cohorts,
        include_prediction = FALSE
      )
    )
  })
  combine <- function(name) {
    out <- do.call(rbind, lapply(rows, `[[`, name))
    rownames(out) <- NULL
    out$fdr_bh <- stats::p.adjust(out$p_value, method = "BH")
    out
  }
  list(knha = combine("knha"), wald = combine("wald"), fixed = combine("fixed"))
}

load_task9_inputs <- function() {
  task7_checksum_path <- file.path("logs", "checksums", "task7_sha256.tsv")
  task8_checksum_path <- file.path("logs", "checksums", "task8_sha256.tsv")
  if (!all(file.exists(c(task7_checksum_path, task8_checksum_path)))) {
    stop("Frozen Task 7/8 checksum manifests are required", call. = FALSE)
  }
  task7_checksum <- read_tsv_task9(task7_checksum_path)
  task8_checksum <- read_tsv_task9(task8_checksum_path)
  target_path <- file.path("data", "clean", "bfff_targets_primary.tsv")
  effect_paths <- file.path("results", "cohort", paste0(task9_cohorts, "_gene_effects.tsv"))
  decision_path <- file.path("results", "cohort", "target_set_primary_decision_partial.tsv")
  vapply(effect_paths, function(path) assert_manifest_hash_task9(task7_checksum, path), character(1))
  assert_manifest_hash_task9(task8_checksum, target_path)
  assert_manifest_hash_task9(task8_checksum, decision_path)

  decision <- read_tsv_task9(decision_path)
  if (nrow(decision) != 1L || !identical(decision$h1_supported, FALSE) ||
      !identical(decision$primary_decision, "H1_not_supported_camera_component_failed")) {
    stop("Task 8 non-rescuing primary decision boundary is not intact", call. = FALSE)
  }
  targets <- read_tsv_task9(target_path)
  required_targets <- c("entrez_id", "gene_symbol")
  if (!all(required_targets %in% names(targets)) || nrow(targets) != 520L ||
      anyDuplicated(targets$entrez_id) || anyNA(targets$entrez_id) ||
      any(!nzchar(as.character(targets$entrez_id)))) {
    stop("Frozen 520-gene Entrez primary target set is invalid", call. = FALSE)
  }
  effects <- lapply(effect_paths, read_tsv_task9)
  names(effects) <- task9_cohorts
  for (cohort in task9_cohorts) {
    x <- effects[[cohort]]
    if (anyDuplicated(x$gene_id) || !identical(unique(as.character(x$cohort)), cohort)) {
      stop("Task 7 effect identity is invalid for ", cohort, call. = FALSE)
    }
  }
  common_ids <- Reduce(
    intersect,
    c(list(as.character(targets$entrez_id)), lapply(effects, function(x) as.character(x$gene_id)))
  )
  common_ids <- as.character(targets$entrez_id[targets$entrez_id %in% common_ids])
  if (!length(common_ids)) stop("No primary targets are estimable in both cohorts", call. = FALSE)
  symbol_map <- setNames(as.character(targets$gene_symbol), as.character(targets$entrez_id))
  meta_input <- do.call(rbind, lapply(task9_cohorts, function(cohort) {
    x <- effects[[cohort]]
    x <- x[x$gene_id %in% common_ids, c(
      "gene_id", "cohort", "direction", "log2FC", "SE", "moderated_t"
    ), drop = FALSE]
    x$gene_symbol <- unname(symbol_map[as.character(x$gene_id)])
    x[, c("gene_id", "gene_symbol", "cohort", "direction", "log2FC", "SE", "moderated_t")]
  }))
  meta_input <- meta_input[order(match(meta_input$gene_id, common_ids),
                                 match(meta_input$cohort, task9_cohorts)), , drop = FALSE]
  rownames(meta_input) <- NULL
  meta_input <- validate_meta_input(meta_input, task9_cohorts)
  list(
    targets = targets,
    effects = effects,
    meta_input = meta_input,
    common_ids = common_ids,
    target_path = target_path,
    effect_paths = effect_paths,
    decision_path = decision_path,
    task7_checksum_path = task7_checksum_path,
    task8_checksum_path = task8_checksum_path
  )
}

make_task9_diagnostics <- function(inputs, fits) {
  wide <- reshape(
    inputs$meta_input[, c("gene_id", "cohort", "log2FC")],
    idvar = "gene_id", timevar = "cohort", direction = "wide"
  )
  direction_consistent <- sign(wide[[paste0("log2FC.", task9_cohorts[[1L]])]]) ==
    sign(wide[[paste0("log2FC.", task9_cohorts[[2L]])]])
  estimates <- fits$knha$estimate
  quantiles <- stats::quantile(estimates, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  data.frame(
    metric = c(
      "frozen_primary_targets", "common_targets_meta_analyzed", "expected_cohorts_per_gene",
      "direction_consistent_genes", "hksj_fdr_lt_0_05", "hksj_positive_estimates",
      "hksj_negative_estimates", "hksj_estimate_min", "hksj_estimate_q1",
      "hksj_estimate_median", "hksj_estimate_q3", "hksj_estimate_max",
      "confirmatory_leave_one_out_status", "primary_H1_rescue_allowed"
    ),
    value = as.character(c(
      nrow(inputs$targets), length(inputs$common_ids), 2L, sum(direction_consistent),
      sum(fits$knha$fdr_bh < 0.05), sum(estimates > 0), sum(estimates < 0),
      quantiles, "not_applicable_k2_confirmatory", FALSE
    )),
    note = c(
      "Frozen before confirmatory outcomes", "Intersection only; no P-value filtering",
      paste(task9_cohorts, collapse = "|"),
      "Same sign of cohort-specific RSV-minus-healthy log2FC; descriptive",
      "BH family is all common frozen primary targets", "Descriptive", "Descriptive",
      rep("REML/HKSJ pooled log2FC distribution", 5L),
      "Leave-one-out is not interpreted with two confirmatory cohorts",
      "FALSE; gene-level meta cannot rescue the failed Task 8 primary camera decision"
    ),
    stringsAsFactors = FALSE
  )
}

main_task9 <- function() {
  inputs <- load_task9_inputs()
  fits <- fit_target_gene_meta(inputs$meta_input, task9_cohorts)
  for (name in names(fits)) {
    x <- fits[[name]]
    if (nrow(x) != length(inputs$common_ids) || any(x$k != 2L) ||
        any(!is.finite(x$estimate)) || any(!is.finite(x$meta_se) | x$meta_se <= 0) ||
        any(!is.finite(x$ci_lb)) || any(!is.finite(x$ci_ub)) ||
        any(x$ci_lb > x$estimate | x$estimate > x$ci_ub) ||
        any(!is.finite(x$p_value) | x$p_value < 0 | x$p_value > 1) ||
        !isTRUE(all.equal(x$fdr_bh, p.adjust(x$p_value, method = "BH"), tolerance = 1e-15))) {
      stop("Invalid meta output for ", name, call. = FALSE)
    }
  }
  if (any(!is.finite(fits$knha$prediction_lb)) || any(!is.finite(fits$knha$prediction_ub)) ||
      any(fits$knha$prediction_lb > fits$knha$estimate) ||
      any(fits$knha$estimate > fits$knha$prediction_ub)) {
    stop("REML/HKSJ prediction intervals are invalid", call. = FALSE)
  }

  output_paths <- c(
    knha = file.path("results", "meta", "target_gene_meta_reml_knha.tsv"),
    wald = file.path("results", "meta", "target_gene_meta_reml_wald.tsv"),
    fixed = file.path("results", "meta", "target_gene_meta_fixed.tsv"),
    loo = file.path("results", "meta", "target_gene_leave_one_out.tsv"),
    diagnostics = file.path("results", "meta", "target_gene_meta_diagnostics.tsv")
  )
  write_tsv_task9(fits$knha, output_paths[["knha"]])
  write_tsv_task9(fits$wald, output_paths[["wald"]])
  write_tsv_task9(fits$fixed, output_paths[["fixed"]])
  leave_one_out <- data.frame(
    analysis_scope = "GSE105450_and_GSE103842_confirmatory_only",
    k = 2L,
    status = "not_applicable_k2_confirmatory",
    reason = "Pre-registration forbids interpretation of leave-one-out with two confirmatory cohorts",
    direction = task9_direction,
    stringsAsFactors = FALSE
  )
  write_tsv_task9(leave_one_out, output_paths[["loo"]])
  diagnostics <- make_task9_diagnostics(inputs, fits)
  write_tsv_task9(diagnostics, output_paths[["diagnostics"]])

  session_path <- file.path("logs", "session_info", "task9_session_info.txt")
  commands_path <- file.path("logs", "session_info", "task9_commands.log")
  dir.create(dirname(session_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())), session_path, useBytes = TRUE)
  writeLines(c(
    "Task 9: pre-registered target-gene meta-analysis",
    "Command: RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_CONFIG_NAMESPACES_CHECK=FALSE Rscript analysis/R/07_meta_targets.R",
    paste0("Confirmatory cohorts only: ", paste(task9_cohorts, collapse = ", ")),
    paste0("Direction: ", task9_direction),
    paste0("Frozen primary targets: ", nrow(inputs$targets)),
    paste0("Common finite target effects meta-analyzed: ", length(inputs$common_ids)),
    "No cohort-level P-value prefiltering; every common target belongs to one BH family",
    "Primary estimator: metafor::rma.uni(yi=log2FC, sei=SE, method=REML, test=knha)",
    "Sensitivities: REML/Wald and fixed inverse-variance",
    "No camera, random gene sets, pathway analysis, GSE38900, GSE188427, or confirmatory leave-one-out",
    "Task 8 primary H1 already failed and cannot be rescued by this secondary gene-level meta-analysis",
    "k=2: heterogeneity estimates and prediction intervals are unstable and descriptive"
  ), commands_path, useBytes = TRUE)

  script_path <- file.path("analysis", "R", "07_meta_targets.R")
  unit_test_path <- file.path("analysis", "tests", "test_meta_targets.R")
  output_test_path <- file.path("analysis", "tests", "test_meta_target_outputs.R")
  checksum_artifacts <- c(
    script_path, unit_test_path, output_test_path,
    inputs$target_path, inputs$effect_paths, inputs$decision_path,
    inputs$task7_checksum_path, inputs$task8_checksum_path,
    unname(output_paths), session_path, commands_path
  )
  if (!all(file.exists(checksum_artifacts))) {
    stop("Task 9 checksum artifact is absent", call. = FALSE)
  }
  checksum <- data.frame(
    artifact = checksum_artifacts,
    role = c(rep("input", 9L), rep("output", length(checksum_artifacts) - 9L)),
    bytes = as.numeric(file.info(checksum_artifacts)$size),
    sha256 = vapply(checksum_artifacts, sha256_file_task9, character(1)),
    stringsAsFactors = FALSE
  )
  write_tsv_task9(checksum, file.path("logs", "checksums", "task9_sha256.tsv"))

  print(diagnostics, row.names = FALSE)
  invisible(list(fits = fits, diagnostics = diagnostics))
}

if (!identical(Sys.getenv("TASK9_SKIP_MAIN"), "1")) main_task9()
