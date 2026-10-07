options(stringsAsFactors = FALSE)

assert_frozen_packages <- function() {
  required <- c("limma", "metafor")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing frozen package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

build_rsv_design <- function(expression, sample_data, cohort_id, direction) {
  if (!identical(direction, "RSV_minus_healthy")) {
    stop("Direction must be RSV_minus_healthy", call. = FALSE)
  }
  if (!is.matrix(expression) || !is.numeric(expression) || !length(expression)) {
    stop("Expression must be a non-empty numeric matrix", call. = FALSE)
  }
  if (is.null(rownames(expression)) || is.null(colnames(expression)) ||
      any(!nzchar(rownames(expression))) || any(!nzchar(colnames(expression))) ||
      anyDuplicated(rownames(expression)) || anyDuplicated(colnames(expression))) {
    stop("Expression requires unique non-empty gene and sample IDs", call. = FALSE)
  }
  if (any(!is.finite(expression))) stop("Expression contains non-finite values", call. = FALSE)

  required <- c("sample_id", "subject_id", "cohort", "platform", "case_status", "age_months", "sex")
  absent <- setdiff(required, names(sample_data))
  if (length(absent)) stop("Sample data missing: ", paste(absent, collapse = ", "), call. = FALSE)
  if (!identical(as.character(sample_data$sample_id), colnames(expression))) {
    stop("Sample data must match expression columns in the same order", call. = FALSE)
  }
  if (anyDuplicated(sample_data$sample_id)) stop("Duplicate sample_id detected", call. = FALSE)
  if (anyDuplicated(sample_data$subject_id)) stop("Duplicate subject_id detected", call. = FALSE)
  if (!identical(unique(as.character(sample_data$cohort)), cohort_id)) {
    stop("Input must contain exactly the requested single cohort", call. = FALSE)
  }
  if (length(unique(as.character(sample_data$platform))) != 1L) {
    stop("Each fit must contain a single platform; cross-platform pooling is forbidden", call. = FALSE)
  }
  if (!setequal(unique(as.character(sample_data$case_status)), c("healthy", "RSV")) ||
      any(!sample_data$case_status %in% c("healthy", "RSV"))) {
    stop("case_status must contain only both healthy and RSV", call. = FALSE)
  }
  if (!is.numeric(sample_data$age_months) || any(!is.finite(sample_data$age_months)) ||
      any(sample_data$age_months < 0)) {
    stop("age_months must be finite, numeric, and non-negative", call. = FALSE)
  }
  if (!setequal(unique(as.character(sample_data$sex)), c("female", "male")) ||
      any(!sample_data$sex %in% c("female", "male"))) {
    stop("sex must contain only both female and male", call. = FALSE)
  }

  model_data <- data.frame(
    case_status = factor(sample_data$case_status, levels = c("healthy", "RSV")),
    age_months = sample_data$age_months,
    sex = factor(sample_data$sex, levels = c("female", "male")),
    stringsAsFactors = FALSE
  )
  design <- stats::model.matrix(~ case_status + age_months + sex, data = model_data)
  if (qr(design)$rank != ncol(design)) stop("Design matrix is not full rank", call. = FALSE)
  if (!"case_statusRSV" %in% colnames(design)) {
    stop("RSV-minus-healthy coefficient is absent", call. = FALSE)
  }
  rownames(design) <- sample_data$sample_id
  list(design = design, coefficient = "case_statusRSV")
}

fit_limma_cohort <- function(expression, sample_data, cohort_id, direction = "RSV_minus_healthy") {
  assert_frozen_packages()
  checked <- build_rsv_design(expression, sample_data, cohort_id, direction)
  fit <- limma::lmFit(expression, checked$design)
  fit <- limma::eBayes(fit)
  coefficient_index <- match(checked$coefficient, colnames(fit$coefficients))
  posterior_se <- fit$stdev.unscaled[, coefficient_index] * sqrt(fit$s2.post)
  if (any(!is.finite(posterior_se) | posterior_se <= 0)) {
    stop("Moderated standard errors must be finite and positive", call. = FALSE)
  }
  gene_effects <- data.frame(
    gene_id = rownames(expression),
    cohort = cohort_id,
    direction = direction,
    log2FC = unname(fit$coefficients[, coefficient_index]),
    SE = unname(posterior_se),
    moderated_t = unname(fit$t[, coefficient_index]),
    p_value = unname(fit$p.value[, coefficient_index]),
    fdr_bh = stats::p.adjust(unname(fit$p.value[, coefficient_index]), method = "BH"),
    average_expression = rowMeans(expression),
    stringsAsFactors = FALSE
  )
  list(
    gene_effects = gene_effects,
    coefficient = checked$coefficient,
    direction = direction,
    design_rank = qr(checked$design)$rank,
    design_columns = ncol(checked$design),
    n_samples = ncol(expression)
  )
}

run_camera_analysis <- function(...) {
  limma::camera(...)
}

test_camera_target_set <- function(expression, sample_data, frozen_targets, cohort_id,
                                   direction = "RSV_minus_healthy") {
  assert_frozen_packages()
  checked <- build_rsv_design(expression, sample_data, cohort_id, direction)
  frozen_targets <- as.character(frozen_targets)
  if (!length(frozen_targets) || anyNA(frozen_targets) || any(!nzchar(frozen_targets))) {
    stop("Frozen target IDs must be non-empty", call. = FALSE)
  }
  if (anyDuplicated(frozen_targets)) stop("Frozen target IDs must be unique", call. = FALSE)
  index <- which(rownames(expression) %in% frozen_targets)
  n_frozen <- length(frozen_targets)
  n_detectable <- length(index)
  coverage <- n_detectable / n_frozen
  if (n_detectable < 10L || coverage < 0.50) {
    reasons <- c(
      if (n_detectable < 10L) "n_detectable_below_10",
      if (coverage < 0.50) "coverage_below_0.50"
    )
    return(data.frame(
      cohort = cohort_id,
      contrast = direction,
      status = "uninterpretable",
      interpretability_reason = paste(reasons, collapse = ";"),
      n_frozen = n_frozen,
      n_detectable = n_detectable,
      coverage = coverage,
      n_genes = n_detectable,
      inter_gene_correlation = NA_real_,
      direction = NA_character_,
      p_value = NA_real_,
      signed_z = NA_real_,
      statistic_definition = paste0(
        "camera not run because frozen-target interpretability gate failed; ",
        "requires n_detectable>=10 and coverage>=0.50"
      ),
      stringsAsFactors = FALSE
    ))
  }
  camera_result <- run_camera_analysis(
    y = expression,
    index = list(frozen_target_set = index),
    design = checked$design,
    contrast = checked$coefficient,
    inter.gene.cor = NA_real_,
    sort = FALSE,
    directional = TRUE
  )
  p_value <- camera_result$PValue[[1L]]
  set_direction <- camera_result$Direction[[1L]]
  if (length(set_direction) != 1L || is.na(set_direction) || !set_direction %in% c("Up", "Down")) {
    stop("camera Direction must be Up or Down", call. = FALSE)
  }
  if (length(p_value) != 1L || !is.finite(p_value) || p_value <= 0 || p_value > 1) {
    stop("camera P value must be finite and in (0,1]", call. = FALSE)
  }
  direction_sign <- if (identical(set_direction, "Up")) 1 else -1
  signed_z <- direction_sign * stats::qnorm(p_value / 2, lower.tail = FALSE)
  data.frame(
    cohort = cohort_id,
    contrast = direction,
    status = "interpretable",
    interpretability_reason = "passed_n_detectable_and_coverage_gate",
    n_frozen = n_frozen,
    n_detectable = n_detectable,
    coverage = coverage,
    n_genes = camera_result$NGenes[[1L]],
    inter_gene_correlation = camera_result$Correlation[[1L]],
    direction = set_direction,
    p_value = p_value,
    signed_z = signed_z,
    statistic_definition = paste0(
      "signed_z=+qnorm(1-P/2) for camera Direction Up and -qnorm(1-P/2) for Down; ",
      "Up means target enrichment toward positive limma t statistics for RSV-minus-healthy"
    ),
    stringsAsFactors = FALSE
  )
}

meta_analyze_gene_effects <- function(gene_effects, expected_cohorts) {
  assert_frozen_packages()
  expected_cohorts <- as.character(expected_cohorts)
  if (!length(expected_cohorts) || anyNA(expected_cohorts) || any(!nzchar(expected_cohorts))) {
    stop("expected_cohorts must be non-empty", call. = FALSE)
  }
  if (anyDuplicated(expected_cohorts)) stop("expected_cohorts must be unique", call. = FALSE)
  required <- c("gene_id", "cohort", "direction", "log2FC", "SE")
  absent <- setdiff(required, names(gene_effects))
  if (length(absent)) stop("Meta input missing: ", paste(absent, collapse = ", "), call. = FALSE)
  if (anyNA(gene_effects$gene_id) || any(!nzchar(as.character(gene_effects$gene_id)))) {
    stop("Meta input requires non-empty gene_id", call. = FALSE)
  }
  if (anyNA(gene_effects$cohort) || any(!nzchar(as.character(gene_effects$cohort)))) {
    stop("Meta input requires non-empty cohort", call. = FALSE)
  }
  if (anyDuplicated(gene_effects[c("gene_id", "cohort")])) {
    stop("Meta input contains duplicate gene-cohort rows", call. = FALSE)
  }
  if (anyNA(gene_effects$direction) || any(gene_effects$direction != "RSV_minus_healthy")) {
    stop("All meta-analysis effects must be RSV_minus_healthy", call. = FALSE)
  }
  if (any(!is.finite(gene_effects$log2FC)) || any(!is.finite(gene_effects$SE) | gene_effects$SE <= 0)) {
    stop("Meta-analysis effects and SEs must be finite with positive SE", call. = FALSE)
  }
  genes <- unique(as.character(gene_effects$gene_id))
  results <- lapply(genes, function(gene_id) {
    one <- gene_effects[gene_effects$gene_id == gene_id, , drop = FALSE]
    if (nrow(one) != length(expected_cohorts) || !setequal(as.character(one$cohort), expected_cohorts)) {
      stop("Each gene cohort set must exactly equal expected_cohorts: ", gene_id, call. = FALSE)
    }
    one <- one[match(expected_cohorts, one$cohort), , drop = FALSE]
    fit <- metafor::rma.uni(
      yi = one$log2FC,
      sei = one$SE,
      method = "REML",
      test = "knha"
    )
    prediction <- stats::predict(fit, level = 95)
    data.frame(
      gene_id = gene_id,
      k = fit$k,
      estimate = unname(fit$b[[1L]]),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      prediction_lb = prediction$pi.lb,
      prediction_ub = prediction$pi.ub,
      p_value = fit$pval,
      tau2 = fit$tau2,
      Q = fit$QE,
      Q_p_value = fit$QEp,
      I2 = fit$I2,
      method = "metafor::rma.uni(method=REML,test=knha)",
      direction = "RSV_minus_healthy",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out$fdr_bh <- stats::p.adjust(out$p_value, method = "BH")
  out
}
