#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

task12_design_spec <- function() {
  list(
    cohorts = c("GSE105450", "GSE103842"),
    manifest_path = file.path("data", "clean", "sample_manifest_frozen.tsv"),
    qc_paths = c(
      GSE105450 = file.path("results", "tables", "GSE105450_qc.tsv"),
      GSE103842 = file.path("results", "tables", "GSE103842_qc.tsv")
    ),
    input_sha256 = c(
      sample_manifest_frozen = "5434ee45d31e80a38db81696118b823e45b5aab2cd2dc628269ab29ee417a859",
      GSE105450_qc = "e4417f3dd2d3fd0c83364a12c15f3b86dab21ee3c0bc117bd703aa0b6cf00e5e",
      GSE103842_qc = "3ed3ef0903ccf962dd0fc327c85e37885519f43e20b9f68eb551e75dd9782f80"
    ),
    expected = data.frame(
      cohort = c("GSE105450", "GSE103842"),
      n = c(122L, 73L),
      n_RSV = c(89L, 61L),
      n_healthy = c(33L, 12L),
      n_qc_excluded = c(0L, 1L),
      design_rank = c(6L, 13L),
      design_columns = c(6L, 13L),
      residual_df = c(116L, 60L),
      design_leverage_sqrt = c(0.22229046974842828, 0.34363953005869768),
      stringsAsFactors = FALSE
    ),
    model_formula = "~ case_status + age_months + sex + technical_batch",
    coefficient = "case_statusRSV",
    direction = "RSV_minus_healthy",
    mu_grid = c(0, 0.10, 0.20, 0.30, 0.50),
    tau_grid = c(0, 0.05, 0.10, 0.20),
    residual_sd_grid = c(0.5, 1.0, 1.5),
    primary_residual_sd = 1.0,
    n_sim = 10000L,
    alpha = 0.05,
    probability_thresholds = c(0.80, 0.90),
    seed = 20260812L,
    rng_kind = "L'Ecuyer-CMRG",
    loop_order = c("residual_sd", "tau", "mu", "replicate"),
    meta_method = "REML",
    meta_test = "knha",
    meta_k = 2L,
    hksj_df = 1L,
    probability_denominator = 10000L,
    null_direction_reporting = "NA_with_denominator_10000",
    instability_caveat = "k=2; tau2 estimates, HKSJ discovery probabilities, and thresholds may be unstable; conditional design sensitivity only",
    nonmonotone_probability_policy = "report_raw_grid_no_isotonic_or_curve_smoothing",
    threshold_rule = "first_nonzero_grid_mu_with_discovery_probability_ge_threshold_no_interpolation",
    not_reached_label = "not_reached_within_prespecified_grid"
  )
}

task12_sha256 <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 failed for ", path, call. = FALSE)
  sub("[[:space:]].*$", "", out[[1L]])
}

task12_read_tsv <- function(path) {
  read.delim(
    path, sep = "\t", quote = "", check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8"
  )
}

task12_verify_safe_inputs <- function(spec = task12_design_spec()) {
  allowed <- c(spec$manifest_path, unname(spec$qc_paths))
  if (!all(file.exists(allowed))) stop("Task 12 safe design input is absent", call. = FALSE)
  observed <- c(
    sample_manifest_frozen = task12_sha256(spec$manifest_path),
    GSE105450_qc = task12_sha256(spec$qc_paths[["GSE105450"]]),
    GSE103842_qc = task12_sha256(spec$qc_paths[["GSE103842"]])
  )
  if (!identical(observed, spec$input_sha256)) {
    stop("Task 12 safe design input SHA-256 mismatch", call. = FALSE)
  }
  invisible(observed)
}

task12_build_cohort_design <- function(cohort, spec = task12_design_spec()) {
  if (!cohort %in% spec$cohorts) stop("Unregistered Task 12 cohort", call. = FALSE)
  task12_verify_safe_inputs(spec)
  manifest <- task12_read_tsv(spec$manifest_path)
  qc <- task12_read_tsv(spec$qc_paths[[cohort]])
  required_manifest <- c(
    "series_accession", "sample_id", "subject_id", "include", "meta_eligible_confirmatory",
    "status", "age_months", "sex", "batch"
  )
  required_qc <- c("sample_id", "exclude_qc")
  if (!all(required_manifest %in% names(manifest)) || !all(required_qc %in% names(qc))) {
    stop("Task 12 design input schema is incomplete", call. = FALSE)
  }
  if (anyDuplicated(manifest$sample_id) || anyDuplicated(qc$sample_id)) {
    stop("Task 12 design input has duplicate sample IDs", call. = FALSE)
  }
  if (anyNA(qc$exclude_qc)) stop("Task 12 QC decision cannot be missing", call. = FALSE)
  retained_ids <- as.character(qc$sample_id[!as.logical(qc$exclude_qc)])
  eligible <- manifest[
    manifest$series_accession == cohort & manifest$include &
      manifest$meta_eligible_confirmatory,
    , drop = FALSE
  ]
  if (!setequal(as.character(eligible$sample_id), as.character(qc$sample_id))) {
    stop("Task 12 QC rows do not equal frozen confirmatory eligibility", call. = FALSE)
  }
  data <- eligible[match(retained_ids, eligible$sample_id), , drop = FALSE]
  if (anyNA(data$sample_id) || !identical(as.character(data$sample_id), retained_ids)) {
    stop("Task 12 retained sample join is not one-to-one", call. = FALSE)
  }
  data$case_status <- factor(
    ifelse(data$status == "RSV_case", "RSV", ifelse(data$status == "healthy_control", "healthy", NA_character_)),
    levels = c("healthy", "RSV")
  )
  data$age_months <- suppressWarnings(as.numeric(data$age_months))
  data$sex <- factor(data$sex, levels = c("female", "male"))
  batch_levels <- sort(unique(as.character(data$batch)), method = "radix")
  data$technical_batch <- factor(as.character(data$batch), levels = batch_levels)
  if (anyNA(data$case_status) || any(!is.finite(data$age_months)) || anyNA(data$sex) || anyNA(data$technical_batch)) {
    stop("Task 12 requires complete frozen covariates", call. = FALSE)
  }
  design <- stats::model.matrix(
    stats::as.formula(spec$model_formula), data = data
  )
  if (qr(design)$rank != ncol(design) || !spec$coefficient %in% colnames(design)) {
    stop("Task 12 frozen design is not full rank or lacks the RSV coefficient", call. = FALSE)
  }
  coefficient_index <- match(spec$coefficient, colnames(design))
  leverage <- sqrt(solve(crossprod(design))[coefficient_index, coefficient_index])
  expected <- spec$expected[spec$expected$cohort == cohort, , drop = FALSE]
  observed <- c(
    n = nrow(data),
    n_RSV = sum(data$case_status == "RSV"),
    n_healthy = sum(data$case_status == "healthy"),
    n_qc_excluded = sum(as.logical(qc$exclude_qc)),
    design_rank = qr(design)$rank,
    design_columns = ncol(design),
    residual_df = nrow(design) - qr(design)$rank
  )
  expected_integer <- unlist(expected[c(
    "n", "n_RSV", "n_healthy", "n_qc_excluded", "design_rank", "design_columns", "residual_df"
  )], use.names = TRUE)
  if (!identical(as.integer(observed), as.integer(expected_integer)) ||
      !isTRUE(all.equal(leverage, expected$design_leverage_sqrt, tolerance = 1e-14))) {
    stop("Task 12 frozen sample/design contract changed", call. = FALSE)
  }
  list(
    cohort = cohort,
    sample_data = data,
    design = design,
    coefficient = spec$coefficient,
    coefficient_index = coefficient_index,
    design_leverage_sqrt = unname(leverage),
    residual_df = unname(observed[["residual_df"]])
  )
}

task12_fixed_cohort_se <- function(residual_sd, designs) {
  if (length(residual_sd) != 1L || !is.finite(residual_sd) || residual_sd <= 0) {
    stop("Task 12 residual SD must be one finite positive value", call. = FALSE)
  }
  leverage <- vapply(designs, `[[`, numeric(1), "design_leverage_sqrt")
  out <- residual_sd * leverage
  if (any(!is.finite(out) | out <= 0)) stop("Task 12 fixed cohort SE is invalid", call. = FALSE)
  out
}

task12_generate_one <- function(mu, tau, cohort_se) {
  if (length(mu) != 1L || !is.finite(mu) || mu < 0 ||
      length(tau) != 1L || !is.finite(tau) || tau < 0 ||
      length(cohort_se) != 2L || any(!is.finite(cohort_se) | cohort_se <= 0)) {
    stop("Task 12 DGP input is invalid", call. = FALSE)
  }
  random_effect <- stats::rnorm(2L, mean = 0, sd = tau)
  true_cohort_effect <- mu + random_effect
  observed_effect <- stats::rnorm(2L, mean = true_cohort_effect, sd = cohort_se)
  data.frame(
    cohort_index = 1:2,
    random_effect = random_effect,
    true_cohort_effect = true_cohort_effect,
    observed_effect = observed_effect,
    SE = cohort_se,
    stringsAsFactors = FALSE
  )
}

task12_fit_registered_meta <- function(observed_effect, cohort_se) {
  if (!requireNamespace("metafor", quietly = TRUE)) stop("Frozen metafor package is unavailable", call. = FALSE)
  if (length(observed_effect) != 2L || any(!is.finite(observed_effect)) ||
      length(cohort_se) != 2L || any(!is.finite(cohort_se) | cohort_se <= 0)) {
    stop("Task 12 meta input must contain two finite effects and positive SEs", call. = FALSE)
  }
  fit <- metafor::rma.uni(
    yi = observed_effect, sei = cohort_se,
    method = "REML", test = "knha"
  )
  values <- c(
    k = 2,
    hksj_df = 1,
    estimate = as.numeric(fit$b[[1L]]),
    meta_se = as.numeric(fit$se),
    ci_lb = as.numeric(fit$ci.lb),
    ci_ub = as.numeric(fit$ci.ub),
    p_value = as.numeric(fit$pval),
    tau2_hat = as.numeric(fit$tau2)
  )
  if (any(!is.finite(values)) || values[["meta_se"]] <= 0 ||
      values[["p_value"]] < 0 || values[["p_value"]] > 1) {
    stop("Task 12 registered meta fit is invalid", call. = FALSE)
  }
  values
}

task12_evaluate_one <- function(meta_values, mu, alpha = task12_design_spec()$alpha) {
  required <- c("estimate", "ci_lb", "ci_ub", "p_value")
  if (!all(required %in% names(meta_values)) || any(!is.finite(meta_values[required])) ||
      length(mu) != 1L || !is.finite(mu) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("Task 12 evaluation input is invalid", call. = FALSE)
  }
  discovered <- unname(meta_values[["p_value"]] < alpha)
  direction_is_defined <- mu > 0
  c(
    discovered = discovered,
    direction_correct_discovery = if (direction_is_defined) discovered && meta_values[["estimate"]] > 0 else NA,
    wrong_direction_discovery = if (direction_is_defined) discovered && meta_values[["estimate"]] < 0 else NA,
    mean_effect_ci_covered = meta_values[["ci_lb"]] <= mu && mu <= meta_values[["ci_ub"]]
  )
}

task12_binomial_summary <- function(successes, n = 10000L, conf_level = 0.95) {
  if (length(successes) != 1L || !is.finite(successes) || successes < 0 ||
      successes != floor(successes) || length(n) != 1L || !is.finite(n) ||
      n != 10000L || n != floor(n) || successes > n ||
      length(conf_level) != 1L || !is.finite(conf_level) || conf_level != 0.95) {
    stop("Task 12 binomial summary requires integer successes in [0,10000], n=10000, and conf_level=0.95", call. = FALSE)
  }
  p <- successes / n
  mcse <- sqrt(p * (1 - p) / n)
  z <- stats::qnorm(0.975)
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
  c(
    successes = successes,
    denominator = n,
    probability = p,
    monte_carlo_se = mcse,
    wilson_lb = max(0, center - half),
    wilson_ub = min(1, center + half),
    conf_level = conf_level,
    z = z
  )
}

task12_summarize_scenario <- function(evaluations, mu, n = 10000L) {
  spec <- task12_design_spec()
  required <- c(
    "discovered", "direction_correct_discovery", "wrong_direction_discovery",
    "mean_effect_ci_covered"
  )
  if (!is.data.frame(evaluations) || !all(required %in% names(evaluations)) ||
      nrow(evaluations) != n || n != 10000L || length(mu) != 1L || !is.finite(mu) || mu < 0) {
    stop("Task 12 scenario summary requires exactly 10000 registered evaluations", call. = FALSE)
  }
  if (anyNA(evaluations$discovered) || anyNA(evaluations$mean_effect_ci_covered)) {
    stop("Task 12 discovery and CI coverage cannot be missing", call. = FALSE)
  }
  if (mu == 0) {
    if (!all(is.na(evaluations$direction_correct_discovery)) ||
        !all(is.na(evaluations$wrong_direction_discovery))) {
      stop("Task 12 null direction metrics must be NA", call. = FALSE)
    }
  } else if (anyNA(evaluations$direction_correct_discovery) ||
             anyNA(evaluations$wrong_direction_discovery)) {
    stop("Task 12 non-null direction metrics cannot be missing", call. = FALSE)
  }
  metrics <- required
  summaries <- lapply(metrics, function(metric) {
    values <- evaluations[[metric]]
    if (all(is.na(values))) {
      out <- c(
        successes = NA_real_, denominator = n, probability = NA_real_,
        monte_carlo_se = NA_real_, wilson_lb = NA_real_, wilson_ub = NA_real_,
        conf_level = 0.95, z = stats::qnorm(0.975)
      )
    } else {
      if (!is.logical(values) || anyNA(values)) stop("Task 12 metric must be complete logical", call. = FALSE)
      out <- task12_binomial_summary(sum(values), n = n, conf_level = 0.95)
    }
    data.frame(
      mu = mu,
      metric = metric,
      interpretation = if (metric == "discovered" && mu == 0) "type_I_error" else metric,
      k = spec$meta_k,
      hksj_df = spec$hksj_df,
      meta_method = spec$meta_method,
      meta_test = spec$meta_test,
      instability_caveat = spec$instability_caveat,
      t(out), check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, summaries)
  rownames(result) <- NULL
  if (any(result$denominator != 10000L)) stop("Task 12 probability denominator changed", call. = FALSE)
  result
}

task12_minimum_grid_effect <- function(mu_grid, discovery_probability, threshold) {
  if (length(mu_grid) != length(discovery_probability) || any(!is.finite(mu_grid)) ||
      any(!is.finite(discovery_probability) | discovery_probability < 0 | discovery_probability > 1) ||
      length(threshold) != 1L || !is.finite(threshold) || threshold <= 0 || threshold >= 1) {
    stop("Task 12 threshold input is invalid", call. = FALSE)
  }
  if (is.unsorted(mu_grid, strictly = TRUE) || mu_grid[[1L]] != 0) {
    stop("Task 12 mu grid must be strictly increasing from zero", call. = FALSE)
  }
  eligible <- which(mu_grid > 0 & discovery_probability >= threshold)
  if (!length(eligible)) return(NA_real_)
  mu_grid[[eligible[[1L]]]]
}
