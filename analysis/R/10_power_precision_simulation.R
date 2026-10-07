#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

source(file.path("analysis", "R", "task12_design_contract.R"), local = FALSE)

task12_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA", fileEncoding = "UTF-8"
  )
}

task12_generate_scenario <- function(mu, tau, cohort_se, n) {
  if (length(n) != 1L || !is.finite(n) || n != floor(n) || n <= 0L) {
    stop("Task 12 scenario size must be one positive integer", call. = FALSE)
  }
  if (length(mu) != 1L || !is.finite(mu) || mu < 0 ||
      length(tau) != 1L || !is.finite(tau) || tau < 0 ||
      length(cohort_se) != 2L || any(!is.finite(cohort_se) | cohort_se <= 0)) {
    stop("Task 12 scenario DGP input is invalid", call. = FALSE)
  }
  # Each row consumes RNG in the exact registered replicate order:
  # u1, u2, sampling error 1, sampling error 2.
  z <- matrix(stats::rnorm(4L * n), nrow = n, ncol = 4L, byrow = TRUE)
  data.frame(
    replicate = seq_len(n),
    observed_effect_1 = mu + tau * z[, 1L] + cohort_se[[1L]] * z[, 3L],
    observed_effect_2 = mu + tau * z[, 2L] + cohort_se[[2L]] * z[, 4L],
    stringsAsFactors = FALSE
  )
}

task12_fast_meta_two_study <- function(observed_effect, cohort_se) {
  if (!is.matrix(observed_effect) || ncol(observed_effect) != 2L ||
      !is.numeric(observed_effect) || any(!is.finite(observed_effect)) ||
      length(cohort_se) != 2L || any(!is.finite(cohort_se) | cohort_se <= 0)) {
    stop("Task 12 fast meta requires an n by 2 finite matrix and two positive SEs", call. = FALSE)
  }
  y1 <- observed_effect[, 1L]
  y2 <- observed_effect[, 2L]
  v <- cohort_se^2
  difference <- y1 - y2

  # For an intercept-only two-study normal-normal model, the REML estimate is
  # max(0, (d^2-v1-v2)/2). The following HKSJ calculation is algebraically
  # identical to metafor::rma.uni(method="REML", test="knha") for k=2.
  tau2_hat <- pmax(0, (difference^2 - sum(v)) / 2)
  w1 <- 1 / (v[[1L]] + tau2_hat)
  w2 <- 1 / (v[[2L]] + tau2_hat)
  weight_sum <- w1 + w2
  estimate <- (w1 * y1 + w2 * y2) / weight_sum
  q_hksj <- w1 * (y1 - estimate)^2 + w2 * (y2 - estimate)^2
  meta_se <- sqrt(q_hksj / weight_sum)
  if (any(!is.finite(meta_se) | meta_se <= 0)) {
    stop("Task 12 optimized HKSJ fit produced a non-positive or non-finite SE", call. = FALSE)
  }
  statistic <- estimate / meta_se
  p_value <- 2 * stats::pt(abs(statistic), df = 1, lower.tail = FALSE)
  critical <- stats::qt(0.975, df = 1)
  ci_lb <- estimate - critical * meta_se
  ci_ub <- estimate + critical * meta_se
  result <- data.frame(
    k = 2L,
    hksj_df = 1L,
    estimate = estimate,
    meta_se = meta_se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    p_value = p_value,
    tau2_hat = tau2_hat,
    stringsAsFactors = FALSE
  )
  numeric_fields <- c("estimate", "meta_se", "ci_lb", "ci_ub", "p_value", "tau2_hat")
  if (any(!is.finite(as.matrix(result[numeric_fields]))) ||
      any(result$p_value < 0 | result$p_value > 1)) {
    stop("Task 12 optimized meta fit is invalid", call. = FALSE)
  }
  result
}

task12_summarize_meta_scenario <- function(meta, mu, spec = task12_design_spec()) {
  required <- c("estimate", "ci_lb", "ci_ub", "p_value")
  if (!is.data.frame(meta) || !all(required %in% names(meta)) ||
      nrow(meta) != spec$n_sim || any(!is.finite(as.matrix(meta[required])))) {
    stop("Task 12 scenario meta output is incomplete", call. = FALSE)
  }
  discovered <- meta$p_value < spec$alpha
  evaluations <- data.frame(
    discovered = discovered,
    direction_correct_discovery = if (mu > 0) discovered & meta$estimate > 0 else rep(NA, spec$n_sim),
    wrong_direction_discovery = if (mu > 0) discovered & meta$estimate < 0 else rep(NA, spec$n_sim),
    mean_effect_ci_covered = meta$ci_lb <= mu & mu <= meta$ci_ub,
    stringsAsFactors = FALSE
  )
  task12_summarize_scenario(evaluations, mu = mu, n = spec$n_sim)
}

task12_direct_equivalence_audit <- function(meta, observed_effect, cohort_se,
                                            indices = c(1L, 2L, 17L, 499L, 9999L)) {
  if (any(indices < 1L | indices > nrow(observed_effect))) {
    stop("Task 12 direct-equivalence indices are outside the scenario", call. = FALSE)
  }
  fields <- c("estimate", "meta_se", "ci_lb", "ci_ub", "p_value", "tau2_hat")
  differences <- numeric(0)
  decision_mismatches <- 0L
  for (i in indices) {
    direct <- task12_fit_registered_meta(observed_effect[i, ], cohort_se)
    fast_values <- unlist(meta[i, fields, drop = FALSE], use.names = TRUE)
    differences <- c(differences, abs(fast_values - direct[fields]))
    decision_mismatches <- decision_mismatches +
      as.integer((meta$p_value[[i]] < 0.05) != (direct[["p_value"]] < 0.05))
  }
  list(
    n = length(indices),
    max_abs_difference = max(differences),
    decision_mismatches = decision_mismatches
  )
}

task12_make_threshold_summary <- function(sensitivity, spec = task12_design_spec()) {
  discovery <- sensitivity[sensitivity$metric == "discovered", , drop = FALSE]
  rows <- list()
  index <- 0L
  for (residual_sd in spec$residual_sd_grid) {
    for (tau in spec$tau_grid) {
      one <- discovery[
        discovery$residual_sd == residual_sd & discovery$tau == tau,
        , drop = FALSE
      ]
      one <- one[match(spec$mu_grid, one$mu), , drop = FALSE]
      if (anyNA(one$mu) || !identical(as.numeric(one$mu), spec$mu_grid)) {
        stop("Task 12 discovery grid is incomplete", call. = FALSE)
      }
      for (threshold in spec$probability_thresholds) {
        index <- index + 1L
        minimum <- task12_minimum_grid_effect(spec$mu_grid, one$probability, threshold)
        rows[[index]] <- data.frame(
          residual_sd = residual_sd,
          tau = tau,
          threshold = threshold,
          minimum_grid_mu = minimum,
          status = if (is.na(minimum)) spec$not_reached_label else "reached_on_prespecified_grid",
          selection_rule = spec$threshold_rule,
          nonmonotone_probability_policy = spec$nonmonotone_probability_policy,
          k = spec$meta_k,
          hksj_df = spec$hksj_df,
          instability_caveat = spec$instability_caveat,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

task12_plot_design_sensitivity <- function(sensitivity, path, spec = task12_design_spec()) {
  discovery <- sensitivity[sensitivity$metric == "discovered", , drop = FALSE]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(
    path, width = 10.2, height = 5.4, useDingbats = FALSE,
    timestamp = FALSE, bg = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(
    mfrow = c(1, 3), mar = c(4.2, 4.2, 3.2, 1.0), oma = c(3.0, 0, 0, 0),
    family = "sans"
  )
  on.exit(graphics::par(old), add = TRUE)
  colors <- c("#1B9E77", "#D95F02", "#7570B3", "#404040")
  symbols <- c(16, 17, 15, 18)
  for (panel in seq_along(spec$residual_sd_grid)) {
    residual_sd <- spec$residual_sd_grid[[panel]]
    one_sd <- discovery[discovery$residual_sd == residual_sd, , drop = FALSE]
    graphics::plot(
      NA_real_, NA_real_, xlim = range(spec$mu_grid), ylim = c(0, 1),
      xlab = expression("Assumed mean effect " * mu * " (log"[2] * "FC)"),
      ylab = if (panel == 1L) "Two-sided discovery probability" else "",
      main = paste0("Residual SD = ", format(residual_sd, nsmall = 1)),
      xaxt = "n", yaxs = "i"
    )
    graphics::axis(1, at = spec$mu_grid, labels = format(spec$mu_grid, nsmall = 1))
    graphics::abline(h = c(0.8, 0.9), col = "#BDBDBD", lty = c(2, 3))
    for (j in seq_along(spec$tau_grid)) {
      tau <- spec$tau_grid[[j]]
      one <- one_sd[one_sd$tau == tau, , drop = FALSE]
      one <- one[match(spec$mu_grid, one$mu), , drop = FALSE]
      graphics::lines(one$mu, one$probability, type = "b", lwd = 1.6,
                      pch = symbols[[j]], col = colors[[j]])
    }
    if (panel == 1L) {
      graphics::legend(
        "topleft", legend = paste0("tau = ", format(spec$tau_grid, nsmall = 2)),
        col = colors, pch = symbols, lty = 1, lwd = 1.6, bty = "n", cex = 0.82
      )
    }
  }
  graphics::mtext(
    "Raw prespecified grid; no interpolation or smoothing. k=2, HKSJ df=1; thresholds may be unstable and are conditional design-sensitivity summaries only.",
    side = 1, outer = TRUE, line = 1.2, cex = 0.75
  )
  invisible(path)
}

task12_verify_checksum_manifest <- function(path) {
  if (!file.exists(path)) stop("Task 12 simulation lock is absent", call. = FALSE)
  manifest <- task12_read_tsv(path)
  if (!all(c("artifact", "sha256") %in% names(manifest)) || !all(file.exists(manifest$artifact))) {
    stop("Task 12 simulation lock is incomplete", call. = FALSE)
  }
  observed <- vapply(manifest$artifact, task12_sha256, character(1))
  if (!identical(unname(observed), unname(manifest$sha256))) {
    stop("Task 12 simulation lock SHA-256 mismatch", call. = FALSE)
  }
  invisible(manifest)
}

task12_run_simulation <- function(spec = task12_design_spec()) {
  task12_verify_safe_inputs(spec)
  designs <- setNames(lapply(spec$cohorts, task12_build_cohort_design, spec = spec), spec$cohorts)
  RNGkind(spec$rng_kind)
  set.seed(spec$seed)
  scenario_rows <- list()
  audit_rows <- list()
  scenario_index <- 0L
  for (residual_sd in spec$residual_sd_grid) {
    cohort_se <- unname(task12_fixed_cohort_se(residual_sd, designs))
    for (tau in spec$tau_grid) {
      for (mu in spec$mu_grid) {
        scenario_index <- scenario_index + 1L
        generated <- task12_generate_scenario(mu, tau, cohort_se, spec$n_sim)
        observed <- cbind(generated$observed_effect_1, generated$observed_effect_2)
        meta <- task12_fast_meta_two_study(observed, cohort_se)
        summary <- task12_summarize_meta_scenario(meta, mu, spec)
        summary <- cbind(
          data.frame(
            scenario_index = scenario_index,
            residual_sd = residual_sd,
            tau = tau,
            cohort_1 = spec$cohorts[[1L]],
            cohort_1_se = cohort_se[[1L]],
            cohort_2 = spec$cohorts[[2L]],
            cohort_2_se = cohort_se[[2L]],
            stringsAsFactors = FALSE
          ),
          summary
        )
        scenario_rows[[scenario_index]] <- summary
        equivalence <- task12_direct_equivalence_audit(meta, observed, cohort_se)
        audit_rows[[scenario_index]] <- data.frame(
          scenario_index = scenario_index,
          residual_sd = residual_sd,
          tau = tau,
          mu = mu,
          n_sim = spec$n_sim,
          seed = spec$seed,
          rng_kind = spec$rng_kind,
          loop_order = paste(spec$loop_order, collapse = ">"),
          parallel = FALSE,
          failed_fits = 0L,
          direct_validation_n = equivalence$n,
          max_abs_difference = equivalence$max_abs_difference,
          direct_decision_mismatches = equivalence$decision_mismatches,
          optimized_formula = "exact_k2_intercept_only_REML_HKSJ",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (scenario_index != 60L) stop("Task 12 did not execute exactly 60 scenarios", call. = FALSE)
  sensitivity <- do.call(rbind, scenario_rows)
  audit <- do.call(rbind, audit_rows)
  rownames(sensitivity) <- NULL
  rownames(audit) <- NULL
  if (any(audit$direct_decision_mismatches != 0L) || any(audit$max_abs_difference > 1e-10)) {
    stop("Task 12 optimized meta calculation failed direct metafor equivalence", call. = FALSE)
  }
  list(
    sensitivity = sensitivity,
    threshold = task12_make_threshold_summary(sensitivity, spec),
    audit = audit
  )
}

task12_lock_simulation <- function(outputs, spec = task12_design_spec()) {
  simulation_log <- file.path("logs", "session_info", "task12_simulation_phase.log")
  dir.create(dirname(simulation_log), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "Task 12 formal simulation phase",
    paste0("Seed: ", spec$seed, "; RNGkind: ", spec$rng_kind),
    paste0("Loop order: ", paste(spec$loop_order, collapse = " > "), "; serial; 60 scenarios x 10000 replicates"),
    "Inputs read before simulation lock: frozen sample manifest and two blind-QC tables only",
    "No target effect, SE, P-value, meta, target-set, random-set, or Task 11 result was read",
    "Estimator: exact vectorized k=2 REML/HKSJ formula, sampled against direct metafor in every scenario",
    "Simulation is fixed-N conditional design sensitivity, not prospective sample-size estimation or post-hoc power",
    "Old MDE filenames superseded before execution by the frozen Task 12 amendment"
  ), simulation_log, useBytes = TRUE)
  lock_artifacts <- c(
    file.path("analysis", "R", "10_power_precision_simulation.R"),
    file.path("analysis", "R", "task12_design_contract.R"),
    file.path("analysis", "tests", "test_task12_power_precision_simulation.R"),
    file.path("analysis", "tests", "test_task12_power_precision_outputs.R"),
    file.path("docs", "science-superpowers", "preregistrations", "amendments",
              "2026-08-05-task12-fixed-n-design-sensitivity.md"),
    file.path("docs", "science-superpowers", "preregistrations", "amendments",
              "2026-08-05-task12-binomial-reporting-clarification.md"),
    spec$manifest_path,
    unname(spec$qc_paths),
    unname(outputs),
    simulation_log
  )
  if (!all(file.exists(lock_artifacts))) stop("Task 12 simulation lock artifact is absent", call. = FALSE)
  lock <- data.frame(
    artifact = lock_artifacts,
    role = c(rep("code_or_contract", 6L), rep("outcome_free_input", 3L),
             rep("simulation_output", length(outputs)), "simulation_log"),
    bytes = as.numeric(file.info(lock_artifacts)$size),
    sha256 = vapply(lock_artifacts, task12_sha256, character(1)),
    stringsAsFactors = FALSE
  )
  lock_path <- file.path("logs", "checksums", "task12_simulation_lock_sha256.tsv")
  task12_write_tsv(lock, lock_path)
  task12_verify_checksum_manifest(lock_path)
  lock_path
}

task12_verified_task9_path <- function(path, checksum_path) {
  checksum <- task12_read_tsv(checksum_path)
  row <- checksum[checksum$artifact == path, , drop = FALSE]
  if (nrow(row) != 1L || !file.exists(path) || task12_sha256(path) != row$sha256[[1L]]) {
    stop("Task 12 observed precision input is not uniquely verified by Task 9 SHA", call. = FALSE)
  }
  path
}

task12_observed_precision <- function(simulation_lock_path) {
  # This verification is the hard phase boundary: no observed outcome result is
  # read unless the formal simulation artifacts already exist and match SHA.
  task12_verify_checksum_manifest(simulation_lock_path)
  task9_checksum_path <- file.path("logs", "checksums", "task9_sha256.tsv")
  if (!file.exists(task9_checksum_path)) stop("Task 9 checksum manifest is absent", call. = FALSE)
  meta_path <- task12_verified_task9_path(
    file.path("results", "meta", "target_gene_meta_reml_knha.tsv"), task9_checksum_path
  )
  effect_paths <- file.path(
    "results", "cohort", paste0(c("GSE105450", "GSE103842"), "_gene_effects.tsv")
  )
  effect_paths <- vapply(effect_paths, task12_verified_task9_path, character(1),
                         checksum_path = task9_checksum_path)
  meta <- task12_read_tsv(meta_path)
  required_meta <- c(
    "gene_id", "k", "ci_lb", "ci_ub", "prediction_lb", "prediction_ub",
    "method", "direction"
  )
  if (!all(required_meta %in% names(meta)) || !nrow(meta) || anyDuplicated(meta$gene_id) ||
      any(meta$k != 2L) || any(meta$method != "REML_HKSJ") ||
      any(meta$direction != "RSV_minus_healthy") ||
      any(!is.finite(as.matrix(meta[c("ci_lb", "ci_ub", "prediction_lb", "prediction_ub")]))) ||
      any(meta$ci_ub < meta$ci_lb) || any(meta$prediction_ub < meta$prediction_lb)) {
    stop("Task 12 observed precision meta input violates the complete primary-family contract", call. = FALSE)
  }
  meta_ids <- as.character(meta$gene_id)
  cohort_effects <- lapply(seq_along(effect_paths), function(i) {
    one <- task12_read_tsv(effect_paths[[i]])
    if (!all(c("gene_id", "log2FC", "direction") %in% names(one)) || anyDuplicated(one$gene_id)) {
      stop("Task 12 cohort effect input schema is invalid", call. = FALSE)
    }
    one$gene_id <- as.character(one$gene_id)
    one <- one[match(meta_ids, one$gene_id), , drop = FALSE]
    if (anyNA(one$gene_id) || !identical(one$gene_id, meta_ids) ||
        any(!is.finite(one$log2FC)) || any(one$direction != "RSV_minus_healthy")) {
      stop("Task 12 cohort effects do not cover the complete unscreened meta family", call. = FALSE)
    }
    one
  })
  ci_width <- meta$ci_ub - meta$ci_lb
  prediction_width <- meta$prediction_ub - meta$prediction_lb
  ci_q <- stats::quantile(ci_width, c(0, 0.25, 0.5, 0.75, 1), names = FALSE, type = 7)
  prediction_q <- stats::quantile(prediction_width, c(0.25, 0.5, 0.75, 1), names = FALSE, type = 7)
  agreement <- mean(sign(cohort_effects[[1L]]$log2FC) == sign(cohort_effects[[2L]]$log2FC))
  data.frame(
    scope = "all_estimable_primary_targets_no_p_filtering",
    n_targets = nrow(meta),
    ci_width_min = ci_q[[1L]],
    ci_width_q1 = ci_q[[2L]],
    ci_width_median = ci_q[[3L]],
    ci_width_q3 = ci_q[[4L]],
    ci_width_max = ci_q[[5L]],
    prediction_width_q1 = prediction_q[[1L]],
    prediction_width_median = prediction_q[[2L]],
    prediction_width_q3 = prediction_q[[3L]],
    prediction_width_max = prediction_q[[4L]],
    cohort_direction_agreement_n = sum(
      sign(cohort_effects[[1L]]$log2FC) == sign(cohort_effects[[2L]]$log2FC)
    ),
    cohort_direction_agreement_fraction = agreement,
    k = 2L,
    hksj_df = 1L,
    interpretation = "descriptive precision of complete unscreened primary target meta family",
    post_hoc_power_computed = FALSE,
    simulation_parameters_updated_from_observed_results = FALSE,
    primary_H1_rescue_allowed = FALSE,
    caveat = "k=2 prediction intervals and heterogeneity are unstable; CI precision does not diagnose power",
    stringsAsFactors = FALSE
  )
}

task12_write_final_archive <- function(paths, simulation_lock_path, spec = task12_design_spec()) {
  session_path <- file.path("logs", "session_info", "task12_session_info.txt")
  commands_path <- file.path("logs", "session_info", "task12_commands.log")
  dir.create(dirname(session_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())), session_path, useBytes = TRUE)
  writeLines(c(
    "Task 12: fixed-N conditional design sensitivity and observed precision",
    "Command: RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_CONFIG_NAMESPACES_CHECK=FALSE Rscript analysis/R/10_power_precision_simulation.R",
    paste0("60 scenarios: residual SD ", paste(spec$residual_sd_grid, collapse = ", "),
           "; tau ", paste(spec$tau_grid, collapse = ", "),
           "; mu ", paste(spec$mu_grid, collapse = ", "), "; 10000 each"),
    paste0("Seed ", spec$seed, "; ", spec$rng_kind, "; serial loop order ",
           paste(spec$loop_order, collapse = " > ")),
    "simulation locked before observed precision read",
    "Observed precision is a separate descriptive step and is not post-hoc power",
    "Task 12 cannot rescue or modify the failed primary H1 joint decision",
    "Only two cohorts: tau2, HKSJ discovery probabilities, thresholds, and prediction intervals may be unstable",
    "Raw grid reported without interpolation, extrapolation, isotonic regression, or curve smoothing",
    "old MDE filenames superseded before execution by frozen amendment; no legacy MDE outputs generated"
  ), commands_path, useBytes = TRUE)
  artifacts <- unique(c(
    file.path("analysis", "R", "10_power_precision_simulation.R"),
    file.path("analysis", "R", "task12_design_contract.R"),
    file.path("analysis", "tests", "test_task12_power_precision_simulation.R"),
    file.path("analysis", "tests", "test_task12_power_precision_outputs.R"),
    file.path("docs", "science-superpowers", "preregistrations", "amendments",
              "2026-08-05-task12-fixed-n-design-sensitivity.md"),
    file.path("docs", "science-superpowers", "preregistrations", "amendments",
              "2026-08-05-task12-binomial-reporting-clarification.md"),
    simulation_lock_path, unname(paths), session_path, commands_path,
    file.path("logs", "checksums", "task9_sha256.tsv"),
    file.path("results", "meta", "target_gene_meta_reml_knha.tsv"),
    file.path("results", "cohort", paste0(spec$cohorts, "_gene_effects.tsv"))
  ))
  if (!all(file.exists(artifacts))) stop("Task 12 final archive artifact is absent", call. = FALSE)
  checksum <- data.frame(
    artifact = artifacts,
    bytes = as.numeric(file.info(artifacts)$size),
    sha256 = vapply(artifacts, task12_sha256, character(1)),
    stringsAsFactors = FALSE
  )
  task12_write_tsv(checksum, file.path("logs", "checksums", "task12_sha256.tsv"))
}

main_task12 <- function() {
  spec <- task12_design_spec()
  if (file.exists(file.path("results", "power", "mde_simulation.tsv")) ||
      file.exists(file.path("results", "figures", "mde_by_tau.pdf"))) {
    stop("Superseded Task 12 MDE output exists; remove it only with explicit authorization", call. = FALSE)
  }
  paths <- c(
    sensitivity = file.path("results", "power", "fixed_n_design_sensitivity.tsv"),
    threshold = file.path("results", "power", "fixed_n_threshold_summary.tsv"),
    audit = file.path("results", "power", "fixed_n_simulation_audit.tsv"),
    figure = file.path("results", "figures", "fixed_n_design_sensitivity_by_tau.pdf"),
    precision = file.path("results", "power", "observed_precision_summary.tsv")
  )
  simulation <- task12_run_simulation(spec)
  task12_write_tsv(simulation$sensitivity, paths[["sensitivity"]])
  task12_write_tsv(simulation$threshold, paths[["threshold"]])
  task12_write_tsv(simulation$audit, paths[["audit"]])
  task12_plot_design_sensitivity(simulation$sensitivity, paths[["figure"]], spec)
  simulation_lock_path <- task12_lock_simulation(paths[c("sensitivity", "threshold", "audit", "figure")], spec)
  precision <- task12_observed_precision(simulation_lock_path)
  task12_write_tsv(precision, paths[["precision"]])
  task12_write_final_archive(paths, simulation_lock_path, spec)
  cat("Task 12 fixed-N design sensitivity and post-lock observed precision completed\n")
}

if (!identical(Sys.getenv("TASK12_SKIP_MAIN"), "1")) {
  main_task12()
}
