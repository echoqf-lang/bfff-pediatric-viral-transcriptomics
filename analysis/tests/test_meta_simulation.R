#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

pipeline_path <- file.path("analysis", "R", "statistical_pipeline.R")
if (!file.exists(pipeline_path)) stop("Task 6 statistical pipeline is absent", call. = FALSE)
source(pipeline_path, local = FALSE)

required_packages <- c("limma", "metafor")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing frozen package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)

expect_error <- function(code, pattern) {
  message <- tryCatch({
    force(code)
    NA_character_
  }, error = function(e) conditionMessage(e))
  if (is.na(message) || !grepl(pattern, message, perl = TRUE)) {
    stop("Expected error matching '", pattern, "' but got: ", message, call. = FALSE)
  }
  invisible(message)
}

make_sample_data <- function(cohort, n_case = 60L, n_control = 30L) {
  n <- n_case + n_control
  status <- c(rep("healthy", n_control), rep("RSV", n_case))
  data.frame(
    sample_id = sprintf("%s_S%03d", cohort, seq_len(n)),
    subject_id = sprintf("%s_P%03d", cohort, seq_len(n)),
    cohort = cohort,
    platform = paste0("SIM_PLATFORM_", cohort),
    case_status = status,
    age_months = pmax(0, rnorm(n, mean = ifelse(status == "RSV", 10, 9), sd = 3)),
    sex = sample(c("female", "male"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

simulate_three_cohorts <- function(seed = 20260804L) {
  set.seed(seed)
  gene_ids <- sprintf("SIM%05d", seq_len(5000L))
  targets <- gene_ids[seq_len(100L)]
  injected <- targets[seq_len(30L)]
  age_beta <- rnorm(length(gene_ids), 0, 0.012)
  sex_beta <- rnorm(length(gene_ids), 0, 0.08)
  baseline <- rnorm(length(gene_ids), 7, 1.2)
  platform_scales <- c(C1 = 0.80, C2 = 1.00, C3 = 1.20)
  cohorts <- lapply(names(platform_scales), function(cohort) {
    sample_data <- make_sample_data(cohort)
    n <- nrow(sample_data)
    true_mean_effect <- as.numeric(gene_ids %in% injected) * 0.4
    cohort_effect <- true_mean_effect + rnorm(length(gene_ids), 0, 0.10)
    centered_age <- sample_data$age_months - mean(sample_data$age_months)
    male <- as.numeric(sample_data$sex == "male")
    case <- as.numeric(sample_data$case_status == "RSV")
    platform_scale <- platform_scales[[cohort]]
    platform_baseline <- 7 + platform_scale * (baseline - 7)
    expression <- matrix(
      platform_baseline +
        platform_scale * outer(age_beta, centered_age) +
        platform_scale * outer(sex_beta, male) +
        outer(cohort_effect, case) +
        matrix(rnorm(length(gene_ids) * n, 0, 0.65 * platform_scale), nrow = length(gene_ids)),
      nrow = length(gene_ids),
      dimnames = list(gene_ids, sample_data$sample_id)
    )
    list(
      cohort = cohort,
      expression = expression,
      sample_data = sample_data,
      target_ids = targets,
      injected_ids = injected,
      true_effect = setNames(true_mean_effect, gene_ids),
      true_cohort_effect = setNames(cohort_effect, gene_ids),
      simulated_tau = 0.10,
      platform_scale = platform_scale
    )
  })
  names(cohorts) <- names(platform_scales)
  cohorts
}

# Small known-effect gate: the production coefficient must be RSV minus healthy.
set.seed(17L)
small_samples <- make_sample_data("KNOWN", n_case = 6L, n_control = 6L)
small_expression <- matrix(rnorm(20L * 12L, 0, 0.20), nrow = 20L,
                           dimnames = list(sprintf("K%02d", seq_len(20L)), small_samples$sample_id))
small_expression[1L, small_samples$case_status == "RSV"] <-
  small_expression[1L, small_samples$case_status == "RSV"] + 1
small_fit <- fit_limma_cohort(
  expression = small_expression,
  sample_data = small_samples,
  cohort_id = "KNOWN",
  direction = "RSV_minus_healthy"
)
stopifnot(small_fit$gene_effects$log2FC[small_fit$gene_effects$gene_id == "K01"] > 0.5)

# Deliberate boundary violations must fail before any model fit.
duplicate_samples <- small_samples
duplicate_samples$subject_id[2L] <- duplicate_samples$subject_id[1L]
expect_error(
  fit_limma_cohort(small_expression, duplicate_samples, "KNOWN", "RSV_minus_healthy"),
  "Duplicate subject_id"
)
expect_error(
  fit_limma_cohort(small_expression, small_samples, "KNOWN", "healthy_minus_RSV"),
  "Direction must be RSV_minus_healthy"
)
pooled_samples <- small_samples
pooled_samples$platform[seq_len(6L)] <- "SIM_PLATFORM_OTHER"
expect_error(
  fit_limma_cohort(small_expression, pooled_samples, "KNOWN", "RSV_minus_healthy"),
  "single platform"
)
rank_deficient <- small_samples
rank_deficient$age_months <- as.numeric(rank_deficient$case_status == "RSV")
expect_error(
  fit_limma_cohort(small_expression, rank_deficient, "KNOWN", "RSV_minus_healthy"),
  "full rank"
)

# camera interpretability is determined from the complete frozen set before camera runs.
make_frozen_targets <- function(n_detectable, n_frozen) {
  c(rownames(small_expression)[seq_len(n_detectable)],
    sprintf("ABSENT%03d", seq_len(n_frozen - n_detectable)))
}
original_camera_runner <- run_camera_analysis
camera_call_count <- 0L
run_camera_analysis <- function(...) {
  camera_call_count <<- camera_call_count + 1L
  original_camera_runner(...)
}
camera_2_of_100 <- test_camera_target_set(
  small_expression, small_samples, make_frozen_targets(2L, 100L), "KNOWN", "RSV_minus_healthy"
)
camera_9_of_100 <- test_camera_target_set(
  small_expression, small_samples, make_frozen_targets(9L, 100L), "KNOWN", "RSV_minus_healthy"
)
camera_10_of_100 <- test_camera_target_set(
  small_expression, small_samples, make_frozen_targets(10L, 100L), "KNOWN", "RSV_minus_healthy"
)
camera_10_of_20 <- test_camera_target_set(
  small_expression, small_samples, make_frozen_targets(10L, 20L), "KNOWN", "RSV_minus_healthy"
)
run_camera_analysis <- original_camera_runner
camera_call_boundary_gate <- camera_call_count == 1L
camera_2_of_100_gate <- identical(camera_2_of_100$status, "uninterpretable") &&
  camera_2_of_100$n_detectable == 2L && all(is.na(camera_2_of_100[, c("direction", "p_value", "signed_z")]))
camera_9_of_100_gate <- identical(camera_9_of_100$status, "uninterpretable") &&
  camera_9_of_100$n_detectable == 9L && all(is.na(camera_9_of_100[, c("direction", "p_value", "signed_z")]))
camera_10_of_100_gate <- identical(camera_10_of_100$status, "uninterpretable") &&
  camera_10_of_100$n_detectable == 10L && camera_10_of_100$coverage == 0.10 &&
  all(is.na(camera_10_of_100[, c("direction", "p_value", "signed_z")]))
camera_10_of_20_gate <- identical(camera_10_of_20$status, "interpretable") &&
  camera_10_of_20$n_detectable == 10L && camera_10_of_20$coverage == 0.50 &&
  camera_10_of_20$direction %in% c("Up", "Down") &&
  is.finite(camera_10_of_20$p_value) && camera_10_of_20$p_value > 0 && camera_10_of_20$p_value <= 1
stopifnot(
  all(vapply(list(camera_2_of_100, camera_9_of_100, camera_10_of_100),
             function(x) identical(x$status, "uninterpretable"), logical(1))),
  all(vapply(list(camera_2_of_100, camera_9_of_100, camera_10_of_100),
             function(x) all(is.na(x[, c("direction", "p_value", "signed_z")])), logical(1))),
  camera_2_of_100$n_detectable == 2L,
  camera_9_of_100$n_detectable == 9L,
  camera_10_of_100$n_detectable == 10L,
  camera_10_of_100$coverage == 0.10,
  camera_10_of_20$status == "interpretable",
  camera_10_of_20$n_detectable == 10L,
  camera_10_of_20$coverage == 0.50,
  camera_2_of_100_gate,
  camera_9_of_100_gate,
  camera_10_of_100_gate,
  camera_10_of_20_gate,
  camera_call_boundary_gate
)
expect_error(
  test_camera_target_set(small_expression, small_samples, c("K01", "K01"), "KNOWN", "RSV_minus_healthy"),
  "unique"
)
expect_error(
  test_camera_target_set(small_expression, small_samples, c("K01", ""), "KNOWN", "RSV_minus_healthy"),
  "non-empty"
)

cohorts <- simulate_three_cohorts(seed = 20260804L)
cohort_fits <- lapply(cohorts, function(x) {
  fit <- fit_limma_cohort(
    expression = x$expression,
    sample_data = x$sample_data,
    cohort_id = x$cohort,
    direction = "RSV_minus_healthy"
  )
  camera <- test_camera_target_set(
    expression = x$expression,
    sample_data = x$sample_data,
    frozen_targets = x$target_ids,
    cohort_id = x$cohort,
    direction = "RSV_minus_healthy"
  )
  list(fit = fit, camera = camera)
})

stopifnot(
  length(cohort_fits) == 3L,
  all(vapply(cohort_fits, function(x) x$fit$design_rank == x$fit$design_columns, logical(1))),
  all(vapply(cohort_fits, function(x) all(is.finite(x$fit$gene_effects$SE) & x$fit$gene_effects$SE > 0), logical(1))),
  all(vapply(cohort_fits, function(x) identical(x$fit$coefficient, "case_statusRSV"), logical(1)))
)

target_ids <- cohorts[[1L]]$target_ids
effects <- do.call(rbind, lapply(names(cohort_fits), function(cohort) {
  out <- cohort_fits[[cohort]]$fit$gene_effects
  out[out$gene_id %in% target_ids, c("gene_id", "cohort", "direction", "log2FC", "SE"), drop = FALSE]
}))
expected_cohorts <- c("C1", "C2", "C3")
meta <- meta_analyze_gene_effects(effects, expected_cohorts = expected_cohorts)

# Meta completeness and direction gates must reject every malformed cohort set.
expect_error(meta_analyze_gene_effects(effects, expected_cohorts = character()), "non-empty")
expect_error(meta_analyze_gene_effects(effects, expected_cohorts = c("C1", "C1")), "unique")
missing_one <- effects[-which(effects$gene_id == target_ids[[1L]] & effects$cohort == "C3"), , drop = FALSE]
expect_error(meta_analyze_gene_effects(missing_one, expected_cohorts), "cohort set")
extra_one <- rbind(effects, transform(effects[1L, , drop = FALSE], cohort = "C4"))
expect_error(meta_analyze_gene_effects(extra_one, expected_cohorts), "cohort set")
duplicate_one <- rbind(effects, effects[1L, , drop = FALSE])
expect_error(meta_analyze_gene_effects(duplicate_one, expected_cohorts), "duplicate gene-cohort")
wrong_direction <- effects
wrong_direction$direction[1L] <- "healthy_minus_RSV"
expect_error(meta_analyze_gene_effects(wrong_direction, expected_cohorts), "RSV_minus_healthy")
empty_gene <- effects
empty_gene$gene_id[1L] <- ""
expect_error(meta_analyze_gene_effects(empty_gene, expected_cohorts), "non-empty gene_id")
empty_cohort <- effects
empty_cohort$cohort[1L] <- ""
expect_error(meta_analyze_gene_effects(empty_cohort, expected_cohorts), "non-empty cohort")

injected_ids <- cohorts[[1L]]$injected_ids
injected_median <- median(meta$estimate[meta$gene_id %in% injected_ids])
null_median <- median(meta$estimate[!meta$gene_id %in% injected_ids])
camera_rows <- do.call(rbind, lapply(cohort_fits, `[[`, "camera"))

# Direct one-gene reference fit proves the wrapper calls the registered estimator.
reference_gene <- target_ids[[1L]]
reference_input <- effects[effects$gene_id == reference_gene, , drop = FALSE]
reference_fit <- metafor::rma.uni(
  yi = reference_input$log2FC,
  sei = reference_input$SE,
  method = "REML",
  test = "knha"
)
reference_prediction <- stats::predict(reference_fit, level = 95)
wrapper_reference <- meta[meta$gene_id == reference_gene, , drop = FALSE]
meta_exact_gate <- isTRUE(all.equal(
  unname(unlist(wrapper_reference[c(
    "estimate", "se", "ci_lb", "ci_ub", "prediction_lb", "prediction_ub",
    "p_value", "tau2", "Q", "Q_p_value", "I2"
  )])),
  unname(c(
    reference_fit$b[[1L]], reference_fit$se, reference_fit$ci.lb, reference_fit$ci.ub,
    reference_prediction$pi.lb, reference_prediction$pi.ub, reference_fit$pval,
    reference_fit$tau2, reference_fit$QE, reference_fit$QEp, reference_fit$I2
  )),
  tolerance = 1e-12
)) && all(meta$method == "metafor::rma.uni(method=REML,test=knha)") && all(meta$k == 3L)
platform_gate <- length(unique(vapply(cohorts, `[[`, numeric(1), "platform_scale"))) == 3L &&
  all(vapply(cohort_fits, function(x) x$fit$n_samples == 90L, logical(1)))
dgp_deviations <- unlist(lapply(cohorts, function(x) x$true_cohort_effect - x$true_effect), use.names = FALSE)
dgp_empirical_tau_sd <- stats::sd(dgp_deviations)
dgp_tau_gate <- dgp_empirical_tau_sd >= 0.095 && dgp_empirical_tau_sd <= 0.105
pi_ordered_gate <- all(meta$prediction_lb <= meta$estimate & meta$estimate <= meta$prediction_ub)

stopifnot(
  nrow(meta) == 100L,
  all(meta$k == 3L),
  all(is.finite(meta$estimate)),
  all(is.finite(meta$se) & meta$se > 0),
  all(is.finite(meta$ci_lb)),
  all(is.finite(meta$ci_ub)),
  all(is.finite(meta$prediction_lb)),
  all(is.finite(meta$prediction_ub)),
  all(is.finite(camera_rows$inter_gene_correlation)),
  all(camera_rows$status == "interpretable"),
  all(camera_rows$n_frozen == 100L),
  all(camera_rows$n_detectable == 100L),
  all(camera_rows$coverage == 1),
  injected_median >= 0.30,
  injected_median <= 0.50,
  abs(null_median) < 0.05,
  all(camera_rows$direction == "Up"),
  all(camera_rows$p_value < 0.05),
  meta_exact_gate,
  platform_gate,
  dgp_tau_gate,
  pi_ordered_gate
)

recovery <- rbind(
  data.frame(
    artifact_type = "recovery_gate", cohort = "META_3_COHORT",
    metric = c(
      "injected_gene_meta_effect_median", "null_gene_meta_effect_median",
      "null_gene_meta_effect_median_absolute", "dgp_cohort_effect_deviation_empirical_sd"
    ),
    value = c(injected_median, null_median, abs(null_median), dgp_empirical_tau_sd),
    criterion = c("0.30_to_0.50", "descriptive_signed_median", "less_than_0.05", "0.095_to_0.105"),
    pass = c(
      injected_median >= 0.30 && injected_median <= 0.50, TRUE,
      abs(null_median) < 0.05, dgp_tau_gate
    ),
    note = "Simulation-only pipeline validation; not a biological result.",
    stringsAsFactors = FALSE
  ),
  data.frame(
    artifact_type = "camera_gate", cohort = camera_rows$cohort,
    metric = "camera_signed_z",
    value = camera_rows$signed_z,
    criterion = "direction_Up_and_two_sided_P_less_than_0.05",
    pass = camera_rows$direction == "Up" & camera_rows$p_value < 0.05,
    note = paste0(
      "camera Direction is relative enrichment for the RSV-minus-healthy contrast; ",
      "signed_z=sign(Direction)*qnorm(1-P/2), derived from camera's two-sided P value."
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    artifact_type = "implementation_gate", cohort = "SIMULATION_ONLY",
    metric = c(
      "design_matrix_full_rank", "moderated_se_finite_positive",
      "duplicate_subject_rejected", "reverse_direction_rejected",
      "cross_platform_pooling_rejected", "rank_deficient_design_rejected",
      "meta_exact_reml_knha", "platforms_modeled_separately",
      "prediction_intervals_ordered_and_contain_estimate",
      "camera_2_of_100_uninterpretable_without_test",
      "camera_9_of_100_uninterpretable_without_test",
      "camera_10_of_100_low_coverage_uninterpretable_without_test",
      "camera_10_of_20_boundary_interpretable",
      "camera_not_called_below_interpretability_threshold"
    ),
    value = c(
      rep(1, 6L), as.numeric(meta_exact_gate), as.numeric(platform_gate), as.numeric(pi_ordered_gate),
      as.numeric(camera_2_of_100_gate), as.numeric(camera_9_of_100_gate),
      as.numeric(camera_10_of_100_gate), as.numeric(camera_10_of_20_gate),
      as.numeric(camera_call_boundary_gate)
    ),
    criterion = rep("must_equal_1", 14L),
    pass = c(
      rep(TRUE, 6L), meta_exact_gate, platform_gate, pi_ordered_gate,
      camera_2_of_100_gate, camera_9_of_100_gate, camera_10_of_100_gate, camera_10_of_20_gate,
      camera_call_boundary_gate
    ),
    note = "Simulation-only implementation gate; not a biological result.",
    stringsAsFactors = FALSE
  )
)

output_path <- file.path("results", "power", "pipeline_recovery_test.tsv")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.table(recovery, output_path, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

session_path <- file.path("logs", "session_info", "task6_simulation_session_info.txt")
dir.create(dirname(session_path), recursive = TRUE, showWarnings = FALSE)
session_lines <- c(
  "Task 6 simulation-only statistical pipeline validation",
  "Seed: 20260804",
  "Simulation: 3 cohorts; each 60 RSV / 30 healthy; 5000 genes; 100 frozen targets; 30 target log2FC=0.4; tau=0.10",
  sprintf("DGP cohort-effect deviation empirical SD across 15000 gene-cohort draws: %.9f", dgp_empirical_tau_sd),
  "camera interpretability gate: n_detectable>=10 and coverage>=0.50; otherwise no camera call and inferential fields are NA",
  "Meta completeness gate: every gene must have exactly expected cohorts C1, C2, and C3",
  "Boundaries: no real sample manifest case_status; no real expression RDS; no biological inference",
  capture.output(sessionInfo())
)
writeLines(session_lines, session_path, useBytes = TRUE)

commands_path <- file.path("logs", "session_info", "task6_commands.log")
writeLines(c(
  "Task 6: simulation-only limma, camera, and REML/HKSJ pipeline validation",
  "Command: RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_CONFIG_NAMESPACES_CHECK=FALSE Rscript analysis/tests/test_meta_simulation.R",
  "Seed: 20260804",
  "No network access; no downloads; no real phenotype-bearing manifest or real expression object read",
  "Each cohort is modeled separately; mixed-platform input is a hard failure",
  "Meta estimator: metafor::rma.uni(yi=log2FC, sei=SE, method=REML, test=knha)",
  "Meta validation: one gene compared field-by-field with direct metafor; all genes require method label and k=3",
  "camera interpretability: complete unique frozen set; n_detectable>=10 and coverage>=0.50; failure returns status=uninterpretable without camera",
  "camera boundary tests: 2/100, 9/100, and 10/100 uninterpretable; 10/20 interpretable",
  "camera signed_z: sign(Direction) * qnorm(1 - two-sided P/2); Up is positive RSV-minus-healthy enrichment",
  "All outputs are simulation-only validation artifacts and are not biological results"
), commands_path, useBytes = TRUE)

sha256_file <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 failed for ", path, call. = FALSE)
  sub("[[:space:]].*$", "", output[[1L]])
}
checksum_artifacts <- c(
  pipeline_path,
  file.path("analysis", "tests", "test_meta_simulation.R"),
  output_path,
  session_path,
  commands_path
)
checksum <- data.frame(
  artifact = checksum_artifacts,
  role = c("input", "input", "output", "output", "output"),
  bytes = as.numeric(file.info(checksum_artifacts)$size),
  sha256 = vapply(checksum_artifacts, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
checksum_path <- file.path("logs", "checksums", "task6_sha256.tsv")
dir.create(dirname(checksum_path), recursive = TRUE, showWarnings = FALSE)
write.table(checksum, checksum_path, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

cat(sprintf("injected median meta effect: %.6f\n", injected_median))
cat(sprintf("null signed median meta effect: %.6f (absolute %.6f)\n", null_median, abs(null_median)))
print(camera_rows[, c("cohort", "n_genes", "direction", "p_value", "signed_z")], row.names = FALSE)
cat("task6 simulation gates passed\n")
