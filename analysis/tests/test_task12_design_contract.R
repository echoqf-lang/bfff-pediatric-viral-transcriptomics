#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

amendment_path <- file.path(
  "docs", "science-superpowers", "preregistrations", "amendments",
  "2026-08-05-task12-fixed-n-design-sensitivity.md"
)
contract_path <- file.path("analysis", "R", "task12_design_contract.R")
clarification_path <- file.path(
  "docs", "science-superpowers", "preregistrations", "amendments",
  "2026-08-05-task12-binomial-reporting-clarification.md"
)
amendment <- paste(readLines(amendment_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
clarification <- paste(readLines(clarification_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
contract_text <- paste(readLines(contract_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
source(contract_path)

required_amendment_text <- c(
  "Status:** FROZEN",
  "fixed-N design-sensitivity assessment",
  "GSE105450：122 人，RSV 89、健康 33",
  "GSE103842：73 人，RSV 61、健康 12",
  "~ case_status + age_months + sex + technical_batch",
  "SE_i = sigma * h_i",
  "sigma={0.5, 1.0, 1.5}",
  "mu` 固定为 log2FC `{0, 0.10, 0.20, 0.30, 0.50}",
  "tau` 固定为 `{0, 0.05, 0.10, 0.20}",
  "10,000 次",
  "method=\"REML\", test=\"knha\"",
  "不模拟 camera 竞争性集合检验",
  "mean_effect_CI_coverage",
  "not_reached_within_prespecified_grid",
  "不报一个伪 MDE",
  "不标注事后功效",
  "不能救回或改写已经失败的 H1 联合判定"
)
stopifnot(all(vapply(required_amendment_text, grepl, logical(1), x = amendment, fixed = TRUE)))
required_clarification_text <- c(
  "分母全部固定为 10,000",
  "I 类错误率",
  "概率固定报告为 `NA`",
  "z=qnorm(0.975)",
  "无连续性校正",
  "`x=0`、`x=5000` 和 `x=10000`",
  "`k=2` 和 `hksj_df=1`",
  "τ²估计、HKSJ 发现率和阈值可能不稳定",
  "不得因为 k=2 的结果不稳定而切换至 Wald、固定效应",
  "等调/曲线平滑"
)
stopifnot(all(vapply(required_clarification_text, grepl, logical(1), x = clarification, fixed = TRUE)))

# The executable design contract has exactly three outcome-free input paths.
spec <- task12_design_spec()
stopifnot(
  identical(spec$cohorts, c("GSE105450", "GSE103842")),
  identical(spec$mu_grid, c(0, 0.10, 0.20, 0.30, 0.50)),
  identical(spec$tau_grid, c(0, 0.05, 0.10, 0.20)),
  identical(spec$residual_sd_grid, c(0.5, 1.0, 1.5)),
  identical(spec$primary_residual_sd, 1.0),
  identical(spec$n_sim, 10000L),
  identical(spec$alpha, 0.05),
  identical(spec$probability_thresholds, c(0.80, 0.90)),
  identical(spec$seed, 20260812L),
  identical(spec$rng_kind, "L'Ecuyer-CMRG"),
  identical(spec$loop_order, c("residual_sd", "tau", "mu", "replicate")),
  identical(spec$meta_method, "REML"),
  identical(spec$meta_test, "knha"),
  identical(spec$meta_k, 2L),
  identical(spec$hksj_df, 1L),
  identical(spec$probability_denominator, 10000L),
  identical(spec$null_direction_reporting, "NA_with_denominator_10000"),
  grepl("k=2", spec$instability_caveat, fixed = TRUE),
  grepl("may be unstable", spec$instability_caveat, fixed = TRUE),
  identical(spec$nonmonotone_probability_policy, "report_raw_grid_no_isotonic_or_curve_smoothing")
)
task12_verify_safe_inputs(spec)
stopifnot(
  grepl("sample_manifest_frozen.tsv", contract_text, fixed = TRUE),
  grepl("GSE105450_qc.tsv", contract_text, fixed = TRUE),
  grepl("GSE103842_qc.tsv", contract_text, fixed = TRUE),
  !grepl("gene_effects.tsv", contract_text, fixed = TRUE),
  !grepl("target_set_camera.tsv", contract_text, fixed = TRUE),
  !grepl("results/meta", contract_text, fixed = TRUE),
  !grepl("results/sensitivity", contract_text, fixed = TRUE)
)

# Frozen samples, complete cases, design ranks, and coefficient leverage must
# be reconstructed from the manifest and blind-QC tables without expression.
designs <- setNames(lapply(spec$cohorts, task12_build_cohort_design, spec = spec), spec$cohorts)
for (cohort in spec$cohorts) {
  one <- designs[[cohort]]
  expected <- spec$expected[spec$expected$cohort == cohort, , drop = FALSE]
  stopifnot(
    nrow(one$sample_data) == expected$n,
    sum(one$sample_data$case_status == "RSV") == expected$n_RSV,
    sum(one$sample_data$case_status == "healthy") == expected$n_healthy,
    qr(one$design)$rank == expected$design_rank,
    ncol(one$design) == expected$design_columns,
    one$residual_df == expected$residual_df,
    identical(one$coefficient, "case_statusRSV"),
    isTRUE(all.equal(one$design_leverage_sqrt, expected$design_leverage_sqrt, tolerance = 1e-14))
  )
}

# Conditional cohort SEs are design leverage times the pre-set residual SD;
# no empirical target-gene SE is accepted by this API.
se_half <- task12_fixed_cohort_se(0.5, designs)
se_one <- task12_fixed_cohort_se(1.0, designs)
se_one_half <- task12_fixed_cohort_se(1.5, designs)
stopifnot(
  identical(names(se_one), spec$cohorts),
  isTRUE(all.equal(unname(se_one), spec$expected$design_leverage_sqrt, tolerance = 1e-14)),
  isTRUE(all.equal(unname(se_half), unname(0.5 * se_one), tolerance = 0)),
  isTRUE(all.equal(unname(se_one_half), unname(1.5 * se_one), tolerance = 0))
)

# DGP behavior is tested on one synthetic fixture only; this is not the formal
# 10,000-replicate simulation and writes no result artifact.
RNGkind(spec$rng_kind)
set.seed(spec$seed)
actual_draw <- task12_generate_one(mu = 0.2, tau = 0.1, cohort_se = se_one)
RNGkind(spec$rng_kind)
set.seed(spec$seed)
expected_u <- rnorm(2L, 0, 0.1)
expected_theta <- 0.2 + expected_u
expected_y <- rnorm(2L, expected_theta, se_one)
stopifnot(
  identical(actual_draw$random_effect, expected_u),
  identical(actual_draw$true_cohort_effect, expected_theta),
  identical(actual_draw$observed_effect, expected_y),
  identical(actual_draw$SE, unname(se_one))
)

# The wrapper must be exactly equivalent to the registered two-study
# REML/HKSJ fit, and the decision rule must remain two-sided alpha 0.05.
fixture_y <- c(0.15, 0.28)
wrapped <- task12_fit_registered_meta(fixture_y, se_one)
direct <- metafor::rma.uni(yi = fixture_y, sei = se_one, method = "REML", test = "knha")
expected_meta <- c(
  k = 2,
  hksj_df = 1,
  estimate = as.numeric(direct$b[[1L]]),
  meta_se = as.numeric(direct$se),
  ci_lb = as.numeric(direct$ci.lb),
  ci_ub = as.numeric(direct$ci.ub),
  p_value = as.numeric(direct$pval),
  tau2_hat = as.numeric(direct$tau2)
)
stopifnot(isTRUE(all.equal(wrapped, expected_meta, tolerance = 0)))
evaluated <- task12_evaluate_one(wrapped, mu = 0.2, alpha = spec$alpha)
stopifnot(
  identical(unname(evaluated[["discovered"]]), wrapped[["p_value"]] < 0.05),
  identical(
    unname(evaluated[["mean_effect_ci_covered"]]),
    wrapped[["ci_lb"]] <= 0.2 && 0.2 <= wrapped[["ci_ub"]]
  )
)

# The binomial reporting contract uses p=x/10000, its binomial MCSE, and the
# uncorrected standard Wilson interval with z=qnorm(.975), including boundaries.
z <- qnorm(0.975)
expected_wilson <- function(x, n = 10000) {
  p <- x / n
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
  c(
    successes = x, denominator = n, probability = p,
    monte_carlo_se = sqrt(p * (1 - p) / n),
    wilson_lb = max(0, center - half), wilson_ub = min(1, center + half),
    conf_level = 0.95, z = z
  )
}
for (x in c(0, 5000, 10000)) {
  stopifnot(isTRUE(all.equal(
    task12_binomial_summary(x, n = 10000, conf_level = 0.95),
    expected_wilson(x), tolerance = 0
  )))
}
stopifnot(
  task12_binomial_summary(0)[["wilson_lb"]] == 0,
  task12_binomial_summary(10000)[["wilson_ub"]] == 1,
  task12_binomial_summary(5000)[["monte_carlo_se"]] == 0.005
)

# All four metrics retain denominator 10000. At mu=0 discovery is explicitly
# type-I error and the two direction metrics are fixed NA, not chosen post hoc.
null_evaluations <- data.frame(
  discovered = rep(c(TRUE, FALSE), c(500L, 9500L)),
  direction_correct_discovery = rep(NA, 10000L),
  wrong_direction_discovery = rep(NA, 10000L),
  mean_effect_ci_covered = rep(c(TRUE, FALSE), c(9500L, 500L))
)
null_summary <- task12_summarize_scenario(null_evaluations, mu = 0)
stopifnot(
  identical(null_summary$metric, c(
    "discovered", "direction_correct_discovery", "wrong_direction_discovery",
    "mean_effect_ci_covered"
  )),
  all(null_summary$denominator == 10000L),
  all(null_summary$k == 2L),
  all(null_summary$hksj_df == 1L),
  all(null_summary$meta_method == "REML"),
  all(null_summary$meta_test == "knha"),
  all(grepl("may be unstable", null_summary$instability_caveat, fixed = TRUE)),
  null_summary$interpretation[[1L]] == "type_I_error",
  null_summary$successes[[1L]] == 500,
  null_summary$probability[[1L]] == 0.05,
  all(is.na(null_summary$probability[2:3])),
  null_summary$probability[[4L]] == 0.95
)
positive_evaluations <- data.frame(
  discovered = rep(c(TRUE, FALSE), c(5000L, 5000L)),
  direction_correct_discovery = rep(c(TRUE, FALSE), c(4500L, 5500L)),
  wrong_direction_discovery = c(rep(FALSE, 4500L), rep(TRUE, 500L), rep(FALSE, 5000L)),
  mean_effect_ci_covered = rep(c(TRUE, FALSE), c(9400L, 600L))
)
positive_summary <- task12_summarize_scenario(positive_evaluations, mu = 0.2)
stopifnot(
  all(positive_summary$denominator == 10000L),
  all(positive_summary$k == 2L),
  all(positive_summary$hksj_df == 1L),
  identical(positive_summary$successes, c(5000, 4500, 500, 9400)),
  identical(positive_summary$probability, c(0.50, 0.45, 0.05, 0.94))
)

# Null per-replicate direction fields are also frozen to NA at evaluation time.
null_evaluated <- task12_evaluate_one(wrapped, mu = 0, alpha = spec$alpha)
stopifnot(
  is.na(null_evaluated[["direction_correct_discovery"]]),
  is.na(null_evaluated[["wrong_direction_discovery"]])
)

# The only meta implementation remains k=2 REML/HKSJ with df=1; no fallback
# estimator or smoothing implementation is present in the executable contract.
meta_body <- paste(deparse(body(task12_fit_registered_meta)), collapse = " ")
stopifnot(
  grepl("method = \"REML\"", meta_body, fixed = TRUE),
  grepl("test = \"knha\"", meta_body, fixed = TRUE),
  !grepl("test = \"z\"", meta_body, fixed = TRUE),
  !grepl("method = \"FE\"", meta_body, fixed = TRUE),
  !grepl("isoreg", contract_text, fixed = TRUE),
  !grepl("smooth.spline", contract_text, fixed = TRUE)
)

# Threshold selection is discrete, has no interpolation, and explicitly
# returns NA when the pre-set grid does not reach the requested probability.
stopifnot(
  identical(task12_minimum_grid_effect(spec$mu_grid, c(0.05, 0.22, 0.81, 0.91, 0.99), 0.80), 0.20),
  identical(task12_minimum_grid_effect(spec$mu_grid, c(0.05, 0.22, 0.79, 0.89, 0.95), 0.90), 0.50),
  is.na(task12_minimum_grid_effect(spec$mu_grid, c(0.05, 0.22, 0.40, 0.60, 0.79), 0.80))
)

cat("Task 12 fixed-N design contract passed; no target effect/result file read and no formal simulation run\n")
