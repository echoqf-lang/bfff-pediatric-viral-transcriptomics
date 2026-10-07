#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(TASK12_SKIP_MAIN = "1")
source(file.path("analysis", "R", "10_power_precision_simulation.R"), local = FALSE)

spec <- task12_design_spec()
designs <- setNames(
  lapply(spec$cohorts, task12_build_cohort_design, spec = spec),
  spec$cohorts
)
cohort_se <- unname(task12_fixed_cohort_se(1, designs))

# Vectorized generation must preserve the registered replicate-level RNG order:
# u1, u2, sampling error 1, sampling error 2 for every replicate.
RNGkind(spec$rng_kind)
set.seed(101L)
actual <- task12_generate_scenario(mu = 0.2, tau = 0.1, cohort_se = cohort_se, n = 7L)
RNGkind(spec$rng_kind)
set.seed(101L)
expected <- do.call(rbind, lapply(seq_len(7L), function(i) {
  task12_generate_one(mu = 0.2, tau = 0.1, cohort_se = cohort_se)
}))
stopifnot(
  isTRUE(all.equal(actual$observed_effect_1, expected$observed_effect[expected$cohort_index == 1L], tolerance = 1e-15)),
  isTRUE(all.equal(actual$observed_effect_2, expected$observed_effect[expected$cohort_index == 2L], tolerance = 1e-15))
)

# The optimized k=2 calculation must be numerically equivalent to direct
# metafor REML/HKSJ, including p-value decisions and HKSJ confidence limits.
RNGkind(spec$rng_kind)
set.seed(808L)
y <- matrix(rnorm(80L), ncol = 2L, byrow = TRUE)
fast <- task12_fast_meta_two_study(y, cohort_se)
for (i in seq_len(nrow(y))) {
  direct <- task12_fit_registered_meta(y[i, ], cohort_se)
  stopifnot(
    isTRUE(all.equal(fast$tau2_hat[[i]], unname(direct[["tau2_hat"]]), tolerance = 1e-12)),
    isTRUE(all.equal(fast$estimate[[i]], unname(direct[["estimate"]]), tolerance = 1e-12)),
    isTRUE(all.equal(fast$meta_se[[i]], unname(direct[["meta_se"]]), tolerance = 1e-12)),
    isTRUE(all.equal(fast$ci_lb[[i]], unname(direct[["ci_lb"]]), tolerance = 1e-12)),
    isTRUE(all.equal(fast$ci_ub[[i]], unname(direct[["ci_ub"]]), tolerance = 1e-12)),
    isTRUE(all.equal(fast$p_value[[i]], unname(direct[["p_value"]]), tolerance = 1e-12)),
    identical(fast$p_value[[i]] < spec$alpha, direct[["p_value"]] < spec$alpha)
  )
}

# A scenario summary has exactly four registered binomial outcomes, fixed
# denominator 10,000, and undefined direction metrics under mu=0.
RNGkind(spec$rng_kind)
set.seed(303L)
null_draw <- task12_generate_scenario(0, 0.1, cohort_se, spec$n_sim)
null_meta <- task12_fast_meta_two_study(
  cbind(null_draw$observed_effect_1, null_draw$observed_effect_2), cohort_se
)
null_summary <- task12_summarize_meta_scenario(null_meta, mu = 0, spec = spec)
stopifnot(
  nrow(null_summary) == 4L,
  all(null_summary$denominator == 10000L),
  all(is.na(null_summary$probability[null_summary$metric %in% c(
    "direction_correct_discovery", "wrong_direction_discovery"
  )])),
  null_summary$interpretation[null_summary$metric == "discovered"] == "type_I_error"
)

cat("Task 12 power/precision simulation unit tests passed\n")
