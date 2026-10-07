#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
source(file.path("analysis", "R", "task12_design_contract.R"), local = FALSE)

read_tsv <- function(path) read.delim(
  path, sep = "\t", quote = "", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)
sha256 <- function(path) sub(
  "[[:space:]].*$", "",
  system2("shasum", c("-a", "256", path), stdout = TRUE)[[1L]]
)

spec <- task12_design_spec()
paths <- c(
  sensitivity = file.path("results", "power", "fixed_n_design_sensitivity.tsv"),
  threshold = file.path("results", "power", "fixed_n_threshold_summary.tsv"),
  audit = file.path("results", "power", "fixed_n_simulation_audit.tsv"),
  precision = file.path("results", "power", "observed_precision_summary.tsv"),
  figure = file.path("results", "figures", "fixed_n_design_sensitivity_by_tau.pdf"),
  simulation_lock = file.path("logs", "checksums", "task12_simulation_lock_sha256.tsv"),
  final_lock = file.path("logs", "checksums", "task12_sha256.tsv"),
  commands = file.path("logs", "session_info", "task12_commands.log"),
  session = file.path("logs", "session_info", "task12_session_info.txt")
)
stopifnot(all(file.exists(paths)), !file.exists(file.path("results", "power", "mde_simulation.tsv")))

sensitivity <- read_tsv(paths[["sensitivity"]])
stopifnot(
  nrow(sensitivity) == 60L * 4L,
  setequal(unique(sensitivity$residual_sd), spec$residual_sd_grid),
  setequal(unique(sensitivity$tau), spec$tau_grid),
  setequal(unique(sensitivity$mu), spec$mu_grid),
  setequal(unique(sensitivity$metric), c(
    "discovered", "direction_correct_discovery", "wrong_direction_discovery",
    "mean_effect_ci_covered"
  )),
  all(table(sensitivity$residual_sd, sensitivity$tau, sensitivity$mu) == 4L),
  all(sensitivity$denominator == spec$n_sim),
  all(sensitivity$k == 2L),
  all(sensitivity$hksj_df == 1L),
  all(sensitivity$meta_method == "REML"),
  all(sensitivity$meta_test == "knha"),
  all(grepl("k=2", sensitivity$instability_caveat, fixed = TRUE)),
  all(sensitivity$probability[!is.na(sensitivity$probability)] >= 0),
  all(sensitivity$probability[!is.na(sensitivity$probability)] <= 1),
  all(sensitivity$wilson_lb[!is.na(sensitivity$wilson_lb)] <= sensitivity$probability[!is.na(sensitivity$probability)]),
  all(sensitivity$probability[!is.na(sensitivity$probability)] <= sensitivity$wilson_ub[!is.na(sensitivity$wilson_ub)])
)
null_direction <- sensitivity$mu == 0 & sensitivity$metric %in% c(
  "direction_correct_discovery", "wrong_direction_discovery"
)
stopifnot(
  all(is.na(sensitivity$successes[null_direction])),
  all(is.na(sensitivity$probability[null_direction])),
  all(sensitivity$denominator[null_direction] == 10000L),
  all(sensitivity$interpretation[sensitivity$mu == 0 & sensitivity$metric == "discovered"] == "type_I_error")
)
defined <- !is.na(sensitivity$probability)
stopifnot(isTRUE(all.equal(
  sensitivity$monte_carlo_se[defined],
  sqrt(sensitivity$probability[defined] * (1 - sensitivity$probability[defined]) / 10000),
  tolerance = 1e-15
)))

threshold <- read_tsv(paths[["threshold"]])
stopifnot(
  nrow(threshold) == 3L * 4L * 2L,
  setequal(threshold$threshold, spec$probability_thresholds),
  all(table(threshold$residual_sd, threshold$tau) == 2L),
  all(threshold$selection_rule == spec$threshold_rule),
  all(threshold$nonmonotone_probability_policy == spec$nonmonotone_probability_policy),
  all(threshold$status %in% c("reached_on_prespecified_grid", spec$not_reached_label)),
  all(is.na(threshold$minimum_grid_mu[threshold$status == spec$not_reached_label]))
)

audit <- read_tsv(paths[["audit"]])
stopifnot(
  nrow(audit) == 60L,
  identical(audit$scenario_index, seq_len(60L)),
  all(audit$n_sim == 10000L),
  all(audit$direct_validation_n == 5L),
  all(audit$direct_decision_mismatches == 0L),
  all(audit$max_abs_difference <= 1e-10),
  all(audit$rng_kind == spec$rng_kind),
  all(audit$seed == spec$seed),
  all(audit$loop_order == paste(spec$loop_order, collapse = ">")),
  all(audit$parallel == FALSE),
  all(audit$failed_fits == 0L)
)

precision <- read_tsv(paths[["precision"]])
stopifnot(
  nrow(precision) == 1L,
  precision$scope == "all_estimable_primary_targets_no_p_filtering",
  precision$n_targets > 0,
  precision$ci_width_min <= precision$ci_width_q1,
  precision$ci_width_q1 <= precision$ci_width_median,
  precision$ci_width_median <= precision$ci_width_q3,
  precision$ci_width_q3 <= precision$ci_width_max,
  precision$prediction_width_q1 <= precision$prediction_width_median,
  precision$prediction_width_median <= precision$prediction_width_q3,
  precision$prediction_width_q3 <= precision$prediction_width_max,
  precision$cohort_direction_agreement_fraction >= 0,
  precision$cohort_direction_agreement_fraction <= 1,
  precision$post_hoc_power_computed == FALSE,
  precision$simulation_parameters_updated_from_observed_results == FALSE,
  precision$primary_H1_rescue_allowed == FALSE
)

sim_lock <- read_tsv(paths[["simulation_lock"]])
stopifnot(
  all(file.exists(sim_lock$artifact)),
  identical(unname(vapply(sim_lock$artifact, sha256, character(1))), unname(sim_lock$sha256)),
  all(!grepl("target_gene_meta|gene_effects|observed_precision", sim_lock$artifact))
)
final_lock <- read_tsv(paths[["final_lock"]])
stopifnot(
  all(file.exists(final_lock$artifact)),
  identical(unname(vapply(final_lock$artifact, sha256, character(1))), unname(final_lock$sha256))
)

commands <- paste(readLines(paths[["commands"]], warn = FALSE), collapse = "\n")
stopifnot(
  grepl("simulation locked before observed precision read", commands, fixed = TRUE),
  grepl("not post-hoc power", commands, fixed = TRUE),
  grepl("cannot rescue", commands, fixed = TRUE),
  grepl("old MDE filenames superseded before execution", commands, fixed = TRUE)
)

cat("Task 12 power/precision simulation output tests passed\n")
