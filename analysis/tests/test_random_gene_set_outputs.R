#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 2)

read_tsv_contract <- function(path) {
  if (!file.exists(path)) stop("Missing Task 10 output: ", path, call. = FALSE)
  read.delim(path, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
}

calibration <- read_tsv_contract(file.path("results", "meta", "random_set_empirical_calibration.tsv"))
null <- read_tsv_contract(file.path("results", "meta", "random_set_null_statistics.tsv"))
audit <- read_tsv_contract(file.path("results", "meta", "random_set_matching_audit.tsv"))
plan <- read_tsv_contract(file.path("results", "meta", "random_set_relaxation_plan.tsv"))
equivalence <- read_tsv_contract(file.path("results", "meta", "random_set_camera_equivalence.tsv"))
diagnostics <- read_tsv_contract(file.path("results", "meta", "random_set_null_diagnostics.tsv"))
extreme_check <- read_tsv_contract(file.path("results", "meta", "random_set_extreme_direct_camera_check.tsv"))
checksum <- read_tsv_contract(file.path("logs", "checksums", "task10_sha256.tsv"))

if (nrow(calibration) != 1L || calibration$random_seed != 20260804L ||
    calibration$random_iterations != 10000L) stop("Calibration header contract failed", call. = FALSE)
if (nrow(null) != 10000L || !identical(null$iteration, seq_len(10000L))) {
  stop("Null distribution must contain exactly iterations 1..10000", call. = FALSE)
}
if (any(!is.finite(null$combined_stouffer_z)) || any(null$GSE105450_p_value <= 0 | null$GSE105450_p_value > 1) ||
    any(null$GSE103842_p_value <= 0 | null$GSE103842_p_value > 1)) {
  stop("Null statistics contain invalid values", call. = FALSE)
}
expected_p <- (1 + sum(abs(null$combined_stouffer_z) >= abs(calibration$observed_combined_stouffer_z))) / 10001
if (!isTRUE(all.equal(calibration$empirical_two_sided_p, expected_p, tolerance = 1e-15))) {
  stop("Empirical P does not recompute from all 10,000 null statistics", call. = FALSE)
}
if (!identical(as.character(calibration$primary_decision), "H1_not_supported_camera_component_failed") ||
    calibration$h1_supported || calibration$random_calibration_can_rescue) {
  stop("Task 10 changed the frozen H1 decision", call. = FALSE)
}
if (!setequal(audit$cohort, c("GSE105450", "GSE103842")) ||
    !identical(audit$detectable_targets[match(c("GSE105450", "GSE103842"), audit$cohort)], c(391L, 374L)) ||
    any(!audit$every_set_size_valid) || any(audit$any_within_set_duplicate) ||
    any(audit$any_target_contamination) || any(audit$expression_bin_relaxed) ||
    any(audit$forbidden_matching_fields_used)) {
  stop("Matching audit contract failed", call. = FALSE)
}
if (sum(plan$n_allocate[plan$cohort == "GSE105450"]) != 391L ||
    sum(plan$n_allocate[plan$cohort == "GSE103842"]) != 374L ||
    any(plan$target_expression_bin != plan$source_expression_bin)) {
  stop("Relaxation plan contract failed", call. = FALSE)
}
if (nrow(equivalence) < 14L || any(table(equivalence$cohort) < 7L) ||
    any(equivalence$status != "passed_exact_camera_contract") ||
    any(equivalence$abs_p_difference > equivalence$tolerance) ||
    any(equivalence$abs_correlation_difference > equivalence$tolerance) ||
    any(equivalence$abs_signed_z_difference > equivalence$tolerance) ||
    any(!equivalence$direction_match)) {
  stop("Cached/direct camera equivalence contract failed", call. = FALSE)
}
if (nrow(diagnostics) != 1L || diagnostics$GSE105450_unique_set_hashes != 10000L ||
    diagnostics$GSE103842_unique_set_hashes != 10000L || diagnostics$unique_paired_set_hashes != 10000L ||
    diagnostics$max_abs_stouffer_formula_difference > 1e-12 ||
    diagnostics$max_abs_combined_z >= abs(calibration$observed_combined_stouffer_z)) {
  stop("Random-set diagnostic contract failed", call. = FALSE)
}
if (nrow(extreme_check) != 20L || length(unique(extreme_check$iteration)) != 10L ||
    any(table(extreme_check$cohort) != 10L) || any(!extreme_check$direction_match) ||
    any(extreme_check$abs_p_difference > 1e-12) ||
    any(extreme_check$abs_correlation_difference > 1e-12) ||
    any(extreme_check$abs_signed_z_difference > 1e-12) ||
    any(extreme_check$abs_combined_difference > 1e-12)) {
  stop("Extreme-set direct camera verification failed", call. = FALSE)
}
figure <- file.path("results", "figures", "random_set_null_distribution.pdf")
if (!file.exists(figure) || file.info(figure)$size < 1000) stop("Null PDF is absent or empty", call. = FALSE)
if (!all(file.exists(checksum$artifact)) || any(vapply(seq_len(nrow(checksum)), function(i) {
  observed <- system2("shasum", c("-a", "256", checksum$artifact[[i]]), stdout = TRUE)
  !identical(sub("[[:space:]].*$", "", observed[[1L]]), checksum$sha256[[i]])
}, logical(1)))) stop("Task 10 checksum manifest failed", call. = FALSE)

cat("Task 10 output contracts: PASS\n")
