#!/usr/bin/env Rscript

candidate_paths <- c(
  file.path(
    "manuscript", "v6_bmc_submission_assets", "supplementary_tables",
    "Supplementary_Table_1_complete_blood_effects_meta.tsv"
  ),
  file.path(
    "results", "supplementary_tables",
    "Supplementary_Table_1_complete_blood_effects_meta.tsv"
  )
)
path <- candidate_paths[file.exists(candidate_paths)][1]

stopifnot(file.exists(path))
x <- read.delim(path, check.names = FALSE)

required <- c(
  "meta_estimate",
  "meta_conventional_hksj_se",
  "meta_conventional_hksj_ci_low",
  "meta_conventional_hksj_ci_high",
  "meta_conventional_hksj_p_value",
  "meta_conventional_hksj_fdr_bh",
  "meta_modified_hksj_se",
  "meta_modified_hksj_ci_low",
  "meta_modified_hksj_ci_high",
  "meta_modified_hksj_p_value",
  "meta_modified_hksj_fdr_bh",
  "meta_tau2",
  "meta_i2",
  "meta_direction_pattern"
)

stopifnot(nrow(x) == 86L)
stopifnot(ncol(x) == 37L)
stopifnot(all(required %in% names(x)))
stopifnot(sum(is.finite(x$meta_conventional_hksj_p_value)) == 81L)
stopifnot(sum(is.finite(x$meta_modified_hksj_p_value)) == 81L)
stopifnot(sum(x$meta_conventional_hksj_fdr_bh < 0.05, na.rm = TRUE) == 5L)
stopifnot(sum(x$meta_modified_hksj_fdr_bh < 0.05, na.rm = TRUE) == 0L)
stopifnot(all(
  x$meta_modified_hksj_ci_low <= x$meta_estimate |
    is.na(x$meta_modified_hksj_ci_low)
))
stopifnot(all(
  x$meta_modified_hksj_ci_high >= x$meta_estimate |
    is.na(x$meta_modified_hksj_ci_high)
))

message("V6 Supplementary Table 1 dual-HKSJ reporting contract passed")
