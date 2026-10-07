#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(TASK5_FUNCTIONS_ONLY = "true")
source(file.path("analysis", "R", "04_prepare_microarray.R"))

# Phenotype/outcome columns must never enter the blind QC object.
allowed <- data.frame(
  sample_id = c("S1", "S2"), platform = "GPL10558", technical_batch = c("A", "B"),
  case_status = c("case", "control"), severity = c("high", "healthy"),
  stringsAsFactors = FALSE
)
blind <- make_blind_sample_data(allowed)
stopifnot(
  identical(names(blind), c("sample_id", "platform", "technical_batch")),
  !any(c("case_status", "condition", "status", "severity") %in% names(blind))
)

# The expression component can be explicitly removed from a parsed container
# without using rm() on a list-member expression.
cleaned_container <- remove_expression_component(list(expression = matrix(1), keep = "yes"))
stopifnot(is.null(cleaned_container$expression), identical(cleaned_container$keep, "yes"))

stopifnot(identical(strip_trailing_whitespace(c("default", "value ", "value\t")), c("default", "value", "value")))

# Already log2-scale data are unchanged; positive intensity-scale data get only
# a log2 transform (never a second between-array normalization).
logged <- matrix(c(6, 7, 8, 9, 10, 11), nrow = 3)
logged_result <- transform_processed_expression(logged)
stopifnot(
  identical(logged_result$expression, logged),
  identical(logged_result$action, "already_log2_no_transform"),
  !logged_result$between_array_normalization_applied
)
intensity <- matrix(c(1, 16, 256, 2, 32, 512), nrow = 3)
intensity_result <- transform_processed_expression(intensity)
stopifnot(
  identical(intensity_result$expression, log2(pmax(intensity, 1))),
  identical(intensity_result$action, "log2_pmax1_only"),
  !intensity_result$between_array_normalization_applied
)

# A sample is excluded only after at least two independent blind metrics exceed
# their median +/- 5 MAD limits. A one-metric outlier remains included.
fixture_metrics <- data.frame(
  sample_id = paste0("S", 1:9),
  metric_a = c(rep(0, 8), 100),
  metric_b = c(rep(0, 7), 100, 100),
  metric_c = rep(0, 9),
  stringsAsFactors = FALSE
)
gate <- apply_two_metric_qc_gate(fixture_metrics, c("metric_a", "metric_b", "metric_c"))
stopifnot(
  gate$n_blind_flags[8] == 1L,
  !gate$exclude_qc[8],
  gate$n_blind_flags[9] == 2L,
  gate$exclude_qc[9]
)

# The frozen QC scale is the standard 1.4826-scaled MAD, not raw MAD.
scaled_fixture <- c(0:8, 18)
raw_center <- stats::median(scaled_fixture)
raw_mad <- stats::mad(scaled_fixture, center = raw_center, constant = 1)
stopifnot(
  abs(scaled_fixture[10] - raw_center) > 5 * raw_mad,
  !robust_flag(scaled_fixture)[10]
)

# MAD=0 is deterministic: all-equal values have no flags, while any nonzero
# deviation from the median is flagged. Non-finite inputs stop the cohort.
stopifnot(
  !any(robust_flag(rep(3, 7))),
  identical(robust_flag(c(rep(0, 6), 1)), c(rep(FALSE, 6), TRUE))
)
expect_error <- function(expression) inherits(try(force(expression), silent = TRUE), "try-error")
stopifnot(
  expect_error(robust_flag(c(1, NA_real_))),
  expect_error(apply_two_metric_qc_gate(
    data.frame(sample_id = c("S1", "S2"), a = c(1, Inf), b = c(1, 2)), c("a", "b")
  ))
)

# Only unique current Entrez mappings survive. When multiple probes map to one
# Entrez gene, the probe with highest all-sample mean is selected without labels.
probe_expression <- matrix(
  c(1, 1, 1, 5, 5, 5, 9, 9, 9, 3, 3, 3),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("P1", "P2", "P3", "P4"), paste0("S", 1:3))
)
probe_map <- data.frame(
  probe_id = c("P1", "P2", "P3", "P4"),
  entrez_id = c("10", "10", NA, "20///21"),
  official_symbol = c("G10", "G10", NA, NA),
  mapping_status = c("unique_current", "unique_current", "empty", "multi_entrez"),
  stringsAsFactors = FALSE
)
collapsed <- collapse_probes_blind(probe_expression, probe_map)
stopifnot(
  identical(rownames(collapsed$expression), "10"),
  identical(collapsed$selected_probes$probe_id, "P2"),
  identical(unname(collapsed$expression[1, ]), c(5, 5, 5))
)

# QC PDFs must be deterministic and must not leave an implicit Rplots device.
local({
  figure_dir <- tempfile("task5_figure_fixture_")
  dir.create(figure_dir)
  old_wd <- setwd(figure_dir)
  on.exit(setwd(old_wd), add = TRUE)
  fixture_qc <- data.frame(
    sample_id = paste0("S", 1:3), median_sample_correlation = c(0.8, 0.9, 0.85),
    detection_rate_p_lt_0_05 = c(0.4, 0.5, 0.6), log_array_weight = c(0, 0.1, -0.1),
    exclude_qc = c(FALSE, FALSE, TRUE), n_blind_flags = c(0L, 0L, 2L)
  )
  first_pdf <- file.path(figure_dir, "fixture_1.pdf")
  second_pdf <- file.path(figure_dir, "fixture_2.pdf")
  write_qc_figure(
    "FIXTURE", matrix(1:12, nrow = 4, dimnames = list(NULL, paste0("S", 1:3))),
    fixture_qc, matrix(1:6, nrow = 3), first_pdf
  )
  Sys.sleep(1.1)
  write_qc_figure(
    "FIXTURE", matrix(1:12, nrow = 4, dimnames = list(NULL, paste0("S", 1:3))),
    fixture_qc, matrix(1:6, nrow = 3), second_pdf
  )
  sha <- function(path) sub("[[:space:]].*$", "", system2("shasum", c("-a", "256", path), stdout = TRUE)[[1L]])
  stopifnot(
    file.exists(first_pdf), file.exists(second_pdf),
    identical(sha(first_pdf), sha(second_pdf)),
    !file.exists(file.path(figure_dir, "Rplots.pdf"))
  )
})

cat("task5 fixture tests passed\n")
