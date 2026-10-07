#!/usr/bin/env Rscript

library(testthat)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_root, "analysis", "R", "05_fit_cohort_models.R"))) {
  project_root <- normalizePath(file.path(project_root, "..", ".."), mustWork = TRUE)
}
old_wd <- getwd()
setwd(project_root)
on.exit(setwd(old_wd), add = TRUE)

old_skip <- Sys.getenv("TASK7_SKIP_MAIN", unset = NA_character_)
Sys.setenv(TASK7_SKIP_MAIN = "1")
source(file.path("analysis", "R", "05_fit_cohort_models.R"), local = FALSE)
if (is.na(old_skip)) Sys.unsetenv("TASK7_SKIP_MAIN") else Sys.setenv(TASK7_SKIP_MAIN = old_skip)

make_fixture <- function(confounded_batch = FALSE) {
  sample_id <- sprintf("S%02d", seq_len(16L))
  case_status <- rep(c("healthy", "RSV"), each = 8L)
  batch <- if (confounded_batch) {
    ifelse(case_status == "healthy", "B1", "B2")
  } else {
    rep(c("B1", "B2", "B2", "B1", "B2", "B1", "B1", "B2"), times = 2L)
  }
  sample_data <- data.frame(
    sample_id = sample_id,
    subject_id = paste0("P", sample_id),
    cohort = "FIXTURE",
    platform = "GPL_FIXTURE",
    case_status = case_status,
    age_months = rep(c(2, 4, 6, 8), times = 4L),
    sex = rep(c("female", "male"), times = 8L),
    technical_batch = batch,
    stringsAsFactors = FALSE
  )
  set.seed(20260805L)
  expression <- matrix(
    rnorm(40L * 16L, sd = 0.15), nrow = 40L,
    dimnames = list(as.character(seq_len(40L)), sample_id)
  )
  expression[1L, case_status == "RSV"] <- expression[1L, case_status == "RSV"] + 1.5
  list(expression = expression, sample_data = sample_data)
}

test_that("fixture recovers RSV-minus-healthy direction and full output schema", {
  fixture <- make_fixture()
  fit <- fit_task7_cohort(fixture$expression, fixture$sample_data, "FIXTURE")
  expected <- c(
    "gene_id", "cohort", "direction", "log2FC", "SE", "moderated_t",
    "p_value", "fdr_bh", "average_expression"
  )
  expect_identical(names(fit$gene_effects), expected)
  expect_gt(fit$gene_effects$log2FC[fit$gene_effects$gene_id == "1"], 1)
  expect_true(all(is.finite(fit$gene_effects$SE) & fit$gene_effects$SE > 0))
  expect_identical(fit$coefficient, "case_statusRSV")
  expect_true(fit$batch_included)
  expect_equal(fit$design_rank, fit$design_columns)
})

test_that("complete-case rule excludes missing age without imputation", {
  fixture <- make_fixture()
  fixture$sample_data$age_months[[1L]] <- NA_real_
  prepared <- prepare_model_data(fixture$expression, fixture$sample_data, "FIXTURE")
  expect_equal(prepared$n_input, 16L)
  expect_equal(prepared$n_complete, 15L)
  expect_equal(prepared$n_excluded_missing_covariate, 1L)
  expect_false("S01" %in% prepared$sample_data$sample_id)
  expect_equal(ncol(prepared$expression), 15L)
})

test_that("batch is omitted when it makes the design rank deficient", {
  fixture <- make_fixture(confounded_batch = TRUE)
  fit <- fit_task7_cohort(fixture$expression, fixture$sample_data, "FIXTURE")
  expect_false(fit$batch_included)
  expect_identical(fit$batch_reason, "omitted_rank_deficient_or_status_collinear")
  expect_equal(fit$design_rank, fit$design_columns)
})

test_that("mapping and analysis boundaries fail closed", {
  fixture <- make_fixture()
  bad_order <- fixture$sample_data[rev(seq_len(nrow(fixture$sample_data))), ]
  expect_error(prepare_model_data(fixture$expression, bad_order, "FIXTURE"), "same order")
  duplicate <- fixture$sample_data
  duplicate$subject_id[[2L]] <- duplicate$subject_id[[1L]]
  expect_error(prepare_model_data(fixture$expression, duplicate, "FIXTURE"), "Duplicate subject_id")
  one_group <- fixture$sample_data
  one_group$case_status <- "RSV"
  expect_error(fit_task7_cohort(fixture$expression, one_group, "FIXTURE"), "both healthy and RSV")
  expect_error(fit_task7_cohort(fixture$expression, fixture$sample_data, "FIXTURE", "healthy_minus_RSV"), "RSV_minus_healthy")
})

cat("task7 unit tests passed\n")
