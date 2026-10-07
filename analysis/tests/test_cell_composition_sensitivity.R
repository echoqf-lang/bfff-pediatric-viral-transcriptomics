#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(testthat))
if (!file.exists(file.path("analysis", "R", "18_cell_composition_sensitivity.R"))) setwd(file.path("..", ".."))
source(file.path("analysis", "R", "18_cell_composition_sensitivity.R"))

test_that("NNLS estimator recovers known mixtures", {
  result <- validate_fraction_estimator()
  expect_gt(result$spearman, 0.80)
  expect_true(all(result$estimates >= 0))
  expect_equal(colSums(result$estimates), rep(1, ncol(result$estimates)), tolerance = 1e-7)
})

test_that("missing or inapplicable signature produces an explicit stop", {
  config <- read_cell_config()
  qc <- build_stopped_qc(config, 0.95)
  expect_true(all(qc$analysis_status == "stopped_by_qc"))
  expect_false(any(qc$real_data_qc_pass))
  expect_true(all(nzchar(qc$stop_reason)))
})

test_that("real QC output records all three cohorts and no adjusted effects", {
  path <- file.path("results", "single_gene_upgrade_v1", "cell_sensitivity", "cell_fraction_qc.tsv")
  skip_if_not(file.exists(path), "real QC output not generated yet")
  qc <- read_upgrade_tsv(path)
  expect_equal(sort(qc$cohort), sort(c("GSE38900", "GSE77087", "GSE103842")))
  expect_true(all(qc$analysis_status == "stopped_by_qc"))
  expect_false(file.exists(file.path(dirname(path), "cell_fraction_group_effects.tsv")))
  expect_false(file.exists(file.path(dirname(path), "86_gene_effect_attenuation.tsv")))
})
