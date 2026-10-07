#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(testthat))
source(file.path("analysis", "R", "single_gene_upgrade_pipeline.R"))

test_that("effect-table assertions reject invalid inputs", {
  good <- data.frame(gene_id = c("1", "2"), cohort = "A", direction = "RSV_minus_healthy", log2FC = c(0.1, 0.2), SE = c(0.1, 0.1))
  expect_invisible(assert_unique_gene_effects(good, "A"))
  expect_error(assert_unique_gene_effects(rbind(good, good[1, ]), "A"), "unique")
  bad <- good
  bad$SE[1] <- 0
  expect_error(assert_unique_gene_effects(bad, "A"), "positive")
})

test_that("HKSJ pipeline recovers a known simulated signal", {
  set.seed(20260805L)
  genes <- as.character(seq_len(86L))
  cohorts <- c("D", "R1", "R2")
  truth <- c(rep(0.8, 10L), rep(0, 76L))
  simulated <- do.call(rbind, lapply(seq_along(cohorts), function(j) {
    data.frame(
      gene_id = genes,
      cohort = cohorts[j],
      direction = "RSV_minus_healthy",
      log2FC = truth + stats::rnorm(86L, sd = 0.08),
      SE = rep(0.10, 86L),
      stringsAsFactors = FALSE
    )
  }))
  result <- fit_hksj_meta(simulated, cohorts)
  result <- result[match(genes, result$gene_id), ]
  expect_equal(nrow(result), 86L)
  expect_true(all(c("gene_id", "estimate", "se", "ci_low", "ci_high", "adhoc_ci_low", "adhoc_ci_high", "tau2", "i2", "p_value", "fdr_bh", "adhoc_fdr_bh") %in% names(result)))
  expect_gt(mean(result$estimate[1:10]), 0.5)
  expect_lt(abs(mean(result$estimate[11:86])), 0.2)
  expect_true(all(result$k == 3L))
})

test_that("direction patterns are deterministic", {
  expect_identical(classify_direction_pattern(1, 2, 3), "all_up")
  expect_identical(classify_direction_pattern(-1, -2, -3), "all_down")
  expect_identical(classify_direction_pattern(1, 2, -3), "two_up_one_down")
})

test_that("real three-cohort outputs honor the frozen contract when present", {
  root <- file.path("results", "single_gene_upgrade_v1")
  paths <- file.path(root, "cohort", c("GSE38900_86_gene_effects.tsv", "GSE77087_86_gene_effects.tsv", "GSE103842_86_gene_effects.tsv"))
  skip_if_not(all(file.exists(paths)), "real outputs not generated yet")
  x <- lapply(paths, read.delim, check.names = FALSE)
  expect_true(all(vapply(x, nrow, integer(1)) == 86L))
  expect_true(all(vapply(x, function(z) !anyDuplicated(z$gene_id), logical(1))))
  summary <- read.delim(file.path(root, "meta", "three_cohort_meta_summary.tsv"), check.names = FALSE)
  expect_equal(summary$gse77087_mapped_n, 86L)
  expect_equal(summary$gse38900_gse77087_direction_concordant_n, 77L)
  expect_equal(summary$gse77087_candidate_bh_and_direction_n, 46L)
  expect_equal(summary$gse77087_genome_fdr_and_direction_n, 41L)
  expect_equal(summary$gse77087_genome_fdr_effect_and_direction_n, 39L)
  meta <- read.delim(file.path(root, "meta", "three_cohort_86_gene_reml_hksj.tsv"), check.names = FALSE)
  expect_true(all(meta$k == 3L))
  expect_equal(nrow(meta), summary$meta_complete_n)
})
