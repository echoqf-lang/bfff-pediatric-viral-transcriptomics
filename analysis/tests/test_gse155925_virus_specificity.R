#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(testthat))

test_that("GSE155925 specificity outputs preserve all frozen genes", {
  root <- file.path("results", "single_gene_upgrade_v1", "specificity")
  effect_path <- file.path(root, "GSE155925_86_gene_rsv_vs_other_virus.tsv")
  class_path <- file.path(root, "cross_context_gene_classification.tsv")
  summary_path <- file.path(root, "specificity_summary.tsv")
  expect_true(file.exists(effect_path))
  expect_true(file.exists(class_path))
  expect_true(file.exists(summary_path))
  effects <- read.delim(effect_path, check.names = FALSE)
  classes <- read.delim(class_path, check.names = FALSE)
  summary <- read.delim(summary_path, check.names = FALSE)
  expect_equal(nrow(effects), 86L)
  expect_equal(nrow(classes), 86L)
  expect_equal(length(unique(effects$gene_id)), 86L)
  expect_true(all(effects$direction == "single_RSV_minus_other_single_virus"))
  expect_equal(sum(summary$n_genes), 86L)
  expect_true(all(summary$n_rsv == 31L))
  expect_true(all(summary$n_other_single_virus == 17L))
  allowed <- c(
    "relative_rsv_enriched",
    "relative_rsv_depleted_or_context_dependent",
    "not_distinguished_from_other_viruses_at_current_precision",
    "not_estimable"
  )
  expect_true(all(unique(effects$specificity_class) %in% allowed))
  expect_true("not_distinguished_from_other_viruses_at_current_precision" %in% effects$specificity_class)
})
