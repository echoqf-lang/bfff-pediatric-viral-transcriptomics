#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(testthat))
if (!file.exists(file.path("analysis", "config", "single_gene_upgrade.yml"))) setwd(file.path("..", ".."))

root <- file.path("results", "single_gene_upgrade_v1")

test_that("main and supplement tables preserve the full analysis families", {
  main_path <- file.path(root, "tables", "Table_main_three_cohort_meta.tsv")
  supp_path <- file.path(root, "tables", "Table_supplement_all_86_genes.tsv")
  expect_true(file.exists(main_path), info = "main table is missing")
  expect_true(file.exists(supp_path), info = "supplement table is missing")
  main <- read.delim(main_path, check.names = FALSE)
  supp <- read.delim(supp_path, check.names = FALSE)
  expect_equal(nrow(main), 81L)
  expect_equal(nrow(supp), 86L)
  expect_equal(anyDuplicated(main$gene_id), 0L)
  expect_equal(anyDuplicated(supp$gene_id), 0L)
  expect_true(all(c("estimate", "adhoc_ci_low", "adhoc_ci_high", "direction_pattern") %in% names(main)))
  expect_true(all(c("GSE38900_log2FC", "GSE77087_log2FC", "GSE103842_log2FC", "severity_trend_log2FC", "specificity_log2FC") %in% names(supp)))
})

test_that("every planned figure has PDF SVG and 600 dpi PNG outputs", {
  stems <- c("Fig_multicohort_effects", "Fig_severity_gradient", "Fig_virus_specificity", "Fig_blood_modules")
  paths <- unlist(lapply(stems, function(stem) file.path(root, "figures", paste0(stem, c(".pdf", ".svg", ".png")))))
  expect_true(all(file.exists(paths)), info = paste("missing:", paste(paths[!file.exists(paths)], collapse = ", ")))
  expect_true(all(file.info(paths)$size > 1000L))
})

test_that("manuscript count statements match frozen summaries", {
  manuscript <- paste(readLines(file.path("manuscript", "draft_v3_单基因跨队列证据边界版.md"), warn = FALSE), collapse = "\n")
  meta <- read.delim(file.path(root, "meta", "three_cohort_meta_summary.tsv"), check.names = FALSE)
  severity <- read.delim(file.path(root, "severity", "severity_model_diagnostics.tsv"), check.names = FALSE)
  specificity <- read.delim(file.path(root, "specificity", "specificity_summary.tsv"), check.names = FALSE)
  modules <- read.delim(file.path(root, "modules", "module_replication_summary.tsv"), check.names = FALSE)
  expect_match(manuscript, paste0(meta$meta_complete_n, "个基因完成三队列"), fixed = TRUE)
  expect_match(manuscript, paste0(meta$all_three_same_direction_n, "个三队列方向一致"), fixed = TRUE)
  expect_match(manuscript, paste0(severity$trend_bh_lt_0_05_n, "个基因呈有序趋势"), fixed = TRUE)
  expect_match(manuscript, "DHFR", fixed = TRUE)
  expect_match(manuscript, paste0(sum(modules$evaluable_all_three), "个模块在三个队列均可评价"), fixed = TRUE)
  expect_match(manuscript, "stopped_by_qc", fixed = TRUE)
  expect_equal(sum(specificity$n_genes), 86L)
})
