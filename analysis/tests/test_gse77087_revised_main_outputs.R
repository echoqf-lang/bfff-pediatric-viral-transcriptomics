#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(testthat)

root <- file.path("results", "revised_main_gse77087")
required <- file.path(root, c(
  "analysis_summary.tsv",
  "cohort/GSE77087_sample_manifest.tsv",
  "cohort/GSE77087_blind_qc.tsv",
  "cohort/model_diagnostics.tsv",
  "target_set/camera.tsv",
  "target_set/weighted_stouffer.tsv",
  "target_set/roast.tsv",
  "target_set/fgsea_exploratory.tsv",
  "meta/target_gene_reml_hksj.tsv",
  "sensitivity/GSE77087_hospitalized_only_camera.tsv",
  "tables/Table2_revised_target_set.tsv",
  "figures/Fig2_revised_target_set.pdf"
))
expect_true(all(file.exists(required)))

read_tsv <- function(path) read.delim(path, sep = "\t", quote = "", check.names = FALSE)

test_that("修订结果样本数、模型和证据标签一致", {
  summary <- read_tsv(file.path(root, "analysis_summary.tsv"))
  diagnostics <- read_tsv(file.path(root, "cohort/model_diagnostics.tsv"))
  expect_equal(summary$gse77087_pre_qc_n, 104L)
  expect_equal(summary$gse77087_qc_excluded_n, 0L)
  expect_equal(summary$gse77087_model_n, 104L)
  expect_equal(summary$gse103842_model_n, 73L)
  expect_false(summary$original_preregistered_h1_changed)
  expect_true(all(diagnostics$batch_included))
  expect_true(all(grepl("technical_batch", diagnostics$model_formula)))
})

test_that("CAMERA、Stouffer和Meta结果完整且不可被写成显著", {
  camera <- read_tsv(file.path(root, "target_set/camera.tsv"))
  stouffer <- read_tsv(file.path(root, "target_set/weighted_stouffer.tsv"))
  meta <- read_tsv(file.path(root, "meta/target_gene_reml_hksj.tsv"))
  expect_identical(as.character(camera$cohort), c("GSE77087", "GSE103842"))
  expect_true(all(camera$direction == "Up"))
  expect_true(all(camera$p_value >= 0.05))
  expect_gt(stouffer$weighted_stouffer_two_sided_p, 0.05)
  expect_equal(nrow(meta), 374L)
  expect_equal(sum(meta$fdr_bh < 0.05), 0L)
})

test_that("随机检验记录的种子与实际修订种子一致", {
  roast <- read_tsv(file.path(root, "target_set/roast.tsv"))
  fgsea <- read_tsv(file.path(root, "target_set/fgsea_exploratory.tsv"))
  expect_true(all(roast$seed == 20260805L))
  expect_true(all(fgsea$seed == 20260805L))
  expect_true(all(roast$nrot == 99999L))
})

test_that("没有下载GSE77087原始归档", {
  files <- list.files(file.path("data", "raw", "GSE77087"), recursive = TRUE, full.names = FALSE)
  expect_false(any(grepl("RAW|[.]CEL|FASTQ|SRA", files, ignore.case = TRUE)))
  expect_setequal(files, c("filelist.txt", "GSE77087_quick.soft.txt", "GSE77087_series_matrix.txt.gz"))
})

cat("GSE77087 revised-main output tests passed\n")
