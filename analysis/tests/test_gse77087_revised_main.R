#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(testthat)

project_root <- normalizePath(getwd(), mustWork = TRUE)
implementation <- file.path(project_root, "analysis", "R", "12_gse77087_revised_main.R")
expect_true(file.exists(implementation), info = "修订主分析脚本尚未实现")

old_skip <- Sys.getenv("GSE77087_REVISED_SKIP_MAIN", unset = NA_character_)
Sys.setenv(GSE77087_REVISED_SKIP_MAIN = "1")
source(implementation, local = FALSE)
if (is.na(old_skip)) Sys.unsetenv("GSE77087_REVISED_SKIP_MAIN") else Sys.setenv(GSE77087_REVISED_SKIP_MAIN = old_skip)

test_that("元数据解析固定为104个独立全血样本", {
  metadata <- read_gse77087_metadata(file.path("data", "raw", "GSE77087", "GSE77087_series_matrix.txt.gz"))
  expect_equal(nrow(metadata), 104L)
  expect_equal(sum(metadata$case_status == "healthy"), 23L)
  expect_equal(sum(metadata$case_status == "RSV"), 81L)
  expect_equal(sum(metadata$clinical_group == "inpatient"), 61L)
  expect_equal(sum(metadata$clinical_group == "outpatient"), 20L)
  expect_setequal(unique(metadata$technical_batch), c("2011", "2015"))
  expect_false(anyDuplicated(metadata$sample_id) > 0L)
  expect_true(all(metadata$tissue == "whole blood"))
  expect_true(all(is.finite(metadata$age_months)))
  expect_true(all(metadata$sex %in% c("female", "male")))
})

test_that("修订主分析只允许GSE77087与GSE103842", {
  expect_identical(revised_main_cohorts, c("GSE77087", "GSE103842"))
  expect_false("GSE105450" %in% revised_main_cohorts)
  expect_identical(revised_analysis_label, "post_outcome_investigator_selected_revised_main")
  expect_identical(task8_seed, revised_seed)
})

test_that("处理后矩阵过滤完全盲于病例标签", {
  x <- rbind(
    invariant = rep(5, 8),
    low = seq(1, 1.7, length.out = 8),
    informative = c(5, 6, 7, 8, 5, 6, 7, 8)
  )
  first <- filter_processed_features_blind(x)
  second <- filter_processed_features_blind(x[, 8:1, drop = FALSE])
  expect_setequal(rownames(first$expression), rownames(second$expression))
  expect_false("invariant" %in% rownames(first$expression))
  expect_true("informative" %in% rownames(first$expression))
})

test_that("模型方向和协变量规则被锁定", {
  sample_data <- data.frame(
    sample_id = paste0("S", 1:12), subject_id = paste0("P", 1:12),
    cohort = "GSE77087", platform = "GPL10558",
    case_status = rep(c("healthy", "RSV"), each = 6),
    age_months = c(1, 2, 4, 7, 8, 11, 1.5, 3, 5, 6, 9, 12),
    sex = rep(c("female", "male"), 6),
    technical_batch = rep(rep(c("2011", "2015"), each = 3), 2)
  )
  expression <- matrix(rnorm(60), nrow = 5, dimnames = list(as.character(1:5), sample_data$sample_id))
  prepared <- prepare_revised_model_data(expression, sample_data, "GSE77087")
  expect_identical(prepared$coefficient, "case_statusRSV")
  expect_match(prepared$formula, "age_months")
  expect_match(prepared$formula, "sex")
  expect_match(prepared$formula, "technical_batch")
  expect_equal(qr(prepared$design)$rank, ncol(prepared$design))
})

cat("GSE77087 revised-main unit tests passed\n")
