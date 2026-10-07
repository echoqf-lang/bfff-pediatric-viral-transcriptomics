#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(testthat))
if (!file.exists(file.path("analysis", "R", "17_fixed_blood_modules.R"))) setwd(file.path("..", ".."))
source(file.path("analysis", "R", "17_fixed_blood_modules.R"))

test_that("fixed module dictionary is complete and unique", {
  path <- file.path("data", "clean", "fixed_blood_modules.tsv")
  skip_if_not(file.exists(path), "fixed module dictionary not generated yet")
  modules <- read_upgrade_tsv(path)
  expect_invisible(assert_module_dictionary(modules))
  expect_equal(length(unique(modules$module_id)), 382L)
  expect_equal(nrow(modules), 14168L)
})

test_that("module score recovers simulated direction", {
  modules <- data.frame(module_id = rep("M.test", 10), gene_symbol = paste0("G", 1:10),
                        function_name = "simulation", position = "X", module_color = "#FFFFFF", cluster = "X")
  names(modules)[names(modules) == "function_name"] <- "function"
  positive <- data.frame(gene_symbol = c(paste0("G", 1:10), paste0("B", 1:30)),
                         log2FC = c(rep(1, 10), rep(0, 30)), statistic = c(rep(5, 10), rep(0, 30)))
  negative <- positive
  negative$log2FC <- -negative$log2FC
  negative$statistic <- -negative$statistic
  up <- score_fixed_modules(positive, modules, "sim_up")
  down <- score_fixed_modules(negative, modules, "sim_down")
  expect_gt(up$median_log2FC, 0)
  expect_lt(down$median_log2FC, 0)
  expect_identical(up$camera_direction, "Up")
  expect_identical(down$camera_direction, "Down")
})

test_that("low coverage modules remain explicit and are not tested", {
  modules <- data.frame(module_id = rep("M.low", 10), gene_symbol = paste0("G", 1:10),
                        function_name = "simulation", position = "X", module_color = "#FFFFFF", cluster = "X")
  names(modules)[names(modules) == "function_name"] <- "function"
  effects <- data.frame(gene_symbol = paste0("G", 1:4), log2FC = rep(1, 4), statistic = rep(3, 4))
  result <- score_fixed_modules(effects, modules, "sim")
  expect_false(result$evaluable)
  expect_equal(result$coverage, 0.4)
  expect_true(is.na(result$p_value))
})

test_that("real module outputs retain every module for every cohort", {
  path <- file.path("results", "single_gene_upgrade_v1", "modules", "cohort_module_effects.tsv")
  skip_if_not(file.exists(path), "real module outputs not generated yet")
  x <- read_upgrade_tsv(path)
  expect_equal(nrow(x), 382L * 3L)
  expect_equal(sort(unique(x$cohort)), sort(c("GSE38900", "GSE77087", "GSE103842")))
  expect_true(all(table(x$cohort) == 382L))
  expect_true(all(x$coverage >= 0 & x$coverage <= 1))
  expect_true(all(is.na(x$p_value[!x$evaluable])))
})
