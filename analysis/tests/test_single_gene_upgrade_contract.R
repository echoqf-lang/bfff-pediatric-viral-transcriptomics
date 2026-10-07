suppressPackageStartupMessages(library(testthat))

test_that("single-gene upgrade contract is frozen exactly", {
  expect_true(file.exists("analysis/config/single_gene_upgrade.yml"))
  expect_true(file.exists("data/clean/bfff_86_gene_set_frozen.tsv"))
  expect_true(file.exists("logs/checksums/single_gene_upgrade_inputs_sha256.tsv"))
  cfg <- yaml::read_yaml("analysis/config/single_gene_upgrade.yml")
  expect_identical(cfg$analysis_label, "post_outcome_exploratory_multicohort_upgrade_v1")
  expect_identical(cfg$meta$method, "REML_HKSJ")
  expect_equal(cfg$families$single_gene_n, 86L)
  expect_equal(unlist(cfg$severity$codes, use.names = FALSE), c(0L, 1L, 2L))
  expect_false(cfg$docking$run_by_default)

  genes <- read.delim("data/clean/bfff_86_gene_set_frozen.tsv", colClasses = "character")
  expect_equal(nrow(genes), 86L)
  expect_equal(length(unique(genes$entrez_id)), 86L)
  expect_equal(length(unique(genes$gene_symbol)), 86L)
  expect_false(anyNA(genes$entrez_id))
  expect_identical(genes$gene_symbol, sort(genes$gene_symbol))
})

