#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

library(testthat)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_root, "analysis", "R", "06_test_target_set.R"))) {
  project_root <- normalizePath(file.path(project_root, "..", ".."), mustWork = TRUE)
}
old_wd <- getwd()
setwd(project_root)
on.exit(setwd(old_wd), add = TRUE)

old_skip <- Sys.getenv("TASK8_SKIP_MAIN", unset = NA_character_)
Sys.setenv(TASK8_SKIP_MAIN = "1")
source(file.path("analysis", "R", "06_test_target_set.R"), local = FALSE)
if (is.na(old_skip)) Sys.unsetenv("TASK8_SKIP_MAIN") else Sys.setenv(TASK8_SKIP_MAIN = old_skip)

test_that("weighted Stouffer uses signed two-sided camera P values and frozen sample weights", {
  camera <- data.frame(
    cohort = c("GSE105450", "GSE103842"),
    status = "interpretable",
    direction = c("Up", "Up"),
    p_value = c(0.01, 0.02),
    signed_z = qnorm(c(0.01, 0.02) / 2, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
  out <- combine_camera_stouffer(camera)
  weights <- sqrt(c(122, 73))
  expected_z <- sum(weights * camera$signed_z) / sqrt(sum(weights^2))
  expect_equal(out$combined_z, expected_z, tolerance = 1e-14)
  expect_equal(out$combined_p_value, 2 * pnorm(-abs(expected_z)), tolerance = 1e-14)
  expect_identical(out$direction_consistent, TRUE)
  expect_identical(out$camera_replication_component_met, TRUE)
})

test_that("directional disagreement is retained and cannot pass the camera component", {
  camera <- data.frame(
    cohort = c("GSE105450", "GSE103842"),
    status = "interpretable",
    direction = c("Up", "Down"),
    p_value = c(0.001, 0.001),
    signed_z = c(1, -1) * qnorm(0.001 / 2, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
  out <- combine_camera_stouffer(camera)
  expect_identical(out$direction_consistent, FALSE)
  expect_identical(out$camera_replication_component_met, FALSE)
})

test_that("camera interpretability uses the complete 520-gene frozen set", {
  frozen <- as.character(seq_len(520L))
  expression_259 <- matrix(
    rnorm(259L * 12L), nrow = 259L,
    dimnames = list(frozen[seq_len(259L)], paste0("S", seq_len(12L)))
  )
  design <- model.matrix(~ rep(c("healthy", "RSV"), each = 6L))
  colnames(design)[2L] <- "case_statusRSV"
  out <- run_camera_prepared(expression_259, design, "case_statusRSV", frozen, "FIXTURE")
  expect_identical(out$status, "uninterpretable")
  expect_equal(out$n_frozen, 520L)
  expect_equal(out$n_detectable, 259L)
  expect_equal(out$coverage, 259 / 520)
  expect_true(is.na(out$p_value))
})

test_that("roast API extraction returns exactly Down, Up, and Mixed with Holm adjustment", {
  api_result <- list(
    p.value = matrix(
      c(0.04, 0.10, 0.07, 0.02, 0.25, 0.30, 0.20, 0.15),
      nrow = 4L, byrow = FALSE,
      dimnames = list(c("Down", "Up", "UpOrDown", "Mixed"), c("Active.Prop", "P.Value"))
    )
  )
  out <- extract_roast_directions(api_result, "FIXTURE", 520L, 300L, 99999L)
  expect_identical(out$test_direction, c("Down", "Up", "Mixed"))
  expect_equal(out$p_value, c(0.25, 0.30, 0.15))
  expect_equal(out$p_holm, p.adjust(c(0.25, 0.30, 0.15), method = "holm"))
  expect_true(all(out$nrot == 99999L))
  data_frame_api <- api_result
  data_frame_api$p.value <- as.data.frame(data_frame_api$p.value)
  out_data_frame <- extract_roast_directions(data_frame_api, "FIXTURE", 520L, 300L, 99999L)
  expect_equal(out_data_frame, out)
})

test_that("partial decision remains pending until random-set calibration", {
  camera <- data.frame(
    cohort = c("GSE105450", "GSE103842"),
    status = "interpretable",
    direction = c("Up", "Up"),
    p_value = c(0.01, 0.02),
    signed_z = qnorm(c(0.01, 0.02) / 2, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
  decision <- make_partial_decision(camera, combine_camera_stouffer(camera))
  expect_identical(decision$primary_decision, "pending_random_set_calibration")
  expect_true(is.na(decision$h1_supported))
  expect_true(is.na(decision$random_set_empirical_p_value))
  expect_identical(decision$camera_replication_component_met, TRUE)
  expect_identical(decision$random_set_calibration_status, "pending_required_for_final_support_decision")
})

test_that("failed camera component fixes H1 as unsupported before non-rescuing random calibration", {
  camera <- data.frame(
    cohort = c("GSE105450", "GSE103842"),
    status = "interpretable",
    direction = c("Up", "Up"),
    p_value = c(0.20, 0.13),
    signed_z = qnorm(c(0.20, 0.13) / 2, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
  stouffer <- combine_camera_stouffer(camera)
  expect_false(stouffer$camera_replication_component_met)
  decision <- make_partial_decision(camera, stouffer)
  expect_identical(decision$primary_decision, "H1_not_supported_camera_component_failed")
  expect_identical(decision$h1_supported, FALSE)
  expect_identical(decision$random_set_calibration_status, "pending_scheduled_nonrescuing")
  expect_true(is.na(decision$random_set_empirical_p_value))
})

test_that("gene-effect reproduction normalizes Entrez storage type but not values", {
  saved <- data.frame(
    gene_id = 1:2, cohort = "FIXTURE", direction = "RSV_minus_healthy",
    log2FC = c(0.1, -0.2), SE = c(0.05, 0.08), moderated_t = c(2, -2.5),
    p_value = c(0.04, 0.02), fdr_bh = c(0.04, 0.04), average_expression = c(5, 6)
  )
  reconstructed <- saved
  reconstructed$gene_id <- as.character(reconstructed$gene_id)
  expect_true(same_gene_effects(saved, reconstructed))
  reconstructed$log2FC[[1L]] <- reconstructed$log2FC[[1L]] + 0.01
  expect_false(same_gene_effects(saved, reconstructed))
})

cat("task8 unit tests passed\n")
