#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
if (!file.exists(file.path("manuscript", "draft_v4_english_submission.md"))) setwd(file.path("..", ".."))

path <- file.path("manuscript", "draft_v4_english_submission.md")
text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

required <- c(
  "Eighty-one genes had complete estimates", "69 of 81", "Five genes", "no gene met BH FDR <0.05",
  "41 genes", "ADORA2B, CACNA1E, CHRNA5, and HBE1", "only DHFR",
  "345 modules", "213 had concordant median effects", "68 were significant",
  "stopped_by_qc", "520 frozen targets"
)
stopifnot(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))

boundary <- c(
  "None of the cohorts included exposure",
  "do not establish BFFF-mediated regulation",
  "not validation of a drug mechanism"
)
stopifnot(all(vapply(boundary, grepl, logical(1), x = text, fixed = TRUE)))

forbidden <- c("proves the mechanism", "validated therapeutic target", "precision treatment", "selectively activates")
stopifnot(!any(vapply(forbidden, grepl, logical(1), x = text, fixed = TRUE)))

figure_labels <- paste0("**Figure ", 1:4, ".")
stopifnot(all(figure_labels %in% regmatches(text, gregexpr("\\*\\*Figure [1-4]\\.", text))[[1]]))
stopifnot(grepl("Repository link to be added", text, fixed = TRUE))
stopifnot(grepl("To be completed by the authors", text, fixed = TRUE))

message("English submission contract checks passed.")
