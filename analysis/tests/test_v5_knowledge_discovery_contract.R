#!/usr/bin/env Rscript

path <- "manuscript/draft_v5_english_knowledge_discovery.md"
stopifnot(file.exists(path))
text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

required <- c(
  "## 1. Introduction",
  "## 2. Materials and methods",
  "## 3. Results",
  "## 4. Discussion",
  "## 5. Conclusions",
  "ciliary impairment",
  "inflammatory activation",
  "phase-dependent repair",
  "85 of 86",
  "24 of 86",
  "16 candidates",
  "GSE97742",
  "GSE41374",
  "None of the analyzed cohorts included BFFF exposure"
)
stopifnot(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))

forbidden <- c(
  "BFFF regulated",
  "BFFF reversed",
  "BFFF inhibited",
  "BFFF activated",
  "confirmed therapeutic targets",
  "established mechanism of BFFF"
)
stopifnot(!any(vapply(forbidden, grepl, logical(1), x = text, fixed = TRUE)))

message("V5 knowledge-discovery manuscript contract passed")
