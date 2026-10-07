#!/usr/bin/env Rscript

root <- "manuscript/v5_submission_assets"
figures <- file.path(root, "figures", sprintf("Figure_%d", 1:5))
for (stem in figures) {
  stopifnot(file.exists(paste0(stem, ".png")))
  stopifnot(file.exists(paste0(stem, ".pdf")))
  stopifnot(file.info(paste0(stem, ".png"))$size > 10000)
  stopifnot(file.info(paste0(stem, ".pdf"))$size > 5000)
}

tables <- file.path(root, "tables", c(
  "Table_1_cohort_characteristics.tsv",
  "Table_2_86_gene_evidence_summary.tsv",
  "Table_3_16_bridge_candidates.tsv"
))
stopifnot(all(file.exists(tables)))
stopifnot(nrow(read.delim(tables[1], check.names = FALSE)) == 5L)
stopifnot(nrow(read.delim(tables[2], check.names = FALSE)) >= 8L)
stopifnot(nrow(read.delim(tables[3], check.names = FALSE)) == 16L)
stopifnot(file.exists(file.path(root, "tables", "V5_main_tables.xlsx")))

message("V5 figure/table output contract passed")
