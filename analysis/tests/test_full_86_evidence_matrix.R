#!/usr/bin/env Rscript

source("analysis/R/25_build_full_86_evidence_matrix.R")

contract <- load_target_contract(
  "data/clean/bfff_86_gene_set_frozen.tsv",
  "analysis/config/bfff_barrier_inflammation_repair_panel.tsv"
)
matrix <- build_evidence_matrix(contract)

stopifnot(nrow(matrix) == 86L)
stopifnot(!anyDuplicated(matrix$gene_symbol))
stopifnot(sum(matrix$is_representative) == 20L)
stopifnot(all(c(
  "GSE38900_log2FC", "GSE77087_log2FC", "GSE103842_log2FC",
  "GSE97742_RSV_minus_hRV_delta", "GSE41374_logFC_shifted_log2",
  "blood_GSE77087_direction_concordant", "airway_any_fdr_support"
) %in% names(matrix)))

message("Full-86 evidence matrix test passed")
