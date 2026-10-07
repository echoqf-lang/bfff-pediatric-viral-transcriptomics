#!/usr/bin/env Rscript

source("analysis/R/24_bfff_full_86_airway_validation.R")

contract <- load_target_contract(
  "data/clean/bfff_86_gene_set_frozen.tsv",
  "analysis/config/bfff_barrier_inflammation_repair_panel.tsv"
)

stopifnot(nrow(contract) == 86L)
stopifnot(!anyDuplicated(contract$gene_symbol))
stopifnot(sum(contract$is_representative) == 20L)
stopifnot(all(contract$axis[!contract$is_representative] == "Other_unclassified"))
stopifnot(all(c("discovery_source", "axis", "role", "priority", "is_representative") %in% names(contract)))

output_root <- "results/bfff_full_86_airway_v2"
if (dir.exists(output_root)) {
  gse97742 <- read.delim(file.path(output_root, "GSE97742_full_86_RSV_vs_hRV_recovery_interaction.tsv"))
  gse41374 <- read.delim(file.path(output_root, "GSE41374_full_86_robust_summary.tsv"))
  stopifnot(nrow(gse97742) == 86L, nrow(gse41374) == 86L)
  stopifnot(sum(!is.na(gse97742$p_value)) == 85L)
  mapped <- !is.na(gse97742$p_value)
  stopifnot(isTRUE(all.equal(
    gse97742$fdr_bh_86_family[mapped],
    p.adjust(gse97742$p_value[mapped], method = "BH")
  )))
  stopifnot(sum(!is.na(gse41374$logFC.shifted_log2)) == 86L)
}

message("Full-86 airway contract test passed")
