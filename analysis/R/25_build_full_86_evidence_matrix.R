#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
source("analysis/R/24_bfff_full_86_airway_validation.R")

append_columns <- function(base, addition, columns, prefix) {
  index <- match(base$gene_symbol, addition$gene_symbol)
  values <- addition[index, columns, drop = FALSE]
  names(values) <- paste0(prefix, gsub("\\.", "_", columns))
  cbind(base, values)
}

build_evidence_matrix <- function(contract) {
  blood_root <- file.path("results", "single_gene_upgrade_v1", "cohort")
  airway_root <- file.path("results", "bfff_full_86_airway_v2")
  output <- contract

  for (cohort in c("GSE38900", "GSE77087", "GSE103842")) {
    effects <- read.delim(file.path(blood_root, paste0(cohort, "_86_gene_effects.tsv")))
    output <- append_columns(
      output, effects, c("log2FC", "SE", "p_value", "fdr_bh"), paste0(cohort, "_")
    )
  }

  recovery <- read.delim(file.path(
    airway_root, "GSE97742_full_86_RSV_vs_hRV_recovery_interaction.tsv"
  ))
  names(recovery)[names(recovery) == "rsv_minus_hrv_delta"] <- "RSV_minus_hRV_delta"
  output <- append_columns(
    output, recovery,
    c("RSV_minus_hRV_delta", "se", "p_value", "fdr_bh_86_family"),
    "GSE97742_"
  )

  cross_sectional <- read.delim(file.path(
    airway_root, "GSE41374_full_86_robust_summary.tsv"
  ))
  output <- append_columns(
    output, cross_sectional,
    c(
      "logFC.shifted_log2", "p_value.shifted_log2", "fdr_bh_86_family.shifted_log2",
      "logFC.asinh_robust", "p_value.asinh_robust", "fdr_bh_86_family.asinh_robust",
      "direction_concordant", "robust_86_family_fdr"
    ),
    "GSE41374_"
  )

  output$blood_GSE77087_direction_concordant <- with(
    output, sign(GSE38900_log2FC) == sign(GSE77087_log2FC)
  )
  output$blood_GSE103842_direction_concordant <- with(
    output, sign(GSE38900_log2FC) == sign(GSE103842_log2FC)
  )
  output$blood_GSE77087_genome_fdr_support <- with(
    output, blood_GSE77087_direction_concordant & GSE77087_fdr_bh < 0.05
  )
  output$blood_GSE103842_genome_fdr_support <- with(
    output, blood_GSE103842_direction_concordant & GSE103842_fdr_bh < 0.05
  )
  output$airway_GSE97742_fdr_support <- output$GSE97742_fdr_bh_86_family < 0.05
  output$airway_GSE41374_fdr_support <- output$GSE41374_robust_86_family_fdr
  output$airway_any_fdr_support <- with(
    output, airway_GSE97742_fdr_support | airway_GSE41374_fdr_support
  )
  support_columns <- c(
    "blood_GSE77087_genome_fdr_support", "blood_GSE103842_genome_fdr_support",
    "airway_GSE97742_fdr_support", "airway_GSE41374_fdr_support"
  )
  output$n_external_contexts_fdr <- rowSums(output[, support_columns], na.rm = TRUE)
  output$evidence_class <- ifelse(
    output$blood_GSE77087_genome_fdr_support & output$airway_any_fdr_support,
    "blood_external_and_airway",
    ifelse(
      output$blood_GSE77087_genome_fdr_support,
      "blood_external_only",
      ifelse(output$airway_any_fdr_support, "airway_only", "no_selected_external_fdr")
    )
  )
  output
}

write_evidence_outputs <- function(matrix) {
  output_root <- file.path("results", "bfff_full_86_evidence_matrix_v1")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  write.table(
    matrix, file.path(output_root, "Table_full_86_cross_context_evidence.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )

  summary <- data.frame(
    metric = c(
      "target_universe", "representative_nodes", "GSE77087_concordant_genome_fdr_support",
      "GSE103842_concordant_genome_fdr_support", "GSE97742_airway_fdr_support",
      "GSE41374_airway_robust_fdr_support", "any_airway_fdr_support",
      "blood_GSE77087_and_any_airway_support"
    ),
    value = c(
      nrow(matrix), sum(matrix$is_representative),
      sum(matrix$blood_GSE77087_genome_fdr_support, na.rm = TRUE),
      sum(matrix$blood_GSE103842_genome_fdr_support, na.rm = TRUE),
      sum(matrix$airway_GSE97742_fdr_support, na.rm = TRUE),
      sum(matrix$airway_GSE41374_fdr_support, na.rm = TRUE),
      sum(matrix$airway_any_fdr_support, na.rm = TRUE),
      sum(matrix$blood_GSE77087_genome_fdr_support & matrix$airway_any_fdr_support, na.rm = TRUE)
    )
  )
  write.table(summary, file.path(output_root, "evidence_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  axis_levels <- unique(matrix$axis)
  axis_summary <- do.call(rbind, lapply(axis_levels, function(axis_name) {
    part <- matrix[matrix$axis == axis_name, , drop = FALSE]
    data.frame(
      axis = axis_name,
      n_genes = nrow(part),
      n_GSE77087_concordant_genome_fdr = sum(part$blood_GSE77087_genome_fdr_support, na.rm = TRUE),
      n_any_airway_fdr = sum(part$airway_any_fdr_support, na.rm = TRUE),
      n_blood_and_airway = sum(
        part$blood_GSE77087_genome_fdr_support & part$airway_any_fdr_support, na.rm = TRUE
      )
    )
  }))
  write.table(axis_summary, file.path(output_root, "axis_evidence_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(output_root, "session_info.txt"))
  invisible(summary)
}

if (sys.nframe() == 0L) {
  contract <- load_target_contract(
    "data/clean/bfff_86_gene_set_frozen.tsv",
    "analysis/config/bfff_barrier_inflammation_repair_panel.tsv"
  )
  summary <- write_evidence_outputs(build_evidence_matrix(contract))
  message(
    "Full-86 evidence matrix complete: genes=", summary$value[summary$metric == "target_universe"],
    " any_airway_fdr=", summary$value[summary$metric == "any_airway_fdr_support"]
  )
}
