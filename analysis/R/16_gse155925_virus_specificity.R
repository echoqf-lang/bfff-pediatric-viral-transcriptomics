#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
source(file.path("analysis", "R", "single_gene_upgrade_pipeline.R"))

main_gse155925_specificity <- function() {
  frozen <- read_upgrade_tsv(file.path("data", "clean", "bfff_86_gene_set_frozen.tsv"), col_classes = "character")
  existing <- read_upgrade_tsv(file.path("results", "sensitivity", "gene_effects_GSE155925.tsv"))
  existing$gene_id <- as.character(existing$gene_id)
  idx <- match(frozen$ensembl_gene_id, existing$gene_id)
  out <- data.frame(
    gene_id = frozen$entrez_id,
    gene_symbol = frozen$gene_symbol,
    ensembl_gene_id = frozen$ensembl_gene_id,
    cohort = "GSE155925",
    direction = "single_RSV_minus_other_single_virus",
    log2FC = existing$log2FC[idx],
    SE = existing$SE[idx],
    p_value = existing$p_value[idx],
    genome_fdr = existing$fdr_bh[idx],
    stringsAsFactors = FALSE
  )
  out$candidate_fdr <- bh_family(out$p_value)
  out$ci_low <- out$log2FC - 1.96 * out$SE
  out$ci_high <- out$log2FC + 1.96 * out$SE

  meta <- read_upgrade_tsv(file.path("results", "single_gene_upgrade_v1", "meta", "three_cohort_86_gene_reml_hksj.tsv"))
  out$three_cohort_meta_effect <- meta$estimate[match(out$gene_id, as.character(meta$gene_id))]
  out$direction_matches_rsv_vs_healthy <- sign(out$log2FC) == sign(out$three_cohort_meta_effect)
  estimable <- complete.cases(out[c("log2FC", "SE", "candidate_fdr", "three_cohort_meta_effect")])
  out$specificity_class <- "not_estimable"
  out$specificity_class[estimable & out$candidate_fdr < 0.05 & out$direction_matches_rsv_vs_healthy] <- "relative_rsv_enriched"
  out$specificity_class[estimable & out$candidate_fdr < 0.05 & !out$direction_matches_rsv_vs_healthy] <- "relative_rsv_depleted_or_context_dependent"
  out$specificity_class[estimable & out$candidate_fdr >= 0.05] <- "not_distinguished_from_other_viruses_at_current_precision"

  root <- file.path("results", "single_gene_upgrade_v1", "specificity")
  write_upgrade_tsv(out, file.path(root, "GSE155925_86_gene_rsv_vs_other_virus.tsv"))
  classification <- out[, c("gene_id", "gene_symbol", "three_cohort_meta_effect", "log2FC", "SE", "ci_low", "ci_high", "candidate_fdr", "direction_matches_rsv_vs_healthy", "specificity_class")]
  write_upgrade_tsv(classification, file.path(root, "cross_context_gene_classification.tsv"))
  summary <- as.data.frame(table(out$specificity_class), stringsAsFactors = FALSE)
  names(summary) <- c("specificity_class", "n_genes")
  summary$analysis_label <- "post_outcome_exploratory_multicohort_upgrade_v1"
  summary$n_rsv <- 31L
  summary$n_other_single_virus <- 17L
  write_upgrade_tsv(summary, file.path(root, "specificity_summary.tsv"))
  message("SPECIFICITY mapped=", sum(estimable), " candidate_BH<0.05=", sum(out$candidate_fdr < 0.05, na.rm = TRUE))
  invisible(list(effects = out, summary = summary))
}

if (sys.nframe() == 0L) main_gse155925_specificity()

