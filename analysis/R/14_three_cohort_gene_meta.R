#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
source(file.path("analysis", "R", "single_gene_upgrade_pipeline.R"))

extract_discovery_effects <- function(frozen) {
  detail <- utils::read.csv(file.path("validation", "validated_targets_detail.csv"), row.names = 1L, check.names = FALSE)
  deg <- utils::read.csv(file.path("GEO_data", "GSE38900_DEG_v2.csv"), row.names = 1L, check.names = FALSE)
  idx <- match(frozen$gene_symbol, rownames(detail))
  idx_deg <- match(frozen$gene_symbol, rownames(deg))
  out <- data.frame(
    gene_id = frozen$entrez_id,
    gene_symbol = frozen$gene_symbol,
    cohort = "GSE38900",
    direction = "RSV_minus_healthy",
    log2FC = detail$logFC[idx],
    SE = abs(deg$logFC[idx_deg] / deg$t[idx_deg]),
    p_value = deg$P.Value[idx_deg],
    fdr_bh = deg$adj.P.Val[idx_deg],
    stringsAsFactors = FALSE
  )
  assert_unique_gene_effects(out, "GSE38900")
  out
}

extract_existing_effects <- function(path, cohort, frozen) {
  all <- read_upgrade_tsv(path)
  all$gene_id <- as.character(all$gene_id)
  idx <- match(frozen$entrez_id, all$gene_id)
  out <- data.frame(
    gene_id = frozen$entrez_id,
    gene_symbol = frozen$gene_symbol,
    cohort = cohort,
    direction = "RSV_minus_healthy",
    log2FC = all$log2FC[idx],
    SE = all$SE[idx],
    p_value = all$p_value[idx],
    fdr_bh = all$fdr_bh[idx],
    stringsAsFactors = FALSE
  )
  out
}

main_three_cohort_meta <- function() {
  frozen <- read_upgrade_tsv(file.path("data", "clean", "bfff_86_gene_set_frozen.tsv"), col_classes = "character")
  g38900 <- extract_discovery_effects(frozen)
  g77087 <- extract_existing_effects(
    file.path("results", "revised_main_gse77087", "cohort", "GSE77087_gene_effects.tsv"),
    "GSE77087", frozen
  )
  g103842 <- extract_existing_effects(
    file.path("results", "revised_main_gse77087", "cohort", "GSE103842_gene_effects.tsv"),
    "GSE103842", frozen
  )

  coverage <- mean(stats::complete.cases(g103842[c("log2FC", "SE")]))
  if (coverage < 0.70) stop("insufficient_mapping_coverage", call. = FALSE)
  for (x in list(g77087[complete.cases(g77087[c("log2FC", "SE")]), ], g103842[complete.cases(g103842[c("log2FC", "SE")]), ])) {
    assert_unique_gene_effects(x, unique(x$cohort))
  }

  root <- file.path("results", "single_gene_upgrade_v1")
  write_upgrade_tsv(g38900, file.path(root, "cohort", "GSE38900_86_gene_effects.tsv"))
  write_upgrade_tsv(g77087, file.path(root, "cohort", "GSE77087_86_gene_effects.tsv"))
  write_upgrade_tsv(g103842, file.path(root, "cohort", "GSE103842_86_gene_effects.tsv"))

  complete_ids <- Reduce(intersect, lapply(list(g38900, g77087, g103842), function(x) as.character(x$gene_id[complete.cases(x[c("log2FC", "SE")])])))
  combined <- do.call(rbind, lapply(list(g38900, g77087, g103842), function(x) x[x$gene_id %in% complete_ids, c("gene_id", "cohort", "direction", "log2FC", "SE")]))
  meta <- fit_hksj_meta(combined, c("GSE38900", "GSE77087", "GSE103842"))
  meta$gene_symbol <- frozen$gene_symbol[match(meta$gene_id, frozen$entrez_id)]
  meta <- meta[, c("gene_id", "gene_symbol", setdiff(names(meta), c("gene_id", "gene_symbol")))]
  write_upgrade_tsv(meta, file.path(root, "meta", "three_cohort_86_gene_reml_hksj.tsv"))

  same_389_770 <- sign(g38900$log2FC) == sign(g77087$log2FC)
  g77087$candidate_fdr <- bh_family(g77087$p_value)
  summary <- data.frame(
    analysis_label = "post_outcome_exploratory_multicohort_upgrade_v1",
    frozen_n = 86L,
    gse77087_mapped_n = sum(complete.cases(g77087[c("log2FC", "SE")])),
    gse103842_mapped_n = sum(complete.cases(g103842[c("log2FC", "SE")])),
    meta_complete_n = nrow(meta),
    gse38900_gse77087_direction_concordant_n = sum(same_389_770, na.rm = TRUE),
    gse77087_candidate_bh_and_direction_n = sum(same_389_770 & g77087$candidate_fdr < 0.05, na.rm = TRUE),
    gse77087_genome_fdr_and_direction_n = sum(same_389_770 & g77087$fdr_bh < 0.05, na.rm = TRUE),
    gse77087_genome_fdr_effect_and_direction_n = sum(same_389_770 & g77087$fdr_bh < 0.05 & abs(g77087$log2FC) > 0.5, na.rm = TRUE),
    all_three_same_direction_n = sum(meta$direction_pattern %in% c("all_up", "all_down")),
    meta_bh_lt_0_05_n = sum(meta$fdr_bh < 0.05),
    meta_adhoc_bh_lt_0_05_n = sum(meta$adhoc_fdr_bh < 0.05),
    stringsAsFactors = FALSE
  )
  write_upgrade_tsv(summary, file.path(root, "meta", "three_cohort_meta_summary.tsv"))
  message("META complete=", nrow(meta), " all_three_same_direction=", summary$all_three_same_direction_n,
          " HKSJ_BH<0.05=", summary$meta_bh_lt_0_05_n,
          " modified_HKSJ_BH<0.05=", summary$meta_adhoc_bh_lt_0_05_n)
  invisible(list(meta = meta, summary = summary))
}

if (sys.nframe() == 0L) main_three_cohort_meta()
