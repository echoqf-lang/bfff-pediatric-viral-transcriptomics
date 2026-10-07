#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
}

expected <- list(
  GSE105450 = list(n = 122L, healthy = 33L, RSV = 89L, genes = 13715L, batches = 3L, rank = 6L),
  GSE103842 = list(n = 73L, healthy = 12L, RSV = 61L, genes = 12475L, batches = 10L, rank = 13L)
)
expected_schema <- c(
  "gene_id", "cohort", "direction", "log2FC", "SE", "moderated_t",
  "p_value", "fdr_bh", "average_expression"
)

for (cohort in names(expected)) {
  wanted <- expected[[cohort]]
  effects_path <- file.path("results", "cohort", paste0(cohort, "_gene_effects.tsv"))
  diagnostics_path <- file.path("results", "cohort", paste0(cohort, "_model_diagnostics.tsv"))
  effects <- read_tsv(effects_path)
  diagnostics <- read_tsv(diagnostics_path)
  stopifnot(
    identical(names(effects), expected_schema),
    nrow(effects) == wanted$genes,
    !anyDuplicated(effects$gene_id),
    all(effects$cohort == cohort),
    all(effects$direction == "RSV_minus_healthy"),
    all(is.finite(effects$log2FC)),
    all(is.finite(effects$SE) & effects$SE > 0),
    all(is.finite(effects$moderated_t)),
    all(is.finite(effects$p_value) & effects$p_value >= 0 & effects$p_value <= 1),
    all(is.finite(effects$fdr_bh) & effects$fdr_bh >= 0 & effects$fdr_bh <= 1),
    isTRUE(all.equal(effects$moderated_t, effects$log2FC / effects$SE, tolerance = 1e-10)),
    isTRUE(all.equal(effects$fdr_bh, p.adjust(effects$p_value, method = "BH"), tolerance = 1e-15)),
    nrow(diagnostics) == 1L,
    diagnostics$n_rds_samples == wanted$n,
    diagnostics$n_complete_cases == wanted$n,
    diagnostics$n_excluded_missing_covariate == 0L,
    diagnostics$n_healthy == wanted$healthy,
    diagnostics$n_RSV == wanted$RSV,
    diagnostics$technical_batch_included,
    diagnostics$technical_batch_n_levels == wanted$batches,
    diagnostics$design_rank == wanted$rank,
    diagnostics$design_columns == wanted$rank,
    diagnostics$design_full_rank,
    diagnostics$genes_tested_unfiltered == wanted$genes,
    !diagnostics$severity_in_primary_model
  )
}

checksum <- read_tsv(file.path("logs", "checksums", "task7_sha256.tsv"))
stopifnot(!anyDuplicated(checksum$artifact))
for (i in seq_len(nrow(checksum))) {
  stopifnot(file.exists(checksum$artifact[[i]]))
  actual <- system2("shasum", c("-a", "256", checksum$artifact[[i]]), stdout = TRUE)
  actual <- sub("[[:space:]].*$", "", actual[[1L]])
  stopifnot(identical(actual, checksum$sha256[[i]]))
}

script <- paste(readLines(file.path("analysis", "R", "05_fit_cohort_models.R"), warn = FALSE), collapse = "\n")
forbidden_calls <- c("run_camera_analysis\\(", "test_camera_target_set\\(", "meta_analyze_gene_effects\\(", "fgsea::")
stopifnot(!any(vapply(forbidden_calls, grepl, logical(1), x = script, perl = TRUE)))

cat("task7 real-output verification passed\n")
