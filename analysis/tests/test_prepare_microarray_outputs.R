#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

expected <- list(
  GSE105450 = list(frozen = 122L, excluded = 0L, retained = 122L, genes = 13715L, targets = 391L, excluded_ids = character()),
  GSE103842 = list(frozen = 74L, excluded = 1L, retained = 73L, genes = 12475L, targets = 374L, excluded_ids = "GSM2784347")
)
forbidden <- c("case_status", "condition", "status", "severity", "pathogen", "coinfection", "treatment")

summary <- read.delim(
  file.path("results", "tables", "confirmatory_qc_mapping_summary.tsv"),
  sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE
)
stopifnot(identical(summary$cohort, names(expected)))

for (cohort in names(expected)) {
  wanted <- expected[[cohort]]
  row <- summary[summary$cohort == cohort, , drop = FALSE]
  qc <- read.delim(file.path("results", "tables", paste0(cohort, "_qc.tsv")), sep = "\t", quote = "", check.names = FALSE)
  coverage <- read.delim(file.path("results", "tables", paste0(cohort, "_target_coverage.tsv")), sep = "\t", quote = "", check.names = FALSE)
  derived <- readRDS(file.path("data", "derived", paste0(cohort, "_expression.rds")))
  stopifnot(
    nrow(qc) == wanted$frozen,
    sum(qc$exclude_qc) == wanted$excluded,
    identical(qc$sample_id[qc$exclude_qc], wanted$excluded_ids),
    all(qc$mad_constant == 1.4826),
    all(qc$mad_rule == "two_sided_abs_x_minus_median_gt_5_scaled_MAD"),
    all(qc$nonfinite_policy == "stop_entire_cohort"),
    ncol(derived$expression) == wanted$retained,
    nrow(derived$expression) == wanted$genes,
    identical(colnames(derived$expression), derived$sample_data$sample_id),
    !any(names(derived$sample_data) %in% forbidden),
    identical(names(derived$sample_data), c("sample_id", "platform", "technical_batch")),
    !anyDuplicated(rownames(derived$expression)),
    all(grepl("^[0-9]+$", rownames(derived$expression))),
    nrow(coverage) == 520L,
    sum(coverage$detectable) == wanted$targets,
    unique(coverage$interpretability) == "interpretable",
    isFALSE(derived$preprocessing$between_array_normalization_applied),
    isFALSE(derived$preprocessing$phenotype_fields_available_in_qc_object),
    identical(derived$preprocessing$qc_mad_constant, 1.4826),
    identical(derived$preprocessing$qc_mad_rule, "two_sided_abs_x_minus_median_gt_5_scaled_MAD"),
    grepl("No RAW/CEL/FASTQ/SRA read", derived$preprocessing$processed_only_boundary, fixed = TRUE),
    row$samples_frozen == wanted$frozen,
    row$samples_retained == wanted$retained,
    row$primary_targets_detectable == wanted$targets
  )
}

annotation_path <- file.path("data", "raw", "GPL10558", "GPL10558.annot.gz")
annotation_sha <- system2("shasum", c("-a", "256", annotation_path), stdout = TRUE)
annotation_sha <- sub("[[:space:]].*$", "", annotation_sha[[1L]])
stopifnot(
  as.numeric(file.info(annotation_path)$size) == 7290886,
  annotation_sha == "c914fdbe1130906ce3b9c97f5a75280591c89c617c168fb687c641f683e76b45"
)

checksum <- read.delim(file.path("logs", "checksums", "task5_sha256.tsv"), sep = "\t", quote = "", check.names = FALSE)
required_provenance <- c(
  file.path("analysis", "R", "04_prepare_microarray.R"),
  file.path("analysis", "tests", "test_prepare_microarray.R"),
  file.path("analysis", "tests", "test_prepare_microarray_outputs.R"),
  file.path("analysis", "tests", "integration_prepare_microarray_determinism.R"),
  file.path("docs", "science-superpowers", "preregistrations", "amendments", "2026-08-05-processed-qc-detection-filter.md"),
  file.path("docs", "science-superpowers", "preregistrations", "amendments", "2026-08-05-scaled-mad-blind-qc.md"),
  file.path("logs", "session_info", "task5_supersession.md")
)
stopifnot(all(required_provenance %in% checksum$artifact))
for (i in seq_len(nrow(checksum))) {
  stopifnot(file.exists(checksum$artifact[[i]]))
  actual <- system2("shasum", c("-a", "256", checksum$artifact[[i]]), stdout = TRUE)
  actual <- sub("[[:space:]].*$", "", actual[[1L]])
  stopifnot(identical(actual, checksum$sha256[[i]]))
}

cat("task5 output verification passed\n")
