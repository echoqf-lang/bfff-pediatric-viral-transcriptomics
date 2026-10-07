#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

read_tsv <- function(path) read.delim(
  path, sep = "\t", quote = "", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

base <- file.path("results", "sensitivity")
required <- file.path(base, c(
  "all_sensitivity_results.tsv",
  "camera_and_stouffer_results.tsv",
  "target_coverage.tsv",
  "cohort_model_diagnostics.tsv",
  "cohort_qc_audit.tsv",
  "gene_effects_GSE105450_hospital.tsv",
  "gene_effects_GSE188427.tsv",
  "gene_effects_GSE103119.tsv",
  "gene_effects_GSE155925.tsv",
  "meta_hospital_reml_hksj.tsv",
  "meta_replacement_reml_hksj.tsv",
  "rnaseq_target_selection.tsv",
  "rnaseq_alias_audit.tsv",
  "not_implemented.tsv",
  "input_access_audit.tsv"
))
stopifnot(all(file.exists(required)))

set_results <- read_tsv(file.path(base, "camera_and_stouffer_results.tsv"))
stopifnot(
  nrow(set_results) == 11L,
  setequal(unique(set_results$family), c("S1", "S2", "E1")),
  all(table(set_results$family) == c(E1 = 4L, S1 = 4L, S2 = 3L)),
  all(is.finite(set_results$p_value)),
  all(is.finite(set_results$p_holm)),
  all(set_results$p_holm >= set_results$p_value),
  all(set_results$primary_H1_rescue_allowed == FALSE)
)
for (family in c("S1", "S2", "E1")) {
  x <- set_results[set_results$family == family, , drop = FALSE]
  stopifnot(isTRUE(all.equal(x$p_holm, p.adjust(x$p_value, method = "holm"), tolerance = 1e-14)))
}

coverage <- read_tsv(file.path(base, "target_coverage.tsv"))
stopifnot(
  all(coverage$n_detectable == coverage$camera_index_length),
  all(coverage$n_detectable >= 10L),
  all(coverage$coverage_fraction >= 0.50),
  all(coverage$interpretability == "interpretable")
)

diagnostics <- read_tsv(file.path(base, "cohort_model_diagnostics.tsv"))
stopifnot(
  all(diagnostics$design_rank == diagnostics$design_columns),
  all(diagnostics$primary_H1_rescue_allowed == FALSE),
  diagnostics$n_model[diagnostics$cohort == "GSE105450_hospital"] == 89L,
  diagnostics$n_model[diagnostics$cohort == "GSE188427"] <= 198L,
  diagnostics$n_model[diagnostics$cohort == "GSE103119"] <= 31L,
  diagnostics$n_model[diagnostics$cohort == "GSE155925"] == 48L
)

selection <- read_tsv(file.path(base, "rnaseq_target_selection.tsv"))
stopifnot(!anyDuplicated(selection$entrez_id), !anyDuplicated(selection$ensembl_id))
aliases <- read_tsv(file.path(base, "rnaseq_alias_audit.tsv"))
stopifnot(setequal(
  aliases$ensembl_id[aliases$mapping_status == "globally_ambiguous_removed"],
  c("ENSG00000223572", "ENSG00000237289")
))

not_implemented <- read_tsv(file.path(base, "not_implemented.tsv"))
stopifnot(
  setequal(not_implemented$analysis, c("CBC_adjustment", "cell_proportion_adjustment", "GSE155925_severity")),
  all(not_implemented$status == "not_implemented_by_frozen_contract")
)

access <- read_tsv(file.path(base, "input_access_audit.tsv"))
stopifnot(
  !any(grepl("RAW[.]tar|[.]CEL|FASTQ|SRA", access$path, ignore.case = TRUE)),
  all(access$network_download == FALSE)
)

decision <- read_tsv(file.path("results", "cohort", "target_set_primary_decision_partial.tsv"))
stopifnot(identical(decision$h1_supported, FALSE))

checksum <- read_tsv(file.path("logs", "checksums", "task11_sha256.tsv"))
stopifnot(all(file.exists(checksum$artifact)))
observed <- vapply(checksum$artifact, function(path) {
  sub("[[:space:]].*$", "", system2("shasum", c("-a", "256", path), stdout = TRUE)[[1L]])
}, character(1))
stopifnot(identical(unname(observed), unname(checksum$sha256)))

cat("Task 11 sensitivity-analysis output tests passed\n")
