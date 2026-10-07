#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

fail <- function(message) stop(message, call. = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) fail("Unable to determine validator script location")
project_root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", script_arg))), "..", ".."))
project_path <- function(...) file.path(project_root, ...)

read_tsv_required <- function(path, required_columns) {
  if (!file.exists(path)) fail(paste0("Required input is absent: ", path))
  x <- read.delim(path, sep = "\t", check.names = FALSE, na.strings = "", stringsAsFactors = FALSE)
  missing <- setdiff(required_columns, names(x))
  if (length(missing)) fail(paste0("Invalid schema in ", path, ": missing ", paste(missing, collapse = ", ")))
  x
}

read_csv_required <- function(path, required_columns) {
  if (!file.exists(path)) fail(paste0("Required input is absent: ", path))
  x <- read.csv(path, check.names = FALSE, na.strings = "", stringsAsFactors = FALSE)
  missing <- setdiff(required_columns, names(x))
  if (length(missing)) fail(paste0("Invalid schema in ", path, ": missing ", paste(missing, collapse = ", ")))
  x
}

assert_equal <- function(observed, expected, label) {
  if (!identical(as.integer(observed), as.integer(expected))) {
    fail(sprintf("%s: expected %d, observed %d", label, expected, observed))
  }
}

frozen_path <- project_path("data", "clean", "bfff_86_gene_set_frozen.tsv")
summary_path <- project_path("results", "bfff_full_86_evidence_matrix_v1", "evidence_summary.tsv")
matrix_path <- project_path("results", "bfff_full_86_evidence_matrix_v1", "Table_full_86_cross_context_evidence.tsv")
cohort_path <- project_path("manuscript", "v6_bmc_submission_assets", "tables", "Table_1_cohort_characteristics.tsv")
manuscript_path <- project_path("manuscript", "draft_v6_bmc_genomics.md")
audit_path <- project_path("results", "tables", "target_provenance_audit.tsv")
drug_target_path <- project_path("manuscript", "supplementary", "Table_S1_drug_compound_targets.csv")
disease_target_path <- project_path("manuscript", "supplementary", "Table_S2_disease_targets.csv")
intersection_path <- project_path("manuscript", "supplementary", "Table_S3_intersection_targets.csv")

frozen <- read_tsv_required(frozen_path, c("entrez_id", "gene_symbol"))
matrix <- read_tsv_required(matrix_path, c("entrez_id", "gene_symbol", "blood_GSE77087_genome_fdr_support", "airway_any_fdr_support"))
summary <- read_tsv_required(summary_path, c("metric", "value"))
cohorts <- read_tsv_required(cohort_path, "GEO_accession")
audit <- read_tsv_required(audit_path, "herb_cn")
drug_targets <- read_csv_required(drug_target_path, "gene")
disease_targets <- read_csv_required(disease_target_path, "gene_symbol")
intersection_targets <- read_csv_required(intersection_path, "gene_symbol")

unique_required_values <- function(x, label) {
  values <- trimws(as.character(x))
  if (anyNA(values) || any(!nzchar(values))) {
    fail(paste0(label, " contains empty or missing values"))
  }
  unique(values)
}

provenance_sets <- list(
  herbs = unique_required_values(audit$herb_cn, "Target provenance herb_cn"),
  putative = unique_required_values(drug_targets$gene, "Table S1 gene"),
  disease = unique_required_values(disease_targets$gene_symbol, "Table S2 gene_symbol"),
  search_space = unique_required_values(intersection_targets$gene_symbol, "Table S3 gene_symbol")
)
provenance_counts <- c(
  herbs = length(provenance_sets$herbs),
  putative = length(provenance_sets$putative),
  disease = length(provenance_sets$disease),
  search_space = length(provenance_sets$search_space)
)
assert_equal(provenance_counts[["herbs"]], 14L, "Unique herb count")
assert_equal(provenance_counts[["putative"]], 2092L, "Unique putative-target count")
assert_equal(provenance_counts[["disease"]], 3053L, "Unique disease-target count")
assert_equal(provenance_counts[["search_space"]], 734L, "Unique intersection-target count")

expected_intersection <- intersect(provenance_sets$putative, provenance_sets$disease)
if (!setequal(provenance_sets$search_space, expected_intersection)) {
  fail("Table S3 gene_symbol set is not the exact intersection of Table S1 gene and Table S2 gene_symbol")
}

validate_gene_keys <- function(x, label) {
  identifiers <- lapply(c("entrez_id", "gene_symbol"), function(column) trimws(as.character(x[[column]])))
  if (any(vapply(identifiers, function(value) anyNA(value) || any(!nzchar(value)), logical(1)))) {
    fail(paste0(label, " contains empty or missing Entrez IDs or gene symbols"))
  }
  if (anyDuplicated(identifiers[[1L]]) || anyDuplicated(identifiers[[2L]])) {
    fail(paste0(label, " contains duplicate Entrez IDs or gene symbols"))
  }
  key <- paste(identifiers[[1L]], identifiers[[2L]], sep = "\r")
  if (anyDuplicated(key)) fail(paste0(label, " contains duplicate (entrez_id, gene_symbol) keys"))
  key
}
frozen_keys <- validate_gene_keys(frozen, "Frozen 86-gene set")
matrix_keys <- validate_gene_keys(matrix, "Evidence matrix")
assert_equal(nrow(frozen), 86L, "Frozen gene count")
assert_equal(nrow(matrix), 86L, "Evidence matrix row count")
if (!setequal(frozen_keys, matrix_keys)) fail("Frozen gene set and evidence matrix (entrez_id, gene_symbol) keys differ")

summary_value <- function(metric) {
  value <- summary$value[summary$metric == metric]
  numeric_value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(numeric_value) || numeric_value != floor(numeric_value)) {
    fail(paste0("Invalid finite integer summary metric: ", metric))
  }
  as.integer(numeric_value)
}
assert_equal(summary_value("any_airway_fdr_support"), 24L, "Summary airway support")
assert_equal(summary_value("blood_GSE77087_and_any_airway_support"), 16L, "Summary blood-to-airway support")

to_flag <- function(x, label) {
  if (is.logical(x)) return(x)
  normalized <- toupper(trimws(as.character(x)))
  if (any(!is.na(normalized) & !normalized %in% c("TRUE", "FALSE"))) fail(paste0("Invalid boolean values in ", label))
  normalized == "TRUE"
}
airway_support <- to_flag(matrix$airway_any_fdr_support, "airway_any_fdr_support")
blood_support <- to_flag(matrix$blood_GSE77087_genome_fdr_support, "blood_GSE77087_genome_fdr_support")
missing_support <- is.na(airway_support) | is.na(blood_support)
if (any(missing_support)) {
  missing_genes <- as.character(matrix$gene_symbol[missing_support])
  cat(sprintf("FIGURE_INPUT_CONTRACT_NOTE missing_support_genes=%d genes=%s; NA treated as unsupported\n",
              sum(missing_support), paste(missing_genes, collapse = ",")))
  airway_support[is.na(airway_support)] <- FALSE
  blood_support[is.na(blood_support)] <- FALSE
}
airway <- sum(airway_support)
bridge <- sum(blood_support & airway_support)
assert_equal(airway, 24L, "Recomputed airway support")
assert_equal(bridge, 16L, "Recomputed blood-to-airway support")

required_cohorts <- c("GSE38900", "GSE77087", "GSE103842", "GSE97742", "GSE41374", "GSE155925")
observed_cohorts <- as.character(cohorts$GEO_accession)
if (anyDuplicated(observed_cohorts) || !setequal(observed_cohorts, required_cohorts)) {
  fail("Cohort characteristics table does not contain exactly the six required GSE cohorts")
}

if (!file.exists(manuscript_path)) fail(paste0("Required input is absent: ", manuscript_path))
manuscript <- paste(readLines(manuscript_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
manuscript_count_labels <- c(
  sprintf("%d herbs", provenance_counts[["herbs"]]),
  sprintf("%s putative formula targets", format(provenance_counts[["putative"]], big.mark = ",")),
  sprintf("%s disease-associated database genes", format(provenance_counts[["disease"]], big.mark = ",")),
  sprintf("%s-gene search space", format(provenance_counts[["search_space"]], big.mark = ","))
)
manuscript_count_present <- vapply(
  manuscript_count_labels,
  function(label) grepl(tolower(label), tolower(manuscript), fixed = TRUE),
  logical(1)
)
if (!all(manuscript_count_present)) {
  fail(paste0(
    "Manuscript provenance counts are inconsistent or absent: ",
    paste(manuscript_count_labels[!manuscript_count_present], collapse = "; ")
  ))
}
if (!grepl("no(\\s+analyzed)?\\s+cohort included bfff exposure", manuscript, ignore.case = TRUE, perl = TRUE)) {
  fail("Manuscript does not state that no cohort included BFFF exposure")
}

cat("FIGURE_INPUT_CONTRACT_PASS frozen=86 airway=24 bridge=16 search_space=734 cohorts=6\n")
