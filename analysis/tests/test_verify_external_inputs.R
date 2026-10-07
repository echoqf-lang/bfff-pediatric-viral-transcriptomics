#!/usr/bin/env Rscript

fail <- function(message) stop(message, call. = FALSE)
assert_true <- function(value, message) if (!isTRUE(value)) fail(message)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script <- file.path(project_root, "analysis/R/verify_external_inputs.R")
assert_true(file.exists(script), "Verifier script is missing")

expected_targets <- c(
  "1.TCMSP_data/ob大于等于30_dl大于等于0.18_filtered_targets_with_uniprot.csv" = "ef80d7b6f18444ad6ea643da90abe5fa0057eae0f8af9126d50250c04b78193e",
  "4.BATMAN_data/BATMAN_All_data.csv" = "4950e91dd067814ef42241bdc1186bc2c4e9e775b0aa0c80ef89aba323cd83cf",
  "herbs.txt" = "d1f6eca8f350bacabb4615a0db8f8b3e4efd7f0ef12d98aa5f6c40885cbe5405"
)
target_manifest <- read.delim(
  file.path(project_root, "logs/checksums/task2_sha256.tsv"),
  sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE
)
target_inputs <- target_manifest[target_manifest$role == "input", c("path", "sha256"), drop = FALSE]
actual_targets <- setNames(target_inputs$sha256, target_inputs$path)
assert_true(identical(actual_targets[names(expected_targets)], expected_targets), "Frozen target inputs changed")

gpl_manifest <- read.delim(
  file.path(project_root, "logs/checksums/gpl10558_annotation_sha256.tsv"),
  sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE
)
assert_true(nrow(gpl_manifest) == 1L, "GPL10558 manifest must contain one row")
assert_true(
  identical(gpl_manifest$official_url[[1L]], "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz") &&
    identical(gpl_manifest$relative_path[[1L]], "data/raw/GPL10558/GPL10558.annot.gz") &&
    identical(as.numeric(gpl_manifest$bytes[[1L]]), 7290886) &&
    identical(gpl_manifest$sha256[[1L]], "c914fdbe1130906ce3b9c97f5a75280591c89c617c168fb687c641f683e76b45"),
  "Frozen GPL10558 annotation identity changed"
)

fixture_root <- tempfile("verify-external-inputs-")
fixture_project <- file.path(fixture_root, "project")
fixture_external <- file.path(fixture_root, "external")
dir.create(file.path(fixture_project, "logs/checksums"), recursive = TRUE)
dir.create(file.path(fixture_project, "analysis/config"), recursive = TRUE)
dir.create(fixture_external, recursive = TRUE)
on.exit(unlink(fixture_root, recursive = TRUE, force = TRUE), add = TRUE)

write_fixture <- function(relative_path, text) {
  path <- file.path(fixture_external, relative_path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, path, useBytes = TRUE)
  path
}

sha256_file <- function(path) {
  shasum <- Sys.which("shasum")
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(shasum)) {
    output <- system2(shasum, c("-a", "256", shQuote(path)), stdout = TRUE)
  } else if (nzchar(sha256sum)) {
    output <- system2(sha256sum, shQuote(path), stdout = TRUE)
  } else {
    fail("Neither shasum nor sha256sum is available for the test")
  }
  substr(output[[1L]], 1L, 64L)
}

target_path <- write_fixture("private/targets.tsv", "target fixture")
geo_path <- write_fixture("data/raw/GSETEST/processed.txt.gz", "processed fixture")
gpl_path <- write_fixture("data/raw/GPLTEST/annot.gz", "annotation fixture")

write.table(
  data.frame(role = "input", path = "private/targets.tsv", sha256 = sha256_file(target_path)),
  file.path(fixture_project, "logs/checksums/task2_sha256.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
geo_url <- "https://example.org/processed.txt.gz"
write.table(
  data.frame(
    file_role = c("series_matrix", "raw_archive"),
    relative_path = c("data/raw/GSETEST/processed.txt.gz", "data/raw/GSETEST/never-read_RAW.tar"),
    official_url = c(geo_url, "https://example.org/never-read_RAW.tar"),
    bytes = c(file.info(geo_path)$size, 999L),
    sha256 = c(sha256_file(geo_path), paste(rep("0", 64L), collapse = ""))
  ),
  file.path(fixture_project, "logs/checksums/geo_sha256.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(
    platform = "GPLTEST",
    official_url = "https://example.org/annot.gz",
    relative_path = "data/raw/GPLTEST/annot.gz",
    bytes = file.info(gpl_path)$size,
    sha256 = sha256_file(gpl_path)
  ),
  file.path(fixture_project, "logs/checksums/gpl10558_annotation_sha256.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(c(
  'automatic_download_allowed_roles: ["official_quick_metadata", "series_matrix", "non_normalized_expression", "raw_counts", "official_supplement_filelist"]',
  'automatic_download_forbidden_types: ["RAW.tar", "CEL", "CEL.gz", "FASTQ", "SRA"]',
  paste0('url: "', geo_url, '"')
), file.path(fixture_project, "analysis/config/cohorts.yml"))

run_verifier <- function(extra_args = character()) {
  output_file <- tempfile("verify-output-")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script), extra_args),
    stdout = output_file,
    stderr = output_file,
    env = c(
      paste0("VERIFY_PROJECT_ROOT=", shQuote(fixture_project)),
      paste0("VERIFY_EXTERNAL_INPUT_ROOT=", shQuote(fixture_external))
    )
  )
  list(status = status, output = readLines(output_file, warn = FALSE))
}

before <- sort(list.files(fixture_external, recursive = TRUE, all.files = TRUE))
verified <- run_verifier()
after <- sort(list.files(fixture_external, recursive = TRUE, all.files = TRUE))
assert_true(identical(verified$status, 0L), paste(verified$output, collapse = "\n"))
assert_true(any(grepl("# read_only=TRUE", verified$output, fixed = TRUE)), "Read-only marker missing")
assert_true(any(grepl("# network_attempted=FALSE", verified$output, fixed = TRUE)), "Network marker missing")
assert_true(any(grepl("verified=3 failed=0", verified$output, fixed = TRUE)), "Success summary is wrong")
assert_true(!any(grepl("never-read_RAW.tar", verified$output, fixed = TRUE)), "raw_archive must be excluded")
assert_true(identical(before, after), "Verifier changed the external input tree")

listed <- run_verifier("--list")
assert_true(identical(listed$status, 0L), paste(listed$output, collapse = "\n"))
assert_true(any(grepl("not_checked", listed$output, fixed = TRUE)), "List mode did not report not_checked")

unlink(geo_path)
missing <- run_verifier()
assert_true(!identical(missing$status, 0L), "Missing processed input should fail")
assert_true(any(grepl("missing", missing$output, fixed = TRUE)), "Missing status was not printed")
assert_true(any(grepl(geo_url, missing$output, fixed = TRUE)), "Recovery URL was not printed")
assert_true(!file.exists(geo_path), "Verifier attempted to restore a missing input")

writeLines("changed target fixture", target_path, useBytes = TRUE)
mismatch <- run_verifier()
assert_true(!identical(mismatch$status, 0L), "Changed target input should fail")
assert_true(any(grepl("sha256_mismatch", mismatch$output, fixed = TRUE)), "SHA mismatch was not printed")

cat("PASS: read-only external-input verifier\n")
