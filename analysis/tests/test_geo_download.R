library(testthat)
library(yaml)

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    marker <- file.path(current, "analysis", "config", "cohorts.yml")
    if (file.exists(marker)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Cannot locate project root containing analysis/config/cohorts.yml", call. = FALSE)
    }
    current <- parent
  }
}

project_root <- find_project_root()
config <- read_yaml(file.path(project_root, "analysis/config/cohorts.yml"))
source(file.path(project_root, "analysis/R/geo_download_policy.R"), local = TRUE)
cohorts <- config$cohorts
accessions <- vapply(cohorts, `[[`, character(1), "accession")
roles <- vapply(cohorts, `[[`, character(1), "analysis_role")

test_that("project root discovery does not depend on a testthat package layout", {
  expect_identical(find_project_root(project_root), project_root)
  expect_identical(find_project_root(file.path(project_root, "analysis", "tests")), project_root)
})

test_that("confirmatory and substitution roles are frozen and mutually exclusive", {
  expect_setequal(accessions[roles == "confirmatory"], c("GSE105450", "GSE103842"))
  expect_false("GSE188427" %in% accessions[roles == "confirmatory"])
  substitute <- cohorts[[match("GSE188427", accessions)]]
  expect_identical(substitute$analysis_role, "substitution_sensitivity")
  expect_identical(substitute$replaces, "GSE105450")
  expect_false(all(c("GSE105450", "GSE188427") %in% accessions[roles == "confirmatory"]))
})

test_that("exploratory evidence boundaries are explicit", {
  observed <- cohorts[[match("GSE38900", accessions)]]
  expect_identical(observed$analysis_role, "exploratory_previously_observed")
  expect_true(observed$outcome_previously_observed)
  expect_identical(cohorts[[match("GSE155925", accessions)]]$analysis_role, "exploratory")
  expect_identical(cohorts[[match("GSE103119", accessions)]]$analysis_role, "exploratory")
})

test_that("download destinations and URLs are unique", {
  keys <- unlist(lapply(cohorts, function(x) {
    vapply(x$downloads, function(y) paste(x$accession, y$filename, sep = "/"), character(1))
  }))
  urls <- unlist(lapply(cohorts, function(x) vapply(x$downloads, `[[`, character(1), "url")))
  expect_equal(anyDuplicated(keys), 0L)
  expect_equal(anyDuplicated(urls), 0L)
  expect_true(all(startsWith(urls, "https://")))
  expect_true(all(grepl("^'.*&.*'$", shQuote(urls[grepl("&", urls)], type = "sh"))))
})

test_that("core acquisition is implemented through actual GEOquery calls", {
  script <- readLines(file.path(project_root, "analysis/R/02_download_geo_metadata.R"), warn = FALSE)
  expect_true(any(grepl("GEOquery::getGEO\\(", script)))
  expect_true(any(grepl("GEOquery::getGEOSuppFiles\\(", script)))
  expect_true(any(grepl("GEOquery::getGEOfile\\(", script)))
  items <- unlist(lapply(cohorts, `[[`, "downloads"), recursive = FALSE)
  methods <- vapply(items, function(x) geo_acquisition_method_for_role(x$role), character(1))
  expect_equal(sum(startsWith(methods, "GEOquery::")), 19L)
  expect_equal(sum(methods == "curl_auxiliary_filelist"), 5L)
  expect_equal(sum(methods == "retain_existing_source_no_auto_download"), 5L)
  expect_true(all(vapply(items[methods == "curl_auxiliary_filelist"], function(x) {
    identical(x$role, "official_supplement_filelist")
  }, logical(1))))
})

test_that("lightweight policy forbids automatic source-file acquisition", {
  expect_identical(config$policies$default_network_mode, "reuse_only")
  expect_true(all(c("RAW.tar", "CEL", "CEL.gz", "FASTQ", "SRA") %in%
                    config$policies$automatic_download_forbidden_types))
  expect_false("raw_archive" %in% config$policies$automatic_download_allowed_roles)
  script <- readLines(file.path(project_root, "analysis/R/02_download_geo_metadata.R"), warn = FALSE)
  expect_true(any(grepl("retain_existing_source_no_auto_download", script, fixed = TRUE)))
  expect_true(any(grepl("GEOQUERY_ALLOW_PROCESSED_DOWNLOAD", script, fixed = TRUE)))
  fake_dir <- tempfile("geo_lightweight_policy_")
  dir.create(fake_dir)
  on.exit(unlink(fake_dir, recursive = TRUE, force = TRUE), add = TRUE)
  fake <- list(downloads = list(
    list(filename = "analysis_matrix.txt.gz", role = "series_matrix"),
    list(filename = "source_RAW.tar", role = "raw_archive")
  ))
  empty_status <- lightweight_availability(fake, fake_dir)
  expect_false(empty_status$analysis_matrix_present)
  expect_equal(length(empty_status$missing_required_non_source), 1L)
  expect_equal(length(empty_status$missing_retained_source), 1L)
  expect_true(file.create(file.path(fake_dir, "analysis_matrix.txt.gz")))
  ready_status <- lightweight_availability(fake, fake_dir)
  expect_true(ready_status$analysis_matrix_present)
  expect_equal(length(ready_status$missing_required_non_source), 0L)
  expect_equal(length(ready_status$missing_retained_source), 1L)
})

test_that("download manifest matches immutable files", {
  manifest_path <- file.path(project_root, config$checksum_output)
  skip_if_not(file.exists(manifest_path), "downloads not completed")
  manifest <- read.delim(manifest_path, sep = "\t", quote = "", check.names = FALSE)
  expect_false(file.exists(paste0(manifest_path, ".partial")))
  expect_equal(anyDuplicated(manifest$relative_path), 0L)
  expect_true(all(startsWith(
    manifest$current_script_acquisition_path[
      !manifest$file_role %in% c("official_supplement_filelist", "raw_archive")
    ],
    "GEOquery::"
  )))
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(project_root, manifest$relative_path[[i]])
    expect_true(file.exists(path), info = path)
    expect_equal(as.numeric(file.info(path)$size), as.numeric(manifest$bytes[[i]]), info = path)
    actual <- sub("[[:space:]].*$", "", system2("shasum", c("-a", "256", path), stdout = TRUE)[[1L]])
    expect_identical(actual, manifest$sha256[[i]], info = path)
  }
})

test_that("source cleanup list excludes processed analysis matrices", {
  cleanup_path <- file.path(project_root, "logs/checksums/source_file_cleanup_candidates.tsv")
  skip_if_not(file.exists(cleanup_path), "cleanup candidate list not generated")
  cleanup <- read.delim(cleanup_path, sep = "\t", quote = "", check.names = FALSE)
  expect_true(all(cleanup$cleanup_status == "candidate_only_not_deleted"))
  expect_true(all(cleanup$file_role == "raw_archive"))
  expect_false(any(grepl("non-normalized|Raw_counts_matrix", cleanup$relative_path, ignore.case = TRUE)))
  for (i in seq_len(nrow(cleanup))) {
    path <- file.path(project_root, cleanup$relative_path[[i]])
    actual <- sub("[[:space:]].*$", "", system2("shasum", c("-a", "256", path), stdout = TRUE)[[1L]])
    expect_identical(actual, cleanup$sha256[[i]])
  }
})

test_that("corrective GEOquery redownload matches every frozen file", {
  validation_path <- file.path(project_root, "logs/checksums/geoquery_redownload_validation.tsv")
  skip_if_not(file.exists(validation_path), "full GEOquery redownload validation not completed")
  validation <- read.delim(validation_path, sep = "\t", quote = "", check.names = FALSE)
  expect_equal(nrow(validation), 29L)
  expect_equal(anyDuplicated(paste(validation$accession, validation$filename)), 0L)
  expect_false(any(validation$comparison == "mismatch"))
  expect_true(all(validation$comparison %in% c("exact_bytes", "same_uncompressed_content")))
  expect_equal(sum(startsWith(validation$acquisition_function, "GEOquery::")), 24L)
  expect_equal(sum(validation$acquisition_function == "curl_auxiliary_filelist"), 5L)
})

test_that("raw tar archives match official filelist sizes and are readable", {
  for (cohort in cohorts) {
    cohort_dir <- file.path(project_root, config$raw_root, cohort$accession)
    filelist_path <- file.path(cohort_dir, "filelist.txt")
    if (!file.exists(filelist_path)) next
    filelist <- read.delim(filelist_path, sep = "\t", quote = "", check.names = FALSE)
    archive <- filelist[filelist[[1L]] == "Archive", , drop = FALSE]
    expect_equal(nrow(archive), 1L)
    archive_path <- file.path(cohort_dir, archive$Name[[1L]])
    expect_equal(as.numeric(file.info(archive_path)$size), as.numeric(archive$Size[[1L]]))
    expect_equal(system2("tar", c("-tf", shQuote(archive_path)), stdout = FALSE), 0L)
  }
})

test_that("metadata samples are unique and match frozen official counts/platforms", {
  metadata_path <- file.path(project_root, config$metadata_output)
  skip_if_not(file.exists(metadata_path), "metadata freeze not completed")
  metadata <- read.delim(metadata_path, sep = "\t", quote = "\"", check.names = FALSE)
  expect_equal(anyDuplicated(metadata$geo_accession), 0L)
  for (cohort in cohorts) {
    subset <- metadata[metadata$series_accession == cohort$accession, , drop = FALSE]
    expect_equal(nrow(subset), cohort$expected_sample_count, info = cohort$accession)
    expect_setequal(unique(subset$platform_id), cohort$expected_platforms)
    expect_true(all(subset$task3_inclusion_status == "all_official_samples_unfiltered"))
  }
})
