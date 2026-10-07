#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("yaml", quietly = TRUE)) stop("Package yaml is required", call. = FALSE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
base_config <- yaml::read_yaml(file.path(project_root, "analysis/config/cohorts.yml"))
accessions <- vapply(base_config$cohorts, `[[`, character(1), "accession")
cohort <- base_config$cohorts[[match("GSE105450", accessions)]]
test_root <- tempfile("geo_processed_only_integration_")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

run_main <- function(config, label) {
  case_root <- file.path(test_root, label)
  dir.create(case_root, recursive = TRUE, showWarnings = FALSE)
  config_path <- file.path(case_root, "cohorts.yml")
  yaml::write_yaml(config, config_path)
  output_path <- file.path(case_root, "run.log")
  expression <- paste0(
    "source(\"renv/activate.R\"); ",
    "source(\"analysis/R/02_download_geo_metadata.R\")"
  )
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("-e", shQuote(expression)),
    stdout = output_path,
    stderr = output_path,
    env = c(
      "R_PROFILE_USER=/dev/null",
      "R_ENVIRON_USER=/dev/null",
      "RENV_CONFIG_SANDBOX_ENABLED=FALSE",
      "RENV_CONFIG_NAMESPACES_CHECK=FALSE",
      paste0("GEO_COHORTS_CONFIG=", shQuote(config_path)),
      paste0("GEO_DOWNLOAD_POLICY=", shQuote(file.path(project_root, "analysis/R/geo_download_policy.R")))
    )
  )
  list(status = status, output = readLines(output_path, warn = FALSE))
}

make_case_config <- function(label) {
  case_root <- file.path(test_root, label)
  config <- base_config
  config$cohorts <- list(cohort)
  config$raw_root <- file.path(case_root, "raw")
  config$metadata_output <- file.path(case_root, "metadata", "all_samples_raw.tsv")
  config$checksum_output <- file.path(case_root, "checksums", "geo_sha256.tsv")
  config
}

processed_label <- "processed_without_raw"
processed_config <- make_case_config(processed_label)
processed_dir <- file.path(processed_config$raw_root, "GSE105450")
dir.create(processed_dir, recursive = TRUE)
copy_names <- vapply(
  Filter(function(x) x$role != "raw_archive", cohort$downloads),
  `[[`,
  character(1),
  "filename"
)
stopifnot(all(file.copy(
  file.path(base_config$raw_root, "GSE105450", copy_names),
  file.path(processed_dir, copy_names)
)))
processed_run <- run_main(processed_config, processed_label)
stopifnot(processed_run$status == 0L)
processed_manifest <- read.delim(
  processed_config$checksum_output,
  sep = "\t",
  quote = "",
  check.names = FALSE
)
stopifnot(nrow(processed_manifest) == 5L)
raw_row <- processed_manifest[processed_manifest$file_role == "raw_archive", , drop = FALSE]
stopifnot(
  nrow(raw_row) == 1L,
  raw_row$action == "optional_source_absent",
  raw_row$bytes == 0,
  raw_row$sha256 == ""
)
stopifnot(nrow(read.delim(processed_config$metadata_output, sep = "\t", quote = "\"")) == 127L)

missing_label <- "missing_processed_matrix"
missing_config <- make_case_config(missing_label)
missing_dir <- file.path(missing_config$raw_root, "GSE105450")
dir.create(missing_dir, recursive = TRUE)
metadata_names <- vapply(
  Filter(function(x) x$role %in% c("official_quick_metadata", "official_supplement_filelist"), cohort$downloads),
  `[[`,
  character(1),
  "filename"
)
stopifnot(all(file.copy(
  file.path(base_config$raw_root, "GSE105450", metadata_names),
  file.path(missing_dir, metadata_names)
)))
files_before <- sort(list.files(missing_dir, recursive = TRUE, all.files = TRUE))
missing_run <- run_main(missing_config, missing_label)
stopifnot(missing_run$status != 0L)
stopifnot(any(grepl("No directly analyzable processed matrix is present", missing_run$output, fixed = TRUE)))
stopifnot(identical(files_before, sort(list.files(missing_dir, recursive = TRUE, all.files = TRUE))))
stopifnot(!file.exists(missing_config$checksum_output))

cat("processed_without_raw=exit0 optional_source_absent recorded\n")
cat("missing_processed_matrix=nonzero no_files_added\n")
