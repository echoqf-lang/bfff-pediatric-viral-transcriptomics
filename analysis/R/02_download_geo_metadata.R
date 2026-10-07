#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

config_path <- Sys.getenv("GEO_COHORTS_CONFIG", unset = "analysis/config/cohorts.yml")
policy_path <- Sys.getenv("GEO_DOWNLOAD_POLICY", unset = "analysis/R/geo_download_policy.R")
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package yaml is required; restore the frozen renv environment", call. = FALSE)
}
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("Package GEOquery is required by the frozen environment", call. = FALSE)
}
if (!file.exists(config_path)) stop("Missing cohorts config", call. = FALSE)
if (!file.exists(policy_path)) stop("Missing GEO download policy helpers", call. = FALSE)
source(policy_path, local = TRUE)

config <- yaml::read_yaml(config_path)
raw_root <- config$raw_root
metadata_output <- config$metadata_output
checksum_output <- config$checksum_output
checksum_partial <- paste0(checksum_output, ".partial")
dir.create(raw_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(metadata_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(checksum_output), recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop(sprintf("SHA-256 failed for %s: %s", path, paste(out, collapse = " ")), call. = FALSE)
  }
  sub("[[:space:]].*$", "", out[[1L]])
}

read_previous_manifest <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.delim(path, sep = "\t", quote = "", check.names = FALSE)
}

previous_manifest <- read_previous_manifest(checksum_output)
if (is.null(previous_manifest)) previous_manifest <- read_previous_manifest(checksum_partial)

write_manifest_atomic <- function(rows, path) {
  manifest <- do.call(rbind, rows)
  temporary <- paste0(path, ".tmp")
  write.table(
    manifest,
    temporary,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  if (!file.rename(temporary, path)) stop(sprintf("Atomic manifest write failed: %s", path), call. = FALSE)
  invisible(manifest)
}

validate_existing <- function(path, url) {
  if (is.null(previous_manifest)) {
    stop(sprintf("Refusing to adopt existing unmanifested final file: %s", path), call. = FALSE)
  }
  relative_path <- sub(paste0("^", normalizePath(".", winslash = "/"), "/?"), "", normalizePath(path, winslash = "/"))
  row <- previous_manifest[
    previous_manifest$relative_path == relative_path & previous_manifest$official_url == url,
    , drop = FALSE
  ]
  if (nrow(row) != 1L) stop(sprintf("No unique prior manifest row for %s", path), call. = FALSE)
  actual_size <- as.numeric(file.info(path)$size)
  actual_sha <- sha256_file(path)
  if (actual_size != as.numeric(row$bytes) || actual_sha != row$sha256) {
    stop(sprintf("Immutable raw file changed: %s", path), call. = FALSE)
  }
  invisible(TRUE)
}

download_auxiliary_filelist <- function(url, path) {
  if (file.exists(path)) {
    validate_existing(path, url)
    return("reused_verified")
  }
  part_path <- paste0(path, ".part")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  args <- c(
    "-L", "--fail", "--show-error", "--retry", "8", "--retry-all-errors",
    "--retry-delay", "3", "--connect-timeout", "30", "--continue-at", "-",
    "--output", shQuote(part_path), shQuote(url)
  )
  status <- system2("curl", args)
  if (status != 0L) stop(sprintf("Download failed (%d): %s", status, url), call. = FALSE)
  if (!file.exists(part_path) || file.info(part_path)$size <= 0) {
    stop(sprintf("Downloaded file is empty: %s", url), call. = FALSE)
  }
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzip_status <- system2("gzip", c("-t", part_path))
    if (gzip_status != 0L) stop(sprintf("Gzip integrity check failed: %s", part_path), call. = FALSE)
  }
  if (!file.rename(part_path, path)) stop(sprintf("Atomic finalize failed: %s", path), call. = FALSE)
  Sys.chmod(path, mode = "0444", use_umask = TRUE)
  "downloaded"
}

network_download_authorized <- identical(
  Sys.getenv("GEOQUERY_ALLOW_PROCESSED_DOWNLOAD", unset = ""),
  "explicit_user_authorization"
)

acquire_quick_soft <- function(cohort, cohort_dir) {
  accession <- cohort$accession
  item <- Filter(function(x) identical(x$role, "official_quick_metadata"), cohort$downloads)[[1L]]
  final_path <- file.path(cohort_dir, item$filename)
  stage_dir <- tempfile(pattern = paste0(accession, "_quick_"))
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  geoquery_path <- file.path(stage_dir, paste0(accession, ".soft"))
  if (file.exists(final_path) && !file.copy(final_path, geoquery_path, overwrite = FALSE)) {
    stop(sprintf("Unable to prepare GEOquery quick-SOFT cache for %s", accession), call. = FALSE)
  }
  returned_path <- GEOquery::getGEOfile(accession, destdir = stage_dir, amount = "quick")
  if (!file.exists(returned_path)) stop(sprintf("GEOquery quick SOFT absent for %s", accession), call. = FALSE)
  if (file.exists(final_path)) {
    if (file.info(returned_path)$size != file.info(final_path)$size ||
        sha256_file(returned_path) != sha256_file(final_path)) {
      stop(sprintf("GEOquery quick SOFT differs from immutable file for %s", accession), call. = FALSE)
    }
  } else if (!file.rename(returned_path, final_path)) {
    stop(sprintf("Unable to finalize GEOquery quick SOFT for %s", accession), call. = FALSE)
  }
  invisible(final_path)
}

acquire_geoquery_core <- function(cohort, raw_root) {
  accession <- cohort$accession
  cohort_dir <- file.path(raw_root, accession)
  dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)

  acquire_quick_soft(cohort, cohort_dir)

  geo_objects <- GEOquery::getGEO(
    accession,
    GSEMatrix = TRUE,
    getGPL = FALSE,
    destdir = cohort_dir
  )
  rm(geo_objects)
  invisible(gc(verbose = FALSE))

  expected_supp <- vapply(
    Filter(function(x) geo_acquisition_method_for_role(x$role) == "GEOquery::getGEOSuppFiles", cohort$downloads),
    `[[`,
    character(1),
    "filename"
  )
  returned_supp <- character()
  if (length(expected_supp) > 0L) {
    escaped_supp <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", expected_supp)
    supp_filter <- paste0("^(", paste(escaped_supp, collapse = "|"), ")$")
    supp_info <- GEOquery::getGEOSuppFiles(
      accession,
      makeDirectory = TRUE,
      baseDir = raw_root,
      fetch_files = TRUE,
      filter_regex = supp_filter
    )
    returned_supp <- if (is.null(supp_info)) character() else supp_info$fname
  }
  if (!setequal(returned_supp, expected_supp)) {
    stop(sprintf("GEOquery supplementary-file set mismatch for %s", accession), call. = FALSE)
  }

  filelist_items <- Filter(
    function(x) identical(geo_acquisition_method_for_role(x$role), "curl_auxiliary_filelist"),
    cohort$downloads
  )
  for (item in filelist_items) {
    download_auxiliary_filelist(item$url, file.path(cohort_dir, item$filename))
  }
  invisible(TRUE)
}

strip_quotes <- function(x) {
  x <- sub('^"', "", x)
  sub('"$', "", x)
}

read_series_matrix_metadata <- function(path) {
  con <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  records <- list()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L || identical(line, "!series_matrix_table_begin")) break
    if (!startsWith(line, "!Sample_")) next
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(fields) < 2L) next
    key <- sub("^!Sample_", "", fields[[1L]])
    records[[length(records) + 1L]] <- list(key = key, values = strip_quotes(fields[-1L]))
  }
  if (length(records) == 0L) stop(sprintf("No sample metadata found in %s", path), call. = FALSE)
  keys <- vapply(records, `[[`, character(1), "key")
  keys <- make.unique(keys, sep = "__repeat_")
  n_samples <- length(records[[which(sub("__repeat_.*$", "", keys) == "geo_accession")[[1L]]]]$values)
  if (any(vapply(records, function(x) length(x$values), integer(1)) != n_samples)) {
    stop(sprintf("Inconsistent sample metadata width in %s", path), call. = FALSE)
  }
  out <- as.data.frame(
    setNames(lapply(records, `[[`, "values"), keys),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out
}

parse_quick_soft <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  accession <- sub("^!Series_geo_accession = ", "", grep("^!Series_geo_accession = ", lines, value = TRUE))
  platforms <- sub("^!Series_platform_id = ", "", grep("^!Series_platform_id = ", lines, value = TRUE))
  sample_ids <- sub("^!Series_sample_id = ", "", grep("^!Series_sample_id = ", lines, value = TRUE))
  if (length(accession) != 1L || length(platforms) < 1L || length(sample_ids) < 1L) {
    stop(sprintf("Incomplete official quick SOFT: %s", path), call. = FALSE)
  }
  list(accession = accession, platforms = platforms, sample_ids = sample_ids)
}

validate_raw_archive <- function(cohort_dir) {
  filelist_path <- file.path(cohort_dir, "filelist.txt")
  if (!file.exists(filelist_path)) return("no_official_filelist")
  filelist <- read.delim(
    filelist_path,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  archive_rows <- filelist[filelist[[1L]] == "Archive", , drop = FALSE]
  if (nrow(archive_rows) != 1L) {
    stop(sprintf("Expected one Archive row in %s", filelist_path), call. = FALSE)
  }
  archive_path <- file.path(cohort_dir, archive_rows$Name[[1L]])
  if (!file.exists(archive_path)) return("optional_source_absent")
  if (as.numeric(file.info(archive_path)$size) != as.numeric(archive_rows$Size[[1L]])) {
    stop(sprintf("Archive size differs from official filelist: %s", archive_path), call. = FALSE)
  }
  tar_status <- system2("tar", c("-tf", shQuote(archive_path)), stdout = FALSE, stderr = FALSE)
  if (tar_status != 0L) stop(sprintf("Tar integrity check failed: %s", archive_path), call. = FALSE)
  "present_verified"
}

manifest_rows <- list()
metadata_rows <- list()
for (cohort in config$cohorts) {
  accession <- cohort$accession
  cohort_dir <- file.path(raw_root, accession)
  dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)
  availability <- lightweight_availability(cohort, cohort_dir)
  if (!availability$analysis_matrix_present) {
    if (!network_download_authorized) {
      stop(
        sprintf(
          paste0(
            "No directly analyzable processed matrix is present for %s. ",
            "Automatic RAW/CEL/FASTQ/SRA acquisition is forbidden; obtain explicit authorization for processed-data download."
          ),
          accession
        ),
        call. = FALSE
      )
    }
  }
  if (length(availability$missing_required_non_source) > 0L) {
    if (!network_download_authorized) {
      stop(
        sprintf(
          "Missing processed/metadata files for %s in reuse-only mode; no network download was attempted",
          accession
        ),
        call. = FALSE
      )
    }
    acquire_geoquery_core(cohort, raw_root)
  }
  archive_validation_status <- validate_raw_archive(cohort_dir)

  for (item in cohort$downloads) {
    final_path <- file.path(cohort_dir, item$filename)
    if (identical(geo_acquisition_method_for_role(item$role), "retain_existing_source_no_auto_download") &&
        !file.exists(final_path)) {
      manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
        accession = accession,
        analysis_role = cohort$analysis_role,
        platform_expected = paste(cohort$expected_platforms, collapse = ";"),
        file_role = item$role,
        relative_path = file.path(raw_root, accession, item$filename),
        official_url = item$url,
        downloaded_at = "",
        bytes = 0,
        sha256 = "",
        action = archive_validation_status,
        initial_acquisition_method = "not_acquired_processed_only_policy",
        current_script_acquisition_path = geo_acquisition_method_for_role(item$role),
        stringsAsFactors = FALSE
      )
      write_manifest_atomic(manifest_rows, checksum_partial)
      next
    }
    if (!file.exists(final_path) || file.info(final_path)$size <= 0) {
      stop(sprintf("Expected acquisition output is absent or empty: %s", final_path), call. = FALSE)
    }
    if (!is.null(previous_manifest)) validate_existing(final_path, item$url)
    if (grepl("\\.gz$", final_path, ignore.case = TRUE) &&
        system2("gzip", c("-t", shQuote(final_path))) != 0L) {
      stop(sprintf("Gzip integrity check failed: %s", final_path), call. = FALSE)
    }
    Sys.chmod(final_path, mode = "0444", use_umask = TRUE)
    relative_path <- file.path(raw_root, accession, item$filename)
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      accession = accession,
      analysis_role = cohort$analysis_role,
      platform_expected = paste(cohort$expected_platforms, collapse = ";"),
      file_role = item$role,
      relative_path = relative_path,
      official_url = item$url,
      downloaded_at = format(file.info(final_path)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
      bytes = as.numeric(file.info(final_path)$size),
      sha256 = sha256_file(final_path),
      action = "present_verified",
      initial_acquisition_method = "curl_first_execution_documented_deviation",
      current_script_acquisition_path = geo_acquisition_method_for_role(item$role),
      stringsAsFactors = FALSE
    )
    write_manifest_atomic(manifest_rows, checksum_partial)
  }
  quick_path <- file.path(cohort_dir, paste0(accession, "_quick.soft.txt"))
  official <- parse_quick_soft(quick_path)
  if (official$accession != accession) stop(sprintf("Accession mismatch for %s", accession), call. = FALSE)
  if (!setequal(official$platforms, cohort$expected_platforms)) {
    stop(sprintf("Official platform mismatch for %s", accession), call. = FALSE)
  }
  if (length(official$sample_ids) != cohort$expected_sample_count) {
    stop(sprintf("Official sample count mismatch for %s", accession), call. = FALSE)
  }

  matrices <- vapply(
    Filter(function(x) identical(x$role, "series_matrix"), cohort$downloads),
    function(x) file.path(cohort_dir, x$filename),
    character(1)
  )
  cohort_metadata <- do.call(rbind, lapply(matrices, read_series_matrix_metadata))
  if (!"geo_accession" %in% names(cohort_metadata) || !"platform_id" %in% names(cohort_metadata)) {
    stop(sprintf("Required sample fields absent for %s", accession), call. = FALSE)
  }
  if (anyDuplicated(cohort_metadata$geo_accession)) {
    stop(sprintf("Duplicate GSM across series matrices for %s", accession), call. = FALSE)
  }
  if (!setequal(cohort_metadata$geo_accession, official$sample_ids)) {
    stop(sprintf("Series matrix GSM set differs from official quick SOFT for %s", accession), call. = FALSE)
  }
  if (!setequal(unique(cohort_metadata$platform_id), cohort$expected_platforms)) {
    stop(sprintf("Sample platform mismatch for %s", accession), call. = FALSE)
  }
  cohort_metadata$series_accession <- accession
  cohort_metadata$analysis_role <- cohort$analysis_role
  cohort_metadata$outcome_previously_observed <- cohort$outcome_previously_observed
  cohort_metadata$task3_inclusion_status <- "all_official_samples_unfiltered"
  cohort_metadata$task3_note <- cohort$notes
  metadata_rows[[length(metadata_rows) + 1L]] <- cohort_metadata
}

all_names <- unique(unlist(lapply(metadata_rows, names)))
align_columns <- function(x) {
  missing <- setdiff(all_names, names(x))
  for (name in missing) x[[name]] <- NA_character_
  x[all_names]
}
all_metadata <- do.call(rbind, lapply(metadata_rows, align_columns))
if (anyDuplicated(all_metadata$geo_accession)) stop("GSM identifiers are not globally unique", call. = FALSE)

write.table(
  all_metadata,
  metadata_output,
  sep = "\t",
  quote = TRUE,
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
manifest <- write_manifest_atomic(manifest_rows, checksum_output)
if (file.exists(checksum_partial) && unlink(checksum_partial) != 0L) {
  stop(sprintf("Unable to remove completed partial manifest: %s", checksum_partial), call. = FALSE)
}

cat(sprintf("Frozen %d official samples across %d GEO series\n", nrow(all_metadata), length(config$cohorts)))
cat(sprintf("Recorded %d immutable downloads\n", nrow(manifest)))
