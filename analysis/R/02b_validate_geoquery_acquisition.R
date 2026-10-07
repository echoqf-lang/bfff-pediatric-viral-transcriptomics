#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("yaml", quietly = TRUE) || !requireNamespace("GEOquery", quietly = TRUE)) {
  stop("Frozen yaml and GEOquery packages are required", call. = FALSE)
}

config <- yaml::read_yaml("analysis/config/cohorts.yml")
if (!identical(
  Sys.getenv("GEOQUERY_ALLOW_FULL_REDOWNLOAD", unset = ""),
  "explicit_user_authorization"
)) {
  stop(
    paste0(
      "Full GEOquery re-download is disabled by default. ",
      "It may run only after new explicit user authorization by setting ",
      "GEOQUERY_ALLOW_FULL_REDOWNLOAD=explicit_user_authorization."
    ),
    call. = FALSE
  )
}
validation_output <- "logs/checksums/geoquery_redownload_validation.tsv"
validation_root <- Sys.getenv("GEOQUERY_VALIDATION_ROOT", unset = "")
if (!nzchar(validation_root)) validation_root <- tempfile("bfff_geoquery_validation_")
if (dir.exists(validation_root) && length(list.files(validation_root, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("GEOQUERY_VALIDATION_ROOT must be absent or empty", call. = FALSE)
}
dir.create(validation_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(validation_output), recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop(sprintf("SHA-256 failed: %s", path), call. = FALSE)
  sub("[[:space:]].*$", "", out[[1L]])
}

uncompressed_sha256 <- function(path) {
  command <- sprintf("gzip -cd %s | shasum -a 256", shQuote(path))
  out <- system2("sh", c("-c", shQuote(command)), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop(sprintf("Uncompressed SHA-256 failed: %s", path), call. = FALSE)
  sub("[[:space:]].*$", "", out[[1L]])
}

method_for_role <- function(role) {
  if (identical(role, "official_quick_metadata")) return("GEOquery::getGEOfile")
  if (identical(role, "series_matrix")) return("GEOquery::getGEO")
  if (identical(role, "official_supplement_filelist")) return("curl_auxiliary_filelist")
  "GEOquery::getGEOSuppFiles"
}

download_auxiliary_filelist <- function(url, path) {
  status <- system2(
    "curl",
    c("-L", "--fail", "--show-error", "--retry", "8", "--output", shQuote(path), shQuote(url))
  )
  if (status != 0L || !file.exists(path) || file.info(path)$size <= 0) {
    stop(sprintf("Auxiliary filelist download failed: %s", url), call. = FALSE)
  }
}

comparison_rows <- list()
for (cohort in config$cohorts) {
  accession <- cohort$accession
  cohort_stage <- file.path(validation_root, accession)
  dir.create(cohort_stage, recursive = TRUE, showWarnings = FALSE)

  geo_objects <- GEOquery::getGEO(
    accession,
    GSEMatrix = TRUE,
    getGPL = FALSE,
    destdir = cohort_stage
  )
  rm(geo_objects)
  invisible(gc(verbose = FALSE))

  GEOquery::getGEOSuppFiles(
    accession,
    makeDirectory = TRUE,
    baseDir = validation_root,
    fetch_files = TRUE
  )
  quick_path <- GEOquery::getGEOfile(accession, destdir = cohort_stage, amount = "quick")

  for (item in cohort$downloads) {
    method <- method_for_role(item$role)
    staged_path <- if (identical(item$role, "official_quick_metadata")) {
      quick_path
    } else {
      file.path(cohort_stage, item$filename)
    }
    if (identical(method, "curl_auxiliary_filelist")) {
      download_auxiliary_filelist(item$url, staged_path)
    }
    frozen_path <- file.path(config$raw_root, accession, item$filename)
    if (!file.exists(staged_path) || !file.exists(frozen_path)) {
      stop(sprintf("Comparison file absent for %s/%s", accession, item$filename), call. = FALSE)
    }
    staged_sha <- sha256_file(staged_path)
    frozen_sha <- sha256_file(frozen_path)
    comparison <- if (identical(staged_sha, frozen_sha) &&
                      file.info(staged_path)$size == file.info(frozen_path)$size) {
      "exact_bytes"
    } else if (grepl("\\.gz$", item$filename, ignore.case = TRUE) &&
               identical(uncompressed_sha256(staged_path), uncompressed_sha256(frozen_path))) {
      "same_uncompressed_content"
    } else {
      "mismatch"
    }
    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
      accession = accession,
      file_role = item$role,
      filename = item$filename,
      acquisition_function = method,
      staged_bytes = as.numeric(file.info(staged_path)$size),
      frozen_bytes = as.numeric(file.info(frozen_path)$size),
      staged_sha256 = staged_sha,
      frozen_sha256 = frozen_sha,
      comparison = comparison,
      validated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      stringsAsFactors = FALSE
    )
  }
}

validation <- do.call(rbind, comparison_rows)
temporary_output <- paste0(validation_output, ".tmp")
write.table(
  validation,
  temporary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
if (!file.rename(temporary_output, validation_output)) stop("Validation log finalize failed", call. = FALSE)
if (any(validation$comparison == "mismatch")) stop("GEOquery redownload differs from frozen files", call. = FALSE)

keep_validation <- identical(tolower(Sys.getenv("GEOQUERY_KEEP_VALIDATION", unset = "false")), "true")
if (!keep_validation && unlink(validation_root, recursive = TRUE, force = TRUE) != 0L) {
  stop(sprintf("Unable to remove temporary validation root: %s", validation_root), call. = FALSE)
}
cat(sprintf("Validated %d files via GEOquery core paths; exact=%d, uncompressed_equal=%d\n",
            nrow(validation), sum(validation$comparison == "exact_bytes"),
            sum(validation$comparison == "same_uncompressed_content")))
