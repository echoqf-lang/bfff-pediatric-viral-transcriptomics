#!/usr/bin/env Rscript

# Read-only verifier for external inputs needed by the processed-only workflow.
# This script never opens a network connection and never creates or changes files.

args <- commandArgs(trailingOnly = TRUE)
list_only <- identical(args, "--list")
if (length(args) > 1L || (length(args) == 1L && !list_only)) {
  stop("Usage: Rscript --vanilla analysis/R/verify_external_inputs.R [--list]", call. = FALSE)
}

project_root <- normalizePath(
  Sys.getenv("VERIFY_PROJECT_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)
external_root <- normalizePath(
  Sys.getenv("VERIFY_EXTERNAL_INPUT_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

allowed_geo_roles <- c(
  "official_quick_metadata",
  "series_matrix",
  "non_normalized_expression",
  "raw_counts",
  "official_supplement_filelist"
)
forbidden_source_types <- c("RAW.tar", "CEL", "CEL.gz", "FASTQ", "SRA")

read_tsv <- function(relative_path, required_columns) {
  path <- file.path(project_root, relative_path)
  if (!file.exists(path)) {
    stop(sprintf("Required manifest is missing: %s", relative_path), call. = FALSE)
  }
  value <- read.delim(
    path,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
  missing_columns <- setdiff(required_columns, names(value))
  if (length(missing_columns) > 0L) {
    stop(sprintf(
      "Manifest %s is missing columns: %s",
      relative_path,
      paste(missing_columns, collapse = ", ")
    ), call. = FALSE)
  }
  value
}

target_manifest <- read_tsv(
  "logs/checksums/task2_sha256.tsv",
  c("role", "path", "sha256")
)
target_manifest <- target_manifest[target_manifest$role == "input", , drop = FALSE]
if (nrow(target_manifest) == 0L) {
  stop("No role=input entries found in logs/checksums/task2_sha256.tsv", call. = FALSE)
}

geo_manifest <- read_tsv(
  "logs/checksums/geo_sha256.tsv",
  c("file_role", "relative_path", "official_url", "bytes", "sha256")
)
geo_manifest <- geo_manifest[
  geo_manifest$file_role %in% allowed_geo_roles,
  ,
  drop = FALSE
]
if (nrow(geo_manifest) == 0L) {
  stop("No permitted processed GEO entries found in logs/checksums/geo_sha256.tsv", call. = FALSE)
}

gpl_manifest <- read_tsv(
  "logs/checksums/gpl10558_annotation_sha256.tsv",
  c("platform", "official_url", "relative_path", "bytes", "sha256")
)
if (nrow(gpl_manifest) == 0L) {
  stop("GPL10558 annotation checksum manifest is empty", call. = FALSE)
}

config_path <- file.path(project_root, "analysis/config/cohorts.yml")
if (!file.exists(config_path)) {
  stop("Required config is missing: analysis/config/cohorts.yml", call. = FALSE)
}
config_text <- paste(readLines(config_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
missing_roles <- allowed_geo_roles[!vapply(
  allowed_geo_roles,
  grepl,
  logical(1),
  x = config_text,
  fixed = TRUE
)]
missing_forbidden <- forbidden_source_types[!vapply(
  forbidden_source_types,
  grepl,
  logical(1),
  x = config_text,
  fixed = TRUE
)]
missing_urls <- unique(geo_manifest$official_url[!vapply(
  geo_manifest$official_url,
  grepl,
  logical(1),
  x = config_text,
  fixed = TRUE
)])
if (length(missing_roles) > 0L || length(missing_forbidden) > 0L || length(missing_urls) > 0L) {
  stop(paste(
    "cohorts.yml does not agree with the processed-only checksum inventory.",
    if (length(missing_roles) > 0L) paste("Missing roles:", paste(missing_roles, collapse = ", ")) else "",
    if (length(missing_forbidden) > 0L) paste("Missing forbidden types:", paste(missing_forbidden, collapse = ", ")) else "",
    if (length(missing_urls) > 0L) paste("Missing URLs:", paste(missing_urls, collapse = ", ")) else ""
  ), call. = FALSE)
}

inventory <- rbind(
  data.frame(
    category = "user_provided_target",
    role = target_manifest$role,
    relative_path = target_manifest$path,
    official_url = "USER_PROVIDED_NO_PUBLIC_DOWNLOAD_CLAIMED",
    expected_bytes = NA_character_,
    expected_sha256 = target_manifest$sha256,
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "geo_processed",
    role = geo_manifest$file_role,
    relative_path = geo_manifest$relative_path,
    official_url = geo_manifest$official_url,
    expected_bytes = as.character(geo_manifest$bytes),
    expected_sha256 = geo_manifest$sha256,
    stringsAsFactors = FALSE
  ),
  data.frame(
    category = "geo_platform_annotation",
    role = gpl_manifest$platform,
    relative_path = gpl_manifest$relative_path,
    official_url = gpl_manifest$official_url,
    expected_bytes = as.character(gpl_manifest$bytes),
    expected_sha256 = gpl_manifest$sha256,
    stringsAsFactors = FALSE
  )
)

sha256_file <- function(path) {
  shasum <- Sys.which("shasum")
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(shasum)) {
    output <- system2(shasum, c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  } else if (nzchar(sha256sum)) {
    output <- system2(sha256sum, shQuote(path), stdout = TRUE, stderr = TRUE)
  } else {
    stop("Neither shasum nor sha256sum is available", call. = FALSE)
  }
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(sprintf("SHA-256 command failed for %s", path), call. = FALSE)
  }
  match <- regmatches(output[[1L]], regexpr("^[0-9a-fA-F]{64}", output[[1L]]))
  if (length(match) != 1L || !nzchar(match)) {
    stop(sprintf("Could not parse SHA-256 for %s", path), call. = FALSE)
  }
  tolower(match)
}

verify_one <- function(relative_path, expected_bytes, expected_sha256) {
  if (list_only) {
    return(c(status = "not_checked", actual_bytes = NA_character_, actual_sha256 = NA_character_))
  }
  path <- file.path(external_root, relative_path)
  if (!file.exists(path)) {
    return(c(status = "missing", actual_bytes = NA_character_, actual_sha256 = NA_character_))
  }
  actual_bytes <- as.character(file.info(path)$size)
  if (!is.na(expected_bytes) && nzchar(expected_bytes) && actual_bytes != expected_bytes) {
    return(c(status = "byte_mismatch", actual_bytes = actual_bytes, actual_sha256 = NA_character_))
  }
  actual_sha256 <- sha256_file(path)
  status <- if (identical(tolower(actual_sha256), tolower(expected_sha256))) {
    "verified"
  } else {
    "sha256_mismatch"
  }
  c(status = status, actual_bytes = actual_bytes, actual_sha256 = actual_sha256)
}

checks <- t(vapply(
  seq_len(nrow(inventory)),
  function(i) verify_one(
    inventory$relative_path[[i]],
    inventory$expected_bytes[[i]],
    inventory$expected_sha256[[i]]
  ),
  character(3)
))
result <- cbind(
  inventory[c("category", "role", "relative_path")],
  data.frame(
    status = checks[, "status"],
    expected_bytes = inventory$expected_bytes,
    actual_bytes = checks[, "actual_bytes"],
    expected_sha256 = inventory$expected_sha256,
    actual_sha256 = checks[, "actual_sha256"],
    official_url = inventory$official_url,
    stringsAsFactors = FALSE
  )
)

cat("# external_input_verifier_version=1\n")
cat("# read_only=TRUE\n")
cat("# network_attempted=FALSE\n")
cat("# mode=", if (list_only) "list" else "verify", "\n", sep = "")
cat("# project_root=", project_root, "\n", sep = "")
cat("# external_input_root=", external_root, "\n", sep = "")
write.table(result, stdout(), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

if (!list_only) {
  failed <- result$status != "verified"
  cat(sprintf(
    "# summary verified=%d failed=%d total=%d\n",
    sum(!failed),
    sum(failed),
    nrow(result)
  ))
  if (any(failed)) {
    quit(save = "no", status = 1L)
  }
}
