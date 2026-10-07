#!/usr/bin/env Rscript

expected_seed <- 20260804L
project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
expected_raw_root <- file.path(project_root, "data", "raw")
dir.create(expected_raw_root, recursive = TRUE, showWarnings = FALSE)
expected_raw_root <- normalizePath(
  expected_raw_root,
  winslash = "/",
  mustWork = TRUE
)
configured_raw_root <- Sys.getenv(
  "ANALYSIS_RAW_ROOT",
  unset = "data/raw"
)
if (!grepl("^/", configured_raw_root)) {
  configured_raw_root <- file.path(project_root, configured_raw_root)
}
configured_raw_root <- normalizePath(
  configured_raw_root,
  winslash = "/",
  mustWork = FALSE
)

if (!identical(configured_raw_root, expected_raw_root)) {
  stop(
    sprintf(
      "ANALYSIS_RAW_ROOT must resolve to %s; received %s",
      expected_raw_root,
      configured_raw_root
    ),
    call. = FALSE
  )
}

tracked_files <- system2("git", c("ls-files"), stdout = TRUE, stderr = TRUE)
git_status <- attr(tracked_files, "status")
if (!is.null(git_status) && git_status != 0L) {
  stop("Unable to audit Git-tracked files", call. = FALSE)
}

raw_file_patterns <- c(
  "^data/raw(/|$)",
  "_raw\\.rds$",
  "series_matrix.*\\.txt($|\\.)",
  "non-normalized",
  "\\.CEL\\.gz$",
  "expression.*\\.(rds|RData|txt)(\\.gz)?$"
)
tracked_raw_files <- tracked_files[
  Reduce(
    `|`,
    lapply(raw_file_patterns, grepl, x = tracked_files, ignore.case = TRUE)
  )
]
if (length(tracked_raw_files) > 0L) {
  stop(
    sprintf(
      "Raw expression/download files must not be tracked by Git: %s",
      paste(tracked_raw_files, collapse = ", ")
    ),
    call. = FALSE
  )
}

configured_seed <- Sys.getenv("ANALYSIS_SEED", unset = as.character(expected_seed))
configured_seed <- suppressWarnings(as.integer(configured_seed))

if (is.na(configured_seed) || configured_seed != expected_seed) {
  stop(
    sprintf(
      "ANALYSIS_SEED must be exactly %d; received %s",
      expected_seed,
      Sys.getenv("ANALYSIS_SEED", unset = "<unset>")
    ),
    call. = FALSE
  )
}

set.seed(expected_seed)
if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
  stop("Global random seed was not initialized", call. = FALSE)
}

required_packages <- c(
  "GEOquery", "limma", "Biobase", "AnnotationDbi", "fgsea",
  "metafor", "msigdbr", "org.Hs.eg.db", "sva", "variancePartition", "yaml", "testthat"
)

if (!requireNamespace("renv", quietly = TRUE)) {
  stop("The project renv package is not available", call. = FALSE)
}

project_library <- normalizePath(
  renv::paths$library(project = project_root),
  winslash = "/",
  mustWork = TRUE
)
active_libraries <- normalizePath(
  .libPaths(),
  winslash = "/",
  mustWork = TRUE
)
if (!project_library %in% active_libraries) {
  stop(
    sprintf("Project renv library is not active: %s", project_library),
    call. = FALSE
  )
}

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    sprintf("Missing required packages: %s", paste(missing_packages, collapse = ", ")),
    call. = FALSE
  )
}

lockfile <- renv::lockfile_read(file.path(project_root, "renv.lock"))
locked_versions <- vapply(
  required_packages,
  function(package) {
    record <- lockfile$Packages[[package]]
    if (is.null(record) || is.null(record$Version)) NA_character_ else record$Version
  },
  character(1)
)
installed_versions <- vapply(
  required_packages,
  function(package) as.character(packageVersion(package)),
  character(1)
)
version_mismatches <- required_packages[
  is.na(locked_versions) |
    as.character(package_version(installed_versions)) !=
      as.character(package_version(locked_versions))
]
if (length(version_mismatches) > 0L) {
  mismatch_details <- sprintf(
    "%s(installed=%s, locked=%s)",
    version_mismatches,
    installed_versions[version_mismatches],
    locked_versions[version_mismatches]
  )
  stop(
    sprintf("Package versions differ from renv.lock: %s", paste(mismatch_details, collapse = ", ")),
    call. = FALSE
  )
}

invisible(lapply(required_packages, function(package) {
  suppressPackageStartupMessages(
    library(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  )
}))

dir.create("logs/session_info", recursive = TRUE, showWarnings = FALSE)

bioconductor_version <- if (requireNamespace("BiocManager", quietly = TRUE)) {
  as.character(BiocManager::version())
} else {
  NA_character_
}

package_versions <- data.frame(
  package = c(required_packages, "BiocManager", "renv", "PLATFORM_ANNOTATION_PACKAGE"),
  version = c(
    vapply(required_packages, function(package) as.character(packageVersion(package)), character(1)),
    if (requireNamespace("BiocManager", quietly = TRUE)) as.character(packageVersion("BiocManager")) else NA_character_,
    if (requireNamespace("renv", quietly = TRUE)) as.character(packageVersion("renv")) else NA_character_,
    NA_character_
  ),
  status = c(
    rep("installed_and_loadable", length(required_packages)),
    if (requireNamespace("BiocManager", quietly = TRUE)) "installed" else "not_installed",
    if (requireNamespace("renv", quietly = TRUE)) "installed" else "not_installed",
    "pending_platform_metadata"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  package_versions,
  "logs/session_info/package_versions.csv",
  row.names = FALSE,
  na = ""
)

session_lines <- c(
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  sprintf("Analysis seed: %d", expected_seed),
  sprintf("R version: %s", R.version.string),
  sprintf("R platform: %s", R.version$platform),
  sprintf("Bioconductor version: %s", bioconductor_version),
  sprintf("Project renv library: %s", project_library),
  sprintf("Confirmatory GEO download root: %s", expected_raw_root),
  "Git raw-file policy audit: passed",
  "Required package versions match renv.lock: passed",
  "Platform annotation package: pending until platform metadata are determined (Task 3)",
  "Data flow: data/raw -> data/metadata and data/clean -> data/derived -> results/cohort -> results/meta -> results/tables and results/figures",
  "Raw-data policy: data/raw is immutable; transformations must write new artifacts outside data/raw",
  "",
  capture.output(sessionInfo())
)
strip_trailing_whitespace <- function(x) sub("[[:space:]]+$", "", x)
stopifnot(
  identical(strip_trailing_whitespace("Matrix products: default"), "Matrix products: default"),
  identical(strip_trailing_whitespace(paste0("value", " ", "\t")), "value")
)
session_lines <- strip_trailing_whitespace(session_lines)

writeLines(session_lines, "logs/session_info/session_info.txt", useBytes = TRUE)
cat("Environment check passed with seed", expected_seed, "\n")
