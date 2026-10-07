#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
if (!requireNamespace("GEOquery", quietly = TRUE)) stop("GEOquery is required", call. = FALSE)

accessions <- c("GSE97742", "GSE41374")
raw_root <- file.path("data", "raw", "geo_processed_only")
derived_root <- file.path("data", "derived", "barrier_inflammation_repair")
dir.create(raw_root, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_root, recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  sub("[[:space:]].*$", "", out[[1L]])
}

manifest <- list()
for (accession in accessions) {
  cohort_dir <- file.path(raw_root, accession)
  dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)
  gse <- GEOquery::getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE, destdir = cohort_dir)
  if (length(gse) != 1L) stop(sprintf("Expected one platform for %s", accession), call. = FALSE)
  eset <- gse[[1L]]
  saveRDS(eset, file.path(derived_root, paste0(accession, "_processed_eset.rds")), compress = "xz")
  matrix_files <- list.files(cohort_dir, pattern = "series_matrix.*\\.txt\\.gz$", full.names = TRUE)
  if (length(matrix_files) != 1L) stop(sprintf("Expected one Series Matrix for %s", accession), call. = FALSE)
  manifest[[accession]] <- data.frame(
    accession = accession,
    relative_path = matrix_files,
    bytes = file.info(matrix_files)$size,
    sha256 = sha256_file(matrix_files),
    acquisition = "GEOquery Series Matrix only; no raw IDAT/CEL",
    stringsAsFactors = FALSE
  )
  message(accession, " samples=", ncol(Biobase::exprs(eset)), " features=", nrow(Biobase::exprs(eset)))
}

write.table(
  do.call(rbind, manifest),
  file.path(derived_root, "processed_input_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8"
)
