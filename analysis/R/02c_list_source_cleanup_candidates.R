#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("yaml", quietly = TRUE)) stop("Package yaml is required", call. = FALSE)
config <- yaml::read_yaml("analysis/config/cohorts.yml")
manifest <- read.delim(config$checksum_output, sep = "\t", quote = "", check.names = FALSE)
output <- "logs/checksums/source_file_cleanup_candidates.tsv"

source_pattern <- "(\\.CEL(\\.gz)?$|\\.fastq(\\.gz)?$|\\.fq(\\.gz)?$|\\.sra$|_RAW\\.tar$)"
candidates <- manifest[
  (manifest$file_role == "raw_archive" |
    grepl(source_pattern, manifest$relative_path, ignore.case = TRUE)) &
    manifest$action != "optional_source_absent" &
    as.numeric(manifest$bytes) > 0,
  c("accession", "file_role", "relative_path", "bytes", "sha256"),
  drop = FALSE
]
candidates$source_type <- ifelse(
  grepl("_RAW\\.tar$", candidates$relative_path, ignore.case = TRUE),
  "RAW.tar",
  "individual_source_file"
)
candidates$cleanup_status <- "candidate_only_not_deleted"
candidates <- candidates[
  , c("accession", "source_type", "file_role", "relative_path", "bytes", "sha256", "cleanup_status")
]

if (any(candidates$file_role %in% c("non_normalized_expression", "raw_counts"))) {
  stop("Processed/non-normalized analysis inputs must not be cleanup candidates", call. = FALSE)
}
if (any(!file.exists(candidates$relative_path))) stop("Cleanup candidate file is absent", call. = FALSE)

write.table(
  candidates,
  output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
cat(sprintf(
  "Listed %d source-file cleanup candidates; releasable bytes=%d; deleted=0\n",
  nrow(candidates),
  sum(candidates$bytes)
))
