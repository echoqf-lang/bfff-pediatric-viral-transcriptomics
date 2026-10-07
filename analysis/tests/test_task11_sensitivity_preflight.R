#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

read_series_structure <- function(path) {
  con <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  metadata <- list()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("Missing series_matrix_table_begin: ", path, call. = FALSE)
    if (identical(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_geo_accession\t")) metadata$sample_id <- line
    if (startsWith(line, "!Sample_title\t")) metadata$title <- line
    if (startsWith(line, "!Sample_description\t")) metadata$description <- line
    if (startsWith(line, "!Sample_data_processing\t")) metadata$processing <- line
    if (startsWith(line, "!Sample_characteristics_ch1\t")) {
      metadata$characteristics <- c(metadata$characteristics, line)
    }
  }
  parse_values <- function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
    gsub('^"|"$', "", fields)
  }
  header <- readLines(con, n = 1L, warn = FALSE)
  columns <- parse_values(header)
  ids <- list()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line) || identical(line, "!series_matrix_table_end")) break
    ids[[length(ids) + 1L]] <- gsub('^"|"$', "", sub("\t.*$", "", line))
  }
  list(
    sample_id = parse_values(metadata$sample_id),
    title = parse_values(metadata$title),
    description = if (is.null(metadata$description)) character() else parse_values(metadata$description),
    processing = unique(parse_values(metadata$processing)),
    characteristic_labels = unique(vapply(metadata$characteristics, function(x) {
      sub(":.*$", "", parse_values(x)[[1L]])
    }, character(1))),
    table_columns = columns,
    feature_id = unlist(ids, use.names = FALSE)
  )
}

manifest <- read.delim(
  file.path("data", "clean", "sample_manifest_frozen.tsv"),
  sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
targets <- read.delim(
  file.path("data", "clean", "bfff_targets_sensitivity.tsv"),
  sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE
)

stopifnot(
  nrow(targets) == 1270L,
  !anyDuplicated(targets$entrez_id),
  all(grepl("^[0-9]+$", targets$entrez_id)),
  sum(nzchar(targets$ensembl_gene_ids)) == 1267L,
  length(unique(unlist(strsplit(
    targets$ensembl_gene_ids[nzchar(targets$ensembl_gene_ids)], ";", fixed = TRUE
  )))) == 1387L
)

for (cohort in c("GSE105450", "GSE103842")) {
  derived <- readRDS(file.path("data", "derived", paste0(cohort, "_expression.rds")))
  expected <- if (cohort == "GSE105450") 885L else 816L
  stopifnot(sum(targets$entrez_id %in% rownames(derived$expression)) == expected)
}

gse105450 <- manifest[manifest$series_accession == "GSE105450" & manifest$include, , drop = FALSE]
hospital_only <- gse105450[gse105450$severity %in% c("hospitalized", "healthy"), , drop = FALSE]
hospital_only$case_status <- factor(
  ifelse(hospital_only$status == "RSV_case", "RSV", "healthy"),
  levels = c("healthy", "RSV")
)
hospital_design <- model.matrix(~ case_status + age_months + sex + batch, data = hospital_only)
gse105450_derived <- readRDS(file.path("data", "derived", "GSE105450_expression.rds"))
stopifnot(
  nrow(hospital_only) == 89L,
  sum(hospital_only$case_status == "RSV") == 56L,
  sum(hospital_only$case_status == "healthy") == 33L,
  all(hospital_only$sample_id %in% colnames(gse105450_derived$expression)),
  qr(hospital_design)$rank == ncol(hospital_design)
)

gse188427 <- read_series_structure(file.path("data", "raw", "GSE188427", "GSE188427_series_matrix.txt.gz"))
gse188427_manifest <- manifest[manifest$series_accession == "GSE188427", , drop = FALSE]
gse188427_included <- gse188427_manifest[gse188427_manifest$include, , drop = FALSE]
gse188427_entrez <- sub("_at$", "", gse188427$feature_id)
stopifnot(
  identical(gse188427$sample_id, gse188427$table_columns),
  setequal(gse188427$sample_id, gse188427_manifest$sample_id),
  nrow(gse188427_included) == 198L,
  sum(gse188427_included$status == "RSV_case") == 147L,
  sum(gse188427_included$status == "healthy_control") == 51L,
  all(is.na(gse188427_included$age_months)),
  all(is.na(gse188427_included$sex)),
  length(unique(gse188427_included$batch)) == nrow(gse188427_included),
  length(gse188427_entrez) == 18604L,
  !anyDuplicated(gse188427_entrez),
  all(grepl("^[0-9]+$", gse188427_entrez)),
  sum(targets$entrez_id %in% gse188427_entrez) == 1201L,
  any(grepl("Robust Multi-array Average", gse188427$processing, fixed = TRUE)),
  any(grepl("alternative CDF mappings (ENTREZG, version 20)", gse188427$processing, fixed = TRUE))
)

gse155925 <- read_series_structure(file.path("data", "raw", "GSE155925", "GSE155925_series_matrix.txt.gz"))
count_path <- file.path("data", "raw", "GSE155925", "GSE155925_Raw_counts_matrix.txt.gz")
count_header <- names(read.delim(count_path, sep = "\t", nrows = 0L, check.names = FALSE, quote = ""))[-1L]
count_lines <- readLines(gzfile(count_path), warn = FALSE)
count_ensembl <- sub("^.*:", "", sub("\t.*$", "", count_lines[-1L]))
target_ensembl <- strsplit(targets$ensembl_gene_ids, ";", fixed = TRUE)
gse155925_manifest <- manifest[manifest$series_accession == "GSE155925" & manifest$include, , drop = FALSE]
gse155925_manifest$virus_group <- factor(
  ifelse(gse155925_manifest$status == "RSV_case", "RSV", "other"),
  levels = c("other", "RSV")
)
batch_parts <- do.call(rbind, strsplit(gse155925_manifest$batch, ";", fixed = TRUE))
gse155925_manifest$hospital_batch <- batch_parts[, 1L]
gse155925_manifest$enrollment_batch <- batch_parts[, 2L]
gse155925_design <- model.matrix(
  ~ virus_group + age_months + sex + hospital_batch + enrollment_batch,
  data = gse155925_manifest
)
stopifnot(
  identical(count_header, gse155925$title),
  nrow(gse155925_manifest) == 48L,
  sum(gse155925_manifest$virus_group == "RSV") == 31L,
  sum(gse155925_manifest$virus_group == "other") == 17L,
  all(gse155925_manifest$coinfection == "none"),
  all(!is.na(gse155925_manifest$age_months)),
  all(!is.na(gse155925_manifest$sex)),
  qr(gse155925_design)$rank == ncol(gse155925_design),
  length(count_ensembl) == 19919L,
  !anyDuplicated(count_ensembl),
  sum(vapply(target_ensembl, function(x) any(x %in% count_ensembl), logical(1))) == 1254L
)

gse103119 <- read_series_structure(file.path("data", "raw", "GSE103119", "GSE103119_series_matrix.txt.gz"))
gse103119_manifest <- manifest[manifest$series_accession == "GSE103119" & manifest$include, , drop = FALSE]
gse103119_manifest$case_status <- factor(
  ifelse(gse103119_manifest$status == "RSV_case", "RSV", "healthy"),
  levels = c("healthy", "RSV")
)
gse103119_design <- model.matrix(~ case_status + age_months + sex, data = gse103119_manifest)
detection_header <- strsplit(
  readLines(gzfile(file.path("data", "raw", "GSE103119", "GSE103119_non-normalized.txt.gz")), n = 1L),
  "\t", fixed = TRUE
)[[1L]]
detection_array_ids <- detection_header[seq.int(2L, length(detection_header), by = 2L)]
included_array_ids <- sub("^.* ", "", gse103119_manifest$title)
stopifnot(
  identical(gse103119$sample_id, gse103119$table_columns),
  setequal(gse103119$sample_id, manifest$sample_id[manifest$series_accession == "GSE103119"]),
  nrow(gse103119_manifest) == 31L,
  sum(gse103119_manifest$case_status == "RSV") == 11L,
  sum(gse103119_manifest$case_status == "healthy") == 20L,
  all(gse103119_manifest$pathogen[gse103119_manifest$case_status == "RSV"] == "RSV"),
  all(gse103119_manifest$coinfection == "none"),
  all(included_array_ids %in% detection_array_ids),
  file.exists(file.path("data", "raw", "GPL10558", "GPL10558.annot.gz")),
  qr(gse103119_design)$rank == ncol(gse103119_design)
)

all_characteristic_labels <- unique(c(
  read_series_structure(file.path("data", "raw", "GSE105450", "GSE105450_series_matrix.txt.gz"))$characteristic_labels,
  read_series_structure(file.path("data", "raw", "GSE103842", "GSE103842_series_matrix.txt.gz"))$characteristic_labels,
  gse188427$characteristic_labels,
  gse155925$characteristic_labels,
  gse103119$characteristic_labels
))
stopifnot(!any(grepl(
  "CBC|complete blood|white blood|WBC|neutroph|lymph|monocyte|cell proportion|cell count|leukocyte",
  all_characteristic_labels, ignore.case = TRUE
)))

cat("Task 11 zero-download sensitivity preflight checks passed\n")
