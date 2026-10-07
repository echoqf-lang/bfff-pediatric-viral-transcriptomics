#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

amendment_path <- file.path(
  "docs", "science-superpowers", "preregistrations", "amendments",
  "2026-08-05-task11-sensitivity-operationalization.md"
)
amendment <- paste(readLines(amendment_path, warn = FALSE), collapse = "\n")
source(file.path("analysis", "R", "task11_design_contract.R"))
required_contract_text <- c(
  "Status:** FROZEN",
  "573c5f577dc071af45289ff68e051d762f7b4403288dc78267e62dab8755abcd",
  "~ case_status",
  "~ case_status + age_months + sex + technical_batch",
  "~ case_status + age_months + sex",
  "~ virus_group + age_months + sex + hospital_batch + enrollment_batch",
  "RSV减健康",
  "单RSV减其他单呼吸道病毒",
  "少于该完整靶点集50%",
  "GSE188427与GSE105450绝不进入同一个",
  "Holm校正",
  "不能救回已经失败的H1联合判定",
  "不实施CBC或细胞比例调整"
)
stopifnot(all(vapply(required_contract_text, grepl, logical(1), x = amendment, fixed = TRUE)))

sha256 <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  sub("[[:space:]].*$", "", out[[1L]])
}

manifest_path <- file.path("data", "clean", "sample_manifest_frozen.tsv")
target_path <- file.path("data", "clean", "bfff_targets_sensitivity.tsv")
stopifnot(
  sha256(manifest_path) == "5434ee45d31e80a38db81696118b823e45b5aab2cd2dc628269ab29ee417a859",
  sha256(target_path) == "573c5f577dc071af45289ff68e051d762f7b4403288dc78267e62dab8755abcd"
)

manifest <- read.delim(
  manifest_path, sep = "\t", quote = "", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)
targets <- read.delim(target_path, sep = "\t", quote = "", check.names = FALSE)
stopifnot(nrow(targets) == 1270L, !anyDuplicated(targets$entrez_id))

gse188427_all <- manifest[manifest$series_accession == "GSE188427", , drop = FALSE]
gse188427 <- gse188427_all[gse188427_all$include, , drop = FALSE]
stopifnot(
  nrow(gse188427) == 198L,
  sum(gse188427$status == "RSV_case") == 147L,
  sum(gse188427$status == "healthy_control") == 51L,
  sum(gse188427_all$needs_manual_review) == 7L,
  sum(gse188427$needs_manual_review) == 3L,
  all(is.na(gse188427$age_months)),
  all(is.na(gse188427$sex)),
  length(unique(gse188427$batch)) == nrow(gse188427)
)

gse105450 <- manifest[
  manifest$series_accession == "GSE105450" & manifest$include &
    manifest$severity %in% c("hospitalized", "healthy"),
  , drop = FALSE
]
gse105450$case_status <- factor(
  ifelse(gse105450$status == "RSV_case", "RSV", "healthy"),
  levels = c("healthy", "RSV")
)
gse105450_design <- model.matrix(~ case_status + age_months + sex + batch, gse105450)
gse105450_expression <- readRDS(file.path("data", "derived", "GSE105450_expression.rds"))
stopifnot(
  nrow(gse105450) == 89L,
  sum(gse105450$case_status == "RSV") == 56L,
  sum(gse105450$case_status == "healthy") == 33L,
  all(gse105450$sample_id %in% colnames(gse105450_expression$expression)),
  qr(gse105450_design)$rank == ncol(gse105450_design)
)

gse103119 <- manifest[manifest$series_accession == "GSE103119" & manifest$include, , drop = FALSE]
gse103119$case_status <- factor(
  ifelse(gse103119$status == "RSV_case", "RSV", "healthy"),
  levels = c("healthy", "RSV")
)
gse103119_design <- model.matrix(~ case_status + age_months + sex, gse103119)
stopifnot(
  nrow(gse103119) == 31L,
  sum(gse103119$case_status == "RSV") == 11L,
  sum(gse103119$case_status == "healthy") == 20L,
  all(gse103119$pathogen[gse103119$case_status == "RSV"] == "RSV"),
  qr(gse103119_design)$rank == ncol(gse103119_design),
  file.exists(file.path("data", "raw", "GSE103119", "GSE103119_series_matrix.txt.gz")),
  file.exists(file.path("data", "raw", "GSE103119", "GSE103119_non-normalized.txt.gz")),
  file.exists(file.path("data", "raw", "GPL10558", "GPL10558.annot.gz"))
)

gse155925 <- manifest[manifest$series_accession == "GSE155925" & manifest$include, , drop = FALSE]
gse155925$virus_group <- factor(
  ifelse(gse155925$status == "RSV_case", "RSV", "other"),
  levels = c("other", "RSV")
)
batch_parts <- do.call(rbind, strsplit(gse155925$batch, ";", fixed = TRUE))
gse155925$hospital_batch <- factor(batch_parts[, 1L])
gse155925$enrollment_batch <- factor(batch_parts[, 2L])
gse155925_design <- model.matrix(
  ~ virus_group + age_months + sex + hospital_batch + enrollment_batch,
  gse155925
)
stopifnot(
  nrow(gse155925) == 48L,
  sum(gse155925$virus_group == "RSV") == 31L,
  sum(gse155925$virus_group == "other") == 17L,
  all(table(gse155925$virus_group, gse155925$hospital_batch) > 0L),
  all(table(gse155925$virus_group, gse155925$enrollment_batch) > 0L),
  qr(gse155925_design)$rank == ncol(gse155925_design)
)

count_path <- file.path("data", "raw", "GSE155925", "GSE155925_Raw_counts_matrix.txt.gz")
count_header <- names(read.delim(count_path, nrows = 0L, sep = "\t", quote = "", check.names = FALSE))[-1L]
count_lines <- readLines(gzfile(count_path), warn = FALSE)
normalized_ensembl <- sub("[.][0-9]+$", "", sub("^.*:", "", sub("\t.*$", "", count_lines[-1L])))
gse155925_all <- manifest[manifest$series_accession == "GSE155925", , drop = FALSE]
stopifnot(
  identical(count_header, gse155925_all$title),
  all(grepl("^ENSG[0-9]+$", normalized_ensembl)),
  !anyDuplicated(normalized_ensembl)
)

# RNA-seq target identity contract: ambiguous aliases are removed globally,
# the highest blind mean logCPM wins, ties use Ensembl lexical order, and
# primary/expanded indices are subsets of one frozen selection table.
target_fixture <- data.frame(
  entrez_id = c("1", "2", "3", "4"),
  ensembl_gene_ids = c("ENSG00000000001.5;ENSG00000000002", "ENSG00000000003", "ENSG00000000004", "ENSG00000000004"),
  stringsAsFactors = FALSE
)
background_fixture <- sprintf("ENSG%011d", 1:4)
selection_fixture <- task11_select_rnaseq_target_units(
  filtered_ensembl = background_fixture,
  mean_logcpm = c(7, 7, 4, 20),
  targets = target_fixture
)
stopifnot(
  identical(selection_fixture$ambiguous_ensembl, "ENSG00000000004"),
  identical(selection_fixture$selection$entrez_id, c("1", "2")),
  identical(selection_fixture$selection$ensembl_id, c("ENSG00000000001", "ENSG00000000003"))
)
primary_index <- task11_target_index(selection_fixture$selection, c("1"), background_fixture)
expanded_index <- task11_target_index(selection_fixture$selection, c("1", "2", "3", "4"), background_fixture)
stopifnot(
  primary_index$coverage_unique_entrez == length(primary_index$index),
  expanded_index$coverage_unique_entrez == length(expanded_index$index),
  all(primary_index$index %in% expanded_index$index),
  !anyDuplicated(expanded_index$selected$entrez_id)
)
actual_alias_map <- task11_build_alias_map(targets)
stopifnot(identical(
  actual_alias_map$ambiguous_ensembl,
  c("ENSG00000223572", "ENSG00000237289")
))
rnaseq_background_body <- paste(deparse(body(task11_prepare_rnaseq_background)), collapse = " ")
stopifnot(
  grepl("filterByExpr(dge, design = design)", rnaseq_background_body, fixed = TRUE),
  grepl("keep.lib.sizes = FALSE", rnaseq_background_body, fixed = TRUE),
  grepl("calcNormFactors(dge, method = \"TMM\")", rnaseq_background_body, fixed = TRUE),
  grepl("prior.count = 0.25", rnaseq_background_body, fixed = TRUE),
  grepl("normalized.lib.sizes = TRUE", rnaseq_background_body, fixed = TRUE)
)

# camera contract: execute the wrapper on synthetic data and require exact
# equality to the only allowed full-design camera call.
set.seed(20260805L)
camera_y <- matrix(rnorm(80L * 12L), nrow = 80L, dimnames = list(as.character(1:80), paste0("S", 1:12)))
camera_group <- factor(rep(c("healthy", "RSV"), each = 6L), levels = c("healthy", "RSV"))
camera_age <- rep(seq_len(6L), 2L)
camera_design <- model.matrix(~ camera_group + camera_age)
camera_contrast <- c(0, 1, 0)
camera_index <- 1:12
wrapped_camera <- task11_camera(camera_y, camera_index, camera_design, camera_contrast)
direct_camera <- limma::camera(
  y = camera_y, index = camera_index, design = camera_design, contrast = camera_contrast,
  directional = TRUE, inter.gene.cor = NA_real_
)
stopifnot(identical(wrapped_camera, direct_camera))
camera_body <- paste(deparse(body(task11_camera)), collapse = " ")
stopifnot(
  grepl("directional = TRUE", camera_body, fixed = TRUE),
  grepl("inter.gene.cor = NA_real_", camera_body, fixed = TRUE),
  !grepl("cameraPR", camera_body, fixed = TRUE),
  !grepl("0.01", camera_body, fixed = TRUE)
)

# Signed Stouffer contract: signs, sqrt(final n) weights, Z, and two-sided P.
stouffer <- task11_signed_stouffer(c(0.04, 0.10), c("Up", "Down"), c(89, 73))
expected_z_i <- c(
  stats::qnorm(0.04 / 2, lower.tail = FALSE),
  -stats::qnorm(0.10 / 2, lower.tail = FALSE)
)
expected_weight <- sqrt(c(89, 73))
expected_z <- sum(expected_weight * expected_z_i) / sqrt(sum(expected_weight^2))
expected_p <- 2 * stats::pnorm(abs(expected_z), lower.tail = FALSE)
stopifnot(
  isTRUE(all.equal(stouffer$z_i, expected_z_i, tolerance = 0)),
  isTRUE(all.equal(stouffer$weight, expected_weight, tolerance = 0)),
  isTRUE(all.equal(stouffer$z_combined, expected_z, tolerance = 0)),
  isTRUE(all.equal(stouffer$p_two_sided, expected_p, tolerance = 0))
)

# GSE188427 blind-QC contract: deterministic variance/tie ordering, exactly
# 5,000 genes, intercept-only arrayWeights, exact PC distance, no phenotype.
tie_expression <- rbind(
  ENSG00000000002 = c(0, 1, 2, 3),
  ENSG00000000001 = c(0, 1, 2, 3),
  ENSG00000000003 = c(0, 0, 0, 1)
)
stopifnot(identical(
  task11_top_variable_gene_ids(tie_expression, 2L),
  c("ENSG00000000001", "ENSG00000000002")
))
set.seed(20260806L)
qc_expression <- matrix(
  rnorm(5005L * 16L), nrow = 5005L,
  dimnames = list(as.character(seq_len(5005L)), paste0("Q", seq_len(16L)))
)
qc_result <- task11_gse188427_blind_qc(qc_expression)
expected_top <- task11_top_variable_gene_ids(qc_expression, 5000L)
expected_weights <- limma::arrayWeights(qc_expression, design = matrix(1, nrow = 16L, ncol = 1L))
expected_pca <- stats::prcomp(t(qc_expression[expected_top, , drop = FALSE]), center = TRUE, scale. = FALSE)
expected_standardized <- sweep(expected_pca$x[, 1:5, drop = FALSE], 2L, expected_pca$sdev[1:5], "/")
expected_distance <- sqrt(rowSums(expected_standardized^2))
stopifnot(
  length(qc_result$top_gene_ids) == 5000L,
  identical(qc_result$top_gene_ids, expected_top),
  isTRUE(all.equal(qc_result$qc$log_array_weight, log(as.numeric(expected_weights)), tolerance = 1e-12)),
  isTRUE(all.equal(qc_result$qc$pca_distance_5pc, unname(expected_distance[qc_result$qc$sample_id]), tolerance = 1e-12)),
  !any(names(qc_result$qc) %in% c("case_status", "status", "severity", "arm", "pathogen"))
)

cat("Task 11 operational design contract passed; no effect/result files read\n")
