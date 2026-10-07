#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

read_tsv_local <- function(path) {
  read.delim(path, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
}

paths <- c(
  camera = file.path("results", "cohort", "target_set_camera.tsv"),
  roast = file.path("results", "cohort", "target_set_roast.tsv"),
  fgsea = file.path("results", "cohort", "target_set_fgsea_exploratory.tsv"),
  decision = file.path("results", "cohort", "target_set_primary_decision_partial.tsv")
)
stopifnot(all(file.exists(paths)))
camera <- read_tsv_local(paths[["camera"]])
roast <- read_tsv_local(paths[["roast"]])
fgsea <- read_tsv_local(paths[["fgsea"]])
decision <- read_tsv_local(paths[["decision"]])
cohorts <- c("GSE105450", "GSE103842")
targets <- read_tsv_local(file.path("data", "clean", "bfff_targets_primary.tsv"))
target_ids <- as.character(targets$entrez_id)

stopifnot(
  nrow(targets) == 520L,
  !anyDuplicated(target_ids),
  identical(as.character(camera$cohort), cohorts),
  all(camera$status == "interpretable"),
  all(camera$n_frozen == 520L),
  all(camera$n_detectable >= 10L),
  all(camera$coverage >= 0.50),
  all(camera$direction %in% c("Up", "Down")),
  all(is.finite(camera$p_value) & camera$p_value > 0 & camera$p_value <= 1),
  all(is.finite(camera$inter_gene_correlation)),
  all(camera$camera_p_sidedness == "two_sided"),
  all(camera$background == "all_detectable_unique_Entrez_genes_in_cohort"),
  all(camera$model_coefficient == "case_statusRSV"),
  identical(as.integer(camera$n_independent_model), c(122L, 73L)),
  isTRUE(all.equal(camera$stouffer_weight_sqrt_n, sqrt(c(122, 73)), tolerance = 1e-14)),
  all(camera$task7_effect_hash_match),
  all(camera$task7_effects_sha256 == camera$reconstructed_effects_sha256)
)
expected_signed_z <- ifelse(
  camera$direction == "Up",
  qnorm(camera$p_value / 2, lower.tail = FALSE),
  -qnorm(camera$p_value / 2, lower.tail = FALSE)
)
stopifnot(isTRUE(all.equal(camera$signed_z, expected_signed_z, tolerance = 1e-14)))
for (i in seq_along(cohorts)) {
  object <- readRDS(file.path("data", "derived", paste0(cohorts[[i]], "_expression.rds")))
  expected_detectable <- sum(target_ids %in% rownames(object$expression))
  stopifnot(
    camera$n_detectable[[i]] == expected_detectable,
    isTRUE(all.equal(camera$coverage[[i]], expected_detectable / 520, tolerance = 1e-15))
  )
}

stopifnot(
  nrow(roast) == 6L,
  identical(unique(as.character(roast$cohort)), cohorts),
  all(roast$n_frozen == 520L),
  all(roast$nrot == 99999L),
  all(roast$seed == 20260804L),
  all(roast$role == "key_supplementary_self_contained"),
  all(is.finite(roast$p_value) & roast$p_value >= 0 & roast$p_value <= 1),
  all(is.finite(roast$p_holm) & roast$p_holm >= 0 & roast$p_holm <= 1)
)
for (cohort in cohorts) {
  one <- roast[roast$cohort == cohort, , drop = FALSE]
  stopifnot(
    identical(as.character(one$test_direction), c("Down", "Up", "Mixed")),
    isTRUE(all.equal(one$p_holm, p.adjust(one$p_value, method = "holm"), tolerance = 1e-15))
  )
}

stopifnot(
  nrow(fgsea) == 2L,
  identical(as.character(fgsea$cohort), cohorts),
  all(fgsea$ranking_statistic == "limma_moderated_t"),
  all(fgsea$n_frozen == 520L),
  all(fgsea$seed == 20260804L),
  all(fgsea$role == "exploratory_ranked_sensitivity_not_independent_validation"),
  all(is.finite(fgsea$p_value) & fgsea$p_value >= 0 & fgsea$p_value <= 1),
  all(is.finite(fgsea$ES)),
  all(is.finite(fgsea$NES))
)
for (i in seq_len(nrow(fgsea))) {
  leading <- strsplit(fgsea$leading_edge_entrez[[i]], "|", fixed = TRUE)[[1L]]
  leading <- leading[nzchar(leading)]
  stopifnot(all(leading %in% target_ids))
}

weights <- sqrt(c(122, 73))
combined_z <- sum(weights * camera$signed_z) / sqrt(sum(weights^2))
combined_p <- 2 * pnorm(-abs(combined_z))
direction_consistent <- length(unique(camera$direction)) == 1L
camera_component <- all(camera$p_value < 0.05) && direction_consistent && combined_p < 0.05
stopifnot(
  nrow(decision) == 1L,
  decision$both_cohorts_interpretable,
  decision$GSE105450_camera_p_lt_0_05 == (camera$p_value[[1L]] < 0.05),
  decision$GSE103842_camera_p_lt_0_05 == (camera$p_value[[2L]] < 0.05),
  decision$camera_direction_consistent == direction_consistent,
  isTRUE(all.equal(decision$weighted_stouffer_z, combined_z, tolerance = 1e-14)),
  isTRUE(all.equal(decision$weighted_stouffer_two_sided_p, combined_p, tolerance = 1e-14)),
  decision$camera_replication_component_met == camera_component,
  identical(decision$primary_decision, "H1_not_supported_camera_component_failed"),
  identical(decision$random_set_calibration_status, "pending_scheduled_nonrescuing"),
  is.na(decision$random_set_empirical_p_value),
  is.na(decision$random_set_empirical_p_lt_0_05),
  identical(decision$h1_supported, FALSE),
  identical(decision$pending_requirement, "Task10_matched_random_gene_set_empirical_p_value")
)

all_text <- paste(capture.output(lapply(c(camera, roast, fgsea, decision), print)), collapse = "\n")
stopifnot(!grepl("GSE38900|GSE188427", all_text, perl = TRUE))

checksum <- read_tsv_local(file.path("logs", "checksums", "task8_sha256.tsv"))
stopifnot(!anyDuplicated(checksum$artifact), !any(grepl("RAW|CEL|FASTQ|SRA", checksum$artifact, ignore.case = TRUE)))
for (i in seq_len(nrow(checksum))) {
  stopifnot(file.exists(checksum$artifact[[i]]))
  actual <- system2("shasum", c("-a", "256", checksum$artifact[[i]]), stdout = TRUE)
  actual <- sub("[[:space:]].*$", "", actual[[1L]])
  stopifnot(identical(actual, checksum$sha256[[i]]))
}

cat("task8 real-output verification passed\n")
