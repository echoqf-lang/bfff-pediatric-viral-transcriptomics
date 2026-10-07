library(testthat)

test_file_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(test_file_dir, "..", ".."), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_root, "analysis", "R", "03_audit_samples.R"))) {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
old_wd <- getwd()
setwd(project_root)
old_skip <- Sys.getenv("TASK4_SKIP_MAIN", unset = NA_character_)
Sys.setenv(TASK4_SKIP_MAIN = "1")
source(file.path("analysis", "R", "03_audit_samples.R"), local = FALSE)
if (is.na(old_skip)) Sys.unsetenv("TASK4_SKIP_MAIN") else Sys.setenv(TASK4_SKIP_MAIN = old_skip)
setwd(old_wd)

raw <- read.delim(
  file.path(project_root, "data", "metadata", "all_samples_raw.tsv"),
  check.names = FALSE, quote = "\"", stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
manifest <- audit_samples(raw)

test_that("all official samples are represented exactly once", {
  expect_equal(nrow(manifest), 1031L)
  expect_equal(anyDuplicated(manifest$sample_id), 0L)
  expect_setequal(manifest$sample_id, raw$geo_accession)
})

test_that("confirmatory cohorts contain only independent acute RSV and healthy children", {
  g1 <- subset(manifest, series_accession == "GSE105450" & include)
  g2 <- subset(manifest, series_accession == "GSE103842" & include)
  expect_equal(as.integer(table(g1$status)), c(33L, 89L))
  expect_equal(as.integer(table(g1$severity[g1$status == "RSV_case"])), c(56L, 33L))
  expect_equal(as.integer(table(g2$status)), c(12L, 62L))
  expect_true(all(g1$timepoint == "acute_primary"))
  expect_true(all(g2$timepoint == "acute_primary"))
  expect_true(all(g1$tissue == "whole_blood"))
  expect_true(all(g2$tissue == "whole_blood"))
})

test_that("five GSE105450 healthy technical repeats are identified without expression", {
  z <- subset(manifest, series_accession == "GSE105450")
  excluded <- subset(z, exclusion_reason == "technical_replicate")
  expect_equal(nrow(excluded), 5L)
  expect_true(all(excluded$status == "healthy_control"))
  paired <- table(z$subject_id[z$replicate_type == "technical"])
  expect_equal(length(paired), 5L)
  expect_true(all(paired == 2L))
})

test_that("longitudinal and exploratory restrictions are enforced", {
  g188 <- subset(manifest, series_accession == "GSE188427")
  expect_equal(sum(g188$include), 198L)
  expect_equal(as.integer(table(g188$status[g188$include])), c(51L, 147L))
  expect_equal(as.integer(table(g188$severity[g188$include & g188$status == "RSV_case"])), c(110L, 34L))
  expect_equal(sum(g188$include & g188$status == "RSV_case" & is.na(g188$severity)), 3L)
  expect_true(all(!g188$include[g188$timepoint %in% c("day_30", "day_180")]))
  expect_true(all(table(g188$subject_id[g188$status == "RSV_case"]) <= 3L))
  expect_equal(sum(!is.na(g188$metadata_discrepancy)), 7L)
  expect_equal(sum(g188$include & !is.na(g188$metadata_discrepancy)), 3L)
  expect_true(all(g188$needs_manual_review[!is.na(g188$metadata_discrepancy)]))
  expect_true(all(is.na(g188$severity[!is.na(g188$metadata_discrepancy)])))
  expect_false(any(g188$severity_analysis_eligible[!is.na(g188$metadata_discrepancy)]))
  expect_false(any(g188$severity_analysis_eligible & g188$needs_manual_review))
  expect_true(all(g188$metadata_discrepancy[!is.na(g188$metadata_discrepancy)] ==
                    "source_name arm IN conflicts with title and characteristics arm OUT; severity unresolved"))
  subject_n <- table(g188$subject_id)
  repeated <- unname(subject_n[g188$subject_id]) > 1L
  expect_true(all(g188$replicate_type[repeated] == "longitudinal"))
  expect_true(all(g188$replicate_type[!repeated] == "none"))
  expect_true(all(!is.na(g188$replicate_group[repeated])))
  expect_true(all(is.na(g188$replicate_group[!repeated])))

  g389 <- subset(manifest, series_accession == "GSE38900")
  expect_equal(sum(g389$exclusion_reason == "recovery_timepoint", na.rm = TRUE), 21L)
  expect_true(all(g389$center[g389$include] %in% c("Dallas", "Turku")))
  expect_true(all(g389$status[g389$include] %in% c("RSV_case", "healthy_control")))
  expect_true(all(g389$outcome_previously_observed[g389$include]))
  expect_false(any(g389$center == "Columbus_NCH" & g389$include))

  g389_raw <- raw[match(g389$sample_id, raw$geo_accession), , drop = FALSE]
  local_token <- gsub("[^0-9]", "", sub("RSV.*$", "", g389_raw$description))
  token_centers <- split(g389$center, local_token)
  cross_center_tokens <- names(token_centers)[vapply(
    token_centers, function(x) length(unique(na.omit(x))) > 1L, logical(1)
  )]
  cross_rows <- local_token %in% cross_center_tokens
  expect_equal(length(cross_center_tokens), 6L)
  expect_equal(sum(cross_rows), 12L)
  expect_true(all(vapply(split(g389$subject_id[cross_rows], local_token[cross_rows]),
                         function(x) length(unique(x)) == 2L, logical(1))))
  expect_equal(sum(cross_rows & g389$center == "Dallas" & g389$include), 6L)

  subject_n389 <- table(g389$subject_id)
  paired_ids <- names(subject_n389[subject_n389 > 1L])
  expect_equal(length(paired_ids), 12L)
  expect_true(all(vapply(paired_ids, function(id) {
    z <- g389[g389$subject_id == id, ]
    setequal(z$status, c("RSV_case", "recovery")) && nrow(z) == 2L && all(z$center == "Dallas")
  }, logical(1))))
  expect_true(all(g389$replicate_type[g389$subject_id %in% paired_ids] == "longitudinal"))
})

test_that("exploratory cohorts do not manufacture healthy controls", {
  g155 <- subset(manifest, series_accession == "GSE155925")
  expect_false(any(g155$status == "healthy_control"))
  expect_false(any(g155$include & g155$status == "virus_negative_symptomatic"))
  expect_false(any(g155$include & g155$coinfection != "none"))

  g103 <- subset(manifest, series_accession == "GSE103119")
  expect_equal(sum(g103$include & g103$status == "RSV_case"), 11L)
  expect_equal(sum(g103$include & g103$status == "healthy_control"), 20L)
  expect_true(all(g103$pathogen[g103$include & g103$status == "RSV_case"] == "RSV"))
})

test_that("confirmatory meta membership cannot leak substitution or observed cohorts", {
  eligible <- subset(manifest, meta_eligible_confirmatory)
  expect_setequal(unique(eligible$series_accession), c("GSE105450", "GSE103842"))
  expect_false(any(eligible$outcome_previously_observed))
  expect_false(any(eligible$series_accession %in% c("GSE188427", "GSE38900")))
})

test_that("sex is normalized to the frozen vocabulary", {
  expect_true(all(is.na(manifest$sex) | manifest$sex %in% c("male", "female")))
  expect_setequal(unique(na.omit(manifest$sex)), c("male", "female"))
})

test_that("independent review packet is evidence-complete and unsigned by the analyst", {
  packet <- build_review_packet(raw, manifest)
  required <- c(
    "sample_id", "series_accession", "platform", "raw_title", "raw_source_name",
    "raw_description", paste0("raw_characteristic_", 1:9), "rule_id",
    "machine_include", "machine_exclusion_reason", "machine_status",
    "machine_timepoint", "machine_tissue", "machine_virus", "machine_severity",
    "metadata_discrepancy", "subject_id", "subject_id_evidence", "review_status",
    "reviewer", "reviewed_at", "review_decision", "final_include",
    "final_exclusion_reason", "review_notes", "packet_schema_version"
  )
  expect_identical(names(packet), required)
  expect_equal(nrow(packet), 1031L)
  expect_equal(anyDuplicated(packet$sample_id), 0L)
  expect_true(all(packet$review_status == "pending_independent_review"))
  expect_true(all(is.na(packet$reviewer)))
  expect_true(all(is.na(packet$reviewed_at)))
  expect_true(all(is.na(packet$review_decision)))
  expect_true(all(is.na(packet$final_include)))
  expect_true(all(!is.na(packet$rule_id) & !is.na(packet$subject_id_evidence)))
  expect_false(validate_independent_review(packet)$complete)
  blank_status <- packet
  blank_status$review_status[[1L]] <- NA_character_
  expect_false(validate_independent_review(blank_status)$complete)

  packet_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(packet_path), add = TRUE)
  write_tsv(packet, packet_path)
  expect_identical(preserve_existing_reviews(packet, packet_path), packet)
  changed <- packet
  changed$raw_title[[1L]] <- paste0(changed$raw_title[[1L]], "_changed")
  expect_identical(preserve_existing_reviews(changed, packet_path), changed)

  partially_signed <- packet
  partially_signed$review_status[[1L]] <- "completed_independent_review"
  partially_signed$reviewer[[1L]] <- "independent_reviewer_01"
  partially_signed$reviewed_at[[1L]] <- "2026-08-05T12:00:00+0800"
  partially_signed$review_decision[[1L]] <- "agree"
  partially_signed$final_include[[1L]] <- partially_signed$machine_include[[1L]]
  partially_signed$final_exclusion_reason[[1L]] <- partially_signed$machine_exclusion_reason[[1L]]
  write_tsv(partially_signed, packet_path)
  expect_error(preserve_existing_reviews(changed, packet_path), "Machine evidence changed")

  signed <- packet
  signed$review_status <- "completed_independent_review"
  signed$reviewer <- "independent_reviewer_01"
  signed$reviewed_at <- "2026-08-05T12:00:00+0800"
  signed$review_decision <- "agree"
  signed$final_include <- signed$machine_include
  signed$final_exclusion_reason <- signed$machine_exclusion_reason
  expect_true(validate_independent_review(signed)$complete)
  expect_equal(nrow(apply_completed_review(manifest, signed)), 1031L)

  signed$reviewer[[1L]] <- "sample_audit_analyst"
  expect_false(validate_independent_review(signed)$complete)
  signed$reviewer[[1L]] <- "independent_reviewer_01"
  excluded_row <- which(!signed$machine_include)[[1L]]
  signed$review_decision[[excluded_row]] <- "disagree"
  signed$final_include[[excluded_row]] <- TRUE
  signed$final_exclusion_reason[[excluded_row]] <- NA_character_
  expect_false(validate_independent_review(signed)$complete)
})

test_that("series-matrix validation reads header only", {
  files <- Sys.glob(file.path(project_root, "data", "raw", "*", "*series_matrix*.gz"))
  expect_length(files, 7L)
  headers <- lapply(files, read_series_matrix_header)
  expect_true(all(vapply(headers, `[[`, logical(1), "stopped_before_table")))
  expect_true(all(vapply(headers, function(x) length(x$sample_ids) > 0L, logical(1))))
  for (series in unique(raw$series_accession)) {
    relevant <- vapply(headers, function(h) identical(h$accession, series), logical(1))
    ids <- unique(unlist(lapply(headers[relevant], `[[`, "sample_ids"), use.names = FALSE))
    expect_setequal(ids, raw$geo_accession[raw$series_accession == series])
    quick <- file.path(project_root, "data", "raw", series, paste0(series, "_quick.soft.txt"))
    expect_identical(read_quick_soft_accession(quick), series)
  }
})

test_that("committed independent review freezes provenance without weakening restrictions", {
  packet_path <- file.path(project_root, "data", "clean", "sample_review_packet.tsv")
  signed <- read.delim(packet_path, check.names = FALSE, quote = "\"", stringsAsFactors = FALSE,
                       na.strings = "", fileEncoding = "UTF-8")
  signed[setdiff(review_columns, "final_include")] <- lapply(
    signed[setdiff(review_columns, "final_include")], as.character
  )
  signed$final_include <- as.logical(signed$final_include)
  expect_identical(sha256_file(packet_path),
                   "140239a08915dc88872e7dedfaebe366a00a104ca3c24b936c0861a9229ffdd6")
  expect_true(validate_independent_review(signed)$complete)
  expect_equal(nrow(signed), 1031L)
  expect_true(all(signed$reviewer == "independent_sample_reviewer"))
  expect_true(all(signed$review_decision == "agree"))
  expect_identical(signed$final_include, signed$machine_include)

  final <- apply_completed_review(manifest, signed)
  expect_equal(nrow(final), 1031L)
  expect_equal(sum(final$include), 611L)
  expect_equal(sum(final$meta_eligible_confirmatory), 196L)
  expect_equal(length(unique(final$subject_id[final$meta_eligible_confirmatory])), 196L)
  expect_true(all(final$reviewer == "independent_sample_reviewer"))
  expect_true(all(final$freeze_status == "frozen_after_complete_independent_review"))
  conflicts <- final$series_accession == "GSE188427" & final$needs_manual_review
  expect_equal(sum(conflicts), 7L)
  expect_true(all(is.na(final$severity[conflicts])))
  expect_false(any(final$severity_analysis_eligible[conflicts]))
  expect_true(all(grepl("^GSE38900:[^:]+:[^:]+$",
                        final$subject_id[final$series_accession == "GSE38900"])))
  expect_true(all(final$processed_only_qc_boundary ==
                    "metadata_frozen; CEL/FASTQ-level QC not claimed"))
})
