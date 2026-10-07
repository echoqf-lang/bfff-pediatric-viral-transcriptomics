#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

read_tsv <- function(path) read.delim(
  path, sep = "\t", quote = "", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

sha256 <- function(path) sub(
  "[[:space:]].*$", "",
  system2("shasum", c("-a", "256", path), stdout = TRUE)[[1L]]
)

paths <- c(
  script = file.path("analysis", "R", "11_make_tables_figures.R"),
  table1 = file.path("results", "tables", "Table1_cohorts.tsv"),
  table2 = file.path("results", "tables", "Table2_target_set_tests.tsv"),
  table3 = file.path("results", "tables", "Table3_target_meta.tsv"),
  table_s = file.path("results", "tables", "TableS_all_effects.tsv"),
  fig1 = file.path("results", "figures", "Fig1_workflow.pdf"),
  fig2 = file.path("results", "figures", "Fig2_target_set_replication.pdf"),
  fig_s1 = file.path("results", "figures", "FigS1_all_target_effects.pdf"),
  fig4 = file.path("results", "figures", "Fig4_heterogeneity.pdf"),
  fig5 = file.path("results", "figures", "Fig5_sensitivity.pdf"),
  manifest = file.path("results", "reproducibility_manifest.tsv")
)
stopifnot(all(file.exists(paths)), all(file.info(paths)$size > 0))
stopifnot(!file.exists(file.path("results", "figures", "Fig3_meta_forest.pdf")))

table1 <- read_tsv(paths[["table1"]])
required_table1 <- c(
  "cohort", "analysis_role", "platform", "tissue", "timepoint", "contrast",
  "official_n", "included_n", "analysis_n", "reference_n", "case_n",
  "summary_population_n", "summary_population_definition",
  "age_months_summary", "age_missing_n", "sex_summary", "sex_missing_n",
  "severity_summary", "severity_missing_n", "covariates", "inclusion_summary",
  "exclusion_summary", "analysis_boundary"
)
stopifnot(
  nrow(table1) == 6L,
  setequal(table1$cohort, c("GSE103119", "GSE103842", "GSE105450", "GSE155925", "GSE188427", "GSE38900")),
  all(required_table1 %in% names(table1)),
  table1$analysis_n[table1$cohort == "GSE103842"] == 73L,
  table1$analysis_n[table1$cohort == "GSE188427"] == 195L,
  table1$summary_population_n[table1$cohort == "GSE103842"] == 73L,
  table1$summary_population_n[table1$cohort == "GSE188427"] == 195L,
  is.na(table1$analysis_n[table1$cohort == "GSE38900"]),
  table1$summary_population_n[table1$cohort == "GSE38900"] == 138L,
  grepl("included_pre-QC", table1$summary_population_definition[table1$cohort == "GSE38900"], fixed = TRUE),
  all(grepl("analysis", table1$summary_population_definition[table1$cohort != "GSE38900"], fixed = TRUE)),
  grepl("previously observed", table1$analysis_boundary[table1$cohort == "GSE38900"], fixed = TRUE),
  grepl("unadjusted", table1$analysis_boundary[table1$cohort == "GSE188427"], fixed = TRUE)
)

table2 <- read_tsv(paths[["table2"]])
stopifnot(
  setequal(table2$analysis_id, c(
    "GSE105450_camera", "GSE103842_camera", "weighted_signed_stouffer",
    "matched_random_empirical_calibration"
  )),
  all(table2$primary_H1_supported == FALSE),
  all(table2$primary_H1_rescue_allowed == FALSE),
  all(table2$direction[table2$analysis_id %in% c(
    "GSE105450_camera", "GSE103842_camera", "weighted_signed_stouffer"
  )] == "Up"),
  all(table2$p_value[table2$analysis_id %in% c("GSE105450_camera", "GSE103842_camera")] >= 0.05),
  table2$p_value[table2$analysis_id == "weighted_signed_stouffer"] >= 0.05,
  table2$p_value[table2$analysis_id == "matched_random_empirical_calibration"] < 0.05,
  !any(c("ci_lb", "ci_ub") %in% names(table2)),
  all(table2$ci_status %in% c("not_estimable", "not_defined")),
  all(nzchar(table2$ci_reason)),
  all(table2$primary_decision == "H1_not_supported_camera_component_failed"),
  grepl("conditional membership-randomization sensitivity", table2$method[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE),
  grepl("after primary outcomes were known", table2$evidence_timing[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE),
  grepl("before any random draws", table2$evidence_timing[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE),
  grepl("different from the camera competitive null", table2$null_hypothesis_boundary[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE),
  grepl("narrow conditional random-null distribution", table2$interpretation[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE),
  grepl("cannot validate camera", table2$interpretation[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE),
  grepl("independent replication", table2$interpretation[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE),
  grepl("roast", table2$supplementary_signal_context[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE),
  grepl("fgsea", table2$supplementary_signal_context[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE),
  grepl("broad RSV whole-blood response or cell-composition shifts", table2$supplementary_signal_context[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE),
  grepl("do not establish competitive target-set enrichment or formula mechanism", table2$supplementary_signal_context[table2$analysis_id == "matched_random_empirical_calibration"], fixed = TRUE)
)

table3 <- read_tsv(paths[["table3"]])
stopifnot(
  nrow(table3) == 358L,
  all(table3$k == 2L),
  all(table3$method_primary == "REML_HKSJ"),
  sum(table3$hksj_fdr_bh < 0.05) == 0L,
  all(grepl("unstable", table3$heterogeneity_caveat, fixed = TRUE)),
  all(table3$primary_H1_rescue_allowed == FALSE)
)

table_s <- read_tsv(paths[["table_s"]])
stopifnot(
  nrow(table_s) == 520L,
  any(!table_s$meta_estimable),
  sum(table_s$meta_estimable) == 358L,
  identical(table_s$gene_symbol, sort(table_s$gene_symbol, na.last = TRUE)),
  all(c(
    "GSE105450_log2FC", "GSE105450_SE", "GSE103842_log2FC", "GSE103842_SE",
    "HKSJ_estimate", "HKSJ_ci_lb", "HKSJ_ci_ub", "HKSJ_prediction_lb",
    "HKSJ_prediction_ub", "HKSJ_tau2", "HKSJ_I2", "HKSJ_fdr_bh"
  ) %in% names(table_s))
)

pdf_info <- lapply(paths[grep("^fig", names(paths))], function(path) {
  output <- system2("pdfinfo", path, stdout = TRUE, stderr = TRUE)
  stopifnot(is.null(attr(output, "status")))
  output
})
stopifnot(
  all(vapply(pdf_info, function(x) any(grepl("^Pages:", x)), logical(1))),
  as.integer(sub("^Pages:[[:space:]]+", "", grep("^Pages:", pdf_info[["fig_s1"]], value = TRUE))) >= 9L
)

manifest <- read_tsv(paths[["manifest"]])
stopifnot(
  all(c("record_type", "name", "path", "sha256", "version", "seed", "note") %in% names(manifest)),
  all(c("input", "archived_contextual_dependency", "script", "package", "seed", "output") %in% unique(manifest$record_type)),
  any(manifest$record_type == "seed" & manifest$seed == 20260804L),
  any(manifest$record_type == "seed" & manifest$seed == 20260812L),
  any(manifest$record_type == "output" & manifest$path == paths[["manifest"]] & is.na(manifest$sha256)),
  any(manifest$record_type == "output" & manifest$name == "Fig3_meta_forest.pdf" & is.na(manifest$path) &
        manifest$note == "cancelled_not_generated: zero HKSJ BH-FDR targets; no non-significant top-20 filling"),
  sum(manifest$record_type == "archived_contextual_dependency") == 3L
)
task10_timing_rows <- manifest$record_type == "input" & manifest$name %in% c(
  "2026-08-05-random-set-cross-cohort-operationalization.md",
  "2026-08-05-random-set-quintile-ties.md"
)
stopifnot(
  sum(task10_timing_rows) == 2L,
  all(grepl("post-outcome operational clarification", manifest$note[task10_timing_rows], fixed = TRUE)),
  all(grepl("before random draws", manifest$note[task10_timing_rows], fixed = TRUE))
)
locked <- manifest$record_type %in% c("input", "archived_contextual_dependency", "script", "output") & !is.na(manifest$sha256)
stopifnot(
  all(file.exists(manifest$path[locked])),
  identical(
    unname(vapply(manifest$path[locked], sha256, character(1))),
    unname(manifest$sha256[locked])
  )
)

# Regression guard: every plotting helper must close its explicit device and
# must not fall through to the default Rplots.pdf device while restoring par().
old_skip <- Sys.getenv("TASK13_SKIP_MAIN", unset = NA_character_)
Sys.setenv(TASK13_SKIP_MAIN = "1")
source(paths[["script"]], local = FALSE)
if (is.na(old_skip)) Sys.unsetenv("TASK13_SKIP_MAIN") else Sys.setenv(TASK13_SKIP_MAIN = old_skip)
stopifnot(grDevices::dev.cur() == 1L)
probe_pdf <- tempfile(fileext = ".pdf")
task13_plot_workflow(probe_pdf)
stopifnot(grDevices::dev.cur() == 1L, file.exists(probe_pdf), file.info(probe_pdf)$size > 0)
unlink(probe_pdf)

cat("Task 13 tables, figures, and reproducibility manifest tests passed\n")
