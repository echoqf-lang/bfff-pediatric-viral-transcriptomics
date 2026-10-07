#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

final_v6_scripts <- c(
  "analysis/R/12_gse77087_revised_main.R",
  "analysis/R/13_freeze_single_gene_upgrade.R",
  "analysis/R/14_three_cohort_gene_meta.R",
  "analysis/R/15_gse77087_severity_gradient.R",
  "analysis/R/16_gse155925_virus_specificity.R",
  "analysis/R/17_fixed_blood_modules.R",
  "analysis/R/18_cell_composition_sensitivity.R",
  "analysis/R/19_make_upgrade_tables_figures.R",
  "analysis/R/20_finalize_upgrade_audit.R",
  "analysis/R/21_acquire_airway_processed_geo.R",
  "analysis/R/22_airway_barrier_inflammation_repair.R",
  "analysis/R/23_bfff_target_panel_airway_validation.R",
  "analysis/R/24_bfff_full_86_airway_validation.R",
  "analysis/R/25_build_full_86_evidence_matrix.R",
  "analysis/R/27_make_v6_bmc_figures_tables.R",
  "analysis/R/28_validate_v6_biological_figure_inputs.R",
  "analysis/R/29_make_v6_biological_figures.R",
  "analysis/R/30_validate_v6_biological_figure_outputs.R"
)

missing_scripts <- final_v6_scripts[!file.exists(final_v6_scripts)]
if (length(missing_scripts)) {
  stop("Archive is missing required scripts: ", paste(missing_scripts, collapse = ", "), call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (!identical(args, "--execute")) {
  cat("BFFF V6 analysis entry point\n")
  cat("The DOI archive excludes downloaded GEO expression payloads.\n")
  cat("Read metadata/data_provenance.md and derived_data/model_inputs/README.md first.\n")
  cat("Scripts in frozen order:\n")
  cat(paste0("- ", final_v6_scripts, collapse = "\n"), "\n")
  cat("After restoring inputs and renv, run: Rscript run_analysis.R --execute\n")
  quit(status = 0L)
}

required_inputs <- c(
  "data/raw/GSE77087/GSE77087_series_matrix.txt.gz",
  "data/derived/GSE103842_expression.rds",
  "data/derived/barrier_inflammation_repair/GSE97742_processed_eset.rds",
  "data/derived/barrier_inflammation_repair/GSE41374_processed_eset.rds"
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop(
    "Processed GEO inputs are absent. Restore the documented inputs before execution: ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

for (script in final_v6_scripts) {
  message("Running ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (!identical(status, 0L)) stop("Analysis script failed: ", script, call. = FALSE)
}
