#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

pipeline_path <- file.path("analysis", "R", "statistical_pipeline.R")
if (!file.exists(pipeline_path)) stop("Frozen Task 6 statistical pipeline is absent", call. = FALSE)
source(pipeline_path, local = FALSE)

task7_cohorts <- c("GSE105450", "GSE103842")
task7_direction <- "RSV_minus_healthy"

sha256_file <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 failed for ", path, call. = FALSE)
  sub("[[:space:]].*$", "", output[[1L]])
}

read_tsv <- function(path) {
  read.delim(
    path, sep = "\t", quote = "", check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8"
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    na = "NA", fileEncoding = "UTF-8"
  )
}

assert_checksum_row <- function(checksum, artifact, expected_role = NULL) {
  row <- checksum[checksum$artifact == artifact, , drop = FALSE]
  if (nrow(row) != 1L) stop("Frozen checksum row absent or duplicated: ", artifact, call. = FALSE)
  if (!is.null(expected_role) && !identical(as.character(row$role), expected_role)) {
    stop("Unexpected frozen checksum role: ", artifact, call. = FALSE)
  }
  if (!file.exists(artifact)) stop("Frozen artifact absent: ", artifact, call. = FALSE)
  observed <- sha256_file(artifact)
  if (!identical(observed, as.character(row$sha256))) {
    stop("Frozen SHA-256 mismatch: ", artifact, call. = FALSE)
  }
  observed
}

verify_frozen_inputs <- function() {
  task5_path <- file.path("logs", "checksums", "task5_sha256.tsv")
  task6_path <- file.path("logs", "checksums", "task6_sha256.tsv")
  supersession_path <- file.path("logs", "session_info", "task5_supersession.md")
  if (!all(file.exists(c(task5_path, task6_path, supersession_path)))) {
    stop("Frozen Task 5/6 provenance is incomplete", call. = FALSE)
  }
  supersession <- paste(readLines(supersession_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (!grepl("only artifacts regenerated under the scaled-MAD amendment may be used", supersession, fixed = TRUE)) {
    stop("Scaled-MAD supersession boundary is not frozen", call. = FALSE)
  }
  task5 <- read_tsv(task5_path)
  task6 <- read_tsv(task6_path)
  artifacts <- c(
    file.path("data", "clean", "sample_manifest_frozen.tsv"),
    file.path("data", "derived", "GSE105450_expression.rds"),
    file.path("data", "derived", "GSE103842_expression.rds")
  )
  hashes <- setNames(vapply(
    artifacts, function(path) assert_checksum_row(task5, path), character(1)
  ), artifacts)
  hashes[[pipeline_path]] <- assert_checksum_row(task6, pipeline_path, "input")
  list(
    hashes = hashes,
    task5_checksum_sha256 = sha256_file(task5_path),
    task6_checksum_sha256 = sha256_file(task6_path),
    task5_supersession_sha256 = sha256_file(supersession_path)
  )
}

prepare_model_data <- function(expression, sample_data, cohort_id) {
  if (!is.matrix(expression) || !is.numeric(expression) || !length(expression)) {
    stop("Expression must be a non-empty numeric matrix", call. = FALSE)
  }
  if (is.null(rownames(expression)) || is.null(colnames(expression)) ||
      any(!nzchar(rownames(expression))) || any(!nzchar(colnames(expression))) ||
      anyDuplicated(rownames(expression)) || anyDuplicated(colnames(expression))) {
    stop("Expression requires unique non-empty gene and sample IDs", call. = FALSE)
  }
  if (any(!is.finite(expression))) stop("Expression contains non-finite values", call. = FALSE)
  required <- c(
    "sample_id", "subject_id", "cohort", "platform", "case_status",
    "age_months", "sex", "technical_batch"
  )
  absent <- setdiff(required, names(sample_data))
  if (length(absent)) stop("Sample data missing: ", paste(absent, collapse = ", "), call. = FALSE)
  if (!identical(as.character(sample_data$sample_id), colnames(expression))) {
    stop("Sample data must match expression columns in the same order", call. = FALSE)
  }
  if (anyDuplicated(sample_data$sample_id)) stop("Duplicate sample_id detected", call. = FALSE)
  if (anyDuplicated(sample_data$subject_id)) stop("Duplicate subject_id detected", call. = FALSE)
  if (!identical(unique(as.character(sample_data$cohort)), cohort_id)) {
    stop("Input must contain exactly the requested single cohort", call. = FALSE)
  }
  if (length(unique(as.character(sample_data$platform))) != 1L) {
    stop("Each fit must contain a single platform", call. = FALSE)
  }
  if (!setequal(unique(as.character(sample_data$case_status)), c("healthy", "RSV")) ||
      any(!sample_data$case_status %in% c("healthy", "RSV"))) {
    stop("case_status must contain only both healthy and RSV", call. = FALSE)
  }

  age <- suppressWarnings(as.numeric(sample_data$age_months))
  valid_age <- !is.na(age) & is.finite(age) & age >= 0
  valid_sex <- !is.na(sample_data$sex) & sample_data$sex %in% c("female", "male")
  complete <- valid_age & valid_sex
  if (!any(complete)) stop("No complete cases for age and sex", call. = FALSE)
  kept <- sample_data[complete, , drop = FALSE]
  kept$age_months <- age[complete]
  kept$case_status <- factor(kept$case_status, levels = c("healthy", "RSV"))
  kept$sex <- factor(kept$sex, levels = c("female", "male"))
  if (!setequal(unique(as.character(kept$case_status)), c("healthy", "RSV"))) {
    stop("Complete cases must retain both healthy and RSV", call. = FALSE)
  }
  if (nlevels(droplevels(kept$sex)) < 2L) stop("Sex is not estimable in complete cases", call. = FALSE)
  expr_kept <- expression[, complete, drop = FALSE]

  base_formula <- stats::as.formula("~ case_status + age_months + sex")
  base_design <- stats::model.matrix(base_formula, data = kept)
  if (qr(base_design)$rank != ncol(base_design)) {
    stop("Base age/sex-adjusted design matrix is not full rank", call. = FALSE)
  }

  batch_values <- as.character(kept$technical_batch)
  batch_available <- all(!is.na(batch_values) & nzchar(batch_values)) && length(unique(batch_values)) >= 2L
  batch_included <- FALSE
  batch_reason <- if (batch_available) "candidate_not_yet_evaluated" else "omitted_unavailable_or_single_level"
  design <- base_design
  model_formula <- base_formula
  batch_levels <- sort(unique(batch_values[!is.na(batch_values) & nzchar(batch_values)]), method = "radix")
  if (batch_available) {
    kept$technical_batch <- factor(batch_values, levels = batch_levels)
    candidate_formula <- stats::as.formula("~ case_status + age_months + sex + technical_batch")
    candidate <- stats::model.matrix(candidate_formula, data = kept)
    if (qr(candidate)$rank == ncol(candidate)) {
      design <- candidate
      model_formula <- candidate_formula
      batch_included <- TRUE
      batch_reason <- "included_available_estimable_not_status_collinear"
    } else {
      batch_reason <- "omitted_rank_deficient_or_status_collinear"
    }
  }
  if (qr(design)$rank != ncol(design)) stop("Selected design matrix is not full rank", call. = FALSE)
  if (!"case_statusRSV" %in% colnames(design)) stop("RSV-minus-healthy coefficient is absent", call. = FALSE)
  rownames(design) <- kept$sample_id

  list(
    expression = expr_kept,
    sample_data = kept,
    design = design,
    coefficient = "case_statusRSV",
    formula = paste(deparse(model_formula), collapse = ""),
    batch_included = batch_included,
    batch_reason = batch_reason,
    batch_levels = batch_levels,
    n_input = nrow(sample_data),
    n_complete = nrow(kept),
    n_excluded_missing_covariate = sum(!complete),
    n_missing_age = sum(!valid_age),
    n_missing_sex = sum(!valid_sex)
  )
}

fit_task7_cohort <- function(expression, sample_data, cohort_id,
                             direction = task7_direction) {
  if (!identical(direction, task7_direction)) stop("Direction must be RSV_minus_healthy", call. = FALSE)
  assert_frozen_packages()
  prepared <- prepare_model_data(expression, sample_data, cohort_id)
  fit <- limma::lmFit(prepared$expression, prepared$design)
  fit <- limma::eBayes(fit)
  coefficient_index <- match(prepared$coefficient, colnames(fit$coefficients))
  posterior_se <- fit$stdev.unscaled[, coefficient_index] * sqrt(fit$s2.post)
  if (any(!is.finite(posterior_se) | posterior_se <= 0)) {
    stop("Moderated standard errors must be finite and positive", call. = FALSE)
  }
  p_values <- unname(fit$p.value[, coefficient_index])
  gene_effects <- data.frame(
    gene_id = rownames(prepared$expression),
    cohort = cohort_id,
    direction = direction,
    log2FC = unname(fit$coefficients[, coefficient_index]),
    SE = unname(posterior_se),
    moderated_t = unname(fit$t[, coefficient_index]),
    p_value = p_values,
    fdr_bh = stats::p.adjust(p_values, method = "BH"),
    average_expression = rowMeans(prepared$expression),
    stringsAsFactors = FALSE
  )
  if (nrow(gene_effects) != nrow(expression) || anyDuplicated(gene_effects$gene_id)) {
    stop("Every detectable gene must be emitted exactly once", call. = FALSE)
  }
  c(
    list(gene_effects = gene_effects),
    prepared[c(
      "coefficient", "formula", "batch_included", "batch_reason", "batch_levels",
      "n_input", "n_complete", "n_excluded_missing_covariate", "n_missing_age", "n_missing_sex"
    )],
    list(
      design_rank = qr(prepared$design)$rank,
      design_columns = ncol(prepared$design),
      n_healthy = sum(prepared$sample_data$case_status == "healthy"),
      n_RSV = sum(prepared$sample_data$case_status == "RSV")
    )
  )
}

load_cohort_inputs <- function(cohort, manifest) {
  rds_path <- file.path("data", "derived", paste0(cohort, "_expression.rds"))
  object <- readRDS(rds_path)
  required_object <- c("cohort", "expression", "sample_data", "qc_excluded_sample_ids", "preprocessing")
  if (!all(required_object %in% names(object))) stop(cohort, ": incomplete derived object", call. = FALSE)
  if (!identical(as.character(object$cohort), cohort)) stop(cohort, ": cohort identity mismatch", call. = FALSE)
  expression <- object$expression
  blind <- object$sample_data
  if (!identical(colnames(expression), as.character(blind$sample_id))) {
    stop(cohort, ": blind sample IDs do not match expression columns", call. = FALSE)
  }
  if (anyDuplicated(manifest$sample_id)) stop("Frozen manifest sample_id is not unique", call. = FALSE)
  eligible <- manifest[
    manifest$series_accession == cohort & manifest$include & manifest$meta_eligible_confirmatory,
    , drop = FALSE
  ]
  modeled_ids <- as.character(blind$sample_id)
  excluded_ids <- as.character(object$qc_excluded_sample_ids)
  if (anyDuplicated(modeled_ids) || anyDuplicated(excluded_ids) || length(intersect(modeled_ids, excluded_ids))) {
    stop(cohort, ": retained/excluded sample identities are invalid", call. = FALSE)
  }
  if (!setequal(eligible$sample_id, c(modeled_ids, excluded_ids))) {
    stop(cohort, ": RDS retained plus QC-excluded samples do not equal frozen eligible manifest", call. = FALSE)
  }
  joined <- eligible[match(modeled_ids, eligible$sample_id), , drop = FALSE]
  if (anyNA(joined$sample_id) || !identical(as.character(joined$sample_id), modeled_ids)) {
    stop(cohort, ": manifest join is not strict one-to-one", call. = FALSE)
  }
  if (!identical(as.character(joined$platform), as.character(blind$platform))) {
    stop(cohort, ": platform differs between manifest and blind RDS", call. = FALSE)
  }
  if (!identical(as.character(joined$batch), as.character(blind$technical_batch))) {
    stop(cohort, ": technical batch differs between manifest and blind RDS", call. = FALSE)
  }
  if (any(!joined$status %in% c("healthy_control", "RSV_case"))) {
    stop(cohort, ": unexpected frozen case label", call. = FALSE)
  }
  case_status <- ifelse(joined$status == "RSV_case", "RSV", "healthy")
  model_data <- data.frame(
    sample_id = modeled_ids,
    subject_id = as.character(joined$subject_id),
    cohort = cohort,
    platform = as.character(joined$platform),
    case_status = case_status,
    age_months = suppressWarnings(as.numeric(joined$age_months)),
    sex = as.character(joined$sex),
    technical_batch = as.character(blind$technical_batch),
    stringsAsFactors = FALSE
  )
  list(
    expression = expression,
    sample_data = model_data,
    rds_path = rds_path,
    n_manifest_eligible = nrow(eligible),
    n_qc_excluded = length(excluded_ids),
    qc_excluded_ids = excluded_ids,
    preprocessing = object$preprocessing
  )
}

main <- function() {
  frozen <- verify_frozen_inputs()
  manifest_path <- file.path("data", "clean", "sample_manifest_frozen.tsv")
  manifest <- read_tsv(manifest_path)
  if (!identical(sort(unique(manifest$series_accession[manifest$meta_eligible_confirmatory])), sort(task7_cohorts))) {
    stop("Confirmatory manifest scope differs from the two frozen cohorts", call. = FALSE)
  }
  outputs <- character()
  summaries <- list()
  for (cohort in task7_cohorts) {
    input <- load_cohort_inputs(cohort, manifest)
    fit <- fit_task7_cohort(input$expression, input$sample_data, cohort, task7_direction)
    effects_path <- file.path("results", "cohort", paste0(cohort, "_gene_effects.tsv"))
    diagnostics_path <- file.path("results", "cohort", paste0(cohort, "_model_diagnostics.tsv"))
    rds_sha <- frozen$hashes[[input$rds_path]]
    diagnostics <- data.frame(
      cohort = cohort,
      direction = task7_direction,
      coefficient = fit$coefficient,
      model_formula = fit$formula,
      n_manifest_eligible = input$n_manifest_eligible,
      n_qc_excluded = input$n_qc_excluded,
      n_rds_samples = ncol(input$expression),
      n_complete_cases = fit$n_complete,
      n_excluded_missing_covariate = fit$n_excluded_missing_covariate,
      n_missing_age = fit$n_missing_age,
      n_missing_sex = fit$n_missing_sex,
      n_healthy = fit$n_healthy,
      n_RSV = fit$n_RSV,
      technical_batch_included = fit$batch_included,
      technical_batch_reason = fit$batch_reason,
      technical_batch_n_levels = length(fit$batch_levels),
      technical_batch_levels = paste(fit$batch_levels, collapse = "|"),
      design_rank = fit$design_rank,
      design_columns = fit$design_columns,
      design_full_rank = fit$design_rank == fit$design_columns,
      genes_tested_unfiltered = nrow(fit$gene_effects),
      manifest_sha256 = frozen$hashes[[manifest_path]],
      expression_rds_sha256 = rds_sha,
      statistical_pipeline_sha256 = frozen$hashes[[pipeline_path]],
      limma_version = as.character(utils::packageVersion("limma")),
      severity_in_primary_model = FALSE,
      stringsAsFactors = FALSE
    )
    write_tsv(fit$gene_effects, effects_path)
    write_tsv(diagnostics, diagnostics_path)
    outputs <- c(outputs, effects_path, diagnostics_path)
    summaries[[cohort]] <- diagnostics
  }

  session_path <- file.path("logs", "session_info", "task7_session_info.txt")
  commands_path <- file.path("logs", "session_info", "task7_commands.log")
  dir.create(dirname(session_path), recursive = TRUE, showWarnings = FALSE)
  session_lines <- sub("[[:space:]]+$", "", capture.output(sessionInfo()))
  writeLines(session_lines, session_path, useBytes = TRUE)
  writeLines(c(
    "Task 7: confirmatory cohort-level limma models only",
    "Command: RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_CONFIG_NAMESPACES_CHECK=FALSE Rscript analysis/R/05_fit_cohort_models.R",
    paste0("Direction: ", task7_direction),
    "Cohorts: GSE105450 and GSE103842 only; no camera, target-set test, meta-analysis, pathway analysis, or exploratory cohort",
    "Labels: frozen phenotype manifest joined one-to-one by sample_id only after blind QC RDS load",
    "Covariates: continuous age_months and sex; complete cases without imputation",
    "Technical batch: included only when available and selected design is full rank",
    "Severity: excluded from primary case-control model",
    paste0("Task5 checksum manifest SHA-256: ", frozen$task5_checksum_sha256),
    paste0("Task5 supersession SHA-256: ", frozen$task5_supersession_sha256),
    paste0("Task6 checksum manifest SHA-256: ", frozen$task6_checksum_sha256),
    vapply(summaries, function(x) paste0(
      x$cohort, ": n=", x$n_complete_cases, ", healthy=", x$n_healthy,
      ", RSV=", x$n_RSV, ", formula=", x$model_formula,
      ", rank=", x$design_rank, "/", x$design_columns,
      ", genes=", x$genes_tested_unfiltered
    ), character(1))
  ), commands_path, useBytes = TRUE)
  outputs <- c(outputs, session_path, commands_path)

  script_path <- file.path("analysis", "R", "05_fit_cohort_models.R")
  test_path <- file.path("analysis", "tests", "test_fit_cohort_models.R")
  output_test_path <- file.path("analysis", "tests", "test_fit_cohort_model_outputs.R")
  checksum_artifacts <- c(
    script_path, test_path, output_test_path, pipeline_path,
    manifest_path,
    file.path("data", "derived", paste0(task7_cohorts, "_expression.rds")),
    file.path("logs", "checksums", "task5_sha256.tsv"),
    file.path("logs", "checksums", "task6_sha256.tsv"),
    file.path("logs", "session_info", "task5_supersession.md"),
    outputs
  )
  roles <- c(
    rep("input", 10L),
    rep("output", length(outputs))
  )
  checksum <- data.frame(
    artifact = checksum_artifacts,
    role = roles,
    bytes = as.numeric(file.info(checksum_artifacts)$size),
    sha256 = vapply(checksum_artifacts, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_tsv(checksum, file.path("logs", "checksums", "task7_sha256.tsv"))

  combined <- do.call(rbind, summaries)
  rownames(combined) <- NULL
  print(combined[, c(
    "cohort", "n_complete_cases", "n_excluded_missing_covariate", "n_healthy", "n_RSV",
    "technical_batch_included", "technical_batch_n_levels", "design_rank", "design_columns",
    "genes_tested_unfiltered"
  )], row.names = FALSE)
  invisible(combined)
}

if (!identical(Sys.getenv("TASK7_SKIP_MAIN"), "1")) main()
