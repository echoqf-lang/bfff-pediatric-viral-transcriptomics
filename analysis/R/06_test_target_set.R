#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

task7_script <- file.path("analysis", "R", "05_fit_cohort_models.R")
if (!file.exists(task7_script)) stop("Frozen Task 7 model script is absent", call. = FALSE)
old_task7_skip <- Sys.getenv("TASK7_SKIP_MAIN", unset = NA_character_)
Sys.setenv(TASK7_SKIP_MAIN = "1")
source(task7_script, local = FALSE)
if (is.na(old_task7_skip)) Sys.unsetenv("TASK7_SKIP_MAIN") else Sys.setenv(TASK7_SKIP_MAIN = old_task7_skip)

task8_cohorts <- c("GSE105450", "GSE103842")
task8_independent_n <- c(GSE105450 = 122L, GSE103842 = 73L)
task8_seed <- 20260804L
task8_nrot <- 99999L

assert_task8_packages <- function() {
  required <- c("limma", "fgsea")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing frozen package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

hash_table <- function(x) {
  path <- tempfile("task8_hash_", fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    na = "NA", fileEncoding = "UTF-8"
  )
  sha256_file(path)
}

assert_task7_artifact <- function(checksums, artifact, expected_role) {
  row <- checksums[checksums$artifact == artifact, , drop = FALSE]
  if (nrow(row) != 1L || !identical(as.character(row$role), expected_role)) {
    stop("Task 7 checksum contract failed for ", artifact, call. = FALSE)
  }
  if (!file.exists(artifact)) stop("Task 7 artifact is absent: ", artifact, call. = FALSE)
  observed <- sha256_file(artifact)
  if (!identical(observed, as.character(row$sha256))) {
    stop("Task 7 artifact SHA-256 mismatch: ", artifact, call. = FALSE)
  }
  observed
}

read_frozen_targets <- function(path = file.path("data", "clean", "bfff_targets_primary.tsv")) {
  targets <- read_tsv(path)
  required <- c("entrez_id", "gene_symbol")
  if (!all(required %in% names(targets))) stop("Frozen primary target schema is invalid", call. = FALSE)
  entrez <- as.character(targets$entrez_id)
  if (nrow(targets) != 520L || anyNA(entrez) || any(!nzchar(entrez)) || anyDuplicated(entrez)) {
    stop("Expected exactly 520 unique non-empty frozen primary Entrez targets", call. = FALSE)
  }
  entrez
}

verify_task8_inputs <- function() {
  assert_task8_packages()
  task7_checksum_path <- file.path("logs", "checksums", "task7_sha256.tsv")
  task5_checksum_path <- file.path("logs", "checksums", "task5_sha256.tsv")
  if (!all(file.exists(c(task7_checksum_path, task5_checksum_path)))) {
    stop("Frozen Task 5/7 checksum manifests are incomplete", call. = FALSE)
  }
  task7 <- read_tsv(task7_checksum_path)
  task5 <- read_tsv(task5_checksum_path)
  target_path <- file.path("data", "clean", "bfff_targets_primary.tsv")
  target_sha <- assert_checksum_row(task5, target_path, "input")
  paths <- c(
    task7_script,
    file.path("data", "clean", "sample_manifest_frozen.tsv"),
    file.path("data", "derived", paste0(task8_cohorts, "_expression.rds")),
    file.path("results", "cohort", paste0(task8_cohorts, "_gene_effects.tsv")),
    file.path("results", "cohort", paste0(task8_cohorts, "_model_diagnostics.tsv"))
  )
  roles <- c("input", rep("input", 3L), rep("output", 4L))
  hashes <- setNames(vapply(
    seq_along(paths), function(i) assert_task7_artifact(task7, paths[[i]], roles[[i]]), character(1)
  ), paths)
  list(
    hashes = hashes,
    target_sha256 = target_sha,
    task7_checksum_sha256 = sha256_file(task7_checksum_path),
    task5_checksum_sha256 = sha256_file(task5_checksum_path)
  )
}

same_gene_effects <- function(saved, reconstructed, tolerance = 1e-13) {
  numeric_fields <- c("log2FC", "SE", "moderated_t", "p_value", "fdr_bh", "average_expression")
  required <- c("gene_id", "cohort", "direction", numeric_fields)
  if (!identical(names(saved), names(reconstructed)) || !all(required %in% names(saved)) ||
      nrow(saved) != nrow(reconstructed)) {
    return(FALSE)
  }
  identities_match <-
    identical(as.character(saved$gene_id), as.character(reconstructed$gene_id)) &&
    identical(as.character(saved$cohort), as.character(reconstructed$cohort)) &&
    identical(as.character(saved$direction), as.character(reconstructed$direction))
  numerics_match <- all(vapply(numeric_fields, function(field) {
    isTRUE(all.equal(saved[[field]], reconstructed[[field]], tolerance = tolerance))
  }, logical(1)))
  identities_match && numerics_match
}

verify_reconstructed_model <- function(cohort, manifest, frozen) {
  input <- load_cohort_inputs(cohort, manifest)
  prepared <- prepare_model_data(input$expression, input$sample_data, cohort)
  fit <- fit_task7_cohort(input$expression, input$sample_data, cohort, task7_direction)
  effects_path <- file.path("results", "cohort", paste0(cohort, "_gene_effects.tsv"))
  diagnostics_path <- file.path("results", "cohort", paste0(cohort, "_model_diagnostics.tsv"))
  saved_effects <- read_tsv(effects_path)
  diagnostics <- read_tsv(diagnostics_path)
  if (nrow(diagnostics) != 1L || !identical(as.character(diagnostics$coefficient), fit$coefficient)) {
    stop(cohort, ": Task 7 coefficient does not reproduce", call. = FALSE)
  }
  if (!identical(as.character(diagnostics$model_formula), fit$formula) ||
      diagnostics$design_rank != qr(prepared$design)$rank ||
      diagnostics$design_columns != ncol(prepared$design)) {
    stop(cohort, ": Task 7 design does not reproduce", call. = FALSE)
  }
  if (!same_gene_effects(saved_effects, fit$gene_effects)) {
    stop(cohort, ": Task 7 gene effects do not reproduce", call. = FALSE)
  }
  reconstructed_effect_sha <- hash_table(fit$gene_effects)
  expected_effect_sha <- frozen$hashes[[effects_path]]
  if (!identical(reconstructed_effect_sha, expected_effect_sha)) {
    stop(cohort, ": reconstructed gene-effect byte hash differs from frozen Task 7 output", call. = FALSE)
  }
  if (nrow(prepared$sample_data) != unname(task8_independent_n[[cohort]])) {
    stop(cohort, ": complete-case model N differs from the frozen Stouffer weight N", call. = FALSE)
  }
  list(
    expression = prepared$expression,
    sample_data = prepared$sample_data,
    design = prepared$design,
    coefficient = prepared$coefficient,
    gene_effects = fit$gene_effects,
    design_sha256 = hash_table(data.frame(sample_id = rownames(prepared$design), prepared$design, check.names = FALSE)),
    task7_effects_sha256 = expected_effect_sha,
    reconstructed_effects_sha256 = reconstructed_effect_sha,
    model_formula = fit$formula,
    n_model = nrow(prepared$sample_data)
  )
}

run_camera_prepared <- function(expression, design, coefficient, frozen_targets, cohort_id) {
  if (!is.matrix(expression) || !is.numeric(expression) || any(!is.finite(expression))) {
    stop("camera expression must be a finite numeric matrix", call. = FALSE)
  }
  if (!is.matrix(design) || !is.numeric(design) || nrow(design) != ncol(expression) ||
      qr(design)$rank != ncol(design)) {
    stop("camera design must be a conformable full-rank numeric matrix", call. = FALSE)
  }
  if (!coefficient %in% colnames(design)) stop("camera coefficient is absent", call. = FALSE)
  frozen_targets <- as.character(frozen_targets)
  if (!length(frozen_targets) || anyNA(frozen_targets) || any(!nzchar(frozen_targets)) || anyDuplicated(frozen_targets)) {
    stop("Frozen target IDs must be unique and non-empty", call. = FALSE)
  }
  index <- which(rownames(expression) %in% frozen_targets)
  n_frozen <- length(frozen_targets)
  n_detectable <- length(index)
  coverage <- n_detectable / n_frozen
  if (n_detectable < 10L || coverage < 0.50) {
    reason <- c(if (n_detectable < 10L) "n_detectable_below_10", if (coverage < 0.50) "coverage_below_0.50")
    return(data.frame(
      cohort = cohort_id, contrast = task7_direction, status = "uninterpretable",
      interpretability_reason = paste(reason, collapse = ";"), n_frozen = n_frozen,
      n_detectable = n_detectable, coverage = coverage, n_genes = n_detectable,
      inter_gene_correlation = NA_real_, direction = NA_character_, p_value = NA_real_,
      signed_z = NA_real_, stringsAsFactors = FALSE
    ))
  }
  result <- run_camera_analysis(
    y = expression,
    index = list(frozen_primary_target_set = index),
    design = design,
    contrast = coefficient,
    inter.gene.cor = NA_real_,
    sort = FALSE,
    directional = TRUE
  )
  required <- c("NGenes", "Correlation", "Direction", "PValue")
  if (nrow(result) != 1L || !all(required %in% names(result))) {
    stop("Unexpected limma::camera API result", call. = FALSE)
  }
  p <- result$PValue[[1L]]
  direction <- as.character(result$Direction[[1L]])
  correlation <- result$Correlation[[1L]]
  if (!direction %in% c("Up", "Down") || !is.finite(p) || p <= 0 || p > 1 ||
      !is.finite(correlation)) {
    stop("Invalid limma::camera inferential result", call. = FALSE)
  }
  signed_z <- if (direction == "Up") qnorm(p / 2, lower.tail = FALSE) else -qnorm(p / 2, lower.tail = FALSE)
  data.frame(
    cohort = cohort_id, contrast = task7_direction, status = "interpretable",
    interpretability_reason = "passed_n_detectable_and_coverage_gate",
    n_frozen = n_frozen, n_detectable = n_detectable, coverage = coverage,
    n_genes = result$NGenes[[1L]], inter_gene_correlation = correlation,
    direction = direction, p_value = p, signed_z = signed_z,
    stringsAsFactors = FALSE
  )
}

extract_roast_directions <- function(roast_result, cohort_id, n_frozen, n_detectable, nrot) {
  pmat <- roast_result$p.value
  required_rows <- c("Down", "Up", "Mixed")
  if (!(is.matrix(pmat) || is.data.frame(pmat)) || !all(required_rows %in% rownames(pmat)) ||
      !all(c("Active.Prop", "P.Value") %in% colnames(pmat))) {
    stop("Unexpected limma::roast API result", call. = FALSE)
  }
  p <- as.numeric(pmat[required_rows, "P.Value"])
  active <- as.numeric(pmat[required_rows, "Active.Prop"])
  if (any(!is.finite(p) | p < 0 | p > 1) || any(!is.finite(active) | active < 0 | active > 1)) {
    stop("Invalid limma::roast result", call. = FALSE)
  }
  data.frame(
    cohort = cohort_id,
    contrast = task7_direction,
    test_direction = required_rows,
    n_frozen = as.integer(n_frozen),
    n_detectable = as.integer(n_detectable),
    coverage = n_detectable / n_frozen,
    active_proportion = active,
    p_value = p,
    p_holm = p.adjust(p, method = "holm"),
    nrot = as.integer(nrot),
    seed = task8_seed,
    role = "key_supplementary_self_contained",
    stringsAsFactors = FALSE
  )
}

run_roast_prepared <- function(expression, design, coefficient, frozen_targets, cohort_id, nrot = task8_nrot) {
  index <- which(rownames(expression) %in% frozen_targets)
  if (length(index) < 10L || length(index) / length(frozen_targets) < 0.50) {
    stop(cohort_id, ": roast not run because target set is uninterpretable", call. = FALSE)
  }
  result <- limma::roast(
    y = expression, index = index, design = design, contrast = coefficient,
    set.statistic = "mean", nrot = nrot
  )
  extract_roast_directions(result, cohort_id, length(frozen_targets), length(index), nrot)
}

run_fgsea_exploratory <- function(gene_effects, frozen_targets, cohort_id) {
  stats <- gene_effects$moderated_t
  names(stats) <- as.character(gene_effects$gene_id)
  if (any(!is.finite(stats)) || anyDuplicated(names(stats)) || any(!nzchar(names(stats)))) {
    stop(cohort_id, ": invalid moderated-t ranking", call. = FALSE)
  }
  pathway <- intersect(as.character(frozen_targets), names(stats))
  result <- fgsea::fgsea(
    pathways = list(frozen_primary_target_set = pathway),
    stats = stats,
    minSize = 1L,
    maxSize = length(stats) - 1L,
    eps = 0
  )
  required <- c("pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge")
  if (nrow(result) != 1L || !all(required %in% names(result))) {
    stop("Unexpected fgsea API result", call. = FALSE)
  }
  data.frame(
    cohort = cohort_id,
    contrast = task7_direction,
    ranking_statistic = "limma_moderated_t",
    pathway = as.character(result$pathway[[1L]]),
    n_frozen = length(frozen_targets),
    n_detectable = length(pathway),
    coverage = length(pathway) / length(frozen_targets),
    ES = result$ES[[1L]],
    NES = result$NES[[1L]],
    p_value = result$pval[[1L]],
    fdr_bh = result$padj[[1L]],
    leading_edge_entrez = paste(as.character(result$leadingEdge[[1L]]), collapse = "|"),
    seed = task8_seed,
    role = "exploratory_ranked_sensitivity_not_independent_validation",
    stringsAsFactors = FALSE
  )
}

combine_camera_stouffer <- function(camera) {
  if (!identical(as.character(camera$cohort), task8_cohorts)) {
    stop("Stouffer input must contain the two confirmatory cohorts in frozen order", call. = FALSE)
  }
  interpretable <- all(camera$status == "interpretable")
  if (!interpretable) {
    return(list(
      both_interpretable = FALSE, direction_consistent = FALSE,
      combined_z = NA_real_, combined_p_value = NA_real_,
      camera_replication_component_met = FALSE
    ))
  }
  expected_z <- ifelse(
    camera$direction == "Up",
    qnorm(camera$p_value / 2, lower.tail = FALSE),
    -qnorm(camera$p_value / 2, lower.tail = FALSE)
  )
  if (!isTRUE(all.equal(camera$signed_z, expected_z, tolerance = 1e-14))) {
    stop("camera signed Z does not match the two-sided directional definition", call. = FALSE)
  }
  weights <- sqrt(unname(task8_independent_n[task8_cohorts]))
  combined_z <- sum(weights * camera$signed_z) / sqrt(sum(weights^2))
  combined_p <- 2 * pnorm(-abs(combined_z))
  direction_consistent <- length(unique(camera$direction)) == 1L
  component <- all(camera$p_value < 0.05) && direction_consistent && combined_p < 0.05
  list(
    both_interpretable = TRUE,
    direction_consistent = direction_consistent,
    combined_z = combined_z,
    combined_p_value = combined_p,
    camera_replication_component_met = component
  )
}

make_partial_decision <- function(camera, stouffer) {
  camera_failed <- !isTRUE(stouffer$camera_replication_component_met)
  data.frame(
    analysis_scope = "GSE105450_and_GSE103842_confirmatory_only",
    both_cohorts_interpretable = stouffer$both_interpretable,
    GSE105450_camera_p_lt_0_05 = isTRUE(camera$p_value[camera$cohort == "GSE105450"] < 0.05),
    GSE103842_camera_p_lt_0_05 = isTRUE(camera$p_value[camera$cohort == "GSE103842"] < 0.05),
    camera_direction_consistent = stouffer$direction_consistent,
    weighted_stouffer_z = stouffer$combined_z,
    weighted_stouffer_two_sided_p = stouffer$combined_p_value,
    weighted_stouffer_p_lt_0_05 = isTRUE(stouffer$combined_p_value < 0.05),
    weight_GSE105450_sqrt_n = sqrt(task8_independent_n[["GSE105450"]]),
    weight_GSE103842_sqrt_n = sqrt(task8_independent_n[["GSE103842"]]),
    random_set_empirical_p_value = NA_real_,
    random_set_empirical_p_lt_0_05 = NA,
    camera_replication_component_met = stouffer$camera_replication_component_met,
    random_set_calibration_status = if (camera_failed) {
      "pending_scheduled_nonrescuing"
    } else {
      "pending_required_for_final_support_decision"
    },
    primary_decision = if (camera_failed) {
      "H1_not_supported_camera_component_failed"
    } else {
      "pending_random_set_calibration"
    },
    h1_supported = if (camera_failed) FALSE else NA,
    pending_requirement = "Task10_matched_random_gene_set_empirical_p_value",
    stringsAsFactors = FALSE
  )
}

main_task8 <- function() {
  frozen <- verify_task8_inputs()
  targets <- read_frozen_targets()
  manifest <- read_tsv(file.path("data", "clean", "sample_manifest_frozen.tsv"))
  if (!identical(sort(unique(manifest$series_accession[manifest$meta_eligible_confirmatory])), sort(task8_cohorts))) {
    stop("Task 8 scope differs from the two frozen confirmatory cohorts", call. = FALSE)
  }
  models <- lapply(task8_cohorts, verify_reconstructed_model, manifest = manifest, frozen = frozen)
  names(models) <- task8_cohorts

  camera_rows <- lapply(task8_cohorts, function(cohort) {
    model <- models[[cohort]]
    out <- run_camera_prepared(model$expression, model$design, model$coefficient, targets, cohort)
    out$model_coefficient <- model$coefficient
    out$model_formula <- model$model_formula
    out$n_independent_model <- model$n_model
    out$stouffer_weight_sqrt_n <- sqrt(model$n_model)
    out$design_sha256 <- model$design_sha256
    out$task7_effects_sha256 <- model$task7_effects_sha256
    out$reconstructed_effects_sha256 <- model$reconstructed_effects_sha256
    out$task7_effect_hash_match <- identical(model$task7_effects_sha256, model$reconstructed_effects_sha256)
    out$camera_p_sidedness <- "two_sided"
    out$background <- "all_detectable_unique_Entrez_genes_in_cohort"
    out
  })
  camera <- do.call(rbind, camera_rows)
  rownames(camera) <- NULL
  stouffer <- combine_camera_stouffer(camera)
  decision <- make_partial_decision(camera, stouffer)

  set.seed(task8_seed)
  roast <- do.call(rbind, lapply(task8_cohorts, function(cohort) {
    model <- models[[cohort]]
    run_roast_prepared(model$expression, model$design, model$coefficient, targets, cohort, task8_nrot)
  }))
  rownames(roast) <- NULL

  set.seed(task8_seed)
  fgsea <- do.call(rbind, lapply(task8_cohorts, function(cohort) {
    run_fgsea_exploratory(models[[cohort]]$gene_effects, targets, cohort)
  }))
  rownames(fgsea) <- NULL

  output_paths <- c(
    camera = file.path("results", "cohort", "target_set_camera.tsv"),
    roast = file.path("results", "cohort", "target_set_roast.tsv"),
    fgsea = file.path("results", "cohort", "target_set_fgsea_exploratory.tsv"),
    decision = file.path("results", "cohort", "target_set_primary_decision_partial.tsv")
  )
  write_tsv(camera, output_paths[["camera"]])
  write_tsv(roast, output_paths[["roast"]])
  write_tsv(fgsea, output_paths[["fgsea"]])
  write_tsv(decision, output_paths[["decision"]])

  session_path <- file.path("logs", "session_info", "task8_session_info.txt")
  commands_path <- file.path("logs", "session_info", "task8_commands.log")
  dir.create(dirname(session_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())), session_path, useBytes = TRUE)
  writeLines(c(
    "Task 8: frozen primary target-set testing",
    "Command: RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_CONFIG_NAMESPACES_CHECK=FALSE Rscript analysis/R/06_test_target_set.R",
    paste0("Confirmatory cohorts only: ", paste(task8_cohorts, collapse = ", ")),
    "Excluded from Task 8: GSE38900, GSE188427, random-set calibration, gene-level meta-analysis, sensitivity target set",
    "Expression input: frozen processed-only scaled-MAD Task 5 RDS; no RAW/CEL/FASTQ/SRA input or download",
    "Model: exact Task 7 formula, design, coefficient, complete-case sample set, and batch rule reconstructed and byte-hash checked",
    paste0("Frozen primary targets: ", length(targets), " unique Entrez IDs"),
    paste0("camera: competitive, directional=TRUE, reported P is two-sided, estimated inter-gene correlation"),
    paste0("roast: set.statistic=mean, nrot=", task8_nrot, ", seed=", task8_seed, ", within-cohort Holm across Down/Up/Mixed"),
    paste0("fgsea: moderated-t ranking, seed=", task8_seed, ", exploratory only; leading edge is not independent validation"),
    paste0("Weighted Stouffer N: GSE105450=", task8_independent_n[["GSE105450"]], ", GSE103842=", task8_independent_n[["GSE103842"]]),
    paste0(
      "Primary decision: ", decision$primary_decision,
      "; Task 10 empirical random-set P remains scheduled but cannot rescue a failed camera component"
    ),
    paste0("Task 5 checksum manifest SHA-256: ", frozen$task5_checksum_sha256),
    paste0("Task 7 checksum manifest SHA-256: ", frozen$task7_checksum_sha256),
    paste0("Frozen target SHA-256: ", frozen$target_sha256)
  ), commands_path, useBytes = TRUE)

  checksum_inputs <- c(
    file.path("analysis", "R", "06_test_target_set.R"),
    file.path("analysis", "tests", "test_target_set_analysis.R"),
    file.path("analysis", "tests", "test_target_set_outputs.R"),
    task7_script,
    file.path("data", "clean", "bfff_targets_primary.tsv"),
    file.path("data", "clean", "sample_manifest_frozen.tsv"),
    file.path("data", "derived", paste0(task8_cohorts, "_expression.rds")),
    file.path("results", "cohort", paste0(task8_cohorts, "_gene_effects.tsv")),
    file.path("results", "cohort", paste0(task8_cohorts, "_model_diagnostics.tsv")),
    file.path("logs", "checksums", "task5_sha256.tsv"),
    file.path("logs", "checksums", "task7_sha256.tsv")
  )
  checksum_outputs <- c(unname(output_paths), session_path, commands_path)
  checksum_artifacts <- c(checksum_inputs, checksum_outputs)
  checksum <- data.frame(
    artifact = checksum_artifacts,
    role = c(rep("input", length(checksum_inputs)), rep("output", length(checksum_outputs))),
    bytes = as.numeric(file.info(checksum_artifacts)$size),
    sha256 = vapply(checksum_artifacts, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_tsv(checksum, file.path("logs", "checksums", "task8_sha256.tsv"))

  print(camera[, c("cohort", "n_frozen", "n_detectable", "coverage", "direction", "p_value", "inter_gene_correlation")], row.names = FALSE)
  print(decision, row.names = FALSE)
  invisible(list(camera = camera, roast = roast, fgsea = fgsea, decision = decision))
}

if (!identical(Sys.getenv("TASK8_SKIP_MAIN"), "1")) main_task8()
