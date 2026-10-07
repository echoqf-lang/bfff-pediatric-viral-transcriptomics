#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
source(file.path("analysis", "R", "single_gene_upgrade_pipeline.R"))

read_cell_config <- function(path = file.path("analysis", "config", "cell_signature_source.yml")) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml package is required", call. = FALSE)
  yaml::read_yaml(path)
}

estimate_nnls_fractions <- function(signature, mixtures) {
  signature <- as.matrix(signature)
  mixtures <- as.matrix(mixtures)
  if (!identical(rownames(signature), rownames(mixtures))) stop("Signature and mixture genes must be aligned", call. = FALSE)
  out <- apply(mixtures, 2L, function(y) {
    fit <- stats::optim(rep(1 / ncol(signature), ncol(signature)),
                        fn = function(p) sum((y - as.numeric(signature %*% p))^2),
                        method = "L-BFGS-B", lower = rep(0, ncol(signature)))
    p <- fit$par
    if (!is.finite(sum(p)) || sum(p) <= 0) return(rep(NA_real_, length(p)))
    p / sum(p)
  })
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  rownames(out) <- colnames(signature)
  colnames(out) <- colnames(mixtures)
  out
}

validate_fraction_estimator <- function(seed = 20260806L) {
  set.seed(seed)
  signature <- matrix(c(10, 1, 8, 1, 6, 1, 1, 10, 1, 8, 1, 6), ncol = 2L)
  rownames(signature) <- paste0("marker", seq_len(nrow(signature)))
  colnames(signature) <- c("cell_A", "cell_B")
  true_a <- seq(0.05, 0.95, length.out = 30L)
  truth <- rbind(cell_A = true_a, cell_B = 1 - true_a)
  mixtures <- signature %*% truth + matrix(stats::rnorm(nrow(signature) * ncol(truth), sd = 0.05), nrow = nrow(signature))
  rownames(mixtures) <- rownames(signature)
  estimates <- estimate_nnls_fractions(signature, mixtures)
  rho <- stats::cor(as.numeric(truth), as.numeric(estimates), method = "spearman", use = "complete.obs")
  list(spearman = rho, truth = truth, estimates = estimates)
}

build_stopped_qc <- function(config, simulation_rho) {
  cohorts <- c("GSE38900", "GSE77087", "GSE103842")
  data.frame(
    cohort = cohorts,
    analysis_status = "stopped_by_qc",
    signature_available = FALSE,
    signature_population_applicable = FALSE,
    marker_coverage = NA_real_,
    boundary_fraction = NA_real_,
    stable_major_cell_class = FALSE,
    simulation_spearman = simulation_rho,
    simulation_pass = simulation_rho > as.numeric(config$minimum_simulation_spearman),
    real_data_qc_pass = FALSE,
    stop_reason = config$stop_reason,
    fallback = config$fallback,
    stringsAsFactors = FALSE
  )
}

main_cell_composition_sensitivity <- function() {
  config <- read_cell_config()
  validation <- validate_fraction_estimator()
  if (!is.finite(validation$spearman) || validation$spearman <= as.numeric(config$minimum_simulation_spearman)) {
    stop("Cell-fraction estimator failed the preregistered simulation gate", call. = FALSE)
  }
  if (is.null(config$signature_file) || !nzchar(config$signature_file) || !file.exists(config$signature_file) ||
      !identical(config$status, "ready_for_blind_qc")) {
    qc <- build_stopped_qc(config, validation$spearman)
    path <- file.path("results", "single_gene_upgrade_v1", "cell_sensitivity", "cell_fraction_qc.tsv")
    write_upgrade_tsv(qc, path)
    message("CELL_SENSITIVITY status=stopped_by_qc simulation_spearman=", signif(validation$spearman, 4),
            " reason=reference_signature_unavailable_or_not_applicable")
    return(invisible(list(status = "stopped_by_qc", qc = qc)))
  }
  stop("A real-data signature was declared ready but the locked implementation has not been supplied", call. = FALSE)
}

if (sys.nframe() == 0L) main_cell_composition_sensitivity()
