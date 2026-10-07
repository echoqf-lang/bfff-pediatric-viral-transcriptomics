#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
source(file.path("analysis", "R", "single_gene_upgrade_pipeline.R"))

load_gse77087_functions <- function() {
  old <- Sys.getenv("GSE77087_REVISED_SKIP_MAIN", unset = NA_character_)
  Sys.setenv(GSE77087_REVISED_SKIP_MAIN = "1")
  source(file.path("analysis", "R", "12_gse77087_revised_main.R"), local = FALSE)
  if (is.na(old)) Sys.unsetenv("GSE77087_REVISED_SKIP_MAIN") else Sys.setenv(GSE77087_REVISED_SKIP_MAIN = old)
}

extract_limma_coefficient <- function(fit, coefficient, gene_ids, family_n = 86L) {
  idx <- match(gene_ids, rownames(fit$coefficients))
  coef_idx <- match(coefficient, colnames(fit$coefficients))
  se <- fit$stdev.unscaled[, coef_idx] * sqrt(fit$s2.post)
  out <- data.frame(
    gene_id = gene_ids,
    log2FC = fit$coefficients[idx, coef_idx],
    SE = se[idx],
    moderated_t = fit$t[idx, coef_idx],
    p_value = fit$p.value[idx, coef_idx],
    stringsAsFactors = FALSE
  )
  out$fdr_bh <- bh_family(out$p_value)
  if (nrow(out) != family_n) stop("Severity output must preserve all 86 frozen genes", call. = FALSE)
  out
}

main_gse77087_severity <- function() {
  if (!requireNamespace("limma", quietly = TRUE)) stop("limma package is required", call. = FALSE)
  load_gse77087_functions()
  prepared <- prepare_gse77087()
  frozen <- read_upgrade_tsv(file.path("data", "clean", "bfff_86_gene_set_frozen.tsv"), col_classes = "character")
  expression <- prepared$expression
  sample_data <- prepared$sample_data
  clinical <- as.character(sample_data$clinical_group)
  if (!identical(c(sum(clinical == "healthy"), sum(clinical == "outpatient"), sum(clinical == "inpatient")), c(23L, 20L, 61L))) {
    stop("GSE77087 severity group counts changed", call. = FALSE)
  }
  model_data <- data.frame(
    severity_ordinal = c(healthy = 0, outpatient = 1, inpatient = 2)[clinical],
    severity_group = factor(clinical, levels = c("healthy", "outpatient", "inpatient")),
    age_months = sample_data$age_months,
    sex = factor(sample_data$sex),
    technical_batch = factor(sample_data$technical_batch),
    stringsAsFactors = FALSE
  )
  design_trend <- stats::model.matrix(~ severity_ordinal + age_months + sex + technical_batch, data = model_data)
  design_group <- stats::model.matrix(~ severity_group + age_months + sex + technical_batch, data = model_data)
  if (qr(design_trend)$rank != ncol(design_trend) || qr(design_group)$rank != ncol(design_group)) stop("Severity design is rank deficient", call. = FALSE)

  fit_trend <- limma::eBayes(limma::lmFit(expression, design_trend))
  trend <- extract_limma_coefficient(fit_trend, "severity_ordinal", frozen$entrez_id)
  trend$gene_symbol <- frozen$gene_symbol
  trend$contrast <- "one_level_increase_healthy_outpatient_hospitalized"

  fit_group <- limma::eBayes(limma::lmFit(expression, design_group))
  outpatient <- extract_limma_coefficient(fit_group, "severity_groupoutpatient", frozen$entrez_id)
  outpatient$gene_symbol <- frozen$gene_symbol
  outpatient$contrast <- "outpatient_minus_healthy"
  inpatient <- extract_limma_coefficient(fit_group, "severity_groupinpatient", frozen$entrez_id)
  inpatient$gene_symbol <- frozen$gene_symbol
  inpatient$contrast <- "hospitalized_minus_healthy"
  direct_matrix <- limma::makeContrasts(
    hospitalized_minus_outpatient = severity_groupinpatient - severity_groupoutpatient,
    levels = design_group
  )
  fit_direct <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expression, design_group), direct_matrix))
  direct <- extract_limma_coefficient(fit_direct, "hospitalized_minus_outpatient", frozen$entrez_id)
  direct$gene_symbol <- frozen$gene_symbol
  direct$contrast <- "hospitalized_minus_outpatient"
  contrasts <- rbind(outpatient, inpatient, direct)

  root <- file.path("results", "single_gene_upgrade_v1", "severity")
  write_upgrade_tsv(trend, file.path(root, "GSE77087_86_gene_ordinal_trend.tsv"))
  write_upgrade_tsv(contrasts, file.path(root, "GSE77087_86_gene_group_contrasts.tsv"))
  diagnostics <- data.frame(
    n_total = nrow(sample_data), n_healthy = sum(clinical == "healthy"),
    n_outpatient = sum(clinical == "outpatient"), n_hospitalized = sum(clinical == "inpatient"),
    trend_design_rank = qr(design_trend)$rank, trend_design_columns = ncol(design_trend),
    group_design_rank = qr(design_group)$rank, group_design_columns = ncol(design_group),
    mapped_genes = sum(complete.cases(trend[c("log2FC", "SE")])),
    trend_bh_lt_0_05_n = sum(trend$fdr_bh < 0.05, na.rm = TRUE),
    hospitalized_vs_outpatient_bh_lt_0_05_n = sum(direct$fdr_bh < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  write_upgrade_tsv(diagnostics, file.path(root, "severity_model_diagnostics.tsv"))
  message("SEVERITY n=104 trend_BH<0.05=", diagnostics$trend_bh_lt_0_05_n)
  invisible(list(trend = trend, contrasts = contrasts, diagnostics = diagnostics))
}

if (sys.nframe() == 0L) main_gse77087_severity()
