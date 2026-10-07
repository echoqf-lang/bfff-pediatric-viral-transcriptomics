#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

task13_read_tsv <- function(path) read.delim(
  path, sep = "\t", quote = "", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8"
)

task13_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA", fileEncoding = "UTF-8"
  )
}

task13_sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 failed for ", path, call. = FALSE)
  sub("[[:space:]].*$", "", output[[1L]])
}

task13_open_pdf <- function(path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(
    path, width = width, height = height, useDingbats = FALSE,
    timestamp = FALSE, bg = "white", family = "sans"
  )
}

task13_compact_table <- function(x) {
  if (!length(x)) return("none")
  tab <- table(x, useNA = "ifany")
  paste(paste0(names(tab), "=", as.integer(tab)), collapse = "; ")
}

task13_age_summary <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("not available")
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE, type = 7)
  sprintf("median %.2f [IQR %.2f-%.2f]", q[[2L]], q[[1L]], q[[3L]])
}

task13_collapse_flow <- function(flow, included) {
  one <- flow[startsWith(flow$decision_reason, if (included) "included:" else "excluded:"), , drop = FALSE]
  if (!nrow(one)) return("none")
  paste(paste0(sub("^[^:]+:", "", one$decision_reason), "=", one$n_samples), collapse = "; ")
}

task13_build_table1 <- function(manifest, flow, sensitivity_diagnostics, qc_audit) {
  confirmatory_paths <- c(
    GSE105450 = file.path("results", "cohort", "GSE105450_model_diagnostics.tsv"),
    GSE103842 = file.path("results", "cohort", "GSE103842_model_diagnostics.tsv")
  )
  confirmatory <- do.call(rbind, lapply(confirmatory_paths, task13_read_tsv))
  model <- data.frame(
    cohort = c(confirmatory$cohort, sensitivity_diagnostics$cohort, "GSE38900"),
    analysis_n = c(confirmatory$n_complete_cases, sensitivity_diagnostics$n_model, NA_integer_),
    reference_n = c(confirmatory$n_healthy, sensitivity_diagnostics$n_reference, 31L),
    case_n = c(confirmatory$n_RSV, sensitivity_diagnostics$n_case, 107L),
    stringsAsFactors = FALSE
  )
  model <- model[!duplicated(model$cohort), , drop = FALSE]

  role_map <- c(
    GSE105450 = "confirmatory_primary",
    GSE103842 = "confirmatory_primary",
    GSE188427 = "exploratory_substitution_sensitivity_overlap_risk",
    GSE38900 = "exploratory_previously_observed",
    GSE103119 = "small_sample_exploratory",
    GSE155925 = "cross_virus_exploratory"
  )
  contrast_map <- c(
    GSE105450 = "RSV minus healthy",
    GSE103842 = "RSV minus healthy",
    GSE188427 = "RSV minus healthy",
    GSE38900 = "RSV minus healthy",
    GSE103119 = "RSV minus healthy",
    GSE155925 = "single RSV minus other single virus"
  )
  covariate_map <- c(
    GSE105450 = "age_months + sex + technical_batch",
    GSE103842 = "age_months + sex + technical_batch",
    GSE188427 = "none; age and sex unavailable",
    GSE38900 = "legacy exploratory analysis; not re-fitted in Tasks 2-12",
    GSE103119 = "age_months + sex",
    GSE155925 = "age_months + sex + hospital_batch + enrollment_batch"
  )
  boundary_map <- c(
    GSE105450 = "adjusted confirmatory model; outpatient and hospitalized RSV cases",
    GSE103842 = "adjusted confirmatory model; one blind-QC exclusion",
    GSE188427 = "exploratory non-H1 substitution; unadjusted model; three blind-QC exclusions; overlap risk with GSE105450",
    GSE38900 = "previously observed exploratory cohort; retained for context and not reanalysed as confirmation",
    GSE103119 = "small-sample exploratory model; not independent confirmation",
    GSE155925 = "cross-virus exploratory contrast; no healthy controls; not independent confirmation"
  )

  cohorts <- sort(unique(manifest$series_accession))
  rows <- lapply(cohorts, function(cohort) {
    all <- manifest[manifest$series_accession == cohort, , drop = FALSE]
    included <- all[all$include, , drop = FALSE]
    if (cohort == "GSE38900") {
      summary_population <- included
      summary_definition <- "included_pre-QC manifest population; no Task2-12 model"
    } else {
      excluded_ids <- qc_audit$sample_id[
        qc_audit$cohort == cohort & !is.na(qc_audit$exclude_qc) & qc_audit$exclude_qc
      ]
      summary_population <- included[!included$sample_id %in% excluded_ids, , drop = FALSE]
      summary_definition <- "analysis population used in model after applicable blind-QC exclusions"
    }
    flow_one <- flow[flow$series_accession == cohort, , drop = FALSE]
    official <- unique(flow_one$official_samples)
    if (length(official) != 1L) stop(cohort, ": official sample count is ambiguous", call. = FALSE)
    m <- model[model$cohort == cohort, , drop = FALSE]
    if (nrow(m) != 1L) stop(cohort, ": analysis counts are absent", call. = FALSE)
    if (cohort != "GSE38900" && nrow(summary_population) != m$analysis_n) {
      stop(cohort, ": summary population does not equal analyzed population", call. = FALSE)
    }
    data.frame(
      cohort = cohort,
      analysis_role = unname(role_map[[cohort]]),
      platform = paste(sort(unique(included$platform)), collapse = ";"),
      tissue = paste(sort(unique(included$tissue)), collapse = ";"),
      timepoint = paste(sort(unique(included$timepoint)), collapse = ";"),
      contrast = unname(contrast_map[[cohort]]),
      official_n = official,
      included_n = nrow(included),
      analysis_n = m$analysis_n,
      reference_n = m$reference_n,
      case_n = m$case_n,
      summary_population_n = nrow(summary_population),
      summary_population_definition = summary_definition,
      age_months_summary = task13_age_summary(summary_population$age_months),
      age_missing_n = sum(is.na(summary_population$age_months)),
      sex_summary = task13_compact_table(summary_population$sex),
      sex_missing_n = sum(is.na(summary_population$sex)),
      severity_summary = task13_compact_table(summary_population$severity),
      severity_missing_n = sum(is.na(summary_population$severity)),
      covariates = unname(covariate_map[[cohort]]),
      inclusion_summary = task13_collapse_flow(flow_one, TRUE),
      exclusion_summary = task13_collapse_flow(flow_one, FALSE),
      analysis_boundary = unname(boundary_map[[cohort]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

task13_signed_z <- function(p, direction) {
  sign <- ifelse(direction == "Down", -1, 1)
  sign * stats::qnorm(p / 2, lower.tail = FALSE)
}

task13_supplementary_signal_context <- function(roast, fgsea) {
  cohorts <- c("GSE105450", "GSE103842")
  roast_up <- roast[roast$test_direction == "Up", , drop = FALSE]
  roast_up <- roast_up[match(cohorts, roast_up$cohort), , drop = FALSE]
  fgsea <- fgsea[match(cohorts, fgsea$cohort), , drop = FALSE]
  if (anyNA(roast_up$cohort) || anyNA(fgsea$cohort) || nrow(roast_up) != 2L || nrow(fgsea) != 2L) {
    stop("Reviewed roast/fgsea context is incomplete", call. = FALSE)
  }
  sprintf(
    paste0(
      "roast self-contained Up Holm P: GSE105450=%.3g, GSE103842=%.3g; ",
      "fgsea ranked NES/FDR: GSE105450=%.3f/%.3g, GSE103842=%.3f/%.3g. ",
      "These strong positive signals are compatible with a broad RSV whole-blood response or cell-composition shifts; ",
      "they do not establish competitive target-set enrichment or formula mechanism."
    ),
    roast_up$p_holm[[1L]], roast_up$p_holm[[2L]],
    fgsea$NES[[1L]], fgsea$fdr_bh[[1L]], fgsea$NES[[2L]], fgsea$fdr_bh[[2L]]
  )
}

task13_build_table2 <- function(camera, decision, empirical, task10_timing, supplementary_signal_context) {
  camera_rows <- data.frame(
    analysis_id = paste0(camera$cohort, "_camera"),
    scope = camera$cohort,
    method = "camera competitive gene-set test",
    direction = camera$direction,
    statistic = camera$signed_z,
    statistic_scale = "signed standard-normal score reconstructed from two-sided camera P",
    ci_status = "not_estimable",
    ci_reason = "camera provides a competitive test statistic and P value, not a prespecified target-set effect estimate with CI",
    p_value = camera$p_value,
    coverage = camera$coverage,
    coverage_detail = paste0(camera$n_detectable, "/", camera$n_frozen),
    n_independent = camera$n_independent_model,
    primary_H1_supported = FALSE,
    primary_H1_rescue_allowed = FALSE,
    primary_decision = decision$primary_decision,
    evidence_timing = "primary camera outcomes were generated under the frozen confirmatory rule",
    null_hypothesis_boundary = "competitive camera null: the target set is no more differentially expressed than the detectable background",
    supplementary_signal_context = NA_character_,
    interpretation = "direction shown; cohort-level camera P did not meet the frozen P<0.05 component",
    stringsAsFactors = FALSE
  )
  combined <- data.frame(
    analysis_id = "weighted_signed_stouffer",
    scope = "GSE105450+GSE103842",
    method = "sqrt(n)-weighted signed Stouffer",
    direction = ifelse(decision$weighted_stouffer_z >= 0, "Up", "Down"),
    statistic = decision$weighted_stouffer_z,
    statistic_scale = "combined standard-normal score",
    ci_status = "not_estimable",
    ci_reason = "signed Stouffer combines test evidence; no target-set effect estimand or CI was prespecified",
    p_value = decision$weighted_stouffer_two_sided_p,
    coverage = NA_real_,
    coverage_detail = paste(camera$n_detectable, camera$n_frozen, sep = "/", collapse = "; "),
    n_independent = sum(camera$n_independent_model),
    primary_H1_supported = FALSE,
    primary_H1_rescue_allowed = FALSE,
    primary_decision = decision$primary_decision,
    evidence_timing = "primary signed Stouffer combination was generated under the frozen confirmatory rule",
    null_hypothesis_boundary = "combines signed camera evidence and therefore inherits the camera competitive testing target",
    supplementary_signal_context = NA_character_,
    interpretation = "combined P did not meet the frozen P<0.05 component",
    stringsAsFactors = FALSE
  )
  random <- data.frame(
    analysis_id = "matched_random_empirical_calibration",
    scope = "post-outcome-clarified conditional membership-randomization sensitivity; paired cohort-specific matched sets",
    method = "10,000 matched random-set conditional membership-randomization sensitivity",
    direction = "Up",
    statistic = empirical$observed_combined_stouffer_z,
    statistic_scale = "observed combined Stouffer score compared with conditional matched-set null",
    ci_status = "not_defined",
    ci_reason = "empirical matched-set calibration is a conditional rank-tail test without an effect-scale CI",
    p_value = empirical$empirical_two_sided_p,
    coverage = NA_real_,
    coverage_detail = paste0(empirical$GSE105450_set_size, "/520; ", empirical$GSE103842_set_size, "/520"),
    n_independent = sum(camera$n_independent_model),
    primary_H1_supported = FALSE,
    primary_H1_rescue_allowed = empirical$random_calibration_can_rescue,
    primary_decision = empirical$primary_decision,
    evidence_timing = task10_timing,
    null_hypothesis_boundary = "conditional null over matched non-target gene membership; different from the camera competitive null",
    supplementary_signal_context = supplementary_signal_context,
    interpretation = "empirical P=0.00010 reflects a narrow conditional random-null distribution; it cannot validate camera, serve as independent replication, or rescue the failed primary H1 rule",
    stringsAsFactors = FALSE
  )
  out <- rbind(camera_rows, combined, random)
  rownames(out) <- NULL
  out
}

task13_build_table3 <- function(hksj, wald, fixed) {
  wald <- wald[match(hksj$gene_id, wald$gene_id), , drop = FALSE]
  fixed <- fixed[match(hksj$gene_id, fixed$gene_id), , drop = FALSE]
  if (anyNA(wald$gene_id) || anyNA(fixed$gene_id)) stop("Meta sensitivity methods do not align", call. = FALSE)
  data.frame(
    gene_id = as.character(hksj$gene_id), gene_symbol = hksj$gene_symbol,
    k = hksj$k, method_primary = hksj$method,
    hksj_estimate = hksj$estimate, hksj_meta_se = hksj$meta_se,
    hksj_ci_lb = hksj$ci_lb, hksj_ci_ub = hksj$ci_ub,
    hksj_prediction_lb = hksj$prediction_lb, hksj_prediction_ub = hksj$prediction_ub,
    hksj_p_value = hksj$p_value, hksj_fdr_bh = hksj$fdr_bh,
    tau2 = hksj$tau2, Q = hksj$Q, Q_p_value = hksj$Q_p_value, I2 = hksj$I2,
    wald_p_value = wald$p_value, wald_fdr_bh = wald$fdr_bh,
    fixed_p_value = fixed$p_value, fixed_fdr_bh = fixed$fdr_bh,
    direction = hksj$direction,
    heterogeneity_caveat = hksj$heterogeneity_caveat,
    primary_H1_rescue_allowed = FALSE,
    stringsAsFactors = FALSE
  )
}

task13_build_table_s <- function(targets, effects1, effects2, hksj) {
  targets$entrez_id <- as.character(targets$entrez_id)
  effects1$gene_id <- as.character(effects1$gene_id)
  effects2$gene_id <- as.character(effects2$gene_id)
  hksj$gene_id <- as.character(hksj$gene_id)
  e1 <- effects1[match(targets$entrez_id, effects1$gene_id), , drop = FALSE]
  e2 <- effects2[match(targets$entrez_id, effects2$gene_id), , drop = FALSE]
  meta <- hksj[match(targets$entrez_id, hksj$gene_id), , drop = FALSE]
  out <- data.frame(
    gene_id = targets$entrez_id,
    gene_symbol = targets$gene_symbol,
    GSE105450_log2FC = e1$log2FC, GSE105450_SE = e1$SE,
    GSE105450_moderated_t = e1$moderated_t, GSE105450_p_value = e1$p_value,
    GSE105450_fdr_bh = e1$fdr_bh,
    GSE103842_log2FC = e2$log2FC, GSE103842_SE = e2$SE,
    GSE103842_moderated_t = e2$moderated_t, GSE103842_p_value = e2$p_value,
    GSE103842_fdr_bh = e2$fdr_bh,
    meta_estimable = !is.na(meta$gene_id),
    HKSJ_estimate = meta$estimate, HKSJ_meta_SE = meta$meta_se,
    HKSJ_ci_lb = meta$ci_lb, HKSJ_ci_ub = meta$ci_ub,
    HKSJ_prediction_lb = meta$prediction_lb, HKSJ_prediction_ub = meta$prediction_ub,
    HKSJ_p_value = meta$p_value, HKSJ_fdr_bh = meta$fdr_bh,
    HKSJ_tau2 = meta$tau2, HKSJ_I2 = meta$I2,
    evidence_boundary = ifelse(
      !is.na(meta$gene_id),
      "complete two-cohort primary target meta; HKSJ-BH family includes all 358 estimable targets",
      "not estimable in both confirmatory cohorts; no imputation"
    ),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$gene_symbol, out$gene_id, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

task13_draw_box <- function(x, y, w, h, label, fill = "#F4F7FB", border = "#355C7D", cex = 0.85) {
  graphics::rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = border, lwd = 1.5)
  graphics::text(x, y, label, cex = cex)
}

task13_plot_workflow <- function(path) {
  task13_open_pdf(path, 10.5, 6.2)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(1, 1, 1, 1), family = "sans")
  graphics::plot.new(); graphics::plot.window(xlim = c(0, 10), ylim = c(0, 6))
  task13_draw_box(1.2, 4.7, 1.8, 0.9, "TCMSP + BATMAN\nprovenance audit")
  task13_draw_box(3.5, 4.7, 1.8, 0.9, "Frozen primary set\n520 Entrez genes", fill = "#E7F3E8", border = "#2E7D32")
  task13_draw_box(6.0, 5.2, 2.2, 0.9, "GSE105450\n122 modeled")
  task13_draw_box(6.0, 4.0, 2.2, 0.9, "GSE103842\n73 modeled")
  task13_draw_box(8.8, 4.6, 1.9, 1.3, "Competitive camera\n+ signed Stouffer\n+ conditional random sets*", fill = "#FFF4E5", border = "#B26A00", cex = 0.78)
  task13_draw_box(3.0, 2.2, 2.2, 1.0, "358 two-cohort targets\nREML/HKSJ meta\n0 BH-FDR < 0.05")
  task13_draw_box(6.0, 2.2, 2.4, 1.2, "roast/fgsea positive\n(broad RSV response or\ncell-composition signal;\nnot competitive/mechanistic)", cex = 0.72)
  task13_draw_box(8.8, 2.2, 1.8, 1.1, "Frozen decision\nH1 not supported", fill = "#FDECEC", border = "#B71C1C")
  arrows <- rbind(c(2.1,4.7,2.6,4.7), c(4.4,4.7,4.9,5.1), c(4.4,4.6,4.9,4.1),
                  c(7.1,5.2,7.8,4.8), c(7.1,4.0,7.8,4.4), c(3.5,4.2,3.1,2.8),
                  c(4.4,4.5,5.1,2.7), c(7.1,2.2,7.8,2.2), c(8.8,4.0,8.8,2.8))
  apply(arrows, 1, function(a) graphics::arrows(a[1], a[2], a[3], a[4], length = 0.08, col = "#555555", lwd = 1.3))
  graphics::mtext("Figure 1. Frozen analysis workflow and evidence boundaries", side = 3, line = -1.2, font = 2, cex = 1.05)
  graphics::text(5, 0.70, "* Random-set operation was finalized after primary outcomes were known but before random draws; its membership null differs from camera and cannot rescue H1.", cex = 0.67, col = "#7A1E1E")
  graphics::text(5, 0.35, "All matrices were analyzed separately; GSE188427 was not combined with overlapping GSE105450.", cex = 0.74, col = "#444444")
}

task13_plot_replication <- function(table2, path) {
  rows <- table2[table2$analysis_id != "matched_random_empirical_calibration", , drop = FALSE]
  labels <- c("GSE105450 camera", "GSE103842 camera", "Combined Stouffer")
  task13_open_pdf(path, 11.2, 6.6)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(10.5, 9.8, 3.8, 2.0), family = "sans")
  y <- rev(seq_len(nrow(rows)))
  statistic_range <- range(c(0, rows$statistic))
  xr <- c(statistic_range[[1L]] - 0.15, statistic_range[[2L]] + 5.3)
  graphics::plot(NA, xlim = xr, ylim = c(0.5, nrow(rows) + 0.7), yaxt = "n",
                 xlab = "Signed normal-score statistic (positive = Up)", ylab = "",
                 main = "Figure 2. Primary target-set replication")
  graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE)
  graphics::abline(v = 0, col = "#666666", lty = 2)
  graphics::points(rows$statistic, y, pch = 19, cex = 1.2, col = "#B23A48")
  info <- ifelse(
    is.na(rows$coverage),
    sprintf("P=%.3f | coverage %s", rows$p_value, rows$coverage_detail),
    sprintf("P=%.3f | coverage %.1f%% (%s)", rows$p_value, 100 * rows$coverage, rows$coverage_detail)
  )
  info <- paste0(info, " | CI: not estimable")
  graphics::text(max(rows$statistic) + 0.25, y, labels = info, adj = c(0, 0.5), cex = 0.70)
  random <- table2[table2$analysis_id == "matched_random_empirical_calibration", , drop = FALSE]
  graphics::mtext(
    "Conditional membership-randomization operation finalized after primary outcomes, but before random draws; its null differs from camera.",
    side = 1, line = 5.0, cex = 0.69, col = "#7A1E1E"
  )
  graphics::mtext(
    sprintf("Empirical P=%.5f arose against a narrow conditional random null; it cannot validate camera, independently replicate, or rescue H1.", random$p_value),
    side = 1, line = 6.3, cex = 0.69, col = "#7A1E1E"
  )
  graphics::mtext(
    "roast/fgsea strongly positive: compatible with broad RSV whole-blood or cell-composition signal, not competitive enrichment or formula mechanism.",
    side = 1, line = 7.6, cex = 0.67, col = "#6A3D00"
  )
  graphics::mtext(
    "camera/Stouffer/empirical rank tests have no prespecified target-set effect estimand; an effect CI is not estimable.",
    side = 3, line = 0.4, cex = 0.76, col = "#444444"
  )
}

task13_plot_all_target_effects <- function(hksj, path, per_page = 40L) {
  hksj <- hksj[order(hksj$gene_symbol, hksj$gene_id, na.last = TRUE), , drop = FALSE]
  pages <- split(seq_len(nrow(hksj)), ceiling(seq_len(nrow(hksj)) / per_page))
  xlim <- range(c(hksj$ci_lb, hksj$ci_ub), finite = TRUE)
  task13_open_pdf(path, 9.2, 11.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in seq_along(pages)) {
    one <- hksj[pages[[p]], , drop = FALSE]
    old <- graphics::par(mar = c(4.5, 8.8, 3.2, 1.2), family = "sans")
    y <- rev(seq_len(nrow(one)))
    graphics::plot(NA, xlim = xlim, ylim = c(0.5, nrow(one) + 0.5), yaxt = "n",
                   xlab = "REML/HKSJ pooled log2FC (95% CI)", ylab = "",
                   main = sprintf("Supplementary Figure S1. All estimable targets (%d/%d)", p, length(pages)))
    graphics::axis(2, at = y, labels = paste0(one$gene_symbol, " [", one$gene_id, "]"), las = 1, tick = FALSE, cex.axis = 0.64)
    graphics::abline(v = 0, lty = 2, col = "#777777")
    graphics::segments(one$ci_lb, y, one$ci_ub, y, col = "#567189", lwd = 1)
    graphics::points(one$estimate, y, pch = 19, cex = 0.55, col = "#B23A48")
    graphics::mtext("Alphabetical display only; no gene passed HKSJ BH-FDR < 0.05; k=2 intervals are unstable.", side = 3, line = 0.3, cex = 0.72)
    graphics::par(old)
  }
}

task13_plot_heterogeneity <- function(hksj, effects1, effects2, path) {
  effects1$gene_id <- as.character(effects1$gene_id); effects2$gene_id <- as.character(effects2$gene_id)
  e1 <- effects1[match(as.character(hksj$gene_id), effects1$gene_id), , drop = FALSE]
  e2 <- effects2[match(as.character(hksj$gene_id), effects2$gene_id), , drop = FALSE]
  task13_open_pdf(path, 10.5, 8.2)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.3, 3.4, 1.0), oma = c(1.5, 0, 1.5, 0), family = "sans")
  graphics::hist(hksj$tau2, breaks = 30, col = "#B8D8D8", border = "white", main = "A. REML tau-squared", xlab = "tau-squared")
  graphics::mtext(sprintf("Mass at zero: %d/%d", sum(hksj$tau2 == 0), nrow(hksj)), side = 3, line = 0.2, cex = 0.75)
  graphics::hist(hksj$I2, breaks = seq(0, 100, by = 5), col = "#F4C7AB", border = "white", main = "B. I-squared (descriptive)", xlab = "I-squared (%)", xlim = c(0, 100))
  graphics::mtext("I-squared=0 is not evidence of homogeneity when k=2", side = 3, line = 0.2, cex = 0.72, col = "#7A1E1E")
  pi_order <- order(hksj$estimate, hksj$gene_symbol)
  pi <- hksj[pi_order, , drop = FALSE]
  y_pi <- seq_len(nrow(pi))
  pi_xlim <- range(c(pi$prediction_lb, pi$prediction_ub), finite = TRUE)
  graphics::plot(NA, xlim = pi_xlim, ylim = c(1, nrow(pi)), yaxt = "n",
                 xlab = "REML/HKSJ log2FC scale", ylab = "Targets sorted by pooled effect",
                 main = "C. Actual 95% prediction intervals")
  graphics::abline(v = 0, lty = 2, col = "#777777")
  graphics::segments(pi$prediction_lb, y_pi, pi$prediction_ub, y_pi,
                     col = grDevices::adjustcolor("#355C7D", 0.34), lwd = 0.55)
  graphics::points(pi$estimate, y_pi, pch = 16, cex = 0.22,
                   col = grDevices::adjustcolor("#B23A48", 0.65))
  graphics::mtext("All 358 intervals shown; k=2 prediction intervals are unstable", side = 3, line = 0.2, cex = 0.70)
  lim <- range(c(e1$log2FC, e2$log2FC), finite = TRUE)
  graphics::plot(e1$log2FC, e2$log2FC, pch = 16, cex = 0.55, col = grDevices::adjustcolor("#2E7D32", 0.55),
                 xlim = lim, ylim = lim, xlab = "GSE105450 log2FC", ylab = "GSE103842 log2FC", main = "D. Cohort-specific target effects")
  graphics::abline(a = 0, b = 1, lty = 2, col = "#777777")
  graphics::abline(h = 0, v = 0, lty = 3, col = "#BBBBBB")
  graphics::mtext(sprintf("Same direction: %d/%d", sum(sign(e1$log2FC) == sign(e2$log2FC)), nrow(hksj)), side = 3, line = 0.2, cex = 0.75)
  graphics::mtext("Figure 4. Heterogeneity and precision summaries (all 358 targets; k=2 caveat applies)", outer = TRUE, side = 3, line = 0.2, font = 2)
}

task13_plot_sensitivity <- function(sensitivity, unavailable, path) {
  z <- ifelse(
    sensitivity$analysis_type == "signed_stouffer",
    sensitivity$z_combined,
    task13_signed_z(sensitivity$p_value, sensitivity$direction)
  )
  display_id <- sensitivity$analysis_id
  display_id <- gsub("__primary_520", " | primary", display_id, fixed = TRUE)
  display_id <- gsub("__expanded_1270", " | expanded", display_id, fixed = TRUE)
  display_id <- gsub("_plus_", " + ", display_id, fixed = TRUE)
  display_id <- gsub("GSE105450_hospital", "GSE105450 hospital", display_id, fixed = TRUE)
  display_id <- gsub("hospital +", "hospital +", display_id, fixed = TRUE)
  display_id <- gsub("replacement +", "GSE188427 +", display_id, fixed = TRUE)
  labels <- paste0(sensitivity$family, " | ", display_id)
  role_colors <- c(
    sensitivity = "#355C7D", exploratory_nonH1_unadjusted = "#B26A00",
    exploratory_nonH1 = "#B26A00", small_sample_exploratory = "#7B5EA7",
    cross_virus_exploratory = "#7B5EA7"
  )
  colors <- unname(role_colors[sensitivity$role])
  colors[is.na(colors)] <- "#555555"
  task13_open_pdf(path, 12.4, 8.2)
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(2.25, 1.35))
  old <- graphics::par(mar = c(7.0, 14.5, 3.5, 1.3), family = "sans")
  y <- rev(seq_len(nrow(sensitivity)))
  xr <- c(min(c(z, -1.96)) - 0.15, max(c(z, 1.96)) + 2.4)
  graphics::plot(NA, xlim = xr, ylim = c(0.4, nrow(sensitivity) + 0.7), yaxt = "n",
                 xlab = "Signed normal-score statistic", ylab = "", main = "Figure 5. Prespecified sensitivity results")
  graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.66)
  graphics::abline(v = 0, lty = 2, col = "#777777")
  graphics::abline(v = c(-1.96, 1.96), lty = 3, col = "#BBBBBB")
  graphics::points(z, y, pch = ifelse(sensitivity$p_holm < 0.05, 19, 1), cex = 1.05, col = colors)
  info <- sprintf("P=%.3g; Holm=%.3g%s", sensitivity$p_value, sensitivity$p_holm,
                  ifelse(is.na(sensitivity$coverage_fraction), "", sprintf("; cov=%.1f%%", 100 * sensitivity$coverage_fraction)))
  graphics::text(z, y, labels = info, pos = 4, cex = 0.59)
  graphics::mtext("Filled points: Holm P<0.05. All rows are non-rescuing; GSE188427 is unadjusted exploratory.", side = 1, line = 5.2, cex = 0.70, col = "#7A1E1E")
  graphics::par(mar = c(2, 1.0, 3.5, 1.0))
  graphics::plot.new()
  graphics::text(0.02, 0.95, "Unavailable by frozen contract", adj = c(0, 1), font = 2, cex = 0.95)
  y0 <- 0.80
  for (i in seq_len(nrow(unavailable))) {
    graphics::text(0.02, y0, unavailable$analysis[[i]], adj = c(0, 1), font = 2, cex = 0.72)
    reason <- paste(strwrap(unavailable$reason[[i]], width = 42), collapse = "\n")
    graphics::text(0.02, y0 - 0.055, reason, adj = c(0, 1), cex = 0.59)
    y0 <- y0 - 0.24
  }
  boundary <- paste(strwrap(
    "No sensitivity or exploratory result may rescue the failed primary H1 decision.", width = 42
  ), collapse = "\n")
  graphics::text(0.02, 0.04, boundary, adj = c(0, 0), cex = 0.64, col = "#7A1E1E")
  graphics::par(old)
}

task13_manifest <- function(inputs, contextual_dependencies, script, outputs) {
  input_rows <- data.frame(
    record_type = "input", name = basename(inputs), path = inputs,
    sha256 = vapply(inputs, task13_sha256, character(1)),
    version = NA_character_, seed = NA_integer_,
    note = "immutable reviewed input consumed by Task 13", stringsAsFactors = FALSE
  )
  task10_amendment <- input_rows$name %in% c(
    "2026-08-05-random-set-cross-cohort-operationalization.md",
    "2026-08-05-random-set-quintile-ties.md"
  )
  input_rows$note[task10_amendment] <- paste0(
    "post-outcome operational clarification: primary outcomes were known; ",
    "cross-cohort random-set rules were frozen before random draws"
  )
  script_rows <- data.frame(
    record_type = "script", name = basename(script), path = script,
    sha256 = task13_sha256(script), version = NA_character_, seed = NA_integer_,
    note = "Task 13 table and figure generator", stringsAsFactors = FALSE
  )
  contextual_rows <- data.frame(
    record_type = "archived_contextual_dependency",
    name = basename(contextual_dependencies), path = contextual_dependencies,
    sha256 = vapply(contextual_dependencies, task13_sha256, character(1)),
    version = NA_character_, seed = NA_integer_,
    note = c(
      "archived cohort-level missingness table; Task13 summaries are recomputed from the frozen manifest using explicit denominators",
      "archived Task12 design-sensitivity result; not read to construct Task13 tables or figures",
      "archived Task12 observed-precision result; not read to construct Task13 tables or figures"
    ),
    stringsAsFactors = FALSE
  )
  packages <- c("R", "graphics", "grDevices", "stats", "utils")
  versions <- c(R.version.string, vapply(packages[-1L], function(x) as.character(utils::packageVersion(x)), character(1)))
  package_rows <- data.frame(
    record_type = "package", name = packages, path = NA_character_, sha256 = NA_character_,
    version = versions, seed = NA_integer_, note = "runtime used for deterministic base-R rendering",
    stringsAsFactors = FALSE
  )
  seed_rows <- data.frame(
    record_type = "seed", name = c("Task10 matched random sets", "Task12 fixed-N design sensitivity", "Task13 rendering"),
    path = NA_character_, sha256 = NA_character_, version = NA_character_,
    seed = c(20260804L, 20260812L, NA_integer_),
    note = c(
      "post-outcome operational clarification frozen before random draws; conditional membership sensitivity only; cannot validate camera, independently replicate, or rescue H1",
      "inherited reviewed fixed-N conditional design sensitivity; not post-hoc power; cannot rescue primary H1",
      "no random number generation"
    ),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path("results", "reproducibility_manifest.tsv")
  output_sha <- vapply(outputs, task13_sha256, character(1))
  output_rows <- data.frame(
    record_type = "output", name = basename(outputs), path = outputs, sha256 = output_sha,
    version = NA_character_, seed = NA_integer_, note = ifelse(
      basename(outputs) == "FigS1_all_target_effects.pdf",
      "nine-page alphabetical display of all 358 estimable targets; no significance-based filling",
      "deterministic Task 13 output"
    ),
    stringsAsFactors = FALSE
  )
  self <- data.frame(
    record_type = "output", name = basename(manifest_path), path = manifest_path,
    sha256 = NA_character_, version = NA_character_, seed = NA_integer_,
    note = "self checksum omitted to avoid recursive identity", stringsAsFactors = FALSE
  )
  cancelled <- data.frame(
    record_type = "output", name = "Fig3_meta_forest.pdf", path = NA_character_,
    sha256 = NA_character_, version = NA_character_, seed = NA_integer_,
    note = "cancelled_not_generated: zero HKSJ BH-FDR targets; no non-significant top-20 filling",
    stringsAsFactors = FALSE
  )
  rbind(input_rows, contextual_rows, script_rows, package_rows, seed_rows, output_rows, cancelled, self)
}

main_task13 <- function() {
  inputs <- c(
    file.path("data", "clean", "sample_manifest_frozen.tsv"),
    file.path("data", "clean", "bfff_targets_primary.tsv"),
    file.path("results", "tables", "cohort_inclusion_flow.tsv"),
    file.path("results", "cohort", "GSE105450_gene_effects.tsv"),
    file.path("results", "cohort", "GSE103842_gene_effects.tsv"),
    file.path("results", "cohort", "GSE105450_model_diagnostics.tsv"),
    file.path("results", "cohort", "GSE103842_model_diagnostics.tsv"),
    file.path("results", "tables", "GSE105450_qc.tsv"),
    file.path("results", "tables", "GSE103842_qc.tsv"),
    file.path("results", "cohort", "target_set_camera.tsv"),
    file.path("results", "cohort", "target_set_primary_decision_partial.tsv"),
    file.path("results", "meta", "random_set_empirical_calibration.tsv"),
    file.path("results", "meta", "target_gene_meta_reml_knha.tsv"),
    file.path("results", "meta", "target_gene_meta_reml_wald.tsv"),
    file.path("results", "meta", "target_gene_meta_fixed.tsv"),
    file.path("results", "sensitivity", "all_sensitivity_results.tsv"),
    file.path("results", "sensitivity", "cohort_model_diagnostics.tsv"),
    file.path("results", "sensitivity", "cohort_qc_audit.tsv"),
    file.path("results", "sensitivity", "not_implemented.tsv"),
    file.path("results", "cohort", "target_set_roast.tsv"),
    file.path("results", "cohort", "target_set_fgsea_exploratory.tsv"),
    file.path("docs", "science-superpowers", "preregistrations", "amendments", "2026-08-05-random-set-cross-cohort-operationalization.md"),
    file.path("docs", "science-superpowers", "preregistrations", "amendments", "2026-08-05-random-set-quintile-ties.md")
  )
  contextual_dependencies <- c(
    file.path("results", "tables", "sample_covariate_missingness.tsv"),
    file.path("results", "power", "fixed_n_design_sensitivity.tsv"),
    file.path("results", "power", "observed_precision_summary.tsv")
  )
  if (!all(file.exists(c(inputs, contextual_dependencies)))) stop("Task 13 reviewed inputs are incomplete", call. = FALSE)

  manifest <- task13_read_tsv(inputs[[1L]])
  targets <- task13_read_tsv(inputs[[2L]])
  flow <- task13_read_tsv(inputs[[3L]])
  effects1 <- task13_read_tsv(inputs[[4L]])
  effects2 <- task13_read_tsv(inputs[[5L]])
  qc_gse105450 <- task13_read_tsv(inputs[[8L]])
  qc_gse103842 <- task13_read_tsv(inputs[[9L]])
  qc_gse105450$cohort <- "GSE105450"
  qc_gse103842$cohort <- "GSE103842"
  qc_confirmatory <- rbind(qc_gse105450, qc_gse103842)
  camera <- task13_read_tsv(inputs[[10L]])
  decision <- task13_read_tsv(inputs[[11L]])
  empirical <- task13_read_tsv(inputs[[12L]])
  hksj <- task13_read_tsv(inputs[[13L]])
  wald <- task13_read_tsv(inputs[[14L]])
  fixed <- task13_read_tsv(inputs[[15L]])
  sensitivity <- task13_read_tsv(inputs[[16L]])
  sensitivity_diagnostics <- task13_read_tsv(inputs[[17L]])
  qc_sensitivity <- task13_read_tsv(inputs[[18L]])
  unavailable <- task13_read_tsv(inputs[[19L]])
  roast <- task13_read_tsv(inputs[[20L]])
  fgsea <- task13_read_tsv(inputs[[21L]])
  task10_operationalization <- paste(readLines(inputs[[22L]], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  task10_ties <- paste(readLines(inputs[[23L]], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (!grepl("Task 8 已判定", task10_operationalization, fixed = TRUE) ||
      !grepl("在生成任何随机集合前冻结", task10_operationalization, fixed = TRUE) ||
      !grepl("在生成任何随机集合前冻结", task10_ties, fixed = TRUE)) {
    stop("Task 10 post-outcome/before-draw timing evidence is absent", call. = FALSE)
  }
  task10_timing <- paste0(
    "complete cross-cohort operationalization was clarified after primary outcomes were known, ",
    "but frozen before any random draws"
  )
  supplementary_signal_context <- task13_supplementary_signal_context(roast, fgsea)
  qc_audit <- rbind(qc_confirmatory[, c("sample_id", "exclude_qc", "cohort")],
                    qc_sensitivity[, c("sample_id", "exclude_qc", "cohort")])

  if (nrow(targets) != 520L || nrow(hksj) != 358L || sum(hksj$fdr_bh < 0.05) != 0L ||
      !identical(decision$h1_supported, FALSE) || !identical(empirical$random_calibration_can_rescue, FALSE)) {
    stop("Task 13 evidence boundary check failed", call. = FALSE)
  }

  table1 <- task13_build_table1(manifest, flow, sensitivity_diagnostics, qc_audit)
  table2 <- task13_build_table2(
    camera, decision, empirical, task10_timing, supplementary_signal_context
  )
  table3 <- task13_build_table3(hksj, wald, fixed)
  table_s <- task13_build_table_s(targets, effects1, effects2, hksj)

  table_paths <- c(
    file.path("results", "tables", "Table1_cohorts.tsv"),
    file.path("results", "tables", "Table2_target_set_tests.tsv"),
    file.path("results", "tables", "Table3_target_meta.tsv"),
    file.path("results", "tables", "TableS_all_effects.tsv")
  )
  task13_write_tsv(table1, table_paths[[1L]])
  task13_write_tsv(table2, table_paths[[2L]])
  task13_write_tsv(table3, table_paths[[3L]])
  task13_write_tsv(table_s, table_paths[[4L]])

  figure_paths <- c(
    file.path("results", "figures", "Fig1_workflow.pdf"),
    file.path("results", "figures", "Fig2_target_set_replication.pdf"),
    file.path("results", "figures", "FigS1_all_target_effects.pdf"),
    file.path("results", "figures", "Fig4_heterogeneity.pdf"),
    file.path("results", "figures", "Fig5_sensitivity.pdf")
  )
  task13_plot_workflow(figure_paths[[1L]])
  task13_plot_replication(table2, figure_paths[[2L]])
  task13_plot_all_target_effects(hksj, figure_paths[[3L]])
  task13_plot_heterogeneity(hksj, effects1, effects2, figure_paths[[4L]])
  task13_plot_sensitivity(sensitivity, unavailable, figure_paths[[5L]])

  script <- file.path("analysis", "R", "11_make_tables_figures.R")
  reproducibility <- task13_manifest(inputs, contextual_dependencies, script, c(table_paths, figure_paths))
  task13_write_tsv(reproducibility, file.path("results", "reproducibility_manifest.tsv"))
  invisible(list(tables = table_paths, figures = figure_paths))
}

if (!identical(Sys.getenv("TASK13_SKIP_MAIN"), "1")) main_task13()
