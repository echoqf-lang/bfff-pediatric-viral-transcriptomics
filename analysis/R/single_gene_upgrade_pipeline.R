options(stringsAsFactors = FALSE)

read_upgrade_tsv <- function(path, col_classes = NA) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = col_classes)
}

write_upgrade_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

bh_family <- function(p) stats::p.adjust(p, method = "BH")

assert_unique_gene_effects <- function(x, cohort_name, expected_direction = "RSV_minus_healthy") {
  required <- c("gene_id", "cohort", "direction", "log2FC", "SE")
  absent <- setdiff(required, names(x))
  if (length(absent)) stop("Effect table missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  if (nrow(x) == 0L || anyNA(x$gene_id) || any(!nzchar(as.character(x$gene_id)))) stop("gene_id is incomplete", call. = FALSE)
  if (anyDuplicated(as.character(x$gene_id))) stop("gene_id must be unique within cohort", call. = FALSE)
  if (!identical(unique(as.character(x$cohort)), cohort_name)) stop("Unexpected cohort label", call. = FALSE)
  if (any(as.character(x$direction) != expected_direction)) stop("Unexpected effect direction", call. = FALSE)
  if (any(!is.finite(x$log2FC)) || any(!is.finite(x$SE) | x$SE <= 0)) stop("Effects require finite log2FC and positive finite SE", call. = FALSE)
  invisible(TRUE)
}

classify_direction_pattern <- function(beta_1, beta_2, beta_3) {
  signs <- c(sign(beta_1), sign(beta_2), sign(beta_3))
  if (any(!is.finite(signs))) return("not_estimable")
  if (all(signs > 0)) return("all_up")
  if (all(signs < 0)) return("all_down")
  if (sum(signs > 0) == 2L) return("two_up_one_down")
  if (sum(signs < 0) == 2L) return("two_down_one_up")
  "contains_zero"
}

fit_hksj_meta <- function(effect_table, expected_cohorts) {
  if (!requireNamespace("metafor", quietly = TRUE)) stop("metafor package is required", call. = FALSE)
  expected_cohorts <- as.character(expected_cohorts)
  if (length(expected_cohorts) != 3L || anyDuplicated(expected_cohorts)) stop("Exactly three unique cohorts are required", call. = FALSE)
  if (anyDuplicated(effect_table[c("gene_id", "cohort")])) stop("Duplicate gene-cohort rows", call. = FALSE)
  genes <- unique(as.character(effect_table$gene_id))
  out <- lapply(genes, function(id) {
    one <- effect_table[as.character(effect_table$gene_id) == id, , drop = FALSE]
    if (nrow(one) != 3L || !setequal(as.character(one$cohort), expected_cohorts)) return(NULL)
    one <- one[match(expected_cohorts, one$cohort), , drop = FALSE]
    if (any(!is.finite(one$log2FC)) || any(!is.finite(one$SE) | one$SE <= 0)) return(NULL)
    fit <- metafor::rma.uni(yi = one$log2FC, sei = one$SE, method = "REML", test = "knha")
    pred <- tryCatch(stats::predict(fit, level = 95), error = function(e) NULL)
    fit_adhoc <- metafor::rma.uni(yi = one$log2FC, sei = one$SE, method = "REML", test = "adhoc")
    data.frame(
      gene_id = id,
      k = fit$k,
      estimate = as.numeric(fit$b),
      se = fit$se,
      ci_low = fit$ci.lb,
      ci_high = fit$ci.ub,
      prediction_low = if (is.null(pred)) NA_real_ else pred$pi.lb,
      prediction_high = if (is.null(pred)) NA_real_ else pred$pi.ub,
      p_value = fit$pval,
      adhoc_se = fit_adhoc$se,
      adhoc_ci_low = fit_adhoc$ci.lb,
      adhoc_ci_high = fit_adhoc$ci.ub,
      adhoc_p_value = fit_adhoc$pval,
      tau2 = fit$tau2,
      i2 = fit$I2,
      direction_pattern = classify_direction_pattern(one$log2FC[1], one$log2FC[2], one$log2FC[3]),
      method = "REML_HKSJ",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0L) stop("No genes have complete three-cohort effects", call. = FALSE)
  out$fdr_bh <- bh_family(out$p_value)
  out$adhoc_fdr_bh <- bh_family(out$adhoc_p_value)
  out
}
