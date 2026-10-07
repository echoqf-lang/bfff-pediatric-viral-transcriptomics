#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

read_tsv <- function(path) {
  read.delim(
    path, sep = "\t", quote = "", check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8"
  )
}

cohorts <- c("GSE105450", "GSE103842")
targets <- read_tsv(file.path("data", "clean", "bfff_targets_primary.tsv"))
effects <- setNames(lapply(cohorts, function(cohort) {
  read_tsv(file.path("results", "cohort", paste0(cohort, "_gene_effects.tsv")))
}), cohorts)
common <- Reduce(
  intersect,
  c(list(as.character(targets$entrez_id)), lapply(effects, function(x) as.character(x$gene_id)))
)
common <- as.character(targets$entrez_id[targets$entrez_id %in% common])

paths <- c(
  knha = file.path("results", "meta", "target_gene_meta_reml_knha.tsv"),
  wald = file.path("results", "meta", "target_gene_meta_reml_wald.tsv"),
  fixed = file.path("results", "meta", "target_gene_meta_fixed.tsv"),
  loo = file.path("results", "meta", "target_gene_leave_one_out.tsv"),
  diagnostics = file.path("results", "meta", "target_gene_meta_diagnostics.tsv")
)
stopifnot(all(file.exists(paths)))
out <- lapply(paths, read_tsv)
expected_schema <- c(
  "gene_id", "gene_symbol", "k", "estimate", "meta_se", "ci_lb", "ci_ub",
  "prediction_lb", "prediction_ub", "p_value", "tau2", "Q", "Q_p_value",
  "I2", "method", "direction", "expected_cohorts", "heterogeneity_caveat", "fdr_bh"
)
for (name in c("knha", "wald", "fixed")) {
  x <- out[[name]]
  stopifnot(
    identical(names(x), expected_schema),
    nrow(x) == length(common),
    identical(as.character(x$gene_id), common),
    !anyDuplicated(x$gene_id),
    all(x$k == 2L),
    all(x$direction == "RSV_minus_healthy"),
    all(x$expected_cohorts == paste(cohorts, collapse = "|")),
    all(is.finite(x$estimate)),
    all(is.finite(x$meta_se) & x$meta_se > 0),
    all(is.finite(x$ci_lb) & is.finite(x$ci_ub)),
    all(x$ci_lb <= x$estimate & x$estimate <= x$ci_ub),
    all(is.finite(x$p_value) & x$p_value >= 0 & x$p_value <= 1),
    all(is.finite(x$fdr_bh) & x$fdr_bh >= 0 & x$fdr_bh <= 1),
    isTRUE(all.equal(x$fdr_bh, p.adjust(x$p_value, method = "BH"), tolerance = 1e-15)),
    all(grepl("k=2", x$heterogeneity_caveat, fixed = TRUE))
  )
}
stopifnot(
  all(out$knha$method == "REML_HKSJ"),
  all(out$wald$method == "REML_Wald"),
  all(out$fixed$method == "fixed_inverse_variance"),
  all(is.finite(out$knha$prediction_lb)),
  all(is.finite(out$knha$prediction_ub)),
  all(out$knha$prediction_lb <= out$knha$estimate),
  all(out$knha$estimate <= out$knha$prediction_ub),
  all(is.na(out$wald$prediction_lb)), all(is.na(out$wald$prediction_ub)),
  all(is.na(out$fixed$prediction_lb)), all(is.na(out$fixed$prediction_ub))
)

# Independent real-input reference fits at both ends of the frozen target order.
for (gene_id in common[c(1L, length(common))]) {
  one <- do.call(rbind, lapply(cohorts, function(cohort) {
    x <- effects[[cohort]]
    x[x$gene_id == gene_id, c("gene_id", "cohort", "direction", "log2FC", "SE", "moderated_t")]
  }))
  stopifnot(
    nrow(one) == 2L,
    identical(as.character(one$cohort), cohorts),
    all(one$direction == "RSV_minus_healthy"),
    all(is.finite(one$SE) & one$SE > 0),
    isTRUE(all.equal(one$moderated_t, one$log2FC / one$SE, tolerance = 1e-10))
  )
  direct_knha <- metafor::rma.uni(yi = one$log2FC, sei = one$SE, method = "REML", test = "knha")
  direct_wald <- metafor::rma.uni(yi = one$log2FC, sei = one$SE, method = "REML", test = "z")
  direct_fixed <- metafor::rma.uni(yi = one$log2FC, sei = one$SE, method = "FE", test = "z")
  direct_prediction <- predict(direct_knha, level = 95)
  got_knha <- out$knha[out$knha$gene_id == gene_id, , drop = FALSE]
  got_wald <- out$wald[out$wald$gene_id == gene_id, , drop = FALSE]
  got_fixed <- out$fixed[out$fixed$gene_id == gene_id, , drop = FALSE]
  stopifnot(
    isTRUE(all.equal(got_knha$estimate, as.numeric(direct_knha$b[[1L]]), tolerance = 1e-12)),
    isTRUE(all.equal(got_knha$p_value, direct_knha$pval, tolerance = 1e-12)),
    isTRUE(all.equal(got_knha$prediction_lb, as.numeric(direct_prediction$pi.lb), tolerance = 1e-12)),
    isTRUE(all.equal(got_knha$prediction_ub, as.numeric(direct_prediction$pi.ub), tolerance = 1e-12)),
    isTRUE(all.equal(got_wald$p_value, direct_wald$pval, tolerance = 1e-12)),
    isTRUE(all.equal(got_fixed$estimate, as.numeric(direct_fixed$b[[1L]]), tolerance = 1e-12)),
    isTRUE(all.equal(got_fixed$p_value, direct_fixed$pval, tolerance = 1e-12))
  )
}

stopifnot(
  nrow(out$loo) == 1L,
  out$loo$k == 2L,
  identical(out$loo$status, "not_applicable_k2_confirmatory"),
  identical(out$loo$direction, "RSV_minus_healthy")
)
diag <- setNames(as.character(out$diagnostics$value), out$diagnostics$metric)
wide <- reshape(
  do.call(rbind, lapply(cohorts, function(cohort) {
    x <- effects[[cohort]]
    x[x$gene_id %in% common, c("gene_id", "cohort", "log2FC")]
  })),
  idvar = "gene_id", timevar = "cohort", direction = "wide"
)
same_sign <- sign(wide[[paste0("log2FC.", cohorts[[1L]])]]) ==
  sign(wide[[paste0("log2FC.", cohorts[[2L]])]])
stopifnot(
  as.integer(diag[["frozen_primary_targets"]]) == 520L,
  as.integer(diag[["common_targets_meta_analyzed"]]) == length(common),
  as.integer(diag[["direction_consistent_genes"]]) == sum(same_sign),
  as.integer(diag[["hksj_fdr_lt_0_05"]]) == sum(out$knha$fdr_bh < 0.05),
  identical(diag[["confirmatory_leave_one_out_status"]], "not_applicable_k2_confirmatory"),
  identical(diag[["primary_H1_rescue_allowed"]], "FALSE")
)

checksum <- read_tsv(file.path("logs", "checksums", "task9_sha256.tsv"))
stopifnot(!anyDuplicated(checksum$artifact))
for (i in seq_len(nrow(checksum))) {
  stopifnot(file.exists(checksum$artifact[[i]]))
  actual <- system2("shasum", c("-a", "256", checksum$artifact[[i]]), stdout = TRUE)
  actual <- sub("[[:space:]].*$", "", actual[[1L]])
  stopifnot(identical(actual, checksum$sha256[[i]]))
}

script <- paste(readLines(file.path("analysis", "R", "07_meta_targets.R"), warn = FALSE), collapse = "\n")
stopifnot(
  grepl('task9_cohorts <- c\\("GSE105450", "GSE103842"\\)', script, perl = TRUE),
  !grepl("run_camera_analysis\\(|test_camera_target_set\\(|fgsea::|random_gene", script, perl = TRUE)
)

cat("task9 real-output verification passed\n")
