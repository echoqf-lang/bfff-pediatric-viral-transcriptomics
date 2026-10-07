#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(TASK9_SKIP_MAIN = "1")
source(file.path("analysis", "R", "07_meta_targets.R"), local = FALSE)

expect_error <- function(code, pattern) {
  message <- tryCatch({
    force(code)
    NA_character_
  }, error = function(e) conditionMessage(e))
  if (is.na(message) || !grepl(pattern, message, perl = TRUE)) {
    stop("Expected error matching '", pattern, "' but got: ", message, call. = FALSE)
  }
}

fixture <- data.frame(
  gene_id = rep(c("101", "202"), each = 2L),
  gene_symbol = rep(c("GENEA", "GENEB"), each = 2L),
  cohort = rep(task9_cohorts, 2L),
  direction = "RSV_minus_healthy",
  log2FC = c(0.20, 0.35, -0.30, -0.10),
  SE = c(0.10, 0.15, 0.12, 0.20),
  moderated_t = c(2.00, 0.35 / 0.15, -2.50, -0.50),
  stringsAsFactors = FALSE
)

validated <- validate_meta_input(fixture, task9_cohorts)
fits <- fit_target_gene_meta(validated, task9_cohorts)
stopifnot(
  nrow(fits$knha) == 2L,
  nrow(fits$wald) == 2L,
  nrow(fits$fixed) == 2L,
  all(fits$knha$k == 2L),
  all(fits$wald$k == 2L),
  all(fits$fixed$k == 2L),
  all(fits$knha$direction == "RSV_minus_healthy"),
  all(fits$knha$expected_cohorts == paste(task9_cohorts, collapse = "|")),
  all(fits$knha$fdr_bh == p.adjust(fits$knha$p_value, method = "BH")),
  all(fits$wald$fdr_bh == p.adjust(fits$wald$p_value, method = "BH")),
  all(fits$fixed$fdr_bh == p.adjust(fits$fixed$p_value, method = "BH")),
  all(fits$knha$prediction_lb <= fits$knha$estimate),
  all(fits$knha$estimate <= fits$knha$prediction_ub)
)

one <- fixture[fixture$gene_id == "101", , drop = FALSE]
direct_knha <- metafor::rma.uni(
  yi = one$log2FC, sei = one$SE, method = "REML", test = "knha"
)
direct_wald <- metafor::rma.uni(
  yi = one$log2FC, sei = one$SE, method = "REML", test = "z"
)
direct_fixed <- metafor::rma.uni(
  yi = one$log2FC, sei = one$SE, method = "FE", test = "z"
)
direct_prediction <- predict(direct_knha, level = 95)
got_knha <- fits$knha[fits$knha$gene_id == "101", , drop = FALSE]
got_wald <- fits$wald[fits$wald$gene_id == "101", , drop = FALSE]
got_fixed <- fits$fixed[fits$fixed$gene_id == "101", , drop = FALSE]
stopifnot(
  isTRUE(all.equal(
    unname(unlist(got_knha[c(
      "estimate", "meta_se", "ci_lb", "ci_ub", "prediction_lb", "prediction_ub",
      "p_value", "tau2", "Q", "Q_p_value", "I2"
    )])),
    unname(c(
      direct_knha$b[[1L]], direct_knha$se, direct_knha$ci.lb, direct_knha$ci.ub,
      direct_prediction$pi.lb, direct_prediction$pi.ub, direct_knha$pval,
      direct_knha$tau2, direct_knha$QE, direct_knha$QEp, direct_knha$I2
    )), tolerance = 1e-12
  )),
  isTRUE(all.equal(
    unname(unlist(got_wald[c("estimate", "meta_se", "ci_lb", "ci_ub", "p_value", "tau2")])),
    unname(c(direct_wald$b[[1L]], direct_wald$se, direct_wald$ci.lb,
             direct_wald$ci.ub, direct_wald$pval, direct_wald$tau2)),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    unname(unlist(got_fixed[c("estimate", "meta_se", "ci_lb", "ci_ub", "p_value")])),
    unname(c(direct_fixed$b[[1L]], direct_fixed$se, direct_fixed$ci.lb,
             direct_fixed$ci.ub, direct_fixed$pval)),
    tolerance = 1e-12
  ))
)

duplicate <- rbind(fixture, fixture[1L, , drop = FALSE])
expect_error(validate_meta_input(duplicate, task9_cohorts), "duplicate gene-cohort")
missing <- fixture[-1L, , drop = FALSE]
expect_error(validate_meta_input(missing, task9_cohorts), "exactly the two expected cohorts")
extra <- rbind(fixture, transform(fixture[1L, , drop = FALSE], cohort = "GSE38900"))
expect_error(validate_meta_input(extra, task9_cohorts), "confirmatory cohorts only")
wrong_direction <- fixture
wrong_direction$direction[1L] <- "healthy_minus_RSV"
expect_error(validate_meta_input(wrong_direction, task9_cohorts), "RSV_minus_healthy")
bad_se <- fixture
bad_se$SE[1L] <- 0
expect_error(validate_meta_input(bad_se, task9_cohorts), "positive finite SE")
mismatch <- fixture
mismatch$moderated_t[1L] <- 999
expect_error(validate_meta_input(mismatch, task9_cohorts), "posterior SE does not correspond")
expect_error(validate_meta_input(fixture, c("GSE105450")), "exactly GSE105450 and GSE103842")

cat("task9 unit tests passed\n")
