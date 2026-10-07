#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(TASK11_SKIP_MAIN = "1")
source(file.path("analysis", "R", "09_sensitivity_analyses.R"), local = FALSE)

set.seed(20260805L)
y <- matrix(
  rnorm(100L * 20L), nrow = 100L,
  dimnames = list(as.character(seq_len(100L)), paste0("S", seq_len(20L)))
)
group <- factor(rep(c("healthy", "RSV"), each = 10L), levels = c("healthy", "RSV"))
design <- model.matrix(~ group)
contrast <- c(0, 1)
index <- seq_len(20L)
wrapped <- task11_camera(y, index, design, contrast)
direct <- limma::camera(
  y = y, index = index, design = design, contrast = contrast,
  directional = TRUE, inter.gene.cor = NA_real_
)
stopifnot(identical(wrapped, direct))

fit <- task11_fit_limma(y, design, "groupRSV", "fixture", "RSV_minus_healthy")
stopifnot(
  nrow(fit$gene_effects) == nrow(y),
  all(fit$gene_effects$direction == "RSV_minus_healthy"),
  isTRUE(all.equal(
    fit$gene_effects$moderated_t,
    fit$gene_effects$log2FC / fit$gene_effects$SE,
    tolerance = 1e-10
  ))
)

stouffer <- task11_signed_stouffer(c(0.04, 0.10), c("Up", "Down"), c(89, 73))
expected_z <- sum(sqrt(c(89, 73)) * c(
  qnorm(0.02, lower.tail = FALSE),
  -qnorm(0.05, lower.tail = FALSE)
)) / sqrt(89 + 73)
stopifnot(
  isTRUE(all.equal(stouffer$z_combined, expected_z, tolerance = 0)),
  isTRUE(all.equal(stouffer$p_two_sided, 2 * pnorm(abs(expected_z), lower.tail = FALSE), tolerance = 0))
)

targets_fixture <- data.frame(
  entrez_id = c("1", "2", "3", "4"),
  ensembl_gene_ids = c(
    "ENSG00000000001;ENSG00000000002",
    "ENSG00000000003",
    "ENSG00000000004",
    "ENSG00000000004"
  ),
  stringsAsFactors = FALSE
)
background <- sprintf("ENSG%011d", 1:4)
selection <- task11_select_rnaseq_target_units(background, c(7, 7, 4, 20), targets_fixture)
stopifnot(
  identical(selection$ambiguous_ensembl, "ENSG00000000004"),
  identical(selection$selection$entrez_id, c("1", "2")),
  identical(selection$selection$ensembl_id, c("ENSG00000000001", "ENSG00000000003"))
)

actual_targets <- task11_read_tsv(file.path("data", "clean", "bfff_targets_sensitivity.tsv"))
actual_targets$ensembl_gene_ids[is.na(actual_targets$ensembl_gene_ids)] <- ""
actual_aliases <- task11_build_alias_map(actual_targets)
stopifnot(identical(
  actual_aliases$ambiguous_ensembl,
  c("ENSG00000223572", "ENSG00000237289")
))
primary <- task11_target_index(selection$selection, "1", background)
expanded <- task11_target_index(selection$selection, c("1", "2", "3", "4"), background)
stopifnot(
  primary$coverage_unique_entrez == length(primary$index),
  expanded$coverage_unique_entrez == length(expanded$index),
  all(primary$index %in% expanded$index)
)

bad_combo <- try(
  task11_assert_meta_pair(c("GSE105450", "GSE188427")),
  silent = TRUE
)
stopifnot(inherits(bad_combo, "try-error"))
stopifnot(
  identical(
    task11_assert_meta_pair(c("GSE105450_hospital", "GSE103842")),
    c("GSE105450_hospital", "GSE103842")
  ),
  identical(
    task11_assert_meta_pair(c("GSE188427", "GSE103842")),
    c("GSE188427", "GSE103842")
  )
)

cat("Task 11 sensitivity-analysis unit tests passed\n")
