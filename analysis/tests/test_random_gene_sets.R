#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 2)

script <- file.path("analysis", "R", "08_random_gene_sets.R")
if (!file.exists(script)) stop("Task 10 implementation is missing", call. = FALSE)
old_skip <- Sys.getenv("TASK10_SKIP_MAIN", unset = NA_character_)
Sys.setenv(TASK10_SKIP_MAIN = "1")
source(script, local = FALSE)
if (is.na(old_skip)) Sys.unsetenv("TASK10_SKIP_MAIN") else Sys.setenv(TASK10_SKIP_MAIN = old_skip)

expect_true <- function(x, label) {
  if (!isTRUE(x)) stop(label, call. = FALSE)
}

expect_equal <- function(x, y, tolerance = 1e-12, label = "values differ") {
  if (!isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))) {
    stop(label, ": ", paste(capture.output(all.equal(x, y, tolerance = tolerance)), collapse = " "), call. = FALSE)
  }
}

# Tied values stay together; bins are deterministic and remain in 1..5.
bins <- empirical_quintile(c(1, 1, 1, 2, 3, 4, 5, 6, 7, 8))
expect_true(length(bins) == 10L && all(bins %in% 1:5), "quintile bins invalid")
expect_true(length(unique(bins[1:3])) == 1L, "equal values split across quintile bins")
expect_true(identical(bins, empirical_quintile(c(1, 1, 1, 2, 3, 4, 5, 6, 7, 8))), "quintiles not deterministic")

# Synthetic allocation requires exact, then probe-only, then adjacent-detection relaxation.
synthetic <- data.frame(
  gene_id = as.character(1:12),
  expression_bin = 1L,
  detection_bin = c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L),
  probe_bin = c(1L, 1L, 2L, 2L, 3L, 1L, 1L, 2L, 2L, 3L, 1L, 2L),
  is_target = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
plan <- build_relaxation_plan(synthetic)
expect_true(sum(plan$n_allocate) == sum(synthetic$is_target), "allocation size differs from target size")
expect_true(all(plan$source_expression_bin == plan$target_expression_bin), "expression bin was relaxed")
expect_true(all(plan$relaxation_stage %in% c("exact", "probe_only", "adjacent_detection")), "unknown relaxation stage")
expect_true(any(plan$relaxation_stage == "probe_only"), "probe-only relaxation not exercised")
expect_true(any(plan$relaxation_stage == "adjacent_detection"), "detection relaxation not exercised")

set.seed(20260804L)
draws <- draw_matched_sets(synthetic, plan, n_sets = 20L)
target_ids <- synthetic$gene_id[synthetic$is_target]
expect_true(length(draws) == 20L, "wrong number of draws")
expect_true(all(lengths(draws) == length(target_ids)), "draw size differs from target size")
expect_true(all(vapply(draws, function(x) !anyDuplicated(x), logical(1))), "within-set duplicate")
expect_true(all(vapply(draws, function(x) !length(intersect(x, target_ids)), logical(1))), "target contamination")

# Cached evaluator must reproduce direct limma::camera for fixed sets.
set.seed(718L)
y <- matrix(rnorm(600 * 30), nrow = 600, ncol = 30)
rownames(y) <- as.character(seq_len(nrow(y)))
group <- factor(rep(c("healthy", "RSV"), each = 15L), levels = c("healthy", "RSV"))
age <- rep(seq_len(15L), 2L)
sex <- factor(rep(c("female", "male"), 15L), levels = c("female", "male"))
design <- model.matrix(~ group + age + sex)
sets <- list(
  as.character(1:40),
  as.character(seq.int(51L, 110L)),
  as.character(seq.int(201L, 350L))
)
cache <- prepare_camera_cache(y, design, "groupRSV")
cached <- evaluate_camera_sets(cache, sets)
direct <- do.call(rbind, lapply(sets, function(ids) {
  z <- limma::camera(
    y = y, index = list(test = match(ids, rownames(y))), design = design,
    contrast = "groupRSV", inter.gene.cor = NA_real_, sort = FALSE,
    directional = TRUE
  )
  p <- z$PValue[[1L]]
  direction <- as.character(z$Direction[[1L]])
  data.frame(
    n_genes = z$NGenes[[1L]], correlation = z$Correlation[[1L]],
    direction = direction, p_value = p,
    signed_z = if (direction == "Up") qnorm(p / 2, lower.tail = FALSE) else -qnorm(p / 2, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}))
expect_equal(cached$n_genes, direct$n_genes, tolerance = 0, label = "cached camera set sizes differ")
expect_equal(cached$correlation, direct$correlation, tolerance = 1e-12, label = "cached camera correlations differ")
expect_true(identical(cached$direction, direct$direction), "cached camera directions differ")
expect_equal(cached$p_value, direct$p_value, tolerance = 1e-12, label = "cached camera P values differ")
expect_equal(cached$signed_z, direct$signed_z, tolerance = 1e-12, label = "cached camera signed Z differs")

fixture_frame <- data.frame(
  gene_id = rownames(y),
  mean_expression = rowMeans(y),
  detection_rate = seq(0, 1, length.out = nrow(y)),
  probe_count = rep(1:5, length.out = nrow(y)),
  is_target = seq_len(nrow(y)) <= 40L,
  stringsAsFactors = FALSE
)
fixture_sets <- make_deterministic_camera_fixtures(fixture_frame, set_size = 40L)
expect_true(length(fixture_sets) >= 5L, "fewer than five deterministic camera fixtures")
expect_true(all(lengths(fixture_sets) == 40L), "fixture set size differs")
expect_true(all(vapply(fixture_sets, function(x) !anyDuplicated(x), logical(1))), "fixture contains duplicate genes")
expect_true(all(vapply(fixture_sets, function(x) !length(intersect(x, fixture_frame$gene_id[fixture_frame$is_target])), logical(1))), "fixture contains target genes")

combined <- combine_signed_z(c(1, -0.5), c(sqrt(122), sqrt(73)))
expected <- (sqrt(122) - 0.5 * sqrt(73)) / sqrt(122 + 73)
expect_equal(combined, expected, tolerance = 1e-15, label = "weighted Stouffer differs")

plot_range <- null_plot_xlim(c(-0.9, -0.2, 0.4, 1.4), observed_z = 1.91)
expect_true(plot_range[[1L]] < -1.91 && plot_range[[2L]] > 1.91, "null plot xlim clips data or two-sided observed Z")
expect_true(diff(plot_range) > 0 && all(is.finite(plot_range)), "null plot xlim is invalid")

cat("Task 10 unit contracts: PASS\n")
