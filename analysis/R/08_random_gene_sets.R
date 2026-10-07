#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

task8_script <- file.path("analysis", "R", "06_test_target_set.R")
if (!file.exists(task8_script)) stop("Frozen Task 8 script is absent", call. = FALSE)
old_task8_skip <- Sys.getenv("TASK8_SKIP_MAIN", unset = NA_character_)
Sys.setenv(TASK8_SKIP_MAIN = "1")
source(task8_script, local = FALSE)
if (is.na(old_task8_skip)) Sys.unsetenv("TASK8_SKIP_MAIN") else Sys.setenv(TASK8_SKIP_MAIN = old_task8_skip)

task10_seed <- 20260804L
task10_iterations <- 10000L
task10_block_size <- 250L
task10_expected_detectable <- c(GSE105450 = 391L, GSE103842 = 374L)
task10_amendments <- c(
  file.path("docs", "science-superpowers", "preregistrations", "amendments", "2026-08-05-random-set-cross-cohort-operationalization.md"),
  file.path("docs", "science-superpowers", "preregistrations", "amendments", "2026-08-05-random-set-quintile-ties.md")
)

empirical_quintile <- function(x) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x))) {
    stop("Quintile input must be finite numeric", call. = FALSE)
  }
  out <- ceiling(5 * rank(x, ties.method = "average") / length(x))
  as.integer(pmax(1L, pmin(5L, out)))
}

stratum_key <- function(expression_bin, detection_bin, probe_bin) {
  paste(expression_bin, detection_bin, probe_bin, sep = ":")
}

append_allocation <- function(rows, target, source, n_allocate, stage) {
  if (n_allocate <= 0L) return(rows)
  rows[[length(rows) + 1L]] <- data.frame(
    target_expression_bin = target[["expression_bin"]],
    target_detection_bin = target[["detection_bin"]],
    target_probe_bin = target[["probe_bin"]],
    source_expression_bin = source[["expression_bin"]],
    source_detection_bin = source[["detection_bin"]],
    source_probe_bin = source[["probe_bin"]],
    n_allocate = as.integer(n_allocate),
    relaxation_stage = stage,
    stringsAsFactors = FALSE
  )
  rows
}

build_relaxation_plan <- function(matching_frame) {
  required <- c("gene_id", "expression_bin", "detection_bin", "probe_bin", "is_target")
  if (!all(required %in% names(matching_frame))) stop("Matching frame is incomplete", call. = FALSE)
  if (anyDuplicated(matching_frame$gene_id) || anyNA(matching_frame$gene_id) ||
      any(!nzchar(as.character(matching_frame$gene_id)))) {
    stop("Matching frame gene IDs must be unique and non-empty", call. = FALSE)
  }
  bins <- matching_frame[c("expression_bin", "detection_bin", "probe_bin")]
  if (any(!vapply(bins, function(x) all(x %in% 1:5), logical(1)))) {
    stop("Matching bins must be integers in 1..5", call. = FALSE)
  }
  key <- stratum_key(bins$expression_bin, bins$detection_bin, bins$probe_bin)
  target_tab <- aggregate(
    rep(1L, sum(matching_frame$is_target)),
    by = matching_frame[matching_frame$is_target, c("expression_bin", "detection_bin", "probe_bin"), drop = FALSE],
    FUN = sum
  )
  names(target_tab)[4L] <- "remaining"
  source_tab <- aggregate(
    rep(1L, sum(!matching_frame$is_target)),
    by = matching_frame[!matching_frame$is_target, c("expression_bin", "detection_bin", "probe_bin"), drop = FALSE],
    FUN = sum
  )
  names(source_tab)[4L] <- "remaining"
  target_tab$key <- stratum_key(target_tab$expression_bin, target_tab$detection_bin, target_tab$probe_bin)
  source_tab$key <- stratum_key(source_tab$expression_bin, source_tab$detection_bin, source_tab$probe_bin)
  target_tab <- target_tab[order(target_tab$expression_bin, target_tab$detection_bin, target_tab$probe_bin), ]
  source_tab <- source_tab[order(source_tab$expression_bin, source_tab$detection_bin, source_tab$probe_bin), ]
  rows <- list()

  # Exact matching has disjoint source strata and therefore cannot create competition.
  for (ti in seq_len(nrow(target_tab))) {
    si <- match(target_tab$key[[ti]], source_tab$key)
    if (!is.na(si)) {
      take <- min(target_tab$remaining[[ti]], source_tab$remaining[[si]])
      rows <- append_allocation(rows, target_tab[ti, ], source_tab[si, ], take, "exact")
      target_tab$remaining[[ti]] <- target_tab$remaining[[ti]] - take
      source_tab$remaining[[si]] <- source_tab$remaining[[si]] - take
    }
  }

  allocate_pairs <- function(target_indices, source_indices, stage, pair_order) {
    if (!length(target_indices) || !length(source_indices)) return(invisible(NULL))
    pairs <- expand.grid(ti = target_indices, si = source_indices, KEEP.OUT.ATTRS = FALSE)
    pairs <- pairs[target_tab$remaining[pairs$ti] > 0L & source_tab$remaining[pairs$si] > 0L, , drop = FALSE]
    if (!nrow(pairs)) return(invisible(NULL))
    order_columns <- pair_order(pairs)
    ord <- do.call(order, c(order_columns, list(method = "radix")))
    pairs <- pairs[ord, , drop = FALSE]
    for (k in seq_len(nrow(pairs))) {
      ti <- pairs$ti[[k]]
      si <- pairs$si[[k]]
      take <- min(target_tab$remaining[[ti]], source_tab$remaining[[si]])
      if (take > 0L) {
        rows <<- append_allocation(rows, target_tab[ti, ], source_tab[si, ], take, stage)
        target_tab$remaining[[ti]] <<- target_tab$remaining[[ti]] - take
        source_tab$remaining[[si]] <<- source_tab$remaining[[si]] - take
      }
    }
    invisible(NULL)
  }

  # First relaxation: probe-count bin only, expression and detection fixed.
  for (expression_bin in 1:5) for (detection_bin in 1:5) {
    ti <- which(target_tab$remaining > 0L & target_tab$expression_bin == expression_bin &
                  target_tab$detection_bin == detection_bin)
    si <- which(source_tab$remaining > 0L & source_tab$expression_bin == expression_bin &
                  source_tab$detection_bin == detection_bin)
    allocate_pairs(ti, si, "probe_only", function(pairs) {
      distance <- abs(target_tab$probe_bin[pairs$ti] - source_tab$probe_bin[pairs$si])
      list(distance, source_tab$probe_bin[pairs$si], target_tab$probe_bin[pairs$ti])
    })
  }

  # Second relaxation: adjacent detection bins progressively; expression remains exact.
  for (expression_bin in 1:5) {
    ti <- which(target_tab$remaining > 0L & target_tab$expression_bin == expression_bin)
    si <- which(source_tab$remaining > 0L & source_tab$expression_bin == expression_bin)
    allocate_pairs(ti, si, "adjacent_detection", function(pairs) {
      detection_distance <- abs(target_tab$detection_bin[pairs$ti] - source_tab$detection_bin[pairs$si])
      probe_distance <- abs(target_tab$probe_bin[pairs$ti] - source_tab$probe_bin[pairs$si])
      list(
        detection_distance, probe_distance,
        source_tab$detection_bin[pairs$si], source_tab$probe_bin[pairs$si],
        target_tab$detection_bin[pairs$ti], target_tab$probe_bin[pairs$ti]
      )
    })
  }

  if (any(target_tab$remaining > 0L)) {
    failed <- target_tab[target_tab$remaining > 0L, c("key", "remaining")]
    stop(
      "Matching is infeasible without relaxing expression bin: ",
      paste(paste0(failed$key, " n=", failed$remaining), collapse = "; "),
      call. = FALSE
    )
  }
  if (!length(rows)) stop("No allocations were produced", call. = FALSE)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$target_key <- stratum_key(out$target_expression_bin, out$target_detection_bin, out$target_probe_bin)
  out$source_key <- stratum_key(out$source_expression_bin, out$source_detection_bin, out$source_probe_bin)
  if (sum(out$n_allocate) != sum(matching_frame$is_target)) stop("Allocation total mismatch", call. = FALSE)
  if (any(out$target_expression_bin != out$source_expression_bin)) stop("Expression bin was relaxed", call. = FALSE)
  source_capacity <- table(key[!matching_frame$is_target])
  source_used <- tapply(out$n_allocate, out$source_key, sum)
  if (any(source_used > source_capacity[names(source_used)])) stop("Allocation exceeds source capacity", call. = FALSE)
  out
}

draw_one_matched_set <- function(matching_frame, plan) {
  candidates <- split(as.character(matching_frame$gene_id[!matching_frame$is_target]),
                      stratum_key(
                        matching_frame$expression_bin[!matching_frame$is_target],
                        matching_frame$detection_bin[!matching_frame$is_target],
                        matching_frame$probe_bin[!matching_frame$is_target]
                      ))
  selected <- character()
  for (source_key in unique(plan$source_key)) {
    rows <- which(plan$source_key == source_key)
    total <- sum(plan$n_allocate[rows])
    pool <- candidates[[source_key]]
    if (length(pool) < total) stop("Source stratum depleted: ", source_key, call. = FALSE)
    sampled <- sample(pool, size = total, replace = FALSE)
    selected <- c(selected, sampled)
  }
  target_ids <- as.character(matching_frame$gene_id[matching_frame$is_target])
  if (length(selected) != length(target_ids) || anyDuplicated(selected) || length(intersect(selected, target_ids))) {
    stop("Invalid matched random set", call. = FALSE)
  }
  selected
}

draw_matched_sets <- function(matching_frame, plan, n_sets) {
  if (length(n_sets) != 1L || !is.finite(n_sets) || n_sets < 1 || n_sets != as.integer(n_sets)) {
    stop("n_sets must be a positive integer", call. = FALSE)
  }
  replicate(as.integer(n_sets), draw_one_matched_set(matching_frame, plan), simplify = FALSE)
}

prepare_camera_cache <- function(expression, design, coefficient) {
  y <- as.matrix(expression)
  storage.mode(y) <- "double"
  design <- as.matrix(design)
  storage.mode(design) <- "double"
  if (any(!is.finite(y)) || any(!is.finite(design)) || nrow(design) != ncol(y)) {
    stop("Invalid camera cache input", call. = FALSE)
  }
  if (is.character(coefficient)) coefficient <- match(coefficient, colnames(design))
  if (length(coefficient) != 1L || is.na(coefficient)) stop("camera coefficient is absent", call. = FALSE)
  p <- ncol(design)
  G <- nrow(y)
  n <- ncol(y)
  df_residual <- n - p
  if (df_residual < 1L || G < 3L || qr(design)$rank != p) stop("Invalid camera dimensions/design", call. = FALSE)
  j <- c((seq_len(p))[-coefficient], coefficient)
  if (coefficient < p) design <- design[, j, drop = FALSE]
  QR <- qr(design)
  effects <- qr.qty(QR, t(y))
  unscaled_t <- effects[p, ]
  if (QR$qr[p, p] < 0) unscaled_t <- -unscaled_t
  U <- effects[-seq_len(p), , drop = FALSE]
  sigma2 <- colMeans(U^2)
  U <- t(U) / sqrt(pmax(sigma2, 1e-8))
  sv <- limma::squeezeVar(sigma2, df = df_residual, covariate = NULL)
  mod_t <- unscaled_t / sqrt(sv$var.post)
  df_total <- min(df_residual + sv$df.prior, G * df_residual)
  Stat <- limma::zscoreT(mod_t, df = df_total, approx = TRUE, method = "hill")
  list(
    gene_ids = rownames(y), G = G, df_camera = min(df_residual, G - 2L),
    U = U, residual_dimensions = ncol(U), Stat = Stat,
    mean_stat = mean(Stat), var_stat = stats::var(Stat)
  )
}

evaluate_camera_sets <- function(cache, sets) {
  if (!length(sets)) stop("camera set list is empty", call. = FALSE)
  sizes <- lengths(sets)
  if (any(sizes < 1L) || any(sizes >= cache$G)) stop("Invalid camera set size", call. = FALSE)
  indices <- lapply(sets, function(ids) {
    idx <- match(as.character(ids), cache$gene_ids)
    if (anyNA(idx) || anyDuplicated(idx)) stop("camera set has absent or duplicate genes", call. = FALSE)
    idx
  })
  membership <- Matrix::sparseMatrix(
    i = rep.int(seq_along(indices), sizes),
    j = unlist(indices, use.names = FALSE), x = 1,
    dims = c(length(indices), cache$G)
  )
  mean_in <- as.numeric(membership %*% cache$Stat) / sizes
  sum_u <- as.matrix(membership %*% cache$U)
  raw_vif <- rowSums(sum_u^2) / (sizes * cache$residual_dimensions)
  correlation <- (raw_vif - 1) / (sizes - 1)
  vif <- pmax(1, raw_vif)
  m2 <- cache$G - sizes
  delta <- cache$G / m2 * (mean_in - cache$mean_stat)
  pooled <- ((cache$G - 1) * cache$var_stat - delta^2 * sizes * m2 / cache$G) / (cache$G - 2)
  statistic <- delta / sqrt(pooled * (vif / sizes + 1 / m2))
  down <- stats::pt(statistic, df = cache$df_camera)
  up <- stats::pt(statistic, df = cache$df_camera, lower.tail = FALSE)
  p <- 2 * pmin(down, up)
  direction <- ifelse(down < up, "Down", "Up")
  signed_z <- ifelse(direction == "Up", qnorm(p / 2, lower.tail = FALSE), -qnorm(p / 2, lower.tail = FALSE))
  data.frame(
    n_genes = as.integer(sizes), correlation = correlation,
    direction = direction, p_value = p, signed_z = signed_z,
    stringsAsFactors = FALSE
  )
}

combine_signed_z <- function(signed_z, weights) {
  if (length(signed_z) != length(weights) || !length(signed_z) ||
      any(!is.finite(signed_z)) || any(!is.finite(weights) | weights <= 0)) {
    stop("Invalid Stouffer input", call. = FALSE)
  }
  sum(weights * signed_z) / sqrt(sum(weights^2))
}

make_deterministic_camera_fixtures <- function(matching_frame, set_size) {
  required <- c("gene_id", "mean_expression", "detection_rate", "probe_count", "is_target")
  if (!all(required %in% names(matching_frame))) stop("Fixture frame is incomplete", call. = FALSE)
  candidates <- matching_frame[!matching_frame$is_target, , drop = FALSE]
  if (nrow(candidates) < set_size || set_size < 2L) stop("Insufficient fixture candidates", call. = FALSE)
  take <- function(order_index) as.character(candidates$gene_id[order_index[seq_len(set_size)]])
  lexicographic <- order(as.character(candidates$gene_id), method = "radix")
  list(
    lexicographic_first = take(lexicographic),
    lexicographic_last = take(rev(lexicographic)),
    low_mean_expression = take(order(candidates$mean_expression, candidates$gene_id, method = "radix")),
    high_mean_expression = take(order(-candidates$mean_expression, candidates$gene_id, method = "radix")),
    low_detection_probe_order = take(order(candidates$detection_rate, candidates$probe_count, candidates$gene_id, method = "radix")),
    high_detection_probe_order = take(order(-candidates$detection_rate, -candidates$probe_count, candidates$gene_id, method = "radix"))
  )
}

build_matching_frame <- function(cohort, model, targets) {
  object <- readRDS(file.path("data", "derived", paste0(cohort, "_expression.rds")))
  feature <- object$feature_data
  required <- c("entrez_id", "all_sample_mean_log2", "detection_rate_p_lt_0_05", "detectable_blind")
  if (!all(required %in% names(feature))) stop(cohort, ": feature_data lacks matching fields", call. = FALSE)
  feature$entrez_id <- as.character(feature$entrez_id)
  feature <- feature[match(rownames(model$expression), feature$entrez_id), , drop = FALSE]
  if (anyNA(feature$entrez_id) || !identical(feature$entrez_id, rownames(model$expression)) ||
      any(!feature$detectable_blind)) {
    stop(cohort, ": detectable universe does not match feature_data", call. = FALSE)
  }
  mapping_path <- file.path("results", "tables", paste0(cohort, "_probe_entrez_mapping_audit.tsv"))
  mapping <- read_tsv(mapping_path)
  mapped <- mapping[mapping$mapping_status == "unique_current" & grepl("^[0-9]+$", mapping$entrez_id), , drop = FALSE]
  probe_counts <- table(as.character(mapped$entrez_id))
  probe_count <- as.integer(probe_counts[feature$entrez_id])
  if (anyNA(probe_count) || any(probe_count < 1L)) stop(cohort, ": invalid unique-current probe count", call. = FALSE)
  out <- data.frame(
    gene_id = feature$entrez_id,
    mean_expression = as.numeric(feature$all_sample_mean_log2),
    detection_rate = as.numeric(feature$detection_rate_p_lt_0_05),
    probe_count = probe_count,
    is_target = feature$entrez_id %in% targets,
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(out$mean_expression)) || any(!is.finite(out$detection_rate)) ||
      any(out$detection_rate < 0 | out$detection_rate > 1) || anyDuplicated(out$gene_id)) {
    stop(cohort, ": invalid blind matching frame", call. = FALSE)
  }
  out$expression_bin <- empirical_quintile(out$mean_expression)
  out$detection_bin <- empirical_quintile(out$detection_rate)
  out$probe_bin <- empirical_quintile(out$probe_count)
  if (sum(out$is_target) != unname(task10_expected_detectable[[cohort]])) {
    stop(cohort, ": detectable target count differs from frozen Task 8", call. = FALSE)
  }
  out
}

validate_observed_task8 <- function(models, targets) {
  saved_camera <- read_tsv(file.path("results", "cohort", "target_set_camera.tsv"))
  saved_decision <- read_tsv(file.path("results", "cohort", "target_set_primary_decision_partial.tsv"))
  if (!identical(as.character(saved_camera$cohort), task8_cohorts) || nrow(saved_decision) != 1L) {
    stop("Frozen Task 8 outputs have invalid scope", call. = FALSE)
  }
  reconstructed <- do.call(rbind, lapply(task8_cohorts, function(cohort) {
    model <- models[[cohort]]
    run_camera_prepared(model$expression, model$design, model$coefficient, targets, cohort)
  }))
  numeric_fields <- c("p_value", "signed_z", "inter_gene_correlation")
  if (!identical(as.character(reconstructed$direction), as.character(saved_camera$direction)) ||
      !all(vapply(numeric_fields, function(field) {
        isTRUE(all.equal(reconstructed[[field]], saved_camera[[field]], tolerance = 1e-13))
      }, logical(1)))) {
    stop("Task 8 observed camera result does not reconstruct", call. = FALSE)
  }
  stouffer <- combine_camera_stouffer(reconstructed)
  if (!isTRUE(all.equal(stouffer$combined_z, saved_decision$weighted_stouffer_z[[1L]], tolerance = 1e-13))) {
    stop("Task 8 observed Stouffer Z does not reconstruct", call. = FALSE)
  }
  list(camera = reconstructed, combined_z = stouffer$combined_z, decision = saved_decision)
}

validate_cached_camera_real <- function(models, matching, observed) {
  rows <- lapply(task8_cohorts, function(cohort) {
    model <- models[[cohort]]
    cache <- prepare_camera_cache(model$expression, model$design, model$coefficient)
    target_ids <- matching[[cohort]]$gene_id[matching[[cohort]]$is_target]
    cached <- evaluate_camera_sets(cache, list(target_ids))
    direct <- observed$camera[observed$camera$cohort == cohort, ]
    numeric_ok <-
      isTRUE(all.equal(cached$p_value, direct$p_value, tolerance = 1e-12)) &&
      isTRUE(all.equal(cached$signed_z, direct$signed_z, tolerance = 1e-12)) &&
      isTRUE(all.equal(cached$correlation, direct$inter_gene_correlation, tolerance = 1e-12))
    if (!numeric_ok || !identical(cached$direction, as.character(direct$direction))) {
      stop(cohort, ": cached evaluator differs from direct limma::camera", call. = FALSE)
    }
    observed_row <- data.frame(
      cohort = cohort, fixture = "observed_target_set",
      direct_p = direct$p_value, cached_p = cached$p_value,
      abs_p_difference = abs(direct$p_value - cached$p_value),
      direct_correlation = direct$inter_gene_correlation,
      cached_correlation = cached$correlation,
      abs_correlation_difference = abs(direct$inter_gene_correlation - cached$correlation),
      direct_signed_z = direct$signed_z, cached_signed_z = cached$signed_z,
      abs_signed_z_difference = abs(direct$signed_z - cached$signed_z),
      direction_match = identical(cached$direction, as.character(direct$direction)),
      tolerance = 1e-12, status = "passed_exact_camera_contract",
      stringsAsFactors = FALSE
    )
    fixtures <- make_deterministic_camera_fixtures(matching[[cohort]], cached$n_genes)
    cached_fixture <- evaluate_camera_sets(cache, fixtures)
    direct_fixture <- do.call(rbind, lapply(fixtures, function(ids) {
      result <- limma::camera(
        y = model$expression, index = list(fixture = match(ids, rownames(model$expression))),
        design = model$design, contrast = model$coefficient, inter.gene.cor = NA_real_,
        sort = FALSE, directional = TRUE
      )
      p <- result$PValue[[1L]]
      direction <- as.character(result$Direction[[1L]])
      data.frame(
        direct_p = p, direct_correlation = result$Correlation[[1L]],
        direct_signed_z = if (direction == "Up") qnorm(p / 2, lower.tail = FALSE) else -qnorm(p / 2, lower.tail = FALSE),
        direct_direction = direction, stringsAsFactors = FALSE
      )
    }))
    fixture_rows <- data.frame(
      cohort = cohort, fixture = names(fixtures),
      direct_p = direct_fixture$direct_p, cached_p = cached_fixture$p_value,
      abs_p_difference = abs(direct_fixture$direct_p - cached_fixture$p_value),
      direct_correlation = direct_fixture$direct_correlation,
      cached_correlation = cached_fixture$correlation,
      abs_correlation_difference = abs(direct_fixture$direct_correlation - cached_fixture$correlation),
      direct_signed_z = direct_fixture$direct_signed_z,
      cached_signed_z = cached_fixture$signed_z,
      abs_signed_z_difference = abs(direct_fixture$direct_signed_z - cached_fixture$signed_z),
      direction_match = direct_fixture$direct_direction == cached_fixture$direction,
      tolerance = 1e-12, status = "passed_exact_camera_contract",
      stringsAsFactors = FALSE
    )
    if (any(fixture_rows$abs_p_difference > 1e-12) ||
        any(fixture_rows$abs_correlation_difference > 1e-12) ||
        any(fixture_rows$abs_signed_z_difference > 1e-12) || any(!fixture_rows$direction_match)) {
      stop(cohort, ": deterministic cached/direct camera fixture failed", call. = FALSE)
    }
    rbind(observed_row, fixture_rows)
  })
  do.call(rbind, rows)
}

null_plot_xlim <- function(random_z, observed_z, padding_fraction = 0.04) {
  if (!is.numeric(random_z) || !length(random_z) || any(!is.finite(random_z)) ||
      length(observed_z) != 1L || !is.finite(observed_z) ||
      length(padding_fraction) != 1L || !is.finite(padding_fraction) || padding_fraction <= 0) {
    stop("Invalid null-plot range input", call. = FALSE)
  }
  limits <- range(c(random_z, observed_z, -abs(observed_z)))
  span <- diff(limits)
  if (span <= 0) span <- max(1, abs(limits[[1L]]))
  limits + c(-1, 1) * padding_fraction * span
}

write_null_figure <- function(null, observed_z, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = 7.2, height = 5.4, useDingbats = FALSE, timestamp = FALSE, bg = "white")
  old <- graphics::par(mar = c(4.5, 4.5, 2.2, 1))
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  graphics::hist(
    null$combined_stouffer_z, breaks = 60, col = "#B8D8E8", border = "white",
    main = "Matched random-set null distribution", xlab = "Weighted directional Stouffer Z",
    ylab = "Random-set count", xlim = null_plot_xlim(null$combined_stouffer_z, observed_z)
  )
  graphics::abline(v = observed_z, col = "#B2182B", lwd = 2.2)
  graphics::abline(v = -abs(observed_z), col = "#B2182B", lwd = 1.2, lty = 2)
  graphics::legend(
    "topright", legend = c("Observed Z", "Two-sided |observed Z|"),
    col = "#B2182B", lwd = c(2.2, 1.2), lty = c(1, 2), bty = "n"
  )
}

verify_formal_random_sets <- function(null, matching, plans, caches, models, observed_z) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Frozen digest package is required", call. = FALSE)
  top_iterations <- null$iteration[order(abs(null$combined_stouffer_z), decreasing = TRUE)[seq_len(10L)]]
  retained <- setNames(lapply(task8_cohorts, function(x) vector("list", length(top_iterations))), task8_cohorts)
  hashes <- setNames(lapply(task8_cohorts, function(x) character(task10_iterations)), task8_cohorts)
  pair_hash <- character(task10_iterations)
  stratum_mismatch <- setNames(logical(length(task8_cohorts)), task8_cohorts)
  target_strata <- lapply(matching, function(frame) {
    table(stratum_key(
      frame$expression_bin[frame$is_target], frame$detection_bin[frame$is_target], frame$probe_bin[frame$is_target]
    ))
  })
  set.seed(task10_seed)
  for (iteration in seq_len(task10_iterations)) {
    current_hashes <- character(length(task8_cohorts))
    names(current_hashes) <- task8_cohorts
    for (cohort in task8_cohorts) {
      ids <- draw_one_matched_set(matching[[cohort]], plans[[cohort]])
      frame <- matching[[cohort]]
      drawn_rows <- match(ids, frame$gene_id)
      observed_strata <- table(stratum_key(
        frame$expression_bin[drawn_rows], frame$detection_bin[drawn_rows], frame$probe_bin[drawn_rows]
      ))
      all_keys <- union(names(target_strata[[cohort]]), names(observed_strata))
      expected_counts <- as.integer(target_strata[[cohort]][all_keys]); expected_counts[is.na(expected_counts)] <- 0L
      observed_counts <- as.integer(observed_strata[all_keys]); observed_counts[is.na(observed_counts)] <- 0L
      stratum_mismatch[[cohort]] <- stratum_mismatch[[cohort]] || !identical(expected_counts, observed_counts)
      hash <- digest::digest(paste(sort(ids), collapse = ","), algo = "sha256", serialize = FALSE)
      hashes[[cohort]][[iteration]] <- hash
      current_hashes[[cohort]] <- hash
      position <- match(iteration, top_iterations)
      if (!is.na(position)) retained[[cohort]][[position]] <- ids
    }
    pair_hash[[iteration]] <- digest::digest(paste(current_hashes, collapse = ":"), algo = "sha256", serialize = FALSE)
  }
  if (any(stratum_mismatch) || any(vapply(retained, function(x) any(lengths(x) == 0L), logical(1)))) {
    stop("Formal random-set replay failed strata or retention checks", call. = FALSE)
  }

  extreme_rows <- list()
  row_index <- 0L
  direct_z <- matrix(NA_real_, nrow = length(top_iterations), ncol = length(task8_cohorts),
                     dimnames = list(as.character(top_iterations), task8_cohorts))
  for (cohort in task8_cohorts) {
    cached <- evaluate_camera_sets(caches[[cohort]], retained[[cohort]])
    for (j in seq_along(top_iterations)) {
      iteration <- top_iterations[[j]]
      ids <- retained[[cohort]][[j]]
      direct <- limma::camera(
        y = models[[cohort]]$expression,
        index = list(extreme_random_set = match(ids, rownames(models[[cohort]]$expression))),
        design = models[[cohort]]$design, contrast = models[[cohort]]$coefficient,
        inter.gene.cor = NA_real_, sort = FALSE, directional = TRUE
      )
      p <- direct$PValue[[1L]]
      direction <- as.character(direct$Direction[[1L]])
      signed_z <- if (direction == "Up") qnorm(p / 2, lower.tail = FALSE) else -qnorm(p / 2, lower.tail = FALSE)
      direct_z[as.character(iteration), cohort] <- signed_z
      saved_p <- null[[paste0(cohort, "_p_value")]][[iteration]]
      saved_z <- null[[paste0(cohort, "_signed_z")]][[iteration]]
      saved_direction <- null[[paste0(cohort, "_direction")]][[iteration]]
      row_index <- row_index + 1L
      extreme_rows[[row_index]] <- data.frame(
        iteration = iteration, cohort = cohort, set_sha256 = hashes[[cohort]][[iteration]],
        direct_p = p, cached_p = cached$p_value[[j]], saved_p = saved_p,
        abs_p_difference = max(abs(p - cached$p_value[[j]]), abs(p - saved_p)),
        direct_correlation = direct$Correlation[[1L]], cached_correlation = cached$correlation[[j]],
        abs_correlation_difference = abs(direct$Correlation[[1L]] - cached$correlation[[j]]),
        direct_signed_z = signed_z, cached_signed_z = cached$signed_z[[j]], saved_signed_z = saved_z,
        abs_signed_z_difference = max(abs(signed_z - cached$signed_z[[j]]), abs(signed_z - saved_z)),
        direct_direction = direction, cached_direction = cached$direction[[j]], saved_direction = saved_direction,
        direction_match = direction == cached$direction[[j]] && direction == saved_direction,
        stringsAsFactors = FALSE
      )
    }
  }
  extreme <- do.call(rbind, extreme_rows)
  direct_combined <- vapply(top_iterations, function(iteration) {
    combine_signed_z(direct_z[as.character(iteration), ], sqrt(unname(task8_independent_n[task8_cohorts])))
  }, numeric(1))
  saved_combined <- null$combined_stouffer_z[top_iterations]
  combined_map <- setNames(direct_combined, top_iterations)
  saved_map <- setNames(saved_combined, top_iterations)
  extreme$direct_combined_z <- combined_map[as.character(extreme$iteration)]
  extreme$saved_combined_z <- saved_map[as.character(extreme$iteration)]
  extreme$abs_combined_difference <- abs(extreme$direct_combined_z - extreme$saved_combined_z)
  if (any(extreme$abs_p_difference > 1e-12) || any(extreme$abs_correlation_difference > 1e-12) ||
      any(extreme$abs_signed_z_difference > 1e-12) || any(extreme$abs_combined_difference > 1e-12) ||
      any(!extreme$direction_match)) {
    stop("Direct limma::camera verification of extreme random sets failed", call. = FALSE)
  }

  weights <- sqrt(unname(task8_independent_n[task8_cohorts]))
  manual_combined <- (weights[[1L]] * null$GSE105450_signed_z + weights[[2L]] * null$GSE103842_signed_z) /
    sqrt(sum(weights^2))
  ks_105450 <- suppressWarnings(stats::ks.test(null$GSE105450_p_value, "punif"))
  ks_103842 <- suppressWarnings(stats::ks.test(null$GSE103842_p_value, "punif"))
  q <- stats::quantile(null$combined_stouffer_z, c(0.001, 0.01, 0.025, 0.5, 0.975, 0.99, 0.999), names = FALSE)
  diagnostics <- data.frame(
    iterations = nrow(null), observed_combined_z = observed_z,
    combined_z_min = min(null$combined_stouffer_z), combined_z_max = max(null$combined_stouffer_z),
    max_abs_combined_z = max(abs(null$combined_stouffer_z)),
    combined_z_mean = mean(null$combined_stouffer_z), combined_z_sd = stats::sd(null$combined_stouffer_z),
    combined_z_q001 = q[[1L]], combined_z_q01 = q[[2L]], combined_z_q025 = q[[3L]],
    combined_z_median = q[[4L]], combined_z_q975 = q[[5L]], combined_z_q99 = q[[6L]], combined_z_q999 = q[[7L]],
    GSE105450_signed_z_min = min(null$GSE105450_signed_z), GSE105450_signed_z_max = max(null$GSE105450_signed_z),
    GSE105450_signed_z_mean = mean(null$GSE105450_signed_z), GSE105450_signed_z_sd = stats::sd(null$GSE105450_signed_z),
    GSE103842_signed_z_min = min(null$GSE103842_signed_z), GSE103842_signed_z_max = max(null$GSE103842_signed_z),
    GSE103842_signed_z_mean = mean(null$GSE103842_signed_z), GSE103842_signed_z_sd = stats::sd(null$GSE103842_signed_z),
    cohort_signed_z_correlation = stats::cor(null$GSE105450_signed_z, null$GSE103842_signed_z),
    cohort_signed_z_same_sign_fraction = mean(sign(null$GSE105450_signed_z) == sign(null$GSE103842_signed_z)),
    GSE105450_p_min = min(null$GSE105450_p_value), GSE105450_p_mean = mean(null$GSE105450_p_value),
    GSE105450_p_median = stats::median(null$GSE105450_p_value), GSE105450_p_lt_0_05_fraction = mean(null$GSE105450_p_value < 0.05),
    GSE105450_uniform_KS_D = unname(ks_105450$statistic), GSE105450_uniform_KS_p = ks_105450$p.value,
    GSE103842_p_min = min(null$GSE103842_p_value), GSE103842_p_mean = mean(null$GSE103842_p_value),
    GSE103842_p_median = stats::median(null$GSE103842_p_value), GSE103842_p_lt_0_05_fraction = mean(null$GSE103842_p_value < 0.05),
    GSE103842_uniform_KS_D = unname(ks_103842$statistic), GSE103842_uniform_KS_p = ks_103842$p.value,
    random_camera_p_expected_uniform = FALSE,
    random_camera_p_reference = "conditional_gene_membership_resampling_in_one_fixed_experiment_not_repeated_experiment_null",
    max_abs_stouffer_formula_difference = max(abs(manual_combined - null$combined_stouffer_z)),
    GSE105450_unique_set_hashes = length(unique(hashes$GSE105450)),
    GSE103842_unique_set_hashes = length(unique(hashes$GSE103842)),
    unique_paired_set_hashes = length(unique(pair_hash)),
    any_stratum_count_mismatch = any(stratum_mismatch),
    extreme_iterations_direct_camera_checked = length(top_iterations),
    extreme_direct_camera_tolerance = 1e-12,
    stringsAsFactors = FALSE
  )
  list(diagnostics = diagnostics, extreme = extreme)
}

main_task10 <- function() {
  started <- Sys.time()
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Frozen Matrix package is required", call. = FALSE)
  if (!all(file.exists(task10_amendments))) stop("Task 10 pre-draw amendments are absent", call. = FALSE)
  frozen <- verify_task8_inputs()
  targets <- read_frozen_targets()
  manifest <- read_tsv(file.path("data", "clean", "sample_manifest_frozen.tsv"))
  models <- lapply(task8_cohorts, verify_reconstructed_model, manifest = manifest, frozen = frozen)
  names(models) <- task8_cohorts
  observed <- validate_observed_task8(models, targets)
  matching <- lapply(task8_cohorts, function(cohort) build_matching_frame(cohort, models[[cohort]], targets))
  names(matching) <- task8_cohorts
  plans <- lapply(matching, build_relaxation_plan)
  caches <- lapply(task8_cohorts, function(cohort) {
    prepare_camera_cache(models[[cohort]]$expression, models[[cohort]]$design, models[[cohort]]$coefficient)
  })
  names(caches) <- task8_cohorts
  camera_contract <- validate_cached_camera_real(models, matching, observed)

  # This is the first random draw in Task 10.
  set.seed(task10_seed)
  null_parts <- vector("list", ceiling(task10_iterations / task10_block_size))
  all_sizes_valid <- TRUE
  any_duplicate <- FALSE
  any_target_contamination <- FALSE
  part <- 0L
  for (block_start in seq.int(1L, task10_iterations, by = task10_block_size)) {
    part <- part + 1L
    block_end <- min(task10_iterations, block_start + task10_block_size - 1L)
    n_block <- block_end - block_start + 1L
    sets <- setNames(vector("list", length(task8_cohorts)), task8_cohorts)
    sets[[task8_cohorts[[1L]]]] <- vector("list", n_block)
    sets[[task8_cohorts[[2L]]]] <- vector("list", n_block)
    for (j in seq_len(n_block)) {
      for (cohort in task8_cohorts) {
        ids <- draw_one_matched_set(matching[[cohort]], plans[[cohort]])
        sets[[cohort]][[j]] <- ids
        all_sizes_valid <- all_sizes_valid && length(ids) == unname(task10_expected_detectable[[cohort]])
        any_duplicate <- any_duplicate || anyDuplicated(ids) > 0L
        any_target_contamination <- any_target_contamination ||
          length(intersect(ids, matching[[cohort]]$gene_id[matching[[cohort]]$is_target])) > 0L
      }
    }
    camera <- lapply(task8_cohorts, function(cohort) evaluate_camera_sets(caches[[cohort]], sets[[cohort]]))
    names(camera) <- task8_cohorts
    combined <- vapply(seq_len(n_block), function(j) {
      combine_signed_z(
        c(camera[[task8_cohorts[[1L]]]]$signed_z[[j]], camera[[task8_cohorts[[2L]]]]$signed_z[[j]]),
        sqrt(unname(task8_independent_n[task8_cohorts]))
      )
    }, numeric(1))
    null_parts[[part]] <- data.frame(
      iteration = block_start:block_end,
      GSE105450_direction = camera$GSE105450$direction,
      GSE105450_p_value = camera$GSE105450$p_value,
      GSE105450_signed_z = camera$GSE105450$signed_z,
      GSE103842_direction = camera$GSE103842$direction,
      GSE103842_p_value = camera$GSE103842$p_value,
      GSE103842_signed_z = camera$GSE103842$signed_z,
      combined_stouffer_z = combined,
      stringsAsFactors = FALSE
    )
    if (block_end %in% c(2500L, 5000L, 7500L, 10000L)) {
      cat("Task 10 progress: ", block_end, "/", task10_iterations, " iterations completed\n", sep = "")
      flush.console()
    }
  }
  null <- do.call(rbind, null_parts)
  if (nrow(null) != task10_iterations || !identical(null$iteration, seq_len(task10_iterations)) ||
      !all_sizes_valid || any_duplicate || any_target_contamination || any(!is.finite(null$combined_stouffer_z))) {
    stop("10,000-set validation failed", call. = FALSE)
  }

  observed_z <- observed$combined_z
  extreme_count <- sum(abs(null$combined_stouffer_z) >= abs(observed_z))
  empirical_p <- (1 + extreme_count) / (task10_iterations + 1)
  signed_percentile <- 100 * mean(null$combined_stouffer_z <= observed_z)
  absolute_percentile <- 100 * mean(abs(null$combined_stouffer_z) <= abs(observed_z))
  calibration <- data.frame(
    analysis_scope = "GSE105450_and_GSE103842_confirmatory_only",
    random_seed = task10_seed, random_iterations = task10_iterations,
    observed_statistic = "Task8_weighted_directional_Stouffer_Z",
    observed_combined_stouffer_z = observed_z,
    random_extreme_abs_count = extreme_count,
    empirical_two_sided_p = empirical_p,
    observed_signed_percentile = signed_percentile,
    observed_absolute_extremeness_percentile = absolute_percentile,
    empirical_p_lt_0_05 = empirical_p < 0.05,
    prior_camera_component_met = FALSE,
    primary_decision = "H1_not_supported_camera_component_failed",
    h1_supported = FALSE,
    random_calibration_can_rescue = FALSE,
    matching_identity = "independent_cohort_sets_paired_by_iteration",
    GSE105450_set_size = task10_expected_detectable[["GSE105450"]],
    GSE103842_set_size = task10_expected_detectable[["GSE103842"]],
    stouffer_weight_GSE105450 = sqrt(task8_independent_n[["GSE105450"]]),
    stouffer_weight_GSE103842 = sqrt(task8_independent_n[["GSE103842"]]),
    empirical_formula = "(1+sum(abs(T_random)>=abs(T_observed)))/(10000+1)",
    stringsAsFactors = FALSE
  )
  formal_verification <- verify_formal_random_sets(null, matching, plans, caches, models, observed_z)

  plan_output <- do.call(rbind, lapply(task8_cohorts, function(cohort) {
    x <- plans[[cohort]]
    x$cohort <- cohort
    x[, c("cohort", setdiff(names(x), "cohort")), drop = FALSE]
  }))
  rownames(plan_output) <- NULL
  matching_audit <- do.call(rbind, lapply(task8_cohorts, function(cohort) {
    frame <- matching[[cohort]]
    plan <- plans[[cohort]]
    by_stage <- tapply(plan$n_allocate, plan$relaxation_stage, sum)
    data.frame(
      cohort = cohort, detectable_universe = nrow(frame), detectable_targets = sum(frame$is_target),
      non_target_candidates = sum(!frame$is_target), requested_iterations = task10_iterations,
      completed_iterations = nrow(null), every_set_size_valid = all_sizes_valid,
      any_within_set_duplicate = any_duplicate, any_target_contamination = any_target_contamination,
      exact_allocations_per_set = unname(by_stage["exact"] %||% 0L),
      probe_only_allocations_per_set = unname(by_stage["probe_only"] %||% 0L),
      adjacent_detection_allocations_per_set = unname(by_stage["adjacent_detection"] %||% 0L),
      expression_bin_relaxed = FALSE,
      expression_bins_present = paste(sort(unique(frame$expression_bin)), collapse = "|"),
      detection_bins_present = paste(sort(unique(frame$detection_bin)), collapse = "|"),
      probe_bins_present = paste(sort(unique(frame$probe_bin)), collapse = "|"),
      mean_expression_source = "Task5_feature_data_QC_retained_all_sample_mean_log2",
      detection_source = "Task5_feature_data_Detection_P_rate_QC_retained",
      probe_count_source = "Task5_mapping_audit_unique_current_probe_count_per_Entrez",
      forbidden_matching_fields_used = FALSE,
      stringsAsFactors = FALSE
    )
  }))

  output_paths <- c(
    calibration = file.path("results", "meta", "random_set_empirical_calibration.tsv"),
    null = file.path("results", "meta", "random_set_null_statistics.tsv"),
    matching_audit = file.path("results", "meta", "random_set_matching_audit.tsv"),
    relaxation_plan = file.path("results", "meta", "random_set_relaxation_plan.tsv"),
    camera_contract = file.path("results", "meta", "random_set_camera_equivalence.tsv"),
    diagnostics = file.path("results", "meta", "random_set_null_diagnostics.tsv"),
    extreme_check = file.path("results", "meta", "random_set_extreme_direct_camera_check.tsv"),
    figure = file.path("results", "figures", "random_set_null_distribution.pdf")
  )
  write_tsv(calibration, output_paths[["calibration"]])
  write_tsv(null, output_paths[["null"]])
  write_tsv(matching_audit, output_paths[["matching_audit"]])
  write_tsv(plan_output, output_paths[["relaxation_plan"]])
  write_tsv(camera_contract, output_paths[["camera_contract"]])
  write_tsv(formal_verification$diagnostics, output_paths[["diagnostics"]])
  write_tsv(formal_verification$extreme, output_paths[["extreme_check"]])
  write_null_figure(null, observed_z, output_paths[["figure"]])

  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  session_path <- file.path("logs", "session_info", "task10_session_info.txt")
  commands_path <- file.path("logs", "session_info", "task10_commands.log")
  anomaly_path <- file.path("logs", "session_info", "task10_anomaly_investigation.md")
  if (!file.exists(anomaly_path)) stop("Task 10 anomaly investigation record is absent", call. = FALSE)
  dir.create(dirname(session_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())), session_path, useBytes = TRUE)
  writeLines(c(
    "Task 10: 10,000 matched random gene-set empirical calibration",
    "Command: RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_CONFIG_NAMESPACES_CHECK=FALSE Rscript analysis/R/08_random_gene_sets.R",
    paste0("Seed: ", task10_seed, "; iterations completed: ", nrow(null)),
    "Scope: GSE105450 and GSE103842 confirmatory cohorts only; no sensitivity/exploratory cohort and no meta-analysis",
    "Sets: independently matched within each cohort universe, paired by iteration number",
    paste0("Set sizes: GSE105450=", task10_expected_detectable[["GSE105450"]], "; GSE103842=", task10_expected_detectable[["GSE103842"]]),
    "Matching: QC-retained mean expression quintile, Detection-P detectability quintile, unique-current probe-count quintile",
    "Relaxation: probe bin first, adjacent detection bin second, expression bin never relaxed",
    "Candidates: detectable non-target genes only; without replacement within every set",
    "camera: exact cached sufficient-statistic implementation validated against direct limma::camera at tolerance 1e-12; cameraPR not used",
    paste0("Observed Task 8 Stouffer Z reconstructed before draws: ", format(observed_z, digits = 16)),
    paste0("Elapsed seconds: ", format(elapsed, digits = 12)),
    "Primary decision remains H1_not_supported_camera_component_failed; random calibration is nonrescuing"
  ), commands_path, useBytes = TRUE)

  checksum_inputs <- c(
    file.path("analysis", "R", "08_random_gene_sets.R"),
    file.path("analysis", "tests", "test_random_gene_sets.R"),
    file.path("analysis", "tests", "test_random_gene_set_outputs.R"),
    task10_amendments,
    file.path("logs", "checksums", "task5_sha256.tsv"),
    file.path("logs", "checksums", "task7_sha256.tsv"),
    file.path("logs", "checksums", "task8_sha256.tsv"),
    file.path("data", "clean", "bfff_targets_primary.tsv"),
    file.path("data", "derived", paste0(task8_cohorts, "_expression.rds")),
    file.path("results", "tables", paste0(task8_cohorts, "_probe_entrez_mapping_audit.tsv")),
    file.path("results", "cohort", "target_set_camera.tsv"),
    file.path("results", "cohort", "target_set_primary_decision_partial.tsv")
  )
  checksum_outputs <- c(unname(output_paths), session_path, commands_path, anomaly_path)
  checksum_artifacts <- c(checksum_inputs, checksum_outputs)
  checksum <- data.frame(
    artifact = checksum_artifacts,
    role = c(rep("input", length(checksum_inputs)), rep("output", length(checksum_outputs))),
    bytes = as.numeric(file.info(checksum_artifacts)$size),
    sha256 = vapply(checksum_artifacts, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_tsv(checksum, file.path("logs", "checksums", "task10_sha256.tsv"))
  print(calibration, row.names = FALSE)
  print(matching_audit, row.names = FALSE)
  invisible(list(calibration = calibration, null = null, matching_audit = matching_audit))
}

`%||%` <- function(x, y) if (length(x) && !is.na(x)) x else y

if (!identical(Sys.getenv("TASK10_SKIP_MAIN"), "1")) main_task10()
