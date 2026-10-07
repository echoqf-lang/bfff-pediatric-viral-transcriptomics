#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
source(file.path("analysis", "R", "single_gene_upgrade_pipeline.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package is required", call. = FALSE)

root <- file.path("results", "single_gene_upgrade_v1")
palette_cb <- c(blue = "#0072B2", orange = "#D55E00", green = "#009E73", purple = "#CC79A7", grey = "#7A7A7A")

theme_publication <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "#222222", size = 10),
      axis.text = ggplot2::element_text(colour = "#222222", size = 8),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(colour = "#444444", size = 9),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

save_publication_plot <- function(plot, stem, width, height) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(paste0(stem, ".pdf"), width = width, height = height, useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()
  grDevices::svg(paste0(stem, ".svg"), width = width, height = height, pointsize = 10)
  print(plot)
  grDevices::dev.off()
  grDevices::png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 600, type = "cairo")
  print(plot)
  grDevices::dev.off()
}

split_ordered_rows <- function(x, n_panels = 2L) {
  n <- nrow(x)
  x$display_panel <- factor(
    pmin(n_panels, ceiling(seq_len(n) / ceiling(n / n_panels))),
    levels = seq_len(n_panels), labels = letters[seq_len(n_panels)]
  )
  x
}

prefix_columns <- function(x, prefix, keys = c("gene_id", "gene_symbol")) {
  names(x)[!names(x) %in% keys] <- paste0(prefix, names(x)[!names(x) %in% keys])
  x
}

make_tables <- function() {
  meta <- read_upgrade_tsv(file.path(root, "meta", "three_cohort_86_gene_reml_hksj.tsv"))
  meta <- meta[order(meta$adhoc_p_value, meta$gene_symbol), ]
  write_upgrade_tsv(meta, file.path(root, "tables", "Table_main_three_cohort_meta.tsv"))

  frozen <- read_upgrade_tsv(file.path("data", "clean", "bfff_86_gene_set_frozen.tsv"), col_classes = "character")
  base <- frozen[c("entrez_id", "gene_symbol")]
  names(base)[1] <- "gene_id"
  cohorts <- lapply(c("GSE38900", "GSE77087", "GSE103842"), function(cohort) {
    x <- read_upgrade_tsv(file.path(root, "cohort", paste0(cohort, "_86_gene_effects.tsv")))
    x <- x[c("gene_id", "log2FC", "SE", "p_value", "fdr_bh")]
    prefix_columns(x, paste0(cohort, "_"), keys = "gene_id")
  })
  supp <- Reduce(function(x, y) merge(x, y, by = "gene_id", all.x = TRUE, sort = FALSE), c(list(base), cohorts))
  meta_supp <- meta[c(
    "gene_id", "estimate",
    "se", "ci_low", "ci_high", "p_value", "fdr_bh",
    "adhoc_se", "adhoc_ci_low", "adhoc_ci_high", "adhoc_p_value", "adhoc_fdr_bh",
    "tau2", "i2", "direction_pattern"
  )]
  names(meta_supp) <- c(
    "gene_id", "meta_estimate",
    "meta_conventional_hksj_se", "meta_conventional_hksj_ci_low",
    "meta_conventional_hksj_ci_high", "meta_conventional_hksj_p_value",
    "meta_conventional_hksj_fdr_bh",
    "meta_modified_hksj_se", "meta_modified_hksj_ci_low",
    "meta_modified_hksj_ci_high", "meta_modified_hksj_p_value",
    "meta_modified_hksj_fdr_bh",
    "meta_tau2", "meta_i2", "meta_direction_pattern"
  )
  supp <- merge(supp, meta_supp, by = "gene_id", all.x = TRUE, sort = FALSE)
  severity <- read_upgrade_tsv(file.path(root, "severity", "GSE77087_86_gene_ordinal_trend.tsv"))
  severity <- prefix_columns(severity[c("gene_id", "log2FC", "SE", "p_value", "fdr_bh")], "severity_trend_", keys = "gene_id")
  supp <- merge(supp, severity, by = "gene_id", all.x = TRUE, sort = FALSE)
  specificity <- read_upgrade_tsv(file.path(root, "specificity", "GSE155925_86_gene_rsv_vs_other_virus.tsv"))
  specificity <- prefix_columns(specificity[c("gene_id", "log2FC", "SE", "p_value", "candidate_fdr", "specificity_class")], "specificity_", keys = "gene_id")
  supp <- merge(supp, specificity, by = "gene_id", all.x = TRUE, sort = FALSE)
  supp <- supp[match(base$gene_id, supp$gene_id), ]
  write_upgrade_tsv(supp, file.path(root, "tables", "Table_supplement_all_86_genes.tsv"))
  invisible(list(main = meta, supplement = supp))
}

plot_multicohort_effects <- function() {
  meta <- read_upgrade_tsv(file.path(root, "meta", "three_cohort_86_gene_reml_hksj.tsv"))
  meta <- meta[order(meta$estimate), ]
  meta <- split_ordered_rows(meta)
  meta$gene_symbol <- factor(meta$gene_symbol, levels = meta$gene_symbol)
  meta$direction <- ifelse(meta$estimate >= 0, "Positive", "Negative")
  ggplot2::ggplot(meta, ggplot2::aes(x = estimate, y = gene_symbol, colour = direction)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#555555", linewidth = 0.35, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = adhoc_ci_low, xmax = adhoc_ci_high), width = 0, linewidth = 0.45) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(~display_panel, nrow = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c(Positive = palette_cb[["orange"]], Negative = palette_cb[["blue"]])) +
    ggplot2::labs(title = "Three-cohort random-effects estimates for the frozen candidate genes",
                  subtitle = "All 81 genes with complete effects; conservative modified HKSJ 95% CIs",
                  x = "Pooled log2 fold change (RSV - healthy)", y = NULL) +
    theme_publication()
}

plot_severity_gradient <- function() {
  x <- read_upgrade_tsv(file.path(root, "severity", "GSE77087_86_gene_ordinal_trend.tsv"))
  x$ci_low <- x$log2FC - 1.96 * x$SE
  x$ci_high <- x$log2FC + 1.96 * x$SE
  x <- x[order(x$log2FC), ]
  x <- split_ordered_rows(x)
  x$gene_symbol <- factor(x$gene_symbol, levels = x$gene_symbol)
  x$trend_fdr <- ifelse(x$fdr_bh < 0.05, "BH FDR < 0.05", "BH FDR >= 0.05")
  ggplot2::ggplot(x, ggplot2::aes(x = log2FC, y = gene_symbol, colour = trend_fdr)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#555555", linewidth = 0.35, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.45) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(~display_panel, nrow = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c("BH FDR < 0.05" = palette_cb[["orange"]], "BH FDR >= 0.05" = palette_cb[["grey"]])) +
    ggplot2::labs(title = "Ordered severity gradient in GSE77087",
                  subtitle = "All 86 genes; one-level increase from healthy to outpatient to hospitalized",
                  x = "Adjusted log2 fold change per severity level", y = NULL) +
    theme_publication()
}

plot_virus_specificity <- function() {
  x <- read_upgrade_tsv(file.path(root, "specificity", "GSE155925_86_gene_rsv_vs_other_virus.tsv"))
  x <- x[is.finite(x$log2FC) & is.finite(x$SE), ]
  x <- x[order(x$log2FC), ]
  x <- split_ordered_rows(x)
  x$gene_symbol <- factor(x$gene_symbol, levels = x$gene_symbol)
  x$classification <- ifelse(x$candidate_fdr < 0.05, "Candidate-family BH FDR < 0.05", "Not distinguished at current precision")
  ggplot2::ggplot(x, ggplot2::aes(x = log2FC, y = gene_symbol, colour = classification)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#555555", linewidth = 0.35, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.45) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(~display_panel, nrow = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c("Candidate-family BH FDR < 0.05" = palette_cb[["purple"]],
                                            "Not distinguished at current precision" = palette_cb[["grey"]])) +
    ggplot2::labs(title = "RSV versus other single-virus infections in GSE155925",
                  subtitle = "All 79 estimable frozen genes; positive values indicate higher expression in RSV",
                  x = "Adjusted log2 fold change (RSV - other virus)", y = NULL) +
    theme_publication()
}

plot_blood_modules <- function() {
  effects <- read_upgrade_tsv(file.path(root, "modules", "cohort_module_effects.tsv"))
  replication <- read_upgrade_tsv(file.path(root, "modules", "module_replication_summary.tsv"))
  ids <- replication$module_id[replication$evaluable_all_three & replication$significant_cohort_n == 3L &
                                 replication$direction_pattern %in% c("all_up", "all_down")]
  x <- effects[effects$module_id %in% ids, ]
  x$label <- paste0(x$module_id, " | ", x[["function"]])
  ordering <- aggregate(median_log2FC ~ label, x, mean)
  ordering <- ordering$label[order(ordering$median_log2FC)]
  x$label <- factor(x$label, levels = ordering)
  panel_lookup <- split_ordered_rows(data.frame(label = ordering))
  x$display_panel <- panel_lookup$display_panel[match(as.character(x$label), panel_lookup$label)]
  x$cohort <- factor(x$cohort, levels = c("GSE38900", "GSE77087", "GSE103842"))
  limit <- max(abs(x$median_log2FC), na.rm = TRUE)
  ggplot2::ggplot(x, ggplot2::aes(x = cohort, y = label, fill = median_log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::facet_wrap(~display_panel, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_gradient2(low = palette_cb[["blue"]], mid = "white", high = palette_cb[["orange"]],
                                  midpoint = 0, limits = c(-limit, limit), name = "Median\nlog2FC") +
    ggplot2::labs(title = "Replicated BloodGen3 transcriptional modules",
                  subtitle = "68 modules with BH FDR < 0.05 and concordant median effects in all three cohorts",
                  x = NULL, y = NULL) +
    theme_publication() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      axis.text.x = ggplot2::element_text(face = "bold", size = 8, angle = 45, hjust = 1, vjust = 1),
      panel.spacing.x = grid::unit(8, "mm"),
      legend.position = "right"
    )
}

main_make_upgrade_outputs <- function() {
  tables <- make_tables()
  figure_root <- file.path(root, "figures")
  save_publication_plot(plot_multicohort_effects(), file.path(figure_root, "Fig_multicohort_effects"), 6.69, 8.25)
  save_publication_plot(plot_severity_gradient(), file.path(figure_root, "Fig_severity_gradient"), 6.69, 8.25)
  save_publication_plot(plot_virus_specificity(), file.path(figure_root, "Fig_virus_specificity"), 6.69, 8.25)
  save_publication_plot(plot_blood_modules(), file.path(figure_root, "Fig_blood_modules"), 6.69, 8.25)
  message("OUTPUTS main_table=", nrow(tables$main), " supplement=", nrow(tables$supplement), " figures=4x3_formats")
  invisible(tables)
}

if (sys.nframe() == 0L) main_make_upgrade_outputs()
