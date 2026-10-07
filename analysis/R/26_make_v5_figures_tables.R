#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
required <- c("ggplot2", "patchwork", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)

library(ggplot2)
library(patchwork)

asset_root <- file.path("manuscript", "v5_submission_assets")
figure_root <- file.path(asset_root, "figures")
table_root <- file.path(asset_root, "tables")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_root, recursive = TRUE, showWarnings = FALSE)

palette <- list(
  navy = "#17365D", blue = "#3C78A8", cyan = "#68A9C4", red = "#C9574E",
  orange = "#E59A52", green = "#5B8E71", purple = "#8064A2", grey = "#6B7280",
  light_blue = "#EAF2F8", light_red = "#FBE9E7", light_green = "#EAF4ED",
  light_purple = "#F0ECF7", light_grey = "#F3F4F6", dark = "#1F2937"
)

theme_v5 <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = palette$dark),
      plot.title = element_text(size = base_size + 3, face = "bold", margin = margin(b = 8)),
      plot.subtitle = element_text(size = base_size, colour = palette$grey, margin = margin(b = 8)),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, colour = palette$dark),
      strip.text = element_text(size = base_size + 1, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "#E5E7EB"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(12, 14, 12, 14)
    )
}

save_figure <- function(plot, stem, width, height) {
  png_path <- file.path(figure_root, paste0(stem, ".png"))
  ggsave(png_path, plot, width = width, height = height,
         units = "in", dpi = 360, bg = "white", limitsize = FALSE)
  ggsave(file.path(figure_root, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
  if (nzchar(Sys.which("sips"))) {
    system2("sips", c("-s", "dpiWidth", "360", "-s", "dpiHeight", "360", png_path),
            stdout = FALSE, stderr = FALSE)
  }
}

write_tsv <- function(x, name) {
  write.table(x, file.path(table_root, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

make_figure_1 <- function() {
  boxes <- data.frame(
    xmin = c(0.4, 3.1, 5.8, 8.8, 8.8, 12.2), xmax = c(2.6, 5.3, 8.2, 11.5, 11.5, 15.2),
    ymin = c(2.4, 2.4, 2.4, 3.8, 1.2, 2.4), ymax = c(4.0, 4.0, 4.0, 5.3, 2.7, 4.0),
    label = c(
      "Database space\n2,092 BFFF targets\n3,053 disease genes",
      "Candidate overlap\n734 genes\nsearch space",
      "Outcome-selected\ndiscovery set\n86 genes (60 up, 26 down)",
      "Blood replication\nGSE77087 + GSE103842\n77/86 two-cohort concordance",
      "Airway localization\nGSE97742 + GSE41374\nfixed 3-process framework",
      "Cross-context synthesis\n24 airway-supported\n16 blood-to-airway bridges"
    ),
    fill = c(palette$light_grey, palette$light_blue, palette$light_purple,
             palette$light_red, palette$light_green, "#FFF4E6")
  )
  arrows <- data.frame(
    x = c(2.6, 5.3, 8.2, 8.2, 11.5, 11.5), y = c(3.2, 3.2, 3.2, 3.2, 4.55, 1.95),
    xend = c(3.1, 5.8, 8.8, 8.8, 12.2, 12.2), yend = c(3.2, 3.2, 4.55, 1.95, 3.35, 3.05)
  )
  ggplot() +
    geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), linewidth = 0.8,
                 colour = palette$navy, arrow = grid::arrow(length = unit(0.14, "in"))) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              colour = palette$navy, linewidth = 0.8) +
    geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
              size = 3.35, lineheight = 1.1, colour = palette$dark) +
    annotate("rect", xmin = 5.8, xmax = 15.2, ymin = 0.25, ymax = 0.85,
             fill = "#FFF8E1", colour = "#B7791F", linewidth = 0.6) +
    annotate("text", x = 10.5, y = 0.55,
             label = "Evidence boundary: no cohort included BFFF exposure; disease association is not formula perturbation",
             size = 3.35, colour = "#7C4A03", fontface = "bold") +
    scale_fill_identity() + coord_cartesian(xlim = c(0.1, 15.4), ylim = c(0.05, 5.55), clip = "off") +
    labs(title = "V5 study design: from candidate space to cross-context knowledge",
         subtitle = "Blood cohorts establish reproducibility; airway cohorts localize the host-response program") +
    theme_void(base_family = "sans", base_size = 10) +
    theme(plot.title = element_text(size = 14, face = "bold", colour = palette$dark),
          plot.subtitle = element_text(size = 10, colour = palette$grey),
          plot.margin = margin(14, 16, 14, 16))
}

make_figure_2 <- function() {
  root <- file.path("results", "single_gene_upgrade_v1", "cohort")
  cohorts <- c("GSE38900", "GSE77087", "GSE103842")
  dat <- do.call(rbind, lapply(cohorts, function(cohort) {
    x <- read.delim(file.path(root, paste0(cohort, "_86_gene_effects.tsv")))
    x[, c("gene_symbol", "cohort", "log2FC", "fdr_bh")]
  }))
  gene_order <- read.delim("data/clean/bfff_86_gene_set_frozen.tsv")$gene_symbol
  blocks <- setNames(rep(c("Genes 1–43", "Genes 44–86"), each = 43), gene_order)
  dat$block <- unname(blocks[dat$gene_symbol])
  dat$gene_symbol <- factor(dat$gene_symbol, levels = rev(gene_order))
  dat$cohort <- factor(dat$cohort, levels = cohorts, labels = c("GSE38900\nDiscovery", "GSE77087\nReplication", "GSE103842\nReplication"))
  dat$mark <- ifelse(dat$fdr_bh < 0.05, "•", "")
  ggplot(dat, aes(cohort, gene_symbol, fill = pmax(pmin(log2FC, 2.5), -2.5))) +
    geom_tile(colour = "white", linewidth = 0.22) +
    geom_text(aes(label = mark), size = 2.4, colour = "black", na.rm = TRUE) +
    facet_wrap(~block, ncol = 2, scales = "free_y") +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
                         limits = c(-2.5, 2.5), oob = scales::squish,
                         name = "log2 fold change\n(clipped at ±2.5)") +
    labs(title = "Cross-cohort blood effects for the complete 86-gene candidate universe",
         subtitle = "Dots indicate cohort-specific transcriptome-wide BH FDR < 0.05; GSE38900 is outcome-selected discovery",
         x = NULL, y = NULL) +
    theme_v5(9.5) +
    theme(axis.text.y = element_text(size = 7.8), axis.text.x = element_text(size = 8.5),
          panel.grid = element_blank(), legend.position = "bottom",
          strip.background = element_rect(fill = palette$light_grey, colour = NA))
}

make_figure_3 <- function() {
  dat <- read.delim("results/barrier_inflammation_repair_v1/three_axis_camera_summary.tsv")
  comparison_labels <- c(
    GSE97742_RSVsi_acute_vs_discharge = "GSE97742 RSV\nAcute − discharge",
    GSE97742_hRV_acute_vs_discharge = "GSE97742 hRV\nAcute − discharge",
    GSE41374_RSV_vs_healthy = "GSE41374\nRSV − healthy"
  )
  axis_labels <- c(Barrier_cilia = "Barrier / cilia", Inflammation = "Inflammation",
                   Repair_turnover = "Repair / turnover")
  term_order <- unique(dat$term)
  dat$comparison_label <- factor(comparison_labels[dat$comparison], levels = comparison_labels)
  term_labels <- tools::toTitleCase(gsub("_", " ", dat$term))
  ordered_labels <- tools::toTitleCase(gsub("_", " ", term_order))
  dat$term_label <- factor(term_labels, levels = rev(ordered_labels))
  dat$axis_label <- factor(axis_labels[dat$axis], levels = axis_labels)
  direction_sign <- ifelse(dat$Direction == "Up", 1, -1)
  dat$signed_score <- direction_sign * pmin(-log10(dat$fdr_bh_fixed_family), 12)
  dat$mark <- ifelse(dat$fdr_bh_fixed_family < 0.05 & dat$direction_matches_expected, "*", "")
  ggplot(dat, aes(comparison_label, term_label, fill = signed_score)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = mark), size = 4.0, fontface = "bold") +
    facet_grid(axis_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
                         limits = c(-12, 12), oob = scales::squish,
                         name = "Signed −log10(FDR)\nBlue: down; red: up") +
    labs(title = "Airway transcriptomes reveal a stable ciliary–inflammatory core",
         subtitle = "Repair/turnover support is confined to paired acute−discharge comparisons; stars mark expected-direction FDR < 0.05",
         x = NULL, y = NULL) +
    theme_v5(10) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8.4),
          axis.text.x = element_text(size = 9), strip.placement = "outside",
          strip.background = element_rect(fill = palette$light_grey, colour = NA),
          legend.position = "bottom")
}

make_figure_4 <- function() {
  x <- read.delim("results/bfff_full_86_evidence_matrix_v1/Table_full_86_cross_context_evidence.tsv")
  sources <- list(
    "GSE38900\nBlood discovery" = list(effect = x$GSE38900_log2FC, fdr = x$GSE38900_fdr_bh,
                                        support = rep(TRUE, nrow(x))),
    "GSE77087\nBlood replication" = list(effect = x$GSE77087_log2FC, fdr = x$GSE77087_fdr_bh,
                                           support = x$blood_GSE77087_genome_fdr_support),
    "GSE103842\nBlood replication" = list(effect = x$GSE103842_log2FC, fdr = x$GSE103842_fdr_bh,
                                            support = x$blood_GSE103842_genome_fdr_support),
    "GSE97742\nAirway recovery interaction" = list(effect = x$GSE97742_RSV_minus_hRV_delta,
                                                     fdr = x$GSE97742_fdr_bh_86_family,
                                                     support = x$airway_GSE97742_fdr_support),
    "GSE41374\nAirway RSV−healthy" = list(effect = x$GSE41374_logFC_shifted_log2,
                                              fdr = x$GSE41374_fdr_bh_86_family_shifted_log2,
                                              support = x$airway_GSE41374_fdr_support)
  )
  dat <- do.call(rbind, lapply(names(sources), function(source) {
    z <- sources[[source]]
    data.frame(gene_symbol = x$gene_symbol, source = source, effect = z$effect, fdr = z$fdr,
               support = z$support, representative = x$is_representative)
  }))
  gene_order <- x$gene_symbol
  blocks <- setNames(rep(c("Genes 1–43", "Genes 44–86"), each = 43), gene_order)
  dat$block <- unname(blocks[dat$gene_symbol])
  dat$gene_symbol <- factor(dat$gene_symbol, levels = rev(gene_order))
  dat$source <- factor(dat$source, levels = names(sources))
  dat$signed_score <- sign(dat$effect) * pmin(-log10(dat$fdr), 6)
  dat$mark <- ifelse(dat$support %in% TRUE, "•", "")
  ggplot(dat, aes(source, gene_symbol, fill = signed_score)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = mark), size = 2.2, colour = "black", na.rm = TRUE) +
    facet_wrap(~block, ncol = 2, scales = "free_y") +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
                         limits = c(-6, 6), oob = scales::squish,
                         na.value = "#E5E7EB", name = "Signed −log10(FDR)\n(clipped at ±6)") +
    labs(title = "Complete 86-gene blood-to-airway evidence matrix",
         subtitle = "Dots denote prespecified support in each context; effect scales are context-specific and are not pooled",
         x = NULL, y = NULL) +
    theme_v5(9.5) +
    theme(axis.text.y = element_text(size = 7.8), axis.text.x = element_text(size = 7.9, angle = 25, hjust = 1),
          panel.grid = element_blank(), legend.position = "bottom",
          strip.background = element_rect(fill = palette$light_grey, colour = NA))
}

make_figure_5 <- function() {
  boxes <- data.frame(
    xmin = c(0.4, 3.0, 3.0, 5.8, 5.8, 8.7), xmax = c(2.2, 5.1, 5.1, 8.0, 8.0, 11.1),
    ymin = c(2.4, 3.6, 1.2, 3.6, 1.2, 2.4), ymax = c(3.8, 5.0, 2.6, 5.0, 2.6, 3.8),
    label = c(
      "Acute pediatric\nviral airway state",
      "Ciliary / barrier impairment\nCilium organization: decreased\nCilium movement: decreased",
      "Inflammatory activation\nInterferon, myeloid and\nneutrophil programs: increased",
      "Phase-dependent repair / turnover\nWound healing and apoptosis regulation\nEpithelial proliferation: increased",
      "Systemic blood counterpart\nInterferon–myeloid: increased\nLymphocyte modules: decreased",
      "16 blood-to-airway\nbridge candidates\nCross-context priorities"
    ),
    fill = c("#FFF4E6", palette$light_blue, palette$light_red,
             palette$light_green, palette$light_purple, "#E8F1F8")
  )
  solid <- data.frame(
    x = c(2.2, 2.2, 5.1, 5.1, 8.0, 8.0), y = c(3.1, 3.1, 4.3, 1.9, 4.3, 1.9),
    xend = c(3.0, 3.0, 5.8, 5.8, 8.7, 8.7), yend = c(4.3, 1.9, 4.3, 1.9, 3.3, 2.9)
  )
  ggplot() +
    geom_segment(data = solid, aes(x, y, xend = xend, yend = yend), colour = palette$navy,
                 linewidth = 0.85, arrow = grid::arrow(length = unit(0.13, "in"))) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              colour = palette$navy, linewidth = 0.8) +
    geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
              size = 2.9, lineheight = 1.08) +
    annotate("rect", xmin = 8.2, xmax = 11.5, ymin = 0.35, ymax = 1.35,
             fill = "#FFF8E1", colour = "#B7791F", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = 9.9, y = 0.85,
             label = "BFFF perturbation hypothesis\nUntested: exposure, direction, target engagement",
             size = 2.7, colour = "#7C4A03") +
    annotate("segment", x = 9.9, y = 1.35, xend = 9.9, yend = 2.4,
             colour = "#B7791F", linetype = "dashed", linewidth = 0.8,
             arrow = grid::arrow(length = unit(0.12, "in"))) +
    annotate("text", x = 6.9, y = 5.4, label = "Data-derived host-response model",
             size = 3.6, fontface = "bold", colour = palette$navy) +
    annotate("text", x = 6.9, y = 0.25,
             label = "Solid arrows: observed cross-context organization   |   Dashed arrow: future BFFF intervention hypothesis",
             size = 3.0, colour = palette$grey) +
    scale_fill_identity() + coord_cartesian(xlim = c(0.15, 11.35), ylim = c(0.05, 5.7), clip = "off") +
    labs(title = "A ciliary–inflammatory core with phase-dependent compensatory repair",
         subtitle = "The model organizes disease-associated evidence without implying a causal sequence or formula action") +
    theme_void(base_family = "sans", base_size = 10) +
    theme(plot.title = element_text(size = 14, face = "bold", colour = palette$dark),
          plot.subtitle = element_text(size = 10, colour = palette$grey),
          plot.margin = margin(14, 16, 14, 16))
}

make_tables <- function() {
  table1 <- data.frame(
    GEO_accession = c("GSE38900", "GSE77087", "GSE103842", "GSE97742", "GSE41374"),
    Platform = c("GPL6884 (predominant)", "GPL10558", "GPL10558", "GPL10558", "GPL10558"),
    Compartment = c("Whole blood", "Whole blood", "Whole blood", "Nasopharyngeal swab", "Nasal wash"),
    Design = c("Cross-sectional", "Cross-sectional", "Cross-sectional", "Paired longitudinal", "Cross-sectional"),
    Main_comparison = c("RSV vs healthy", "RSV vs healthy", "RSV vs healthy",
                        "Acute vs discharge; RSV vs hRV change", "RSV vs healthy"),
    Analysis_n = c(138, 104, 73, 136, 86),
    Group_detail = c("107 RSV; 31 healthy", "81 RSV; 23 healthy", "61 RSV; 12 healthy",
                     "38 RSV pairs; 30 hRV pairs", "76 RSV; 10 healthy"),
    Evidence_role = c("Outcome-selected discovery", "Primary exploratory blood replication",
                      "Second blood replication", "Longitudinal airway localization",
                      "Independent cross-sectional airway localization"),
    Model_or_adjustment = c("Existing discovery workflow", "Age + sex + technical batch",
                            "Age + sex + technical batch", "Subject blocking; interaction adjusted for age + sex",
                            "Age + sex; two-transformation robustness"),
    Interpretation_boundary = c(
      "Not independent validation", "No BFFF exposure", "No BFFF exposure",
      "Discharge is not a healthy baseline", "Small healthy group; processed-matrix transformation sensitivity"
    )
  )
  write_tsv(table1, "Table_1_cohort_characteristics.tsv")

  evidence <- read.delim("results/bfff_full_86_evidence_matrix_v1/evidence_summary.tsv")
  axis <- read.delim("results/bfff_full_86_evidence_matrix_v1/axis_evidence_summary.tsv")
  evidence_labels <- c(
    target_universe = "Frozen candidate universe",
    representative_nodes = "Representative annotated nodes",
    GSE77087_concordant_genome_fdr_support = "GSE77087 concordant genome-wide FDR support",
    GSE103842_concordant_genome_fdr_support = "GSE103842 concordant genome-wide FDR support",
    GSE97742_airway_fdr_support = "GSE97742 recovery-interaction support",
    GSE41374_airway_robust_fdr_support = "GSE41374 two-transformation robust support",
    any_airway_fdr_support = "Support in at least one airway analysis",
    blood_GSE77087_and_any_airway_support = "Blood-to-airway bridge candidates"
  )
  axis_labels <- c(
    Other_unclassified = "Other/unclassified: total / blood / airway / bridge",
    Inflammation = "Inflammation: total / blood / airway / bridge",
    Repair_stress = "Repair/stress: total / blood / airway / bridge",
    Barrier_remodeling = "Barrier/remodeling: total / blood / airway / bridge"
  )
  table2 <- rbind(
    data.frame(Section = "Overall evidence", Metric = unname(evidence_labels[evidence$metric]), Value = evidence$value,
               Interpretation = c(
                 "Formal frozen candidate universe", "Representative annotations only",
                 "Direction-concordant and transcriptome-wide FDR < 0.05",
                 "Direction-concordant and transcriptome-wide FDR < 0.05",
                 "RSV-vs-hRV recovery interaction; 86-gene-family FDR < 0.05",
                 "Concordant and FDR < 0.05 under both transformations",
                 "Supported in at least one airway analysis",
                 "GSE77087 blood support plus at least one airway analysis"
               )),
    data.frame(Section = "Representative-axis distribution",
               Metric = unname(axis_labels[axis$axis]),
               Value = paste(axis$n_genes, axis$n_GSE77087_concordant_genome_fdr,
                             axis$n_any_airway_fdr, axis$n_blood_and_airway, sep = "/"),
               Interpretation = "Counts: total genes / GSE77087 blood support / any airway support / blood-to-airway support")
  )
  write_tsv(table2, "Table_2_86_gene_evidence_summary.tsv")

  x <- read.delim("results/bfff_full_86_evidence_matrix_v1/Table_full_86_cross_context_evidence.tsv")
  bridge <- x$blood_GSE77087_genome_fdr_support & x$airway_any_fdr_support
  table3 <- x[bridge, c(
    "gene_symbol", "is_representative", "axis",
    "GSE38900_log2FC", "GSE38900_fdr_bh",
    "GSE77087_log2FC", "GSE77087_fdr_bh",
    "GSE103842_log2FC", "GSE103842_fdr_bh",
    "GSE97742_RSV_minus_hRV_delta", "GSE97742_fdr_bh_86_family",
    "GSE41374_logFC_shifted_log2", "GSE41374_fdr_bh_86_family_shifted_log2",
    "airway_GSE97742_fdr_support", "airway_GSE41374_fdr_support", "n_external_contexts_fdr"
  )]
  names(table3) <- c(
    "Gene", "Representative_node", "Representative_axis",
    "GSE38900_log2FC", "GSE38900_genome_FDR",
    "GSE77087_log2FC", "GSE77087_genome_FDR",
    "GSE103842_log2FC", "GSE103842_genome_FDR",
    "GSE97742_RSV_minus_hRV_recovery_delta", "GSE97742_86_gene_FDR",
    "GSE41374_logFC_shifted_log2", "GSE41374_86_gene_FDR_shifted_log2",
    "GSE97742_support", "GSE41374_robust_support", "External_context_support_count"
  )
  write_tsv(table3, "Table_3_16_bridge_candidates.tsv")
}

main <- function() {
  save_figure(make_figure_1(), "Figure_1", 15.5, 6.6)
  save_figure(make_figure_2(), "Figure_2", 11.8, 12.4)
  save_figure(make_figure_3(), "Figure_3", 10.5, 9.2)
  save_figure(make_figure_4(), "Figure_4", 13.2, 12.4)
  save_figure(make_figure_5(), "Figure_5", 11.8, 6.8)
  make_tables()
  capture.output(sessionInfo(), file = file.path(asset_root, "figure_table_session_info.txt"))
  message("V5 figures and TSV tables created under ", asset_root)
}

if (sys.nframe() == 0L) main()
