#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

# Reuse the validated V5 data assembly and plotting helpers, while keeping all
# V6 outputs isolated and redefining only the journal-facing composition.
source(file.path("analysis", "R", "26_make_v5_figures_tables.R"), local = FALSE)

asset_root <- file.path("manuscript", "v6_bmc_submission_assets")
figure_root <- file.path(asset_root, "figures")
table_root <- file.path(asset_root, "tables")
supplement_root <- file.path(asset_root, "supplementary_tables")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_root, recursive = TRUE, showWarnings = FALSE)
dir.create(supplement_root, recursive = TRUE, showWarnings = FALSE)

theme_v5 <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = palette$dark),
      plot.title = element_text(size = base_size + 3.5, face = "bold", margin = margin(b = 8)),
      plot.subtitle = element_text(size = base_size, colour = palette$grey, margin = margin(b = 9)),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 0.5, colour = palette$dark),
      strip.text = element_text(size = base_size + 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "#E5E7EB"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      plot.margin = margin(14, 16, 14, 16)
    )
}

make_figure_1 <- function() {
  boxes <- data.frame(
    xmin = c(0.4, 3.2, 6.3, 6.3, 9.5), xmax = c(2.7, 5.8, 8.9, 8.9, 12.6),
    ymin = c(2.2, 2.2, 3.55, 0.85, 2.2), ymax = c(3.9, 3.9, 5.15, 2.45, 3.9),
    label = c(
      "Frozen candidate universe\n86 outcome-selected genes\n(60 up, 26 down)",
      "Blood reproducibility\nGSE38900 + GSE77087 + GSE103842\ncomplete cohort-specific effects",
      "Airway functional localization\nGSE97742 + GSE41374\nfixed 16-set framework",
      "Viral-context analyses\nseverity, BloodGen3 and GSE155925\nnegative results retained",
      "Cross-context synthesis\n24 airway-supported genes\n16 blood-to-airway candidates"
    ),
    fill = c(palette$light_purple, palette$light_red, palette$light_green,
             palette$light_blue, "#FFF4E6")
  )
  arrows <- data.frame(
    x = c(2.7, 5.8, 5.8, 8.9, 8.9), y = c(3.05, 3.05, 3.05, 4.35, 1.65),
    xend = c(3.2, 6.3, 6.3, 9.5, 9.5), yend = c(3.05, 4.35, 1.65, 3.35, 2.75)
  )
  ggplot() +
    geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), linewidth = 0.9,
                 colour = palette$navy, arrow = grid::arrow(length = unit(0.14, "in"))) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              colour = palette$navy, linewidth = 0.9) +
    geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
              size = 3.35, lineheight = 1.1, colour = palette$dark) +
    annotate("rect", xmin = 0.4, xmax = 5.8, ymin = 0.10, ymax = 0.78,
             fill = palette$light_grey, colour = palette$grey, linewidth = 0.55) +
    annotate("text", x = 3.1, y = 0.44,
             label = "Candidate provenance: 2,092 BFFF-related targets intersected with\n3,053 disease genes to define a 734-gene search space",
             size = 2.85, lineheight = 1.05, colour = palette$grey) +
    annotate("rect", xmin = 6.3, xmax = 12.6, ymin = 0.10, ymax = 0.78,
             fill = "#FFF8E1", colour = "#B7791F", linewidth = 0.6) +
    annotate("text", x = 9.45, y = 0.44,
             label = "Evidence boundary: no cohort included BFFF exposure;\nintervention effects remain untested",
             size = 2.9, lineheight = 1.05, colour = "#7C4A03", fontface = "bold") +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.15, 12.85), ylim = c(0.02, 5.45), clip = "off") +
    labs(title = "Multi-cohort transcriptomic integration across blood and airway",
         subtitle = "Candidate provenance is separated from disease-response evidence and the untested intervention layer") +
    theme_void(base_family = "sans", base_size = 11) +
    theme(plot.title = element_text(size = 15, face = "bold", colour = palette$dark),
          plot.subtitle = element_text(size = 11, colour = palette$grey),
          plot.margin = margin(16, 18, 16, 18))
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
  dat$cohort <- factor(dat$cohort, levels = cohorts,
                       labels = c("GSE38900\nDiscovery", "GSE77087\nReplication", "GSE103842\nReplication"))
  dat$mark <- ifelse(dat$fdr_bh < 0.05, "•", "")
  heat <- ggplot(dat, aes(cohort, gene_symbol, fill = pmax(pmin(log2FC, 2.5), -2.5))) +
    geom_tile(colour = "white", linewidth = 0.24) +
    geom_text(aes(label = mark), size = 2.55, colour = "black", na.rm = TRUE) +
    facet_wrap(~block, ncol = 2, scales = "free_y") +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
                         limits = c(-2.5, 2.5), oob = scales::squish,
                         name = "log2 fold change\n(clipped at ±2.5)") +
    labs(title = "A  |  Complete cohort-specific effects",
         subtitle = "Dots: cohort-specific transcriptome-wide BH FDR < 0.05",
         x = NULL, y = NULL) +
    theme_v5(10.5) +
    theme(axis.text.y = element_text(size = 8.6), axis.text.x = element_text(size = 9.5),
          panel.grid = element_blank(), legend.position = "bottom",
          strip.background = element_rect(fill = palette$light_grey, colour = NA))

  summary <- read.delim("results/single_gene_upgrade_v1/meta/three_cohort_meta_summary.tsv")
  metrics <- data.frame(
    label = c("GSE38900–GSE77087\ndirection concordance",
              "All three cohorts\nsame direction",
              "Conventional HKSJ\nBH FDR < 0.05",
              "Modified HKSJ\nBH FDR < 0.05"),
    value = c(summary$gse38900_gse77087_direction_concordant_n,
              summary$all_three_same_direction_n,
              summary$meta_bh_lt_0_05_n,
              summary$meta_adhoc_bh_lt_0_05_n),
    denominator = c(86, summary$meta_complete_n, summary$meta_complete_n, summary$meta_complete_n),
    order = 4:1
  )
  metrics$label <- factor(metrics$label, levels = metrics$label[order(metrics$order)])
  metrics$ratio <- metrics$value / metrics$denominator
  metrics$value_label <- paste0(metrics$value, "/", metrics$denominator)
  robust <- ggplot(metrics, aes(ratio, label)) +
    geom_col(width = 0.58, fill = c(palette$grey, palette$orange, palette$blue, palette$navy)) +
    geom_text(aes(label = value_label), hjust = -0.12, size = 4.0, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1.13), breaks = c(0, 0.5, 1), labels = scales::percent) +
    labs(title = "B  |  Reproducibility and small-k robustness",
         subtitle = "Conservative variance handling removes nominal meta-analytic discoveries",
         x = "Proportion of evaluable genes", y = NULL) +
    theme_v5(11) +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 10),
          legend.position = "none")

  heat / robust + patchwork::plot_layout(heights = c(3.2, 1)) +
    patchwork::plot_annotation(
      title = "Blood effects across the frozen 86-gene universe",
      subtitle = "Direction reproducibility is broad, whereas meta-analytic significance is sensitive to conservative small-k inference",
      theme = theme(plot.title = element_text(size = 15, face = "bold", colour = palette$dark),
                    plot.subtitle = element_text(size = 11, colour = palette$grey),
                    plot.margin = margin(10, 12, 4, 12)))
}

make_figure_5 <- function() {
  boxes <- data.frame(
    xmin = c(0.5, 3.2, 3.2, 6.1, 6.1, 9.1), xmax = c(2.4, 5.4, 5.4, 8.3, 8.3, 11.6),
    ymin = c(2.35, 3.65, 1.05, 3.65, 1.05, 2.35), ymax = c(3.95, 5.25, 2.65, 5.25, 2.65, 3.95),
    label = c(
      "Acute pediatric\nviral airway state",
      "Stable airway core\nCiliary programs decreased",
      "Stable airway core\nInflammatory programs increased",
      "Phase-dependent response\nRepair / turnover increased\nduring acute-to-discharge change",
      "Systemic counterpart\nInterferon–myeloid increased\nLymphocyte modules decreased",
      "Cross-context priorities\n16 blood-to-airway\ncandidate genes"
    ),
    fill = c("#FFF4E6", palette$light_blue, palette$light_red,
             palette$light_green, palette$light_purple, "#E8F1F8")
  )
  links <- data.frame(
    x = c(2.4, 2.4, 5.4, 5.4, 8.3, 8.3), y = c(3.15, 3.15, 4.45, 1.85, 4.45, 1.85),
    xend = c(3.2, 3.2, 6.1, 6.1, 9.1, 9.1), yend = c(4.45, 1.85, 4.45, 1.85, 3.45, 2.85)
  )
  ggplot() +
    geom_segment(data = links, aes(x, y, xend = xend, yend = yend),
                 colour = palette$navy, linewidth = 0.9) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              colour = palette$navy, linewidth = 0.9) +
    geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
              size = 3.25, lineheight = 1.1) +
    annotate("rect", xmin = 8.55, xmax = 12.0, ymin = 0.03, ymax = 0.92,
             fill = "#FFF8E1", colour = "#B7791F", linetype = "dashed", linewidth = 0.85) +
    annotate("text", x = 10.275, y = 0.475,
             label = "BFFF intervention layer (untested)\nExposure, direction and target engagement require experiments",
             size = 2.55, lineheight = 1.05, colour = "#7C4A03") +
    annotate("segment", x = 10.25, y = 0.92, xend = 10.25, yend = 2.35,
             colour = "#B7791F", linetype = "dashed", linewidth = 0.85) +
    annotate("text", x = 6.3, y = 5.65, label = "Observed transcriptomic organization",
             size = 4.0, fontface = "bold", colour = palette$navy) +
    annotate("text", x = 4.45, y = 0.38,
             label = "Solid connectors organize co-occurring evidence;\nno causal or temporal sequence is implied",
             size = 2.7, lineheight = 1.05, colour = palette$grey) +
    scale_fill_identity() + coord_cartesian(xlim = c(0.2, 12.1), ylim = c(0, 5.95), clip = "off") +
    labs(title = "A ciliary–inflammatory core with phase-dependent repair",
         subtitle = "The disease-response model is data-derived; formula action remains a separate hypothesis") +
    theme_void(base_family = "sans", base_size = 11) +
    theme(plot.title = element_text(size = 15, face = "bold", colour = palette$dark),
          plot.subtitle = element_text(size = 11, colour = palette$grey),
          plot.margin = margin(16, 18, 16, 18))
}

make_tables_v6 <- function() {
  table1 <- data.frame(
    GEO_accession = c("GSE38900", "GSE77087", "GSE103842", "GSE97742", "GSE41374", "GSE155925"),
    Platform = c("GPL6884 (predominant)", "GPL10558", "GPL10558", "GPL10558", "GPL10558", "GPL16791"),
    Compartment = c("Whole blood", "Whole blood", "Whole blood", "Nasopharyngeal swab", "Nasal wash", "Whole blood"),
    Design = c("Cross-sectional", "Cross-sectional", "Cross-sectional", "Paired longitudinal", "Cross-sectional", "Cross-sectional"),
    Main_comparison = c("RSV vs healthy", "RSV vs healthy", "RSV vs healthy",
                        "Acute vs discharge; RSV vs hRV change", "RSV vs healthy",
                        "Single RSV vs other single-virus infections"),
    Analysis_n = c(138, 104, 73, 136, 86, 48),
    Group_detail = c("107 RSV; 31 healthy", "81 RSV; 23 healthy", "61 RSV; 12 healthy",
                     "38 RSV pairs; 30 hRV pairs", "76 RSV; 10 healthy", "31 RSV; 17 other-virus"),
    Evidence_role = c("Outcome-selected discovery", "Primary exploratory blood replication",
                      "Second blood replication", "Longitudinal airway localization",
                      "Independent cross-sectional airway localization", "Exploratory viral-context comparison"),
    Model_or_adjustment = c("Existing discovery workflow", "Age + sex + technical batch",
                            "Age + sex + technical batch", "Subject blocking; interaction adjusted for age + sex",
                            "Age + sex; two-transformation robustness",
                            "Age + sex + hospital batch + enrollment batch"),
    Interpretation_boundary = c(
      "Not independent validation", "No BFFF exposure", "No BFFF exposure",
      "Discharge is not a healthy baseline", "Small healthy group; transformation sensitivity",
      "No healthy controls; not an independent confirmation"
    )
  )
  write_tsv(table1, "Table_1_cohort_characteristics.tsv")

  evidence <- read.delim("results/bfff_full_86_evidence_matrix_v1/evidence_summary.tsv")
  axis <- read.delim("results/bfff_full_86_evidence_matrix_v1/axis_evidence_summary.tsv")
  summary <- read.delim("results/single_gene_upgrade_v1/meta/three_cohort_meta_summary.tsv")
  table2 <- data.frame(
    Evidence_layer = c("Frozen candidate universe", "Blood direction reproducibility",
                       "Three-cohort direction reproducibility", "Conventional HKSJ meta-analysis",
                       "Modified HKSJ sensitivity", "GSE77087 genome-wide blood support",
                       "GSE103842 genome-wide blood support", "Any airway support", "Blood-to-airway candidates"),
    Numerator = c(86, summary$gse38900_gse77087_direction_concordant_n,
                  summary$all_three_same_direction_n, summary$meta_bh_lt_0_05_n,
                  summary$meta_adhoc_bh_lt_0_05_n,
                  evidence$value[evidence$metric == "GSE77087_concordant_genome_fdr_support"],
                  evidence$value[evidence$metric == "GSE103842_concordant_genome_fdr_support"],
                  evidence$value[evidence$metric == "any_airway_fdr_support"],
                  evidence$value[evidence$metric == "blood_GSE77087_and_any_airway_support"]),
    Denominator = c(86, 86, summary$meta_complete_n, summary$meta_complete_n,
                    summary$meta_complete_n, 86, 81, 86, 86),
    Interpretation = c(
      "Formal analysis universe", "GSE38900 and GSE77087 same effect direction",
      "Same effect direction across all estimable blood cohorts",
      "REML-HKSJ with conventional variance estimate; exploratory",
      "Ad hoc modified variance estimate; conservative small-k diagnostic",
      "Direction-concordant and transcriptome-wide FDR < 0.05",
      "Direction-concordant and transcriptome-wide FDR < 0.05",
      "Support in at least one prespecified airway analysis",
      "GSE77087 blood support plus at least one airway analysis"
    )
  )
  write_tsv(table2, "Table_2_86_gene_evidence_summary.tsv")

  x <- read.delim("results/bfff_full_86_evidence_matrix_v1/Table_full_86_cross_context_evidence.tsv")
  bridge <- x$blood_GSE77087_genome_fdr_support & x$airway_any_fdr_support
  table3 <- x[bridge, c(
    "gene_symbol", "is_representative", "axis", "GSE38900_log2FC", "GSE38900_fdr_bh",
    "GSE77087_log2FC", "GSE77087_fdr_bh", "GSE103842_log2FC", "GSE103842_fdr_bh",
    "GSE97742_RSV_minus_hRV_delta", "GSE97742_fdr_bh_86_family",
    "GSE41374_logFC_shifted_log2", "GSE41374_fdr_bh_86_family_shifted_log2",
    "airway_GSE97742_fdr_support", "airway_GSE41374_fdr_support", "n_external_contexts_fdr")]
  names(table3) <- c(
    "Gene", "Representative_node", "Representative_axis", "GSE38900_log2FC", "GSE38900_genome_FDR",
    "GSE77087_log2FC", "GSE77087_genome_FDR", "GSE103842_log2FC", "GSE103842_genome_FDR",
    "GSE97742_RSV_minus_hRV_recovery_delta", "GSE97742_86_gene_FDR",
    "GSE41374_logFC_shifted_log2", "GSE41374_86_gene_FDR_shifted_log2",
    "GSE97742_support", "GSE41374_robust_support", "External_context_support_count")
  write_tsv(table3, "Table_3_16_bridge_candidates.tsv")

  write.table(read.delim("results/single_gene_upgrade_v1/tables/Table_supplement_all_86_genes.tsv"),
              file.path(supplement_root, "Supplementary_Table_1_complete_blood_effects_meta.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write.table(read.delim("results/barrier_inflammation_repair_v1/three_axis_camera_summary.tsv"),
              file.path(supplement_root, "Supplementary_Table_2_fixed_airway_CAMERA.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write.table(x, file.path(supplement_root, "Supplementary_Table_3_complete_cross_context_matrix.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write.table(read.delim("results/single_gene_upgrade_v1/severity/GSE77087_86_gene_ordinal_trend.tsv"),
              file.path(supplement_root, "Supplementary_Table_4A_severity_ordinal_trend.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write.table(read.delim("results/single_gene_upgrade_v1/specificity/GSE155925_86_gene_rsv_vs_other_virus.tsv"),
              file.path(supplement_root, "Supplementary_Table_4B_viral_context.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write.table(read.delim("results/single_gene_upgrade_v1/modules/cohort_module_effects.tsv"),
              file.path(supplement_root, "Supplementary_Table_4C_BloodGen3_modules.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write.table(read.delim("results/single_gene_upgrade_v1/cell_sensitivity/cell_fraction_qc.tsv"),
              file.path(supplement_root, "Supplementary_Table_4D_cell_fraction_QC.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write.table(read.delim("results/revised_main_gse77087/tables/Table2_revised_target_set.tsv"),
              file.path(supplement_root, "Supplementary_Table_4E_complete_target_set_tests.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

main_v6 <- function() {
  save_figure(make_figure_1(), "Figure_1", 13.2, 6.5)
  save_figure(make_figure_2(), "Figure_2", 12.5, 16.2)
  save_figure(make_figure_3(), "Figure_3", 11.2, 9.5)
  save_figure(make_figure_4(), "Figure_4", 14.2, 13.0)
  save_figure(make_figure_5(), "Figure_5", 12.5, 7.2)
  make_tables_v6()
  capture.output(sessionInfo(), file = file.path(asset_root, "figure_table_session_info.txt"))
  message("V6 figures and TSV tables created under ", asset_root)
}

if (sys.nframe() == 0L) main_v6()
