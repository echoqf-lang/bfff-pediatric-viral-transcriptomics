#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

# Reuse the validated V6 data assembly and plotting functions. This script
# changes typography only and writes separate review copies so the current
# submission figures remain recoverable.
source(file.path("analysis", "R", "27_make_v6_bmc_figures_tables.R"), local = FALSE)

pdf_root <- file.path("output", "pdf")
preview_root <- file.path("output", "figures")
dir.create(pdf_root, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_root, recursive = TRUE, showWarnings = FALSE)

save_submission_figure <- function(plot, stem, width, height) {
  ggsave(
    file.path(pdf_root, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    file.path(preview_root, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 360,
    bg = "white",
    limitsize = FALSE
  )
}

figure_2 <- make_figure_2()
# The lower panel has no long gene-label column. Add a deliberate left inset so
# its plotting region sits beneath the upper heatmap body rather than extending
# into the upper panel's label gutter.
figure_2[[2]] <- figure_2[[2]] +
  scale_x_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.5, 1),
    labels = scales::percent,
    expand = expansion(mult = c(0, 0))
  ) +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 92))
figure_2 <- figure_2 &
  theme(
    axis.text.y = element_text(size = 11.2),
    axis.text.x = element_text(size = 11.5),
    strip.text = element_text(size = 12.5, face = "bold"),
    legend.title = element_text(size = 11.5, face = "bold"),
    legend.text = element_text(size = 10.8)
  )

figure_3 <- strip_figure_header(make_figure_3()) +
  theme(
    axis.text.y = element_text(size = 11.2),
    axis.text.x = element_text(size = 11.5),
    strip.text = element_text(size = 12.0, face = "bold"),
    legend.title = element_text(size = 11.5, face = "bold"),
    legend.text = element_text(size = 10.8)
  )

figure_4 <- strip_figure_header(make_figure_4()) +
  theme(
    axis.text.y = element_text(size = 10.8),
    axis.text.x = element_text(size = 10.8, angle = 25, hjust = 1),
    strip.text = element_text(size = 12.0, face = "bold"),
    legend.title = element_text(size = 11.2, face = "bold"),
    legend.text = element_text(size = 10.5)
  )

save_submission_figure(
  figure_2,
  "Figure_2_submission_larger_text",
  12.5,
  15.7
)
save_submission_figure(
  figure_3,
  "Figure_3_submission_larger_text",
  11.2,
  9.2
)
save_submission_figure(
  figure_4,
  "Figure_4_submission_larger_text",
  14.2,
  12.6
)

message("Larger-text Figure 2-4 review copies created under output/.")
