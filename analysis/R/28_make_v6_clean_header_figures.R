#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

# Reuse the validated V6 data assembly and plotting functions. This script
# writes clean-header submission copies without overwriting the existing assets.
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

save_submission_figure(
  make_figure_2(),
  "Figure_2_submission_revised",
  12.5,
  15.7
)
save_submission_figure(
  strip_figure_header(make_figure_3()),
  "Figure_3_submission_revised",
  11.2,
  9.2
)
save_submission_figure(
  strip_figure_header(make_figure_4()),
  "Figure_4_submission_revised",
  14.2,
  12.6
)

message("Clean-header Figure 2-4 submission copies created under output/.")
