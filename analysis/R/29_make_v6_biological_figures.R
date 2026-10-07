#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

fail <- function(message) stop(message, call. = FALSE)

parse_figure_request <- function(args) {
  unknown <- args[!grepl("^--figure=", args)]
  if (length(unknown)) fail(paste0("Unknown argument(s): ", paste(unknown, collapse = ", ")))
  figure_arg <- grep("^--figure=", args, value = TRUE)
  if (length(figure_arg) != 1L) fail("Specify exactly one of --figure=1, --figure=5, or --figure=all")
  requested <- sub("^--figure=", "", figure_arg)
  if (!requested %in% c("1", "5", "all")) {
    fail("--figure must be one of: 1, 5, all")
  }
  requested
}

# Fail unsupported requests before reading any study input or creating output.
requested_figure <- parse_figure_request(commandArgs(trailingOnly = TRUE))
if (!capabilities("cairo")) fail("Cairo graphics support is required")

# Fontconfig needs a writable cache in restricted execution environments.
font_cache_root <- file.path(tempdir(), "v6-biological-figures-font-cache")
fontconfig_cache <- file.path(font_cache_root, "fontconfig")
if (!dir.exists(fontconfig_cache) &&
    !dir.create(fontconfig_cache, recursive = TRUE, showWarnings = FALSE)) {
  fail(paste0("Unable to create writable Fontconfig cache: ", fontconfig_cache))
}
Sys.setenv(XDG_CACHE_HOME = font_cache_root)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) fail("Unable to determine figure script location")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
project_path <- function(...) file.path(project_root, ...)

# Run the frozen input contract before creating an output directory or opening
# a graphics device. The validated objects are retained in a private
# environment so key figure counts come from the real inputs, not plot labels.
contract <- new.env(parent = globalenv())
sys.source(
  project_path("analysis", "R", "28_validate_v6_biological_figure_inputs.R"),
  envir = contract
)

figure_counts <- list(
  frozen = as.integer(nrow(contract$frozen)),
  airway = as.integer(contract$airway),
  bridge = as.integer(contract$bridge)
)

required_provenance_counts <- c("herbs", "putative", "disease", "search_space")
if (!exists("provenance_counts", envir = contract, inherits = FALSE) ||
    !identical(names(contract$provenance_counts), required_provenance_counts)) {
  fail("Input contract did not provide the required table-driven provenance counts")
}
provenance_counts <- contract$provenance_counts

source(project_path("analysis", "R", "v6_biological_figure_primitives.R"), local = FALSE)

output_dir <- project_path(
  "manuscript", "v6_bmc_submission_assets", "figures_biological"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# The 13.2-inch source is scaled to approximately 177 mm at publication.
# A 13.5-point source label therefore remains just above 7 points after scale.
figure_type <- list(
  title = 24,
  subtitle = 14,
  panel = 17,
  primary = 14.5,
  body = 13.5,
  count = 18,
  boundary = 14.5
)

render_device <- function(open_device, draw_fun) {
  open_device()
  device_id <- grDevices::dev.cur()
  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && device_id %in% open_devices) {
      grDevices::dev.off(device_id)
    }
  }, add = TRUE)
  draw_fun()
  grDevices::dev.off(device_id)
  invisible(TRUE)
}

render_pair <- function(stem, draw_fun, width = 13.2, height = 7.2, dpi = 360) {
  if (!is.character(stem) || length(stem) != 1L || is.na(stem) || !nzchar(stem)) {
    fail("stem must be one non-empty character value")
  }
  if (!is.function(draw_fun)) fail("draw_fun must be a function")
  for (item in list(width = width, height = height, dpi = dpi)) {
    if (!is.numeric(item) || length(item) != 1L || is.na(item) ||
        !is.finite(item) || item <= 0) {
      fail("width, height and dpi must be positive finite numeric values")
    }
  }
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  png_path <- file.path(output_dir, paste0(stem, ".png"))

  render_device(
    function() grDevices::cairo_pdf(
      filename = pdf_path, width = width, height = height, family = "sans"
    ),
    draw_fun
  )
  render_device(
    function() grDevices::png(
      filename = png_path, width = width, height = height, units = "in",
      res = dpi, type = "cairo", bg = "white"
    ),
    draw_fun
  )

  if (nzchar(Sys.which("sips"))) {
    sips_output <- system2(
      "sips",
      c("-s", "dpiWidth", as.character(dpi),
        "-s", "dpiHeight", as.character(dpi), png_path),
      stdout = TRUE,
      stderr = TRUE
    )
    sips_status <- attr(sips_output, "status")
    if (!is.null(sips_status) && sips_status != 0L) {
      fail(paste0("Unable to set PNG DPI metadata: ", paste(sips_output, collapse = " ")))
    }
  }

  generated <- c(pdf_path, png_path)
  if (!all(file.exists(generated)) || any(file.info(generated)$size <= 0)) {
    fail(paste0("Figure export failed for ", stem))
  }
  cat(sprintf(
    "BIOLOGICAL_FIGURE_RENDERED stem=%s pdf=%s png=%s\n",
    stem, pdf_path, png_path
  ))
  invisible(generated)
}

draw_panel <- function(x, y, width, height, fill) {
  grid::grid.roundrect(
    x = unit(x + width / 2, "npc"),
    y = unit(y + height / 2, "npc"),
    width = unit(width, "npc"),
    height = unit(height, "npc"),
    r = unit(2.5, "mm"),
    gp = grid::gpar(fill = fill, col = "#D7E0E7", lwd = 0.9)
  )
  invisible(TRUE)
}

draw_metric_badge <- function(value, label, x, y, fill, border) {
  grid::grid.roundrect(
    unit(x, "npc"), unit(y, "npc"),
    width = unit(0.093, "npc"), height = unit(0.130, "npc"),
    r = unit(2, "mm"),
    gp = grid::gpar(fill = fill, col = border, lwd = 1.2)
  )
  draw_label(as.character(value), x, y + 0.028, size = figure_type$count,
             colour = border, fontface = "bold")
  draw_label(label, x, y - 0.025, size = figure_type$body,
             colour = v6_palette$grey)
  invisible(TRUE)
}

draw_figure_boundary <- function(label, x, y, width, height) {
  grid::grid.roundrect(
    unit(x + width / 2, "npc"), unit(y + height / 2, "npc"),
    width = unit(width, "npc"), height = unit(height, "npc"),
    r = unit(2, "mm"),
    gp = grid::gpar(fill = v6_palette$light_gold,
                    col = v6_palette$evidence_gold, lwd = 1.2, lty = 2)
  )
  draw_label(label, x + width / 2, y + height / 2,
             size = figure_type$boundary,
             colour = v6_palette$evidence_gold, fontface = "bold")
  invisible(TRUE)
}

draw_virus_particle <- function(x, y, scale = 1, colour = v6_palette$purple) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1 ||
      !is.numeric(y) || length(y) != 1L || !is.finite(y) || y < 0 || y > 1 ||
      !is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    fail("Virus-particle coordinates and scale must be valid finite scalars")
  }
  radius <- 0.012 * scale
  for (angle in seq(0, 2 * pi - pi / 4, by = pi / 4)) {
    x0 <- x + cos(angle) * radius
    y0 <- y + sin(angle) * radius
    x1 <- x + cos(angle) * radius * 1.55
    y1 <- y + sin(angle) * radius * 1.55
    grid::grid.lines(unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
                     gp = grid::gpar(col = colour, lwd = 1.0))
    grid::grid.circle(unit(x1, "npc"), unit(y1, "npc"),
                      r = unit(0.0026 * scale, "snpc"),
                      gp = grid::gpar(fill = colour, col = NA))
  }
  grid::grid.circle(unit(x, "npc"), unit(y, "npc"),
                    r = unit(radius, "snpc"),
                    gp = grid::gpar(fill = "#EEEAF7", col = colour, lwd = 1.2))
  grid::grid.circle(unit(x - radius * 0.30, "npc"), unit(y + radius * 0.18, "npc"),
                    r = unit(radius * 0.15, "snpc"),
                    gp = grid::gpar(fill = colour, col = NA))
  grid::grid.circle(unit(x + radius * 0.28, "npc"), unit(y - radius * 0.20, "npc"),
                    r = unit(radius * 0.12, "snpc"),
                    gp = grid::gpar(fill = colour, col = NA))
  invisible(TRUE)
}

draw_program_tag <- function(label, x, y, width, fill, border) {
  grid::grid.roundrect(
    unit(x + width / 2, "npc"), unit(y, "npc"),
    width = unit(width, "npc"), height = unit(0.060, "npc"),
    r = unit(1.8, "mm"),
    gp = grid::gpar(fill = fill, col = border, lwd = 1.0)
  )
  draw_label(label, x + width / 2, y, size = figure_type$body,
             colour = border, fontface = "bold")
  invisible(TRUE)
}

draw_figure_1 <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, 1), yscale = c(0, 1)))
  on.exit(grid::popViewport(), add = TRUE)
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  draw_label(
    "Multi-cohort integration across pediatric blood and airway",
    0.025, 0.956, size = figure_type$title, fontface = "bold", just = "left"
  )
  draw_label(
    "Candidate provenance, compartment-stratified evidence, and cross-context synthesis",
    0.025, 0.908, size = figure_type$subtitle,
    colour = v6_palette$grey, just = "left"
  )

  draw_panel(0.020, 0.160, 0.390, 0.725, "#F8FBF7")
  draw_panel(0.425, 0.160, 0.295, 0.725, "#F7FAFC")
  draw_panel(0.735, 0.160, 0.245, 0.725, "#FCF9FD")

  draw_label("A  Candidate-space provenance", 0.038, 0.848,
             size = figure_type$panel, fontface = "bold", just = "left")
  draw_label("B  Cohorts by compartment", 0.442, 0.848,
             size = figure_type$panel, fontface = "bold", just = "left")
  draw_label("C  Cross-context synthesis", 0.752, 0.848,
             size = figure_type$panel, fontface = "bold", just = "left")

  # Panel A: connectors denote provenance and analytic selection, not binding.
  draw_connector(0.105, 0.655, 0.118, 0.710, arrow = TRUE,
                 colour = v6_palette$grey, width = 1.1)
  draw_connector(0.212, 0.700, 0.219, 0.650, arrow = TRUE,
                 colour = v6_palette$grey, width = 1.1)
  draw_connector(0.212, 0.550, 0.219, 0.610, arrow = TRUE,
                 colour = v6_palette$grey, width = 1.1)
  draw_connector(0.312, 0.625, 0.315, 0.625, arrow = TRUE,
                 colour = v6_palette$grey, width = 1.1)

  draw_herb_cluster(0.073, 0.640, scale = 1.05)
  draw_label(sprintf("%d herbs", provenance_counts[["herbs"]]),
             0.073, 0.535, size = figure_type$body, colour = v6_palette$green,
             fontface = "bold")
  draw_metric_badge(
    format(provenance_counts[["putative"]], big.mark = ","), "putative\ntargets",
    0.165, 0.710, v6_palette$light_blue, v6_palette$blue
  )
  draw_metric_badge(
    format(provenance_counts[["disease"]], big.mark = ","), "disease\ngenes",
    0.165, 0.535, "#F3EAF1", v6_palette$purple
  )
  draw_metric_badge(
    format(provenance_counts[["search_space"]], big.mark = ","), "search\nspace",
    0.265, 0.625, v6_palette$light_gold, v6_palette$evidence_gold
  )
  draw_metric_badge(
    format(figure_counts$frozen, big.mark = ","), "outcome\nselected",
    0.362, 0.625, "#EEEAF7", v6_palette$purple
  )
  draw_label("Database\nintersection", 0.238, 0.425,
             size = figure_type$body, colour = v6_palette$grey)
  draw_label("GSE38900\nselection", 0.330, 0.345,
             size = figure_type$body, colour = v6_palette$grey)
  draw_label("Not confirmed\npharmacologic targets", 0.215, 0.225,
             size = figure_type$body, colour = v6_palette$evidence_gold,
             fontface = "bold")

  # Panel B: blood and airway are visibly separate analytic compartments.
  draw_child_torso(0.472, 0.620, scale = 1.08)
  draw_label("Pediatric\ncohorts", 0.472, 0.485, size = figure_type$body,
             colour = v6_palette$navy, fontface = "bold")
  draw_connector(0.518, 0.655, 0.532, 0.710, colour = "#A8B5C0", width = 1.0)
  draw_connector(0.518, 0.585, 0.532, 0.455, colour = "#A8B5C0", width = 1.0)
  draw_blood_drop(0.545, 0.720, scale = 0.68)
  draw_upper_airway_icon(0.545, 0.445, scale = 0.74)

  draw_label("Blood · whole blood", 0.575, 0.755, size = figure_type$primary,
             colour = v6_palette$red, fontface = "bold", just = "left")
  draw_label("GSE38900\nGSE77087\nGSE103842", 0.575, 0.680,
             size = figure_type$body,
             colour = v6_palette$navy, just = "left")
  draw_label("Viral context\nGSE155925", 0.575, 0.585,
             size = figure_type$body,
             colour = v6_palette$grey, just = "left")

  grid::grid.lines(unit(c(0.530, 0.700), "npc"), unit(c(0.525, 0.525), "npc"),
                   gp = grid::gpar(col = "#D7E0E7", lwd = 0.9))

  draw_label("Upper airway", 0.575, 0.495, size = figure_type$primary,
             colour = v6_palette$blue, fontface = "bold", just = "left")
  draw_label("GSE97742\nGSE41374", 0.575, 0.425, size = figure_type$body,
             colour = v6_palette$navy, just = "left")
  draw_label("NP / nasal samples", 0.575, 0.345, size = figure_type$body,
             colour = v6_palette$grey, just = "left")
  draw_label("Compartments\nanalyzed separately", 0.573, 0.250,
             size = figure_type$body,
             colour = v6_palette$evidence_gold, fontface = "bold")

  # Neutral workflow connection: evidence is integrated without merging matrices.
  draw_connector(0.704, 0.355, 0.744, 0.355, arrow = TRUE,
                 colour = v6_palette$grey, width = 1.1)

  # Panel C: biological programs are transcriptomic summaries, not causal paths.
  draw_airway_epithelium(0.755, 0.680, 0.060, 0.055, damaged_fraction = 1)
  draw_label("Ciliary\nprograms ↓", 0.835, 0.735, size = figure_type$primary,
             colour = v6_palette$blue, fontface = "bold", just = "left")
  draw_label("airway-localized", 0.835, 0.650, size = figure_type$body,
             colour = v6_palette$grey, just = "left")

  draw_immune_cell(0.770, 0.545, scale = 0.80, type = "myeloid",
                   show_label = FALSE)
  draw_label("M", 0.770, 0.545, size = figure_type$body,
             colour = v6_palette$navy, fontface = "bold")
  draw_immune_cell(0.807, 0.545, scale = 0.80, type = "neutrophil")
  draw_label("Inflammatory\nprograms ↑", 0.835, 0.565,
             size = figure_type$primary,
             colour = v6_palette$red, fontface = "bold", just = "left")
  draw_label("shared viral\nresponse", 0.835, 0.475, size = figure_type$body,
             colour = v6_palette$grey, just = "left")

  draw_repair_cells(0.785, 0.380, scale = 0.82)
  draw_label("Repair / turnover", 0.835, 0.405, size = figure_type$primary,
             colour = v6_palette$green, fontface = "bold", just = "left")
  draw_label("phase-dependent", 0.835, 0.345, size = figure_type$body,
             colour = v6_palette$green, fontface = "bold", just = "left")
  draw_metric_badge(
    format(figure_counts$airway, big.mark = ","), "airway-\nsupported",
    0.807, 0.245, v6_palette$light_gold, v6_palette$evidence_gold
  )
  draw_metric_badge(
    format(figure_counts$bridge, big.mark = ","), "blood–airway",
    0.925, 0.245, v6_palette$light_blue, v6_palette$blue
  )

  draw_figure_boundary(
    "No cohort included BFFF exposure — intervention effects remain untested",
    0.030, 0.048, 0.940, 0.073
  )
  invisible(TRUE)
}

draw_figure_5 <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, 1), yscale = c(0, 1)))
  on.exit(grid::popViewport(), add = TRUE)
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  draw_label(
    "A ciliary–inflammatory core with phase-dependent repair",
    0.025, 0.956, size = figure_type$title, fontface = "bold", just = "left"
  )
  draw_label(
    "Cross-compartment organization of pediatric RSV/hRV host-response evidence",
    0.025, 0.908, size = figure_type$subtitle,
    colour = v6_palette$grey, just = "left"
  )

  draw_panel(0.020, 0.145, 0.700, 0.735, "#F7FAFC")
  draw_panel(0.735, 0.145, 0.245, 0.735, "#FCF9FD")
  draw_label("Observed upper-airway transcriptomic organization", 0.038, 0.846,
             size = figure_type$panel, fontface = "bold", just = "left")
  draw_label("Cross-compartment\nprioritization", 0.752, 0.828,
             size = figure_type$panel, fontface = "bold", just = "left")

  # Luminal context: generic particles represent viral context only, without
  # receptor binding or virus-specific entry mechanisms.
  draw_upper_airway_icon(0.065, 0.788, scale = 0.70)
  draw_label("Nasal / upper-airway context", 0.100, 0.805,
             size = figure_type$primary, colour = v6_palette$blue,
             fontface = "bold", just = "left")
  draw_virus_particle(0.500, 0.790, scale = 0.92, colour = v6_palette$purple)
  draw_virus_particle(0.610, 0.790, scale = 0.92, colour = v6_palette$red)
  draw_label("RSV", 0.500, 0.748, size = figure_type$body,
             colour = v6_palette$purple, fontface = "bold")
  draw_label("hRV", 0.610, 0.748, size = figure_type$body,
             colour = v6_palette$red, fontface = "bold")

  # The epithelial drawing is a visual mapping of transcriptomic programs,
  # not microscopy. Reduced/shortened cilia occupy only part of the layer.
  draw_airway_epithelium(0.065, 0.575, 0.610, 0.135, damaged_fraction = 0.45)
  draw_program_tag("Cilium organization / movement ↓", 0.235, 0.535, 0.270,
                   v6_palette$light_blue, v6_palette$blue)
  draw_label("Transcriptomic program mapped to ciliated epithelium",
             0.370, 0.494, size = figure_type$body,
             colour = v6_palette$grey)

  # Submucosal programs share a field but are deliberately not connected into
  # a directional cilium -> inflammation -> repair causal chain.
  grid::grid.roundrect(
    unit(0.370, "npc"), unit(0.382, "npc"),
    width = unit(0.610, "npc"), height = unit(0.185, "npc"),
    r = unit(2, "mm"),
    gp = grid::gpar(fill = "#FBF7F6", col = "#E7DAD7", lwd = 0.8)
  )
  draw_label("Inflammatory programs ↑", 0.105, 0.444,
             size = figure_type$primary, colour = v6_palette$red,
             fontface = "bold", just = "left")
  draw_immune_cell(0.145, 0.382, scale = 0.92, type = "myeloid",
                   show_label = FALSE)
  draw_immune_cell(0.265, 0.382, scale = 0.92, type = "neutrophil",
                   show_label = FALSE)
  draw_immune_cell(0.385, 0.382, scale = 0.92, type = "interferon",
                   show_label = FALSE)
  draw_label("Myeloid", 0.145, 0.328, size = figure_type$body,
             colour = v6_palette$orange, fontface = "bold")
  draw_label("Neutrophil", 0.265, 0.328, size = figure_type$body,
             colour = v6_palette$red, fontface = "bold")
  draw_label("Interferon", 0.385, 0.328, size = figure_type$body,
             colour = v6_palette$purple, fontface = "bold")

  draw_label("Repair / turnover", 0.535, 0.444,
             size = figure_type$primary, colour = v6_palette$green,
             fontface = "bold")
  draw_repair_cells(0.535, 0.382, scale = 1.02)
  draw_label("Phase-dependent", 0.535, 0.328, size = figure_type$body,
             colour = v6_palette$green, fontface = "bold")
  grid::grid.lines(unit(c(0.500, 0.570), "npc"), unit(c(0.292, 0.292), "npc"),
                   arrow = grid::arrow(ends = "last", length = unit(2.2, "mm")),
                   gp = grid::gpar(col = v6_palette$green, lwd = 1.1))
  draw_label("Acute", 0.485, 0.292, size = figure_type$body,
             colour = v6_palette$green, just = "right")
  draw_label("Discharge", 0.585, 0.292, size = figure_type$body,
             colour = v6_palette$green, just = "left")

  # Systemic blood remains a separate compartment. The two expression-module
  # directions are descriptive and are not linked to a treatment effect.
  draw_label("Systemic blood compartment", 0.065, 0.255,
             size = figure_type$primary, colour = v6_palette$red,
             fontface = "bold", just = "left")
  draw_blood_vessel(0.065, 0.168, 0.315, 0.065)
  draw_label("Interferon–myeloid modules ↑", 0.405, 0.216,
             size = figure_type$body, colour = v6_palette$purple,
             fontface = "bold", just = "left")
  draw_label("Lymphocyte-related modules ↓", 0.405, 0.174,
             size = figure_type$body, colour = v6_palette$blue,
             fontface = "bold", just = "left")

  # Dashed neutral connections organize airway and blood evidence around the
  # cross-compartment set; arrows are intentionally omitted.
  draw_connector(0.675, 0.635, 0.808, 0.620, dashed = TRUE,
                 colour = v6_palette$grey, width = 1.0)
  draw_connector(0.675, 0.205, 0.808, 0.590, dashed = TRUE,
                 colour = v6_palette$grey, width = 1.0)
  draw_metric_badge(
    format(figure_counts$bridge, big.mark = ","), "blood–airway\ncandidates",
    0.855, 0.610, v6_palette$light_blue, v6_palette$blue
  )
  draw_label("Cross-compartment\npriority candidates", 0.855, 0.493,
             size = figure_type$primary, colour = v6_palette$navy,
             fontface = "bold")
  draw_label("Not confirmed\npharmacologic targets", 0.855, 0.410,
             size = figure_type$body, colour = v6_palette$evidence_gold,
             fontface = "bold")

  # BFFF is visually isolated: no line or arrow connects it to any biological
  # process or candidate node.
  grid::grid.roundrect(
    unit(0.855, "npc"), unit(0.258, "npc"),
    width = unit(0.205, "npc"), height = unit(0.195, "npc"),
    r = unit(2.5, "mm"),
    gp = grid::gpar(fill = v6_palette$light_gold,
                    col = v6_palette$evidence_gold, lwd = 1.2, lty = 2)
  )
  draw_label("BFFF intervention layer", 0.855, 0.315,
             size = figure_type$primary, colour = v6_palette$evidence_gold,
             fontface = "bold")
  draw_label("UNTESTED", 0.855, 0.269, size = figure_type$body,
             colour = v6_palette$evidence_gold, fontface = "bold")
  draw_label("Exposure and direction", 0.855, 0.240,
             size = figure_type$body, colour = v6_palette$grey)
  draw_label("Target binding", 0.855, 0.212,
             size = figure_type$body, colour = v6_palette$grey)
  draw_label("require experimental testing", 0.855, 0.184,
             size = figure_type$body, colour = v6_palette$grey)

  draw_label(
    "Connections organize co-occurring evidence; no causal or temporal sequence is implied",
    0.500, 0.067, size = figure_type$body, colour = v6_palette$grey
  )
  invisible(TRUE)
}

if (requested_figure %in% c("1", "all")) {
  render_pair("Figure_1_biological", draw_figure_1)
}
if (requested_figure %in% c("5", "all")) {
  render_pair("Figure_5_biological", draw_figure_5)
}
