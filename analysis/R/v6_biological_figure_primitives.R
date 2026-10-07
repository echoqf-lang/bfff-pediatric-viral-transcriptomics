#!/usr/bin/env Rscript

# Reusable vector primitives for the V6 biological figures. Anchors use the
# active parent viewport's npc coordinates; compact icons use square snpc
# viewports so their local geometry is not distorted by the page aspect ratio.

library(grid)

v6_palette <- list(
  navy = "#17365D", blue = "#4C87AA", cyan = "#79AFC2",
  red = "#C96561", orange = "#D99555", green = "#82A97E",
  purple = "#8E72A3", skin = "#EFC5A6", grey = "#667482",
  light_blue = "#E6F0F5", light_red = "#F8E5E2",
  light_green = "#E8F1E6", light_gold = "#FFF6DF",
  herb_fill = "#82A97E", herb_light = "#E8F1E6",
  epithelium_fill = "#F7FBFD", epithelium_border = "#79AFC2",
  epithelial_nucleus = "#8E72A3", blood_cell_border = "#B64D4A",
  interferon_fill = "#EEEAF7", evidence_gold = "#B7791F"
)

npc_unit <- function(value) unit(value, "npc")

assert_finite_scalar <- function(value, name, lower = -Inf, upper = Inf,
                                 lower_inclusive = TRUE) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop(sprintf("%s must be one finite, non-missing numeric value", name), call. = FALSE)
  }
  lower_ok <- if (lower_inclusive) value >= lower else value > lower
  if (!lower_ok || value > upper) {
    interval <- paste0(if (lower_inclusive) "[" else "(", lower, ", ", upper, "]")
    stop(sprintf("%s must lie in %s", name, interval), call. = FALSE)
  }
  invisible(TRUE)
}

assert_text <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("%s must be one non-missing character value", name), call. = FALSE)
  }
  invisible(TRUE)
}

assert_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("%s must be one non-missing logical value", name), call. = FALSE)
  }
  invisible(TRUE)
}

assert_anchor <- function(x, y) {
  assert_finite_scalar(x, "x")
  assert_finite_scalar(y, "y")
  invisible(TRUE)
}

assert_scale <- function(scale) {
  assert_finite_scalar(scale, "scale", lower = 0, lower_inclusive = FALSE)
  invisible(TRUE)
}

assert_size <- function(value, name) {
  assert_finite_scalar(value, name, lower = 0, lower_inclusive = FALSE)
  invisible(TRUE)
}

push_icon_viewport <- function(x, y, scale, base_size) {
  pushViewport(viewport(
    x = npc_unit(x), y = npc_unit(y),
    width = unit(base_size * scale, "snpc"),
    height = unit(base_size * scale, "snpc"),
    just = "centre"
  ))
}

icon_lwd <- function(base, scale) base * scale
icon_radius <- function(base, scale) unit(base * scale, "mm")

draw_label <- function(label, x, y, size = 10, colour = v6_palette$navy,
                       fontface = "plain", just = "centre") {
  assert_text(label, "label")
  assert_anchor(x, y)
  assert_size(size, "size")
  assert_text(colour, "colour")
  assert_text(fontface, "fontface")
  assert_text(just, "just")
  grid.text(label, x = npc_unit(x), y = npc_unit(y), just = just,
            gp = gpar(fontfamily = "sans", fontsize = size, col = colour,
                      fontface = fontface))
  invisible(TRUE)
}

draw_connector <- function(x0, y0, x1, y1, dashed = FALSE, arrow = FALSE,
                           colour = v6_palette$navy, width = 1.5) {
  assert_anchor(x0, y0)
  assert_anchor(x1, y1)
  assert_flag(dashed, "dashed")
  assert_flag(arrow, "arrow")
  assert_text(colour, "colour")
  assert_size(width, "width")
  grid.lines(npc_unit(c(x0, x1)), npc_unit(c(y0, y1)),
             arrow = if (arrow) grid::arrow(length = unit(2.2, "mm")) else NULL,
             gp = gpar(col = colour, lwd = width, lty = if (dashed) 2 else 1))
  invisible(TRUE)
}

draw_herb_cluster <- function(x, y, scale = 1) {
  assert_anchor(x, y)
  assert_scale(scale)
  push_icon_viewport(x, y, scale, base_size = 0.105)
  on.exit(popViewport(), add = TRUE)
  offsets <- matrix(c(-0.23, 0.08, 0.00, 0.28, 0.23, 0.08,
                      -0.13, -0.20, 0.13, -0.20), ncol = 2, byrow = TRUE)
  grid.lines(npc_unit(c(0.5, 0.5)), npc_unit(c(0.12, 0.83)),
             gp = gpar(col = v6_palette$herb_fill, lwd = icon_lwd(1.7, scale)))
  for (i in seq_len(nrow(offsets))) {
    leaf_x <- 0.5 + offsets[i, 1]
    leaf_y <- 0.5 + offsets[i, 2]
    grid.roundrect(npc_unit(leaf_x), npc_unit(leaf_y),
                   width = npc_unit(0.27), height = npc_unit(0.15),
                   r = icon_radius(1.5, scale),
                   gp = gpar(fill = if (i %% 2L) v6_palette$herb_light else v6_palette$herb_fill,
                             col = v6_palette$herb_fill, lwd = icon_lwd(0.65, scale)))
    grid.lines(npc_unit(c(0.5, leaf_x)), npc_unit(c(0.41, leaf_y)),
               gp = gpar(col = v6_palette$herb_fill, lwd = icon_lwd(0.65, scale)))
  }
  invisible(TRUE)
}

draw_child_torso <- function(x, y, scale = 1) {
  assert_anchor(x, y)
  assert_scale(scale)
  push_icon_viewport(x, y, scale, base_size = 0.115)
  on.exit(popViewport(), add = TRUE)
  grid.circle(npc_unit(0.5), npc_unit(0.79), r = npc_unit(0.16),
              gp = gpar(fill = v6_palette$skin, col = v6_palette$navy,
                        lwd = icon_lwd(1, scale)))
  grid.roundrect(npc_unit(0.5), npc_unit(0.43), width = npc_unit(0.57), height = npc_unit(0.54),
                 r = icon_radius(2, scale),
                 gp = gpar(fill = v6_palette$light_blue, col = v6_palette$navy,
                           lwd = icon_lwd(1.2, scale)))
  grid.lines(npc_unit(c(0.22, 0.05)), npc_unit(c(0.53, 0.31)),
             gp = gpar(col = v6_palette$navy, lwd = icon_lwd(1.2, scale)))
  grid.lines(npc_unit(c(0.78, 0.95)), npc_unit(c(0.53, 0.31)),
             gp = gpar(col = v6_palette$navy, lwd = icon_lwd(1.2, scale)))
  grid.lines(npc_unit(c(0.36, 0.64)), npc_unit(c(0.46, 0.46)),
             gp = gpar(col = v6_palette$cyan, lwd = icon_lwd(1.2, scale)))
  invisible(TRUE)
}

draw_blood_drop <- function(x, y, scale = 1) {
  assert_anchor(x, y)
  assert_scale(scale)
  push_icon_viewport(x, y, scale, base_size = 0.090)
  on.exit(popViewport(), add = TRUE)
  grid.polygon(npc_unit(c(0.5, 0.18, 0.16, 0.5, 0.84, 0.82)),
               npc_unit(c(0.93, 0.56, 0.30, 0.08, 0.30, 0.56)),
               gp = gpar(fill = v6_palette$light_red, col = v6_palette$red,
                         lwd = icon_lwd(1.3, scale)))
  grid.circle(npc_unit(0.5), npc_unit(0.37), r = npc_unit(0.13),
              gp = gpar(fill = v6_palette$red, col = NA))
  invisible(TRUE)
}

draw_upper_airway_icon <- function(x, y, scale = 1) {
  assert_anchor(x, y)
  assert_scale(scale)
  push_icon_viewport(x, y, scale, base_size = 0.100)
  on.exit(popViewport(), add = TRUE)
  grid.lines(npc_unit(c(0.18, 0.40, 0.70, 0.76, 0.56)),
             npc_unit(c(0.73, 0.81, 0.66, 0.42, 0.12)),
             gp = gpar(col = v6_palette$navy, lwd = icon_lwd(1.5, scale)))
  grid.lines(npc_unit(c(0.56, 0.40, 0.57)), npc_unit(c(0.54, 0.43, 0.30)),
             gp = gpar(col = v6_palette$cyan, lwd = icon_lwd(2, scale)))
  grid.lines(npc_unit(c(0.57, 0.57)), npc_unit(c(0.30, 0.08)),
             gp = gpar(col = v6_palette$cyan, lwd = icon_lwd(2, scale)))
  grid.lines(npc_unit(c(0.43, 0.72)), npc_unit(c(0.21, 0.21)),
             gp = gpar(col = v6_palette$cyan, lwd = icon_lwd(1, scale)))
  invisible(TRUE)
}

draw_airway_epithelium <- function(x, y, width, height, damaged_fraction = 0.45) {
  assert_anchor(x, y)
  assert_size(width, "width")
  assert_size(height, "height")
  assert_finite_scalar(damaged_fraction, "damaged_fraction", lower = 0, upper = 1)
  cell_count <- 12L
  cell_width <- width / cell_count
  damaged_count <- as.integer(round(cell_count * damaged_fraction))
  damaged_start <- cell_count - damaged_count + 1L
  grid.rect(npc_unit(x + width / 2), npc_unit(y + height / 2),
            width = npc_unit(width), height = npc_unit(height),
            gp = gpar(fill = v6_palette$light_blue, col = v6_palette$blue, lwd = 1))
  for (i in seq_len(cell_count)) {
    centre_x <- x + (i - 0.5) * cell_width
    damaged <- damaged_count > 0L && i >= damaged_start
    grid.roundrect(npc_unit(centre_x), npc_unit(y + height * 0.47),
                   width = npc_unit(cell_width * 0.88), height = npc_unit(height * 0.82),
                   r = unit(1.1, "mm"),
                   gp = gpar(fill = v6_palette$epithelium_fill,
                             col = v6_palette$epithelium_border, lwd = 0.65))
    grid.circle(npc_unit(centre_x), npc_unit(y + height * 0.36),
                r = unit(min(cell_width * 0.17, height * 0.085), "snpc"),
                gp = gpar(fill = v6_palette$epithelial_nucleus, col = NA))
    cilia_count <- if (damaged) 1L else 3L
    cilia_height <- if (damaged) height * 0.10 else height * 0.23
    cilia_x <- if (cilia_count == 1L) centre_x else {
      centre_x + seq(-0.22, 0.22, length.out = cilia_count) * cell_width
    }
    for (tip_x in cilia_x) {
      grid.lines(npc_unit(c(tip_x, tip_x)),
                 npc_unit(c(y + height * 0.88, y + height * 0.88 + cilia_height)),
                 gp = gpar(col = v6_palette$navy, lwd = 0.65))
    }
  }
  grid.lines(npc_unit(c(x, x + width)), npc_unit(c(y + height * 0.06, y + height * 0.06)),
             gp = gpar(col = v6_palette$navy, lwd = 1.2))
  invisible(TRUE)
}

draw_immune_cell <- function(x, y, scale = 1,
                             type = c("myeloid", "neutrophil", "interferon"),
                             show_label = TRUE) {
  assert_anchor(x, y)
  assert_scale(scale)
  if (!is.character(type) || length(type) < 1L || anyNA(type)) {
    stop("type must be a non-missing immune-cell type", call. = FALSE)
  }
  if (!is.logical(show_label) || length(show_label) != 1L || is.na(show_label)) {
    stop("show_label must be one non-missing logical value", call. = FALSE)
  }
  type <- match.arg(type)
  cell_style <- switch(
    type,
    myeloid = list(fill = v6_palette$orange, nucleus = v6_palette$navy, symbol = "M"),
    neutrophil = list(fill = v6_palette$light_red, nucleus = v6_palette$red, symbol = "N"),
    interferon = list(fill = v6_palette$interferon_fill, nucleus = v6_palette$purple, symbol = "IFN")
  )
  push_icon_viewport(x, y, scale, base_size = 0.062)
  on.exit(popViewport(), add = TRUE)
  grid.circle(npc_unit(0.5), npc_unit(0.5), r = npc_unit(0.38),
              gp = gpar(fill = cell_style$fill, col = cell_style$nucleus,
                        lwd = icon_lwd(1.1, scale)))
  if (type == "neutrophil") {
    for (offset in c(-0.12, 0, 0.12)) {
      grid.circle(npc_unit(0.5 + offset), npc_unit(0.5), r = npc_unit(0.09),
                  gp = gpar(fill = cell_style$nucleus, col = NA))
    }
  } else if (show_label) {
    draw_label(cell_style$symbol, 0.5, 0.5,
               size = (if (type == "interferon") 6 else 8) * scale,
               colour = cell_style$nucleus, fontface = "bold")
  }
  invisible(TRUE)
}

draw_repair_cells <- function(x, y, scale = 1) {
  assert_anchor(x, y)
  assert_scale(scale)
  push_icon_viewport(x, y, scale, base_size = 0.090)
  on.exit(popViewport(), add = TRUE)
  for (centre_x in c(0.30, 0.70)) {
    grid.roundrect(npc_unit(centre_x), npc_unit(0.45), width = npc_unit(0.30), height = npc_unit(0.43),
                   r = icon_radius(1.1, scale),
                   gp = gpar(fill = v6_palette$light_green, col = v6_palette$green,
                             lwd = icon_lwd(1, scale)))
    grid.circle(npc_unit(centre_x), npc_unit(0.38), r = npc_unit(0.06),
                gp = gpar(fill = v6_palette$epithelial_nucleus, col = NA))
  }
  grid.lines(npc_unit(c(0.43, 0.57)), npc_unit(c(0.82, 0.82)),
             arrow = grid::arrow(ends = "both", length = icon_radius(1.5, scale)),
             gp = gpar(col = v6_palette$green, lwd = icon_lwd(1, scale)))
  invisible(TRUE)
}

draw_blood_vessel <- function(x, y, width, height) {
  assert_anchor(x, y)
  assert_size(width, "width")
  assert_size(height, "height")
  grid.roundrect(npc_unit(x + width / 2), npc_unit(y + height / 2),
                 width = npc_unit(width), height = npc_unit(height), r = unit(2, "mm"),
                 gp = gpar(fill = v6_palette$light_red, col = v6_palette$red, lwd = 1.2))
  cell_x <- seq(x + width * 0.12, x + width * 0.88, length.out = 6)
  for (centre_x in cell_x) {
    grid.circle(npc_unit(centre_x), npc_unit(y + height / 2),
                r = unit(min(width * 0.0375, height * 0.24), "snpc"),
                gp = gpar(fill = v6_palette$red, col = v6_palette$blood_cell_border, lwd = 0.4))
  }
  invisible(TRUE)
}

draw_count_badge <- function(value, label, x, y, fill, border) {
  assert_text(value, "value")
  assert_text(label, "label")
  assert_anchor(x, y)
  assert_text(fill, "fill")
  assert_text(border, "border")
  grid.roundrect(npc_unit(x), npc_unit(y), width = npc_unit(0.093), height = npc_unit(0.062),
                 r = unit(2, "mm"), gp = gpar(fill = fill, col = border, lwd = 1.2))
  draw_label(value, x, y + 0.010, size = 11, colour = border, fontface = "bold")
  draw_label(label, x, y - 0.013, size = 6.8, colour = v6_palette$grey)
  invisible(TRUE)
}

draw_evidence_boundary <- function(label, x, y, width, height) {
  assert_text(label, "label")
  assert_anchor(x, y)
  assert_size(width, "width")
  assert_size(height, "height")
  grid.roundrect(npc_unit(x + width / 2), npc_unit(y + height / 2),
                 width = npc_unit(width), height = npc_unit(height), r = unit(2, "mm"),
                 gp = gpar(fill = v6_palette$light_gold, col = v6_palette$evidence_gold,
                           lwd = 1.2, lty = 2))
  draw_label(label, x + width / 2, y + height / 2, size = 8.5,
             colour = v6_palette$evidence_gold, fontface = "bold")
  invisible(TRUE)
}
