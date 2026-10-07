#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

fail <- function(message) stop(message, call. = FALSE)

validator_args <- commandArgs(trailingOnly = TRUE)
if (length(validator_args) > 1L ||
    (length(validator_args) == 1L && !identical(validator_args, "--quick"))) {
  fail("Validator accepts no argument or exactly --quick")
}
quick_mode <- identical(validator_args, "--quick")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) fail("Unable to determine validator script location")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
project_path <- function(...) file.path(project_root, ...)

generator_path <- project_path("analysis", "R", "29_make_v6_biological_figures.R")
output_dir <- project_path(
  "manuscript", "v6_bmc_submission_assets", "figures_biological"
)
stems <- c("Figure_1_biological", "Figure_5_biological")
png_paths <- file.path(output_dir, paste0(stems, ".png"))
pdf_paths <- file.path(output_dir, paste0(stems, ".pdf"))
all_paths <- c(png_paths, pdf_paths)

if (!file.exists(generator_path)) fail(paste0("Figure generator is absent: ", generator_path))

run_generator <- function(argument) {
  rscript <- file.path(R.home("bin"), "Rscript")
  output <- suppressWarnings(system2(
    rscript,
    # Avoid recursively activating the project's renv lock while this
    # validator is already running under that project profile.
    c("--vanilla", shQuote(generator_path), argument),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), output = output)
}

snapshot_outputs <- function(paths) {
  present <- file.exists(paths)
  hashes <- rep(NA_character_, length(paths))
  hashes[present] <- unname(tools::md5sum(paths[present]))
  setNames(hashes, paths)
}

# Full mode checks the generator CLI and rendering determinism. Quick mode is a
# read-only validation of already-rendered outputs for fast review cycles.
if (!quick_mode) {
  # An invalid request must fail before any figure output is touched.
  before_invalid <- snapshot_outputs(all_paths)
  invalid_run <- run_generator("--figure=invalid")
  after_invalid <- snapshot_outputs(all_paths)
  if (invalid_run$status == 0L) fail("Invalid --figure argument unexpectedly succeeded")
  if (!identical(before_invalid, after_invalid)) {
    fail("Invalid --figure argument changed one or more figure outputs")
  }
  if (!any(grepl("--figure must be one of", invalid_run$output, fixed = TRUE))) {
    fail("Invalid --figure argument did not return the expected fail-fast diagnostic")
  }

  # Generate both figures twice. PNG bytes should be deterministic; PDF
  # creation metadata are intentionally excluded from hash comparison.
  first_run <- run_generator("--figure=all")
  if (first_run$status != 0L) {
    fail(paste0("--figure=all failed: ", paste(first_run$output, collapse = " | ")))
  }
  if (!all(file.exists(all_paths)) || any(file.info(all_paths)$size <= 0L)) {
    fail("--figure=all did not create four non-empty figure outputs")
  }
  first_png_hash <- unname(tools::md5sum(png_paths))

  second_run <- run_generator("--figure=all")
  if (second_run$status != 0L) {
    fail(paste0("Second --figure=all run failed: ", paste(second_run$output, collapse = " | ")))
  }
  second_png_hash <- unname(tools::md5sum(png_paths))
  if (!identical(first_png_hash, second_png_hash)) {
    fail("PNG outputs are not byte-for-byte deterministic across two runs")
  }

  for (stem in stems) {
    expected_line <- paste0("BIOLOGICAL_FIGURE_RENDERED stem=", stem)
    if (!any(grepl(expected_line, second_run$output, fixed = TRUE))) {
      fail(paste0("--figure=all did not report rendering ", stem))
    }
  }
} else {
  if (!all(file.exists(all_paths)) || any(file.info(all_paths)$size <= 0L)) {
    fail("Quick validation requires four existing non-empty figure outputs")
  }
  cat("BIOLOGICAL_FIGURE_VALIDATION_NOTE mode=quick rendering_skipped=true\n")
}

read_u32_be <- function(bytes) {
  if (length(bytes) != 4L) fail("Internal PNG parser received a non-4-byte integer")
  values <- as.numeric(as.integer(bytes))
  sum(values * c(256^3, 256^2, 256, 1))
}

read_png_metadata <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  signature <- readBin(con, what = "raw", n = 8L)
  expected_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (!identical(signature, expected_signature)) fail(paste0("Invalid PNG signature: ", path))

  width <- height <- dpi_x <- dpi_y <- NA_real_
  repeat {
    length_raw <- readBin(con, what = "raw", n = 4L)
    if (length(length_raw) == 0L) break
    if (length(length_raw) != 4L) fail(paste0("Truncated PNG chunk length: ", path))
    chunk_length <- read_u32_be(length_raw)
    type_raw <- readBin(con, what = "raw", n = 4L)
    if (length(type_raw) != 4L) fail(paste0("Truncated PNG chunk type: ", path))
    chunk_type <- rawToChar(type_raw)
    chunk_data <- readBin(con, what = "raw", n = chunk_length)
    if (length(chunk_data) != chunk_length) fail(paste0("Truncated PNG chunk data: ", path))
    crc <- readBin(con, what = "raw", n = 4L)
    if (length(crc) != 4L) fail(paste0("Truncated PNG CRC: ", path))

    if (identical(chunk_type, "IHDR")) {
      if (chunk_length != 13L) fail(paste0("Invalid PNG IHDR length: ", path))
      width <- read_u32_be(chunk_data[1:4])
      height <- read_u32_be(chunk_data[5:8])
    } else if (identical(chunk_type, "pHYs")) {
      if (chunk_length != 9L) fail(paste0("Invalid PNG pHYs length: ", path))
      if (as.integer(chunk_data[9]) == 1L) {
        dpi_x <- read_u32_be(chunk_data[1:4]) * 0.0254
        dpi_y <- read_u32_be(chunk_data[5:8]) * 0.0254
      }
    } else if (identical(chunk_type, "IEND")) {
      break
    }
  }
  c(width = width, height = height, dpi_x = dpi_x, dpi_y = dpi_y)
}

for (path in png_paths) {
  metadata <- read_png_metadata(path)
  if (!identical(unname(metadata[c("width", "height")]), c(4752, 2592))) {
    fail(sprintf(
      "Unexpected PNG dimensions for %s: %.0f x %.0f",
      basename(path), metadata[["width"]], metadata[["height"]]
    ))
  }
  if (any(!is.finite(metadata[c("dpi_x", "dpi_y")])) ||
      any(abs(metadata[c("dpi_x", "dpi_y")] - 360) > 0.1)) {
    fail(sprintf(
      "Unexpected or absent PNG DPI metadata for %s: %.3f x %.3f",
      basename(path), metadata[["dpi_x"]], metadata[["dpi_y"]]
    ))
  }
}

read_binary_text <- function(path) {
  size <- file.info(path)$size
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  bytes <- readBin(con, what = "raw", n = size)
  # PDF structure tokens are ASCII. Replacing high bytes prevents locale/UTF-8
  # warnings when scanning a binary stream that also contains compressed data.
  bytes[as.integer(bytes) > 127L] <- as.raw(32L)
  paste0(rawToChar(bytes, multiple = TRUE), collapse = "")
}

pdf_page_count <- function(path) {
  pdfinfo <- Sys.which("pdfinfo")
  if (nzchar(pdfinfo)) {
    info <- suppressWarnings(system2(pdfinfo, shQuote(path), stdout = TRUE, stderr = TRUE))
    status <- attr(info, "status")
    if (is.null(status) || status == 0L) {
      page_line <- grep("^Pages:[[:space:]]+", info, value = TRUE)
      if (length(page_line) == 1L) {
        return(as.integer(sub("^Pages:[[:space:]]+", "", page_line)))
      }
    }
  }
  content <- read_binary_text(path)
  matches <- gregexpr("/Type[[:space:]]*/Page([^sA-Za-z]|$)", content, perl = TRUE)[[1L]]
  count <- if (identical(matches, -1L)) 0L else length(matches)
  cat(sprintf(
    "BIOLOGICAL_FIGURE_VALIDATION_NOTE pdf=%s page_check=fallback_object_scan\n",
    basename(path)
  ))
  count
}

pdf_has_raster_image <- function(path) {
  pdfimages <- Sys.which("pdfimages")
  if (nzchar(pdfimages)) {
    listing <- suppressWarnings(system2(
      pdfimages, c("-list", shQuote(path)), stdout = TRUE, stderr = TRUE
    ))
    status <- attr(listing, "status")
    if (!is.null(status) && status != 0L) {
      fail(paste0("pdfimages failed for ", basename(path)))
    }
    data_lines <- listing[grepl("^[[:space:]]*[0-9]+[[:space:]]+[0-9]+", listing)]
    return(length(data_lines) > 0L)
  }
  content <- read_binary_text(path)
  cat(sprintf(
    paste0(
      "BIOLOGICAL_FIGURE_VALIDATION_NOTE pdf=%s ",
      "vector_check=fallback_no_PDF_image_subtype_object\n"
    ),
    basename(path)
  ))
  grepl("/Subtype[[:space:]]*/Image([^A-Za-z]|$)", content, perl = TRUE)
}

for (path in pdf_paths) {
  pages <- pdf_page_count(path)
  if (!identical(pages, 1L)) {
    fail(sprintf("Expected a single-page PDF for %s; observed %s", basename(path), pages))
  }
  if (pdf_has_raster_image(path)) {
    fail(paste0("Raster image object detected in vector PDF: ", basename(path)))
  }
}

normalize_text <- function(x) {
  x <- gsub("ﬁ", "fi", x, fixed = TRUE)
  x <- gsub("ﬂ", "fl", x, fixed = TRUE)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

extract_pdf_text <- function(path) {
  pdftotext <- Sys.which("pdftotext")
  if (nzchar(pdftotext)) {
    destination <- tempfile(fileext = ".txt")
    on.exit(unlink(destination), add = TRUE)
    output <- suppressWarnings(system2(
      pdftotext, c(shQuote(path), shQuote(destination)), stdout = TRUE, stderr = TRUE
    ))
    status <- attr(output, "status")
    if (is.null(status) || status == 0L) {
      return(paste(readLines(destination, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
    }
  }

  python <- Sys.which("python3")
  if (nzchar(python)) {
    code <- paste0(
      "import sys; from PyPDF2 import PdfReader; ",
      "print('\\n'.join((p.extract_text() or '') for p in PdfReader(sys.argv[1]).pages))"
    )
    output <- suppressWarnings(system2(
      python, c("-c", shQuote(code), shQuote(path)), stdout = TRUE, stderr = TRUE
    ))
    status <- attr(output, "status")
    if (is.null(status) || status == 0L) return(paste(output, collapse = "\n"))
  }
  NULL
}

generator_lines <- readLines(generator_path, warn = FALSE, encoding = "UTF-8")
generator_text <- paste(generator_lines, collapse = "\n")

# Reuse the frozen input contract as the authoritative source of displayed
# counts; do not duplicate study counts inside this output validator.
contract <- new.env(parent = globalenv())
sys.source(
  project_path("analysis", "R", "28_validate_v6_biological_figure_inputs.R"),
  envir = contract
)
display_counts <- list(
  herbs = sprintf("%d herbs", contract$provenance_counts[["herbs"]]),
  putative = format(contract$provenance_counts[["putative"]], big.mark = ","),
  disease = format(contract$provenance_counts[["disease"]], big.mark = ","),
  search_space = format(contract$provenance_counts[["search_space"]], big.mark = ","),
  frozen = format(nrow(contract$frozen), big.mark = ","),
  airway = format(contract$airway, big.mark = ","),
  bridge = format(contract$bridge, big.mark = ",")
)

required_source_fragments <- c(
  'provenance_counts[["herbs"]]',
  'provenance_counts[["putative"]]',
  'provenance_counts[["disease"]]',
  'provenance_counts[["search_space"]]',
  "figure_counts$frozen",
  "figure_counts$airway",
  "figure_counts$bridge",
  "Multi-cohort integration across pediatric blood and airway",
  "Candidate-space provenance",
  "Blood · whole blood",
  "Upper airway",
  "No cohort included BFFF exposure — intervention effects remain untested",
  "A ciliary–inflammatory core with phase-dependent repair",
  "Observed upper-airway transcriptomic organization",
  "Cilium organization / movement ↓",
  "Inflammatory programs ↑",
  "Repair / turnover",
  "Systemic blood compartment",
  "BFFF intervention layer",
  "UNTESTED",
  "no causal or temporal sequence is implied"
)
missing_source <- required_source_fragments[!vapply(
  required_source_fragments,
  function(fragment) grepl(fragment, generator_text, fixed = TRUE),
  logical(1)
)]
if (length(missing_source)) {
  fail(paste0(
    "Generator is missing required label/count bindings: ",
    paste(missing_source, collapse = "; ")
  ))
}

expected_pdf_content <- list(
  Figure_1_biological = c(
    "Multi-cohort integration across pediatric blood and airway",
    "Candidate-space provenance", display_counts$herbs, display_counts$putative,
    display_counts$disease, display_counts$search_space, display_counts$frozen,
    "GSE38900", "GSE77087", "GSE103842", "GSE155925", "GSE97742", "GSE41374",
    "Upper airway", display_counts$airway, display_counts$bridge,
    "No cohort included BFFF exposure"
  ),
  Figure_5_biological = c(
    "A ciliary", "upper-airway transcriptomic organization",
    "Cilium organization / movement", "Inflammatory programs", "Repair / turnover",
    "Systemic blood compartment", "Interferon", "Myeloid", "Neutrophil",
    display_counts$bridge, "BFFF intervention layer", "UNTESTED",
    "no causal or temporal sequence is implied"
  )
)

for (index in seq_along(pdf_paths)) {
  path <- pdf_paths[[index]]
  stem <- stems[[index]]
  extracted <- extract_pdf_text(path)
  if (is.null(extracted)) {
    cat(sprintf(
      paste0(
        "BIOLOGICAL_FIGURE_VALIDATION_NOTE pdf=%s text_check=generator_source_only ",
        "reason=no_pdftotext_or_PyPDF2\n"
      ),
      basename(path)
    ))
    next
  }
  normalized <- normalize_text(extracted)
  required <- expected_pdf_content[[stem]]
  missing <- required[!vapply(
    required,
    function(label) grepl(normalize_text(label), normalized, fixed = TRUE),
    logical(1)
  )]
  if (length(missing)) {
    fail(paste0(
      "Extracted PDF text is missing required content in ", basename(path), ": ",
      paste(missing, collapse = "; ")
    ))
  }
}

sizes <- file.info(all_paths)$size
determinism_status <- if (quick_mode) "not_rechecked" else "true"
cat(sprintf(
  paste0(
    "BIOLOGICAL_FIGURE_OUTPUT_PASS figures=2 files=4 ",
    "png=4752x2592 dpi=360 pdf_pages=1 vector=true png_deterministic=%s ",
    "bytes=%s\n"
  ),
  determinism_status,
  paste(sizes, collapse = ",")
))
