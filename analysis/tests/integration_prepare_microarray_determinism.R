#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

run_task5 <- function(label) {
  log_path <- file.path(tempdir(), paste0("task5_determinism_", label, ".log"))
  expression <- 'source("renv/activate.R"); source("analysis/R/04_prepare_microarray.R")'
  status <- system2(
    file.path(R.home("bin"), "Rscript"), c("-e", shQuote(expression)),
    stdout = log_path, stderr = log_path,
    env = c(
      "R_PROFILE_USER=/dev/null", "R_ENVIRON_USER=/dev/null",
      "RENV_CONFIG_SANDBOX_ENABLED=FALSE", "RENV_CONFIG_NAMESPACES_CHECK=FALSE",
      "TASK5_FUNCTIONS_ONLY=false"
    )
  )
  if (status != 0L) stop("Task 5 determinism run failed: ", paste(readLines(log_path, warn = FALSE), collapse = "\n"), call. = FALSE)
  checksum_path <- file.path("logs", "checksums", "task5_sha256.tsv")
  manifest <- read.delim(checksum_path, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
  manifest_sha <- sub(
    "[[:space:]].*$", "",
    system2("shasum", c("-a", "256", checksum_path), stdout = TRUE)[[1L]]
  )
  list(manifest = manifest, manifest_sha = manifest_sha)
}

if (file.exists("Rplots.pdf")) stop("Unexpected project-root Rplots.pdf before determinism test", call. = FALSE)
first <- run_task5("first")
if (file.exists("Rplots.pdf")) stop("First Task 5 run created Rplots.pdf", call. = FALSE)
second <- run_task5("second")
if (file.exists("Rplots.pdf")) stop("Second Task 5 run created Rplots.pdf", call. = FALSE)

stopifnot(
  identical(first$manifest, second$manifest),
  identical(first$manifest_sha, second$manifest_sha),
  all(c(
    file.path("results", "figures", "qc", "GSE105450_qc.pdf"),
    file.path("results", "figures", "qc", "GSE103842_qc.pdf")
  ) %in% first$manifest$artifact)
)

cat("task5 full-pipeline two-run SHA determinism passed\n")
