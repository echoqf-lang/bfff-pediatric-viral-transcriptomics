#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

read_tsv_upgrade <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv_upgrade <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest package is required", call. = FALSE)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

main_freeze_single_gene_upgrade <- function() {
  config_path <- file.path("analysis", "config", "single_gene_upgrade.yml")
  discovery_path <- file.path("validation", "validated_targets_detail.csv")
  alias_path <- file.path("data", "clean", "bfff_target_gene_aliases.tsv")
  output_path <- file.path("data", "clean", "bfff_86_gene_set_frozen.tsv")
  checksum_path <- file.path("logs", "checksums", "single_gene_upgrade_inputs_sha256.tsv")
  deviation_path <- file.path("logs", "session_info", "single_gene_upgrade_deviations.md")

  if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml package is required", call. = FALSE)
  cfg <- yaml::read_yaml(config_path)
  discovery <- utils::read.csv(discovery_path, row.names = 1L, check.names = FALSE)
  aliases <- read_tsv_upgrade(alias_path)
  aliases$entrez_id <- as.character(aliases$entrez_id)

  symbols <- rownames(discovery)
  if (length(symbols) != cfg$families$single_gene_n || anyDuplicated(symbols)) {
    stop("Discovery set must contain exactly 86 unique symbols", call. = FALSE)
  }
  matched <- aliases[match(symbols, aliases$official_symbol), c("entrez_id", "official_symbol", "ensembl_gene_id")]
  if (anyNA(matched$entrez_id) || anyDuplicated(matched$entrez_id)) {
    stop("Frozen discovery genes require unique complete Entrez mappings", call. = FALSE)
  }
  names(matched)[names(matched) == "official_symbol"] <- "gene_symbol"
  matched$discovery_source <- "GSE38900_Dallas_Turku_post_outcome_discovery"
  matched <- matched[order(matched$gene_symbol), ]
  rownames(matched) <- NULL
  write_tsv_upgrade(matched, output_path)

  inputs <- c(config_path, discovery_path, alias_path)
  checksums <- data.frame(path = inputs, sha256 = vapply(inputs, sha256_file, character(1)))
  write_tsv_upgrade(checksums, checksum_path)

  dir.create(dirname(deviation_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(deviation_path)) {
    writeLines(c("# Single-gene upgrade deviations", "", "- Initialization: no deviations."), deviation_path)
  }
  message("FROZEN genes=", nrow(matched), " unique_entrez=", length(unique(matched$entrez_id)), " missing=", sum(is.na(matched$entrez_id)))
  invisible(matched)
}

if (sys.nframe() == 0L) main_freeze_single_gene_upgrade()

