#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setlocale("LC_CTYPE", "C.UTF-8")
set.seed(20260804L)

required <- c("AnnotationDbi", "org.Hs.eg.db", "yaml")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)

rules <- yaml::read_yaml("analysis/config/target_rules.yml")
trim <- function(x) trimws(enc2utf8(as.character(x)))
empty_to_na <- function(x) { x <- trim(x); x[x == ""] <- NA_character_; x }
collapse_unique <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)]), method = "radix"), collapse = ";")
standardize_name <- function(x) {
  x <- tolower(trim(x))
  x <- gsub("[^[:alnum:]]+", " ", x)
  trimws(gsub("[[:space:]]+", " ", x))
}
sha256 <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  sub("[[:space:]].*$", "", out[[1]])
}
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

input_specs <- rules$inputs
input_paths <- vapply(input_specs, `[[`, character(1), "path")
expected_hashes <- vapply(input_specs, `[[`, character(1), "sha256")
observed_hashes <- vapply(input_paths, sha256, character(1))
if (!identical(unname(observed_hashes), unname(expected_hashes))) {
  stop("Input SHA-256 mismatch", call. = FALSE)
}

herbs <- trim(readLines(input_paths[["herbs_file"]], encoding = "UTF-8"))
configured_herbs <- vapply(rules$herbs, `[[`, character(1), "cn")
if (!identical(herbs, configured_herbs) || length(unique(herbs)) != 14L) {
  stop("herbs.txt and configured 14-herb order differ", call. = FALSE)
}

tc <- read.csv(input_paths[["tcms_target_file"]], check.names = FALSE, fileEncoding = "UTF-8")
ba <- read.csv(input_paths[["batman_target_file"]], check.names = FALSE, fileEncoding = "UTF-8")
tc$source_row_number <- seq_len(nrow(tc)) + 1L
ba$source_row_number <- seq_len(nrow(ba)) + 1L
tc$source_row_id <- paste0("TCMSP:line_", tc$source_row_number)
ba$source_row_id <- paste0("BATMAN:line_", ba$source_row_number)

tc_human <- trim(tc$`来源`) == "Homo sapiens (Human)"
tc_required <- tc_human & !is.na(empty_to_na(tc$Entry)) & !is.na(empty_to_na(tc$`GENE name`))
tc_keys <- unique(trim(tc$Entry[tc_required]))
tc_link <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db, keys = tc_keys, keytype = "UNIPROT", columns = "ENTREZID"
))
tc_link <- unique(tc_link[!is.na(tc_link$ENTREZID) & nzchar(tc_link$ENTREZID), ])
names(tc_link)[names(tc_link) == "UNIPROT"] <- "uniprot_key"

ba_keys <- unique(trim(ba$Gene_ID))
all_entrez <- unique(c(tc_link$ENTREZID, ba_keys))
gene_raw <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db, keys = all_entrez, keytype = "ENTREZID", columns = c("SYMBOL", "ENSEMBL")
))
gene_raw <- gene_raw[!is.na(gene_raw$SYMBOL) & nzchar(gene_raw$SYMBOL), ]
gene_raw$ENTREZID <- trim(gene_raw$ENTREZID)
gene_raw$SYMBOL <- trim(gene_raw$SYMBOL)
official_counts <- tapply(gene_raw$SYMBOL, gene_raw$ENTREZID, function(z) length(unique(z)))
if (any(official_counts > 1L)) stop("OrgDb returned multiple official symbols for one Entrez ID", call. = FALSE)

gene_ids <- sort(unique(gene_raw$ENTREZID))
gene_info <- data.frame(
  entrez_id = gene_ids,
  official_symbol = vapply(gene_ids, function(k) collapse_unique(gene_raw$SYMBOL[gene_raw$ENTREZID == k]), character(1)),
  ensembl_gene_ids = vapply(gene_ids, function(k) collapse_unique(gene_raw$ENSEMBL[gene_raw$ENTREZID == k]), character(1)),
  n_ensembl_aliases = vapply(gene_ids, function(k) length(unique(na.omit(gene_raw$ENSEMBL[gene_raw$ENTREZID == k]))), integer(1)),
  stringsAsFactors = FALSE
)
gene_info <- gene_info[nzchar(gene_info$official_symbol), ]
gene_aliases <- unique(gene_raw[!is.na(gene_raw$ENSEMBL) & nzchar(gene_raw$ENSEMBL), c("ENTREZID", "SYMBOL", "ENSEMBL")])
names(gene_aliases) <- c("entrez_id", "official_symbol", "ensembl_gene_id")
gene_aliases <- gene_aliases[order(as.integer(gene_aliases$entrez_id), gene_aliases$ensembl_gene_id), ]

tc_valid <- tc[tc_required, ]
tc_valid$uniprot_key <- trim(tc_valid$Entry)
tc_ev <- merge(tc_valid, tc_link, by = "uniprot_key", all = FALSE, sort = FALSE)
tc_ev <- merge(tc_ev, gene_info, by.x = "ENTREZID", by.y = "entrez_id", all = FALSE, sort = FALSE)

score <- trim(ba$Score)
score_num <- suppressWarnings(as.numeric(score))
ba$evidence_class <- ifelse(
  score == "known target in HERB", "batman_known_herb_exact",
  ifelse(!is.na(score_num) & score_num >= 0.90, "batman_predicted_ge_0_90",
    ifelse(!is.na(score_num) & score_num >= 0.84 & score_num < 0.90,
      "batman_predicted_0_84_0_89",
      ifelse(grepl("^known target in ", score), "batman_known_other", "batman_other")
    )
  )
)
ba$entrez_id <- trim(ba$Gene_ID)
ba_ev <- merge(ba, gene_info, by = "entrez_id", all = FALSE, sort = FALSE)

make_evidence <- function(source) {
  if (source == "TCMSP") {
    data.frame(
      source = source, source_row_id = tc_ev$source_row_id, source_row_number = tc_ev$source_row_number,
      herb_cn = trim(tc_ev$`中药名`), compound_id_type = "TCMSP_MOL_ID",
      compound_id = trim(tc_ev$MOL_ID), compound_name_raw = trim(tc_ev$molecule_name),
      compound_name_standardized = standardize_name(tc_ev$molecule_name),
      source_gene_id_type = "UniProt", source_gene_id = tc_ev$uniprot_key,
      source_symbol_raw = trim(tc_ev$`GENE name`), official_symbol = tc_ev$official_symbol,
      gene_symbol = tc_ev$official_symbol, entrez_id = tc_ev$ENTREZID,
      ensembl_gene_ids = tc_ev$ensembl_gene_ids, uniprot_id = tc_ev$uniprot_key,
      evidence_class = "tcms_human_uniprot_hgnc", score_raw = "",
      included_primary = TRUE, included_sensitivity = TRUE,
      input_file_sha256 = observed_hashes[["tcms_target_file"]],
      mapping_method = "UNIPROT_to_ENTREZID_org.Hs.eg.db",
      symbol_conflict = trim(tc_ev$`GENE name`) != tc_ev$official_symbol,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      source = source, source_row_id = ba_ev$source_row_id, source_row_number = ba_ev$source_row_number,
      herb_cn = trim(ba_ev$Chinese_Name), compound_id_type = "PubChem_CID",
      compound_id = trim(ba_ev$CID), compound_name_raw = trim(ba_ev$Compound_Name),
      compound_name_standardized = standardize_name(ba_ev$Compound_Name),
      source_gene_id_type = "Entrez", source_gene_id = trim(ba_ev$Gene_ID),
      source_symbol_raw = trim(ba_ev$Gene_Name), official_symbol = ba_ev$official_symbol,
      gene_symbol = ba_ev$official_symbol, entrez_id = ba_ev$entrez_id,
      ensembl_gene_ids = ba_ev$ensembl_gene_ids, uniprot_id = NA_character_,
      evidence_class = ba_ev$evidence_class, score_raw = trim(ba_ev$Score),
      included_primary = ba_ev$evidence_class == "batman_known_herb_exact",
      included_sensitivity = ba_ev$evidence_class %in% c("batman_known_herb_exact", "batman_predicted_ge_0_90"),
      input_file_sha256 = observed_hashes[["batman_target_file"]],
      mapping_method = "ENTREZID_to_org.Hs.eg.db",
      symbol_conflict = !is.na(empty_to_na(ba_ev$Gene_Name)) & trim(ba_ev$Gene_Name) != ba_ev$official_symbol,
      stringsAsFactors = FALSE
    )
  }
}
all_evidence <- rbind(make_evidence("TCMSP"), make_evidence("BATMAN"))
required_nonempty <- c("source", "source_row_id", "herb_cn", "compound_id_type", "compound_id",
  "compound_name_standardized", "official_symbol", "entrez_id", "evidence_class")
if (any(vapply(all_evidence[required_nonempty], function(x) any(is.na(x) | !nzchar(x)), logical(1)))) {
  stop("Mapped evidence has an empty required key", call. = FALSE)
}
if (any(!all_evidence$herb_cn %in% herbs)) stop("Evidence contains an unconfigured herb", call. = FALSE)
lineage_key <- paste(all_evidence$source_row_id, all_evidence$entrez_id, sep = "\r")
if (anyDuplicated(lineage_key)) stop("Duplicate source-row/Entrez lineage", call. = FALSE)
all_evidence <- all_evidence[order(match(all_evidence$herb_cn, herbs), all_evidence$source,
  all_evidence$source_row_number, as.integer(all_evidence$entrez_id)), ]

make_gene_set <- function(evidence, include_column) {
  x <- evidence[evidence[[include_column]], ]
  ids <- sort(unique(x$entrez_id), method = "radix")
  out <- lapply(ids, function(id) {
    z <- x[x$entrez_id == id, ]
    symbols <- unique(z$official_symbol)
    if (length(symbols) != 1L) stop("Non-unique official symbol for Entrez ", id, call. = FALSE)
    data.frame(
      entrez_id = id, gene_symbol = symbols,
      ensembl_gene_ids = collapse_unique(z$ensembl_gene_ids),
      n_ensembl_aliases = length(unique(unlist(strsplit(z$ensembl_gene_ids[nzchar(z$ensembl_gene_ids)], ";", fixed = TRUE)))),
      supporting_herbs = collapse_unique(z$herb_cn), n_supporting_herbs = length(unique(z$herb_cn)),
      supporting_sources = collapse_unique(z$source), n_supporting_sources = length(unique(z$source)),
      supporting_source_rows = collapse_unique(z$source_row_id), n_supporting_source_rows = length(unique(z$source_row_id)),
      supporting_compounds = collapse_unique(paste(z$compound_id_type, z$compound_id, z$compound_name_standardized, sep = ":")),
      n_supporting_compounds = length(unique(paste(z$compound_id_type, z$compound_id, z$compound_name_standardized, sep = "\r"))),
      evidence_classes = collapse_unique(z$evidence_class), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
primary <- make_gene_set(all_evidence, "included_primary")
sensitivity <- make_gene_set(all_evidence, "included_sensitivity")

tc_mapped <- unique(tc_link$uniprot_key[tc_link$ENTREZID %in% gene_info$entrez_id])
tc_entrez_by_key <- function(k) collapse_unique(tc_link$ENTREZID[tc_link$uniprot_key == k & tc_link$ENTREZID %in% gene_info$entrez_id])
tc_symbols_by_key <- function(k) collapse_unique(gene_info$official_symbol[gene_info$entrez_id %in% tc_link$ENTREZID[tc_link$uniprot_key == k]])
tc_audit <- data.frame(
  source = "TCMSP", source_row_id = tc$source_row_id, source_row_number = tc$source_row_number,
  herb_cn = trim(tc$`中药名`), compound_id_type = "TCMSP_MOL_ID", compound_id = trim(tc$MOL_ID),
  compound_name_raw = trim(tc$molecule_name), evidence_class = "tcms_human_uniprot_hgnc",
  source_gene_id_type = "UniProt", source_gene_id = trim(tc$Entry), source_symbol_raw = trim(tc$`GENE name`),
  species_raw = trim(tc$`来源`),
  mapping_status = ifelse(!tc_human, "excluded_nonhuman", ifelse(!tc_required, "excluded_missing_required_id",
    ifelse(trim(tc$Entry) %in% tc_mapped, "mapped", "unmapped"))),
  mapped_entrez_ids = vapply(trim(tc$Entry), tc_entrez_by_key, character(1)),
  official_symbols = vapply(trim(tc$Entry), tc_symbols_by_key, character(1)),
  multi_entrez_mapping = vapply(trim(tc$Entry), function(k) length(unique(tc_link$ENTREZID[tc_link$uniprot_key == k])) > 1L, logical(1)),
  symbol_conflict = vapply(seq_len(nrow(tc)), function(i) {
    sy <- tc_symbols_by_key(trim(tc$Entry[[i]])); nzchar(sy) && !trim(tc$`GENE name`[[i]]) %in% strsplit(sy, ";", fixed = TRUE)[[1]]
  }, logical(1)), stringsAsFactors = FALSE
)
ba_audit <- data.frame(
  source = "BATMAN", source_row_id = ba$source_row_id, source_row_number = ba$source_row_number,
  herb_cn = trim(ba$Chinese_Name), compound_id_type = "PubChem_CID", compound_id = trim(ba$CID),
  compound_name_raw = trim(ba$Compound_Name), evidence_class = ba$evidence_class,
  source_gene_id_type = "Entrez", source_gene_id = trim(ba$Gene_ID), source_symbol_raw = trim(ba$Gene_Name),
  species_raw = "Homo sapiens (mapped against taxon 9606 OrgDb)",
  mapping_status = ifelse(trim(ba$Gene_ID) %in% gene_info$entrez_id, "mapped", "unmapped"),
  mapped_entrez_ids = ifelse(trim(ba$Gene_ID) %in% gene_info$entrez_id, trim(ba$Gene_ID), ""),
  official_symbols = vapply(trim(ba$Gene_ID), function(k) collapse_unique(gene_info$official_symbol[gene_info$entrez_id == k]), character(1)),
  multi_entrez_mapping = FALSE,
  symbol_conflict = vapply(seq_len(nrow(ba)), function(i) {
    sy <- gene_info$official_symbol[gene_info$entrez_id == trim(ba$Gene_ID[[i]])]
    length(sy) == 1L && !is.na(empty_to_na(ba$Gene_Name[[i]])) && trim(ba$Gene_Name[[i]]) != sy
  }, logical(1)), stringsAsFactors = FALSE
)
mapping_audit <- rbind(tc_audit, ba_audit)
mapping_audit$mapping_conflict <- mapping_audit$multi_entrez_mapping | mapping_audit$symbol_conflict
unmapped <- mapping_audit[mapping_audit$mapping_status != "mapped", ]

classes <- c("tcms_human_uniprot_hgnc", "batman_known_herb_exact", "batman_predicted_ge_0_90",
  "batman_predicted_0_84_0_89", "batman_known_other", "batman_other")
grid <- expand.grid(herb_cn = herbs, source = c("TCMSP", "BATMAN"), evidence_class = classes, stringsAsFactors = FALSE)
provenance_audit <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  z <- all_evidence[all_evidence$herb_cn == g$herb_cn & all_evidence$source == g$source & all_evidence$evidence_class == g$evidence_class, ]
  raw <- mapping_audit[mapping_audit$herb_cn == g$herb_cn & mapping_audit$source == g$source & mapping_audit$evidence_class == g$evidence_class, ]
  data.frame(g, n_input_rows = nrow(raw), n_mapped_input_rows = sum(raw$mapping_status == "mapped"),
    n_unmapped_or_excluded_input_rows = sum(raw$mapping_status != "mapped"), n_evidence_lineages = nrow(z),
    n_unique_entrez_genes = length(unique(z$entrez_id)), n_unique_compounds = length(unique(paste(z$compound_id_type, z$compound_id, z$compound_name_standardized))),
    n_primary_genes = length(unique(z$entrez_id[z$included_primary])), n_sensitivity_genes = length(unique(z$entrez_id[z$included_sensitivity])),
    stringsAsFactors = FALSE)
}))

write_tsv(all_evidence, "data/clean/bfff_targets_all_evidence.tsv")
write_tsv(primary, "data/clean/bfff_targets_primary.tsv")
write_tsv(sensitivity, "data/clean/bfff_targets_sensitivity.tsv")
write_tsv(gene_aliases, "data/clean/bfff_target_gene_aliases.tsv")
write_tsv(provenance_audit, "results/tables/target_provenance_audit.tsv")
write_tsv(mapping_audit, "results/tables/target_mapping_audit.tsv")
write_tsv(unmapped, "results/tables/target_unmapped_or_excluded.tsv")

cat(sprintf("Task 2 target build complete: evidence_lineages=%d primary_entrez=%d sensitivity_entrez=%d\n",
  nrow(all_evidence), nrow(primary), nrow(sensitivity)))
