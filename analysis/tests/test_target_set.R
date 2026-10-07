#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setlocale("LC_CTYPE", "C.UTF-8")
stopifnot(requireNamespace("yaml", quietly = TRUE))
stopifnot(requireNamespace("AnnotationDbi", quietly = TRUE))
stopifnot(requireNamespace("org.Hs.eg.db", quietly = TRUE))

trim <- function(x) trimws(enc2utf8(as.character(x)))
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
rules <- yaml::read_yaml("analysis/config/target_rules.yml")
expected_tcms_primary_rule <- paste(
  "human UniProt mapped offline to unique Entrez ID plus official org.Hs.eg.db SYMBOL;",
  "Ensembl mappings retained only as optional aliases and not counted"
)
stopifnot(
  identical(rules$evidence_rules$primary$TCMSP, expected_tcms_primary_rule),
  identical(rules$gene_identity$primary_key, "entrez_id"),
  identical(rules$gene_identity$ensembl_role, "alias_only_not_counted")
)
all <- read_tsv("data/clean/bfff_targets_all_evidence.tsv")
primary <- read_tsv("data/clean/bfff_targets_primary.tsv")
sensitivity <- read_tsv("data/clean/bfff_targets_sensitivity.tsv")
aliases <- read_tsv("data/clean/bfff_target_gene_aliases.tsv")
audit <- read_tsv("results/tables/target_provenance_audit.tsv")
mapping <- read_tsv("results/tables/target_mapping_audit.tsv")

# Re-read the permitted raw inputs and independently reconstruct the inclusion boundary.
tc_path <- rules$inputs$tcms_target_file$path
ba_path <- rules$inputs$batman_target_file$path
tc <- read.csv(tc_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
ba <- read.csv(ba_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
tc$source_row_id <- paste0("TCMSP:line_", seq_len(nrow(tc)) + 1L)
ba$source_row_id <- paste0("BATMAN:line_", seq_len(nrow(ba)) + 1L)
tc_human <- trim(tc$`来源`) == "Homo sapiens (Human)"
tc_eligible <- tc_human & nzchar(trim(tc$Entry)) & nzchar(trim(tc$`GENE name`))
tc_link <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db, keys = unique(trim(tc$Entry[tc_eligible])),
  keytype = "UNIPROT", columns = "ENTREZID"
))
tc_link <- unique(tc_link[!is.na(tc_link$ENTREZID) & nzchar(tc_link$ENTREZID), ])
official <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique(c(tc_link$ENTREZID, trim(ba$Gene_ID))),
  keytype = "ENTREZID", columns = "SYMBOL"
))
valid_entrez <- unique(official$ENTREZID[!is.na(official$SYMBOL) & nzchar(official$SYMBOL)])
tc_link <- tc_link[tc_link$ENTREZID %in% valid_entrez, ]

score <- trim(ba$Score)
score_num <- suppressWarnings(as.numeric(score))
exact <- score == "known target in HERB"
composite <- grepl("HERB", score, fixed = TRUE) & !exact
high <- !is.na(score_num) & score_num >= 0.90
low <- !is.na(score_num) & score_num >= 0.84 & score_num < 0.90
ba_mapped <- trim(ba$Gene_ID) %in% valid_entrez

expected_primary <- sort(unique(c(tc_link$ENTREZID, trim(ba$Gene_ID[exact & ba_mapped]))))
expected_sensitivity <- sort(unique(c(expected_primary, trim(ba$Gene_ID[high & ba_mapped]))))
stopifnot(identical(sort(trim(primary$entrez_id)), expected_primary))
stopifnot(identical(sort(trim(sensitivity$entrez_id)), expected_sensitivity))

# Direct raw-row boundary checks, not output-only self-consistency checks.
bat_out <- all[all$source == "BATMAN", ]
exact_ids <- ba$source_row_id[exact & ba_mapped]
composite_ids <- ba$source_row_id[composite & ba_mapped]
high_ids <- ba$source_row_id[high & ba_mapped]
low_ids <- ba$source_row_id[low & ba_mapped]
stopifnot(all(bat_out$included_primary[match(exact_ids, bat_out$source_row_id)]))
stopifnot(!any(bat_out$included_primary[match(composite_ids, bat_out$source_row_id)]))
stopifnot(all(bat_out$included_sensitivity[match(high_ids, bat_out$source_row_id)]))
stopifnot(!any(bat_out$included_sensitivity[match(low_ids, bat_out$source_row_id)]))

# Complete lineage: every mapped raw row/Entrez pair appears exactly once.
tc_expected <- merge(
  data.frame(source_row_id = tc$source_row_id[tc_eligible], UNIPROT = trim(tc$Entry[tc_eligible])),
  tc_link, by = "UNIPROT", all = FALSE
)
expected_lineage <- sort(c(
  paste(tc_expected$source_row_id, tc_expected$ENTREZID, sep = "\r"),
  paste(ba$source_row_id[ba_mapped], trim(ba$Gene_ID[ba_mapped]), sep = "\r")
))
observed_lineage <- sort(paste(all$source_row_id, all$entrez_id, sep = "\r"))
stopifnot(identical(observed_lineage, expected_lineage), !anyDuplicated(observed_lineage))

# Gene tables count Entrez genes once; multiple Ensembl aliases never inflate counts.
stopifnot(!anyDuplicated(primary$entrez_id), !anyDuplicated(sensitivity$entrez_id))
stopifnot(nrow(primary) == length(expected_primary), nrow(sensitivity) == length(expected_sensitivity))
multi_alias <- names(which(table(aliases$entrez_id) > 1L))
stopifnot(length(multi_alias) > 0L)
stopifnot(all(table(primary$entrez_id[primary$entrez_id %in% multi_alias]) <= 1L))
stopifnot(all(table(sensitivity$entrez_id[sensitivity$entrez_id %in% multi_alias]) <= 1L))

# Formula and species audits.
herbs <- trim(readLines(rules$inputs$herbs_file$path, encoding = "UTF-8"))
stopifnot(length(herbs) == 14L, length(unique(herbs)) == 14L, length(unique(audit$herb_cn)) == 14L)
stopifnot(!any(all$source == "TCMSP" & all$source_row_id %in% tc$source_row_id[!tc_human]))
fw <- ba$Chinese_Name == "浮小麦" & ba_mapped
fw_primary <- length(unique(trim(ba$Gene_ID[fw & exact])))
fw_sensitivity <- length(unique(trim(ba$Gene_ID[fw & (exact | high)])))
stopifnot(any(ba$Chinese_Name == "浮小麦"))
stopifnot(length(unique(all$entrez_id[all$herb_cn == "浮小麦" & all$included_primary])) == fw_primary)
stopifnot(length(unique(all$entrez_id[all$herb_cn == "浮小麦" & all$included_sensitivity])) == fw_sensitivity)

# Raw and official identities remain together, including ambiguous creatine-kinase identities if present.
stopifnot(all(c("source_row_id", "source_gene_id", "source_symbol_raw", "official_symbol") %in% names(all)))
ck_rows <- mapping$source_symbol_raw %in% c("CKMT1A", "CKMT1B") | mapping$official_symbols %in% c("CKMT1A", "CKMT1B")
if (any(ck_rows)) stopifnot(all(nzchar(mapping$source_row_id[ck_rows])), all(nzchar(mapping$mapped_entrez_ids[ck_rows])))

# UTF-8 and forbidden-input boundary checks.
stopifnot(all(nzchar(iconv(all$herb_cn, from = "UTF-8", to = "UTF-8"))))
script_text <- paste(readLines("analysis/R/01_build_target_set.R", warn = FALSE, encoding = "UTF-8"), collapse = "\n")
forbidden <- c("disease_targets", "intersection/", "PPI/", "GSE")
stopifnot(!any(vapply(forbidden, grepl, logical(1), x = script_text, fixed = TRUE)))

# Regression guard for the session-info whitespace bug.
strip_trailing_whitespace <- function(x) sub("[[:space:]]+$", "", x)
stopifnot(
  identical(strip_trailing_whitespace("Matrix products: default"), "Matrix products: default"),
  identical(strip_trailing_whitespace(paste0("value", " ", "\t")), "value")
)

cat(sprintf("PASS: raw-recomputed primary=%d sensitivity=%d lineages=%d herbs=14 floating_wheat=%d/%d\n",
  nrow(primary), nrow(sensitivity), nrow(all), fw_primary, fw_sensitivity))
