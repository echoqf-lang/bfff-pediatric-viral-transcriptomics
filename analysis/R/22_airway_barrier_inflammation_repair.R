#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
required <- c("Biobase", "limma", "AnnotationDbi", "illuminaHumanv4.db", "org.Hs.eg.db", "GO.db")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)

input_root <- file.path("data", "derived", "barrier_inflammation_repair")
output_root <- file.path("results", "barrier_inflammation_repair_v1")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

go_spec <- data.frame(
  axis = c(rep("Barrier_cilia", 6), rep("Inflammation", 5), rep("Repair_turnover", 5)),
  term = c(
    "cilium organization", "cilium movement", "epithelial cell differentiation",
    "cell-cell junction organization", "tight junction organization", "extracellular matrix organization",
    "inflammatory response", "response to type I interferon", "neutrophil activation",
    "myeloid leukocyte activation", "cytokine-mediated signaling pathway",
    "wound healing", "tissue regeneration", "epithelial cell proliferation",
    "regulation of apoptotic process", "cell cycle"
  ),
  expected = c(rep("Down", 5), "Up", rep("Up", 10)),
  stringsAsFactors = FALSE
)

go_terms <- AnnotationDbi::Term(GO.db::GOTERM)
go_spec$go_id <- vapply(go_spec$term, function(term) {
  ids <- names(go_terms)[tolower(go_terms) == tolower(term)]
  if (length(ids) != 1L) stop("GO term did not resolve uniquely: ", term, call. = FALSE)
  ids
}, character(1))

go_map <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = go_spec$go_id,
  keytype = "GOALL",
  columns = c("GOALL", "SYMBOL", "ONTOLOGYALL")
)
go_map <- go_map[go_map$ONTOLOGYALL == "BP" & !is.na(go_map$SYMBOL), c("GOALL", "SYMBOL")]
go_sets <- split(go_map$SYMBOL, factor(go_map$GOALL, levels = go_spec$go_id))
names(go_sets) <- paste(go_spec$axis, go_spec$term, sep = " | ")
write.table(go_spec, file.path(output_root, "frozen_go_terms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  do.call(rbind, lapply(seq_along(go_sets), function(i) data.frame(set = names(go_sets)[i], gene_symbol = unique(go_sets[[i]])))),
  file.path(output_root, "frozen_go_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

collapse_to_symbols <- function(eset) {
  expression <- Biobase::exprs(eset)
  if (as.numeric(stats::quantile(expression, 0.99, na.rm = TRUE)) > 100) expression <- log2(pmax(expression, 0) + 1)
  symbols <- AnnotationDbi::mapIds(
    illuminaHumanv4.db::illuminaHumanv4.db,
    keys = rownames(expression), keytype = "PROBEID", column = "SYMBOL", multiVals = "first"
  )
  keep <- !is.na(symbols) & nzchar(symbols)
  expression <- expression[keep, , drop = FALSE]
  symbols <- unname(symbols[keep])
  variances <- apply(expression, 1L, stats::var, na.rm = TRUE)
  ordering <- order(symbols, -variances, rownames(expression))
  expression <- expression[ordering, , drop = FALSE]
  symbols <- symbols[ordering]
  expression <- expression[!duplicated(symbols), , drop = FALSE]
  rownames(expression) <- symbols[!duplicated(symbols)]
  expression
}

run_analysis <- function(expression, design, coefficient, comparison, sample_table) {
  fit <- limma::eBayes(limma::lmFit(expression, design))
  gene <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
  gene$gene_symbol <- rownames(gene)
  gene$comparison <- comparison
  gene <- gene[, c("comparison", "gene_symbol", setdiff(names(gene), c("comparison", "gene_symbol")))]

  index <- lapply(go_sets, function(genes) which(rownames(expression) %in% genes))
  camera <- limma::camera(expression, index = index, design = design, contrast = coefficient)
  camera$set <- rownames(camera)
  camera$comparison <- comparison
  camera$axis <- sub(" \\|.*$", "", camera$set)
  camera$term <- sub("^[^|]+\\| ", "", camera$set)
  camera$expected <- go_spec$expected[match(camera$term, go_spec$term)]
  camera$direction_matches_expected <- camera$Direction == camera$expected
  camera$fdr_bh_fixed_family <- p.adjust(camera$PValue, method = "BH")
  camera <- camera[, c("comparison", "axis", "term", "NGenes", "Direction", "expected",
                       "direction_matches_expected", "PValue", "fdr_bh_fixed_family")]

  write.table(gene, file.path(output_root, paste0(comparison, "_gene_effects.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(camera, file.path(output_root, paste0(comparison, "_camera.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(sample_table, file.path(output_root, paste0(comparison, "_samples.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  camera
}

gse97742 <- readRDS(file.path(input_root, "GSE97742_processed_eset.rds"))
expression_97742 <- collapse_to_symbols(gse97742)
pd_97742 <- Biobase::pData(gse97742)
pd_97742$phase <- sub("phase: ", "", pd_97742[["phase:ch1"]], fixed = TRUE)
pd_97742$status <- sub("subject status: ", "", pd_97742[["subject status:ch1"]], fixed = TRUE)
pd_97742$infection <- ifelse(grepl("RSVsi", pd_97742$status, fixed = TRUE), "RSVsi",
                             ifelse(grepl("RSVco", pd_97742$status, fixed = TRUE), "RSVco", "hRV"))
pd_97742$subject <- sub(".* (RSVsi|RSVco|hRV) ", "", pd_97742$title)
pd_97742$subject_key <- paste(pd_97742$infection, pd_97742$subject, sep = "_")

camera_results <- list()
for (infection in c("RSVsi", "hRV")) {
  keep <- pd_97742$infection == infection
  meta <- pd_97742[keep, c("geo_accession", "title", "phase", "infection", "subject_key"), drop = FALSE]
  counts <- table(meta$subject_key, meta$phase)
  paired <- rownames(counts)[counts[, "acute"] == 1L & counts[, "discharge"] == 1L]
  meta <- meta[meta$subject_key %in% paired, , drop = FALSE]
  meta <- meta[order(meta$subject_key, factor(meta$phase, levels = c("discharge", "acute"))), ]
  x <- expression_97742[, rownames(meta), drop = FALSE]
  meta$subject_key <- factor(meta$subject_key)
  meta$phase <- factor(meta$phase, levels = c("discharge", "acute"))
  design <- stats::model.matrix(~ subject_key + phase, data = meta)
  coefficient <- which(colnames(design) == "phaseacute")
  comparison <- paste0("GSE97742_", infection, "_acute_vs_discharge")
  camera_results[[comparison]] <- run_analysis(x, design, coefficient, comparison, meta)
}

gse41374 <- readRDS(file.path(input_root, "GSE41374_processed_eset.rds"))
expression_41374 <- collapse_to_symbols(gse41374)
pd_41374 <- Biobase::pData(gse41374)
meta_41374 <- data.frame(
  geo_accession = pd_41374$geo_accession,
  title = pd_41374$title,
  age_months = as.numeric(pd_41374[["age (month):ch1"]]),
  sex = factor(pd_41374[["gender:ch1"]]),
  group = factor(pd_41374[["infection:ch1"]], levels = c("healthy", "RSV")),
  row.names = rownames(pd_41374), stringsAsFactors = FALSE
)
complete <- stats::complete.cases(meta_41374[, c("age_months", "sex", "group")])
meta_41374 <- meta_41374[complete, , drop = FALSE]
x_41374 <- expression_41374[, rownames(meta_41374), drop = FALSE]
design_41374 <- stats::model.matrix(~ age_months + sex + group, data = meta_41374)
coefficient_41374 <- which(colnames(design_41374) == "groupRSV")
comparison_41374 <- "GSE41374_RSV_vs_healthy"
camera_results[[comparison_41374]] <- run_analysis(
  x_41374, design_41374, coefficient_41374, comparison_41374, meta_41374
)

summary <- do.call(rbind, camera_results)
rownames(summary) <- NULL
write.table(summary, file.path(output_root, "three_axis_camera_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

decision <- aggregate(
  cbind(significant = summary$fdr_bh_fixed_family < 0.05,
        expected_and_significant = summary$fdr_bh_fixed_family < 0.05 & summary$direction_matches_expected) ~ comparison + axis,
  data = summary, FUN = sum
)
write.table(decision, file.path(output_root, "axis_support_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Completed three comparisons; output rows=", nrow(summary))
