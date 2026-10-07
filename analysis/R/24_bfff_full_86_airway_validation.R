#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

load_target_contract <- function(target_file, representative_file) {
  targets <- read.delim(target_file, check.names = FALSE)
  panel <- read.delim(representative_file, check.names = FALSE)
  stopifnot(nrow(targets) == 86L, !anyDuplicated(targets$gene_symbol))
  stopifnot(nrow(panel) == 20L, !anyDuplicated(panel$gene_symbol))
  stopifnot(all(panel$gene_symbol %in% targets$gene_symbol))

  panel_index <- match(targets$gene_symbol, panel$gene_symbol)
  targets$is_representative <- !is.na(panel_index)
  targets$axis <- ifelse(
    targets$is_representative,
    panel$axis[panel_index],
    "Other_unclassified"
  )
  targets$role <- ifelse(targets$is_representative, panel$role[panel_index], NA_character_)
  targets$priority <- ifelse(targets$is_representative, panel$priority[panel_index], NA_character_)
  targets
}

require_analysis_packages <- function() {
  required <- c("Biobase", "limma", "AnnotationDbi", "illuminaHumanv4.db")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

map_and_collapse <- function(expression, target_symbols) {
  symbols <- AnnotationDbi::mapIds(
    illuminaHumanv4.db::illuminaHumanv4.db,
    keys = rownames(expression), keytype = "PROBEID", column = "SYMBOL", multiVals = "first"
  )
  keep <- !is.na(symbols) & symbols %in% target_symbols
  expression <- expression[keep, , drop = FALSE]
  symbols <- unname(symbols[keep])
  variances <- apply(expression, 1L, stats::var, na.rm = TRUE)
  ordering <- order(symbols, -variances, rownames(expression))
  expression <- expression[ordering, , drop = FALSE]
  symbols <- symbols[ordering]
  retain <- !duplicated(symbols)
  expression <- expression[retain, , drop = FALSE]
  rownames(expression) <- symbols[retain]
  expression
}

annotate_results <- function(contract, results) {
  result_index <- match(contract$gene_symbol, results$gene_symbol)
  extra_columns <- setdiff(names(results), "gene_symbol")
  cbind(contract, results[result_index, extra_columns, drop = FALSE])
}

run_gse97742_interaction <- function(eset, contract) {
  expression <- map_and_collapse(Biobase::exprs(eset), contract$gene_symbol)
  pd <- Biobase::pData(eset)
  pd$phase <- sub("phase: ", "", pd[["phase:ch1"]], fixed = TRUE)
  status <- sub("subject status: ", "", pd[["subject status:ch1"]], fixed = TRUE)
  pd$infection <- ifelse(grepl("RSVsi", status, fixed = TRUE), "RSVsi",
                         ifelse(grepl("RSVco", status, fixed = TRUE), "RSVco", "hRV"))
  pd$subject <- sub(".* (RSVsi|RSVco|hRV) ", "", pd$title)
  pd$subject_key <- paste(pd$infection, pd$subject, sep = "_")
  pd$age_months <- as.numeric(pd[["age:ch1"]])
  pd$sex <- factor(toupper(pd[["Sex:ch1"]]))
  pd <- pd[pd$infection %in% c("RSVsi", "hRV"), , drop = FALSE]

  rows <- lapply(rownames(expression), function(gene) {
    values <- expression[gene, rownames(pd)]
    wide <- reshape(
      data.frame(subject_key = pd$subject_key, infection = pd$infection,
                 age_months = pd$age_months, sex = pd$sex, phase = pd$phase, value = values),
      idvar = c("subject_key", "infection", "age_months", "sex"),
      timevar = "phase", direction = "wide"
    )
    wide <- wide[stats::complete.cases(wide), , drop = FALSE]
    wide$delta <- wide$value.acute - wide$value.discharge
    wide$infection <- factor(wide$infection, levels = c("hRV", "RSVsi"))
    fit <- stats::lm(delta ~ infection + age_months + sex, data = wide)
    coefficient <- summary(fit)$coefficients["infectionRSVsi", ]
    data.frame(
      gene_symbol = gene,
      n_subjects = nrow(wide),
      rsv_minus_hrv_delta = unname(coefficient["Estimate"]),
      se = unname(coefficient["Std. Error"]),
      p_value = unname(coefficient["Pr(>|t|)"]),
      rsv_mean_delta = mean(wide$delta[wide$infection == "RSVsi"]),
      hrv_mean_delta = mean(wide$delta[wide$infection == "hRV"])
    )
  })
  results <- do.call(rbind, rows)
  results$fdr_bh_86_family <- p.adjust(results$p_value, method = "BH")
  annotate_results(contract, results)
}

fit_gse41374 <- function(eset, contract, transform_name) {
  raw <- Biobase::exprs(eset)
  transformed <- switch(
    transform_name,
    shifted_log2 = log2(raw - min(raw, na.rm = TRUE) + 1),
    asinh_robust = asinh(raw / stats::median(abs(raw), na.rm = TRUE)),
    stop("Unknown transform", call. = FALSE)
  )
  expression <- map_and_collapse(transformed, contract$gene_symbol)
  pd <- Biobase::pData(eset)
  meta <- data.frame(
    age_months = as.numeric(pd[["age (month):ch1"]]),
    sex = factor(pd[["gender:ch1"]]),
    group = factor(pd[["infection:ch1"]], levels = c("healthy", "RSV")),
    row.names = rownames(pd)
  )
  meta <- meta[stats::complete.cases(meta), , drop = FALSE]
  expression <- expression[, rownames(meta), drop = FALSE]
  design <- stats::model.matrix(~ age_months + sex + group, data = meta)
  fit <- limma::eBayes(limma::lmFit(expression, design))
  table <- limma::topTable(fit, coef = "groupRSV", number = Inf, sort.by = "none")
  results <- data.frame(
    gene_symbol = rownames(table),
    transform = transform_name,
    n_samples = nrow(meta),
    logFC = table$logFC,
    se = table$logFC / table$t,
    t = table$t,
    p_value = table$P.Value,
    stringsAsFactors = FALSE
  )
  results$fdr_bh_86_family <- p.adjust(results$p_value, method = "BH")
  annotate_results(contract, results)
}

run_full_analysis <- function() {
  require_analysis_packages()
  contract <- load_target_contract(
    "data/clean/bfff_86_gene_set_frozen.tsv",
    "analysis/config/bfff_barrier_inflammation_repair_panel.tsv"
  )
  input_root <- file.path("data", "derived", "barrier_inflammation_repair")
  output_root <- file.path("results", "bfff_full_86_airway_v2")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  gse97742 <- run_gse97742_interaction(
    readRDS(file.path(input_root, "GSE97742_processed_eset.rds")), contract
  )
  write.table(gse97742, file.path(output_root, "GSE97742_full_86_RSV_vs_hRV_recovery_interaction.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  eset41374 <- readRDS(file.path(input_root, "GSE41374_processed_eset.rds"))
  gse41374 <- rbind(
    fit_gse41374(eset41374, contract, "shifted_log2"),
    fit_gse41374(eset41374, contract, "asinh_robust")
  )
  write.table(gse41374, file.path(output_root, "GSE41374_full_86_two_transform_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  compact <- gse41374[, c("gene_symbol", "transform", "logFC", "p_value", "fdr_bh_86_family")]
  wide <- reshape(compact, idvar = "gene_symbol", timevar = "transform", direction = "wide")
  wide$direction_concordant <- with(wide, sign(logFC.shifted_log2) == sign(logFC.asinh_robust))
  wide$robust_86_family_fdr <- with(
    wide,
    direction_concordant & fdr_bh_86_family.shifted_log2 < 0.05 &
      fdr_bh_86_family.asinh_robust < 0.05
  )
  wide <- annotate_results(contract, wide)
  write.table(wide, file.path(output_root, "GSE41374_full_86_robust_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  summary <- data.frame(
    metric = c(
      "target_universe", "representative_nodes", "GSE97742_mapped",
      "GSE97742_86_family_fdr_lt_0.05", "GSE41374_mapped",
      "GSE41374_robust_86_family_fdr_lt_0.05"
    ),
    value = c(
      nrow(contract), sum(contract$is_representative),
      sum(!is.na(gse97742$p_value)), sum(gse97742$fdr_bh_86_family < 0.05, na.rm = TRUE),
      sum(!is.na(wide$logFC.shifted_log2)), sum(wide$robust_86_family_fdr, na.rm = TRUE)
    )
  )
  write.table(summary, file.path(output_root, "analysis_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(output_root, "session_info.txt"))

  message(
    "Full-86 airway analysis complete: target_universe=", nrow(contract),
    " GSE97742_mapped=", summary$value[summary$metric == "GSE97742_mapped"],
    " GSE41374_mapped=", summary$value[summary$metric == "GSE41374_mapped"]
  )
}

if (sys.nframe() == 0L) run_full_analysis()
