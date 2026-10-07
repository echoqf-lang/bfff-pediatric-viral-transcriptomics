#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
required <- c("Biobase", "limma", "AnnotationDbi", "illuminaHumanv4.db")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)

panel <- read.delim("analysis/config/bfff_barrier_inflammation_repair_panel.tsv", check.names = FALSE)
input_root <- file.path("data", "derived", "barrier_inflammation_repair")
output_root <- file.path("results", "bfff_target_panel_airway_v1")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

map_and_collapse <- function(expression) {
  symbols <- AnnotationDbi::mapIds(
    illuminaHumanv4.db::illuminaHumanv4.db,
    keys = rownames(expression), keytype = "PROBEID", column = "SYMBOL", multiVals = "first"
  )
  keep <- !is.na(symbols) & symbols %in% panel$gene_symbol
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

run_97742_interaction <- function() {
  eset <- readRDS(file.path(input_root, "GSE97742_processed_eset.rds"))
  expression <- map_and_collapse(Biobase::exprs(eset))
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

  rows <- list()
  for (gene in rownames(expression)) {
    values <- expression[gene, rownames(pd)]
    wide <- reshape(
      data.frame(subject_key = pd$subject_key, infection = pd$infection, age_months = pd$age_months,
                 sex = pd$sex, phase = pd$phase, value = values),
      idvar = c("subject_key", "infection", "age_months", "sex"), timevar = "phase", direction = "wide"
    )
    wide <- wide[stats::complete.cases(wide), ]
    wide$delta <- wide$value.acute - wide$value.discharge
    wide$infection <- factor(wide$infection, levels = c("hRV", "RSVsi"))
    fit <- stats::lm(delta ~ infection + age_months + sex, data = wide)
    tab <- summary(fit)$coefficients
    rows[[gene]] <- data.frame(
      gene_symbol = gene, n_subjects = nrow(wide),
      rsv_minus_hrv_delta = tab["infectionRSVsi", "Estimate"],
      se = tab["infectionRSVsi", "Std. Error"], p_value = tab["infectionRSVsi", "Pr(>|t|)"],
      rsv_mean_delta = mean(wide$delta[wide$infection == "RSVsi"]),
      hrv_mean_delta = mean(wide$delta[wide$infection == "hRV"])
    )
  }
  out <- do.call(rbind, rows)
  out$fdr_bh_panel <- p.adjust(out$p_value, method = "BH")
  merge(panel, out, by = "gene_symbol", all.x = TRUE, sort = FALSE)
}

fit_41374 <- function(transform_name) {
  eset <- readRDS(file.path(input_root, "GSE41374_processed_eset.rds"))
  raw <- Biobase::exprs(eset)
  transformed <- switch(
    transform_name,
    shifted_log2 = log2(raw - min(raw, na.rm = TRUE) + 1),
    asinh_robust = asinh(raw / stats::median(abs(raw), na.rm = TRUE)),
    stop("Unknown transform", call. = FALSE)
  )
  expression <- map_and_collapse(transformed)
  pd <- Biobase::pData(eset)
  meta <- data.frame(
    age_months = as.numeric(pd[["age (month):ch1"]]),
    sex = factor(pd[["gender:ch1"]]),
    group = factor(pd[["infection:ch1"]], levels = c("healthy", "RSV")),
    row.names = rownames(pd)
  )
  keep <- stats::complete.cases(meta)
  meta <- meta[keep, , drop = FALSE]
  expression <- expression[, rownames(meta), drop = FALSE]
  design <- stats::model.matrix(~ age_months + sex + group, data = meta)
  fit <- limma::eBayes(limma::lmFit(expression, design))
  tab <- limma::topTable(fit, coef = "groupRSV", number = Inf, sort.by = "none")
  tab$gene_symbol <- rownames(tab)
  tab$transform <- transform_name
  tab$fdr_bh_panel <- p.adjust(tab$P.Value, method = "BH")
  tab[, c("gene_symbol", "transform", "logFC", "P.Value", "fdr_bh_panel")]
}

interaction <- run_97742_interaction()
write.table(interaction, file.path(output_root, "GSE97742_RSV_vs_hRV_recovery_interaction.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

gse41374 <- rbind(fit_41374("shifted_log2"), fit_41374("asinh_robust"))
write.table(gse41374, file.path(output_root, "GSE41374_two_transform_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

wide <- reshape(gse41374, idvar = "gene_symbol", timevar = "transform", direction = "wide")
wide$direction_concordant <- sign(wide$logFC.shifted_log2) == sign(wide$logFC.asinh_robust)
wide$robust_panel_fdr <- wide$direction_concordant & wide$fdr_bh_panel.shifted_log2 < 0.05 &
  wide$fdr_bh_panel.asinh_robust < 0.05
wide <- merge(panel, wide, by = "gene_symbol", all.x = TRUE, sort = FALSE)
write.table(wide, file.path(output_root, "GSE41374_robust_target_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Target panel complete: genes=", nrow(panel),
        " interaction_fdr=", sum(interaction$fdr_bh_panel < 0.05, na.rm = TRUE),
        " gse41374_robust=", sum(wide$robust_panel_fdr, na.rm = TRUE))
