#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
source(file.path("analysis", "R", "single_gene_upgrade_pipeline.R"))

read_module_config <- function(path = file.path("analysis", "config", "fixed_blood_module_source.yml")) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml package is required", call. = FALSE)
  yaml::read_yaml(path)
}

extract_fixed_modules <- function(config, output_path = file.path("data", "clean", "fixed_blood_modules.tsv")) {
  archive <- config$source_file
  if (!file.exists(archive)) stop("BloodGen3 source archive is missing", call. = FALSE)
  actual_sha <- system2("shasum", c("-a", "256", shQuote(archive)), stdout = TRUE)
  actual_sha <- strsplit(actual_sha[[1]], "[[:space:]]+")[[1]][1]
  if (!identical(actual_sha, config$sha256)) stop("BloodGen3 source archive checksum mismatch", call. = FALSE)

  extract_dir <- tempfile("bloodgen3_")
  dir.create(extract_dir)
  on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)
  utils::untar(archive, files = "BloodGen3Module/R/sysdata.rda", exdir = extract_dir)
  env <- new.env(parent = emptyenv())
  load(file.path(extract_dir, "BloodGen3Module", "R", "sysdata.rda"), envir = env)
  members <- env$Module_listGen3
  ann <- env$Gen3_ann
  out <- merge(members[c("Module", "Gene", "Function", "position")],
               ann[c("Module", "Module_color", "Cluster")], by = "Module", all.x = TRUE, sort = FALSE)
  names(out) <- c("module_id", "gene_symbol", "function", "position", "module_color", "cluster")
  out <- unique(out)
  out <- out[order(out$module_id, out$gene_symbol), ]
  rownames(out) <- NULL
  if (length(unique(out$module_id)) != as.integer(config$expected_modules)) stop("Unexpected number of BloodGen3 modules", call. = FALSE)
  if (nrow(out) != as.integer(config$expected_memberships)) stop("Unexpected number of BloodGen3 memberships", call. = FALSE)
  write_upgrade_tsv(out, output_path)
  out
}

assert_module_dictionary <- function(modules) {
  required <- c("module_id", "gene_symbol", "function", "position", "module_color", "cluster")
  if (length(setdiff(required, names(modules)))) stop("Module dictionary is missing columns", call. = FALSE)
  if (anyNA(modules$module_id) || anyNA(modules$gene_symbol) || any(!nzchar(modules$gene_symbol))) stop("Module membership is incomplete", call. = FALSE)
  if (anyDuplicated(modules[c("module_id", "gene_symbol")])) stop("Duplicate module-gene memberships", call. = FALSE)
  invisible(TRUE)
}

collapse_gene_symbols <- function(x) {
  x <- x[is.finite(x$log2FC) & is.finite(x$statistic) & !is.na(x$gene_symbol) & nzchar(x$gene_symbol), ]
  x <- x[order(x$gene_symbol, -abs(x$statistic)), ]
  x[!duplicated(x$gene_symbol), ]
}

read_gse38900_genome_effects <- function() {
  x <- utils::read.csv(file.path("GEO_data", "GSE38900_DEG_v2.csv"), row.names = 1L, check.names = FALSE)
  collapse_gene_symbols(data.frame(gene_symbol = rownames(x), log2FC = x$logFC, statistic = x$t,
                                   p_value = x$P.Value, fdr_bh = x$adj.P.Val, stringsAsFactors = FALSE))
}

read_entrez_genome_effects <- function(path) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("AnnotationDbi and org.Hs.eg.db are required", call. = FALSE)
  }
  x <- read_upgrade_tsv(path)
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(as.character(x$gene_id)),
                               keytype = "ENTREZID", columns = "SYMBOL")
  x$gene_symbol <- map$SYMBOL[match(as.character(x$gene_id), map$ENTREZID)]
  collapse_gene_symbols(data.frame(gene_symbol = x$gene_symbol, log2FC = x$log2FC,
                                   statistic = x$moderated_t, p_value = x$p_value,
                                   fdr_bh = x$fdr_bh, stringsAsFactors = FALSE))
}

score_fixed_modules <- function(effects, modules, cohort, minimum_coverage = 0.50, inter_gene_cor = 0.01) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("limma package is required", call. = FALSE)
  effects <- collapse_gene_symbols(effects)
  stats_vec <- effects$statistic
  names(stats_vec) <- effects$gene_symbol
  split_modules <- split(modules$gene_symbol, modules$module_id)
  module_ids <- names(split_modules)
  total_n <- lengths(split_modules)
  observed_n <- vapply(split_modules, function(g) sum(unique(g) %in% names(stats_vec)), integer(1))
  coverage <- observed_n / total_n
  evaluable <- coverage >= minimum_coverage & observed_n >= 2L

  out <- data.frame(
    cohort = cohort,
    direction = "RSV_minus_healthy",
    module_id = module_ids,
    function_name = modules[["function"]][match(module_ids, modules$module_id)],
    total_genes = total_n,
    observed_genes = observed_n,
    coverage = coverage,
    evaluable = evaluable,
    median_log2FC = NA_real_,
    camera_direction = NA_character_,
    p_value = NA_real_,
    fdr_bh = NA_real_,
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "function_name"] <- "function"
  for (i in seq_along(module_ids)) {
    genes <- intersect(unique(split_modules[[i]]), effects$gene_symbol)
    out$median_log2FC[i] <- if (length(genes)) stats::median(effects$log2FC[match(genes, effects$gene_symbol)], na.rm = TRUE) else NA_real_
    if (!evaluable[i]) next
    camera <- limma::cameraPR(stats_vec, index = which(names(stats_vec) %in% genes), inter.gene.cor = inter_gene_cor, sort = FALSE)
    out$camera_direction[i] <- as.character(camera$Direction[1])
    out$p_value[i] <- camera$PValue[1]
  }
  out$fdr_bh[evaluable] <- bh_family(out$p_value[evaluable])
  out
}

module_replication_summary <- function(cohort_results) {
  cohorts <- unique(cohort_results$cohort)
  wide <- reshape(cohort_results[c("module_id", "cohort", "evaluable", "median_log2FC", "fdr_bh")],
                  idvar = "module_id", timevar = "cohort", direction = "wide")
  eval_cols <- paste0("evaluable.", cohorts)
  effect_cols <- paste0("median_log2FC.", cohorts)
  fdr_cols <- paste0("fdr_bh.", cohorts)
  wide$evaluable_all_three <- apply(wide[eval_cols], 1L, function(z) all(z %in% TRUE))
  wide$direction_pattern <- apply(wide[effect_cols], 1L, function(z) {
    if (any(!is.finite(z))) return("not_estimable")
    if (all(z > 0)) return("all_up")
    if (all(z < 0)) return("all_down")
    "mixed"
  })
  wide$significant_cohort_n <- rowSums(sapply(wide[fdr_cols], function(z) is.finite(z) & z < 0.05))
  wide
}

map_frozen_genes_to_modules <- function(frozen, modules) {
  merged <- merge(frozen[c("entrez_id", "gene_symbol")], modules, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  merged <- merged[order(match(merged$gene_symbol, frozen$gene_symbol), merged$module_id), ]
  rownames(merged) <- NULL
  merged
}

main_fixed_blood_modules <- function() {
  config <- read_module_config()
  modules <- extract_fixed_modules(config)
  assert_module_dictionary(modules)
  effects <- list(
    GSE38900 = read_gse38900_genome_effects(),
    GSE77087 = read_entrez_genome_effects(file.path("results", "revised_main_gse77087", "cohort", "GSE77087_gene_effects.tsv")),
    GSE103842 = read_entrez_genome_effects(file.path("results", "revised_main_gse77087", "cohort", "GSE103842_gene_effects.tsv"))
  )
  scored <- do.call(rbind, lapply(names(effects), function(cohort) {
    score_fixed_modules(effects[[cohort]], modules, cohort,
                        minimum_coverage = config$minimum_module_coverage,
                        inter_gene_cor = config$inter_gene_correlation)
  }))
  rownames(scored) <- NULL
  root <- file.path("results", "single_gene_upgrade_v1", "modules")
  write_upgrade_tsv(scored, file.path(root, "cohort_module_effects.tsv"))
  replication <- module_replication_summary(scored)
  write_upgrade_tsv(replication, file.path(root, "module_replication_summary.tsv"))
  frozen <- read_upgrade_tsv(file.path("data", "clean", "bfff_86_gene_set_frozen.tsv"), col_classes = "character")
  mapping <- map_frozen_genes_to_modules(frozen, modules)
  write_upgrade_tsv(mapping, file.path(root, "frozen_86_gene_module_mapping.tsv"))
  message("MODULES total=", length(unique(modules$module_id)),
          " evaluable_all_three=", sum(replication$evaluable_all_three),
          " same_direction_all_three=", sum(replication$evaluable_all_three & replication$direction_pattern %in% c("all_up", "all_down")))
  invisible(list(scored = scored, replication = replication, mapping = mapping))
}

if (sys.nframe() == 0L) main_fixed_blood_modules()
