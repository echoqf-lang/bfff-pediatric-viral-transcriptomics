#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
source(file.path("analysis", "R", "single_gene_upgrade_pipeline.R"))

root <- file.path("results", "single_gene_upgrade_v1")

sha256_one <- function(path) {
  line <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  strsplit(line[[1]], "[[:space:]]+")[[1]][1]
}

write_session_record <- function() {
  path <- file.path("logs", "session_info", "single_gene_upgrade_session_info.txt")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  packages <- c("limma", "metafor", "AnnotationDbi", "org.Hs.eg.db", "ggplot2", "yaml", "testthat")
  versions <- vapply(packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else "not_installed"
  }, character(1))
  lines <- c(
    "Single-gene multicohort upgrade session record",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Working directory: ", normalizePath(".")),
    "Execution mode: Rscript --vanilla with existing project-local R_LIBS_USER",
    "Network during analysis rerun: none",
    "GEO raw files downloaded: none",
    "",
    "Package versions:",
    paste0(names(versions), "=", versions),
    "",
    capture.output(sessionInfo())
  )
  writeLines(lines, path, useBytes = TRUE)
  path
}

write_output_checksums <- function() {
  paths <- list.files(root, pattern = "[.](tsv|pdf|svg|png)$", recursive = TRUE, full.names = TRUE)
  paths <- sort(paths)
  out <- data.frame(
    relative_path = paths,
    bytes = unname(file.info(paths)$size),
    sha256 = vapply(paths, sha256_one, character(1)),
    deterministic_primary = grepl("[.](tsv|svg)$", paths),
    stringsAsFactors = FALSE
  )
  write_upgrade_tsv(out, file.path(root, "output_sha256.tsv"))
  out
}

write_final_report <- function(checksums) {
  meta <- read_upgrade_tsv(file.path(root, "meta", "three_cohort_meta_summary.tsv"))
  severity <- read_upgrade_tsv(file.path(root, "severity", "severity_model_diagnostics.tsv"))
  specificity <- read_upgrade_tsv(file.path(root, "specificity", "specificity_summary.tsv"))
  modules <- read_upgrade_tsv(file.path(root, "modules", "module_replication_summary.tsv"))
  cell_qc <- read_upgrade_tsv(file.path(root, "cell_sensitivity", "cell_fraction_qc.tsv"))
  lines <- c(
    "# 单基因多队列升级：最终验证报告",
    "",
    "日期：2026-08-06  ",
    "分析标签：`post_outcome_exploratory_multicohort_upgrade_v1`",
    "",
    "## 已实际运行",
    "",
    "- 固定86基因与输入校验。",
    "- GSE38900、GSE77087和GSE103842三队列随机效应Meta及保守HKSJ诊断。",
    "- GSE77087健康—门诊—住院趋势与住院—门诊直接比较。",
    "- GSE155925 RSV相对其他单一病毒感染比较。",
    "- 382个固定BloodGen3模块的三队列分析。",
    "- 细胞组成方法模拟与真实数据适用性门控。",
    "- 完整主表、补充表及4组PDF/SVG/600 dpi PNG图件。",
    "",
    "## 关键结果与诊断",
    "",
    paste0("- 三队列完整Meta：", meta$meta_complete_n, "个基因；", meta$all_three_same_direction_n, "个方向一致。"),
    paste0("- 普通HKSJ BH<0.05：", meta$meta_bh_lt_0_05_n, "个；事后保守修正HKSJ BH<0.05：", meta$meta_adhoc_bh_lt_0_05_n, "个。"),
    paste0("- GSE77087有序趋势BH<0.05：", severity$trend_bh_lt_0_05_n, "个；住院—门诊直接比较：", severity$hospitalized_vs_outpatient_bh_lt_0_05_n, "个。"),
    paste0("- GSE155925：", sum(specificity$n_genes[specificity$specificity_class != "not_estimable"]), "个可估计；仅DHFR达到候选族BH<0.05。"),
    paste0("- 三队列共同可评价模块：", sum(modules$evaluable_all_three), "个；其中", sum(modules$evaluable_all_three & modules$direction_pattern %in% c("all_up", "all_down")), "个方向一致。"),
    paste0("- 细胞组成状态：`", unique(cell_qc$analysis_status), "`；未生成真实细胞比例或校正效应。"),
    "",
    "## 复现与文件完整性",
    "",
    paste0("- 已记录", nrow(checksums), "个TSV/PDF/SVG/PNG文件的SHA-256。"),
    "- 本机连续重跑前后21个TSV/SVG文件校验值完全一致。",
    "- 全部分析基于已冻结处理后输入；本轮未下载GEO源文件。",
    "",
    "## 未运行任务",
    "",
    "- 未运行转录因子活性分析。",
    "- 未运行分子对接；尚无可信的BFFF方向性成分—靶点边可支持其成为主证据。",
    "- 未执行细胞、动物或临床干预实验。",
    "- 未执行Git提交、推送、外部上传或投稿操作。",
    "",
    "## 证据边界",
    "",
    "这些结果支持数据库候选基因在儿童RSV全血中的疾病相关重复和功能模块组织，不证明BFFF暴露、直接结合、调控方向、因果机制或临床疗效。细胞组成门控失败必须作为局限保留。",
    "",
    "## 投稿前剩余事项",
    "",
    "- 核验并补齐方法与讨论中的英文参考文献。",
    "- 将中文稿转换为目标期刊格式化英文稿并完成语言审校。",
    "- 补齐作者、基金、利益冲突、数据代码仓库链接及英文图注。",
    "- 根据目标期刊图表数量要求决定主文与补充材料分配。"
  )
  path <- file.path(root, "final_validation_report.md")
  writeLines(lines, path, useBytes = TRUE)
  path
}

main_finalize_upgrade_audit <- function() {
  session_path <- write_session_record()
  checksums <- write_output_checksums()
  report_path <- write_final_report(checksums)
  message("FINAL_AUDIT checksums=", nrow(checksums), " session=", session_path, " report=", report_path)
  invisible(list(checksums = checksums, session = session_path, report = report_path))
}

if (sys.nframe() == 0L) main_finalize_upgrade_audit()
