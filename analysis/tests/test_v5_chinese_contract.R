target <- file.path("manuscript", "draft_v5_中文知识发现版.md")

stopifnot(file.exists(target))
txt <- paste(readLines(target, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

required <- c(
  "# 补肺防感方覆盖候选基因的跨情境转录组映射",
  "## 摘要",
  "## 1 引言",
  "## 2 材料与方法",
  "## 3 结果",
  "## 4 讨论",
  "## 5 结论",
  "纤毛功能受损",
  "炎症激活",
  "时相依赖性修复",
  "85/86",
  "24/86",
  "16个",
  "GSE97742",
  "GSE41374",
  "所有纳入分析的队列均未包含补肺防感方暴露"
)

for (item in required) {
  stopifnot(grepl(item, txt, fixed = TRUE))
}

forbidden <- c(
  "补肺防感方调控了",
  "补肺防感方逆转了",
  "补肺防感方抑制了",
  "证实了补肺防感方的作用机制",
  "验证了补肺防感方的疗效"
)

for (item in forbidden) {
  stopifnot(!grepl(item, txt, fixed = TRUE))
}

message("V5 Chinese manuscript contract passed")
