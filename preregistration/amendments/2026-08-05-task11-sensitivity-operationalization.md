# Operational amendment: Task 11 敏感性与探索分析

**Status:** FROZEN upon the initial Git commit containing this file，且早于本文件所述任何效应模型、P值或集合结果的计算或查看

**Parent preregistration:** `docs/science-superpowers/preregistrations/2026-08-04-bufei-fanggan-rsv-validation.md`

**Reason:** 原预注册已指定分析方向，但 processed-only 输入需要在效应建模前冻结平台特异处理、缺失协变量、对比、多重性和不可实施项。本修订不改变核心 H1、确认性队列或主要判定。

## 输入与禁止事项

- 冻结样本清单：`data/clean/sample_manifest_frozen.tsv`，SHA-256 `5434ee45d31e80a38db81696118b823e45b5aab2cd2dc628269ab29ee417a859`。
- 扩展靶点集：`data/clean/bfff_targets_sensitivity.tsv`，1,270个唯一Entrez基因，SHA-256 `573c5f577dc071af45289ff68e051d762f7b4403288dc78267e62dab8755abcd`。
- 本任务不读取CEL、FASTQ、SRA或GEO `RAW.tar`档案，不联网下载。GSE155925允许读取作者提供的本地基因级未标准化计数矩阵；它不是FASTQ或SRA源文件。
- 在本修订冻结前不拟合或查看Task 11效应、P值、camera、Meta或富集结果。不得按任何结果改变样本、靶点、协变量、背景或分析版本。

## 1. GSE188427替换敏感性

- 固定纳入198个独立Day 1/健康样本：147个RSV病例、51个健康对照。7个arm冲突样本仅影响住院/门诊身份；其中3个Day 1样本仍可按一致的disease字段进入overall病例—健康比较，但不得进入严重程度分层。
- series matrix已由作者使用RMA并按`ENTREZG, version 20` alternative CDF汇总。18,604个行ID按`<EntrezID>_at`去除后缀直接映射；不再次RMA或分位数标准化，不使用标准GPL25336注释重新映射。
- 只使用198个预定样本做盲态processed-only QC，传入QC函数的对象只能含表达矩阵和样本ID，不得含case status、severity、arm或其他表型。任何非有限表达停止队列；其余盲态指标固定为缺失率、样本中位数、IQR、前5,000个高变基因的样本中位Spearman相关、前5 PC标准化距离和`limma::arrayWeights`对数。基因按全样本方差降序、方差并列时按规范化Entrez ID字典序选前5,000；PCA固定为`t(top5000)`、`center=TRUE`、`scale.=FALSE`，取最多前5 PC并将各PC得分除以相应`sdev`后计算欧氏距离。`arrayWeights`固定使用全表达矩阵和仅截距设计`matrix(1,n_samples,1)`。每项按`abs(x-median)>5*scaled MAD`标记，`stats::mad(constant=1.4826)`；MAD=0时仅非零偏差标记；至少2项标记才排除。指标构造和排除门复用Task 5实现，除无Detection P指标外不提供替代版本。
- 月龄和性别对198人全部缺失；`A3733_...`为逐样本唯一array描述而非技术批次。模型固定为未调整`~ case_status`，病例系数方向为RSV减健康。
- 该分析是探索性的、非H1、非完全规格替换；即使结果支持，也不能称为按主模型协变量调整后的独立重复。

## 2. GSE105450仅住院敏感性

- 直接复用已冻结的`GSE105450_expression.rds`盲态QC矩阵，不重新QC或标准化。
- 固定89人：56个住院RSV病例和33个独立健康对照；33个门诊RSV不进入该敏感性。
- 模型固定为`~ case_status + age_months + sex + technical_batch`；技术批次仅在设计满秩时保留，预检为6/6满秩。病例系数方向为RSV减健康。

## 3. GSE103119探索分析

- 固定31人：11个逐样本病原标签严格等于RSV的单感染病例、20个健康对照。方向为RSV减健康。
- 作者背景扣除和平均信号缩放后的series matrix只作`log2(pmax(intensity,1))`，不做第二次阵列间标准化。
- 用本地non-normalized processed表的Detection P；探针在保留样本中`Detection P<0.05`比例至少10%才视为可检测。
- 盲态QC、scaled-MAD常数、两指标排除和非有限值停止规则完全沿用Task 5。用本地官方GPL10558注释映射当前唯一Entrez；多重、空或过期映射删除；同一基因选择全部保留样本平均表达最高的探针，不使用状态、效应或P值。
- 模型固定为`~ case_status + age_months + sex`，预检4/4满秩；无可靠技术批次字段，不从表达结果推断批次。

## 4. GSE155925跨病毒探索分析

- 固定48个单病毒住院儿童：31个RSV、17个其他单病毒；10个病毒阴性住院儿童和6个混合病毒样本不进入模型。当前无可用严重程度梯度，不运行严重程度模型。
- count matrix列按`Case 1`至`Case 64`与series matrix title/GSM顺序严格匹配。基因行按最后一个冒号后的Ensembl ID解析，先去除可选版本后缀；必须满足`^ENSG[0-9]+$`。规范化后重复Ensembl行在过滤前按样本求和；显示symbol不参与身份或筛选。
- 仅对48个预定样本建立`edgeR::DGEList`；用完整设计执行`filterByExpr(dge, design=design)`，背景固定为过滤后唯一Ensembl行。随后`dge[keep,,keep.lib.sizes=FALSE]`、`calcNormFactors(method="TMM")`、`limma::voom`和经验贝叶斯拟合；不对计数作log2伪计数后直接线性回归。
- RNA-seq靶点单位固定为唯一Entrez。先在完整1,270集alias表中全局删除任何映射到多于1个Entrez的Ensembl alias；当前已知共享alias `ENSG00000223572`和`ENSG00000237289`均同时指向Entrez 1159与548596，必须删除，不能重复计权或按主/扩展集分别选择。
- 对每个Entrez，仅在去歧义且通过`filterByExpr`的候选Ensembl中选择一行。盲态排序量固定为过滤后TMM对象的`edgeR::cpm(dge, log=TRUE, prior.count=0.25, normalized.lib.sizes=TRUE)`在全部48个预定样本的行均值；均值最高者入选，并列按规范化Ensembl ID字典序。一个Entrez最多进入camera一次。
- 先用完整1,270集生成一张冻结选择表，520主集只能按Entrez从同一表取子集，不得另行选择较有利Ensembl。覆盖按有入选行的唯一Entrez计数，且必须严格等于传给camera的去重index长度。
- `batch`固定拆分为GEO原始字段`hospital_batch`与`enrollment_batch`。两批次的每个水平在RSV和其他病毒组内均有样本，完整设计预检为6/6满秩，故模型固定为`~ virus_group + age_months + sex + hospital_batch + enrollment_batch`。
- 系数方向固定为单RSV减其他单呼吸道病毒；它不是RSV减健康效应。

## 5. 靶点覆盖、背景和集合检验

- 主靶点集仍为520个唯一Entrez基因；扩展敏感性集固定为上述1,270个Entrez基因。RNA-seq以冻结Ensembl aliases连接靶点，不能映射的Entrez靶点保持未检测，不按symbol补配。
- 每个队列、每个靶点集分别计算实际可检测覆盖。可检测目标少于10个或少于该完整靶点集50%时标记为不可解释，不运行或解释camera；不得用其他队列靶点数补齐。
- camera背景固定为该模型实际保留、唯一映射且通过检测/表达过滤的全部基因；传入完整表达/voom对象、完整设计和预定病例系数。所有Task 11集合检验只允许`limma::camera(..., directional=TRUE, inter.gene.cor=NA_real_)`，由camera估计集合内相关；禁止固定`inter.gene.cor=0.01`、`cameraPR`或按结果选择camera版本。检验双侧；方向由病例/指定组减参照组的moderated统计量确定。
- 扩展靶点集在GSE105450与GSE103842的既定主模型中重复camera，并用模型实际独立样本数平方根加权的带方向Stouffer合并；不重估或改变确认性主靶点判定。

## 6. 两队列敏感性组合与Meta

- 仅住院组合固定为GSE105450住院模型与GSE103842主模型；替换组合固定为GSE188427 overall未调整模型与GSE103842主模型。
- GSE188427与GSE105450绝不进入同一个“独立”Meta或Stouffer组合。权重均为各模型QC和完整案例后独立受试者数平方根。
- 两个组合各自对主靶点运行双侧camera和带方向Stouffer。基因级Meta只纳入组合内两个队列均可估的全部主靶点，不按队列P值筛选；沿用REML+Hartung–Knapp，另报REML+Wald和固定效应，方向统一为RSV减健康。k=2不运行留一法，异质性谨慎描述。
- signed Stouffer固定如下。队列`i`的双侧camera P为`p_i`：Direction=`Up`时`z_i=+qnorm(p_i/2, lower.tail=FALSE)`，Direction=`Down`时取其负值；其他方向值停止。权重`w_i=sqrt(n_i)`，其中`n_i`为该模型经QC和协变量完整案例后最终独立受试者数。综合`Z=sum(w_i*z_i)/sqrt(sum(w_i^2))`，双侧`P=2*pnorm(abs(Z), lower.tail=FALSE)`；不得改用原始计划样本数、单侧P或无符号P。

## 7. 多重性与完整报告

- 集合层分为三个预设分析族并分别Holm校正：S1为GSE105450住院camera、GSE188427 camera及两个相应的两队列Stouffer共4项；S2为1,270集在GSE105450、GSE103842的camera及其Stouffer共3项；E1为GSE103119和GSE155925各自的520集与1,270集camera共4项。
- 基因级GSE105450住院组合和GSE188427替换组合分别构成完整主靶点分析族，各自在全部可估基因内作BH；两族不合并，也不选择较有利版本。
- 同时报告原始精确P值与校正P值、效应、区间、覆盖率、样本数和不可解释项。所有预设版本均进入结果，不因方向或显著性省略。

## 8. CBC、细胞比例与主结论边界

- 本地五个相关series matrix没有CBC、WBC、绝对白细胞分类或实测细胞比例字段，也没有冻结的儿科全血反卷积参考。本Task 11不实施CBC或细胞比例调整。
- 不因缺少细胞调整从表达结果选择替代代理变量；以后若引入成人参考，必须另行事前冻结并仅标记探索性。
- 本Task 11全部为敏感性或探索性结果，不能救回已经失败的H1联合判定，不能证明补肺防感方疗效、表达逆转、直接靶向或因果机制。
