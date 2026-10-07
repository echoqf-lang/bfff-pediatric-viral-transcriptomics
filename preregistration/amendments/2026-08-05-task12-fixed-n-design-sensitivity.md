# Task 12 操作化修订：固定样本下的设计敏感度与精度

**Status:** FROZEN by the Git commit first containing this amendment, before the Task 12 formal simulation

**Scope:** 本修订只闭合原预注册 Task 12 尚未操作化的方差、队列内标准误、随机生成、判定和报告规则。不修改确认性主要判定，不执行正式模拟，不生成 MDE 或设计敏感度结果。

## 1. 设计定位与输入边界

- 这是公开数据已给定样本的 **fixed-N design-sensitivity assessment**，不是前瞻性样本量估算，也不是事后功效。
- 只允许读取以下与结局无关的设计输入：
  - `data/clean/sample_manifest_frozen.tsv`，SHA-256 `5434ee45d31e80a38db81696118b823e45b5aab2cd2dc628269ab29ee417a859`；
  - `results/tables/GSE105450_qc.tsv`，SHA-256 `e4417f3dd2d3fd0c83364a12c15f3b86dab21ee3c0bc117bd703aa0b6cf00e5e`；
  - `results/tables/GSE103842_qc.tsv`，SHA-256 `3ed3ef0903ccf962dd0fc327c85e37885519f43e20b9f68eb551e75dd9782f80`。
- 模拟代码不得读取任何基因效应、SE、P值、集合检验、Meta 结果、随机集合结果或 Task 11 结果；不得用观察到的目标基因结果选择残差 SD、效应或 τ。
- 确认队列与最终完整案例固定为：
  - GSE105450：122 人，RSV 89、健康 33；全部年龄、性别完整，无 QC 剔除；
  - GSE103842：73 人，RSV 61、健康 12；全部年龄、性别完整，1 人已按盲态 QC 剔除。
- 两队列使用与主分析一致的设计 `~ case_status + age_months + sex + technical_batch`；病例系数固定为 `case_statusRSV`，方向为 RSV 减健康。GSE105450 设计阶数为 6/6，GSE103842 为 13/13。

## 2. 队列内标准误的结果无关设定

- 不使用真实表达矩阵估计基因特异方差，也不使用 limma 观察到的后验 SE。
- 对队列 `i` 的病例系数，令 `h_i = sqrt([(X_i'X_i)^(-1)]_case,case)`，其中 `X_i` 只由上述冻结样本与协变量构建。队列内已知标准误设为 `SE_i = sigma * h_i`。
- 冻结设计杠杆项：GSE105450 `h=0.22229046974842828`；GSE103842 `h=0.34363953005869768`。
- 由于无结果独立的基因特异残差 SD 依据，不给出伪精确的单一方差。预先固定 `sigma={0.5, 1.0, 1.5}` 个 log2 表达单位：`1.0` 为主设计标尺，`0.5` 和 `1.5` 分别为低方差和高方差设计情景。这三者是结果无关的尺度敏感度网格，属于情景假设，不声称是经验估计。
- 模拟直接在队列系数层面运行，把上述 SE 视为固定已知设计参数；不模拟个体表达、残差方差估计或 limma 经验贝叶斯收缩。因此结果只是条件于这些 SD 情景的设计敏感度，不是精确的基因功效。

## 3. 数据生成过程与分析

- 综合真效应 `mu` 固定为 log2FC `{0, 0.10, 0.20, 0.30, 0.50}`；研究间标准差 `tau` 固定为 `{0, 0.05, 0.10, 0.20}`，单位均为 log2FC。
- 每个 `sigma × tau × mu` 组合顺序模拟 10,000 次，使用 R `RNGkind("L'Ecuyer-CMRG")`、单一起始种子 `20260812`，循环顺序固定为 `sigma` 升序、`tau` 升序、`mu` 升序、重复编号升序，不并行。
- 对每次重复和队列 `i`：先抽取 `u_i ~ Normal(0, tau^2)`，得队列真效应 `theta_i=mu+u_i`；再抽取 `y_i ~ Normal(theta_i, SE_i^2)`。两队列的 `u_i` 和采样误差相互独立。
- 每次重复对两个 `y_i, SE_i` 运行与主基因 Meta 一致的 `metafor::rma.uni(method="REML", test="knha")`，双侧 `alpha=0.05`，不运行 Wald 或固定效应替代分析。
- 该模拟的分析对象是“单个假定基因在两个确认队列的 REML/HKSJ Meta”；它不模拟 camera 竞争性集合检验、随机集合校准或 H1 联合判定，因此“发现概率”不得改写主集合结论。
- 任何拟合错误、非有限输出或非正标准误都使该情景失败并停止正式输出；不丢弃、不补抽、不换种子。

## 4. 结局、阈值与蒙特卡罗不确定性

- `discovery_probability = mean(P_REML_HKSJ < 0.05)`。`mu=0` 时该量称为 I 类错误率，不称为功效。
- 对 `mu>0`，同时报告 `direction_correct_discovery_probability = mean(P<0.05 and pooled_estimate>0)` 和错误方向显著比例；80%/90% 阈值仍按预注册的双侧发现概率判定。
- `mean_effect_CI_coverage = mean(CI_lower <= mu and mu <= CI_upper)`，覆盖对象是超总体平均真效应 `mu`，不是本次抽取的两个 `theta_i`。
- 每个比例同时报告基于 10,000 次的二项分布蒙特卡罗标准误 `sqrt(p*(1-p)/10000)` 和 95% Wilson 区间。
- 对每个 `sigma × tau` 和阈值 0.80/0.90，“最小情景效应”定义为预设的非零 `mu` 网格中第一个 `discovery_probability >= threshold` 的值。不插值、不外推、不做等调平滑；若网格内未达到，记为 `not_reached_within_prespecified_grid`，不报一个伪 MDE。

## 5. 输出命名与真实结果精度的隔离

- 原计划中 `mde_simulation.tsv` 和 `mde_by_tau.pdf` 的命名被本修订取代，因为离散情景网格不是连续 MDE 估计。正式执行时写入：
  - `results/power/fixed_n_design_sensitivity.tsv`；
  - `results/power/fixed_n_threshold_summary.tsv`；
  - `results/figures/fixed_n_design_sensitivity_by_tau.pdf`。
- 模拟完成并锁定后，才可以在独立步骤读取真实基因 Meta 结果，写入 `results/power/observed_precision_summary.tsv`。该表对未经 P 值筛选的全部可估主靶点汇总：95% CI 宽度的中位数、IQR、最小值和最大值；预测区间宽度的中位数、IQR和最大值；两队列 log2FC 同号的比例。
- 真实精度表不计算、不倒推、不标注事后功效，不把 CI 包含零解释为“功效不足”，不向模拟回填任何真实效应或 SE。
- Task 12 只描述固定设计在预设情景下的灵敏度与估计精度；不能救回或改写已经失败的 H1 联合判定。

## 6. 假设质量审计

| 输入 | 质量标签 | 用途/边界 |
|---|---|---|
| 最终样本数、病例/对照、完整案例、协变量设计 | Known / directly provided | 由冻结 manifest 和盲态 QC 确定 |
| log2FC 与 τ 网格 | Known / pre-registered | 不根据真实结果修改 |
| residual SD `sigma={0.5,1.0,1.5}` | Guessed / scenario-only | 仅供尺度敏感度，不当作经验方差 |
| 队列内 SE | Deterministic conditional calculation | 仅由 `sigma` 和冻结 `X` 计算 |
| 基因特异方差及 limma 收缩后 SE 分布 | Missing and intentionally unused | 禁止从目标结果估计；因此不声称精确功效 |
