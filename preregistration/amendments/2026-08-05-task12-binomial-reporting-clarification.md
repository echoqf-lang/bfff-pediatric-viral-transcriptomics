# Task 12 协议审查澄清：概率分母、Wilson 区间与 k=2 报告

**Status:** FROZEN by the Git commit first containing this clarification, before the Task 12 formal simulation

**Scope:** 本澄清只补齐 Task 12 已冻结 DGP 的概率汇总和 k=2 输出合同；不改变样本、效应网格、τ网格、残差 SD 情景、随机生成或主模型，不运行正式模拟。

## 固定概率分母

- `discovery`、`direction_correct_discovery`、`wrong_direction_discovery` 和 `mean_effect_ci_coverage` 的场景级分母全部固定为 10,000；不丢弃重复、不使用非缺失个数作可变分母。
- `mu=0` 时，`discovery` 的解释固定为 I 类错误率；此时没有预设真方向，因此方向正确/错误发现的成功数和概率固定报告为 `NA`，但 `denominator` 字段仍固定写 10,000。不允许根据结果改成其他规则。

## 二项比例汇总

- 对成功数 `x` 和 `n=10000`，固定 `p=x/n`，蒙特卡罗标准误为 `sqrt(p*(1-p)/n)`。
- 95% Wilson 区间使用无连续性校正的标准公式，固定 `z=qnorm(0.975)`：
  - `center=(p+z^2/(2*n))/(1+z^2/n)`；
  - `half=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n)`；
  - 下界 `max(0, center-half)`，上界 `min(1, center+half)`。
- `x=0`、`x=5000` 和 `x=10000` 均使用同一公式，不加 0.5、不使用正态 Wald 比例区间。

## k=2 的强制输出边界

- 每次 Meta 输出必须显式写入 `k=2` 和 `hksj_df=1`。
- 所有场景表和图注必须声明：只有两个队列，τ²估计、HKSJ 发现率和阈值可能不稳定，仅用于条件设计敏感度描述。
- 不得因为 k=2 的结果不稳定而切换至 Wald、固定效应、另一个 τ² 估计量，或对非单调的发现概率做等调/曲线平滑。按预设网格原样报告。
