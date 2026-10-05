# 数据来源与字段说明

这些是已有检查的精简结果，没有为撰写交接 note 重新计算。程序、解析分布和独立参考计算的结果分别标注，不能互相替代。

| 文件 | 来源及用途 |
|---|---|
| `free-limit-summary.csv` | 无耦合小晶格生产程序检查；每个 seed 单独一行。`energy_real` 是整个晶格的 Majorana 能量 |
| `free-limit-sampling-moments.csv` | 同一无耦合检查的位移均值和位移平方均值；解析独立样本标准误差仅作为参考尺度 |
| `free-limit-structure-moments.csv` | 同一检查的非零允许动量声子结构因子，归一化为 $S_{zz}/N_c$；不是自旋 $S^zS^z$ |
| `fixed-path-fock-reference.csv` | 固定指定声子路径后，费米子投影计算与独立多体直接计算的差异；不验证路径采样是否充分 |
| `coupled-ed-qmc-comparison.csv` | 小型费米子–声子模型与独立声子基底参考的逐项比较。不是完整生产蜂窝格点扫描 |
| `forced-oscillator-frequency-scan.csv` | `../validation/forced-oscillator/analytic_mixing.py` 产生的独立 Python 高斯分布采样实验；没有调用生产代码 |
| `adiabatic-L6-summary.json` | $L=6$、$\lambda=2.4$ 固定规范场的静态最低已找到解；能量为整个 36 原胞系统的值 |

## 受力谐振子频率扫描

`forced-oscillator-frequency-scan.csv` 中 `tau` 是整条路径长度 $2\Theta$，`intervals` 是虚时间区间数，离散步长为 0.1。`exact_center_mean` 和 `exact_center_variance` 对应相同离散路径分布。`stationary_independence_acceptance` 是此解析高斯模型中自由路径独立提议的平稳平均接受概率，不能当作大晶格的测量接受率。

`analytic_free_start_*` 和 `analytic_equilibrium_start_*` 均为独立 Python Metropolis–Hastings 实验；每个实验运行 5000 次提议。该 CSV 在首次发布时将这些列错误命名为 `core_*`，本次已更正来源标签，数值保持原值。

使用生产采样类的另一组受力谐振子检查是每次 1000 个提议。原始摘要在 `../validation/forced-oscillator/core-mix-free-seed40956.txt`、`core-mix-equilibrium-seed40955.txt` 和 `core-guided-seed40956.txt`。最后一个程序分支单独统计非零均值提议的接受次数，应读取 `guided_accept`；其中 `accepted=0 proposed=0` 是未调用的另一采样分支的计数。

## 独立小模型对照

`coupled-ed-qmc-comparison.csv` 中 `E` 表示 Majorana 能量；`X0`、`X1`、`X2` 是三个对应键组的平均位移；`X20`、`X21`、`X22` 分别是这三个键组的平均位移平方，不是第 20、21、22 个变量，也不是平均位移再平方。`stderr` 取 100 和 200 个测量值分块所得的两个标准误差中较大者。`difference_in_stderr` 是与参考值的差除以这一误差；预先规定的比较门槛是 6 个标准误差，不能把通过门槛解释为所有参数的正确性保证。
