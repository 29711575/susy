# 统一计算工具

入口为 `python3 scripts/qmc.py`，在 quick1010 根目录运行。需要 Python 3.9+、
NumPy；诊断 PNG 另需 Matplotlib。全部操作读取既有数据或 PBS 状态，
不会提交、取消、停止计算或修改 seed。无需重新编译现有计算程序。
新输出默认存放在上海日期的 `audit/MMDD/qmc-tools/` 下，也可指定 `--out`；
输出目录必须不存在，避免覆盖历史分析。

## 从原始块重做最终分析

```bash
python3 scripts/qmc.py analyze PATH/mc_results.json --out audit/1023/my-analysis
```

输入也可以是 `mc_blocks.csv`、输出前缀或含 `mc_blocks.csv` 的目录。
程序独立读取原始复数 `sum_phase` 和 `sum_phase_O`，不启动蒙特卡洛，
也不依赖原先最终汇总是否成功。输出新的 `analysis.json` 与 `estimates.csv`，
原文件保留。新 JSON 是离线分析格式，不是原生产 JSON 的字段完全相同的副本。

支持全部原始观测量及按列名重建的全动量 VBS 谱、实空间关联、背景连接减法、
复数连通 susceptibility、幅度 VA/B4/chiA 和有限场派生量。
selected/both 模式的谱及投影谱可恢复；候选最大值需要原 JSON 的参与者元数据，
缺失时明确列入 `derived_omissions`，不猜测候选集合。
原始块未保存的构型、RNG、测量量或运行计数不会由统计重建出来。

新增总能量六列可用相同命令重做联合统计，自动输出每晶胞分项、
`energy_matter_link`、`energy_phonon_link`及`energy_sum_residual`。
历史块若保存`X2_x/y/z`及实际k，可重建
`energy_phonon_potential_center_reconstructed`；另有逐键`BondK/BondXK`及明确
`gauge_u`时可重建中心裸matter/耦合项。缺少动能信息的旧块不能重建总能量。
新能量与历史中心量的测量位置不同，见[TOTAL_ENERGY.md](TOTAL_ENERGY.md)。

对运行中、停止、失败或缺少完成证明的任务，显式加 `--allow-prefix`：

```bash
python3 scripts/qmc.py analyze PATH/mc_blocks.csv --allow-prefix \
  --cut-sweeps 1000 --reblock 1 2 4 8 16 --out audit/1023/my-prefix-analysis
```

读取文件打开时的固定字节上限，仅取完整换行记录，保存使用字节数和 SHA256。
末尾尚未写完的记录只在 `--allow-prefix` 下排除并记录；完整记录损坏、缺块、重复块、
非有限数、错误 phase 分类、完成结果与块计数不一致均拒绝。
所有弃样/终点必须对齐原始块边界；不从聚合块猜测其中的逐 sweep 数据。
`--end-sweeps N` 可固定共同终点。

身份优先使用实际保存的结果、argv、point 和 matter-trial 记录；不把今天的
默认参数套到旧任务。确实缺少的参数可通过 `--metadata FILE.json` 补充；
与已记录身份冲突会拒绝。只需原始量时可用 `--raw-only`。

多链按原始样本数合并复数分子/分母，然后作非线性变换；误差在各独立链内删组，
保持其原来的总样本分配。相邻块只在链内合并；粗分块的短尾仍保留，使用实际
样本比例的加权 pseudovalue 方差（与既有 g0 分析同一方法）。相等组时严格退化为
原联合 delete-one-block jackknife；不等组是有限样本删组估计。
组数不足、相位分母不定或某派生量不定，会明确标注并扣留对应误差。
禁止重复文件、重复 seed、祖先块与已包含祖先的续算块一起合并。

## 任务总表

```bash
python3 scripts/qmc.py registry refresh --out audit/1023/my-registry
python3 scripts/qmc.py registry query --registry audit/1023/my-registry/registry.json --L 12 --state R
python3 scripts/qmc.py registry query --registry audit/1023/my-registry/registry.json --job 90057 --lambda 2.7
```

默认只读采集 20022 和 20015，可重复 `--port` 指定端口。SSH 主机和密钥可用
`--host`、`--ssh-key` 设置。输出包含机器可读 `registry.json`、逐链 `tasks.csv`、
按物理参数/来源/控制组划分的 `cells.csv` 和两端原始任务快照。
每条记录保留端口、PBS 子任务 ID、局部/全局下标、参数、seed、输出前缀、
完整落盘块的进度、身份映射来源和状态采集时间。尚未产生输出的任务保留空值；
缺失身份逐项列在 `identity_missing`。

刷新时用 `--previous OLD/registry.json` 延续之前的记录。此次没有列出的任务保留
上次观察时间并标为历史观察，不推断其已完成或已取消。PBS 数组中列出的 X/F
终止子任务也保留，查询 `--state R` 可只取本次观察到的运行任务。
运行中的启动占位 JSON 不用作实时数值失败计数，最近接受率记录单独列出。

历史包装脚本支持 PBS start 记录及唯一 point 目录；Γ 数组的 `CHAIN_OFFSET`
从与 PBS 作业名匹配的脚本中读取明确整数字面量，并记录映射依据，避免把两个
数组的局部下标误当成同一个参数点。无法确定的映射保留错误，不能用于取数。

提交记录也可导入，显式保持其历史身份：

```bash
python3 scripts/qmc.py registry import PATH/submission.json --points PATH/points.json \
  --previous audit/1023/my-registry/registry.json --out audit/1023/merged-registry
```

当前导入格式是既有 `port/root/groups[{job,indices,exit}]/utc` 提交回执及
带唯一 `index` 的 points 列表。只导入实际成功提交的下标，不能把整个 points 文件
都认作已提交。其他历史格式需先转换为这个明确格式。

## 从任务表取数并统一诊断

```bash
python3 scripts/qmc.py registry fetch --registry audit/1023/my-registry/registry.json \
  --port 20022 --job 89979 --lambda 1 --state R --out audit/1023/my-snapshot
python3 scripts/qmc.py analyze --selection audit/1023/my-snapshot/selection.json \
  --allow-prefix --out audit/1023/my-offline-analysis
python3 scripts/qmc.py diagnose --selection audit/1023/my-snapshot/selection.json \
  --allow-prefix --compare-cuts 100 200 400 --out audit/1023/my-diagnostics
```

`fetch` 只读冻结完整 CSV 前缀，保留远端文件 SHA256、采集起止时间、压缩归档，
下载后核对校验和。运行中数据会继续增长；这是一份带时间范围的快照。

`diagnose` 与 `analyze` 共用同一个读取器和统计实现，生成逐 seed 原始块比值及
累积比值 PNG、共同窗口的 rank/folded split Rhat、块记录尺度 bulk ESS、
联合重分块、逐 seed 统计以及显式指定弃样窗口的敏感性结果。
默认检查已有的 Xx/Xy/Xz、Xminus、matter energy、总能量和声子动能/势能，并在已记录固定 q 时加入
Szz 和复 Delta；可重复 `--observable NAME` 选择原始列或 Xplus/Xminus。

Rhat 使用共同相位比的影响序列，以整个共同窗口为中心，不逐链扣除自身均值。
它基于存下来的聚合块，不是逐 sweep 自相关诊断。每条链取共同长度，排除的
尾部长度记录在报告中；重分块统计与图仍保留请求窗口的全部数据。
链间 Rhat>1.01 时扣留 ESS，常数序列与不足数据单独标注。
没有任何输出自动认证平衡、选择最佳弃样窗口、停止任务或减少 seed。

## 验证与本次示例

```bash
python3 -m unittest discover -s scripts/tests -p test_qmc_tools.py -v
```

本次实际任务表、L12 四链快照/分析/PNG、旧最终汇总失败的恢复，以及与原 C++
结果的交叉验证见 `audit/1023/calculation-tools/REPORT_CN.md`。
