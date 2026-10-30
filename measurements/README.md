# 测量材料与独立复现

本目录包含不同日期、版本和截断口径的冻结测量；不会随 HPC 运行自动更新。
当前入口为 [2026-10-30 序参量快照](20261030/README.md)：L6/L9、频率比较、
全 q / 四 q 定义，以及可离线复算的原始矩与代码。

## 历史 2026-10-10 快照

以下四组是2026-10-10已采集的数据。普通生产链包含已停止但过截断的 partial
数据；complete不是热化证明。历史极化率解读见
[归档交接报告](../docs/history/HANDOFF_REPORT_BEFORE_20261030.md)。

- `L9_tau8_lambda/`：203条有效链、27个λ，τ8，cut500；源279条计划链中76条未过截断。`allowed-grid.csv`是81个允许q；峰曲线搜索80非零q。gzip JSON含每链 `phase_sum` 和 `chi_numerator_sums`。
- `L9_L6_tau_scan/`：L9原始340链与80条有效追加链，85点；L6原始40条与18条有效细扫链，13点。L9τ>6 cut1000、其他0；L6 cut3000；本目录总478有效链、2572730保留sweep。`grid-chi.csv`包括全部q；`L9-original-curves.csv`保留不加补充链的参考。
- `L9_tau8_heatmaps/`：前一τ8λ快照的六边形第一BZ，统一色阶、包含Γ，星号为裸Dirac差参考。每个色块对应最近允许动量单元，不是额外off-grid数据。
- `connected_decomposition/`：分别读取两个数据集做Cfull/Cdisc/χ分解，125个参数点。重叠seed属于不同cut/范围的两次诊断，不能合并成新增独立样本。

`chi_numerator_sums`每q按 `ChiD_z, ChiIconj_z, ChiRaw_zz`排列，各字段是实部/虚部两项；整数组shape为(Nc,3,2)。它们是保留样本的sum(phase×O)，不是已经phase重加权的均值。`phase_sum`是相同样本sum(phase)。`count`是保留sweep，`raw_count`含截断前部分。

对一个物理参数点，先把各链上述和逐元素累加，再共同除以phase_sum；`chi=<Raw>−Nc<D><Iconj>`。不能再按count乘这些“和”，否则重复加权；从逐链均值重构才使用count×平均phase作为权重。

CSV误差是删独立seed jackknife1SE。峰曲线每副本重选max；同q分解保持选定q不变。保留样本不是独立有效样本数，gzip聚合矩不支持重新改变cut、链内重分块或IAT计算。那些操作需要完整HPC blocks/原始轨迹。

源码生产版本不同于当前quick1010；每链manifest、run/source信息保留在gzip JSON中。绝对路径仅为出处，不是离线复现依赖。没有上传HPC认证信息或大型完整逐sweep日志。

复现仅需Python3、NumPy、Matplotlib，从仓库根运行README中的四个命令。各原始README/metadata保存采集窗口与源文件SHA256；根`SHA256SUMS`校验本次导出的全部文件。
