# 2026-10-30：数值核心与近期结果补充

本次基于 GitHub `5f4c357`，补充 10 月 26 日导出后的已验证实现与近期分析。
生产核心逐文件匹配 2026-10-28 冻结归档，SHA256：
`9b7c24398cf4d650480bd89396ab62898efe394859dfe4b39f64589e4b04fdcf`。
根目录继续采用原仓库布局；历史测量与验证文件保持其原始日期。

## 代码

| 更新 | 入口及边界 |
|---|---|
| QR/Green 工作区 | [ADR-030](../decisions/ADR-030-qr-green-workspace-reuse.md)：链私有缓冲，保持原数值认证与参考恢复。 |
| 路径因子复用 | [ADR-031](../decisions/ADR-031-path-factor-storage-and-real-companion.md)：不可变发布与唯一所有权约束，避免重复分配。 |
| 测量因子复用 | [ADR-032](../decisions/ADR-032-measurement-factor-reuse.md)：密封快照复用完整谱信息，保持原测量定义。 |
| 中心传播及谱安全 | [ADR-033](../decisions/ADR-033-endpoint-dense-propagation.md)：直接读取稠密因子；任一谱模的相关指数非有限时，安全判据失败。 |
| 独立验证 | 新增工作区、路径/测量因子、Fock 参考、endpoint 及谱安全验证器，保留各自检查范围。 |

本次没有改变 Hamiltonian、trial 默认、物理测量定义、MH/RNG 或原认证门槛。
历史配对计时未显示跨输入一致的整链或中心重建提速，因此不作整体性能承诺。
原记录与本次导出检查分开保存，见
[验证记录](../../evidence/source-sync-20261030/README.md)。

部分新任务使用的多尺度 Gamma±/z 是另一份已验证冻结实现；该实验 CLI
**没有并入本次普通核心**。本次源码仍支持已有可选 `--gamma-every`，
不能直接使用实验任务中的 `--collective-beta-scales` 等选项。

## 结果与复现

新增 [2026-10-30 测量摘要](../../measurements/20261030/README.md)，包括：

- L6、ω=0.1：原 3 seed、每链 cut4000，全非零 q 的 m² 峰，以及只取裸 Dirac
  点差最近 q={15,20,22,27} 的 m² 峰；参考斜率 2/3 明确标为非拟合。
- L9、ω=1：λ=6/7/8，τ=2–8，cut5000；τ10/12 未达到截断，数值失败链排除。
- L6、λ=2：ω=0.1/0.2/0.3/0.5，cut4000，显式区分 dt 及采样控制。

每组附联合复相位原始矩 gzip JSON、CSV 和不依赖本地 audit/HPC 的作图代码。
它们是指定时间的测量快照，不是实时任务清单；聚合矩不能重建逐 sweep
自相关，也不能仅凭 seed jackknife 和指定 cut 证明热化或遍历。

L6 裸 Dirac 点差最近网格点按真实倒空间距离计算，不能逐坐标四舍五入。
[物理约定](../PHYSICS_CONVENTIONS.md#phy-05-键连接与-dirac-候选)已补充该区别。

本次仓库更新不提交或操作 HPC 作业；既有作业保留其原程序、配置和检查点。
