# 2026-10-26 主线源码同步

本次从本地 quick1010 同步生产主线源码、验证程序和分析工具，基于 GitHub
`5493696301052bd9011b4110c964977d1811e11d`（2026-10-13）。历史 measurements、
evidence 和交接数据保留原样；本次没有提交、停止或修改 HPC 计算。

| 改动 | 当前行为与入口 |
|---|---|
| 物理 matter trial | 默认 `--matter-trial physical`，按几何/规范选择相容奇偶性；典型 L6/++ 改为最低能 blocked trial，L9/++ 原态不变。`unrestricted` 仅作历史控制。见 [ADR-023](../decisions/ADR-023-physical-matter-trial.md)。 |
| 数值正确性 | Gaussian 框架微小结构维护及幅相补偿、signed Green 结构认证、病态中间步骤的增广 Pfaffian 和高精度残差恢复。原 rank/rcond/残差门槛保留。见 [ADR-025](../decisions/ADR-025-signed-prefix-structure-guard.md)、[ADR-026](../decisions/ADR-026-gaussian-frame-phase-repair.md)。 |
| 可选采样 | 独立 xyz Γ 开放路径 pCN、额外 q 池、初始路径与混合诊断；`--gamma-every 0` 仍为普通 CLI 默认。见 [ADR-024](../decisions/ADR-024-opt-in-gamma-path-updates.md)。 |
| 保存/续算 | 信号与定期检查点保存路径、RNG、进度及原始统计；必须用原可执行文件和新输出前缀续算。见 [操作说明](../CHECKPOINT_RESUME.md)。 |
| 能量/诊断 | 新增五项分量及总能量六列，`--measure-total-energy 0` 默认关闭，旧 matter energy 保留；新增逐键及声子诊断，修复 VBS 汇总切片维度。见 [能量定义](../TOTAL_ENERGY.md)。 |
| 离线工具 | `scripts/qmc.py analyze/diagnose/registry` 提供联合复相位统计、重分块、诊断及只读任务表。见 [工具说明](../CALCULATION_TOOLS.md)。 |
| 独立验证 | 增加 physical/Gamma/energy 固定输入验证、有限耦合 Fock/ED 和 pCN 网格验证工具。导出时实际执行的检查见 [导出记录](../../evidence/source-sync-20261026/README.md)。 |

推荐的新任务测量选项显式写为：

```bash
--fermions 1 --vbs 1 --vbs-mode all --vbs-susceptibility 1 \
--vbs-amplitude 1 --vbs-amplitude-response 1 --measure-energy 1 \
--measure-total-energy 0
```

`production-v1`、guided 整键、固定 q 池、中心 Green 加速、VBS 响应/幅度及静态工具
在 2026-10-13 版本已存在。Γ±、Γ–q、多尺度及副本交换仍在本地隔离实验中，
没有并入本次普通生产入口。当前源码的总能量省略默认是0；历史已验证冻结
二进制仍可能默认跟随 matter energy，新任务始终显式传入0。

`docs/provenance/PRODUCTION_SOURCE.local-20261023.json` 是原本地冻结生产版本记录，
其中路径和构建身份是历史来源，不能用作当前 GitHub checkout 的可执行文件选择器。
当前导出文件身份由根目录 `SOURCE_MANIFEST.json` 和 `SHA256SUMS` 给出；本次导出
没有获得新的 HPC 目标运行时验证。目标生产构建仍须从冻结源在目标运行时构建。

同步的技术文档保留本地历史 audit/build 路径和引用；这些材料不全随仓库导出，
不能把历史报告当作本次重跑。已有作业和数据不能追溯获得新功能；不混合不同
trial/sector 或续算祖先数据。本次源码同步和确定性验证不认证热化、混合、投影收敛或 SUSY。
