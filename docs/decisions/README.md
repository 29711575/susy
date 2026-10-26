# 决策记录索引

这些记录保存“为什么这样做、适用条件、不能推出什么、应由什么证据保护”。公式以各专门规范文件为准，ADR不复制长篇推导。

| 记录 | 状态 | 内容 |
|---|---|---|
| [ADR-001](ADR-001-vbs-wavevector.md) | recorded | Dirac动量差作为VBS候选，而非固定所有参数峰位 |
| [ADR-002](ADR-002-background.md) | recorded | connected与projected分开，零模保留 |
| [ADR-003](ADR-003-offgrid.md) | recorded | 非网格平移平均不改变Hamiltonian |
| [ADR-004](ADR-004-driven-protocol.md) | proposed | SSH声子驱动接入、双支时间和最低验证 |
| [ADR-005](ADR-005-matrix-statistics.md) | recorded | 矩阵先平均、联合统计和采样证据边界 |
| [ADR-006](ADR-006-dimensionless-observables.md) | definition-recorded / extension-proposed | 当前关联比、ξ/L命名与Binder扩展边界 |
| [ADR-007](ADR-007-bipartite-spectrum.md) | implemented / validated | 二分格完整 SVD 矩阵函数、持久 overlap 工作区及参考回退 |
| [ADR-008](ADR-008-real-magnitude.md) | implemented / local-validated; HPC pending | 等价实数 magnitude、精确回退、小尺寸分派和按需中心 Green |
| [ADR-009](ADR-009-vbs-static-response.md) | accepted; validation recorded in STATUS | 精确 Mehler 场插入、固定 trial 的 VBS 有限投影响应、联合连接统计 |
| [ADR-010](ADR-010-green-solve-recovery.md) | implemented; validation recorded in STATUS | 病态满秩重叠的 LU 精化与一致扩展精度 Green 求解，保留原认证门槛 |
| [ADR-011](ADR-011-selective-green-precision.md) | fast-1006 implemented / local validated | 生产 Green 先用已有相对检查、values-only SVD + LU，按需恢复扩展精度 |
| [ADR-018](ADR-018-fixed-background-static-curvature.md) | quick1010 implemented / local validated | 共同均匀背景的静态软模筛选、k*I、独立差分与同动量几何对照；区别于有限投影响应 |
| [ADR-019](ADR-019-vbs-amplitude-peaks.md) | quick1010 implemented / focused validation | 更正 note 的模长矩、梯形模长源响应、联合误差和固定模式 Binder 斜率峰；保真估计器留作下一阶段 |
| [ADR-021](ADR-021-center-green-acceleration.md) | quick1010 implemented / native MKL validated | 中心 endpoint LAPACK thin QR、values-only SVD/LU，保留严格数值门槛和完整参考恢复 |
| [ADR-023](ADR-023-physical-matter-trial.md) | implemented; execution evidence in audit/1019 | 默认生产 trial 的物理奇偶选择、简并 Bloch 模式阻塞与历史控制兼容 |
| [ADR-024](ADR-024-opt-in-gamma-path-updates.md) | implemented; optional, default off | 独立 xyz Γ 开放路径 pCN、额外 q 池与混合诊断 |
| [ADR-025](ADR-025-signed-prefix-structure-guard.md) | implemented; extended by ADR-026 | signed prefix 在 Pfaffian 前强制 Green 结构认证 |
| [ADR-026](ADR-026-gaussian-frame-phase-repair.md) | implemented; scoped deterministic/reference validation | Gaussian 框架维护、病态中间定号恢复及相位误差预算 |
| [ADR-027](ADR-027-signal-checkpoint-resume.md) | implemented; local regressions and target build verified | 完整sweep边界信号检查点、RNG/统计恢复、新前缀及PBS实际PID核验 |
| [ADR-028](ADR-028-offline-analysis-and-task-registry.md) | implemented; historical C++ parity and read-only live validation | 原始块独立分析、端口限定任务索引、统一复相位统计与块尺度诊断 |
| [ADR-029](ADR-029-total-energy-at-phonon-link-midpoints.md) | implemented; independent fixed-input and target build verified | 相邻中心声子中点上的完整能量/分项、严格切面重建、联合离线统计 |

`recorded`是对已用约定的记录，不是新的执行授权；`proposed`未定稿，不能当作已实现或用户已批准。已有ADR被替代时标明superseded和替代编号，保留历史原因，不直接删掉旧论证。

新ADR最少写：问题、状态/来源、选择、条件/限制、影响的代码、已有或所需验证。若当前仓库有编号冲突，合并时重编号并更新链接，不覆盖既有决策。
