# 代码地图

生产可执行文件由 `kitaev_ssh/coupled_mc.cpp` 构成，依赖同目录模型/采样代码以及 `inc/` 中的 Majorana Pfaffian 投影实现。

系统正确性检查入口为 `scripts/check_pfqmc.sh` 或
`scripts/run_pfqmc_correctness.py`。新增 `kitaev_ssh/pfqmc_*_validate.cpp`
分别覆盖自由声子 Gaussian/统计参考、小系统独立 Fock 参考、缓存及
Monte Carlo 事务；物理约定、测试范围和运行方式见
[PFQMC_CORRECTNESS.md](PFQMC_CORRECTNESS.md)。

| 路径 | 职责 |
|---|---|
| `kitaev_ssh/honeycomb_lattice.h` | 蜂窝晶格、键方向、边界条件与规范扇区几何 |
| `kitaev_ssh/kitaev_ssh_model.h` | SSH–Kitaev 键耦合与 Majorana 生成元 |
| `kitaev_ssh/phonon_path.h` | 开放虚时间路径、自由声子高斯分布和声子作用量 |
| `kitaev_ssh/phonon_q_mode.h` | 固定允许 q 的实正交模式、实际自由均值中心化、集体 pCN 候选及提议密度；见 ADR-014 |
| `kitaev_ssh/q_mode_validate.cpp` | 自由 Gaussian 详细平衡、Fock/延迟 phase/Green、集体拒绝/恢复事务与默认关闭随机流回归 |
| `kitaev_ssh/phonon_projector.h` | 给定声子路径构造费米子权重和测量量；延迟合并提议只保存 weight_order，canonical factors 在测量时构造，完整 uncached 恢复路径保留（ADR-016） |
| `kitaev_ssh/phonon_mc.h` | 声子路径提议、接受/拒绝与链状态管理；`sample` 可选 sweep 完成回调提供分阶段累计 MoveCounts，不额外测量或抽样 |
| `kitaev_ssh/coupled_mc.cpp` | 命令行生产程序与结果输出；`--acceptance-every N` 输出并刷新运行中接受率 CSV（默认 0 关闭） |
| `kitaev_ssh/production_profile.cpp` | 复用完整生产 CLI 的 profiler，支持实际尺寸、预热、q 更新与全部测量；额外保存 `PREFIX_profile.json`。计时为嵌套 inclusive 时间，不可直接相加；普通生产程序无这些计时开销 |
| `kitaev_ssh/vbs_*.h`、`fermion_observables.h` | VBS 与费米子关联测量 |
| `kitaev_ssh/vbs_susceptibility.h` | 固定 trial 端点的有限投影 VBS 静态场响应原始矩；`--vbs-susceptibility 1`，仅允许动量网格 |
| `kitaev_ssh/vbs_susceptibility_validate.cpp` | Mehler 外场 score、自由声子解析响应和 Fourier/实源归一化确定性验证 |
| `kitaev_ssh/vbs_amplitude.h` | 零场逐构型 A/A2/A4 及可选先模后积分 IA/AIA；对称梯形离散模长源，CLI vbs-amplitude/vbs-amplitude-response，见 ADR-019 |
| `kitaev_ssh/vbs_amplitude_validate.cpp` | 独立 DFT/剪切、模长积分、复权重离散源导数、自由 Rayleigh/半正态矩；CMake/nodeps 目标 kitaev_ssh_vbs_amplitude_validate |
| `kitaev_ssh/analyze_vbs_amplitude.py` | completed 原始块、固定非零 q 的 VA/chiA/B4、每链联合 JK/bootstrap、局部 Binder 斜率窗口与二倍块长检查；JSON/CSV/PNG |
| `kitaev_ssh/test_vbs_amplitude.py`、`scripts/validate_vbs_amplitude.py` | 离线统计/拟合/积分检查与原生新旧开关兼容性、固定场矩对照及 L6 smoke |
| `kitaev_ssh/static_vbs_response.py` | 固定共同 Xbar/Jeff 的静态三键混合 Hessian；全部非零 q、同 q Bloch/实空间对照、实际三档能量差分和多档 BZ 积分；输出 JSON/CSV、候选测量 q-file 与可选 PNG，见 ADR-018 |
| `kitaev_ssh/test_static_vbs_response.py` | 静态能量归一化、k*I、全部模式正交、混合差分、q/-q/自共轭/非网格与零模边界、CLI 和实际积分记录验证 |
| `inc/pure_projector_*.h` | Majorana 投影、Pfaffian 权重、Green 函数及稳定化 |
| `inc/pure_projector_factor.h`、`inc/pure_projector_weight.h` | 不可变传播因子与权重分片；新生成的已认证 dense 整步矩阵可移动到 owner，lvalue 和 generic 构造语义保持 |
| `inc/pure_projector_green.h`、`kitaev_ssh/green_refinement_validate.cpp` | Green 求解失败后的 LU/扩展精度恢复；合法 Gaussian 子空间的独立 MP50 参考和失败门槛回归，见 ADR-010 |
| `inc/pure_projector_svd.h`、`kitaev_ssh/green_precision_benchmark.cpp` | 生产 Green 的仅奇异值分解及 LU 求解、局部精度/耗时基准；按需恢复策略见 ADR-011 |
| `inc/pure_projector_fast.h`、`inc/pure_projector_green.h` | signed LAPACK QR/SVD与基于实测rcond的生产求解容差；`pureProjectorSolveTolerance`及values-only SVD边界参数见 ADR-013 |
| `kitaev_ssh/*validate.cpp` | 固定输入的数值参考和回归目标 |
| `reports/physics-check-20261005/validation/` | 本次新增的独立解析/多体交叉检查 |

主要可调参数和定义见 [`PHYSICS_CONVENTIONS.md`](PHYSICS_CONVENTIONS.md) 与 [`OBSERVABLES.md`](OBSERVABLES.md)。当前声子整条路径更新及其低频采样问题见 [`物理与问题说明`](../reports/physics-check-20261005/PHYSICS_AND_ISSUES_CN.md)。

2026-10-08：`PhononSampler` 增加可选第五参数 `QModeOptions` 和 `qCounts()`，
总 `counts()` 包括原提议与集体提议。CLI `--q-mode-every/n1/n2/beta` 默认关闭，
`--q-diagnostic-every` 输出独立原始 phase/固定 q 物理量/CPU 时间；
`_q_acceptance.csv` 与 JSON `move_counts` 分别记录原提议和集体提议。
算法及统计边界见 [ADR-014](decisions/ADR-014-fixed-q-collective-update.md)。

2026-10-10 quick1010 整合固定参考整键 pCN 与固定动量池，见
[ADR-017](decisions/ADR-017-guided-bond-and-q-pool.md)。
`phonon_guided_proposal.h` 缓存参考均值、精确 Mehler 源权重并生成整键候选；
`PhononSampler` 新增第六参数 `GuidedBondOptions`、第七参数固定 q 池和
`qPoolCounts()`，旧调用保持兼容。`phonon_q_mode.h` 支持多键正交模式及
`readQModePool()` 的 Gamma/共轭重复检查。
`coupled_mc.cpp` 新增 `--proposal bond-guided --bond-guide --bond-beta`、
手工参考力与 `--q-mode-pool --q-mode-bonds`，记录实际参考力和池内累计计数。
`guided_proposal_validate.cpp` 对应 CMake 和免 CMake 构建目标
`kitaev_ssh_guided_validate`。Hamiltonian、物理声子权重、trial、相位重加权、
原测量与 Pfaffian/Green 内核保持原实现。生产 profiler 复用同一 CLI。

VBS 静态响应的场、归一化和复相位连接减法见 [ADR-009](decisions/ADR-009-vbs-static-response.md)。
准备包 `audit/production-L9-omega0p3-vbs-response-20261005/` 包含参数表、
PBS arrays、目标主机源码构建脚本和 `summarize_chi.py` 联合连接/重分块汇总。

`kitaev_ssh/deferred_phase_validate.cpp` 还包含满秩、小归一化行列式的独立Majorana旋转振幅回归；`audit/tau12-numerical-repair-20261005/replay.cpp` 从保存路径重建合并/未合并权重和完整中心测量，审计材料与正式重交记录在同一目录。

2026-10-06 有限 VBS 场：`kitaev_ssh/vbs_field.h` 定义固定允许 q 的实 z 键源、
原始纵向 score 矩与逐构型复四阶矩；`phonon_path.h` 新增可选物理 forces、
强迫 Mehler 权重/高斯均值/version3 snapshot；`phonon_mc.h` 的 conditional
和整键 proposal 使用相同强迫高斯。`coupled_mc.cpp` 的 `--vbs-field-*` 接口、
CSV Field_* 和 JSON vbs_field/blocking 见 ADR-012。
`vbs_field_validate.cpp` 为对应数学与详细平衡验证，CMake 和 nodeps 均可构建。
