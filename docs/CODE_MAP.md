# 代码地图

2026-10-23：完整能量估计器见[能量定义与操作](TOTAL_ENERGY.md)和
[ADR-029](decisions/ADR-029-total-energy-at-phonon-link-midpoints.md)。
`total_energy.h`实现强迫Mehler链路中点插入，`phonon_projector.h::energyLinkGreens`
重建相邻中心两切面，`PhononSampler::totalEnergy`复用同构型一次参考恢复。
CLI追加六个能量列并记录位置/源/零点能约定；旧`majorana_energy`定义和列位置保留。
`total_energy_validate.cpp`及correctness入口`--group energy`提供固定输入独立谱/Fock验证，
不运行MC。离线工具识别新分项，并把历史可重建势能明确标为中心切面。

2026-10-23：原始块离线最终分析、统一任务总表与统计诊断见
[统一工具说明](CALCULATION_TOOLS.md)和[ADR-028](decisions/ADR-028-offline-analysis-and-task-registry.md)。
`scripts/qmc.py`为入口；`scripts/qmc_tools/blocks.py`冻结完整记录及身份核验，
`statistics.py`共享复相位联合删组，`observables.py`按列名恢复已有派生量，
`diagnostics.py`提供块尺度诊断，`registry.py`及`remote_*.py`只读索引/取数。
`scripts/tests/test_qmc_tools.py`包含原C++结果和失败最终汇总恢复的回归。
生产采样器与PRODUCTION_SOURCE不变，无新计算/取消/减seed。

2026-10-23：信号保存/续算入口见[操作说明](CHECKPOINT_RESUME.md)及
[ADR-027](decisions/ADR-027-signal-checkpoint-resume.md)。
`kitaev_ssh/mc_checkpoint.h`负责异步安全信号标记、完整sweep边界检查点、
校验/原子发布及恢复选项锁定；`checkpoint_io.h`提供序列化辅助。
`phonon_mc.h`、`phase_statistics.h`保存RNG/normal缓存、随机排列、阶段/cadence、
全部计数及未满原始统计块；投影缓存重建并通过原严格数值门槛。
`scripts/checkpoint_job.py`核对PBS节点/PID身份，只向指定MC进程发信号并等保存回执。
既有旧二进制作业不变；续算原始块包含祖先记录，禁止再相加。

2026-10-22 隔离实验入口（未合入活动源码/生产默认）：
[`audit/1022/gamma-next-kernels/`](../audit/1022/gamma-next-kernels/README_CN.txt)
保存 `phonon_collective_block.h`（Γ±/Γ–q基）、独立collective sampler/CLI候选、
`path_snapshot.h`（完整路径与目标身份快照）、`replica_exchange.h`（四交叉权重）、
三个固定输入验证器和真实L9固定候选重放工具。准确接口与范围见该报告/USAGE。


2026-10-20 后续运行入口：`audit/1020/precision-max-budget-20022-workq/`
保存独立32倍预算计划、PBS提交/启动核验、`run_chain.py`与按需
`check_live.py`；后者只给完整落盘块的共同长度诊断，不停止任务、不作正式PASS。
冻结历史trial控制及提前停止统计边界见该目录README.txt；生产接口未改。

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
| `kitaev_ssh/phonon_gamma.h` | 默认关闭的三键独立Γ归一化实模式开放pCN路径与提议抽样；Γ与非零q池独立，见ADR-024 |
| `kitaev_ssh/gamma_validate.cpp` | 复用q模式测试并独立检查Γ g0协方差、16维物理trial Fock MH/phase/缓存/恢复和关闭兼容 |
| `kitaev_ssh/physical_matter_trial.h` | 默认物理 matter 奇偶选择、规范修饰 Bloch 最低模式阻塞、trial身份输出；CLI `--matter-trial physical|unrestricted`，见ADR-023 |
| `kitaev_ssh/physical_matter_trial_validate.cpp` | 独立 Clifford 重排、四Wilson/剪切/规范/平移与最低相容能量、16维 Fock 有符号路径/Green 对照；correctness `--group physical` |
| `kitaev_ssh/phonon_path.h` | 开放虚时间路径、自由声子高斯分布和声子作用量 |
| `kitaev_ssh/phonon_q_mode.h` | 固定允许 q 的实正交模式、实际自由均值中心化、集体 pCN 候选及提议密度；见 ADR-014 |
| `kitaev_ssh/q_mode_validate.cpp` | 自由 Gaussian 详细平衡、Fock/延迟 phase/Green、集体拒绝/恢复事务与默认关闭随机流回归 |
| `kitaev_ssh/phonon_projector.h` | 给定声子路径构造费米子权重和测量量；延迟合并提议只保存 weight_order，canonical factors 在测量时构造，完整 uncached 恢复路径保留（ADR-016） |
| `inc/pure_projector_endpoint.h`、`kitaev_ssh/center_green_validate.cpp` | 中心 endpoint 的无主元 LAPACK thin QR、严格尺度/正交性认证与完整 Green 对照；CLI `--center-green-fast 1`，uncached 参考保留原实现，见 ADR-021 |
| `kitaev_ssh/phonon_mc.h` | 声子路径提议、接受/拒绝与链状态管理；`sample` 可选 sweep 完成回调提供分阶段累计 MoveCounts，不额外测量或抽样 |
| `kitaev_ssh/coupled_mc.cpp` | 命令行生产程序与结果输出；`--acceptance-every N` 输出并刷新运行中接受率 CSV（默认 0 关闭） |
| `kitaev_ssh/numerical_policy.h` | 已验收数值优化的固定生产配置 `production-v1`；单一 reference 诊断配置，采样/物理参数保持独立，见 ADR-022 |
| `kitaev_ssh/production_profile.cpp` | 复用完整生产 CLI 的 profiler，支持实际尺寸、预热、q 更新与全部测量；额外保存 `PREFIX_profile.json`。计时为嵌套 inclusive 时间，不可直接相加；普通生产程序无这些计时开销 |
| `kitaev_ssh/vbs_*.h`、`fermion_observables.h` | VBS 与费米子关联测量 |
| `kitaev_ssh/vbs_susceptibility.h` | 固定 trial 端点的有限投影 VBS 静态场响应原始矩；`--vbs-susceptibility 1`，仅允许动量网格 |
| `kitaev_ssh/vbs_susceptibility_validate.cpp` | Mehler 外场 score、自由声子解析响应和 Fourier/实源归一化确定性验证 |
| `kitaev_ssh/vbs_amplitude.h` | 零场逐构型 A/A2/A4 及可选先模后积分 IA/AIA；对称梯形离散模长源，CLI vbs-amplitude/vbs-amplitude-response，见 ADR-019 |
| `kitaev_ssh/vbs_amplitude_validate.cpp` | 独立 DFT/剪切、模长积分、复权重离散源导数、自由 Rayleigh/半正态矩；CMake/nodeps 目标 kitaev_ssh_vbs_amplitude_validate |
| `kitaev_ssh/analyze_vbs_amplitude.py` | completed 原始块、固定非零 q 的 VA/chiA/B4、每链联合 JK/bootstrap、局部 Binder 斜率窗口与二倍块长检查；JSON/CSV/PNG |
| `kitaev_ssh/test_vbs_amplitude.py`、`scripts/validate_vbs_amplitude.py` | 离线统计/拟合/积分检查与原生新旧开关兼容性、固定场矩对照及 L6 smoke |
| `kitaev_ssh/static_vbs_response.py` | 固定共同 Xbar/Jeff 的静态三键混合 Hessian；全部非零 q、同 q Bloch/实空间对照、实际三档能量差分和多档 BZ 积分；输出 JSON/CSV、候选测量 q-file 与可选 PNG，见 ADR-018 |
| `kitaev_ssh/static_limit_response.py` | 严格静态有限环面全键多初值松弛，默认保留 trial matter 奇偶；解析全键 Hessian、单畴三键全 q 源响应及简并平移畴范围；JSON/CSV/NPY/NPZ/PNG，见 ADR-020 与 STRICT_STATIC_README_CN.txt |
| `kitaev_ssh/test_static_vbs_response.py` | 静态能量归一化、k*I、全部模式正交、混合差分、q/-q/自共轭/非网格与零模边界、CLI 和实际积分记录验证 |
| `inc/pure_projector_*.h` | Majorana 投影、Pfaffian 权重、Green 函数及稳定化 |
| `inc/pure_projector_factor.h`、`inc/pure_projector_weight.h` | 不可变传播因子与权重分片；新生成的已认证 dense 整步矩阵可移动到 owner，lvalue 和 generic 构造语义保持 |
| `inc/pure_projector_green.h`、`kitaev_ssh/green_refinement_validate.cpp` | Green 求解失败后的 LU/扩展精度恢复；合法 Gaussian 子空间的独立 MP50 参考和失败门槛回归，见 ADR-010 |
| `inc/pure_projector_svd.h`、`kitaev_ssh/green_precision_benchmark.cpp` | 生产 Green 的仅奇异值分解及 LU 求解、局部精度/耗时基准；按需恢复策略见 ADR-011 |
| `inc/pure_projector_fast.h`、`inc/pure_projector_green.h` | signed LAPACK QR/SVD与基于实测rcond的生产求解容差；`pureProjectorSolveTolerance`及values-only SVD边界参数见 ADR-013 |
| `inc/pure_projector_fast.h::pureProjectorStableReferenceWeight` | signed partial/final Green在Pfaffian前保留原结构门槛；ADR-025护栏由ADR-026的结构维护与增广定号恢复扩充 |
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

2026-10-17：g=0 解析极限验证：
`pfqmc_free_phonon_validate --deterministic` 覆盖独立 Mehler 精度矩阵、
端点/内部 conditional/action ratio、所有方向 DFT、响应和逐构型 A4；
`pfqmc_fock_validate --g0` 用独立 16 维 Fock Hamiltonian 对照同扇区
自由基态的 phase、matter norm、中心 Green 和能量。
生产 CLI 可选 `--phonon-diagnostics 1`（默认关闭）追加各节点跨键平均的
X/X2、中心同键/邻接键编号相关和 matter log-weight 原始矩，沿用共同
相位分母，不改变已有测量或采样。
`scripts/g0_validation_plan.py` 生成固定 seed/参数预算的 CLI 测试计划；
`scripts/analyze_g0_validation.py` 恢复原始块，链内重分块后联合相位比
jackknife，核对有限 tau 的解析极限、Binder 和精度/误差平台门槛。
本轮执行记录及运行前资源修订见 `audit/1017/g0-analytic-validation/`。

2026-10-17：有限 g 独立验证接口：
`kitaev_ssh/bond_diagnostics.h` 提供验证专用中心各键 Ki=i*cA*cB 和 XKi，
生产 CLI 用 `--bond-diagnostics 0|1` 控制相位加权原始矩（2026-10-22起默认1，显式0关闭）；
Ki 不包含 u，majorana_energy 仍仅为 matter energy。
`kitaev_ssh/interacting_fock_validate.cpp` 保存固定声子路径，独立全 Fock
传播/中心插入及含 proposal density 的 MH 比较，并覆盖接受/拒绝与恢复。
`validation/interacting_ed/` 包含独立全部12量子声子的总占据截断 ED、
有限 tau 与对称离散传播、固定 trial 源响应、配置生成和生产块联合
jackknife 比较。入口和复现说明见该目录 README.txt。

2026-10-19：pCN 自相关、固定物理网格效率与投影收敛验证：
生产 CLI 的 `--validation-series 1` 可选输出逐测量未加权复数原始矩和共同
相位，`--path-series-bond B` 可选输出该键的全部时间节点，均默认关闭。
两个开关不抽取随机数；HPC 上开关前后原有 blocks 逐字节相同。
`validation/pcn_grid_projection/` 复用生产 CLI 和独立 full12-phonon ED：
`generate.py` 固定试跑及纯pCN计划；`plan_formal.py` 只依试跑预算一次，
正式计划哈希后才提交；`run.py`/`run_formal.py` 保存命令、seed、退出码、
线程及二进制哈希；`analyze_pcn.py` 核对完整 Mehler 噪声协方差、端点、
线性量和平方量的 ACF/IAT；`analyze_ed.py` 分列截断、dt 与投影偏差；
`report.py` 以共同相位联合分块 jackknife 汇总 Binder、逐链 ACF、均值、
精度和同误差时间估计。`maybe_submit_report.py` 待全部正式链终止后提交
汇总，包括失败标记，避免已完成 PBS 数组的依赖挂起。物理及默认采样
政策未变；本验证明确 qoff。协议见目录 README.txt，证据在
`audit/1019/pcn-grid-projection-validation/`。

2026-10-21：`audit/1021/precision-progress-review-89640-89643/delivery/`
新增冻结逐seed原始块归档及离线`reproduce/analyze_frozen.py`；重现92量、
联合和逐seed完整blocking，并独立输出数值相容/统计精度/有限投影状态。
`REPRODUCE.txt`给出实际命令、数据列和中途诊断边界；不改变生产接口。

## 2026-10-21 数值修复接口（ADR-026）

- `inc/pure_projector_canonical.h`：受独立舍入界限制的 Gaussian 框架微小维护与行列式幅相补偿，附 process-lifetime 统计。
- `inc/pure_projector_augmented_pf.h`：病态中间定号的等价增广 Pfaffian、重构/逆残差和分支预算；原矩阵预算 F 与通用 skew 目标界 2F 分开解释。
- `inc/pure_projector_green.h`、`inc/pure_projector_endpoint.h`：罕见长双精度框架恢复，以及对同一已保存双精度 Green 元素重新累加残差；不投影/裁剪 Green，保持原数值门槛。
- `kitaev_ssh/coupled_mc.cpp`：输出 `canonical_frame_roundoff`、`green_precision_recovery`、`augmented_phase_recovery` 的 process-lifetime 统计，含初始化/恢复，不作为热化或接受率；工作区既有 opt-in CLI 保留。
