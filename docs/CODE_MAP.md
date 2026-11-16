# QMC 代码阅读顺序

本页按一次计算的执行过程组织，不要求了解开发历史。模型和算法的整体说明见[审阅说明](QMC_REVIEW_GUIDE_CN.md)。下列文件都属于此仓库；历史报告中指向原工作站 `audit/` 或构建目录的路径不一定随仓库提供。

## 1. 从输入到测量的关系

```text
coupled_mc.cpp：读取参数，建立晶格、模型、试探态、路径和统计对象
    ↓
PhononSampler：生成路径提议并决定接受或拒绝
    ↓
CoupledPhononPath：为旧路径和候选路径计算权重，管理传播缓存
    ↓
inc/pure_projector_*.h：费米子传播、Pfaffian 权重、中心 Green 函数
    ↓
coupled_mc.cpp 中的测量回调：计算声子、费米子及响应量
    ↓
PhaseStatistics：累计相位及相位加权测量，保存统计块
    ↓
scripts/qmc.py 或指定数据集的分析脚本：重加权、派生量、误差与图表
```

Green 函数是计算费米子关联所需的二点函数矩阵；其源码约定为 Gij=δij−〈ci cj〉，详见[验证说明](PFQMC_CORRECTNESS.md)。

## 2. 建议依次打开的文件

| 顺序 | 文件与主要对象 | 审阅重点 |
|---|---|---|
| 1 | [coupled_mc.cpp](../kitaev_ssh/coupled_mc.cpp)，`main` | 参数如何传入各模块、测量开关、输出和异常处理；先理解主流程即可 |
| 2 | [honeycomb_lattice.h](../kitaev_ssh/honeycomb_lattice.h) | 原胞和格点索引、三种键、周期和剪切边界、规范变量 |
| 3 | [kitaev_ssh_model.h](../kitaev_ssh/kitaev_ssh_model.h) | 有效键耦合及 Majorana 矩阵的符号和倍数 |
| 4 | [physical_matter_trial.h](../kitaev_ssh/physical_matter_trial.h) | 所需宇称、实际试探态和最低能模式选择 |
| 5 | [phonon_path.h](../kitaev_ssh/phonon_path.h) | 时间节点、端点波函数、声子作用量与条件分布 |
| 6 | [phonon_mc.h](../kitaev_ssh/phonon_mc.h)，`PhononSampler` | `sweep`、`attempt`、`attemptBond`、`attemptGuidedBond`、`attemptQMode`、`attemptGamma`和`sample`；接受、拒绝、异常恢复 |
| 7 | [phonon_projector.h](../kitaev_ssh/phonon_projector.h)，`CoupledPhononPath`、`evaluatePath` | 路径到传播因子的映射、因子顺序、权重绝对值与相位、测量中心及缓存 |
| 8 | [phase_statistics.h](../kitaev_ssh/phase_statistics.h)，`PhaseStatistics`、`summarizeBlocks`、`summarizeDerivedBlocks` | 相位共同分母、统计块、派生量和联合误差 |

读完上述流程，再进入通用投影内核。`inc/` 还包含其他模型用到的通用代码，初次阅读不必逐一展开。

## 3. 按问题深入

| 问题 | 对应文件 |
|---|---|
| 整键路径提议的参考分布和概率修正 | [phonon_guided_proposal.h](../kitaev_ssh/phonon_guided_proposal.h) |
| 指定非零动量或空间均匀模式的更新 | [phonon_q_mode.h](../kitaev_ssh/phonon_q_mode.h)、[phonon_gamma.h](../kitaev_ssh/phonon_gamma.h) |
| Pfaffian 权重与相位 | [pure_projector_weight.h](../inc/pure_projector_weight.h)、[pure_projector_fast.h](../inc/pure_projector_fast.h)、[pure_projector_augmented_pf.h](../inc/pure_projector_augmented_pf.h) |
| 传播稳定化、中心 Green 函数 | [pure_projector_canonical.h](../inc/pure_projector_canonical.h)、[pure_projector_endpoint.h](../inc/pure_projector_endpoint.h)、[pure_projector_green.h](../inc/pure_projector_green.h) |
| 声子 Fourier 结构因子 | [vbs_observables.h](../kitaev_ssh/vbs_observables.h) |
| 静态外场响应和有限场 | [vbs_susceptibility.h](../kitaev_ssh/vbs_susceptibility.h)、[vbs_field.h](../kitaev_ssh/vbs_field.h) |
| 模长及模长响应 | [vbs_amplitude.h](../kitaev_ssh/vbs_amplitude.h) |
| 费米子关联 | [fermion_observables.h](../kitaev_ssh/fermion_observables.h) |
| 总能量 | [total_energy.h](../kitaev_ssh/total_energy.h)、[定义说明](TOTAL_ENERGY.md) |
| 保存与恢复计算 | [mc_checkpoint.h](../kitaev_ssh/mc_checkpoint.h)、[操作说明](CHECKPOINT_RESUME.md) |
| 离线统计 | [scripts/qmc.py](../scripts/qmc.py)、[scripts/qmc_tools/](../scripts/qmc_tools/)、[工具说明](CALCULATION_TOOLS.md) |

稳定化指避免长时间矩阵乘积导致有效数字丢失的数值处理；主要使用 QR 分解、奇异值分解及必要的高精度重新计算。相关误差判据集中在[数值约束](NUMERICAL_INVARIANTS.md)。不能通过放宽误差判据掩盖计算失败。

## 4. 输入、输出和验证入口

主程序输出前缀由 `--output PREFIX`指定。`PREFIX_blocks.csv`保存相位和各测量的原始统计块，`PREFIX_results.json`保存参数、计数和汇总状态，`PREFIX_numerical.jsonl`保存代表性的数值诊断；其他输出随测量或检查点选项启用。CSV 应按列名读取，不能假设所有开关组合的列位置相同。

检查完整调用流程时，可对照 `coupled_mc.cpp` 中建立观测量名称、测量回调、保存块和最终汇总四部分。检查传播和权重时，优先对照独立参考程序，而非只比较两条共享同一内核的代码路径。

[验证说明](PFQMC_CORRECTNESS.md)介绍自由声子、独立多体参考和缓存测试；[报告索引](REPORT_INDEX_CN.md)区分验证程序与实际执行日志。旧版代码地图中的逐日开发记录保留在[原文](history/CODE_MAP_BEFORE_20261116.md)，仅用于追溯。
