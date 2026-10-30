# SSH–Kitaev 模块入口

本目录包含耦合声子投影 Monte Carlo、自由模型参考、观测量和独立验证。
当前生产入口为 [coupled_mc.cpp](coupled_mc.cpp)，构建方法见[仓库 README](../README.md)。

| 模块 | 入口 |
|---|---|
| 声子开放路径与更新 | [phonon_path.h](phonon_path.h)、[phonon_mc.h](phonon_mc.h) |
| 投影权重、trial 与缓存 | [phonon_projector.h](phonon_projector.h) |
| guided 整键及非零 q 更新 | [phonon_guided_proposal.h](phonon_guided_proposal.h)、[phonon_q_mode.h](phonon_q_mode.h) |
| 等时全动量关联 | [vbs_observables.h](vbs_observables.h) |
| 线性 connected 极化率 | [vbs_susceptibility.h](vbs_susceptibility.h) |
| 模长矩与响应 | [vbs_amplitude.h](vbs_amplitude.h)、[用法](../VBS_AMPLITUDE_README_CN.txt) |
| 自由模型确定性参考 | [free_projector.cpp](free_projector.cpp)、[validate.cpp](validate.cpp) |
| 动量选择与候选诊断 | [说明](VBS_SELECTED_MOMENTA.md)、[vbs_candidates.h](vbs_candidates.h) |

[物理约定](../docs/PHYSICS_CONVENTIONS.md)、[数值约束](../docs/NUMERICAL_INVARIANTS.md)
和[观测量定义](../docs/OBSERVABLES.md)分别解释模型、认证及统计口径。
离线处理新任务推荐 [scripts/qmc.py](../scripts/qmc.py) 和
[任务工具说明](../docs/CALCULATION_TOOLS.md)；现成测量可按[测量索引](../measurements/README.md)复现。

早期 P0–P1 自由基线说明已移至[历史归档](../docs/history/KITAEV_FREE_BASELINE_README_CN.md)。
旧 PBS 文件与无当前调用者的计划已收进[历史脚本压缩包](../docs/history/legacy-hpc-scripts-before-20261030.zip)，
其中硬编码路径与参数不是新任务模板。

仍被旧工具或回归使用的 `fermion_gscan_plan.json`、`phase_tau_scan_plan.json`、
`t32_replication_plan.json`、`t32_replication8_plan.json` 和 `p5_plan.json` 保留在原位。
它们属于历史计划或测试样本，不选择当前生产参数，也不授权运行旧计算。
