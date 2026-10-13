# ADR-021：保留严格认证的中心 Green 重建加速

日期：2026-10-10。状态：quick1010 已实现，原生 MKL 确定性验收通过。
来源：用户要求评估并开启尚未接入的优化；允许源码审查和集中数值对照验收，
不要求额外 QMC 长链。

## 决定

新增 `--center-green-fast 1` 和 `PhononWeightOptions::center_green_fast`，
只选择中心 endpoint 重建的数值实现。占据维数 >=32 的 endpoint 传播
使用 LAPACK 无主元 thin QR（zgeqrf/zungqr），覆盖 R 前保存全部对角
绝对值的 log 尺度；不截断模式。endpoint overlap 使用仅奇异值 SVD
和 pivoted LU 求解。小尺寸 QR 保留 Eigen。

矩阵乘积顺序、左端 adjoint、物理 Green 的最终 adjoint、signed Pfaffian
相位、Hamiltonian、trial、声子作用量、提议分布和完整 MH 校正保持原定义。
rank/minimum-rcond 仍为 1e-12，solve/Green/QR 正交性门槛仍为 1e-10。
不启用 oblique Green 放宽，也不启用基于 rcond 放宽 solve 的选项。
已有精化和扩展精度恢复保持；`uncachedReference` 刻意不传新开关，
同一候选只重试一次完整原参考，持续失败仍终止。

`reuse_rebuild_workspace` 随快速选项设置，但当前 endpoint 实现不消费
该 stack 选项；本次收益不能归因于 stack 工作区复用。

新 CLI 省略时仍为 0，保留旧调用的算术行为；未来普通生产提交显式写
`--center-green-fast 1`。有意对照显式写 0。现有已提交任务的冻结程序
不支持此新参数，不热切换、不改写现有作业目录。

## 验收及范围

新 validator 对 n=4/32/36/81 的复矩阵 QR 子空间及 log 尺度、秩边界、
Green overlap 拒绝边界作对照。L6/M12、L6/M100/omega.3/lambda2.4、
L9/M200/omega1/lambda3.7、L9/M240/omega.3/lambda3 的完整重建、
eager/deferred、局部接受/拒绝及完整参考比较均通过，最大缩放差
9.28e-15，验收阈值 2e-9 未调整。

三组固定 seed、guided bare/beta1、q 池/xyz/q-beta.5 的旧原生程序、
新关闭、新开启短程对照通过：接受/拒绝计数及最终路径完全一致。
新关闭与旧程序原始块逐数完全相同；新开启原始块最大缩放差 3.93e-15。
恢复/终止、计数、RNG、事务和延迟测量注入检查通过。故障注入使用
原小晶格，检验恢复控制流；LAPACK 分支的数值边界由新 validator 检验。

交错四次 endpoint 计时中 L6/M100、L9/M200、L9/M240 分别降低约
7.2%、8.4%、6.1%。这是人工平滑路径、局部原生 MKL 的中心重建段计时，
不是整链或 HPC 效率保证，也不认证热化/平衡。

当前整键/q 更新不接入共享库的 HS 局部低秩 walker。单 SSH 键在
测试路径上改变生成矩阵的秩为2，改变矩阵指数的数值秩为8（阈值1e-12）；
不能把生成矩阵的秩当作传播矩阵差的秩。本次不启用该独立引擎。
`PFQMC_UDT_DISABLE_RANK_LOSS_GUARD` 会删除秩损失、指数范围、有限性
和部分/最终 Q 正交性保护，不是本次验收的生产配置；继续不定义它。

代码：`inc/pure_projector_endpoint.h`、`kitaev_ssh/phonon_projector.h`、
`kitaev_ssh/coupled_mc.cpp`、`kitaev_ssh/center_green_validate.cpp`、
`kitaev_ssh/numerical_failure_validate.cpp` 及 CMake/nodeps 目标。
完整证据：[验收报告](../../evidence/recorded-validation/center-green/REPORT.txt)。

2026-10-11 补充适用范围：上文关于QR关闭宏删除保护的描述针对共享库
UDT引擎。当前SSH生产入口不调用该引擎，实际投影QR/Green检查由各自代码
控制，不受该宏影响。只增加宏的原生新构建与已验收程序逐字节一致，
故当前长链通过这个宏获得的算法耗时收益为0%。它不能算作当前生产
遗漏的有效加速开关；见原工作区 audit/1010/qr-guard-cost/REPORT.txt（历史记录，未随本次导出）。

2026-10-11 局部Green适配可行性补充：独立完整支撑稠密原型在L6/L9
固定路径上通过权重比/完整同cut Green对照。L9/M200,M240热cache的
单节点权重kernel约快7.5–8.3倍，整键/xyz q逐节点替换约慢18–24倍。
这不是生产集成或整链效率验收；缓存建立/失效与提议数量需要另计。
现有精确小支撑低秩内核仍不适用，不用数值秩阈值删模式，也不替换当前
非对易矩阵指数的离散方式。当前不新增生产开关，现有guided/q保持。
见原工作区 audit/1011/local-green-feasibility/REPORT.txt（历史记录，未随本次导出）。
