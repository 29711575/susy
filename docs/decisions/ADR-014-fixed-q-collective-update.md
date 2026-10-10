# ADR-014：固定允许动量的声子集体更新

日期：2026-10-08。用户批准实现及新旧算法各四 seed 的 HPC 对照；
正式提交状态以本轮 audit 的 submission.json 为准。

原 bond-refresh 每次刷新一条键的完整时间路径，空间上仍逐键移动。
新增可选 z 键固定允许 q 的空间模式更新，保留整键 sweep 后每 N sweep
增加一次集体提议。默认 N=0，旧随机流与数值结果保持。

将 X−freeMean 投影到实正交空间基：一般 q 的 cos/sin 各一条，
自共轭 q 仅一条。允许网格包含剪切几何，禁止开启 Gamma 集体更新。
每条模式的完整时间系数 a 采用 a'=sqrt(1−beta²)a+beta*xi，
xi 的协方差为实际自由声子路径的 exp(−omega*|t−t'|)/(2m omega)。
保留正交补、实际 trial、force、Hamiltonian 和相位重加权。
该提议对自由声子 Gaussian 可逆，接受率仅需 matter log magnitude 差。

候选全部时间节点一次重建，通过一次不可变快照提交或拒绝。
数值失败时保留候选与 u，使用无缓存参考重试一次，仍无法认证则中止。
不能把无法认证的候选当作正常拒绝。

CLI：`--q-mode-every/n1/n2/beta`；`--q-diagnostic-every` 在测量点
输出原始 phase、D、|D|²、已有能量/Xz 及 CPU/墙钟时间，不添加测量或 RNG。
诊断 interval 必须是测量 stride 的倍数。独立 q/base 接受与失败计数另存
`_q_acceptance.csv` 和 JSON `move_counts`；旧 acceptance CSV 仍记录总数。

本轮本地验证：自由 Mehler 与一般 MH 抵消、小 Fock 权重比、
L3 剪切/非对易 honeycomb 的延迟 phase/Green、同候选恢复和拒绝事务通过。
g=0 实际采样均值/方差/时间协方差及默认关闭的 RNG 回归通过。
与之前冻结生产源码同编译器重编译的 local/bond CLI 原始 blocks 逐字节一致，
最终路径、全部物理观测量和接受计数一致；开启诊断也不改变 blocks。
详细记录见 `audit/qmode-pilot-L9-lambda2p4-tau14-20261008/`。

计划对照：L9、omega .3、lambda 2.4、tau14，warmup1000+
measurement10000，各四独立 seed；beta1、每 sweep 一次 q=(5/9,4/9) 更新。
两组使用同一冻结新程序，baseline 关闭 q 更新，逐 sweep 保留原始诊断。
只用20015、每链1CPU、数学库1线程、720h、无显式内存，原任务版本保留。
短正确性检查和接受率不能证明收敛、遍历性或效率优势。

## 2026-10-10 新生产提交的默认配置

用户要求集体更新等适用优化在正确性验证后默认开启。未来普通生产提交
显式设置`--q-mode-every 1 --q-mode-beta 1`及适配几何/目标的非零允许
动量；当前J=(1,1,0.5)/shear0固定q研究，L6索引(3,3)，L9索引(5,4)。
几何或目标改变时重新选择索引，不能启用Gamma或复制其他尺寸的索引。
结合后续已验证计算优化，使用最新冻结生产程序并记录参数和版本。

CLI省略参数的默认N=0仍保留；这是生产提交的默认偏好调整，提交者
不能再依赖该旧CLI默认。基线对照显式N=0，用户明确设置优先，已有
提交和冻结快照保持。正确性通过允许使用，不等于认证长链收敛或普遍
加速；暂定性能结论仍按原统计边界报告。
