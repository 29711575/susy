# quick1010 优化交接

2026-10-10 用户更正 note 后，当前工作转为零场模长峰测量，见 ADR-019。
新 vbs_amplitude.h 与两个可选 CLI 保存 A/A2/A4、IA/AIA；离散梯形源独立
定义，不套用线性 Mehler score。analyze_vbs_amplitude.py 固定非零 q，
联合相位比与每链内 jackknife/bootstrap，输出 VA/chiA/B4、局部 Binder
斜率窗口依赖、二倍块长诊断与候选峰。原生专项、12 Python测试、3组新旧
关闭/开启原始块逐列回归、L6维度与5点短扫描全流程 PASS；尚无物理峰或
临界点结论。新本机入口 build-amplitude-20261010/run_mc.sh；未来生产需
冻结新版并原生构建。既有 HPC 作业保持旧冻结代码，不能从旧块补出新矩。
保真估计器未实现，note 明确列为后续可选；两半 Pfaffian log 权重的混合
导数不能省略。详细用法 VBS_AMPLITUDE_README_CN.txt，证据在
../fast-1006/audit/1010/vbs-amplitude-peaks-corrected-note-20261010/REPORT.txt。
下段静态工具来自被更换的旧 note，保留为独立参考。

2026-10-10 已按 `new/note.docx` 实现/修正固定共同均匀背景的静态三键
软模筛选；docx 没有外部超链接，其旧路径由用户指定的 quick1010 覆盖。
CLI 采用 note 的 J=(1,1,.5)、g=.3、mass1.7/omega.8/k1.088，
实际差分、两档 BZ 参考、零模和同 q 几何对照均保存；13测试与105个非零 q
检查 PASS。三种36原胞几何曲率均正；剪切明显改善 Q 覆盖，并在一个
共有 q 上改善有限 k 偏差，但不能据此确定有限频率量子相变或普遍收益。
入口 `kitaev_ssh/static_vbs_response.py --output-dir DIR --plot`，
详见 ADR-018 和 `../fast-1006/audit/1010/static-vbs-note-20261010/REPORT.txt`。
本次静态工具不改采样 C++、优化开关、已有作业或冻结程序。

2026-10-10 当前代码修改/优化基线已由用户指定为 quick1010。
本目录从用户选定的 fast-1006 核心包解压，并整合 new/ 内的 guided 整键
pCN 和固定 q 池源码；本机原生 MKL 的新提议与旧模式兼容性验证通过。
新整键生产配置按 AGENTS 显式记录 guided bare/beta1，单 q 规定沿用；
动量池保持按目标选用。CLI 省略 proposal 仍为 local。
测量、Hamiltonian、trial 和 Pfaffian/Green 核心保持原实现。
本次不宣称 L9 混合或收敛改善；完整验证范围及源码/二进制出处见
../fast-1006/audit/1010/quick1010-guided-integration/REPORT.txt。
以下 fast-1006 状态为继承历史，不是本次重跑结果。

2026-10-10 用户更新新计算默认配置：集体q更新等适用的优化开关在正确性
验证通过后默认开启，无需再次询问。新提交使用已验收冻结新版，显式写入
`--q-mode-every 1 --q-mode-beta 1`及适配当前几何/目标的非零允许动量；
现有J=(1,1,0.5)/shear0固定q研究，L6用n=(3,3)，L9用n=(5,4)。
combined/deferred phase/prefix+Cayley缓存沿用已验证生产配置。
对照组显式关闭待比较开关；用户明确设置优先；已有作业保持原配置。
这是未来生产提交的默认偏好，CLI省略q参数仍为关闭，提交脚本必须显式
开启。诊断、额外测量及物理外场/驱动参数另按实际任务选择，不一概开启。
本轮只记录偏好，未修改采样源码、重编译或提交任务；详见AGENTS.md。

2026-10-10 本次长链检查：156968/156969的8对均R，首版本测量2281–2803，
完整配对0/8；累计数值失败/恢复/未认证拒绝0。长链验收及配对CPU收益
仍未完成，不能标为通过。此前正确性已验收的半步修改可用于新生产
提交，集体q继续默认开启，提交说明保留长链对照进行中；核心源码与
验收冻结版一致。本次未提交或修改任务。见audit/1010/lazy-half-long-review-20261010/REPORT.txt。

2026-10-10 长链对照已提交并按用户要求停止：20015/workq，156968[]为
q关闭四对，156969[]为q开启四对，共8任务/16链；L6/tau14/M140、
预热1000+测量10000、chi+fermions+energy。每对同seed顺序跑新旧冻结
fast-1006程序，交替先运行版本，自动在结束时比较物理输出与CPU。
提交快照8R；没有等待或收集结果，长链一致性/速度仍待验证。
完整参数、版本和脚本见audit/1010/lazy-half-long-L6-20261010/submission.json。
已有生产任务保持原版本。

2026-10-09 启动、10-10 完成的半步矩阵按需构造已验收：延迟合并提议
只保存接受判据所需的完整weight_order，内部canonical半步在测量或指数
溢出fallback时构造；局部/整键/q更新及完整uncached恢复路径保持安全。
20015构建156965和数组156966/156967共5子任务全部exit0、PASS。
L6状态切换检查、既有624接受/240拒绝缓存事务和失败恢复验证通过；
配对最终路径/计数精确一致，原始blocks逐字节一致、物理JSON/CSV差0。
相对上一轮3b10b149版本，L6/M140/warmup0四次正反配对CPU减少
8.55%–12.16%，中位10.48%；短工程测速，不外推L9/长链收敛收益。
profiling半步构造1430次，原逐节点方式对应61617次，嵌套计时不可相加。
当前冻结源码manifest f4689c4f…、HPC binary2fd634ae…；本机新程序
build-lazy-half-20261009/kitaev_ssh_mc（2a56ce40…），同目录run_mc.sh
设置MKL和单线程。本机仅构建/--help，数值验证均在20015。
证据见audit/1009/lazy-half-L6-20261009/REPORT.txt，设计见ADR-016。
既有生产任务未改版本；新计算使用该冻结快照或其核验后的更新版本。

上一轮2026-10-09代码优化：CLI去掉无用预热phase/Green，
生产整步矩阵直接移动到不可变存储，新增复用完整CLI的production profiler。
20015两批共13个L6测试子任务全部PASS；配对路径/接受计数一致、原始blocks
逐字节一致，缓存624接受/240拒绝事务通过。L6/M140/warmup0四次正反配对
CPU节省3.97%–5.19%，中位4.46%；不外推L9或长链收敛收益。
最终冻结源码manifest c82f83a6…、HPC binary3b10b149…；本机新版程序
build-code-speed-20261009/kitaev_ssh_mc（c5db92c7…），可通过同目录
run_mc.sh设置MKL路径和单线程环境。本机完成构建/--help，数值测试均在20015。
源码快照、diff、全部运行/速度/身份记录见audit/1009/code-speed-L6-20261009/REPORT.txt。
已提交生产任务保留原版本；新计算使用本轮冻结源码或其核验后的新快照。

2026-10-07 用户确认：fast-1006 为以后代码修改、优化和新计算的默认代码基线。
本轮逐文件/逐差异核对当前 fast-1002，未发现其独有、尚未迁入的核心改进。
核对记录在 audit/baseline-selection-20261007/comparison.json；这是源码比较，
本轮未重跑计算。fast-1002 与 fast-0929 保留作历史参照。
新计算从 fast-1006 冻结源码构建，现有已提交计算保留其记录版本。
2026-10-07程序保留在 build-signed-svd-precision-20261007/kitaev_ssh_mc，
对应数值对照与速度证据见 audit/signed-svd-precision-20261007/REPORT.md。

2026-10-06 从当前 fast-1002 完整核心复制开始，复制清单和原始代码快照见 audit/precision-speed-20261006。fast-1002 与 fast-0929 保留。当前未初始化 Git。

目标：减少不必要的高精度恢复和重复运行时检查，提高大规模 QMC 速度，同时保持权重、相位、Metropolis 和测量定义。主要改动在 inc/pure_projector_green.h 与 pure_projector_fast.h。

主模型负责分析、核心实现和最终审查。简单代码修改可以使用 gpt-6.1-sol low；验证、数据和 HPC 操作由主模型完成。用户要求不做繁冗验证；针对本次改动执行保存构型、独立 Green 参考和短链必要对照，记录实际速度。

最新状态见 docs/STATUS.md，代码地图见 docs/CODE_MAP.md。旧历史记录来自 fast-1002，引用旧 audit 时到 ../fast-1002 查找。
