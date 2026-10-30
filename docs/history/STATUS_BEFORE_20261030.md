> 历史归档：以下“当前”、默认值、作业状态和命令属于原记录日期。当前入口见 [仓库 README](../../README.md)。

# 2026-10-30 仓库更新

本次同步 2026-10-28 已验证冻结源码中的工作区/因子复用及参考传播检查，
新增近期测量与可离线复算的紧凑原始矩。更新范围、版本限制和实际检查见
[更新说明](../updates/20261030-source-and-results.md)、
[结果摘要](../../measurements/20261030/README.md) 与
[验证记录](../../evidence/source-sync-20261030/README.md)。
以下内容保留为对应日期的历史状态，不能据此推断当前队列或采样已平衡。

# 2026-10-26 当前源码同步

主线行为和范围见 [更新说明](../updates/20261026-source-sync.md)，实际验证见 [导出记录](../../evidence/source-sync-20261026/README.md)。physical trial 默认启用；总能量与 Γ 更新默认关闭。历史测量和作业保留各自冻结版本，以下旧状态保留为来源记录。

# 2026-10-13 当前更新：核心数值默认

production-v1 已统一 combined/deferred/prefix/Cayley/center-fast；普通脚本可以省略逐项数值开关。--numerics reference 及旧显式对照有效。guided、q 池与物理参数没有被硬编码。

12次固定seed短程CLI运行、8组原始块/路径/计数一致性对照和故障注入检查通过；证据见 evidence/numerics-policy-20261013/。现有任务及历史测量不改写；不认证平衡或统计效率。

下文保留原 2026-10-10 状态记录。

# quick1010 当前状态：2026-10-10

## 当前更正 note：零场模长及 Binder 峰探针已实现

用户更换了 `new/note.docx`，当前 SHA256 为
`6c78e4ceec66453972f66b0c64736c7c98d66f6012ffcf66a4fb4980adfa94c2`；
以下旧静态 note 的工作保留为独立参考，当前任务以 ADR-019 为准。
正确 docx 的7个网站链接位于 Word 字段，已提取并核对主要原论文/ALPS教程。

新 `--vbs-amplitude 1` 保存全允许 q、各键的逐构型 A/A2/A4；可选
`--vbs-amplitude-response 1` 保存先取模再梯形积分的 IA 和中心 A*IA。
默认均关闭，需 VBS/all，测量不改变零场权重、trial 或提议随机流。
新分析入口 `analyze_vbs_amplitude.py` 固定非零模式，计算 VA、chiA、B4，
逐链重分块、联合 jackknife/bootstrap、局部一/二阶 Binder 斜率和窗口依赖，
保留复数与相位/M2诊断、二倍块长误差及网格候选峰。
离散模长源独立定义与验证，不使用旧线性 Mehler score。

本次原生 MKL/GCC15 专项 PASS：独立 DFT 最大偏差6.94e-17，复权重
离散源差分最终偏差2.71e-10；6万自由声子配置的复/实模式 Binder 分别
2.00286/3.02164。12项 Python 专项测试 PASS；解析自由实模长响应中
dt=.2/.1/.05/.025 的连续积分偏差分别7.10e-4/1.81e-4/4.60e-5/1.16e-5。
三组原生新旧关闭/开启对照：旧 phase/原始测量列逐列一致、旧派生量一致，
与 h=0 场测量的 A2/A4 一致；L6 全维度 smoke 与参数依赖拒绝检查 PASS。
五个 lambda 的短 MC 及 JSON/CSV/PNG 后处理全流程已运行；它们仅验证接口，
不认证平衡、投影收敛、物理峰、临界点或 SUSY。

本机新入口 `build-amplitude-20261010/run_mc.sh`；既有 HPC 冻结作业及旧
二进制不产生新矩，保持原设置。未为此提交新的 HPC 生产扫描。
保真易感率按 note 的可选下一阶段保留推导方向，尚未实现；Pfaffian 两半
导数需包括混合 log 权重项，不能直接复制展开计数公式。
用法见 `../VBS_AMPLITUDE_README_CN.txt`；完整原始证据见
`../../fast-1006/audit/1010/vbs-amplitude-peaks-corrected-note-20261010/REPORT.txt`。

## 旧 note 静态软模工作（更正后保留为独立参考）

2026-10-10 读取 `new/note.docx`；文件没有外部超链接，旧工作区路径由
用户指定的 quick1010 覆盖。固定共同均匀背景的静态三键混合曲率筛选合理，
范围限于静态能量近似，不直接确定有限频率量子相变点。
修正既有 `static_vbs_response.py` 的 k=1 硬编码、未实际执行的差分和
仅占位的参考网格检查；采用 note 的 mass1.7/omega.8/k1.088/g.3，
普通 q 六个 cos/sin/键模式，自共轭 q 三个模式，零模显式保留异常与能量扫描。
200/400 参考改为分批解析顶角，保留完整两个 Fourier 顶角；新增同 q
几何比较、实际三档完整能量差分、测量候选 q-file 和 PNG。

本次 13 项测试 PASS；独立生产 C++ 模型导出的三种几何 A 矩阵与工具
逐元素一致，Majorana 能量归一化差<=2.85e-14。105 个非零允许 q 的
实空间与有限 k Bloch 矩阵最大差2.23e-16；三档差分 (.001,.0005,.00025)
完整 Hessian 最大偏差3.85e-8。共同背景使用 800 BZ 网格，400→800
Xbar 最大变化1.83e-6；12个所选 q 的 200→400 参考最小曲率最大变化
3.22e-5。这些网格变化不是严格误差界。全套静态计算约7秒，单线程本地完成。

各几何静态最小曲率：L6 1.0011679668、12x3/shear0 0.9895910824、
12x3/shear3 0.9990249471，均正。L6 最软四点包含 Q 附近允许点及其
共轭；未剪切12x3的最软点为 ±(.25,1/3)，偏离预测 Q；剪切最软点为
±(7/12,5/12)，与 Q 附近候选相同。剪切使最近 Q mismatch 从
0.34961 降到0.01956；共有 q=(2/3,1/3) 的有限-400参考偏差，
未剪切12x3为0.003753，剪切为0.002643。该同 q 点改善不外推全网格或
所有尺寸；尚未搜索连续 BZ 的全局最小值。没有发现此参数背景的非零 q
局部静态失稳，不据此否定其他参数、背景或量子转变。

设计见 ADR-018；完整 JSON/CSV/PNG、测试、独立矩阵、实际命令和评估见
`../../fast-1006/audit/1010/static-vbs-note-20261010/REPORT.txt`。
本次仅更改静态 Python 工具及文档，生产采样 C++、优化设置、trial、
既有作业与冻结源码保持原记录；本次未提交新的 HPC 作业。

## 新增 L6、tau10 guided 长链对照（2026-10-10）

按用户要求追加 156986/156987 同设计的 L6、tau10 测试：
157002[] 为 lambda2.4 四对，157003[] 为 lambda1.3 四对，
seed106101600–106101607，两臂 bond-refresh / bond-guided bare、beta1。
每链 warmup500、measurement2000、measure1/block1，dt0.1/M100；
固定 q=(3,3)、every1/beta1，combined/defer-phase/prefix/Cayley 全显式开启，
沿用原冻结程序中无需开关的 lazy half-step 优化；动量池是替代方案，未叠加。
保持原 linear chi/energy，fermions0；未改物理态或增加新版模长观测。
沿用 156986/156987 同一验证通过的 frozen candidate（binary37cd275e…、
source fbe588d0…），两个提议臂使用同一程序；157001 为新 L6/tau10 smoke 门禁，
长数组 afterok 依赖并检查门禁状态。157001 smoke 已 PASS/exit0，两提议均 complete/数值失败0；
末查157002四子任务R，157003两R两Q，长计算已启动；
所有包/程序/依赖校验通过，未改既有作业，未认证平衡或最终等误差加速。
完整参数/来源/提交记录见
[审计](../../../fast-1006/audit/1010/quick1010-guided-L6-tau10-20261010/REPORT.txt)。

## HPC guided 原生验证通过，对照链已运行（统计验收待完成）

2026-10-10 用户授权验证计算，已在 20015/workq 提交：156985 原生构建及
五项专项验证；156986[] 为 L9/lambda2.4 四对、156987[] 为 L9/lambda1.3
四对，两臂分别 bond-refresh 与 guided bare/beta1，均 tau8/M80、
warmup500+measurement2000、固定单 q(5,4)；156988[] 为 L6 四组新旧源码
短回归。合计24条数组链，另有构建阶段专项验证与 smoke。
每个数组子任务2CPU并发两条链，数学库1线程，walltime720h，无显式内存。
初始核查156985为R、三个数组父任务H/子任务Q，均依赖156985成功。
随后156985完成、exit0，五项原生专项验证及L6新旧/guided smoke均PASS；
三个数组12个子任务全部R，24条链均启动，初始数值失败记录为0。
这不表示长链统计、收敛或效率验收通过。参数、源码身份、
实际作业号与状态见
`../../fast-1006/audit/1010/quick1010-guided-validation-20261010/REPORT.txt`。
使用本次整合时冻结的源码，tar SHA256 fbe588d0…；既有任务不变。

## 当前源码与本次整合

用户于 2026-10-10 指定 quick1010 为后续代码修改、优化基线。
原始来源为 fast-1006-core-current-20261010.zip，本次审查并选择性整合
`new/fast1006_guided_sampling_source_20261010.zip`。
新增固定参考整键 pCN、固定动量池及多键模式、CLI 元数据与验证目标；
补齐免 CMake 构建入口。Hamiltonian、trial、物理声子权重、相位重加权、
Pfaffian/Green 内核和已有科学估计器保持原实现。

本次在原生 MKL、GCC15、-O3、数学库单线程环境完成新提议验证及
旧 local/bond-refresh、q 开/关四组兼容性检查；原始 blocks/预热 CSV
逐字节相同，共有物理 JSON、最终路径与计数一致。
零参考力/beta1 与旧整键刷新一致，guided+xyz池及物理源 CLI 测试通过。
guided、q-mode、MC、numerical-failure、VBS susceptibility 五个原生专项
验证全部通过，生产程序与这些验证目标均成功编译。L9/M8/tau.8 的
warmup1+sweeps2 尺寸检查完成、数值失败为零；它不是 tau8 生产链。
完整结果以 `../../fast-1006/audit/1010/quick1010-guided-integration/REPORT.txt`
为准。本机入口为 `../build-guided-20261010/run_mc.sh`；同目录构建记录
及本次审计源码快照用于追踪版本。

未来普通整键生产配置按 AGENTS 显式选 guided bare/beta1 与原单 q 配置。
CLI 默认 proposal=local 仍保留；动量池另按目的选用。
本次没有重跑上传包的性能 benchmark，未认证 L9 加速、热化、长链或投影
收敛，也未替换既有作业的冻结程序。前后源码、diff、命令与原始验证记录
保存在上述新审计目录。原始及整合出处见 `../BASELINE_ORIGIN.txt`。

## 继承的历史状态记录（以下保持原文）

# fast-1002 当前状态：2026-10-02

## 版本和工作区

工作目录为 `/home/asus/newsusy/fast-1002`，完整复制自 fast-0929。
来源包标记 `2e53aee`；本目录不是该提交的 Git checkout。
用户已明确授权逐项审查并迁入 new_0923 的适用后续改进；HPC 生产二进制
的来源和本轮源码验证分别记录。

## 源码范围

生产程序是固定 g、固定规范 matter 扇区的 SSH–Kitaev 声子 PfQMC，
支持开放声子路径、剪切周期边界、VBS 测量及同候选无缓存重试。
CLI 默认 `mass=4`、`omega-ph=0.5`、`J=(1,1,0.5)`、`g=0.3`，
`proposal=local`、`weight-engine=combined`，延迟相位和 prefix/Cayley
缓存开关均为 1；实际缓存路径随引擎及提议变化。
生产 honeycombGenerator 捕获固定模型；通用 density 驱动程序
不能作为 SSH 声子生产时间驱动已经接入或验证的证据。

## 本工作区已经验证的内容

初始复制后核心哈希和权限与来源一致；下列为首轮已验收记录。
再次核对 fast-0929 的 349 个来源文件，全部哈希仍保持原值。

- Nc>=32 且精确二分格结构检查通过的 honeycomb 节点，以实 Nc×Nc 完整 SVD
  替代 2Nc×2Nc 复谱构建；magnitude overlap 使用持久 zgesvd 工作区。
  不截断模式，不放宽原稳定性/物理门槛；signed Pfaffian、测量 Green、
  generic/canonical、eager local h*h 与无缓存重试路径保留。
- 新增矩阵函数/Fock/SVD 边界测试通过。原包缺失的旧测试 oracle 补齐；
  两个旧测试改为认证后读取 lazy Green，原阈值保留。baseline CTest 23/23、
  candidate 25/25 通过，原始失败与修正 diff 保留。
- 四组固定 MC 回归路径/计数一致，最大区块误差 2.84e-14。
  HPC L6 200 samples、L9 60 samples 的全 VBS/费米子输出比较通过，
  最大区块误差 1.71e-13，数值失败零，最终路径/计数完全一致。
- HPC20022、单 CPU、数学库线程1、M80 的 L12 完整一轮采样
  1475.96 → 644.84 秒（2.29×）；三次同 CPU 配对采样中位加速 2.42×。
  含较长回放全部测量输出的 L6/L9 总时间加速分别为 1.68× / 1.99×。

以上为固定参数的效率和回归证据，不提供采样/投影收敛、独立 block、
相位优势、相变或 SUSY 证据。优化后 L12 约 51.5% 采样时间在 magnitude QR/SVD，
后续替换需要新的数值验收。

- [复制检查](../../audit/core-copy-20261002/REPORT.md)
- [当前优化与验收](../../audit/optimization-20261002/REPORT.md)、
  [ADR-007](../decisions/ADR-007-bipartite-spectrum.md)
- [历史报告汇编](REPORTS.md)：保留原版本的结果、参数与证据边界。
- [代码地图](../CODE_MAP.md)、[物理约定](../PHYSICS_CONVENTIONS.md)、
  [测量约定](../OBSERVABLES.md)、[数值不变量](../NUMERICAL_INVARIANTS.md)。

后续优化应记录具体变更、基准条件、数值回归和实际耗时，
仅在新证据产生时更新状态。

## 第二轮当前实现与验收

已审查 new_0923 当前未提交代码、20260930 完整源包及后续运行来源记录。
迁入小尺寸完整 SVD 分派及 measure-energy 按需测量，保留 fast-1002 物理默认值。
本轮新增认证实数 magnitude 传播/QR/SVD与有界缓存；signed phase、所需Green
和同候选无缓存恢复保持。精确结构不符或原数值门槛附近回放复数候选。

本地完整编译成功、CTest 27/27、8/8 默认及energy0输出对照通过；
最大物理JSON误差1.07e-14、CSV1.42e-14，路径计数一致、失败零。
新增实数测试46组直接比较、40次提议；最大log ratio误差3.98e-13。
HPC 配对测速和较长全观测量回放正在验收，尚无本轮实际加速结论。
详见 [第二轮审计](../../audit/optimization-round2-20261002/REPORT.md)、
[ADR-008](../decisions/ADR-008-real-magnitude.md)。

## 最新生产提交

HPC20022 最新优化二进制 c1132b5f 已通过节点 smoke 和四项专项验证。
按用户要求无PBS验证依赖提交60条L9/tau10链，每seed15000测量：
workq 86065[] 为lambda2.4、2.5各10seed；仅2.4中5seed测量费米子关联。
testq 86064[] 为lambda2、2.2、2.6、3各10seed，费米子关闭。
提交不代表计算完成；测速86063[]仍单独运行。
见[提交记录](../../audit/production-L9-tau10-m15000-10seed-20261002/submission.json)。

## 2026-10-03 进度核验

生产60/60链仍R，写出区块含260–470测量/seed（目标15000），尚无完整链，60份运行状态均numerical_failures=0。
第二轮12组配对短基准通过既定比较：最终路径和接受计数一致；log绝对误差最大1.14e-13，phase距离最大1.69e-15。
相对首轮优化基线，12次bond提议/M80的三次配对加速中位值：L12 1.282倍，L9 1.214倍，L6 1.110倍，L4 2.193倍。
这是有限提议测速，不是完整sweep或15000测量生产链总耗时。
[验证原始记录](../../audit/optimization-round2-20261002/hpc/results/paired-profile/paired-profile-verification.json)。

2026-10-03追加testq 86066[]：lambda2/2.2/2.4/2.5/2.6/3各8seed，共48链，全部fermions1，L9/tau10/每seed15000测量；旧60链保留，无PBS依赖。记录见[追加提交](../../audit/production-L9-tau10-m15000-fermions-8seed-20261003/submission.json)。

用户随后取消86066[]全部48链，已确认结束。替代testq 86067[]仅lambda2.4/2.5各10seed全部fermions1，其他参数不变。原始86064[]/86065[]共60链继续运行。

2026-10-03当前热图已直接在HPC读取并汇总80条有效链（原60+替代费米子20），取消86066排除。
仅绘制Szz_connected(k)/Nc，Nc81，标记裸J=(1,1,0.5)狄拉克点差。六个lambda共用线性色标，无插值；当前每链260–1100测量，warmup0/no cut，仍为阶段性快照。
完整矩阵Hermitian/反演残差检查通过；汇总、原60、费米子20分别留图。
数据及可复现脚本见audit/heatmaps-current-20261003；最终仅用hpc_reduced.json，未完成raw传输不纳入。

2026-10-03用户要求跨版本10000sweep对照：20022 testq 86073[]，fast0929/fast1002首轮/new0923-20260930每版2seed，共6新链；lambda2.4/L9/tau10，seed10024010/11，VBSall+fermions1。最新版复用86065[] index10/11前10000测量，不重复提交。尚待QMC完成和固定阈值比较。记录audit/qmc-version-consistency-L9-lambda2p4-tau10-10000-20261003。

2026-10-03用户要求准备L9 lambda2/2.4/2.5/2.6/3 x14tau扫描。250条链/70参数点，fermions0；tau<4每点3seed各3000测量，tau>4每点5seed各2000测量。参数与PBS脚本静态核验通过。用户随后明确只准备不提交，未上传或qsub，无新job ID；计划20022workq后迁实际Q子任务到testq。见audit/production-L9-lambda-tau-scan-20261003/preparation.json。

2026-10-04扫描准备更新：tau<=4.2各4seed/3000测量（220链），tau>4.2各6seed/2000测量（90链），共310链/70点/840000测量。新准备包audit/production-L9-lambda-tau-scan-20261004替代Oct3旧准备；旧包留档，未上传、未提交。

2026-10-04新增L6/lambda2.4绝热静态参照：J=(1,1,0.5), k=1, g=sqrt(2.4)，零通量Wilson(+,+)。独立全108键变量优化，每种matter奇偶约束128初态； unrestricted与QMC自由试探态同奇偶两种搜索均得到均匀最低已找到解。Etotal=-84.9942873901，Ematter=-120.8288913304，Eelastic=35.8346039403；Xbar=(0.9302793735,0.9302793696,0.5098744078)。非零允许动量Szz/Nc为数值零；不能推断无限尺寸无VBS。最小Hessian特征值0.01440556，gradient最大4.82e-10；独立实空间Hermitian谱和Bloch能量核验误差1.43e-14；300个随机扰动均增能。多初态结果不是全局最优数学证明，也不是完整物理自旋投影。见audit/adiabatic-L6-lambda2p4-20261004/results.json及verification.json。

## 2026-10-04 L4 lambda/tau scan submitted

User revised L6/8000 to L4/10000 before production submission and explicitly selected port20022 workq. Latest validated c1132b5f binary; 2 lambda ×15 tau ×4 seeds=120 chains, m=omega=k=1, dt=.1, warmup0, energy+VBSall, fermions0. Node smoke86125 passed10 samples, zero numerical failures and PBS exit0; first smoke86122 failed before sampling from4 sweeps/block10 mismatch, retained diagnostics. Production arrays86146[],86147[],86148[] submitted longest tau group first,40 chains each. Initial states {'R': 120}. Audit: audit/production-L4-lambda2-4-tau-scan-m10000-4seed-20261004/submission.json. Computations ongoing; no convergence claim.

用户要求缩短omega扫描tau后重新提交20015：156772[] workq，L6/lambda2.4，omega0.4/0.2/0.1/0.05分别tau5/10/20/40，各2seed、10000测量、warmup0、dt0.1、VBSall、fermions0、Majorana能量1。最新二进制哈希核验通过；8/8子任务初查R。每链先做节点2样本smoke，成功后开始生产；初查smoke尚未写出，不声称已完成。原20022长tau扫描86075[]保留。详见audit/omega-scan-L6-lambda2p4-tau5-40-20015-20261004/submission.json。

2026-10-05短tau扫描复查与修复：原20015 156772[]全部在smoke前退出127，缺GCC8.3运行库。156773[]遇oneAPI初始化非零返回；156780[]节点smoke本身complete/零失败，但默认measure-every2仅1样本，与2样本门槛不符，未进生产。均保留失败目录。补齐四项GCC运行库、保留环境初始化日志并用ldd检查实际依赖、明确smoke measure-every1后，重试156781[]在20015workq已8/8节点smoke通过(samples2,failures0)、8/8生产R。实际argv核对8/8通过：L6/lambda2.4, omega0.4/0.2/0.1/0.05对应tau5/10/20/40，各2seed、10000sweep、warmup0、dt0.1。已写100/50/20/10样本每seed。审计audit/omega-scan-L6-lambda2p4-tau5-40-20015-retry-20261005/submission.json。

## 2026-10-05 L4 old-image comparison

Collected120 production result JSONs on20022;119 complete with numerical_failures0; lambda4/tau14 seed10043001 aborted at uncertified zero_weight after canonical retry, excluded (3/4 chains at cell). Digitizedlambda2/4 on15 matchingtau points from user image, natural-log e^-n axis calibration. Raw Szz/Nc maxima inclusive/exclusive q0; independent-chain means/SEM and separately pooled phase-reweighted spectrum maxima exported. Curves broadly close at lowtau; longtau deviations and incompletecell prevent an equivalence claim from image alone. Evidenceaudit/old-image-vs-L4-fast1002-20261005.

2026-10-05新增高omega扫描：20022 workq 86680[]，L6/lambda2.4，omega=1,2,4,8,16，每点2seed、10000 sweeps、warmup0；固定omega*tau=2与omega*dt=0.1（20 intervals），故tau=2/omega、dt=0.1/omega。10/10节点smoke完成(samples2,numerical_failures0)，10/10正式结果进入running且零数值失败。首个提交86673[]因smoke路径误用M=1而被程序拒绝，未产生生产数据；成功重交86680[]。见audit/omega-scan-L6-lambda2p4-highomega-20022-20261005/submission.json。

高频扫描86680[]完成复核：10/10链PBS exit_code=0，均status=complete、samples=10000、numerical_failures=0；最终日志时间2026-10-05T07:33:49Z UTC。结果同上目录远端workq/runs；submission.json更新完成状态。

2026-10-05绝热参照复核：重新读取SSH–Kitaev原论文正文式(1)–(3)及补充材料B/C，核对固定k,g,J、m→∞下最小化费米子能量与弹性能之和的定义。采用不导入QMC或静态能量代码的标准库动量求和，L6/Wilson(+,+)已保存均匀解的Majorana能量为-120.82889133041127，与原参照相差1.42e-14；弹性能35.83460394032415、总能量-84.99428739008712。此计算验证该构型的能量，不认证其全局最优性。原128起点优化及局部稳定性检查仍是静态最低解的现有证据。有限频率的Majorana能量不受该分量数值的下界约束；基态omega→0比较需分别确认采样、投影和离散收敛，固定omega*总轮廓长度不保证投影收敛。高频CSV各链分段均值已复核，但未据此认证平衡。复现与来源哈希：audit/adiabatic-reference-review-20261005/independent_check.py及results.json。

2026-10-05统一高omega tau12重算：20022workq初交86807[]，omega16/8/4/2/1对应index0-9，各2seed、10000sweep、warmup0，dt0.1/omega。节点smoke10/10通过。初交index2-7权重认证异常退出；完整reference/eager、qr-block1、关闭prefix/cayley缓存保持物理参数重交86808[] index2-7。最新初交0/1/8/9及重交2/3/6共7条R；重交4/5/7仍zero_weight/green_residual=inf异常abort，失败数据不可作正常测量。尚未完成全10链；见audit/omega-scan-L6-lambda2p4-highomega-tau12-20022-20261005/submission.json与startup-check.txt。

2026-10-05用户要求准备omega0.3相变峰扫描，明确不提交。已在audit/prepare-omega0p3-susceptibility-20261005生成未确认的L9/三lambda/四tau/每点4seed参数草案和不可提交的PBS草案，k1映射mass100/9、参数及bash语法静态核验通过。现场CLI检查无极化率/虚时间积分响应开关；susceptibility_enabled=false，尚需确认chi定义与扫描范围，再实现和验证测量。未上传、未提交、未运行MC、未替换生产二进制。

2026-10-05 VBS静态响应准备完成：用户确认lambda1.3/2/2.4/2.6/3和17个总tau，10000sweep；沿用L9、4seed、warmup0。新增默认关闭的--vbs-susceptibility，all允许动量网格上保存精确Mehler外场score与中心D的24个原始矩/q；固定trial端点，响应连接减法使用共同复phase分母、独立保存conj(I)和实源一阶矩。GNU/MKL本地验证通过：外场有限差分最大误差8.56e-11，自由有限投影解析响应误差6.67e-16，直接Fourier误差1.19e-16；M2/M6每例32测量的off/default/on与旧修复binary的采样路径、接受/phase、全部旧原始区块和VBS派生字段逐项一致，5项非法配置拒绝。合成复phase汇总核验误差1.80e-15，另用两条32测量短链验证汇总接口；不认证生产收敛或自关联独立。340链/85点的本地准备包冻结，17个tau数组每元素4CPU并发4seed，无mem，walltime720h；完成时尚未上传/提交。见audit/production-L9-omega0p3-vbs-response-20261005/preparation.json，定义见ADR-009。用户随后授权：tau<4改5000，只先提交tau>4到20015，排队再20022workq；这属于后续执行阶段，原准备包保留。

2026-10-06 tau12假零权重修复：86808 index2–7最终均因full-rank QR归一化重叠行列式小于绝对阈值而abort。inc/pure_projector_fast.h改用SVD对数幅值与LU逐枢轴单位相位，不以已通过rank/rcond检查的奇异值乘积大小判物理零权重。原rank/rcond/residual/Pfaffian容差及失败即停止政策保持；Hamiltonian、随机流、MH公式、物理参数不变。新增18对Majorana旋转精确振幅回归，本地五组验证及六个保存候选的合并/未合并/完整测量通过；HPC验证86810 exit0同样通过。编译86809的GCC8 filesystem验证链接问题在验证重试显式-lstdc++fs解决，原日志保留。新生产二进制f03bfc35…；20022workq重交86811[] index2–7，每omega8/4/2各2seed、L6/lambda2.4/tau12/10000sweep/warmup0，原种子不变。原86807的omega16/1共4链仍运行。本次使用明确冻结的数值修复快照，不覆盖后续并行新增VBS静态响应功能，该功能未纳入此批二进制。见audit/tau12-numerical-repair-20261005/REPAIR.md、hpc-validation.json、submission.json。验证通过不代表生产完成或投影/热化收敛。

2026-10-06 omega0.3扫描提交与队列迁移完成：L9，lambda1.3/2/2.4/2.6/3 ×17总tau ×4seed=340独立链。tau>4共140链各10000sweep，每cell4CPU并发4seed；tau<4共200链各5000sweep，初始单一大数组156824[] index140–339，每元素1seed/1CPU。VBS静态外场响应开启，warmup0、dt0.1、energy1、fermions0，无mem、walltime720h。先20015提交，按新鲜Q状态逐子任务取消/核验再20022 workq提交；最后按用户要求将86820[]仍Q的211–339共129链转交86821[] testq，保留已R子任务。最终快照20015 workq59R，20022 workq152R，20022 testq63R+66Q；274/274已启动链节点smoke通过，当前wrapper/numerical failures均0，340seed唯一分配无重复活动子任务。初次长tau156805[]–156811[]因TSV CRLF尾部导致CLI smoke退出，140链均未进生产；LF参数表与防御CR剥离修复后两端完整表驱动wrapper smoke通过，失败尝试保留排除。全部提交完成不等于QMC完成或热化/自关联/投影收敛。权威记录与逐seed分配见audit/production-L9-omega0p3-submission-summary-20261006/submission.json和seed-assignments.tsv；各端完整快照见短tau数组审计目录。

2026-10-06长链复查：tau12共10链尚无完成；86811重交index2/3/4/5/7因signed measurement partial Green的solve_residual_exceeded再次abort，combined与uncached reference均不通过，残差约1.0e-10–8.9e-10，原1e-10门槛未放宽。原假零权重修复与保存构型验证成立，但未解决后续长链的全部数值稳定性。86811 index6约1800/10000测量仍R；86807 index8/9约4410/4400测量仍R，index0/1调度R但未见生产区块。失败候选另需诊断，不能把这些部分结果算为完整计算。见audit/tau12-numerical-repair-20261005/progress-20261006.json。

2026-10-06后续队列再分配：按用户要求从86821[] testq新鲜Q子任务取25seed到20015workq 156834[] index274–298、25seed到20022workq 86823[] index299–323，均25R且节点smoke通过；随后将剩余16个Q全部移20015workq 156835[] index324–339。最后快照20015workq共100链(93R+6Q+1complete)，20022workq177R，testq保留63R、无Q。已完成index329为5000samples、数值失败0、PBS/wrapper exit0；全部340seed唯一分配，取消源子任务均不再活动，原运行链保留；334条已启动链全部smoke通过、当前无wrapper/numerical失败。参数不变，每seed1CPU，5000/10000sweep规则保留。收集时务必加入两个rebalance25根目录和rebalance-rest根目录，完整八根路径与逐seed分配见audit/omega0p3-queued-redistribute50-20261006/submission.json、seed-assignments.tsv。先前记录保留为历史快照。

2026-10-06 omega0.3当前曲线：直接汇总两端八个活动结果根目录340/340 seed、85/85参数点，955930个完整区块测量，77条complete、18个4/4完成cell。每链540–5000测量，warmup0/no cut；长tau14各seed540–600。运行JSON仅含argv/状态，使用实际argv与manifest核验参数，并按CSV文件打开时字节长度上限只读取完整有限区块；最终无排除/无重复。共同phase原始矩先合并再按OBS-03只扣均匀背景，非零允许q上保留完整Szz/Nc；VBS静态chi按ADR-009逐q连接减法，独立保存conj(I)均值，不减自由声子响应。主图五lambda的非零允许81网格中80点峰值随总tau变化；另固定允许q=(5/9,4/9)作比较。1SE为四独立seed删整链jackknife，峰值在每次删除后重算；不认证自关联/热化/投影收敛或相变点。复phase合成独立直接平均验证最大误差1.69e-16，85cell图值/网格导出一致，无非网格独立QMC或插值。PNG/PDF/逐q数据/逐seed样本数与可复现脚本见audit/omega0p3-curves-current-20261006。

2026-10-06 lambda2.6逐seed sweep诊断：从两端当前CSV完整区块重新读取tau5.2/6.4/8/10/12/14各4seed共24链，每链580–1650测量，零排除、独立seed不混合。主图为每seed累计共同phase重加权完整Szz谱之后取非零允许q最大/Nc；另有从sweep100显示的放大图(累计仍从1开始)和固定q=(5/9,4/9)的10sweep区块图。末点累计谱与直接总和误差0。图显示长tau仍有seed差异/漂移，不认证热化或自关联。原始区块抽取、CSV、PNG/PDF及复现脚本见audit/omega0p3-lambda2p6-seed-sweeps-20261006。

2026-10-06 tau12 Green 异常退出修复：86811 index2/3/4/5/7保存路径因满秩病态重叠的double SVD求解残差超限停止。pure_projector_green.h增加失败后LU精化和由原子空间一致形成overlap/solve/Green的long double恢复，原rank/rcond/solve/Green门槛及signed层已有oblique相对检查保留；正常通过案例原运算逐元素保持。无Hamiltonian/Metropolis/相位/测量/物理参数改动。主模型本地七组验证通过，48个合法病态Gaussian案例独立MP50最大相对误差8.44e-14，五退出路径合并/未合并完整测量10/10通过。HPC88009构建成功，Bash4.2空数组nounset包装器故障在88010显式命令重试；88010 exit0且11保存构型两模式完整测量22/22通过，HPC MP50最大误差8.13e-14。节点88011 smoke complete/2samples/零失败/exit0。新binary e1c3ee8c…，冻结源清单b7c63b38…，保留新增VBS静态响应功能。20022workq重交88012[] index2–5及88013单链index7，原seed不变，L6/lambda2.4/tau12/10000sweep/warmup0；提交后5Q。原86811[6]及86807[0/1/8/9]仍R且保留。验证与重交不代表长链完成或热化/投影/离散收敛。详见ADR-010及audit/tau12-green-repair-20261006/REPAIR.md、hpc-validation.json、submission.json。

## 2026-10-06：运行中接受率输出（仅代码，本轮未提交 HPC）

生产 CLI 新增 `--acceptance-every N`，默认 0 关闭；开启后每 N 个
已完成 sweep 和阶段末尾刷新 `PREFIX_acceptance.csv`。分开记录预热和
测量累计接受率及数值恢复计数，测量计数扣除预热。保留现有 LU 精化
和扩展精度恢复实现。

本轮本地编译成功；L=2 固定 seed、预热 3 / 测量 7 sweep 验证：
原基线、新版关闭输出、新版开启输出的原始 blocks 文件逐字节一致，
最终路径及接受/提议计数一致；CSV 阶段末尾与结果 JSON 计数一致。
另一次 200 sweep 测试确认进程仍在运行时可读已刷新行；非法负间隔
及既有接受率文件覆盖均被拒绝。证据在
`audit/live-acceptance-code-20261006/validation.json` 和 `validate.py`。
本轮没有提交任务，也没有更改已有 HPC 计算。

2026-10-06追加同一快照tau10五lambda静态极化率六边形图：逐q Re chi_zz保留全部81允许点含Gamma，用有限9×9盒子归一化|W|²核重建连续展示，统一线性色标0.845604–1.529540。直接动量点/裸J狄拉克差/非零允许点峰代表分别标点/星/圈；非网格像素不是新QMC数据。网格复现残差0、周期权重残差3.78e-15、权重和1且无越界过冲；各点每4seed630–940测量，未认证热化/自关联/投影收敛。PNG/PDF/网格表/核脚本见audit/omega0p3-chi-heatmaps-tau10-20261006。

2026-10-06追加L6 omega0.3极化率扫描：lambda1.3/2/2.4/2.6/3 ×总tau14/10 ×4seed=40独立链，各15000sweep、1CPU，dt0.1、warmup0、VBSall+静态外场响应、energy1、fermions0，walltime720h、无显式mem。使用已在HPC验证的e1c3ee8c Green修复冻结binary，保留chi；不替换为当前fast1006优化候选。20015实际L6/tau14节点preflight156855 exit0且2samples/schema通过；先提交长tau数组156856[] index0–19，再156857[] index20–39，均workq。用户明确覆写本批fallback：即使Q也全部保留20015，不转20022。生产未取消、未迁移；仅20022预备检查88014在Q时取消并确认F，该端未提交生产。最后快照30R+10Q，30/30已启动链smoke通过、当前wrapper/numerical失败0；不认证热化或投影/自关联收敛。完整记录audit/production-L6-omega0p3-chi-tau10-14-m15000-20261006/submission.json及seed-assignments.tsv。

2026-10-06 fast-1006按需精度优化：从当前fast-1002复制271核心/构建文件，逐文件一致且原fast-1002保留。生产Green先执行已有oblique相对检查，values-only SVD保留rank/rcond后LU求解，失败才精化/扩展精度，复用LU与延迟分配；原容差/物理参数/MH/测量保持，strict默认路径保持。独立MP50生产8病态案例最大相对误差6.04e-11，strict/权重相位/失败事务/deferred既有回归通过。五退出构型的权重、相位、完整中心Green、能量与重新编译参照逐项一致，完整回放56.71→51.33秒。L6/M64与L9/M16各8测量短链采样52.18→49.05、180.58→163.23秒，约减少6.0%/9.6%；路径/接受拒绝次数及115498数值字段一致、零失败；涵盖fermion、VBSall和静态响应。单次短链不认证长链收益/收敛。局部900 Green重算31.0→1.87秒，扩展精度案例4→0，不是整链倍数。并行修改phonon_mc/path新增静态源接口已保留，本轮仅认证F=0原生产路径，不认证非零源功能。当前优化主目录fast-1006，未提交或替换HPC；见ADR-011、audit/precision-speed-20261006/REPORT.md及summary.json。

2026-10-06 有限 VBS 共轭场对照准备：用户授权 20022testq。实现固定允许 q 的
cos 源、精确强迫 Mehler 权重与一致的局部/整键自由提议，trial 随 h 固定。
新增 Field_m/score/乘积及逐构型复 Fourier 二/四阶矩、纵向 connected response
和 phase joint blocking；非零场 snapshots 保存 forces(version3)。本地强迫核、
非零 h 有限差分 score、解析自由响应/归一化/norm、两种 proposal 详细平衡及
快照回放通过。两种 proposal 默认关闭旧 blocks 逐字节一致；显式 h=0 的最终
路径和全部旧观测量一致。共8短链 complete/8samples/零数值失败。本地证据
在 audit/vbs-finite-field-controls-20261006/validation-local.json；数学检查 binary
kitaev_ssh_vbs_field_validate。78条对照计划已冻结，远端构建/节点检查及提交
状态另行追加；此记录不声称生产提交或临界点已确定。见 ADR-012。

2026-10-07用户修订有限场尺寸对照：不要L12，改成L6。现行78链中66条L9、
12条L6；主扫描48、L6尺寸12、tau14投影6、有限差分8、相位2、dt2。
L6固定允许q=(3/6,2/6)，为蜂窝倒格度规到裸Dirac差最近并列点之一；不与L9
q=(5/9,4/9)混称同一动量。未提交的L12旧计划保留作记录，已被现行points.tsv取代。
20022testq构建/实际L9-M100-h.05节点预检88015 exit0，2samples且零数值失败，
强迫核/score/自由解析与详细平衡验证通过；binary153b848ac17c…、冻结源c4bc201e…。
本条只记录构建证据；生产任务提交情况由下一条追加。

2026-10-07有限场78链实际提交20022testq：88016[] index0–7(较长L9)，
88017[] index8–65(L9/tau10)，88018[] index66–77(L6)，每链1CPU、
walltime720h、数学单线程、warmup1000/measurement10000。最终没有L12。
初查64R+14Q，64已启动链节点smoke通过、生产argv逐项与冻结计划一致，
PBS非零退出/当前数值及wrapper失败均0。未声称生产完成或临界点已确定。
运行/任务/参数/seed/源与binary校验记录在audit/vbs-finite-field-controls-20261006/
submission.json、live-check.json、seed-assignments.tsv。新纵向chi与既有
cos响应的同相位归一化短链核对误差2.08e-17。

## 2026-10-07：现有 omega=0.3、L=9 序参量快照

从原 fast-1002 生产任务两端口汇总全部 340 条链、85 个参数点，
196 条链完成、48 个参数点 4/4 seed 完成，共 1,348,180 个测量 sweep。
沿用 OBS-03 的均匀背景扣除、共同相位加权后取非零允许动量最大
Szz/Nc；误差为删除整条独立 seed 的 jackknife。未做预热截断，
长 tau 为未完成链快照，不能据此认证收敛或相变。
已完成链采用最终 JSON 的相位加权均值和样本数；运行链采用抓取时
已经落盘的完整块。两端服务器时钟不同，抓取本地时间约 20:06–20:08
（北京时间）。原始汇总、数值表和图保存在
`audit/omega0p3-order-current-20261007/`；本次没有提交或变更 HPC 任务。

2026-10-07 L6 omega0.3 tau10/14极化率批次只读进度复查：156856[]/156857[]全部40R，无Q、无complete，全部保留20015workq，当前数值失败0。按实际CSV完整10测量区块统计累计369960/600000测量(61.66%)，tau10各seed6450–12530、tau14各seed8050–8850，warmup0；运行JSON的samples字段尚无累计数，进度据CSV而非将缺字段当零。未改动/迁移任务。快照和逐cell seed范围见audit/production-L6-omega0p3-progress-20261007。

2026-10-07 L6 omega0.3序参量六边形热图：20015该批40seed/10cell全部汇总，共371380完整区块测量，各seed6500–12580，仍0complete，warmup0/no cut。绘制lambda1.3/2/2.4/2.6/3、tau10/14的Szz_connected(k)/Nc，Nc36；共同phase矩先跨seed合并，再只扣均匀系综均值(OBS-03)，保留Gamma连接涨落。全部36允许动量为独立QMC数据，连续显示使用6×6有限盒子归一化|W|²核重建；共有线性色标0.004617826–0.007459616，裸J狄拉克点差标星、非零允许点峰标圈。允许网格复现/反演残差检查通过，权重和残差2.11e-15，无核重建过冲。逐q四独立seed删整链jackknife1SE另存；未认证热化/自关联/投影收敛。组合图、两tau单图、PNG/PDF、逐q表和冻结原始矩及复现脚本见audit/L6-omega0p3-order-heatmaps-current-20261007。


2026-10-07 fast-1006后续加速与降精度：本轮从当前fast1006基线继续，
生产signed QR/values-only SVD使用LAPACK并复用SVD工作区，Green工厂传递
加速设置。生产solve残差由1e-10放宽到1e-9，条件好时按实测rcond最高1e-8；
Green残差仍1e-8，rank/rcond及独立strict默认保持，矩阵仍double，失败仍恢复。
最终重编译L6-M64/L9-M16各4测量短链24.97→22.16、85.54→78.39秒，
约减少11.3%/8.4%；最终路径与接受计数完全一致、零失败，81286个JSON/CSV
数值字段最大缩放差8.73e-13。12个合法Gaussian MP50对照最大Green误差
2.96e-10；既有SVD边界/deferred/real回归通过。五保存退出构型回放45.07→44.16秒，
权重/相位/Green/能量最大缩放差1.94e-14。仅L9两个由舍入级phase stderr
产生的衍生SNR单独记录，原phase/误差/物理阈值保持。实数提议DGESDD候选
无稳定收益已撤回，最终实数driver保持基线。专项验证器构建在撤回前，保留
Green/signed源码未变；最终MC与回放在撤回后重编译。单次短链不是长链收益
保证，未提交/替换HPC。最新本机程序build-signed-svd-precision-20261007/kitaev_ssh_mc，
SHA0889df75750d…；见ADR-013及audit/signed-svd-precision-20261007/REPORT.md。
指定6.1-sol low代码代理未及时返回补丁后停止，主模型完成最终实现与验证。

2026-10-07 L6极化率随lambda曲线：复用刚才热图的同一40seed/371380测量冻结快照，tau10/14各5lambda，逐q连接响应按共同phase分母计算ChiRaw−Nc<D><conjI>，不扣自由声子响应；主图为35个非零允许动量上最大Re chi_zz，四独立seed删整链jackknife1SE且每次删除重新连接/取峰。另存固定允许q=(1/2,1/2)曲线对照。tau10峰值1.08857至1.63288；tau14为1.25175/1.37628/1.54037/1.45454/1.66316，lambda2.4处误差0.12925。未据五lambda当前快照认证相变峰。PNG/PDF/逐点数据/全部36q及复现代码见audit/L6-omega0p3-chi-vs-lambda-20261007。


2026-10-07默认代码基线确认：逐文件比较当前fast-1002/fast-1006的inc/src/
kitaev_ssh/scripts/reproduction/cmake及CMake，共271个对应核心/构建文件，
262个字节一致，9个差异均为fast-1006的生产精度/速度改动、有限VBS场扩展或
相应验证/构建注册；fast-1002没有独有核心文件，fast-1006另有3个新增文件。
fast-1002/perf_work的34个旧测速/构建/补丁辅助文件留作历史工具，不含独立
更新的生产算法。未发现fast-1002独有、未迁入的核心改进。按用户要求，今后
代码修改、优化及新计算默认以fast-1006为基线；旧源码和已提交任务版本保留。
根及两目录AGENTS/README/CODEX_HANDOFF已同步，避免旧入口继续指定fast-1002。
本轮仅源码审查与规则更新，未重跑编译/QMC/HPC；核对数据和完整差异见
 audit/baseline-selection-20261007/comparison.json与core-differences.patch。

2026-10-07 L9 omega0.3 加密lambda生产提交：按用户改为tau8，lambda1.9至3.0步长0.1，每点目标10独立seed，复用原tau8扫描16条链，仅新增104条。新增链全部使用fast-1006冻结源码SHA d3ae252e2ee3…，新编译程序SHA87786fc8d90e…；已有链保留原fast-1002版本及fermions0。20015 workq数组156874[].mgt，索引0–103，每seed1CPU、数学库1线程、10000sweep、walltime720h，无显式内存。全允许动量VBS静态极化率与fermions开启，接受率每10sweep输出。计算节点L9/M80两sweep预检156873.mgt退出0，极化率及648个费米关联通道输出检查通过；提交后104个子任务均R。这仅认证启动与输出接口，不认证收敛。tau10准备包已被tau8替代，未提交tau10生产。冻结来源、覆盖及提交证据见audit/production-L9-omega0p3-chi-dense-tau8-20261007/submission.json。

2026-10-08绝热/omega结果汇总与迁移：156781低频omega0.05/0.1/0.2/0.4对应tau40/20/10/5共8链均complete、10000samples、数值失败0、wrapper exit0；86680高频短tau10链同样完整；固定tau12仅原86807 index8/9两条omega1完整。原86075已终止，仅部分区块，排除正常点。原修复重交88012[2–5]及88013仍Q，按用户要求新鲜核对Q后取消并确认终止，迁20015workq156879[2–5]及156880(index7)。作为原任务迁移保留e1c3ee8c… binary/源清单b7c63b38…，非新fast1006代码计算；GCC8运行库随包，smoke156876 exit0，5/5生产节点smoke通过且当前5R/数值失败0，无重复活动seed。20完整链原phase/energy区块重建与JSON一致，输出未切热化区间的能量随sweep、能量随omega及360行允许网格Szz/Nc，误差为2seed删整链jackknife，不认证平衡/投影收敛或omega0外推。证据、原始JSON、参数、图和CSV见audit/adiabatic-omega-output-migration-20261008/。

2026-10-08低omega远离绝热参照诊断：复核156781测量为中心Majorana能量分量，与静态-120.82889133同口径，k1/m=1/omega²/零trial shift。总tau=2/omega使单侧omega*Theta=1保持固定，不认证声子基态投影；恒力谐振子同trial的精确中心均值Xeq(1-exp(-omega*Theta))仅作为边界偏差类比，不能当作相互作用模型精确预测。低omega0.05两链接受率约4.2%，0.4约43.7%；0.05能量phase合并从前2000 sweep -81.1626到后2000 -86.8882，全程均值-84.9342，已有初始漂移证据。不能把距静态差距全归咎热化；尚需独立tau与采样收敛检查。先前tau=2/omega建议未认证基态omega外推。数值诊断见audit/omega-adiabatic-diagnosis-20261008/diagnostics.json。

2026-10-08 L9 omega0.3序参量双对数图更新：按20261007 order-current.pdf右图相同定义，重新读取两端八根目录原340seed扫描，85/85cell、340/340seed无排除重复，共1536610完整区块测量(前次1348180)，200complete链、50个4/4完成cell；各seed2020–6900测量。共同phase加权后取80非零允许q上最大Szz/Nc，独立seed删整链jackknife1SE且每次重新取峰；warmup0/no cut，运行点空心。旧快照新汇总公式复现原曲线残差低于1e-14，PNG/PDF已检查。新加密tau8批次为另批，未混入原五lambda曲线；未认证热化、自关联或投影收敛。新鲜冻结矩、覆盖、CSV及脚本见audit/omega0p3-order-loglog-current-20261008/。


2026-10-08 fast-1006 系统正确性验证：新增独立 Gaussian 自由声子、4/8 Majorana Fock、快速更新/缓存三个验证器和统一 runner；复用八个已有回归，本地11/11 PASS。HPC20015 workq冻结源码重编译，156891与数组156892/156893全部exit0；L2/4/6 local/bond六条自由MC均PASS，627个统计检查最大3.303SE、预设6SE门槛，数值失败0。独立有符号复权重误差2.15e-14、中心Green3.57e-15；864缓存事务含624接受/240拒绝，拒绝深层缓存精确不变，PBC/shear及24设置通过。有限投影自由静态响应按实际Mehler score为(1-exp(-omega_ph*tau/2))/k，未强制短tau等于1/k。本轮未发现生产实现错误，未改生产模型/API/容差；非等时/键四点是测试重建，未认证完整自旋投影或相互作用长链收敛。三条旧慢参考测试主动取消后按生产deferred模式重跑，旧记录保留。完整结果、源码/二进制SHA及运行方式见audit/pfqmc-correctness-20261008/REPORT.md与docs/PFQMC_CORRECTNESS.md。

2026-10-08 L9 omega0.3逐seed sweep图：重新读取tau8/10/12/14、lambda1.3/2/2.4/2.6/3的80条原扫描链，无排除/重复，共200140测量、每链2020–3380，均仍running。每seed独立累计共同phase原始Szz谱后在80个非零允许q取最大/Nc，每10sweep一点，累计从sweep1开始、warmup0/no cut。各tau五lambda分面PNG及显示从sweep100开始的PNG放大版已生成，准确seed编号见图例；累计前缀复现上一冻结快照残差低于1e-12，不能据此认证热化/自关联。图、CSV、覆盖与可复现脚本见audit/omega0p3-seed-sweeps-tau8-14-20261008/。用户要求以后仅输出PNG，已记入根AGENTS.md。

2026-10-08按用户要求覆盖上一逐seed sweep图：保持相同冻结文件长度/样本末点，80链全部丢弃sweep1–1000的原始phase分子分母，从1001重新累计Szz谱后取80非零允许q最大/Nc；首导出点1010。保留120140测量，每seed1020–2380。原tau8/10/12/14及zoom共8张PNG、CSV、README和metadata已覆盖，zoom显示从1100开始但累计仍1001开始。扣除后累计与全量减去丢弃矩直接复算残差8.24e-18；冻结原图snapshot保留。此用户指定热化cut不等于热化/自关联收敛认证。见audit/omega0p3-seed-sweeps-tau8-14-20261008/。

2026-10-08 L9 omega0.3序参量随tau双对数图热化cut更新：统一对原340seed、85cell的同一20261008冻结总矩扣除sweep1–1000原始phase分子/分母，再跨四seed共同重加权、取80非零允许q最大Szz/Nc；保留1196610测量，每seed1020–5900。全部85cell四seed覆盖完整，无排除重复；seed jackknife1SE每次重算峰。tau8/10/12/14的80条逐seed截断累计前缀独立核对最大残差7.81e-18，PNG双对数轴与误差棒已检查。只生成PNG，不认证热化/自关联/投影收敛。新图、CSV、冻结丢弃矩、原总矩hash与脚本见audit/omega0p3-order-loglog-cut1000-20261008/。

2026-10-08 L9 tau8极化率lambda密集图：合并20022作业86815原16seed与20015作业156874新增104seed，严格核对提交清单，lambda1.9–3.0步长0.1、12点各10seed，无缺失/重复/数值失败。全部仍running，旧链2780–2840样本、新链1010–1210；本次按用户要求合并可用完整区块，不切warmup，保留不等链长并按样本汇总共同phase原始矩后做chi连接减法，再取80非零允许q实部最大。误差为删除整seed jackknife1SE且每replicate重算最大；属于当前快照，不认证热化或收敛。仅PNG及CSV/冻结汇总/脚本/快照时间和SHA见audit/L9-tau8-chi-lambda-86815-156874-20261008/。SSH大数据停滞通过TCP_MAXSEG512代理恢复，生产作业未改动。

2026-10-08 tau8各seed极化率随sweep：按用户指定lambda2.4/2.5/2.9重新读取30seed，均running且无缺失/重复/数值失败；每seed从sweep1累计phase原始矩再连接减法、取80非零允许q的Re chi_zz最大，每10sweep一点、未切热化。各点末端1050–2860/1110–1140/1090–1210；三分面总PNG与三张独立PNG、CSV及冻结快照/时间/sha见audit/L9-tau8-chi-seed-sweeps-20261008/。

2026-10-08补充lambda1.3–1.8 tau8生产：按用户要求6点步长0.1各10seed共60链，20015workq阵列156895[0–59]全部R、无Q；节点预检156894两样本通过exit0，60/60每链小测试通过、生产argv匹配、当前数值失败0。L9/omega0.3/J1,1,.5/dt0.1/interval80/warmup0/sweeps10000，chi+fermions+energy开启，每链1CPU/数学库单线程、720h、不指定内存。使用冻结fast1006源d3ae252e…与binary87786fc8…；核对当前生产源码一致，仅构建文件新增测试目标。独立seed106081300–106081359，任务表lambda降序，独立目录和所有清单见audit/production-L9-omega0p3-chi-lambda1p3-1p8-tau8-20261008/。

2026-10-08 lambda2.8–3 step0.025加密：按用户要求跳过已有2.8/2.9/3.0各10seed，仅六新点2.825/2.85/2.875/2.925/2.95/2.975各10seed，共60链。20015workq156896[0–59]最初50R/10Q；用户明确撤回20022迁移要求，迁移未执行（仅准备包）。按用户要求取消上一批lambda1.5的156895[30–39]、seed106081330–339，qdel全部成功并确认10X，其余50R；随后156896全部60R/0Q。tau8/L9/omega0.3/sweeps10000/warmup0/chi+fermions+energy、冻结fast1006源和binary延续，各链1CPU/单线程/720h不指定内存。参数清单与覆盖记录见audit/production-L9-omega0p3-chi-lambda2p8-3-step0p025-tau8-20261008/；取消证据保留在lambda1p3-1p8原audit目录。

## 2026-10-08 L6 omega0.3 tau10/14 完整链热图

HPC20015 156856[]/156857[] 的40/40链均complete、numerical_failures=0，
每链15000测量，共600000；lambda1.3/2/2.4/2.6/3，每参数4seed。
共同相位合并原始矩后生成connected Szz/Nc与Re chi_zz各tau一张PNG，
每类跨tau共用色标，沿用参考有限盒|W|²重构；允许动量CSV保留独立seed jackknife1SE。
warmup0、不丢弃sweep；完成不表示热化、自关联或投影收敛。
见[audit及验证](../../audit/L6-omega0p3-order-chi-heatmaps-final-20261008/README.txt)。

2026-10-08按用户取消并新交L6：已逐一取消L9tau8的lambda1.8(156895[0–9])、2.5(156874[42–51])、2.2(156874[68–77])共30条链，qdel全部成功且确认30X/无活动匹配链；部分输出保留，其他lambda未取消。新L6/omega0.3/lambda2.85,2.9,2.95/tau10/dt0.1/interval100/sweeps20000/warmup0各10seed共30链，20015workq156899[0–29]全部R无Q。预检156898两样本chi+fermions通过exit0，30/30节点小测试通过、生产参数匹配、数值失败0。冻结fast1006生产源码与当前一致，binary87786fc8…，每链1CPU/单线程/720h、不指定内存；延续用户不迁20022要求。完整参数、冻结源/运行库、取消和启动证据见audit/production-L6-omega0p3-chi-lambda2p85-2p95-tau10-m20000-20261008/。

2026-10-08任务缩减与高lambda新增：按用户要求完整取消156895的lambda1.4/1.7共20链，剩余1.3/1.6各取消高编号6seed保留4seed，共12链；156896六点及156899三点每点取消高编号4seed，共36链。合计68条qdel成功，最终无活动匹配任务，156895剩8R、156896剩36R、156899剩18R。新L9/omega0.3/tau8/lambda3.1,3.2,3.4,3.6各5seed共20链，20015workq156903全部20R/0Q，预检156902两样本通过exit0、20/20节点小测试通过/数值失败0，10000sweeps/warmup0/chi+fermions+energy、冻结fast1006、1CPU单线程720h无内存请求。完成后北京时间21:46查询20015总760CPU占713/空47(workq0Q)，cn18空1/cn19空20/cn20空22/cn3空4。完整清单、取消/启动证据和资源快照见audit/production-L9-omega0p3-chi-lambda3p1-3p6-tau8-5seed-20261008/。

2026-10-08缩减156874并新增tau扫描：按用户要求156874当前10个活动lambda点各保留最低编号4seed，共40R，取消44条多余链全部成功并确认无活动取消目标；此前2.2/2.5取消保持不变，20022既有链未动。新L9/omega0.3/tau14,12,10/lambda1.3,2,2.4,2.6,3每组合新增6独立seed共90链，10000sweeps/warmup0/dt0.1/chi+fermions+energy，冻结fast1006生产源码与当前一致。最长tau先交：20015workq156905(tau14)、156906(tau12)、156907(tau10)各30R/0Q；预检156904两样本通过exit0，90/90每链节点小测试和argv核对通过/数值失败0。每链1CPU/数学库单线程/720h/no memory，保留用户不迁20022要求。所有输入、冻结源/运行库、保留seed及取消/提交/启动证据见audit/production-L9-omega0p3-chi-tau10-14-add6-20261008/。

2026-10-08 omega0.3大lambda序参量混合诊断：以最近裸狄拉克点差的允许q=(5/9,4/9)为主，从此前约19:00北京时冻结CSV上限只读重抽原tau8/14、lambda1.3/2.4/3各4seed共24链的10-sweep原始区块。丢弃1-1000后，各tau统一窗口1001-2600/1001-2000；rank/fold split Rhat在大lambda的m2为1.007-1.029，复D实/虚较大值为1.046-1.181；200对10-sweep分箱SE倍率中位1.47-1.92。支持相关采样和复序参量慢混合，不证明结构性非遍历、精确tau_int或物理平衡；tau14仅5个200箱，Rhat近1也不是认证，相位锁定本身不改变m2。跨tau链长不同，另存共同1001-2000窗口诊断。与此前截断谱逐q复核残差0、反演6.58e-18、phase偏差1.56e-15、零排除重复；补充新鲜tau8 lambda2.9短链另存不混入主表。固定q累计、复平面及ACF/分箱PNG已检查，源码显示free整键路径提议、空间仍单键更新，作为慢混合机制推断；未变更生产任务。数值、冻结区块、来源及复现脚本见audit/omega0p3-highlambda-mixing-diagnosis-20261008/README.txt。

2026-10-08采样收敛优化研究：源码确认单键bond-refresh已是自由Gaussian pCN的beta1独立提议，且每提议full rebuild全部时间节点。提出第一优先级为保留整键更新、少量混入q=(5/9,4/9)实共轭模式集体更新；第二为不改物理trial/force的参考均值移位并加入精确Hastings修正，尚未实现生产C++。本机实际运行解析受力Gaussian surrogate：tau8/14每方法4链、warmup2000+measure20000，原m2相关因子中位4.64/15.86，移位1.02/1.01，q混合1.00/1.01；beta0.5虽将接受率升至77.3%/67.2%，m2相关因子11.31/13.54，说明高接受率不必改善混合。模型无真实matter/Pfaffian、不能作为生产加速倍率；32toy链解析m2核对通过。另用相同Mehler作用量60例检查移位/集体一般MH与简化公式，残差4.39e-13/1.54e-12。方案、PNG、脚本/数值及真实代码验证计划见audit/sampling-optimization-research-20261008/REPORT.txt；未提交或改变HPC任务。

## 2026-10-08 88016/88017 现有数据检查

截至主机UTC 2026-10-08T14:29:18.109286Z，64条已启动链全部仍预热，
最近记录warmup400–800/1000、完整测量区块0；另88017[64/65]排队。
数值失败和未认证拒绝均0；没有m/chi/Binder可分析结果。
dt=.05链最近100sweep约229秒/sweep，按预热速度外推11000sweep约700h，
尚未包含测量开销，不能承诺在720h内完成。未改任务。
见[现有数据检查](../../audit/vbs-field-88016-88017-analysis-20261008/REPORT.txt)。

## 2026-10-08 按用户要求停止88016/88017

20022的88016[]和88017[]已qdel成功；复查qstat无这两个数组的子任务，退出清理完成。
88018[]保留，释放资源后12条链开始运行。原始文件和分析快照保留。
见[取消记录](../../audit/cancel-88016-88017-20261008/cancellation.json)。

## 2026-10-08 86814/86815/86818现有测量分析

主机UTC14:52–14:53直接收集20022三个数组的48条有效链（四seed/cell），
共223200测量、全部running、数值失败0。86814仅tau10/lambda1.3和2；
86815 tau8、86818 tau4.2各含lambda1.3/2/2.4/2.6/3。
原始矩共同相位合并后生成固定q=(5/9,4/9)及非零动量最大chi与Szz_conn/Nc，
并保留独立seed jackknife、100–1000sweep组诊断、前后半段与discard1000敏感性。
tau4.2到8响应明显增大；tau8/10仍有链差和窗口漂移，当前不能认证投影收敛。
tau8 chi_max最大采样点在lambda3右端；固定q在lambda2.6略高不足定位临界。
未修改、取消或提交任务。PNG图及证据见
[分析报告](../../audit/L9-86814-86815-86818-analysis-20261008/REPORT.txt)。


2026-10-08 固定q集体更新：可选实cos/sin空间pCN，默认关闭；保留实际freeMean、trial/force、Hamiltonian、phase重加权及同候选/u一次参考重试。新q验证器与既有mc回归本地通过，local L2/bond L3与先前冻结源码同工具链重编译的blocks逐字节一致、最终路径/全部观测量/计数一致；诊断开启不改变blocks。新增逐sweep原始phase、D/m2、energy/Xz、CPU时间与分开q/base计数，见ADR-014。冻结源manifest c7c3f1fe…；20015节点构建/预检156922.mgt尚未核对完成。按用户要求直接提交正式数组156923[].mgt，L9/omega0.3/lambda2.4/tau14，warmup1000+measurement10000，beta1/every1 q=(5/9,4/9)与关闭集体更新各4独立seed，共8链；同一冻结新程序。1CPU/链、BLAS/OMP1、720h、无mem，不迁移20022。提交快照8Q/父H，依赖afterok:156922.mgt，预检成功后PBS自动放行。用户要求提交后停止，本轮未等待预检结束或收集生产结果，不宣称加速/收敛；记录见audit/qmode-pilot-L9-lambda2p4-tau14-20261008/submission.json。

2026-10-09：取消156874[]的lambda2/2.4/2.6/3中16条运行链（另8条此前已退出）；新提交L9 omega0.3 tau8 lambda4/4.5/5各5seed，20015 workq 156926[]，10000sweep、warmup0、当前冻结fast1006源码（18edf30c），复查15R。审计目录audit/1009/production-L9-omega0p3-chi-lambda4-4p5-5-tau8-5seed-20261009。未来audit按MMDD归档。

2026-10-09 L6集体更新对照新增及sweep修订：沿用omega0.3/lambda2.4/tau14，q=(1/2,1/2)、自共轭单实模式、beta1/every1，与关闭集体更新各4独立seed。复用156923已验证冻结fast1006 source manifest c7c3f1fe…/binary a8b9c10b…，节点156941完成L6/M140两组2sweep预检exit0。原10000sweep数组156942刚启动后，按用户追加30000sweep要求取消，8/8确认F；保留原seed106099000–7在新独立目录重交156943[].mgt，warmup仍1000、measurement30000、measure1/block10/chi+fermions+energy。20015workq、1CPU/链/数学库1线程、720h/no memory、不迁20022；提交快照8Q，无重复活动原链。按用户要求提交后停止，未收集/分析生产结果，不认证收敛或加速。新audit：audit/1009/qmode-pilot-L6-lambda2p4-tau14-m30000/submission.json，原准备与取消记录另保留在同级qmode-pilot-L6-lambda2p4-tau14。

## 2026-10-09 fast-1006 代码优化与 L6 验收

用户授权优化并指定L6/20015。CLI预热raw路径诊断移至onSweep，关闭时也不再
触发无用的signed phase/Green；sample其他调用者和正式测量认证保持。新增
kitaev_ssh_production_profile复用完整生产CLI，并区分整键/q候选计时，普通程序
无新增计时器开销。按L6 profiling证据，再将已认证整步dense矩阵直接移入
不可变owner，保留lvalue/generic行为、矩阵系数、检查、缓存身份与数值门槛。

20015/workq首批156945/156946/156947共8子任务，最终156963/156964共5子任务
均exit0、PASS、结束；覆盖L6/M140/q开关、预热诊断、shear1/local、noGreen及
profiler对照。配对路径与接受/拒绝计数一致，全部原始blocks逐字节一致、数值
失败0、比较的物理JSON字段最大缩放差0。最终构建156962 exit0，既有缓存
事务验证624接受/240拒绝通过，log ratio误差2.49e-14。首构建156944仅末尾
比较器误要求signed工作量不变而退出1，修正计数口径后smoke通过，原日志保留。

L6/omega.3/lambda2.4/M140/dt.1、chi+fermions+energy、1CPU/数学库1线程。
首批4预热+4测量正反顺序配对：预热wall节省中位2.64%，整段CPU节省中位0.78%。
最终warmup0/4测量四次正反配对CPU节省3.97%–5.19%，中位4.46%；这是短工程
测速，未认证长链/L9加速或等误差采样收益。profiling仍约65%用于传播矩阵
生成/存储，后续优先研究半步矩阵按需构造。既有生产作业未改动。

冻结最终源码manifest c82f83a6…，HPC binary3b10b149…；本机新构建
build-code-speed-20261009/kitaev_ssh_mc（c5db92c7…），本机只验证构建/--help。
源码、diff、全部参数/CPU计时/调度/二进制身份见
[本轮审计报告](../../audit/1009/code-speed-L6-20261009/REPORT.txt)，设计见ADR-015。

## 2026-10-10 半步矩阵按需构造验收（10-09启动）

用户继续授权优化，沿用L6/20015。延迟phase且合并节点传播的提议
不再提前生成内部canonical半步；快照只保存完整认证的weight_order。
两端半步保留，正式测量继续原完整signed/Green重建。测量及uncached
恢复后继续局部更新可安全回到空canonical表示；指数溢出仍按需生成
并验证半步，再尝试原h*h。RNG、接受/相位重加权/测量定义、所有
实际传播矩阵finite/Hermitian门槛及同候选/u一次重试中止策略保持。

20015/workq构建156965 exit0；数组156966/156967共5个L6子任务全部
F/exit0/PASS。覆盖M140 q开启/关闭、shear1局部更新及预热诊断、
noGreen和完整production profiler。最终路径和接受/拒绝计数精确一致，
所有raw blocks逐字节一致，物理JSON/CSV最大缩放差0，数值失败0。
新增L6检查覆盖lazy/canonical系数、普通/二分谱分解、拒绝深值不变、
测量/局部/恢复/集体状态切换及不完整输入；既有缓存624接受/240拒绝
和失败恢复检查通过。既有生产任务保留版本，未修改。

相对上一轮3b10b149程序，L6/omega.3/lambda2.4/M140/dt.1，warmup0、
4测量sweep、chi+fermions+energy，四次正反顺序配对CPU减少
8.55%–12.16%，中位10.48%。局部M4和noGreen M16用例也通过。
M140/q0 profiler半步构造1430次，旧逐节点方式对应61617次；嵌套
计时不可相加。短测速不认证L9/长链/收敛或等误差采样收益。

冻结源码manifest f4689c4f…、HPC binary2fd634ae…；本机
build-lazy-half-20261009/kitaev_ssh_mc（2a56ce40…）构建/--help通过，
同目录run_mc.sh设置本机MKL和单线程；数值验证全部在20015。
本任务按启动日保留在audit/1009/lazy-half-L6-20261009，记录见
[验收报告](../../audit/1009/lazy-half-L6-20261009/REPORT.txt)，设计见ADR-016。


## 2026-10-10 L6 半步优化长链对照已提交并停止

用户要求长链测试、提交HPC后停止。20015/workq新交156968[]（q关闭，
索引0–3）和156969[]（q开启，索引4–7），共8个配对任务、16条链。
L6/shear0/omega.3/lambda2.4/tau14/dt.1/M140，bond-refresh，预热1000、
测量10000、measure1/block10、chi+fermions+energy；q=(1/2,1/2)，
beta1/every1或0。每组4个独立seed106110000–7；同seed新旧程序顺序
配对，按replica交替先运行版本。预计较长q0组先提交，1CPU/任务、
数学库1线程、720h、无mem，现有生产任务未改。

复用上一轮已验收冻结fast-1006程序：新source f4689c4f…/binary2fd634ae…，
显式对照source c82f83a6…/binary3b10b149…；新目录保存两份source归档与
二进制，提交前校验源码归档、程序、库及ldd。本轮不重编译、不重跑
既有数值验收。配对结束由任务自动核对路径/计数/raw blocks和物理
JSON/CSV，并记录CPU耗时。提交快照8R；按用户要求停止，未等待、
收集或分析长链结果，不认证长链一致性/速度/收敛。
记录见[提交审计](../../audit/1010/lazy-half-long-L6-20261010/submission.json)。


## 2026-10-10 156943 L6 q集体更新中期对照

20015主机UTC00:02快照，8/8链R、数值失败0；测量2550–2609/30000。
只读收集冻结c7c3f1fe…/a8b9c10b…批次，q=(1/2,1/2)开启/关闭各4seed，
统一前2540测量sweep，逐sweep phase/m2/ReD/Xz/能量与raw blocks核对PASS。
m2相位重加权组值0.004801±0.000148与0.004956±0.000121（4seed jackknife）。
m2单链tau_int中位3.48→0.78sweep，ReD 36.05→0.87；测量CPU/sweep
9.277→9.343s（+0.71%），m2 CPU/暂定ESS中位64.34→14.60s，
短链成本比约4.41，额外cut500/1000约4.75/4.15。新增q接受85.55%，
base两组均约14.9%。固定q去相关有初步收益，不认证长链等误差加速。
connected chi组值仍1.036与1.194，累计值漂移，Xz/能量仍有慢变量；
不能认证所有物理量或seed收敛。原任务未改；未混入后续代码优化批次。
数据快照、窗口/重分块、每seed与PNG见
[中期报告](../../audit/1010/qmode-L6-156943-interim/REPORT.txt)。

## 2026-10-10 已验证优化默认开启的用户偏好

用户要求以后集体更新等适用优化在正确性验证后默认开启。已写入根目录
和fast-1006的AGENTS.md、交接及ADR-014；未来新生产提交显式q-every1/
beta1与适配几何/目标的非零允许动量，使用已验收冻结新版及生产加速配置。
用户明确设置和基线对照优先，已有任务保持原版本/开关。CLI省略参数的
q默认仍0，未来提交脚本须显式开启；本轮未修改计算源码或提交新任务。


## 2026-10-10 半步优化长链中期验收检查

用户要求核对长链是否成功、是否可提交新修改，并确认集体更新已默认
用于新生产提交。本轮只读检查20015的156968/156969；远端UTC00:20
快照8/8配对任务R，均仍跑第一版本，测量2281–2803/10000，完整链
0/16、完整配对0/8。最新累计acceptance/q计数中的first failures、
reference recovered/failed、untrusted rejected均0，numerical事件文件
不存在，stderr只有正常启动说明。running结果JSON不代表已完成验收。

半步按需构造的长链输出一致性和CPU收益仍待完整同seed配对，不能
标注长链通过，不能用不同seed的首条链耗时推断配对加速。此前针对性
正确性已通过，当前核心源码仍与f4689c4f冻结验收版一致（只有AGENTS
生产偏好文档更新）；可以纳入后续新生产提交，说明保留“长链对照
进行中”。集体q更新继续默认显式every1/beta1和适配动量；q0为专门
对照组。未修改核心代码、提交/取消新HPC任务或执行Git commit/push。
记录见[中期审查](../../audit/1010/lazy-half-long-review-20261010/REPORT.txt)。

2026-10-10 L9/omega0.3序参量随tau当前图：只读重采原340独立seed扫描，
lambda1.3/2/2.4/2.6/3，各17tau(0.2–14)、4seed，无缺失/重复/数值失败。
两端HPC快照时间分别北京时间08:55:16–08:55:50(20015)、
08:49:41–08:51:09(20022)，220/340链、55/85cell完成，保留1829760样本。
沿用OBS03 m_z²=max非零允许q Re(Szz)/81；先joint phase合并，再取最大，
独立seed删除jackknife1SE且每副本重选峰。warmup0、未额外切样本；
完成链用最终JSON、部分链只读文件打开时字节上限内完整CSV区块。
新增tau8/长tau补充及qmode试验不混入原扫描；未认证热化、自相关或投影收敛。
85点与保留80非零q网格最大值一致，反演残差1.30e-18；PNG已视觉核查。
线性/双对数PNG、CSV、时间/SHA与复现脚本见
[audit](../../audit/1010/L9-omega0p3-order-vs-tau-current/README.txt)。未改生产任务。

2026-10-10 tau12/14序参量随sweep：只读重采原tau扫描40独立seed，
lambda1.3/2/2.4/2.6/3每tau各4seed，均running，无缺失/重复/数值失败。
20015快照北京时间08:58:41–08:59:00；tau12保留3650–3960、tau14保留3220–3350
测量sweep。从测量sweep1按每10样本完整块逐seed累计共同phase原始Szz矩，
再取80非零q最大并除81；各seed独立显示，warmup0且无额外cut，未估计轨迹误差。
固定q49累计值也保留CSV，补充add6/qmode链未混入。连续完整块及40seed覆盖检查
通过，两张PNG视觉核查完成；不据累计曲线认证热化或收敛。见
[audit](../../audit/1010/L9-tau12-14-order-vs-sweep/README.txt)。未修改生产任务。

2026-10-10按用户要求对原L9/omega0.3 tau扫描的tau>6截去前1000测量sweep：
适用tau6.4/8/10/12/14，共100链；保持08:49–08:55原快照各seed样本上限，
两端只读重采完整CSV块并验证截断前phase/Szz全部原始和与冻结快照一致（残差0）。
保留1729760样本，tau<=6的60cell主值/误差/峰q/样本数完全不变。
每seed从sweep1001重新累计phase原始矩；随tau图合并4seed后重取80非零q最大，
重新做删除seed jackknife1SE；tau12/14逐seed轨迹从1001累计并保留原sweep横坐标。
轨迹末端和截断汇总序参量差<=8.68e-19；线性/双对数及两张轨迹PNG已核查。
此为用户指定热化截断，不认证剩余数据热化或收敛；补充批次仍不混入，生产任务未改。
结果、CSV/冻结快照与复现脚本见
[audit](../../audit/1010/L9-order-cut1000-tau-gt6/README.txt)。

2026-10-10双对数序参量图加入后续生产补充：只读重采20015六批349链
（原tau8加密104、低lambda60、2.825–2.975加密60、高lambda3.1–3.6共20、
tau10/12/14 add6共90、lambda4/4.5/5共15）。新HPC快照北京时间09:06:37–09:07:52，
263链在截去1–1000 sweep后保有样本，86链raw_count<=1000排除；无其他收集失败。
原340链保持此前冻结样本上限/截断口径；同lambda,tau在每seed截断后合并共同phase
原始和再取80非零q最大，重算删除独立seed jackknife1SE，不平均逐seed峰值。
总603有效独立链、107cell/27lambda、1890830保留样本；其他lambda单点加在tau8，
主图为tau双对数，右图展开tau8全部lambda便于辨别。65无补充cell与原截断表一致，
107点逐q最大值核查通过、非零q反演一致；PNG视觉核查完成。
包括已停止但有截断后样本的部分生产链，未认证热化/收敛；qmode试验未混入。
fermions旗标为附加观测量，未移除matter权重；原/补充物理参数核验相同。
早前dense tau10计划已改tau8且未生产提交，现场该根无生产结果，不加入虚构点。
来源/seed覆盖/排除/CSV/SHA及复现脚本见
[audit](../../audit/1010/L9-order-loglog-supplemented-cut1000/README.txt)。未修改生产任务。

2026-10-10按用户要求tau8极化率随lambda图改为前800测量sweep截断：
原20seed+后续tau8五批259seed共279候选链，保留此前各链冻结样本上限，
203链有截断后样本、76链raw_count<=800排除，无其他收集失败；27lambda，
4–10seed/点，共251470保留样本。重抽ChiD/ChiIconj/ChiRaw的共同phase原始矩，
每seed去除1–800后同lambda合并phase比值，chi=Raw-81*D*Iconj，再取80非零q
Re chi最大；删除独立seed jackknife1SE重算连接减法和峰。无自由声子减法。
原tau扫描和补充批次同口径合并；对照其他tau图未修改截断规则。
27点逐q最大值核验通过，独立相位比值重构差2.84e-14；203seed及800样本
逐链截断核验通过。全lambda/1.3–3.1放大PNG视觉核查完成；此为指定截断，
不认证热化、相关采样或收敛，包含有保留样本的停止部分链。见
[audit](../../audit/1010/L9-tau8-chi-vs-lambda-cut800/README.txt)。未改生产任务。
