# 有限虚时间驱动与有限时间标度

状态：**方法依据 + 向 SSH–Kitaev 映射的设计草案，不是生产实现说明或成功结果。**
来源：[FTS](SOURCES.md#fts) 讲义第2、3、5、8页；[C0](SOURCES.md#c0) 通用驱动/声子接口；[DISCUSS](SOURCES.md#discuss)。


## 证据边界

- **论文事实**：FTS 讲义及 SSH–Kitaev 论文提供方法背景和物理动机。
- **数学推导**：本页的协议、双支轮廓和标度公式是条件推导。
- **源码现状**：通用驱动头文件和当前 SSH 入口的描述来自源码核对。
- **报告结果**：目前没有 SSH–Kitaev 线性驱动生产结果；历史报告不等于本轮验证。
- **未验证计划**：控制参数、初态、离散规则和驱动收益均待实现与验证。


## DRV-01 已有方法材料实际研究什么

FTS 讲义研究方格上的单 Dirac 费米子模型、U=U0±Rτ 的虚时间驱动和无序初态。第3页依次提出：用 m²(2R)/m²(R) 初定位临界点；用无量纲关联量的 RL^r 标度确定 r；固定 RL^r 再确定临界点和 ν，用 r=z+1/ν 求 z。

讲义的模型、初态、铁磁观测量、具体 Uc/指数/尺寸均不是 SSH–Kitaev 结果。本文件保留这一路线，并将其应用于本项目的部分明确标为草案。没有重新验证讲义引用的其他论文。

## DRV-02 拟议的驱动对象与双支态

为避免 g 同时表示微观耦合和到临界点距离，控制参数记 p；首个候选是 p=g，是否最终选择须确认。候选升驱动 R>0：

$$
p(\tau)=p_i+R\tau,\quad0\le\tau\le\tau_f,\quad
\tau_f=(p_f-p_i)/R.
$$

下降协议需使用一致的有符号速率/正长度约定；p_i=p_f 的常量协议单独给定长度，不能用上式的 0/0 定义。

拟计算的态为

$$
|\psi_R\rangle=U_R|\psi_i\rangle,\quad
U_R=\mathcal T_\tau\exp[-\int_0^{\tau_f}H(p(\tau))d\tau],
\quad\langle O\rangle_R=\frac{\langle\psi_i|U_R^\dagger O U_R|\psi_i\rangle}
{\langle\psi_i|U_R^\dagger U_R|\psi_i\rangle}.
$$

在对称双支、中心测量约定下，总轮廓长度 2τf，确定性参数沿轮廓是 p_i→p_f→p_i，可表示为 p_contour(s)=p_i+R min(s,2τf−s)。这是本次讨论的拟议映射，不是 FTS 讲义规定的 SSH 代码接口。

**参数折返不等于声子构型镜像。** 当前优先评估保留一个总开放轮廓：两侧的内部节点是不同积分变量，中心节点是唯一共享坐标，两外端各有 trial 因子。一个数组即可保存 `X_0,...,X_{2N}`；这不要求两个完整 `PhononPath` 对象，也不表示两侧在耦合后的分布中统计独立。逐构型正权重不由总范数正性推出；不能强制相同场来消掉相位。

对声子坐标对角的 `X`、`XX` 型观测量，中心测量必须使用这个唯一共享坐标和同一 matter cut。若直接相乘两个现成开放路径权重，会重复中心 trial 因子；双对象实现只有在共享中心、单次中心积分和正确端点条件都明确后才与总轮廓等价。

单支 τf、总 tau_total、Δτ 和 Monte Carlo sweeps 必须分别记录。常量 tau_total=0.2 的旧数据不是 R=某值的线性驱动数据。[H, §1.4；R-SAMPLING；DISCUSS]

## DRV-03 初态与控制参数待定

g_i=0 的固定扇区自由 Majorana 态 × 自由声子真空是讨论中的优先候选，仍应核对有限系统零模、奇偶及初态实现。不把任意 Gaussian 试探态当作非零耦合的精确基态。

若驱动 g，λ=g²/k 一般非线性；若驱动 λ，必须重新定义其速率及 g_j，不混称 R。改变 m/ω 则会影响自由声子核和提议，超出当前“只改变 g”的最小方案。

需要验证本初态与 VBS 观测量是否满足所选 FTS 形式；“都叫无序”不足以证明与讲义初态具有同样标度。

## DRV-04 当前代码可复用与缺口

C0 中 `PureImaginaryTimeProtocol` 有 constant/suddenQuench/linearRamp 和分片中点值。`pureProjectorIndependentBraKetContour` 有独立 HS bra/ket 的反序伴随逻辑。

当前 HEAD 已现场核对：`MajoranaGenerator` 仍只接收 X，`honeycombGenerator` 捕获固定模型，`nodeHalfStepInto` 使用 `path.node(node)` 和 `exp(-dt K/2)`；模型新增的显式 `g` 重载只用于 validation，`coupled_mc.cpp` 仍没有有限速率驱动 CLI。当前树已有 [`reproduction/pure_state_driven_stage1/`](../reproduction/pure_state_driven_stage1) 和 `pure_state_driven_qmc` 构建目标，但该基线仍未证明适配 SSH 声子节点时间、双支总轮廓与唯一中心测量。[CODE_MAP.md](CODE_MAP.md)

当前已增加 `KitaevSSHModel::frameworkGeneratorWith(X,g)` 的显式 `g` 组装重载，并在 [中点 canonical 回归](../audit/midpoint-driven-regression-20260924) 中以 9 组小系统协议与 Fock 有序乘积核对。这个接口和回归只属于适配前的验证层；`MajoranaGenerator`、采样器、缓存和生产 CLI 仍未获得时间协议上下文。

以上是 C0 源码快照事实，不等于对用户最新工作树的完整审查。[CODE_MAP.md](CODE_MAP.md)

## DRV-05 最小实现约束：尚待接入

- 明确将节点/分片时间和协议值传入生成元，或采用等价不可变上下文。不能按调用次数计时。
- 节点声子半步与通用协议中点值不是同一数组；选定离散规则并用小系统验证后才能接入。
- 若采用区间中点 `g_{j+1/2}`，同一内部节点的左右半步一般分别是 `exp[-dt K(X_j,g_{j-1/2})/2]` 与 `exp[-dt K(X_j,g_{j+1/2})/2]`；它们通常不相等，不能继续按固定 g 把一个 node half-step 平方，也不能直接复用旧 node 谱/Cayley 缓存。
- 第一版可先使用未合并、canonical、按 `(node, adjacent interval, protocol)` 生成的参考路径；node-merge、prefix/Cayley 和恢复后的优化 bookkeeping 逐项证明前应明确禁用或拒绝，而不是静默沿用。
- 中心能量等显含耦合的测量用终点参数；声子测量保持相同终点语义。
- 权重、缓存、局部候选、参考重算和诊断共享同一协议。不可恢复成固定 g 参考问题。
- 保留常量入口及其回归。离散后保存实际 R，不把近似 R 的数组标成严格 R 与 2R 配对。

这里规定的是正确性条件，不锁死类名、CLI 参数或缓存数据结构；这些细节留给完整工作树的专项核对。[ADR-004](decisions/ADR-004-driven-protocol.md)

## DRV-06 本次复核的路径修正

当前 `PhononPath` 已能表示总长度 `tau_total=2*tau_f` 的开放轮廓：`M+1` 个节点、唯一 `center()`、`action()` 只在两个外端加入 trial。保留该总轮廓并附加分支视图/协议表是优先方案；“两个完整路径对象”不是物理必需条件。

现有固定 g 的 `pathFactors`、`prepareOwnedNode`、`prepareOwnedPath`、`localMergedNode` 和恢复重建都利用同一节点半步在相邻位置重复出现并可平方。在线性中点 g 下，该复用一般改变离散对象；直到非对易小系统验证通过前，驱动参考必须保留相邻区间身份并禁用这些不适用优化。该结论是本次源码复核与条件推导，不是已完成实现。

## FTS-01 讲义形式与本项目符号映射

讲义 p2 中 g=U−Uc；本文件把到临界点的距离记 δ=p_f−p_c，序参量指数记 β_op。对讲义给出的从无序初态出发的形式，写为

$$
m^2(\delta,R,L)=L^{-d}R^{\alpha}
 f(\delta R^{-1/(\nu r)},L^{-1}R^{-1/r}),\quad
\alpha=\frac{2\beta_{\rm op}/\nu-d}{r},\quad r=z+1/\nu.
$$

这是符号转写。将左侧替换为本项目 m_z²=S_T,c;zz(Q)/Nc 是候选应用，尚无 SSH–Kitaev 驱动数据验证。

讲义称“大速率下”忽略第二变量。此前讨论把它解释为驱动长度 ξ_R~R^(−1/r) 相对 L 较小，同时还需检查是否已经离开纯微观瞬态；不能解释为 R 任意增大都有效。[FTS p2；DISCUSS]

## FTS-02 速率比：由上式直接相除

忽略第二变量并令 x=δ R^(−1/(νr))，有

$$
\frac{m^2(\delta,2R,L)}{m^2(\delta,R,L)}
=2^\alpha\frac{f(2^{-1/(\nu r)}x,0)}{f(x,0)}.
$$

在 δ=0、相应前提成立时，共同值为 2^α。这是基于讲义形式的代数展开，不是已经观测到的交点。

速率对必须使用同一 L、同一初态定义、同一终点、相同 Q 和通道规则。均匀背景在各个驱动系综分别估计但规则相同。不得按每条链样本最高峰改变 Q，也不从数据中任意扣自由声子参考来促成拟合。

## FTS-03 无量纲量与后续标度

讲义 p5 所称“关联长度比”的公式是 1−S(Δk)/S(0)。本项目候选用 OBS-05 的中心 Q 和正负邻点平均；名称不能代替公式。

候选形式是 $\mathcal R_i=\mathcal F_i(\delta L^{1/\nu},RL^r)$。在临界附近研究 RL^r，再固定 RL^r 改进 pc、ν，并计算 z=r−1/ν。非公度窗口、Q 随参数变化、形状与初态都是迁移时要检查的限制，不在此宣称已消除。

这里不沿用讲义的 Uc、r、ηb、ηf 数值。若检验费米子反常维数，需使用本模型合适的规范弦相关定义，而非把 SLAC 局域电子量机械替换。[FTS p3/5/8；P1 p5]

## FTS-04 方法验收不是每个终点先到基态

对每个有限 (R,τf) 协议要充分采样并检查离散误差；不以终点长投影基态收敛作为第一步驱动研究的前置条件。

“可用标度区”和“相位仍可承受的区间”是否重叠是要验证的问题。需要记录相位、统计精度和墙钟成本；比较的是提取同类临界信息的可信度与成本，不是要求驱动态结构因子等于平衡基态结果。

截至本包，没有 SSH–Kitaev 线性驱动生产结果、速率交点或符号优势的数值证据。下一项具体任务和运行预算仍需用户授权。
