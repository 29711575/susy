# VBS 观测量及统计定义

状态：已有定义整理；Binder 与物理 ξ/L 的扩展不是已完成的生产结果。
来源：[H](SOURCES.md#h) 第7–11节、[P1/P2](SOURCES.md)、[C0](SOURCES.md#c0)、[C1](SOURCES.md#c1) 及 [U-CORRECTION](SOURCES.md#u-correction)。


## 证据边界

- **论文事实**：VBS 候选、场展开和物理解释来自 P1/P2。
- **数学推导**：Fourier 窗口、connected/projected、矩阵和关联比公式。
- **源码现状**：导出器、窗口函数和统计接口按当前代码核对。
- **报告结果**：数值表和历史测试只代表其 audit 或报告。
- **未验证计划**：Binder、物理 ξ/L 和驱动标度观测量仍未认证。


## OBS-01 一个波矢不等于一个序参量

以 QD=K+−K−≡2K̃ mod G 为候选波矢，用三个键场的 Fourier 振幅定义序参量。区分 Q_bare、Q_uniform、最近允许 q、静态最软 q 与 QMC 峰。P2 §B 说明 Dirac 点会随耦合移动；不能把裸 Q 固定为所有参数和速率的真实峰。

采用列向量 F=(Fx,Fy,Fz)^T：

$$
F_\mu(Q;X)=\sum_{r\in\Lambda_L}e^{-iQ\cdot r}X_\mu(r),\quad
\Delta(Q;X)=F(Q;X)/N_c,\quad
S_{\rm dir}(Q;X)=F(Q;X)F(Q;X)^\dagger/N_c.
$$

二阶结构因子除以 Nc 后对应序参量二阶矩。有限对称系综不以 ⟨Δ(Q)⟩ 是否非零作为唯一有序判据。固定 z 键通道是合法选择；临界声子是三键模式的线性组合，不是三个独立临界复玻色场。[P1, p4 场展开；P2, §B；H, §7]

## OBS-02 非网格测量与平移平均

有限允许网格 K_L 上定义

$$
W_L(p)=\frac1{N_c}\sum_{r\in\Lambda_L}e^{-ip\cdot r},\qquad
w_k(Q)=|W_L(Q-k)|^2,\qquad\sum_kw_k=1.
$$

逐构型一阶恒等式：F(Q)=Σk W_L(Q−k)F(k)。固定原点的二阶外积含 k≠k′ 交叉项；不能逐构型直接扔掉。

定义 T_aX 为同一有限环面上的空间平移，则逐构型严格有

$$
S_T(Q;X)=\frac1{N_c}\sum_aS_{\rm dir}(Q;T_aX)
=\sum_kw_k(Q)\,F(k;X)F(k;X)^\dagger/N_c.
$$

后文 S_T(Q)=⟨S_T(Q;X)⟩。把它视为固定原点测量的同一期望，需要实际 Hamiltonian、边界、初态/试探态密度矩阵和规范处理保持相应平移对称性。平移平均不增加独立样本数。[H, §7–8]

完整网格二阶矩可重建任意 Q 的平移平均二阶矩，不会给有限 Hamiltonian 增加新动量模式。路径加密是展示和估计器重组，不是改善物理动量分辨率。

## OBS-03 raw、connected、projected

主要报告定义为

$$
S_{T,c}(Q)=S_T(Q)-N_c|W_L(Q)|^2
\langle\Delta(0)\rangle\langle\Delta(0)\rangle^\dagger.
$$

逐构型去掉整个均匀模式后再平移平均：

$$
S_{T,\rm prj}(Q)=S_T(Q)-|W_L(Q)|^2S_{\rm raw}(0).
$$

两者之差：S_T,c−S_T,prj=|W_L(Q)|² S_c(0)。前者只扣均匀系综均值，后者还删除零模涨落。[H, §9]

非零允许 k 的 W_L(k)=0，因此这里三种平移平均定义重合；这不是非网格点也无需扣背景的理由。`W_L(Q)=0` 本身亦不证明 Q 在允许网格上。

q=0 及相差倒格矢的等价点不进入非零动量绘图和峰搜索，但零模一阶、二阶统计保留。u+v=1 路径的 (0,1)、(1,0) 两端都等价于零动量。不得删掉用于背景重构的原始零模数据。

不要用有限样本所有非零网格均值的外积，悄悄替代上述平移平均 connected 口径；一般固定原点 connected 与此在有限样本上不是相同估计器。[H, §9；U-CORRECTION]

## OBS-04 标量通道与归一化

输出 Szz、Tr S、λmax(S)，以及每个量除以 Nc 的版本。Szz 是 z 键声子通道，不是自旋 SzSz。

$$
m_z^2(Q)=S_{T,c;zz}(Q)/N_c.
$$

先平均完整 3×3 矩阵并扣背景，再求特征值。不同动量各自的 λmax 可能对应不同模式；需要固定通道比较时优先明确 z 键或事先规定的同一个线性组合。不能把两 seed 的 λmax 平均成“合并结果”。[H, §10]

峰搜索输出要保留代表规则、对称/共轭集合及统计限制。驱动 R 与 2R 比较应使用相同终点、几何和同一 Q/通道规则，不按每条链的样本最高峰任意换观测量。[DISCUSS]

## OBS-05 无量纲关联比

在普通 L×L 周期几何中，令 C_L(q)=S_T,c;zz(q)，定义

$$
\mathcal R_i^{\rm win}(Q)=1-
\frac{C_L(Q+b_i/L)+C_L(Q-b_i/L)}{2C_L(Q)},\qquad i=1,2.
$$

这是项目现有窗口关联比的定义。中心和邻点共用原始区块及相位分母，并联合 jackknife。分母不能可靠估计时标记不可用；合法负值保留。[H, §11.1]

“差一个 2π/L”在 q·a_i 的无量纲坐标中成立；分数坐标位移为 (1/L,0) 或 (0,1/L)。笛卡尔方向和剪切几何改用 PHYSICS_CONVENTIONS 的实际倒格矢。

沿 u+v=1 的位移是 (1/L,−1/L)。当前 L=6 时 Q+−(1/6,−1/6)≈(0.413764,0.586236)，接近 Q−。所以这条路径上的一个网格间隔未必是在“离开”共轭双峰，不能机械解释为单峰宽度。[DISCUSS，由 H §3 坐标直接计算]

本定义不自动等于 ξ/L。非网格窗口、有限几何和中心失配仍影响尺寸函数；L=6/9/12 的相对网格偏移不同。文献中的二阶矩 ξ/L 需要其自身小动量、长度单位及峰形约定，不在本项目中自动套用。[H, §11；FTS, p5]

## OBS-06 自由声子参考不是要减掉的背景

在 g=0、采用对应自由声子真空的参考问题中，σ²=1/(2mω)，独立键场给 S_c(Q)=σ²I。m=ω=1、Nc=36 时：Szz/Nc=λmax/Nc=1/72，Tr S/Nc=1/24。

这是已有讨论中的独立参考；不是 g=0.3 的精确答案，不从相互作用数据中扣掉。projected 参考还受 (1−w0) 因子影响。去均值后保留局域量子涨落，不表示代码未扣背景。[H, §9.4；DISCUSS；U-CORRECTION]

## OBS-07 Binder：候选扩展，非当前前置条件

若为固定复标量 M(Q) 真正采集二阶、四阶矩，可定义 B4=⟨|M|⁴⟩/⟨|M|²⟩²。已有二阶矩区块不足以恢复四阶矩；不能把区块 S 或 λmax 平方当成 ⟨|M|⁴⟩。

若使用空间平移平均，四阶量必须对 |M(T_aX)|⁴ 本身平均，不能用平移平均的二阶量的平方代替。均值扣除、复通道和相位重加权须先确定。

此前讨论的圆对称复高斯基准 B4=2 有条件；一般零均值复高斯变量有 ⟨|M|⁴⟩=2⟨|M|²⟩²+|⟨M²⟩|²。非网格模式的 ⟨M²⟩ 未必为零。这里仅保留讨论边界，不据此声明当前 Binder 已实现或已认证。[DISCUSS]

## 按需 Majorana 能量与中心 Green（2026-10-02）

`--measure-energy` 默认 1，保留 Majorana 能量估计器；显式 0 只省略其
输出字段及区块列，不改变声子/VBS 的测量定义和复数相位重加权。
fermions=0 时可省去测量中心 Green；fermions=1 始终强制它，费米子
矩阵、相位和统计处理保持原定义。输出 measure_energy 和
measurement_green_enabled 区分测量选择。参见 [ADR-008](decisions/ADR-008-real-magnitude.md)。
# VBS 静态极化率扩展（2026-10-05）

`--vbs-susceptibility 1` 在 `--vbs 1 --vbs-mode all` 下保存允许动量上的
中心 VBS 序参量对静态外场的有限投影响应所需原始矩；默认关闭。
采用精确 Mehler 核外场导数，trial 端点固定，连接减法在共同复相位
重加权之后进行。字段命名、实 cos/sin 与复 Fourier 归一化、统计及
适用范围见 [ADR-009](decisions/ADR-009-vbs-static-response.md)。
短 tau 响应不自动代表基态极化率；等时结构因子不是本开关的响应估计量。

## 运行中接受率输出（2026-10-06）

`--acceptance-every N`（非负整数，默认 0 关闭）在每 N 个完成的 sweep
及各阶段最后一个 sweep 后刷新 `PREFIX_acceptance.csv`。
`stage` 为 `warmup` 或 `measurement`，`sweep` 是阶段内已完成的 sweep 数。
`proposed`、`accepted`、拒绝及恢复计数都是该阶段累计值；测量阶段扣除了
预热计数，`acceptance=accepted/proposed`。`elapsed_seconds` 从程序开始计时。
该回调只读取计数，不产生提议、随机数或额外 Green 测量。
运行中可读取完整 CSV 行；终止后至多保留最后一次刷新时的进度，
该文件不是恢复计算的 checkpoint。既有接受率文件也触发输出覆盖保护。

## 有限静态场（2026-10-06）

`--vbs-field-h` 实际改变 Hamiltonian 和声子路径权重，trial 固定；
`Field_m` 为固定 cos(q.R+phi) 的中心投影，`Field_score` 是有限 h 处完整
log 权重导数。`chi_parallel=<Field_m_score>-<Field_m><Field_score>`，
相位比值合并之后连接减法。`Field_absD4` 为逐构型固定复 Fourier 四阶矩。
非零场破坏平移对称，不能把窗口平移平均冒充固定相位二阶外积。
定义、归一化及证据边界见 [ADR-012](decisions/ADR-012-finite-vbs-field.md)。

## 预热路径诊断（2026-10-09）

`--warmup-observable-every N` 每 N 个预热 sweep 保存中心声子路径的
`raw_all_q_max` 与 q 索引。这是未相位重加权的逐构型路径诊断，
不是物理系综估计量。CLI 从 sweep 完成回调读取路径，不认证相位或
构造 Green；关闭诊断时也不安装需要这些测量的预热回调。
正式测量仍按 stride 完整认证相位及所需 Green，提议和测量的
同候选无缓存恢复策略保持。见 [ADR-015](decisions/ADR-015-warmup-path-diagnostics.md)。

## 零场模长探针（更正 note，2026-10-10）

`--vbs-amplitude 1` 保存 A=|F/Nc|、A²及逐构型A⁴（全允许 q、各键）。
`--vbs-amplitude-response 1` 另存 IA=Σt_j A(j)、中心 A*IA；
端点 t=dt/2、内部 t=dt，先取模再积分。对应明确的离散模长源，
固定 trial 下的连续极限是 H_epsilon=H0-epsilon*Nc*|D|，需要 dt/tau 收敛。
共同相位比后计算 VA=Nc*(<A2>-<A>²)、chiA=Nc*(<AIA>-<A><IA>)，
B4=<A4>/<A2>²；不使用复均值绝对值替换平方，不按负号剪裁。
VA是等时方差，chiA不等于线性chi_zz或有限场chi_parallel。
`analyze_vbs_amplitude.py` 固定模式并做逐链联合误差、Binder局部斜率峰与
窗口/块长检查；峰只作伪临界候选。详细离散定义及保真边界见
[ADR-019](decisions/ADR-019-vbs-amplitude-peaks.md) 和
[用法](../VBS_AMPLITUDE_README_CN.txt)。
