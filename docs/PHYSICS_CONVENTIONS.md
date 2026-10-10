# 物理、几何与时间约定

状态：既有约定整理；驱动新增约定见草案 [DRIVEN_FTS_METHOD.md](DRIVEN_FTS_METHOD.md)。
来源：[P1](SOURCES.md#p1) 式(1)–(3)、[H](SOURCES.md#h) 第1–6节、[C0](SOURCES.md#c0) 的模型/几何/路径源码。


## 证据边界

- **论文事实**：Hamiltonian、场论和几何约定的外部来源。
- **数学推导**：坐标、倒格、Majorana 矩阵和时间关系的推导。
- **源码现状**：代码中的索引、生成元、路径和默认参数，以当前接口为准。
- **报告结果**：运行报告中的参数只描述对应运行，不改写模型定义。
- **未验证计划**：驱动协议的新约定在验证前保持草案状态。


## PHY-01 Hamiltonian 与规范范围

取 ℏ=1、Pauli 算符 τ，而不是未经转换的自旋 S=τ/2：

$$
H=\sum_{\langle ij\rangle_\mu}(J_\mu+gX_{ij})\tau_i^\mu\tau_j^\mu
+\sum_{\langle ij\rangle}\left[\frac{P_{ij}^2}{2m}+\frac{kX_{ij}^2}{2}\right].
$$

固定静态规范变量后：

$$
H_{\mathrm M+ph}=i\sum_bu_b(J_{\mu(b)}+gX_b)c_{A(b)}c_{B(b)}+H_{\mathrm{ph}},
\quad\{c_i,c_j\}=2\delta_{ij},\quad\omega_{\mathrm{ph}}=\sqrt{k/m}.
$$

每条物理键只存一次，方向 A→B；反向规范边带相应负号。当前生产记录主要是零通量、Wilson=(+1,+1)。四扇区几何验证不等于生产扇区扫描，也不等于局域约束和物理自旋投影已完成。[P1, 式(2)；H, §1.2；C0, CONVENTIONS §5]

## PHY-02 Majorana 因子与代码生成元

$$
H_M=\frac i4 c^TAc,\qquad A_{AB}=2u_bt_b,\quad A_{BA}=-2u_bt_b,
\qquad t_b=J_\mu+gX_b.
$$

C0 中 `K_code=-iA`，`frameworkGenerator()` 与 `realAntisymmetricA()` 不是同一个矩阵。`B(dt)=exp(-dt*K_code)` 的约定须与现有多体/Fock 作用顺序一起使用，不凭符号名称更改正负号。

若 T 为实际 A–B hopping 矩阵，静态能量基准为

$$
E_M=-\sum_j\sigma_j(T)=-\tfrac14\operatorname{Tr}|iA|.
$$

单键 `H=it c_A c_B` 的基态能是 −|t|；复费米子配对定义不同会改变“占据态”标签，不改变这条能量基准。[H, §2.3；C0, CONVENTIONS §2–3]

## PHY-03 参数不能混用

`a=Jz/J` 是交换各向异性；`a0=bond_length` 是最近邻键长。`m` 为声子质量；`m_VBS^2` 是另一种量。弹性常数 k 不是动量向量。J=1 单位下 λ=g²/k；代码的 X 未经过论文式(3)中的 g 重标度。

因此比较 X 序参量与论文重标度序参量时要转换；比较 R 与 2R 时必须使用同一种定义。`g(τ)` 线性不等于 `λ(τ)` 线性。

最近实际运行是 J=(1,1,0.5)、g=0.3、m=ω=k=1。C0 的 CLI 默认 mass=1.7、omega-ph=0.8；运行参数不能从默认值猜测。ω/J=1 不满足 P2 §C 的小频率条件 ω≪J。[H, §1.3、§4；R-SAMPLING；C0, coupled_mc.cpp]

## PHY-04 晶格、Fourier 与度规

$$
 a_1=(\sqrt3a_0/2,3a_0/2),\quad a_2=(-\sqrt3a_0/2,3a_0/2),
$$
$$
 b_1=(2\pi/(\sqrt3a_0),2\pi/(3a_0)),\quad
 b_2=(-2\pi/(\sqrt3a_0),2\pi/(3a_0)).
$$

R=x a1+y a2，q=u b1+v b2，q·R=2π(ux+vy)。三个键场都用 A 端原胞坐标标记；键中心坐标需要额外基底相位，不可直接替换。

$$
 q_x=\frac{2\pi}{\sqrt3a_0}(u-v),\quad q_y=\frac{2\pi}{3a_0}(u+v),
\quad |\delta q|^2=\frac{16\pi^2}{9a_0^2}(du^2+dv^2-du\,dv).
$$

Nc=LxLy、Ns=2Nc、Nb=3Nc。C0 中先排全部 A 再排全部 B，不是每个原胞交错 AB。索引以公开转换接口为准。[H, §2；C0, honeycomb_lattice.h、CONVENTIONS §1]

## PHY-05 键连接与 Dirac 候选

z: A(R)→B(R)；x: A(R)→B(R−a1)；y: A(R)→B(R−a2)。由这些连接得到

$$
f(q)=\widetilde J_z+\widetilde J_xe^{-2\pi iu}+\widetilde J_ye^{-2\pi iv}.
$$

当正的 $\widetilde J_x=\widetilde J_y$ 且 a_eff=J̃z/J̃x 处于相应无隙范围时，选择 [0,1)² 代表有

$$
\kappa=\frac{\arccos(-a_{\mathrm{eff}}/2)}{2\pi},\quad
K_+=(\kappa,1-\kappa),\quad K_-=(1-\kappa,\kappa).
$$

a=0.5 裸候选：Q+=(0.5804306232551663,0.4195693767448337)，Q−为共轭，等价于 K+−K− 模倒格矢。均匀自洽 J̃ 的候选、有限网格点和实际 QMC 峰分开记录。[H, §2–4；P1, 场展开]

## PHY-06 周期几何与反演

整数剪切：T1=Lx a1，T2=s a1+Ly a2。允许动量满足

$$
L_xu\in\mathbb Z,\quad su+L_yv\in\mathbb Z,\quad
u=n_1/L_x,\quad v=n_2/L_y-sn_1/(L_xL_y).
$$

等价有限倒格基矢为 G1=b1/Lx−s b2/(LxLy)、G2=b2/Ly。非零剪切改变有限 Hamiltonian，不是仅改变测量；固定形状序列可以做标度，但改善效果须实测。

反演：q→−q 模倒格矢，即 (u,v)→(−u,−v) mod 1。最近四点的反演对为

- (1/2,1/3) ↔ (1/2,2/3)；
- (1/3,1/2) ↔ (2/3,1/2)。

坐标交换是另一种操作。复矩阵的非对角虚部不应为使标量相同而被抹掉。[H, §5–6；U-CORRECTION]

## PHY-07 三种时间

Monte Carlo sweep 是采样步数；Δτ 是物理虚时间离散步长；驱动 R 是外部参数的虚时间变化率。三者不可混同。

C0 `PhononPath` 为含端点的开放路径：`totalTau()=intervals*dt`，`center()=intervals/2`，中间节点测量。不要直接把它叫 β 或文献中的单侧 Θ。对于本项目拟议的对称双支驱动态定义，τ_total=2τ_f；这个映射以对应协议和离散实现为前提。[H, §1.4；C0, phonon_path.h；DRV-02]

当前优先评估保留一个总开放 `PhononPath`：两外端带 trial，中间只有一个共享中心节点，两侧其余节点不强制镜像。两个完整路径对象不是物理必需条件；若采用双对象表示，必须证明共享中心只积分一次且不重复中心 trial 因子。该路径表示仍是设计选择，不是已实现的驱动接口。
