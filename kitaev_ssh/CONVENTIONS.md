# SSH–Kitaev P0–P1 物理与数值约定

本文记录 P0–P1 时点的固定 gauge 自由 matter-sector 基线；该阶段没有
声子路径采样。本 fast-0929 快照另有声子路径与 Monte Carlo 实现，见
`phonon_path.h`、`phonon_mc.h` 和 `coupled_mc.cpp`，仍不是完整自旋 Hilbert
空间的投影。以下测试核验描述属于原记录时点，不表示 2026-10-02 重跑。取 `{c_i,c_j}=2 delta_ij`、`hbar=1`。

## 1. 自由度、几何、索引

最近邻距离为 1。实空间原始基矢和基底为

```
a1 = (sqrt(3)/2, 3/2),  a2 = (-sqrt(3)/2, 3/2)
rA = (0,0),             rB = (0,1)
b1 = (2*pi/sqrt(3), 2*pi/3), b2 = (-2*pi/sqrt(3), 2*pi/3)
ai . bj = 2*pi*delta_ij
R(x,y) = x*a1 + y*a2
```

`0 <= x < Lx`、`0 <= y < Ly`，空间 PBC；本实现要求 `Lx,Ly >= 2`，
更小的退化环面明确拒绝，而不合并重键。

```
Nc = Lx*Ly
Ns = NM = 2*Nc                 # 一个物理格点一个 matter Majorana
Nb = 3*Nc
cell(x,y) = x*Ly+y
site(x,y,A) = cell(x,y)
site(x,y,B) = Nc+cell(x,y)
majorana index = site index
bond(x,y,mu) = 3*cell(x,y)+mu   # mu: x=0, y=1, z=2
Phi.shape = (NM, NM/2)
```

排列是整块 A 后整块 B，不是逐原胞 AB 交错。这个排列允许复用原数值层
的 parity 配对约定，但**不能**沿用旧模型 `NM=2*physical_sites` 的维数规则。
所有坐标/索引的公开转换做范围校验，内部构造跨边界键时才取模。

每条物理键只存储一次，统一 A→B：

| 类型 | A 端 | 未取模的 B 端 | 实空间位移 |
| --- | --- | --- | --- |
| x | A(R) | B(R-a1) | (-sqrt(3)/2,-1/2) |
| y | A(R) | B(R-a2) | (+sqrt(3)/2,-1/2) |
| z | A(R) | B(R) | (0,1) |

`Bond` 同时保存类型、原胞差 `cell_delta`、边界绕行整数 `boundary_wrap`、
未折叠的实空间位移 `displacement`。其中
`B_unwrapped = B_wrapped + wrap_x*Lx*a1 + wrap_y*Ly*a2`。
矩阵元使用 `+=`/`-=` 累加；不由折叠坐标反推键类型或方向。

## 2. 多体 H、实矩阵 A、框架 K_code、传播 B

定义每条 A→B 键的 `u_b in {-1,+1}`、`t_b = J_mu + g*X_b`，反向
`u_BA=-u_AB`。默认 `Jx=Jy=1`、`Jz=0.6`、`g=0`、`X_b=0`。

```
H_matter = sum_b i*u_b*t_b*c_A(b)*c_B(b)
         = (i/4) c^T A c
A_AB += 2*u_b*t_b;  A_BA -= 2*u_b*t_b
A is real, A^T = -A

K_code = -i*A
H_matter = -(1/4) c^T K_code c
K_code_AB += -2*i*u_b*t_b
K_code^T = -K_code; K_code^dagger = K_code

B(dt) = exp(-dt*K_code)
```

`KitaevSSHModel::realAntisymmetricA()` 与 `frameworkGenerator()` 明确区分
这两个矩阵；`bondHopping()` 返回实系数 `t_b`，不是多体算符。没有含义
混淆的 `getHamiltonian()` 别名。模型输入的长度、有限性、gauge 值、通量、
矩阵维数及两种对称性均检查。

原 `inc/skewMatUtils.h::expm(H,x)` 内部用 `SelfAdjointEigenSolver`，
只适用于 Hermitian 输入。这里传 **K_code** 和 `x=-dt`；不能传实反对称 A。
传播结果有限性亦检查。静态同一个 K 的指数因子没有 Trotter 分裂误差。

`U=exp(-dt*H_matter)` 满足逐行指标约定
`U*c_i*U^(-1) = sum_j B_ij*c_j`。这是从多体对易关系推导并用单键直接
验证的等式；转置排列系数会反转符号，不能仅凭矩阵名称判断传播方向。

固定 gauge 的 Majorana 表示、A 的因子 2、基态能量参见
[Kitaev 原论文第 2–4 节，式 (13)、(18)、(22)、(25)](https://arxiv.org/pdf/cond-mat/0506438)。

## 3. 单键 Fock 锚点、Gaussian 子空间和 Green

按 `f=(c_A-i*c_B)/2` 定义复费米子，在 `(|0>,|1>)` 中

```
c_A = [[0,1],[1,0]], c_B = [[0,i],[-i,0]]
i*t*c_A*c_B = diag(t,-t)
K_code = [[0,-2*i*t],[2*i*t,0]]
```

对 `t>0`，多体基态是占据态 `|1>`，能量 `-t`；K 的负本征向量是
`Phi=(1,-i)^T/sqrt(2)`。它是数值层称为 occupied subspace 的空间，
**不是**原来某个 number-conserving 模型的粒子占据列表。
`beta=Phi^dagger*c/sqrt(2)=f^dagger` 湮灭此 `|1>`。
一般情形，`GaussianTrialState::fromMajoranaHamiltonian(K)` 选 K 的全部
负本征向量，`Phi^dagger*Phi=I`、`Phi^T*Phi=0`，确定被
`beta_a=Phi_a^dagger*c/sqrt(2)` 湮灭的纯 Gaussian 多体态。

```
G_code = I - 2*Phi*Phi^dagger       # 未传播的试探态
G_code[i,j] = -<c_i c_j>           # i != j
G_code[i,i] = 0; <c_i*c_i> = 1
<c_i c_j> = delta_ij - G_code[i,j]
```

单键基态 `G_AB=-i`。能量估计是
`E = -i*sum_b u_b*t_b*G_code[A(b),B(b)] = -Tr(K_code*G_code)/4`。
不能把 G 的零对角直接当作 Majorana 自关联。

指定路径 p 的请求中 string 估计量为

```
<c_i * product_(directed edges in p)(u) * c_j>
    = product(u) * (delta_ij - G_code[i,j])
```

每个反向步骤额外带负号。`gaugeStringGreen()` 实现这个指定路径的
估计器，不是生产测量累加器；没有裸 correlator 代替 string，也没有固定 VBS 波矢。

## 4. 投影范数、复用接口与 Pfaffian 分支

`theta` 是**每一侧**的投影长度，`M=theta/dt` 必须为正整数。
归一化试探态的被测状态是 `|psi_theta>=exp(-theta*H)|psi_T>`，

```
Z_theta = <psi_T| exp(-2*theta*H) |psi_T> = ||psi_theta||^2
E_theta = <psi_theta|H|psi_theta>/Z_theta
Z_theta^2 = det[Phi^dagger * B^(2*M) * Phi]  # 静态模型, eta=1
```

单键一次因子物理振幅为 `exp(dt*t)`，不是 `exp(2*dt*t)`；对基态
`log(Z_theta)=-2*theta*E0`。输出名 `log_norm_squared` 表示 `log(Z_theta)`，
不是再把 determinant 当作单份 Majorana 范数。

`free_projector.h` 仅构造 `PureProjectorSlice` 的显式 bra/ket 两半，
调用 `GaussianTrialState`、`PureProjectorStackManager` 和
`pureProjectorEndpointRebuild`。Green 的重建沿用原 thin-QR、SVD overlap
求解；范数沿用 endpoint 的 QR 尺度累计和 `1/2 log|det|`。
没有复制第二套稳定化算法，没有借用 `singleMaj` 平方权重。

当前轮廓严格是 `U^dagger U`，正范数来自这个事实，不是旧 HS1 正性假设。
该适配器不能直接作为未来非静态声子路径的权重/接受率实现。

独立单键测试暴露并修复两处共享权重代码的排序相位：
`pfaffianForSignOfProductWithStatus(G1,G2)` 返回块矩阵
`[[G1,-I],[I,G2]]` 的原始 Pfaffian 相位。在 `G1=0` 时为
`(-1)^(NM/2)`，而恒等算符的连续 lift 应为正。
`PureProjectorWeightEvaluator` 与 `pureProjectorStableReferenceWeight`
现补乘这个常数，与原局部比值例程已有的 `extraSign` 一致。
原始 Pfaffian 接口本身不改。单键 0/1/2/3 因子均对照 Fock 振幅，
2×2 的非对易试探态再核验有限长度范数及相位；修改后重跑旧回归。
这不是一般复杂轮廓或任意大时间步的分支完备性证明。

## 5. Plaquette、环面扇区、matter parity

以 A(R) 为起点的有向六边形是

```
A(R) -> B(R) -> A(R+a1) -> B(R+a1-a2)
     -> A(R+a1-a2) -> B(R-a2) -> A(R)
types/directions: z+, x-, y+, z-, x+, y-
```

沿遍历方向的原始乘积有三个反向符号，因此
`W_p = - product_directed_hexagon(u) = product_six_stored_A_to_B(u)`。
`W_p=+1` 是零通量；默认全 A→B `u=+1` 时有向原始乘积为 `-1`。
这与 [Kitaev 式 (15) 的偶子格到奇子格乘积定义](https://arxiv.org/pdf/cond-mat/0506438)
一致，不能将原始有向乘积 `+1` 误认作这里的零通量。

两条非收缩闭路从 A(0,0) 出发，重复 `x+ z-` 或 `y+ z-`，分别沿
`-Lx*a1`、`-Ly*a2` 绕行。定义
`W_axis = (-1)^Laxis * product_directed_loop(u)`，即偶长闭路的
`product(-i*u)`，与 [Kitaev 式 (16) 的 Wilson 路径约定](https://arxiv.org/pdf/cond-mat/0506438)
相容。所有 u=+1 时 `(Wx,Wy)=(+1,+1)`；原始乘积是
`((-1)^Lx,(-1)^Ly)`。这是约定相关的 Wilson 标签，故程序同时报告二者。
本任务要求的开放 `c_i product(u) c_j` 不额外乘 `(-i)^length`，
不要将其与自旋 Wilson 算符混淆。

`zeroFluxGauge(wx,wy)` 允许显式在跨界缝上的 x/y 键翻号，产生另外三个
环面扇区；空间几何仍为 PBC。生产基线入口固定 (+1,+1)，不会为避零模
偷偷换扇区。验证中的静态非零 g 小系统显式用 (-1,+1)；几何测试遍历四扇区。

以配对 `f_r=(c_A(r)-i*c_B(r))/2` 定义
`P_matter=product_r (i*c_A(r)*c_B(r)) = (-1)^N_f`。
`GaussianTrialState::fermionParity()` 的现有 Pfaffian 约定与这个定义经
单键和 2×2 Fock 检验。程序报告真实试探态 parity；二次演化保持它，
没有继承旧链固定的 `parity=-1`。

局域 Z2 变换 `c'=D*c`、`u'_ij=d_i*u_ij*d_j` 下，
`K'=D*K*D`、`Phi'=D*Phi`、`G'=D*G*D`。
固定配对下 `P'_matter=det(D)*P_matter`，所以 matter parity 本身不是
gauge-invariant 自旋物理量；能量和上述带 string 的估计量不变。

**尚未实施**扩展 Majorana Hilbert 空间的局域约束/自旋物理态投影。
特定几何和 gauge 对物理 matter parity 的筛选尚未加入，不能把这里的
最低 matter 态一概声称为物理自旋基态。

## 6. 独立 Bloch 谱、Dirac 节点和零模政策

用 cell-periodic Fourier 基 `exp[i(k1*x+k2*y)]`（不含 B 基底位置相位），
其中 `k1=k.a1`、`k2=k.a2`。从三条位移独立得到

```
f(k) = Jz + Jx*exp(-i*k1) + Jy*exp(-i*k2)
K_code(k) = [[0,-2*i*f(k)],[2*i*conj(f(k)),0]]
epsilon(k) = +/- 2*|f(k)|
E0 = -sum_discrete_k |f(k)| = -Tr(|K_code|)/4
G_AB(k) = -i*f(k)/|f(k)|
```

周期均匀 gauge 时 `k1=2*pi*nx/Lx`、`k2=2*pi*ny/Ly`。
有限尺寸对照使用这个解析谱和逆 Fourier Green，不再次调用实空间模型
生成器作参考。能量因子与 [Kitaev 式 (32)–(35)](https://arxiv.org/pdf/cond-mat/0506438)
一致：本文 f 的定义不含原论文 f 的前因子 2，几何/Fourier 方向也已明示。

对 `Jx=Jy=1,Jz=0.6`，令 `q=acos(-0.3)`，两节点为
`(k1,k2)=(q,-q),(-q,q)`，即
`k_cartesian=(+/-2*q/sqrt(3),0)`，模倒格矢等价。
节点邻域
`epsilon_+ ~= sqrt(vx^2*delta_kx^2 + vy^2*delta_ky^2)`，
`vx=2*sqrt(3)*sin(q)`、`vy=3*Jz=1.8`。
测试沿两个独立 Cartesian 方向的对称差商与解析速度比较。
L=4、6 的默认 PBC 网格未必命中各向异性节点；不要求有精确零本征值。
这些节点不被指定为未来 VBS 通道或测量波矢。

现沿用 Gaussian factory 的绝对阈值 `|epsilon| <= 1e-10`（J=1 单位），
发现零/近零模即返回 `zero_mode_ambiguous`。不加小质量、不改 PBC、不换
gauge，也不随意选零模占据。L=6、Jx=Jy=Jz=1 的显式 PBC 测试覆盖该错误。
未来支持零模需要显式定义简并子空间占据及 parity 政策；当前绝对阈值
也不是任意能量单位缩放后的通用误差政策。

## 7. P2 前保留的声子约定（P0–P1 历史记录）

更新：以下实现状态是 P0–P1 时点；当前路径、MC 和移位端点功能分别见
`P2.md`、`P3.md`、`P4.md`、`SHIFTED_TRIAL.md`。默认 d=0 的物理约定未变。

未重标度坐标下，`m>0`、`omega_ph>0`：

```
k = m*omega_ph^2
lambda = g^2/k
psi_T(X) proportional to exp[-(m*omega_ph/2)*sum_b X_b^2]
<X_b^2>_vacuum = 1/(2*m*omega_ph)
```

这是 [MIT 8.04 Lecture 8，振子长度与式 (0.27) 的基态 Gaussian](https://ocw.mit.edu/courses/8-04-quantum-physics-i-spring-2013/808334d09369e3a726f6f20d82315c20_MIT8_04S13_Lec08.pdf)
给出的未重标度宽度。最初提示词中 `1/(2*m*omega_ph)` 放在指数系数的
写法不适用于此坐标；方差不能与波函数指数系数混淆。

后续采用带试探波函数端点的**开放虚时间**轮廓，不能套用有限温度声子的
时间 PBC。至少需要 `X[0..2M][bond]`：两个端点有试探波函数权重，bra/ket
在测量切面共用同一个 `X[M][bond]`，不是两份互不相关的坐标。
仍需在 P2 推导：离散动能/势能作用量、端点因子、时间节点到传播因子的
映射（端点/中点等取值规则）、一个坐标改变影响哪些因子及多因子事务提交。
现有 `PureProjectorSlice` 能存矩阵，旧离散 HS configuration 不能直接
表达以上连续路径及 proposal 关联。本轮的静态 `X_b` 不是这种数据结构，
不提供假声子 Monte Carlo、ramp、VBS 或接受率入口。
