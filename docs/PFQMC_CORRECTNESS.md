# PfQMC 正确性验证

所有新测试使用 `tau = Ntau * dtau` 表示总投影长度，
`tau_mid = tau / 2` 表示中心。`Ntau` 为偶数，声子路径有
`Ntau + 1` 个积分节点。保留生产 API 的 `intervals`、`dt` 等名称，
不改命令行接口。

## 与生产计算相同的约定

Hamiltonian 为

$$H_M=i\sum_b u_b(J_{\mu(b)}+gX_b)c_{A(b)}c_{B(b)},\qquad
H_{\rm ph}=\sum_b[P_b^2/(2m)+kX_b^2/2],\quad k=m\omega_{\rm ph}^2.$$

Majorana 算符满足 `{c_i,c_j}=2 delta_ij`。每条 A 到 B 物理键只计一次。
生产生成元是 `K=-iA`，`A_AB=2u_b(J_mu+gX_b)`。
每个声子区间的物质传播在 ket 作用次序中先应用左节点半步、再应用右
节点半步；中间节点出现两次，两端各一次。中心 cut 是 `Ntau` 个
半步因子之后。两侧声子路径独立，不能强制镜像以制造正权重。

物质 trial 是零声子场固定规范模型的 Gaussian 基态。独立 Fock 参考从
多体 Hamiltonian 对角化得到同一个态，不调用 Pfaffian 权重函数。
声子 trial 是宽度 `m*omega_ph` 的谐振子真空，可带位移；自由极限测试
使用零位移。声子作用量为精确 Mehler 核加两个端点 trial 波函数。

`PathWeight.center_green` 的约定为 `G_ij=delta_ij-<c_i c_j>`。
逐构型观测量通过共同复相位分母重加权，再做连接减法。
固定规范 matter 测试不等于完整物理自旋投影测试。

## 新增测试

| 程序 | 参考及覆盖 |
|---|---|
| `kitaev_ssh_pfqmc_free_phonon_validate` | 独立 Gaussian 精度矩阵求逆，作用量、条件分布、局部比值、L=2/4/6、多 omega_ph/tau/dtau、全允许动量 Szz/chi、实际 g=0 局部与整键采样 |
| `kitaev_ssh_pfqmc_fock_validate` | 4/8 Majorana Fock 空间显式算符、有序多体指数、带符号复权重和比值、中心 Green、键能及四点关联、非对易路径、顺序敏感性控制 |
| `kitaev_ssh_pfqmc_cache_validate` | 同 proposal 对规范半步的无缓存完整稳定重算，PBC/shear、24 组缓存/稳定化设置、连续局部与整键更新、实际 MH 接受概率与决策、接受提交和拒绝后的缓存深层值 |

负实权重与复相位见证使用额外的小 Majorana 参考模型，指数因子均有
真实 Hermitian Hamiltonian 来源，没有人为乘 `eta=-1`。这些见证检验
Pfaffian 符号和相位，不声称指定蜂窝参数已经出现符号问题。

生产程序目前没有非等时 Majorana 关联及键能四点关联的测量 API。
新 Fock 测试用生产 Green/传播矩阵重建这些量，再与独立多体插入比较；
这是传播及 Wick 收缩验证，不能记为生产测量接口验证。

## 自由声子参考与统计门槛

零位移自由真空 trial 下，精确 Mehler 离散化在节点上满足

$$C_{ij}=\langle X_iX_j\rangle=
\frac{1}{2m\omega_{\rm ph}}e^{-\omega_{\rm ph}|i-j|dtau}.$$

因此有限 tau 下中心均值为零，方差与允许动量 Szz 为
`omega_ph/(2k)`。独立精度矩阵先由 `coth`/`sinh` 核系数和端点 trial
组装再求逆，也与实际生产作用量及条件分布核对。

现有零场静态响应用 Mehler 核的外场 score。离散积分权重为

$$w_j=\frac{\tanh(\omega_{\rm ph}dtau/2)}{\omega_{\rm ph}}
\begin{cases}1&j=0,Ntau,\\2&0<j<Ntau.\end{cases}$$

用该权重定义 `I_z` 时，有限投影参考为

$$\chi_{zz}=N_c[\langle D_z(\tau_{\rm mid})I_z^*\rangle
-\langle D_z(\tau_{\rm mid})\rangle\langle I_z^*\rangle]
=\frac{1-e^{-\omega_{\rm ph}tau/2}}{k}.$$

它在长 tau 极限趋于 `1/k`。有限 dtau 的梯形积分具有额外离散偏差，
测试单独输出这项差异，不修改生产估计器使其符合另一种积分定义。

统计测试使用固定随机种子、40 个等长 block、共同相位分母和联合删块
jackknife，预先固定 `6*stderr + 2e-10` 门槛。确定性测试使用舍入误差
门槛。短统计测试不是相互作用模型热化、长投影或自关联收敛证明。
逐构型原始 Szz 与独立 Fourier 外积严格比较；自由统计检查额外用每个
q 的样本均值形成连接协方差，因此保留有限样本均值扣除的偏差。
这只是零均值自由 Gaussian 参考的统计诊断，不替代 OBS-03 的生产
均匀背景扣除定义，也不改变相互作用数据的处理。

## 运行

构建依赖与现有项目相同：C++17、gfortran、Eigen、Boost、MKL。
设置 `EIGEN3_INCLUDE_DIR`、`BOOST_INCLUDE_DIR` 和 `MKLROOT` 后：

```bash
bash scripts/check_pfqmc.sh build-pfqmc-correctness audit/my-new-check
```

输出目录必须不存在。每项保存 stdout、stderr、退出码、耗时及 binary
SHA256，汇总为 `summary.json`。已有二进制无需重编译：

```bash
export LD_LIBRARY_PATH="$MKLROOT/lib:$MKLROOT/lib/intel64:${LD_LIBRARY_PATH:-}"
python3 scripts/run_pfqmc_correctness.py \
  --build-dir build-pfqmc-correctness --output-dir audit/my-new-check
```

`--group free|fock|cache|regression|optimized` 可以单独运行。整套包含三个新测试和
已有的 Mehler 静态响应、数值失败注入、病态 Green 独立 MP50 参考、
费米子 Fourier/规范测量，以及 L=2/4/6/9/12 实数权重和二分格谱优化
验证。后两项覆盖小系统缓存测试不会进入的大尺寸优化分派。
已有 `phonon_validate` 另外验证 Gaussian 积分的绝对归一化与自由路径
范数，因此整套同时检查声子权重常数和物质 Pfaffian 振幅的绝对大小。
`mc_validate` 覆盖生产 `PhaseStatistics` 的复相位重加权、共同分母、
联合 jackknife 和不可解析相位标记，补足自由极限相位恒定时不能检查的
统计路径。
失败注入还运行生产 CLI 并检查六种
complete/aborted 状态，验证同候选恢复、失败中止、计数和 RNG 不变。

CMake 注册三个新测试，可用 `ctest -L pfqmc_correctness --output-on-failure`。
默认 free 测试有 11 组独立直接采样及 L=2 两种实际 MC 更新。
`--direct-only` 仅跳过实际 MC；不代表 MC 已验证。

扩展实际 MC 用 `kitaev_ssh_pfqmc_free_phonon_validate --mc-case INDEX`，
INDEX=0/1 为 L=2 local/bond，2/3 为 L=4，4/5 为 L=6。
每例有 100 预热 sweep、800 次测量、stride=4。
自由统计链使用生产默认的仅奇异值、节点合并与延迟相位模式，每次测量
都认证签名 Pfaffian；缓存测试另外检查 eager 与 deferred 两条路径。
预计超过 5 分钟的计算按根目录规则提交 HPC20015；PBS 示例及本次
冻结源码和执行日志保存在 `audit/pfqmc-correctness-20261008/`。
数学库线程均为 1，每个数组任务运行独立随机链；长组先提交，不依赖
数组索引控制调度次序。

实际通过情况只以本次 audit 报告及退出码为准，不沿用历史运行结论。
