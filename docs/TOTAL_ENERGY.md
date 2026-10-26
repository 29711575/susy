# 总能量与分项测量

接口：`--measure-total-energy 0|1`，活动源码默认0（2026-10-23按用户要求改为关闭）。
新任务脚本显式设置`--measure-total-energy 0`，需要时才设为1。
已验证的历史冻结二进制省略时仍跟随`--measure-energy`，因此使用它提交新任务时
也必须显式传0；冻结源/二进制不原地改写，后续构建再冻结新默认。
旧 `majorana_energy` 的定义、中心测量位置与列位置保持原样。
总能量开启时在全部旧原始列之后追加六列，共同使用原来的复相位分母：

| 原始列 | 定义 |
|---|---|
| `energy_bare_matter` | `sum_b J_mu u_b <i c_A c_B>` |
| `energy_matter_phonon` | `sum_b g u_b <X_b i c_A c_B>` |
| `energy_phonon_kinetic` | `sum_b <P_b^2/(2m)>` |
| `energy_phonon_potential` | `sum_b <k X_b^2/2>` |
| `energy_external_source` | `-sum_b F_b <X_b>`；无外场时为0 |
| `energy_total` | 上述五项之和，包含零点能 |

这些都是整个有限晶格的广延量。除以 `Lx*Ly` 得到每晶胞值。
最终 JSON 的 `total_energy` 记录 schema、原始偏移、切面、时间、归一化与
有限步长约定；`gauge_u` 记录实际规范符号。`timing` 单列两次额外 Green
及能量计算耗时；该时间已包含在 measurement callback 中，不能再相加。

## 为什么新分项使用声子传播中点

当前离散传播为 `T_h=F exp(-h H_ph(F_source)) F`，其中
`F=exp(-h H_matter(X)/2)`；开放路径使用偶数 M 个区间、M+1个节点。
旧中心测量切面位于两个 matter 半步之间。位移与这些半步对易，声子动能
则不对易，因此把相邻 Mehler 核的导数直接加到旧中心 matter 能量上，
不能得到同一切面上的有限 dt 总能量。

新估计器选中心相邻的两条声子链路，分别在各自传播中点插入所有能量分项，
再等权平均。时间是 `tau/2-dt/2` 和 `tau/2+dt/2`，canonical matter
factor cut 是 `M-1` 与 `M+1`。matter 能量使用相应切面的 Green，
线性 X 耦合使用下面的精确条件均值。

这是指定离散投影下的对称混合插入估计器，与旧 `majorana_energy` 的有限
dt 测量位置不同。不能把两种位置的分项混合相加，也不宣称其有限 dt/有限
投影值满足变分上界。连续时间与长投影比较仍需相应收敛检查。

## 对项目 Mehler 核的推导

本节为对实际离散测度的推导。谐振子传播核的一般基础参见
[Moriconi, Am. J. Phys. 72, 1258 (2004)](https://arxiv.org/abs/physics/0402069)
和 [NIST Mehler 公式](https://dlmf.nist.gov/18.18.E28)。

单键 `H_F=P^2/(2m)+kX^2/2-FX`，`omega=sqrt(k/m)`，`q=F/k`，
一条链路两端坐标 x,y，长度 h。把传播核平分后，中点条件分布正比于
`K_F(x,z;h/2) K_F(z,y;h/2)`，均值与方差为

```
mu = q + (x+y-2q)/(2 cosh(omega*h/2))
v  = tanh(omega*h/2)/(2m*omega)
```

因此中点 `X^2=mu^2+v`，不是简单使用端点平方平均。外场势能为 `-F*mu`。

由于同一个 H_F 与自身的传播对易，链路中的 H_F 插入为

```
E_F = -d_h log K_F(x,y;h)
    = omega/2*coth(omega*h)
      + k/8*sech^2(omega*h/2)
        * [(x+y-2q)^2 - (x-y)^2/tanh^2(omega*h/2)]
      - F^2/(2k)

E_kinetic = E_F - k/2*(mu^2+v) + F*mu
```

实现使用 `expm1`、负指数和 long double 中间量，避免极小 h 的不必要抵消和
大 `omega*h` 的双曲函数溢出。动能构型估计值可以为负；没有截断、投影或
替换。该估计器仍会在小 dt 下有较大统计方差，长双精度不消除采样噪声。

给定链路端点，传播中点没有 matter 因子，故相应 Green 不依赖插入的中点 z。
`<X i c_A c_B>` 的条件估计正好是 `mu * <i c_A c_B>_cut`。
没有抽取额外中点坐标，也不消耗提议或测量随机数。

## 数值策略与兼容

两次目标切面 Green 复用已测量的 canonical half-step factors，调用原 endpoint
认证与阈值。失败后仅在同一构型上重建无缓存参考一次，核对 signed overlap，
再失败则保存原数值诊断并中止；不作为普通 Metropolis 拒绝，也不增加提议计数。
权重、Hamiltonian、接受率、提议 cadence、相位和原始测量不改变。

`--measure-total-energy 0` 可关闭全部新增测量。既有运行任务保留原可执行文件。
旧检查点继续使用其原可执行文件恢复；不能跨新旧二进制改变测量 schema。

## 离线分析

`scripts/qmc.py analyze` 自动读取六个原始列，输出相同联合相位统计以及
`energy_matter_link`、`energy_phonon_link`、分项每晶胞值和 `energy_sum_residual`。
总能量误差保留分项协方差，不独立相加各分项标准误差。

历史块包含 X2 时，可恢复 `energy_phonon_potential_center_reconstructed`；
另有逐键 BondK/BondXK 和显式 `gauge_u` 时，可恢复中心裸 matter/耦合能。
它们明确带 center 标签。没有保存新动能相关信息的旧块不能补造总能量。
统一诊断默认纳入已经保存的总能量、声子动能和势能。

## 独立验证

`kitaev_ssh_total_energy_validate` 只用确定性输入，不运行 MC 链：

- 独立谐振子谱展开和真正的无限空间 X/X2/P2 ladder 插入，检查强迫 Mehler 中点。
- 自由开放路径的解析 Gaussian 二阶矩积分，检查每振子动能/势能各 omega/4。
- 两Majorana的两个固定 matter 扇区与全部三个量子声子：对完整路径高斯精确积分，
  与独立48/64级振子矩阵的实际有限 dt Trotter 传播比较全部六项；含力及位移 trial。
- L2 非对易路径的16维 full-Fock 任意切面 Green/能量比较，物理/非受限 trial、
  g=0/.4/.9、生产/参考两套设置；检查无随机流/状态改变，以及失败重试/中止。

原生、Python 和目标20022构建证据见 `audit/1023/total-energy-estimator/`。
