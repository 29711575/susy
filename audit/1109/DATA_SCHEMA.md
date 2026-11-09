# 离线数据与统计口径

## 文件

- `data/seed-moments.json.gz`：174条选定 seed 链的联合复相位矩，限定 λ=2、L6/L7、τ≤24、Jz=0.6。不是逐 sweep 原始数据。
- `data/reference-results.csv`：59个点的冻结结果，包含不足3 seed 的点。
- `provenance/`：生产 manifest、L7 自由谱与 trial 验证、两端快照身份。manifest 的完整参数范围可能大于本包导出范围。
- `derived/recomputed-results.csv`：脚本从每 seed 矩重算的峰值、误差、相位。
- `derived/momentum-spectra.csv`：逐 q 的 m²(q)、χ(q)。CSV 中 `m2_q` 是数值，与 reference CSV 中整数峰位置字段 `m2_q` 的含义不同；本表另用 `q` 表示整数索引。
- `derived/matched-window-fits.csv`：两尺寸同点集、3 seed、相位通过判据的拟合；`window` 是请求窗口，`taus` 是实际使用点。
- `derived/bare-background-reference.csv`：自由参考/实测总量的比值，不是已分离的真实背景比例。
- `derived/bare-shells.json`：X=0 的 Bloch 裸谱；单准粒子能量不是物理固定宇称多体激发隙。

## 每条链

`point` 保留 L、J、omega、lambda_value、g、mass、tau、dt、intervals、seed、replica、生产 sweeps 目标。`source` 和 `port` 是远端来源索引，不是运行分析脚本所需的路径。

`saved` 是冻结时累计保存的 measurement sweep；`retained` 是扣除前4000 measurement sweep 后累积的样本数。保留独立链边界，不把已含祖先的完整续跑与祖先相加。

`raw` 长度为 `2 + 8*(L*L-1)`。前两项为共同相位分母的实部、虚部。对每个 q=1,…,L²−1：

```text
i = 2 + 8*(q-1)
raw[i:i+8] = complex-pair(Szz), complex-pair(ChiRaw),
             complex-pair(ChiD), complex-pair(ChiIconj)
```

每个 complex-pair 按实部、虚部排列，所有量都是带相位的未归一化和。q=n1*L+n2；分数倒格坐标为(n1/L,n2/L)，Γ=q0没有包含。

先跨同参数独立 seed 累加 raw；共同分母为复数 P。令 Z_j=raw_j/P，则逐 q：

```text
m2(q) = Re(Z_Szz)/L²
chi(q) = Re(Z_ChiRaw - L² * Z_ChiD * Z_ChiIconj)
phase_real = Re(P) / sum(retained)
```

不要将相位加权后的 `<I*>` 替换为 `<I>*`，也不要先分别将各个 seed 的比值取平均。m²、χ分别在非零 q 中选峰。

## 误差及拟合

n>1 时删整条 seed，重算相位比、连接减法及峰。误差为
`sqrt((n-1)/n * sum((estimate_delete_i - mean_delete)^2))`。
只有1 seed时误差为null，不是0。

相位判据只是 `abs(phase_real)>5*phase_real_error`。拟合另外要求每点恰好3 seed，并在每个频率内取两尺寸 τ 点的交集。拟合使用原总量，不减背景；ln y 方差近似(SE/y)²，输出协方差为(XᵀWX)⁻¹，不按残差重新缩放。因此当 χ²/dof 很大时，形式SE不包含幂律模型失配。

不能用本包分析逐 sweep 自相关、热化漂移、block cut 敏感性或集体更新接受率。需要这些证据时，应向数据持有人索取必要的远端小型分段统计，而不是声称此处已认证平衡。
