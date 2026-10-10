# 均匀移位声子试探态：实现与有界对照（2026-09-16）

已实现并在 HPC 用新二进制完成对照。Xz 的均匀投影漂移明显减小；Btrace 有减小迹象但统计较弱。两点结果不构成投影收敛证明。

## 参数与兼容

新增 `--trial-dx/--trial-dy/--trial-dz`，均默认 **0**：各方向在空间均匀、bra/ket 两端相同、每次运行固定的波函数中心 d；不是静态 Hamiltonian 位移，也不改变声子核或 matter 试探态。

`PhononPath` 保存逐键 d，只改变两个端点作用量。`freeMean(j,b)` 给出随时间节点变化的均值，`freeSample` 在原 AR 零均值路径上加该均值。`freeLogNorm()` 是积分后的自由路径范数；原 `logNormalization()` 仍只是核及归一化端点前因子。

`bond-refresh` 在 X−mu 上计算正反提议密度，内部单键临时路径继承当前键的 d，因此仍只用真实 matter 权重绝对值比接受。`local` 同样支持非零 d：只改端点条件均值和端点 localLogRatio。既有路径/因子/缓存事务、Hamiltonian、物理 Pfaffian 与相位均未改动。

路径格式升为 `KITAEV_SSH_PHONON_PATH 2`，参数行后保存逐键 d；v1 按零位移读取。示例入口 `--load` 与显式时间/声子参数或任意 `--trial-d*` 冲突时拒绝；MC 本身没有新增续跑接口。MC JSON 增加 `trial_shift`、`phonon_trial` 和 `matter_trial`。已知旧结果缺失这些字段时按原未位移真空处理。不同 d 禁止合并；同 Hamiltonian/T/h 可比较新旧态；投影、步长比较禁止同时改变 d。

## 检查

专用 `kitaev_ssh_shifted_validate` 在本地和 HPC 均通过：自由 Gaussian 均值最大误差 8.5e−15、协方差 5.5e−15、积分 log 范数 1.7e−14；端点/内部局部变化与条件提议抵消误差 ≤2.2e−15；整键提议抵消 ≤5.6e−15。固定 eta=±1 单键与严格有限离散公式比较：log 范数 ≤8.0e−15、均值 ≤2.1e−15、协方差 ≤5.4e−15（容差 5e−12）。接受/拒绝事务、单键 sweep 继承 d、v1/v2 存档检查通过。

HPC 同编译环境下，d=0 新旧二进制的短链原始块最大差 **0**。CLI 非有限值/重复参数拒绝、非零 d 保存加载及 load 参数冲突检查通过；11 项分析测试含非零 d 分组与比较冲突检查通过；受构建配置影响的旧 density_sign smoke 通过。两个新输出的独立 NumPy 联合统计检查最大误差约 1.8e−15。没有重复全套 P0–P5。

## 实际运行

预定预算见 `shifted_trial_plan.json`，没有结果驱动的加样本。PBS **155138.mgt**，cn1、**4 CPU**：先并行编译新 MC/验证程序，验证通过后四条链各单线程并行。新 MC SHA256 为 `75440fcce5f9781c2f89eb776d78d4209cc46509b07cdf8def34727d1be2edf4`，不是旧 d=0 二进制。

物理参数 L=2、J=(1,1,.5)、g=.3、mass=1.7、omega=.8、Wilson++；d=(.1378676470588235,...), h=.2。该 d 来自本点 −gB/k，不是程序默认值。每链热化400、测量1600、stride1、原始块25、bond-refresh。

| T | 种子 | 每链耗时/s | 接受率 |
| --- | --- | --- | --- |
| 1.6 | 6701 / 6702 | 15.70 / 16.07 | .9300 / .9319 |
| 3.2 | 6801 / 6802 | 26.88 / 27.30 | .8824 / .8871 |

PBS 墙钟 **81 s**（含编译验证），CPU 时间179 s；四链数值失败均0。两组复平均相位为 `1+3.20e−19 i`、`1−6.17e−19 i`，虚部 SE 8.50e−19、1.28e−18；未强制设为1。原始数据/验证日志位于 `build-kitaev-p01/shifted_trial/results/`，PBS记录为 `shifted_trial/pbs_record.txt`；远端目录 `/home/sunxr/kitaev_shift_trial_20260916.1TSYvf`。

## T=1.6→3.2 的漂移

旧 d=0 数据复用 P5 的 t16_a/b/c、t32_a/b/c，各组24000测量；新态各组3200。主块100 sweeps，新链各16块，未跨链拼块。保持 `S=<sFF†>/(Nc<s>)`；Btrace=`Nc*sum_mu(<X_mu>²)` 是声子均匀背景，不是 Majorana 键算符。

下表 dO 为后者减前者；SE 来自原联合相位/分子 jackknife，独立组差值误差平方相加；区间为近似正态边际95%，依赖足够分块。完整 abs(dO)、z、复数估计和保守区间在 `analysis/comparisons.csv`、`groups.json`。

| 观测量 | 旧态 dO±SE | 移位态 dO±SE | 移位态95%区间 |
| --- | --- | --- | --- |
| Xz | +.03433±.00317 | −.00269±.00893 | [−.02020,.01482] |
| Btrace | +.06276±.00375 | +.02173±.01755 | [−.01267,.05614] |
| Szz(q1) | +.00063±.00508 | +.02162±.01214 | [−.00217,.04541] |
| Tr S(q1) | +.00196±.00924 | +.01112±.02236 | [−.03270,.05495] |
| Szz(q3) | +.01016±.00494 | −.01406±.01321 | [−.03994,.01183] |
| Tr S(q3) | +.01253±.00891 | +.02070±.02456 | [−.02745,.06884] |

q1=(0,1)、q3=(1,1) 固定，未逐配置选峰。Xz 的“新漂移−旧漂移”为 −.03702±.00948，边际区间 [−.05560,−.01845]；Btrace 为 −.04103±.01795，区间 [−.07621,−.00585]。这些差来自四个独立分布；不把同 T 的新旧态差直接解释为收敛。

分块限制重要：新 T=3.2 的 Xz SE 随100→400 sweep 块从 .00705 增至 .00971，Btrace 从 .01375 增至 .01957；400-sweep 层每链只有4块，既不能忽略增长，也不能用它宣称精确误差平台。现有≥16块敏感性包络下，新漂移 Xz 区间 [−.02163,.01625]、Btrace [−.01496,.05842]。Xz 改善有证据，Btrace 的改善结论较弱；四个结构因子均未满足 .02/.04 精度。不能从“不显著”宣布投影收敛或 VBS。

图：[均匀背景](../build-kitaev-p01/shifted_trial/analysis/uniform_background.png)、[固定 q 差值](../build-kitaev-p01/shifted_trial/analysis/differences.png)、[分块](../build-kitaev-p01/shifted_trial/analysis/blocking.png)。U16/U32 表示旧零位移，D16/D32 表示本次移位态。

## 命令与下一步

新参数实际调用示例（另选空输出前缀）：

```bash
./kitaev_ssh_mc --Lx 2 --Ly 2 --Jx 1 --Jy 1 --Jz .5 --g .3 --mass 1.7 --omega-ph .8 \
 --trial-dx .1378676470588235 --trial-dy .1378676470588235 --trial-dz .1378676470588235 \
 --intervals 8 --dt .2 --proposal bond-refresh --seed 6701 --warmup 400 --sweeps 1600 \
 --measure-every 1 --block-samples 25 --qr-block 4 --vbs 1 --output NEW_PREFIX
```

完整 HPC 编译和四链命令在 `shifted_trial.pbs` / `run_p5.sh`；通过20015端口上传至新目录后 `qsub kitaev_ssh/shifted_trial.pbs`。本地 `cmake --build build-kitaev-p01 --target kitaev_ssh_mc kitaev_ssh_shifted_validate -j 2`，然后 `ctest --test-dir build-kitaev-p01 -R '^(kitaev_ssh_shifted_trial|density_sign)$' --output-on-failure`。分析：`P5_PYTHON=/home/asus/.cache/pfqmc_plot_env/bin/python bash kitaev_ssh/analyze_shifted_trial.sh`，输出目录拒绝覆盖。

本轮按预算结束。下一步优先固定新态参数改善独立统计与分块稳定性，再决定是否补一个更长 T；不自动倍增投影或展开尺寸/耦合扫描。
