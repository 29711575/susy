# 全动量网格与非公度路径结构因子导出

本增量只新增 `export_vbs_momentum.py`、其确定性测试和本说明。
不改 C++ Monte Carlo、投影权重、随机数、边界、现有测量布局或既有后处理。
依赖项目已有 `predict_vbs_momenta.py` 和 NumPy；`--plot` 另需 Matplotlib。

## 1. 现有生产程序已经支持所需的测量

在原运行命令中开启：

```text
--vbs 1 --vbs-mode all
```

这会保存全部 `Nc=Lx*Ly` 个允许动量的完整 3×3 矩阵和一阶量。
已有 `all` 或 `both` 输出可以直接使用，不必重新采样。

需要同时直接测量非公度点以作固定原点诊断时，使用：

```text
--vbs 1 --vbs-mode both --vbs-q-file selected_q.txt
```

对 Jx=Jy=1、Jz=0.5，最小的 `selected_q.txt` 可以是：

```text
Qbare_plus  0.5804306232551663  0.4195693767448337
Qbare_minus 0.4195693767448337  0.5804306232551663
```

这是裸模型的 ±2 K_tilde 候选，不是已测得的相互作用峰位。
其它耦合应重新生成，不要复用这两个数字。
已有预测器还可生成均匀自洽候选和邻近网格点：

```bash
python kitaev_ssh/predict_vbs_momenta.py \
  --Lx 6 --Ly 6 --shear 0 --Jx 1 --Jy 1 --Jz 0.5 \
  --uniform-effective --g 0.3 --mass 1 --omega-ph 1 \
  --uniform-grid 400 --neighbors 1 \
  --output-dir q_candidates_new --prefix ivbs
```

此命令只生成候选，不运行 QMC。输出是 `q_candidates_new/ivbs_qfile.txt` 等文件。
400 是均匀自洽积分网格，不是 QMC 尺寸，也不意味着积分误差已收敛。
生产运行使用原来经过核对的其它参数；上述片段不替代完整运行命令。

## 2. 从一条已完成的链导出

`run_L6_results.json` 和 `run_L6_blocks.csv` 是示例文件名，请替换为真实输出。
脚本要求 `status=complete`、完整网格、相等且非空的区块，且 CSV 计数与 JSON 一致。

```bash
python kitaev_ssh/export_vbs_momentum.py run_L6_results.json \
  --output-prefix spectra/L6 \
  --vertices '0,1;1,0' --points-per-segment 201 --plot
```

`--blocks` 未指定时，自动寻找同前缀的 `_blocks.csv`。
默认路径 `(u,v)=(t,1-t), 0<=t<=1` 经过当前两个裸模型 VBS 候选。
路径使用未折回的倒格分数坐标；横轴是实际倒格基矢定义的累计笛卡尔距离。
标签标明输入分数坐标，不把移动后的 Dirac 点混称为几何高对称 K 点。
脚本自动把路径上精确的允许动量和目标候选插入路径采样点，不用最近点替代候选。

输出：

| 文件 | 内容 |
|---|---|
| `_grid.csv` | 全部允许动量、分数及笛卡尔坐标、三类矩阵的标量通道与区块误差 |
| `_path.csv` | 指定路径的平移平均窗口结果，包括路径上精确的网格点和候选 |
| `_targets.csv` | 精确的裸 ±2 K_tilde、生产 qpoints，以及可选 `--q-file` 中的其它候选；不在路径上的也保留 |
| `_association.csv` | 固定 `Q_D`、`-Q_D` 及四个 `±delta` 邻点的 connected `Szz`、`m2_z`、`R_b1/R_b2` 和联合删块误差 |
| `_matrices.npz` | 上述三表依次拼接后的完整 3×3 矩阵中心值和实、虚部标准误差 |
| `_metadata.json` | 实际参数、文件来源、背景均值、相位诊断、估计器定义及限制 |
| `_path.png`, `_grid.png` | 加 `--plot` 时生成的路径图及允许动量散点图 |

NPZ 中 `matrix_center` 维度是 `(点数, 3, 3, 3)`：第二维依次为 raw、connected、projected；
最后两维为 x/y/z 键方向。`matrix_se_real`、`matrix_se_imag` 有同样维度。
`kind` 和 `uv` 可定位每个点。原 CSV 区块仍应保留：NPZ 的边际误差不能替代跨动量协方差。

每张 CSV 同时保存 `szz`、`trace`、`lambda_max`，以及各自 `/Nc` 和名义区块误差。
例如 `connected_szz`、`connected_szz_se`、`connected_szz_per_cell`。
量纲使用项目未重标度 X，不是 gX；不同尺寸比较时不要混用 S 与 S/Nc。
没有负值截断、白噪声补偿、曲线平滑或对本征值的插值。

## 3a. 固定 VBS 窗口关联比

导出器默认从运行元数据中的实际裸 `Jx,Jy,Jz` 调用预测器，取
`Q_D=K_plus-K_minus`（模倒格矢）并同时输出 `-Q_D`。也可以显式固定预先确定的中心：

```bash
python kitaev_ssh/export_vbs_momentum.py run_L6_results.json \
  --output-prefix spectra/L6_fixed --fixed-q 0.5804306232551663,0.4195693767448337
```

普通 `shear=0` 周期边界使用精确分数步长
`delta1=(1/Lx,0)`、`delta2=(0,1/Ly)`；已验证的剪切环面使用元数据记录的有限倒格基矢，
未支持的几何会直接报错。坐标只做模倒格矢归一化，不吸附到最近网格。

`_association.csv` 的主通道是同一个固定 z 键 connected 通道：
`C(q)=Re[S_T,c(q)_zz]`、`m2_z=C/Nc`，以及
`R_b=1-[C(Q+delta_b)+C(Q-delta_b)]/[2C(Q)]`。中心值是全样本 plug-in；
每个删块副本共同重算相位分母、完整矩阵窗口、零模背景和两个比值，因此保留跨动量协方差。
合法负值原样输出。中心或邻点等价于 `q=0` 时会显式标记，整组 `R_win` 标为不可用；原始 q=0
仍用于背景计算。`raw`、平移平均 `connected` 和逐构型 `projected` 仍在原矩阵表中分别保留。

这个量命名为窗口结构因子关联比 `R_win`，不等同于已认证的 `xi/L`。非公度窗口、公度失配、
采样和投影收敛问题仍需单独验证；本轮没有实现 Binder ratio，也不拟合临界交点。

## 3. 背景扣除的两个口径

采用项目定义 F_mu(q)=sum_R exp(-iq.R)X_mu(R)，S=<FF†>/Nc。

- `raw`：S_T(Q)=sum_k |W(Q-k)|² S_raw(k)，即全部空间平移后的平均。
- `connected`：S_T,c(Q)=S_T(Q)-Nc |W(Q)|² <Delta(0)><Delta(0)>†。
  只扣除三种键各自的均匀系综均值背景；这是默认绘图口径。
- `projected`：S_T,prj(Q)=S_T(Q)-|W(Q)|² S_raw(0)。
  等价于先对每个构型分别扣掉三种键的空间均值，再做空间平移平均；整个零模涨落也被删除。

两种扣除相差 |W(Q)|² S_connected(0)，不能混称 connected。
在每个非零允许网格点上 W(Q)=0，三类结果严格相同。
在非公度点上通常 W(Q)不为0，所以“不画 q=0”不等于已经扣除背景。

q=0 仍保存在原统计、导出表及矩阵中，标记 `is_q0=1`，不画在两张图中。
这样既满足非零波矢观察需求，也保留非公度点背景扣除所需信息。
不应在数据保存阶段删掉它。共轭到任意整数倒格矢的点也同样排除。

## 4. 更改图、路径和区块长度

全部三类估计器总会导出；下列选项只控制绘图：

```bash
# 完全删除逐构型零模；画归一化的 z 键通道。
python kitaev_ssh/export_vbs_momentum.py run_L6_results.json \
  --output-prefix spectra/L6_projected \
  --estimator projected --channel szz --per-cell --plot

# 一条可比较于 M–Gamma–K–M 的路径；使用项目倒格分数坐标。
python kitaev_ssh/export_vbs_momentum.py run_L6_results.json \
  --output-prefix spectra/L6_MGKM \
  --vertices '0.5,0.5;0,0;0.6666666666666666,0.3333333333333333;0.5,0.5' \
  --plot

# 只在同一条链内部合并相邻区块，重新做共同相位重加权和联合 jackknife。
python kitaev_ssh/export_vbs_momentum.py run_L6_results.json \
  --output-prefix spectra/L6_reblock2 --reblock 2 --plot
```

`--reblock` 必须整除原区块数，且至少剩两块，不丢弃余数。
程序不跨 seed 合块，不把多条链伪装为一条连续时间序列。
块独立性、热化和投影收敛仍要另行检验。

为了跨尺寸比较，分别用各尺寸真实结果及同一路径导出，先比较网格点和 S/Nc；
可将 `_path.csv` 的同名列画在同一图里。不同尺寸路径点数可能不同，因为精确网格点被额外插入，
应按各自 `distance` 列绘图，不能按行号对齐。

## 5. 计算成本与物理限制

主循环的完整网格 Fourier 为三个键方向的 FFT；指定 Nq 个非公度点的直接求和为 O(Nq Nc)。
因此无需把绘图所需几百个点全放进生产 q-file：保存全网格，在后处理中重建路径即可。
本导出器对 M 个绘图点、B 个区块的主要工作是 O(B M Nc) 的小矩阵线性组合和 3×3 本征值。
它不调用费米子投影、不重新求大矩阵 Pfaffian，也不生成随机数。

现有 C++ `all/both` 派生输出还含完整实空间关联的双重求和和多档 blocking；
从 selected 改成 both 时，并非所有附加工作都只有 FFT，I/O 也会增加。
这里没有实测完整 QMC 的新增墙钟比例。

离网格路径是平移平均后的有限窗口测量，不是同一短链的固定原点 `Ssel`。
与固定原点的物理系综期望相等，需要 Hamiltonian、边界/规范处理及试探态密度矩阵的平移对称性。
更密的绘图点不意味着更高的物理动量分辨率，也不能去除原环面的公度效应。

有尖峰或峰随尺寸变窄，尚不足以单独证明热力学极限长程序。
在控制公度、窗口及采样/投影误差的可比较序列上，还需检查峰强度及 S_peak/Nc 的尺寸行为。
本增量没有将峰宽强行转成 xi/L，没有新增未经验证的 correlation-length ratio。

## 6. 本轮验证范围

新增 14 项确定性测试通过；连同原窗口、预测器、静态响应测试，共 40 项通过。
验证包括矩形/剪切几何的窗口权重、显式遍历平移轨道、恒定和涨落背景、矩阵先平均、
共同复相位重加权、联合删块、重分块、实际生产字段格式的合成夹具、异常输入拒绝和绘图冒烟测试。

没有真实 QMC 区块随本源码包提供。本轮没有编译 C++，没有新开 QMC 链，
没有输出新的物理结构因子曲线；绘图测试使用了明确的合成数据，不作为研究结果。

运行测试：

```bash
cd kitaev_ssh
python -m unittest test_export_vbs_momentum -v
python -m unittest test_vbs_window test_predict_vbs_momenta test_static_vbs_response test_export_vbs_momentum -v
```
