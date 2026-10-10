更正后的 note：零场模长峰观测量（2026-10-10）

当前依据是 new/note.docx（SHA256 6c78e4ceec66453972f66b0c64736c7c98d66f6012ffcf66a4fb4980adfa94c2）。
旧 note 的静态曲率工具保留作独立参考；当前任务与设计见 docs/decisions/ADR-019-vbs-amplitude-peaks.md。

已实现
1. 每个允许 q、x/y/z 键分别保存逐构型 A=|D|、A2、A4；D=F/Nc，X 未经 g 重标度。
2. 可选完整路径模长积分 IA，以及中心 A*IA。先对每个节点 Fourier 场取模，再积分。
   离散权重 dt/2（两端）、dt（内部）。响应是这个离散源的零场导数，需检查 dt/tau 收敛。
3. 联合相位比计算 VA=Nc*(<A2>-<A>^2)、B4=<A4>/<A2>^2、chiA=Nc*(<AIA>-<A><IA>)。
4. 固定显式非零 q 的 lambda 扫描：局部一/二阶 Binder 拟合及 |dB4/dlambda|，
   每链内联合 block bootstrap、联合 jackknife、两倍块长检查、拟合窗口依赖与网格候选峰。

在已确定物理参数、采样参数及优化开关的 MC 命令中加入：
  --vbs 1 --vbs-mode all --vbs-amplitude 1 --vbs-amplitude-response 1
仅等时矩时省略 response 开关。两个新增开关均默认0，不能与 selected/both 模式同时使用。
无需 vbs-susceptibility 或 vbs-field；也可与它们同时测量，原始列按旧列之后追加。
可运行的本机程序：build-amplitude-20261010/run_mc.sh；不能用旧 build-guided 程序测新量。
批量生产前需冻结当前新源码并在目标机器原生构建；已有 HPC 冻结作业不会产生新矩。
IA 开启后每次测量约增加 M+1 次全网格 Fourier（中心只计算一次），
以及 15*Nc 个原始复相位加权列（只测等时则为9*Nc）。不增加 Green 重建。

从 /home/asus/newsusy 使用后处理（同尺寸、几何、trial、dt/tau、裸 J/k/m；仅 g 改变）：
  MPLCONFIGDIR=/tmp/quick1010-amplitude-mpl OPENBLAS_NUM_THREADS=1 \
  python3 quick1010/kitaev_ssh/analyze_vbs_amplitude.py \
    --n1 3 --n2 3 --bond z --windows 3 5 7 --degree 2 --replicas 1000 \
    --output fast-1006/audit/1010/YOUR_NEW_SCAN_ANALYSIS \
    path/to/chain1_results.json path/to/chain2_results.json ...
例子的 n=(3,3) 仅适用于既定 L6/shear0 模式；L9 既定研究用 n=(5,4)。
不同尺寸分别分析，不能混入同一条扫描。lambda=g^2/(mass*omega_ph^2)，取 g>=0。
所有输入必须 complete、有匹配原始块、独立 seed；旧 S/chi 块不足以重建 A 与 A4。
--block-sweeps N 可显式重分块：每条链原始块长须整除 N，链尾不得丢弃。
默认采用输入块长最小公倍数，并记录二倍块长检查。块长仍需按自相关检查选择。

输出 amplitude_scan.json、Binder 斜率 CSV、PNG。保留复数均值、联合实虚协方差、
相位/M2/虚部诊断；未解析的 bootstrap 比值不筛掉后继续给误差区间。
候选峰区间是网格分辨率下的统计区间，不含拟合窗口系统误差或投影/尺寸/模式漂移。
另外保存峰所在网格单元；若某窗口仅剩不到3个内部拟合中心，不输出峰位区间。
不足16块、相位或M2未解析、虚部异常会标记；即使这些检查通过也未证明平衡或块独立。
Gamma 的原始矩保留，但该峰分析入口不选 Gamma；全程不逐构型或逐lambda换q找最大值。
普通复高斯模式 B4->2，自共轭实高斯模式 B4->3；不统一指定一个基准。

保真易感率：当前未实现。note 将其列为下一阶段可选项；先完成上述三个探针。
现有 Pfaffian 路径没有 SSE 顶点计数，不能复制计数协方差或从 chi_zz 重建 chi_F。
ADR-019 记录两半耦合重叠的推导方向及必要的混合导数项。

验证及科学来源：../fast-1006/audit/1010/vbs-amplitude-peaks-corrected-note-20261010/REPORT.txt
