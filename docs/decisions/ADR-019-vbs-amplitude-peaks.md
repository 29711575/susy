# ADR-019：更正 note 的零场模长峰探针

状态：implemented / focused native and offline validation，2026-10-10。
当前 `new/note.docx` SHA256 为
`6c78e4ceec66453972f66b0c64736c7c98d66f6012ffcf66a4fb4980adfa94c2`。
用户说明先前 note 错误，要求重新工作；ADR-018 和历史静态结果保留为独立参考。
正确文档的七个链接存于 Word `w:instrText` 中，已提取、核查主要来源。

## 判断与范围

模长涨落、模长源响应和 Binder 斜率是合理的候选探针，可以从零场增加测量。
其峰能否出现、是否对应目标临界模式，需真实扫描检验。经典绝对值序参量及
Binder 矩比实现可参考 [ALPS 官方教程](https://alps.comp-phys.org/tutorials/mcs/mc07/)，
但其 Ising 例子不证明当前量子模型的峰形。
[Jaumà 等的原论文](https://arxiv.org/pdf/2307.13065) 展示 spin-glass Binder 导数
峰定位，方法动机可借鉴，模型与插值选择不直接移植。

先实现 note 优先的三个探针，保真易感率另作可选下一阶段。采样分布、trial、
裸物理参数、优化提议及 Pfaffian/Green 核心不因新增观测量而改变。
当前接口仅支持完整允许网格；各键分别保存，分析默认固定 z 键。
没有实现逐构型取最强键/最大 q，也没有自动改变目标模式。

## 实现定义（本项目的离散选择）

`D_mu(q,j)=sum_R X_mu(R,j)*exp(-iqR)/Nc`，`A_mu=abs(D_mu)`。
中心节点逐构型保存 `A,A2,A4`，其中 `A4=(A2)^2` 在构型层面计算，
不是块平均 `A2` 的平方。可选 `IA=sum_j t_j A(j)` 和 `A(center)*IA`，
`t_0=t_M=dt/2`、其他节点 `t_j=dt`。每个节点取模在时间积分之前完成。
全部原始量交给现有 PhaseStatistics 保存共同 `sum phase` 和 `sum phase*O`。
CSV 名称为 `AmpA/AmpA2/AmpA4/AmpIA/AmpAIA_q<index>_<bond>`；9或15列/q。

有限 dt 的源**定义**为
`W_epsilon=W_0 exp(epsilon*Nc*IA)`，即各时间节点插入对角模长源的对称离散。
直接求归一化比的导数得
`d_epsilon <A(center)>|0 = Nc*(<A(center)*IA>-<A(center)><IA>)`。
固定 trial 下连续极限对应 `H_epsilon=H_0-epsilon*Nc*|D_mu(q)|`。
这是新的非局域源；没有把旧精确线性 Mehler `tanh(omega*dt/2)/omega`
权重替换进模长积分，没有声称当前有限 dt 估计就是连续精确响应。
新增 CLI 仅测 epsilon=0，没有运行有限模长场或改变 Hamiltonian。

`VA=Nc*(<A2>-<A>^2)` 是等时模长方差，不能自动乘 tau 得 chiA。
`B4=<A4>/<A2>^2`。先计算共同相位比，再做连接、乘法和比值；
即使比值有统计虚部，也不把 `<A>^2` 改成 `abs(<A>)^2`。
保留虚部、负估计及误差诊断，不做强制实数/正数修补。
自由平稳高斯场的非自共轭 q 有两个实分量，Rayleigh 模长给 B4=2；
自共轭 q 的实场给 B4=3。此基准是高斯矩推导，不是对临界分布的假设。

## 扫描及误差

`analyze_vbs_amplitude.py` 重用已验证的 completed/JSON/CSV 一致性读取器与
独立链联合 jackknife；支持不等链长度，但保留原链样本权重。
每条链内重分块，不跨链连接，不丢余数。bootstrap 联合重采样完整块行，
每次同时重算 phase、A、A2、A4、IA、AIA，再做派生量和曲线拟合。
统计方法假设重分块之后相关性足够小；16块/5sigma 门槛不认证自相关或平衡。
报告二倍块长的点误差及签名斜率误差，供检查块长敏感性。

同一扫描固定尺寸、剪切、裸 J/m/omega、trial、dt、投影时间和一个显式非零 q。
lambda=g^2/k，g>=0；不同 lambda 必须使用独立 seed，此接口不处理共同随机数扫描。
Binder 在每个固定窗口内用 degree1或2最小二乘，求带符号导数后取绝对值；
所有 bootstrap 曲线用相同规则重复拟合。不外推端点，不采用高阶全局样条。
3/5/7点窗口分别报告网格候选峰及统计区间；窗口差异是系统依赖，不混进单一误差棒。
同时保存峰所在的网格单元；某窗口不足3个内部拟合中心时，峰位显式未定义，
避免只有一个中心也输出零宽峰位区间。
边界最大值提示扩大扫描范围，短块/相位/M2未解析提示不足。
bootstrap 比值失败则保留失败计数并停止相关峰区间，不能筛掉失败副本伪造精度。
VA/chiA 同样报告网格候选峰；全部候选均不自动等于 lambda_c。

通常量子 FSS 的适用条件与修正见
[Campostrini、Pelissetto、Vicari 原论文](https://arxiv.org/html/1401.0788v1)。
由 `B4=B((lambda-lambda_c)L^(1/nu))+...` 求导，可得峰高 `L^(1/nu)` 和
通常峰位漂移 `L^(-1/nu)`；这是有条件的推导。当前非公度模式还需控制 q 失配、
形状与投影时间，不能只用 L6/L9 的一个峰证明标度或 SUSY。

## 保真估计器评估（推导方向，未实现）

[Wang、Liu、Imriška、Ma、Troyer 的原论文](https://arxiv.org/html/1502.06969)
及 [PRX 正式版本](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.5.031007)
给出投影 QMC 两半展开计数的协方差公式，其前提是对应的展开权重幂次结构。
本项目连续声子路径的费米子振幅为 Pfaffian，没有这些顶点计数。
固定 J/k/m 且 lambda>0 时，由 PHY-01 直接推导
`V_lambda=dH/dlambda=g/(2*lambda)*sum_b X_b*u_b*i*c_A*c_B`。
此算符保持平移对称，但仍可能带无隙 QSL 背景。

可用固定 trial 的两半独立耦合重叠
`Z(lambda_L,lambda_R)=<trial|exp(-Theta H(lambda_L))*exp(-Theta H(lambda_R))|trial>`，
归一化重叠的二阶展开给 `chiF_Theta=partial_L partial_R log Z`（对称点）。
若以路径权重 W 实现，必须包括
`<partial_L partial_R log W>+Cov(partial_L log W,partial_R log W)`。
Pfaffian 取对数后一般存在混合导数，不能仅复制 SSE 的计数协方差。
需专用两半耦合/振幅导数实现、固定中心只积分一次、独立小系统重叠验证，
然后检查 Theta/dt 与有限差分步长。当前已保存 chi_zz 不提供这些信息。

基态谱式及额外时间权重的标度依据见
[Albuquerque 等原论文](https://arxiv.org/html/0912.2689) 和
[PRB 正式版本](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.064418)：
`chiF=integral_0^infinity dt*t*<V_lambda(t)V_lambda(0)>_connected`。
[Kitaev–SSH 原文](https://arxiv.org/html/2510.07383v2) 的
`Delta_|phi|2=3-1/nu≈1.9098` 给 nu≈.917；结合常用临界标度可推得
`(chiF/Nc)_sing ~ L^(2/nu-2)≈L^.18`、`alpha=2-3nu≈-.75`。
这两项是理论指数推论，不是当前 QMC 结果，不保证小尺寸尖峰。

## 本次证据

原生 GCC15/MKL：独立 DFT（含剪切）、真正四阶矩、先模后积分、复权重
离散源有限差分，以及6万自由声子配置的 Rayleigh/半正态基准通过。
Python：12项相位/联合误差/独立链/局部斜率/自由模长响应 dt 收敛检查通过。
三个固定种子旧程序/新程序关闭/新程序开启对照，旧相位和测量块列逐列一致；
固定场 h=0 的 absD2/absD4 与新矩一致。L6 短维度 smoke 通过。
五个 lambda 点的短 MC 仅用于全流程验证，不能用于临界点或收敛结论。
复现命令、失败的测试夹具及修正记录、最终原始结果、源码/二进制哈希在
`../../../fast-1006/audit/1010/vbs-amplitude-peaks-corrected-note-20261010/REPORT.txt`。
