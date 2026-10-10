# 驱动源码核对报告复核：路径表示、中心拼接和半步复用

日期：2026-09-24。性质：报告复核 + 局部源码复读 + 明确列出前提的推导。

本文件不是 SSH 驱动实现报告，不是新增测试日志，也不授权 QMC、HPC 或代码修改。

## 1. 来源和路径身份

- 用户本轮上传的 [源码核对报告](../../audit/driven-source-review-20260924/REPORT.md)，简称 R-DRIVEN-AUDIT。原文保存不改；其 L4 报告工作树为 `/home/asus/newsusy/new_0923`，L5 报告 HEAD 为 `4130a2e20e4223f629c9cfbc048d8e3719dc055b`。这是原报告的记录；本次合并又对当前本地工作树做了现场复核，结果见第 8 节。
- `C0 局部复读摘录`（附件中引用的历史摘录未作为独立本地文件复制） 来自早先上传的 `core-new_0923-gpt-analysis.tar(1).gz`。该包不包含完整 Stage 1；其内容只作为历史 C0 证据，不覆盖当前工作树现场复核。
- 从 `U_R† O U_R` 推导的中心拼接属于 `derivation`，不是 R-DRIVEN-AUDIT 已经给出的实现。
- 外部原始文献 Cuervo、Roy、Boninsegni，arXiv:cond-mat/0412707，§III 式(8)–(9)，仅核对普通坐标路径积分中“两外端试探波函数 + 全路径核乘积 + 中心对角观测量”的通用结构。不将其氦体系、基态极限要求或正概率解释移植成 SSH 驱动结论。下述驱动映射由目标态定义单独推导。

工作树迁移只更新当前路径登记；旧 audit 中的 `/home/asus/susy/...` 是历史来源，不能批量改成新地址来假装原实验发生在新路径。执行命令使用当前仓库根目录与相对路径，并现场读取 Git 状态。

## 2. 可以吸收的源码现状

R-DRIVEN-AUDIT L10–18、L36–54 补足了 C0 的资料缺口：Stage 1 是 spinless tV 链的 HS/纯态试验程序，有 protocol、独立 HS bra/ket contour、dense/Fock 对照测试源码和 CMake 注册；它不是连续声子的 SSH–Kitaev 驱动模块。

SSH 入口仍是固定 `--g`；`MajoranaGenerator` 只接收 X，模型捕获固定 g；没有 SSH 驱动参数和相应测试目标。显式时间上下文、实际速率记录、同协议失败恢复、相应验证确需接入。

该报告明确未运行测试。其“可运行”“这验证了”等措辞只能按已有代码/注册事实理解，不能在 STATUS 中写成“本轮编译和测试通过”。

## 3. 关键更正：两个物理分支不要求两个完整 PhononPath 对象

R-DRIVEN-AUDIT L58 已指出现有 `PhononPath` 有独立节点坐标；但 L60、L102–106、L140、L146 把“必须两条独立 PhononPath”视为最小方案，缺少推导支持。

设单支有 N 个区间，tau_f=N*dt。一个完整开放轮廓可以存储

$$
X_0,X_1,\ldots,X_N,\ldots,X_{2N},\qquad
\tau_{\rm total}=2N\Delta\tau.
$$

定义两个 canonical 分支（各自从初态走到驱动终点）：

$$
X^{\rm ket}_j=X_j,\qquad X^{\rm bra}_j=X_{2N-j},\qquad 0\le j\le N.
$$

于是 `X_ket[N]=X_bra[N]=X_N`，其他节点没有镜像相等约束。**两侧是不同积分变量，不表示它们在耦合后的采样分布中统计独立。**

一个数组已经能储存这些自由度；对称协议只约束确定性 g 的取值，不约束声子构型镜像。是否新增容器是软件设计问题，不是物理上必须加倍路径。

C0 的 `PhononPath` 恰好要求总区间数为偶数，`center()=M/2`，`action()` 只在两个外端加 trial term，局部中心提议连接左右邻居。这支持优先评估保留单个总轮廓 + 分支视图/协议表。尚未由此验证驱动生成元、传播顺序、相位或 Green。

还有一个具体兼容性问题：若把现有类直接用作每个分支，它会拒绝奇数 N；而完整轮廓 M=2N 合法。两对象方案不是无需改动即可复用。

## 4. 中心坐标不是从 bra/ket 任意选一个

此处限于声子坐标表象对角的观测量（当前 X、XX 型 VBS，以及声子空间中的单位算符所配的 matter 观测量）。若以后测量含声子动量/非对角核的对象，需要另行推导。

记终态在声子坐标上的 matter 向量为 `|chi(X_c)> = <X_c|U_R|psi_i>`，则

$$
\langle O(X)\rangle_R=
\frac{\int dX_c\,O(X_c)\langle\chi(X_c)|\chi(X_c)\rangle}
{\int dX_c\,\langle\chi(X_c)|\chi(X_c)\rangle}.
$$

因为

$$
\langle X_b|O(\hat X)|X_k\rangle
=O(X_k)\delta(X_b-X_k),
$$

bra 与 ket 必须在同一个中心坐标拼接。两分支端点若重复存储，也必须通过共享/约束实现同一变量和单次积分，不能把它们当无约束独立端点，再任意挑一个用于测量。

对 factorized 初态，声子权重有结构

$$
W_{\rm ph}[X]=\Phi_{\rm ph}(X_0)\Phi_{\rm ph}^*(X_{2N})
\prod_{\ell=0}^{2N-1} D(X_{\ell+1},X_\ell;\Delta\tau),
$$

另乘保留复相位的 matter 路径权重。**中心没有额外的 trial 波函数。** 若直接相乘两个现成 `PhononPath.logWeight()`，各自都会在中心加 trial term；中心共享后也会多出 `|Phi_ph(X_N)|^2` 及相应归一化。要选双对象方案，必须先处理这些端点/提议条件，不能把两个完整开放路径权重直接相乘。

所以 R-DRIVEN-AUDIT L106 的“定义中心终点使用哪一个分支节点”应改为“定义唯一中心共享坐标，以及同一 cut 的 matter Green”。保留总轮廓时，现有 `path.center()` 可以继续表示该中心。

## 5. 区间中点方案需要更改半步因子归属，不能只加一个 g 参数

以下分析以报告 L88–98 提出的“区间中点 g + 两端节点 X”为前提；不是认定此离散方案已经通过验证。

同一个内部节点 X_j 位于两个相邻区间，因此有

$$
F_j^-=\exp[-\Delta\tau K(X_j,g_{j-1/2})/2],\qquad
F_j^+=\exp[-\Delta\tau K(X_j,g_{j+1/2})/2].
$$

线性非恒定协议下一般 `F_j^- != F_j^+`。既有代码按 node 生成一个半步并写入两个 factor slots，又把内部节点权重合成同一个半步的平方。C0 中涉及 `pathFactors`、`prepareOwnedNode`、`prepareOwnedPath`、`localMergedNode`、候选更新及无缓存恢复后的优化记录重建。

若 `K(X,g)=K_J+g K_X`，有

$$
[K(X,g_-),K(X,g_+)]=(g_+-g_-)[K_J,K_X].
$$

一般非零，所以两个半步的有序乘积一般也不再是 Hermitian 矩阵；不能复用同一个节点谱分解、平方或 tanh-Cayley 公式，更不能强制 Hermitian 化。

第一版参考路径可按 `(node, adjacent interval, protocol)` 明确生成半步并保留完整有序序列，先禁用或拒绝不适用的 node merge、prefix/Cayley 等优化组合。**只关闭 merge 开关但仍重复使用一个 node 半步，也不正确。**

保留全部优化不是第一版驱动正确性的前置条件。以后恢复优化应逐项证明同一离散对象下等价，并公开记录未支持的组合；不得静默使用旧常量结果。

另一种节点取 g 的方案可以单独比较，但属于不同离散规定，不能与当前中点方案混称为同一个已经验证的 Trotter 实现。

## 6. 最小待定稿方案，而非新增实现任务

优先核实以下方案的可行性：现有总轮廓、唯一中心、两外端 trial、显式折返协议与半步区间身份；复用既有声子自由核和提议；先建立无缓存/未合并的正确驱动参考，再按需恢复优化。

确定性参数可写为

$$
g_c(s)=g_i+R\min(s,2\tau_f-s),\quad 0\le s\le2\tau_f.
$$

这不是实物回程驱动的新实验；第二半是 `U_R†` 在完整 norm 轮廓中的表示。实际存储/action 顺序仍需按 Majorana 约定用非对易多体参考核对。

必要的后续验证包括：常量协议逐对象回归；共同中心与外端因子；两支只改一侧的局部候选不会复制到另一侧；非对易分片顺序；时间步长与实际速率；同协议恢复；明确拒绝暂不支持的优化。当前没有新增或执行这些测试。

## 7. 文档状态和执行边界

`CODE_MAP` 和 `STATUS` 吸收新路径、tV/HS 模块与固定 g 现状。`ADR-004` 保持 proposed，不把 R-DRIVEN-AUDIT 的两对象设计升级为 accepted。历史运行结论与 correction 不变。

本轮已经完成的是报告复核与文档增量，没有修改任何 .h/.cpp/.py 生产文件，没有新 QMC/HPC；源码现场确认使用当前用户工作树的只读命令。


## 8. 当前工作树现场确认

本次合并前实际读取当前 `/home/asus/newsusy/new_0923`，HEAD 为 `4130a2e20e4223f629c9cfbc048d8e3719dc055b`。现场源码与上述复核在关键点上一致：

- `PhononPath` 为偶数区间的单个开放总轮廓，`center()` 是唯一中心，`action()` 只在两个外端加入 trial；节点 `j` 与 `M-j` 没有镜像约束。
- `pathFactors`、`prepareOwnedNode`、`prepareOwnedPath` 和 `localMergedNode` 在固定 g 下复用同一节点半步，内部节点可平方或进入谱/Cayley 优化。
- `MajoranaGenerator`/`honeycombGenerator` 的现有生成接口仍只接收 X，模型捕获固定 `g`；`coupled_mc.cpp` 没有有限速率驱动 CLI。
- `reproduction/pure_state_driven_stage1/` 的测试注册和源码确实存在，但它们验证的是 spinless tV/HS 驱动基线，不是 SSH 声子驱动。

因此，本地文档采用“总 `PhononPath` + 协议映射 + 唯一中心”的优先设计；双对象表示保留为待比较方案。中点 g 下相邻区间的半步不能沿用固定 g 的复制/平方/旧缓存，第一版可暂不支持相关优化，但必须显式拒绝或记录限制。上述确认是本次源码阅读结果，不是编译、测试或 QMC 结果。
