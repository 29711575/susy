# ADR-007：固定 g honeycomb 节点的二分格谱分解

状态：implemented / validated，2026-10-02。来源：当前代码的二分格结构、下述代数推导，以及 [本次优化审计](../../audit/optimization-20261002/REPORT.md)。

## 问题和选择

大尺寸 bond-refresh 提议会重建全部时间节点。L12、M80 的原始采样中，因子构建占约 79% 时间，全尺寸复 Hermitian 谱分解及半步构建占约 63%。将固定 g honeycomb 的二分格结构用于矩阵函数，可保留全部模式并减少计算规模。

在现有 A/B 排序下，生成元为

\[
K=\begin{pmatrix}0&-iT\\iT^T&0\end{pmatrix},\quad T\in\mathbb R^{N_c\times N_c}.
\]

对实矩阵做完整 SVD，\(T=U\Sigma V^T\)，则

\[
e^{-aK}=\begin{pmatrix}
U\cosh(a\Sigma)U^T&iU\sinh(a\Sigma)V^T\\
-iV\sinh(a\Sigma)U^T&V\cosh(a\Sigma)V^T
\end{pmatrix}.
\]

\[
2(I+e^{-aK})^{-1}-I=\tanh(aK/2)=\begin{pmatrix}
0&-iU\tanh(a\Sigma/2)V^T\\
iV\tanh(a\Sigma/2)U^T&0
\end{pmatrix}.
\]

这些是矩阵函数恒等式，不要求 T 对称、满秩或非负，也不依赖特定 shear/gauge。半步取 a=dt/2，内部合并因子取 a=dt；时间排序、中心 cut 和 Pfaffian 相位处理保持原有规则。

## 适用条件与回退

`honeycombGenerator` 显式标记候选结构；每次仍检查 K 的形状、有限性、同子格块精确为零、交叉块纯虚，以及两块互为 adjoint。任意非零额外项都会保留完整复谱路径，不使用小量截断。生产分派只在 Nc>=32 时开启；不截断奇异值或零模式。

LAPACK 失败会回到完整谱分解。半步和合并矩阵保留原有 Hermiticity/有限性门槛；不安全的 Cayley 条件保留原 FullPivLU 认证。通用生成元、canonical `pathFactors`、无缓存同候选重试和 eager local 的原 h*h 运算保留。未来时间驱动的相邻区间身份仍受 ADR-004 约束。

`PFQMC_NO_BIPARTITE_SPECTRUM=1` 可禁用这个分派。复 overlap 的 magnitude-only SVD 另使用链私有的持久 zgesvd 工作区，算法仍为原 Eigen/MKL 的 zgesvd；小尺寸、rank/rcond/zero 临界条件保留 Eigen。`PFQMC_NO_OVERLAP_SVD=1` 可禁用该工作区路径。signed Pfaffian 和测量 Green 的求解器没有替换。

## 代码与验证

代码为 `kitaev_ssh/bipartite_spectrum.h`、`phonon_projector.h` 和 `inc/pure_projector_svd.h` / `pure_projector_fast.h`。新增矩阵及 Fock oracle、overlap 边界测试；原回归阈值不变。旧两个测试补齐测量认证再读取 Green，原始失败及修正后的结果均留档。完整采样输出、配置和计时见本次审计。

计时只支持给定固定 g、bond-refresh 参数下的效率结论；短回放与同种子吻合不能证明热化、相位优势、投影收敛或物理相变。
