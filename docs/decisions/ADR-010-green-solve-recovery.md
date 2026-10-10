# ADR-010：病态重叠的 Green 求解恢复

日期：2026-10-06。状态：实现；验证与生产状态以 STATUS 和审计记录为准。

## 问题

长 τ 的带符号权重需要中间 Green。已保存的五个退出构型在重叠满秩、
rcond 约 1.6e-7–1.1e-6 时，原 SVD 逆作用的残差超过 1e-10。
仅精化 double 舍入后的重叠不足以保证 Green 的纯态恒等式：重叠形成误差
可被大逆矩阵放大。不能通过普通拒绝、跳过测量、裁剪相位或放宽门槛处理。

## 决定

正常通过检查的计算保持原算法。只有求解或 Green 认证失败时：

1. 对相同 double 系统尝试 pivoted LU；必要时用 long double 累加残差，
   最多四次候选迭代，并以原 double 残差公式检查候选。
2. 若求解仍不通过，或由舍入重叠构造的 Green 不通过纯态检查，
   从输入左右列空间用扩展精度一致重算重叠、求解和 Green。
   求解的残差在其实际精度计算，归一化及阈值仍为原值。
3. 输出 Green 舍回 double 后，继续原 Green 纯态/结构检查。
   signed 层已有的 oblique-Green 相对残差检查保持原定义和门槛。

原 double SVD 的 rank/rcond 检查先于任何恢复。真实秩缺失、条件数不可信、
非有限值以及恢复后仍未认证的结果继续失败，不参与正常拒绝或测量。
没有修改 Hamiltonian、路径排序、随机流、Metropolis、相位重加权或测量定义。

诊断新增 `initial_solve_residual`、`solve_refinement_steps` 和
`extended_precision_solve`；steps=5 表示扩展精度重算，而非第五次 double 迭代。
long double 与 double 精度相同的平台不启用扩展恢复，也不降低门槛。

## 验证

`green_refinement_validate.cpp` 用合法 pure Gaussian 基与独立 MP50 解比较。
48 个病态案例覆盖初始求解超限和 Green 认证失败；8 个正常案例逐元素保持；
rank/rcond、NaN/Inf 继续失败。MP50 使用输入列空间的真实重叠，另记录与
double 舍入重叠的差异，避免把实现自己的舍入系统误当独立物理参考。
验证器保留原 Green 门槛失败，不强制把病态案例的默认严格检查改成成功。

五个已退出的完整测量路径分别回放合并和未合并因子，检查幅值、相位及中心
Green；同时保留前一轮假零权重修复的六个路径回归。
详见 `audit/tau12-green-repair-20261006/`。通过回归不证明所有长链完成或收敛。
