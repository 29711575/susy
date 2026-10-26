# ADR-026：Gaussian 结构维护与病态中间定号恢复

日期：2026-10-21。状态：本轮有限范围验收；覆盖所列确定性、独立多精度与保存路径验证，不认证热化、混合、收敛或普适前向误差界。

普通 QR 的 `Q†Q=I` 不保证 Majorana Gaussian 子空间的 `QᵀQ=0`。长投影中，微小的非 Gaussian 分量可能增长；此外，中间 trial-overlap 图册可病态，而最终中心端点仍良好。这两类问题需要分别处理。

Gaussian 传播 QR 后执行 `δQ=−conj(Q)(QᵀQ)/2`，再 QR。只允许独立于用户容差的 `128*d*epsilon(double)` 相对步长，并检查正交性与输出 isotropy。权重路径保留附加 `det(R2)` 的幅度和相位；实数捷径采用对应的 `diag(I,−I)` 度量。通用 QR API 默认不进行此修正。缓存恢复同样检查 Gaussian 框架；超界不能变成接受。

Green 的罕见恢复只对已通过独立结构检查的框架进行长双精度微小维护，并对维护后的重叠重新检查原来的秩与条件数门槛，然后求解并构造 Green。若双精度 `G*G` 检查失败，用长双精度累加重新检查**同一组已保存双精度元素**。归一化与阈值不变，不裁剪、反对称化或投影 Green。Gaussian 的 thin-QR、prepared、endpoint 与 stack 消费路径强制执行结构门槛；通用直接求解 API 保留显式选项的语义。

有符号权重仅在合格的中间 Green 求解／结构／纯度失败时尝试增广 Pfaffian。令 `M=P†Q`，`U=[Q,P*]`，`E=[[0,Mᵀ],[-M,0]]`，使用

```
K = [ C  -I   0 ]
    [ I   0   U ]
    [ 0  -Uᵀ  E ]
```

其 Pfaffian 与 `Pf(E)` 的比等价于原来的乘积定号公式，避免显式形成巨大 oblique Green。保留原始 rank/rcond 门槛，独立检查实际传播因子与 Cayley 关系、Cayley 所代表因子的平方行列式一致性，以及框架维护大小。完整主元 Pfaffian 分解与独立行主元 LU 的重构残差、带舍入误差估计的逆残差，给出分支误差预算；只有预算小于正负平方根分支角度余量的一半才接受。普通 Pfaffian 路径与最终 Green 门槛保留。算法不强制实相位或正权重。

输出字段 `maximum_phase_bound_radians` 及内部 `ratio_posterior_phase_bound` 报告原矩阵参考的预算 F，并以 `F < branch_margin/2` 验收。若原 Cayley 有舍入造成的微小非反对称分量，真正的 Pfaffian 目标为 `Ks=(K−Kᵀ)/2`；通用保守相位上界是 **2F**。现有半余量判据已经保证此上界小于实际分支余量，无需修改矩阵或阈值。证明使用 `||Khat−Ks||F≤||Khat−K||F` 及逆矩阵扰动界；完整推导与浮点模型假设见本轮 `augmented/MATH_REVIEW.txt`。对于非严格反对称的原 Cayley，原始平方行列式残差是一项额外一致性检查，不能单独称为精确 Pfaffian 恒等式。

这个误差界约束的是已构造矩阵上的算术，不是全部传播相对于精确物理路径的前向误差证明；因此需要独立多精度实际路径对照。中间单步对数差的前向误差不是单独的拒绝依据，因为最终幅度仍由原来的望远镜求和给出；相关偏差继续记录。

新内核默认启用，无新增采样参数。输出记录 `canonical_frame_roundoff`、`green_precision_recovery`、`augmented_phase_recovery` 三组 process-lifetime 统计，含初始化和恢复尝试，不作为热化或接受率指标。测试强制使用增广路径的开关仅供编译期测试，不保留隐藏的生产环境变量。

验证记录位于 `audit/1021/qmc-numerics-round2/`：独立 Fock 正负／复相位、完整 MH、缓存事务、原秩零点及超界拒绝、完整 A/B/C 路径与多精度中心矩阵、真实 L9 保存路径。最终状态、源与程序 SHA、HPC 构建及重放回执以该目录的最终报告和 `PRODUCTION_SOURCE.json` 为准。

Hamiltonian、trial 选择、离散时间、声子提议、观测量与已有作业不变。选择的生产冻结源以 ADR-025 的已验证 guard 快照为起点；工作区原有 opt-in Gamma 等功能保留，输出差异须单独合入，不能用较早快照整文件覆盖。本决定不声称平衡、遍历性、SUSY、临界指数或计算加速已经得到验证。
