# ADR-008：二分格实数 magnitude 与按需中心 Green

日期：2026-10-02。基线为 fast-1002 首轮已验收版本；本轮来源、验收和
耗时记录见 [第二轮审计](../../audit/optimization-round2-20261002/REPORT.md)。

用户明确授权先审查 new_0923 后续改进再选择迁入。其当前生产核心与
20260930 的完整源包一致；实二分格完整 SVD 与 ADR-007 重合。迁入小尺寸
分派和 `--measure-energy`，保留 fast-1002 的物理参数默认值；没有替换
整个源码基线或使用 new_0923 生产二进制。

对 K=[[0,-iT],[iT^t,0]]，D=diag(I,iI) 使 D^dag K D 实对称。
因此 D^dag exp(-aK) D 为实矩阵。自由 trial 的占据子空间在该表示下
有实正交基 Q0；Q0 与 D^dag Phi 只相差占据空间内的酉变换，
QR 的 |det R| 和 overlap 的奇异值给出原来的 log|W|。

`freeGroundState` 保留原复数工厂、零模拒绝、Phi、G_T 和诊断，额外构建
并认证实占据基及不可变 source Phi。Nc>=16 的预验证 magnitude 路径
可使用实 GEMM、dgeqrf_work/dorgqr_work、dgesvd_work。工作区为链私有，
因子实表示缓存按 slice 数有界，并保留不可变原矩阵 owner。

每个因子必须精确满足实同子格/纯虚交叉块；不截断小系数。缺少明确 trial
来源、Phi 修改、结构不符、非有限值、LAPACK 失败，或临近原 rank/rcond/
zero 门槛，均重放整个同一候选的原复数 magnitude。回退不是新的接受或
拒绝标准。signed Pfaffian、partial/final Green、canonical half steps、
无缓存同候选恢复以及 eager local h*h 运算保持原实现。

新增 `PFQMC_NO_REAL_MAGNITUDE=1` 可禁用实数 magnitude。原有
`PFQMC_NO_BIPARTITE_SPECTRUM` 和 `PFQMC_NO_OVERLAP_SVD` 保留。

`--measure-energy` 默认 1。显式设为 0 时只省去 majorana_energy 估计器；
若 fermions=0，可不重建测量中心 Green，复数相位仍完整认证和重加权。
fermions=1 始终强制中心 Green。输出记录 measure_energy 和
measurement_green_enabled，省略的能量不输出占位数值；公共估计量及
链/block 边界保持原定义。此开关不省去 signed phase 内部必需的 Green。

缓存的实数组约 M*(2Nc)^2*8 字节；保留原矩阵 owner 还可能延长被拒绝
候选的复矩阵寿命，额外内存高于仅实数组的估计。整体随路径长度和矩阵尺寸增长，
不随采样提议数增长。短回放和速度测试仅验证实现，不提供物理或统计收敛证据。
