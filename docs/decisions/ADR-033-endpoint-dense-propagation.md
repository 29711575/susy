# ADR-033：中心投影直接读取稠密因子与谱安全检查

日期：2026-10-28。实现及验证证据位于
`audit/1028/reference-core-optimization/`。2026-10-29 仅迁移当前核心：最终冻结源
现位于 `artifacts/production/20261028-reference-core/source/`，归档及二进制
字节未变；迁移记录见 `audit/1029/core-artifact-relocation/`。

## 决定

参考有限温 PfQMC 对已知矩阵结构直接分派的做法，在纯态中心重建中消除
不可变稠密因子的重复物化。仅当 `accelerated_linear_algebra` 开启、调用方
提供私有工作区且因子确有稠密存储时，直接计算 `B*right` 或 `B.adjoint()*left`，
将输出写入工作区的传播矩阵，再与局部 frame 交换所有权。
null 工作区、非加速和压缩因子仍走原表达式。返回的 Q 和 Green 拥有独立存储。

保持因子元素、复共轭、相乘顺序、QR 位置、标度累积、canonical 维护、
秩/条件数/残差门槛及失败恢复。没有改动 Hamiltonian、trial、MH、RNG、
观测量定义或测量默认开关。诊断计数不参与数值决策。

另外修正 `PhononNodeSpectrum::cayleySafe`：完整谱中的任一
`1+exp(-a*E)` 非有限时返回 false，不能跳过该模态后计算条件比。
配对谱以 `abs(a)*sigma_max` 确定最大指数，使正负尺度的安全判据一致。
原生产尺度为正，配对谱的正常生产算术不变。后续因子有限性检查继续保留。

## 验证

- `endpoint_propagation_validate.cpp` 对复数、signed zero、稠密/压缩分派、
  工作区独立性和尺寸变化，以及 L2/L6/L9 多 cut、失败后恢复进行确定性验证。
  相同后端要求逐字节一致，另与原 Eigen 路径比较。
- `spectral_safety_validate.cpp` 用已知配对谱的解析条件比检查溢出、零模、
  正负尺度及 NaN；相同测试对旧基线复现失败，对新实现通过。
- `reference_core_validate.cpp` 用独立 Clifford/Fock 构造覆盖 1/2/3/5 模、
  两种 trial 占据和标量相位分支，并对两个物理 cut 直接计算四点函数。
  标量 eta 控制仅检验分支；真实负权重仍由已有 Fock 验证覆盖。
- 核心权重使用反序因子，其返回 Green 对应物理 cut0；原序 endpoint 的
  Green 取伴随后对应物理 cutN。测试显式记录该约定。新增 oracle 的早期
  符号及顺序错误、独立 NumPy 诊断和失败记录均保留，没有据此更改生产核。

构建、已有确定性验证、保存路径回放及目标运行时验证的实际退出码和哈希
以审计回执为准。最终源码与中间冻结版仅有上述新 Fock 测试文件的差异，
其他245个文件逐项相同，已编译的相关验证保留其实际来源。

本轮8项确定性验证、7组L6/L9保存路径对照及本地/目标最终构建均通过，
保存路径物理量逐项完全一致。整段测量和中心专项各做配对计时后，
没有跨输入一致的稳定时间收益；保留的是经过逐字节验证的复制/分配消除，
不宣称整体加速。完整时间、失败审计与数值误差见本轮报告。

## 范围

这些检查不认证热化、混合、投影或时间离散收敛，也不认证符号问题消失。
保存路径计时仅说明该输入与运行时的行为，不代表完整 Monte Carlo 或同误差吞吐。
本次不运行预检、smoke 或新 Monte Carlo 链。冻结源覆盖 SSH nodeps 生产构建；
完整 CMake 注册检查使用含历史 reproduction 文件的工作区。
