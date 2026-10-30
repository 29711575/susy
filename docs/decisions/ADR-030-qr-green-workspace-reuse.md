# ADR-030：QR、Gaussian 框架与 Green 工作区复用

日期：2026-10-27。状态：implemented / deterministic and target-build verified。目标GCC8.3/MKL2021.2构建与666项固定输入验证已通过；完整生产链加速尚未测量。

## 问题与选择

生产重建反复执行同尺寸 QR、Gaussian 框架维护和 Green 求解，局部临时矩阵及 LAPACK 工作缓冲重复分配。选择显式、每链独占的工作区复用，不改变分解算法或物理估计器。

`PureComplexQRWorkspace` 缓存 LAPACK QR 与显式薄 Q 的工作长度查询结果、tau 和 work 缓冲。缓存按矩阵尺寸区分；开始尺寸变更查询时先失效，全部查询成功后提交新尺寸。`formQ` 校验同一 packed 矩阵的尺寸和数据地址，每次调用消费待完成状态，防止失败后误用旧分解。

`PureCanonicalWorkspace` 与实数 magnitude 工作区复用 Gram、isotropy、修正量和 QR 临时存储；输入与乘积输出分离。Green 工作区复用 overlap、右端、解、LU 和主 SVD 存储。SVD 重算时重置后端与向量模式选择，保留原 Eigen / LAPACKE 分派及显式阈值。返回的 Q 和 Green 保有独立矩阵所有权，后续工作区复用不会修改已返回结果。

## 接口和生命周期

可选工作区指针默认为 `nullptr`，保留独立参考调用；小尺寸 Eigen 分派保持。生产 signed/magnitude、center-fast endpoint 和原 local rebuild 复用路径显式传入工作区。prepared Green 失败后的 uncached 双侧 QR 参考路径仍不传工作区。复杂的 extended solve、混合精度恢复临时对象及 corrected SVD 保持独立，不能覆盖主 SVD。

一个工作区归一个链/调用者所有，不能被并发或重入计算共享；工作区 scratch 不可作为该调用的输入帧。Green 工作区复制时生成空缓存，赋值重置缓存，移动转移所有权；复制链不会共享分解存储。`noalias` 仅用于目标与输入独立的 scratch 乘积，不构成对外部输入别名的自动保护。

## 不变项和边界

原 rank、overlap rcond、solve / Green / structure 残差门槛及 canonical correction bound 保持；行列式幅相补偿、Pfaffian 相位、MH 接受公式、RNG 使用和失败后参考恢复政策保持。测量频率与默认开关保持，包括 fermion、全动量 VBS、susceptibility、amplitude / response、matter energy 开启，以及显式 `--measure-total-energy 0`。既有任务和检查点继续使用记录的源与程序。

本决策不引入新的 pivoting、截秩、近似传播或门槛放宽，也不提供热化、混合、统计收敛或全程性能保证。工作区改变可能影响浮点运算派发；须用固定输入与冻结基线比较认证结果和物理输出，不能仅凭减少分配推断整链收益。

## 代码与证据

主要实现位于 `inc/pure_projector_qr_workspace.h`、`pure_projector_canonical.h`、`pure_projector_real.h`、`pure_projector_green.h`、`pure_projector_svd.h`；调用接入为 `pure_projector_fast.h`、`pure_projector_endpoint.h`、`pure_projector_stack.h` 与 `kitaev_ssh/phonon_projector.h`。确定性检查入口为 `kitaev_ssh/workspace_reuse_validate.cpp`，在 CMake 与 `scripts/build_nodeps.sh` 注册。

本轮冻结源共 241 文件，归档 SHA256 为 `a643f1e2ff48a7e8acdfbb64af2f1faae0f8a4a5d870477a295bbf07760b97e1`。固定输入、冻结基线、生命周期与数值门槛检查及目标构建的实际状态统一见 `audit/1027/qr-green-workspace-reuse/REPORT_CN.md`（原工作区验证报告）；本文不替代报告中的证据与限制，不记录尚未获得的性能数值。
