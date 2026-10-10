# ADR-013：signed 传播加速与生产求解容差

日期：2026-10-07。状态：fast-1006 实现；最终短链/保存构型对照通过。
来源：用户要求继续加速，并允许继续降低不必要的计算精度。

## 决定

生产 `reference_singular_values_only` 路径在占据维数 >=32 时，使用
LAPACK 无主元 QR 生成传播 Q。覆盖 R 之前保留复对角行列式及其相位；
QR 基底重相位与 overlap determinant 相位抵消，原 Pfaffian 分支选择保持。
其他 signed QR 分支保留原计算。

signed overlap 的 values-only 分解复用 `PureOverlapSingularValues` 工作区。
rank/rcond 边界仍回退 Eigen；signed 不以 QR 归一化 determinant 的小乘积
决定物理零权重。接口尾部参数默认保留原 magnitude 调用行为。
`pureFastGreenOptions` 传递已有 accelerated_linear_algebra 设置。

生产 signed Green 求解残差门槛改为
`max(调用方原容差, 1e-9, 1e-8 * 实测 overlap rcond)`。
修正前次使用配置 minimum_overlap_rcond 判断、导致默认不生效的条件。
这降低的是恢复/验收精度，矩阵仍是 complex<double>。
Green 残差门槛保持 1e-8，rank/minimum-rcond 保持 1e-12；失败仍精化并
按需用 long double 一致重算。独立 PureProjectorOptions 默认及中心测量
路径的严格容差保持。没有用 stochastic error 代替实际数值对照。

这项生产容差选择细化 ADR-010/011；旧条目的原容差记录保留为历史。
Hamiltonian、随机抽样、Metropolis 公式、相位重加权和测量定义保持。

## 证据与范围

最终重编译程序的 L6/M64 与 L9/M16 两条4测量短链，采样
24.97→22.16秒、85.54→78.39秒，约减少11.3%/8.4%；
最终路径、接受/拒绝计数完全一致，原始输出最大缩放差8.73e-13，零失败。
单次短链配对不是大统计量长链耗时保证。

五个保存退出构型回放均成功，权重、相位、中心 Green 和能量最大缩放差
1.94e-14；合计45.07→44.16秒。12个n=12/36、rcond约1e-4/1e-6/1e-7
合法Gaussian案例对独立MP50参考，Green最大相对误差2.96e-10。
既有SVD边界、deferred phase与实数传播回归通过；deferred的小det精确
旋转振幅案例补充values-only开/关覆盖，原验收阈值保持。

比较器最初对L9的两个衍生相位SNR报差异：平均phase为1，stderr约
7.36e-17/2.34e-16，导致除法放大舍入噪声。仅在两边误差均低于
64 eps、phase幅值均在64 eps内接近1时单独记录该SNR，不作相对相等要求。
相位本身、误差、原始CSV及物理输出仍按原2e-9阈值比较；两份SNR原值留档。

实数提议SVD的DGESDD替代已尝试，目标n=81/144无稳定收益，已撤回。
最终实数提议实现逐字节保持本轮基线，最终MC/回放程序在撤回后重编译。
专项验证器构建于撤回之前，所测的保留Green/signed源码未再改变。
本轮未提交或替换HPC任务。

记录：`../../audit/signed-svd-precision-20261007/REPORT.md`、`summary.json`。
