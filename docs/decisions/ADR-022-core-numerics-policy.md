# ADR-022：把已验收数值优化统一到核心默认配置

日期：2026-10-13。用户要求减少计算脚本中的优化开关。

## 决定

生产 CLI 从 `kitaev_ssh/numerical_policy.h` 读取固定的 `production-v1`：
combined 权重（仅奇异值 SVD、同节点半步合并与适用的按需构造）、
延迟 phase、适用的 prefix/Cayley 缓存、center Green 加速全部启用。
普通生产脚本不再需要逐项传入这五个数值选项，也不必写 `--numerics`。
数学内核依然保留根据矩阵结构/规模选择实现及数值失败的完整参考恢复。
没有删除严格秩、rcond、残差、正交性检查或相位重加权。

新增单一诊断入口 `--numerics reference`，选择未合并/eager、关闭缓存和
center-fast 的数值配置。它不宣称关闭共享库的所有优化，也不等同于独立
`uncachedReference`。后者继续用于同候选、同随机数的一次认证恢复。

历史五个数值 CLI 选项继续作为诊断/兼容覆盖入口接受，不再逐项占用
普通 `--help`；显式参数覆盖配置的对应字段。结果 JSON 写出
`numerics_requested`、`numerics_profile` 及每项实际数值，兼容覆盖标为
`compatibility-overrides`。不能仅凭命令行省略项推断旧输出设置。

## 有意保留的配置

guided/bond-refresh/local、q 开关/池/方向/幅度/频率影响采样提议与混合，
继续作为运行配置。几何、物理参数、trial、场、投影长度、warmup、seed、
测量/观测量及诊断间隔也保持独立。不能把某些短链上的接受率优势当作
所有参数下的独立样本效率证明。已有冻结 HPC 作业及历史脚本不改写。

## 兼容性边界

普通 CLI 的 center-green-fast 省略默认从 0 改为 1；显式 0 仍有效。
此前其他四项生产 CLI 默认已为 combined/1/1/1。直接数学/测试 API 的
`PhononWeightOptions` 默认不改，以保护独立参考和既有验证语义。
此前中心 Green 加速证据仍见 ADR-021；本次不把新默认包装当作新增
算法加速，亦不声称“绝对正确”或“所有参数都更优”。

本轮构建、固定 seed 的旧二进制显式参数/新隐式配置/新显式参数对照、
关闭 center 对照、参考配置和故障恢复检查保存在
`../../audit/1013/core-numerics-defaults/`。通过结果以该目录报告为准。
