# 2026-10-13 quick1010 源码更新

本次将当前已验证的本地源码导出到 GitHub main；旧版 `27f37df0e1d13bf35fbc03b99e14f7c5ec634742` 完整保留在标签 `quick1010-20261010-before-center-green`，不重写历史。

## 改动

- 新增 `--center-green-fast 1`：中心 endpoint 使用 LAPACK 无主元 thin QR，以及已有的仅奇异值 SVD/LU 求解；保留严格秩、rcond、残差和正交检查及完整 uncached 参考恢复。省略参数仍为 0，普通新生产显式写 1。
- 增加中心 Green 验证器、故障恢复测试和 CMake/nodeps 构建目标。
- 增加独立严格静态工具 `kitaev_ssh/static_limit_response.py`，逐条键松弛并计算单畴响应；它不修改 QMC 采样器，使用时需要 SciPy。

guided 整键、固定 q 池、Hamiltonian、trial、声子作用量及主要 VBS 观测实现与旧 GitHub 快照相同。q 池文件选择运行配置；扩池或改变尝试频率不属于本次源码发布。

## 验证与数据来源

历史原生 MKL 的中心 Green 数值对照、CLI 兼容性及故障恢复检查记录在 [历史验证](../../evidence/recorded-validation/center-green/README.md)，设计与限制见 [ADR-021](../decisions/ADR-021-center-green-acceleration.md)。本次导出检查单独保存到 `evidence/source-update-20261013/`；历史记录不冒充本次重跑。

已有 measurements 仍为 2026-10-10 快照；每个生产任务保留自己的冻结源码、二进制及参数。较新的生产使用 source archive `1a37f776e24445ed7a86422e6fd827fcd3e538539d66a0467d72b2495865d003`；此前长 tau 批次使用 `805d57952de6597088b650d00014f2ed4b9fba9b392e86532e848ec75939b55c`，未自动更新。构建或短检查通过不证明热化、平衡、投影收敛或单位 CPU 的统计效率。

## 查看旧版

```bash
git fetch --tags
git switch --detach quick1010-20261010-before-center-green
```

旧版完整仓库包括当时的源码、说明、证据和测量文件。
