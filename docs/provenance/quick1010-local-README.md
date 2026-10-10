# fast-1006 核心代码与计算效率优化

系统正确性测试见 [docs/PFQMC_CORRECTNESS.md](docs/PFQMC_CORRECTNESS.md)：
设置现有构建依赖后，运行 `bash scripts/check_pfqmc.sh`，或使用独立
测试程序/CTest。测试结果与冻结源码保存在新的 audit 目录中。

2026-10-07 已确认为以后代码修改、优化和新计算的默认基线。
当前 fast-1002 的核心改进已包含在这里；它与 fast-0929 保留作历史参照。
源码核对见 `audit/baseline-selection-20261007/comparison.json`。
最新版包含 signed QR/SVD 加速、生产求解容差调整及有限 VBS 外场功能。
当前本机程序为 `build-signed-svd-precision-20261007/kitaev_ssh_mc`；
最近的对照与速度记录见 `audit/signed-svd-precision-20261007/REPORT.md`。

2026-10-06 从当前 `/home/asus/newsusy/fast-1002` 复制完整核心代码开始。
fast-1002 与 fast-0929 保留为来源。当前目录未初始化 Git。

核心包含 `inc/`、`src/`、`kitaev_ssh/`、`reproduction/`、`cmake/`、
`scripts/` 和顶层构建文件，含 PFAPACK 源码及许可证。生产入口是
`kitaev_ssh/coupled_mc.cpp`。VBS 全动量关联、费米子关联、能量及
默认关闭的 VBS 静态响应测量均保留。

2026-10-06 首次优化减少已满足生产 Green 相对误差检查时的扩展精度重算及后续
重复 QR。常规双精度、物理模型、相位重加权、Metropolis 和测量定义
保持；真正无法通过求解检查的构型仍可恢复到扩展精度。
实现、验证和实际耗时见 `audit/precision-speed-20261006/REPORT.md`。

`baseline.sha256` 和 `baseline-core.tar.gz` 保存复制时的271个核心文件，
与 fast-1002 逐文件一致；新审计材料写在本目录的 `audit/`。
旧审计/构建目录未复制；历史文档提到的旧报告应到 `../fast-1002/audit/` 查找。
最新状态以本目录 `docs/STATUS.md` 中 fast-1006 条目为准。

构建入口为 CMake 或 `scripts/build_nodeps.sh`。后者需要设置
`EIGEN3_INCLUDE_DIR`、`BOOST_INCLUDE_DIR` 和 `MKLROOT`，例如：

```bash
bash scripts/build_nodeps.sh build kitaev_ssh_mc kitaev_ssh_green_precision_benchmark
```

`kitaev_ssh_green_precision_benchmark [repeats]` 是局部 Green 求解基准，
其加速倍数不代表整条 Monte Carlo 链的加速倍数。
物理约定见 `docs/PHYSICS_CONVENTIONS.md`，测量定义见 `docs/OBSERVABLES.md`，
工作规则见 `AGENTS.md`。
