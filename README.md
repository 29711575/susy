# 蜂窝晶格 Kitaev–声子模型的量子蒙特卡洛计算

本仓库实现固定规范场下、键声子与 Majorana 费米子耦合模型的投影量子蒙特卡洛方法。程序对声子的虚时间路径进行抽样，用 Pfaffian（反对称矩阵的 Pfaffian）计算费米子权重，并测量路径中心的关联和响应。

**初次阅读请从 [QMC 代码审阅说明](docs/QMC_REVIEW_GUIDE_CN.md) 开始。** 该说明面向没有参与开发、没有阅读过仓库的研究者，介绍模型、抽样方法、统计处理、代码入口和需要进一步检查的问题。

## 阅读入口

| 想了解的问题 | 文档 |
|---|---|
| 程序计算什么，算法如何连接起来 | [QMC 代码审阅说明](docs/QMC_REVIEW_GUIDE_CN.md) |
| 具体实现在哪个文件、哪些函数 | [代码阅读顺序](docs/CODE_MAP.md) |
| 报告应该读哪份，哪些已经过时 | [报告与验证材料索引](docs/REPORT_INDEX_CN.md) |
| Hamiltonian、边界、时间、宇称和动量的精确定义 | [物理约定](docs/PHYSICS_CONVENTIONS.md) |
| 结构因子、响应、关联比及误差的精确定义 | [观测量定义](docs/OBSERVABLES.md) |
| 已经做过哪些独立检查 | [正确性验证方法](docs/PFQMC_CORRECTNESS.md)、[2026-10-30 实际执行记录](evidence/source-sync-20261030/README.md) |
| 如何从原始统计文件重新分析 | [计算与分析工具](docs/CALCULATION_TOOLS.md) |
| 绝热极限的另一套计算 | [静态代码与结果](limit/README.md) |

## 版本和结果日期

本次文档整理以 GitHub 主分支提交 `b085974c18b5db4935536aa4fe84b5f778428cc5` 为基础，日期为2026年11月16日。

- QMC 数值核心对应2026年10月28日保存的源码，10月30日有重新构建和确定性验证记录。
- 仓库包含10月10日、10月30日的 QMC 结果，以及11月9日的 L=6、7 对比分析；它们是各自日期采集的数据，不是实时结果。
- 11月10日补充的是 `limit/` 中的静态结果，不表示 QMC 核心更新到了11月10日。
- 部分后续计算使用的多尺度集体更新和副本交换实现**没有包含在本仓库 QMC 核心中**。不能仅靠这里的源码复现所有后来计算；须核对具体结果附带的源码和运行参数。

已有测试不能单独证明实际计算已充分抽样、已达到基态或已找到临界点。报告中的统计误差也不包含全部系统误差。

## 构建与查看接口

需要 C++17、GNU C/C++/Fortran、Eigen3、Boost 头文件和 Intel MKL。PFAPACK 源码及其许可证随仓库提供，编译结果和外部运行库不随仓库提供。

下面只构建主程序和三个确定性验证程序，不启动蒙特卡洛计算。在仓库根目录执行；前三个路径应改为本机安装位置。

```bash
export EIGEN3_INCLUDE_DIR=/path/to/directory/containing/Eigen
export BOOST_INCLUDE_DIR=/path/to/directory/containing/boost
export MKLROOT=/path/to/oneapi/mkl/version
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
bash scripts/build_nodeps.sh build kitaev_ssh_mc \
  kitaev_ssh_reference_core_validate \
  kitaev_ssh_endpoint_propagation_validate \
  kitaev_ssh_spectral_safety_validate
export LD_LIBRARY_PATH="$MKLROOT/lib:$MKLROOT/lib/intel64:${LD_LIBRARY_PATH:-}"
./build/kitaev_ssh_mc --help
./build/kitaev_ssh_reference_core_validate
./build/kitaev_ssh_endpoint_propagation_validate
./build/kitaev_ssh_spectral_safety_validate
```

这些命令的历史执行结果和适用范围见[验证记录](evidence/source-sync-20261030/README.md)。本次整理没有重新编译或运行 QMC。

## 目录用途

`kitaev_ssh/` 保存本模型、声子路径、抽样和测量实现；`inc/` 保存通用 Majorana 投影与线性代数实现；`scripts/` 保存构建和分析工具；`validation/` 保存专门验证程序与方法；`evidence/` 保存已经执行的检查记录；`measurements/` 和 `audit/1109/` 保存指定日期的结果及分析；`limit/` 保存绝热极限的静态计算。

`docs/history/` 是历史报告，`docs/provenance/` 是原工作区来源记录。它们保留原始背景，初次审阅不需要通读。`AGENTS.md` 是自动化工具的操作约束，不是物理或算法说明。

文件来源见 [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json)，文件校验清单见 [SHA256SUMS](SHA256SUMS)。旧数据保留原参数与原来源，不能从当前程序默认值反推历史运行设置。
