# SSH–Kitaev 声子投影 PfQMC

本仓库整理 `quick1010` 的生产核心、独立验证和冻结测量快照。
当前数值核心对应 2026-10-28 已验证源码，最近整理日期为 2026-10-30。
研究对象是固定规范与指定 matter 扇区的声子耦合 Majorana 投影问题。

## 从哪里开始

- [当前研究摘要](HANDOFF_REPORT_CN.md)：结果、仍待检查的问题和结论边界。
- [当前状态](docs/STATUS.md)：源码版本、默认设置和可用验证。
- [近期结果与离线作图](measurements/20261030/README.md)：L6/L9 序参量及频率比较。
- [代码地图](docs/CODE_MAP.md)、[物理约定](docs/PHYSICS_CONVENTIONS.md)、[观测量定义](docs/OBSERVABLES.md)。

## 代码与工具

| 用途 | 入口 |
|---|---|
| 生产 CLI | [kitaev_ssh/coupled_mc.cpp](kitaev_ssh/coupled_mc.cpp) |
| 更新与权重 | [phonon_mc.h](kitaev_ssh/phonon_mc.h)、[phonon_projector.h](kitaev_ssh/phonon_projector.h) |
| VBS 关联、线性响应、模长响应 | [模块说明](kitaev_ssh/README.md)、[模长用法](VBS_AMPLITUDE_README_CN.txt) |
| 检查点保存与续算 | [说明](docs/CHECKPOINT_RESUME.md)、[工具](scripts/checkpoint_job.py) |
| 离线分析和任务索引 | [用法](docs/CALCULATION_TOOLS.md)、[scripts/qmc.py](scripts/qmc.py) |
| 独立确定性检查 | [验证记录](evidence/source-sync-20261030/README.md) |
| 静态工具 | [全键松弛](STRICT_STATIC_README_CN.txt)、[固定背景曲率](STATIC_VBS_README_CN.txt) |

生产设置显式选择 physical matter trial、guided bare/beta1 和非零 q 更新；
参数及测量选项见 [AGENTS.md](AGENTS.md)。总能量默认关闭，matter energy 保留。
现有独立 Γ 更新默认关闭；部分数据使用的多尺度 Gamma± 是另一冻结实现，
其 `--collective-beta-scales` 等实验选项没有包含在本仓库普通核心中。

## 构建

需要 C++17、GNU C/C++/Fortran、Eigen3、Boost headers 和 Intel MKL。
PFAPACK 源码及许可证随仓库提供；构建依赖和二进制不随仓库发布。

```bash
export EIGEN3_INCLUDE_DIR=/path/to/directory/containing/Eigen
export BOOST_INCLUDE_DIR=/path/to/directory/containing/boost
export MKLROOT=/path/to/oneapi/mkl/version
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
bash scripts/build_nodeps.sh build kitaev_ssh_mc \
  kitaev_ssh_reference_core_validate kitaev_ssh_spectral_safety_validate
export LD_LIBRARY_PATH="$MKLROOT/lib:${LD_LIBRARY_PATH:-}"
./build/kitaev_ssh_mc --help
./build/kitaev_ssh_reference_core_validate
./build/kitaev_ssh_spectral_safety_validate
```

已有计算保留原可执行文件、运行库和物理参数；续算不能直接换用新构建。
完整构建与确定性验证范围见 [本次导出验证](evidence/source-sync-20261030/README.md)。
通过验证不表示采样已经热化、混合或达到投影收敛。

## 数据和历史

[测量索引](measurements/README.md) 区分 10 月 30 日结果和 10 月 10 日历史快照。
每组保存原始矩、采集时间、cut 和独立 seed 误差口径；它们不是实时任务状态。
Python 复现依赖见 [requirements-analysis.txt](requirements-analysis.txt)。

旧状态、交接报告及硬编码集群路径的提交脚本已归入 [历史材料](docs/history/README.md)，
不作为当前操作指南。源码与导出文件身份见 [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json)
和 [SHA256SUMS](SHA256SUMS)。
