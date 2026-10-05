# fast-1002：SSH–Kitaev 声子 PfQMC 代码与研究交接

本目录整理了截至 2026-10-05 的 `fast-1002` 核心源码、构建所需的源码依赖、物理与数值约定、效率优化记录，以及针对低频声子采样问题的验证摘要。它面向第一次接触本项目的合作者。

## 从哪里开始

1. 阅读 [`reports/physics-check-20261005/PHYSICS_AND_ISSUES_CN.md`](reports/physics-check-20261005/PHYSICS_AND_ISSUES_CN.md)：模型、当前疑点、证据和建议验证步骤。
2. 阅读 [`reports/physics-check-20261005/VALIDATION_CN.md`](reports/physics-check-20261005/VALIDATION_CN.md)：验证范围、结果和局限。
3. 查看 [`docs/CODE_MAP.md`](docs/CODE_MAP.md) 了解源码模块划分，以及 [`reports/OPTIMIZATION_HISTORY_CN.md`](reports/OPTIMIZATION_HISTORY_CN.md) 了解两轮效率优化的范围和计时结果。
4. 生产程序入口是 [`kitaev_ssh/coupled_mc.cpp`](kitaev_ssh/coupled_mc.cpp)；声子路径、投影权重和蒙特卡洛更新分别主要位于 `kitaev_ssh/phonon_path.h`、`kitaev_ssh/phonon_projector.h` 与 `kitaev_ssh/phonon_mc.h`。

## 源码与范围

`inc/`、`src/`、`cmake/`、`kitaev_ssh/` 和 `reproduction/` 保留了构建和运行生产程序及现有回归目标所需的源码。`fast-1002/` 的 C++ 生产代码没有为这份交接包另作物理或算法修改。`reports/physics-check-20261005/validation/` 中的额外程序是审计用小模型或参考计算，明确标为验证代码，不是生产程序功能。

程序实现固定耦合常数的基态投影计算，并在指定规范场及费米子宇称扇区中处理 Majorana 费米子。它不是完整的物理自旋投影计算。PfQMC 原始论文验证了有限温度 Pfaffian 权重方法；本文档另外记录本项目将该权重用于声子路径和基态投影后的具体检查。

## 构建

构建需要 CMake、C++17、Eigen、Boost、Intel MKL，以及兼容的 C/Fortran 编译器。常用环境为 GNU C++ 与 Intel oneAPI C/Fortran。依赖目录由环境提供，不打包进仓库。示例：

```bash
source /opt/intel/oneapi/setvars.sh
cmake -S fast-1002 -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DEIGEN3_INCLUDE_DIR=/path/to/eigen \
  -DBOOST_INCLUDE_DIR=/path/to/boost
cmake --build build --parallel 1
ctest --test-dir build --output-on-failure --parallel 1
```

上面的绝对路径仅为依赖示例，请替换成自己的安装位置。2026-10-05 的独立全新构建和 27 项 CTest 结果记录在验证报告中；它们不表示长投影或低频蒙特卡洛已收敛。

## 主要文献

- Z.-Y. Han, Z.-Q. Wan, H. Yao, “Pfaffian quantum Monte Carlo: solution to Majorana sign ambiguity and applications,” [arXiv:2408.10311](https://arxiv.org/abs/2408.10311)。
- Z. Wu, Z.-Q. Wan, S.-K. Jian, H. Yao, “Emergent spacetime supersymmetry at 2D fractionalized quantum criticality,” [arXiv:2510.07383](https://arxiv.org/abs/2510.07383)。

## 文件校验

`SOURCE_MANIFEST.sha256` 列出这份交接包中的源码文件校验值，方便接收者核对下载是否完整。代码许可证及 PFAPACK 随附许可请分别查看相应目录中的许可证文件。
