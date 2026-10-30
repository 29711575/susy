> 历史归档：以下“当前”、默认值、作业状态和命令属于原记录日期。当前入口见 [仓库 README](../../README.md)。

# SSH–Kitaev 自由模型基线（P0–P1）

实际耦合 Monte Carlo、相位重加权与分块统计入口见 [P3.md](../../kitaev_ssh/P3.md)。

当前新增的连续声子路径、物理权重适配及示例见 [P2.md](../../kitaev_ssh/P2.md)。下文保留
P0–P1 基线说明；其中“声子路径未实现”的阶段性描述不适用于新增 P2 模块。

物理定义见 [CONVENTIONS.md](../../kitaev_ssh/CONVENTIONS.md)，本机实测误差与限制见
[VALIDATION.md](REPORTS.md#report-21)。入口仅执行 `g=0`、固定零通量 gauge、
空间 PBC 的 Gaussian 纯态投影，不执行 Monte Carlo。

## 文件与改动原因

| 文件 | 职责/修改原因 |
| --- | --- |
| `honeycomb_lattice.h` | 二维几何、一个 site 一个 Majorana 的索引、键、跨界信息、六边形和有向路径 |
| `kitaev_ssh_model.h` | 校验静态 J/u/X，生成实 A、Hermitian K_code 和传播矩阵 |
| `free_projector.h` | 薄适配：复用 GaussianTrialState、stack、endpoint、Green、QR/SVD；能量与指定路径 string 估计 |
| `free_projector.cpp` | 可复现的自由模型运行入口，无 HS、phonon MC 或伪占位参数 |
| `validate.cpp` | 单键 Fock、2×2 多体 ED、L=4/6 Bloch/投影/gauge、连续 Dirac 和输入负例 |
| `CONVENTIONS.md` / `VALIDATION.md` / `README.md` | 约定、验证证据、构建与交付说明 |
| `../CMakeLists.txt` | 沿用依赖/编译选项，新增两个 target 和一项 CTest |
| `../inc/pure_projector_weight.h` | 补齐全局 Pfaffian 分支的恒等锚点 `(-1)^(NM/2)` |
| `../inc/pure_projector_fast.h` | 同步修正 stable reference weight 的相同排序符号 |

未删除模块。旧 `spinless_tV.h`、`kitaevChain.h`、density 入口、脚本、
回归源码、数值库/稳定化算法均保留。只有表中两处共享权重符号作了
数学上必要的修正，并由单键奇数模式测试与旧 density 回归覆盖。
旧发布包的 `SOURCE.json`、`SHA256SUMS` 保留为原始基线记录；打 patch 后
其 CMake 与上述两头文件的原始校验和预期不再匹配，不代表这些变化未经记录。

没有为类名清单增加无职责的 `KitaevMajorana` 包装。
没有 `phonon_field.h` / VBS 生产模块：连续路径数据、声子初态端点权重和更新
均属于 P2 以后，当前静态 X 数组不能冒充声子 path integral。

## 编译和测试

在发布目录运行，依赖保持原仓库的 GCC C++、Intel C/Fortran、MKL、
Eigen 和 Boost、随包 PFAPACK。不引入新的运行时依赖。

```bash
source /opt/intel/oneapi/setvars.sh
export EIGEN3_INCLUDE_DIR=/path/to/eigen3
export BOOST_INCLUDE_DIR=/path/to/boost_1_70_0
cmake -S . -B build-kitaev-p01 -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=ifx \
  -DEIGEN3_INCLUDE_DIR="$EIGEN3_INCLUDE_DIR" -DBOOST_INCLUDE_DIR="$BOOST_INCLUDE_DIR"
cmake --build build-kitaev-p01 --parallel 1
ctest --test-dir build-kitaev-p01 --output-on-failure --parallel 1

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE
./build-kitaev-p01/kitaev_ssh_validate
./build-kitaev-p01/kitaev_ssh_free --Lx 4 --Ly 4 --theta 0.4 --dt 0.05 --block 4
./build-kitaev-p01/kitaev_ssh_free --Lx 6 --Ly 6 --theta 0.4 --dt 0.05 --block 4
```

本机 `cmake` 不在 PATH，实际验证使用以下路径（版本与完整结果见验证报告）：

```bash
export PATH=/home/asus/code/driven-new-validation-20260908/green-residual/tools/cmake/data/bin:$PATH
export EIGEN3_INCLUDE_DIR=/home/asus/old-p/kitaev_ssh_pfqmc_hpc/third_party/eigen3
export BOOST_INCLUDE_DIR=/home/asus/.cache/pfqmc_boost_1_70_0
```

设置这些路径后，上述完整命令可原样执行。CTest 串行运行，并设置数值库
单线程；直接运行入口时也请保留这些环境变量。

## 参数与输出

仅支持 `--Lx`、`--Ly`、`--Jx`、`--Jy`、`--Jz`、`--theta`、`--dt`、
`--block`、`--help`。默认是 4×4、1/1/0.6、theta=0.4、dt=0.05、block=4。
`theta` 为每侧长度，总轮廓长 `2*theta`。block 单位是传播 factor。
未知/重复参数、非整数尺寸、不整除时间步、非有限输入返回非零状态。
不接受本轮未实现的 `--lambda`、`--mass`、`--omega-ph`、`--ramp-rate`。

成功输出单行 JSON：site/Majorana/bond 数、耦合、局域通量、两个 Wilson
标签及原始绕行乘积、matter parity、归一化能量、每 site 能量、
`log_norm_squared` 和数值残差。`phonon_sampling=false`、
`physical_spin_projection=false` 明示边界。

零模错误示例（**预期 exit 1**，不调整模型或 gauge）：

```bash
./build-kitaev-p01/kitaev_ssh_free --Lx 6 --Ly 6 --Jz 1
```

## 旧模型 smoke test

无需修改原脚本，指定本轮重编译的旧入口即可：

```bash
export QMC_BINARY="$PWD/build-kitaev-p01/pure_state_driven_qmc"
KITAEV_SMOKE_DIR=$(mktemp -d /tmp/kitaev-ssh-density-smoke.XXXXXX)
for trial in free cdw-even; do
  bash scripts/run-density.sh --L 18 --boundary pbc --parity -1 \
    --protocol quench --Vi 0 --Vf 4 --tau 0.1 --dt 0.01 \
    --initial-state "$trial" --block 16 --burn 1 --measurements 1 \
    --seed 18150915 --output "$KITAEV_SMOKE_DIR/$trial.csv"
done
```

此少量 sweep 仅检查旧入口仍可运行，不用于物理统计精度声明。

## Patch 交付

交付的 `kitaev-ssh-p01.patch` 是相对原发布目录的 unified diff，包含本目录
及三处已有文件的改动，不包含构建产物、临时 PDF 或 smoke CSV。
在未打补丁的原始源码副本中运行：

```bash
git apply --check /path/to/kitaev-ssh-p01.patch
git apply /path/to/kitaev-ssh-p01.patch
```

`git apply` 不要求初始化 Git 仓库。当前工作目录已经实施这些变更，不应
重复应用；共享头文件的原始基线及 patch 可用于审查/恢复这两处局部修改。

## 下一阶段边界

P2 前需推导开放声子时间轮廓和节点到 factor 的映射，并设计
`X[0..2M][bond]`、端点试探权重、关联多因子 proposal/提交与拒绝恢复。
还需验证该轮廓下单份 Majorana 权重的相位，不能从当前静态 `U^dagger U`
正性推出全部声子 configuration 正权重。零模占据政策和物理自旋 parity
投影若要纳入目标，也须显式设计。这里不自动扩展这些工作。
