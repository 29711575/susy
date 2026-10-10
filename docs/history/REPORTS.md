# 历史报告汇编

整理日期：2026-10-02。以下为来源版本的报告与交接原文，按原文件独立保留；
只调整标题层级和相对 Markdown 链接。不同报告的参数、默认值、工作路径与
“当前”“本轮”表述属于各自记录时点，不代表 fast-1002 的当前验收结果。
历史 audit、构建目录和结果通常未随源码包分发，相关引用可能不可访问。
当前状态以 [STATUS.md](../STATUS.md) 为准；当前构建入口以 [README.md](../../README.md) 为准。

逐字原文件与 SHA256SUMS 保存在 [原文压缩包](original-reports-20261002.zip)。
本次未编译、重跑历史测试或查询 HPC。

## 目录

- [CORE_PACKAGE_README.md](#report-01)
- [CORE_README.md](#report-02)
- [FAST_0929_README.md](#report-03)
- [IVBS_HANDOFF.md](#report-04)
- [IVBS_ROUND3_FIXES.md](#report-05)
- [MERGE_HANDOFF.md](#report-06)
- [PERFORMANCE_HANDOFF.md](#report-07)
- [PERFORMANCE_OPTIMIZATION.md](#report-08)
- [kitaev_ssh/FERMION_GSCAN.md](#report-09)
- [kitaev_ssh/L4_ED_FEASIBILITY.md](#report-10)
- [kitaev_ssh/L4_MIXING.md](#report-11)
- [kitaev_ssh/L4_PERFORMANCE.md](#report-12)
- [kitaev_ssh/L4_PILOT.md](#report-13)
- [kitaev_ssh/PHASE_TAU_SCAN.md](#report-14)
- [kitaev_ssh/SHIFTED_LONG.md](#report-15)
- [kitaev_ssh/SHIFTED_T32_STATIONARITY.md](#report-16)
- [kitaev_ssh/T32_INDEPENDENT_REPLICATION.md](#report-17)
- [kitaev_ssh/T32_REPLICATION8.md](#report-18)
- [kitaev_ssh/T32_T64_PROJECTION.md](#report-19)
- [kitaev_ssh/T64_STEP_HALVING.md](#report-20)
- [kitaev_ssh/VALIDATION.md](#report-21)
- [kitaev_ssh/WEIGHT_OPTIMIZATION.md](#report-22)
- [kitaev_ssh/kitaev-ssh-p5-convergence-and-prompt.md](#report-23)
- [kitaev_ssh/kitaev-ssh-p5-followup-projection-and-prompt.md](#report-24)
- [CODEX_HANDOFF.md](#report-25)
- [docs/STATUS.md](#report-26)

<a id="report-01"></a>

## CORE_PACKAGE_README.md

### SSH–Kitaev core source snapshot

> 历史文档：下文的版本、目录、默认值和测试结果属于原记录时点，不代表 fast-0929 当前状态。当前快照入口见 [FAST_0929_README.md](REPORTS.md#report-03)，源码核对状态见 [docs/STATUS.md](REPORTS.md#report-26)。历史 audit 未随本包分发；本次未重跑这些测试。

Snapshot date: 2026-09-16 (Asia/Shanghai).

This source-only package contains the current pure-state projector framework,
the independent `kitaev_ssh/` model and phonon modules, numerical validation
programs, analysis scripts, build configuration, and the retained legacy
density-model sources and regression tests. Build products, raw Monte Carlo
blocks, PBS output, Python bytecode, editor metadata, and credentials are not
included.

Current physics scope: fixed-gauge matter sector on the periodic honeycomb
lattice; finite-frequency continuous bond phonons on an open projection
contour; physical signed/complex Pfaffian weight; shifted equal endpoint
phonon trial state; bond-refresh sampling; phase-reweighted joint jackknife;
VBS structure factors on the complete momentum grid. Physical spin projection
and a general zero-mode occupation policy remain outside the implemented
scope.

Canonical imaginary-time notation in current output is `tau`. Readers retain
compatibility with historical `theta_per_branch` and `total_time` fields.

Recent L=4 work:

- `kitaev_ssh/L4_PERFORMANCE.md` records the measured performance profile.
- `kitaev_ssh/convergence_analysis.py` supports even L, maps reciprocal
  fractions `(0,1/2)` and `(1/2,1/2)` to the actual grid, uses `Nc=Lx*Ly`, and
  retains full complex momentum data.
- `kitaev_ssh/l4_mixing_plan.json` and `l4_mixing.pbs` define the currently
  running fixed-budget mixing check. No results from that unfinished job are
  present in this source-only snapshot.

The profiling hooks in `kitaev_ssh/path_profile.h` are inactive in production
unless `KITAEV_SSH_PROFILE` is defined. The dedicated target
`kitaev_ssh_profile` enables them in `profile_chain.cpp`.

Build requirements and paths are described by `CMakeLists.txt` and
`scripts/build.sh`. Eigen, Boost, Intel MKL, and the bundled PFAPACK sources are
required by the existing build. Do not infer physics convergence from smoke,
profiling, or short-chain validation programs.

Suggested review order:

1. `kitaev_ssh/CONVENTIONS.md`
2. `kitaev_ssh/honeycomb_lattice.h`
3. `kitaev_ssh/kitaev_ssh_model.h`
4. `kitaev_ssh/phonon_path.h`
5. `kitaev_ssh/phonon_projector.h`
6. `kitaev_ssh/phonon_mc.h`
7. `kitaev_ssh/coupled_mc.cpp`
8. `kitaev_ssh/vbs_observables.h`
9. `kitaev_ssh/convergence_analysis.py`
10. `inc/pure_projector_*` and `inc/gaussian_trial_state.h`


<a id="report-02"></a>

## CORE_README.md

### Core Source Snapshot

> 历史文档：下文的版本、目录、默认值和测试结果属于原记录时点，不代表 fast-0929 当前状态。当前快照入口见 [FAST_0929_README.md](REPORTS.md#report-03)，源码核对状态见 [docs/STATUS.md](REPORTS.md#report-26)。历史 audit 未随本包分发；本次未重跑这些测试。

This directory is the clean, buildable core-source snapshot extracted from
`/home/asus/susy/pure-new-ivbs-perf-merged-20260920`.

Included:

- `CMakeLists.txt`, `cmake/`, `inc/`, and `src/`;
- the top-level `kitaev_ssh` runtime, VBS measurement, candidate diagnostics,
  deferred-phase implementation, validators, and analysis scripts;
- the `pure_state_driven_stage1` sources and fixtures required by CMake tests;
- reusable build/regression scripts and performance profiling source/tools;
- the handoff and implementation documents needed to understand the merged
  semantics.

Excluded intentionally:

- `build-*` directories and compiled binaries;
- `audit/` logs and source snapshots;
- Python caches and generated result files;
- `perf_work/results/`, historical patch stacks, scan output directories,
  and generated predicted-momentum reports.

The build still expects the external Eigen and Boost include directories used
by the merged tree. MKLROOT selection and cache-isolation logic remain in
`CMakeLists.txt`.

Configuration example:

```bash
cmake=/home/asus/code/driven-new-validation-20260908/green-residual/tools/cmake/data/bin/cmake
export CPATH=/opt/intel/oneapi/mkl/2026.1/include
export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2026.1/lib:${LD_LIBRARY_PATH:-}
export EIGEN3_INCLUDE_DIR=/home/asus/old-p/kitaev_ssh_pfqmc_hpc/third_party/eigen3
export BOOST_INCLUDE_DIR=/home/asus/.cache/pfqmc_boost_1_70_0
"$cmake" -S . -B /tmp/pure-new-ivbs-perf-core-build
cmake --build /tmp/pure-new-ivbs-perf-core-build --target kitaev_ssh_mc -j8
```

The source-only file hash list is `CORE_MANIFEST.sha256`.


<a id="report-03"></a>

## FAST_0929_README.md

### newsusy fast-0929 core snapshot

This is the accelerated source snapshot from commit `2e53aee` on branch
`optimize-deferred-phase-default`.

Included source areas are `inc/`, `kitaev_ssh/`, `src/`, `cmake/`,
`reproduction/`, `scripts/`, build files, and validation documentation.
Build products, `audit/` results, raw QMC output, and temporary files are not
included.

The current default phonon parameters are `k=1` and `omega-ph=0.5`, which
implies `mass=4` through `k=mass*omega-ph^2`. The coupling scan uses
`lambda=g^2/k`.

Current production accelerations include:

- deferred phase-on-demand magnitude updates;
- validated-factor scan reuse;
- MKL QR for large magnitude-weight workspaces;
- default local workspace and factor-cache reuse;
- cached validated-factor identity residuals.

The last two validation runs on HPC preserved byte-identical `L=6` block
outputs while reducing the short benchmark by about 1--2% relative to the
previous optimized build. Configure with the existing `CMakeLists.txt`; Eigen,
Boost and Intel MKL remain external build dependencies; PFAPACK C/Fortran
sources and licenses are bundled under `inc/pfapack/`.

#### 当前快照入口（2026-10-02 源码核对）

本目录 `/home/asus/newsusy/fast-0929` 是源码快照；上述提交和分支是包的来源标记，
不是对本目录 Git HEAD 的重新认证。后续开发工作树是 `/home/asus/newsusy/new_0923`；
其后续优化二进制与 HPC 作业不能自动归属于本快照。

- 生产入口：`kitaev_ssh/coupled_mc.cpp`，CMake 目标 `kitaev_ssh_mc`。
- CLI 默认：`mass=4`、`omega-ph=0.5`、`g=0.3`、`J=(1,1,0.5)`；
  `proposal=local`、`weight-engine=combined`、`defer-phase=1`、
  `prefix-cache=1`、`cayley-cache=1`。底层 options 的默认值可能不同；
  缓存是否实际使用取决于所选更新和权重路径。
- 固定 g 的声子路径采样、剪切周期边界及 VBS 测量已存在。
  生产 `honeycombGenerator` 捕获固定模型，未提供节点时间/协议参数；
  通用 density 驱动程序不等于 SSH 声子生产驱动已完成。
- 当前状态见 [docs/STATUS.md](REPORTS.md#report-26)，模块关系见
  [docs/CODE_MAP.md](../CODE_MAP.md)。其历史段落按原日期保留。

本次只核对源码及修正文档，未编译、运行 QMC、重跑回归或查询 HPC。
`FAST_0929_MANIFEST.sha256` 和 `CORE_MANIFEST.sha256` 保留为原包校验记录；
文档修改后不能用于验证整个当前目录。此次修改文件的校验值单独保存在
`../new_0923/audit/fast0929-doc-sync-20261002/after.sha256`，原文保存在同目录 `before/`。


<a id="report-04"></a>

## IVBS_HANDOFF.md

### iVBS 指定动量测量任务 — 第二轮完成交接（2026-09-20）

工作目录：`/home/asus/susy/pure-new`。该目录不是 Git 工作树；第二轮基线由
`/home/asus/susy/pure-state-density-work-optimized-20260918` 的
`45de6b2`（第一轮 iVBS 合并提交）重建，基线文件哈希记录在
`/home/asus/susy/tmp/round2/baseline_sha256.txt`。`pure-new` 中已有的 HPC 侧
`kitaev_ssh/phonon_mc.h` 修改属于本轮基线，未纳入第二轮增量补丁。

#### 0. 任务边界

- 只修改 VBS 指定动量的测量、候选诊断、后处理、验证和文档。
- 不修改哈密顿量、声子空间 PBC、QMC 权重、接受规则、试探态、更新算法或 RNG
  调用顺序。
- 不实现 continuous twist；任意 Q 仍是原 PBC 系综上的 finite-window estimator。
- JSON 持续写入 `"spatial_boundary":"pbc"` 和 `"twist_applied":false`，不改变
  VBS 允许动量网格。
- 不提交 HPC 作业，不运行生产规模物理计算。

#### 1. 第一轮基线状态

第一轮已提供 `--vbs-mode all|selected|both`、`--vbs-q-file`、`vbs_momenta.h`、
finite-window raw/projected 观测量、动量预测工具和 Task-B twist 可行性审查。
`all` 布局及 `both` 的 all 前缀保持原有名称和顺序。

#### 2. 第二轮修改

1. 修复 selected 模式 derived 名称和值的错位：selected 现在严格输出并映射
   `8*Nq+3` 个 derived 项；旧错位结果不可按位置迁移，必须重新生成。
2. 新增 `kitaev_ssh/vbs_candidates.h`：严格分开的 `grid_raw`、`offgrid_raw`、
   `offgrid_projected` 组；排除 Q=0；保留实部评分、虚部诊断、负值、窗口泄漏、
   +/-Q 类、自共轭类、并列集合和 finest-block 选择稳定性。
3. `coupled_mc.cpp` 将候选诊断写入独立的 `vbs.candidates` 对象；分类 winner
   index/count 的 stderr 和协方差写为 `null`，不把类别标签当连续物理量。
4. 扩展 `vbs_q_validate.cpp` 到 43 项，覆盖网格/非网格 Fourier、均匀背景、调制、
   共轭关系、布局和复相位 jackknife、q-file 错误、RNG 轨迹及候选诊断。
5. 新增 `summarize_vbs_candidates.py`：逐文件读取 `vbs.candidates`，输出 CSV/JSON；
   缺少候选对象或 finest-block 稳定性时返回 `status=unavailable`，不跨文件合并。
6. 更新 `VBS_SELECTED_MOMENTA.md` 和 `TWIST_FEASIBILITY_REVIEW.md`，修正章节、
   路径、候选 schema 及连续 twist 的理论表述。

#### 3. 已完成验证

##### 构建

```bash
export EIGEN3_INCLUDE_DIR=/home/asus/old-p/kitaev_ssh_pfqmc_hpc/third_party/eigen3
export BOOST_INCLUDE_DIR=/home/asus/.cache/pfqmc_boost_1_70_0
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export LD_LIBRARY_PATH=$MKLROOT/lib:$LD_LIBRARY_PATH
cd /home/asus/susy/pure-new
bash scripts/build_nodeps.sh /home/asus/susy/tmp/round2/final_build
```

结果：13/13 目标构建成功，退出码 0，耗时约 799.39 秒。完整日志：
`/home/asus/susy/tmp/round2/build_all.log`。

##### Validation

以下 8 个二进制均退出码 0；逐个日志保存在
`/home/asus/susy/tmp/round2/validate/`：

`kitaev_ssh_validate`、`kitaev_ssh_phonon_validate`、`kitaev_ssh_mc_validate`、
`kitaev_ssh_shifted_validate`、`kitaev_ssh_vbs_validate`、
`kitaev_ssh_vbs_q_validate`、`kitaev_ssh_fermion_validate`、
`kitaev_ssh_weight_validate`。

`kitaev_ssh_vbs_q_validate` 报告 43 项 PASS，最大报告误差约 `3.0e-14`，测试
容差为 `5e-12`。原有 P4 `kitaev_ssh_vbs_validate` 及其余验证也全部 PASS。

##### 三模式回归

```bash
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export LD_LIBRARY_PATH=$MKLROOT/lib:$LD_LIBRARY_PATH
export PYTHON=/home/asus/.venv/bin/python
bash scripts/regress_vbs_modes.sh \
  /home/asus/susy/tmp/round2/final_build \
  /home/asus/susy/tmp/round2/final_regress
```

结果：`REGRESS_VBS_MODES PASS`。L=4x4（Nc=16）acceptance `0.99599`，
L=5x3（Nc=15）acceptance `0.99606`；all/selected/both 的相位、计数和 7 个
基本观测量列逐位一致。脚本同时检查 selected/both 维数、all 前缀、候选三组、
categorical null error、`W_L_abs2` 以及 CSV/JSON 后处理。

##### 后处理异常分支

对 all-mode JSON 运行 `summarize_vbs_candidates.py` 得到 `status=unavailable`、
零 rows，未从中心值伪造频率。结果保存在
`/home/asus/susy/tmp/round2/unavailable/`。

#### 4. 交付物

- `vbs_round2_incremental.patch`：相对 pure-new 第二轮基线的增量补丁；预测工具
  的 L=4 生成物和预先存在的 HPC 侧 `phonon_mc.h` 修改排除在外。
- `kitaev_ssh/vbs_candidates.h`
- `kitaev_ssh/summarize_vbs_candidates.py`
- 更新后的 `kitaev_ssh/vbs_observables.h`、`coupled_mc.cpp`、`vbs_q_validate.cpp`
  和 `scripts/regress_vbs_modes.sh`
- 更新后的 `kitaev_ssh/VBS_SELECTED_MOMENTA.md`、`TWIST_FEASIBILITY_REVIEW.md`、
  `IVBS_HANDOFF.md`

补丁验证已完成：在由 `45de6b2` 重建并叠加基线 `phonon_mc.h` 的干净副本上，
`patch --dry-run -p1` 和实际应用均成功；9 个第二轮受影响文件逐个 `cmp` 和
SHA-256 比对一致。

#### 5. 红线和未验证项

- 没有生产规模计算或物理图；所有数值均为合成场、单元验证或小链回归证据。
- selected off-grid 峰不消除非公度失配，也不证明 iVBS 长程序；应结合
  `W_L_abs2`、raw/projected 诊断和选择稳定性解释。
- `candidate_max_*` 仍是兼容性保留的混合 raw 最大值；新分析应使用
  `vbs.candidates.groups`。
- 实空间 connected correlator 只在 all/both 网格模式存在。
- continuous twist 仍是“尚未实现 / 尚未证明”；没有 twist CLI，也没有把 off-grid
  Fourier 测量重新命名为 twisted boundary condition。

#### 6. 第三轮小修补（2026-09-20）

第三轮增量说明见 `IVBS_ROUND3_FIXES.md`，补丁为
`vbs_round3_incremental.patch`。本轮相对本交接所描述的第二轮完成树工作，
不从原始压缩包重建，也不覆盖独立的性能目录。

- `vbs.candidates` 版本为 `ivbs-candidates-v2`；中心、完整普通 jackknife
  和条件性选峰稳定性是三个独立状态，并携带原因和原始/有效/无效副本计数。
- 原始 B 个删一副本中任一副本无效时，中心值可保留，但普通 stderr/covariance
  为 `null`；有效副本至少两个时的频率标记为
  `conditional_on_valid_replicas=true`，不足两个时为 `null`。
- 后处理默认名保留参数名中的小数点；默认拒绝覆盖，`--force` 只能用于明确
  指定的输出。输入、CSV、JSON 的路径别名和 symlink 别名在写出前拒绝。
- v1 候选对象不再被当作具有新有效性语义的数据；汇总状态为 unavailable。
  中心存在但统计稳定性或普通误差缺失时为 partial，CSV 不可用数值留空。

本轮新增测试与日志：`kitaev_ssh/vbs_q_validate.cpp`、
`kitaev_ssh/test_summarize_vbs_candidates.py`、
`/home/asus/susy/tmp/round3/vbs_q_validate.log`、
`/home/asus/susy/tmp/round3/python_unit.log`。端到端示例位于
`/home/asus/susy/tmp/round3/e2e/`。没有修改 QMC 采样或物理测量定义，没有提交
HPC 作业，也没有运行生产规模计算。


<a id="report-05"></a>

## IVBS_ROUND3_FIXES.md

### iVBS 第三轮小修补（2026-09-20）

本轮基线是 `/home/asus/susy/pure-new` 的第二轮完成树。该目录不是 Git
工作树；基线逐文件哈希保存在
`/home/asus/susy/tmp/round3/baseline_sha256.txt`，第二轮补丁哈希保存在
`/home/asus/susy/tmp/round3/round2_patch_sha256.txt`。第三轮补丁
`vbs_round3_incremental.patch` 相对这棵第二轮树生成，不从原始压缩包重建。

#### 修复范围

1. `summarize_vbs_candidates.py` 只移除末尾约定的 `_results.json`，再追加
   `_candidates_summary.csv/.json`，因此 `scan_lambda_3.0_results.json` 和
   `scan_lambda_3.5_results.json` 不再碰撞。CSV/JSON 同文件、输入别名、相对
   路径和符号链接别名均在写出前检查。默认不覆盖已有文件；`--force` 只允许
   覆盖明确用 `--csv`/`--json` 指定的输出。JSON/CSV 内容先完整构造，再分别
   通过临时文件替换；这不是跨文件事务。

2. 候选诊断升级为 `ivbs-candidates-v2`。每组和每个估计量分别记录
   `center_defined/reason`、`jackknife_defined/reason` 和
   `stability_defined/reason`，并记录 `total_replicas`、`valid_replicas`、
   `invalid_replicas`。完整普通 jackknife 只在原始 B 个删一副本全部有效、
   `B>=2` 且误差有限时定义；不能用筛选后的 R 个副本代替 B。中心有效而副本
   不完整时保留中心值，stderr 和 covariance 为 JSON `null`。至少两个有效
   副本时仍可输出条件性选择频率，并设置
   `conditional_on_valid_replicas=true`；不足两个时频率为 `null`。

3. Python 汇总器严格验证 v2 的 groups、points、qpoint 索引、参与者/标签/频率
   数组、计数、有限数值和状态一致性。v1 被标记为 unavailable，不补默认有效
   状态或零频率。中心可用但误差/稳定性不可用时输出 partial；CSV 对不可用数值
   留空并保留状态及原因。结构损坏、类型错误、非有限 JSON 或 malformed JSON
   会在写出前失败，已有结果不被截断。

#### 验证产物

- `kitaev_ssh/vbs_q_validate.cpp` 新增完整 B、B=4/3 有效、0/1/2 有效、中心
  分母无效、非有限值和空 block 测试；转置关系和有限样本共轭残差的测试语义
  已分开。
- `kitaev_ssh/test_summarize_vbs_candidates.py` 覆盖小数/多输出名、默认拒绝覆盖、
  显式 `--force`、输入/同文件/符号链接冲突、partial/unavailable、v1 和错误
  频率长度。机器可读运行结果见
  `/home/asus/susy/tmp/round3/python_unit.log`。
- 真实 C++ 验证日志见 `/home/asus/susy/tmp/round3/vbs_q_validate.log`。
- 最小真实端到端样例在 `/home/asus/susy/tmp/round3/e2e/`：
  `coupled_mc` 生成 `scan_lambda_3.0_results.json`，随后生成同名参数保留的
  CSV/汇总 JSON。单元测试中的 synthetic v2 fixture 同时覆盖 partial 和
  unavailable 汇总状态。

本轮没有修改 QMC 采样内核、哈密顿量、权重、接受规则、提议分布、RNG 顺序、
raw/projected 定义、all/selected/both 布局或候选分组物理规则；没有提交 HPC
作业，也没有运行生产规模物理计算。


<a id="report-06"></a>

## MERGE_HANDOFF.md

### iVBS + Performance Merge Handoff

#### Scope and source state

Output tree:

`/home/asus/susy/pure-new-ivbs-perf-merged-20260920`

The output tree was copied from the complete working tree of
`/home/asus/susy/pure-new-faster-perf`, not from a clean checkout. Its HEAD was
`a672fcc4608307c3c324d6225e925635e4220697`, and the dirty status and binary
working-tree diff are frozen under `audit/merge-20260920/source_status/`.
No reset, checkout, clean, or restore operation was used. The source manifest
is `performance_source_manifest.sha256`.

`/home/asus/susy/pure-new` is not a Git worktree. Its source manifest is frozen
as `ivbs_source_manifest.sha256`. It already contained Round 3 before this
merge (`ivbs-candidates-v2`, schema revision 2, the Round 3 Python test, and
the Round 3 validity fields), so it was used as the Round 3 reference tree.
The Round 3 patch hash is recorded in `round3_patch.sha256` and is:

`e77ad2644a1cd173993743b086a934e99aa57ca612e47bf5c10dc20b58a24bc1`

The two input trees were not used as output workspaces. Their source manifests
match again at handoff when generated-file Python caches are excluded, and the
performance Git dirty status matches its frozen status exactly.

#### Merge decisions

- `kitaev_ssh/coupled_mc.cpp` was merged by hand. The performance
  `--weight-engine`, A1/A2/A3/deferred construction, metadata, validation, and
  `--defer-phase 1` default remain; eager mode uses `--defer-phase 0`. The iVBS `--vbs`, `--vbs-mode`,
  `--vbs-q-file`, selected/both measurement, v2 candidate JSON, and Round 3
  candidate-width fix are present. A final phonon-path and center-coordinate
  snapshot was added after sampling solely for the fixed-seed cross-regression;
  it does not change proposal, acceptance, measurement, or RNG calls.
- `kitaev_ssh/vbs_observables.h` keeps the Round 3 selected derived layout
  (`8*Nq+3`) and the existing both layout (`22*Nc+8*Nq+3`); both keeps the all
  prefix.
- `kitaev_ssh/vbs_candidates.h` and
  `kitaev_ssh/summarize_vbs_candidates.py` keep schema v2, independent center,
  ordinary-jackknife, and conditional stability states, invalid replica
  accounting, real/imag covariance, fail-closed input handling, decimal
  filename preservation, and no-overwrite/path-alias checks.
- `inc/pure_projector_fast.h`, `kitaev_ssh/kitaev_ssh_model.h`,
  `kitaev_ssh/phonon_projector.h`, and `kitaev_ssh/phonon_mc.h` remain the
  performance versions. A1 prepared trial QR, A2 static spectral reuse, A3
  geometry/gauge cache, output-buffer generation, deferred phase-on-demand, and
  cache validity/failure semantics were not replaced by iVBS files.
- `CMakeLists.txt` remains the performance CMake with MKLROOT cache identity
  handling and all 17 original CTest registrations. The iVBS validator is
  retained as one of those targets.
- `scripts/regress_vbs_deferred_modes.sh` is the new six-combination audit
  regression. It compares fixed-seed all/selected/both trajectories, final
  phonon paths, centers, seven basic observables, phase/block data, and v2
  candidate diagnostics for eager and deferred modes.

#### Verification

Configuration used the supplied CMake 4.4.3 and MKL root
`/opt/intel/oneapi/mkl/2026.1`, with the specified one-thread runtime
environment. Build directory: `build-merge-20260920`.

- Full build: passed; only existing Eigen/Boost compiler warnings were emitted.
- Incremental rebuild after the final `coupled_mc.cpp` metadata addition:
  passed.
- CTest: `17/17` passed, final log `audit/merge-20260920/test_logs/ctest_final.log`.
- `kitaev_ssh_vbs_q_validate`: passed.
- `kitaev_ssh/test_summarize_vbs_candidates.py`: passed; decimal filenames,
  v1 unavailable, partial status, undefined CSV values, and overwrite/path
  safety cases are covered.
- `scripts/regress_vbs_modes.sh`: passed for 4x4 and 5x3 fixed-seed
  all/selected/both runs, including layout, candidate, and postprocessing
  checks. The log is `audit/merge-20260920/test_logs/regress_vbs_modes.log`.
- `scripts/regress_vbs_deferred_modes.sh`: passed for all/selected/both x
  defer 0/1 on 3x2 with seed 314159. The machine-readable result and all six
  raw outputs are under `audit/merge-20260920/cross_modes/`.
- The deferred-phase CTest validator passed. Its cache-validity and cadence
  checks, together with `PhononSampler::ensureMeasured()` being called at the
  measurement stride, verify that VBS integration did not turn every proposal
  into a phase/center-Green materialization.

Every VBS integration result reports `spatial_boundary = pbc` and
`twist_applied = false`. No HPC job or production-scale physical calculation
was submitted. The short integration is a deterministic software regression,
not a claim of statistical convergence or a new speedup measurement.

#### Delta and artifacts

The actual source/document delta relative to the performance working-tree
snapshot is recorded in:

- `audit/merge-20260920/performance_snapshot_merge_delta.json`
- `audit/merge-20260920/performance_snapshot_merge_delta.patch`

The complete machine-readable merge record is:

`audit/merge-20260920/merge_manifest.json`

There are no unresolved merge or test failures. Generated Python cache files
were excluded from the final source tree; all test logs and regression outputs
remain under the audit directory.


<a id="report-07"></a>

## PERFORMANCE_HANDOFF.md

### Performance Acceptance Handoff

> 历史文档：下文的版本、目录、默认值和测试结果属于原记录时点，不代表 fast-0929 当前状态。当前快照入口见 [FAST_0929_README.md](REPORTS.md#report-03)，源码核对状态见 [docs/STATUS.md](REPORTS.md#report-26)。历史 audit 未随本包分发；本次未重跑这些测试。

Repository: `/home/asus/susy/pure-new-faster-perf`

Branch: `master`

Acceptance head: `a672fcc` plus the current uncommitted A3/deferred and
acceptance changes. The worktree is intentionally dirty; preserve unrelated
existing changes.

#### Result

The implementation is ready for handoff within the stated scope. Build and all
17 registered CTest tests pass when the MKL compile/runtime environment is
present. The combined eager/deferred HPC pair is bitwise identical under the
comparison tool and shows a 1.8365x wall-time speedup for deferred mode.

#### Rebuild and test

```bash
cd /home/asus/susy/pure-new-faster-perf
cmake=/home/asus/code/driven-new-validation-20260908/green-residual/tools/cmake/data/bin/cmake
ctest=/home/asus/code/driven-new-validation-20260908/green-residual/tools/cmake/data/bin/ctest
export CPATH=/opt/intel/oneapi/mkl/2026.1/include
export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2026.1/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE
"$cmake" -S . -B build-acceptance-20260920
gmake -C build-acceptance-20260920 -j8
"$ctest" --test-dir build-acceptance-20260920 --output-on-failure
```

The exact acceptance logs are under `audit/acceptance-20260920/`.

#### Evidence map

- `audit/acceptance-20260920/ACCEPTANCE_REPORT.md`: A-G acceptance decisions.
- `audit/acceptance-20260920/manifest.json`: machine-readable configuration,
  hashes, HPC completion states, and local evidence.
- `audit/acceptance-20260920/combined_eager_deferred_comparison.json`: exact
  eager/deferred output comparison.
- `audit/acceptance-20260920/raw_hpc/`: copied job scripts, outputs, timing,
  exit status, blocks, and binary hashes.
- `audit/acceptance-20260920/deferred_phase_validate.stdout`: detailed cache,
  cadence, and failure-recovery validation.
- `audit/acceptance-20260920/mklroot-cache-regression.log`: MKL cache isolation
  regression result.

#### Operational boundaries

- `--defer-phase 1` is the default; use `--defer-phase 0` for eager mode.
- No new HPC job was submitted, cancelled, or restarted for this handoff.
- The reference eager HPC run has only a blocks snapshot. It lacks the final
  `.time`, `.exit`, and `results.json`, so it is not used as a final timing or
  correctness result.
- HPC status was not re-queried successfully during this pass because the
  connection returned `Permission denied (publickey...)`; the copied artifacts
  remain available and hashed.
- Negative/complex phase counts are finite-sample observations. They do not
  prove positivity or the absence of a sign problem for other parameters.

#### Remaining follow-up

The next performance experiment should be a separately authorized production
design with explicit statistical goals. It should not silently make deferred
mode the default or treat the incomplete reference eager snapshot as a result.


<a id="report-08"></a>

## PERFORMANCE_OPTIMIZATION.md

### Performance Optimization Record

Status: acceptance completed on 2026-09-20. This document records the existing
A1/A2/A3 work and the deferred phase-on-demand phase. It does not introduce a
new optimization direction.

#### Scope

- A1 uses a prepared trial QR path for repeated projector evaluations.
- A2 reuses the static spectral decomposition of the immutable node geometry.
- A3 separates static lattice/coupling/gauge geometry from the changing phonon
  coordinates and provides an output-buffer generator path.
- Deferred phase-on-demand computes the certified magnitude during proposals and
  computes the Pfaffian phase and center Green only when a measurement requires
  them. It is now the command-line default; eager mode remains available with
  `--defer-phase 0`.

The relevant implementation is in `inc/pure_projector_fast.h`,
`kitaev_ssh/kitaev_ssh_model.h`, `kitaev_ssh/phonon_projector.h`,
`kitaev_ssh/phonon_mc.h`, and `kitaev_ssh/coupled_mc.cpp`.

#### Reconstruction workspace reuse

Local updates now retain the left/right reconstruction frames and factor
application workspaces between proposals. This removes repeated temporary
matrix allocations while preserving the existing QR checkpoint boundaries,
residual checks, and fail-closed recovery path. The setting is controlled by
`PureFastOptions::reuse_rebuild_workspace` and can still be disabled for an
ablation comparison.

#### HPC comparison

The combined eager/deferred comparison used the same L=4 by L=4 configuration,
seed, engine, proposal schedule, and 800 measured samples (20 blocks of 40).
The comparison tool reports zero numeric differences in both JSON results and
blocks CSV, including phase numerator/denominator and complex real/imaginary
components.

| configuration | mode | wall time (s) | result |
| --- | ---: | ---: | --- |
| combined, g=0.3, seed=8201 | eager | 2709.29 | complete |
| combined, g=0.3, seed=8201 | deferred | 1475.30 | complete |
| reference, g=0.6, seed=8301 | deferred | 1696.92 | complete |

For the paired combined run, deferred mode is `1.8365x` faster and reduces
wall time by `45.55%` (`1 - 1475.30 / 2709.29`). The reference row is a
separate engine/coupling point and is not an eager/deferred speedup pair.

Raw files, timing data, hashes, configuration, and completion states are in
`audit/acceptance-20260920/manifest.json` and its `raw_hpc/` directory.

#### Correctness evidence

- Weight replay: maximum log-weight error `1.4210854715202004e-13`, maximum
  phase error `1.6123219951621876e-15`, maximum Green error
  `4.9311717106680608e-15`.
- A3 generator validation: 21 cases, generator/`generateInto`/half-step errors
  all zero.
- Deferred replay: maximum merged log error
  `8.5265128291212022e-14`; phase and Green errors zero; both accept and reject
  transactions covered.
- Deferred measurement cadence: `measure-every=1` produced 6 measurements and
  `measure-every=3` produced 2 measurements.
- Failure injection: a measurement exception leaves the accepted snapshot and
  invalid-cache state unchanged; measurement succeeds after recovery.
- Synthetic scalar phases include two nontrivial complex-phase cases. HPC
  phase diagnostics observed 800 positive, zero negative, and zero nonreal
  samples per complete combined chain. This is finite-run evidence, not a
  positivity theorem.

#### Reproduction environment

The acceptance build used CMake 4.4.3 from the bundled path, `/usr/bin/c++`,
and MKL 2026.1. Runtime tests require:

```text
CPATH=/opt/intel/oneapi/mkl/2026.1/include
LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2026.1/lib:$LD_LIBRARY_PATH
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE
```

The first CTest invocation without the runtime library path failed during
dynamic loading; the same 17 tests passed after the path was supplied. The
failure is retained in `ctest.log`, and the passing run is
`ctest-with-mkl-runtime.log`.

The deferred path is now the CLI default after the acceptance checks above.
Use `--defer-phase 0` for the eager reference path. Existing HPC jobs retain
the binary and command line with which they were submitted.


<a id="report-09"></a>

## kitaev_ssh/FERMION_GSCAN.md

### Bounded L=4 coupling scan with matter fermions

Approved: (g, intervals) = (0,64), (.3,64), (.5,64), (.7,64), (.5,80), two independent seeds per group. L=4, J=(1,1,.5), mass=4, omega_ph=.5, k=1, dt=.1; tau_total=intervals*dt. Combined engine, bond-refresh, qr-block=4, trial shifts=0. 800 warmup + 1600 measurement sweeps, stride1, block25, vbs1, fermions1. Ten single-threaded chains. No adaptive extension.

The optional --fermions 1 flag only reads center_green at the existing measurement times, including retained states after rejection. It measures C_ij=delta_ij-G_ij for matter Majoranas (not spin correlations) and the same correlation times product of directed u along a specified path. Each channel averages over all 16 source cells, preserving source/target sublattice and periodic cell displacement. All 16 displacements and AA/AB/BA/BB are retained: 128 complex channels (bare + string). Full s and sO sums enter the existing PhaseStatistics; no rounding to signs, extra RNG calls, or Green rebuild.

Path: travel the nonnegative cell displacement along +a1, then +a2, on the initial sublattice via two-bond steps; if sublattices differ, finish with z. u_BA=-u_AB. No (-i)^length prefactor. JSON records origin-zero directed bonds and displacement; translations use identical paths with periodic wrapping. In the stored u=+1 gauge, string channels can differ from bare by deterministic path signs; this is not independent physical evidence.

Validation: before/after ctest, explicit weight and phonon tests, independent free Bloch Fourier comparison, nonuniform complex translation and gauge-invariance checks, exact 468-transaction trace with measurement off/on, and old/new production CLI field comparison. Production sampling/weight/Green/proposal files are untouched in this step.

Run directory: runs/fermion-gscan-20260917. Manifest captures binary hash, argv, seeds, compiler/MKL and source hashes. Per-chain process.json, completion.json and stderr.log retain provenance and failures. Source snapshots include new files not captured by git diff. Local WSL processes, no PBS submission.

Environment: source /opt/intel/oneapi/setvars.sh; OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE. Analysis Python: /home/asus/.cache/pfqmc_plot_env/bin/python.

Analysis: kitaev_ssh/analyze_fermion_gscan.py --run-dir runs/fermion-gscan-20260917. Uses existing convergence_analysis.joint_estimate; only equal (g,tau,other distribution parameters) are pooled. CSV records Re/Im, SE components/covariance, magnitude and jackknife magnitude SE. NPZ retains joint covariance, source blocks retain raw sO. Blocks 25/50/100/200 and first/second halves screen uncertainty and drift; neither proves equilibrium. g=0 compared with independent analytic Fourier values. Denominator failure retains phase-only diagnostics, not fabricated physical ratios.

No sign observed is not a strict positivity proof. No L expansion, dt extrapolation or new Pfaffian optimization. Numerical failures stop the affected chain and are not treated as rejections. Completion analysis never restarts sampling.


<a id="report-10"></a>

## kitaev_ssh/L4_ED_FEASIBILITY.md

### L=4 ED feasibility and exact matter benchmark

#### Scope

The existing `kitaev_ssh_free` calculation is an exact Gaussian ED of the
fixed-gauge quadratic matter Hamiltonian. It diagonalizes the 32x32 Majorana
generator for `Lx=Ly=4` (16 unit cells, 32 Majoranas, 16 complex matter
modes). It is not an exact diagonalization of the coupled quantum
electron--phonon Hamiltonian.

#### Reproducible L=4 matter result

For the same zero-flux/Wilson-(+1,+1) sector used by the scan, with
`J=(1,1,0.5)` and `X=0`:

```text
build/kitaev_ssh_free --Lx 4 --Ly 4 --Jx 1 --Jy 1 --Jz .5 \
  --tau-per-branch .4 --dt .05 --block 4
```

The exact matter ground-state energy is

```text
E0 = -21.806344131545224
E0 / N_site = -0.68144825411078824
```

The independent Gaussian/projector consistency residuals are
`green_residual=2.56e-15` and `solve_residual=9.21e-16`.

As a same-size cross-check, the existing `lambda=0`, `tau=8` L=4 QMC
records `majorana_energy=-21.806344131545206 +/- 5.22e-14`; its difference
from the ED value is about `1.8e-14`. This validates the L=4 matter sector at
zero coupling, but it does not validate the interacting phonon structure
factor.

#### Why full L=4 electron--phonon ED is not practical

At `L=4` there are 48 bond phonon modes and 16 complex matter modes. A
per-mode phonon cutoff `nmax` gives the Hilbert-space dimension

```text
2^16 * (nmax+1)^48.
```

Even `nmax=1` is `2^64 ~= 1.84e19` states; `nmax=2` is about `5.2e27`.
No dense or ordinary sparse ED optimization turns this into a practical
exact L=4 calculation. A truncated/selected-mode calculation would be an
approximation and must be labeled as such.

#### Efficient options

1. For fixed phonon coordinates, use the 32x32 Gaussian ED above (exact and
   essentially instantaneous).
2. If an explicit matter Fock spectrum is required, use parity-resolved,
   matrix-free Lanczos in the 2^16 matter basis rather than a dense
   65536x65536 matrix.
3. For the full coupled problem, keep the QMC formulation, or restrict ED to
   L=2 / a deliberately small collective-phonon truncation.


<a id="report-11"></a>

## kitaev_ssh/L4_MIXING.md

### L=4 mixing and stationarity batch

PBS job 155149.mgt completed on port 20015 with exit status 0. All four
single-threaded chains (seeds 7801--7804) exited normally. The fixed parameters
and budget are recorded in `l4_mixing_plan.json`: Lx=Ly=4, tau=6.4, dt=.1,
J=(1,1,.5), g=.3, mass=1.7, omega_ph=.8, inherited trial shift
d=(.1378676470588235,.1378676470588235,.1378676470588235), 800 warmup and
1600 measured sweeps per chain, stride 1 and raw blocks of 25 sweeps.

The PBS wall time was 1:53:40 and CPU time 7:31:13. Individual elapsed times
were 6761.87--6775.96 seconds including warmup. Acceptance ranged from
.81350 to .81461. All chains had zero numerical failures, resolved phase
1+O(1e-17)i, and CSV/JSON ratio discrepancies no larger than 1.43e-14.

The chain-preserving joint complex jackknife at the planned 100-sweep block
size gives:

| observable | estimate | SE |
|---|---:|---:|
| Xz | .106135 | .002411 |
| uniform Btrace | .982904 | .019787 |
| Szz(0,2) | .382047 | .006587 |
| Tr S(0,2) | 1.115469 | .012618 |
| Szz(2,2) | .366669 | .006645 |
| Tr S(2,2) | 1.116216 | .011183 |

For 50/100/200-sweep blocks the SE envelopes (maximum/minimum) are 1.032,
1.090, 1.190, 1.256, 1.087 and 1.227 in the table order. The 50 and 100
layers are close. The 200 layer has only eight blocks per chain and usually
smaller errors, so it is a sensitivity diagnostic rather than evidence of an
error plateau.

Time stationarity fails in the q1 channel. Treating each chain's paired time
difference as one independent repeat and using a Student-t interval with
three degrees of freedom gives:

| observable | second half - first half | t95 interval |
|---|---:|---:|
| Xz | -.00270 | [-.02243,.01703] |
| uniform Btrace | -.02486 | [-.10622,.05650] |
| Szz(0,2) | +.03774 | [.02696,.04853] |
| Tr S(0,2) | +.07991 | [.04133,.11848] |
| Szz(2,2) | -.01060 | [-.06326,.04206] |
| Tr S(2,2) | +.00548 | [-.05122,.06218] |

The independent first/last-quarter check agrees: Szz(0,2) changes by
+.04064 with t95 [.03123,.05005], and Tr S(0,2) by +.07733 with t95
[.03322,.12144]. All four chains have positive q1 changes. Quarter means of
Szz(0,2) are .36056, .36579, .40063 and .40120; the common shift occurs
between the second and third quarters. Pairwise full-chain means remain
mutually consistent (largest planned-observable pairwise |z| is 1.84), so the
problem is a common time drift rather than one visibly discrepant chain.

This batch therefore does not establish stationarity for the q1 structure
factor channel and is not ready for a production L=4 projection comparison.
The most informative next calculation is a fresh, fixed-budget batch with a
longer warmup while retaining tau, dt, trial shift, estimator definitions and
the now-validated combined weight engine. It should repeat the same paired
half/quarter q1 diagnostic. No VBS conclusion follows from this batch.

Raw data and analysis are in `build-kitaev-p01/l4_mixing/`. Reproduce the
joint analysis without plots with:

```bash
bash kitaev_ssh/analyze_l4_mixing.sh
```


<a id="report-12"></a>

## kitaev_ssh/L4_PERFORMANCE.md

### L=4 performance and joint analysis

One new short chain, seed 7701, 30 sweeps (1440 whole-bond proposals), tau=6.4, h=.1, 64 intervals, J=(1,1,.5), g=.3, mass=1.7, omega=.8, fixed inherited trial shifts .1378676470588235. Initialization and measurement/output are outside the timed region. Each math library and chain uses one thread. This is a performance trajectory, not equilibrium data. No images or long chains were generated.

| operation | seconds | fraction of 83.7674 s |
|---|---:|---:|
| physical signed Pfaffian weight routine |69.2522|82.672%|
| factor generation and storage |11.0359|13.174%|
| center Green endpoint rebuild |2.80084|3.344%|
| factor validation including temporary matrix copy |.50884|.607%|
| proposal path copy |.00654|.0078%|
| reversed factor list copy |.00416|.0050%|
| snapshot path copy |.00851|.0102%|
| remaining proposal/action/commit/timer overhead |.15046|.180%|

Each measured operation ran 1440 times. The top-level categories do not overlap. Numerical-core internal matrix allocation/copy costs remain charged to the Pfaffian/Green routine; this is not a global allocator profile. Factor generation includes model construction, generator checks, matrix exponentials and factor storage. The signed weight routine includes its QR/SVD, intermediate Green functions, LU transforms and Pfaffian sign/phase tracking, not just the bare Pfaffian kernel.

The main cost is the complete signed overlap propagation through 128 factors for every bond proposal. Intermediate subspaces change between equal adjacent half-steps, so neither their overlaps nor phase decisions can simply be reused. Node exponentials are already generated once per node. Explicit external copies are negligible. A possible future cache of slice-only Cayley transforms would need a deeper profile and a validated immutable-factor cache; the current evidence does not quantify its benefit. No local optimization was implemented and no speedup is claimed. All physical phase and endpoint consistency checks remain active. Instrumentation is enabled only by KITAEV_SSH_PROFILE in the profiling executable; the production code compiles through direct calls.

Measured throughput: 2.79225 s/sweep, acceptance .834028, final matter phase 1-8.29e-16 i. The short profile agrees with the earlier 2.77–2.90 s/sweep pilot. It does not establish absence of failures on long contours or equilibrium sampling.

HPC: job 155147 failed at compilation due to an unavailable C.UTF-8 locale and performed no sweeps. Setting LC_ALL=C allowed job 155148.mgt to complete (exit 0). Its PBS walltime was 2m16s and CPU time 2m50s including two concurrent compilations and the affected MC validation. The source directory is `/home/sunxr/kitaev_l4_profile.AtQTige5`; full build/run command is `qsub kitaev_ssh/l4_profile.pbs`. Local target: `kitaev_ssh_profile`. Logs and profile.json are in `build-kitaev-p01/l4_profile/`.

The affected existing mc_validate passed, including Fock acceptance ratios, accept/reject transaction checks, phase accumulation and free/coupled single-bond checks. The existing Python analysis tests plus L=4 normalization/momentum tests passed 15/15. No full historical regression suite was rerun.

#### Joint analysis

The existing chain-weighted complex-ratio jackknife now supports even grids including L=4. Historical q1/q3 labels denote fixed reciprocal fractions (0,1/2) and (1/2,1/2). Output `fixed_momentum_aliases` maps these to actual integers (0,2)/(2,2), flattened indices 2/10 at L=4. Grid/Nc mismatches and differing grid order across chains are rejected. Odd grids lacking the prescribed momenta are rejected rather than rounded.

Nc is Lx*Ly=16, including the uniform background Nc*sum_mu(<X_mu>^2) and its correlated jackknife error. Input S already contains 1/Nc and is not divided again. Other-nonzero-momentum backgrounds now average all Nc-2 eligible momenta rather than the two L=2 points. All 16x9 complex S entries have primary-block joint estimates, real/imag errors and covariance in `full_grid_estimates`; original complete phase/numerator blocks remain untouched. Historical L=2 behavior is retained. Images now require explicit `--plots`.

The existing pilot's four chains were analyzed independently of L=2, primary block 5, sensitivity 10/20/25 sweeps (20/10/5/4 blocks per chain). Primary interface-check results:

| quantity | mean | joint SE |
|---|---:|---:|
| Xz |.114033|.007388|
| uniform Btrace |1.023594|.067534|
| Szz (0,2) |.380350|.029541|
| Tr S (0,2) |1.119433|.054052|
| Szz (2,2) |.343968|.026201|
| Tr S (2,2) |1.067454|.045218|

These errors are short-block diagnostics only. The 50-sweep warmup and available block counts cannot establish thermalization or an error plateau. Reproduce, choosing an absent output directory:

```bash
MPLCONFIGDIR=/tmp/kitaev-profile-mpl /home/asus/.venv/bin/python kitaev_ssh/convergence_analysis.py \
 --group L4_pilot=build-kitaev-p01/l4_pilot/results/l4_7601_results.json,build-kitaev-p01/l4_pilot/results/l4_7602_results.json,build-kitaev-p01/l4_pilot/results/l4_7603_results.json,build-kitaev-p01/l4_pilot/results/l4_7604_results.json \
 --block-sweeps 5 --output /tmp/kitaev-l4-joint-reproduction
```

Saved output is `build-kitaev-p01/l4_pilot/joint_analysis/`, containing groups.json, metadata, half/quarter diagnostics, all without plots.

#### Proposed next budget (not submitted)

Four fresh independent chains, each 800 warmup plus 1600 measured sweeps, raw block 25, stride 1. At 2.79 s/sweep this is about 1.86 hours/chain or 7.45 CPU-hours, with four concurrent single-threaded processes; request a 2.5-hour walltime allowance. Keep the present physics and trial shifts fixed. Primary block 100 with 50/200 sensitivity; examine chain means, half/quarter evolution, phase and numerical failures. The 200-sweep layer has only eight blocks/chain and is diagnostic, not an error plateau certificate. Lock this budget before observing results and do not extend based on significance. This mixing pilot precedes any expensive L=4 projection comparison; it does not aim to attain the L=2 .02/.04 equivalence precision.


<a id="report-13"></a>

## kitaev_ssh/L4_PILOT.md

### L=4 bounded pilot and L=2 error budget

The existing L=2 checks support the following measured changes. Intervals are the existing four-observable corrected, blocking-sensitivity intervals; they are not bounds on the remaining infinite-projection bias.

| Observable | tau 3.2 to 6.4, h=.2 | h=.2 to .1, tau=6.4 | target |
|---|---|---|---|
| Szz (0,1/2) | [-.006766,.014037] | [-.014793,.011918] | +/-.02 |
| Tr S (0,1/2) | [-.016906,.019157] | [-.022004,.020502] | +/-.04 |
| Szz (1/2,1/2) | [-.006270,.014682] | [-.015405,.008558] | +/-.02 |
| Tr S (1/2,1/2) | [-.008110,.029937] | [-.023004,.021250] | +/-.04 |

Both measured changes satisfy their separate criteria. These intervals must not be interpreted as an additive rigorous error budget: the two comparisons share the coarse tau=6.4 reference, and the projection check was performed only at h=.2. Fine-step residual estimates require the explicitly conditional second-order assumption in T64_STEP_HALVING.md. Seed 7503's quarter drift remains a stationarity limitation.

The next fixed budget is stored in l4_pilot_plan.json, with full commands in l4_pilot.pbs: four independent seeds 7601–7604, 50 warmup and 100 measurement sweeps each, raw blocks 5, tau=6.4, h=.1, intervals=64. Couplings, mass, frequency and trial shift are fixed as recorded there. This is only a timing, acceptance, phase and numerical stability pilot. It cannot establish equilibrium, size convergence or VBS. No sampling extension is authorized by the plan.

HPC submission: `qsub l4_pilot.pbs` in `/home/sunxr/kitaev_l4_pilot.tjdwUTsw`, job `155146.mgt`, four single-threaded chains. The validated HPC binary is reused with the SHA256 pinned in the PBS script. It writes legacy imaginary-time metadata, which the current reader converts to tau.

L=4 fixed physical momenta are integer (0,2) and (2,2), corresponding to reciprocal fractions (0,1/2) and (1/2,1/2). The existing P5 fixed-q multi-chain analyzer explicitly supports L=2 only, so it must not be applied to these L=4 files; the existing full-grid P4 reader supports the pilot outputs.

#### Completed pilot

PBS job 155146.mgt completed with exit 0, walltime 7m19s and CPU time 28m38s. Every chain completed 100 measurements and 20 raw blocks; numerical failures were zero. Mean phases were 1 with imaginary parts below 1.4e-16. CSV/JSON discrepancies were at most 1.07e-14.

| seed | elapsed seconds | acceptance | Xz | Szz(0,2) | Szz(2,2) |
|---|---:|---:|---:|---:|---:|
|7601|434.47|.81542|.13086|.34717|.36263|
|7602|432.84|.82104|.09908|.45633|.35343|
|7603|433.34|.81646|.09535|.37712|.31353|
|7604|415.58|.81896|.13084|.34078|.34629|

These are short-run diagnostic means. At raw block 5, individual-chain Szz standard errors span .0376–.0769 and Tr S errors .0714–.1247; neither thermalization nor a block-error plateau is established. The trial shift was inherited from L=2 and is not asserted to be optimal for L=4.

Measured cost is 2.77–2.90 seconds per sweep including amortized output overhead, about 57 times the previous L=2 fine-step cost per sweep. A chain of 3200 warmup plus 8000 measured sweeps would take approximately 8.6–9.0 hours if this throughput persists; four parallel chains would consume about 35–36 CPU-hours. This extrapolates cost, not attainable error.

The next informative calculation is a separately budgeted L=4 mixing pilot with substantially longer warmup and enough measurements for within-chain blocks, before an expensive projection comparison. No further job was submitted. Raw blocks, result JSON, PBS output and exit file are in `build-kitaev-p01/l4_pilot/`. Reproduce the numerical summary using `bash kitaev_ssh/analyze_l4_pilot.sh`; it writes summary.csv with all saved blocking levels. No images were generated, following the user's preference.


<a id="report-14"></a>

## kitaev_ssh/PHASE_TAU_SCAN.md

### L=4 tau 扫描：预定方案与验收

本轮只测完整权重的平均符号/复 phase 随总投影时间的变化，不改 Hamiltonian、weight、Green、proposal、RNG 消费顺序，不扩大 L，不做 Pfaffian 性能优化。

#### 固定输入

Lx=Ly=4，J=(1,1,0.5)，用户已确认 g=0.5；m=4，omega_ph=0.5，k=m*omega_ph²=1，lambda=0.25。dtau=0.1，combined，bond-refresh，trial shifts=(0,0,0)，固定 u=+1 / Wilson(+1,+1)，free X=0 matter trial，qr-block=4，VBS=1。

只通过 --intervals 8,16,32,48,64,80 和 --dt 0.1 实现 tau_total=0.8,1.6,3.2,4.8,6.4,8.0；不存在新 --tau/--dtau CLI。

所有链统一 200 warmup + 800 measurement sweeps；stride=1；原始block=25；额外50/100块长度敏感性检查。不根据观测显著性延长或提前结束。每个tau两个不同seed，跨tau也不复用seed。计划完整保存在 phase_tau_scan_plan.json。

12条本地WSL链并行，数学库每链1线程。不是HPC/PBS提交，故无HPC作业号。每链argv、PID、时间、exit code保存在run目录，不覆盖旧结果。

#### 最小诊断改动

- phase_statistics.h：PhaseBlock/PhaseSummary增加positive/negative/nonreal测量计数，mergeAdjacentBlocks保持计数守恒；pending块清零恢复默认0。仍先验证单位复相位，保留原phase及phase*O求和。
- 分类阈值abs(Im phase)>1e-10时为明显非实；否则Re<0为负，其余为正。互斥分类，仅production测量计数，不计warmup或接受提案次数。
- 复均值SE=hypot(SE_Re,SE_Im)，SNR=abs(mean phase)/SE；保留Re/Im协方差。模长SE使用原jackknifeErrors和相同删块replica，未重写统计框架。
- coupled_mc.cpp只新增输出helper和JSON字段，原命令行、采样和Green调用原样保留。SE=0时SNR输出null并标注zero_observed_variance，避免非法JSON infinity。
- mc_validate.cpp增加正/负/非实分类、计数合并守恒、非rounding和模长误差验证。

#### 轨迹与测试验收

补丁前后13/13 ctest通过，另显式执行weight_validate、phonon_validate。weight输出逐字节相同，phonon日志逐字节相同。

隔离test-only sampler observer记录每事务accept、log ratio、完整log weight、phase，以及声子场/Green/RNG engine/normal_distribution cache/shuffle order摘要。L4 intervals=8/80各3 sweeps的bond-refresh，加L2小系统3 sweeps local，共468事务；补丁前后trace逐字节一致。观察器只存在于runs/.../validation/trace_before和trace_after，不进入生产代码。

真实生产CLI对照：原二进制mc_before与新binary，同输入L4 intervals8，seed940008，warmup1+sweeps4，所有旧JSON字段（耗时除外）和旧CSV字段严格相等，新增计数守恒。

回归误差：log weight=1.1368683772161603e-13；phase=2.0466243424411105e-15；Green=0；log ratio=1.4210854715202004e-13；接受概率=1.1246559239452836e-13。weight_validate的48/24是强制事务覆盖，不是自然MH接受率。

#### 运行与输出

```bash
cd /home/asus/susy/pure-state-density-work
source /opt/intel/oneapi/setvars.sh
python3 -u kitaev_ssh/run_phase_tau_scan.py --run-dir runs/phase-tau-20260917
```

上面是已提交命令的记录，不要对同目录重跑；runner会拒绝覆盖manifest。运行前有ctest、weight、phonon、trace、CLI对照门槛。扫描结束后自动调用analyze_phase_tau_scan.py，后者复用convergence_analysis.read_chain/joint_estimate/estimates_json，不混合不同tau，不跨独立链合并block。

最终输出：results/phase-tau-20260917/PHASE_TAU_REPORT.md、per_chain.csv、pooled_by_tau.csv、blocking_sensitivity.csv、stationarity.csv、phase_tau_summary.json。

#### 解释边界

“未观测到负号”只表示这些已记录的有限采样中negative_count=0；不等同于“严格证明无符号”。复相位虚部可能仅有浮点噪声，此时很大的分母SNR和零实部误差不能排除稀有负号或未访问的构型区。尚未建立平稳性、可靠tau_int、dtau收敛、投影收敛或尺寸收敛。半链比较和block敏感性是诊断，不是严格平稳性证明。

源commit：4c02dc60e1ffcd5d142d34f588f15ced5424eee2（有诊断补丁）；binary SHA256：ee5eff9ef63c25562119828a8f6549fe6bb4adb018668d955f31cfd50762cc15。
编译器：g++ (Ubuntu 15.2.0-16ubuntu1) 15.2.0；MKL：/opt/intel/oneapi/mkl/2026.1。


<a id="report-15"></a>

## kitaev_ssh/SHIFTED_LONG.md

### Shifted-trial bounded long-chain comparison (2026-09-16)

This batch used the previously validated shifted-trial binary, unchanged physics and an advance-fixed budget. Parameters: Lx=Ly=2, J=(1,1,0.5), g=0.3, mass=1.7, omega_ph=0.8, zero-flux Wilson (++), equal endpoint trial shifts d=(0.1378676470588235,0.1378676470588235,0.1378676470588235), dt=0.2, free whole-bond refresh. T=1.6/M=8 uses seeds 6901,6902; T=3.2/M=16 uses seeds 7001,7002. Each chain: 800 warmup sweeps, 8000 measured sweeps, stride 1, raw blocks of 25 sweeps. Four CPU slots; each process and math library limited to one thread. Planned parameters and binary hash are in `shifted_long_plan.json`; precise command lines are in `shifted_long.pbs`. Submitted as PBS 155139.mgt on cn1; exit code 0 for job and all chains; walltime 00:01:56, CPU time 00:05:48. Per-chain sampling elapsed times were 60.43, 60.23, 112.56, 112.90 s. Acceptances: 0.9332, 0.9335, 0.8869, 0.8855; zero numerical failures; mean complex phase 1+O(10^-18)i; CSV/JSON observable ratio discrepancies <8e-15.

New long chains **alone** determine the main results below, using chain-preserving joint phase/numerator jackknife with 200-sweep primary blocks (40 blocks per chain). Old 1600-measurement shifted chains appear only in `history_comparison/` as separate groups. All four fixed q observables retain their predeclared |difference| targets (.02 for Szz; .04 for Tr S). `Btrace=4*sum_mu(<X_mu>^2)` is the uniform phonon background, not a Majorana bond expectation.

| Observable | T=1.6 mean ± SE | T=3.2 mean ± SE | T3.2−T1.6 ± SE | marginal 95% interval | four-q/blocking envelope 95% (where relevant) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Xz | .136796 ± .003010 | .132205 ± .002817 | −.004591 ± .004123 | [−.012672,.003490] | — |
| Btrace | .221583 ± .005702 | .221892 ± .005341 | +.000310 ± .007813 | [−.015003,.015622] | — |
| Szz q1=(0,1) | .376138 ± .003790 | .384299 ± .005175 | +.008161 ± .006414 | [−.004411,.020733] | [−.009019,.025341] |
| Tr S q1 | 1.119034 ± .008008 | 1.129856 ± .007398 | +.010821 ± .010902 | [−.010547,.032190] | [−.017423,.039066] |
| Szz q3=(1,1) | .366490 ± .004520 | .370110 ± .004783 | +.003620 ± .006581 | [−.009279,.016520] | [−.013331,.020572] |
| Tr S q3 | 1.122089 ± .007578 | 1.129173 ± .008152 | +.007084 ± .011130 | [−.014730,.028898] | [−.022745,.036913] |

The last interval uses the worst eligible SE among the 100/200/400-sweep levels and the existing four-q Bonferroni safeguard. Thus Tr S at both q satisfies the **finite T-pair** .04 tolerance by this diagnostic, but Szz at both q fails the complete-interval .02 criterion (q3 exceeds by only .000572). These are approximate intervals, not infinite-projection or VBS evidence. The new Szz(q1) drift +.008161 is smaller than the short-chain +.021622; their difference −.013461 has independent-chain SE about .01373, so the earlier drift is not independently contradicted or established.

Block SE at 100/200/400 sweeps (each chain has 80/40/20 blocks at these levels):

| Observable | T=1.6 | T=3.2 |
| --- | --- | --- |
| Xz | .002749 / .003010 / .002896 | .002866 / .002817 / .002695 |
| Btrace | .004950 / .005702 / .005478 | .005696 / .005341 / .005342 |
| Szz q1 | .004230 / .003790 / .003385 | .004823 / .005175 / .004209 |
| Tr S q1 | .007372 / .008008 / .008535 | .007181 / .007398 / .006045 |
| Szz q3 | .004324 / .004520 / .004224 | .004837 / .004783 / .005063 |
| Tr S q3 | .007632 / .007578 / .006604 | .008437 / .008152 / .009054 |

Blocking appears reasonably level over this range, though stationarity remains a limitation. T=1.6 chains differ in Szz(q1) by −.01570 ± .00758 (marginal 95% [−.03055,−.00085], but four-observable simultaneous interval crosses zero). T=3.2 seed 7001 shows second-half Xz−first-half +.01779 ± .00739; seed 7002 shows second-half Tr S(q3)−first-half −.06547 ± .01943. These are screens over many comparisons; the half-chain SE neglects covariance across its shared boundary, and independent-chain agreement cannot prove mixing. Old short-chain values are compatible with new ones within their large independent-chain error; the short T=3.2 Btrace .24448±.01375 versus the long .22189±.00534 differs −.02258±.01476 (95% [−.05150,.00634]), so short chains are not used to tighten the primary intervals.

Reproduce analysis locally from the raw `build-kitaev-p01/shifted_long/results/` files with `P5_PYTHON=/home/asus/.cache/pfqmc_plot_env/bin/python bash kitaev_ssh/analyze_shifted_long.sh` from the repository root. Original blocks, JSON, stdout/stderr and PBS stdout/exit file are in `build-kitaev-p01/shifted_long/`; main CSV, grouped jackknife, blocking and fixed-q/background/difference plots in `build-kitaev-p01/shifted_long/analysis/`; separate old-chain checks in `build-kitaev-p01/shifted_long/history_comparison/`. The only statistical code adjustment was adding Btrace to existing per-chain half diagnostics; its targeted regression and all 12 affected statistical tests passed.

Most informative next calculation: first resolve the T=3.2 half-chain stationarity screen using prespecified longer warmup and additional **independent** chains at the same d,T,h, without merging old short chains to force tighter errors. Then, if stable, compare a single shifted T=6.4 point against T=3.2 for fixed q1/q3 Szz; this batch did not add any samples or T points after inspecting results.


<a id="report-16"></a>

## kitaev_ssh/SHIFTED_T32_STATIONARITY.md

### T=3.2 shifted-trial chain stationarity screen (2026-09-16)

Cross-chain supplement (existing data only): joint phase/numerator pooling with 200-sweep blocks gives Szz(q1) Q4−Q1 = +.024871 ± .009418, normal approximate 95% [.006411,.043331], nominal two-sided p=.008273. Leave out seeds 7101/7102/7103/7104 respectively: +.030896±.010208 [.010889,.050903]; +.025644±.010651 [.004769,.046520]; +.016070±.011410 [−.006294,.038434]; +.026875±.011192 [.004938,.048812]. All chains remain in the formal dataset. These are chain-stratified joint ratio jackknives, keeping the original screen's omission of same-chain cross-quarter covariance; quarters are separated by 4000 measured sweeps. The individual contrasts reproduce the original JSON; equal allocation and numerically unit phases also recover the pooled value/SE from independent chain contributions to <3e-17. Seed 7103 supplies 51.5% of the summed signed contrast; omitting it reduces the pooled point estimate 35.4%, retains a positive point estimate, and makes its marginal interval cross zero. Positive direction is shared, while strength of evidence is sensitive to 7103.

Multiplicity supplement: the six prespecified observables generate 24 half contrasts (4×6), 24 Q4−Q1 contrasts (4×6), and 36 between-chain contrasts (6 pairs×6), total 84. Same-direction/repeated-trigger rules reuse those contrasts; six 400/200 SE ratios, phase/numerical guards and 96 descriptive quarter means do not supply additional calibrated p-values. The general script also saves 70 nonprimary contrasts (4 X2 half contrasts and 66 nonprimary between-chain contrasts); they are outside the prespecified readiness decision. Normal two-sided p-values with Holm adjustment within each family have minima .789861 (halves), .026824 (quarters: seed 7103 Szz(q1)), and .485290 (chain pairs). The smallest unadjusted p=.00111768 remains significant within the 24 quarter comparisons, but Holm/Bonferroni across the supplemental 84-comparison union both give .093885, with no rejection at .05. The original six-observable correction for each chain pair did not cover the entire six-pair or 84-comparison family. The 84-union correction is a transparent supplementary scope, not an alteration of the original alert rule; neither its nonsignificance nor the retrospective pooled/leave-one-out intervals proves stationarity. Normal approximation and only ten blocks per quarter limit the p-values' accuracy.

Independent replication remains useful to distinguish reproducible cross-chain positive Szz(q1) early/late change from random trajectory fluctuations whose evidence is amplified by 7103. If commissioned, the pooled cross-chain contrast should be specified as the primary replication statistic before seeing new data; no such sampling was added here. Reproduce this supplement using `MPLCONFIGDIR=/tmp/kitaev-p5-mpl /home/asus/.cache/pfqmc_plot_env/bin/python kitaev_ssh/summarize_stationarity_quarters.py`. Exact individual/pooled/leave-one-out values and all p-value adjustments are saved in `build-kitaev-p01/shifted_t32_stationarity/analysis/quarter_summary_multiplicity.json`.

**Practical decision: not yet ready to launch T=3.2 versus T=6.4.** The prespecified *quarter* trigger fired in independent seed 7103 for Szz(q1); none of the new chains reproduced old seed 7002's Tr S(q3) half-chain decline. An isolated quarter trigger, among many screens, does not prove nonstationarity or a model/estimator error. It does prevent a claim of sufficient stationarity under the rules fixed in `shifted_t32_stationarity_plan.json` before any of these results were viewed. No new samples, additional T values, altered thresholds, Hamiltonian, trial state or estimator were introduced after that trigger.

Four newly seeded chains: 7101,7102,7103,7104. Fixed L=2, J=(1,1,.5), g=.3, mass=1.7, omega_ph=.8, Wilson (++), d=(.1378676470588235,.1378676470588235,.1378676470588235), intervals=16, dt=.2, bond-refresh. Each has 3200 thermalization sweeps followed by 8000 measured sweeps at stride 1, raw blocks of 25. The exact invocation is in `shifted_t32_stationarity.pbs`; the binary SHA256 was `75440fcce5f9781c2f89eb776d78d4209cc46509b07cdf8def34727d1be2edf4`. To submit on the original HPC account: `cd /home/sunxr/kitaev_t32_stationarity_20260916.zNKJNeyd && /opt/pbs/bin/qsub shifted_t32_stationarity.pbs`. PBS job **155140.mgt** ran four single-threaded chains on cn1 with four CPU slots. All four exited normally, all 320 raw blocks/chain are present, zero numerical failures; job exit 0, measured PBS walltime 2m44s and CPU time 10m36s. Chain elapsed sampling times 157.55,157.43,160.13,161.97 seconds; acceptance fractions .88406,.88543,.88473,.88471. The retained physical complex average phase is 1+O(10^-18)i (numerical phase error only); CSV/JSON reweighting ratio differences <7e-15.

All results below use **only the four new chains** for the group mean and a primary 200-sweep block, which leaves 40 blocks per chain. `Btrace=4*sum_mu(<X_mu>^2)` is the uniform phonon background, not the Majorana bond operator. The sequential quarter bounds are measured-sweep indices [0,2000), [2000,4000), [4000,6000), [6000,8000); each quarter contains ten independent-of-other-chains 200-sweep blocks. Quarter/half comparisons are screening diagnostics; their quoted difference SEs omit covariance between portions of the *same* chain and their many comparisons are not a calibrated familywise test.

| Observable | Four-chain mean ± SE (200) | SE at 100 / 200 / 400 |
| --- | ---: | --- |
| Xz | .139659 ± .001987 | .001947 / .001987 / .001977 |
| Btrace | .231780 ± .003591 | .003640 / .003591 / .003936 |
| Szz q1=(0,1) | .377620 ± .003490 | .003270 / .003490 / .003642 |
| Tr S q1 | 1.130932 ± .005649 | .005344 / .005649 / .005776 |
| Szz q3=(1,1) | .370173 ± .003213 | .002988 / .003213 / .003529 |
| Tr S q3 | 1.127949 ± .006002 | .005877 / .006002 / .005761 |

There are 80/40/20 blocks per chain at the 100/200/400 layers; the largest 400/200 SE ratio across the six observables is 1.10, so the prespecified blocking trigger (>1.5 for at least two) did not fire. Across the six pairs of independent chains, the largest standardized difference is 2.47 for Szz(q1); no targeted quantity in any pair has a six-observable simultaneous 95% interval excluding zero. Neither result establishes mixing on its own.

New-chain Tr S(q3) second-half minus first-half values, ordered by seeds 7101–7104: −.02525±.02522, −.00543±.02520, +.01240±.02205, +.01303±.02415. Their fourth-minus-first-quarter differences are −.02646±.04650, +.04176±.03430, −.00765±.02956, +.03258±.03349. Old *separate*, 800-warmup seed 7002 gave −.06547±.01943 in halves and −.10609±.02511 in first/fourth quarters. This old effect does **not** recur in four new chains, though four null screens cannot prove the old trajectory was harmless.

The specified quarter test nevertheless triggers at **seed 7103, Szz(q1)**: first quarter .34303±.01273 and fourth .39430±.00925, difference +.05128±.01573 (3.26 diagnostic SE). Its second and third quarters were .35997±.01706 and .36195±.01215. Szz(q1) fourth-minus-first was positive in *all* four new chains (+.00680±.02194, +.02255±.01996, +.05128±.01573, +.01886±.01709), but only 7103 exceeds two SE; intermediate quarters are not uniformly monotonic. None of the 24 new-chain half diagnostics exceeds 3 SE, and the half differences above two SE in Tr S(q1) have opposite signs (seeds 7103 and 7104). Thus the **specific failed readiness condition** is one quarter contrast, not reproduced seed-7002 Tr S(q3) drift, a widespread half-chain trend, blocking growth, phase instability, or numerical failure. Because that contrast was among 24 prespecified first/last diagnostics with only ten blocks/quarter, classify it as *unresolved stationarity* rather than a proven drift; do not waive the prespecified threshold after viewing it.

Full complex means, joint block jackknife covariance, every chain half and quarter contrast, plus blocking and between-chain checks: `build-kitaev-p01/shifted_t32_stationarity/analysis/groups.json`. Direct tables: `halves.csv`, `quarters.csv`, `quarter_first_last.csv`, `metadata.csv` in that analysis directory; `quarters_NEW.png` and `blocking.png` show the trajectories. Old seed 7002 is analyzed separately in `old_seed7002_diagnostic/`, never merged with the four new longer-warmup chains. Original `_blocks.csv`, `_results.json`, per-chain logs and PBS stdout/exit code: `build-kitaev-p01/shifted_t32_stationarity/`. Reproduce from repository root using `STAT_ANALYSIS_DIR=/tmp/kitaev-t32-repro-new STAT_OLD_ANALYSIS_DIR=/tmp/kitaev-t32-repro-old P5_PYTHON=/home/asus/.cache/pfqmc_plot_env/bin/python bash kitaev_ssh/analyze_shifted_t32_stationarity.sh` (choose absent, task-specific output paths to preserve original analyses). A run into a separate `/tmp/kitaev-t32-repro.90p5ve5T/` produced byte-identical new-chain `groups.json` (verified by `cmp`). The only statistical additions are 200-sweep four-quarter screens and plots/CSVs plus six-observable simultaneous chain comparisons; 13 targeted regression tests pass. No T=6.4 run or physical VBS conclusion is made here.


<a id="report-17"></a>

## kitaev_ssh/T32_INDEPENDENT_REPLICATION.md

### T=3.2 independent paired-chain replication

Fixed plan `t32_replication_plan.json` was saved and uploaded before job submission. Seeds 7201–7204 were absent from existing run metadata. Each chain: 3200 warmup, 8000 measured sweeps, stride 1, raw blocks 25. L=2, J=(1,1,.5), g=.3, mass=1.7, omega=.8, d_x=d_y=d_z=.1378676470588235, T=3.2, h=.2, M=16, Wilson (++), bond-refresh; original validated binary SHA256 `75440fcce5f9781c2f89eb776d78d4209cc46509b07cdf8def34727d1be2edf4`. No outcome-dependent extension was made.

The sole primary endpoint is each chain's phase-reweighted Szz(q1) mean over measured sweeps 6001–8000 minus its mean over sweeps 1–2000. The four paired contrasts are the independent units. The mean is their equal-weight average, SE is their sample SD divided by sqrt(4), and the 95% interval uses Student t with 3 degrees of freedom (critical value 3.1824463). The within-chain first/last relationship remains paired: Var(mean difference) = [sample Var(last)+sample Var(first)−2 sample Cov(first,last)]/4, not a sum assuming zero covariance. Each quarter still uses the original complex phase/numerator estimator; the mean contrast agrees exactly here with the pooled complex-ratio contrast. Approximate normality of chain contrasts is an assumption with only four replicates, so this is not a distribution-free guarantee.

| Seed | First quarter | Last quarter | Paired difference |
| --- | ---: | ---: | ---: |
| 7201 | .367365 | .356131 | −.011234 |
| 7202 | .364188 | .385638 | +.021450 |
| 7203 | .375634 | .381076 | +.005442 |
| 7204 | .364020 | .368702 | +.004681 |

Primary mean difference **+.005085**, SE **.006673**, chain SD .013347, 95% t(3) interval **[−.016153,+.026323]**, nominal two-sided t p=.5015. The sampled first/last chain-mean covariance is +1.37326e-5. The prespecified practical tolerance is ±.02; the complete interval does not fit inside it. Thus the old positive change is not independently supported, while small-change equivalence is also not established. In fact the interval still includes the old +.02487 point estimate. This is an uncertainty limit rather than evidence establishing either persistent drift or stationarity.

Historical seeds 7101–7104 were not pooled into the replication. Applying the same paired-chain method to them gives +.024871 ± .009422, 95% t(3) [−.005115,+.054858]. Earlier normal/quarter-block intervals should not be substituted for this small-replicate interval in this comparison.

Decision after this fixed batch: **do not yet proceed to T=6.4 under the predefined ±.02 interval criterion**. The reason is insufficient precision in the sole primary interval, not an added auxiliary 3-sigma veto. If the tolerance remains required, a separately planned independent-chain budget would be needed to resolve it; no additional chains were commissioned here.

Auxiliary 100/200/400-sweep blocking is retained in `auxiliary/`; the maximum 400/200 SE ratio among the six observables is 1.045. The existing screen does display seed 7203 Tr S(q3) half difference +.06731±.02148 and first/last-quarter difference +.08436±.02571; these remain auxiliary observations and do not alter the primary decision rule. All numerical failures were zero, quarter phases resolved, mean phase 1+O(1e-18)i; all chains have 320 raw blocks. No assumption of squared fermion weight or discarded phase was introduced.

PBS job **155141.mgt**, directory `/home/sunxr/kitaev_t32_replication_20260916.jOW7GVYu`, four single-threaded processes on four CPU slots. Job and all chains exited 0. PBS walltime **2m32s**, CPU time **9m59s**. Individual elapsed times 149.18,149.20,150.08,149.57 s; acceptance fractions .88573,.88661,.88698,.88541. Actual commands are in `t32_replication.pbs`; submission was `cd /home/sunxr/kitaev_t32_replication_20260916.jOW7GVYu && /opt/pbs/bin/qsub t32_replication.pbs`.

Raw blocks, run metadata and logs: `build-kitaev-p01/t32_replication/results/`. Original PBS stdout and exit code: `build-kitaev-p01/t32_replication/kitaev_t32_repl.o155141` and `job_exit.txt`. Main machine-readable result and plot: `paired_replication.json`, `paired_replication.png` in that directory. Repeat the primary calculation with `MPLCONFIGDIR=/tmp/kitaev-p5-mpl /home/asus/.cache/pfqmc_plot_env/bin/python kitaev_ssh/paired_t32_replication.py`; the full analysis invocation was `P5_PYTHON=/home/asus/.cache/pfqmc_plot_env/bin/python bash kitaev_ssh/analyze_t32_replication.sh` (set `REPLICATION_AUX_DIR` to a fresh directory to preserve prior auxiliary outputs). Only the Student-t critical value and paired-variance identity were directly checked for the narrow new analysis; no new testing framework or repeated historical suite was added.


<a id="report-18"></a>

## kitaev_ssh/T32_REPLICATION8.md

### T=3.2 eight-chain independent stationarity replication

The fixed plan was saved in `t32_replication8_plan.json` before submission. New seeds 7301–7308 were verified absent from existing result metadata. Physics, shifted endpoint state, T=3.2, h=.2, bond refresh and binary SHA256 `75440fcce5f9781c2f89eb776d78d4209cc46509b07cdf8def34727d1be2edf4` were unchanged. Each chain used 3200 warmup, 8000 measured sweeps, stride 1 and raw blocks of 25. Exactly eight chains were run; no result-dependent extension was made.

The sole primary statistic was fixed as Szz(q1) over measured sweeps 6001–8000 minus measured sweeps 1–2000. Each chain supplies one paired difference. The reported SE is the sample SD of the eight chain differences divided by sqrt(8); the 95% interval uses Student t with df=7, critical value 2.364624251. Blocks are not treated as independent repeats. The t interval assumes approximately normal independent chain differences, and n=8 remains a small sample.

| Seed | Q1 | Q4 | Q4−Q1 |
| --- | ---: | ---: | ---: |
| 7301 | .358981 | .369347 | +.010365 |
| 7302 | .376945 | .373942 | −.003004 |
| 7303 | .367759 | .380253 | +.012495 |
| 7304 | .387458 | .378600 | −.008857 |
| 7305 | .372021 | .370269 | −.001751 |
| 7306 | .368163 | .382696 | +.014533 |
| 7307 | .370084 | .368203 | −.001881 |
| 7308 | .377848 | .364394 | −.013454 |

Primary result: mean **+.001056**, chain SD .010288, SE **.003638**, Student-t(7) 95% interval **[−.007546,+.009657]**. The complete interval lies inside the prespecified ±.02 practical tolerance. Therefore this particular Szz(q1) early/late change satisfies the planned tolerance in the new independent batch. This does not establish stationarity of every observable or projection convergence.

Historical batches were recomputed with their appropriate t(3) intervals and were not pooled: seeds 7101–7104 give +.024871±.009422, 95% [−.005115,+.054858]; seeds 7201–7204 give +.005085±.006673, 95% [−.016153,+.026323]. The new eight-chain result is the sole decision batch.

Auxiliary joint-jackknife SEs at 100/200/400 sweeps were: Xz .001432/.001466/.001461; Btrace .002795/.002910/.003056; Szz(q1) .002205/.002272/.002116; TrS(q1) .003762/.003841/.003922; Szz(q3) .002273/.002211/.002320; TrS(q3) .004004/.003889/.003980. Every chain has 320 raw blocks, phase is 1+O(1e-18)i, and all numerical-failure counts are zero. Existing auxiliary half/quarter diagnostics are preserved in `auxiliary/`; they add no decision veto in this batch.

Runs used two simultaneous PBS jobs with four single-threaded chains each: 155142.mgt and 155143.mgt on cn1. Both exited 0. Each walltime was 2m31s; CPU times were 9m48s and 9m54s, total 19m42s. Per-chain elapsed times were 145.14–149.42 s and acceptance fractions .88337–.88746. The locked PBS command is `t32_replication8.pbs`, submitted from `/home/sunxr/kitaev_t32_replication8_20260916.4jf3XKnO` with `qsub -v BATCH=1 t32_replication8.pbs` and `qsub -v BATCH=2 t32_replication8.pbs`.

Raw blocks, metadata and logs are under `build-kitaev-p01/t32_replication8/results/`; PBS stdout files are in the parent result directory. Machine-readable primary and historical results are `paired_replication8.json`; plot `paired_replication8.png`; auxiliary blocking and diagnostics are in `auxiliary/`. Reproduce with `P5_PYTHON=/home/asus/.cache/pfqmc_plot_env/bin/python bash kitaev_ssh/analyze_t32_replication8.sh`, setting `R8_AUX_DIR` to a fresh path if preserving the existing auxiliary output.

Decision: the predefined Szz(q1) stationarity gate now permits planning the T=3.2 versus T=6.4 projection comparison. That future comparison still needs its own fixed chain budget and paired observables; no T=6.4 sampling was started here.


<a id="report-19"></a>

## kitaev_ssh/T32_T64_PROJECTION.md

### Shifted-trial T=3.2 to 6.4 bounded projection comparison

The plan in `t32_t64_projection_plan.json` was locked before reading new results. T=3.2 uses the independent seeds 7301–7308 already generated; T=6.4 uses new seeds 7401–7404, intervals=32 and dt=.2. Each new chain had 3200 warmup and 8000 measured sweeps, stride 1 and raw blocks 25. All other physics, Wilson sector, trial shift and bond-refresh algorithm were unchanged. The validated binary SHA256 was `75440fcce5f9781c2f89eb776d78d4209cc46509b07cdf8def34727d1be2edf4`. Exactly four new chains were run and no result-dependent samples were added.

The comparison uses the existing chain-preserving joint complex phase/numerator jackknife. T=3.2 and T=6.4 retain separate chain boundaries, sample counts and distribution groups; the independent group covariance contributions are summed only for the difference. Primary block length is 200 sweeps, leaving 40 blocks in every chain. Sensitivity uses 100/200/400 sweeps, with 80/40/20 blocks per chain. Marginal intervals below are normal approximate 95% intervals conditional on adequate blocking. The conservative fixed-q interval is the existing four-observable Bonferroni interval enlarged to the largest eligible SE across block levels (`blocking_look_envelope`, one planned look).

| Observable | T=3.2 | T=6.4 | T6.4−T3.2 ± SE | Marginal 95% | Conservative fixed-q/block interval |
| --- | ---: | ---: | ---: | ---: | ---: |
| Xz | .134639±.001466 | .135718±.002243 | +.001079±.002680 | [−.004173,+.006331] | block envelope [−.004342,+.006500] |
| Btrace | .223885±.002910 | .220138±.004054 | −.003747±.004990 | [−.013528,+.006034] | block envelope [−.014111,+.006617] |
| Szz q1=(0,1) | .374211±.002272 | .377846±.003385 | +.003635±.004077 | [−.004355,+.011626] | **[−.006766,+.014037]** |
| Tr S q1 | 1.120362±.003841 | 1.121488±.006043 | +.001125±.007160 | [−.012909,+.015160] | **[−.016906,+.019157]** |
| Szz q3=(1,1) | .369301±.002211 | .373507±.003494 | +.004206±.004135 | [−.003898,+.012310] | **[−.006270,+.014682]** |
| Tr S q3 | 1.124484±.003889 | 1.135397±.006479 | +.010914±.007556 | [−.003896,+.025724] | **[−.008110,+.029937]** |

All four complete conservative intervals lie inside the prospective targets: ±.02 for Szz and ±.04 for Tr S. Thus at L=2 and this parameter point, residual changes from T=3.2 to 6.4 are small enough under the predefined fixed-q criterion; this conclusion uses interval containment rather than lack of significance. The tightest margins are Szz(q3), whose conservative upper endpoint is .01468 versus .02, and Tr S(q3), .02994 versus .04. Xz and Btrace have no predeclared equivalence target here; their intervals contain zero and are reported without declaring convergence from nonsignificance.

Blocking is stable. T=6.4 SEs at 100/200/400 are: Xz .002292/.002243/.002345; Btrace .004316/.004054/.003955; Szz(q1) .003419/.003385/.003285; Tr S(q1) .005955/.006043/.005986; Szz(q3) .003464/.003494/.003204; Tr S(q3) .006236/.006479/.006377. Across T=6.4 auxiliary diagnostics, the largest targeted half-chain standardized difference is 1.42, largest chain-pair difference is 1.48, and largest first/fourth-quarter diagnostic is 2.26 (seed 7404, Tr S(q3)); these are descriptive checks and add no new veto. No clear collective time trend or chain inconsistency appears, but absence of an auxiliary excursion is not proof of mixing.

All four new chains completed with 320 raw blocks and zero numerical failures. Mean phase is 1+O(1e-18)i, CSV/JSON estimator discrepancies are below 5e-15. Acceptance fractions were .81768–.82284. Per-chain elapsed times were 288.84–289.73 s. PBS job **155144.mgt** ran four single-threaded processes on four CPU slots at cn1; exit 0, walltime **5m01s**, CPU time **19m16s**. The command is in `t64_projection.pbs`, submitted from `/home/sunxr/kitaev_t64_projection_20260916.W4Ya7ubK` using `qsub t64_projection.pbs`.

Raw T=6.4 blocks, JSON and chain logs are in `build-kitaev-p01/t32_t64_projection/results/`; PBS stdout and exit status are in its parent directory. Full group covariance, half/quarter and chain checks are in `analysis/groups.json`; comparison rows in `analysis/comparisons.csv`; plots include `differences.png`, `blocking.png`, `fixed_observables.png` and `uniform_background.png`. Reproduce with `PROJECTION_ANALYSIS_DIR=/tmp/kitaev-t32-t64-repro P5_PYTHON=/home/asus/.cache/pfqmc_plot_env/bin/python bash kitaev_ssh/analyze_t32_t64_projection.sh` using an absent output directory.

This finite-system projection check supports using T=6.4 as the current projection baseline for these fixed observables. It does not prove all observables or larger systems have converged. Before a size/coupling scan, the most informative next bounded calculation is a fixed-T=6.4 step-size comparison at h=.1 (M=64), because projection-length uncertainty is now below the chosen targets while the h=.2 discretization error remains unquantified. No T=12.8, step-size run or scan was started here.


<a id="report-20"></a>

## kitaev_ssh/T64_STEP_HALVING.md

### T=6.4 shifted-trial h=.2 to .1 bounded step comparison

The plan in `t64_step_plan.json` was locked before reading new results. The coarse group reuses independent seeds 7401–7404 at h=.2, M=32. The fine group uses new seeds 7501–7504 at h=.1, M=64. T=6.4, all physical parameters, Wilson sector, trial displacement and bond-refresh algorithm are identical; trial d was not retuned. Each fine chain used 3200 warmup and 8000 measured sweeps, stride 1 and raw blocks 25. The validated binary SHA256 was `75440fcce5f9781c2f89eb776d78d4209cc46509b07cdf8def34727d1be2edf4`. Exactly four new chains were run, with no outcome-dependent extension.

The comparison uses chain-preserving joint complex phase/numerator jackknife with separate group samples and weights. Primary blocks are 200 sweeps (40 per chain); 100 and 400 sweeps leave 80 and 20 blocks per chain. Marginal intervals are normal approximate 95% intervals conditional on adequate blocking. The conservative fixed-q interval is the existing four-observable Bonferroni interval enlarged using the maximum eligible block-level SE.

| Observable | h=.2 | h=.1 | h=.1−h=.2 ± SE | Marginal 95% | Conservative fixed-q/block interval |
| --- | ---: | ---: | ---: | ---: | ---: |
| Xz | .135718±.002243 | .135087±.002245 | −.000630±.003173 | [−.006850,+.005589] | block [−.007461,+.006200] |
| Btrace | .220138±.004054 | .224986±.004326 | +.004848±.005929 | [−.006773,+.016468] | block [−.007129,+.016825] |
| Szz q1=(0,1) | .377846±.003385 | .376409±.003918 | −.001438±.005178 | [−.011586,+.008711] | **[−.014793,+.011918]** |
| Tr S q1 | 1.121488±.006043 | 1.120737±.005914 | −.000751±.008455 | [−.017323,+.015821] | **[−.022004,+.020502]** |
| Szz q3=(1,1) | .373507±.003494 | .370084±.003087 | −.003423±.004662 | [−.012561,+.005714] | **[−.015405,+.008558]** |
| Tr S q3 | 1.135397±.006479 | 1.134520±.005823 | −.000877±.008711 | [−.017951,+.016196] | **[−.023004,+.021250]** |

All four complete conservative intervals fit within the prospective targets: ±.02 for Szz and ±.04 for Tr S. The step change therefore satisfies the predefined fixed-q criterion. This is interval containment, not a conclusion from nonsignificance. Xz and Btrace have no predefined tolerance in this comparison.

Fine-step SEs at 100/200/400 sweeps were: Xz .002152/.002245/.002578; Btrace .004314/.004326/.004165; Szz(q1) .003543/.003918/.004061; Tr S(q1) .005991/.005914/.005989; Szz(q3) .003211/.003087/.002924; Tr S(q3) .005874/.005823/.005331. The block dependence is modest and included in the conservative intervals.

Auxiliary diagnostics show no broad chain separation: largest targeted chain-pair standardized difference is 2.11 and largest half-chain difference is 1.74. One isolated time diagnostic is visible in seed 7503: Tr S(q3) Q4−Q1 = −.09175±.02708 (3.39 diagnostic SE), accompanied by Szz(q3) −.04607±.01656. These quarter diagnostics use only ten 200-sweep blocks per quarter and omit their same-chain covariance. Per the locked plan, they are retained as limitations and do not create a new single-item veto; they mean this four-chain fine group should not be described as proving universal stationarity.

Conditional second-order estimate: only if `O(h)=O(0)+c h^2+O(h^4)` applies at h=.2 and .1, then `O(.1)-O(0)=-[O(.1)-O(.2)]/3`. Linear propagation of the conservative block/multiple-observable intervals gives:

| Observable | Conditional fine−continuum estimate ± SE | Conditional conservative interval |
| --- | ---: | ---: |
| Xz | +.000210±.001058 | [−.002067,+.002487] |
| Btrace | −.001616±.001976 | [−.005608,+.002376] |
| Szz q1 | +.000479±.001726 | [−.003973,+.004931] |
| Tr S q1 | +.000250±.002818 | [−.006834,+.007335] |
| Szz q3 | +.001141±.001554 | [−.002853,+.005135] |
| Tr S q3 | +.000292±.002904 | [−.007083,+.007668] |

These are conditional extrapolation estimates, not measured continuum values. Two step sizes cannot demonstrate second-order asymptotics or exclude higher-order/nonasymptotic behavior.

All four fine chains completed with 320 raw blocks, zero numerical failures and phase 1+O(1e-18)i. Acceptance fractions were .81876–.82253. Per-chain elapsed times were 569.09–570.92 s. PBS job **155145.mgt** used four single-threaded processes on cn1; exit 0, walltime **9m34s**, CPU time **38m00s**. It was submitted as `qsub t64_h01.pbs` from `/home/sunxr/kitaev_t64_h01_20260916.rIe569eq`.

Raw fine-step blocks, JSON and logs are under `build-kitaev-p01/t64_step/results/`; PBS stdout/exit files are in the parent directory. Full group diagnostics and covariance are in `analysis/groups.json`, measured differences in `analysis/comparisons.csv`, conditional estimates in `analysis/conditional_step_residual.csv`, with difference/blocking/quarter plots in the same directory. Reproduce using an absent output path, e.g. `STEP_ANALYSIS_DIR=/tmp/kitaev-t64-step-repro P5_PYTHON=/home/asus/.cache/pfqmc_plot_env/bin/python bash kitaev_ssh/analyze_t64_step.sh`.

At L=2, the chosen fixed q observables now meet both the T=3.2→6.4 projection target and the h=.2→.1 step target. A useful next bounded computation is an L=4 timing and stability pilot at T=6.4, h=.1 with prespecified fixed observables and chain budget; convergence must be rechecked at that size. No h=.05, T=12.8 or size/coupling job was started here.


<a id="report-21"></a>

## kitaev_ssh/VALIDATION.md

### P0–P1 验证记录（2026-09-15）

本记录对应当前 patch，所有误差来自实际运行；不是预计结果。
完整构建、测试和运行命令见 [README.md](../../kitaev_ssh/README.md)。

**验收：P0–P1 在本轮限定的静态 fixed-gauge matter-sector 范围内通过。**
配置、全部 target 重编译均 exit 0。共享符号修正后的最终 CTest 为
8/8 通过，总耗时 20.68 秒：新 kitaev_ssh_p01，以及原 density_sign、
cdw_initial、channel_negative、channel_positive、channel_eta、
channel_unknown、channel_anomalies。另有旧入口两种初态 smoke 成功。
这不宣称 P2 的有限频率声子 QMC 已实现。

#### 环境与可复现性

工作目录：`pure-state-density-release-20260915`。独立构建目录：
`build-kitaev-p01`；GCC 15.2.0、Intel icx/ifx 2026.1、MKL、
Eigen 3.4.0、Boost 1.70、CMake 4.4.3。使用 Release、C++17、原仓库
`-fno-fast-math -ffp-contract=off`，串行编译、MKL 顺序库和单线程运行。

没有缺失数值依赖。CMake 不在默认 PATH，已使用 README 中记录的现存
工具绝对路径。构建中出现原 Boost 的 `std::unary_function` deprecated
以及 Eigen/MKL 的 enum comparison 警告，不影响编译和测试。

#### 独立对照及容差依据

单键直接使用 2×2 Pauli/Fock 矩阵，验证多体能级、算符共轭演化、Green
符号、parity 及 0/1/2/3 个传播因子的振幅。2×2 蜂窝用占据比特/Jordan–Wigner
直接构造 8 个 Majorana 算符和 16×16 多体 H，再指数演化。不从被测
K_code/Green/权重函数反构造参考。两边仅共享基础线性代数库。

额外试探 Hamiltonian 是目标 H 加 `0.17*i*c_0*c_1 + 0.11*i*c_4*c_6`；
其与 H 的对易子范数为 4.1439691118539956，目标基态重叠平方
0.99545068781388346。对两个试探态，取 theta=0.05、0.2、0.8、2，
dt=0.025，block=1、3、8，比较归一化能量、Green、Z_theta。
短轮廓还用既有全局 Pfaffian 例程独立检查范数和相位。

L=4、6 **不做完整多体 ED**。参考是从三条键位移独立推导的 2×2 Bloch
矩阵、逐离散动量谱和 Fourier Green。投影长度 theta=0.1、0.8、2，
dt=0.05、block=4。随机数固定种子 20260915+L，每个尺寸测试三组局域 Z2
变换，包含正反单键、多步和绕行开放路径。g=0 的 X 独立性要求逐元素零误差。

双精度机器 epsilon 约 2.22e-16。小矩阵代数比较容差 1e-14–2e-14；
最大 72 维谱和 Fourier 求和用 2e-12（维数乘 epsilon 再留求解/求和裕量）。
重复至 160 个 factor 的投影误差阈值 2e-10；粗略 `factors*NM*epsilon`
累计尺度在 1e-12 级，阈值预留 QR/SVD、范数对数累计及数学库差异，
不以最坏指数传播的未稳定化条件数为借口放宽阈值。实际残差远低于阈值。
结构和计数的整数检查容差为 0，g=0 矩阵/可观测量独立性亦为 0。
下面所有阈值均在首次测试前设定，单键符号失败通过修代码解决，未调宽容差。

Dirac 速度用 h=1e-3、5e-4、1e-4 的对称差商；解析截断项分别为
`h^2/8` 和 `3*h^2/32` 量级，所以最大相对容差 3e-7 是 O(h²) 截断界，
不是机器精度比较；节点能量容差 2e-15。

#### 关键最大误差

| 比较 | 最大误差 | 容差 |
| --- | ---: | ---: |
| 单键能级 | 0 | 1e-14 |
| 单键算符演化 / 因子 2 / Green | 2.23e-16 | 2e-14 |
| 单键两种 Pfaffian 权重 | 4.45e-16 | 2e-14 |
| 2×2 ED 投影能量 | 1.78e-15 | 2e-10 |
| 2×2 ED 投影 Green | 8.00e-15 | 2e-10 |
| 2×2 ED log(Z) / Z 相对误差 | 1.43e-14 | 2e-10 |
| 2×2 ED Pfaffian 权重相对误差 | 7.78e-16 | 2e-10 |
| 2×2 ED Pfaffian 相位 | 2.60e-17 | 2e-12 |
| L=4/6 Bloch 谱 / 正负配对 | 1.29e-14 / 9.33e-15 | 2e-12 |
| L=4/6 Bloch 基态 Green | 1.34e-15 | 2e-12 |
| L=4/6 基态投影能量 / Green | 7.82e-14 / 7.11e-15 | 2e-10 |
| L=4/6 基态投影 log(Z) | 2.28e-13 | 2e-10 |
| gauge 变换 K / 能量 / Green | 0 / 2.14e-14 / 5.89e-15 | 0 / 2e-10 / 2e-10 |
| gauge 变换 log(Z) / string | 1.64e-13 / 3.06e-15 | 2e-10 |
| g=0 更换 X 的 K / 能量 / Green / log(Z) | 0 | 0 |
| Green 结构 / overlap 求解残差 | 3.48e-15 / 1.27e-15 | 1e-10 |
| 静态非零 g、显式 Wilson(-1,+1) 的 ED Green | 1.34e-15 | 2e-12 |
| 连续 Dirac 节点能量 | 0 | 2e-15 |
| Dirac x/y 速度相对误差 | 1.251e-7 / 9.375e-8 | 3e-7 |

表内非零误差向上取简短位数；完整数值及每项比较次数见下方原始输出。
几何测试覆盖 2×2、3×4、4×4、6×6，包括全部四个环面扇区、六边形
闭合/唯一顶点/每键属于两个 plaquette、每 site 各一条 x/y/z 键、PBC 位移。

#### 自由入口实测

参数 theta=0.4、dt=0.05、block=4，J=(1,1,0.6)，g=0，Wilson(+1,+1)：

| L | site=NM | bond | E0（Bloch） | projector E | log(Z) | matter parity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 32 | 48 | -22.431439354917106 | -22.431439354917110 | 17.945151483933635 | -1 |
| 6 | 72 | 108 | -50.372700234707914 | -50.372700234707978 | 40.298160187766321 | -1 |

两有限网格的最小正 K 本征值分别为 1.2、0.8，没有精确零模。
连续节点 `(k1,k2)=+/- (1.8754889808102941,-1.8754889808102941)`，
Cartesian `kx=+/-2.1656281358660006, ky=0`，
速度 `vx=3.3045423283716611, vy=1.8`。

显式 `--Lx 6 --Ly 6 --Jz 1` 返回 exit 1 和 `zero_mode_ambiguous`，
无自动质量/边界/gauge 改动。显式 `--lambda 0.1` 返回 exit 1 的未知参数
错误，确认未把未实现功能做成可运行假入口。

#### 发现并修复的问题

第一次单键物理权重比较失败：旧代码给 `-exp(0.19*0.73)`，应为正，
绝对误差 2.2975588289580355（阈值 2e-14）。定位到两处全局权重
缺少原始块 Pfaffian 的 `(-1)^(NM/2)` 排序常数。仅补这个相位因子，
不改旧模型或低层 Pfaffian；修正后奇数/偶数因子单键振幅均通过。
详见 CONVENTIONS 第 4 节和实际 patch。

对原始 SHA256SUMS 的再检查恰有三个预期差异：CMakeLists.txt、
inc/pure_projector_weight.h、inc/pure_projector_fast.h；其余 121 个原始
清单文件校验通过。原始模型、density 入口、脚本和回归源码未修改。

#### 旧模型 smoke 实测

使用 README 中命令，L=18、Vf=4、PBC、tau=0.1、dt=0.01、block=16、
burn=1、measurements=1、seed=18150915，两种初态均正常完成并生成 CSV：

| 初态 | 接受率 | 测量能量 | 最大 endpoint Green 残差 | 配置哈希 |
| --- | ---: | ---: | ---: | --- |
| free | 0.9125 | -18.301717082150436 | 6.56e-15 | 14311318817813035757 |
| cdw-even | 0.77916666666666667 | -20.733065282564631 | 5.95e-15 | 1078904481023143977 |

两次 sign mismatch、MP fallback failure 均为 0。CSV 暂存在
`/tmp/kitaev-ssh-density-smoke.nTPiaA/`；不依赖这些临时文件重跑验证。
没有把 smoke 的单次测量当作统计收敛结果。

#### 边界与进入 P2 的条件

本阶段只验收静态固定 gauge matter sector。未实现物理自旋投影、零模
占据策略、连续声子路径、ramp 或 VBS 生产测量。未进行大系统长投影稳定性
或一般复杂轮廓分支完备性验收。

P2 必须先推导开放虚时间的离散声子作用量和端点试探权重，确定共享
X[M] 与时间节点/factor 映射，设计连续 proposal 的多因子原子提交及
拒绝恢复，然后为整个路径建立单份 Majorana 权重的独立小系统验证。
不能把当前静态正范数或旧 HS1 正性扩展成声子正权重假设。

#### 新验证器完整标准输出

```text
DETAIL ed_noncommuting_ground_overlap=0.99545068781388346 commutator_norm=4.1439691118539956
DETAIL L=4 sites=32 bonds=48 matter_parity=-1 E0=-22.431439354917106 min_positive_energy=1.1999999999999984
DETAIL L=6 sites=72 bonds=108 matter_parity=-1 E0=-50.372700234707914 min_positive_energy=0.80000000000000038
DETAIL dirac_k1=1.8754889808102941 dirac_k2=-1.8754889808102941 cartesian_kx=2.1656281358660006 cartesian_ky=0 vx=3.3045423283716611 vy=1.8
PASS bloch_finite_grid_nonzero max_error=0 tolerance=0 comparisons=52
PASS bloch_ground_green max_error=1.3317876976328526e-15 tolerance=2e-12 comparisons=2
PASS bloch_spectrum max_error=1.2878587085651816e-14 tolerance=2e-12 comparisons=104
PASS bond_length max_error=0 tolerance=1.0000000000000001e-15 comparisons=204
PASS dirac_node_energy max_error=0 tolerance=2.0000000000000002e-15 comparisons=2
PASS dirac_x_slope_relative max_error=1.2500011381533227e-07 tolerance=2.9999999999999999e-07 comparisons=6
PASS dirac_y_slope_relative max_error=9.3749934149123249e-08 tolerance=2.9999999999999999e-07 comparisons=6
PASS ed_noncommuting_trial max_error=0 tolerance=0 comparisons=1
PASS ed_pfaffian_phase max_error=2.5991574317133729e-17 tolerance=2e-12 comparisons=4
PASS ed_pfaffian_status max_error=0 tolerance=0 comparisons=4
PASS ed_pfaffian_weight_relative max_error=7.7759063297812442e-16 tolerance=2.0000000000000001e-10 comparisons=4
PASS ed_projection_energy_descent max_error=0 tolerance=0 comparisons=8
PASS ed_projector_energy max_error=1.7763568394002505e-15 tolerance=2.0000000000000001e-10 comparisons=24
PASS ed_projector_green max_error=7.9936059655306692e-15 tolerance=2.0000000000000001e-10 comparisons=24
PASS ed_projector_log_norm max_error=1.4210854715202004e-14 tolerance=2.0000000000000001e-10 comparisons=24
PASS ed_projector_norm_relative max_error=1.4210854715202105e-14 tolerance=2.0000000000000001e-10 comparisons=24
PASS ed_trial_green max_error=1.9984014443252818e-15 tolerance=9.9999999999999998e-13 comparisons=2
PASS ed_trial_overlap max_error=0 tolerance=0 comparisons=1
PASS ed_trial_parity max_error=2.2204460492503131e-16 tolerance=9.9999999999999998e-13 comparisons=2
PASS exact_zero_mode_diagnostic max_error=0 tolerance=0 comparisons=1
PASS fock_majorana_algebra max_error=0 tolerance=0 comparisons=64
PASS g0_X_matrix_independence max_error=0 tolerance=0 comparisons=2
PASS g0_X_observable_independence max_error=0 tolerance=0 comparisons=6
PASS gauge_energy max_error=2.1316282072803006e-14 tolerance=2.0000000000000001e-10 comparisons=6
PASS gauge_generator max_error=0 tolerance=0 comparisons=6
PASS gauge_green max_error=5.8882402522229788e-15 tolerance=2.0000000000000001e-10 comparisons=6
PASS gauge_log_norm max_error=1.6342482922482304e-13 tolerance=2.0000000000000001e-10 comparisons=6
PASS gauge_matter_parity max_error=0 tolerance=0 comparisons=6
PASS gauge_string max_error=3.0501155054132858e-15 tolerance=2.0000000000000001e-10 comparisons=1098
PASS gauge_wilson_invariance max_error=0 tolerance=0 comparisons=12
PASS generator_hermitian max_error=0 tolerance=0 comparisons=2
PASS generator_skew max_error=0 tolerance=0 comparisons=2
PASS geometry_boundary_counts max_error=0 tolerance=0 comparisons=4
PASS geometry_coordinates max_error=6.6613381477509392e-16 tolerance=1e-14 comparisons=408
PASS geometry_counts max_error=0 tolerance=0 comparisons=4
PASS geometry_degree max_error=0 tolerance=0 comparisons=408
PASS geometry_hexagon_vertices max_error=0 tolerance=0 comparisons=68
PASS geometry_indices max_error=0 tolerance=0 comparisons=408
PASS geometry_plaquette_closure max_error=0 tolerance=0 comparisons=68
PASS geometry_plaquette_count max_error=0 tolerance=0 comparisons=4
PASS geometry_plaquette_displacement max_error=0 tolerance=1e-14 comparisons=68
PASS geometry_plaquette_incidence max_error=0 tolerance=0 comparisons=204
PASS geometry_type_counts max_error=0 tolerance=0 comparisons=12
PASS geometry_unique_bonds max_error=0 tolerance=0 comparisons=4
PASS ground_projector_energy max_error=7.815970093361102e-14 tolerance=2.0000000000000001e-10 comparisons=6
PASS ground_projector_green max_error=7.1089496906354422e-15 tolerance=2.0000000000000001e-10 comparisons=6
PASS ground_projector_log_norm max_error=2.2737367544323206e-13 tolerance=2.0000000000000001e-10 comparisons=6
PASS invalid_input_rejection max_error=0 tolerance=0 comparisons=13
PASS model_eigensolver max_error=0 tolerance=0 comparisons=2
PASS projector_green_residual max_error=3.4733196087061919e-15 tolerance=1e-10 comparisons=6
PASS projector_solve_residual max_error=1.2613743984373174e-15 tolerance=1e-10 comparisons=6
PASS real_A_skew max_error=0 tolerance=0 comparisons=2
PASS reciprocal_vectors max_error=0 tolerance=2.9999999999999998e-15 comparisons=16
PASS reference_eigensolver max_error=0 tolerance=0 comparisons=17
PASS single_bond_factor_two max_error=2.2204460492503131e-16 tolerance=2e-14 comparisons=1
PASS single_bond_ground_green max_error=2.2204460492503131e-16 tolerance=2e-14 comparisons=1
PASS single_bond_levels max_error=0 tolerance=1e-14 comparisons=1
PASS single_bond_log_weight max_error=2.2204460492503131e-16 tolerance=2e-14 comparisons=4
PASS single_bond_matter_parity max_error=0 tolerance=0 comparisons=1
PASS single_bond_operator_evolution max_error=2.2204460492503131e-16 tolerance=2e-14 comparisons=2
PASS single_bond_pfaffian_status max_error=0 tolerance=0 comparisons=5
PASS single_bond_raw_pfaffian_anchor max_error=0 tolerance=0 comparisons=1
PASS single_bond_stable_weight max_error=4.4408920985006262e-16 tolerance=2e-14 comparisons=4
PASS single_bond_weight max_error=4.4408920985006262e-16 tolerance=2e-14 comparisons=4
PASS spectrum_pairing max_error=9.3258734068513149e-15 tolerance=2e-12 comparisons=104
PASS static_X_generator max_error=0 tolerance=0 comparisons=1
PASS static_X_trial_green max_error=1.3353649518931756e-15 tolerance=2e-12 comparisons=1
PASS trial_dimensions max_error=0 tolerance=0 comparisons=2
PASS wilson_sectors max_error=0 tolerance=0 comparisons=32
PASS winding_closure max_error=0 tolerance=0 comparisons=32
PASS winding_displacement max_error=8.8817841970012523e-16 tolerance=2e-14 comparisons=32
PASS zero_flux_sectors max_error=0 tolerance=0 comparisons=272
PASS P0-P1 fixed-gauge matter-sector baseline; no phonon sampling or spin projection
```


<a id="report-22"></a>

## kitaev_ssh/WEIGHT_OPTIMIZATION.md

### L=4 weight-engine benchmark

Parameters were Lx=Ly=4, tau=6.4, dt=0.1 (64 intervals), J=(1,1,.5),
g=.3, mass=1.7, omega_ph=.8, and the inherited trial shift
0.1378676470588235 on every bond. Each run used one math-library thread and
30 sweeps with the same seed 7701 and proposal trajectory. Initialization and
JSON serialization were outside the timed region.

The 20022 HPC array used two repetitions per engine. Mean wall time per sweep:

| engine | seconds/sweep | speedup vs reference |
|---|---:|---:|
| reference | 5.78047929 | 1.000 |
| svd-values | 5.66364755 | 1.0206 |
| node-merge | 3.57863558 | 1.6153 |
| combined | 3.38812236 | 1.7061 |

The combined profile measured 3.39681232 s/sweep. Its timed fractions were
physical Pfaffian 70.20%, factor generation/storage 20.59%, center Green
5.24%, merge 2.59%, and factor validation 1.10%; explicit path/factor copies
were each below 0.02%.

Across the fixed trajectory, all engines accepted 1174 proposals. The largest
cross-engine differences were 0 in final log absolute weight,
5.15e-16 in final phase, and 4.81e-16 in smoke block observables. Independent
replay at L=4 gave maximum log-absolute error 1.42e-13, phase error
1.63e-15, Green error 0, log-ratio error 1.71e-13, and acceptance-probability
error 1.14e-13 (48 accepted and 24 rejected transactions).

The production default for the new SSH--Kitaev CLI is now `combined`; the
`reference`, `svd-values`, and `node-merge` engines remain selectable with
`--weight-engine`. Green reconstruction always retains the original half-step
ordering. A two-sweep default smoke on the rebuilt Intel/MKL binary recorded
`weight_engine=combined` in its JSON metadata.

For eager local proposals, validated prefix checkpoints and Cayley LU records
are now enabled by default (`--prefix-cache 1 --cayley-cache 1`). They are
inactive in the default deferred phase-on-demand path and can be disabled
individually for an ablation comparison. An L=4 same-seed short run reduced
wall time from 73.73 s to 30.33 s with identical blocks and zero numerical
failures.

Reproducibility:

* raw benchmark and smoke files: `build-kitaev-p01/hpc_20022_results/results/`;
* summary: `build-kitaev-p01/hpc_20022_results/results/benchmark_summary.json`;
* plan and scripts: `kitaev_ssh/weight_optimization_plan.json`,
  `weight_build.pbs`, `weight_array.pbs`, and `weight_replay.pbs`;
* HPC jobs: build 84616 (binaries), array 84619, replay 84621, default rebuild
  and smoke 84622, on SSH port 20022.

The local CMake rebuild was attempted but cannot find `mkl_types.h`; the
production compile and default smoke were completed in the HPC Intel/MKL
environment. No long sampling job was submitted.


<a id="report-23"></a>

## kitaev_ssh/kitaev-ssh-p5-convergence-and-prompt.md

### P5：采样与投影收敛（2026-09-15）

本轮得到可复现的有限 `(L,T,h)` 基准，**尚未得到整体投影收敛或可信 VBS 峰**。原始任务文件 `/home/asus/susy/kitaev-ssh-p5-convergence-and-prompt.md` 保留不变。本记录只汇总运行与结论，不重复物理推导。

#### 选择、代码与统计约定

采用已有 `bond-refresh`：按 P2 的相关 Gaussian AR 分布独立刷新一整条键的开放路径，`q_b=Wph_b*exp(omega*T/2)`；正反提议比抵消声子权重，接受率为 `min(1,abs(Wf_new)/abs(Wf_old))`。这不是独立刷新各节点，也不把两侧路径镜像复制。复用既有物理 Pfaffian、半步映射、事务与 Fock 检查，未更改 C++ 数值核心、旧模型或任何观测量归一化。

新增 `convergence_analysis.py`、`test_convergence_analysis.py`、`run_p5.sh`、`analyze_p5.sh`、`p5_plan.json`。分析从已有 `sum(s)`、`sum(sO)` 和计数开始；不同 T/h 不合并，重分块不跨链，不丢尾块。独立链按样本数保持原权重，使用分层联合 delete-one-block jackknife，相加各链协方差；不是平均各链比值。复数及实虚协方差、背景扣除均保留。不同分布或重复种子输入报错。

`S=<F F†>/Nc`、`m_z²=Szz/Nc`，`q=n_x*b1/Lx+n_y*b2/Ly`。本轮只固定 L=2 的 `q1=(0,1)`、`q3=(1,1)`，不重新搜索峰。背景 `Btrace=Nc*sum_mu <X_mu>²` 联合重采样；它不是调制序参量。

新增运行前选择绝对变化尺度：Szz 为 0.02，Tr S 为 0.04。主块长 100 sweeps；诊断 25/50/100/200/400，始终保留链边界。表内 `d=B-A`、`SE_d=sqrt(SE_A²+SE_B²)`、`z=d/SE_d`、近似边际 95% 区间 `d±1.96*SE_d`，依赖足够分块的假设。CSV 另有四观测量 Bonferroni 区间、两次查看数据的修正及最大分块 SE 包络；后者是保守敏感性检查，不是长自相关尾或序贯覆盖率的严格保证。没有从块文件臆造独立样本数。

#### P4 实际数据审计

以下历史比较均 L=2、g=0.3、lambda=0.08272058823529412、mass=1.7、omega=0.8、J=(1,1,0.5)。`q2=(1,0)`。表中 |d| 是绝对差；符号见 d。

| P4 比较 B−A | 观测量/动量 | d | abs(d) | SE_d | z | 边际 95% 区间 |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| local half_step−seed_a；T=.8，h .2→.1 | Szz(q2) | −.148883 | .148883 | .048057 | −3.098 | [−.243075,−.054691] |
| refresh_half−refresh_a；T=.8，h .2→.1 | Szz(q3) | −.028246 | .028246 | .020174 | −1.400 | [−.067787,.011295] |
| refresh_long−refresh_a；h=.2，T .8→1.6 | Tr S(q1) | +.070998 | .070998 | .031476 | +2.256 | [.009305,.132692] |

这些是 P4 多动量/通道比较后挑出的差异，不是三个事先规定的检验；2.3 sigma 对应非零 q1 的 Tr S，不是 q=0。历史重现使用每链 16 块：local 200 sweeps/块，refresh 100 sweeps/块。

完整实际元数据（含 L、g/lambda、m、omega、T、h、种子、热化、测量、块长、复相位及误差、耗时）见 [p4_metadata.csv](../../build-kitaev-p01/p5_runs/p4_audit/p4_metadata.csv)，历史差值见 [p4_discrepancies.csv](../../build-kitaev-p01/p5_runs/p4_audit/p4_discrepancies.csv)。local seed_a/b 为 4201/4202，warmup=400、sweeps=3200、stride=2，分别 81.22/81.68 s；refresh_a/b 为 5201/5202，warmup=100、sweeps=1600、stride=1，分别 9.73/9.92 s。两组都是 M=4、h=.2、T=.8，故算法比较没有夹带步长变化。local half_step/long_projection 热化各 800；全部详细条目保存在上述 CSV。

在同目标的两组独立种子上，`SE²*sum(chain_seconds)` 的 local/refresh 比依次为 Szz(q1) 24.2、Tr S(q1) 66.8、Szz(q3) 41.7、Tr S(q3) 36.6。旧 local 分块误差尚不稳定，这只是效率诊断，不能据此证明 local 有偏或百分比精度等价。完整绝对差和区间在 [algorithm_benchmark/comparisons.csv](../../build-kitaev-p01/p5_runs/algorithm_benchmark/comparisons.csv)。结合既有单键/Fock/提议密度检查，后续选 refresh；没有在昂贵点重复低效 local。

#### 新运行与比较

所有新链：Lx=Ly=2、J=(1,1,.5)、g=.3、mass=1.7、omega=.8、k=1.088、lambda=.08272058823529412；空间 PBC、零通量 Wilson++、固定 gauge matter-sector；开放轮廓总长 `T=M*h=2*Theta`。每链 warmup=400、测量 8000 sweeps、stride=1、原始块=25 samples，主分析每链 80 块。每个 refresh sweep 更新每条物理键路径一次。

| T | h | M | 独立种子 | 总测量数 | 各链耗时之和/s | 接受率范围 |
| ---: | ---: | ---: | --- | ---: | ---: | --- |
| .8 | .2 | 4 | 6101,6102 | 16000 | 92.87 | .9630–.9636 |
| 1.6 | .2 | 8 | 6201,6202,6203 | 24000 | 262.61 | .9321–.9347 |
| 3.2 | .2 | 16 | 6301,6302,6303 | 24000 | 509.66 | .8840–.8861 |
| 3.2 | .1 | 32 | 6401,6402 | 16000 | 681.30 | .8854–.8866 |

总计 10 链、80000 次测量，链耗时之和 **1546.44 s（25.77 分钟）**，不是并行作业总墙钟时间。全部完成，数值失败 0。平均相位依次为 `1+i*(4.083e-18,6.964e-18,1.349e-17,1.606e-17)`，虚部 SE 为 `(5.48e-19,6.76e-19,9.03e-19,1.60e-18)`，实部 SE 数值为零；未把相位强制设为 1。[metadata.csv](../../build-kitaev-p01/p5_runs/final_analysis/metadata.csv) 保存每链数据。

先增加 T=.8/1.6 的独立统计，发现原 .071 差异未重现而均匀背景明显变化，因此仅延长到 T=3.2。之后按固定预算各补一条 T=1.6/3.2 链，并仅在 T=3.2 做 h/2；没有为取得“通过”持续加样本。过程决定保存在 `p5_plan.json`。

| 比较 B−A | 观测量 | d | abs(d) | SE_d | z | 边际 95% 区间 |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| T1.6−T.8，h=.2 | Szz(q1) | .001985 | .001985 | .005565 | .36 | [−.00892,.01289] |
| 同上 | Tr S(q1) | .008427 | .008427 | .009958 | .85 | [−.01109,.02794] |
| 同上 | Szz(q3) | −.005002 | .005002 | .005451 | −.92 | [−.01569,.00568] |
| 同上 | Tr S(q3) | .001558 | .001558 | .009402 | .17 | [−.01687,.01999] |
| T3.2−T1.6，h=.2 | Szz(q1) | .000627 | .000627 | .005077 | .12 | [−.00932,.01058] |
| 同上 | Tr S(q1) | .001962 | .001962 | .009235 | .21 | [−.01614,.02006] |
| 同上 | Szz(q3) | .010159 | .010159 | .004942 | 2.06 | [.00047,.01984] |
| 同上 | Tr S(q3) | .012530 | .012530 | .008909 | 1.41 | [−.00493,.02999] |
| h.1−h.2，T=3.2 | Szz(q1) | −.008206 | .008206 | .005636 | −1.46 | [−.01925,.00284] |
| 同上 | Tr S(q1) | −.015130 | .015130 | .009882 | −1.53 | [−.03450,.00424] |
| 同上 | Szz(q3) | −.002532 | .002532 | .005549 | −.46 | [−.01341,.00834] |
| 同上 | Tr S(q3) | −.011337 | .011337 | .009891 | −1.15 | [−.03072,.00805] |

更多背景比较、原始均值与保守区间见 [comparisons.csv](../../build-kitaev-p01/p5_runs/final_analysis/comparisons.csv)。[groups.json](../../build-kitaev-p01/p5_runs/final_analysis/groups.json) 含完整复数结果、联合协方差、各分块层、独立链及半链检查；原始 `*_blocks.csv`、`*_results.json`、stdout/stderr 位于 `build-kitaev-p01/p5_runs/`。

#### 支持的结论与限制

- **统计**：25→200 sweeps 分块未显示普遍持续增长，400-sweep 层每链仅 20 块，仍有散布；例如 T=.8 的 Szz(q3) SE 从 100-sweep 的 .00423 到 400-sweep 的 .00512。四个固定观测量独立链比较最大 |z|=1.53，只是辅助诊断，不是混合证明。保守包络保留这些分块敏感性，不用最小误差层挑结果。
- **投影**：T1.6→3.2 的 q1 Szz/Tr S 包络分别 [−.01358,.01483]、[−.02377,.02769]，限制了这次投影加倍的变化尺度；不能据此限制所有更长 T。q3 Szz 包络 [−.00358,.02389] 尚未满足 ±.02。均匀背景 Btrace 增加 .062764，SE=.003751，边际区间 [.055411,.070116]；Xz 增加 .034334，SE=.003168。**整体尚未投影收敛**。
- **步长**：h.1−h.2 的 q1 包络仍为 Szz [−.02406,.00765]、Tr S [−.04376,.01350]，不能仅凭不显著宣布步长收敛。若已在二阶渐近区间，细步长残差 `(O_h−O_h/2)/3` 的保守范围：q1 Szz [−.00255,.00802]、Tr S [−.00450,.01459]；q3 Szz [−.00437,.00606]、Tr S [−.00541,.01296]。条件估计小于预定尺度，但两步长不能验证该假设，未做第三步长。
- **峰**：最细步长的固定 q1 Tr S 相对其余非零 q 平均差为 −.00539±.00917（SE）；q3 Szz 为 .00301±.00559。没有可信调制峰，更没有非公度 VBS 结论。
- **下一项**：先小规模测试 T=6.4、h=.2 的成本与稳定性，再在同一物理点比较固定 q 和均匀背景；优先解决投影，不追加尺寸/耦合扫描或盲目缩小步长。一般长轮廓数值恢复和大尺寸混合尚未验证。

实际采样图：[固定观测量](../../build-kitaev-p01/p5_runs/final_analysis/fixed_observables.png)、[分块稳定性](../../build-kitaev-p01/p5_runs/final_analysis/blocking.png)、[绝对差及预定尺度](../../build-kitaev-p01/p5_runs/final_analysis/differences.png)、[均匀背景](../../build-kitaev-p01/p5_runs/final_analysis/uniform_background.png)。没有生成合成物理构型。

#### 可复现命令与检查

在 release 目录执行；运行脚本含完整物理参数和逐链种子。现有结果拒绝覆盖；重跑时使用新的 `P5_OUTPUT_DIR`，同种子重跑结果不能当作新增独立链。

```bash
source /opt/intel/oneapi/setvars.sh
/home/asus/code/driven-new-validation-20260908/green-residual/tools/cmake/data/bin/cmake --build build-kitaev-p01 --target kitaev_ssh_mc -j 1
for stage in short long projection_a projection_b supplement fine_a fine_b; do
  bash kitaev_ssh/run_p5.sh "$stage"
done
MPLCONFIGDIR=/tmp/kitaev-p5-mpl /home/asus/.cache/pfqmc_plot_env/bin/python kitaev_ssh/test_convergence_analysis.py
P5_PYTHON=/home/asus/.cache/pfqmc_plot_env/bin/python bash kitaev_ssh/analyze_p5.sh
```

已有数据无需重跑 MC：设 `P5_ANALYSIS_DIR=build-kitaev-p01/p5_reanalysis` 执行最后一条即可重新生成图表。P4 审计可执行 `convergence_analysis.py --audit-p4 build-kitaev-p01/p4_runs --output NEW_EMPTY_DIRECTORY`。本机系统 Python 没有 NumPy，使用上述已有 NumPy/Matplotlib 环境；无需新依赖。

实际构建检查通过；10 个新增分析测试通过（0.102 s），包括 P4 单链联合 jackknife 数值回归、复相位非线性量、不同链长权重、链边界、重复种子/错误目标、未解析相位以及二阶残差符号。运行脚本语法及重复输出保护检查通过。见 `p5_runs/build.log`、`analysis_tests.log`、`numerical_code_sha256.txt`。未修改公共数值层或构建配置，因此没有重跑无关旧 density/P0–P4 全套；复用了 P4 已完成的相关单键/Fock/事务验证。


<a id="report-24"></a>

## kitaev_ssh/kitaev-ssh-p5-followup-projection-and-prompt.md

### P5 有界长投影检查：2026-09-16

**结论：T=6.4 运行稳定，均匀背景仍有投影漂移；漂移与未位移声子初态松弛的一阶诊断相容。短批不足以确认固定动量结构因子的 .02/.04 收敛，不再追加样本。** 原始任务文件保留不变。

#### 预算与实际执行

沿用 Lx=Ly=2、J=(1,1,.5)、g=.3、mass=1.7、omega=.8、k=1.088、lambda=.0827205882、Wilson++、原自由 matter 试探态和未位移声子真空。开放轮廓 T=6.4、h=.2、M=32；bond-refresh、物理 Pfaffian 和观测量归一化均不变。

在查看新结果前固定两条链，各热化 400、测量 1600 sweeps、stride=1、原始块=25。目的为检查稳定性、成本和约 .028 的 Xz 漂移，不是让全部结构因子通过精度要求。计划见 `p5_followup_plan.json`。

经 SSH **20015** 提交 PBS **155137.mgt**，计算节点 cn1，申请 **2 CPU**；各链一个 CPU、MKL/OMP=1，同时运行。

| 种子 | 测量数 | 链程序耗时/s | 接受率 | 数值失败 |
| ---: | ---: | ---: | ---: | ---: |
| 6601 | 1600 | 49.5358 | .8209375 | 0 |
| 6602 | 1600 | 49.1553 | .8221875 | 0 |

PBS 墙钟 **84 s**，包含自由诊断编译；CPU 时间 132 s；链耗时之和 98.6912 s。两个退出码均 0。相位分别为 `1−5.77e−19 i`、`1−4.33e−18 i`，虚部 SE 为 `3.14e−18`、`2.48e−18`，保留复数而未强制正权重。这里的稳定仅指本次路径、长度和检查没有报错，不外推到更长轮廓。

使用已通过 HPC smoke 的二进制 `/home/sunxr/kitaev_ssh_project_smoke_v2.AEKaLn/mc_smoke`，SHA256 `3c79a5650e39549dc90a989ee3a53aa9b45123c7c74d8364958bb977ad2c38f1`。P5 旧链没有重跑。跨编译平台不保证随机轨迹逐位一致；原数值约定及依赖验证沿用 HPC smoke，未把平台差异解释为物理效应。

#### 实际自由基准与响应诊断

`free_projector.cpp` 仅增加诊断输出，用已验证的符号 `B_b=-i*u_b*G_code(A,B)`，按每方向 Nc 条键平均。得到三个方向 B_mu=-.5，最大误差约 4.5e−16；E0=-5.000000000000001，matter parity=-1，投影前后 Green 最大差 1.04e−15。独立 Bloch 基准四个振幅为 2.5、.5、.5、−1.5，能量 -5、各方向每键导数 -.5，与实际态相符。

MC 使用 `freeGroundState(model)`，模型构造的静态 X 默认全零，因此即使 g=.3，试探生成矩阵仍等于上述 g=0 基准。没有改变试探态、占据或 gauge；这是固定 gauge matter-sector，不是物理自旋投影。

一阶响应只作未拟合诊断，未校正采样结果：

`Xz(T,h) ≈ g/(2k) * [1-exp(-omega*T/2)] * (omega*h/2)/tanh(omega*h/2)`。

使用实际自由 B_mu 计算离散 Cw 与闭式之差最大 2.78e−17。T3.2→6.4 预测 Xz 增量 **.02773354**，实际 **.02731147±.00851758**（SE）；残差 −.00042207，截断的高阶物理修正不包含在该统计 SE 中。Btrace 同阶诊断增量 .07562225，实际 .08156447±.01066552。方向和量级支持初态松弛解释，不证明耦合模型严格等于线性源模型。

#### 与旧 T=3.2 数据比较

基线仅使用 P5 的 t32_a/b/c，保留各链边界。总和之比及联合分块 jackknife 完全复用 `convergence_analysis.py`，没有合并不同 T 的分布。主块长 100 sweeps，新链各 16 块；差值 `d=O(6.4)-O(3.2)`，`SE_d=sqrt(SE_old²+SE_new²)`。下表是近似正态边际 95% 区间 `d±1.96 SE_d`，以足够分块为前提，并非精确覆盖保证。

| 观测量 | T=3.2 | T=6.4 | d | SE_d | z | 边际95%区间 |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Xz | .095406 | .122717 | +.027311 | .008518 | 3.21 | [.010617,.044006] |
| Btrace | .114060 | .195624 | +.081564 | .010666 | 7.65 | [.060660,.102469] |
| Szz(q1) | .373008 | .354900 | −.018108 | .008205 | −2.21 | [−.034190,−.002025] |
| Tr S(q1) | 1.120154 | 1.126253 | +.006098 | .016364 | .37 | [−.025976,.038173] |
| Szz(q3) | .372344 | .368203 | −.004141 | .009711 | −.43 | [−.023175,.014892] |
| Tr S(q3) | 1.131197 | 1.139777 | +.008579 | .018050 | .48 | [−.026799,.043958] |

q1=(0,1)、q3=(1,1) 为固定整数动量索引，未重新找峰。CSV 同时保存 abs(d)、完整复数统计和更保守区间。

分块敏感性仍限制精度：新数据的 Xz SE 在块长25/50/100时为 .00702/.00709/.00817；200/400时只有每链8/4块，不能把此时较小误差当成收敛。Tr S(q1) SE 在100/400时为 .01503/.02157。独立链检查的固定量最大 |z|≈1.48，只是辅助检查，不能证明已混合。

沿用 P5 的分块敏感性包络（只选每链至少16块的层）和四固定观测量修正，q1 Szz/Tr S 的差区间为 [−.04114,.00493]/[−.04277,.05497]，q3 为 [−.02964,.02136]/[−.03990,.05706]。四项均未完全落入原定 ±.02/±.04 尺度；Szz(q1) 的边际约2.2σ也不能单独宣告物理变化。整体投影未收敛，当前数据不支持可信 VBS 峰。

#### 文件与复现

本地原始新数据：`build-kitaev-p01/p5_followup/results/`；实际 PBS 记录：`p5_followup/pbs_record.txt`。远端保留 `/home/sunxr/kitaev_ssh_T64_pilot_20260916.4V2fid`。

- [全部比较及绝对差](../../build-kitaev-p01/p5_followup/analysis/comparisons.csv)、[元数据](../../build-kitaev-p01/p5_followup/analysis/metadata.csv)、[联合统计与分块](../../build-kitaev-p01/p5_followup/analysis/groups.json)。
- [实际漂移与一阶诊断图](../../build-kitaev-p01/p5_followup/analysis/projection_response.png)、[固定动量差值图](../../build-kitaev-p01/p5_followup/analysis/differences.png)、[分块图](../../build-kitaev-p01/p5_followup/analysis/blocking.png)。

完整提交、编译参数在 `p5_followup.pbs`；MC 参数在原 `run_p5.sh` 新增的 `followup_a/b`。在新的远端工作目录上传上述脚本和自由入口后：

```bash
qsub kitaev_ssh/p5_followup.pbs
### 每条链实际调用（PBS 脚本已并行执行）：
P5_BINARY=/home/sunxr/kitaev_ssh_project_smoke_v2.AEKaLn/mc_smoke \
  P5_OUTPUT_DIR=results bash kitaev_ssh/run_p5.sh followup_a
### followup_b 对应 seed 6602；输出拒绝覆盖。
```

下载新数据后，在本地 release 根目录执行：

```bash
P5_PYTHON=/home/asus/.cache/pfqmc_plot_env/bin/python bash kitaev_ssh/analyze_p5_followup.sh
```

分析输出目录拒绝覆盖，重分析需设 `P5_FOLLOWUP_DIR` 指向含同一批 results 的新目录。新增的响应绘图只读取已有联合统计，不另造估计器。相关10项 P5 分析测试通过（.130 s），HPC 自由诊断编译运行通过；未重复 P0–P4 全套，也未改共享数值层。

#### 下一项建议

不自动加到 T=12.8。若优先降低均匀松弛成本，最有信息量的是**均匀移位声子试探态的单键验证，再做原/新试探态的同 T 小对照**：先证明端点权重、自由提议均值与正反提议比一致，再检查能否减少背景漂移；两种试探态的数据不得混合。本轮未实现它。盲目延长该弱耦合 L=2 基准不如上述受控对照；后续也可开展明确标注有限 T 的少量探索，但不能跳过出现实际 VBS 信号后所需的收敛检查。


<a id="report-25"></a>

## CODEX_HANDOFF.md

### fast-1002 optimization handoff

Active optimization workspace: `/home/asus/newsusy/fast-1002`.
Created on 2026-10-02 by copying the complete fast-0929 source tree.
All source files and fixtures were checked against their source SHA256 values.
Only workspace documentation was adapted; no optimization or numerical test
was performed during packaging. No Git repository was initialized.

Start with README.md, AGENTS.md, docs/STATUS.md and docs/CODE_MAP.md.
Keep fast-0929 as the preserved reference. Later new_0923 implementations and
HPC binaries are separate baselines. Build and run evidence for future changes
belongs in a new audit directory here.

The material below is retained historical handoff text. Its paths, defaults,
job states and results refer to their original dates, not this workspace.

---

### Historical `newest_0921` baseline handoff

> 历史文档：下文的版本、目录、默认值和测试结果属于原记录时点，不代表 fast-0929 当前状态。当前快照入口见 [FAST_0929_README.md](REPORTS.md#report-03)，源码核对状态见 [docs/STATUS.md](REPORTS.md#report-26)。历史 audit 未随本包分发；本次未重跑这些测试。

`/home/asus/susy/newest_0921` was the active repository and combined
baseline at the time of this historical handoff (branch `newest-0921`). The input geometry and performance
trees are references only; do not continue development in them.

#### Combined scope

- Accepted Stage 1.1 performance work: prepared trial/generator data,
  persistent numerical workspaces, immutable dense factor sharing, local
  factor reuse, detailed timings, and opt-in deferred phase measurement.
- Integer-sheared honeycomb torus: quotient matrix
  `[[Lx, shear], [0, Ly]]`, shear-aware bonds/gauge paths/momentum grid/FFT and
  geometry metadata. This changes finite-size geometry and is not a continuous
  twisted boundary condition.
- Correct separation of the bare iVBS target from its nearest allowed sheared
  grid momentum.
- Opt-in eager-local signed-prefix checkpoints (`--prefix-cache 1`) and LU
  Cayley reuse (`--cayley-cache 1`). Both remain disabled by default. No
  bidirectional splice, low-rank update, or propagation rewrite is present.
- Numerical candidate recovery: one same-candidate canonical uncached retry,
  then fail closed with replay diagnostics and atomic running/complete/aborted
  result JSON. Deferred measurements retry before adding a sample.

#### Integration decisions

The performance workspace/factor representation is the base implementation.
Geometry changes are layered through the lattice/model/measurement code.
Prefix and Cayley caches use the actual immutable factor storage and preserve
accept/reject ownership. Numerical recovery reconstructs the bookkeeping
required by the next optimized proposal only after its uncached evaluation is
certified.

The `core/` directory from historical commit `dcc101a` is intentionally absent:
it was a duplicate source package frozen at older revision `ba6c1fa`, rather
than another implementation to merge.

#### Build and verification

The clean local verification build is `build-final-check` (ignored by Git).
It used CMake 4.4.3, GCC/gfortran 15.2, Eigen from
`/home/asus/old-p/kitaev_ssh_pfqmc_hpc/third_party/eigen3`, Boost 1.70 from
`/home/asus/.cache/pfqmc_boost_1_70_0`, and oneAPI MKL 2026.1.

Fresh merged-source checks passed:

- `kitaev_ssh_weight_validate`
- `kitaev_ssh_deferred_phase_validate`, including cadence and injected-failure
  checks
- `kitaev_ssh_phonon_validate`
- `kitaev_ssh_mc_validate`
- `kitaev_ssh_sheared_pbc_validate`
- `kitaev_ssh_numerical_failure_validate`, including six CLI recovery/abort
  status cases
- `scripts/test_stage1_tools.py`
- `kitaev_ssh/test_summarize_vbs_candidates.py`
- `audit/numerical-recovery-20260921/check_outputs.py`

Current-binary prefix-on Cayley off/on smokes for shear 0 and 1 produced
identical normalized JSON and byte-identical block CSV. LU calls fell from
2437 to 494 with 1943 hits in each tiny case; prefix resumes and propagated
factor counts were unchanged. Use the retained three-run results in
`audit/cayley-cache-20260921/` for the bounded performance record.

Detailed history and limitations are in
`audit/sheared_pbc/HANDOFF.md`, `audit/sheared_pbc/execution.log`,
`audit/cayley-cache-20260921/REPORT.md`, and
`audit/numerical-recovery-20260921/REPORT.md`.

No HPC work or production-scale physics run was performed for this merge.

#### VBS momentum pipeline follow-up — 2026-09-22

Work continued in the isolated worktree `/home/asus/susy/newest_0921_edit` on
branch `newest-0921-edit`, starting from tagged baseline `6a3745c`. The integer
shear geometry and C++ q-file duplicate policy were unchanged.

`predict_vbs_momenta.py` now separates the complete labeled candidate report
from the parser-clean measurement list. Actual normalized coordinates are
written once in deterministic order; JSON retains all aliases and sources.
Distinct off-grid targets remain distinct when only their nearest grid point
coincides. User `--extra-q` on-grid metadata now uses the actual sheared dual
grid.

An opt-in uniform-background predictor solves
`Xbar=-(g/k)<u i c_A c_B>` using the repository Hamiltonian/f(k) convention and
reports convergence controls, residual, `Xbar` and `Jeff`. It is explicitly a
geometry guide, not a physical iVBS peak. `compare_sheared_geometries.py` can
emit bare/effective nearest-Q mismatch, reduced-torus shape and free `min|f|`
for explicit candidates without assigning a score.

Fresh checks passed: sheared-PBC validator, vbs_q validator, six predictor
tests, three Python-to-C++ q-file end-to-end tests, and selected/both generated
q-file QMC smokes with zero numerical failures. Full commands, raw metrics,
examples and limitations are in
`audit/vbs-momentum-pipeline-20260922/REPORT.md` and its `execution.log`.
No HPC job, production run, cache/projector change or performance benchmark was
performed.

#### Momentum-distance precision follow-up — 2026-09-22

Continued on branch `newest-0921-edit` while preserving all uncommitted work.
The previous report's bare/uniform distance `0.00259669399478971` was a prose
calculation error: it used `(4*pi/sqrt(3))*metric_distance` instead of the
project-basis conversion `(4*pi/3)*metric_distance`. The corrected value is
`0.0014992019768949238`. A shared `momentum_distance()` helper now feeds the
predictor and geometry comparison, and tests cover direct reciprocal vectors,
collinear additivity, periodic equivalence, and Q/-Q handling.

Uniform midpoint integrations at 200, 400, and 800 were recorded with fixed
point residuals below `1e-13`; the 400-to-800 Q change is
`6.99e-7`, while the bare-to-uniform shift is `1.499e-3`. This supports
numerical resolvability of the correction but does not make the uniform
approximation a physical iVBS peak. Fresh `(6,6,0)`, `(12,3,0)`, and
`(12,3,3)` reports and q-files were parser-validated. The detailed append-only
audit is `audit/vbs-momentum-pipeline-20260922-distance-precision/REPORT.md`.
No MC chain, projector/cache path, integer-shear implementation, HPC job, or
performance benchmark was rerun.

#### Static VBS response benchmark — 2026-09-22

Using `J=(1,1,.5)`, `g=.3`, `mass=1`, `omega=1`, `k=1`, the common 800x800
uniform background converged to `Xbar=(.180556937519,.180556937519,.095079563664)`
and `Jeff=(1.054167081256,1.054167081256,.528523869099)`, with
`Q_uniform+=(.580655231756,.419344768244)`, residual `8e-14`, and 42 iterations.
The bare/uniform Cartesian shift is `0.001629579024`.

Added the audit-only `kitaev_ssh/static_vbs_response.py` tool. It scans every
nonzero allowed q on `(6,6,0)`, `(12,3,0)`, and `(12,3,3)`, constructs all
orthonormal x/y/z cosine/sine bond modes, computes matter/phonon/total Hessians,
and compares selected target/soft modes with dense 200x200 and 400x400 Bloch
references. The phonon Hessian is identity; `g=0`, q/-q, self-conjugate-mode,
and finite-energy/Bloch checks pass.

Shear-3 reduces the uniform mismatch to `0.0194302`, versus `0.3496062` for
the two unsheared examples. The softest allowed q's are not universally the
nearest target, and finite/reference eigenvalue differences remain
approximately `0.06–0.08`; momentum matching and static-response finite-size
effects are not yet equivalent. Detailed raw output is in
`audit/static-vbs-response-k1-20260922/`. No QMC or HPC work was performed.

#### Static reference and nearest-tie correction — 2026-09-22

The previous dense response reference was incomplete: its 4x4 Bloch block had
only `A_k -> B_(k+q)` and omitted `A_(k+q) -> B_k`, with an additional
self-conjugate duplicate-counting error. `static_vbs_response.py` now exposes
an explicit-k-list `momentum_reference()` with both real-mode vertices and the
correct generic/self-conjugate factors. On the same finite allowed-k list, all
key total-Hessian matrices agree with the real-space Hessian below `1.2e-16`
elementwise. Corrected 400-grid references are 0.9096558406, 0.9181375452,
and 0.8964976189 for the three highlighted q's.

Modes now use `sum_(r,mu) v_a v_b/Ncell=delta_ab`, and the differentiated
energy is explicitly `E_matter[Jeff+g deltaX] + 1/2 sum_b (Xbar_b+deltaX_b)^2`
per cell. Added regressions for finite-k matrix equality, mixed bond-channel
finite differences, q/-q, self-conjugate q, g=0, and 200/400/800 reference
checks.

Nearest reporting now retains all ties within Cartesian tolerance `1e-10` and
chooses a stable grid-index display representative. On `(6,6,0)`, bare and
uniform targets both have nearest set `{(1/2,1/3),(2/3,1/2)}`; the earlier
single-label switch was not a physical nearest-point change. Full append-only
record: `audit/static-vbs-response-reference-check-20260922/REPORT.md`.

#### Shape-preserving static size comparison — 2026-09-22

Using the common 800x800 uniform background, added the audit-only scan
`audit/static-vbs-size-comparison-k1-20260922/` for A/B/C geometries at l=1,2:
`(6l,6l,0)`, `(12l,3l,0)`, `(12l,3l,3l)`. The scan uses the validated explicit
allowed-k response, shared Cartesian `d_pair(q,{+Q,-Q})`, and 200/400 dense
references, while retaining the prior independently executed 800x800 anchors.
It preserves all q coordinates and nearest/soft ties. Shear-3 keeps pair
mismatch 0.0194302 at both sizes, while unsheared elongated cases are 0.3496062
(l=1) and 0.2930683 (l=2). Same-q response errors are nonuniform, so coverage
and response effects remain distinct; delta decomposition is recorded in JSON.
No QMC, HPC, projector, cache, MC, or shear-geometry changes were made.

#### Stronger-coupling static screen — 2026-09-22

Added `audit/static-vbs-coupling-screen-k1-20260922/`. The parameterized audit
screens g=.3,.6,.9,1.2, refines uniform 400/800 integrations for selected
g=.6 and 1.2, then compares all six l=1,2 geometries using explicit allowed-k
responses and key dense 200/400 references. At g=.6 shear-3 pair mismatch is
0.0146532; at g=1.2 it is 0.00289510, while elongated unsheared controls remain
about .30-.35. Softest q and common fixed q can differ substantially, and
same-q delta_k is not uniformly improved by shear. g=1.2 has positive but
400-to-800 sensitive reference curvature, so it is a stress test rather than a
critical-point claim. No HPC/QMC/projector/cache/MC changes were made.
#### Small PfQMC sign screen — 2026-09-22

Added `audit/pfqmc-sign-screen-k1-20260922/`. Restored the existing C++
executable's oneAPI MKL runtime through `LD_LIBRARY_PATH` without source or
weight-definition changes. Short 2x2 and 3x2 runs at g=.3 and .6, with shear
0/1 and tau=4/8 where applicable, used two seeds per case; all completed with
zero numerical failures and no observed negative or complex samples. This is
only a feasibility result: eight blocks per run cannot establish absence of an
average sign problem. The q-file C++ e2e test passed 3/3 with the restored
runtime. No HPC job or long QMC chain was run.

#### PfQMC lambda=2.5 all-q batch submission — 2026-09-22

Submitted the requested 1000-warmup/3000-measurement, four-seed batch on port
20015 `workq` as PBS array `155565[].mgt`. Parameters are `J=(1,1,.5)`,
`mass=omega_ph=k=1`, `g=sqrt(2.5)`, the requested 15-point tau list, `dt=.1`,
`vbs=1`, `prefix-cache=1`, `cayley-cache=1`, and both `defer-phase=0,1`.
Geometries are `(6,6,0)`, `(12,3,3)`, and the added `(4,4,0)`, for 360 chains.
The production executable is the remote-compatible binary whose
`coupled_mc.cpp` source hash matches the current checkout; failed ABI/build,
header, and CRLF attempts are excluded and documented in
`audit/pfqmc-lambda2p5-submit-20260922/`.

At the latest check, 269 children were running, 43 queued, 46 exiting, and 2
completed; 49 result JSONs had explicit numerical-recovery aborts at large tau.
These failures are retained as diagnostics, not converted into measurements.

The first failure set was audited on 2026-09-22. Its errors are numerical
certification failures in the strong-coupling long-projection path
(`overlap_untrusted`, `zero_weight`, and occasional
`green_residual_exceeded`), with rank/rcond and residual diagnostics preserved;
they are not queue, ABI, or argument-parsing failures. A four-point 20022
probe using `bond-refresh` and `qr-block=8` recovered the short `(4,4,0),
tau=2.0` case, while three longer L6/L12 probes were still active. The 93
aborted rows were resubmitted independently as array `84941[].mgt` on port
20022, with unchanged physical parameters, seeds, caches, defer modes, and
statistics. The retry results are kept separate so they cannot overwrite or
silently replace the original measurements.

At 20:28 CST, the 32 rows already aborted in `84941[].mgt` were resubmitted
with seeds increased by 10000 (`155567[].mgt` on port 20015 `workq`). The
retry keeps `bond-refresh`, `qr-block=8`, both caches, defer mode, and the
1000/3000 schedule unchanged. It uses the isolated directory
`/home/sunxr/pfqmc_lambda2p5_seedretry_20260922_20015`; its initial check had
20 children running and 12 early numerical aborts, so changing the seed is
being treated as a stochastic retry, not as evidence that the certification
problem is fixed.

#### Fixed-QD extended stability run — 2026-09-24

The four-chain L=6 fixed-window audit completed as PBS array `155676[].mgt`
on port 20015 `workq`. Seeds 923005–923008 each have 1600 measurements,
160 ten-sample blocks, zero numerical failures, and phase +1 on every sample.
Parameters remain fixed at `J=(1,1,.5)`, `g=.3`, `tau=.2`, `dt=.1`, and the
same ten off-grid QD/neighbor windows.

The full-chain center `Szz/Nc` values are 0.013522, 0.013638, 0.015070, and
0.013690. `R_b1` is −0.0176, −0.0016, +0.0271, and −0.0126; `R_b2` is
−0.0109, +0.0107, +0.0416, and −0.0094. Only seed 923007 has both ratios
positive, while half-chain and quarter estimates move within/across seeds.
Reblocking errors also remain sensitive, and the sampler does not verify block
independence. The longer data therefore do not yet establish a reproducible
fixed-QD local peak. They test sampling at `tau=.2` only, not projection-length
convergence.

The analysis is in
`audit/vbs-window-stability-20260924/REPORT.md`; collected raw results, checksums,
segment/reblocking tables, and full complex fixed-window matrices are retained
alongside it. Continue at the fixed point and monitor phase, `Delta(0)`,
connected `Szz`, trace, and matrix channels before changing projection settings.


<a id="report-26"></a>

## docs/STATUS.md

### fast-1002 状态与历史记录

#### 当前快照核对：2026-10-02

适用范围：`/home/asus/newsusy/fast-1002` 优化工作区，完整复制自 fast-0929，来源包标记为 `2e53aee`。
2026-10-02 完成逐文件复制一致性检查，尚未进行优化、构建或数值测试。
生产入口仍是固定 g 的 `kitaev_ssh_mc`；CLI 默认 `mass=4`、`omega-ph=0.5`、
`proposal=local`、`weight-engine=combined`，延迟相位、prefix/Cayley 缓存开关均为 1。
缓存实际生效范围由引擎和提议路径决定。支持声子采样、剪切边界、VBS 测量及
同候选无缓存重试；尚未从生产入口看到 SSH 声子时间驱动协议接入。

本次证据为来源源码阅读及复制 SHA256 核对，没有新构建、回归、性能或物理结果。后续 `new_0923`
工作树和 HPC 作业需单独记录；不以其结果更新本快照能力。
当前入口见 [README.md](../../README.md)，复制记录见 [REPORT.md](../../audit/core-copy-20261002/REPORT.md)。
[FAST_0929_README.md](REPORTS.md#report-03) 保留来源快照说明。

#### 历史状态：2026-09-24

以下内容是原工作树历史记录；“当前”“本轮”均指原记录时点，
其中 audit 路径未随本包分发，运行中任务状态也未在此次重新查询。

状态范围：截至 2026-09-24 对 `/home/asus/newsusy/new_0923` 的本地核对；运行结果仍按其原 audit 的证据等级记录。
来源：[H、C0、C1、R-EXPORT、R-SAMPLING、R-CONTOUR-REGRESSION、R-MIDPOINT-REGRESSION、U-CORRECTION、U-GOAL](../SOURCES.md)。


#### 证据边界

- **论文事实**：仅作为研究目标和比较背景。
- **数学推导**：只记录影响当前状态判断的已写明关系。
- **源码现状**：当前工作树核对结果，随版本变化。
- **报告结果**：历史运行和测试按来源引用，绝不改写为本轮重新通过。
- **未验证计划**：优先问题和草案不构成执行授权或完成声明。


#### 当前研究阶段

当前首要阶段是完善代码并验证代码正确性：优先做局部接口核对、基准、回归和数值不变量检查。固定参数 SSH–Kitaev 声子 PfQMC 和 VBS 测量基础仍在整理；有限虚时间驱动/有限时间标度只是长期背景，不要求每项工作都推进到该目标。尚未收到 SSH–Kitaev 线性驱动专项核对完成报告或生产结果。

此前文档复核轮次仅授权写说明；本轮额外完成了固定-g 回归和 validation-only 中点协议回归。仍没有授权生产驱动接入、下一轮采样、HPC 或 commit。

#### 驱动源码复核（2026-09-24）

本地现场核对与 [源码复核](../reviews/DRIVEN_SOURCE_REVIEW_20260924.md) 的关键判断一致：Stage 1 是 spinless tV/HS 基线；生产 SSH 入口仍固定 `g`，`MajoranaGenerator` 只接收 X，模型新增的显式 `g` 重载仅用于 validation；`coupled_mc.cpp` 没有驱动 CLI。现有 `PhononPath` 已能表示总开放轮廓，具有两外端 trial 和唯一中心，内部两侧节点不强制镜像。

优先设计改为“总 `PhononPath` + 协议映射 + 唯一中心”；两个完整路径对象只是待比较的软件表示。若采用区间中点 `g`，相邻区间在同一 X 节点的半步一般不同，固定 g 下的 node 平方、谱/Cayley 缓存和恢复重建不能直接沿用。第一版可暂不支持这些优化，但必须显式限制并先做 canonical 参考验证。[R-DRIVEN-RECHECK]

#### 固定-g 总轮廓回归（2026-09-24）

针对当前源码重新配置、构建并运行了已有 `kitaev_ssh_p2`（`kitaev_ssh/phonon_validate.cpp`），结果 1/1 通过。当前测试确认非镜像节点、唯一中心、两外端路径语义、Fock 有序乘积、中心 Green、固定-g 优化路径和 `dt` 二阶收敛；实际日志和二进制校验值见 [contour-regression-20260924](../../audit/contour-regression-20260924/REPORT.md)。这不是驱动实现验证，也不是 QMC 结果。

#### 中点协议 canonical 回归（2026-09-24）

为验证显式瞬时 `g` 的最小矩阵接口，新增 `frameworkGeneratorWith(X,g)` 重载，并在已有 `kitaev_ssh_p2` 中加入 validation-only canonical 中点路径：9 组协议（8 组线性 ramp、1 组 constant）、非镜像 2×2 小系统，与 Fock 有序乘积比较 matter 对数幅度、相位和中心 Green；恒定协议还逐因子回归固定-g 路径。直接验证和 CTest 均通过，详见 [midpoint-driven-regression-20260924](../../audit/midpoint-driven-regression-20260924)。自由声子抽样统计提高到 100,000 次，最大标准误差单位偏差 3.41，低于 6σ 门槛。

这只确认显式 `g` 组装和 canonical 因子顺序；没有接入 `MajoranaGenerator` 的时间上下文、`coupled_mc` 驱动 CLI、双支采样、缓存/恢复优化或生产测量，因此不能称为 SSH–Kitaev 驱动已实现。

#### 已有证据支持的工作

固定规范模型、声子开放路径、Pfaffian 相位、local/bond-refresh、数值失败参考重算有源码/历史验证记录；完整 3×3 VBS 全网格和指定非网格测量已有实现。窗口平移平均、connected/projected 和联合区块分析已有数学关系与测试记录。[H、C0]

动量导出增量已由用户报告合并，相关 40 项测试通过；真实短链全网格/路径导出已完成。这里引用历史报告，不是本次文档编写重新运行的测试。[C1、U-MERGE、R-EXPORT]

#### 最新实际采样基线

R-SAMPLING：HEAD 4130a2e…，J=(1,1,0.5)、g=0.3、m=ω=k=1、(Lx,Ly,s)=(6,6,0)、dt=0.1、intervals=2、tau_total=0.2、Wilson=(1,1)、local、combined、both。每 seed 100 热化+400 测量，10 samples/block，共40块。

| seed | 完成/样本/块 | 墙钟秒 | 数值失败 |
|---|---|---:|---:|
| 923003 | complete / 400 / 40 | 1202.64 | 0 |
| 923004 | complete / 400 / 40 | 1196.10 | 0 |

这两条确实有 HPC 完成记录（155649[].mgt）。早期“未提交 HPC”的 handoff 状态仅对当时有效，不能覆盖此后发生的运行。两条平均相位约1，不是驱动优势证据。[R-SAMPLING]

#### 最新 correction：中心值未变，误差与配对已纠正

以下 Q± 的 connected 值均除以 Nc=36，SE 统一为 10 samples/block 的删块 jackknife：

| seed | Szz/Nc ± SE | Tr S/Nc ± SE | λmax/Nc ± SE |
|---|---:|---:|---:|
| 923003 | 0.013506 ± 0.000615 | 0.042121 ± 0.001286 | 0.015774 ± 0.001039 |
| 923004 | 0.013817 ± 0.000593 | 0.044504 ± 0.001390 | 0.016767 ± 0.000829 |

seed923004 旧主表误用的 (0.000668,0.001690,0.001048) 属于20 samples/block，不再标作10/block。四个邻近点单点估计没有算错，但旧配对归属交叉；按 PHY-06 的正确反演对使用。[R-SAMPLING + U-CORRECTION]

本地 correction 产物位于 [`audit/vbs-momentum-sampling-extension-L6-20260923/correction_20260923/`](../../audit/vbs-momentum-sampling-extension-L6-20260923/correction_20260923)，其中有报告、重算脚本和 CSV；原始数据与 `analysis_v3` 未修改。本轮核对了报告、脚本、输出结构和来源路径，未重跑 QMC 或历史计算。用户摘要仍作为运行事实的来源标签，具体修正证据以本地 correction 文件为准。

#### 当前不能作出的结论

固定 tau_total=0.2 下采样稳定性尚未证实；不能反向断言已经证明未热化。固定通道、零模和本征值需分开诊断；只有10样本区块，不能恢复逐 sweep 自相关或有效样本数。没有基态投影收敛结论。[U-CORRECTION]

尚无可用证据证明：SSH–Kitaev 有限速率驱动已接入并验证、R 与 2R 交点、驱动标度窗口、符号/成本优势、L=9/12 的驱动临界标度、SUSY 数值验证或完整物理自旋投影。不能把通用驱动头文件、构建目标或旧建议当作完成证据。

#### 当前优先问题，而非自动待执行命令

当前完整工作树已有 [`reproduction/pure_state_driven_stage1/`](../../reproduction/pure_state_driven_stage1) 和 CMake 的 `pure_state_driven_qmc` 目标；其中 `driven_qmc.cpp`、`stage1_common.h`、`stage1_dense.h` 可作为驱动基线材料。它们尚未证明已适配 SSH 声子生成元、双支时间和中心测量。DRV-02 的 g 线性协议仍是草案；生成元时间接口及离散规则待定。更细粒度采样诊断也是缺口，但不应把所有终点长投影基态收敛设为开始驱动研究的门槛。[U-GOAL、DISCUSS]

已有局部文档/脚本和 HPC playbook 只在本地报告过；实际路径、HEAD、未提交修改由接手者读取，不从此快照猜测。

#### 维护规则

新证据到来时更新本页日期、适用版本和来源；历史表不覆盖旧 audit。模型定义见 PHYSICS_CONVENTIONS，估计器见 OBSERVABLES，驱动方案见 DRIVEN_FTS_METHOD，避免在本页复制完整推导。

#### 固定窗口长链波动原因检查（2026-09-24）

四条1600样本长链的十样本基础块存在可测残余相关：C/Nc相邻块相关0.31–0.55，R影响序列0.18–0.44。单节点条件更新的自由三节点解析参考也有慢模；高接受率不能排除此问题。更粗分块明显缓解C和R1链间张力，R2仍有一定张力。符号翻转不能单独证明未热化；当前没有证据认定严重采样偏差。独立求和/比值路径与导出值相差至多2.22e-16，均值背景对R的影响至多2.59e-5。详情及复现见[原因检查](../../audit/vbs-stability-cause-20260924/REPORT.md)。本轮没有新QMC或HPC运行。

#### local 与 bond-refresh 高统计量对照已提交（2026-09-24）

针对上述慢相关诊断，已提交固定 `L=6, tau=0.2` 的提议算法对照：PBS 数组
`155698[].mgt`，20015 `workq`。local 与 bond-refresh 各 4 条独立链，每条
800 热化、6400 测量，逐测量保存，共每种提议 25600 个测量点。初次和
05:17 UTC 的展开查询均显示 8 个子任务运行中；结果尚未收集，不能据此
更新采样或收敛判断。提交、采集和预先固定的统计口径见
[`audit/vbs-local-vs-bond-refresh-L6-20260924/`](../../audit/vbs-local-vs-bond-refresh-L6-20260924)。

#### VBS 汇总计算优化（2026-09-24）

本轮缓存固定几何的逆 Fourier 相位、仅求所需本征值并复用统计临时缓冲。
L=6 固定输入的6400次 derive 局部基准中位耗时从0.611秒降到0.241秒，
该输入输出完全一致；现有 p3、p4、vbs_q_validate 三项回归通过。
尚无整条生产链加速数据；155698数组使用原二进制。
细节见[效率审计](../../audit/vbs-efficiency-20260924/REPORT.md)。

#### local 与 bond-refresh 对照已完成（2026-09-24）

数组 `155698[].mgt` 的 8 条链均完成：两种提议各 4 条、每链 6400 测量、无数值失败。
对照显示 local 影响序列 lag-1 约 0.8，bond-refresh 约 0；local 的 40 测量重分块误差约为逐测量误差的 3 倍，而 bond-refresh 基本不膨胀。bond-refresh 在固定点的单位测量时间约快 2.4 倍。该证据支持此前的大波动主要来自 local 慢混合和误差低估，但尚不提供投影长度收敛证据。详见
[`audit/vbs-local-vs-bond-refresh-L6-20260924/REPORT.md`](../../audit/vbs-local-vs-bond-refresh-L6-20260924/REPORT.md)。

#### tau=0.2 与 tau=0.8 投影扫描已完成（2026-09-24）

数组 `155718[].mgt` 和 `155719[].mgt` 均完成，分别为 bond-refresh 的 L=6、tau=0.2
四链和 tau=0.8 六链，每链 2000 测量。两点中心 `Szz/Nc` 均约 0.014，`R_b1/R_b2`
均与零相容；这是有限的两点稳定性证据，不是 tau 无穷收敛证明。tau=1.6 数组
`155720[].mgt` 当时仍在运行。详见 [`audit/vbs-tau-scan-L6-20260924/REPORT.md`](../../audit/vbs-tau-scan-L6-20260924/REPORT.md)。

