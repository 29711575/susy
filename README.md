# SSH–Kitaev 声子投影 PfQMC：quick1010 源码与极化率交接

本仓库是 2026-10-26 同步的 `quick1010` 主线源码整理版，供合作者检查计算方法、采样和测量。源码放在仓库根目录；旧仓库中的 `fast-1002/` 内容已由本版本替换，原提交历史仍可回溯。


2026-10-26 更新：默认 physical matter trial、长投影 Gaussian/Pfaffian 数值修复、
可选 Γ 更新、信号检查点/续算、完整能量估计器和离线分析工具已同步。
**总能量默认关闭，旧 matter energy 保留；Γ 更新默认关闭。**
见 [本次改动及兼容说明](docs/updates/20261026-source-sync.md) 和
[实际导出验证](evidence/source-sync-20261026/README.md)。
历史测量仍对应各自冻结版本；以下 2026-10-13/10-10 说明保留为历史记录。


2026-10-13 源码更新：增加显式 `--center-green-fast 1` 中心 Green 重建加速及严格静态全键响应工具。guided 整键和固定 q 池实现沿用上一版。旧版完整保存在标签 [`quick1010-20261010-before-center-green`](https://github.com/29711575/susy/tree/quick1010-20261010-before-center-green)；详见 [更新与验证说明](docs/updates/20261013-source-update.md)。下列测量仍是 2026-10-10 的历史快照，各自冻结版本保持原记录。


数值配置已进一步统一到核心 `production-v1`：普通脚本无需再写 `--weight-engine`、`--defer-phase`、`--prefix-cache`、`--cayley-cache` 或 `--center-green-fast`。参考诊断用 `--numerics reference`；旧逐项参数仍支持显式覆盖。guided、q 池与物理/观测参数保持独立。见 [ADR-022](docs/decisions/ADR-022-core-numerics-policy.md)。此次重构前版本保留在标签 [`quick1010-20261013-before-numerics-policy`](https://github.com/29711575/susy/tree/quick1010-20261013-before-numerics-policy)。

**先读 [交接报告](HANDOFF_REPORT_CN.md)**。目前极化率曲线没有已确认的临界峰；尚未确认生产链热化、遍历或投影收敛。这里提供的数据来自较早冻结的生产版本，不能当作当前 quick1010 新更新算法的长链效果。

## 代码入口

| 任务 | 文件 |
|---|---|
| 生产 CLI 与结果写出 | [kitaev_ssh/coupled_mc.cpp](kitaev_ssh/coupled_mc.cpp) |
| MC 更新与接受率 | [kitaev_ssh/phonon_mc.h](kitaev_ssh/phonon_mc.h) |
| guided 整键提议 | [kitaev_ssh/phonon_guided_proposal.h](kitaev_ssh/phonon_guided_proposal.h) |
| 中心 Green 加速与验证 | [inc/pure_projector_endpoint.h](inc/pure_projector_endpoint.h)、[验证器](kitaev_ssh/center_green_validate.cpp) |
| 严格静态全键响应 | [工具](kitaev_ssh/static_limit_response.py)、[用法](STRICT_STATIC_README_CN.txt) |
| 集体动量更新 | [kitaev_ssh/phonon_q_mode.h](kitaev_ssh/phonon_q_mode.h) |
| 线性 connected 极化率 | [kitaev_ssh/vbs_susceptibility.h](kitaev_ssh/vbs_susceptibility.h) |
| 新增模长矩与模长源响应 | [kitaev_ssh/vbs_amplitude.h](kitaev_ssh/vbs_amplitude.h)、[用法](VBS_AMPLITUDE_README_CN.txt) |
| 模长/Binder 分析 | [kitaev_ssh/analyze_vbs_amplitude.py](kitaev_ssh/analyze_vbs_amplitude.py) |
| 几何、物理及统计约定 | [docs/PHYSICS_CONVENTIONS.md](docs/PHYSICS_CONVENTIONS.md)、[docs/OBSERVABLES.md](docs/OBSERVABLES.md) |

## 2026-10-10 测量快照

| 数据 | 范围与口径 | 入口 |
|---|---|---|
| L9，τ=8，χ 随 λ | λ=1.3–5 的27个已有点；203条有效链；各链截500 sweep | [图](measurements/L9_tau8_lambda/chi-vs-lambda-tau8-cut500.png)、[CSV](measurements/L9_tau8_lambda/chi-vs-lambda-tau8-cut500.csv) |
| L9，χ 随 τ | 5个λ、17个τ；τ>6截1000，其他不截；追加τ10/12/14生产链 | [普通图](measurements/L9_L6_tau_scan/chi-vs-tau-L9.png)、[双对数图](measurements/L9_L6_tau_scan/chi-vs-tau-L9-loglog.png)、[固定q](measurements/L9_L6_tau_scan/chi-vs-tau-L9-fixed-q49.png) |
| L6，各τ的χ随λ | τ10/14；各链截3000；τ10含λ2.85/2.9/2.95细扫描 | [τ10](measurements/L9_L6_tau_scan/chi-vs-lambda-L6-tau10.png)、[τ14](measurements/L9_L6_tau_scan/chi-vs-lambda-L6-tau14.png)、[对照](measurements/L9_L6_tau_scan/chi-vs-lambda-L6-all-tau.png) |
| L9，τ8动量热图 | 六边形第一BZ、裸狄拉克点差标星；非零q峰标圈；包含Γ | [λ1.3–5总览](measurements/L9_tau8_heatmaps/chi-heatmaps-lambda1p3-5.png)、[λ2.8–3细扫描](measurements/L9_tau8_heatmaps/chi-heatmaps-lambda2p8-3-step0p025.png) |
| full / disconnected 分解 | 同一q核对χ=Cfull−Cdisc；含逐seed一阶矩 | [随τ](measurements/connected_decomposition/L9-full-disc-vs-tau-fixed-q49.png)、[随λ](measurements/connected_decomposition/L9-tau8-full-disc-vs-lambda.png)、[详细表](measurements/connected_decomposition/decomposition.csv) |

## 构建与离线复现

C++17、GNU C/C++/Fortran、Intel oneAPI MKL、Eigen3、Boost headers；PFAPACK源码及原许可证随仓库提供。依赖包不随仓库上传。

```bash
export EIGEN3_INCLUDE_DIR=/path/to/directory/containing/Eigen
export BOOST_INCLUDE_DIR=/path/to/directory/containing/boost
export MKLROOT=/path/to/oneapi/mkl/version
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
bash scripts/build_nodeps.sh build kitaev_ssh_mc \
  kitaev_ssh_vbs_susceptibility_validate kitaev_ssh_vbs_amplitude_validate
export LD_LIBRARY_PATH="$MKLROOT/lib:${LD_LIBRARY_PATH:-}"
./build/kitaev_ssh_vbs_susceptibility_validate
./build/kitaev_ssh_vbs_amplitude_validate
```

Python 后处理需要 NumPy、Matplotlib；不需要 HPC 登录。冻结的 gzip JSON 是各独立链截断后的共同相位分母与原始矩之和，可以重算均值和删整链误差；它们不包含逐 sweep 时间序列，不能据此估算链内自相关。

```bash
python3 measurements/L9_tau8_lambda/plot.py
python3 measurements/L9_L6_tau_scan/plot.py
python3 measurements/connected_decomposition/analyze.py
python3 measurements/connected_decomposition/seed_check.py
```

普通生产设置显式使用 `bond-guided/bare/beta1` 和集体q；核心数值默认提供 combined/defer-phase/prefix/Cayley，见 [AGENTS.md](AGENTS.md)。2026-10-13 的更新加入数值加速，当时保持原 trial、采样和观测定义；2026-10-26 同步进一步改变默认物理 trial 并增加可选更新及测量，详见顶部更新说明。既有 HPC 任务保留各自冻结版本。普通新生产的核心默认已启用中心 Green 加速；有意对照仍可显式写 `--center-green-fast 0`。

[源码来源清单](SOURCE_MANIFEST.json)、[本次导出验证](evidence/export-validation/README.md)、[历史验证](evidence/recorded-validation/) 分别记录不同证据。`docs/provenance/` 保存原本机文档；其中历史路径不是本仓库的相对路径，不能把历史测试记录误认为本次重跑。
