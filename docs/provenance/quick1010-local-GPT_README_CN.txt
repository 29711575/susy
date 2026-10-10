【优化候选版入口】本目录在原快照上新增可选 guided-bond / q-pool 采样。先读 OPTIMIZATION_README_CN.md 和 docs/SAMPLING_OPTIMIZATION_RESULTS_20261010.md。以下是原快照说明，不是本次优化的测试状态。

fast-1006 当前核心源码快照 — GPT 阅读入口

当前工作区为 quick1010，后续代码修改与优化使用本目录。
2026-10-10 已整合 fixed guided-bond pCN 与固定 q 池扩展；当前原生 MKL
验证及来源以 docs/STATUS.md 和整合审计 REPORT.txt 为准。下文为继承的
原包阅读说明与历史物理分析背景。

打包目的：分析 SSH–Kitaev 声子 Pfaffian projector QMC 的物理、数值和采样实现，尤其是 VBS 序参量与静态极化率。
本包取自当前 fast-1006 工作区，保留原目录结构；不是 fast-1002 旧包，也不是已提交作业使用的统一版本。既有生产作业保留各自冻结源码/二进制；当前代码含后续优化，不能认为所有图中的链均用了当前更新算法。

推荐阅读顺序：
1. docs/CODE_MAP.md、docs/PHYSICS_CONVENTIONS.md、kitaev_ssh/CONVENTIONS.md。
2. kitaev_ssh/coupled_mc.cpp：生产 CLI、参数、初态、测量与输出。
3. kitaev_ssh/kitaev_ssh_model.h、honeycomb_lattice.h、phonon_path.h：Hamiltonian、几何、开放投影路径与 Gaussian trial。
4. kitaev_ssh/phonon_projector.h、phonon_mc.h、phonon_q_mode.h：权重、提议、接受/恢复及集体模式更新。
5. inc/pure_projector_*.h：Majorana 权重、phase、Green 与稳定化。
6. kitaev_ssh/vbs_observables.h、vbs_susceptibility.h、phase_statistics.h；docs/OBSERVABLES.md、NUMERICAL_INVARIANTS.md、decisions/ADR-009-vbs-static-response.md。
7. docs/decisions/ADR-014-fixed-q-collective-update.md、ADR-016-lazy-canonical-half-steps.md；CODEX_HANDOFF.md 用于最新优化状态。README.md 中旧构建名称不代表当前最新二进制。

本次讨论的物理参数：L=9，shear=0，J=(1,1,0.5)，omega=0.3，k=1，mass=100/9，lambda=g^2/k，dt=0.1；tau 为总开放轮廓长度，中心距端点 tau/2。tau=8 对应 intervals=80。X 没有按论文的 gX 重标度。

序参量口径：m_z^2=max_{q!=0} Re Szz(q)/81；先合并相位加权矩，再搜索 80 个非零允许 q。OBS03 connected 均匀背景在这些非零允许 q 上为0。
静态极化率口径：chi_zz(q)=<ChiRaw(q)>-81*<ChiD(q)>*<ChiIconj(q)>；一阶/二阶矩共用 phase 分母。ChiIconj 独立保存，不能用重加权均值的共轭替代。chi_max=max_{q!=0} Re chi_zz(q)。不是 tau*S(q)，也不减自由声子响应。

最新 tau8 极化率图：每条链弃去测量 sweep1–800，原作业 warmup=0；800 是后处理切除，不是已认证的热化时间。27 lambda 点，203条有保留样本的独立链，每点4–10seed，251470保留样本。部分链已停止，所有点仍属部分链快照；跨 seed 删除整链 jackknife1SE，副本重算连接减法和动量峰；未认证热化、自相关或基态投影收敛。

analysis_context/ 保存本轮极化率后处理脚本和27点汇总，用于核对从生产原始矩到绘图的统计逻辑；不包含逐链宽CSV或远端登录资料。
用户希望重点分析：为什么 tau8 下极化率随 lambda 只有缓慢增强、没有清晰临界峰；应区分物理有限尺寸/投影、动量与通道选择、慢混合、统计误差以及潜在代码/归一化问题。不要预设必须有峰，也不要把未出现峰直接判为无相变或代码错误。

构建：CMakeLists.txt 或 scripts/build_nodeps.sh；需要外部 Eigen、Boost、BLAS/LAPACK/MKL 等依赖，按原构建文档配置。PFAPACK 源码/许可证和确定性验证程序已包含。本次只做打包及内容校验，未重新编译或运行数值测试。

排除：audit 历史数据、build-*、编译产物、远端凭据、原始生产CSV。docs 中历史 audit 链接可能指向包外材料。
校验：SHA256SUMS 给出每个打包文件的 SHA256。SNAPSHOT.json 给出打包时间及来源。
