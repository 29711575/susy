# ADR-015：预热路径诊断、矩阵存储与生产负载 profiling

日期：2026-10-09。用户授权代码优化，指定测试 L6、HPC20015。

此前 coupled_mc 无条件传入 sample 的 onWarmup 回调。sample 的这一
接口承诺先 ensureMeasured，因此即使诊断关闭，每个预热 sweep 仍
认证 signed phase 和配置所需的中心 Green。开启诊断实际输出的
rawAllQMax 只依赖中心声子路径，不需要这两种计算。

CLI 改从既有 onSweep 回调读取路径并按诊断 interval 输出，传给
sample 的 onWarmup 为空。sample 的公共接口及其他调用者的语义保持。
正式测量相位重加权、RNG、提议、接受判据、trial/force/Hamiltonian、
精度门槛和需要进行的提议/测量失败恢复均保持。预热不再主动执行
未被观测量使用的 signed 检查；该诊断仍是 raw path，不能按物理平均解释。

新增 kitaev_ssh_production_profile 直接复用 coupled_mc 的完整 CLI；
开启编译时 PathProfile/CoreProfile，额外写 PREFIX_profile.json。
原 kitaev_ssh_profile 的 L4 legacy 接口保留。新计时覆盖初始化、
预热、测量，区分整键/q 候选计算，并保留生产 timing 的阶段耗时。
嵌套计时为 inclusive，不能相加声称总耗时；生产 MC 不开启 profiling。

首批 L6 生产 profiling 定位到 factor_generation_and_storage 约占完整
测试耗时的三分之二。本轮随后为 validated dense factor/slice 增加
rvalue 所有权转移，在 prepareOwnedPath 完成整步矩阵检查后将 weightB
直接移入不可变 owner，减少一次 dense 复制。原 lvalue、generic、
有限性/Hermitian 和 identity residual 检查保持；不改变矩阵系数、
传播顺序、缓存绑定或数值门槛。

第二批 L6 长用例采用 warmup=0，单独测量整步存储改动对正式采样的
影响；复用已验收的同工具链 baseline，另构建最终程序、production
profiler 和既有缓存事务验证器。首批候选和所有历史证据保留。

修改前与候选源码、冻结 manifests、L6 配对测试和运行证据位于
`audit/1009/code-speed-L6-20261009/`；验收与实测速率由该目录报告记录。

本轮20015实际验收：两批13个子任务全部exit0、PASS，配对最终路径与
接受/拒绝计数一致，全部raw blocks逐字节一致；既有缓存事务验证器
覆盖624次接受/240次拒绝通过。第一项预热wall节省中位2.64%，其
4预热+4测量长用例的整段CPU节省中位0.78%；最终warmup0的L6/M140
四次正反顺序配对CPU节省3.97%–5.19%，中位4.46%。这不是长链或L9
速度认证，也不是收敛或等误差采样收益。完整证据见
[本轮报告](../../audit/1009/code-speed-L6-20261009/REPORT.txt)。
