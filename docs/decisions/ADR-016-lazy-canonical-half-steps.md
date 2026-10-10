# ADR-016：延迟测量提议按需构造半步矩阵

日期：2026-10-09。用户继续授权代码优化，沿用 L6、HPC20015 测试。

上一轮生产负载 profiling 显示传播矩阵生成和存储仍占约65%。延迟
phase 且合并节点传播时，接受判据使用内部节点 exp(-dt K) 和两端
exp(-dt K/2)。此前还为每个内部节点生成 exp(-dt K/2)，存入 canonical
factors；正式测量本来就重新生成完整 factors 和 Cayley，因此提议
阶段的内部半步矩阵没有被使用。

prepareOwnedPath 和 prepareOwnedNode 增加 wantHalfSteps 参数，默认 true，
保持直接调用者原接口行为。CoupledPhononPath 的延迟合并提议选择 false。
内部节点只生成原谱分解的整步传播矩阵；两端半步仍生成并认证。提议
快照 factors 为空，weight_order 完整且不可变，phase/Green 无效。
evaluatePathMagnitude 仅在 prevalidated 且合并顺序长度为 nodes 时接受
空 factors；普通、未合并和 signed/Green 调用仍要求 canonical factors。

ensureMeasured 继续在提交新快照前重建并认证完整 canonical factors、
signed 权重和请求的 Green。从测量快照或 uncached 恢复快照继续局部
更新时，候选只更新合并顺序并回到空 canonical 表示，不依赖旧半步
槽位。拒绝或认证失败仍保留原快照；提议、RNG、接受和测量定义保持。

生成元、谱分解和所有实际使用的传播矩阵保留原 finite/Hermitian 检查
及1e-12门槛。未被使用的内部半步矩阵及其检查延迟至测量。若直接
整步指数非有限，仍按原策略生成并验证该节点半步、尝试 h*h；重算
失败仍触发同候选及同 u 的一次 uncached 重试，随后中止。uncached
参考路径完整保留，不采用 lazy 因子、接受态缓存或放宽数值门槛。

production profiler 增加 gen_half_assemble 计数，嵌套计时仍为 inclusive。
针对性检查位于既有 pfqmc_cache_validate：L6 下比较 eager/lazy 合并
系数、普通与非二分谱分解、拒绝深值不变、测量/局部/集体/恢复状态
切换，以及不完整因子输入拒绝。另使用既有 numerical_failure_validate
检查一次恢复、恢复后 Metropolis 拒绝、持久失败中止和测量恢复。

冻结源码与本轮运行证据位于 audit/1009/lazy-half-L6-20261009/。
速度结论须以本轮实际配对结果为准，不从省去矩阵个数推断加速倍率。

2026-10-10 验收：20015构建156965、数组156966/156967的5个L6子任务
全部exit0、PASS。L6状态切换、624接受/240拒绝缓存事务和既有失败
恢复验证通过；配对最终路径/计数精确一致，全部raw blocks逐字节
一致，物理JSON/CSV最大缩放差0。相对上一轮已验收版本，L6/M140
四次正反顺序配对CPU减少8.55%–12.16%，中位10.48%。M140/q0的
4测量sweep profiler记录1430次半步构造，原逐节点方式对应61617次。
这仅是短工程验收，不认证L9、长链速度或等误差采样收益。完整证据
见[本轮报告](../../audit/1009/lazy-half-L6-20261009/REPORT.txt)。
