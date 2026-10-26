# ADR-027：信号触发的 Monte Carlo 检查点与续算

状态：实现及本地回归完成；目标编译与核验记录见 `audit/1023/signal-checkpoint/`。
用户于2026-10-23要求在手动停止前保存当前构型，以便从中断点继续。

## 接口与安全边界

生产CLI默认安装SIGUSR1、SIGUSR2、SIGTERM处理器。信号处理器仅置
`volatile sig_atomic_t` 标记，不分配内存、不写文件、不修改采样状态。
采样器在完整sweep、相应测量与统计更新均完成后检查请求：

- SIGUSR1：原子保存，随后继续运行。
- SIGUSR2：原子保存，随后以75退出；结果状态为`checkpointed`，不能作完整生产结果。
- SIGTERM：与SIGUSR2相同，但外部强制杀死的宽限期可能短于当前sweep；仅作尽力保护。
- SIGKILL、节点失联和进程崩溃不能触发保存；有需要时显式设置定期检查点。

`--checkpoint-file FILE.dat`默认为`PREFIX_checkpoint.dat`；
`--checkpoint-every N`默认为0，表示仅信号触发，N>0时每阶段每N个完整sweep保存。
检查点仅允许.dat后缀，防止与CSV/JSON/path输出冲突；最终汇总期间的请求也会在结果发布前处理。
没有请求时不额外抽随机数、不增加测量、不更改权重或数值门槛。
`--resume FILE --output NEW_PREFIX`继承检查点记录的选项并锁定物理、
数值、更新及测量设置。可以修改`--sweeps`为**全链总目标**，不是追加数量；
必须包含已完成的measurement sweeps并满足原stride/block约束。
允许改变输出/检查点位置及定期保存间隔。恢复必须使用同一可执行文件字节；
还应保留原库运行环境，不保证跨libstdc++/MKL版本迁移。

## 保存内容

检查点是版本化、有长度与FNV-1a损坏校验的单文件，包含：
完整M+1虚时间节点（含端点）、声子参数、trial位移和外力；
mt19937_64、normal_distribution缓存、当前随机排列；
warmup/measurement进度、stride间隙进度、q/Gamma cadence；
全部基础/q池/Gamma的累计与测量计数、统计计时；
全部原始相位分子/分母区块及未满区块；
warmup计数基线、数值失败/恢复计数、原选项、q池/选定VBS动量文件内容。

文件先写同目录独占临时文件，写完fsync并关闭，rename替换，再刷新目录。
无法保存时终止并报告失败，不发布成功receipt；前一完整检查点在rename前保持有效。
校验用于发现损坏，不是防恶意篡改的真实性证明。

缓存、投影器工作矩阵和Green不保存。恢复时从坐标重建并通过原严格数值门槛，
失败只对同一构型做既有参考恢复，不允许重新随机抽一个初态替代。
RNG、随机排列及统计状态恢复不意味着一般相互作用路径在缓存重建后保证
逐比特相同：舍入差异可能在后续MH阈值附近改变接受分支。
本地g=0信号中断/恢复逐字节回归和有限g短链回归提供实施证据，不能推广为
所有历史构型/编译器的逐比特续算证明，也不认证热化或采样收敛。
恢复构型重建期间发生的数值恢复计数不能被继承计数覆盖。

## 输出与重复计数

恢复写入**新输出前缀**，不覆盖祖先输出或检查点。新`_blocks.csv`重写已继承
完整块，然后续写；未满块保存于检查点并在续算中补齐，绝不丢弃或重复测量。
最终results中的samples、原始块、接受计数为全链累计。
诊断/接受率/混合轨迹CSV只包含本次进程新写的行，sweep编号接续；
`elapsed_seconds`、主`timing`及process-lifetime算术计数属于本次进程。
q/Gamma累计计时随采样器恢复，不能与主timing混作同一计时范围。
结果`resume`记录祖先、继承样本数、诊断分段范围和恢复后的有效选项。

后处理必须把最新完整resumed blocks作为该seed的权威全链记录，
**不能再把祖先blocks与它相加**。未满块若尚未补齐，只在检查点中存在；
`checkpointed`输出不能按`complete`分析。

## PBS操作

CLI写`PREFIX_checkpoint.control.json`，记录PID、hostname、Linux进程start_ticks和PBS_JOBID。
成功保存后写`PREFIX_checkpoint.receipt.json`，记录generation及保存进度。
`scripts/checkpoint_job.py`通过指定SSH端口，在登录节点核对qstat执行节点、
control的job ID，再在计算节点核对PID start_ticks和/proc环境中的PBS_JOBID，
仅向该实际MC进程发送信号并等待新generation成功receipt。
它不qdel、不qsub、不向整个PBS进程组发信号，避免先杀死Python/Bash包装器。
此方式需要登录节点可SSH到执行节点；不支持时明确报错，不用其他方式盲目终止。
超时不等于已保存；不要在未确认receipt时qdel。

既有PBS包装器可能将退出75显示为非零退出，这是有意暂停而非生产完成；
续交时包装器须识别该码，并使用新输出前缀。旧二进制没有信号接口，不能在运行中
通过修改磁盘源码追溯启用，不能直接给旧任务发送这些信号。
已有作业的源码、参数及二进制保持原样，未在本任务中发出任何生产停止信号。

## 验证

本地三种提议local/bond-refresh/bond-guided，q/Gamma非平凡cadence，
无请求及定期保存与修改前原始块逐字节一致；USR1继续、USR2退出、TERM热化中退出，
热化/测量进度和未满块恢复；有限g原始块回归；更改物理/测量选项、截断/损坏文件、
覆盖输出拒绝；原q池文件删除后从内嵌内容恢复。
必要目标环境编译仅做源码/可执行文件/库与help检查，无QMC预检、smoke或新生产链。
准确版本、命令、失败尝试与核验范围以审计记录为准。
