# ADR-031：路径传播因子存储池与精确实数伴随表示

日期：2026-10-27。状态：implemented / deterministic and target-build verified。154项专项、cache/Fock、11组22次回放及本地/目标构建已验证；分组性能与限制以`audit/1027/path-factor-reuse/REPORT_CN.md`（原工作区本轮报告）为准。本文不预填收益数字。

## 问题与选择

整键候选重建反复生成同尺寸传播矩阵，再转为实数 magnitude 使用的 `D†BD`。选择每链独占的 dense owner 存储池，并由严格二分格谱装配器在生成复矩阵时同时生成精确实数伴随表示。复矩阵仍是 signed、phase、Green、hash 与独立参考计算的权威输入。

`PureFactorMatrix::DenseStoragePool` 保留已分配的槽位。每次完整路径生成调用 `beginBatch()` 重置游标，一个 batch 对已有槽位作一次向前扫描；只有共享所有权计数 `use_count()==1`、即仅池本身持有的槽位才可写。仍被 factor、Snapshot、缓存或外部 proposal 持有的槽位被跳过，不能覆写。没有空闲槽位时新增槽位，避免逐因子从头扫描带来的二次开销。

writer 完整赋值并通过有限值检查后，才把槽位发布为不可变 `PureFactorMatrix`。writer 抛异常或生成非有限矩阵时，不发布失败结果，也不改变已有输出或其他已发布 owner；下一 batch 可再次使用未发布槽位。pool 禁止复制，移动转移池所有权；factor 的复制继续共享不可变 owner，移动继续转移该 owner。

## 精确实数证书

`certified_real` 的发布入口私有，仅严格二分格谱装配器可用。复传播矩阵与实数伴随矩阵从同一组未经截断的谱块装配；证书不是对任意复矩阵作容差判定后附加的标志。实数 magnitude 工作区可直接借用证书，在本次求值期间由 slice 保活，不再重复转换或额外保留被拒候选的 owner。

generic 构造没有证书，继续执行原精确结构检查与转换缓存；任意微小额外复耦合都必须走原 fallback，不能丢弃。generic owner 释放后，缓存可保留已分配的实数缓冲，但不可把旧缓存当成新矩阵的转换结果。

## 数值与生命周期边界

保留 complex B 的系数及 signed zero、原 identity norm、系数 hash 和 immutable identity 语义；不得通过证书重建或替换权威复矩阵。原 Pfaffian phase、signed weight、QR/rank/残差门槛、参考入口、同候选 uncached 恢复与 overflow/异常处理政策保持。非有限谱装配不能发布证书；原 overflow fallback 继续按原路径处理，不能借证书绕过数值失败。

常规候选在接受构型与当前候选间复用存储，steady 状态通常需要两套存活因子。外部保留的 Snapshot、多份 proposal 或缓存会 pin 旧 owner，使池的高水位增长；这是保活语义的结果，不能承诺固定两套或固定内存上限。释放这些 owner 后才允许复用。每链 pool 与工作区不可跨链并发或重入共享。

本决策不改变 Hamiltonian、trial、MC 接受公式、提议 RNG、测量频率或默认测量开关，不影响已提交任务的源码与可执行文件。少分配、少转换本身不能证明整链加速或同误差采样效率提高。

## 代码与所需验证

owner 与 pool 位于 `inc/pure_projector_factor.h`，slice 所有权接入位于 `inc/pure_projector_weight.h`；私有证书装配位于 `kitaev_ssh/bipartite_spectrum.h`；实数借用与 generic cache 生命周期位于 `inc/pure_projector_real.h`；路径池接入位于 `kitaev_ssh/phonon_projector.h`。

专项入口 `kitaev_ssh/path_factor_reuse_validate.cpp` 检查 live owner、复制/移动、多尺寸、writer 失败、Snapshot 与多 proposal 保活，以及复矩阵/hash/identity 和实数证书的一致性。还需配合已有 cache/Fock 验证、冻结 baseline/new 的 L6 A/B/C 与 L9 历史 unrestricted、固定生成 physical 路径回放，覆盖接受、拒绝、测量及独立参考恢复。性能分别记录 prepare、候选、measurement 与完整确定性流水，编译与计时分开，数学线程为1；profile 归因与主计时分开。实际通过情况、内存统计含义和收益边界统一记录于`audit/1027/path-factor-reuse/REPORT_CN.md`（原工作区报告）。
