# fast-1006 采样优化：实现、验证与性能边界

来源说明：以下为用户提供的优化包报告，保留其原始结论与运行口径。
quick1010 整合另有原生 MKL 验证，见
`../../fast-1006/audit/1010/quick1010-guided-integration/REPORT.txt`。
本次整合未重跑下文的长 benchmark，不将其性能数字记作本次实测。

日期：2026-10-10。基线为用户提供的 `fast-1006-core-current-20261010.zip`，
SHA256 `0a200aea9559becbafd23874166bc1e6da0d27ea25b3b0514aec5328c43732b7`。
原压缩包及其解压参考源码未修改。本报告记录本轮实际运行，不把此前审计当成本轮测试。

## 结论

已实现一版可选启用的采样扩展。新整键提议在小系统上明显改善均匀声子背景的慢相关，
一个 L2 对照中固定动量极化率的每有效样本 CPU 成本也下降。没有证明所有系统或所有
观测量都同样加速，更没有认证 L9 临界计算的热化 sweep 数。新功能默认关闭；旧的
local/bond-refresh 行为保留。

重点区分：平衡态自相关时间、从初态达到平衡的热化时间、每 sweep CPU 时间、物理
投影长度 Theta 是四件不同的事。本轮短链性能估计不能确定最小热化截断。

## 1. 改了什么，以及为什么不改变物理分布

### 1.1 固定参考均值的整键 pCN 提议

旧 `bond-refresh` 从自由声子分布重新抽取一条键的全部虚时间路径。相互作用若把均匀
位移推离这个参考均值，提议可能反复落在低权重位置。新 `bond-guided` 使用更接近
耦合体系均值的固定 Gaussian 参考，但不改变被采样的权重、实际外场或试探态。

目标幅值测度仍是

    p(X) proportional to W_ph(X) |W_f(X)|,

Pfaffian 相位仍按原逻辑重加权。选固定参考力 f_b，定义

    p_ref(X) proportional to W_ph(X) exp(sum_b f_b I_b(X)),
    I_b=sum_j w_j X_jb,
    w_0=w_M=tanh(omega*dt/2)/omega, w_internal=2*w_0.

这里 w 是原 Mehler 核线性源的精确离散权重，不另行改用 dt 梯形权重。
Gaussian 参考的均值为

    mu_ref(t,b)=mu_free(t,b)+(f_b/k)*[1-(exp(-omega*t)+exp(-omega*(Theta-t)))/2],

协方差为 C_jl=exp(-omega*|t_j-t_l|)/(2*m*omega)。提议为

    X'_b = mu_ref_b + sqrt(1-beta_p^2)*(X_b-mu_ref_b) + beta_p*xi,
    xi ~ N(0,C), 0<beta_p<=1.

它关于 p_ref 可逆。因此精确 Metropolis–Hastings 比必须是

    log r = log|W_f(X')|-log|W_f(X)| - f_b*(I_b(X')-I_b(X)).

最后一项不能漏；漏掉会改变目标分布。新代码实现了该项，并另外计算完整前后向 Gaussian
提议密度供独立交叉检查。参考均值在运行前固定；估得不准影响效率，不是物理近似。

`--bond-guide bare` 从同一有限图的裸费米子键关联估计 f_b=-g*<u_b i c_i c_j>，
按键类型平均，初始化一次后冻结。用户也可固定手工参考力。它不是实际的 VBS 共轭场。
`--bond-beta` 是提议步长，不是逆温度。

算法属于保 Gaussian 参考的 pCN 思路；参见 Cotter, Roberts, Stuart, White,
“MCMC Methods for Functions: Modifying Old Algorithms to Make Them Faster”,
arXiv:1202.0709，Statistical Science 28, 424–446 (2013)。以上偏移参考的接受率
由本程序的开放 Mehler 路径测度推导，未照搬有限温周期路径的接受率。

每 sweep 的基础提议仍是一条键一次，不靠增加提议次数来人为减少 sweep 数。
缓存均值和积分权重只在初始化时构造。

### 1.2 固定动量池及多键模式（可选，未证明普遍加速）

旧的单 q、z 键 cos/sin 模式保留。新增固定动量池，每次计划的集体更新均匀随机选一项；
每个 q 项可包括 x、y、z 中指定的键类型及对应实余弦/正弦模式。

- q 选择不依赖当前构型或当前统计峰，不在线追着峰改 proposal。
- 正交化检查、q/-q 重复项检查、Gamma 拒绝、shear 动量均纳入。
- q 更新仍以原自由 Gaussian 为参考，其接受率仍是 matter 幅值比。
- 基础 guided 更新与 q 更新的参考不同，二者不能混用接受率修正。
- 物理测量仍为原 all-q 结构因子和极化率矩，未限制在池内。
- 池模式的 q-diagnostic 仅记录第一项的 z 分量；每池项累计提议/接受数包含
  warmup+measurement，并在 JSON 中标明。

动量池覆盖更多候选模式，但在相同集体提议预算下也会稀释某个固定 q 的更新频率。
本轮测试确实观察到指定 q 极化率比 guided+单 q 更慢，所以不要默认将池视为更快。

## 2. 实际完成的正确性检查

测试环境：单线程 OpenBLAS、Eigen、GCC 14、-O1；没有 MKL，使用隔离的审计兼容层。
该兼容层不在生产源码中，也不是交付的推荐生产后端。原 MKL 构建/生产机器必须再跑。

本轮新验证程序 `guided_proposal_validate` 的最终输出：

    PASS
    worst deterministic residual = 7.780442956573097e-13
    largest finite-sample difference = 0.014818438728930697

第二个数来自有限 Monte Carlo 样本的线性模型热浴检验，不是详细平衡误差。
覆盖开放端点和真实源、参考均值 Cw、完整 Gaussian 前后向密度、2-Majorana Fock
权重、不对易 shear 蜂窝构型、接受/拒绝的完整恢复、同候选同随机数重试、持续故障中止、
多键模式正交性、池输入检查、g=0 接受率和线性模型精确热浴。

已有 q_mode_validate、mc_validate、numerical_failure_validate、
vbs_susceptibility_validate 均在修改后的代码上重新运行通过。后者源项有限差分
误差 8.5546e-11，自由响应误差 6.6613e-16，Fourier 检查误差 7.0763e-17。

另外编译原始与修改后生产入口，分别运行 local/bond-refresh、q 开/关四组相同 seed
对照。四组原始 blocks CSV 字节一致，最终路径、提议/接受计数和相位也一致。
这说明新增功能关闭时，本轮所测旧路径保持不变；不是对所有平台逐 bit 一致的保证。

测试开发中一次扩展检查误把未 materialize 的参考相位当已测量相位比较，产生过一条
测试断言失败。修正独立参考为显式相位测量后重新通过；未为通过该测试修改生产权重。
证据包保留开发日志，最终判据见 guided_validation_release.log/exit。

### 独立的不对易相互作用检验

模型为 4 个 Majorana 加 1 个量子声子，不是蜂窝 VBS 生产模型。总投影 4，omega=.7，
k=1，g=.8，dt=.2。故意使用不准确的固定参考力 .4，beta_p=.75，以检验接受率修正。
QMC warmup1000、采样60000、seed9931；精确侧用独立 Fock×24 级谐振子空间的同 dt
转移矩阵和固定试探态源项差分。

| 量 | 同 dt 独立转移矩阵 | 新 guided QMC（1SE） |
|---|---:|---:|
| <X> | 0 | -0.005880 ± 0.008147 |
| <X²> | 0.423331795 | 0.426238681 ± 0.006230147 |
| chi | 0.990538181 | 0.987590083 ± 0.016166988 |

三者均相符。该链平均相位 1，基础接受率 .7894，无数值故障。这不代表大系统无符号问题。

L9 另做了 CLI/尺寸烟雾测试：M8、dt.1、总投影仅 .8、两个 sweep，新 guided+xyz池
正常完成且无数值故障。它**不是**生产 tau8 的混合或收敛测试。

## 3. 性能：同一目标权重、同一每 sweep 提议次数

蜂窝 benchmark 均为 J=(1,1,.5)、lambda2.4、omega.3、k1、Theta8、M80、dt.1。
基线是已有的 `bond-refresh+单q z键集体更新`，并非只有局域更新的弱基线。
改进组只把整键刷新替换为 guided bare beta_p1，其余单 q 集体更新保持不变。

L2：每组4条独立链（seed7411,7427,7451,7481），每链 warmup500、采样5000。
每 sweep 两组都是12次整键+1次单q更新；q=(1/2,1/2)。

| 指标 | 原 bond-refresh+单q | 新 guided+同一单q |
|---|---:|---:|
| 整键接受率 | 0.3609 | 0.5940 |
| CPU 秒/sweep | 0.01732 | 0.01708 |
| 均匀 Xz 的 tau_int，中位数/sweep | 23.5 | 5.23 |
| 均匀 Xz 的 CPU 秒/有效样本，中位数 | 0.826 | 0.184 |
| 固定 q m² 的 tau_int，中位数/sweep | 0.883 | 0.680 |
| 固定 q 连通 chi 的 tau_int，中位数/sweep | 1.04 | 0.745 |
| 固定 q 连通 chi 的 CPU 秒/有效样本，中位数 | 0.0375 | 0.0251 |

均匀 Xz 有效样本成本的点估计下降约4.5倍，固定 q 极化率约1.5倍；单 sweep 成本没有
显著改变。有符号 D 的相关性未改善。不能将较大的 Xz 改善系数当成 chi 的改善系数。

均值对照：m² 原 .0457328±.00111、新 .0458582±.00104；chi 原 .929543±.0259、新
.924741±.0191。误差为多链保守摘要（链内IAT误差合成与链间SE取较大者），不是已认证
长链置信区间。原始 traces 和重分块结果在证据包内。

L3：每组两条链、warmup300、采样1200，每 sweep27+1次提议。
接受率 .3610→.7957，CPU/sweep .2802→.2797。均匀 x/y 慢相关明显改善，但旧链长度
不足以可靠估算慢模式 IAT，**不报告其确定加速倍数**。固定 q chi 的中位成本 .453→.555
秒/有效样本，未显示改善。样本较短、仅两seed，不能作显著变慢或不变的结论。

额外 L2 guided+动量池xyz pilot 两条链：指定 q chi 的成本 .0453秒/有效样本，反而
高于 guided+单q 的 .0251。此处同时改变了动量分配、键组合和 q 步长，不能单独归因
于某一个因素，只说明这套池配置不是无条件更高效。

共完成14条蜂窝 benchmark，全部正常退出，无数值故障。实际耗时与命令记录在
logs/benchmark_runs.json。CPU/有效样本以 2*tau_int*t_sweep 估计；相位加权连通 chi
使用联合矩的 delta-method 影响函数估计 IAT，不是只拿 ChiRaw 的 IAT 代替。
提供16、32、64、128、256块长检查。短链IAT可能漏掉弱长尾，上述倍数是点估计，
不是L9生产加速承诺，更不是最小warmup证明。

## 4. 如何使用及如何在 L9 判断是否值得换

首先只更改基础提议：

```bash
--proposal bond-guided --bond-guide bare --bond-beta 1
```

其余保持原命令，包括物理 tau、dt、Hamiltonian、Wilson扇区、trial、原固定q与测量。
完整可执行 CLI 示例在源码根目录 OPTIMIZATION_README_CN.md。

新参考会影响效率，若接受率仍偏低可单独对照 beta_p=.5,.8,1，但不要根据测量中的
即时峰值不停调参数。任何选定参数在正式测量前固定并重新检查热化。

第二步才独立试验动量池，L9/shear0示例已含共轭：

```text
5 4
3 5
4 6
```

```bash
--q-mode-every 1 --q-mode-pool examples/qpool_L9_s0.txt --q-mode-bonds xyz --q-mode-beta 0.5
```

不要同时传原 q-mode-n1/n2；不同几何必须重新定义池。池不会删除池外物理测量。

L9 建议先选择一个弱耦合点和一个响应增强点，各比较原提议与guided提议，至少4条
独立链、不同声子初态（但保持同一个trial）和一致每sweep预算。先采用充足初始预算，
根据分段均值、增加弃样长度、跨初态一致性和重分块判断，而不是预先假定 warmup=500
或800足够。观察均匀 Xx/y/z、有符号 D、m²、ChiRaw、ChiD、ChiIconj、phase，成本应
以实际关心的连通极化率统计误差衡量。

当前用户说明中的前800 sweep 后处理弃样并未认证热化，因此不能把800除以4.5后
宣布新热化长度。也不能减小 tau 来节省 sweeps；物理投影误差不会被好更新消除。
新算法不解决符号问题、不改变物理自旋投影范围、也未实现 HMC 或副本交换。

## 5. 文件及可复现性

- 完整可选优化源码：fast1006_guided_sampling_source_20261010.zip。
- 增量补丁：fast1006_guided_sampling.patch，以原快照为基线，可用 patch -p1。
- 测试/benchmark证据包：fast1006_guided_sampling_tests_20261010.zip，含补丁、独立
  测试程序、构建适配器、命令、原始短链CSV、结果和分析脚本，不含二进制及第三方库。
- 新SHA256SUMS是修改后目录的清单；ORIGINAL_SHA256SUMS和SNAPSHOT.original.json
  仅记录原版出处，不要拿它们验证修改后的目录。
- 原科学估计器、Hamiltonian、PhononPath、Pfaffian/Green内核和旧后处理保持逐字节不变。

新目标 `kitaev_ssh_guided_validate` 已加入原CMake测试区。用户原MKL环境需自行重新
编译并运行 README 所列验证；本轮没有访问生产机器或替换旧冻结二进制。
