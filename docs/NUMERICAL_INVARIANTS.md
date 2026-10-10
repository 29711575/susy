# 数值与统计不可混淆的约束

状态：继承已有物理/数值原则；下列建议测试不等于本轮执行结果。
来源：[H](SOURCES.md#h) 第10、14、16节，[C0](SOURCES.md#c0)，[R-SAMPLING](SOURCES.md#r-sampling) 及 [U-CORRECTION](SOURCES.md#u-correction)。


## 证据边界

- **论文事实**：权重、相位和统计原则的外部依据。
- **数学推导**：重加权、背景、矩阵统计和区块估计的公式。
- **源码现状**：失败恢复、缓存和测试函数是查看到的接口事实。
- **报告结果**：历史测试或运行报告按来源标签引用，不称为本轮重跑。
- **未验证计划**：驱动新增验证清单是要求，不是已完成测试。


## NUM-01 权重与相位

对路径 X，w(X)=|w(X)|exp(iθ)。按 |w| 采样后，

$$
\langle O\rangle=\frac{\langle e^{i\theta}O\rangle_{|w|}}
{\langle e^{i\theta}\rangle_{|w|}}.
$$

一阶、二阶、不同动量和矩阵元使用同一批样本和相位分母。相位接近 1 不授权省略重加权；`defer-phase` 是推迟认证，不是设相位为 1。

U†U 的总范数为正，不说明独立 bra/ket 场构型上的 Pfaffian 权重逐构型为正。不能用正平方根、平方权重或镜像声子路径来“修复符号”。[H, §10.2、§14；C0, phonon_projector.h 的路径权重]

## NUM-02 无法认证不是正常拒绝

现有策略名 `uncached_reference_once_then_abort`：首次失败，对同一候选做一次无缓存参考重算；保持原接受随机数与阈值。成功则继续原判据，仍失败则保存诊断并中止。拒绝不污染已接受状态；初始化重抽、提议失败和测量失败分开计数。

`running/complete/aborted` 区分保存。被中断且状态仍 running 的文件，不可按完整链静默分析。故障诊断不等于完整 RNG 可续算检查点。[H, §14.3；C0, phonon_mc.h]

## NUM-03 缓存只用于同一个物理对象

模型、几何、规范、初态、离散选项和声子节点共同决定传播对象。当前固定模型的不可变捕获对常量计算合法；若扩展驱动，节点耦合与协议身份也应纳入正确性设计，或由等价不可变结构确保不混用。

不规定必须使用某种 hash/key 实现，但必须证明缓存、参考重算和局部提议看到的是同一协议。生成元不能依赖调用次数推断物理时间。驱动分片不对易，顺序变化需显式验证。[C0；DRV 草案]

采用区间中点 `g_{j+1/2}` 时，`X_j` 的左右半步属于不同的 `(node, adjacent interval, protocol)` 对象；不能按固定 g 将一个半步复制、平方或直接套用旧谱/Cayley 缓存。第一版可以保留总 `PhononPath`、唯一中心和两外端 trial，但 factor identity 必须保留相邻区间。协议身份至少包括 `g_i,g_f,R,tau_f,dt`、中点规则、分支/中心 cut、模型几何/规范、初态和声子路径身份。

## NUM-04 先矩阵后非线性量

窗口加权、相位重加权和 connected 背景均先在矩阵层面完成。λmax、比值和拟合随后计算。有限样本的 λmax 可有非线性估计偏差；累计下降不能单独证明未热化。

可报告 Hermitian 化的谱，但必须保存/报告非 Hermitian 残差。不得截去负本征值或负关联比以制造“正确形状”。[H, §10；U-CORRECTION；DISCUSS]

## NUM-05 联合区块与误差口径

每次删除同一基础区块，重算相位分母、全部相关矩阵/均值、窗口、背景和非线性量。跨动量协方差不能忽略。不同 seed 不拼成一段连续时间序列，也不跨 seed 合块。

主值默认 full-sample plug-in；删块副本均值、偏差修正值另名输出。每张表/图记录 samples/block、block_count、样本选择、归一化和误差方法。10/20/40 samples/block 是不同误差估计，不能换标签。

只有 10 样本区块时，可以研究更粗重分块和区块累计，不能重建逐 sweep 自相关时间或有效样本数。块数达到软件最低门槛不是独立性的证明。[H, §10；U-CORRECTION]

## NUM-06 建议保留的最小元数据

现有字段优先复用；缺少的字段才考虑增量，不要求为本文件立即改格式。

运行身份应能追到源码 commit/dirty patch、构建/库版本、输入文件、模型/几何/扇区、seed、实际更新模式、试探态、Δτ、总轮廓长度、采样计划、区块配置和完成状态。参数不能从 CLI 默认值反推。

驱动扩展还应记录控制参数名称、起终点、请求和实际速率、单支时间、总长度、节点/中点取值规则及测量位置。原始矩阵/相位统计和文件校验值应允许后处理追溯。[R-SAMPLING；DISCUSS]

## NUM-07 证据对应的测试，不重复发明测试名

C0 中已看到的 Python 测试包括：

| 约束 | 已查看的测试文件/函数 |
|---|---|
| 完整矩阵先平均 | `test_vbs_window.py::test_matrix_is_weighted_before_largest_eigenvalue` |
| 窗口与空间平移 | `test_vbs_window.py::test_translation_invariant_known_covariance_and_all_shifts` |
| 联合动量协方差 | `test_vbs_window.py::test_joint_momentum_covariance_survives_window_jackknife` |
| 均值扣除重采样 | `test_vbs_window.py::test_connected_window_jackknife_recomputes_first_moments` |
| 白噪声与 projected | `test_vbs_window.py::test_white_noise_covariance_and_projected_zero_mode_factor` |
| 静态实空间/动量核 | `test_static_vbs_response.py::test_allowed_k_reference_matches_real_space_elementwise` |

这些是查看到的函数，不是声称本轮重跑通过。导出增量与后续合并的 40 项通过记录来源见 C1/U-MERGE；生产 HPC 运行证据见 R-SAMPLING。

驱动新增验证仍是要求：常量协议回归；非对易分片的小系统有序乘积与相位/Green 核对；单支/双支、节点时间和终点耦合；同协议参考重算；时间步长改变时保持实际协议一致。先检查完整仓库已有测试，再补缺项；不因上面列了要求就宣称功能通过。

本轮 validation-only 中点回归已覆盖常量协议、9 组小系统中点有序乘积、相位和中心 Green；它没有覆盖双支采样、缓存/恢复或生产测量，因此上述不变量清单仍未全部完成。[R-MIDPOINT-REGRESSION]

## NUM-08 新结果更新规则

新运行/后处理写新 audit；不改旧数据以消除勘误。状态变化更新 STATUS 并链接证据；设计改变更新相关 ADR。文档中的目标不是权限：需要实际运行、改采样或提交 HPC 时，以用户本轮授权为准。
