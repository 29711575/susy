# ADR-023：默认生产 trial 的物理 matter 奇偶选择

日期：2026-10-19；用户明确要求修正默认 sector。固定规范、Wilson、Hamiltonian
及已提交计算保持原值；新 CLI 默认改为 `--matter-trial physical`。

采用正号自旋 Hamiltonian、标准 `tau=i b c`、`D=bxyz c=+1`。
代码 A→B 存储规范满足 `u_code=-u_phys`。完整局域约束乘积给出

    P_required = (-1)^[Lx+Ly+s*(Lx-s)] * (-1)^Nc * product(u_code).

此处 `P_code=product(i c_A c_B)=(-1)^N_f`；单独调用名为
`pureProjectorPhysicalParity` 的函数只计算这个 Fock parity，不检查几何。
推导、局域8-Majorana算符锚点和此前独立检验保存在
`audit/1012/physical-matter-parity/REPORT_CN.txt`；不能遗漏 u 的符号映射。
一般物理子空间与边界/规范相关，见 Pedrocchi/Chesi/Loss,
https://arxiv.org/abs/1105.4573 。本项目具体符号由上述局域算符锚点决定。

## 选择规则

`physical_matter_trial.h::selectFreeMatterTrial` 在采样开始前只选择一次：

1. 保留原 `freeGroundState` 数学自由真空；已有奇偶相容时原 Phi 完全不变。
2. 若奇偶不相容，交换一个最低能负能模式及其粒子空穴共轭，构造纯
   blocked Gaussian，重验正交/各向同性/Green/奇偶和自由能量增量。
3. 简并时按 n1、n2 字典序选择最低能的规范修饰 Bloch 模式，避免任意
   本征求解器简并基破坏均匀参考态的平移对称性。支持剪切、四种 Wilson
   以及等价局域规范；从规范 winding 得到扭曲动量网格。
4. L6/s0/++/J=(1,1,.5) 选择 (n1,n2)=(0,3)，P 从 -1 修正至 +1；
   自由能量代价为1。L9同设置已相容，原 trial 不变。
5. 零模继续以原1e-10门槛拒绝。若非均匀参考不存在所需最低 Bloch 模式，
   明确拒绝，要求显式构造物理态，不任意挑态、不改规范或 Hamiltonian。

不随声子构型重新选奇偶、重新阻塞或重算目标态。传播保持固定 parity。
实幅值加速仅在修正后重新认证；不满足条件时保留完整复数路径。

## 接口、兼容与证据范围

生产 `kitaev_ssh_mc` 和无 MC 的 `kitaev_ssh_free` 默认 physical。
`--matter-trial unrestricted` 是显式历史/数值控制，记录其真实相容性；
不得将不相容控制解释为原自旋模型结果。原 `freeGroundState(model)` 保留
无约束数学含义，供独立 Fock、Bloch 与数值回归使用，不代表生产默认。
演示程序同样采用 physical。静态极限历史优化工具的 `--parity trial`
仍表示其原先无约束自由真空的固定奇偶，不因本次改变历史参照。

完整结果新增 `matter_trial_sector`，保存 policy、required/unrestricted/actual
parity、blocked、能量代价、最低模式简并数及选择的动量。
生产启动即写 `PREFIX_matter_trial.json`，输出保护包含该文件。
`matter_trial` 名称区分物理奇偶自由态与无约束控制。

这项选择保证固定规范下奇偶投影不为零，不是显式规范轨道求和，也不
保证该 Wilson 全局最低、不保证与相互作用基态具有非零重叠或已收敛。
有限 tau 仍依赖 trial；本次不能恢复或重解释旧 L6 生产数据。
未来冻结新源码并在目标平台构建时使用新默认；已有冻结源码/二进制及
运行任务不被修改。

## 验证

专用 `physical_matter_trial_validate.cpp`：独立 Clifford 重排与几何条件、
四Wilson/剪切/局域规范、最低相容能量、L6平移对称、L9原 Phi保留、
零模门槛；16维 Fock 独立偶子空间能量及固定 blocked trial 的有符号
路径权重/中心 Green，并核对 reference/production 两种执行路径。
CMake 与 nodeps 均接入，`run_pfqmc_correctness.py --group physical` 可单独执行。
本轮实际执行结果、CLI及兼容检查在 `audit/1019/physical-matter-trial-default/`，
以该目录报告为准；本 ADR 不把尚未完成的检查写成通过。
