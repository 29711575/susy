# quick1010 采样优化使用说明

2026-10-10 已在 quick1010 整合本扩展，并在本机原生 MKL 后端通过新提议
正确性验证及四组旧模式兼容性对照；详细结果见
`../fast-1006/audit/1010/quick1010-guided-integration/REPORT.txt`。
包内 OpenBLAS 性能结果是外部提供的历史证据，本次没有重跑其长 benchmark。
未来普通整键生产链显式使用 guided bare/beta1，单 q 设置沿用原规定；
动量池按目标另选。CLI 省略 proposal 仍是 local，生产配置需显式记录。

本版本基于用户上传的 `fast-1006-core-current-20261010.zip`，SHA256 为
`0a200aea9559becbafd23874166bc1e6da0d27ea25b3b0514aec5328c43732b7`。
这是可选择启用的采样扩展，不是新的物理 Hamiltonian，也不是已认证的
L9 生产版本。原 `local` / `bond-refresh` 默认行为保留。旧作业和冻结二进制
没有被更改。新增优化不会替代投影长度、dt、物理态或符号问题的检查。

## 1. 先只换整键提议，便于做有意义的对照

在原运行命令上，将

    --proposal bond-refresh

替换为

    --proposal bond-guided --bond-guide bare --bond-beta 1

其他参数全部不变，包括 tau、dt、trial、测量频率和原单 q 集体更新。
`bond-beta` 是 pCN 提议步长，不是逆温度；合法范围是 (0,1]。
`bare` 从同一有限图的裸费米子基态键关联计算固定参考力，初始化一次后冻结。
它是用于生成提议的均值近似，接受率包含精确修正，不把目标分布近似为自由态。

参考力可显式设置，例如

    --proposal bond-guided --bond-guide zero --bond-guide-x 0.9 --bond-guide-y 0.9 --bond-guide-z 0.5 --bond-beta 0.8

这些数只是语法示例，不是通用最优参数。不要把 `bond-guide-*` 与
`vbs-field-h` 混淆：前者不加实际外场、不移动 trial、不改变物理测量。
新结果 JSON 会写出实际参考力、beta 和这些约束。

## 2. 再独立测试动量池（非强制、未证明普遍加速）

原单 q 更新仍可使用，例如 L9 普通 PBC 的

    --q-mode-every 1 --q-mode-n1 5 --q-mode-n2 4 --q-mode-beta 1

如果要覆盖多个候选软模，可替换为

    --q-mode-every 1 --q-mode-pool examples/qpool_L9_s0.txt --q-mode-bonds xyz --q-mode-beta 0.5

不能同时传 pool 和 q-mode-n1/n2。池文件只允许整数 n1 n2，支持 # 注释。
每次计划的 q 更新从池中均匀随机选一项，不从当前峰动态选择；每项已经
包含 q/-q 的实余弦/正弦子空间，所以不可再添加共轭重复项。
Gamma 被拒绝。`q-mode-bonds` 也支持 x、y、z、xy 等不同、不重复的键类型。

示例池仅适用于 Lx=Ly=9、shear=0；不同几何必须重选整数索引。
所有键类型同时提议和分散到多个 q 可能降低某一个已充分更新的 q 的效率。
不要把动量池当成无条件比单 q 更快的默认配置。

`q-diagnostic-every` 在池模式下仅记录第一项的 z 键分量。
科学分析继续使用原完整 all-q 结构因子/极化率矩。JSON 新增每个池项的累计
接受/提议计数，范围为 warmup + measurement；原分阶段 q/base 计数不变。

## 3. 本机示例：这是性能/正确性对照，不是收敛认证

以下 warmup/sweeps 为演示预算，不保证任何参数已经平衡或达到基态。
使用新输出目录，不覆盖旧链。

```bash
./build/kitaev_ssh_mc \
  --Lx 9 --Ly 9 --shear 0 --Jx 1 --Jy 1 --Jz 0.5 \
  --g 1.5491933384829668 --omega-ph 0.3 --mass 11.11111111111111 \
  --intervals 80 --dt 0.1 --warmup 500 --sweeps 2000 \
  --measure-every 1 --block-samples 10 --seed 7101 \
  --proposal bond-guided --bond-guide bare --bond-beta 1 \
  --q-mode-every 1 --q-mode-n1 5 --q-mode-n2 4 --q-mode-beta 1 \
  --vbs 1 --vbs-mode all --vbs-susceptibility 1 --measure-energy 0 \
  --acceptance-every 10 --q-diagnostic-every 1 --output test_guided_seed7101
```

基线命令只改为 `--proposal bond-refresh` 并移除 bond-guide/bond-beta，换新
输出前缀。应有多条独立链，并比较不同起始状态（若原主程序不支持外部初态，
不要仅通过修改 trial-shift 冒充同一个有限投影问题的不同初态）。
`--measure-energy 0` 仅是避免不需要的中心 Green 重建，可根据实际观测需求改为1。

## 4. 构建与验证

沿用原 CMake、MKL、Eigen、Boost/PFAPACK 配置；新增目标：

```bash
cmake --build build --target kitaev_ssh_mc kitaev_ssh_guided_validate kitaev_ssh_q_mode_validate kitaev_ssh_mc_validate kitaev_ssh_numerical_failure_validate
ctest --test-dir build --output-on-failure -R 'kitaev_ssh_(guided_proposals|q_mode|p3|numerical_failure)$'
```

上述测试目标需按原项目开启测试选项。不要为提速打开 fast-math、降低数值
认证阈值、关闭相位重加权或把无法计算的候选当作拒绝。
上传包原先在隔离 OpenBLAS 适配环境测试；本次整合另外完成原生 MKL
专项验证，实际重跑范围与命令以整合审计记录为准。
免 CMake 构建脚本现也支持 `kitaev_ssh_guided_validate`。

## 5. 如何决定能否少跑 sweep

比较的目标是相同物理参数、相同统计误差下的 CPU 时间，而不是接受率或
总 sweep 数。建议监测各键的均匀均值、固定 q 的有符号 D、|D|²、ChiRaw、
ChiD、ChiIconj 和 phase。连接极化率的误差仍须从相同块中的联合矩计算。

    每有效样本的时间约为 2 * tau_int(O) * t_sweep

改为新提议后，不要预先把旧的 800-sweep 截断直接除以某个小系统加速因子。
先检查热化截断翻倍、时间序列和跨初态一致性。短链 IAT 不是最低热化长度
的证明；均匀位移变快也不意味着 VBS 幅度/模式转换都等比例变快。

更不能为了减少总时间把 tau 从已收敛值调小：Monte Carlo sweep 和物理
投影长度是两种不同的“时间”。本优化保持后者不变。

## 改动位置

- `kitaev_ssh/phonon_guided_proposal.h`：固定参考 Gaussian、AR 提议密度、精确修正。
- `kitaev_ssh/phonon_mc.h`：guided bond 接受/恢复；固定 q 池调度；旧路径默认不变。
- `kitaev_ssh/phonon_q_mode.h`：多键模式、池读取、共轭去重检查。
- `kitaev_ssh/coupled_mc.cpp`：新增 CLI、一次性参考力计算、元数据。
- `kitaev_ssh/guided_proposal_validate.cpp`、`CMakeLists.txt`：新验证目标。
- `docs/decisions/ADR-017-guided-bond-and-q-pool.md`：详细平衡推导和范围。

Pfaffian/Green 内核、Hamiltonian、声子路径权重、易感率估计器和原后处理均未修改。
