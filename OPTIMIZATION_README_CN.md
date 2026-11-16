# 整键路径提议与多动量更新

本文只说明功能及适用范围。模型与抽样流程见[QMC审阅说明](docs/QMC_REVIEW_GUIDE_CN.md)，详细算法与历史计时见[2026年10月10日报告](docs/SAMPLING_OPTIMIZATION_RESULTS_20261010.md)。

## 整键路径提议

`--proposal bond-guided --bond-guide bare --bond-beta 1`为一条键的全部虚时间坐标生成候选路径。它使用初始化时确定的参考力改善候选分布；接受概率保留精确的提议密度修正，因此参考力不是实际外场，也不改变目标Hamiltonian。

`bond-beta`取值在(0,1]，表示提议幅度，不是逆温度。参考力估计不准确会影响效率，但不能因此省略接受率中的修正。具体实现见[phonon_guided_proposal.h](kitaev_ssh/phonon_guided_proposal.h)与[phonon_mc.h](kitaev_ssh/phonon_mc.h)。

## 多个固定动量

`--q-mode-pool FILE`可替代单个`--q-mode-n1/--q-mode-n2`。文件列出整数动量索引；每次更新均匀随机选择一项，不根据当前测量峰调整。每项已包括q与−q对应的实余弦、正弦模式，不应重复列入共轭动量。零动量不在该列表中。

`--q-mode-bonds`指定参与的键类型。动量整数必须按实际晶格和边界重新选择；示例列表不能直接套用到不同尺寸。多个动量或多个键同时变化，并不保证比单模式更高效。

## 验证与性能限制

[guided_proposal_validate.cpp](kitaev_ssh/guided_proposal_validate.cpp)检查提议概率和接受率；[q_mode_validate.cpp](kitaev_ssh/q_mode_validate.cpp)检查动量模式更新。历史整合报告记录过本机验证，但部分原始日志路径未随仓库提供；可访问的执行记录应从[报告索引](docs/REPORT_INDEX_CN.md)查找。

评价效率应比较相同物理参数、相同误差下的计算时间，并检查相关时间和从初态达到平衡所需的长度。接受率提高或小系统变快，不证明所有观测量都会加速，也不能据此缩短物理投影长度或按比例减少舍弃样本数。

[完整旧说明](docs/history/OPTIMIZATION_BEFORE_20261116.md)保留命令示例、旧验证范围和来源。当前构建入口统一见[首页](README.md)，不再在本页重复列出演示计算预算。
