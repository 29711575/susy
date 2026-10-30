# 当前源码与分析状态

整理日期：2026-10-30。本文件描述仓库内容，不维护实时 HPC 队列。

## 当前版本

- 数值核心与 2026-10-28 已验证冻结源逐文件相符；当前身份由
  [SOURCE_MANIFEST.json](../SOURCE_MANIFEST.json) 和 [SHA256SUMS](../SHA256SUMS) 记录。
- 默认 physical matter trial、production-v1 数值策略；guided、非零 q 和测量开关
  仍须在任务中显式记录。已有任务按原程序与参数继续。
- 总能量默认0，旧 matter energy 保留；独立 Γ 更新为可选项。
  多尺度 Gamma±/z 的实验 CLI 没有合入本仓库普通核心。
- 默认新 HPC 提交端口为20015；用户明确指定的端口和队列优先。
  构建检查、必要目标编译和生产数值认证保留，不设置 preflight/smoke 链门禁。

## 已有证据

[10月30日导出验证](../evidence/source-sync-20261030/README.md)记录从当前源码的本地
CLI/确定性验证器构建、3项确定性检查，以及离线统计回归与原始矩复算。
原工作区目标运行时检查单独保存为历史证据，不冒充此次重跑。
这些检查不认证平衡、混合、投影/dt 收敛或整体性能提升。

## 测量入口

[近期测量摘要](../measurements/20261030/README.md)记录 L6/L9 序参量与频率比较。
每组有明确的采集时间、截断和可用 seed；聚合矩不包含逐 sweep 自相关信息。
历史10月10日快照见[测量索引](../measurements/README.md)。

## 历史边界

旧状态日志已移至 [STATUS_BEFORE_20261030.md](history/STATUS_BEFORE_20261030.md)。
其中旧作业是否运行、旧二进制位置、旧默认值和旧验证门禁均不是当前操作要求。
更多归档见[历史索引](history/README.md)。
