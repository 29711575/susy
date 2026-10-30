# 历史来源记录

本目录保存原工作区文档、源码身份和旧校验清单，用于追溯历史数据。
文件中的“当前”、默认端口、build/audit 路径及 AGENTS 规则属于原记录日期，
**不选择当前构建或当前任务设置**。

当前操作规则见[根 AGENTS.md](../../AGENTS.md)，源码身份见
[根 SOURCE_MANIFEST.json](../../SOURCE_MANIFEST.json)，文件校验见
[根 SHA256SUMS](../../SHA256SUMS)。

- `SOURCE_MANIFEST.*.json`、`SHA256SUMS.*`：对应日期的导出身份，不能拿来校验当前树。
- `PRODUCTION_SOURCE.local-20261023.json`：当时本地冻结生产版本，路径未随仓库迁移。
- `quick1010-local-*`：原工作区文字和清单的来源副本，不是当前配置。
- `STRICT_STATIC_README.local-20261010.txt`：当时的静态工具说明。

这些文件为保证历史来源可核对而保留；当前入口文档已移除对它们的操作依赖。
