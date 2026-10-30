# 2026-10-30 过时内容清理

本次基于提交 `ccafc3ff494c2a99f670bf729059c092bbb8806d`，整理入口与退役任务文件。
数值生产核心、历史测量、原始矩和已有验证日志没有改动；不操作 HPC。

- 根 README 只保留当前构建、代码、测量与历史入口。
- HANDOFF_REPORT_CN.md 改为当前研究摘要；旧版完整内容移至 docs/history。
- docs/STATUS.md 不再混放旧运行状态；原日志另存历史目录并标明时间边界。
- kitaev_ssh/README.md 改为耦合 MC 模块入口；早期 P0–P1 说明移至历史目录。
- 32个退役 PBS/计划文件归档于
  [legacy-hpc-scripts-before-20261030.zip](../../docs/history/legacy-hpc-scripts-before-20261030.zip)，
  包含原路径、原字节和 ARCHIVE_MANIFEST.json。5个仍被工具/回归引用的计划留在原位。
- 测量总索引区分10月10日和10月30日；来源目录补充“仅历史记录”的边界。
- 静态/模长说明中的本机旧目录改为仓库相对路径，旧版本限制保留。
- 修正当前规则中的旧状态入口，以及历史报告中未发布压缩包的失效链接。

检查结果见 [summary.json](summary.json)：核心245文件哈希、原测量/验证文件身份、
历史脚本归档字节、计划依赖、当前入口链接和完整文件校验清单。
