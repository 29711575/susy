# 历史材料索引

这些材料保留来源版本的模型讨论、数值记录和任务参数，不作为当前运行指南。
当前入口是[仓库 README](../../README.md)和[当前状态](../STATUS.md)。

| 材料 | 用途 |
|---|---|
| [旧交接报告](HANDOFF_REPORT_BEFORE_20261030.md) | 历史极化率曲线、connected 分解与未解决问题 |
| [旧状态日志](STATUS_BEFORE_20261030.md) | 旧源码阶段、任务与验证记录；不是实时任务状态 |
| [P0–P1 自由基线说明](KITAEV_FREE_BASELINE_README_CN.md) | 耦合 MC 实现之前的自由模型构建与验证范围 |
| [更早报告汇编](REPORTS.md) | 原工作区历史报告；部分原 audit 不随仓库发布 |
| [旧 HPC 脚本和计划](legacy-hpc-scripts-before-20261030.zip) | 32个退役 PBS/计划文件，保留原路径、原字节与 SHA256 |

压缩包内 `ARCHIVE_MANIFEST.json` 给出原路径和校验值。需要追溯时可在独立临时目录
解压；其中硬编码的程序、库、队列和旧参数只对应原任务，不应直接重提。
仍被维护工具和回归引用的5个计划文件保留在源码目录，未归档。

历史测量与验证证据分别保留在 `measurements/` 和 `evidence/`，
原本地源码身份清单保留在[来源目录](../provenance/README.md)。
