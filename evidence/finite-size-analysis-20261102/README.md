# 不同 L 分析补充：复核记录

日期：2026-11-02，Asia/Shanghai。基线：远端 main
`e5fc359cdb15b09417446e7e56e543c00691f114`，开始时工作树干净。
用户要求简单检查仓库中不同 L 的分析，不全面则补充。

原文已有 L6/L9 physical trial、动量邻近集合及参考斜率说明，但缺少统一的
L6–L9 裸谱表、奇偶与公度效应区分、短时初态/背景条件及比较方法。
本次增加[分析说明](../../docs/FINITE_SIZE_AND_SHORT_TIME.md)，并链接当前入口。
不修改物理/采样核心或历史数值，不启动 QMC、优化或 HPC 作业。

## 实际复算

从仓库根目录运行，依赖 Python 3 和 NumPy：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 validation/finite_size_reference.py
```

[脚本](../../validation/finite_size_reference.py)直接构造 T，对 L=6/7/8/9 比较实空间
奇异值与独立 Bloch 谱，核对 det 符号并按既定 Fock 约定计算物理 trial 修正。
输出保存为 [bare-spectrum.csv](bare-spectrum.csv)。这些是 X=0 的确定性参考，
不是相互作用 QMC 的 gap、临界耦合或平均符号预测。

[检查摘要](checks.json)记录裸谱与已有参考的核对、修改范围、相对链接及校验清单检查。
入口与新分析的相对文件链接均检查存在；所有原核心、测量数值、图和历史证据保持原字节。
SHA256SUMS 与 SOURCE_MANIFEST.json 的当前文件身份同步更新；旧导出与构建记录仍保留
原日期和范围，本次不冒充重新构建。文献支持与本次推导的边界写在分析说明中。
