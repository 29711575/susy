# ADR-005：矩阵统计顺序与采样证据

状态：recorded。依据：H §10/14、R-SAMPLING、U-CORRECTION；见[来源](../SOURCES.md)。

## 问题

矩阵特征值、均值外积和比值是非线性量；跨动量与相位分母相关。局部错误可能不改变单点中心值，却改变误差和比较。

## 记录的选择

先从共同区块重算相位、完整矩阵和背景，再求Szz、迹、本征值、比值。joint delete-one block jackknife每次同步删除相关统计。不同seed不拼连续时间序列、不平均各自λmax。

中心值、block长度、块数、样本段和误差方法必须同表同图可追溯。最新923004的10/block误差以STATUS为准，20/block不再换标签复用。

## 对“漂移”的解释边界

λmax累计下降不能单独判未热化；固定通道和零模需一起看。两个seed误差重叠亦不能证明稳定。当前10样本基础块不可恢复逐sweep相关性或ESS。

反演点标量相同不许可删除复矩阵虚部或强制平均非等价点。修正配对归属不意味着重新采样或改单点估计器。

## 代码与验证

相关：phase_statistics.h、vbs_window.py、导出器及audit后处理。已有协方差、均值扣除和矩阵先平均测试见NUM-07。增加逐测量诊断是候选增量，不是本次已实现功能。

历史audit不回写；修正结果新目录保存并链接证据。[NUMERICAL_INVARIANTS](../NUMERICAL_INVARIANTS.md)
