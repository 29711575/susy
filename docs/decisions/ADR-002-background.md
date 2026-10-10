# ADR-002：均值背景与整个零模不是同一个扣除

状态：recorded。依据：H §9、C0/C1、U-CORRECTION；见[来源](../SOURCES.md)。

## 问题

SSH耦合下可以有均匀声子均值；非网格Fourier窗口会把它贡献到非零Q。只隐藏q=0并不能完成背景处理。

## 记录的选择

主结果采用OBS-03的平移平均connected；同时保存raw和projected对照。connected只扣均匀系综均值的外积；projected逐构型去掉整个零模，包括涨落。

q=0及等价端点不进非零动量图或峰搜索，但原始零模一阶/二阶统计必须保留。自由声子涨落参考不作为可任意减去的“trivial均值背景”。

## 不能采用的替代

不能把Sprj直接命名connected；不能删除零模原始数据后声称以后可恢复背景；不能在最终矩阵扣一次均值而沿用raw误差。

## 代码与验证

相关：vbs_window.py的不同connected辅助函数、导出器、生产区块和联合jackknife。

应检查白噪声、均匀均值、非零允许点以及off-grid的不同定义。每次删块重算零模均值和共同相位分母。已有测试函数及证据状态见NUM-07。

这是一条定义约束，不要求数据本身接近任何预先设定的背景数值。[OBSERVABLES](../OBSERVABLES.md)
