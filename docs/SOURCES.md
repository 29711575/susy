# 来源与参考

## 方法论文

- Z.-Y. Han, Z.-Q. Wan, H. Yao, “Pfaffian quantum Monte Carlo: solution to Majorana sign ambiguity and applications,” [arXiv:2408.10311](https://arxiv.org/abs/2408.10311)。用于 Majorana 约定和 Pfaffian 权重方法背景。原文演示有限温度计算；基态投影是本代码的具体使用场景。
- Z. Wu, Z.-Q. Wan, S.-K. Jian, H. Yao, “Emergent spacetime supersymmetry at 2D fractionalized quantum criticality,” [arXiv:2510.07383](https://arxiv.org/abs/2510.07383)。SSH–Kitaev 声子耦合模型的物理背景。
- F. L. Pedrocchi, S. Chesi, D. Loss, “Physical solutions of the Kitaev honeycomb model,” [arXiv:1105.4573](https://arxiv.org/abs/1105.4573)。有限环面的物理宇称约束。
- F. Zschocke, M. Vojta, “Physical states and finite-size effects in Kitaev's honeycomb model,” [arXiv:1502.03063](https://arxiv.org/abs/1502.03063)。物理态选择与有限尺寸响应。
- Y.-R. Shu, S. Yin, D.-X. Yao, “Universal short time quantum critical dynamics of finite size systems,” [arXiv:1705.05931](https://arxiv.org/abs/1705.05931)。有限尺寸虚时间短时标度；用于本模型的2/3关系是条件性推导。

## 代码来源

这份源码快照来自研究工作区 `fast-1002`，其初始核心以 `fast-0929` 为基线，并经两轮计算效率优化。`SOURCE_MANIFEST.sha256` 可用于核对随包源码。PFAPACK 文件及其独立许可证位于 `inc/pfapack/`；顶层许可证见 `LICENSE`。
