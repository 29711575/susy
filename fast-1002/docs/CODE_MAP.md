# 代码地图

生产可执行文件由 `kitaev_ssh/coupled_mc.cpp` 构成，依赖同目录模型/采样代码以及 `inc/` 中的 Majorana Pfaffian 投影实现。

| 路径 | 职责 |
|---|---|
| `kitaev_ssh/honeycomb_lattice.h` | 蜂窝晶格、键方向、边界条件与规范扇区几何 |
| `kitaev_ssh/kitaev_ssh_model.h` | SSH–Kitaev 键耦合与 Majorana 生成元 |
| `kitaev_ssh/phonon_path.h` | 开放虚时间路径、自由声子高斯分布和声子作用量 |
| `kitaev_ssh/phonon_projector.h` | 给定声子路径时构造费米子投影权重和测量量 |
| `kitaev_ssh/phonon_mc.h` | 声子路径提议、接受/拒绝与链状态管理 |
| `kitaev_ssh/coupled_mc.cpp` | 命令行生产程序与结果输出 |
| `kitaev_ssh/vbs_*.h`、`fermion_observables.h` | VBS 与费米子关联测量 |
| `inc/pure_projector_*.h` | Majorana 投影、Pfaffian 权重、Green 函数及稳定化 |
| `kitaev_ssh/*validate.cpp` | 固定输入的数值参考和回归目标 |
| `reports/physics-check-20261005/validation/` | 本次新增的独立解析/多体交叉检查 |

模型参数和 Majorana 记号见[物理与问题说明](../reports/physics-check-20261005/PHYSICS_AND_ISSUES_CN.md)。优化算法的代数细节见 `decisions/ADR-007-bipartite-spectrum.md` 与 `decisions/ADR-008-real-magnitude.md`。
