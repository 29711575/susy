# ADR-018：固定均匀背景的静态软模筛选

2026-10-10 更正：以下记录依据用户随后说明有误的旧 note；工具及独立验证
仍保留为静态参考，不代表当前 note 的测量任务。当前实现见 ADR-019。

状态：implemented / locally validated，2026-10-10。用户要求审查
`new/note.docx` 的建议，并在 quick1010 基础上实施合理部分。
docx 指向旧工作区，没有外部超链接；本次采用用户指定的 quick1010。

## 问题与选择

裸 Dirac 点或均匀自洽背景给出的 Q=±2K 是候选；有限系统的最软模式
不能由最近网格点直接决定。单独检查 z 键还可能漏掉三键混合通道。
沿用并修正 `kitaev_ssh/static_vbs_response.py`，在同一个 BZ 均匀背景上
比较有限几何 (6,6,0)、(12,3,0)、(12,3,3) 的全部非零允许动量。
工具独立于生产采样程序，实际参数和源码摘要写入输出。

每条键只存一次，A→B，零通量和 Wilson=(+1,+1)：
`A_AB=2T`，`E_matter=-sum singular_values(T)=-Tr|iA|/4`。
实际 hopping 是 `Jeff_mu+g*deltaX_mu`，其中 `Jeff=J+g*Xbar`；
均匀位移不重复加入。能量为
`E_static=E_matter+spring_k/2*sum (Xbar+deltaX)^2`。
每原胞二阶导数使用正交归一实 cos/sin 模式：
`sum v_alpha*v_beta/Ncell=delta_alpha_beta`，所以声子 Hessian 为 `k*I`。
普通 q 保留六个键/相位模式，自共轭 q 只保留三个实际 cosine 模式。

有限响应用真实键表的单粒子占据带二阶微扰；独立能量差分检查完整
混合矩阵和组合方向，实际使用三档步长。Bloch 参考保留实调制的两个
Fourier 顶角，先与相同有限 k 网格逐元素对照，再使用 200/400 midpoint
积分；直接求 k+q，不吸附节点。批处理解析顶角用于控制耗时和内存。
两档参考的矩阵差和最小曲率差都保存为实际数值，不用网格标签冒充检查。
三种几何共有的、最接近 Q 的允许点另作同动量对比。

CLI 默认采用 note 的 J=(1,1,.5)、g=.3、mass=1.7、omega=.8，
`k=mass*omega^2=1.088`；显式 `--spring-k` 只作一致性检查。
旧 Python 函数的 `spring_k` 省略值保留 1 供已有调用使用；CLI 将实际 k
显式传给全部函数。`--finite-step` 实际产生 step、step/2、step/4 差分。
`--uniform-grid` 不再限于原来硬编码的 400/800 两个键。

## 范围与异常

本工具计算静态固定背景的能量 Hessian，包含 matter 反馈和弹性势能。
不包含声子量子动力学，mass/omega 在这里仅通过 k 进入。
它不等于 QMC 等时 S(q)，也不等于 ADR-009 的有限投影外场响应 chi。
负曲率只能表明该背景对该调制局部不稳定；正曲率不能排除一级转变或
其他背景的有序态。K 的逆不得自动充当有限频率 QMC susceptibility。

有限 torus 的非网格 q 明确拒绝；连续 q 只允许用于 Bloch 参考。
零模不跳过、不加展宽；标记 `zero_mode_uncertified`，不输出假 Hessian，
保留能量扫描。缺少认证点时 `soft_search_complete=false`。
`certified` 只表示此响应公式未遇到零模；不认证 BZ 积分误差或量子临界点。
全部简并最软 q 和软子空间保留；不把任意实特征向量称为唯一键模式。

导出的 q_candidates 文件是允许测量波矢及其共轭，可以供现有
`--vbs-q-file` 的 selected/both 等时测量接口读取。
现有 susceptibility 接口仍使用 all 网格；不因候选文件改采样、外场或 trial。
若要投影固定三键通道，应在相同相位分母下先合并完整 QMC 矩阵和误差，
不能逐构型最大化后再平均。静态候选不证明相变点或 QMC 峰。

## 验证与参考

本次 13 项静态专项测试、生产 C++ 模型矩阵独立导出、105 个有限 q 的
Bloch/实空间对照、三档完整差分和 200/400 同 q 积分已完成。
原始证据和命令在
`../../../fast-1006/audit/1010/static-vbs-note-20261010/`；实际结果见 STATUS。
没有为此修改或重跑 QMC 内核、已有 HPC 作业或历史结果。

外部核对来源（note 本身没有链接）：
- [SSH–Kitaev 论文](https://arxiv.org/html/2510.07383v2)，模型及静态/有限频率边界；三键场的临界组合。
- [NumPy eigh 官方文档](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html)，实对称/Hermitian 矩阵的特征值与列特征向量约定。
