# ADR-025：signed prefix在Pfaffian前强制Green结构认证

日期：2026-10-21。状态：最小护栏已应用active源码；生产源码身份及本机native构建已更新；目标HPC运行时仍须重建。不是完整数值修复。

## 证据与问题

真实L6 C-before路径在等价factor拆分下出现signed phase±1且均success。普通thin QR约束Q†Q=I，却不保证Gaussian isotropy QᵀQ=0。oblique projector P²=P可使G²=I，即使G不反对称；rcond和小solve residual亦不能替代结构认证。PFAPACK只读上三角时，非反对称G会破坏用于定号的Pfaffian/determinant身份式前提。

独立低维污染反例中旧Green与Pfaffian均success，开启原有结构门槛后返回green_structure_exceeded。真实B重放的失败位置属于**中间signed prefix**，不能误写成最终measurement-center失败。实际路径和反例是不同证据。

## 决定

仅在`inc/pure_projector_fast.h`的`pureProjectorStableReferenceWeight`使局部greenOptions可变并设置`fail_on_green_structure=true`。继承原数值阈值、完整MH与一次无缓存重算后中止语义；不反对称化G、不投影Q、不强制正号、不修改Hamiltonian/trial/source。覆盖当前phonon signed入口、partial/final Green及其uncached恢复；magnitude-only候选不返回phase，不由此认证长期Gaussian结构。旧public参考Evaluator不在此局部护栏覆盖中，不可用于绕开失败。

active文件SHA256：`c341ddbd6fa15bdbd5c1c5eccfd1fcceb1c9c6e7c624ad31b8e1ea9d75ce3876`，与已验证audit overlay一致。修改前文件保存在`audit/1021/qmc-susy-resolution/active-guard/before/`。主智能体随后完成1019物理trial冻结源码加唯一护栏的归档，356文件，archive SHA256为`ecdd154319d914ec6ec868de730b7cd47a243c0e313c16f1ad22542a75c01c90`。本机native build/help均exit0，binary SHA256为`b2bfe8a14e30a7cda02b22bc73810c4cd83e2c0116b47c53e17a18be4cda4051`；PRODUCTION_SOURCE已指向该guard身份并明确long_projection_numerical_repair_unresolved。归档不含shifted Γ或canonical correction。已有提交保留其冻结版本。

## 验证与后果

audit overlay的污染反例以及原27个Gaussian固定输入、144个独立16维Fock事务通过；没有新QMC链或混合认证。本次文档阶段只核对active/overlay哈希及修改范围，不冒称重跑测试或完成目标HPC构建。护栏可把原silent success转为明确失败，不能恢复正确phase/权重，也不保证所有过门槛结果均正确。不能据此判定所有历史L9数据错误。

独立静态L6 blocked试验展示普通QR放大误差的机制；canonical correction Qc=Q−Q* (QᵀQ)/2最初提供静态可行性证据，未接入active。静态成功本身不是最终修复认证；同日后续动态参考记录如下。

同日后续证据：20022/workq六项独立50/80位固定路径诊断全部exit0，确定C-before及C-candidate符号均为−1，确认旧实现C-before的+1是错误结果。隔离逐步canonical correction恢复了C路径的相位与幅值，C-before完整中心矩阵相对MP误差降至约2–4e−15；C-candidate仍有贴近原门槛的表示误差，A/B仍存在病态prefix。该候选未晋升生产，ADR的实际决定仍仅为上述fail-closed护栏。详细证据见本次audit的`numerical/MP_REPORT_CN.md`、`center-comparison.json`及矩阵方向勘误。

证据：[guard审查与覆盖](../../audit/1021/qmc-susy-resolution/sampler/phase_guard/REPORT_CN.md)、[验证manifest](../../audit/1021/qmc-susy-resolution/sampler/phase_guard/MANIFEST.json)、[真实phase变体](../../audit/1021/qmc-susy-resolution/numerical/phase-variant-summary.csv)、[总报告](../../audit/1021/qmc-susy-resolution/REPORT_CN.md)。
