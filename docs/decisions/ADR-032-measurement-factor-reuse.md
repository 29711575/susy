# ADR-032: 测量阶段因子存储复用

日期：2026-10-27。状态：已实施，确定性与本地/目标冻结构建已验证，生产源已提升。

## 背景与决定

从已验证 ADR-031 冻结源继续优化测量阶段因子重复构建。即时前驱的242文件与当前工作区逐项匹配；归档 SHA-256 为 `1e6ba311974007c2f9e14c781c16e8d2aebe5852de953defbcffc907b010a4bb`。

只允许复用具有明确生命周期和所有权的存储；原复数元素、signed zero、幅相、数值门槛、参考恢复、MH/RNG及测量开关策略保持。`MajoranaGenerator` 的工厂保存私有不可变model及原始wrapper身份；公开flag或callable被替换后不能通过身份检查。`Snapshot` 私有 `all_spectral_` 记录来自直接谱构建的完整路径资格；局部更新采用旧资格与新节点谱状态的AND，不能恢复已失效资格。`OwnedMeasurementReuse` 只能由链从私有不可变状态颁发；公开PreparedPath不能提供复用凭证。测量仍生成half-step和Cayley，只复用已有merged full-step谱因子。

`cacheDiagnostics()` 增加 `measurement_spectral_reuses`、`measurement_full_rebuilds` 与 `measurement_merged_node_reuses`，分别记录所选测量分支尝试和复用节点数；含失败尝试，不是接受率或吞吐指标。2867项专项、cache/Fock、11组22次保存输入回放及本地/目标冻结构建已通过；分组性能及限制见审计报告。

## 验证与边界

审计目录 `audit/1027/measurement-factor-reuse/` 保留即时前驱归档、清单和生产选择文件。确定性独立参考专项、L6五个保存构型、历史及固定生成physical L9 tau20各三对回放全部通过；相同O3、单数学线程，L9串行交替且不与编译并行。L6只用于正确性，其时间与编译重叠不作为性能证据。测量与可选总能量时间分别报告；总能量只用于确定性正确性检查，生产默认仍关闭。

无新MC、预检或smoke chain。保存输入回放不证明热化、混合、完整sweep或同误差吞吐；目标运行时必要构建及源/库/二进制完整性核验全部通过，生产源已提升。现有任务保持原源与设置。

本轮主要收益是RSS：历史/physical配对节省中位80.535/80.742MiB；测量耗时分别增加0.2103%与降低2.1499%，pipeline均未改善。仅固定输入三对，不能声称整体或同误差吞吐加速。
