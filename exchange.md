# exchange

2026-11-07。本文简述本地 quick1010 已验证冻结实现中的 λ 副本交换与交叉权重计算。该功能已有代码，但刚才展示的 L=9、ω=0.1、τ=30/40 共40条链没有启用交换，仍是独立链。

**版本边界：** 本文依据 `20261031-qr-stability` 冻结源码；GitHub 当前普通核心仍是另一导出版本，尚未包含本文的 `replica_exchange.h` / `file_replica_exchange.h`。本次只添加说明，不同步源码、不修改运行任务。仓库版本说明见 [README](README.md) 与 [源码更新记录](docs/updates/20261030-source-and-results.md)。

## 1. 四个权重与接受率

相邻参数槽 A、B 分别持有完整声子路径 X、Y。每次交换需要：

| 权重 | 含义 | 文件交换协议的取得方式 |
|---|---|---|
| W_A(X) | A 参数下自己的路径 | `ensureMeasured()` 后绑定当前状态 |
| W_B(Y) | B 参数下自己的路径 | 对方同样取得，再通过文件传入 |
| W_A(Y) | A 参数下评估对方路径 | A 调用 `exchangeWeight(Y)` 完整重算 |
| W_B(X) | B 参数下评估对方路径 | B 调用 `exchangeWeight(X)` 完整重算 |

`ReplicaFullWeight` 保存声子对数权重、matter 对数模长、复相位与有效性标志。定义

$$
\ell_A(X)=\log |W_A(X)|
=\mathrm{log\_phonon}_A(X)+\mathrm{log\_abs\_matter}_A(X).
$$

`logSwapRatio()` 用四个对数权重计算

$$
\Delta=\ell_A(Y)+\ell_B(X)-\ell_A(X)-\ell_B(Y),\qquad
P_{\rm acc}=\min(1,e^\Delta).
$$

代码以 `long double` 累加差值，并比较 `log(u) < min(0, Δ)`，避免直接相乘大权重。接受率采样绝对权重乘积；相位必须有效并通过检查，但不放入交换接受率，观测量仍按各参数槽的相位重加权。

## 2. 交叉权重实际怎样算

实际文件协议调用链是：

```text
FileReplicaExchange::observe()
  → PhononSampler::exchangeWeight(peerPath)
  → CoupledPhononPath::uncachedReference(peerPath, true, true)
  → pathFactors() → evaluatePath()
```

路径坐标保持不变，使用接收槽自己的模型生成器、耦合 g、trial 和 QR 配置重新构造传播因子，计算完整 signed matter 权重。两个 `true` 表示按测量路径计算完整相位，并计算中心 Green。

**这里没有用当前路径的 prefix/Cayley 缓存、合并权重顺序或 prepared trial 数据来近似交叉权重。** `uncachedReference()` 从新 canonical 半步因子作完整传播；一般 sweep 中开启的缓存加速，不能直接理解为交换权重也享有同样的增量加速。

另有 `ReplicaExchangeTarget::evaluate()`，通过新建目标自己的 `CoupledPhononPath` 再 `ensureMeasured()` 求完整权重，主要供独立验证使用；不要与上面的实际文件协议入口混淆。

## 3. 接受交换以后

协议交换完整路径及 walker 标签，参数槽、trial 与槽内采样器保持所属关系。接受后调用

```text
installExchangePath(peerPath, cross.logAbsTotal())
  → proposeQModePath() → 在本槽参数下重建候选 → accept()
```

安装时重新计算的总对数模长必须与已认证交叉权重相差不超过 `1e-7`，否则报错。它没有直接搬用对方的投影缓存，也没有直接安装刚才交叉评估产生的整个 snapshot。无效权重或未恢复的数值故障会中止交换链组，不能伪装成普通 MH 拒绝。

## 4. 支持范围与额外成本

当前协议只实现相邻 **λ 梯子**：保持声子参数不变，通过 `g=sqrt(λ·m·ω²)` 改变耦合，要求零静态位移背景。两槽必须有一致的几何、J、规范、trial policy、实际宇称和精确 trial 基，以及相同的时间切片数、Δτ、声子质量/频率、trial shift 和源项；因此它不直接支持跨 τ、ω 或不同 trial 的交换。

- 每对有配对的副本需要两个完整交叉重算，双方可以并行；已有自权重如未完成相位/Green，还要补算。
- 接受后双方再各自重建一次候选状态；交叉评估已计算的 Green 也不是免费复用到安装状态。
- 文件协议传输完整路径、权重和确认包，按奇偶轮交替配对，需要等待对端并清理文件；有同步、序列化与共享文件系统开销。
- 辅助 λ 槽占用额外 CPU、内存和热化时间。收益应比较全部槽与等待成本计入后的误差或 ESS/CPU小时，不能只看交换接受率。

CLI 的 `--exchange-every` 默认是0。启用还需指定槽号、λ列表、共享目录和独立交换 seed。交换只在 measurement 阶段按 cadence 触发；当前不支持单槽 checkpoint resume，不能独立续算某一个交换槽来代替整个链组恢复。

## 5. 源码与已有验证

本文源码归档：本地 `quick1010/artifacts/production/20261031-qr-stability/source.tar.gz`，SHA256：

```text
c8f12bb0b57697b1555aa1099f99ef7b76829f3d7dd78ebf9386ff6a9d4c7e99
```

归档内主要入口为 `kitaev_ssh/{replica_exchange.h,file_replica_exchange.h,phonon_mc.h,phonon_projector.h,coupled_mc.cpp}`。目标/完整路径指纹用于绑定权重身份，是实现内部的一致性检查，不是密码学摘要。

已有独立验证源码 `replica_exchange_validate.cpp` 覆盖小系统 Fock 四权重对照、详细平衡、trial/τ/源项不兼容拒绝、缓存重建和无效权重处理；`file_replica_exchange_validate.cpp` 覆盖文件协议的接受/拒绝、walker 与清理。2026-10-25 本地 `audit/1025/L9-lambda3p7-tau10-expanded-controls-20015/native/exchange-validation.log` 留有 PASS 记录。这是历史固定输入正确性证据，不是当前 L9 参数下的混合效率或性能认证；本次未新增蒙特卡洛计算。
