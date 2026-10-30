# 2026-10-30 导出验证

本次将已验证的 2026-10-28 数值核心导出到仓库，并从此检出重新构建主程序、
reference-core、endpoint 和谱安全验证器。实际结果见 [summary.json](summary.json)，
导出完整性、脚本及统计复算见 [export-checks.json](export-checks.json)。
245 个核心文件（不含仓库指导文件）逐字节匹配冻结源，见
[frozen-source-files.json](frozen-source-files.json)。

| 检查 | 记录 |
|---|---|
| 从当前检出原生构建 CLI 及三个确定性验证器 | [构建日志](native-build.log) |
| CLI help、动态依赖与确定性测试 | [summary.json](summary.json)、各 validator 日志 |
| Python 离线统计回归 | [日志](python-regression.log)：19项发现、15通过、4项缺少历史本地样本而跳过 |
| 数据和图重算 | 三组快照均在仓库内离线重现；四 q 的160个均值与jackknife误差另用标量实现复算 |
| 来源与文档 | 核心哈希、Python/Shell语法、新报告相对链接、旧测量文件未改动 |

安装依赖并设置根 README 所列环境后，可运行：

```bash
bash scripts/build_nodeps.sh build kitaev_ssh_mc \
  kitaev_ssh_reference_core_validate \
  kitaev_ssh_endpoint_propagation_validate \
  kitaev_ssh_spectral_safety_validate
./build/kitaev_ssh_mc --help
./build/kitaev_ssh_reference_core_validate
./build/kitaev_ssh_endpoint_propagation_validate
./build/kitaev_ssh_spectral_safety_validate
python3 -m unittest discover -s scripts/tests -p test_qmc_tools.py -v
```

[recorded-final-summary.json](recorded-final-summary.json) 与
[recorded-target-validation.json](recorded-target-validation.json) 是原工作区
10 月 28 日验收记录，包含其原路径、编译器和误差，不是本次重跑日志。
本次原生检查使用 GNU + MKL sequential 单线程，不运行 MC 链，也不操作 HPC。
保存路径回放和历史速度记录的适用范围见 ADR-030–033；通过这些检查不能认证
平衡、混合、投影/时间步长收敛、临界性或整体计算提速。
