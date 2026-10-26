# 2026-10-26 源码导出验证

本记录对应本次导出的主线源码；历史 evidence 和 measurements 保留各自版本。
341个源码、构建及验证文件先冻结，再从冻结目录编译，最后逐文件确认与导出源码一致。
文件身份见 [frozen-source-files.json](frozen-source-files.json)，检查命令/退出码见
[checks.json](checks.json)，本机二进制哈希与范围见 [summary.json](summary.json)。

| 检查 | 结果 |
|---|---|
| Python AST、Shell bash -n | 通过，覆盖冻结源中的全部相应文件 |
| 离线统计 unittest | 19项发现，15项通过；4项缺少历史本地 audit 样本而跳过，见 [日志](python-regression.log) |
| CLI/能量/physical trial 原生构建 | 从本次冻结源用 GNU、Eigen3、Boost、MKL sequential 重新构建，exit0；见 [构建日志](native-build.log) |
| CLI help | exit0；确认 physical 默认和 total energy 默认0，见 [接口](cli-help.log) |
| 总能量确定性验证 | exit0，独立谱/Fock/固定输入验证，见 [日志](total-energy-deterministic.log) |
| physical trial 确定性验证 | exit0，奇偶/几何/规范/固定输入 Fock 检查，见 [日志](physical-trial-deterministic.log) |

可在安装构建依赖后，从仓库根目录重现：

```bash
python3 -m unittest discover -s scripts/tests -p test_qmc_tools.py -v
bash scripts/build_nodeps.sh build kitaev_ssh_mc \
  kitaev_ssh_total_energy_validate kitaev_ssh_physical_trial_validate
./build/kitaev_ssh_mc --help
./build/kitaev_ssh_total_energy_validate
./build/kitaev_ssh_physical_trial_validate
```

构建依赖环境变量和动态库路径按根 README 配置。日志保留实际命令和本地冻结路径，
这些路径是来源记录。本次未运行MC链，未操作HPC，也没有新的目标运行时验证；
不据此认证热化、混合、投影收敛或任意长投影的普适数值正确性。
