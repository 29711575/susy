# 本次导出验证（2026-10-10）

这些检查针对本仓库导出的 quick1010，区别于 `../recorded-validation/` 的历史记录。

- 328个继承源码/文档文件分别与当前quick1010和SOURCE_MANIFEST的SHA256逐项一致；整理没有修改生产物理代码。
- GNU C/C++/Fortran、Eigen3、Boost、MKL 2026.1，单线程构建生产程序及线性极化率、模长观测量两个验证程序成功。外部依赖路径在构建日志中仅用于记录本机环境，不是源码依赖的固定路径；参考根README设置自己的依赖位置。
- `kitaev_ssh_vbs_susceptibility_validate`通过：源微分误差8.55e-11、自由响应误差6.66e-16、直接DFT误差1.19e-16。
- `kitaev_ssh_vbs_amplitude_validate`通过：直接DFT误差6.94e-17、模长源微分误差2.71e-10。
- `python3 -m unittest test_vbs_amplitude test_static_vbs_response`：25项通过。记录环境Python3.14、NumPy2.5.3、Matplotlib3.11.2；Matplotlib缓存目录提示不影响测试通过。
- 根README的四个冻结快照分析命令逐一运行成功，重新生成曲线和分解。与原审计目录同名CSV逐字节一致；列表在summary.json。L9τ8独立加权重构最大残差2.84e-14，L9/L6τ扫描4.26e-14，原始扫描参考8.53e-14；full−disc与既有χ中心值一致。

日志、`summary.json`和measurements中的validation/metadata共同保留证据。构建目录、可执行文件、认证凭据没有上传。

这些检查验证实现与离线加权复现，不认证旧生产链的热化、自相关、相位遍历或基态极限；本次没有提交新的HPC计算。
