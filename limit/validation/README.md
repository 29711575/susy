# 静态数据与源码验证记录

这里收录已有冻结结果的验证摘要，原 JSON 保持原字节：

- [L6/L9 原始摘要](L6-L9-original-verification.json)：138点，独立能量、响应与外场核对。
- [L7/L8 原始摘要](L7-L8-original-verification.json)：70点，后续采用准确响应适配器重新认证。
- [响应源码身份](response-source-manifest.original.json)：解析固定宇称 Hessian 与 MP40 回退来源。
- [响应精度审计](response-precision-audit.original.json)：280个实源中268个按绝对界认证，
  12个按软响应相对界认证；保留逐点绝对不确定性与比较界。

这些原始记录含本机/HPC 旧路径，仅表示出处；离线读取 CSV 或作图不依赖这些路径。
L7/L8 原摘要中的较大绝对外场差值，应结合后续软点的相对误差界解释，不能只看单个绝对门槛。
这里没有重新执行历史208点优化或全量外场认证，也没有改变其结果。

迁移时执行的检查见 [package-checks.json](package-checks.json)：
源文件与冻结哈希一致；合并表逐字段匹配两批原表；所有点满足记录的物理宇称；
各尺寸的范围、分辨率和 m² 归一化核对；原有预测/曲率回归与独立 Fock/导数验证；
全文件清单、兼容入口及文档相对链接检查。

独立物理核验可从仓库根目录运行，不生成 MC 链：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 limit/source/validate_physics.py /tmp/susy-static-independent-checks.json
```

该核验检查小 Fock 矩阵、导数和有限几何，不证明全局最优或热力学相界。
