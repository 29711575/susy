严格静态极限：全键松弛与单畴极化率（2026-10-10）

入口 kitaev_ssh/static_limit_response.py，仅依赖 NumPy/SciPy，绘图另需 Matplotlib。
从 /home/asus/newsusy 执行，输出目录必须新建或为空：

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLCONFIGDIR=/tmp/newsusy-static-limit-mpl \
/home/asus/.venv/bin/python quick1010/kitaev_ssh/static_limit_response.py \
  --Lx 9 --Ly 9 --shear 0 --Jx 1 --Jy 1 --Jz .5 \
  --lambda 2.4 --spring-k 1 --parity trial --starts 80 --seed 106101700 \
  --output-dir quick1010/audit/1010/YOUR_NEW_TASK --plot

固定 k,g 下取 m→∞；能量无声子动能，不需要 mass/omega/tau/dt/MC sweep。
E(X)=-sum singular_values[T(J+gX)]+k/2 sum X²，固定零通量、Wilson(+,+)。
默认保留自由零位移试探态的 matter 奇偶扇区；若 det T 换号，加入最低
单粒子奇偶激发的 2 sigma_min 能量及梯度。--parity unrestricted 为显式诊断。
零模或不可微的奇偶激发简并拒绝输出已认证的局部响应。

允许所有 3LxLy 个实键变量独立变化，以均匀、调制、随机初值找局部最低能。
解析全键 Hessian 保存为 hessian.npy；三键 Fourier 响应为
chi_mu,nu(q)=F_mu(q) H_X^-1 F_nu(q)^*/Nc，D=sum exp(-iq.r)X/Nc。
它是某个稳定静态畴的线性外场导数，包含弹性和 matter 反馈，无自由声子减法。
全允许 q 含 Gamma；逐 q 最大值不是 lambda 扫描的临界峰。

results.json 保存参数、能量、稳定性、畴、三键矩阵及源码身份；
susceptibility.csv 保存全 q 标量；equilibrium_bonds.csv 为静态构型；
multistart.csv 和 multistart_configurations.npz 保留全部优化候选。
validation.json 检查独立能量归一化、梯度/Hessian、自由响应、正负外场。
--plot 只产生 PNG。确定性数值不附 Monte Carlo 统计误差。

若静态构型已破坏平移对称，严格零场有简并平移畴。外场在不同符号下
选择不同畴，因此整体最小化得到的 h→0 左右极限可能不同，不能把
单畴 Hessian 数值冒充对称恢复后的唯一有限极化率。
多初值搜索不证明全局最小；固定规范 matter 也不等于完整物理自旋投影。

本次实际 L9/lambda2.4 报告及独立 Fock/占据带/平移检查：
../fast-1006/audit/1010/L9-strict-static-lambda2p4-20261010/REPORT.txt
设计见 docs/decisions/ADR-020-strict-static-limit-response.md。
旧 static_vbs_response.py 为固定共同背景曲率工具，继续独立保留。

2026-10-10 L9/lambda2–3扫描与PNG：
audit/1010/L9-strict-static-lambda2-3-scan-20261010/
lambda步长.05、21点，每点80初值，另将21个最低找到构型在所有点交叉松弛。
reviewed/为验收后数据；scan_summary.csv为主绘图输入，plot_scan.py可重画。
三个软模点验证使用更小外场；validation()新增可选field_strengths参数，
默认及数值容差不变。分支切换处曲线不连接，不把采样最高点认定为临界点。
