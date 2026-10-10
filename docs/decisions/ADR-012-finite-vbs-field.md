# ADR-012：有限静态 VBS 共轭场与纵向响应

状态：accepted。用户 2026-10-06 要求做有限场对照并提交 20022 testq。
本次只实现常量外场；不实现 ADR-004 的虚时间驱动。

在允许环面动量 q=(n1/Lx,n2/Ly-shear*n1/(Lx*Ly)) 上固定
f_b=delta_{mu,z} cos(q.R+phi)，H_h=H-h sum_b f_b X_b=H-h Nc m。
X 保持原未重标度定义，h 不并入 g；两端 trial 态随 h 固定。
CLI `--vbs-field-h h --vbs-field-n1 n1 --vbs-field-n2 n2
--vbs-field-phase phi` 开启物理场及选定投影矩。显式 h=0 可做同通道零场对照，
要求 VBSall 和 vbs-susceptibility=1。默认无场保持既有输出和 RNG 路径。

单键 F=h*f 的精确声子核为
K_F(x,y)=K_0(x-F/k,y-F/k) exp(dt*F²/(2k))。
路径 log 权重增量为
sum_links F*c*(x+y)+M*F²*(dt-2c)/(2k)，c=tanh(omega*dt/2)/omega。
PhononPath 储存物理 forces，并在 action/localLogRatio 中计入该项。
强迫自由路径的协方差不变、均值为
mu_j=d/2*(exp(-omega*t_j)+exp(-omega*(tau-t_j)))
+F/k*(1-(exp(-omega*t_j)+exp(-omega*(tau-t_j)))/2)。
局部 conditional 与整键 freeSample 同时使用该强迫高斯分布，因此完整
声子目标/提议仍严格抵消，接受率保留 matter 比值。固定 trial 不改为 F/k。
非零场路径快照升级 version3 保存 forces；零场仍写 version2，读入兼容1/2/3。

保存中心 m、完整 score=d log W_h/dh、m*score、m²/m⁴，及固定复 Fourier
D_z 的 |D|²/|D|⁴。完整 score 包含路径独立的 h sum_b f_b² M(dt-2c)/k。
纵向响应 chi_parallel=<m*score>-<m><score>，相位重加权之后连接减法，
每次 jackknife 删除后重算均值乘积。B4_complex=<|D|⁴>/<|D|²>²。
B4 不把区块二阶矩平方冒充逐构型四阶矩；有限场 B4 不作零场 crossing。
复 Fourier chi_zz 与所选单实投影 chi_parallel 归一化不同，不混用。

h 非零显式破坏平移对称。允许 q 的直接 Fourier 二阶矩仍合法，但
窗口重建只能解释为平移平均，不能当作有限场固定相位构型的完整二阶外积。
序参量保留非零允许 q 的 raw |D|²；响应须扣相应通道的均值乘积。
所有 grid 原始矩继续保存以检查错选波矢/三键模式。

数学检查：强迫 Mehler 权重和有限差分 score；解析自由中心响应
chi_free=(1-exp(-omega*tau/2))/k，普通 cos 源归一化乘mean(f²)=1/2；
局部和整键一般 MH 比值一致、g=0 接受率；forces 快照回放。
零场两种 proposal 的已有 blocks 逐字节回归和显式 h=0 路径/既有观测量回归。
本轮日志在 audit/vbs-finite-field-controls-20261006/。

临界判断：固定 q/通道扫描 lambda，不按每个 lambda 的最大动量切换。
有限场峰仅给出 lambda_p(L,h,tau)。不同 h、投影、尺寸及 dt 对照后再圈出
候选区间；不能把一次网格最大值直接当作 lambda_c。两个 seed 是试探扫描，
不认证自关联或平衡，运行中应检查接受率/漂移并按独立链重分块。
