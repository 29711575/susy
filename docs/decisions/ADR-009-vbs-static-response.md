# ADR-009：VBS 对静态外场的有限投影响应

状态：accepted，用户于 2026-10-05 明确要求 omega=0.3 扫描开启 VBS
序参量对外场的静态响应。只准备，不提交。实现/验证证据分别记录在 STATUS
和 `audit/production-L9-omega0p3-vbs-response-20261005/validation/`。

## 外场与有限路径定义

对允许的环面动量 q，采用 `H_h=H-h sum_r f(r) X_z(r)`；
两外端的 trial 态随 h 固定，中心节点测量。tau 是总开放轮廓长度，
中心投影时间为 tau/2。响应是这个有限投影期望值的导数；不假设基态极限。
f=cos(q.r) 对应 Re D_z，f=-sin(q.r) 对应 Im D_z，
其中 `D_mu(q)=sum_r exp(-iq.r)X_mu(r,center)/Nc`。

声子链使用精确 Mehler 核。线性场下的单键核满足
`K_h(x,y)=K_0(x-h*f/k,y-h*f/k)*exp(dt*h*h*f*f/(2*k))`。
在零场处 log 核导数为 `f*(x+y)*tanh(omega*dt/2)/omega`。
因此定义 c=tanh(omega*dt/2)/omega，节点权重两端 c、内部 2c。
试探态固定；matter 乘子没有显式 h 导数。
这给出所用分片路径权重的静态响应，dt 收敛仍须另行验证。

令 `I_mu(q)=sum_j w_j F_mu(q,j)/Nc`，则保存原始矩
`Nc*D_mu*conj(I_nu)`，最后用共同复相位分母重加权并做
`chi_mu,nu=<Nc*D_mu*conj(I_nu)>-Nc*<D_mu>*<conj(I_nu)>`。
不以 `tau*S(q)` 或等时涨落直接代替静态响应。
均值乘积必须在原始矩合并后减去，每个允许 q 都做连接减法；
不额外减去自由声子响应。有限样本不强行 Hermitize 响应矩阵。

## 归一化与复相位

实 cos/sin 源还独立保存 Re/Im D_z、Re/Im I_z 及其乘积。
`chi_cos=Nc*Cov(Re D_z,Re I_z)`，
`chi_sin=Nc*Cov(Im D_z,Im I_z)`；
实重加权统计下，复 Fourier 通道的 `Re chi_zz=chi_cos+chi_sin`。
一般复相位有限样本下，混合 quadrature 的虚部会使实部也有差异；
不强制二者相等。另导出 `quadrature_zz=chi_cos+chi_sin` 及联合误差。
特别是复相位重加权后 `<conj I>` 不一般等于 `conj(<I>)`，
`<Re D>` 也不一般等于 `Re(<D>)`。因此分开保存共轭/实分量的一阶矩，
不在重加权后的均值上取共轭或实部冒充对应均值。

## 软件接口与统计

新增 `--vbs-susceptibility 0|1`，默认 0；开启需 `--vbs 1 --vbs-mode all`。
初版只测有限环面真正允许的动量，不接入 selected/both 非公度通道。
每 q 增加 24 个复原始矩：D[3]、I[3]、Iconj[3]、Raw[9]、
DRe_z/DIm_z/IRe_z/IIm_z[4]、CosRaw_zz/SinRaw_zz[2]。
CSV 保存共同 phase 和 phase*observable 的区块和，JSON 标注
`vbs-static-response-v1` 与精确插入权重。
Fourier 线性性允许先在实空间积分路径再做 FFT；只需中心和积分场两次变换。
测量不抽 RNG，不改提议/接受/权重路径；关闭开关保留旧接口行为。

汇总在原始矩层面合并独立 seed，然后连接减法；跨 seed 删除整链估误差。
单链的联合删区块估计在每次删除后重新计算全部 phase 比值和连接乘积，
逐 seed 重分块，不跨链拼块，不静默舍尾。不以区块误差数值宣称自关联已解决。

确定性检查应包含精确强迫核的有限差分 log 导数、自由真空端点的解析响应
`chi_free=(1-exp(-omega*tau/2))/k`、空间 Fourier/实源归一化及 CLI 开关回归。
lambda 网格上的最大值只用于后续定位峰，不能据此认证临界点或热力学相变。
