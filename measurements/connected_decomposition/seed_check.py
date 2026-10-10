import csv,json,os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/chi-full-disc-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
r=list(csv.DictReader((P/'seed-first-moments.csv').open()))
fig,axes=plt.subplots(1,3,figsize=(14,5),layout='constrained');checks=[]
for ax,tau in zip(axes,[10,12,14]):
 rr=[x for x in r if x['campaign']=='tau_scan' and x['L']=='9' and float(x['lambda_value'])==3 and float(x['tau'])==tau and x['selection']=='fixed_target']
 weights=np.array([complex(float(x['complex_phase_weight_real']),float(x['complex_phase_weight_imag'])) for x in rr])
 D=np.array([complex(float(x['mean_D_real']),float(x['mean_D_imag'])) for x in rr]);Iconj=np.array([complex(float(x['mean_Iconj_real']),float(x['mean_Iconj_imag'])) for x in rr]);meanD=(weights*D).sum()
 n=np.array([int(x['retained_sweeps']) for x in rr]);original=np.array([int(x['seed'])<106000000 for x in rr])
 for select,col,label in [(original,'#2463a0','original seeds'),(~original,'#d58b1c','supplemental seeds')]:
  ax.scatter(D[select].real,D[select].imag,s=30+220*n[select]/n.max(),color=col,alpha=.75,edgecolors='white',label=label)
 ax.scatter([meanD.real],[meanD.imag],marker='*',s=220,color='#a65087',edgecolors='#222222',label='pooled mean',zorder=5)
 ax.axhline(0,color='#aaaaaa',lw=.6);ax.axvline(0,color='#aaaaaa',lw=.6);ax.set_aspect('equal',adjustable='datalim')
 ax.set_xlabel(r'$\mathrm{Re}\langle D\rangle_{\rm seed}$');ax.set_ylabel(r'$\mathrm{Im}\langle D\rangle_{\rm seed}$');ax.set_title(r'$\tau='+str(tau)+'$'+f' | retained {n.min()}–{n.max()} sweeps')
 checks.append(dict(tau=tau,L=9,lambda_value=3,q=49,retained_sweeps_min=int(n.min()),retained_sweeps_max=int(n.max()),pooled_disc=float((81*(weights*D).sum()*(weights*Iconj).sum()).real),phase_weighted_mean_seed_disc_diagnostic=float((81*np.sum(weights*D*Iconj)).real),largest_seed_disc=float(max(float(x['C_disc_seed']) for x in rr)),largest_seed_retained=int(max(rr,key=lambda x:float(x['C_disc_seed']))['retained_sweeps'])))
axes[0].legend(frameon=False,fontsize=9)
fig.suptitle(r'$L=9,\lambda=3,\omega=0.3$ | first moment at fixed $q=(5/9,4/9)$',fontsize=15)
fig.supxlabel('Dot size reflects retained sweeps. Different seed means can cancel in the pooled mean; short-chain fluctuations and VBS-sector trapping are not distinguished.',fontsize=9)
fig.savefig(P/'L9-lambda3-seed-first-moments.png',dpi=200);plt.close(fig)
(P/'seed-cancellation-diagnostic.json').write_text(json.dumps(checks,indent=2)+'\n')
report=json.loads((P/'metadata.json').read_text())
text='''当前数据的 connected 极化率分解（冻结快照复核；未重新查询HPC）

口径：C_full=Re<ChiRaw>=Nc Re<D Iconj>，C_disc=Nc Re(<D><Iconj>)，chi=C_full-C_disc。
所有系综平均先合并保留样本的复phase分子/分母。没有平均各seed已扣除的chi，没有把复相位下<Iconj>替换为conj(<I>)。
固定目标是最靠近裸Dirac点差的允许动量代表：L9 n=(5,4), L6 n=(3,3)，并非直接测量连续裸Q。
125参数点的全部q减法身份与此前绘图chi逐点一致（本轮计算残差均0）。

L9 tau扫描（tau>=8，固定目标）扣除比例绝对值中位数0.685%、最大3.554%；connected峰处最大3.598%。
L6 tau10/14固定目标最大0.694%，connected峰处最大1.950%。
L9 tau8 lambda扫描固定目标最大5.788%（lambda3.6），connected峰处最大7.865%（lambda5）。这些为当前点估计，不是误差上界；部分C_disc误差与点值相当或更大。
因此当前曲线主要形状并非disconnected减法产生；完整相关项本身已呈现同样走势。没有由此确认物理相、低能机制、基态或热化/遍历收敛。

L9 lambda3固定目标 q49：
 tau    C_full       C_disc      chi         disc/full
 8      1.053968     0.004128    1.049840     0.392%
10      1.373947     0.011428    1.362519     0.832%
12      1.222786     0.003656    1.219129     0.299%
14      1.107627     0.010174    1.097453     0.919%
此表tau8是tau扫描的4seed、cut1000，不与另外tau8lambda扫描的10seed、cut500混合。

单seed检查：lambda3,tau10有原始链保留3810sweep时C_disc_seed约0.4764，合并后仅0.01143；逐seed< D >在复平面中分散并在合并后抵消。
补充链中的最大单seed值不能用于证明遍历性：例如tau10某链仅保留10sweep，其C_disc_seed约1.542；tau12/14部分链仅保留约330–460sweep。单seed与合并的disconnected差别可能来自有限样本、自关联、未热化或不同VBS相位采样，当前聚合矩不能区分。

误差：删整条seed的jackknife1SE，每个副本重新合并并求均值乘积；分解在同一个选定q计算（峰值q在副本中固定，便于保持同一q加法关系）。之前峰曲线jackknife每副本重选最大q，两种误差口径不同，中心chi相同。
文件：decomposition.csv为固定q/connected峰，all-q-decomposition.csv为完整动量网格，seed-first-moments.csv为每seed一阶矩及诊断；metadata.json存源文件与SHA256。
'''
(P/'REPORT.txt').write_text(text)
print(json.dumps(checks,indent=2))
