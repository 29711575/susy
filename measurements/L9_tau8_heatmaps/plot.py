#!/usr/bin/env python3
import csv,json,os,hashlib
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/chi-heatmap-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.path import Path as MplPath
P=Path(__file__).resolve().parent
S=P.parent/'L9_tau8_lambda'
LAMBDAS=[1.3, 1.6, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.825, 2.85, 2.875, 2.9, 2.925, 2.95, 2.975, 3.0, 3.1, 3.2, 3.4, 3.6, 4.0, 4.5, 5.0]
rows=list(csv.DictReader((S/'allowed-grid.csv').open()))
summary=list(csv.DictReader((S/'chi-vs-lambda-tau8-cut500.csv').open()))
arrays=[];selected=[];metadata=[]
for lam in LAMBDAS:
 rr=[r for r in rows if float(r['lambda_value'])==lam]
 assert len(rr)==81 and len({int(r['q']) for r in rr})==81
 a=np.full((9,9),np.nan)
 for r in rr:
  n1,n2=int(r['n1']),int(r['n2']);assert int(r['q'])==9*n1+n2
  a[n2,n1]=float(r['chi_real'])
 assert np.isfinite(a).all()
 residual=max(abs(a[n2,n1]-a[(-n2)%9,(-n1)%9]) for n1 in range(9) for n2 in range(9))
 assert residual<1e-10
 point=next(r for r in summary if float(r['lambda_value'])==lam)
 nonzero=a.copy();nonzero[0,0]=-np.inf
 assert abs(nonzero.max()-float(point['chi_max']))<1e-12
 metadata.append({'lambda':lam,'seeds':int(point['seeds']),'retained_samples':int(point['samples']),
                  'max_nonzero_chi':float(nonzero.max()),'q_minus_q_real_residual':residual,
                  'peak_indices_n1_n2':[[int(x),int(y)] for y,x in zip(*np.where(np.isclose(nonzero,nonzero.max(),rtol=0,atol=1e-12)))]})
 arrays.append(a);selected.extend(rr)
with (P/'heatmap-values.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(selected[0]));w.writeheader();w.writerows(selected)
vmin=min(a.min() for a in arrays);vmax=max(a.max() for a in arrays)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
# Cartesian reciprocal geometry: nearest-neighbour bond length a=1.
B=np.column_stack(([2*np.pi/np.sqrt(3),2*np.pi/3],[-2*np.pi/np.sqrt(3),2*np.pi/3]))
hexv=np.array([[4*np.pi/(3*np.sqrt(3)),0],[2*np.pi/(3*np.sqrt(3)),2*np.pi/3],[-2*np.pi/(3*np.sqrt(3)),2*np.pi/3],[-4*np.pi/(3*np.sqrt(3)),0],[-2*np.pi/(3*np.sqrt(3)),-2*np.pi/3],[2*np.pi/(3*np.sqrt(3)),-2*np.pi/3]])
inside=MplPath(np.vstack([hexv,hexv[0]]))
# Clip each exact triangular-grid Voronoi hexagon against the first BZ.
def clip(poly):
 for edge in range(6):
  origin=hexv[edge];direction=hexv[(edge+1)%6]-origin
  def side(v):return direction[0]*(v[1]-origin[1])-direction[1]*(v[0]-origin[0])
  out=[]
  if len(poly)==0:break
  previous=poly[-1];sp=side(previous)
  for current in poly:
   sc=side(current)
   if (sp>=-1e-12)!=(sc>=-1e-12):out.append(previous+(current-previous)*sp/(sp-sc))
   if sc>=-1e-12:out.append(current)
   previous=current;sp=sc
  poly=np.asarray(out)
 return poly
polygons=[];cell_indices=[];dots=[];dot_indices=[]
for n1 in range(9):
 for n2 in range(9):
  for i in range(-2,3):
   for j in range(-2,3):
    center=B@np.array([n1/9+i,n2/9+j])
    polygon=clip(center+hexv/9)
    if len(polygon)>=3:
     area=abs(np.sum(polygon[:,0]*np.roll(polygon[:,1],-1)-polygon[:,1]*np.roll(polygon[:,0],-1)))/2
     if area>1e-12:polygons.append(polygon);cell_indices.append((n1,n2))
    if inside.contains_point(center,radius=1e-9):dots.append(center);dot_indices.append((n1,n2))
assert len(set(cell_indices))==81
area=sum(abs(np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1)))/2 for p in polygons)
bz_area=abs(np.linalg.det(B));assert abs(area-bz_area)<1e-10
# Bare Bloch zero: f(k)=exp(i k1)+exp(i k2)+0.5=0.
theta=np.arccos(-.25)
assert abs(2*np.cos(theta)+.5)<1e-14
Q=np.array([theta/np.pi,1-theta/np.pi])
stars=[];star_signs=[]
for sign in [1,-1]:
 for i in range(-2,3):
  for j in range(-2,3):
   point=B@(sign*Q+np.array([i,j]))
   if inside.contains_point(point,radius=1e-9):stars.append(point);star_signs.append(sign)
stars=np.asarray(stars);dots=np.asarray(dots)
assert len(stars)==4
norm=Normalize(vmin=vmin,vmax=vmax)
def panel(ax,a,m):
 values=[a[n2,n1] for n1,n2 in cell_indices]
 im=PolyCollection(polygons,array=np.array(values),cmap='cividis',norm=norm,edgecolors='none',rasterized=True)
 ax.add_collection(im)
 outline=np.vstack([hexv,hexv[0]])
 ax.plot(outline[:,0],outline[:,1],color='#333333',lw=1.2,zorder=3)
 ax.scatter(dots[:,0],dots[:,1],s=9,color='#252525',alpha=.65,zorder=4)
 for point,index in zip(dots,dot_indices):
  if list(index) in m['peak_indices_n1_n2']:ax.plot(*point,'o',ms=12,mfc='none',mec='white',mew=1.7,zorder=5)
 ax.scatter(stars[:,0],stars[:,1],marker='*',s=160,facecolors='#52e8ff',edgecolors='#222222',linewidths=.8,zorder=6,clip_on=False)
 for point,sign in zip(stars,star_signs):
  if point[1]>0:ax.annotate(r'$Q_{\rm bare}^{'+('+' if sign==1 else '-')+'}$',point,xytext=(0,15),textcoords='offset points',ha='center',va='bottom',fontsize=11)
 ax.text(0,.11,r'$\Gamma$',ha='center',va='bottom',color='white',fontsize=12,zorder=6)
 ax.set_aspect('equal');ax.set_xlim(-2.57,2.57);ax.set_ylim(-2.35,2.6)
 ax.set_xticks([-2,0,2]);ax.set_yticks([-2,0,2])
 ax.set_xlabel(r'$q_x$ ($a^{-1}$)');ax.set_ylabel(r'$q_y$ ($a^{-1}$)')
 for spine in ax.spines.values():spine.set_visible(False)
 ax.set_title(r'$\lambda='+str(m['lambda'])+'$'+f"\n{m['seeds']} seeds, {m['retained_samples']} retained samples",fontsize=12)
 return im
fig,axes=plt.subplots(9,3,figsize=(16,44),layout='constrained')
for ax,a,m in zip(axes.flat,arrays,metadata):im=panel(ax,a,m)
fig.colorbar(im,ax=axes,shrink=.75,label=r'$\mathrm{Re}\,\chi_{zz}(\mathbf{q})$ (connected)')
fig.suptitle(r'$L=9,\ \omega=0.3,\ \tau=8$ | discard first 500 sweeps per chain',fontsize=18)
fig.supxlabel('Dots: measured q; circles: nonzero-q maxima; stars: bare Dirac-point differences (Jx=Jy=1, Jz=0.5).\nFirst hexagonal Brillouin zone; periodic Voronoi cells, no interpolation; Gamma included. Shared color scale.',fontsize=12)
fig.savefig(P/'chi-heatmaps-lambda1p3-5.png',dpi=150,facecolor='white');plt.close(fig)
for start in range(0,len(arrays),9):
 fig,axes=plt.subplots(3,3,figsize=(16,16),layout='constrained')
 for ax,a,m in zip(axes.flat,arrays[start:start+9],metadata[start:start+9]):im=panel(ax,a,m)
 fig.colorbar(im,ax=axes,shrink=.75,label=r'$\mathrm{Re}\,\chi_{zz}(\mathbf{q})$ (connected)')
 fig.suptitle(r'$L=9,\ \omega=0.3,\ \tau=8$ | discard first 500 sweeps per chain',fontsize=18)
 fig.supxlabel('Dots: measured q; circles: nonzero-q maxima; stars: bare Dirac-point differences (Jx=Jy=1, Jz=0.5).\nPeriodic Voronoi cells, no interpolation; Gamma included. Shared color scale across all 27 lambda points.',fontsize=12)
 fig.savefig(P/('chi-heatmaps-page'+str(start//9+1)+'.png'),dpi=180,facecolor='white');plt.close(fig)
for a,m in zip(arrays,metadata):
 fig,ax=plt.subplots(figsize=(6.5,6),layout='constrained');im=panel(ax,a,m)
 fig.colorbar(im,ax=ax,label=r'$\mathrm{Re}\,\chi_{zz}(\mathbf{q})$')
 fig.suptitle(r'$L=9,\ \omega=0.3,\ \tau=8$ | cut 500 sweeps',fontsize=13)
 fig.supxlabel('Dots: measured q; circles: maxima; stars: bare Dirac-point differences.\nPeriodic Voronoi cells, Gamma included; shared color scale.',fontsize=9)
 fig.savefig(P/('chi-heatmap-lambda'+str(m['lambda']).replace('.','p')+'.png'),dpi=220,facecolor='white');plt.close(fig)
report={'source':str(S),'source_grid_sha256':hashlib.sha256((S/'allowed-grid.csv').read_bytes()).hexdigest(),
        'source_collection':json.load(open(S/'metadata.json'))['sources'],'discard_sweeps':500,
        'definition':'Real connected chi_zz on all81 allowed q, joint phase raw moments before connected subtraction; Gamma included',
        'coordinates':'Cartesian q in first hexagonal BZ, bond length a=1; exact periodic nearest-grid Voronoi cells clipped at BZ boundary; all81 q including Gamma; no interpolation',
        'reciprocal_vectors_columns':B.tolist(),'bare_Q_fractional':Q.tolist(),'bare_Q_cartesian_markers':stars.tolist(),
        'bare_Q_definition':'Q_bare=+/-2K modulo reciprocal lattice; theta=acos(-Jz/(2Jx)), Jx=Jy=1,Jz=0.5',
        'geometry_validation':{'represented_q':len(set(cell_indices)),'cell_area_sum':area,'BZ_area':bz_area,'area_residual':abs(area-bz_area)},
        'color_scale':{'shared':True,'min':float(vmin),'max':float(vmax)},'points':metadata,
        'limitations':'Same latest snapshot as the preceding lambda plot; partial/stopped chains included; no convergence certification'}
(P/'metadata.json').write_text(json.dumps(report,indent=2)+'\n')
(P/'README.txt').write_text(json.dumps(report,indent=2)+'\nReproduce: python3 plot.py\n')
print(json.dumps(report,indent=2))
