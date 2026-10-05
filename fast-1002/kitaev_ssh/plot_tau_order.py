#!/usr/bin/env python3
"""Plot finite-size VBS proxies versus tau using the existing joint jackknife."""
import argparse,csv,json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from convergence_analysis import fixed_momentum_mapping,joint_estimate,read_chain

def transform(sums,count,names,meta):
    if abs(sums[0])<=1e-12*count: raise ValueError("unresolved phase denominator")
    r=dict(zip(names,sums[1:]/sums[0])); nc=meta["Lx"]*meta["Ly"]
    q=fixed_momentum_mapping(meta); q1=q["q1"]["index"];q3=q["q3"]["index"]
    s=np.array([r[f"S_q{i}_zz"].real/nc for i in range(nc)])
    d=np.array([abs(r[f"Delta_q{i}_z"]) for i in range(nc)])
    m=max(float(s.max()),0.)
    return {"m2_q1":s[q1],"m2_q3":s[q3],"m2_grid_max":m,
            "rms_q1":np.sqrt(max(s[q1],0.)),"rms_q3":np.sqrt(max(s[q3],0.)),
            "rms_grid_max":np.sqrt(m),"abs_delta_q1":d[q1],
            "abs_delta_q3":d[q3],"abs_delta_grid_max":d.max()}

def calculate(manifest):
    m=json.loads(manifest.read_text())
    chains=[read_chain(Path(x["folder"])/"run_results.json") for x in m["commands"]]
    rows=[];peaks=[]
    for intervals in m["plan"]["intervals"]:
        group=[c.reblock(25) for c in chains if c.metadata["intervals"]==intervals]
        e=joint_estimate(group,transform=transform); tau=round(intervals*m["plan"]["dt"],12)
        for i,name in enumerate(e["names"]):
            rows.append({"tau_total":tau,"intervals":intervals,"observable":name,
                         "value":float(e["values"][i].real),
                         "stderr":float(np.sqrt(e["covariance"][i,i])),
                         "block_sweeps":25,"seeds":";".join(map(str,e["seeds"]))})
        sums=sum((c.sums.sum(axis=0) for c in group))
        r=dict(zip(group[0].names,sums[1:]/sums[0]));nc=group[0].metadata["Lx"]*group[0].metadata["Ly"]
        peaks.append({"tau_total":tau,
                      "m2_peak_q_index":int(np.argmax([r[f"S_q{i}_zz"].real/nc for i in range(nc)])),
                      "abs_delta_peak_q_index":int(np.argmax([abs(r[f"Delta_q{i}_z"]) for i in range(nc)]))})
    return rows,peaks

def save(path,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def series(rows,name):
    a=sorted((x for x in rows if x["observable"]==name),key=lambda x:x["tau_total"])
    return np.array([x["tau_total"] for x in a]),np.array([x["value"] for x in a]),np.array([x["stderr"] for x in a])

def decorate(ax,label):
    ax.set(xlabel=r"$\tau_{\mathrm{total}}$",ylabel=label)
    ax.set_xticks([.8,1.6,3.2,4.8,6.4,8.]);ax.grid(color="#d8dde3",linewidth=.7)
    ax.tick_params(direction="in",top=True,right=True)

def plots(rows,out):
    panels=[("m2",r"$m_z^2(q)=S_{zz}(q)/N_c$"),("rms",r"$m_{z,\mathrm{rms}}(q)$"),
            ("abs_delta",r"$|\langle\Delta_z(q)\rangle|$")]
    curves=[("q1",r"$q_1=(0,\pi)$","#0072B2","o"),("q3",r"$q_3=(\pi,\pi)$","#D55E00","s"),
            ("grid_max","all-q max (q=0 for all tau)","#009E73","^")]
    fig,axes=plt.subplots(1,3,figsize=(15,4.7),constrained_layout=True)
    for ax,(prefix,label) in zip(axes,panels):
        for suffix,legend,color,marker in curves:
            x,y,e=series(rows,f"{prefix}_{suffix}")
            ax.errorbar(x,y,yerr=e,marker=marker,ms=5,capsize=3,lw=1.4,color=color,label=legend)
        decorate(ax,label);ax.legend(frameon=False,fontsize=9)
        one,oa=plt.subplots(figsize=(6.7,4.7),constrained_layout=True)
        for suffix,legend,color,marker in curves:
            x,y,e=series(rows,f"{prefix}_{suffix}")
            oa.errorbar(x,y,yerr=e,marker=marker,ms=5,capsize=3,lw=1.5,color=color,label=legend)
        decorate(oa,label);oa.legend(frameon=False);oa.set_title("L=4 finite-size VBS proxy")
        one.savefig(out/f"vbs_{prefix}_vs_tau.png",dpi=220);plt.close(one)
    fig.suptitle("Finite-size VBS proxies vs projection length (two-seed joint jackknife)",fontsize=13)
    fig.savefig(out/"vbs_order_proxies_vs_tau.png",dpi=220);plt.close(fig)

def main():
    p=argparse.ArgumentParser();p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    rows,peaks=calculate(a.manifest);save(a.output/"vbs_order_proxies_by_tau.csv",rows)
    save(a.output/"vbs_peak_momenta_by_tau.csv",peaks);plots(rows,a.output)
    print("wrote 4 PNG and 2 CSV files",a.output)
if __name__=="__main__":main()
