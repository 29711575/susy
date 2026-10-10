#!/usr/bin/env python3
"""Read P5 joint-JK output; overlay a first-order diagnostic, never refit data.

No new statistical estimator or modification of the raw MC observables.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def response(tau, h, g, mass, omega, bond_mean):
    z=omega*h/2
    return -g*np.asarray(bond_mean)/(mass*omega**2)*(-np.expm1(-omega*tau/2))*z/np.tanh(z)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis',type=Path)
    parser.add_argument('free',type=Path)
    args=parser.parse_args()
    groups=json.loads((args.analysis/'groups.json').read_text())
    free=json.loads(args.free.read_text())
    B=np.array([complex(*z) for z in free['trial_bond_mean']])
    if np.max(np.abs(B.imag))>1e-12 or free['trial_to_projected_green_max_error']>1e-12:
        raise ValueError('free baseline diagnostic failed')
    B=B.real
    rows=[]
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),constrained_layout=True)
    for label,group in groups.items():
        p=group['physics']; m=p['intervals']; h=p['dt']; tau=m*h
        if any(p.get('trial_shift',[0,0,0])):
            raise ValueError('this old diagnostic is for undisplaced endpoints only')
        if any(p[k]!=free[k] for k in ('Lx','Ly','Jx','Jy','Jz','wilson_loops')):
            raise ValueError('free baseline geometry/gauge/couplings differ')
        prediction=response(tau,h,p['g'],p['mass'],p['omega_ph'],B)
        j=np.arange(m+1); w=np.full(m+1,h); w[[0,-1]]*=.5
        Cmid=np.exp(-p['omega_ph']*h*abs(j-m//2))/(2*p['mass']*p['omega_ph'])
        direct=-p['g']*B*np.dot(Cmid,w)
        error=float(np.max(abs(prediction-direct)))
        if error>1e-14:
            raise ValueError('closed form and discrete Cw disagree')
        nc=p['Lx']*p['Ly']
        for ax,name,pred in zip(axes,('X_z','uniform_background_trace'),(prediction[2],nc*np.dot(prediction,prediction))):
            e=group['estimates'][name]
            ax.errorbar(tau,e['value'][0],yerr=1.96*e['stderr_real'],fmt='o',capsize=4,label=label)
            rows.append(dict(group=label,tau=tau,h=h,observable=name,actual=e['value'][0],stderr=e['stderr_real'],
                first_order_diagnostic=pred,actual_minus_diagnostic=e['value'][0]-pred,discrete_Cw_error=error))
    p=next(iter(groups.values()))['physics']
    taus=np.linspace(0,max(row['tau'] for row in rows),200)
    curve=np.array([response(tau,p['dt'],p['g'],p['mass'],p['omega_ph'],B) for tau in taus])
    for ax,values,name in zip(axes,(curve[:,2],p['Lx']*p['Ly']*np.sum(curve**2,axis=1)),('X_z','Btrace')):
        ax.plot(taus,values,'--',color='grey',label='Free-baseline response (not exact)')
        ax.set(xlabel='Total open-contour tau',ylabel=name)
        ax.legend(fontsize=8)
    fig.suptitle(f"Actual MC {p['Lx']}x{p['Ly']}, J=(1,1,.5), g={p['g']}, m={p['mass']}, omega={p['omega_ph']}, h={p['dt']}\n"
        "Undisplaced vacuum; 95% marginal joint-JK bars; short tau=6.4 pilot, not convergence proof",fontsize=10)
    fig.savefig(args.analysis/'projection_response.png',dpi=160)
    plt.close(fig)
    with (args.analysis/'projection_response.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print('Actual free trial B_mu:',B,'max discrete Cw error:',max(r['discrete_Cw_error'] for r in rows))


if __name__=='__main__':
    main()
