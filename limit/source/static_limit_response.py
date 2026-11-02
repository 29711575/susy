#!/usr/bin/env python3
"""Relax every bond in the strict static limit and measure source response.

Fixed zero-flux, Wilson(+,+) matter sector; no phonon kinetic term.
Unlike static_vbs_response.py, the background is a finite-torus stationary
configuration of all bond variables, not a prescribed uniform BZ solution.
Multistart local minimization is not a global optimality certificate.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

import numpy as np
import scipy
from scipy.optimize import minimize

from static_vbs_response import bond_table, hopping_matrix, mode_basis, static_hessian
from predict_vbs_momenta import allowed_grid


class StaticModel:
    def __init__(self, lx, ly, shear, j, g, spring_k, parity="trial"):
        self.lx, self.ly, self.shear = lx, ly, shear
        self.nc, self.nb = lx*ly, 3*lx*ly
        self.j, self.g, self.k, self.parity = np.asarray(j), g, spring_k, parity
        bonds = bond_table(lx, ly, shear)
        self.a = np.array([x*ly+y for x,y,mu,bx,by,u in bonds])
        self.b = np.array([bx*ly+by for x,y,mu,bx,by,u in bonds])
        self.mu = np.array([mu for x,y,mu,bx,by,u in bonds])
        self.coords = np.array([(x,y) for x in range(lx) for y in range(ly)])
        self.trial_sign = int(np.linalg.slogdet(self.matrix(np.zeros(self.nb)))[0])
        if self.trial_sign == 0:
            raise ValueError("Bare matter trial has a zero mode; parity is not fixed by det T")

    def matrix(self, x):
        t = np.zeros((self.nc, self.nc))
        np.add.at(t, (self.a, self.b), self.j[self.mu]+self.g*np.asarray(x))
        return t

    def spectral(self, x):
        t = self.matrix(x)
        u, s, vh = np.linalg.svd(t, full_matrices=False)
        sign = int(np.linalg.slogdet(t)[0])
        correction = self.parity == "trial" and sign != self.trial_sign
        return u, s, vh, sign, correction

    def objective(self, x, force=None):
        u, s, vh, sign, correction = self.spectral(x)
        polar = u@vh
        matter = -s.sum()
        if correction:
            matter += 2*s[-1]
            polar -= 2*np.outer(u[:, -1], vh[-1, :])
        force = np.zeros(self.nb) if force is None else force
        energy = matter+self.k/2*np.dot(x, x)-np.dot(force, x)
        return float(energy), self.k*x-self.g*polar[self.a, self.b]-force

    def hessian(self, x):
        u, s, vh, sign, correction = self.spectral(x)
        if s[-1] < 1e-10*max(1., s[0]):
            raise ValueError("Matter zero mode: smooth static Hessian cannot be certified")
        if correction:
            # A fixed-parity excited branch also differentiates its smallest
            # singular value. Do not use the ground-state polar formula here.
            gap = s[-2]-s[-1]
            if gap < 1e-8*max(1., s[0]):
                raise ValueError("Degenerate parity excitation: smooth response undefined")
            eps = 1e-5
            eye = np.eye(self.nb)
            h = np.column_stack([(self.objective(x+eps*v)[1]-
                                  self.objective(x-eps*v)[1])/(2*eps) for v in eye])
            residual = float(np.max(abs(h-h.T)))
            if residual > 2e-6:
                raise ValueError("Fixed-parity Hessian failed symmetry check")
            return (h+h.T)/2, "central gradient differences on fixed-parity branch", residual
        i, j = np.triu_indices(self.nc, 1)
        v = vh.T
        # d^2(-||T||_*) = -sum_(i<j) (A_ij-A_ji)^2/(s_i+s_j).
        c = self.g*(u[self.a[:, None], i]*v[self.b[:, None], j]-
                    u[self.a[:, None], j]*v[self.b[:, None], i])/np.sqrt(s[i]+s[j])
        h = self.k*np.eye(self.nb)-c@c.T
        return h, "analytic derivative of nuclear-norm polar factor", float(np.max(abs(h-h.T)))

    def relax(self, x, force=None, maxiter=2500):
        # |X_b| <= (g+|force_b|)/k follows from |polar_ab| <= 1.
        f = np.zeros(self.nb) if force is None else force
        bound = (self.g+abs(f))/self.k+1e-8
        return minimize(lambda y: self.objective(y, f), np.clip(x, -bound, bound),
                        jac=True, method="L-BFGS-B", bounds=list(zip(-bound, bound)),
                        options={"gtol": 2e-9, "ftol": 2e-15, "maxiter": maxiter,
                                 "maxls": 50, "maxcor": 30})

    def polish(self, x, force=None):
        x = x.copy()
        for _ in range(12):
            energy, grad = self.objective(x, force)
            if np.max(abs(grad)) < 2e-11:
                break
            h, _, _ = self.hessian(x)
            if np.linalg.eigvalsh(h)[0] <= 1e-8:
                break
            step = np.linalg.solve(h, -grad)
            for exponent in range(22):
                trial = x+step*2.**(-exponent)
                new_energy, new_grad = self.objective(trial, force)
                if new_energy <= energy+2e-12 and np.linalg.norm(new_grad) < np.linalg.norm(grad):
                    x = trial
                    break
            else:
                break
        return x


def initial_states(model, count, seed):
    rng = np.random.default_rng(seed)
    scale = model.g/model.k
    starts = [("zero", np.zeros(model.nb)),
              ("uniform_equal", np.full(model.nb, scale/np.sqrt(3)))]
    for mu in range(3):
        x = np.zeros((model.nc, 3)); x[:, mu] = scale
        starts.append((f"uniform_bond_{mu}", x.ravel()))
    base = model.relax(np.zeros(model.nb)).x
    starts.append(("finite_uniform_stationary", base))
    # Fixed-q production candidate, Kekule candidate and several long/short
    # wavelength patterns; these are only initial guesses, never constraints.
    indices = [(5,4), (3,6), (1,0), (0,1), (1,1), (4,4), (4,0), (2,7)]
    grid = allowed_grid(model.lx, model.ly, model.shear)
    unique = set()
    for n1, n2 in indices:
        n1, n2 = n1 % model.lx, n2 % model.ly
        if (n1,n2) in unique or (n1,n2) == (0,0):
            continue
        unique.add((n1,n2))
        q = grid[n1*model.ly+n2]["fraction"]
        phase = 2*np.pi*(model.coords@np.array(q))
        for mu in range(3):
            x = base.reshape(model.nc,3).copy()
            x[:, mu] += .65*scale*np.cos(phase)
            starts.append((f"cos_n{n1}_{n2}_bond{mu}", x.ravel()))
    while len(starts) < count:
        i = len(starts)
        amplitude = (.15, .5, 1., 2.)[i % 4]*scale
        starts.append((f"random_{i}_amp{amplitude:.5g}",
                       base+rng.normal(0., amplitude, model.nb)))
    return starts[:count]


def response_grid(model, x, hinv):
    rows, matrices = [], []
    fields = x.reshape(model.nc,3)
    for point in allowed_grid(model.lx, model.ly, model.shear):
        q = np.array(point["fraction"])
        phase = np.exp(-2j*np.pi*(model.coords@q))
        f = np.zeros((3, model.nb), complex)
        for mu in range(3):
            f[mu, mu::3] = phase
        chi = f@hinv@f.conj().T/model.nc
        d = phase@fields/model.nc
        eigenvalues = np.linalg.eigvalsh(chi)
        rows.append({"index": point["index"], "n1": point["n1"], "n2": point["n2"],
                     "u": float(q[0]), "v": float(q[1]), "chi_xx": float(chi[0,0].real),
                     "chi_yy": float(chi[1,1].real), "chi_zz": float(chi[2,2].real),
                     "chi_eigenvalue_max": float(eigenvalues[-1]),
                     "D_z_abs2": float(abs(d[2])**2),
                     "chi_hermitian_residual": float(np.max(abs(chi-chi.conj().T)))})
        matrices.append({**point, "chi_real": chi.real.tolist(), "chi_imag": chi.imag.tolist(),
                         "eigenvalues": eigenvalues.tolist(),
                         "D_real": d.real.tolist(), "D_imag": d.imag.tolist()})
    return rows, matrices


def validation(model, x, h, hinv, peak, seed, *,
               field_strengths=(.0001, .00005, .000025, .000005, .0000025)):
    if not field_strengths or any(not np.isfinite(v) or v <= 0 for v in field_strengths):
        raise ValueError("Validation field strengths must be finite and positive")
    finest_field = min(field_strengths)
    rng = np.random.default_rng(seed+100)
    directions = rng.normal(size=(4, model.nb))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    gradient_checks = []
    for d in directions:
        e, grad = model.objective(x)
        ep, gp = model.objective(x+1e-5*d)
        em, gm = model.objective(x-1e-5*d)
        gradient_checks.append({"gradient_direction_error": float(abs((ep-em)/2e-5-grad@d)),
                                "hessian_vector_error": float(np.max(abs((gp-gm)/2e-5-h@d)))})
    # Independent energy normalization uses the inherited real-space helper
    # and the 2Nc Majorana matrix, distinct from the polar derivative formula.
    t = model.matrix(x)
    a = np.block([[np.zeros_like(t), 2*t], [-2*t.T, np.zeros_like(t)]])
    unrestricted = float(-.25*np.sum(abs(np.linalg.eigvalsh(1j*a))))
    u, s, vh, sign, correction = model.spectral(x)
    norm_error = abs(unrestricted+s.sum())
    matrix_error = float(np.max(abs(t-hopping_matrix(model.lx, model.ly, model.shear,
                                      model.j, x.reshape(model.nc,3), model.g))))
    points = {(peak["n1"],peak["n2"]), (5 % model.lx,4 % model.ly)}
    field_checks = []
    for n1, n2 in sorted(points):
        point = allowed_grid(model.lx,model.ly,model.shear)[n1*model.ly+n2]
        phase = 2*np.pi*(model.coords@np.array(point["fraction"]))
        for quadrature, values in (("cos", np.cos(phase)), ("minus_sin", -np.sin(phase))):
            source = np.zeros(model.nb); source[2::3] = values
            expected = float(source@hinv@source/model.nc)
            # Weak lattice pinning in a modulated state can make h=1e-3
            # nonlinear. Resolve the derivative rather than loosen its test.
            for strength in field_strengths:
                dx = hinv@source
                plus = model.polish(model.relax(x+strength*dx, strength*source).x, strength*source)
                minus = model.polish(model.relax(x-strength*dx, -strength*source).x, -strength*source)
                measured = float(source@(plus-minus)/(2*strength*model.nc))
                field_checks.append({"n1":n1, "n2":n2, "quadrature":quadrature,
                                     "h":strength, "chi_hessian":expected,
                                     "chi_finite_field":measured,
                                     "absolute_error":abs(measured-expected),
                                     "gradient_plus_max":float(np.max(abs(model.objective(plus,strength*source)[1]))),
                                     "gradient_minus_max":float(np.max(abs(model.objective(minus,-strength*source)[1])))})
    uniform_crosschecks = []
    if np.max(abs(x.reshape(model.nc,3)-x.reshape(model.nc,3).mean(axis=0))) < 1e-8 and not correction:
        jeff = model.j+model.g*x.reshape(model.nc,3).mean(axis=0)
        for n1,n2 in sorted(points):
            q = allowed_grid(model.lx,model.ly,model.shear)[n1*model.ly+n2]["fraction"]
            result = static_hessian(model.lx,model.ly,model.shear,jeff,model.g,q,spring_k=model.k)
            if result["total_hessian"] is not None:
                modes = np.array([item[2].ravel() for item in mode_basis(model.lx,model.ly,model.shear,q)])
                projected = modes@h@modes.T/model.nc
                uniform_crosschecks.append({"n1":n1,"n2":n2,
                    "occupied_band_hessian_error":float(np.max(abs(projected-result["total_hessian"])))})
    free = StaticModel(model.lx,model.ly,model.shear,model.j,0.,model.k,model.parity)
    hf,_,_ = free.hessian(np.zeros(model.nb))
    free_rows,_ = response_grid(free,np.zeros(model.nb),np.linalg.inv(hf))
    free_error = max(abs(row["chi_zz"]-1/model.k) for row in free_rows)
    residual = float(np.max(abs(h@hinv-np.eye(model.nb))))
    last_field_errors = [r["absolute_error"] for r in field_checks if r["h"] == finest_field]
    passed = (norm_error < 1e-9 and matrix_error < 1e-12 and free_error < 1e-11 and
              residual < 1e-8 and max(r["hessian_vector_error"] for r in gradient_checks) < 2e-6 and
              max(r["gradient_direction_error"] for r in gradient_checks) < 2e-7 and
              max(last_field_errors, default=0.) < 2e-5 and
              max((max(r["gradient_plus_max"],r["gradient_minus_max"])
                   for r in field_checks), default=0.) < 1e-9 and
              max((r["occupied_band_hessian_error"] for r in uniform_crosschecks), default=0.) < 1e-8)
    return {"status":"PASS" if passed else "FAIL", "majorana_energy_normalization_error":float(norm_error),
            "inherited_hopping_matrix_error":matrix_error, "inverse_residual":residual,
            "free_chi_1_over_k_error":float(free_error), "directional_checks":gradient_checks,
            "finest_field":finest_field, "finest_field_max_absolute_error":max(last_field_errors,default=0.),
            "finite_field_checks":field_checks, "uniform_occupied_band_crosschecks":uniform_crosschecks}


def write_csv(path, rows):
    with path.open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def plot_response(rows, output, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lx,ly = max(r["n1"] for r in rows)+1,max(r["n2"] for r in rows)+1
    zz = np.array([r["chi_zz"] for r in rows]).reshape(lx,ly)
    ev = np.array([r["chi_eigenvalue_max"] for r in rows]).reshape(lx,ly)
    fig,axes = plt.subplots(1,2,figsize=(10.5,4.4),constrained_layout=True)
    for ax,z,label in zip(axes,(zz,ev),(r"$\chi_{zz}(q)$",r"$\lambda_{\max}[\chi(q)]$")):
        im = ax.imshow(z.T, origin="lower", extent=(-.5,lx-.5,-.5,ly-.5),cmap="viridis")
        fig.colorbar(im,ax=ax,label=label)
        ax.set(xlabel=r"$n_1$",ylabel=r"$n_2$",xticks=range(lx),yticks=range(ly),title=label)
        peak = max((r for r in rows if r["index"]),key=lambda r:r["chi_zz"] if z is zz else r["chi_eigenvalue_max"])
        ax.plot(peak["n1"],peak["n2"],"o",mfc="none",mec="white",ms=14,mew=1.8)
    fig.suptitle(title+"; Gamma included, nonzero-q maximum circled")
    fig.savefig(output,dpi=180); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--Lx",type=int,default=9); ap.add_argument("--Ly",type=int,default=9)
    ap.add_argument("--shear",type=int,default=0)
    ap.add_argument("--Jx",type=float,default=1.); ap.add_argument("--Jy",type=float,default=1.)
    ap.add_argument("--Jz",type=float,default=.5)
    ap.add_argument("--lambda",dest="lam",type=float,required=True)
    ap.add_argument("--spring-k",type=float,default=1.)
    ap.add_argument("--starts",type=int,default=80); ap.add_argument("--seed",type=int,default=106101700)
    ap.add_argument("--parity",choices=("trial","unrestricted"),default="trial")
    ap.add_argument("--output-dir",type=Path,required=True); ap.add_argument("--plot",action="store_true")
    args = ap.parse_args()
    if args.lam < 0 or args.spring_k <= 0 or args.starts < 6 or not all(np.isfinite(v) for v in
            (args.lam,args.spring_k,args.Jx,args.Jy,args.Jz)):
        ap.error("require finite lambda>=0, spring-k>0, couplings and starts>=6")
    out = args.output_dir
    if out.exists() and any(out.iterdir()):
        ap.error("output-dir must be new or empty; historical runs are preserved")
    out.mkdir(parents=True,exist_ok=True)
    start = time.monotonic()
    model = StaticModel(args.Lx,args.Ly,args.shear,[args.Jx,args.Jy,args.Jz],
                        np.sqrt(args.lam*args.spring_k),args.spring_k,args.parity)
    records, xs, best = [], [], None
    for i,(label,x0) in enumerate(initial_states(model,args.starts,args.seed)):
        r = model.relax(x0)
        e,grad = model.objective(r.x)
        xs.append(r.x)
        records.append({"start":i,"label":label,"energy":e,"gradient_max":float(np.max(abs(grad))),
                        "optimizer_success":bool(r.success),"iterations":int(r.nit),"evaluations":int(r.nfev),
                        "message":str(r.message),"det_T_sign":model.spectral(r.x)[3]})
        if best is None or e < best[0]:
            best = (e,r.x.copy(),i)
        write_csv(out/"multistart.csv",records)
        if (i+1)%10 == 0 or i+1 == args.starts:
            print(json.dumps({"completed_starts":i+1,"best_energy":best[0],
                              "elapsed_seconds":time.monotonic()-start}),flush=True)
    np.savez_compressed(out/"multistart_configurations.npz",x=np.array(xs))
    x = model.polish(best[1]); e,grad = model.objective(x)
    h,method,symmetry = model.hessian(x); eigs = np.linalg.eigvalsh(h)
    np.save(out/"hessian.npy",h)
    if eigs[0] <= 1e-8 or np.max(abs(grad)) > 2e-7:
        (out/"failure.json").write_text(json.dumps({"status":"unstable_or_unconverged",
            "minimum_hessian_eigenvalue":float(eigs[0]),"gradient_max":float(np.max(abs(grad)))},indent=2)+"\n")
        raise RuntimeError("Best-found configuration cannot support a certified local response")
    hinv = np.linalg.inv(h)
    rows,matrices = response_grid(model,x,hinv)
    write_csv(out/"susceptibility.csv",rows)
    bonds = [{"x":int(c[0]),"y":int(c[1]),"Xx":float(v[0]),"Xy":float(v[1]),"Xz":float(v[2])}
             for c,v in zip(model.coords,x.reshape(model.nc,3))]
    write_csv(out/"equilibrium_bonds.csv",bonds)
    peak = max((r for r in rows if r["index"]),key=lambda r:r["chi_zz"])
    mixed_peak = max((r for r in rows if r["index"]),key=lambda r:r["chi_eigenvalue_max"])
    checks = validation(model,x,h,hinv,peak,args.seed)
    (out/"validation.json").write_text(json.dumps(checks,indent=2,allow_nan=False)+"\n")
    u,s,vh,sign,correction = model.spectral(x)
    matter = float(-s.sum()+(2*s[-1] if correction else 0.))
    field = x.reshape(model.nc,3)
    peak_point = matrices[peak["index"]]
    peak_d = complex(peak_point["D_real"][2],peak_point["D_imag"][2])
    # Exact translations are degenerate static minima. Nonzero order makes
    # the field-minimized h->0 response distinct from a single-branch derivative.
    translation_d = peak_d*np.exp(-2j*np.pi*(model.coords@np.array([peak["u"],peak["v"]])))
    domains = {"translations":model.nc,"peak_D_z_abs":float(abs(peak_d)),
        "peak_cos_D_range":[float(translation_d.real.min()),float(translation_d.real.max())],
        "peak_minus_sin_D_range":[float(translation_d.imag.min()),float(translation_d.imag.max())],
        "zero_field_derivative_scope":"If these ranges are nonzero, switching among degenerate translations gives different h->0+ and h->0- order values. A unique finite symmetry-restored zero-field susceptibility is not supplied by the local Hessian.",
        "translation_D_real":translation_d.real.tolist(),"translation_D_imag":translation_d.imag.tolist()}
    result = {"schema":"strict-static-limit-response-v1", "argv":sys.argv,
        "parameters":{"Lx":args.Lx,"Ly":args.Ly,"shear":args.shear,"Nc":model.nc,"Nb":model.nb,
                      "J":model.j.tolist(),"lambda":args.lam,"g":model.g,"k":model.k,
                      "limit":"mass -> infinity at fixed k and g; phonon kinetic term absent",
                      "zero_flux":True,"wilson_sector":[1,1],"parity":args.parity},
        "trial_det_T_sign":model.trial_sign,"equilibrium_det_T_sign":sign,
        "parity_excitation_applied":bool(correction),"energy_total":e,"energy_per_cell":e/model.nc,
        "matter_energy":matter,"elastic_energy":float(model.k/2*x@x),
        "minimum_singular_value":float(s[-1]),"gradient_max":float(np.max(abs(grad))),
        "mean_X":field.mean(axis=0).tolist(),"spatial_std_X":field.std(axis=0).tolist(),
        "mean_X2":(field**2).mean(axis=0).tolist(),
        "hessian_method":method,"hessian_symmetry_residual":symmetry,
        "minimum_hessian_eigenvalue":float(eigs[0]),"hessian_eigenvalues":eigs.tolist(),
        "multistart_count":len(records),"best_start_index":best[2],
        "starts_within_1e_minus_7":sum(abs(r["energy"]-e)<1e-7 for r in records),
        "starts_optimizer_success":sum(r["optimizer_success"] for r in records),
        "starts_gradient_below_1e_minus_6":sum(r["gradient_max"]<1e-6 for r in records),
        "nonzero_q_chi_zz_peak":peak,"nonzero_q_mixed_chi_peak":mixed_peak,
        "nonzero_q_chi_zz_peak_ties":[r for r in rows if r["index"] and abs(r["chi_zz"]-peak["chi_zz"])<1e-9],
        "fixed_production_q":next((r for r in rows if r["n1"]==5 and r["n2"]==4),None),
        "translated_domains":domains,
        "normalization":"D_mu(q)=sum_r exp(-iq.r) X_mu(r)/Nc; chi_mu,nu(q)=F_mu H_X^-1 F_nu^*/Nc; source H_h=H-h sum_r f(r)X_nu(r)",
        "response_scope":"Local branch source derivative after relaxing all bonds. For nonuniform domains, same-q diagonal response; no domain switching or finite-temperature averaging.",
        "limitation":"Multistart local search, no mathematical global-minimum certificate; fixed-gauge matter, no full physical-spin projection. No QMC sampling, finite-tau estimate, or free-phonon subtraction.",
        "validation_status":checks["status"],"software":{"numpy":np.__version__,"scipy":scipy.__version__},
        "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                         (Path(__file__),Path(__file__).with_name("static_vbs_response.py"),
                          Path(__file__).with_name("predict_vbs_momenta.py"))},
        "all_q":matrices,"elapsed_seconds":time.monotonic()-start}
    (out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if args.plot:
        plot_response(rows,out/"susceptibility.png",f"Strict static, L={args.Lx}x{args.Ly}, lambda={args.lam:g}, k={args.spring_k:g}")
    print(json.dumps({k:v for k,v in result.items() if k not in ("all_q","hessian_eigenvalues")},indent=2),flush=True)
    if checks["status"] != "PASS":
        raise RuntimeError("Response validation failed; inspect validation.json")


if __name__ == "__main__":
    main()
