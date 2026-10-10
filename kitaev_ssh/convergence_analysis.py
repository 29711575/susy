#!/usr/bin/env python3
"""P5: chain-preserving extension of the P3/P4 joint block jackknife.

No new MC estimator: read existing s and sO blocks. Reblock ONLY inside a
chain; do not merge different distributions. For independent chains of
unequal length, restore the deleted chain's original weight in each replica
and sum its jackknife covariance contribution. For one chain this is exactly
the P3/P4 delete-one-block ratio estimator, including complex covariance.
"""
import argparse
import csv
from dataclasses import dataclass, replace
import json
from pathlib import Path
from statistics import NormalDist

import numpy as np
from plot_vbs import load_run, write_csv, normalize_trial_metadata
import matplotlib.pyplot as plt

PRIMARY = ["Szz_q1", "TrS_q1", "Szz_q3", "TrS_q3"]
STATIONARITY = ["X_z", "uniform_background_trace"] + PRIMARY
EPSILON = dict(zip(PRIMARY, [.02, .04, .02, .04]))
PHYSICS = ("Lx", "Ly", "Jx", "Jy", "Jz", "g", "mass", "omega_ph", "scope",
           "spatial_boundary", "wilson_loops", "time_contour")
DISTRIBUTION = PHYSICS + ("intervals", "dt", "trial_shift")

def check_trial_comparison(a,b,step=False):
    normalize_trial_metadata(a); normalize_trial_metadata(b)
    if a['dt']!=b['dt'] and not np.isclose(a['intervals']*a['dt'],b['intervals']*b['dt'],rtol=0,atol=1e-12):
        raise ValueError('changing projection length and step together confounds comparison')
    if a['trial_shift']!=b['trial_shift']:
        if step or a['intervals']!=b['intervals'] or a['dt']!=b['dt']:
            raise ValueError('changing trial shift with tau/h confounds comparison')


@dataclass
class Chain:
    path: Path
    metadata: dict
    names: list
    sums: np.ndarray  # columns: sum(s), then every existing sum(sO)
    counts: np.ndarray
    csv_ratio_error: float = 0.

    @property
    def seed(self):
        return self.metadata["seed"]

    @property
    def block_sweeps(self):
        return int(self.counts[0])*self.metadata["measure_every"]

    def reblock(self, sweeps):
        if sweeps < self.block_sweeps or sweeps % self.block_sweeps:
            raise ValueError("block duration must be a multiple of each chain's stored block")
        factor = sweeps // self.block_sweeps
        if len(self.counts) % factor:
            raise ValueError("incomplete chain block: no remainder silently discarded")
        return replace(self, sums=self.sums.reshape(-1, factor, self.sums.shape[1]).sum(axis=1),
                       counts=self.counts.reshape(-1, factor).sum(axis=1))


def read_chain(path):
    path = Path(path)
    meta, difference = load_run(path)  # reuse P4 JSON/CSV consistency checks
    names = list(meta["observables"])
    with path.with_name(path.name.removesuffix("_results.json")+"_blocks.csv").open() as stream:
        blocks = list(csv.DictReader(stream))
    if [int(b["block"]) for b in blocks] != list(range(len(blocks))):
        raise ValueError("unordered or missing blocks")
    columns = ["sum_phase"] + ["sum_phase_"+name for name in names]
    sums = np.array([[complex(float(b[c+"_real"]), float(b[c+"_imag"])) for c in columns] for b in blocks])
    counts = np.array([int(b["count"]) for b in blocks])
    if not np.isfinite(sums).all():
        raise ValueError("nonfinite raw block data")
    if not np.isclose(meta["intervals"]*meta["dt"], meta["tau"], rtol=0, atol=1e-12):
        raise ValueError("inconsistent contour metadata")
    return Chain(path, meta, names, sums, counts, difference)


def checked_group(chains):
    for c in chains:
        normalize_trial_metadata(c.metadata)
    if not chains or len({c.seed for c in chains}) != len(chains):
        raise ValueError("independent chains must have distinct seeds")
    if len({str(c.path.resolve()) for c in chains}) != len(chains):
        raise ValueError("duplicated chain file")
    for c in chains:
        for key in DISTRIBUTION:
            if c.metadata[key] != chains[0].metadata[key]:
                raise ValueError(f"cannot pool different distributions: {key}")
        if c.names != chains[0].names:
            raise ValueError("observable schemas differ")
        if c.metadata.get('vbs',{}).get('momenta') != chains[0].metadata.get('vbs',{}).get('momenta'):
            raise ValueError('momentum grid ordering differs between chains')
        if not len(c.counts) or len(set(c.counts)) != 1:
            raise ValueError("each chain requires equal nonempty blocks")
    # Matter trial remains freeGroundState(X=0); shifted phonon endpoints are
    # part of the distribution and must never be pooled across distinct d.
    return chains


def fixed_momentum_mapping(metadata):
    """q1/q3 are historical channel aliases, not general flattened indices."""
    lx, ly = metadata['Lx'], metadata['Ly']
    if lx < 2 or ly < 2 or lx % 2 or ly % 2:
        raise ValueError('fixed half-reciprocal momenta require even Lx and Ly')
    nc = lx*ly
    grid = metadata.get('vbs', {}).get('momenta', [
        {'index': x*ly+y, 'n': [x,y], 'reciprocal_fraction': [x/lx,y/ly]}
        for x in range(lx) for y in range(ly)])
    if metadata.get('vbs', {}).get('Nc', nc) != nc or len(grid) != nc:
        raise ValueError('inconsistent Nc or momentum grid')
    if sorted(m['index'] for m in grid) != list(range(nc)) or len({tuple(m['n']) for m in grid}) != nc:
        raise ValueError('duplicated or missing momentum indices')
    for m in grid:
        x,y=m['n']
        if not (0<=x<lx and 0<=y<ly) or not np.allclose(m['reciprocal_fraction'],[x/lx,y/ly],rtol=0,atol=1e-12):
            raise ValueError('inconsistent reciprocal coordinates')
    result={}
    for alias, target in [('q0',[0,0]), ('q1',[0,ly//2]), ('q3',[lx//2,ly//2])]:
        result[alias]=next(m for m in grid if m['n']==target)
    return result


def fixed_observables(sums, count, names, metadata):
    if abs(sums[0]) <= 1e-12*count:
        raise ValueError("numerically unresolved phase denominator")
    r = dict(zip(names, sums[1:]/sums[0]))
    mapping=fixed_momentum_mapping(metadata)
    nc=metadata['Lx']*metadata['Ly']
    out = {"phase": sums[0]/count}
    for alias, momentum in mapping.items():
        q=momentum['index']
        out[f"Szz_{alias}"] = r[f"S_q{q}_zz"]
        out[f"TrS_{alias}"] = sum(r[f"S_q{q}_{mu}{mu}"] for mu in "xyz")
    for mu in "xyz":
        out[f"X_{mu}"] = r[f"Delta_q{mapping['q0']['index']}_{mu}"]
    out["X2_z"] = r["X2_z"]
    out["uniform_background_zz"] = nc*out["X_z"]**2
    out["uniform_background_trace"] = nc*sum(out[f"X_{mu}"]**2 for mu in "xyz")
    out["Szz_connected_q0"] = out["Szz_q0"]-out["uniform_background_zz"]
    out["TrS_connected_q0"] = out["TrS_q0"]-out["uniform_background_trace"]
    other1=[q for q in range(nc) if q not in (mapping['q0']['index'],mapping['q1']['index'])]
    other3=[q for q in range(nc) if q not in (mapping['q0']['index'],mapping['q3']['index'])]
    out["TrS_q1_minus_other_nonzero_mean"] = out["TrS_q1"]-sum(r[f"S_q{q}_{mu}{mu}"] for q in other1 for mu in "xyz")/len(other1)
    out["Szz_q3_minus_other_nonzero_mean"] = out["Szz_q3"]-sum(r[f"S_q{q}_zz"] for q in other3)/len(other3)
    return out


def full_grid_observables(sums, count, names, metadata):
    if abs(sums[0]) <= 1e-12*count:
        raise ValueError('numerically unresolved phase denominator')
    r=dict(zip(names,sums[1:]/sums[0]))
    nc=metadata['Lx']*metadata['Ly']
    return {'phase':sums[0]/count, **{f'S_q{q}_{mu}{nu}':r[f'S_q{q}_{mu}{nu}']
            for q in range(nc) for mu in 'xyz' for nu in 'xyz'}}


def joint_estimate(chains, transform=fixed_observables):
    checked_group(chains)
    sums = sum((c.sums.sum(axis=0) for c in chains))
    count = sum(int(c.counts.sum()) for c in chains)
    result = transform(sums, count, chains[0].names, chains[0].metadata)
    labels = list(result)
    values = np.array(list(result.values()))
    covariance = np.zeros((2*len(labels), 2*len(labels)))
    phase_covariance = np.zeros((2,2))
    components, lags = [], []
    for c in chains:
        b = len(c.counts)
        if b < 2:
            raise ValueError("each independent chain needs at least two blocks for joint jackknife")
        chain_sum = c.sums.sum(axis=0)
        replicas, phase_replicas = [], []
        for block in c.sums:
            # Other chains stay fixed; reweight this chain after deleting one
            # equal block. Total sample count/chain allocation stays unchanged.
            replica = sums+(chain_sum-b*block)/(b-1)
            phase_replicas.append(replica[0]/count)
            transformed = transform(replica, count, c.names, c.metadata)
            replicas.append([transformed[label] for label in labels])
        replicas = np.array(replicas)
        centered = replicas-replicas.mean(axis=0)
        real = np.concatenate((centered.real, centered.imag), axis=1)
        component = (b-1)/b * (real.T @ real)
        covariance += component
        pp=np.array(phase_replicas);pp-=pp.mean()
        pp=np.column_stack((pp.real,pp.imag))
        phase_covariance+=(b-1)/b*(pp.T@pp)
        components.append(component)
        lag = []
        for j in range(len(labels)):
            a = centered[:, j].real
            denominator = np.dot(a, a)
            lag.append(float(np.dot(a[:-1], a[1:])/denominator) if denominator > 1e-28 else None)
        lags.append(lag)
    if not np.isfinite(values).all() or not np.isfinite(covariance).all():
        raise ValueError("nonfinite joint estimate")
    phase=sums[0]/count
    phase_resolved=bool(abs(phase)>5*np.sqrt(np.trace(phase_covariance)))
    return {"names": labels, "values": values, "covariance": covariance, "chain_covariances": components,
            "average_phase": phase, "phase_covariance": phase_covariance, "phase_resolved_5sigma": phase_resolved,
            "chain_lag1": lags, "samples": count, "chain_blocks": [len(c.counts) for c in chains],
            "chain_samples": [int(c.counts.sum()) for c in chains], "seeds": [c.seed for c in chains]}


def estimates_json(estimate):
    n = len(estimate["names"])
    out = {}
    for i, name in enumerate(estimate["names"]):
        z = estimate["values"][i]
        c = estimate["covariance"]
        out[name] = {"value": [z.real, z.imag], "stderr_real": float(np.sqrt(c[i, i])),
                     "stderr_imag": float(np.sqrt(c[n+i, n+i])), "cov_real_imag": float(c[i, n+i])}
    return out


def comparison(a, b, label, primary=PRIMARY, epsilon=EPSILON, looks=1):
    if set(a["seeds"]) & set(b["seeds"]):
        raise ValueError("difference SE requires independent seed sets")
    if a["names"] != b["names"]:
        raise ValueError("estimate schemas differ")
    rows = []
    for i, name in enumerate(a["names"]):
        d = float((b["values"][i]-a["values"][i]).real)
        se = float(np.sqrt(a["covariance"][i, i]+b["covariance"][i, i]))
        eps = epsilon.get(name)
        family = len(primary) if name in primary else 1
        critical = NormalDist().inv_cdf(1-.05/(2*family))
        look_critical = NormalDist().inv_cdf(1-.05/(2*family*looks))
        low, high = d-1.96*se, d+1.96*se
        slo, shi = d-critical*se, d+critical*se
        llo,lhi=d-look_critical*se,d+look_critical*se
        rows.append({"comparison": label, "observable": name,
                     "a": float(a["values"][i].real), "b": float(b["values"][i].real),
                     "difference_b_minus_a": d, "absolute_difference": abs(d), "difference_se": se,
                     "z": d/se if se > 1e-16 else None, "ci95_low": low, "ci95_high": high,
                     "simultaneous95_low": slo, "simultaneous95_high": shi, "family_size": family,
                     "look_safeguard95_low": llo, "look_safeguard95_high": lhi, "looks": looks,
                     "epsilon": eps, "marginal_interval_inside_epsilon": eps is not None and low > -eps and high < eps,
                     "simultaneous_interval_inside_epsilon": eps is not None and slo > -eps and shi < eps,
                     "look_safeguard_inside_epsilon": eps is not None and llo > -eps and lhi < eps,
                     "phase_resolved_5sigma": a['phase_resolved_5sigma'] and b['phase_resolved_5sigma'],
                     "min_blocks_per_chain": min(a['chain_blocks']+b['chain_blocks']),
                     "interval_method": "normal approximation; simultaneous uses Bonferroni; conditional on adequate blocking"})
    return rows


def second_order_residual(rows):
    """If Oh=O0+c*h^2+O(h^4), fine residual = (Oh-Ofine)/3.

    NOT a fit and NOT a check that this model/step is in the asymptotic regime.
    """
    out=[]
    for r in rows:
        out.append({"observable":r['observable'],"fine_minus_continuum_conditional":-r['difference_b_minus_a']/3,
                    "stderr_conditional":r['difference_se']/3,
                    "ci95_low_conditional":-r['ci95_high']/3,"ci95_high_conditional":-r['ci95_low']/3,
                    "look_safeguard95_low_conditional":-r['look_safeguard95_high']/3,
                    "look_safeguard95_high_conditional":-r['look_safeguard95_low']/3,
                    "block_look_envelope_low_conditional":-r.get('blocking_look_envelope_high',r['look_safeguard95_high'])/3,
                    "block_look_envelope_high_conditional":-r.get('blocking_look_envelope_low',r['look_safeguard95_low'])/3,
                    "assumption":"second-order asymptotic regime; NOT established by two steps alone"})
    return out


def add_block_sensitivity(rows, a, b):
    # Diagnostic envelope, NOT a rigorous bound on all possible long-time tails.
    # Different samplers may require different block durations. Keep only
    # durations leaving >=16 blocks in EACH chain and take each group's max SE.
    for r in rows:
        name=r['observable']
        errors=[]
        for group in (a,b):
            eligible=[level['estimates'][name]['stderr_real'] for level in group['blocking'] if min(level['chain_blocks'])>=16]
            errors.append(max([group['estimates'][name]['stderr_real']]+eligible))
        se=float(np.hypot(*errors));d=r['difference_b_minus_a']
        critical=NormalDist().inv_cdf(1-.05/(2*r['family_size']*r['looks']))
        r.update(max_block_sensitivity_se=se, blocking_envelope95_low=d-1.96*se,
                 blocking_envelope95_high=d+1.96*se,
                 blocking_look_envelope_low=d-critical*se,blocking_look_envelope_high=d+critical*se,
                 blocking_look_inside_epsilon=r['epsilon'] is not None and abs(d)+critical*se<r['epsilon'])
    return rows


def metadata_row(c):
    d = c.metadata
    row = {"path": str(c.path), **{key: d[key] for key in ("Lx", "Ly", "Jx", "Jy", "Jz", "g", "lambda", "mass", "omega_ph",
            "intervals", "dt", "proposal", "seed", "warmup_sweeps", "measurement_sweeps", "measure_every", "samples", "block_samples", "blocks", "acceptance", "elapsed_seconds")}}
    row.update(tau=d["intervals"]*d["dt"], stored_block_sweeps=c.block_sweeps,
               phase_real=d["average_phase"]["value"][0], phase_imag=d["average_phase"]["value"][1],
               phase_stderr_real=d["average_phase"]["stderr_real"], phase_stderr_imag=d["average_phase"]["stderr_imag"],
               csv_ratio_error=c.csv_ratio_error,
               trial_shift=d['trial_shift'],
               trial="free matter ground state at X=0; equal shifted phonon endpoints; fixed zero flux")
    return row


def historical_audit(directory, output):
    all_chains = [read_chain(p) for p in sorted(Path(directory).glob("*_results.json"))
                  if json.loads(p.read_text()).get("vbs_enabled")]
    write_csv(output/"p4_metadata.csv", [metadata_row(c) for c in all_chains])
    by_name = {c.path.name.removesuffix("_results.json"): c for c in all_chains}
    rows = []
    for left, right, sweeps, quantity in (("seed_a", "half_step", 200, "S_q2_zz"),
            ("refresh_a", "refresh_half", 100, "S_q3_zz"), ("refresh_a", "refresh_long", 100, "TrS_q1")):
        def observed(s, count, names, metadata):
            r = dict(zip(names, s[1:]/s[0]))
            return {quantity: sum(r[f"S_q1_{mu}{mu}"] for mu in "xyz") if quantity == "TrS_q1" else r[quantity]}
        a, b = (joint_estimate([by_name[n].reblock(sweeps)], observed) for n in (left, right))
        row = comparison(a, b, right+" - "+left, primary=[quantity]*8, epsilon={})[0]
        q = int(quantity.split("q")[1][0])
        row.update(q_index=q, nx=q//2, ny=q%2, blocks_per_chain=16, block_sweeps=sweeps,
                   caveat="P4 maximum selected over momentum/channel comparisons; NOT a prespecified hypothesis")
        rows.append(row)
    write_csv(output/"p4_discrepancies.csv", rows)
    diagnostics = []
    for name in ("seed_a", "seed_b", "refresh_a", "refresh_b", "refresh_long", "half_step", "refresh_half"):
        c = by_name[name]
        for sweeps in (25, 50, 100, 200, 400):
            try:
                reb = c.reblock(sweeps)
                if len(reb.counts) < 8:
                    continue
                e = joint_estimate([reb])
            except ValueError:
                continue
            for quantity in PRIMARY:
                j = e["names"].index(quantity)
                var = float(e["covariance"][j, j])
                diagnostics.append({"chain": name, "seed": c.seed, "block_sweeps": sweeps, "blocks": len(reb.counts),
                    "observable": quantity, "mean": e["values"][j].real, "stderr": np.sqrt(var),
                    "lag1_of_block_influence": e["chain_lag1"][0][j], "variance_times_seconds": var*c.metadata["elapsed_seconds"]})
    write_csv(output/"p4_blocking.csv", diagnostics)


def group_analysis(label, chains, block_sweeps=100):
    checked_group(chains)
    primary = joint_estimate([c.reblock(block_sweeps) for c in chains])
    full_grid = joint_estimate([c.reblock(block_sweeps) for c in chains], full_grid_observables)
    metadata = chains[0].metadata
    summary = {"label": label, "physics": {k: metadata[k] for k in DISTRIBUTION},
               "fixed_momentum_aliases": fixed_momentum_mapping(metadata),
               "Nc": metadata['Lx']*metadata['Ly'],
               "full_momentum_grid": metadata.get('vbs',{}).get('momenta'),
               "full_momentum_data": "All complex S_mu,nu numerators retained in each source block CSV; no momentum discarded",
               "full_grid_estimates": estimates_json(full_grid),
               "normalization": "Nc=Lx*Ly; S=<s F_mu conj(F_nu)>/(Nc <s>); background=Nc*sum_mu(<X_mu>^2)",
               "phase_convention": "sample |Wf|Wph; retain complex s and sO",
               "primary_block_sweeps": block_sweeps, "chains": [metadata_row(c) for c in chains],
               "estimates": estimates_json(primary), "covariance_order": primary["names"]+['Im '+n for n in primary["names"]],
               "covariance": primary["covariance"].tolist(), "samples": primary["samples"],
               "chain_blocks": primary["chain_blocks"], "chain_samples": primary["chain_samples"],
               "phase_resolved_5sigma": primary['phase_resolved_5sigma'],
               "sum_chain_seconds": sum(c.metadata["elapsed_seconds"] for c in chains),
               "blocking": [], "independence_proved": False, "interval_method": "approximate normal, not an exact coverage guarantee"}
    for sweeps in sorted(set((block_sweeps, 5, 10, 20, 25, 50, 100, 200, 400))):
        try:
            rb = [c.reblock(sweeps) for c in chains]
            if min(len(c.counts) for c in rb) < 4:
                continue
            e = joint_estimate(rb)
        except ValueError:
            continue
        summary["blocking"].append({"block_sweeps": sweeps, "chain_blocks": e["chain_blocks"],
                                    "estimates": estimates_json(e), "chain_lag1": e["chain_lag1"]})
    chain_estimates = [joint_estimate([c.reblock(block_sweeps)]) for c in chains]
    chain_checks = []
    for i in range(len(chains)):
        for j in range(i+1, len(chains)):
            chain_checks.extend(comparison(chain_estimates[i], chain_estimates[j],
                                           f"{label}: seed {chains[j].seed}-{chains[i].seed}",
                                           primary=STATIONARITY))
    summary["between_chain_checks"] = chain_checks
    splits = []
    for c in chains:
        rb = c.reblock(block_sweeps)
        mid = len(rb.counts)//2
        if len(rb.counts)%2 or mid < 4:
            continue
        halves = [replace(rb, sums=rb.sums[sl], counts=rb.counts[sl]) for sl in (slice(0, mid), slice(mid, None))]
        e1, e2 = (joint_estimate([half]) for half in halves)
        for name in PRIMARY+["X_z", "X2_z", "uniform_background_trace"]:
            j = e1['names'].index(name)
            delta = float((e2['values'][j]-e1['values'][j]).real)
            se = float(np.sqrt(e1['covariance'][j, j]+e2['covariance'][j, j]))
            splits.append({"seed": c.seed, "observable": name, "second_half_minus_first": delta,
                           "difference_se": se, "z_diagnostic": delta/se if se>1e-16 else None,
                           "caveat": "stationarity screen only; neglects covariance across the single half-chain boundary"})
    summary["split_chain_diagnostics"] = splits
    # Four nonoverlapping measurement quarters per *individual* chain. Keep
    # complex phase/observable block sums and the joint ratio jackknife intact;
    # neither quarters nor chain boundaries are merged into artificial blocks.
    quarters, quarter_contrasts = [], []
    for c in chains:
        rb = c.reblock(block_sweeps)
        b = len(rb.counts)
        if b % 4 or b // 4 < 4:
            continue
        q = b // 4
        quarter_estimates = []
        for part in range(4):
            sl = slice(part*q, (part+1)*q)
            e = joint_estimate([replace(rb, sums=rb.sums[sl], counts=rb.counts[sl])])
            quarter_estimates.append(e)
            z = estimates_json(e)
            quarters.append({"seed": c.seed, "quarter": part+1,
                             "measurement_sweep_start": part*q*block_sweeps,
                             "measurement_sweep_end": (part+1)*q*block_sweeps,
                             "blocks": q, "samples": int(rb.counts[sl].sum()),
                             "estimates": {name: z[name] for name in STATIONARITY},
                             "caveat": "within-quarter delete-one-block joint jackknife; screening, not independent draws"})
        first, last = quarter_estimates[0], quarter_estimates[-1]
        for name in STATIONARITY:
            j = first['names'].index(name)
            delta = float((last['values'][j]-first['values'][j]).real)
            se = float(np.sqrt(first['covariance'][j,j]+last['covariance'][j,j]))
            quarter_contrasts.append({"seed": c.seed, "observable": name,
                                      "fourth_minus_first": delta, "difference_se": se,
                                      "z_diagnostic": delta/se if se>1e-16 else None,
                                      "caveat": "screen only; first/last quarter covariance not estimated across their shared chain"})
    summary["quarter_chain_diagnostics"] = quarters
    summary["quarter_first_last_diagnostics"] = quarter_contrasts
    return summary, primary


def plot_quarter_diagnostics(groups, output):
    for label, summary in groups.items():
        rows = summary["quarter_chain_diagnostics"]
        if not rows:
            continue
        fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
        for ax, name in zip(axes.flat, STATIONARITY):
            for seed in dict.fromkeys(r["seed"] for r in rows):
                parts = [r for r in rows if r["seed"] == seed]
                ax.errorbar([r["quarter"] for r in parts],
                            [r["estimates"][name]["value"][0] for r in parts],
                            yerr=[r["estimates"][name]["stderr_real"] for r in parts],
                            fmt="o-", capsize=3, label=str(seed))
            ax.set(xticks=(1, 2, 3, 4), title=name, xlabel="measurement quarter",
                   ylabel="mean +/- 1 joint-JK SE")
            ax.legend(title="seed", fontsize=8)
        p = summary["physics"]
        fig.suptitle(f"{label}: tau={p['dt']*p['intervals']:.3g}; shifted trial d={p['trial_shift']}; "
                     f"{summary['primary_block_sweeps']}-sweep block\n"
                     "Four consecutive quarters of each measured chain; screening only; no between-quarter covariance",
                     fontsize=10)
        fig.savefig(output/f"quarters_{label}.png", dpi=160, bbox_inches="tight")
        plt.close(fig)


def plot_groups(groups, comparisons, output):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    labels = list(groups)
    for ax, name in zip(axes.flat, PRIMARY):
        for i, label in enumerate(labels):
            e = groups[label]["estimates"][name]
            ax.errorbar(i, e["value"][0], yerr=1.96*e["stderr_real"], fmt="o", capsize=4)
        ax.set(xticks=range(len(labels)), xticklabels=labels, ylabel=name, title="95% marginal intervals (fixed q)")
        ax.tick_params(axis='x', labelrotation=20)
    p=next(iter(groups.values()))['physics']
    caption=f"Actual MC {p['Lx']}x{p['Ly']}, J=({p['Jx']},{p['Jy']},{p['Jz']}), g={p['g']}, m={p['mass']}, omega={p['omega_ph']}"
    caption+='\nEndpoint d by group: '+ '; '.join(x+'=('+','.join(f'{v:.5g}' for v in groups[x]['physics']['trial_shift'])+')' for x in labels)
    status=f"samples={[groups[x]['samples'] for x in labels]}, phase resolved={[groups[x]['phase_resolved_5sigma'] for x in labels]}"
    fig.suptitle(caption+"; independent-chain joint JK\n"+status+"\n"
                 "No momentum re-selection; intervals require adequate blocking; not a tau->infinity proof", fontsize=9)
    fig.savefig(output/"fixed_observables.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, name in zip(axes.flat, PRIMARY):
        for label, summary in groups.items():
            ax.plot([b["block_sweeps"] for b in summary["blocking"]],
                    [b["estimates"][name]["stderr_real"] for b in summary["blocking"]], 'o-', label=label)
        ax.set(xscale="log", xlabel="sweeps per block, inside each chain", ylabel="joint JK SE", title=name)
        ax.set_xticks([25, 50, 100, 200, 400], labels=['25', '50', '100', '200', '400'])
        ax.minorticks_off()
        ax.legend(fontsize=8)
    fig.suptitle(caption+"\nBlock-size stability; independent chain boundaries retained (no cross-chain blocks)", fontsize=10)
    fig.savefig(output/"blocking.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    for ax,name in zip(axes,('X_z','uniform_background_trace')):
        for i,label in enumerate(labels):
            e=groups[label]['estimates'][name]
            ax.errorbar(i,e['value'][0],yerr=1.96*e['stderr_real'],fmt='o',capsize=4)
        ax.set(xticks=range(len(labels)),xticklabels=labels,ylabel=name)
        ax.tick_params(axis='x',labelrotation=20)
    fig.suptitle(caption+"\nUniform background: 95% joint-JK intervals; not a modulation order parameter",fontsize=10)
    fig.savefig(output/'uniform_background.png',dpi=160,bbox_inches='tight')
    plt.close(fig)
    if comparisons:
        comparison_labels = list(dict.fromkeys(r["comparison"] for r in comparisons))
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        for ax, name in zip(axes.flat, PRIMARY):
            rows = [r for r in comparisons if r["observable"] == name]
            for j, r in enumerate(rows):
                d = r["difference_b_minus_a"]
                ax.errorbar(j, d, yerr=[[d-r['blocking_look_envelope_low']], [r['blocking_look_envelope_high']-d]], fmt='o', capsize=4)
            eps = EPSILON[name]
            ax.axhspan(-eps, eps, color="grey", alpha=.15, label=f"target +/-{eps}")
            ax.axhline(0, color="grey", lw=.7)
            ax.set(xticks=range(len(rows)), xticklabels=[r['comparison'] for r in rows], title=name, ylabel="B-A")
            ax.tick_params(axis='x', labelrotation=25)
            ax.legend(fontsize=8)
        fig.suptitle(caption+"\nAbsolute differences: block-sensitivity envelope, four-observable / look safeguard\n"
                     "Shaded region = prospective tolerance; conservative diagnostic, not an exact tail bound", fontsize=10)
        fig.savefig(output/"differences.png", dpi=160, bbox_inches="tight")
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", action="append", default=[], help="LABEL=path1,path2,...")
    parser.add_argument("--compare", action="append", default=[], help="A:B, report B-A")
    parser.add_argument("--block-sweeps", type=int, default=100)
    parser.add_argument("--group-block", action="append",default=[],help="LABEL=SWEEPS overrides the group's primary block size")
    parser.add_argument("--looks", type=int, default=1, help="additional Bonferroni look safeguard; not exact sequential coverage")
    parser.add_argument("--step-pair", help="coarse:fine labels, same tau and h ratio exactly 2; conditional residual /3")
    parser.add_argument("--audit-p4", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plots", action="store_true", help="explicitly request diagnostic images (default: numerical output only)")
    args = parser.parse_args()
    if args.looks<1:
        raise ValueError("looks must be positive")
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("output directory must be absent or empty")
    args.output.mkdir(parents=True, exist_ok=True)
    if args.audit_p4:
        historical_audit(args.audit_p4, args.output)
    groups, estimates, chains = {}, {}, {}
    overrides={key:int(val) for key,val in (arg.split('=') for arg in args.group_block)}
    for option in args.group:
        label, paths = option.split("=", 1)
        if label in groups:
            raise ValueError("duplicate group label")
        chains[label] = [read_chain(p) for p in paths.split(",")]
        groups[label], estimates[label] = group_analysis(label, chains[label], overrides.get(label,args.block_sweeps))
    comparisons = []
    for option in args.compare:
        left, right = option.split(":")
        a, b = chains[left][0].metadata, chains[right][0].metadata
        check_trial_comparison(a,b)
        if any(a[k] != b[k] for k in PHYSICS):
            raise ValueError("this benchmark comparison must hold physical/gauge parameters fixed")
        if a['proposal'] != b['proposal'] and (a['dt'] != b['dt'] or a['intervals'] != b['intervals']):
            raise ValueError("algorithm and discretization changed together: confounded comparison")
        rows=comparison(estimates[left], estimates[right], right+" - "+left,looks=args.looks)
        comparisons.extend(add_block_sensitivity(rows,groups[left],groups[right]))
    if args.step_pair:
        left,right=args.step_pair.split(':')
        a,b=chains[left][0].metadata,chains[right][0].metadata
        check_trial_comparison(a,b,step=True)
        if any(a[k]!=b[k] for k in PHYSICS) or a['proposal']!=b['proposal']:
            raise ValueError("step comparison must keep physical/gauge parameters and algorithm fixed")
        if not np.isclose(a['intervals']*a['dt'],b['intervals']*b['dt'],rtol=0,atol=1e-12) or not np.isclose(a['dt'],2*b['dt'],rtol=0,atol=1e-12):
            raise ValueError("step pair must have fixed total tau and h/2")
        step_rows=add_block_sensitivity(comparison(estimates[left],estimates[right],right+' - '+left,looks=args.looks),groups[left],groups[right])
        residual=second_order_residual(step_rows)
        write_csv(args.output/'conditional_step_residual.csv',residual)
    if groups:
        (args.output/"groups.json").write_text(json.dumps(groups, indent=2, allow_nan=False)+"\n")
        write_csv(args.output/"metadata.csv", [metadata_row(c) for cc in chains.values() for c in cc])
        quarter_rows = []
        half_rows, first_last_rows = [], []
        for label, summary in groups.items():
            half_rows.extend({"group": label, **r} for r in summary["split_chain_diagnostics"])
            first_last_rows.extend({"group": label, **r} for r in summary["quarter_first_last_diagnostics"])
            for row in summary["quarter_chain_diagnostics"]:
                for name, e in row["estimates"].items():
                    quarter_rows.append({"group": label, "seed": row["seed"], "quarter": row["quarter"],
                                         "measurement_sweep_start": row["measurement_sweep_start"],
                                         "measurement_sweep_end": row["measurement_sweep_end"],
                                         "blocks": row["blocks"], "samples": row["samples"],
                                         "observable": name, "mean_real": e["value"][0],
                                         "mean_imag": e["value"][1], "stderr_real": e["stderr_real"],
                                         "stderr_imag": e["stderr_imag"]})
        if quarter_rows:
            write_csv(args.output/"quarters.csv", quarter_rows)
        if half_rows:
            write_csv(args.output/"halves.csv", half_rows)
        if first_last_rows:
            write_csv(args.output/"quarter_first_last.csv", first_last_rows)
        if args.plots:
            plot_groups(groups, comparisons, args.output)
            plot_quarter_diagnostics(groups, args.output)
    if comparisons:
        write_csv(args.output/"comparisons.csv", comparisons)
        for row in comparisons:
            if row["observable"] in PRIMARY:
                print(row)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
