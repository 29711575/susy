"""Descriptive rank/folded split Rhat on explicitly labeled stored blocks."""
from statistics import NormalDist
import numpy as np


def rank_normalize(x):
    flat = x.ravel()
    order = np.argsort(flat, kind="stable")
    ranked = np.empty(len(flat), float)
    start = 0
    while start < len(flat):
        end = start+1
        while end < len(flat) and flat[order[end]] == flat[order[start]]:
            end += 1
        ranked[order[start:end]] = (start+end-1)/2+1
        start = end
    normal = NormalDist()
    return np.array([normal.inv_cdf((r-.375)/(len(flat)+.25)) for r in ranked]).reshape(x.shape)


def basic_rhat(x):
    n = x.shape[1]
    within = np.mean(np.var(x, axis=1, ddof=1))
    between = n*np.var(np.mean(x, axis=1), ddof=1)
    if within == 0:
        return 1. if between == 0 else float("inf")
    return float(np.sqrt(((n-1)/n*within+between/n)/within))


def bulk_ess(x):
    m, n = x.shape
    centered = x-x.mean(1, keepdims=True)
    fft = np.fft.rfft(centered, n=1 << (2*n-1).bit_length(), axis=1)
    acov = np.fft.irfft(fft*fft.conj(), axis=1)[:, :n]/n
    within = np.var(x, axis=1, ddof=1).mean()
    var_plus = (n-1)/n*within+np.var(x.mean(1), ddof=1)
    if var_plus == 0:
        return None
    rho = 1-(within-acov.mean(0))/var_plus
    rho[0] = 1.
    pairs = []
    for lag in range(0, n-1, 2):
        pair = rho[lag]+rho[lag+1]
        if pair < 0:
            break
        pairs.append(min(pair, pairs[-1]) if pairs else pair)
    tau = max(1., -1+2*sum(pairs))
    return float(min(m*n, m*n/tau))


def rank_diagnostic(values):
    x = np.asarray(values, float)
    if x.shape[0] < 2 or x.shape[1] < 8 or not np.isfinite(x).all():
        return dict(status="insufficient_data", rhat=None, ess_bulk_blocks=None)
    n = x.shape[1]//2
    split = np.concatenate([x[:, :n], x[:, -n:]], axis=0)
    z = rank_normalize(split)
    rhat = max(basic_rhat(z), basic_rhat(rank_normalize(np.abs(split-np.median(split)))))
    constant = bool(np.ptp(split) == 0)
    return dict(status="constant_series" if constant else "descriptive_only", rhat=None if constant else rhat,
                ess_bulk_blocks=None if constant or rhat > 1.01 else bulk_ess(z),
                ess_withheld_for_chain_disagreement=bool(rhat > 1.01),
                block_records_per_chain=x.shape[1], split_records_used_per_chain=2*n,
                equilibration_certified=False)


def diagnose(chains, names=None):
    common = min(len(c.counts) for c in chains)
    sizes = {int(c.counts[0])*c.stride for c in chains}
    if len(sizes) != 1:
        raise ValueError("diagnostic chains need matching stored block durations")
    default = [n for n in ["X_x", "X_y", "X_z", "majorana_energy", "energy_total", "energy_phonon_kinetic", "energy_phonon_potential"] if n in chains[0].names]
    if {'X_x', 'X_y'}.issubset(chains[0].names):
        default += ['Xminus']
    meta = chains[0].metadata
    entries = meta.get('q_mode_pool_entries') or [
        [e['n1'], e['n2']] for e in meta.get('q_mode_pool', {}).get('entries', [])]
    opts = meta.get('recorded_options', {})
    if not entries and 'q-mode-n1' in opts and 'q-mode-n2' in opts:
        entries = [[int(opts['q-mode-n1']), int(opts['q-mode-n2'])]]
    q = entries[0][0]*meta['Ly']+entries[0][1] if entries and 'Ly' in meta else None
    if q is not None:
        default += [n for n in [f'S_q{q}_zz', f'Delta_q{q}_z'] if n in chains[0].names]
    selected = names or default
    if not selected:
        raise ValueError("no default diagnostic observables; specify --observable")
    lookup = {n:i+1 for i, n in enumerate(chains[0].names)}
    total = sum(c.sums[:common].sum(0) for c in chains)
    count = sum(int(c.counts[:common].sum()) for c in chains)
    if abs(total[0]) <= 1e-12*count:
        raise ValueError("phase denominator unresolved; influence diagnostic undefined")
    phase_mean = total[0]/count
    rows, traces = [], []
    for name in selected:
        if name not in lookup and name not in {"Xplus", "Xminus"}:
            raise ValueError("unknown raw diagnostic observable: " + name)
        if name in {"Xplus", "Xminus"} and not {"X_x", "X_y"}.issubset(lookup):
            raise ValueError("Xplus/Xminus require X_x and X_y")
        def numerator(a):
            if name in {"Xplus", "Xminus"}:
                return (a[..., lookup["X_x"]]+(1 if name == "Xplus" else -1)*a[..., lookup["X_y"]])/(2 if name == "Xplus" else 1)
            return a[..., lookup[name]]
        theta = numerator(total)/total[0]
        # Influence series retains ONE common ensemble centering; never
        # subtract each chain's own mean, which would erase disagreement.
        series = [(numerator(c.sums[:common])-theta*c.sums[:common, 0])/(c.counts[:common]*phase_mean) for c in chains]
        for component in ["real", "imag"]:
            values = np.array([getattr(a, component) for a in series])
            rows.append(dict(observable=name, component=component, **rank_diagnostic(values)))
        for chain, values in zip(chains, series):
            for end, value in zip(chain.ends[:common], values):
                traces.append(dict(seed=chain.metadata.get("seed"), path=str(chain.path), sweep=int(end),
                                   observable=name, influence_real=float(value.real), influence_imag=float(value.imag)))
    return dict(common_block_records=common, block_sweeps=next(iter(sizes)),
                windows=[dict(path=str(c.path), start_sweep=c.start_sweep, end_sweep=int(c.ends[common-1]),
                              excluded_tail_blocks=len(c.counts)-common) for c in chains],
                default_fixed_q_index=q,
                series="joint complex phase-ratio influence; stored block means, not per-sweep samples",
                equilibration_certified=False, diagnostics=rows), traces
