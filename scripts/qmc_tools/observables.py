"""Existing ensemble estimators expressed by column name, never by raw offset."""
import re
import numpy as np


def transform_for(names, meta):
    available = set(names)
    nc = meta.get("Lx", 0)*meta.get("Ly", 0)
    grid = sorted({int(m[1]) for n in names if (m := re.fullmatch(r"S_q(\d+)_zz", n))})
    selected = sorted({int(m[1]) for n in names if (m := re.fullmatch(r"Ssel_s(\d+)_zz", n))})
    amp = sorted({(int(m[1]), m[2]) for n in names if (m := re.fullmatch(r"AmpA_q(\d+)_([xyz])", n))})
    chi = sorted({int(m[1]) for n in names if (m := re.fullmatch(r"ChiRaw_q(\d+)_zz", n))})
    if (grid or selected or amp or chi) and not nc:
        raise ValueError("Lx/Ly required for normalized derived observables; provide --metadata")
    if grid and grid != list(range(nc)):
        raise ValueError("incomplete all-momentum structure-factor schema")
    if grid and "shear" not in meta:
        raise ValueError("recorded shear required for real-space reconstruction")
    axes = "xyz"
    energy_names = ['energy_bare_matter', 'energy_matter_phonon', 'energy_phonon_kinetic',
                    'energy_phonon_potential', 'energy_external_source', 'energy_total']
    if available.intersection(energy_names) and not set(energy_names).issubset(available):
        raise ValueError('incomplete total-energy component schema')
    if available.intersection(energy_names) and meta.get('total_energy', {}).get('schema', 'phonon-link-total-energy-v1') != 'phonon-link-total-energy-v1':
        raise ValueError('unknown total-energy measurement convention')
    if 'gauge_u' in meta and (not nc or len(meta['gauge_u']) != 3*nc or any(u not in (-1, 1) for u in meta['gauge_u'])):
        raise ValueError('invalid recorded bond gauge signs')
    matrix_keys = lambda stem, i: [f"{stem}{i}_{a}{b}" for a in axes for b in axes]
    required = []
    for q in grid:
        required += matrix_keys("S_q", q) + [f"Delta_q{q}_{a}" for a in axes]
    for j in selected:
        required += matrix_keys("Ssel_s", j)+matrix_keys("Sprj_s", j)
    for q in chi:
        required += matrix_keys("ChiRaw_q", q)+[f"{stem}_q{q}_{a}" for stem in ["ChiD", "ChiIconj"] for a in axes]
    for q, a in amp:
        required += [f"{stem}_q{q}_{a}" for stem in ["AmpA", "AmpA2", "AmpA4"]]
    missing = set(required)-available
    if missing:
        raise ValueError("incomplete derived observable inputs: " + ",".join(sorted(missing)[:5]))
    inverse = None
    if grid:
        lx, ly, shear = meta["Lx"], meta["Ly"], meta["shear"]
        coords = np.array([divmod(d, ly) for d in range(nc)])
        q = np.array([(n1/lx, n2/ly-shear*n1/(lx*ly)) for n1 in range(lx) for n2 in range(ly)])
        inverse = np.exp(-2j*np.pi*(coords @ q.T))/nc
    participants = meta.get("vbs", {}).get("candidate_set", {}).get("participant_indices")
    if selected and participants is None:
        opts = meta.get("recorded_options", {})
        path = opts.get("vbs-q-file")
        # Raw selected observables and their eigenvalues remain recoverable;
        # candidate membership must come from recorded metadata, never a guess.
        participants = None

    def eigen(out, matrix, suffix, projected=False):
        herm = (matrix+matrix.conj().T)*.5
        stems = ["TrS", "mz2", "leading_eigenvalue", "hermitian_residual"]
        if projected:
            stems = [s+"_proj" for s in stems]
        values = [np.trace(matrix), matrix[2, 2]/nc, np.linalg.eigvalsh(herm)[-1], np.linalg.norm(matrix-matrix.conj().T)]
        out.update((s+suffix, z) for s, z in zip(stems, values))

    def derive(r):
        out = {}
        if grid:
            matrices = np.array([[r[k] for k in matrix_keys("S_q", q)] for q in grid]).reshape(nc, 3, 3)
            for q, matrix in enumerate(matrices):
                eigen(out, matrix, f"_q{q}")
                out[f"m2_z_q{q}"] = matrix[2, 2]/nc
            c = (inverse @ matrices.reshape(nc, 9)).reshape(nc, 3, 3)
            background = np.outer([r[f"Delta_q0_{a}"] for a in axes], [r[f"Delta_q0_{a}"] for a in axes])
            for d in range(nc):
                for i, a in enumerate(axes):
                    for j, b in enumerate(axes):
                        out[f"C_q{d}_{a}{b}"] = c[d, i, j]
                        out[f"C_connected_q{d}_{a}{b}"] = c[d, i, j]-background[i, j]
        for j in selected:
            for stem, proj in [("Ssel_s", False), ("Sprj_s", True)]:
                s = np.array([r[k] for k in matrix_keys(stem, j)]).reshape(3, 3)
                eigen(out, s, f"_s{j}", proj)
        if selected and participants is not None:
            if not set(participants).issubset(selected):
                raise ValueError("invalid recorded candidate membership")
            winner = max(participants, key=lambda j: r[f"Ssel_s{j}_zz"].real) if participants else -1
            out.update(candidate_max_mz2=r[f"Ssel_s{winner}_zz"].real/nc if participants else 0.,
                       candidate_max_winner_index=winner, candidate_participant_count=len(participants))
        for q in chi:
            for a in axes:
                for b in axes:
                    out[f"chi_q{q}_{a}{b}"] = r[f"ChiRaw_q{q}_{a}{b}"]-nc*r[f"ChiD_q{q}_{a}"]*r[f"ChiIconj_q{q}_{b}"]
            for stem, part in [("Cos", "Re"), ("Sin", "Im")]:
                key = f"Chi{stem}Raw_q{q}_zz"
                if key in r:
                    out[f"chi_{stem.lower()}_q{q}_zz"] = r[key]-nc*r[f"ChiD{part}_q{q}_z"]*r[f"ChiI{part}_q{q}_z"]
        for q, a in amp:
            mean, a2, a4 = (r[f"{stem}_q{q}_{a}"] for stem in ["AmpA", "AmpA2", "AmpA4"])
            out[f"VA_q{q}_{a}"] = nc*(a2-mean*mean)
            out[f"B4_q{q}_{a}"] = a4/(a2*a2) if a2 != 0 else complex(np.nan, np.nan)
            if f"AmpAIA_q{q}_{a}" in r:
                out[f"chiA_q{q}_{a}"] = nc*(r[f"AmpAIA_q{q}_{a}"]-mean*r[f"AmpIA_q{q}_{a}"])
        if "Field_m" in r:
            m, score, ms = (r[k] for k in ["Field_m", "Field_score", "Field_m_score"])
            a2, a4 = r["Field_absD2"], r["Field_absD4"]
            out.update(m_parallel=m, chi_parallel=ms-m*score, absD2=a2,
                       B4_complex=a4/(a2*a2) if a2 != 0 else complex(np.nan, np.nan))
        if "X_x" in r and "X_y" in r:
            out["Xplus"] = (r["X_x"]+r["X_y"])/2
            out["Xminus"] = (r["X_x"]-r["X_y"])/2
        if set(energy_names).issubset(r):
            out['energy_matter_link'] = r['energy_bare_matter']+r['energy_matter_phonon']
            out['energy_phonon_link'] = r['energy_phonon_kinetic']+r['energy_phonon_potential']
            out['energy_sum_residual'] = r['energy_total']-sum(r[n] for n in energy_names[:-1])
            if nc:
                for name in energy_names:
                    out[name+'_per_cell'] = r[name]/nc
        # Existing center-node X2 remains useful, but is a DIFFERENT time cut
        # from the new link energy. Never add it to link kinetic/matter energy.
        if nc and 'k' in meta and all('X2_'+a in r for a in axes):
            out['energy_phonon_potential_center_reconstructed'] = nc*meta['k']/2*sum(r['X2_'+a] for a in axes)
        gauge = meta.get('gauge_u')
        if nc and gauge is not None and len(gauge) == 3*nc and all(k in meta for k in ['Jx', 'Jy', 'Jz', 'g']):
            if all(f'BondK_b{b}' in r and f'BondXK_b{b}' in r for b in range(3*nc)):
                js = [meta['Jx'], meta['Jy'], meta['Jz']]
                out['energy_bare_matter_center_reconstructed'] = sum(gauge[b]*js[b%3]*r[f'BondK_b{b}'] for b in range(3*nc))
                out['energy_matter_phonon_center_reconstructed'] = meta['g']*sum(gauge[b]*r[f'BondXK_b{b}'] for b in range(3*nc))
        return out

    derive.omissions = (["selected candidate maxima: participant metadata unavailable"]
                        if selected and participants is None else [])
    return derive
