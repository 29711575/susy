"""Joint common-phase ratios; delete blocks within each independent chain."""
import numpy as np


def evaluate(total, count, names, transform):
    if abs(total[0]) <= 1e-12*count:
        raise ValueError("unresolved complex phase denominator")
    ratios = total[1:]/total[0]
    out = {"phase": total[0]/count}
    out.update(zip(names, ratios))
    out.update(transform(dict(zip(names, ratios))))
    return out


def joint(chains, transform=lambda r: {}):
    total = sum((c.sums.sum(0) for c in chains))
    count = sum(int(c.counts.sum()) for c in chains)
    result = evaluate(total, count, chains[0].names, transform)
    labels = list(result)
    values = np.array(list(result.values()), dtype=complex)
    valid = np.isfinite(values)
    variance = np.zeros((3, len(values)))
    reasons = []
    for c in chains:
        n = len(c.counts)
        if n < 2:
            reasons.append("each chain needs at least two groups for a variance estimate")
            break
        chain_sum = c.sums.sum(0)
        chain_count = int(c.counts.sum())
        # Restore the deleted chain's original allocation. Unequal tail groups
        # use weighted pseudovalues, as in analyze_g0_validation.py. Equal
        # groups reduce exactly to the original stratified block jackknife.
        mean = np.zeros(len(values), complex)
        m2r = np.zeros(len(values)); m2i = m2r.copy(); cov = m2r.copy()
        weight_total = 0.
        try:
            for block, block_count in zip(c.sums, c.counts):
                weight = float(block_count)/chain_count
                replica = evaluate(total+(weight*chain_sum-block)/(1-weight), count, c.names, transform)
                if list(replica) != labels:
                    raise ValueError("derived schema changed within jackknife")
                x = np.array(list(replica.values()), dtype=complex)
                valid &= np.isfinite(x)
                x = np.where(np.isfinite(x), x, 0.)
                x = np.where(np.isfinite(values), values, 0.)+(1-weight)/weight*(np.where(np.isfinite(values), values, 0.)-x)
                weight_total += weight
                delta = x-mean; mean += weight/weight_total*delta; after = x-mean
                m2r += weight*delta.real*after.real; m2i += weight*delta.imag*after.imag
                cov += weight*delta.real*after.imag
        except (ValueError, ZeroDivisionError, FloatingPointError) as error:
            reasons.append(str(error)); break
        variance += np.array([m2r, m2i, cov])/(n-1)
    estimates = {}
    for i, label in enumerate(labels):
        estimates[label] = dict(real=float(values[i].real), imag=float(values[i].imag),
                               error_real=None if reasons or not valid[i] else float(np.sqrt(max(0., variance[0, i]))),
                               error_imag=None if reasons or not valid[i] else float(np.sqrt(max(0., variance[1, i]))),
                               covariance_real_imag=None if reasons or not valid[i] else float(variance[2, i]),
                               defined=bool(valid[i]))
    phase = estimates["phase"]
    resolved = not reasons and abs(values[0]) > 5*np.hypot(phase["error_real"], phase["error_imag"])
    return dict(samples=count, chain_blocks=[len(c.counts) for c in chains],
                chain_group_sample_counts=[c.counts.tolist() for c in chains],
                unequal_tail_groups=any(len(set(c.counts)) != 1 for c in chains),
                too_few_groups=any(len(c.counts) < 16 for c in chains),
                average_phase=phase, phase_resolved_5sigma=bool(resolved),
                ratio_defined=True, jackknife_defined=not reasons,
                limitations=reasons, undefined_observables=[n for i, n in enumerate(labels) if not valid[i]], estimates=estimates)
