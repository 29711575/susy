"""Separate response adapter; the recorded frozen source is never modified.

Ordinary responses retain the 2e-5 absolute accuracy target.  A soft response
may instead be certified to 1e-6 relative precision, with its absolute error
budget explicitly retained.  Finite fields must converge on the same positive
curvature branch.  A comparison at one accidentally favourable field is not a
certificate.  High precision is used to measure/refine endpoint residuals when
the conservative double-precision residual floor would prevent certification.
"""
import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

from pathlib import Path
import sys
import json
import time
import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

TASK = Path(__file__).resolve().parent
sys.path.insert(0, str(TASK / "source"))
from physical_model import PhysicalStaticModel


ABSOLUTE_TARGET = 2e-5
SOFT_RELATIVE_TARGET = 1e-6


class ResponseCertificationError(RuntimeError):
    def __init__(self, report):
        self.report = report
        super().__init__("SOURCE_RESPONSE_NOT_CERTIFIED: " + report.get("reason", ""))


def signed_hessian(model, x):
    """Return (H, method, symmetry_error) for the fixed physical branch.

    For F=sum(d_i*s_i), A=U.T*dT*V, its second derivative is
      sum_i<j [(d_i-d_j)(Aij+Aji)^2/(2(si-sj))
              +(d_i+d_j)(Aij-Aji)^2/(2(si+sj))].
    Thus ++ pairs use skew/sqrt(si+sj); a +,- lowest pair uses
    symmetric/sqrt(si-sj).  The energy is elastic-F, hence H=kI-C*C.T.
    """
    u, s, vh, sign, correction = model.spectral(x)
    if s[-1] < 1e-10 * max(1., s[0]):
        raise ValueError("Matter zero mode: smooth static Hessian cannot be certified")
    if correction and s[-2] - s[-1] < 1e-8 * max(1., s[0]):
        raise ValueError("Degenerate parity excitation: smooth response undefined")
    i, j = np.triu_indices(model.nc, 1)
    opposite = (j == model.nc - 1) if correction else np.zeros(i.size, bool)
    denominator = np.where(opposite, s[i] - s[j], s[i] + s[j])
    v = vh.T
    numerator = (u[model.a[:, None], i] * v[model.b[:, None], j]
                 - np.where(opposite, -1., 1.)
                 * u[model.a[:, None], j] * v[model.b[:, None], i])
    c = model.g * numerator / np.sqrt(denominator)
    h = model.k * np.eye(model.nb) - c @ c.T
    symmetry_error = float(np.max(abs(h - h.T)))
    return h, "analytic signed-singular-value Hessian", symmetry_error


def newton_refine(model, x, force=None, *, maxiter=28, return_diagnostics=False):
    """Newton refinement using gradient merit, avoiding unresolved energy changes."""
    x = np.asarray(x, float).copy()
    f = np.zeros(model.nb) if force is None else np.asarray(force, float)
    best = None
    history = []
    reason = "iteration_limit"
    for iteration in range(maxiter):
        energy, gradient = model.objective(x, f)
        norm = float(np.linalg.norm(gradient))
        h, _, _ = signed_hessian(model, x)
        minimum = float(np.linalg.eigvalsh(h)[0])
        history.append(dict(iteration=iteration, gradient_max=float(abs(gradient).max()),
                            gradient_norm=norm, hessian_minimum=minimum))
        if minimum <= 1e-10:
            reason = "nonpositive_curvature"
            break
        if best is None or norm < best[0]:
            best = norm, x.copy(), gradient.copy(), h.copy()
        if abs(gradient).max() < 6e-16:
            reason = "floating_precision_reached"
            break
        step = np.linalg.solve(h, -gradient)
        for exponent in range(22):
            trial = x + step * 2.**(-exponent)
            trial_gradient = model.objective(trial, f)[1]
            if np.linalg.norm(trial_gradient) < norm:
                x = trial
                break
        else:
            reason = "gradient_precision_stagnation"
            break
    if best is None:
        raise ValueError("Newton initial point is outside a positive-curvature branch")
    info = dict(reason=reason, iterations=len(history), history=history,
                gradient_max=float(abs(best[2]).max()),
                hessian_minimum=float(np.linalg.eigvalsh(best[3])[0]))
    return (best[1], info) if return_diagnostics else best[1]


class AccuratePhysicalStaticModel(PhysicalStaticModel):
    def hessian(self, x):
        return signed_hessian(self, x)

    def polish(self, x, force=None):
        return newton_refine(self, x, force)


def _mp_gradient(model, values, source, field, digits):
    """Independent high-precision SVD gradient at an endpoint, no mp Hessian."""
    import mpmath as mp
    with mp.workdps(digits):
        xx = [mp.mpf(z) for z in values]
        ff = [mp.mpf(float(z)) * mp.mpf(float(field)) for z in source]
        gg, kk = mp.mpf(float(model.g)), mp.mpf(float(model.k))
        jj = [mp.mpf(float(z)) for z in model.j]
        t = mp.matrix(model.nc)
        for b in range(model.nb):
            t[int(model.a[b]), int(model.b[b])] += jj[int(model.mu[b])] + gg * xx[b]
        u, singular, vh = mp.svd_r(t)
        double_x = np.array([float(z) for z in xx])
        correction = model.spectral(double_x)[4]
        gradient = []
        for b in range(model.nb):
            a, bb = int(model.a[b]), int(model.b[b])
            polar = mp.fsum(u[a, i] * vh[i, bb]
                           * (-1 if correction and i == model.nc - 1 else 1)
                           for i in range(model.nc))
            gradient.append(kk * xx[b] - gg * polar - ff[b])
        floor = float(100 * mp.eps * (1 + abs(gg)) * mp.sqrt(model.nb))
    return xx, gradient, floor


def _endpoint(model, x, source, field, dx, *, high_precision, digits,
              precision_target, reference_correction):
    force = field * source
    predicted = x + field * dx
    y, newton = newton_refine(model, predicted, force, return_diagnostics=True)
    h, _, _ = signed_hessian(model, y)
    eigmin = float(np.linalg.eigvalsh(h)[0])
    gradient = model.objective(y, force)[1]
    ep_dx = np.linalg.solve(h, source)
    floor = 5 * np.finfo(float).eps * (1 + model.g) * np.sqrt(model.nb)
    mp_iterations = 0
    values = list(map(float, y))
    if high_precision:
        # A double Hessian is an excellent quasi-Newton preconditioner here.
        # Retain sub-double corrections in mp values; rounding each update back
        # to float would reintroduce the endpoint precision floor.
        import mpmath as mp
        for mp_iterations in range(1, 4):
            values, mp_gradient, floor = _mp_gradient(model, values, source, field, digits)
            gradient = np.array([float(z) for z in mp_gradient])
            y = np.array([float(z) for z in values])
            h, _, _ = signed_hessian(model, y)
            eigmin = float(np.linalg.eigvalsh(h)[0])
            ep_dx = np.linalg.solve(h, source)
            projected = abs(float(ep_dx @ gradient))
            conservative = float(np.linalg.norm(ep_dx) * np.linalg.norm(gradient))
            endpoint_allowance = precision_target * abs(field) * model.nc / 64
            if max(projected, conservative) < endpoint_allowance:
                break
            if mp_iterations == 3:
                # Every returned observable must have a residual evaluated at
                # exactly those values, including an unsuccessful final step.
                break
            correction = np.linalg.solve(h, -gradient)
            with mp.workdps(digits):
                values = [z + mp.mpf(float(dz)) for z, dz in zip(values, correction)]
    projected = abs(float(ep_dx @ gradient))
    rounding = float(np.linalg.norm(ep_dx) * floor)
    # Factor two reserves room for local nonlinear residual correction and
    # the double preconditioner's rounding.  Diagnostics also expose the norm
    # bound and Newton correction so this remains inspectable.
    root_numerator_bound = 2 * (projected + rounding)
    correction_norm = float(np.linalg.norm(np.linalg.solve(h, gradient)))
    displacement = float(np.linalg.norm(y - x))
    predicted_displacement = float(abs(field) * np.linalg.norm(dx))
    ratio = displacement / max(predicted_displacement, np.finfo(float).tiny)
    correction = bool(model.spectral(y)[4])
    same_branch = (eigmin > 1e-10 and correction == reference_correction
                   and ratio < 4. and correction_norm < max(1e-6, displacement * .01))
    with __import__('mpmath').workdps(digits) if high_precision else _null_context():
        if high_precision:
            import mpmath as mp
            observable = mp.fsum(mp.mpf(float(a)) * b for a, b in zip(source, values))
        else:
            observable = float(source @ y)
    return dict(field=float(field), observable=observable,
                gradient_max=float(abs(gradient).max()), hessian_minimum=eigmin,
                root_numerator_bound=root_numerator_bound,
                projected_gradient=float(projected), gradient_roundoff_floor_l2=float(floor),
                norm_gradient_projection_bound=float(np.linalg.norm(ep_dx) * np.linalg.norm(gradient)),
                newton_correction_norm=correction_norm, displacement=displacement,
                displacement_over_linear=ratio, same_branch=bool(same_branch),
                precision="mp%d" % digits if high_precision else "float64",
                mp_gradient_evaluations=mp_iterations, newton=newton)


class _null_context:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def certify_source(model, x, source, expected=None, *, absolute_target=ABSOLUTE_TARGET,
                   relative_target=SOFT_RELATIVE_TARGET, max_fields=14,
                   minimum_field=1e-12, high_precision_fallback=True, mp_digits=40,
                   raise_on_fail=True, start_h_factor=1., **unused):
    """Certify one real source derivative, returning a JSON-serializable record.

    A field sequence first runs in double precision.  When necessary it is
    re-evaluated with mp endpoint residuals/refinement.  Consecutive Richardson
    windows must converge; the entire difference between overlapping windows
    is charged to the truncation budget (rather than divided by an assumed
    asymptotic factor).  Endpoint residual and rounding bounds are propagated
    through the local inverse Hessian and the extrapolation coefficients.
    """
    start = time.time()
    source = np.asarray(source, float)
    x = newton_refine(model, x)
    h, _, _ = signed_hessian(model, x)
    minimum = float(np.linalg.eigvalsh(h)[0])
    if minimum <= 1e-10:
        raise ValueError("Cannot certify response outside a positive-curvature branch")
    dx = np.linalg.solve(h, source)
    analytic = float(source @ dx / model.nc)
    expected = analytic if expected is None else float(expected)
    if abs(expected - analytic) > max(1e-7, relative_target * .1 * abs(analytic)):
        raise ValueError("Expected response does not use the analytic Hessian/refined background")
    target = max(float(absolute_target), float(relative_target) * abs(expected))
    reference_correction = bool(model.spectral(x)[4])
    # Limit predictor motion before entering the nonlinear soft-mode region.
    initial = min(1e-5, .0005 * (1 + float(abs(x).max())) / max(float(abs(dx).max()), 1.))
    initial = max(initial * float(start_h_factor), minimum_field * 4)
    attempts = []
    windows = []
    best = None
    for precise in ([False, True] if high_precision_fallback else [False]):
        fields = []
        previous_window = None
        for index in range(max_fields):
            strength = initial * 2.**(-index)
            if strength < minimum_field:
                break
            try:
                plus = _endpoint(model, x, source, strength, dx, high_precision=precise,
                                 digits=mp_digits, precision_target=target,
                                 reference_correction=reference_correction)
                minus = _endpoint(model, x, source, -strength, dx, high_precision=precise,
                                  digits=mp_digits, precision_target=target,
                                  reference_correction=reference_correction)
                if precise:
                    import mpmath as mp
                    with mp.workdps(mp_digits):
                        measured = float((plus['observable'] - minus['observable'])
                                         / (2 * mp.mpf(float(strength)) * model.nc))
                else:
                    measured = (plus['observable'] - minus['observable']) / (2 * strength * model.nc)
                bound = (plus['root_numerator_bound'] + minus['root_numerator_bound']) / (2 * strength * model.nc)
                observable_rounding = (0. if precise else
                                      4 * np.finfo(float).eps * np.linalg.norm(source)
                                      * np.linalg.norm(x) / (strength * model.nc))
                bound += observable_rounding
                for endpoint in (plus, minus):
                    endpoint['observable'] = str(endpoint['observable']) if precise else endpoint['observable']
                record = dict(h=float(strength), expected=expected, measured=float(measured),
                              error=float(abs(measured - expected)), root_and_roundoff_uncertainty=float(bound),
                              observable_roundoff_uncertainty=float(observable_rounding),
                              same_branch=plus['same_branch'] and minus['same_branch'],
                              plus=plus, minus=minus)
                fields.append(record)
                attempts.append(record)
            except (ValueError, np.linalg.LinAlgError) as err:
                attempts.append(dict(h=float(strength), precision="mp" if precise else "float64", error_message=str(err)))
                continue
            if len(fields) < 3:
                continue
            a, b, c = fields[-3:]
            if not (abs(a['h'] / b['h'] - 2) < 1e-10 and abs(b['h'] / c['h'] - 2) < 1e-10):
                continue
            coarse = (4 * b['measured'] - a['measured']) / 3
            fine = (4 * c['measured'] - b['measured']) / 3
            coarse_bound = (4 * b['root_and_roundoff_uncertainty'] + a['root_and_roundoff_uncertainty']) / 3
            fine_bound = (4 * c['root_and_roundoff_uncertainty'] + b['root_and_roundoff_uncertainty']) / 3
            difference = abs(fine - coarse)
            truncation = difference + coarse_bound + fine_bound
            uncertainty = truncation + fine_bound
            mismatch = abs(fine - expected)
            raw_coarse = abs(b['measured'] - a['measured'])
            raw_fine = abs(c['measured'] - b['measured'])
            contraction = raw_fine <= .6 * raw_coarse + 4 * (coarse_bound + fine_bound)
            stable_windows = (previous_window is not None and
                              difference <= .6 * previous_window['window_difference']
                              + 4 * (coarse_bound + fine_bound + previous_window['endpoint_uncertainty']))
            window = dict(h=[a['h'], b['h'], c['h']], measured=float(fine), coarse_extrapolation=float(coarse),
                          expected=expected, error=float(mismatch), absolute_uncertainty=float(uncertainty),
                          total_comparison_bound=float(mismatch + uncertainty),
                          nonlinear_uncertainty=float(truncation), endpoint_uncertainty=float(fine_bound),
                          window_difference=float(difference),
                          precision="mp%d" % mp_digits if precise else "float64",
                          contracting=bool(contraction), overlapping_windows_agree=bool(stable_windows),
                          same_branch=all(z['same_branch'] for z in (a, b, c)))
            windows.append(window)
            if best is None or window['total_comparison_bound'] < best['total_comparison_bound']:
                best = window
            if (window['same_branch'] and contraction and stable_windows
                    and window['total_comparison_bound'] <= target):
                chosen = [a, b, c]
                report = dict(status="PASS", certification=("ABSOLUTE_BOUNDED" if mismatch + uncertainty <= absolute_target else "RELATIVE_BOUNDED"),
                              expected=expected, measured=float(fine), error=float(mismatch),
                              absolute_uncertainty=float(uncertainty), relative_uncertainty=float(uncertainty / max(abs(expected), 1e-300)),
                              comparison_bound=float(mismatch + uncertainty), absolute_target=float(absolute_target),
                              relative_target=float(relative_target), precision_target=float(target),
                              h=float(c['h']), finest_field=float(c['h']), chosen_h=float(c['h']),
                              gradient_plus=float(c['plus']['gradient_max']), gradient_minus=float(c['minus']['gradient_max']),
                              maximum_endpoint_gradient=max(max(z['plus']['gradient_max'], z['minus']['gradient_max']) for z in chosen),
                              hessian_minimum=minimum, error_budget=dict(nonlinear=float(truncation), endpoint_and_roundoff=float(fine_bound)),
                              chosen_window=window, field_checks=attempts, extrapolation_windows=windows,
                              seconds=time.time() - start,
                              scope="Derivative in one positive-curvature branch; explicit residual and extrapolation error budget, not a global/domain-switching susceptibility")
                return report
            previous_window = window
            # Remaining smaller double fields only enlarge the rounding floor.
            if not precise and fine_bound > target * 2 and index >= 3:
                break
    report = dict(status="FAIL", reason="No converged same-branch extrapolation with a sufficient residual/roundoff budget",
                  expected=expected, precision_target=target, absolute_target=float(absolute_target),
                  relative_target=float(relative_target), best_window=best,
                  field_checks=attempts, extrapolation_windows=windows, seconds=time.time() - start)
    if raise_on_fail:
        raise ResponseCertificationError(report)
    return report


def _diagnostic_case(case):
    result_root = TASK.parents[3] / "limit/results/20261029" / TASK.name
    L, lam = case
    name = ("L%d-lambda%.2f" % (L, lam)).replace(".", "p")
    model = AccuratePhysicalStaticModel(L, lam)
    x = np.load(result_root / "cross2" / name / "equilibrium.npy")
    original = result_root / "original-response-certification/certified" / name / "validation-attempts.json"
    attempt_path = original if original.exists() else result_root / "certified" / name / "validation-attempts.json"
    attempts = json.loads(attempt_path.read_text())
    old = max(attempts[-1]['finite_field_checks'], key=lambda z: z['chi_hessian'])
    phase = 2 * np.pi * (model.coords @ np.array([old['n1'] / L, old['n2'] / L]))
    source = np.zeros(model.nb)
    source[2::3] = np.cos(phase) if old['quadrature'] == 'cos' else -np.sin(phase)
    result = certify_source(model, x, source, raise_on_fail=False)
    result.update(point=name, source_q=[old['n1'], old['n2']], quadrature=old['quadrature'])
    output = result_root / "accurate-response-diagnostics"
    output.mkdir(exist_ok=True)
    (output / (name + ".json")).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return {k: result.get(k) for k in ['point', 'status', 'certification', 'expected', 'measured', 'error', 'absolute_uncertainty', 'seconds']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--diagnose', action='store_true')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--point', nargs=2, metavar=('L', 'LAMBDA'))
    args = parser.parse_args()
    cases = [(int(args.point[0]), float(args.point[1]))] if args.point else [(8, z) for z in [.1, .2, .3, .4, .5, .6, .7]] + [(7, 1.8)]
    if args.diagnose:
        with ProcessPoolExecutor(max_workers=min(args.workers, 8), mp_context=multiprocessing.get_context('fork')) as pool:
            for future in as_completed([pool.submit(_diagnostic_case, case) for case in cases]):
                print(json.dumps(future.result(), allow_nan=False), flush=True)


def _cache_identity(model, x, source, expected, start_h_factor):
    """Exact request identity.  Alternate field grids have distinct identities."""
    import hashlib
    x = np.ascontiguousarray(x, dtype=np.float64)
    source = np.ascontiguousarray(source, dtype=np.float64)
    identity = dict(L=int(model.lx), lambda_value=float(model.lambda_value),
                    g=float(model.g), k=float(model.k), target_detT_sign=int(model.trial_sign),
                    background_sha256=hashlib.sha256(x.tobytes()).hexdigest(),
                    source_sha256=hashlib.sha256(source.tobytes()).hexdigest(),
                    expected_hex=float(expected).hex(), start_h_factor_hex=float(start_h_factor).hex(),
                    absolute_target=ABSOLUTE_TARGET, relative_target=SOFT_RELATIVE_TARGET)
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    return key, identity


def save_diagnostic_cache(model, x, source, report, cache_dir, *, start_h_factor=1.):
    """Package already calculated endpoint evidence with exact canonical inputs.

    This helper changes no numerical solver.  It rebases the comparison against
    the symmetrized dense inverse used by the integration wrapper, without
    changing a field measurement or its uncertainty.  The rebased comparison
    must independently satisfy the same declared precision target.
    """
    report = strengthen_residual_bounds(report, model.nc)
    audit = audit_certificate_windows(report, model.nc)
    if not audit['pass']:
        raise ValueError('Stored endpoint evidence fails the strengthened window/residual audit: ' + str(audit))
    x = newton_refine(model, x)
    h, _, _ = signed_hessian(model, x)
    inverse = np.linalg.inv(h)
    inverse = (inverse + inverse.T) / 2
    expected = float(source @ inverse @ source / model.nc)
    rebased = dict(report, strengthened_window_audit=audit)
    measured = float(report['measured'])
    mismatch = abs(measured - expected)
    budget = float(report['absolute_uncertainty'])
    target = max(ABSOLUTE_TARGET, SOFT_RELATIVE_TARGET * abs(expected))
    if report['status'] != 'PASS' or mismatch + budget > target:
        raise ValueError('A failed or insufficient certificate cannot enter the exact cache')
    rebased.update(expected=expected, error=mismatch, comparison_bound=mismatch + budget,
                   precision_target=target,
                   certification='ABSOLUTE_BOUNDED' if mismatch + budget <= ABSOLUTE_TARGET else 'RELATIVE_BOUNDED')
    key, identity = _cache_identity(model, x, source, expected, start_h_factor)
    directory = Path(cache_dir)
    directory.mkdir(parents=True, exist_ok=True)
    for filename, array in ((key + '-background.npy', x),
                            (key + '-source.npy', np.asarray(source, float))):
        target_path = directory / filename
        if not target_path.exists():
            temporary = target_path.with_name(target_path.name + '.tmp.%d' % os.getpid())
            with temporary.open('wb') as stream:
                np.save(stream, array)
            os.replace(temporary, target_path)
    record = dict(identity=identity, certificate=rebased,
                  background_file=key + '-background.npy', source_file=key + '-source.npy')
    record_path = directory / (key + '.json')
    temporary = record_path.with_name(record_path.name + '.tmp.%d' % os.getpid())
    temporary.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    os.replace(temporary, record_path)
    return dict(key=key, expected=expected, background=str(directory / record['background_file']),
                record=str(directory / (key + '.json')))


def load_diagnostic_cache(model, x, source, expected, cache_dir, *, start_h_factor=1.):
    """Read only exact matching requests; return None for any mismatch."""
    key, identity = _cache_identity(model, x, source, expected, start_h_factor)
    path = Path(cache_dir) / (key + '.json')
    if not path.exists():
        return None
    record = json.loads(path.read_text())
    if record['identity'] != identity or record['certificate']['status'] != 'PASS':
        return None
    result = dict(record['certificate'])
    result['exact_input_cache'] = dict(key=key, record=str(path))
    return result


def audit_certificate_windows(report, nc):
    """Strengthen earlier certificates using their unchanged endpoint evidence."""
    if report.get('status') != 'PASS':
        return dict(**{'pass': False}, reason='certificate_not_pass')
    chosen = report['chosen_window']
    previous = None
    contraction = False
    for window in report['extrapolation_windows']:
        difference = abs(window['measured'] - window['coarse_extrapolation'])
        fine_bound = window['endpoint_uncertainty']
        coarse_bound = max(0., window['nonlinear_uncertainty'] - difference - fine_bound)
        same_precision = previous is not None and previous['precision'] == window['precision']
        if window['h'] == chosen['h'] and window['precision'] == chosen['precision']:
            previous_difference = (abs(previous['measured'] - previous['coarse_extrapolation'])
                                   if same_precision else 0.)
            contraction = (same_precision and difference <= .6 * previous_difference
                           + 4 * (coarse_bound + fine_bound + previous['endpoint_uncertainty']))
            break
        previous = window
    aligned = True
    used = []
    for field in report['field_checks']:
        if field.get('h') not in chosen['h'] or 'plus' not in field:
            continue
        if field['plus']['precision'] != chosen['precision']:
            continue
        used.append(field)
        for endpoint in (field['plus'], field['minus']):
            # Initial version could take a final unevaluated MP correction
            # only when its third residual had not reached this exact bound.
            if endpoint.get('mp_gradient_evaluations', 0) == 3:
                allowance = report['precision_target'] * abs(endpoint['field']) * nc / 64
                if max(endpoint['projected_gradient'], endpoint['norm_gradient_projection_bound']) >= allowance:
                    aligned = False
    bounded = report['comparison_bound'] <= report['precision_target']
    passed = bool(contraction and aligned and bounded and len(used) == 3)
    return dict(**{'pass': passed}, successive_richardson_differences_contract=bool(contraction),
                endpoint_observables_match_evaluated_residuals=bool(aligned),
                comparison_within_declared_target=bool(bounded), chosen_fields=len(used))


def strengthen_residual_bounds(report, nc):
    """Recompute a certificate with Cauchy residual bounds, using saved data.

    No energies, gradients, fields or declared precision targets are changed.
    Every stored window is reconsidered; failure is explicit if none survives.
    """
    import copy
    revised = copy.deepcopy(report)
    revised['prior_projected_residual_summary'] = {
        k: report.get(k) for k in ('status', 'certification', 'absolute_uncertainty',
                                  'comparison_bound', 'chosen_h')}
    for field in revised.get('field_checks', []):
        if 'plus' not in field:
            continue
        for endpoint in (field['plus'], field['minus']):
            if endpoint.get('cauchy_norm_bound_used'):
                continue
            projected = endpoint['projected_gradient']
            old = endpoint['root_numerator_bound']
            roundoff = max(0., old / 2 - projected)
            norm_bound = max(projected, endpoint['norm_gradient_projection_bound'])
            endpoint['root_numerator_bound'] = 2 * (norm_bound + roundoff)
            endpoint['cauchy_norm_bound_used'] = True
        field['root_and_roundoff_uncertainty'] = (
            (field['plus']['root_numerator_bound'] + field['minus']['root_numerator_bound'])
            / (2 * field['h'] * nc) + field.get('observable_roundoff_uncertainty', 0.))
    lookup = {(f['h'], f['plus']['precision']): f
              for f in revised.get('field_checks', []) if 'plus' in f}
    eligible = []
    previous = None
    for window in revised.get('extrapolation_windows', []):
        records = [lookup[(z, window['precision'])] for z in window['h']]
        a, b, c = records
        difference = abs(window['measured'] - window['coarse_extrapolation'])
        coarse_bound = (4 * b['root_and_roundoff_uncertainty'] + a['root_and_roundoff_uncertainty']) / 3
        fine_bound = (4 * c['root_and_roundoff_uncertainty'] + b['root_and_roundoff_uncertainty']) / 3
        truncation = difference + coarse_bound + fine_bound
        uncertainty = truncation + fine_bound
        mismatch = abs(window['measured'] - revised['expected'])
        raw_contraction = (abs(c['measured'] - b['measured'])
                           <= .6 * abs(b['measured'] - a['measured'])
                           + 4 * (coarse_bound + fine_bound))
        same_precision = previous is not None and previous['precision'] == window['precision']
        consecutive_contraction = (same_precision and difference <= .6 * previous['window_difference']
                                   + 4 * (coarse_bound + fine_bound + previous['endpoint_uncertainty']))
        residual_alignment = True
        for field in records:
            for endpoint in (field['plus'], field['minus']):
                if endpoint.get('mp_gradient_evaluations', 0) == 3:
                    allowance = revised['precision_target'] * abs(endpoint['field']) * nc / 64
                    if max(endpoint['projected_gradient'], endpoint['norm_gradient_projection_bound']) >= allowance:
                        residual_alignment = False
        window.update(window_difference=float(difference), absolute_uncertainty=float(uncertainty),
                      total_comparison_bound=float(mismatch + uncertainty),
                      nonlinear_uncertainty=float(truncation), endpoint_uncertainty=float(fine_bound),
                      contracting=bool(raw_contraction), overlapping_windows_agree=bool(consecutive_contraction),
                      residuals_match_returned_values=bool(residual_alignment),
                      error=float(mismatch))
        if (raw_contraction and consecutive_contraction and residual_alignment and window['same_branch']
                and mismatch + uncertainty <= revised['precision_target']):
            eligible.append((window['total_comparison_bound'], window, records))
        previous = window
    revised['residual_bound_policy'] = '2*(||H_endpoint^-1 source||_2*||gradient_endpoint||_2 + SVD rounding floor), propagated through both field roots and Richardson coefficients'
    if not eligible:
        revised.update(status='FAIL', reason='No saved window satisfies the original declared target under the Cauchy residual bound')
        return revised
    _, window, records = min(eligible, key=lambda z: z[0])
    a, b, c = records
    uncertainty = window['absolute_uncertainty']
    revised.update(status='PASS', measured=window['measured'], error=window['error'],
                   absolute_uncertainty=uncertainty,
                   relative_uncertainty=uncertainty / max(abs(revised['expected']), 1e-300),
                   comparison_bound=window['total_comparison_bound'],
                   certification='ABSOLUTE_BOUNDED' if window['total_comparison_bound'] <= revised['absolute_target'] else 'RELATIVE_BOUNDED',
                   h=c['h'], finest_field=c['h'], chosen_h=c['h'],
                   gradient_plus=c['plus']['gradient_max'], gradient_minus=c['minus']['gradient_max'],
                   maximum_endpoint_gradient=max(max(z['plus']['gradient_max'], z['minus']['gradient_max']) for z in records),
                   error_budget=dict(nonlinear=window['nonlinear_uncertainty'], endpoint_and_roundoff=window['endpoint_uncertainty']),
                   chosen_window=window)
    return revised
