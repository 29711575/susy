#!/usr/bin/env python
"""Predict iVBS measurement momenta for the selected-mode VBS analysis.

Physics input (kitaev_ssh/CONVENTIONS.md section 6, verified against the
source): the itinerant Majorana spectrum is f(k)=Jz+Jx*exp(-i k1)+Jy*exp(-i k2)
with k_i = k.a_i in fractional cell coordinates. Isolated Dirac nodes K_tilde
exist only inside the triangle |Jx-Jy| <= Jz <= Jx+Jy. The incommensurate VBS
order parameter carries the intravalley scattering momentum Q = +/-2*K_tilde
(mod reciprocal lattice); a Kekule (period-3) VBS instead lives at
K+/- = (1/3,-1/3),(-1/3,1/3).

DISCLAIMER baked into every output: predictions from the BARE couplings are
search candidates for the selected-mode finite-window measurement. They are
NOT established peak positions of the interacting model; effective-coupling
predictions are a better but still approximate guide. Only the measured
Ssel/Dsel outputs can locate the actual maxima.

Nearest grid points use the true reciprocal-space metric
|du*b1+dv*b2|^2 proportional to du^2+dv^2-du*dv (never per-coordinate
rounding), and the mismatch report gives delta_i = (Q-q_L).(L_i a_i)
= 2*pi*L_i*du_i with the metric-minimal equivalent representative.

Example (L=4, Jx=Jy=1, Jz=0.5):
  /home/asus/.venv/bin/python predict_vbs_momenta.py --Lx 4 --Ly 4 \
      --Jx 1 --Jy 1 --Jz 0.5 --output-dir predicted_q_L4_Jz0.5
"""
import argparse
import csv
import json
import math
import os

import numpy as np

LABEL_CHARS = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.+-")
QFILE_DUPLICATE_TOLERANCE = 1e-9
GENERATOR_EQUIVALENCE_TOLERANCE = 1e-12
NEAREST_TIE_TOLERANCE = 1e-10


def frac(x):
    """Normalize a fractional coordinate into [0,1)."""
    return x - math.floor(x)


RECIPROCAL_B1 = (2. * math.pi / math.sqrt(3.), 2. * math.pi / 3.)
RECIPROCAL_B2 = (-2. * math.pi / math.sqrt(3.), 2. * math.pi / 3.)
PRIMITIVE_A1 = (math.sqrt(3.) / 2., 1.5)
PRIMITIVE_A2 = (-math.sqrt(3.) / 2., 1.5)


def metric_distance(du, dv):
    """Dimensionless reciprocal metric, without the common Cartesian factor."""
    return math.sqrt(max(du * du + dv * dv - du * dv, 0.0))


def reciprocal_cartesian(u, v):
    """Convert fractional reciprocal coordinates to project Cartesian units."""
    return (u * RECIPROCAL_B1[0] + v * RECIPROCAL_B2[0],
            u * RECIPROCAL_B1[1] + v * RECIPROCAL_B2[1])


def cartesian_distance(du, dv):
    """Actual Cartesian distance for a fractional reciprocal displacement."""
    qx, qy = reciprocal_cartesian(du, dv)
    return math.hypot(qx, qy)


def wrapped_delta(du, dv):
    """Metric-minimal representative of target-minus-grid modulo reciprocal vectors."""
    best = None
    for m in range(-2, 3):
        for n in range(-2, 3):
            cand = (du + m, dv + n)
            d = metric_distance(*cand)
            if best is None or d < best[0] - 1e-15:
                best = (d, cartesian_distance(*cand), cand)
    return best  # (dimensionless metric, Cartesian distance, representative)


def momentum_distance(left, right):
    """Distance from ``right`` to ``left`` modulo reciprocal lattice vectors.

    Both the dimensionless metric and Cartesian distance use the same
    metric-minimal fractional representative.  Cartesian units are exactly
    those of HoneycombLattice::reciprocal; the conversion factor is 4*pi/3,
    not 4*pi/sqrt(3).
    """
    metric, cartesian, delta = wrapped_delta(left[0] - right[0],
                                             left[1] - right[1])
    return {"fractional_delta": [float(delta[0]), float(delta[1])],
            "metric_distance": float(metric),
            "cartesian_distance": float(cartesian)}


def allowed_grid(Lx, Ly, shear):
    """Return the dual characters satisfying both integer torus constraints."""
    result = []
    for n1 in range(Lx):
        for n2 in range(Ly):
            result.append({
                "n1": n1,
                "n2": n2,
                "index": n1 * Ly + n2,
                "fraction": (frac(n1 / Lx), frac(n2 / Ly - shear * n1 / (Lx * Ly))),
            })
    return result


def grid_point_for(q, Lx, Ly, shear, tolerance=1e-10):
    for point in allowed_grid(Lx, Ly, shear):
        du = min(abs(q[0] - point["fraction"][0]),
                 1. - abs(q[0] - point["fraction"][0]))
        dv = min(abs(q[1] - point["fraction"][1]),
                 1. - abs(q[1] - point["fraction"][1]))
        if du <= tolerance and dv <= tolerance:
            return point
    return None


def momentum_separation(left, right):
    """Euclidean fractional-coordinate separation modulo the reciprocal lattice."""
    du = frac(left[0] - right[0])
    dv = frac(left[1] - right[1])
    if du > .5:
        du -= 1.
    if dv > .5:
        dv -= 1.
    return math.hypot(du, dv)


def dirac_nodes(Jx, Jy, Jz):
    """Fractional (k1/2pi, k2/2pi) of the Dirac nodes of f(k), or None.

    Closed form: with p=(Jy^2-Jx^2-Jz^2)/(2 Jz) and y=+/-sqrt(Jx^2-p^2),
    z1 = Jx e^{-i k1} = p + i y and z2 = Jy e^{-i k2} = -Jz - p - i y solve
    Jz + z1 + z2 = 0. Jz -> 0 is degenerate (nodal line for |Jx|=|Jy|) and is
    reported, never silently faked. Closed-form nodes are verified against
    |f(k)| with a Newton-on-(k1,k2) numerical fallback.
    """
    scale = max(abs(Jx), abs(Jy), abs(Jz))
    if scale <= 0:
        return {"status": "degenerate_zero_couplings", "nodes": []}
    if Jz == 0.:
        if abs(Jx) == abs(Jy):
            return {"status": "degenerate_nodal_line_no_isolated_dirac_points", "nodes": []}
        return {"status": "no_dirac_points_gapped", "nodes": []}
    p = (Jy * Jy - Jx * Jx - Jz * Jz) / (2. * Jz)
    disc = Jx * Jx - p * p
    if disc < -1e-12 * max(scale * scale, 1.):
        return {"status": "no_dirac_points_gapped", "nodes": []}
    disc = max(disc, 0.)
    nodes, seen = [], []
    for sign in (+1., -1.):
        y = sign * math.sqrt(disc)
        z1 = complex(p, y)
        z2 = complex(-Jz - p, -y)
        if abs(z1) <= 0 or abs(z2) <= 0:
            continue
        k1 = -math.atan2(z1.imag, z1.real)
        k2 = -math.atan2(z2.imag, z2.real)
        node = (frac(k1 / (2. * math.pi)), frac(k2 / (2. * math.pi)))
        if any(abs(node[0] - a) + abs(node[1] - b) < 1e-12 for a, b in seen):
            continue
        seen.append(node)
        nodes.append({"k_fraction": node, "closed_form": True})
    # Boundary of the triangle: nodes merge and are not isolated.
    triangle = abs(Jx - Jy) - 1e-12 * scale <= abs(Jz) <= Jx + Jy + 1e-12 * scale
    merged = len(nodes) < 2
    if not triangle:
        return {"status": "no_dirac_points_gapped", "nodes": []}
    if merged:
        return {"status": "dirac_points_merged_not_isolated", "nodes": []}
    # Numerical verification of f(k)=0 with Newton fallback.
    verified = []
    for node in nodes:
        k = 2. * math.pi * np.array(node["k_fraction"])
        ok = abs(f_value(k, Jx, Jy, Jz)) <= 1e-9 * scale
        if not ok:
            k, ok_newton = newton_refine(k, Jx, Jy, Jz)
            ok = ok_newton and abs(f_value(k, Jx, Jy, Jz)) <= 1e-9 * scale
            node["closed_form"] = False
            node["newton_refined"] = True
        node["residual_abs_f"] = float(abs(f_value(k, Jx, Jy, Jz)))
        node["k_fraction"] = (frac(k[0] / (2. * math.pi)), frac(k[1] / (2. * math.pi)))
        if ok:
            verified.append(node)
    if len(verified) < 2:
        return {"status": "numerical_verification_failed_no_isolated_dirac_points", "nodes": []}
    # Deterministic order: the node with k1 fraction in (0,0.5] first, so that
    # "+2 K_tilde" of the first node matches the Qp used in q_points_example.txt.
    verified.sort(key=lambda n: (0 if n["k_fraction"][0] <= 0.5 else 1, n["k_fraction"]))
    return {"status": "ok", "nodes": verified}


def f_value(k, Jx, Jy, Jz):
    return Jz + Jx * np.exp(-1j * k[0]) + Jy * np.exp(-1j * k[1])


def newton_refine(k, Jx, Jy, Jz, iterations=60):
    """Real 2x2 Newton on Re f = 0, Im f = 0 from a coarse-grid seed."""
    best_k, best = np.asarray(k, dtype=float), abs(f_value(k, Jx, Jy, Jz))
    grid = np.linspace(0., 2. * math.pi, 61)
    for a in grid:
        for b in grid:
            r = abs(f_value((a, b), Jx, Jy, Jz))
            if r < best:
                best, best_k = r, np.array([a, b])
    k = best_k
    for _ in range(iterations):
        f = f_value(k, Jx, Jy, Jz)
        if abs(f) < 1e-14:
            break
        j1 = -1j * Jx * np.exp(-1j * k[0])
        j2 = -1j * Jy * np.exp(-1j * k[1])
        jac = np.array([[j1.real, j2.real], [j1.imag, j2.imag]])
        try:
            step = np.linalg.solve(jac, np.array([-f.real, -f.imag]))
        except np.linalg.LinAlgError:
            break
        k = k + step
    return k, abs(f_value(k, Jx, Jy, Jz)) < 1e-9 * max(abs(Jx), abs(Jy), abs(Jz), 1.)


def q_predictions(nodes):
    """Q = +/-2*K_tilde (mod 1), deduplicated across the node pair."""
    out, seen = [], []
    for node in nodes:
        kf = node["k_fraction"]
        for sign in (+1., -1.):
            q = (frac(sign * 2. * kf[0]), frac(sign * 2. * kf[1]))
            if any(abs(q[0] - a) + abs(q[1] - b) < 1e-12 for a, b in seen):
                continue
            seen.append(q)
            out.append({"u": q[0], "v": q[1], "source_node": kf})
    return out


def nearest_grid(q, Lx, Ly, shear, count):
    """The `count` nearest allowed torus characters under both metrics."""
    rows = []
    for point in allowed_grid(Lx, Ly, shear):
        du, dv = q[0] - point["fraction"][0], q[1] - point["fraction"][1]
        metric, cartesian, (du, dv) = wrapped_delta(du, dv)
        rows.append({**point, "metric_distance": metric,
                     "cartesian_distance": cartesian, "delta": (du, dv)})
    rows.sort(key=lambda r: (r["cartesian_distance"], r["index"]))
    return rows[:count]


def nearest_grid_set(q, Lx, Ly, shear, tolerance=NEAREST_TIE_TOLERANCE):
    """All grid momenta tied for nearest under the Cartesian metric.

    The returned list is deterministic (distance, then grid index).  Tied
    points remain distinct physical momenta; this helper only groups them as
    nearest candidates and never feeds q-file coordinate deduplication.
    """
    rows = nearest_grid(q, Lx, Ly, shear, Lx * Ly)
    best = rows[0]["cartesian_distance"]
    tied = [r for r in rows if r["cartesian_distance"] <= best + tolerance]
    return sorted(tied, key=lambda r: r["index"])


def kekule_status(Lx, Ly, shear):
    # K+/- = (1/3,-1/3),(-1/3,1/3) in fractional reciprocal coordinates,
    # normalized into [0,1): (1/3,2/3) and (2/3,1/3).
    ks = [{"label": "kek_Kp", "u": frac(1. / 3.), "v": frac(-1. / 3.)},
          {"label": "kek_Km", "u": frac(-1. / 3.), "v": frac(1. / 3.)}]
    on_grid = (Lx % 3 == 0) and ((shear - Ly) % 3 == 0)
    for k in ks:
        k["nearest_set"] = nearest_grid_set((k["u"], k["v"]), Lx, Ly, shear)
        exact = grid_point_for((k["u"], k["v"]), Lx, Ly, shear)
        k["on_grid"] = on_grid and exact is not None
        k["grid_point"] = exact
    return ks, on_grid


def zero_mode_report(Lx, Ly, shear, couplings, tolerance=1e-10):
    """Check the uniform (++ Wilson) matter grid without altering parameters."""
    Jx, Jy, Jz = couplings
    values = []
    for point in allowed_grid(Lx, Ly, shear):
        u, v = point["fraction"]
        f = Jz + Jx * np.exp(-2j * math.pi * u) + Jy * np.exp(-2j * math.pi * v)
        values.append((float(abs(f)), point))
    minimum = min(values, key=lambda item: (item[0], item[1]["index"]))
    return {
        "evaluated": True,
        "wilson_sector": [1, 1],
        "min_abs_f": minimum[0],
        "zero_mode_tolerance": tolerance,
        "zero_mode_risk": minimum[0] <= tolerance,
        "minimum_points": [item[1] for item in values if item[0] <= minimum[0] + 1e-12],
        "method": "uniform-gauge f(k)=Jz+Jx exp(-2pi i u)+Jy exp(-2pi i v) on the sheared dual grid",
        "action_on_risk": "freeGroundState must reject; this predictor does not change couplings, sector or occupancy",
    }


def uniform_bond_means(couplings, integration_grid=400):
    """Thermodynamic uniform-gauge <u i c_A c_B> for x/y/z bonds.

    The project convention has E0/Nc=-<|f|> and
    H=sum_b u_b i (J_mu+g X_b)c_A c_B.  Hellmann--Feynman therefore gives
    B_mu=<u i c_A c_B>=d(E0/Nc)/d Jeff_mu.  Midpoint quadrature avoids
    evaluating an exactly gapless Dirac point while retaining a deterministic
    uniform Brillouin-zone rule.
    """
    if integration_grid < 16:
        raise ValueError("uniform integration grid must be >= 16")
    values = np.asarray(couplings, dtype=float)
    if values.shape != (3,) or not np.all(np.isfinite(values)):
        raise ValueError("uniform effective couplings must be three finite values")
    points = (np.arange(integration_grid, dtype=float) + .5) / integration_grid
    u, v = np.meshgrid(points, points, indexing="ij")
    phase_x = np.exp(-2j * math.pi * u)
    phase_y = np.exp(-2j * math.pi * v)
    f = values[2] + values[0] * phase_x + values[1] * phase_y
    magnitude = np.abs(f)
    if not np.all(np.isfinite(magnitude)) or np.any(magnitude <= 1e-14 * max(1., np.max(abs(values)))):
        raise RuntimeError("uniform bond quadrature encountered a zero/invalid |f(k)|")
    conjugate_direction = np.conj(f) / magnitude
    means = -np.array([
        np.mean(np.real(phase_x * conjugate_direction)),
        np.mean(np.real(phase_y * conjugate_direction)),
        np.mean(np.real(conjugate_direction)),
    ])
    if not np.all(np.isfinite(means)):
        raise RuntimeError("uniform bond expectation is non-finite")
    return means


def uniform_background_solution(J, g, spring_k, tolerance=1e-10,
                                max_iterations=200, mixing=.5,
                                integration_grid=400, initial_xbar=None):
    """Solve k Xbar_mu + g <u i c_A c_B>_mu = 0 by damped iteration."""
    J = np.asarray(J, dtype=float)
    if J.shape != (3,) or not np.all(np.isfinite(J)):
        raise ValueError("Jx/Jy/Jz must be finite")
    if not math.isfinite(g):
        raise ValueError("g must be finite")
    if not math.isfinite(spring_k) or spring_k <= 0:
        raise ValueError("spring k must be finite and positive")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("uniform tolerance must be finite and positive")
    if max_iterations < 1:
        raise ValueError("uniform max iterations must be positive")
    if not math.isfinite(mixing) or not (0. < mixing <= 1.):
        raise ValueError("uniform mixing must satisfy 0 < mixing <= 1")

    xbar = (np.zeros(3) if initial_xbar is None
            else np.asarray(initial_xbar, dtype=float).copy())
    if xbar.shape != (3,) or not np.all(np.isfinite(xbar)):
        raise ValueError("uniform initial Xbar must be three finite values")
    initial = xbar.copy()
    residual = math.inf
    bond_mean = np.full(3, np.nan)
    converged = False
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        jeff = J + g * xbar
        bond_mean = uniform_bond_means(jeff, integration_grid)
        target = -(g / spring_k) * bond_mean
        residual = float(np.max(np.abs(target - xbar)))
        if residual <= tolerance:
            converged = True
            break
        xbar = xbar + mixing * (target - xbar)
    jeff = J + g * xbar
    # Report the residual at the returned point, including the last damped
    # update if max_iterations was exhausted.
    bond_mean = uniform_bond_means(jeff, integration_grid)
    target = -(g / spring_k) * bond_mean
    residual = float(np.max(np.abs(target - xbar)))
    converged = residual <= tolerance
    return {
        "status": "converged" if converged else "not_converged",
        "converged": converged,
        "iterations": iterations,
        "max_iterations": max_iterations,
        "tolerance": tolerance,
        "final_residual_max_abs_X": residual,
        "mixing": mixing,
        "initial_Xbar": [float(value) for value in initial],
        "integration_grid": integration_grid,
        "integration_rule": "uniform midpoint quadrature on fractional Brillouin torus",
        "Xbar": [float(value) for value in xbar],
        "Jeff": [float(value) for value in jeff],
        "bond_mean_u_i_cA_cB": [float(value) for value in bond_mean],
        "spring_k": spring_k,
        "equation": "Xbar_mu=-(g/k)<u i c_A c_B>_mu; Jeff_mu=J_mu+g Xbar_mu",
        "hamiltonian_convention": "H=sum_<A,B> u i (J_mu+g X_mu)c_A c_B",
        "scope": "uniform translation-invariant free-Majorana self-consistency; not an iVBS saddle or QMC peak",
    }


def measurement_points(rows, Lx, Ly, shear):
    """Collapse only identical actual q-file coordinates, preserving aliases.

    Distinct off-grid targets are never collapsed merely because their nearest
    allowed grid label agrees.  A sub-1e-9 but non-identical pair is rejected:
    the strict C++ parser would reject it too, and silently merging it would
    change the requested measurement set.
    """
    points = []
    for row in rows:
        coordinate = (frac(float(row["u"])), frac(float(row["v"])))
        match = None
        for point in points:
            separation = momentum_separation(coordinate, point["reciprocal_fraction"])
            if separation <= GENERATOR_EQUIVALENCE_TOLERANCE:
                match = point
                break
            if separation < QFILE_DUPLICATE_TOLERANCE:
                raise ValueError(
                    "distinct generated q points are closer than the strict C++ "
                    f"duplicate tolerance: {row['label']} and {point['label']}")
        source = {"label": row["label"], "class": row["class"], "note": row["note"]}
        if match is not None:
            match["aliases"].append(row["label"])
            match["sources"].append(source)
            continue
        grid = grid_point_for(coordinate, Lx, Ly, shear)
        points.append({
            "index": len(points),
            "label": row["label"],
            "reciprocal_fraction": [coordinate[0], coordinate[1]],
            "phase_fraction": list(grid["fraction"] if grid is not None else coordinate),
            "on_grid": grid is not None,
            "n": [grid["n1"], grid["n2"]] if grid is not None else None,
            "grid_index": grid["index"] if grid is not None else None,
            "aliases": [row["label"]],
            "sources": [source],
            "role": "q_file_measurement_momentum",
        })
    return points


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 epilog="Bare-coupling predictions are search candidates, "
                                        "never established interacting peak positions.")
    ap.add_argument("--Lx", type=int, required=True)
    ap.add_argument("--Ly", type=int, required=True)
    ap.add_argument("--shear", type=int, default=0,
                    help="integer torus shear, 0 <= shear < Lx (default 0)")
    ap.add_argument("--Jx", type=float, default=1.)
    ap.add_argument("--Jy", type=float, default=1.)
    ap.add_argument("--Jz", type=float, default=0.5)
    ap.add_argument("--Jeff-x", type=float, dest="Jeff_x")
    ap.add_argument("--Jeff-y", type=float, dest="Jeff_y")
    ap.add_argument("--Jeff-z", type=float, dest="Jeff_z")
    ap.add_argument("--uniform-effective", action="store_true",
                    help="solve the translation-invariant phonon background and add its effective target")
    ap.add_argument("--g", type=float,
                    help="SSH coupling for --uniform-effective")
    ap.add_argument("--mass", type=float,
                    help="phonon mass; use with --omega-ph unless --spring-k is supplied")
    ap.add_argument("--omega-ph", type=float, dest="omega_ph",
                    help="phonon frequency; use with --mass unless --spring-k is supplied")
    ap.add_argument("--spring-k", type=float,
                    help="positive spring constant k; alternative to mass*omega_ph^2")
    ap.add_argument("--uniform-tolerance", type=float, default=1e-10)
    ap.add_argument("--uniform-max-iterations", type=int, default=200)
    ap.add_argument("--uniform-mixing", type=float, default=.5)
    ap.add_argument("--uniform-grid", type=int, default=400,
                    help="midpoint quadrature points per reciprocal direction (default 400)")
    ap.add_argument("--extra-q", action="append", default=[], metavar="LABEL,U,V",
                    help="user-supplied candidate momenta, reported with class 'user'")
    ap.add_argument("--neighbors", type=int, default=4,
                    help="nearest PBC grid points kept per prediction (default 4)")
    ap.add_argument("--output-dir", default=".")
    ap.add_argument("--prefix", default="ivbs")
    args = ap.parse_args()

    if args.Lx < 2 or args.Ly < 2:
        raise SystemExit("--Lx/--Ly must be >= 2 (honeycomb_lattice requirement)")
    if args.shear < 0 or args.shear >= args.Lx:
        raise SystemExit("--shear must satisfy 0 <= shear < Lx")
    classes = []  # (class_name, couplings)
    classes.append(("bare", (args.Jx, args.Jy, args.Jz)))
    if (args.Jeff_x is not None or args.Jeff_y is not None or args.Jeff_z is not None):
        if None in (args.Jeff_x, args.Jeff_y, args.Jeff_z):
            raise SystemExit("--Jeff-x/y/z must be given together or not at all")
        classes.append(("effective", (args.Jeff_x, args.Jeff_y, args.Jeff_z)))
    uniform_solution = None
    if args.uniform_effective:
        if args.g is None:
            raise SystemExit("--uniform-effective requires --g")
        if args.spring_k is None:
            if args.mass is None or args.omega_ph is None:
                raise SystemExit("--uniform-effective requires --spring-k or both --mass and --omega-ph")
            spring_k = args.mass * args.omega_ph * args.omega_ph
            spring_source = "mass*omega_ph^2"
        else:
            spring_k = args.spring_k
            spring_source = "explicit --spring-k"
        try:
            uniform_solution = uniform_background_solution(
                (args.Jx, args.Jy, args.Jz), args.g, spring_k,
                args.uniform_tolerance, args.uniform_max_iterations,
                args.uniform_mixing, args.uniform_grid)
        except (ValueError, RuntimeError) as error:
            raise SystemExit("uniform-effective predictor failed: " + str(error))
        uniform_solution.update({
            "J_bare": [args.Jx, args.Jy, args.Jz],
            "g": args.g,
            "mass": args.mass,
            "omega_ph": args.omega_ph,
            "spring_source": spring_source,
        })
        if not uniform_solution["converged"]:
            raise SystemExit(
                "uniform-effective predictor did not converge: "
                f"iterations={uniform_solution['iterations']} "
                f"residual={uniform_solution['final_residual_max_abs_X']:.17g} "
                f"tolerance={uniform_solution['tolerance']:.17g}")
        classes.append(("uniform_effective", tuple(uniform_solution["Jeff"])))

    disclaimer = ("Q_pred from bare couplings is a search candidate only; it is not the "
                  "established peak position of the interacting model. Effective-coupling "
                  "predictions are labeled separately and are still approximate. Only the "
                  "measured Ssel/Dsel observables locate the actual maxima.")
    period0=(args.Lx*PRIMITIVE_A1[0],args.Lx*PRIMITIVE_A1[1])
    period1=(args.shear*PRIMITIVE_A1[0]+args.Ly*PRIMITIVE_A2[0],
             args.shear*PRIMITIVE_A1[1]+args.Ly*PRIMITIVE_A2[1])
    report = {"Lx": args.Lx, "Ly": args.Ly, "shear": args.shear,
              "Nc": args.Lx * args.Ly, "sites": 2 * args.Lx * args.Ly,
              "pbc_geometry": "axis_aligned" if args.shear == 0 else "sheared",
              "geometry_id": f"honeycomb-pbc-Lx{args.Lx}-Ly{args.Ly}-shear{args.shear}",
              "torus_matrix": [[args.Lx, args.shear], [0, args.Ly]],
              "period_vectors_cartesian": [list(period0), list(period1)],
              "convention": "f(k)=Jz+Jx e^{-i k1}+Jy e^{-i k2}, k_i=k.a_i fractional; "
                            "Q_pred=+/-2*K_tilde mod reciprocal lattice",
              "momentum_grid_convention": "u=n1/Lx, v=n2/Ly-shear*n1/(Lx*Ly), index=n1*Ly+n2",
              "reciprocal_basis_cartesian": [list(RECIPROCAL_B1), list(RECIPROCAL_B2)],
              "metric": "metric_distance=sqrt(du^2+dv^2-du*dv); cartesian_distance=|du*b1+dv*b2|=(4*pi/3)*metric_distance in HoneycombLattice Cartesian units",
              "mismatch_convention": "DeltaQ=Q_target-q_L; delta1=DeltaQ.T1, delta2=DeltaQ.T2, with the metric-minimal reciprocal representative",
              "disclaimer": disclaimer,
              "momentum_roles": {
                  "bare_target": "continuous +/-2 K_tilde from bare J; geometry guide only",
                  "uniform_effective_target": "continuous +/-2 K_tilde from self-consistent uniform Jeff; improved geometry guide, not a final physical momentum",
                  "nearest_allowed_finite_size": "metric-nearest character of the selected finite torus",
                  "q_file_measurement_momentum": "actual coordinate requested from selected/both QMC by the generated q-file",
                  "qmc_structure_factor_peak": "not predicted here; must be inferred later from measured nonzero Fourier components",
              },
              "qmc_structure_factor_peak": {"status": "not_measured", "reciprocal_fraction": None},
              "uniform_background": uniform_solution,
              "coupling_classes": {}, "kekule": {}, "candidates": [], "warnings": []}

    rows = []  # CSV rows
    label_taken = set()

    def fresh(label):
        base, name = label, label
        while name in label_taken:
            name = base + "_" + str(len(label_taken))
        label_taken.add(name)
        return name

    def add_row(cls, label, u, v, note, grid=None):
        u, v = frac(u), frac(v)
        actual_grid = grid_point_for((u, v), args.Lx, args.Ly, args.shear)
        row = {"class": cls, "label": label, "u": u, "v": v, "note": note,
               "on_grid": actual_grid is not None,
               "n1": actual_grid["n1"] if actual_grid is not None else "",
               "n2": actual_grid["n2"] if actual_grid is not None else "",
               "grid_index": actual_grid["index"] if actual_grid is not None else "",
               "grid_fraction": "", "metric_distance": "", "cartesian_distance": "",
               "delta1_rad": "", "delta1_over_2pi": "",
               "delta2_rad": "", "delta2_over_2pi": ""}
        if actual_grid is not None:
            row["grid_fraction"] = actual_grid["fraction"]
        if grid is not None:
            du, dv = grid["delta"]
            row.update(metric_distance=grid["metric_distance"],
                       cartesian_distance=grid["cartesian_distance"],
                       delta1_rad=2. * math.pi * args.Lx * du,
                       delta1_over_2pi=args.Lx * du,
                       delta2_rad=2. * math.pi * (args.shear * du + args.Ly * dv),
                       delta2_over_2pi=args.shear * du + args.Ly * dv)
        report["candidates"].append(row)
        rows.append(row)

    for cls, (Jx, Jy, Jz) in classes:
        info = dirac_nodes(Jx, Jy, Jz)
        entry = {"J": [Jx, Jy, Jz], "status": info["status"],
                 "triangle_condition": "no isolated Dirac point unless |Jx-Jy| <= Jz <= Jx+Jy"}
        preds = []
        if info["status"] == "ok":
            preds = q_predictions(info["nodes"])
            for i, q in enumerate(preds):
                label = fresh(f"{cls}_Q{'p' if i == 0 else 'm'}")
                near = nearest_grid((q["u"], q["v"]), args.Lx, args.Ly, args.shear, args.neighbors)
                exact = near[0] if near and near[0]["metric_distance"] <= 1e-10 else None
                add_row(cls, label, q["u"], q["v"],
                        "continuous +/-2 K_tilde prediction; may be off the sheared PBC grid",
                        grid=exact)
                for g in near:
                    glabel = fresh(f"{label}_g{g['index']}")
                    add_row(cls, glabel, g["fraction"][0], g["fraction"][1],
                            "nearest sheared PBC grid point of " + label, grid=g)
                q["nearest_grid"] = near
        elif info["status"].startswith(("degenerate", "no_dirac", "dirac_points_merged",
                                        "numerical")):
            report["warnings"].append(
                f"{cls}: {info['status']}; no Q_pred is emitted rather than a fake "
                "prediction (Jz->0 with |Jx|=|Jy| is a degenerate nodal line).")
        entry["Q_pred"] = preds
        entry["nodes"] = info["nodes"]
        entry["target_role"] = ("bare_target" if cls == "bare" else
                                "uniform_effective_target" if cls == "uniform_effective" else
                                "externally_supplied_effective_target")
        report["coupling_classes"][cls] = entry

    report["target_comparisons"] = []
    bare_predictions = report["coupling_classes"].get("bare", {}).get("Q_pred", [])
    for other_class in ("effective", "uniform_effective"):
        other_predictions = report["coupling_classes"].get(other_class, {}).get("Q_pred", [])
        for index, (bare, other) in enumerate(zip(bare_predictions, other_predictions)):
            report["target_comparisons"].append({
                "from": f"bare_Q{'p' if index == 0 else 'm'}",
                "to": f"{other_class}_Q{'p' if index == 0 else 'm'}",
                **momentum_distance((other["u"], other["v"]),
                                    (bare["u"], bare["v"])),
                "units": "HoneycombLattice reciprocal Cartesian units",
            })

    ks, kek_on_grid = kekule_status(args.Lx, args.Ly, args.shear)
    report["kekule"] = {"period": 3, "on_grid": kek_on_grid,
                        "on_grid_requires": "Lx divisible by 3 and (shear-Ly) divisible by 3",
                        "points": []}
    for k in ks:
        label = fresh(k["label"])
        exact_rows=nearest_grid((k["u"],k["v"]),args.Lx,args.Ly,args.shear,args.Lx*args.Ly)
        exact=next((p for p in exact_rows if p["index"]==k["grid_point"]["index"]),None) if k["grid_point"] else None
        add_row("kekule", label, k["u"], k["v"],
                "Kekule (period-3) VBS momentum; "
                + ("on the sheared PBC grid" if k["on_grid"] else "off the sheared PBC grid for this L"),
                grid=exact)
        report["kekule"]["points"].append(
            {"label": label, "fraction": [k["u"], k["v"]],
             "nearest_grid_set": k["nearest_set"], "on_grid": k["on_grid"]})

    for spec in args.extra_q:
        parts = spec.split(",")
        if len(parts) != 3:
            raise SystemExit("--extra-q expects LABEL,U,V: " + spec)
        label, us, vs = parts[0].strip(), parts[1].strip(), parts[2].strip()
        if not label or set(label) - LABEL_CHARS:
            raise SystemExit("--extra-q label must use [A-Za-z0-9_.+-]: " + label)
        u, v = float(us), float(vs)
        add_row("user", fresh("user_" + label), frac(u), frac(v),
                f"user-supplied candidate (input {u},{v})")

    # Resolution warning: an iVBS prediction's nearest grid point (ties kept)
    # coincides with a Kekule nearest grid point -> this lattice cannot
    # separate incommensurate from Kekule ordering at that momentum.
    kek_grid = {p["index"] for k in ks for p in k["nearest_set"]}
    ivbs_nearest = {}
    for cls, _ in classes:
        entry = report["coupling_classes"].get(cls, {})
        for q in entry.get("Q_pred", []):
            tied = nearest_grid_set((q["u"], q["v"]), args.Lx, args.Ly, args.shear)
            q["nearest_grid_set"] = tied
            for p in tied:
                ivbs_nearest.setdefault(p["index"], set()).add(
                    f"{cls}:{q['u']:.6f},{q['v']:.6f}")
    clash = sorted(kek_grid & set(ivbs_nearest))
    if clash:
        report["warnings"].append(
            f"grid indices {clash} are simultaneously a nearest grid point of a Kekule "
            "K point and of an iVBS prediction: this finite lattice cannot resolve "
            "incommensurate from Kekule ordering; treat Ssel there as ambiguous.")
    report["kekule"]["coincident_grid_indices_with_ivbs"] = clash
    report["zero_mode_risk"] = zero_mode_report(
        args.Lx, args.Ly, args.shear, (args.Jx, args.Jy, args.Jz))

    try:
        qpoints = measurement_points(rows, args.Lx, args.Ly, args.shear)
    except ValueError as error:
        raise SystemExit("cannot generate a parser-clean q-file: " + str(error))
    report["q_file_measurement_points"] = qpoints
    report["q_file_deduplication"] = {
        "input_candidate_rows": len(rows),
        "written_unique_points": len(qpoints),
        "merged_alias_count": len(rows) - len(qpoints),
        "rule": "merge actual reciprocal coordinates equal modulo the reciprocal lattice; never merge distinct off-grid targets because nearest-grid labels agree",
        "generator_equivalence_tolerance": GENERATOR_EQUIVALENCE_TOLERANCE,
        "cpp_parser_duplicate_tolerance": QFILE_DUPLICATE_TOLERANCE,
        "order": "stable first occurrence in bare, external effective, uniform-effective, Kekule, then user-extra source order; nearest neighbors use grid index tie-breaks",
    }

    os.makedirs(args.output_dir, exist_ok=True)
    qfile = os.path.join(args.output_dir, args.prefix + "_qfile.txt")
    with open(qfile, "w") as f:
        f.write("# Generated by predict_vbs_momenta.py; format 'label u v' with "
                "Q=u*b1+v*b2 (fractional reciprocal coordinates).\n")
        f.write(f"# Geometry: Lx={args.Lx} Ly={args.Ly} shear={args.shear}; "
                "grid labels are n1,n2 with v=n2/Ly-shear*n1/(Lx*Ly).\n")
        f.write("# " + disclaimer + "\n")
        for point in qpoints:
            u, v = point["reciprocal_fraction"]
            aliases = ",".join(point["aliases"])
            f.write(f"{point['label']} {u:.17g} {v:.17g} # aliases={aliases}\n")
    csv_path = os.path.join(args.output_dir, args.prefix + "_candidates.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    json_path = os.path.join(args.output_dir, args.prefix + "_report.json")
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"wrote {qfile}\nwrote {csv_path}\nwrote {json_path}")
    for w in report["warnings"]:
        print("WARNING:", w)


if __name__ == "__main__":
    main()
