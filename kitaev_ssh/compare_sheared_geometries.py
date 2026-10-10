#!/usr/bin/env python3
"""Enumerate integer-sheared torus candidates and report finite-size geometry.

This tool reports candidates; it never selects a ``physically best`` torus and
never runs QMC.  Momentum distance, period shape, Kekule commensurability and
the uniform (++ Wilson) free-spectrum zero-mode check are separate diagnostics.
"""
import argparse
import csv
import json
import math
from pathlib import Path

from predict_vbs_momenta import (PRIMITIVE_A1, PRIMITIVE_A2, dirac_nodes,
                                 kekule_status, momentum_distance, nearest_grid,
                                 nearest_grid_set,
                                 q_predictions, uniform_background_solution,
                                 zero_mode_report)


def cell_cartesian(vector):
    x, y = vector
    return (x * PRIMITIVE_A1[0] + y * PRIMITIVE_A2[0],
            x * PRIMITIVE_A1[1] + y * PRIMITIVE_A2[1])


def dot(left, right):
    return left[0] * right[0] + left[1] * right[1]


def norm(vector):
    return math.hypot(*cell_cartesian(vector))


def subtract(left, right):
    return (left[0] - right[0], left[1] - right[1])


def scale(coefficient, vector):
    return (coefficient * vector[0], coefficient * vector[1])


def reduced_period_basis(Lx, Ly, shear):
    """Gauss-reduce the integer period basis using the physical metric."""
    a = (Lx, 0)
    b = (shear, Ly)
    for _ in range(100):
        if norm(b) < norm(a) - 1e-14:
            a, b = b, a
            continue
        coefficient = int(round(dot(cell_cartesian(a), cell_cartesian(b)) /
                                max(norm(a) ** 2, 1e-300)))
        candidate = subtract(b, scale(coefficient, a))
        if candidate == b:
            break
        b = candidate
    return a, b


def parse_targets(values):
    targets = []
    for index, value in enumerate(values):
        parts = value.split(",")
        if len(parts) != 2:
            raise ValueError("--target expects U,V")
        u, v = float(parts[0]), float(parts[1])
        if not math.isfinite(u) or not math.isfinite(v):
            raise ValueError("target fractions must be finite")
        targets.append({"label": f"explicit_{index}", "role": "explicit_target",
                        "fraction": (u % 1., v % 1.), "coupling_class": None})
    return targets


def model_targets(J, uniform_solution=None):
    result = []
    classes = [("bare", tuple(J), "bare_target")]
    if uniform_solution is not None:
        classes.append(("uniform_effective", tuple(uniform_solution["Jeff"]),
                        "uniform_effective_target"))
    for name, couplings, role in classes:
        info = dirac_nodes(*couplings)
        if info["status"] != "ok":
            raise ValueError(f"{name} Dirac predictor failed: {info['status']}")
        for index, point in enumerate(q_predictions(info["nodes"])):
            result.append({
                "label": f"{name}_Q{'p' if index == 0 else 'm'}",
                "role": role,
                "fraction": (float(point["u"]), float(point["v"])),
                "coupling_class": name,
            })
    return result


def explicit_geometries(values):
    result = []
    for value in values:
        parts = value.split(",")
        if len(parts) != 3:
            raise ValueError("--geometry expects Lx,Ly,shear")
        item = tuple(int(part) for part in parts)
        result.append(item)
    return result


def enumerate_geometries(cell_counts, explicit, lx_min, lx_max, ly_min, ly_max):
    if explicit:
        return explicit
    result = []
    for cells in cell_counts:
        for Lx in range(lx_min, lx_max + 1):
            if cells % Lx:
                continue
            Ly = cells // Lx
            if Ly < ly_min or Ly > ly_max or Lx < 2 or Ly < 2:
                continue
            for shear in range(Lx):
                result.append((Lx, Ly, shear))
    return result


def describe(Lx, Ly, shear, targets, coupling_sets):
    basis_a, basis_b = reduced_period_basis(Lx, Ly, shear)
    length_a, length_b = norm(basis_a), norm(basis_b)
    shortest, longest = sorted((length_a, length_b))
    period0 = cell_cartesian((Lx, 0))
    period1 = cell_cartesian((shear, Ly))
    result = {
        "Lx": Lx, "Ly": Ly, "shear": shear, "Nc": Lx * Ly,
        "sites": 2 * Lx * Ly,
        "geometry_id": f"honeycomb-pbc-Lx{Lx}-Ly{Ly}-shear{shear}",
        "pbc_geometry": "axis_aligned" if shear == 0 else "sheared",
        "torus_matrix": [[Lx, shear], [0, Ly]],
        "period_vectors_cartesian": [list(period0), list(period1)],
        "reduced_period_basis_cell": [list(basis_a), list(basis_b)],
        "reduced_period_basis_cartesian": [list(cell_cartesian(basis_a)),
                                            list(cell_cartesian(basis_b))],
        "shortest_nonzero_period_length": shortest,
        "period_shape_ratio": longest / max(shortest, 1e-300),
        "period_basis_angle_cosine": dot(cell_cartesian(basis_a), cell_cartesian(basis_b)) /
                                     max(length_a * length_b, 1e-300),
        "targets": [],
        "kekule": {},
        "free_fermion_diagnostics": {
            name: zero_mode_report(Lx, Ly, shear, couplings)
            for name, couplings in coupling_sets.items()
        },
    }
    for target in targets:
        fraction = target["fraction"]
        nearest_candidates = nearest_grid_set(fraction, Lx, Ly, shear)
        point = nearest_candidates[0]
        du, dv = point["delta"]
        distance = momentum_distance(fraction, point["fraction"])
        if abs(distance["cartesian_distance"] - point["cartesian_distance"]) > 2e-14:
            raise RuntimeError("nearest-grid Cartesian distance disagrees with reciprocal basis")
        result["targets"].append({
            "label": target["label"],
            "role": target["role"],
            "coupling_class": target["coupling_class"],
            "target_fraction": list(fraction),
            "nearest": point,
            "nearest_candidates": nearest_candidates,
            "cartesian_mismatch": distance["cartesian_distance"],
            "metric_distance": distance["metric_distance"],
            "metric_fractional_delta": distance["fractional_delta"],
            "delta1": 2. * math.pi * Lx * du,
            "delta2": 2. * math.pi * (shear * du + Ly * dv),
            "delta1_over_2pi": Lx * du,
            "delta2_over_2pi": shear * du + Ly * dv,
        })
    ks, on_grid = kekule_status(Lx, Ly, shear)
    result["kekule"] = {
        "on_grid": on_grid,
        "condition": "Lx % 3 == 0 and (shear-Ly) % 3 == 0",
        "points": ks,
    }
    return result


def flatten(report, target_index, target):
    item = report["targets"][target_index]
    nearest = item["nearest"]
    return {
        "geometry_id": report["geometry_id"], "Lx": report["Lx"], "Ly": report["Ly"],
        "shear": report["shear"], "Nc": report["Nc"], "sites": report["sites"],
        "target_index": target_index, "target_label": target["label"],
        "target_role": target["role"], "coupling_class": target["coupling_class"],
        "target_u": target["fraction"][0], "target_v": target["fraction"][1], "n1": nearest["n1"],
        "n2": nearest["n2"], "q_index": nearest["index"],
        "q_u": nearest["fraction"][0], "q_v": nearest["fraction"][1],
        "cartesian_mismatch": item["cartesian_mismatch"],
        "metric_distance": item["metric_distance"], "delta1": item["delta1"],
        "delta2": item["delta2"],
        "shortest_nonzero_period_length": report["shortest_nonzero_period_length"],
        "period_shape_ratio": report["period_shape_ratio"],
        "kekule_on_grid": report["kekule"]["on_grid"],
        "free_min_abs_f": (report["free_fermion_diagnostics"][target["coupling_class"]]["min_abs_f"]
                           if target["coupling_class"] is not None else ""),
        "zero_mode_risk": (report["free_fermion_diagnostics"][target["coupling_class"]]["zero_mode_risk"]
                           if target["coupling_class"] is not None else ""),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--target", action="append", default=[], metavar="U,V")
    parser.add_argument("--cells", default="36", help="comma-separated cell counts")
    parser.add_argument("--geometry", action="append", default=[], metavar="Lx,Ly,S",
                        help="explicit geometry; skips --cells enumeration")
    parser.add_argument("--Jx", type=float, default=1.)
    parser.add_argument("--Jy", type=float, default=1.)
    parser.add_argument("--Jz", type=float, default=.5)
    parser.add_argument("--uniform-effective", action="store_true",
                        help="add bare and self-consistent uniform-effective +/-Q targets")
    parser.add_argument("--g", type=float)
    parser.add_argument("--mass", type=float)
    parser.add_argument("--omega-ph", type=float, dest="omega_ph")
    parser.add_argument("--spring-k", type=float)
    parser.add_argument("--uniform-tolerance", type=float, default=1e-10)
    parser.add_argument("--uniform-max-iterations", type=int, default=200)
    parser.add_argument("--uniform-mixing", type=float, default=.5)
    parser.add_argument("--uniform-grid", type=int, default=400)
    parser.add_argument("--Lx-min", type=int, default=2)
    parser.add_argument("--Lx-max", type=int, default=10_000)
    parser.add_argument("--Ly-min", type=int, default=2)
    parser.add_argument("--Ly-max", type=int, default=10_000)
    parser.add_argument("--max-period-ratio", type=float, default=None,
                        help="mark shape_in_range using this upper bound; no rows are hidden")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--prefix", default="sheared_geometry")
    args = parser.parse_args()
    targets = parse_targets(args.target)
    couplings = (args.Jx, args.Jy, args.Jz)
    uniform_solution = None
    coupling_sets = {"bare": couplings}
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
        uniform_solution = uniform_background_solution(
            couplings, args.g, spring_k, args.uniform_tolerance,
            args.uniform_max_iterations, args.uniform_mixing, args.uniform_grid)
        uniform_solution.update({"J_bare": list(couplings), "g": args.g,
                                 "mass": args.mass, "omega_ph": args.omega_ph,
                                 "spring_source": spring_source})
        if not uniform_solution["converged"]:
            raise SystemExit(
                "uniform-effective predictor did not converge: "
                f"iterations={uniform_solution['iterations']} "
                f"residual={uniform_solution['final_residual_max_abs_X']:.17g}")
        coupling_sets["uniform_effective"] = tuple(uniform_solution["Jeff"])
        targets.extend(model_targets(couplings, uniform_solution))
    if not targets:
        raise SystemExit("provide at least one --target or use --uniform-effective")
    counts = [int(value) for value in args.cells.split(",") if value]
    geometries = enumerate_geometries(counts, explicit_geometries(args.geometry),
                                      args.Lx_min, args.Lx_max, args.Ly_min, args.Ly_max)
    if not geometries:
        raise SystemExit("no legal geometries")
    reports = []
    for Lx, Ly, shear in geometries:
        if Lx < 2 or Ly < 2 or shear < 0 or shear >= Lx:
            raise SystemExit(f"invalid geometry {(Lx, Ly, shear)}")
        report = describe(Lx, Ly, shear, targets, coupling_sets)
        report["shape_in_range"] = (args.max_period_ratio is None or
                                     report["period_shape_ratio"] <= args.max_period_ratio)
        reports.append(report)
    output = {"schema": "sheared-pbc-geometry-comparison-v2",
              "method": "enumerate legal integer quotient geometries; report all candidates",
              "targets": [{**target, "fraction": list(target["fraction"])} for target in targets],
              "couplings": list(couplings),
              "coupling_sets": {name: list(values) for name, values in coupling_sets.items()},
              "uniform_background": uniform_solution,
              "max_period_ratio": args.max_period_ratio,
              "reports": reports,
              "selection": None,
              "target_comparisons": [],
              "warning": "Raw metrics only: no best-geometry score is computed. Momentum distance is not a physical error estimate or a proof of finite-size convergence; no QMC was run."}
    bare = [target for target in targets if target["coupling_class"] == "bare"]
    for other_class in ("uniform_effective",):
        other = [target for target in targets if target["coupling_class"] == other_class]
        for left, right in zip(bare, other):
            distance = momentum_distance(right["fraction"], left["fraction"])
            output["target_comparisons"].append({
                "from": left["label"], "to": right["label"], **distance,
                "units": "HoneycombLattice reciprocal Cartesian units",
            })
    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / f"{args.prefix}_report.json"
    json_path.write_text(json.dumps(output, indent=2) + "\n")
    rows = [flatten(report, index, target) for report in reports for index, target in enumerate(targets)]
    csv_path = args.output_dir / f"{args.prefix}_table.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {json_path}")
    print(f"wrote {csv_path}")


if __name__ == "__main__":
    main()
