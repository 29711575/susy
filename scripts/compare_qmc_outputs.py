#!/usr/bin/env python3
"""Compare two completed coupled-MC outputs without discarding complex parts."""

import argparse
import csv
import hashlib
import json
import math
import pathlib
import sys


class Comparison:
    def __init__(self, tolerance, ignored_suffixes=()):
        self.tolerance = tolerance
        self.ignored_suffixes = tuple(ignored_suffixes)
        self.numeric_difference_count = 0
        self.exact_difference_count = 0
        self.max_abs = 0.0
        self.max_normalized = 0.0
        self.max_relative = 0.0
        self.first_difference = None
        self.differences = []

    def note_numeric(self, path, left, right):
        left = float(left)
        right = float(right)
        absolute = abs(left - right)
        scale = max(1.0, abs(left), abs(right))
        relative_scale = max(abs(left), abs(right), 1.0e-300)
        normalized = absolute / scale
        relative = absolute / relative_scale
        self.numeric_difference_count += absolute != 0.0
        self.max_abs = max(self.max_abs, absolute)
        self.max_normalized = max(self.max_normalized, normalized)
        self.max_relative = max(self.max_relative, relative)
        if absolute != 0.0 and self.first_difference is None:
            self.first_difference = {
                "path": path,
                "left": left,
                "right": right,
                "absolute_error": absolute,
                "normalized_error": normalized,
                "relative_error": relative,
            }
        if absolute > self.tolerance and len(self.differences) < 20:
            self.differences.append({
                "path": path,
                "left": left,
                "right": right,
                "absolute_error": absolute,
                "normalized_error": normalized,
                "relative_error": relative,
            })

    def note_exact(self, path, left, right):
        self.exact_difference_count += 1
        if len(self.differences) < 20:
            self.differences.append({"path": path, "left": left, "right": right})

    def compare(self, left, right, path):
        if any(path == suffix or path.endswith("." + suffix)
               for suffix in self.ignored_suffixes):
            return
        if isinstance(left, bool) or isinstance(right, bool):
            if left != right:
                self.note_exact(path, left, right)
            return
        if isinstance(left, (int, float)) and isinstance(right, (int, float)):
            if isinstance(left, float) and not math.isfinite(left):
                if left != right:
                    self.note_exact(path, left, right)
                return
            if isinstance(right, float) and not math.isfinite(right):
                if left != right:
                    self.note_exact(path, left, right)
                return
            self.note_numeric(path, left, right)
            return
        if type(left) is not type(right):
            self.note_exact(path, left, right)
            return
        if isinstance(left, dict):
            left_keys = set(left)
            right_keys = set(right)
            for key in sorted(left_keys - right_keys):
                self.note_exact(f"{path}.{key}", left[key], "<missing>")
            for key in sorted(right_keys - left_keys):
                self.note_exact(f"{path}.{key}", "<missing>", right[key])
            for key in sorted(left_keys & right_keys):
                self.compare(left[key], right[key], f"{path}.{key}")
            return
        if isinstance(left, list):
            if len(left) != len(right):
                self.note_exact(f"{path}.length", len(left), len(right))
                return
            for index, (item_left, item_right) in enumerate(zip(left, right)):
                self.compare(item_left, item_right, f"{path}[{index}]")
            return
        if left != right:
            self.note_exact(path, left, right)

    def passed(self):
        return not self.differences and self.exact_difference_count == 0


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    with path.open() as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"{path} does not contain a JSON object")
    return value


def geometry_identity(data):
    """Return a fail-closed spatial identity, accepting only known old PBC files."""
    if "shear" not in data:
        if any(key in data for key in ("pbc_geometry", "torus_matrix",
                                       "geometry_id", "momentum_grid_convention")):
            raise ValueError("geometry metadata is incomplete: shear is missing")
        # Historical files predate integer shear and therefore mean rectangular PBC.
        Lx, Ly = data.get("Lx"), data.get("Ly")
        if not isinstance(Lx, int) or not isinstance(Ly, int):
            raise ValueError("cannot infer legacy rectangular geometry")
        return {"Lx": Lx, "Ly": Ly, "shear": 0,
                "pbc_geometry": "axis_aligned",
                "geometry_id": f"honeycomb-pbc-Lx{Lx}-Ly{Ly}-shear0"}
    shear = data["shear"]
    Lx, Ly = data.get("Lx"), data.get("Ly")
    if isinstance(shear, bool) or not isinstance(shear, int) or not isinstance(Lx, int) or not isinstance(Ly, int):
        raise ValueError("invalid geometry metadata types")
    if Lx < 2 or Ly < 2 or shear < 0 or shear >= Lx:
        raise ValueError("invalid sheared geometry dimensions")
    expected_kind = "axis_aligned" if shear == 0 else "sheared"
    expected_id = f"honeycomb-pbc-Lx{Lx}-Ly{Ly}-shear{shear}"
    if data.get("pbc_geometry") != expected_kind or data.get("geometry_id") != expected_id:
        raise ValueError("invalid or incomplete geometry identity")
    return {"Lx": Lx, "Ly": Ly, "shear": shear,
            "pbc_geometry": expected_kind, "geometry_id": expected_id}


def compare_json(left, right, tolerance):
    left_geometry = geometry_identity(left)
    right_geometry = geometry_identity(right)
    if left_geometry != right_geometry:
        raise ValueError("refusing to compare different spatial geometries: "
                         f"{left_geometry} vs {right_geometry}")
    # Invocation provenance and performance counters differ between engines.
    # Physical values, trajectories, acceptance counts and numerical-failure
    # diagnostics remain in the comparison.
    performance_diagnostics = {
        "elapsed_seconds", "timing", "cache_diagnostics",
        "prefix_resumes", "signed_factor_steps",
        "cayley_lu_calls", "cayley_cache_hits",
        "cayley_cache_fallbacks", "cayley_bytes_estimate",
    }
    allowed = {"defer_phase", "argv"} | performance_diagnostics
    left_checked = {key: value for key, value in left.items() if key not in allowed}
    right_checked = {key: value for key, value in right.items() if key not in allowed}
    ignored_suffixes = ("phase_diagnostics.denominator_snr",)
    comparison = Comparison(tolerance, ignored_suffixes)
    comparison.compare(left_checked, right_checked, "json")
    metadata = {
        "geometry": left_geometry,
        "allowed_differences": {
            "defer_phase": {"eager": left.get("defer_phase"), "deferred": right.get("defer_phase")},
            "elapsed_seconds": {"eager": left.get("elapsed_seconds"), "deferred": right.get("elapsed_seconds")},
            "performance_diagnostics": sorted(performance_diagnostics),
            "derived_fields": [
                "phase_diagnostics.denominator_snr",
                "vbs.blocking[*].phase_diagnostics.denominator_snr",
            ],
        },
        "schema_equal": set(left) == set(right),
        "status_complete": left.get("status") == "complete" and right.get("status") == "complete",
    }
    return comparison, metadata


def csv_number(column, value):
    if column in {"block", "count", "positive_count", "negative_count", "complex_count"}:
        return int(value)
    return float(value)


def compare_csv(left_path, right_path, tolerance):
    with left_path.open(newline="") as stream:
        left_reader = csv.DictReader(stream)
        left_header = left_reader.fieldnames
        left_rows = list(left_reader)
    with right_path.open(newline="") as stream:
        right_reader = csv.DictReader(stream)
        right_header = right_reader.fieldnames
        right_rows = list(right_reader)

    comparison = Comparison(tolerance)
    if left_header != right_header:
        comparison.note_exact("blocks.header", left_header, right_header)
    if len(left_rows) != len(right_rows):
        comparison.note_exact("blocks.row_count", len(left_rows), len(right_rows))
    for row_index, (left_row, right_row) in enumerate(zip(left_rows, right_rows)):
        for column in left_header or []:
            if column not in right_row:
                comparison.note_exact(f"blocks[{row_index}].{column}", left_row[column], "<missing>")
                continue
            try:
                comparison.note_numeric(
                    f"blocks[{row_index}].{column}",
                    csv_number(column, left_row[column]),
                    csv_number(column, right_row[column]),
                )
            except (TypeError, ValueError):
                comparison.note_exact(
                    f"blocks[{row_index}].{column}", left_row[column], right_row[column]
                )
    header = left_header or []
    required = {"block", "count", "sum_phase_real", "sum_phase_imag",
                "positive_count", "negative_count", "complex_count"}
    complex_columns = [column for column in header if column.startswith("sum_phase_")]
    metadata = {
        "header_equal": left_header == right_header,
        "row_count_equal": len(left_rows) == len(right_rows),
        "required_columns_present": required.issubset(header),
        "complex_sum_columns": len(complex_columns),
    }
    if not metadata["required_columns_present"]:
        comparison.note_exact("blocks.required_columns", sorted(required), sorted(set(header)))
    return comparison, metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--eager", type=pathlib.Path, required=True)
    parser.add_argument("--deferred", type=pathlib.Path, required=True)
    parser.add_argument("--eager-blocks", type=pathlib.Path, required=True)
    parser.add_argument("--deferred-blocks", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1.0e-12)
    args = parser.parse_args()

    eager = read_json(args.eager)
    deferred = read_json(args.deferred)
    json_comparison, json_metadata = compare_json(eager, deferred, args.tolerance)
    csv_comparison, csv_metadata = compare_csv(
        args.eager_blocks, args.deferred_blocks, args.tolerance
    )
    result = {
        "schema": "qmc_eager_deferred_comparison_v1",
        "inputs": {
            "eager_results": {"path": str(args.eager), "sha256": sha256(args.eager)},
            "deferred_results": {"path": str(args.deferred), "sha256": sha256(args.deferred)},
            "eager_blocks": {"path": str(args.eager_blocks), "sha256": sha256(args.eager_blocks)},
            "deferred_blocks": {"path": str(args.deferred_blocks), "sha256": sha256(args.deferred_blocks)},
        },
        "json": {
            **json_metadata,
            "passed": json_comparison.passed() and json_metadata["schema_equal"] and
            json_metadata["status_complete"],
            "numeric_difference_count": json_comparison.numeric_difference_count,
            "exact_difference_count": json_comparison.exact_difference_count,
            "max_absolute_error": json_comparison.max_abs,
            "max_normalized_error": json_comparison.max_normalized,
            "max_relative_error": json_comparison.max_relative,
            "first_difference": json_comparison.first_difference,
            "differences_over_tolerance": json_comparison.differences,
        },
        "blocks_csv": {
            **csv_metadata,
            "passed": csv_comparison.passed() and csv_metadata["header_equal"] and
            csv_metadata["row_count_equal"] and csv_metadata["required_columns_present"],
            "numeric_difference_count": csv_comparison.numeric_difference_count,
            "exact_difference_count": csv_comparison.exact_difference_count,
            "max_absolute_error": csv_comparison.max_abs,
            "max_normalized_error": csv_comparison.max_normalized,
            "max_relative_error": csv_comparison.max_relative,
            "first_difference": csv_comparison.first_difference,
            "differences_over_tolerance": csv_comparison.differences,
        },
    }
    result["passed"] = result["json"]["passed"] and result["blocks_csv"]["passed"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({
        "passed": result["passed"],
        "json_max_abs": result["json"]["max_absolute_error"],
        "json_max_normalized": result["json"]["max_normalized_error"],
        "blocks_max_abs": result["blocks_csv"]["max_absolute_error"],
        "blocks_max_normalized": result["blocks_csv"]["max_normalized_error"],
        "output": str(args.output),
    }, sort_keys=True))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"compare_qmc_outputs: {error}", file=sys.stderr)
        sys.exit(2)
