#!/usr/bin/env python3
"""Compare a baseline/current coupled-MC pair without hiding trajectories."""

import argparse
import csv
import json
import math
from pathlib import Path


DEFAULT_IGNORED_JSON_PATHS = {
    ("json", "elapsed_seconds"),
    ("json", "timing"),
    ("json", "cache_diagnostics"),
}


class Comparison:
    def __init__(self, ignored_paths=None):
        self.ignored_paths = set(DEFAULT_IGNORED_JSON_PATHS if ignored_paths is None else ignored_paths)
        self.numeric_difference_count = 0
        self.exact_difference_count = 0
        self.max_absolute_error = 0.0
        self.max_relative_error = 0.0
        self.first_difference = None

    def difference(self, path, left, right, absolute=None, relative=None):
        if absolute is None:
            self.exact_difference_count += 1
        else:
            self.numeric_difference_count += absolute != 0.0
            self.max_absolute_error = max(self.max_absolute_error, absolute)
            self.max_relative_error = max(self.max_relative_error, relative)
        if self.first_difference is None:
            entry = {"path": path, "baseline": left, "current": right}
            if type(left) is not type(right):
                entry["baseline_type"] = type(left).__name__
                entry["current_type"] = type(right).__name__
            if absolute is not None:
                entry.update({"absolute_error": absolute, "relative_error": relative})
            self.first_difference = entry


def _numeric_error(left, right):
    absolute = abs(left - right)
    if isinstance(left, int) and isinstance(right, int):
        scale = max(abs(left), abs(right), 1)
        try:
            relative = float(absolute) / float(scale)
        except OverflowError:
            relative = 1.0 if absolute else 0.0
        return absolute, relative
    scale = max(abs(float(left)), abs(float(right)), 1.0e-300)
    return absolute, float(absolute) / scale


def _mixed_numeric_equal(integer, floating):
    """Use exact integer equality for a finite integral binary64 value."""
    if not math.isfinite(floating) or not floating.is_integer():
        return False
    return int(floating) == integer


def compare_json(left, right, comparison, path="json"):
    if isinstance(left, bool) or isinstance(right, bool):
        if type(left) is not type(right) or left != right:
            comparison.difference(path, left, right)
        return
    if isinstance(left, int) and isinstance(right, int):
        if left != right:
            absolute, relative = _numeric_error(left, right)
            comparison.difference(path, left, right, absolute, relative)
        return
    if isinstance(left, float) and isinstance(right, float):
        if math.isnan(left) or math.isnan(right):
            if not (math.isnan(left) and math.isnan(right)):
                comparison.difference(path, left, right)
            return
        if not math.isfinite(left) or not math.isfinite(right):
            if left != right:
                comparison.difference(path, left, right)
            return
        absolute = abs(left - right)
        if absolute:
            comparison.difference(
                path,
                left,
                right,
                absolute,
                absolute / max(abs(float(left)), abs(float(right)), 1.0e-300),
            )
        return
    if ((isinstance(left, int) and isinstance(right, float)) or
            (isinstance(left, float) and isinstance(right, int))):
        # Mixed int/float values pass only when the finite integral binary64
        # value is exactly the same integer. This makes the policy explicit
        # and avoids converting large integer/integer pairs through binary64.
        integer, floating = (left, right) if isinstance(left, int) else (right, left)
        if _mixed_numeric_equal(integer, floating):
            return
        absolute, relative = _numeric_error(left, right)
        comparison.difference(path, left, right, absolute, relative)
        return
    if isinstance(left, dict) and isinstance(right, dict):
        keys = set(left) | set(right)
        for key in sorted(keys):
            child_path = path + "." + key
            if tuple(child_path.split(".")) in comparison.ignored_paths:
                continue
            if key not in left or key not in right:
                comparison.difference(child_path, left.get(key, "<missing>"), right.get(key, "<missing>"))
            else:
                compare_json(left[key], right[key], comparison, child_path)
        return
    if isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            comparison.difference(path + ".length", len(left), len(right))
            return
        for index, (left_item, right_item) in enumerate(zip(left, right)):
            compare_json(left_item, right_item, comparison, f"{path}[{index}]")
        return
    if left != right:
        comparison.difference(path, left, right)


def compare_csv(left_path, right_path, comparison):
    with left_path.open(newline="") as left_stream:
        left_rows = list(csv.reader(left_stream))
    with right_path.open(newline="") as right_stream:
        right_rows = list(csv.reader(right_stream))
    if left_rows == right_rows:
        return
    for row_index, (left_row, right_row) in enumerate(zip(left_rows, right_rows)):
        if left_row != right_row:
            comparison.difference(f"blocks.csv[{row_index}]", left_row, right_row)
            return
    comparison.difference("blocks.csv.length", len(left_rows), len(right_rows))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--baseline-blocks", type=Path, required=True)
    parser.add_argument("--current-blocks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    comparison = Comparison()
    with args.baseline.open() as stream:
        baseline = json.load(stream)
    with args.current.open() as stream:
        current = json.load(stream)
    compare_json(baseline, current, comparison)
    compare_csv(args.baseline_blocks, args.current_blocks, comparison)
    result = {
        "status": "PASS" if comparison.first_difference is None else "FAIL",
        "ignored_json_paths": [".".join(path) for path in sorted(comparison.ignored_paths)],
        "numeric_difference_count": comparison.numeric_difference_count,
        "exact_difference_count": comparison.exact_difference_count,
        "max_absolute_error": comparison.max_absolute_error,
        "max_relative_error": comparison.max_relative_error,
        "first_difference": comparison.first_difference,
        "baseline": str(args.baseline),
        "current": str(args.current),
        "baseline_blocks": str(args.baseline_blocks),
        "current_blocks": str(args.current_blocks),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: result[k] for k in (
        "status", "numeric_difference_count", "exact_difference_count",
        "max_absolute_error", "max_relative_error", "first_difference")}, sort_keys=True))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
