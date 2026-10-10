#!/usr/bin/env python3
"""Unit tests for Stage 1.1's independent timing and output comparators."""

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "perf_work" / "scripts"))
sys.path.insert(0, str(ROOT / "scripts"))

from stage1_time import TimeParseError, parse_time_text
from compare_stage1_outputs import Comparison, compare_json
from compare_qmc_outputs import compare_json as compare_qmc_json


def check_time(value: str, expected: float) -> None:
    parsed = parse_time_text(
        f"Elapsed (wall clock) time: {value}\n"
        "Maximum resident set size (kbytes): 42\n"
    )
    assert parsed["elapsed_seconds"] == expected, (value, parsed)


def check_verbose_time(value: str, expected: float) -> None:
    parsed = parse_time_text(
        f"Elapsed (wall clock) time (h:mm:ss or m:ss): {value}\n"
        "Maximum resident set size (kbytes): 42\n"
    )
    assert parsed["elapsed_seconds"] == expected, (value, parsed)


def check_comparison(left, right, should_pass: bool, ignored=None) -> Comparison:
    comparison = Comparison(ignored_paths=ignored or set())
    compare_json(left, right, comparison)
    assert (comparison.first_difference is None) == should_pass, comparison.first_difference
    return comparison


def main() -> None:
    check_time("0:12.34", 12.34)
    check_time("1:23.13", 83.13)
    check_time("1:02:03.45", 3723.45)
    check_time("1:00.00", 60.0)
    check_time("12.5", 12.5)
    check_verbose_time("0:12.34", 12.34)
    check_verbose_time("1:23.13", 83.13)
    check_verbose_time("1:02:03.45", 3723.45)
    for bad in ("", "1:60.0", "1:02:60.0", "1:2:3:4", "nan", "-1.0"):
        try:
            check_time(bad, 0.0)
        except TimeParseError:
            pass
        else:
            raise AssertionError(f"accepted invalid elapsed time {bad!r}")
    try:
        parse_time_text("Maximum resident set size (kbytes): 42\n")
    except TimeParseError:
        pass
    else:
        raise AssertionError("accepted missing elapsed field")

    check_comparison({"x": 7}, {"x": 7}, True)
    check_comparison({"x": 0.5}, {"x": 0.5}, True)
    check_comparison({"x": -0.0}, {"x": 0.0}, True)
    for left, right in (
        (9007199254740992, 9007199254740993),
        (18446744073709551614, 18446744073709551615),
        (True, 1),
        (False, 0),
    ):
        check_comparison({"x": left}, {"x": right}, False)
    check_comparison({"x": 1}, {"x": 1.0}, True)
    check_comparison({"x": 9007199254740992}, {"x": 9007199254740992.0}, True)
    check_comparison({"x": 1}, {"x": 1.5}, False)
    check_comparison({"x": float("nan")}, {"x": float("nan")}, True)
    check_comparison({"x": float("inf")}, {"x": float("inf")}, True)
    check_comparison({"x": float("inf")}, {"x": float("-inf")}, False)
    check_comparison({"x": [1, 2]}, {"x": [1]}, False)
    check_comparison({"x": {"y": 1}}, {"x": {"y": 2}}, False)
    check_comparison({"timing": {"x": 1}}, {"timing": {"x": 2}}, True,
                     ignored={("json", "timing")})
    qmc_base = {
        "status": "complete", "Lx": 3, "Ly": 2, "shear": 1,
        "pbc_geometry": "sheared",
        "geometry_id": "honeycomb-pbc-Lx3-Ly2-shear1",
        "defer_phase": 0, "elapsed_seconds": 1.0,
        "timing": {"sampling_seconds": 0.8},
        "cache_diagnostics": {"local_proposals": 4},
        "prefix_resumes": 3, "signed_factor_steps": 7,
        "cayley_lu_calls": 7, "cayley_cache_hits": 0,
        "cayley_cache_fallbacks": 0, "cayley_bytes_estimate": 0,
        "argv": ["mc", "--defer-phase", "0"], "observable": 0.25,
    }
    qmc_other = dict(qmc_base)
    qmc_other.update({
        "defer_phase": 1, "elapsed_seconds": 2.0,
        "timing": {"sampling_seconds": 1.7},
        "cache_diagnostics": {"local_proposals": 8},
        "prefix_resumes": 0, "signed_factor_steps": 11,
        "cayley_lu_calls": 1, "cayley_cache_hits": 6,
        "cayley_cache_fallbacks": 1, "cayley_bytes_estimate": 4096,
        "argv": ["mc", "--defer-phase", "1"],
    })
    qmc_comparison, metadata = compare_qmc_json(qmc_base, qmc_other, 1e-12)
    assert qmc_comparison.passed() and metadata["status_complete"]
    qmc_other["observable"] = 0.5
    qmc_comparison, _ = compare_qmc_json(qmc_base, qmc_other, 1e-12)
    assert not qmc_comparison.passed()
    print("STAGE1_TOOLS PASS")


if __name__ == "__main__":
    main()
