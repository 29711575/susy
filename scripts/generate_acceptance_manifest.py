#!/usr/bin/env python3
"""Create a machine-readable manifest for the 2026-09-20 acceptance evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(root: Path, path: Path) -> dict[str, object]:
    return {
        "path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def command_line(*args: str) -> str:
    try:
        return subprocess.run(args, check=True, text=True, capture_output=True).stdout.splitlines()[0]
    except (OSError, subprocess.CalledProcessError, IndexError):
        return "unavailable"


def parse_time(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    if not path.exists() or not path.read_text().strip():
        return values
    for line in path.read_text().splitlines():
        key, separator, value = line.partition("=")
        if separator:
            try:
                values[key] = float(value)
            except ValueError:
                pass
    return values


def result_record(root: Path, raw: Path, result_path: Path) -> dict[str, object]:
    result = json.loads(result_path.read_text())
    prefix = result_path.name.removesuffix("_results.json")
    time_path = raw / f"{prefix}.time"
    exit_path = raw / f"{prefix}.exit"
    blocks_path = raw / f"{prefix}_blocks.csv"
    phase = result.get("phase_diagnostics", {})
    exit_code = None
    if exit_path.exists() and exit_path.read_text().strip():
        exit_code = int(exit_path.read_text().strip())
    return {
        "prefix": prefix,
        "status": result.get("status"),
        "engine": result.get("weight_engine"),
        "defer_phase": result.get("defer_phase"),
        "seed": result.get("seed"),
        "g": result.get("g"),
        "lambda": result.get("lambda"),
        "Lx": result.get("Lx"),
        "Ly": result.get("Ly"),
        "intervals": result.get("intervals"),
        "dt": result.get("dt"),
        "warmup": result.get("warmup"),
        "sweeps": result.get("sweeps"),
        "measure_every": result.get("measure_every"),
        "block_samples": result.get("block_samples"),
        "samples": result.get("samples"),
        "blocks": result.get("blocks"),
        "acceptance": result.get("acceptance"),
        "elapsed_seconds": result.get("elapsed_seconds"),
        "timing": parse_time(time_path),
        "exit_code": exit_code,
        "numerical_failures": result.get("numerical_failures"),
        "ratio_defined": result.get("ratio_defined"),
        "jackknife_defined": result.get("jackknife_defined"),
        "block_independence_verified": result.get("block_independence_verified"),
        "phase_diagnostics": {
            "positive_count": phase.get("positive_count"),
            "negative_count": phase.get("negative_count"),
            "complex_count": phase.get("complex_count"),
            "reality_tolerance": phase.get("reality_tolerance"),
        },
        "evidence_files": [
            file_record(root, result_path),
            file_record(root, blocks_path) if blocks_path.exists() else None,
            file_record(root, time_path) if time_path.exists() else None,
            file_record(root, exit_path) if exit_path.exists() else None,
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    audit = root / "audit" / "acceptance-20260920"
    raw = audit / "raw_hpc"
    output = (args.output or audit / "manifest.json").resolve()

    result_paths = sorted(raw.glob("*_results.json"))
    jobs = [result_record(root, raw, path) for path in result_paths]
    incomplete_snapshot = raw / "long_reference_g.6_seed8301_defer0_blocks.csv"
    ctest_log = audit / "ctest-with-mkl-runtime.log"
    ctest_text = ctest_log.read_text() if ctest_log.exists() else ""
    passed = len(re.findall(r"Test\s+#\d+:.*Passed", ctest_text))
    total_match = re.search(r"100% tests passed out of (\d+)", ctest_text)
    total = int(total_match.group(1)) if total_match else None
    comparison_path = audit / "combined_eager_deferred_comparison.json"
    comparison = json.loads(comparison_path.read_text()) if comparison_path.exists() else None

    files = sorted(path for path in raw.rglob("*") if path.is_file())
    local_names = [
        "configure.log",
        "build.log",
        "ctest.log",
        "ctest-with-mkl-runtime.log",
        "mklroot-cache-regression.log",
        "deferred_phase_validate.stdout",
        "deferred_phase_validate.stderr",
        "a3_validate.stdout",
        "a3_validate.stderr",
        "weight_validate.stdout",
        "weight_validate.stderr",
        "vbs_q_validate.stdout",
        "vbs_q_validate.stderr",
        "combined_eager_deferred_comparison.json",
    ]
    local_files = [
        file_record(root, audit / name)
        for name in local_names
        if (audit / name).exists()
    ]

    comparison_times = {}
    if comparison:
        comparison_times = comparison.get("json", {}).get("allowed_differences", {})
    eager_seconds = comparison_times.get("elapsed_seconds", {}).get("eager")
    deferred_seconds = comparison_times.get("elapsed_seconds", {}).get("deferred")
    speedup = eager_seconds / deferred_seconds if eager_seconds and deferred_seconds else None
    reduction = 1.0 - deferred_seconds / eager_seconds if eager_seconds and deferred_seconds else None

    branch = command_line("git", "branch", "--show-current")
    head = command_line("git", "rev-parse", "HEAD")
    manifest = {
        "schema": "qmc_acceptance_manifest_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repository": {
            "root": str(root),
            "branch": branch,
            "head": head,
            "worktree_is_dirty": bool(subprocess.run(["git", "status", "--porcelain"], cwd=root, text=True, capture_output=True).stdout.strip()),
        },
        "build": {
            "directory": "build-acceptance-20260920",
            "cmake": "/home/asus/code/driven-new-validation-20260908/green-residual/tools/cmake/data/bin/cmake",
            "cmake_version": command_line("/home/asus/code/driven-new-validation-20260908/green-residual/tools/cmake/data/bin/cmake", "--version"),
            "ctest": "/home/asus/code/driven-new-validation-20260908/green-residual/tools/cmake/data/bin/ctest",
            "compiler": "/usr/bin/c++",
            "compiler_version": command_line("/usr/bin/c++", "--version"),
            "mkl_root": "/opt/intel/oneapi/mkl/2026.1",
            "mkl_libraries": [
                "/opt/intel/oneapi/mkl/2026.1/lib/libmkl_intel_lp64.so.3",
                "/opt/intel/oneapi/mkl/2026.1/lib/libmkl_sequential.so.3",
                "/opt/intel/oneapi/mkl/2026.1/lib/libmkl_core.so.3",
            ],
            "required_runtime_environment": {
                "CPATH": "/opt/intel/oneapi/mkl/2026.1/include",
                "LD_LIBRARY_PATH_prefix": "/opt/intel/oneapi/mkl/2026.1/lib",
                "thread_variables": "OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE",
            },
        },
        "local_tests": {
            "ctest_log": "audit/acceptance-20260920/ctest-with-mkl-runtime.log",
            "passed": passed,
            "total": total,
            "all_passed": total is not None and passed == total,
            "first_run_without_runtime_path": "audit/acceptance-20260920/ctest.log",
        },
        "hpc": {
            "connection_port": 20022,
            "job_id": "84857[]",
            "remote_path": "/home/sunxr/pure-new-faster-perf-testq-20260920",
            "remote_query_status": "copied_evidence_available_ssh_requery_denied",
            "no_new_job_submitted_during_acceptance": True,
            "complete_jobs": jobs,
            "incomplete_reference_eager": {
                "status": "snapshot_only",
                "reason": "missing final .time, .exit, and results.json",
                "evidence_files": [file_record(root, incomplete_snapshot)] if incomplete_snapshot.exists() else [],
            },
            "raw_files": [file_record(root, path) for path in files],
        },
        "comparison": {
            "path": "audit/acceptance-20260920/combined_eager_deferred_comparison.json",
            "passed": comparison.get("passed") if comparison else False,
            "json_max_absolute_error": comparison.get("json", {}).get("max_absolute_error") if comparison else None,
            "blocks_max_absolute_error": comparison.get("blocks_csv", {}).get("max_absolute_error") if comparison else None,
            "combined_eager_seconds": eager_seconds,
            "combined_deferred_seconds": deferred_seconds,
            "deferred_speedup": speedup,
            "deferred_wall_time_reduction": reduction,
        },
        "local_evidence_files": local_files,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(output)


if __name__ == "__main__":
    main()
