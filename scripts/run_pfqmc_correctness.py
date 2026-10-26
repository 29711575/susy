#!/usr/bin/env python3
"""Run the focused PfQMC validators and retain all output (stdlib only)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


GROUPS = {
    "energy": ["kitaev_ssh_total_energy_validate"],
    "gamma": ["kitaev_ssh_gamma_validate"],
    "physical": ["kitaev_ssh_physical_trial_validate"],
    "free": ["kitaev_ssh_pfqmc_free_phonon_validate"],
    "fock": ["kitaev_ssh_pfqmc_fock_validate"],
    "cache": ["kitaev_ssh_pfqmc_cache_validate"],
    "regression": ["kitaev_ssh_phonon_validate",
                   "kitaev_ssh_mc_validate",
                   "kitaev_ssh_vbs_susceptibility_validate",
                   "kitaev_ssh_numerical_failure_validate",
                   "kitaev_ssh_green_refinement_validate",
                   "kitaev_ssh_fermion_validate"],
    "optimized": ["kitaev_ssh_real_magnitude_validate", "kitaev_ssh_bipartite_validate"],
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="new directory; existing evidence is never overwritten")
    parser.add_argument("--group", choices=["all"] + list(GROUPS), default="all")
    args = parser.parse_args()
    build = args.build_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ)
    for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    env["MKL_DYNAMIC"] = "FALSE"
    names = [name for group in GROUPS.values() for name in group] if args.group == "all" else GROUPS[args.group]
    results = []
    for name in names:
        binary = build / name
        command = [str(binary)]
        if name == "kitaev_ssh_numerical_failure_validate":
            command.append(str(output / "failure-cli"))
        begin = time.monotonic()
        record = {"name": name, "command": command}
        try:
            record["binary_sha256"] = hashlib.sha256(binary.read_bytes()).hexdigest()
            with (output / (name + ".stdout")).open("w") as out, (output / (name + ".stderr")).open("w") as err:
                proc = subprocess.run(command, stdout=out, stderr=err, env=env)
            record["exit_code"] = proc.returncode
            record["status"] = "PASS" if proc.returncode == 0 else "FAIL"
            if name == "kitaev_ssh_numerical_failure_validate" and proc.returncode == 0:
                expected = {"recover": "complete", "abort": "aborted",
                            "measure_recover": "complete", "measure_abort": "aborted",
                            "initial_recover": "complete", "initial_redraw": "complete"}
                record["cli_statuses"] = {}
                for case, state in expected.items():
                    result = json.loads((output / "failure-cli" / (case + "_results.json")).read_text())
                    record["cli_statuses"][case] = result["status"]
                    if result["status"] != state:
                        raise RuntimeError("unexpected CLI status for " + case)
        except Exception as error:
            record["status"] = "FAIL"
            record["error"] = str(error)
        record["seconds"] = time.monotonic() - begin
        results.append(record)
        print("{} {} ({:.2f}s)".format(record["status"], name, record["seconds"]), flush=True)
        summary = {"status": "PASS" if all(r["status"] == "PASS" for r in results) else "FAIL",
                   "completed": len(results), "planned": len(names), "tests": results}
        (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0 if summary["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
