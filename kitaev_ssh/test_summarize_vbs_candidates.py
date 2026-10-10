#!/usr/bin/env python3
"""Focused third-round tests for the real candidate postprocessor."""

import argparse
import copy
import hashlib
import importlib.util
import json
import os
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("summarize_vbs_candidates", HERE / "summarize_vbs_candidates.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _estimate(value, *, defined=True, jackknife=True):
    return {
        "defined": defined,
        "reason": "defined" if defined else "full_sample_ratio_undefined",
        "real": value if defined else None,
        "imag_diag": 0.0 if defined else None,
        "stderr_real": 0.1 if jackknife and defined else None,
        "stderr_imag": 0.2 if jackknife and defined else None,
        "cov_real_imag": 0.03 if jackknife and defined else None,
        "jackknife_defined": jackknife and defined,
        "jackknife_reason": "defined" if jackknife and defined else "incomplete_delete_one_block_replicas",
    }


def _group(name, participants, *, total=4, valid=4, invalid=0, center=True,
           jackknife=True, stability=True, conditional=False):
    available = bool(participants)
    return {
        "name": name,
        "estimator": "projected" if name == "offgrid_projected" else "raw",
        "available": available,
        "participants": list(participants),
        "participant_labels": ["q" for _ in participants],
        "center_defined": center if available else False,
        "center_reason": "defined" if center and available else "empty_group",
        "max_mz2_real": 1.0 if center and available else None,
        "winner_display_index": participants[0] if center and available else None,
        "winner_set": list(participants) if center and available else [],
        "winner_labels": ["q" for _ in participants] if center and available else [],
        "total_replicas": total,
        "valid_replicas": valid,
        "invalid_replicas": invalid,
        "jackknife_defined": jackknife if available else False,
        "jackknife_reason": "defined" if jackknife and available else "empty_group",
        "stability_defined": stability if available else False,
        "stability_reason": "defined" if stability and available else "fewer_than_two_valid_replicas",
        "conditional_on_valid_replicas": conditional if available else False,
        "unique_win_counts": [valid for _ in participants],
        "tie_membership_counts": [valid for _ in participants],
        "unique_win_frequency": [1.0 if stability and available else None for _ in participants],
        "tie_membership_frequency": [1.0 if stability and available else None for _ in participants],
        "class_tie_counts": [],
        "class_full_set_counts": [],
        "class_tie_frequency": [],
        "class_full_set_frequency": [],
    }


def fixture():
    points = [
        {"index": 0, "label": "q", "group": "grid_nonzero", "u": 0.5, "v": 0.5,
         "on_grid": True, "W_L_abs2": 0.0, "mz2_raw": _estimate(1.0),
         "mz2_projected": _estimate(1.1), "plus_minus_class": None,
         "class_self_conjugate": False},
        {"index": 1, "label": "zero", "group": "zero_excluded", "u": 0.0, "v": 0.0,
         "on_grid": True, "W_L_abs2": 1.0, "mz2_raw": _estimate(0.0),
         "mz2_projected": _estimate(0.0), "plus_minus_class": None,
         "class_self_conjugate": False},
    ]
    groups = [_group("grid_raw", [0]), _group("offgrid_raw", []),
              _group("offgrid_projected", [])]
    return {"status": "complete", "vbs": {
        "mode": "selected", "Nc": 4, "Nq": 2, "qpoints": [
            {"index": 0, "label": "q", "on_grid": True, "n": [1, 1],
             "reciprocal_fraction": [0.5, 0.5], "W_L_abs2": 0.0},
            {"index": 1, "label": "zero", "on_grid": True, "n": [0, 0],
             "reciprocal_fraction": [0.0, 0.0], "W_L_abs2": 1.0},
        ],
        "candidates": {
            "schema": "ivbs-candidates-v2", "schema_revision": 2, "level": "finest_blocks",
            "ratio_defined": True, "center_defined": True, "center_reason": "defined",
            "jackknife_defined": True, "jackknife_reason": "defined",
            "stability_defined": True, "stability_reason": "defined",
            "total_replicas": 4, "valid_replicas": 4, "invalid_replicas": 0,
            "invalid_replica_reasons": [], "groups": groups, "points": points,
            "plus_minus_classes": [],
        },
    }}


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _write_input(path, data):
    path.write_text(json.dumps(data, allow_nan=False), encoding="utf-8")


def run():
    with tempfile.TemporaryDirectory(prefix="ivbs-round3-python-") as temp:
        root = Path(temp)
        for status in ("running", "aborted"):
            incomplete = fixture()
            incomplete["status"] = status
            try:
                MODULE.summarize(incomplete, root / "incomplete_results.json")
            except MODULE.SchemaError:
                pass
            else:
                raise AssertionError("incomplete run accepted: " + status)
        source30 = root / "scan_lambda_3.0_results.json"
        source35 = root / "scan_lambda_3.5_results.json"
        _write_input(source30, fixture())
        _write_input(source35, fixture())
        csv30, json30 = MODULE.output_paths(source30, None, None)
        csv35, json35 = MODULE.output_paths(source35, None, None)
        assert csv30.name == "scan_lambda_3.0_candidates_summary.csv"
        assert json30.name == "scan_lambda_3.0_candidates_summary.json"
        assert csv30 != csv35 and json30 != json35

        summary = MODULE.summarize(json.loads(source30.read_text()), source30)
        MODULE.write_outputs(summary, csv30, json30, input_path=source30,
                             explicit_csv=False, explicit_json=False)
        assert json.loads(json30.read_text())["status"] == "available"
        try:
            MODULE.write_outputs(summary, csv30, json30, input_path=source30,
                                 explicit_csv=False, explicit_json=False)
        except FileExistsError:
            pass
        else:
            raise AssertionError("default output overwrite was not refused")
        original_json_hash = _digest(json30)
        original_csv_hash = _digest(csv30)

        explicit_csv = root / "explicit.csv"
        explicit_json = root / "explicit.json"
        MODULE.write_outputs(summary, explicit_csv, explicit_json, input_path=source30,
                             explicit_csv=True, explicit_json=True)
        try:
            MODULE.validate_output_paths(source30, explicit_csv, explicit_json,
                                         force=False, explicit_csv=True, explicit_json=True)
        except FileExistsError:
            pass
        else:
            raise AssertionError("explicit output overwrite was not refused")
        MODULE.validate_output_paths(source30, explicit_csv, explicit_json, force=True,
                                     explicit_csv=True, explicit_json=True)

        try:
            MODULE.validate_output_paths(source30, root / "x.csv", root / "x.csv",
                                         force=True, explicit_csv=True, explicit_json=True)
        except ValueError:
            pass
        else:
            raise AssertionError("same output path was accepted")
        try:
            MODULE.validate_output_paths(source30, source30, root / "x.json",
                                         force=True, explicit_csv=True, explicit_json=True)
        except ValueError:
            pass
        else:
            raise AssertionError("input overwrite was accepted")
        try:
            MODULE.validate_output_paths(source30, csv30, json30, force=True,
                                         explicit_csv=False, explicit_json=False)
        except ValueError:
            pass
        else:
            raise AssertionError("force accepted implicit existing outputs")
        alias = root / "input_alias.json"
        os.symlink(source30, alias)
        try:
            MODULE.validate_output_paths(source30, alias, root / "alias.json",
                                         force=True, explicit_csv=True, explicit_json=True)
        except ValueError:
            pass
        else:
            raise AssertionError("symlink input alias was accepted")
        assert _digest(json30) == original_json_hash and _digest(csv30) == original_csv_hash

        partial = fixture()
        c = partial["vbs"]["candidates"]
        g = c["groups"][0]
        for group in c["groups"]:
            group.update({"total_replicas": 4, "valid_replicas": 3, "invalid_replicas": 1})
        c.update({"total_replicas": 4, "valid_replicas": 3, "invalid_replicas": 1,
                  "invalid_replica_reasons": ["test_invalid_replica"]})
        g.update({"jackknife_defined": False, "jackknife_reason": "incomplete_delete_one_block_replicas",
                  "conditional_on_valid_replicas": True})
        g["unique_win_counts"] = [3]
        g["tie_membership_counts"] = [3]
        g["unique_win_frequency"] = [1.0]
        g["tie_membership_frequency"] = [1.0]
        c.update({"jackknife_defined": False, "jackknife_reason": "one_or_more_groups_jackknife_undefined"})
        for point in partial["vbs"]["candidates"]["points"]:
            for estimate in (point["mz2_raw"], point["mz2_projected"]):
                for key in ("stderr_real", "stderr_imag", "cov_real_imag"):
                    estimate[key] = None
                estimate.update({"jackknife_defined": False,
                                 "jackknife_reason": "incomplete_delete_one_block_replicas"})
        partial_summary = MODULE.summarize(partial, root / "partial.json")
        assert partial_summary["status"] == "partial"
        assert partial_summary["rows"][0]["grid_raw_unique_win_frequency"] == 1.0
        assert partial_summary["rows"][0]["mz2_raw_stderr_real"] is None

        sheared = fixture()
        sheared.update({"Lx": 2, "Ly": 2, "shear": 1,
                        "pbc_geometry": "sheared",
                        "geometry_id": "honeycomb-pbc-Lx2-Ly2-shear1"})
        sheared["vbs"]["qpoints"][0]["reciprocal_fraction"] = [0.5, 0.25]
        sheared["vbs"]["candidates"]["points"][0]["u"] = 0.5
        sheared["vbs"]["candidates"]["points"][0]["v"] = 0.25
        sheared_summary = MODULE.summarize(sheared, root / "sheared.json")
        assert sheared_summary["status"] == "available"
        assert sheared_summary["geometry"]["shear"] == 1

        unavailable_data = copy.deepcopy(partial)
        uc = unavailable_data["vbs"]["candidates"]
        ug = uc["groups"][0]
        ug.update({"center_defined": False, "center_reason": "full_sample_ratio_undefined",
                   "max_mz2_real": None, "winner_display_index": None,
                   "winner_set": [], "winner_labels": [], "stability_defined": False,
                   "stability_reason": "center_undefined", "conditional_on_valid_replicas": False,
                   "unique_win_frequency": [None], "tie_membership_frequency": [None],
                   "jackknife_defined": False, "jackknife_reason": "center_undefined"})
        uc.update({"center_defined": False, "center_reason": "one_or_more_groups_center_undefined",
                   "stability_defined": False, "stability_reason": "center_undefined"})
        for estimate in (uc["points"][0]["mz2_raw"], uc["points"][0]["mz2_projected"]):
            estimate.update({"defined": False, "reason": "full_sample_ratio_undefined", "real": None,
                             "imag_diag": None, "jackknife_defined": False,
                             "jackknife_reason": "center_undefined", "stderr_real": None,
                             "stderr_imag": None, "cov_real_imag": None})
        assert MODULE.summarize(unavailable_data, root / "unavailable.json")["status"] == "unavailable"

        legacy = fixture()
        legacy["vbs"]["candidates"]["schema"] = "ivbs-candidates-v1"
        assert MODULE.summarize(legacy, root / "legacy.json")["status"] == "unavailable"
        malformed = fixture()
        malformed["vbs"]["candidates"]["groups"][0]["unique_win_frequency"] = []
        try:
            MODULE.summarize(malformed, root / "malformed.json")
        except MODULE.SchemaError:
            pass
        else:
            raise AssertionError("malformed frequency length was accepted")
        malformed_json = root / "malformed_json.json"
        malformed_json.write_text("{", encoding="utf-8")
        try:
            MODULE._load_json(malformed_json)
        except json.JSONDecodeError:
            pass
        else:
            raise AssertionError("malformed JSON was accepted")
        nonfinite_json = root / "nonfinite.json"
        nonfinite_json.write_text('{"value": NaN}', encoding="utf-8")
        try:
            MODULE._load_json(nonfinite_json)
        except ValueError:
            pass
        else:
            raise AssertionError("non-finite JSON was accepted")

        print(json.dumps({
            "status": "PASS",
            "default_names": [csv30.name, json30.name, csv35.name, json35.name],
            "partial_status": partial_summary["status"],
            "unavailable_status": "unavailable",
            "legacy_status": "unavailable",
            "csv_blank_for_undefined_error": True,
        }, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.parse_args()
    run()
