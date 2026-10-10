#!/usr/bin/env python3
"""Validate and summarize one coupled_mc ``vbs.candidates`` object.

The C++ candidate diagnostics use ``ivbs-candidates-v2``. This reader keeps
center, complete-jackknife, and conditional selection-stability states
separate. It never pools input files and never turns an unavailable number
into zero.
"""

import argparse
import csv
import io
import json
import math
import os
import tempfile
from pathlib import Path


SCHEMA = "ivbs-candidates-summary-v2"
CANDIDATE_SCHEMA = "ivbs-candidates-v2"
GROUPS = ("grid_raw", "offgrid_raw", "offgrid_projected")
ESTIMATORS = {"grid_raw": "raw", "offgrid_raw": "raw", "offgrid_projected": "projected"}


class SchemaError(ValueError):
    """The input is JSON, but not a supported candidate diagnostic."""


def number(value):
    """Return a finite JSON number, otherwise None for a legal null field."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def labels(value):
    if not isinstance(value, list) or not all(isinstance(x, str) for x in value):
        return None
    return ";".join(value)


def parameter_subset(data):
    keys = ("Lx", "Ly", "shear", "pbc_geometry", "geometry_id", "Jx", "Jy", "Jz", "g", "mass", "omega_ph",
            "lambda", "intervals", "dt", "tau", "seed", "proposal")
    return {key: data[key] for key in keys if key in data}


def _canonical(path):
    return Path(path).expanduser().resolve(strict=False)


def _lexists(path):
    return os.path.lexists(os.fspath(path))


def _same_file_or_path(left, right):
    if _canonical(left) == _canonical(right):
        return True
    if _lexists(left) and _lexists(right):
        try:
            return os.path.samefile(left, right)
        except OSError:
            return False
    return False


def normalize_geometry(data):
    """Accept known legacy rectangular output and reject incomplete sheared metadata."""
    if "Lx" not in data or "Ly" not in data:
        return None
    if "shear" not in data:
        if any(key in data for key in ("pbc_geometry", "torus_matrix", "geometry_id",
                                       "momentum_grid_convention")):
            raise SchemaError("geometry metadata is incomplete: missing shear")
        data["shear"] = 0
        data["pbc_geometry"] = "axis_aligned"
        data["geometry_id"] = f"honeycomb-pbc-Lx{data['Lx']}-Ly{data['Ly']}-shear0"
        return {key: data[key] for key in ("Lx", "Ly", "shear", "pbc_geometry", "geometry_id")}
    shear = data["shear"]
    if isinstance(shear, bool) or not isinstance(shear, int) or not 0 <= shear < data["Lx"]:
        raise SchemaError("invalid shear metadata")
    kind = "axis_aligned" if shear == 0 else "sheared"
    identity = f"honeycomb-pbc-Lx{data['Lx']}-Ly{data['Ly']}-shear{shear}"
    if data.get("pbc_geometry") != kind or data.get("geometry_id") != identity:
        raise SchemaError("invalid or incomplete geometry identity")
    return {"Lx": data["Lx"], "Ly": data["Ly"], "shear": shear,
            "pbc_geometry": kind, "geometry_id": identity}


def output_paths(input_path, csv_path, json_path):
    """Return output paths without interpreting dots in the input stem."""
    input_path = Path(input_path)
    name = input_path.name
    stem = name[:-len(".json")] if name.endswith(".json") else input_path.stem
    if stem.endswith("_results"):
        stem = stem[:-len("_results")]
    base = input_path.with_name(stem + "_candidates_summary")
    csv_output = Path(csv_path) if csv_path is not None else Path(str(base) + ".csv")
    json_output = Path(json_path) if json_path is not None else Path(str(base) + ".json")
    if _same_file_or_path(csv_output, json_output):
        raise ValueError("CSV and JSON outputs resolve to the same file")
    if _same_file_or_path(input_path, csv_output) or _same_file_or_path(input_path, json_output):
        raise ValueError("summary output cannot be the input file")
    return csv_output, json_output


def validate_output_paths(input_path, csv_path, json_path, *, force=False,
                          explicit_csv=False, explicit_json=False):
    """Check aliases and overwrite policy before creating any output."""
    if _same_file_or_path(csv_path, json_path):
        raise ValueError("CSV and JSON outputs resolve to the same file")
    if _same_file_or_path(input_path, csv_path) or _same_file_or_path(input_path, json_path):
        raise ValueError("summary output cannot be the input file, even with --force")
    for path, explicit in ((csv_path, explicit_csv), (json_path, explicit_json)):
        if not _lexists(path):
            continue
        if not force:
            raise FileExistsError("refusing to overwrite existing summary: %s" % path)
        if not explicit:
            raise ValueError("--force requires every existing output to be explicitly named")


def unavailable(data, source, reason):
    vbs = data.get("vbs") if isinstance(data.get("vbs"), dict) else {}
    return {
        "schema": SCHEMA,
        "status": "unavailable",
        "reason": reason,
        "source": str(source),
        "input_status": data.get("status"),
        "parameters": parameter_subset(data),
        "geometry": {key: data.get(key) for key in
                      ("Lx", "Ly", "shear", "pbc_geometry", "geometry_id")},
        "vbs": {key: vbs.get(key) for key in
                ("mode", "Nc", "Nq", "normalization", "spatial_boundary", "twist_applied")},
        "candidate_diagnostics": {},
        "groups": [],
        "plus_minus_classes": [],
        "rows": [],
    }


def _dict(value, path):
    if not isinstance(value, dict):
        raise SchemaError("%s must be an object" % path)
    return value


def _list(value, path):
    if not isinstance(value, list):
        raise SchemaError("%s must be an array" % path)
    return value


def _bool(value, path):
    if not isinstance(value, bool):
        raise SchemaError("%s must be boolean" % path)
    return value


def _string(value, path):
    if not isinstance(value, str):
        raise SchemaError("%s must be string" % path)
    return value


def _finite(value, path):
    result = number(value)
    if result is None:
        raise SchemaError("%s must be a finite number or null" % path)
    return result


def _nullable_number(value, path):
    if value is None:
        return None
    return _finite(value, path)


def _nonnegative_int(value, path):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise SchemaError("%s must be a non-negative integer" % path)
    return value


def _index_list(value, path, upper):
    values = _list(value, path)
    result = []
    for pos, item in enumerate(values):
        if isinstance(item, bool) or not isinstance(item, int) or not 0 <= item < upper:
            raise SchemaError("%s[%d] is not a valid point index" % (path, pos))
        result.append(item)
    if len(set(result)) != len(result):
        raise SchemaError("%s contains duplicate point indices" % path)
    return result


def _finite_array(value, path, length, *, allow_null):
    values = _list(value, path)
    if len(values) != length:
        raise SchemaError("%s length %d does not match %d" % (path, len(values), length))
    result = []
    for pos, item in enumerate(values):
        if item is None and allow_null:
            result.append(None)
        else:
            result.append(_finite(item, "%s[%d]" % (path, pos)))
    return result


def _count_array(value, path, length):
    values = _list(value, path)
    if len(values) != length:
        raise SchemaError("%s length %d does not match %d" % (path, len(values), length))
    return [_nonnegative_int(item, "%s[%d]" % (path, pos))
            for pos, item in enumerate(values)]


def _frequency_array(value, path, length, defined):
    values = _finite_array(value, path, length, allow_null=True)
    if defined:
        if any(item is None or item < 0 or item > 1 for item in values):
            raise SchemaError("%s must contain finite frequencies in [0,1]" % path)
    elif any(item is not None for item in values):
        raise SchemaError("%s must be null when stability_defined is false" % path)
    return values


def _normalize_estimate(value, path, *, total_replicas, invalid_replicas):
    value = _dict(value, path)
    required = ("defined", "reason", "real", "imag_diag", "stderr_real",
                "stderr_imag", "cov_real_imag", "jackknife_defined",
                "jackknife_reason")
    missing = [key for key in required if key not in value]
    if missing:
        raise SchemaError("%s missing fields: %s" % (path, ",".join(missing)))
    defined = _bool(value["defined"], path + ".defined")
    reason = _string(value["reason"], path + ".reason")
    jackknife_defined = _bool(value["jackknife_defined"], path + ".jackknife_defined")
    jackknife_reason = _string(value["jackknife_reason"], path + ".jackknife_reason")
    result = {
        "defined": defined,
        "reason": reason,
        "real": _nullable_number(value["real"], path + ".real"),
        "imag_diag": _nullable_number(value["imag_diag"], path + ".imag_diag"),
        "stderr_real": _nullable_number(value["stderr_real"], path + ".stderr_real"),
        "stderr_imag": _nullable_number(value["stderr_imag"], path + ".stderr_imag"),
        "cov_real_imag": _nullable_number(value["cov_real_imag"], path + ".cov_real_imag"),
        "jackknife_defined": jackknife_defined,
        "jackknife_reason": jackknife_reason,
    }
    if defined:
        if result["real"] is None or result["imag_diag"] is None:
            raise SchemaError("%s has defined=true but no finite center value" % path)
    else:
        if any(result[key] is not None for key in ("real", "imag_diag")):
            raise SchemaError("%s has defined=false but contains a center value" % path)
        if jackknife_defined:
            raise SchemaError("%s has jackknife_defined=true without a center" % path)
    if jackknife_defined:
        if not defined or total_replicas < 2 or invalid_replicas != 0:
            raise SchemaError("%s has inconsistent jackknife validity" % path)
        if any(result[key] is None for key in ("stderr_real", "stderr_imag", "cov_real_imag")):
            raise SchemaError("%s has defined jackknife errors missing" % path)
    elif any(result[key] is not None for key in ("stderr_real", "stderr_imag", "cov_real_imag")):
        raise SchemaError("%s has errors although jackknife_defined=false" % path)
    return result


def _normalize_qpoints(vbs, nq, Lx=None, Ly=None, shear=0):
    raw = vbs.get("qpoints")
    if not isinstance(raw, list):
        raise SchemaError("vbs.qpoints must be an array for candidate summaries")
    result = {}
    for pos, value in enumerate(raw):
        value = _dict(value, "vbs.qpoints[%d]" % pos)
        index = value.get("index")
        if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < nq:
            raise SchemaError("vbs.qpoints[%d].index is invalid" % pos)
        if index in result:
            raise SchemaError("duplicate qpoint index %d" % index)
        label = _string(value.get("label"), "vbs.qpoints[%d].label" % pos)
        n = value.get("n")
        if n is not None:
            n = _list(n, "vbs.qpoints[%d].n" % pos)
            if len(n) != 2 or any(isinstance(x, bool) or not isinstance(x, int) for x in n):
                raise SchemaError("vbs.qpoints[%d].n must have two integer entries" % pos)
            if Lx is not None and (not 0 <= n[0] < Lx or not 0 <= n[1] < Ly):
                raise SchemaError("vbs.qpoints[%d].n is outside the dual grid" % pos)
        on_grid = _bool(value.get("on_grid"), "vbs.qpoints[%d].on_grid" % pos)
        if on_grid != (n is not None):
            raise SchemaError("vbs.qpoints[%d] on_grid/n fields are inconsistent" % pos)
        fractions = _list(value.get("reciprocal_fraction"),
                          "vbs.qpoints[%d].reciprocal_fraction" % pos)
        if len(fractions) != 2:
            raise SchemaError("vbs.qpoints[%d].reciprocal_fraction must have two entries" % pos)
        fractions = [_finite(fractions[0], "vbs.qpoints[%d].reciprocal_fraction[0]" % pos),
                     _finite(fractions[1], "vbs.qpoints[%d].reciprocal_fraction[1]" % pos)]
        if on_grid and Lx is not None:
            expected = (n[0] / Lx,
                        (n[1] / Ly - shear * n[0] / (Lx * Ly)) % 1.)
            if any(abs(((fractions[i] - expected[i] + .5) % 1.) - .5) > 1e-9 for i in range(2)):
                raise SchemaError("vbs.qpoints[%d] does not use the declared sheared grid" % pos)
        result[index] = {"index": index, "label": label, "n": n,
                         "on_grid": on_grid,
                         "W_L_abs2": _finite(value.get("W_L_abs2"),
                                             "vbs.qpoints[%d].W_L_abs2" % pos)}
    if len(result) != nq or set(result) != set(range(nq)):
        raise SchemaError("vbs.qpoints indices do not cover Nq")
    return result


def _normalize_group(group, path, nq, class_count, warnings):
    group = _dict(group, path)
    for key in ("name", "estimator", "available", "participants", "participant_labels",
                "center_defined", "center_reason", "max_mz2_real", "winner_display_index", "winner_set",
                "winner_labels", "total_replicas", "valid_replicas", "invalid_replicas",
                "unique_win_counts", "tie_membership_counts", "class_tie_counts",
                "class_full_set_counts"):
        if key not in group:
            raise SchemaError("%s missing %s" % (path, key))
    name = _string(group["name"], path + ".name")
    if name not in GROUPS:
        raise SchemaError("unsupported candidate group %s" % name)
    estimator = _string(group["estimator"], path + ".estimator")
    if estimator != ESTIMATORS[name]:
        raise SchemaError("%s has wrong estimator" % path)
    participants = _index_list(group["participants"], path + ".participants", nq)
    participant_labels = _list(group["participant_labels"], path + ".participant_labels")
    if len(participant_labels) != len(participants) or not all(isinstance(x, str) for x in participant_labels):
        raise SchemaError("%s participant labels do not match participants" % path)
    available = _bool(group["available"], path + ".available")
    if available != bool(participants):
        raise SchemaError("%s available does not match participants" % path)
    center_defined = _bool(group["center_defined"], path + ".center_defined")
    center_reason = _string(group["center_reason"], path + ".center_reason")
    max_mz2 = _nullable_number(group["max_mz2_real"], path + ".max_mz2_real")
    winner_set = _index_list(group["winner_set"], path + ".winner_set", nq)
    winner_labels = _list(group["winner_labels"], path + ".winner_labels")
    if len(winner_labels) != len(winner_set) or not all(isinstance(x, str) for x in winner_labels):
        raise SchemaError("%s winner labels do not match winner_set" % path)
    if not set(winner_set).issubset(set(participants)):
        raise SchemaError("%s winner_set is not a participant subset" % path)
    if center_defined and (max_mz2 is None or not winner_set):
        raise SchemaError("%s center_defined requires max and winner_set" % path)
    if not center_defined and (max_mz2 is not None or winner_set or winner_labels):
        raise SchemaError("%s contains a center winner although center_defined=false" % path)
    display = group["winner_display_index"]
    if display is not None:
        if isinstance(display, bool) or not isinstance(display, int) or display not in winner_set:
            raise SchemaError("%s winner_display_index is invalid" % path)
    total = _nonnegative_int(group["total_replicas"], path + ".total_replicas")
    valid = _nonnegative_int(group["valid_replicas"], path + ".valid_replicas")
    invalid = _nonnegative_int(group["invalid_replicas"], path + ".invalid_replicas")
    if valid + invalid != total:
        raise SchemaError("%s replica counts do not add to total_replicas" % path)
    unique_counts = _count_array(group["unique_win_counts"], path + ".unique_win_counts", len(participants))
    tie_counts = _count_array(group["tie_membership_counts"], path + ".tie_membership_counts", len(participants))
    class_tie_counts = _count_array(group["class_tie_counts"], path + ".class_tie_counts", class_count)
    class_full_counts = _count_array(group["class_full_set_counts"], path + ".class_full_set_counts", class_count)
    for label, values in (("unique_win_counts", unique_counts), ("tie_membership_counts", tie_counts),
                          ("class_tie_counts", class_tie_counts), ("class_full_set_counts", class_full_counts)):
        if any(value > valid for value in values):
            raise SchemaError("%s.%s exceeds valid_replicas" % (path, label))

    if "jackknife_defined" not in group:
        jackknife_defined = False
        jackknife_reason = "missing_jackknife_status"
        warnings.append("%s missing jackknife validity" % path)
    else:
        jackknife_defined = _bool(group["jackknife_defined"], path + ".jackknife_defined")
        jackknife_reason = _string(group.get("jackknife_reason"), path + ".jackknife_reason")
    if "stability_defined" not in group:
        stability_defined = False
        stability_reason = "missing_stability_status"
        warnings.append("%s missing stability validity" % path)
    else:
        stability_defined = _bool(group["stability_defined"], path + ".stability_defined")
        stability_reason = _string(group.get("stability_reason"), path + ".stability_reason")
    conditional = _bool(group.get("conditional_on_valid_replicas", False),
                        path + ".conditional_on_valid_replicas")
    if jackknife_defined and (not center_defined or total < 2 or invalid != 0):
        raise SchemaError("%s has inconsistent jackknife_defined" % path)
    if stability_defined and (not center_defined or valid < 2):
        raise SchemaError("%s has inconsistent stability_defined" % path)
    if conditional != (stability_defined and invalid > 0):
        raise SchemaError("%s conditional_on_valid_replicas is inconsistent" % path)
    if not stability_defined:
        conditional = False

    def frequency(key, length):
        if key not in group:
            if stability_defined:
                raise SchemaError("%s missing %s" % (path, key))
            warnings.append("%s missing %s" % (path, key))
            return [None] * length
        return _frequency_array(group[key], path + "." + key, length, stability_defined)

    unique_frequency = frequency("unique_win_frequency", len(participants))
    tie_frequency = frequency("tie_membership_frequency", len(participants))
    class_tie_frequency = frequency("class_tie_frequency", class_count)
    class_full_frequency = frequency("class_full_set_frequency", class_count)
    return {
        "name": name, "estimator": estimator, "available": available,
        "participants": participants, "participant_labels": participant_labels,
        "center_defined": center_defined, "center_reason": center_reason,
        "max_mz2_real": max_mz2, "winner_display_index": display,
        "winner_set": winner_set, "winner_labels": winner_labels,
        "total_replicas": total, "valid_replicas": valid, "invalid_replicas": invalid,
        "jackknife_defined": jackknife_defined, "jackknife_reason": jackknife_reason,
        "stability_defined": stability_defined, "stability_reason": stability_reason,
        "conditional_on_valid_replicas": conditional,
        "unique_win_counts": unique_counts, "tie_membership_counts": tie_counts,
        "unique_win_frequency": unique_frequency, "tie_membership_frequency": tie_frequency,
        "class_tie_counts": class_tie_counts, "class_full_set_counts": class_full_counts,
        "class_tie_frequency": class_tie_frequency,
        "class_full_set_frequency": class_full_frequency,
    }


def _normalize_point(point, path, qpoints, groups, classes, total, invalid):
    point = _dict(point, path)
    for key in ("index", "label", "group", "u", "v", "on_grid", "W_L_abs2",
                "mz2_raw", "mz2_projected", "plus_minus_class", "class_self_conjugate"):
        if key not in point:
            raise SchemaError("%s missing %s" % (path, key))
    index = point["index"]
    if isinstance(index, bool) or not isinstance(index, int) or index not in qpoints:
        raise SchemaError("%s.index is not a qpoint index" % path)
    label = _string(point["label"], path + ".label")
    if label != qpoints[index]["label"]:
        raise SchemaError("%s.label does not match qpoint label" % path)
    group = _string(point["group"], path + ".group")
    if group not in ("grid_nonzero", "offgrid", "zero_excluded"):
        raise SchemaError("%s.group is invalid" % path)
    class_id = point["plus_minus_class"]
    if class_id is not None and (isinstance(class_id, bool) or not isinstance(class_id, int) or class_id < 0):
        raise SchemaError("%s.plus_minus_class is invalid" % path)
    if class_id is not None and class_id >= len(classes):
        raise SchemaError("%s.plus_minus_class is out of range" % path)
    if point["on_grid"] != qpoints[index]["on_grid"]:
        raise SchemaError("%s.on_grid does not match qpoint" % path)
    result = {
        "index": index, "label": label, "group": group,
        "u": _finite(point["u"], path + ".u"), "v": _finite(point["v"], path + ".v"),
        "on_grid": _bool(point["on_grid"], path + ".on_grid"),
        "W_L_abs2": _finite(point["W_L_abs2"], path + ".W_L_abs2"),
        "plus_minus_class": class_id,
        "class_self_conjugate": _bool(point["class_self_conjugate"], path + ".class_self_conjugate"),
        "mz2_raw": _normalize_estimate(point["mz2_raw"], path + ".mz2_raw",
                                        total_replicas=total, invalid_replicas=invalid),
        "mz2_projected": _normalize_estimate(point["mz2_projected"], path + ".mz2_projected",
                                               total_replicas=total, invalid_replicas=invalid),
    }
    expected = {"grid_nonzero": {"grid_raw"}, "offgrid": {"offgrid_raw", "offgrid_projected"},
                "zero_excluded": set()}[group]
    for name, value in groups.items():
        if (index in value["participants"]) != (name in expected):
            raise SchemaError("%s is not partitioned consistently with group %s" % (path, name))
    if class_id is not None and index not in classes[class_id]["members"]:
        raise SchemaError("%s.plus_minus_class does not contain the point" % path)
    return result


def _normalize_classes(values, nq):
    result = []
    for pos, value in enumerate(_list(values, "vbs.candidates.plus_minus_classes")):
        value = _dict(value, "vbs.candidates.plus_minus_classes[%d]" % pos)
        if value.get("id") != pos:
            raise SchemaError("plus/minus class ids are not contiguous")
        members = _index_list(value.get("members"),
                              "vbs.candidates.plus_minus_classes[%d].members" % pos, nq)
        member_labels = _list(value.get("member_labels"),
                              "vbs.candidates.plus_minus_classes[%d].member_labels" % pos)
        if len(member_labels) != len(members) or not all(isinstance(x, str) for x in member_labels):
            raise SchemaError("plus/minus class member labels do not match members")
        canonical = _list(value.get("canonical_fraction"),
                          "vbs.candidates.plus_minus_classes[%d].canonical_fraction" % pos)
        if len(canonical) != 2:
            raise SchemaError("plus/minus class canonical fraction must have two entries")
        result.append({"id": pos, "members": members, "member_labels": member_labels,
                       "self_conjugate": _bool(value.get("self_conjugate"),
                                                "plus/minus class self_conjugate"),
                       "canonical_fraction": [_finite(canonical[0], "canonical_fraction[0]"),
                                              _finite(canonical[1], "canonical_fraction[1]")]})
    return result


def _point_row(point, qpoint, groups):
    row = {
        "status": "unavailable" if point["group"] == "zero_excluded" else "partial",
        "reason": "excluded_Q0" if point["group"] == "zero_excluded" else None,
        "point_index": point["index"], "label": point["label"],
        "u": point["u"], "v": point["v"], "on_grid": point["on_grid"],
        "n1": qpoint.get("n", [None, None])[0] if qpoint.get("n") is not None else None,
        "n2": qpoint.get("n", [None, None])[1] if qpoint.get("n") is not None else None,
        "W_L_abs2": point["W_L_abs2"], "group": point["group"],
        "mz2_raw_real": point["mz2_raw"]["real"], "mz2_raw_imag_diag": point["mz2_raw"]["imag_diag"],
        "mz2_raw_stderr_real": point["mz2_raw"]["stderr_real"], "mz2_raw_stderr_imag": point["mz2_raw"]["stderr_imag"],
        "mz2_raw_cov_real_imag": point["mz2_raw"]["cov_real_imag"],
        "raw_center_defined": point["mz2_raw"]["defined"], "raw_center_reason": point["mz2_raw"]["reason"],
        "raw_jackknife_defined": point["mz2_raw"]["jackknife_defined"], "raw_jackknife_reason": point["mz2_raw"]["jackknife_reason"],
        "mz2_projected_real": point["mz2_projected"]["real"], "mz2_projected_imag_diag": point["mz2_projected"]["imag_diag"],
        "mz2_projected_stderr_real": point["mz2_projected"]["stderr_real"], "mz2_projected_stderr_imag": point["mz2_projected"]["stderr_imag"],
        "mz2_projected_cov_real_imag": point["mz2_projected"]["cov_real_imag"],
        "projected_center_defined": point["mz2_projected"]["defined"], "projected_center_reason": point["mz2_projected"]["reason"],
        "projected_jackknife_defined": point["mz2_projected"]["jackknife_defined"], "projected_jackknife_reason": point["mz2_projected"]["jackknife_reason"],
    }
    participant_groups = []
    for name in GROUPS:
        group = groups[name]
        if point["index"] in group["participants"]:
            position = group["participants"].index(point["index"])
            participant_groups.append(group)
        else:
            position = None
        prefix = name
        for field in ("center_defined", "center_reason", "jackknife_defined", "jackknife_reason",
                      "stability_defined", "stability_reason", "conditional_on_valid_replicas",
                      "total_replicas", "valid_replicas", "invalid_replicas", "max_mz2_real"):
            row[prefix + "_" + field] = group[field] if position is not None else None
        row[prefix + "_winner_labels"] = labels(group["winner_labels"]) if position is not None else None
        row[prefix + "_unique_win_frequency"] = (group["unique_win_frequency"][position]
                                                  if position is not None else None)
        row[prefix + "_tie_membership_frequency"] = (group["tie_membership_frequency"][position]
                                                      if position is not None else None)
    if point["group"] != "zero_excluded":
        center_ok = point["mz2_raw"]["defined"] if point["group"] == "grid_nonzero" else (
            point["mz2_raw"]["defined"] and point["mz2_projected"]["defined"])
        stability_ok = all(group["stability_defined"] for group in participant_groups)
        jackknife_ok = all(group["jackknife_defined"] for group in participant_groups)
        row["status"] = "available" if center_ok and stability_ok and jackknife_ok else (
            "partial" if center_ok else "unavailable")
        if row["status"] == "unavailable":
            row["reason"] = "center_undefined"
        elif row["status"] == "partial":
            row["reason"] = "jackknife_or_stability_unavailable"
    return row


CSV_FIELDS = [
    "status", "reason", "point_index", "label", "u", "v", "on_grid", "n1", "n2", "W_L_abs2", "group",
    "raw_center_defined", "raw_center_reason", "mz2_raw_real", "mz2_raw_imag_diag",
    "raw_jackknife_defined", "raw_jackknife_reason", "mz2_raw_stderr_real", "mz2_raw_stderr_imag",
    "mz2_raw_cov_real_imag", "projected_center_defined", "projected_center_reason",
    "mz2_projected_real", "mz2_projected_imag_diag", "projected_jackknife_defined",
    "projected_jackknife_reason", "mz2_projected_stderr_real", "mz2_projected_stderr_imag",
    "mz2_projected_cov_real_imag",
]
for _group in GROUPS:
    CSV_FIELDS.extend([
        "%s_center_defined" % _group, "%s_center_reason" % _group,
        "%s_jackknife_defined" % _group, "%s_jackknife_reason" % _group,
        "%s_stability_defined" % _group, "%s_stability_reason" % _group,
        "%s_conditional_on_valid_replicas" % _group,
        "%s_total_replicas" % _group, "%s_valid_replicas" % _group,
        "%s_invalid_replicas" % _group, "%s_max_mz2_real" % _group,
        "%s_winner_labels" % _group, "%s_unique_win_frequency" % _group,
        "%s_tie_membership_frequency" % _group,
    ])


def flatten_csv_row(row):
    return {key: row.get(key) for key in CSV_FIELDS}


def summarize(data, source):
    if not isinstance(data, dict):
        raise SchemaError("top-level results JSON must be an object")
    if data.get("status") != "complete":
        raise SchemaError("refusing incomplete or aborted MC results")
    geometry = normalize_geometry(data)
    vbs = data.get("vbs")
    if vbs is None:
        return unavailable(data, source, "missing vbs object")
    vbs = _dict(vbs, "vbs")
    candidates = vbs.get("candidates")
    if candidates is None:
        return unavailable(data, source, "missing vbs.candidates; stability is unavailable")
    candidates = _dict(candidates, "vbs.candidates")
    schema = candidates.get("schema")
    if schema == "ivbs-candidates-v1":
        return unavailable(data, source, "legacy ivbs-candidates-v1 lacks validity semantics")
    if schema != CANDIDATE_SCHEMA:
        return unavailable(data, source, "unsupported candidates schema")
    if candidates.get("level") != "finest_blocks":
        return unavailable(data, source, "candidate stability is not recorded at finest_blocks")
    if "groups" not in candidates or "points" not in candidates:
        return unavailable(data, source, "candidate schema has no groups or points")
    if "Nq" not in vbs or isinstance(vbs["Nq"], bool) or not isinstance(vbs["Nq"], int) or vbs["Nq"] < 0:
        raise SchemaError("vbs.Nq must be a non-negative integer")
    nq = vbs["Nq"]
    qpoints = _normalize_qpoints(vbs, nq, data.get("Lx"), data.get("Ly"), data.get("shear", 0))
    class_values = _normalize_classes(candidates.get("plus_minus_classes"), nq)
    warnings = []
    required_root = ("ratio_defined", "center_defined", "center_reason", "jackknife_defined",
                     "jackknife_reason", "stability_defined", "stability_reason",
                     "total_replicas", "valid_replicas", "invalid_replicas")
    missing_root = [key for key in required_root if key not in candidates]
    if missing_root:
        return unavailable(data, source, "candidate schema missing validity fields: " + ",".join(missing_root))
    root_ratio = _bool(candidates["ratio_defined"], "vbs.candidates.ratio_defined")
    root_center = _bool(candidates["center_defined"], "vbs.candidates.center_defined")
    root_center_reason = _string(candidates["center_reason"], "vbs.candidates.center_reason")
    root_jack = _bool(candidates["jackknife_defined"], "vbs.candidates.jackknife_defined")
    root_jack_reason = _string(candidates["jackknife_reason"], "vbs.candidates.jackknife_reason")
    root_stability = _bool(candidates["stability_defined"], "vbs.candidates.stability_defined")
    root_stability_reason = _string(candidates["stability_reason"], "vbs.candidates.stability_reason")
    total = _nonnegative_int(candidates["total_replicas"], "vbs.candidates.total_replicas")
    valid = _nonnegative_int(candidates["valid_replicas"], "vbs.candidates.valid_replicas")
    invalid = _nonnegative_int(candidates["invalid_replicas"], "vbs.candidates.invalid_replicas")
    if valid + invalid != total:
        raise SchemaError("candidate root replica counts do not add to total_replicas")
    invalid_reasons = candidates.get("invalid_replica_reasons", [])
    invalid_reasons = _list(invalid_reasons, "vbs.candidates.invalid_replica_reasons")
    if len(invalid_reasons) != invalid or not all(isinstance(item, str) for item in invalid_reasons):
        raise SchemaError("candidate invalid_replica_reasons does not match invalid_replicas")
    groups_raw = _list(candidates["groups"], "vbs.candidates.groups")
    groups = {}
    for pos, group in enumerate(groups_raw):
        normalized = _normalize_group(group, "vbs.candidates.groups[%d]" % pos, nq,
                                      len(class_values), warnings)
        if normalized["name"] in groups:
            raise SchemaError("duplicate candidate group %s" % normalized["name"])
        groups[normalized["name"]] = normalized
    missing_groups = [name for name in GROUPS if name not in groups]
    if missing_groups:
        return unavailable(data, source, "candidate groups missing: " + ",".join(missing_groups))
    for name, group in groups.items():
        if (group["total_replicas"], group["valid_replicas"], group["invalid_replicas"]) != (total, valid, invalid):
            raise SchemaError("group %s replica counts disagree with candidate root" % name)
        for index, label in zip(group["participants"], group["participant_labels"]):
            if label != qpoints[index]["label"]:
                raise SchemaError("group %s participant label does not match qpoint" % name)
        for index, label in zip(group["winner_set"], group["winner_labels"]):
            if label != qpoints[index]["label"]:
                raise SchemaError("group %s winner label does not match qpoint" % name)
    points_raw = _list(candidates["points"], "vbs.candidates.points")
    if not points_raw:
        return unavailable(data, source, "candidate points are empty")
    points = {}
    for pos, point in enumerate(points_raw):
        normalized = _normalize_point(point, "vbs.candidates.points[%d]" % pos,
                                      qpoints, groups, class_values, total, invalid)
        if normalized["index"] in points:
            raise SchemaError("duplicate candidate point index %d" % normalized["index"])
        points[normalized["index"]] = normalized
    if len(points) != nq or set(points) != set(range(nq)):
        raise SchemaError("candidate points do not cover Nq")
    if root_jack and (not root_center or total < 2 or invalid != 0):
        raise SchemaError("candidate root jackknife validity is inconsistent")
    if root_stability and (not root_center or valid < 2):
        raise SchemaError("candidate root stability validity is inconsistent")
    available_groups = [group for group in groups.values() if group["available"]]
    expected_root_center = bool(available_groups) and all(group["center_defined"] for group in available_groups)
    if root_center != expected_root_center:
        raise SchemaError("candidate root center_defined is inconsistent with groups")
    if not available_groups:
        status, reason = "unavailable", "no candidate participants"
    elif not any(group["center_defined"] for group in available_groups):
        status, reason = "unavailable", "no candidate center estimate is defined"
    elif (warnings or not all(group["center_defined"] and group["jackknife_defined"] and
                              group["stability_defined"] for group in available_groups)):
        status, reason = "partial", "center exists but jackknife or stability is unavailable"
    else:
        status, reason = "available", None
    rows = [_point_row(points[index], qpoints[index], groups) for index in sorted(points)]
    return {
        "schema": SCHEMA, "status": status, "reason": reason,
        "source": str(source), "input_status": data.get("status"),
        "parameters": parameter_subset(data),
        "geometry": geometry,
        "vbs": {key: vbs.get(key) for key in
                ("mode", "Nc", "Nq", "normalization", "spatial_boundary", "twist_applied")},
        "candidate_diagnostics": {
            "schema": schema, "level": candidates["level"], "ratio_defined": root_ratio,
            "center_defined": root_center, "center_reason": root_center_reason,
            "jackknife_defined": root_jack, "jackknife_reason": root_jack_reason,
            "stability_defined": root_stability, "stability_reason": root_stability_reason,
            "total_replicas": total, "valid_replicas": valid, "invalid_replicas": invalid,
            "invalid_replica_reasons": invalid_reasons,
            "validation_warnings": warnings, "excluded_Q0_rule": candidates.get("excluded_Q0_rule"),
            "scoring": candidates.get("scoring"), "tie_tolerance": candidates.get("tie_tolerance"),
            "winner_set_rule": candidates.get("winner_set_rule"),
            "stability_definition": candidates.get("stability_definition"),
            "groups": [groups[name] for name in GROUPS], "plus_minus_classes": class_values,
        },
        "groups": [groups[name] for name in GROUPS], "plus_minus_classes": class_values,
        "rows": rows,
    }


def _csv_text(summary):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    rows = summary["rows"]
    if not rows:
        marker = {key: None for key in CSV_FIELDS}
        marker["status"] = summary["status"]
        marker["reason"] = summary.get("reason")
        marker["group"] = "unavailable"
        rows = [marker]
    for row in rows:
        writer.writerow(flatten_csv_row(row))
    return stream.getvalue()


def _temporary_write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".%s." % path.name, suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def write_outputs(summary, csv_path, json_path, *, input_path=None, force=False,
                  explicit_csv=True, explicit_json=True):
    """Validate, construct both documents, then replace each output once."""
    if input_path is not None:
        validate_output_paths(input_path, csv_path, json_path, force=force,
                              explicit_csv=explicit_csv, explicit_json=explicit_json)
    elif _same_file_or_path(csv_path, json_path):
        raise ValueError("CSV and JSON outputs resolve to the same file")
    json_text = json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    csv_text = _csv_text(summary)
    if input_path is None:
        for path, explicit in ((csv_path, explicit_csv), (json_path, explicit_json)):
            if _lexists(path) and (not force or not explicit):
                raise FileExistsError("refusing to overwrite existing summary: %s" % path)
    _temporary_write(json_path, json_text)
    _temporary_write(csv_path, csv_text)


def _check_finite_tree(value, path):
    if isinstance(value, float) and not math.isfinite(value):
        raise SchemaError("%s contains a non-finite number" % path)
    if isinstance(value, dict):
        for key, child in value.items():
            _check_finite_tree(child, "%s.%s" % (path, key))
    elif isinstance(value, list):
        for pos, child in enumerate(value):
            _check_finite_tree(child, "%s[%d]" % (path, pos))


def _load_json(path):
    def reject_constant(value):
        raise ValueError("non-finite JSON constant %s is not allowed" % value)
    data = json.loads(path.read_text(encoding="utf-8"), parse_constant=reject_constant)
    _check_finite_tree(data, "results")
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="one coupled_mc *_results.json file")
    parser.add_argument("--csv", type=Path, help="CSV output path")
    parser.add_argument("--json", type=Path, help="JSON output path")
    parser.add_argument("--force", action="store_true",
                        help="overwrite existing explicitly named summary outputs")
    args = parser.parse_args()
    try:
        if not args.results.is_file():
            raise OSError("input results file does not exist or is not a regular file: %s" % args.results)
        data = _load_json(args.results)
        summary = summarize(data, args.results)
        csv_path, json_path = output_paths(args.results, args.csv, args.json)
        validate_output_paths(args.results, csv_path, json_path, force=args.force,
                              explicit_csv=args.csv is not None, explicit_json=args.json is not None)
        write_outputs(summary, csv_path, json_path, input_path=args.results, force=args.force,
                      explicit_csv=args.csv is not None, explicit_json=args.json is not None)
        print(json.dumps({"status": summary["status"], "csv": str(csv_path),
                          "json": str(json_path), "rows": len(summary["rows"])}, sort_keys=True))
        return 0
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
