"""Freeze complete CSV records, recover identity, and validate chain boundaries."""
from dataclasses import dataclass, replace
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np


def frozen_bytes(path):
    path = Path(path)
    with path.open("rb") as stream:
        size = stream.seek(0, 2)
        stream.seek(0)
        data = stream.read(size)
    used = data.rfind(b"\n") + 1
    return data[:used], dict(path=str(path.resolve()), bytes_cap=size,
                            bytes_used=used, incomplete_tail_bytes=size-used,
                            sha256=hashlib.sha256(data[:used]).hexdigest())


def options(argv):
    out = {}
    i = 0
    while i < len(argv):
        word = str(argv[i])
        if word.startswith("--"):
            if i+1 >= len(argv):
                raise ValueError("missing argument: " + word)
            out[word[2:]] = argv[i+1]
            i += 2
        else:
            i += 1
    return out


def identity(meta, argv=None, point=None):
    """Only recover explicitly recorded settings; never apply today's defaults."""
    result = dict(point or {})
    a = options(argv or meta.get("argv") or (point or {}).get("argv", []))
    a.update(meta.get("resume", {}).get("effective_options", {}))
    aliases = {"Lx":"Lx", "Ly":"Ly", "omega-ph":"omega_ph", "dt":"dt",
               "intervals":"intervals", "mass":"mass", "g":"g", "Jx":"Jx",
               "Jy":"Jy", "Jz":"Jz", "seed":"seed", "shear":"shear",
               "measure-every":"measure_every", "block-samples":"block_samples",
               "sweeps":"measurement_sweeps", "matter-trial":"matter_trial_policy"}
    integers = {"Lx", "Ly", "intervals", "seed", "shear", "measure_every", "block_samples", "measurement_sweeps"}
    for key, target in aliases.items():
        if key in a:
            result[target] = str(a[key]) if target == "matter_trial_policy" else (int(a[key]) if target in integers else float(a[key]))
    for key in ["trial-dx", "trial-dy", "trial-dz"]:
        if key in a and all(k in a for k in ["trial-dx", "trial-dy", "trial-dz"]):
            result["trial_shift"] = [float(a[k]) for k in ["trial-dx", "trial-dy", "trial-dz"]]
    for key in ["Lx", "Ly", "seed", "g", "mass", "omega_ph", "dt", "intervals", "measure_every"]:
        if key in meta and key in result and not np.isclose(float(meta[key]), float(result[key]), rtol=1e-12, atol=1e-12):
            raise ValueError("argv/metadata identity mismatch: " + key)
    result.update(meta)
    if "L" in result:
        result.setdefault("Lx", result["L"])
        result.setdefault("Ly", result["L"])
    if "omega" in result:
        result.setdefault("omega_ph", result["omega"])
    if "J" in result:
        for key, value in zip(["Jx", "Jy", "Jz"], result["J"]):
            result.setdefault(key, value)
    if "lambda_value" in result:
        result.setdefault("lambda", result["lambda_value"])
    if "sweeps" in result:
        result.setdefault("measurement_sweeps", result["sweeps"])
    if "dt" in result and "intervals" in result:
        tau = result["dt"]*result["intervals"]
        if "tau" in result and not np.isclose(result["tau"], tau, rtol=0, atol=1e-10):
            raise ValueError("inconsistent tau/intervals/dt")
        result["tau"] = tau
    if "mass" in result and "omega_ph" in result:
        result.setdefault("k", result["mass"]*result["omega_ph"]**2)
    if "g" in result and result.get("k", 0) > 0:
        result.setdefault("lambda", result["g"]**2/result["k"])
    result["recorded_options"] = a
    if "matter_trial_sector" in result:
        result.setdefault("matter_trial_policy", result["matter_trial_sector"].get("policy"))
    return result


@dataclass
class Chain:
    path: Path
    names: list
    counts: np.ndarray
    sums: np.ndarray  # common complex phase followed by complex phase*O
    classifications: np.ndarray
    metadata: dict
    provenance: dict
    start_sweep: int = 0

    @property
    def stride(self):
        if "measure_every" not in self.metadata:
            raise ValueError(f"{self.path}: measure_every missing; provide --metadata")
        value = int(self.metadata["measure_every"])
        if value <= 0:
            raise ValueError("measure_every must be positive")
        return value

    @property
    def ends(self):
        return self.start_sweep + np.cumsum(self.counts)*self.stride

    def window(self, cut=0, end=None):
        starts = self.ends-self.counts*self.stride
        end = int(self.ends[-1]) if end is None else end
        if cut not in set(starts.tolist()+self.ends.tolist()) or end not in set(starts.tolist()+self.ends.tolist()):
            raise ValueError("window must align with complete stored blocks; no partial block can be reconstructed")
        keep = (starts >= cut) & (self.ends <= end)
        if not keep.any():
            raise ValueError(f"{self.path}: no complete blocks in requested window")
        return replace(self, counts=self.counts[keep], sums=self.sums[keep],
                       classifications=self.classifications[keep], start_sweep=cut)

    def reblock(self, factor):
        if factor < 1:
            raise ValueError("reblock factor must be positive")
        starts = np.arange(0, len(self.counts), factor)
        # Retain short tail groups within their own chain. The statistics
        # module uses their actual sample fractions in delete-group errors.
        return replace(self, counts=np.add.reduceat(self.counts, starts),
                       sums=np.add.reduceat(self.sums, starts, axis=0),
                       classifications=np.add.reduceat(self.classifications, starts, axis=0))


def read_chain(path, metadata=None, allow_prefix=False):
    path = Path(path).expanduser().resolve()
    if path.is_dir():
        path = path/"mc_blocks.csv"
    elif path.name.endswith("_results.json"):
        path = path.with_name(path.name[:-13]+"_blocks.csv")
    elif not path.name.endswith("_blocks.csv"):
        path = Path(str(path)+"_blocks.csv")
    prefix = str(path)[:-11]
    meta, evidence = {}, []
    for candidate in [path.parent/"point.json", path.parent/"source.json", Path(prefix+"_matter_trial.json"), Path(prefix+".argv.json"), Path(prefix+"_results.json")]:
        if candidate.exists():
            raw = candidate.read_bytes()
            item = json.loads(raw)
            evidence.append(dict(path=str(candidate), sha256=hashlib.sha256(raw).hexdigest()))
            if candidate.name.endswith(".argv.json"):
                meta = identity(meta, argv=item)
            elif candidate.name == "point.json":
                meta = identity(meta, point=item)
            elif candidate.name.endswith("_matter_trial.json"):
                meta["matter_trial_sector"] = item
            else:
                meta = identity({**meta, **item})
    if metadata:
        item = json.loads(Path(metadata).read_text()) if isinstance(metadata, (str, Path)) else metadata
        for key, value in item.items():
            if key in meta and key in {"seed", "Lx", "Ly", "dt", "intervals", "mass", "g", "omega_ph", "measure_every"} and meta[key] != value:
                raise ValueError("metadata conflicts with recorded identity: " + key)
        meta = identity({**meta, **item})
    raw, provenance = frozen_bytes(path)
    if provenance["incomplete_tail_bytes"] and not allow_prefix:
        raise ValueError("incomplete final CSV record; use --allow-prefix for an explicitly partial snapshot")
    if meta.get("status") != "complete" and not allow_prefix:
        raise ValueError("sampling completion is not established; use --allow-prefix to recover available complete blocks")
    reader = csv.reader(io.StringIO(raw.decode("utf-8")))
    try:
        header = next(reader)
    except StopIteration:
        raise ValueError("empty raw-block file") from None
    fixed = ["block", "count", "sum_phase_real", "sum_phase_imag", "positive_count", "negative_count", "complex_count"]
    if header[:7] != fixed or len(set(header)) != len(header) or (len(header)-7) % 2:
        raise ValueError("invalid or duplicated block CSV columns")
    names = []
    for i in range(7, len(header), 2):
        if not header[i].startswith("sum_phase_") or not header[i].endswith("_real") or header[i+1] != header[i][:-5]+"_imag":
            raise ValueError("unpaired complex observable columns")
        names.append(header[i][10:-5])
    counts, sums, classes = [], [], []
    for i, row in enumerate(reader):
        if len(row) != len(header):
            raise ValueError(f"malformed complete CSV record at block {i}")
        if int(row[0]) != i or int(row[1]) <= 0:
            raise ValueError("missing, duplicated, or unordered block")
        count = int(row[1]); cat = [int(x) for x in row[4:7]]
        if min(cat) < 0 or sum(cat) != count:
            raise ValueError("invalid phase classifications")
        values = [complex(float(row[2]), float(row[3]))] + [complex(float(row[j]), float(row[j+1])) for j in range(7, len(row), 2)]
        if not np.isfinite(values).all() or abs(values[0]) > count*(1+1e-10):
            raise ValueError("nonfinite sums or impossible phase sum")
        counts.append(count); sums.append(values); classes.append(cat)
    if not counts or len(set(counts)) != 1:
        raise ValueError("require nonempty equal-sized stored blocks")
    if meta.get("status") == "complete" and "samples" in meta and sum(counts) != meta["samples"]:
        raise ValueError("complete result/raw-block sample count mismatch")
    provenance.update(metadata_sources=evidence, original_status=meta.get("status", "unknown"),
                      complete_result=meta.get("status") == "complete")
    return Chain(path, names, np.array(counts, dtype=np.int64), np.array(sums),
                 np.array(classes, dtype=np.int64), meta, provenance)


def checked_group(chains):
    if not chains or len({str(c.path) for c in chains}) != len(chains):
        raise ValueError("empty or duplicated input chains")
    seeds = [c.metadata.get("seed") for c in chains]
    if len(chains) > 1 and (None in seeds or len(set(seeds)) != len(seeds)):
        raise ValueError("independent chains require recorded distinct seeds; ancestor and resumed output must not be pooled")
    # An explicit ancestry match is invalid even if a caller edited a seed.
    prefixes = {str(c.path)[:-11] for c in chains}
    if len(chains) > 1:
        required = ["Lx", "Ly", "shear", "Jx", "Jy", "Jz", "g", "mass", "omega_ph", "dt", "intervals", "trial_shift", "matter_trial_policy"]
        for c in chains:
            missing = [k for k in required if k not in c.metadata]
            if missing:
                raise ValueError("independent-chain identity incomplete; provide metadata: " + ",".join(missing))
    for c in chains:
        ancestor = c.metadata.get("resume", {}).get("ancestor_output_prefix")
        if ancestor and str(Path(ancestor).resolve()) in prefixes:
            raise ValueError("resumed blocks already include their ancestor")
    keys = ("Lx", "Ly", "shear", "Jx", "Jy", "Jz", "g", "mass", "omega_ph", "dt", "intervals", "trial_shift", "matter_trial", "matter_trial_policy", "matter_trial_sector", "scope", "wilson_loops", "gauge_u", "time_contour", "vbs_field", "total_energy", "source_archive_sha256", "binary_sha256")
    for c in chains[1:]:
        if c.names != chains[0].names or c.stride != chains[0].stride:
            raise ValueError("observable schema or measurement stride differs")
        for key in keys:
            if c.metadata.get(key) != chains[0].metadata.get(key):
                raise ValueError("different or missing chain identity: " + key)
    return chains
