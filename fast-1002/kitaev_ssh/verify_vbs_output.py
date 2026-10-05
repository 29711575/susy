#!/usr/bin/env python3
"""Independent NumPy inverse-FFT / eigensolve / joint-JK check of actual output."""
import argparse
import csv
from pathlib import Path
import numpy as np
from plot_vbs import load_run, value


def inverse_structure_factor(s, nx, ny, shear):
    """Independent inverse of the sheared character table, C(d)=sum_q e^-iqd S(q)."""
    if shear == 0:
        return np.fft.fft2(s, axes=(0, 1)) / (nx * ny)
    staged = np.fft.fft(s, axis=1) / ny
    n1 = np.arange(nx)[:, None]
    y = np.arange(ny)[None, :]
    staged *= np.exp(2j * np.pi * shear * n1 * y / (nx * ny))[:, :, None, None]
    return np.fft.fft(staged, axis=0) / nx


def verify(path):
    data, ratio_error = load_run(path)
    raw_names = list(data["observables"])
    nc, nx, ny = data["vbs"]["Nc"], data["Lx"], data["Ly"]
    block_file = path.with_name(path.name.removesuffix("_results.json")+"_blocks.csv")
    with block_file.open() as stream:
        blocks = list(csv.DictReader(stream))
    phase = np.array([complex(float(b["sum_phase_real"]), float(b["sum_phase_imag"])) for b in blocks])
    nums = np.array([[complex(float(b[f"sum_phase_{n}_real"]), float(b[f"sum_phase_{n}_imag"])) for n in raw_names] for b in blocks])

    def derived(ratios):
        r = dict(zip(raw_names, ratios))
        s = np.array([[[r[f"S_q{q}_{mu}{nu}"] for nu in "xyz"] for mu in "xyz"] for q in range(nc)])
        output = {}
        for q in range(nc):
            output[f"TrS_q{q}"] = np.trace(s[q])
            output[f"mz2_q{q}"] = s[q, 2, 2]/nc
            output[f"leading_eigenvalue_q{q}"] = np.linalg.eigvalsh((s[q]+s[q].conj().T)/2)[-1]
            output[f"hermitian_residual_q{q}"] = np.linalg.norm(s[q]-s[q].conj().T)
        c = inverse_structure_factor(s.reshape(nx, ny, 3, 3), nx, ny,
                                     data["shear"]).reshape(nc, 3, 3)
        uniform = np.array([r[f"Delta_q0_{mu}"] for mu in "xyz"])
        for prefix, matrix in (("C", c), ("C_connected", c-uniform[:, None]*uniform[None, :])):
            for d in range(nc):
                for mu, axis in enumerate("xyz"):
                    for nu, other in enumerate("xyz"):
                        output[f"{prefix}_q{d}_{axis}{other}"] = matrix[d, mu, nu]
        return output

    max_value, max_error, max_cov = 0., 0., 0.
    for level in data["vbs"]["blocking"]:
        factor = level["block_samples"] // data["block_samples"]
        p = phase.reshape(-1, factor).sum(axis=1)
        n = nums.reshape(-1, factor, len(raw_names)).sum(axis=1)
        total_p, total_n = p.sum(), n.sum(axis=0)
        target = derived(total_n/total_p)
        for name, z in target.items():
            max_value = max(max_value, abs(z-value(level["derived"][name])))
        if level["jackknife_defined"]:
            replicas = [derived((total_n-nb)/(total_p-pb)) for pb, nb in zip(p, n)]
            for name in target:
                v = np.array([r[name] for r in replicas])
                v -= v.mean()
                scale = (len(v)-1)/len(v)
                re, im = np.sqrt(scale*np.sum(v.real**2)), np.sqrt(scale*np.sum(v.imag**2))
                cov = scale*np.sum(v.real*v.imag)
                actual = level["derived"][name]
                max_error = max(max_error, abs(re-actual["stderr_real"]), abs(im-actual["stderr_imag"]))
                max_cov = max(max_cov, abs(cov-actual["cov_real_imag"]))
    tolerance = 5e-11  # accumulation + two independent FFT/eigen implementations
    if max(max_value, max_error, max_cov) > tolerance:
        raise AssertionError("independent derived statistics check failed")
    print(f"{path.name}: block_ratio={ratio_error:.3g}, inverseFFT/eigen/value={max_value:.3g}, "
          f"joint_JK_error={max_error:.3g}, covariance={max_cov:.3g}; tolerance={tolerance:g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, nargs="+")
    for result in parser.parse_args().results:
        verify(result)
