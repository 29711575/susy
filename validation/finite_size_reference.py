"""Deterministic bare-spectrum reference; no QMC sampling or optimization.

Run from the repository root with Python 3 and NumPy. Prints a CSV to stdout.
Conventions: docs/PHYSICS_CONVENTIONS.md, square torus, shear=0, Wilson++.
"""
import csv
import sys

import numpy as np


def reference(L):
    N = L * L
    T = np.zeros((N, N))
    for x in range(L):
        for y in range(L):
            row = L * x + y
            T[row, row] += 0.5
            T[row, L * ((x - 1) % L) + y] += 1.0
            T[row, L * x + (y - 1) % L] += 1.0
    singular = np.linalg.svd(T, compute_uv=False)
    k = 2 * np.pi * np.arange(L) / L
    bloch = np.abs(0.5 + np.exp(-1j * k[:, None]) + np.exp(-1j * k[None, :]))
    np.testing.assert_allclose(np.sort(singular), np.sort(bloch.ravel()), atol=3e-14, rtol=0)
    sign, _ = np.linalg.slogdet(T)
    assert singular.min() > 1e-10, 'Zero modes require a separate parity analysis'
    # For the recorded A/B Fock convention, vacuum parity=(-1)^N sign(det T).
    vacuum = (-1) ** N * int(sign)
    # PHY-01 at Lx=Ly=L, shear=0, and all code gauge links +1.
    required = (-1) ** (2 * L + N)
    needs_excitation = vacuum != required
    cost = 2 * singular.min() if needs_excitation else 0.0
    return dict(L=L,unit_cells=N,matter_sites=2*N,detT_sign=int(sign),
                vacuum_code_parity=vacuum,required_code_parity=required,
                minimum_T_singular=float(singular.min()),
                minimum_single_quasiparticle_energy=float(2*singular.min()),
                lowest_degeneracy=int(np.isclose(singular, singular.min(), atol=1e-12, rtol=0).sum()),
                needs_excitation=needs_excitation,parity_cost_total=float(cost),
                parity_cost_per_cell=float(cost/N))


if __name__ == '__main__':
    rows = [reference(L) for L in (6, 7, 8, 9)]
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
