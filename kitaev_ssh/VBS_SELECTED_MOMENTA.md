# VBS selected-momentum measurement (`--vbs-mode`)

This document covers the all/selected/both measurement modes on the current
integer-shear torus. With `shear=0` the existing rectangular path and
momentum ordering are retained. With `shear!=0`, the spatial quotient changes:
the seam bonds, fixed-gauge matter Hamiltonian, free trial state and QMC
weights generally change. This is an integer geometric boundary change, not a
measurement-only change and not a continuous U(1) twist. See
`SHEARED_PBC.md` for the quotient and `TWIST_FEASIBILITY_REVIEW.md` for the
separate continuous-twist review.

The JSON continues to report `"spatial_boundary":"pbc"` and
`"twist_applied":false`. The latter means that no continuous U(1) phase was
added; it does not mean that a nonzero integer shear has the same seam as
axis-aligned PBC.

## 1. The three modes (`--vbs 1` master switch)

| mode | momenta evaluated | estimator | real-space C |
|---|---|---|---|
| `all` (default) | full PBC grid, `Nc = Lx*Ly` momenta, index `q=n1*Ly+n2` | MKL complex FFT per bond direction; staged x/phase/y FFT for shear | full inverse-FFT `C`, `C_connected` |
| `selected` | only the `Nq` momenta from `--vbs-q-file` | direct finite-window sums, O(Nq*Nc) per configuration, precomputed cell phases | not produced (`real_space_correlation=unavailable_selected_mode`) |
| `both` | full grid block first, then the selected block | both of the above | as in `all` |

- Without `--vbs 1` the binary is bit-identical in behaviour to the previous
  version. `--vbs-mode`/`--vbs-q-file` without `--vbs 1` is an error, `all`
  mode with `--vbs-q-file` is an error, `selected`/`both` without a q-file is
  an error (no silent fallbacks).
- Every selected momentum is a **finite-window Fourier estimator** at an
  arbitrary `Q`. Off-grid points do not become allowed lattice momenta, do not
  remove the incommensurability mismatch, and never receive a fabricated
  integer `(n1,n2)`; the JSON reports `on_grid:false, n:null` for them.

For a legal sheared geometry, `Lx,Ly >= 2` and `0 <= shear < Lx`. The period
matrix has columns `(Lx,0)` and `(shear,Ly)` in cell coordinates, so the
allowed grid is

```
u = n1/Lx
v = n2/Ly - shear*n1/(Lx*Ly),  0 <= n1 < Lx, 0 <= n2 < Ly
q_index = n1*Ly+n2
```

`n1,n2` are labels of the dual character grid, not independent rectangular
fractions when `shear!=0`. The constraints are `Lx*u in Z` and
`shear*u + Ly*v in Z`. The all-mode sheared transform is an x transform,
the phase `exp(+2*pi*i*shear*n1*y/(Lx*Ly))`, then a y transform; it is not a
complex boundary hopping. The selected path always uses direct sums with the
same quotient character table. `minusQPoint` resolves the negated point by
this table instead of negating `n2` independently.

## 2. q-file format (`--vbs-q-file`)

Each non-comment line is `label u v` with `Q = u*b1 + v*b2` in fractional
reciprocal coordinates of the code (`b1=(2π/√3,2π/3)`,
`b2=(−2π/√3,2π/3)`, `a_i·b_j=2πδ_ij`). Blank lines, `#` and `//` comments are
allowed; tabs and commas count as spaces. Labels use `[A-Za-z0-9_.+-]`.
Parsing rejects: non-numeric/extra/missing fields, non-finite values, empty
lists, empty/forbidden labels, duplicate labels and duplicate momenta after
normalization into `[0,1)` (e.g. `0.25 0` vs `1.25 0`). See
`q_points_example.txt` for an L=4, Jz=0.5 example and
`predicted_q_L4_Jz0.5/` for generated candidates (from
`predict_vbs_momenta.py`).

Per-point JSON metadata: `input_fraction`, `reciprocal_fraction` (normalized),
`phase_fraction` (snapped to exact grid fractions when on-grid), `cartesian`,
`on_grid`, `n` (only when on-grid), `grid_snap_residual`, and the window
`W_L` plus `W_L_abs2=|W_L|^2`.

## 3. Definitions and normalization

With `X_mu(R)` the phonon/bond field of direction `mu` in cell `R`
(cell `= x*Ly+y`, bond `= 3*cell+mu`, **no** bond-midpoint or sublattice form
factor):

```
F_mu(Q)   = sum_R exp(-i Q.R) X_mu(R)          (finite-window Fourier sum)
Delta_mu(Q) = F_mu(Q)/Nc
S_mu,nu(Q)  = <F_mu(Q) conj(F_nu(Q))>/Nc       (complex, kept un-symmetrized)
mz2(Q)      = S_zz(Q)/Nc
W_L(Q)      = (1/Nc) sum_R exp(-i Q.R)         (Dirichlet-kernel window)
```

`W_L(Q)` is computed in closed form
`prod_axis sin(pi f L)/sin(pi f) * phase`, exact 0 on every non-zero PBC grid
momentum and exactly 1 at `Q=0`; on-grid points use snapped phase fractions so
this holds to roundoff.

**Projected (`prj`) quantities** remove the uniform background per
configuration:

```
F_proj_mu(Q) = F_mu(Q) - W_L(Q) * F_mu(0)      (per configuration)
```

For a purely uniform field this vanishes for every `Q` including `Q=0`. This
is **not** the ensemble connected correlator: the ensemble connected
`C_connected` subtracts `mean(Delta_mu(0))*mean(Delta_nu(0))` inside the joint
jackknife and exists only in `all`/`both` (grid) mode. In `selected` mode the
JSON states this explicitly.

## 4. Output layout

- `all` (unchanged, verbatim names/order): `12*Nc` raw entries — `9*Nc` `S_q<mu><nu>`
  then `3*Nc` `Delta_q<mu>`; derived `22*Nc` — per `q`
  `TrS_q/mz2_q/leading_eigenvalue_q/hermitian_residual_q`, then `C_d<mu><nu>`
  and `C_connected_d<mu><nu>`. The leading eigenvalue is computed from the
  Hermitian part of the **ensemble-averaged** `S` (negative values are
  reported, never clipped); `hermitian_residual` is `||S-S^dagger||` of the
  ensemble matrix.
- `selected` (`24*Nq+3` raw): raw `S(9Nq)` → raw `Delta(3Nq)` → projected
  `S(9Nq)` → projected `Delta(3Nq)` → background `Dsel0_<mu>` (= `F_mu(0)/Nc`).
  Names: `Ssel_s{j}_{mu}{nu}`, `Dsel_s{j}_{mu}`, `Sprj_s{j}_{mu}{nu}`,
  `Dprj_s{j}_{mu}`, `Dsel0_{mu}` with `j` the q-file order. Derived
  (`8*Nq+3`): raw and projected `TrS_s{j}/mz2_s{j}/leading_eigenvalue_s{j}/
  hermitian_residual_s{j}` each `4*Nq`, then
  `candidate_max_mz2`, `candidate_max_winner_index`, `candidate_participant_count`.
- `both`: the complete `all` block (identical names, order and values) followed
  by the `selected` block.
- Statistics: all outputs go through the existing phase-reweighting
  (`sum_phase * observable / sum_phase`) and the **joint delete-block
  jackknife** (`summarizeDerivedBlocks`); the measurement draws **no** random
  numbers.

### Legacy `candidate_max_*` semantics and migration

`candidate_max_mz2` is the maximum of the **raw** `mz2` over the listed points
with `Q=0` excluded by construction (`candidate_participant_count` reports how
many participated). This is a **legacy mixed grid+off-grid raw maximum** and is
kept for compatibility; it is not a complete Brillouin-zone maximum and is not
an iVBS order parameter. The argmax is still re-evaluated inside every
jackknife replica, but `candidate_max_winner_index` and
`candidate_participant_count` are categorical labels: their value remains for
compatibility while `stderr_real`, `stderr_imag` and covariance are written as
`null`. No non-integer index error is reported.

Older selected-mode output could contain a full-grid-derived key list paired
with a shorter selected-derived value vector. That positional mapping was not
reliable and must not be migrated by index; regenerate selected results. The
current selected `vbs.blocking[].derived` object has exactly `8*Nq+3` named
entries. The all layout and the all prefix of both remain unchanged.

### Structured candidate diagnostics

Selected and both modes add a separate `vbs.candidates` object with schema
`ivbs-candidates-v2` (`schema_revision=2`) and `level=finest_blocks`; it is not inserted into the
numeric raw or derived vectors. The three groups are deliberately disjoint:

| group | participants | score and interpretation |
|---|---|---|
| `grid_raw` | on-grid, nonzero candidates only | raw `Re[mz2]`; a maximum of the listed grid subset, not the full grid unless the list covers it |
| `offgrid_raw` | off-grid candidates only | raw `Re[mz2]`; finite-window search diagnostic that can contain uniform-background leakage `Xbar^2 |W_L|^2` |
| `offgrid_projected` | the same off-grid candidates only | projected `Re[mz2]`; removing the uniform mode does not remove window effects or noncommensurate mismatch |

Each point records raw/projected center values, imaginary diagnostics, separate
center/jackknife validity and reasons, real/imag errors, covariance,
`W_L_abs2`, group, normalized fractions and its `+/-Q` class. Scoring never
takes an absolute value and never clips negative real estimates. Q=0 is
excluded after normalization and the on-grid test; an empty group is
`available=false` with an empty participant/winner set and
`max_mz2_real=null`.

The center winner is a set selected from the full reweighted sample using the
recorded absolute/relative tie tolerances (`1e-12`, `1e-9`). Delete-one-finest-
block replicas reselect candidates and report unique-win and tie-membership
frequencies, with total, valid and invalid replica counts. A conventional
jackknife error is defined only when the complete original B replicas are
valid. If one replica is invalid, the center is retained but point
`stderr_*`/covariance are JSON `null`, `jackknife_defined=false`, and the
reason identifies the incomplete delete-one-block set. With at least two valid
replicas, selection frequencies may remain defined, but an invalid replica
sets `conditional_on_valid_replicas=true` and the frequency denominator is the
valid count. With fewer than two valid replicas, frequencies are `null`, not
zero. These are selection-stability diagnostics, not phase probabilities or
confidence intervals; tie-membership frequencies can sum to more than one.
The `+/-Q` class list retains original members and marks self-conjugate classes.

The three validity layers are independent: `center_defined` describes the
full-sample estimate, `jackknife_defined` describes the complete ordinary
delete-one-block error, and `stability_defined` describes the conditional
winner diagnostic. The JSON and CSV carry their reasons and the same replica
counts. A summary is `available` only when every non-empty group has all three
states; a center with missing jackknife or stability is `partial`; no defined
center is `unavailable`.

The lightweight reader `kitaev_ssh/summarize_vbs_candidates.py` writes one row
per point to CSV and a structured JSON summary without pooling files:

```bash
/home/asus/.venv/bin/python kitaev_ssh/summarize_vbs_candidates.py PREFIX_results.json
# default: PREFIX_candidates_summary.csv and PREFIX_candidates_summary.json
# explicit output and overwrite:
/home/asus/.venv/bin/python kitaev_ssh/summarize_vbs_candidates.py PREFIX_results.json \
  --csv OUT.csv --json OUT.json --force
```

The default name removes only the final `_results.json` marker and appends the
summary suffix, so dots in parameter names are preserved. Existing outputs are
rejected by default; `--force` is accepted only for explicitly named outputs.
CSV/JSON aliases, input aliases (including symlinks), and malformed or
non-finite input are rejected before any output is written. Writes use a
temporary file per output and replace the two files independently. The reader
reports `status=unavailable` for absent candidates, legacy v1, or an undefined
center, and `status=partial` for a defined center with unavailable ordinary
errors or stability. It does not infer frequencies from center values or
combine different `L`, `lambda`, `tau` or estimators. Off-grid raw/projected
values are supplementary diagnostics; neither one alone proves iVBS long-range
order or removes incommensurability.

## 5. Momentum prediction (`predict_vbs_momenta.py`)

`/home/asus/.venv/bin/python kitaev_ssh/predict_vbs_momenta.py --Lx 4 --Ly 4
--Jx 1 --Jy 1 --Jz 0.5 [--Jeff-x .. --Jeff-y .. --Jeff-z ..]
[--extra-q LABEL,U,V] --output-dir DIR --prefix P` writes `P_qfile.txt` (ready
for `--vbs-q-file`), `P_candidates.csv` and `P_report.json`.

The report and q-file are intentionally different layers. The report retains
every target and source label. Before writing the q-file, the generator
normalizes modulo reciprocal lattice vectors and emits each actual measurement
coordinate once. `aliases` and `sources` in `q_file_measurement_points` record
the merged labels. Distinct off-grid targets are not merged merely because
their nearest allowed grid point agrees. The q-file order is stable. The C++
parser retains its strict duplicate-label and duplicate-momentum rejection.
`--extra-q` grid metadata is computed with the same `(Lx,Ly,shear)` character
table as every other candidate.

`--uniform-effective --g G (--mass M --omega-ph W | --spring-k K)` adds the
low-risk uniform-background approximation. In the project convention
`H=sum u i (J_mu+g X_mu)c_A c_B`, midpoint Brillouin-zone quadrature evaluates
the free-Majorana bond mean and damped iteration solves
`Xbar_mu=-(g/k)<u i c_A c_B>_mu`. JSON records tolerance, maximum/actual
iterations, final residual, convergence status, `Xbar` and `Jeff`. A
non-converged solve exits with an error and writes no q-file. This is a
translation-invariant geometry guide, not an iVBS self-consistent saddle and
not the measured QMC structure-factor peak.

- Dirac nodes from `f(k)=Jz+Jx e^{-i k1}+Jy e^{-i k2}=0` (CONVENTIONS §6):
  closed form `p=(Jy²−Jx²−Jz²)/(2Jz)`, `z1=p+iy`, `z2=−Jz−p−iy`, verified
  numerically with a Newton fallback; outside the triangle `|Jx−Jy|≤Jz≤Jx+Jy`
  (or at the `Jz→0` nodal-line degeneracy) it reports **no isolated Dirac
  points and emits no Q**.
- Emits `Q_pred = ±2 K̃ mod 1` (continuous, may be off-grid), its nearest PBC
  grid points under the true reciprocal metric
  `|du b1+dv b2|² ∝ du²+dv²−du dv` (never per-coordinate rounding), the
  Kekulé `K±=(1/3,−1/3),(−1/3,1/3)`, and user extras; every row is labeled
  `bare` / `effective` / `kekule` / `user`. The disclaimer that bare
  predictions are **not** established interacting peak positions is written
  into every output file.
- Reports the mismatch `delta_i=(Q-q_L).(L_i a_i)=2*pi*L_i*Delta_i` (metric-minimal
  representative) in radians and in `2π` units, whether Kekulé is on-grid
  (`Lx % 3 == 0` and `(shear-Ly) % 3 == 0`), and warns when an iVBS nearest grid point
  coincides with a Kekulé nearest grid point (finite-size ambiguity).
- `compare_sheared_geometries.py --uniform-effective ... --geometry Lx,Ly,s`
  reports separate bare and uniform-effective target rows, nearest allowed
  momenta, raw mismatch, reduced-period shortest length/aspect ratio and the
  matching free-fermion `min |f(k)|`. It deliberately emits no best-geometry
  score or selection.
- L=4, Jz=0.5 artifacts are present under `predicted_q_L4_Jz0.5/`; they are
  regenerable and excluded from the incremental source patch. The bare
  prediction `(0.5804306,0.4195694)` and nearest grid point `q index 10`
  (`n1=n2=2` for the shear-0 L=4 example) agree with `q_points_example.txt` and with
  `check_bare_dirac_peak.py`.

## 6. Build and run

No-cmake path (`build_nodeps.sh` covers the 13 core and validation targets in
this tree; the CMake path also handles the reproduction targets):

```bash
export EIGEN3_INCLUDE_DIR=/home/asus/old-p/kitaev_ssh_pfqmc_hpc/third_party/eigen3
export BOOST_INCLUDE_DIR=/home/asus/.cache/pfqmc_boost_1_70_0
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export LD_LIBRARY_PATH=$MKLROOT/lib:$LD_LIBRARY_PATH
cd /home/asus/susy/pure-new
bash scripts/build_nodeps.sh /home/asus/susy/tmp/round2/final_build       # everything
bash scripts/build_nodeps.sh /home/asus/susy/tmp/round2/final_build kitaev_ssh_mc kitaev_ssh_vbs_q_validate
```

cmake path (GNU branch, needs MKLROOT in the environment):

```bash
cmake -S . -B build -DEIGEN3_INCLUDE_DIR=... -DBOOST_INCLUDE_DIR=... \
      -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran
cmake --build build --target kitaev_ssh_mc kitaev_ssh_vbs_q_validate
ctest --test-dir build -R "kitaev_ssh_vbs_q_validate|kitaev_ssh_p4"
```

Selected-mode run:

```bash
/home/asus/susy/tmp/round2/final_build/kitaev_ssh_mc --Lx 4 --Ly 4 --Jx 1 --Jy 1 --Jz 0.5 \
  --intervals 4 --dt 0.2 --mass 1.7 --omega-ph 0.8 --g 0.3 --seed 42 \
  --warmup 20 --sweeps 80 --measure-every 2 --block-samples 10 --vbs 1 \
  --vbs-mode selected --vbs-q-file kitaev_ssh/q_points_example.txt --output PREFIX
```

## 7. Changed and new files

| file | status | content |
|---|---|---|
| `kitaev_ssh/vbs_momenta.h` | new | `VBSMode`, q-file parsing, normalization, grid test, `W_L`, phase precomputation |
| `kitaev_ssh/vbs_observables.h` | rewritten | all/selected/both modes; `all` layout kept verbatim; `derive()` in `both` mode slices the grid block for `deriveAll` (bug fix, see §9) |
| `kitaev_ssh/coupled_mc.cpp` | modified | `--shear`, `--vbs-mode`, `--vbs-q-file`, mode-aware geometry JSON (`spatial_boundary:"pbc"`, `twist_applied:false` for no U(1) twist) |
| `kitaev_ssh/q_points_example.txt` | new | L=4, Jz=0.5 example list |
| `kitaev_ssh/vbs_q_validate.cpp` | new | measurement-layer test binary (7 categories, §8) |
| `kitaev_ssh/vbs_candidates.h` | new | disjoint candidate groups, tie sets, +/-Q classes and finest-block stability diagnostics |
| `kitaev_ssh/summarize_vbs_candidates.py` | new | one-file CSV/JSON reader for `vbs.candidates` |
| `kitaev_ssh/predict_vbs_momenta.py` | new | momentum prediction tool |
| `kitaev_ssh/VBS_SELECTED_MOMENTA.md` | new | this document |
| `kitaev_ssh/TWIST_FEASIBILITY_REVIEW.md` | new | Task-B feasibility review (no code) |
| `scripts/build_nodeps.sh` | new | cmake-free build mirroring the CMake GNU branch |
| `scripts/regress_vbs_modes.sh` | modified | square/non-square trajectory, layout, candidate JSON and postprocessor regression |
| `CMakeLists.txt` | modified | `kitaev_ssh_vbs_q_validate` target + ctest entry |
| `kitaev_ssh/predicted_q_L4_Jz0.5/` | generated | prediction artifacts (regenerable, not part of the source patch) |

## 8. Tests and results (2026-09-20)

`kitaev_ssh_vbs_q_validate` (43 checks, all tolerance `5e-12`, PASS):

1. **PBC grid consistency** (5×6 and 4×6): selected list over all/most grid
   points; raw `F`/`Delta`/`S` and derived `mz2`/leading eigenvalue agree with
   the all-mode FFT and with independent direct sums (max `3.0e-14`); `both`
   prefix `12*Nc` bitwise equals `all`.
2. **Uniform field**: raw `Delta = v_mu*W_L(Q)` and raw `S =
   Nc v_mu conj(v_nu)|W|²`; projected block vanishes for every `Q` including
   `Q=0`; `Dsel0` equals the uniform background (max `7.1e-15`).
3. **Known modulation** `X=X0+A cos(Q0·R+φ)`, on-grid and off-grid `Q0`:
   matches an independent direct sum and the analytic window sum
   `F(Q)=Nc[X0 W(Q)+(A/2)(e^{iφ}W(Q−Q0)+e^{−iφ}W(Q+Q0))]` (max `7.3e-15`);
   on-grid amplitude `Nc*A/2` exact; off-grid amplitudes deliberately not held
   to `A/2` (window + leakage).
4. **Conjugation for real fields**: `F(−Q)=F(Q)*`, while the linear relation
   `S(−Q)=S(Q)ᵀ` is preserved by the same complex reweighting; the
   conjugation/Hermitian residual can remain nonzero at finite sample size.
   `Szz(−Q)=Szz(Q)`, `W(−Q)=W(Q)*` (max `4.2e-15`).
5. **Layout/statistics** (`Nq≠Nc`): names/measure/derive dimensions for the
   three modes; `both` name order; complex-phase joint delete-block jackknife
   with direct raw-ratio and derived-`mz2` references and correct jackknife
   error/covariance; blocking preserves total ratios; candidate count; 12
   q-file error cases all throw (max `1.8e-15`).
6. **RNG trajectory regression**: same-seed chains in `all`/`selected`/`both`
   produce bitwise-identical acceptance counts, center-node coordinates, log
   weights, matter phases and the 7 basic observables (exact equality); the
   selected grid point reproduces the all-mode FFT value.
7. **Compile + QMC smoke**: three-mode smoke at fixed seed (acceptance
   `0.99599` identical across modes, basic observables identical, `both`
   contains the old raw+derived names and the new selected names) plus
   non-square-lattice runs in the regression sweep.

8. **Candidate diagnostics**: a large uniform background plus weak grid
   modulation keeps grid and off-grid groups separate and exposes finite-window
   leakage; q-file reordering preserves label-based maxima, +/-Q classes and
   frequencies; exact ties produce a two-member winner set with zero unique-win
   frequency; self-conjugate and switching cases are covered; empty groups have
   null semantics; hand-built complex blocks verify the reweighted center,
   imaginary diagnostics, complete-jackknife errors/covariance, one invalid
   replica, zero/one/two valid replicas, undefined centers and non-finite
   values; the finite-sample conjugation residual is kept visible rather than
   symmetrized.

9. **Mode regression and postprocessing**: `scripts/regress_vbs_modes.sh` passed
   L=4x4 (`Nc=16`) and L=5x3 (`Nc=15`) all/selected/both runs. Acceptance was
   `0.99599` and `0.99606`; phase/count/basic-observable columns were bitwise
   identical across modes. The script checked selected `8*Nq+3` derived keys,
   both `22*Nc+8*Nq+3` keys and unchanged all prefixes, candidate partitions,
   categorical null errors, v2 validity states, `W_L_abs2`, and readable
   CSV/JSON summaries. The standalone postprocessor tests also cover decimal
   filename collisions, explicit overwrite policy, symlink aliases, malformed
   arrays, legacy v1, and partial/unavailable output states.

An independent all-mode key audit confirms every original key is present
(`missing_original_keys: []`). The validation logs are retained under
`/home/asus/susy/tmp/round2/validate/`; the final build log is
`/home/asus/susy/tmp/round2/build_all.log`.

The postprocessor's missing-candidate path was also exercised on an all-mode
JSON. It returned `status: unavailable`, zero rows, and did not synthesize
candidate frequencies; the outputs are under
`/home/asus/susy/tmp/round2/unavailable/`.

## 9. Bugs found and fixed during this pass

The interrupted handoff smoke left `m*.json` truncated; re-running it showed
`--vbs-mode both` failing with `VBS ratio dimension` after the blocks were
written. Cause: in `both` mode `derive()` passed the concatenated
(all+selected) ratio vector to `deriveAll`, which requires exactly `12*Nc`
entries. Fix in `vbs_observables.h`: the `both` branch now slices off exactly
the grid block for `deriveAll` while `deriveSelected` consumes the tail at
offset `allStride(Nc)`; each helper keeps its strict dimension check. After
the fix the three-mode smoke passes end to end. The smoke script's layout
assertion was also corrected to look for derived names (`C_q0_xx` etc.) in the
`vbs.blocking[].derived` section where they are actually written.

The candidate winner helper also had an index/compact-score mismatch: group
participants are global q-file indices while group scores are compact. It now
maps score positions back to global indices before applying the tie rule; the
targeted order and tie tests cover this case.

## 10. Not verified / out of scope

- No production-scale physics run is part of this delivery; all numbers above
  are synthetic-field or small-chain regression checks (tests, not physical
  plots).
- The selected-mode estimator does **not** remove the incommensurability
  mismatch and does **not** change the allowed momentum grid of the selected
  integer-shear geometry; interpreting
  off-grid maxima requires the finite-window/leakage analysis (§3) and, in
  particular, off-grid amplitudes must not be read as `A/2`.
- `candidate_max_*` keeps the selection bias (documented, by design); use the
  structured `vbs.candidates.groups` object for new analysis.
- The real-space connected correlator exists only in `all`/`both` mode.
- Twisted boundary conditions: **not implemented**; see
  `TWIST_FEASIBILITY_REVIEW.md` for the feasibility review and its explicit
  "not implemented / not proven" conclusion.
