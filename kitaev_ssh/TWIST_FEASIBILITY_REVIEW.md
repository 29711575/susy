# Feasibility review: continuous boundary twist for the iVBS order parameter

**Task-B document — review only, no code.** Status of every construction
discussed below: **尚未实现 / 尚未证明 (not implemented, not proven)**. This
delivery contains **no** continuous-twist switch: the JSON metadata keeps
`"spatial_boundary":"pbc"`, `"twist_applied":false` in every mode, and the
selected-momentum finite-window Fourier estimator must not be mistaken for a
twisted boundary condition.

This review predates neither the integer-shear implementation nor its fixed
quotient: `--shear s` is a discrete change of the spatial torus and can change
the allowed VBS grid. It is deliberately outside the continuous-twist topic;
the statement above remains that no U(1) phase is applied.

## 1. What the question is

The paper's field theory puts the incommensurate VBS transition at a
fractionalized QCP of a single Dirac cone coupled to a complex boson
`φ` carrying momentum `−2K̃`. On a torus, the SUSY* spectrum and the boson
momentum grid are tied to the spatial boundary conditions. A "continuous
twist" would mean a family of torus Hamiltonians `H(θ)`, `θ ∈ [0,2π)`, under
which the low-energy boson (and/or fermion) picks up `e^{iθ}` around a
non-contractible cycle, so that the allowed `φ` momenta shift continuously
off the PBC grid `Q = 2π(nx/Lx, ny/Ly)`. The motivation would be to tune a
finite lattice into exact commensuration with an arbitrary incommensurate
`2K̃`.

## 2. What the current code actually supports (evidence)

- Gauge/Wilson sectors: `HoneycombLattice::zeroFluxGauge(wx,wy)` accepts only
  `wx,wy ∈ {±1}` (it throws otherwise) and `wilsonLoop(direction,u)` evaluates
  the Kitaev Wilson loop in the four sign sectors. The production entry point
  fixes `(+1,+1)` (`"wilson_loops":[1,1]` in the JSON); `mc_validate`
  exercises `(-1,+1)` explicitly. There is no `U(1)` sector label anywhere in
  the gauge bookkeeping — the stored gauge field `u_b` is `±1` per bond
  (`validateGauge`).
- Boundary degrees of freedom in the fermion sector: a sign `±1` on each
  boundary-crossing bond (`boundary_wrap ≠ 0`) is already expressible; that is
  exactly what `zeroFluxGauge` does. Nothing else is parameterized.
- Phonon sector: each bond carries a **real** coordinate `X_b` on a real
  Gaussian path measure (open imaginary-time contour, PBC in space by bond
  identity). The phonon action and the trial endpoint factors are real and
  positive; the AR transition sampler (`freeConditional`, `freeBondLogDensity`)
  is built on real coordinates. No phase, no anti-periodicity, no complexified
  bond field exists in the code path.
- Measurement: `--vbs-mode selected` evaluates `Σ_R e^{-iQ·R} X(R)` at
  arbitrary `Q`. This is an **estimator** on the same PBC ensemble; it changes
  no weight, no boundary, and no allowed momentum of anything.

## 3. Why a continuous twist cannot be bolted on today

### 3.1 A single real Majorana hop carries only a sign

The matter Hamiltonian in the fixed gauge is `H = (i/4) c^T A c` with `A` a
**real antisymmetric** matrix (`A_AB = ±2 t_b`, CONVENTIONS §2). With the code
definition `K_code = -i A`, the matrix properties are
`K_code^T = -K_code` and `K_code^dagger = K_code`: **Hermitian and
transpose-skew**. This is the property checked by the source; it is not
skew-Hermitian.

The generators `c_i` satisfy `{c_i,c_j}=2 delta_ij` with `c_i^dagger=c_i`.
Consequently, multiplying a single Hermitian Majorana bond term by a scalar
phase, `i t exp(i phi) c_A c_B`, preserves Hermiticity only when the scalar
coefficient is real, `phi in {0,pi}`: a sign. This statement concerns one
Majorana degree of freedom and does not rule out a construction that mixes
several microscopic degrees of freedom.

For comparison, two Majoranas can be mapped to a complex fermion only after a
specific choice such as `f=(c_A-i c_B)/2`. For ordinary complex fermions on
sites `i,j`,
`exp(i theta) f_i^dagger f_j + h.c.` equals
`cos(theta)(f_i^dagger f_j+f_j^dagger f_i) +
i sin(theta)(f_i^dagger f_j-f_j^dagger f_i)`. The sine term is a
number-conserving current-like hopping combination, not pairing and not a
chemical-potential term. A claim that this is equivalent to the single-
Majorana construction would require an explicit mapping of the degrees of
freedom and lattice bonds; none is supplied here. The boundary data currently
implemented in this code remains `Z2`: per-cycle signs, i.e. the four Wilson
sectors.

### 3.2 Changing the fermion boundary does not move the phonon momenta

Suppose one flips boundary bond signs (or otherwise moves single-particle
levels within the `±1` sectors). The phonons are still real bond coordinates
with a real, positive Gaussian action, and the existing code still stores a
single-valued PBC configuration. A direct scalar boundary rule
`X(R+L_i a_i)=exp(i theta_i) X(R)` is incompatible with a nonzero real `X`
for a generic `theta_i`; this rules out that particular microscopic rule.

It does **not** rule out every way of changing the phase of a real modulation.
For example,
`X(R)=Xbar+A cos(Q.R+phi)` remains real when `phi -> phi+theta`. This is
different from the generally non-real field transformation `X -> exp(i theta)X`.
Nor does the scalar-real-coordinate observation rule out a future boundary
construction that mixes several legitimate microscopic degrees of freedom.
Such a construction would need to be specified and derived, rather than
asserted equivalent to the present model. For a fixed integer-shear geometry,
an off-grid estimator still has window leakage (`W_L(Q) != 0`) and mismatch
`delta_i=(Q-q_L).(L_i a_i)`; it is not a twisted ensemble. Changing the
integer shear selects a different quotient and can change the allowed grid,
but it is a discrete seam change, not a continuous phase twist.

### 3.3 The `U(1)` twist of `φ` is emergent, with no microscopic knob

The complex order parameter `phi ~ sum_R exp(-2 i Ktilde.R) X_z(R)` is a linear
combination of real bond fields. Its phase is an emergent low-energy
description, with an approximate `U(1)` rotation when the effective potential
depends only on `|phi|`. A boundary condition on that emergent field is not by
itself a legal microscopic boundary condition. The current code has no
microscopic continuous knob, but the distinction above means that the absence
of a scalar complex phase on `X` is not a proof that all multi-component or
phase-of-a-real-modulation constructions are impossible.

## 4. Checklist a real construction would have to satisfy

Any future proposal must demonstrate all of the following; each item is
currently unproven for this codebase:

1. **Boundary Hamiltonian.** An explicit `H(θ)` on the torus in which `θ`
   enters only boundary-crossing terms and `H(0)` is *exactly* the current
   SSH-Kitaev model (same couplings, same gauge bookkeeping, same phonon
   measure).
2. **Majorana algebra.** The single-particle generator must remain compatible
   with the code's real antisymmetric `A` and its **Hermitian and
   transpose-skew** `K_code=-iA`, preserving `{c_i,c_j}=2δ_ij`, Hermiticity of
   the many-body Hamiltonian and the fixed-gauge `Z2` flux/Wilson bookkeeping; Lieb's
   theorems (zero flux, mirror-symmetric phonon ground state) must be
   re-derived under `H(θ)` or their failure quantified.
3. **Phonon consistency across the seam.** Phonon coordinates stay real with
   a positive real action; the AR path measure and trial-state endpoint
   factors on seam bonds must remain well-defined and consistent with the
   twisted fermion hopping; the zero-flux/zero-phonon-configuration ground
   state statement must be re-examined.
4. **Microscopic action on the VBS field.** Give a legal boundary
   transformation of the actual microscopic degrees of freedom and derive its
   action on the complex iVBS amplitude, including the resulting momentum
   quantization. It is not sufficient to impose
   `X(R+L_i a_i)=exp(i m_i(theta))X(R)` on a real coordinate, and a changed
   finite-size peak is not by itself evidence that the momentum grid moved.
5. **Zero-twist and small-twist limits.** At `theta=0`, the original code path
   should be an exact fallback, including weights, acceptance and RNG usage.
   For nonzero small `theta`, require controlled convergence of the relevant
   physical observables and numerical diagnostics; do not demand bitwise
   equality to the theta-zero run.
6. **Trial states and Pfaffian phase.** The pfQMC weight uses a Gaussian trial
   state built from the free ground state of `K_code`. A continuous twist
   needs a continuous family `|ψ_T(θ)⟩`; the Pfaffian sign/branch tracking
   (`pfaffianForSignOfProductWithStatus`, the `(-1)^(NM/2)` correction) and the
   matter-parity policy must remain consistent across the family, and the
   `zero_mode_ambiguous` policy needs an explicit rule for the twisted ensembles
   where a Dirac point can cross zero as `θ` varies.
7. **Sign/weight behavior.** On intervals with no Pfaffian zero, level
   crossing, degeneracy, parity change or trial-state branch change, specify a
   branch convention and test controlled continuity of `log_abs_matter` and
   `matter_phase`. Global continuity through such singular or branch-changing
   points is not automatic and must not be assumed.
8. **Measurement.** Estimators must be re-derived for the twisted ensemble;
   in particular the finite-window Fourier sum is *not* the twisted-object
   measurement and cannot stand in for it.

## 5. Conclusion

- What exists today: `±1` Wilson sectors, real phonons on an integer-shear PBC
  torus, and a boundary-blind Fourier estimator at arbitrary `Q`.
- What does **not** exist: a continuous (`U(1)`) boundary twist in the current
  fermion or phonon code path, or a demonstrated microscopic mechanism that
  moves the VBS momentum quantization off the PBC grid. The present review also
  does not prove that every broader construction is impossible.
- Therefore: **尚未实现 / 尚未证明**. No continuous-twist flag, no
  twist-parameter CLI, and no re-labeling of off-grid Fourier measurements as
  twisted boundary conditions accompanies this delivery. A future proposal
  must provide the legal microscopic boundary map and its action on the iVBS
  amplitude. The JSON markers `"spatial_boundary":"pbc"`,
  `"twist_applied":false` are the enforceable statement of this fact.
