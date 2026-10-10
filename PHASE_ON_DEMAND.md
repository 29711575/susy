# Phase-on-Demand (Deferred Pfaffian Phase) — EXPERIMENTAL Design Note

Status: **validated optimization, default ON**. Enabled explicitly with `--defer-phase 1` on
`kitaev_ssh_mc` (sets `PhononWeightOptions.defer_phase_measurement=true`), or
`PFQMC_DEFER_PHASE=1` in the profiling driver. The eager/full-phase code path
remains available with `--defer-phase 0`; the current CLI default is deferred mode.

This note records the correctness argument, the state machine, the failure
policy and the validation evidence. It is deliberately separate from the
production optimization report to preserve the historical design and validation scope. The original default-selection
discussion belongs to its recording date; current defaults are listed in README.md.

## 1. What the production chain actually samples

The coupled phonon–matter Markov chain has proposal weight

```
W(X) = exp(-S_ph(X)) * w_matter(X),   w_matter(X) = |w_matter(X)| * phase(X)
```

with `phase(X)` the Pfaffian/sign phase of the signed Majorana overlap.

In `PhononSampler::attempt` / `attemptBond` (`kitaev_ssh/phonon_mc.h`) the
Metropolis decision uses only

```
logRatio = candidate.weight.log_abs_matter - current.weight.log_abs_matter
           + (candidate.log_phonon - current.log_phonon)
accept   = log(u) < min(0, logRatio)
```

**The Pfaffian phase never enters the accept/reject decision.** The phonon
action and the symmetric free-phonon proposal ratio are real and identical in
both engines (the local conditional Gaussian proposal cancels the phonon
Gaussian weight; only the changed bond/coordinate contributes
`log_phonon`, which is unchanged). Therefore the Markov chain's stationary
distribution for the phonon coordinates is proportional to

```
exp(-S_ph(X)) * |w_matter(X)|,
```

i.e. it samples the **magnitude** measure. The phase is used only at
measurement: at fixed strides `statistics.add(s.weight.matter_phase, measure(s))`
forms the phase-reweighted estimator

```
<O> = sum_n phase(X_n) O(X_n) / sum_n phase(X_n)
```

with numerator and denominator accumulated jointly (the existing complex
statistics). This is the standard "sample |W|, reweight by W/|W|" estimator
and is already what the eager code implements — the eager path merely
computes the phase earlier than strictly necessary.

## 2. Why deferring the phase is an exact rescheduling, not an approximation

For every proposal the deferred path computes:

* the same QR propagation of the trial subspace and the same `log|det R|`;
* the **same rank / reciprocal-condition / zero-overlap certification gates**
  (JacobiSVD singular values, identical thresholds), so a proposal that is
  `zero_weight` / `overlap_untrusted` / `propagation_failure` in the eager
  path is rejected/certified identically in the deferred path;
* `log|w|` from the same telescoping overlap ratios, with `log|det overlap|`
  obtained as `sum log sigma_i` (singular values) instead of a separate
  determinant LU — mathematically `|det M| = prod sigma_i`.

It omits, per proposal, work that only feeds the phase or the measurement:

* per-slice Cayley solves (`sw_cayley_prepare/solve`);
* intermediate/partial Green solves (`sw_partial_green`);
* the Pfaffian block and sign branch (`sw_pfaffian`);
* the final Green (`sw_final_green`);
* the endpoint center-Green rebuild (`center_green`).

At each fixed measurement stride the **current committed configuration** is
evaluated with the unchanged full reference routine
(`pureProjectorStableReferenceWeight` + `pureProjectorEndpointRebuild`),
yielding exactly the same phase and center Green the eager engine would have
produced for that configuration. Because

1. the magnitude walk is the same Markov chain (same proposal, same
   log-ratio, same uniform RNG, same accept/reject), and
2. the phase at a configuration depends only on that configuration and is
   reconstructed independently from it (no sign branch carried along the
   Markov trajectory),

the sequence of measured `(phase, O)` pairs is identical to the eager chain
up to roundoff. No determinant-square-root principal value is used; no sign
is continued between configurations.

## 3. State validity and transactions

`PathWeight` carries three validity flags:

| flag | meaning |
|---|---|
| `magnitude_valid` | `log_abs_matter`, `log_phonon` certified and usable for Metropolis |
| `phase_valid` | `matter_phase` is the full Pfaffian phase of this snapshot |
| `green_valid` | `center_green` is the rebuilt center Green of this snapshot |

Eager snapshots have all three true. Deferred proposal/committed snapshots
have only `magnitude_valid`; `matter_phase` is set to a zero sentinel and
`center_green` is empty.

* `Proposal::phaseRatio()` throws `std::logic_error` if either snapshot's
  phase is invalid — there is no default `phase=1`.
* **Reject** leaves the committed snapshot untouched; a previously measured
  phase/Green stays valid and is reused (consecutive rejects at one
  measurement point are counted repeatedly per the original fixed-stride
  rule, with the same cached phase — this is the same configuration repeated,
  as required).
* **Accept** atomically swaps in the candidate snapshot, whose phase/Green
  are invalid; the next measurement calls `ensureMeasured()` first.
* `CoupledPhononPath::ensureMeasured()` recomputes the full signed weight and
  center Green of the **current** configuration via the reference path
  (A1/A2 prepared data with all residual/rcond gates) and atomically replaces
  the immutable snapshot only on success. It takes no RNG.
* Measurement cadence is unchanged (`(i+1) % stride == 0`); the only change
  is that `PhononSampler::sample()` calls `ensureMeasured()` immediately
  before adding the sample. The same configuration appearing at several
  strides through consecutive rejects is re-added each time, as before.

Candidate and accepted data remain strictly separated: deferred candidates
are built by `fullReference`/`propose` on private snapshots; a rejected
candidate is destroyed without touching the committed state.

## 4. Failure policy

The magnitude path keeps every numerical gate. A rank loss, rcond below
threshold, non-finite propagation or zero overlap produces the same failure
status as eager (handled by the existing untrusted-candidate logic; it is not
silently turned into a reject-with-good-data).

At measurement, `ensureMeasured()` runs the authoritative reference
evaluation. If it fails (e.g. Pfaffian untrusted), the exception propagates
and the run aborts with the existing diagnostic path — the configuration is
**not** assigned `phase=1`, skipped, zeroed or silently dropped. There is no
second fallback beyond the reference routine because the reference routine
is itself the certified algorithm; an uncertifiable configuration at
measurement is a hard error exactly as it is at startup in the eager engine.

No audit or recomputation consumes the proposal RNG.

## 5. What was validated

`kitaev_ssh/deferred_phase_validate.cpp` (build target
`deferred_phase_validate`) drives an eager and a deferred engine through the
same bond-path/single-node proposals with the same fixed uniform `u`, and
checks per proposal and per committed state:

* `log_abs_matter` agreement (max abs error);
* identical Metropolis accept/reject decisions;
* `phaseRatio()` refusing invalid snapshots;
* after `ensureMeasured()`, phase complex distance and center-Green matrix
  norm agreement with eager;
* accept invalidates phase/Green, reject preserves a measured phase/Green,
  repeated measurement is stable;
* the initial deferred snapshot is magnitude-only.

Run across six configurations (seeds 7911–7916; g = 0.3/0.6/0.9/1.2;
intervals 48/64/96/128; up to 10 rejects per case): all pass,
`max_log_abs_error <= 1.5e-13`, `max_phase_error <= 1.9e-15`,
`max_green_error <= 5.1e-15`, `max_acceptance_error <= 1.5e-13`.

End-to-end (`kitaev_ssh_mc`, 3x3, intervals 32, bond-refresh, same seed,
stride 2): eager and deferred JSON outputs are **bit-identical in every
observable, phase diagnostic and block** (0 numeric differences; same
acceptance and accepted counts).

## 6. Known coverage limitations

* All certifiable configurations in these tests had Pfaffian phase `+1`
  (`negative_phase_samples = 0`, `positive_count = 20` in the end-to-end
  run). The **negative-sign and genuinely complex phase branches were not
  exercised** because no such physical configuration was generated; the
  sign machinery itself is the unchanged eager code, but a constructed
  negative-phase regression sample is left as future work and must be added
  before this mode is considered for a production default.
* The mode was tested with the combined/merged engine and bond-refresh plus
  single-node proposals; other weight engines were not separately measured.
* Only short/medium chains were run in this round; no equilibrium, ESS or
  integrated-autocorrelation claim is made for the deferred estimator (it is
  the same estimator as eager, but production-length evidence is absent).

## 7. Code locations

* `inc/pure_projector_fast.h`: `pureProjectorMagnitudeWeight`.
* `kitaev_ssh/phonon_projector.h`: validity flags in `PathWeight`,
  `defer_phase_measurement` in `PhononWeightOptions`,
  `evaluatePathMagnitude`, `CoupledPhononPath::ensureMeasured`,
  deferred branches in `fullReference`/`propose`, guarded `phaseRatio`.
* `kitaev_ssh/phonon_mc.h`: `PhononSampler::ensureMeasured` and the
  measurement call in `sample()`.
* `kitaev_ssh/coupled_mc.cpp`: `--defer-phase 0|1` (default 1; eager mode is 0)
  banner and `defer_phase` JSON field.
