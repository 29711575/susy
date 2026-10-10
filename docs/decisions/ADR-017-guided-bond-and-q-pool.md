# ADR-017: fixed guided whole-bond pCN and fixed momentum pool

Date: 2026-10-10. Experimental opt-in extension of the uploaded core snapshot.
No old job or default proposal is changed. This is not a certification of L9
warmup, projector convergence, or average-phase statistics.

Integration note (2026-10-10): this implementation was reviewed and integrated
into quick1010. Native MKL guided-proposal and legacy CLI regressions passed;
the fallback build now registers the guided validator. Exact test scope and
source/binary provenance are recorded in
`../../../fast-1006/audit/1010/quick1010-guided-integration/REPORT.txt`.
The previous portable benchmark results remain supplied historical evidence.
New whole-bond production submissions explicitly use guided bare/beta1 per the
user's validated-optimization preference; the CLI default remains local.
Pool configurations are selected for a stated purpose, with no general
effective-sample performance claim. Existing submitted jobs retain their versions.

## Bottleneck and scope

A legacy `bond-refresh` draws the entire time path of one bond from the physical
free-phonon Gaussian. Coupled equilibrium paths can have large mean displacements.
Refreshing from the zero/endpoint-induced free mean then repeatedly proposes
unfavorable paths. A reference close to this mean can improve decorrelation
without increasing the number of bond attempts per sweep.

The existing fixed-q update remains available. It originally spans only one
z-bond q/-q pair. An optional fixed pool and a bond mask now allow systematic
coverage of several candidate subspaces. A pool is NOT chosen from the current
configuration's peak, and need not improve an already well-sampled single q.

## Guided bond update and exact target

The target remains `p(X) proportional to W_ph(X) |W_f(X)|`; the Pfaffian phase is
reweighted exactly as before. Neither phonon forces, mass, spring constant,
trial wavefunctions nor susceptibility definitions are modified.

Choose fixed real proposal-reference forces f_b. Let

    I_b(X) = sum_j w_j X_jb,
    w_0=w_M=tanh(omega*dt/2)/omega, w_internal=2*w_0.

The reference Gaussian is `p_ref(X) proportional to W_ph(X) exp(sum_b f_b I_b)`.
For each bond its mean is

    mu_ref(j,b)=mu_free(j,b)+f_b/k *
        [1-(exp(-omega*t_j)+exp(-omega*(Theta-t_j)))/2].

The covariance remains `C_jl=exp(-omega*|t_j-t_l|)/(2*m*omega)` including the
fixed open Gaussian endpoints. Propose

    X'_b=mu_ref_b+sqrt(1-beta^2)*(X_b-mu_ref_b)+beta*xi, xi~N(0,C).

This kernel is reversible with respect to p_ref. Hence

    log acceptance ratio = log|W_f(X')|-log|W_f(X)| - f_b*(I_b(X')-I_b(X)).

The source correction is mandatory. Accepting only with the matter ratio here
would sample the WRONG target. The full forward/reverse AR Gaussian density is
also computed for testing. The reference and beta are frozen, not adapted during
measurement. Inaccurate reference forces affect efficiency only.

`--proposal bond-guided --bond-guide bare` uses the *bare* finite-geometry free
fermion correlation to initialize the fixed proposal forces `f_b=-g<ui c_i c_j>`
(up to bond orientation conventions), averaged by bond type. This is not a
mean-field approximation of the target. Manual overrides and a zero guide are
available. Source-free/finite-physical-source cases use the same exact correction.

The proposal reference arrays are precomputed once. Each sweep still proposes
one whole path per bond, with the same collective-attempt cadence if enabled.
No Hamiltonian derivatives or intermediate Green functions are needed per attempt.

## Fixed q pool

`--q-mode-pool file --q-mode-bonds xyz` reads integer n1 n2 entries and builds
orthonormal real cos/sin modes for every requested bond type. Self-conjugate q
uses one real mode per bond. Gamma and q/-q duplicate subspaces are rejected.
The actual allowed momenta include the shear term. On a scheduled q update one
entry is chosen uniformly, independently of the current state; its full selected
subspace gets the existing free-reference pCN update. Acceptance remains the
matter ratio for these q proposals. Base guided proposals and q proposals use
different references and therefore DIFFERENT acceptance corrections.

The single-q default z mask preserves the legacy path and RNG stream. All-pool
counters include warmup+measurement and are labeled as such in results JSON;
legacy q/base stage counters are unchanged. q-diagnostic with a pool refers only
to its first entry, not to all pool members. Use all-q physical measurements for
scientific analysis.

## Failure semantics and tests

Same candidate / same uniform variate on one uncached retry, then fail closed;
no numerical failure becomes a normal Metropolis rejection. Deferred signed
weight/Green materialization, trial constraints and physical metadata are retained.
Tests cover full MH cancellation, source/endpoint covariance, Fock ratios,
noncommuting sheared honeycomb, accepted/rejected snapshot integrity, retries,
invalid input, multi-bond orthonormality, g=0 pool acceptance, and a linear model
where the guided beta=1 proposal is an exact independent heatbath.

Legacy local/bond-refresh runs (q on and off) were compared with the original
executable in an isolated OpenBLAS test build; raw blocks, paths and counts were
byte-for-byte identical. The original MKL backend still needs native regression.

## Performance interpretation

Measure both autocorrelation and seconds/sweep. Report fixed-q response moments,
not only the configurationwise maximum or acceptance rate. A small-system speedup
in the uniform phonon coordinate is not an established L9 critical VBS speedup.
Do not shorten the physical projection Theta to save warmup sweeps. Use identical
Theta, dt, Hamiltonian, trial, boundary and measurements in matched A/B runs.

Theory reference: Cotter et al., MCMC Methods for Functions, arXiv:1202.0709.
The shifted-reference acceptance formula above is derived for this source's
Mehler path measure; it is not copied from a finite-temperature periodic algorithm.
