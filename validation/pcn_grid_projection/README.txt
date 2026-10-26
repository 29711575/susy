PCN / mesh efficiency / finite-projection verification

Production remains kitaev_ssh_mc. Optional --validation-series 1 writes compact
UNWEIGHTED per-measurement complex values and their common phase. Multiply each
value by that phase exactly once before summing. --path-series-bond 2 writes the
z bond at cell(0,0), all imaginary-time nodes. Both options default off, introduce
no RNG draws, and leave Hamiltonian, trial, gauge, estimators and sampler defaults
unchanged. On/off production block bytes agree in the executed native-HPC check.

Pure PCN CLI mapping:
 --proposal bond-guided --bond-guide zero --bond-beta BETA
 --q-mode-every 0 --measure-every 1 --path-series-bond 2
(See actual plan.json args for the complete executable command.)
BETA=.3,.6,1 is bond-beta, NOT q-mode-beta. No local, refresh, q or force update
is added. One MC step is one production sweep: each of12 independent g0 bond
paths is updated once, including the selected bond exactly once. Updates to
other independent g0 bonds do not modify its transition. The measurement is
then recorded. tau_int=1/2+sum(k>=1)rho(k), variance inflation2*tau_int; units
are selected-path updates, interval1. Analysis reports the finite analytic lag
window and infinite geometric answer separately. F and F2 use known analytic
zero mean / variance centering. Full path and innovation covariance, endpoint
moments, ACF/IAT and complete blocking curves are preserved.

Grid physics: L2,J(1,1,.5),g.4,m=omega=k1,tau_tot2,dt.2/.1/.05,trial0,
fixed all-uAB+1/Wilson(+1,+1), odd free matter trial tensor phonon vacuum.
pCN is zero-guide beta.6; guided is fixed bare-force beta.6; local has no beta,
and always uses the exact free Mehler conditional (its dt-dependent variance
is the algorithm definition, not a tuning change). q updates are disabled as
an intentional pure-algorithm control. Forward/reverse guided proposal densities
remain in the unchanged production MH implementation. JSON records effective
combined/defer/prefix/Cayley/center-Green settings, explicitly enabled in args.

Projection: same g.4 parameters, guided beta.6/qoff, tau_tot2/4/8/16,
dt.1/.05. Tau2 overlaps the grid cases. Match the identical-split center-cut ED
before interpreting continuous-time or same-parity ground values. matter_energy
is matter only, excludes oscillator kinetic and zero-point energy. Boson X/X2,
S,m2,m4,B4 and fermionic Ki=i*cA*cB / mixed XKi all enter the check. This small
lattice checks ingredients, not critical exponents or full physical-spin physics.

Generate definitions: python generate.py NEW_EXPERIMENT_DIRECTORY
Run explicit PBS arrays with PLAN=pilot-plan.json or pcn-plan.json. run.py calls
only the frozen production CLI. Independent pilot data determine one formal
allocation via plan_formal.py, then the immutable formal-plan.json is hashed
BEFORE qsub. Four fresh seeds/case; no pooling pilot/formal or old/new campaigns.
The cap is an ESTIMATED4CPU-hours/chain using pilot timing, including warmup;
formal sample counts are fixed and are not stopped by measured time or result.
Report resource-capped/insufficient precision as INCONCLUSIVE. No seed replacement.

report.py generates per-observable CSV, every blocking level, per-seed ACF and
execution records, acceptance, seconds/sweep and estimated CPU time to a COMMON
absolute error across meshes. It uses joint common-phase jackknife for Binder.
IAT of nonlinear estimates comes from their common-phase linear influence;
per-configuration m4/m2^2 is not a Binder time series. Time estimates require
resolved phase, blocking and seed diagnostics; they are conditional stationary
variance scaling, not an observed first-passage/stopping time. Total/measurement
cost and warmup are distinguished. Δt bias, projection bias, ED cutoff error,
solver residual and MC SE are separate columns. Physical projection convergence
in ED does not imply QMC can resolve the same small ground-distance bound.

All full12-phonon ED code is independent Python from validation/interacting_ed,
no production generator/Pfaffian/Green imports. Ncut4/6/8 are full total-occupation
cutoffs. Source/ground propagation residuals and observable convergence, not
boundary probability alone, gate comparisons. High moments use enlarged
intermediate boson bases. If cutoff fails, formal allocation stops rather than
blaming QMC. No finite-dt QMC machine-precision demand against continuous ED.

Executed task artifacts: audit/1019/pcn-grid-projection-validation/. Historical
PASS/INCONCLUSIVE reports under audit/1017 remain unchanged. Plots are PNG only.

Final fixed-budget result (2026-10-19)

60/60 production chains,8257536 measurements, all exit0/numerical_failures0.
Pure PCN PASS;13/15 interacting cases PASS, local dt.1/.05 INCONCLUSIVE.
319/345 observable gates PASS,14 unresolved-correlation and12 precision gates.
All means satisfy the preregistered six-SE compatibility gate; this does not
certify unresolved estimates. No added samples, replacement seeds or relaxed gates.
ED tau16 meets all23 ground tolerances, but QMC resolves only15/23(dt.1) and
10/23(dt.05) under the stricter ground-distance diagnostic. Full ground accuracy
is not certified merely by matching finite-tau ED.
Final report: audit/1019/pcn-grid-projection-validation/FINAL-REPORT.txt.
Complete CSV/blocking/ACF and frozen output: experiment/analysis/ under that audit.
closure/ has execution audit, raw SHA manifests, actual commands, per-seed IAT,
hardware/timing limitations, and an additional reference-error diagnostic; adding
continuous/ground ED cutoff differences changes no classification. Frozen outputs
and original1017 reports are retained. User-requested migration moved33 unstarted
chains to20015, while27 running chains stayed20022. Time costs span heterogeneous
hosts; no claim of isolated algorithm speedup follows from those time ratios.
