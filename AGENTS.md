# Repository instructions

This repository exports the quick1010 main source, updated on 2026-10-30.
Read README.md and docs/updates/20261030-source-and-results.md for current behavior;
HANDOFF_REPORT_CN.md and measurements preserve historical frozen results.
Local audit/build paths in technical documents are provenance, not repository destinations.

Preserve the Hamiltonian, gauge, geometry, physical parameters, measurement definitions,
joint complex-phase raw statistics and independent-chain boundaries unless requested.
Do not infer equilibration, mixing, projection convergence or SUSY from builds or acceptance.
Existing jobs retain their recorded source, executable and settings.

New production uses --matter-trial physical; unrestricted is an explicit historical control.
Use production-v1 numerical defaults and retain strict rank/rcond/structure/residual gates.
Ordinary whole-bond sampling explicitly selects --proposal bond-guided --bond-guide bare
--bond-beta 1. Explicitly enable nonzero-q updates with --q-mode-every 1 --q-mode-beta 1;
Legacy fixed proposal targets are L6=(3,3), L9=(5,4) for J=(1,1,.5), shear0.
These are sampling choices, not measurements of the peak or a nearest-grid rule.
The bare L6 valley-transfer nearest set is (2,3), (3,2), (3,4), (4,3).
Gamma updates and extra q pools remain optional; Gamma CLI default is off.
The multiscale Gamma-pm runtime used in some new measurements is a separate frozen
implementation; this exported core does not accept its collective-beta-scales CLI.

For new calculations explicitly set --fermions 1 --vbs 1 --vbs-mode all
--vbs-susceptibility 1 --vbs-amplitude 1 --vbs-amplitude-response 1 --measure-energy 1
--measure-total-energy 0. Total energy is enabled only when requested.
Explicit user settings and intentional controls take precedence.

Use verified frozen source snapshots, built for the target runtime, for production.
The historical local production manifest in docs/provenance does not select a current build.
Do not run preflight jobs or smoke chains, or gate production on them. Keep parameter,
file/executable checks, required target-runtime builds and production numerical gates.
New HPC submissions use port20015 by default; short calculations with at most14 CPUs
run locally. Explicit user port/queue choices override this default. Default walltime720h, no memory request.
Prefer PBS arrays and independent chains, one BLAS/OpenMP thread per chain.
Export plots as PNG. Keep credentials, builds, raw task trees and local audit outputs out of Git.

Checkpoint resume requires the original executable/runtime and a fresh output prefix;
resumed raw blocks already include their ancestor and must not be added to it again.
Make scoped changes and retain provenance. Do not rewrite Git history or publish without
authorization; an explicit user request to synchronize the repository authorizes publication.
