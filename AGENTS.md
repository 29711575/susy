# Repository instructions

This repository exports the quick1010 source baseline selected on 2026-10-10, updated on 2026-10-13. Read README.md and HANDOFF_REPORT_CN.md before changing measurement or sampling behavior. Historical local instructions are retained in docs/provenance/quick1010-local-AGENTS.md for provenance; machine-specific paths in historical documents are not repository destinations.

Preserve physical parameters, Hamiltonian, trial state, gauge/sector and measurement definitions unless explicitly requested. Keep phase-reweighted raw moments and independent-chain boundaries. Do not infer convergence or a physical peak from build success, acceptance or a short validation run.

Future ordinary whole-bond production runs explicitly select `--proposal bond-guided --bond-guide bare --bond-beta 1`, `--q-mode-every 1 --q-mode-beta 1`, combined weight engine and enabled defer-phase/prefix/Cayley caches. For J=(1,1,.5), shear0, the fixed target is n=(3,3) for L6 and n=(5,4) for L9. Recompute q for other geometries. Optional q pools and multi-bond masks are not universal defaults. Controls retain explicit baseline settings; already-submitted jobs retain their frozen versions. New amplitude observables are separate from the historical linear chi and require explicit measurement options.

For HPC operations, prefer independent-chain CPU parallelism with BLAS/OpenMP one thread per chain, PBS arrays, no explicit memory request, and walltime720h unless overridden. Calculations longer than5min use the configured port20015 environment; this GitHub repository contains no credentials. Export plots as PNG.

Make scoped changes and retain source/measurement provenance. Do not rewrite Git history or publish changes without explicit authorization; an explicit user request to upload or replace repository contents provides that authorization.

Future ordinary production builds explicitly pass `--center-green-fast 1`; intentional controls pass 0. This only selects the correctness-validated center endpoint numerical implementation. Keep strict rank/rcond/residual gates and the uncached reference fallback. Existing frozen jobs retain their recorded versions. The previous GitHub snapshot is retained by tag `quick1010-20261010-before-center-green`.
