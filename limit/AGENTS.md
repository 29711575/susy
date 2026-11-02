# Static/adiabatic package

Apply the repository AGENTS.md. This directory contains the separate finite-torus
static implementation and compact L6/L7/L8/L9 data, published on 2026-11-02.
Read README.md and PROVENANCE.json before calculations or interpretation.

Preserve the byte-identical recorded sources in source/ and accurate_response.py.
Changes or portable adapters must be separate and have their own provenance.
The explicit physical model fixes J=(1,1,.5), k=1, square shear-zero geometry,
zero flux and Wilson++. Legacy --parity trial is not the physical selection.
Use the signed-Hessian response adapter when extending physical calculations;
retain its error budgets and explicit failure at unresolved branches.

Store static outputs and task plots in results/; task validation records belong
in validation/. Keep the published historical tables and verification identities.
Do not replace them with new local optimization results or overwrite frozen sources.
Bare single-domain susceptibility is distinct from domain-switching response.
Finite square L is distinct from a periodic Bloch-reference cell period N.
No global-minimum, full-flux or thermodynamic critical-point claim follows from
these multistart local results. Plotting or documentation does not authorize new scans.
