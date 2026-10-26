Independent interacting ED validation (2x2 periodic honeycomb)

A: kitaev_ssh_interacting_fock_validate, registered in scripts/build_nodeps.sh
and CMake BUILD_TESTING. It writes explicit paths and MH.csv. The original
kitaev_ssh_pfqmc_fock_validate must also be run: synthetic Fock phase witnesses
and physical honeycomb paths are separate coverage.

B: ed_reference.py independently constructs full Fock Majoranas, projects EVEN
operators to the odd trial sector, and retains all 12 quantum bond oscillators
with sum(n_b)<=cut. No production Hamiltonian, Green or Pfaffian imports.
Install the versions in requirements.txt in a separate environment. Set all
BLAS/OpenMP thread counts to one. Example:
 python ed_reference.py --cut 8 --g 0.4 --taus 2 4 --output ed-g0p4-cut8.json
Run g=0 at cuts 2,4,6 and g=.2,.4 at cuts 4,6,8; verify_ed.py summarizes
convergence, eigen/propagation residuals, fixed-trial source differentiation,
and the separate finite-projection bias. See --help for each program.

For finite dt, discrete_reference.py independently applies the symmetric
matter-half / exact-phonon / matter-half operator. Response uses the exact
number-basis exponential derivative. Run cuts 4,6,8 at each production point;
check_discrete_ed.py and check_trotter_order.py distinguish truncation and dt.
High moments use an enlarged intermediate space before projection.

C: prepare_qmc.py --output PACKAGE [--smoke] writes plan.json, points.tsv and
explicit CLI command files. Execute those with a freshly built production
kitaev_ssh_mc in PACKAGE/runtime, preserving runs/chainNNN* and chain boundaries.
The commands check runtime/binary.sha256. Smoke uses short chains and cannot
pass the full statistical gates. Full budgets and gates are in full-config.json;
no automatic seed replacement or tolerance relaxation is performed.
 compare_qmc.py PACKAGE --discrete-ed DISCRETE_DIR --continuous-ed ED_DIR --output OUT
uses joint common-phase block jackknife, connected subtraction and true raw
fourth moments. It checks reblocking, precision and phase resolution. Exit 2
means BLOCKED/INCONCLUSIVE; exit 1 means FAIL. Partial outputs remain evidence.

Canonical execution evidence is under audit/1017/interacting-ed-validation/.
The package/source.tar.gz is the frozen compiled source; later analysis files
have their own hashes. Native and HPC compilers/runtime hashes are recorded
separately. All parameters/seeds come from plan.json, never inferred defaults.
majorana_energy is matter energy, excluding oscillator kinetic and zero-point
energy. ED total energy is not a production total-energy estimator.
These tests do not prove physical-spin projection, thermodynamic convergence,
finite-g correctness outside tested coverage, or SUSY.

Final fixed-budget assessment: all28 chains (640000 production samples) complete.
A PASS; B finite-projection/cutoff PASS; C INCONCLUSIVE because four local
absolute-SE targets are unmet. Six cases pass; no six-SE discrepancy found.
Do not replace this status with PASS merely because every process exits zero.
plot_results.py C_DIRECTORY exports QMC-ED.png and local-precision-blocking.png.
See the audit REPORT.txt and C/C-summary.json for the complete assessment.
