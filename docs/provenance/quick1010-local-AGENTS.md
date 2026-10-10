# `quick1010` project instructions

- As requested on 2026-10-10, this directory is the baseline for subsequent
  code changes and optimization. It was extracted directly from
  `../fast-1006/audit/1010/core-code-for-gpt/fast-1006-core-current-20261010.zip`,
  removing the archive's outer `fast-1006/` directory. The inherited documents
  describe the source snapshot; historical fast-1006 paths and evidence remain
  references to the original workspace. See `BASELINE_ORIGIN.txt` for checksums.

This file is the short project entry point. Read the linked document only when
the task needs its details.

## Always apply

- Follow the workspace rules in [`../AGENTS.md`](../AGENTS.md).
  The user selected `/home/asus/newsusy/quick1010` on 2026-10-10 for subsequent
  code changes and optimization. Preserve fast-1006, fast-1002 and fast-0929
  as source and evidence references.
  Preserve existing work; do not clean, reset, overwrite, or auto-commit.
  Keep `/home/asus/newsusy/fast-0929` as the preserved source reference.
  Do not assume this directory is an initialized Git repository.
- The current priority is calculation speed while preserving physical correctness.
  Reduce unnecessary precision recovery and repeated runtime checks; use focused
  weight/phase/observable comparisons rather than expanding routine validation. The long-term driven finite-time-scaling goal is background,
  not a requirement for every task.
- Do not turn a historical test, run report, or user summary into a claim that
  this turn reran it. Keep paper facts, derivations, source observations,
  reported results, and proposed plans separate; see [`docs/SOURCES.md`](docs/SOURCES.md).
- Do not infer convergence, equilibrium, positivity, absence of a sign
  problem, a phase transition, SUSY, or computational advantage from a short
  run, phase near one, high acceptance, or a successful build.

## Read before acting

- Restoring project context: [`CODEX_HANDOFF.md`](CODEX_HANDOFF.md),
  [`docs/RESEARCH_CHARTER.md`](docs/RESEARCH_CHARTER.md), and
  [`docs/STATUS.md`](docs/STATUS.md).
- Locating code: [`docs/CODE_MAP.md`](docs/CODE_MAP.md).
- Measurement or statistics: [`docs/OBSERVABLES.md`](docs/OBSERVABLES.md),
  [`docs/NUMERICAL_INVARIANTS.md`](docs/NUMERICAL_INVARIANTS.md), and relevant
  files in [`docs/decisions/`](docs/decisions).
- Hamiltonian, geometry, or imaginary-time protocol:
  [`docs/PHYSICS_CONVENTIONS.md`](docs/PHYSICS_CONVENTIONS.md),
  [`docs/DRIVEN_FTS_METHOD.md`](docs/DRIVEN_FTS_METHOD.md), and ADR-004.
- If an HPC submission or subsequent job operation fails, read
  `/home/asus/newsusy/HPC_SUBMISSION_PLAYBOOK.md` for diagnosis or recovery;
  normal submissions do not require reading it first. Follow the root-level
  HPC rules in `/home/asus/newsusy/AGENTS.md`.

## Production optimization defaults (2026-10-10)

- The user requests that applicable optimization switches be enabled by
  default after correctness validation, including collective q updates.
  Use the latest verified frozen production binary; explicitly record all
  selected switches in new job scripts, parameter tables and manifests.
  Do not request another opt-in or wait for a long-chain convergence claim
  before using an already correctness-validated optimization.
- For new ordinary production chains, enable collective q updates with
  `--q-mode-every 1 --q-mode-beta 1`. Set `--q-mode-n1` and `--q-mode-n2`
  explicitly to the allowed nonzero momentum appropriate to the current
  lattice, shear and target observable. For the existing J=(1,1,0.5), shear0
  fixed-q studies, use (n1,n2)=(3,3) for L6 and (5,4) for L9. Recompute the
  allowed indices when geometry or the target momentum changes; never use
  Gamma or blindly reuse another lattice's indices.
- Continue the verified production acceleration settings (combined weight
  engine, deferred phase, prefix/Cayley caching) unless explicitly overridden
  or inapplicable. Diagnostic/measurement switches are selected for the
  requested output; physical-field/driven switches are not optimizations.
- The guided whole-bond proposal is integrated with native MKL correctness
  and legacy-compatibility validation on 2026-10-10. For new ordinary
  whole-bond production chains, explicitly select
  `--proposal bond-guided --bond-guide bare --bond-beta 1`, together with the
  existing single-q settings above. User-specified local/bond-refresh runs
  and baseline/control runs retain their requested proposal.
  The fixed q pool and multi-bond masks remain optional configurations;
  select them for a stated sampling purpose rather than assume universal
  improvement. L9 mixing, warmup length and effective-sample CPU gain have
  not been certified by this integration. Evidence is in
  `../fast-1006/audit/1010/quick1010-guided-integration/REPORT.txt`.
- Intentional baseline/control chains explicitly disable the feature being
  compared, e.g. `--q-mode-every 0`. Explicit user settings take precedence.
  Existing submitted jobs and historical frozen snapshots remain unchanged.
- The executable's current omitted-option default for q updates is still 0
  for interface compatibility. Future production submissions must explicitly
  pass the enabled settings above, so this older CLI default does not decide
  the production configuration. Keep provisional efficiency and convergence
  claims qualified as required under "Always apply".

## Project-specific boundaries

- `connected` and `projected` are distinct observables.
- Average the complete matrix and apply the common phase/background treatment
  before eigenvalues, ratios, or fits.
- The QMC fixed-gauge matter calculation is not full physical-spin projection.
- A proposed driven protocol does not authorize implementation, QMC, HPC, or a
  commit. Keep constant-path behavior and correctness checks separate from
  future driven work.
- Preserve independent-chain boundaries in errors and reblocking.

## After a task

- Update [`docs/STATUS.md`](docs/STATUS.md) only when new evidence changes the
  recorded state.
- Update [`docs/CODE_MAP.md`](docs/CODE_MAP.md) when interfaces or file roles
  change, and the relevant ADR when a design decision changes.
- Put new run or correction material in a new audit directory; keep historical
  audits immutable.

Detailed physics definitions, statistical invariants, protocol design, and
scientific reporting limits belong in the linked documents, not in this entry
file.

# Audit output organization

- As requested on 2026-10-09, create future audit task directories under `fast-1006/audit/MMDD/` by the task date (Asia/Shanghai), for example `fast-1006/audit/1009/<task-name>/`. Keep existing historical directories in place unless the user requests migration.
