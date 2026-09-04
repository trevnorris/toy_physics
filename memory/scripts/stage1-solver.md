---
title: Stage-1 solver script catalog
type: script_catalog
status: current
sources:
  - software/stage1_solver/src/stage1_solver/harness.py
  - software/stage1_solver/src/stage1_solver/branch_harness.py
  - software/stage1_solver/src/stage1_solver/convergence_harness.py
  - software/stage1_solver/README.md
  - software/stage1_solver/STAGE1_VERDICT.md
last_updated: 2026-09-03
---

# Stage-1 solver scripts

Run these from the repository root with
`PYTHONPATH=software/stage1_solver/src`. Generated arrays, checkpoints,
manifests, and convergence tables belong under ignored run directories.

## `python3 -m stage1_solver.harness`

- **Entry point:** `software/stage1_solver/src/stage1_solver/harness.py`.
- **Role:** Step-1/2 numerical validation gate.
- **Runs:** known-answer GPE benchmarks and manufactured-solution tests for
  the five operator blocks.
- **Output:** the configured Markdown validation report; exits nonzero if
  either benchmark group fails.
- **Meaning:** validates the discrete machinery against known continuum
  operators. It does not validate the toy model's physical assumptions.

Source: `software/stage1_solver/src/stage1_solver/harness.py` — function
`main`; `software/stage1_solver/README.md` — heading `What the numbers mean — and what they do not`.

## `python3 -m stage1_solver.branch_harness`

- **Entry point:** `software/stage1_solver/src/stage1_solver/branch_harness.py`.
- **Role:** Step-3 coupled matter/localized-Maxwell engineering smoke.
- **Inputs:** a deterministic `HarnessConfig` with placeholder, target-blind
  engineering parameters.
- **Output:** `software/stage1_solver/reports/step3_coupled_branch_smoke.md`;
  exits nonzero when the engineering gate fails.
- **Meaning:** exercises coupled residual/continuation machinery, not a
  physical branch verdict.

## `python3 -m stage1_solver.convergence_harness`

- **Entry point:**
  `software/stage1_solver/src/stage1_solver/convergence_harness.py`.
- **Role:** Step-4 coupled-branch grid-convergence study.
- **Outputs:**
  `software/stage1_solver/reports/step4_convergence_study.md` plus an ignored
  machine-readable convergence table under `runs/`; exits nonzero if the
  convergence gate fails.
- **Meaning:** tests numerical self-convergence of the engineering branch.

## Current-result boundary

The three entry points above describe the original validation and engineering
sequence. They should not be used to infer the current scientific outcome.
The later frozen effective-closure run is summarized in
`software/stage1_solver/STAGE1_VERDICT.md`: that particular branch has a
provisional robust magnitude miss, while the full Path-A model remains open.

Related sources:

- `memory/sources/software/stage1-solver.md`
- `memory/sources/software/stage1-verdict.md`
