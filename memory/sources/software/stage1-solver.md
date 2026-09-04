---
title: Stage-1 solver and branch-realization groundwork
type: source
status: current
sources:
  - software/stage1_solver/README.md
  - software/stage1_solver/src/stage1_solver/harness.py
  - software/stage1_solver/src/stage1_solver/branch_harness.py
  - software/stage1_solver/src/stage1_solver/convergence_harness.py
  - software/stage1_solver/reports/pathA_23_stage0_action_and_contracts.md
  - software/stage1_solver/derivations/pathA_01_return_source_and_balance.md
last_updated: 2026-09-03
---

# Stage-1 solver and branch-realization groundwork

## Purpose and status boundary

This unit combines three related pieces: the original Stage-1 numerical
validation harness, a specification of a candidate brane-elastic parent
action, and a derivation of return-source/static-balance structure. It is
groundwork for a target-blind branch-realization test, not one uniform result.

The README's build table is an early snapshot: it records validation Steps
1–4 as complete while saying no physical branch had yet been solved. Later
work in this repository did perform a frozen branch test, so current outcome
questions must use `memory/sources/software/stage1-verdict.md`; the README is
still useful for the numerical architecture and validation history.

Sources:

- `software/stage1_solver/README.md` — headings `Build status`, `What to expect`, and `What the numbers mean — and what they do not`
- `software/stage1_solver/STAGE1_VERDICT.md` — heading `The result — a robust MISS`

## Numerical harness

The solver uses float64 PyTorch with a custom conservative finite-volume/
finite-difference discretization, full radial measure, and a matrix-free
Newton–Krylov solve with autodifferentiated Jacobian-vector products and an
Armijo line search.

- `harness.py` runs known-answer GPE benchmarks plus five manufactured-solution
  operator tests and writes the Step-1/2 validation report.
- `branch_harness.py` runs the Step-3 coupled matter/localized-Maxwell
  engineering smoke.
- `convergence_harness.py` runs the Step-4 coupled-branch grid-convergence
  study and writes a report plus an ignored machine-readable table.

The README reports second-order convergence for the five manufactured
operator blocks. This certifies the discretization of the encoded continuum
operators, not the physical correctness of those operators. Generated meshes,
arrays, and checkpoints belong in ignored run directories.

Sources:

- `software/stage1_solver/README.md` — headings `How to run`, `What to expect`, `Numerical approach`, and `Repo hygiene (firewall)`
- `software/stage1_solver/src/stage1_solver/harness.py` — function `main`
- `software/stage1_solver/src/stage1_solver/branch_harness.py` — function `main`
- `software/stage1_solver/src/stage1_solver/convergence_harness.py` — function `main`

## Candidate action and coupling contracts

The Stage-0 report specifies a `NEW_PARENT_ACTION` family containing the kept
bulk scalar sector, an off-brane displacement, and one of three in-plane
constitutive candidates: Cauchy/Navier, rotational/MacCullagh, or
Cosserat/micropolar. It deliberately does not select a winner or claim a
microphysical derivation. Its normal-work and conserved-current couplings are
classified as symmetry-allowed postulates.

The current-admissibility pre-check succeeds only for a conserved defect
current with the scalar-potential completion. That is a structural gauge
check, not a derivation of charge or Maxwell dynamics. Keeping this carrier
beside the already localized Maxwell sector would also require an explicit
double-counting resolution.

Sources:

- `software/stage1_solver/reports/pathA_23_stage0_action_and_contracts.md` — headings `VERDICT`, `Candidate in-plane constitutive menu`, `Microstructure contract`, `Coupling provenance`, and `Current-admissibility pre-check`

## Return source and static balance

The return-source note derives the quadratic matter-to-wall cross source
$S_\eta^{(\psi)}=-k_1\delta\rho$. The accompanying
$k_2\rho_0\eta$ term is a wall-stiffness renormalization, not a new return
channel. Equality of the forward and return cross-Hessians is the key
reciprocity statement.

The declared fixed-localization Maxwell action has no direct wall variation,
so gauge return is matter-mediated and left as an on-shell effective-action
variation rather than assigned an invented local formula. Closing the return
loop can change the self-consistent background and denominator-side
coefficients, but it introduces no independent numerator magnitude in this
quadratic derivation.

Sources:

- `software/stage1_solver/derivations/pathA_01_return_source_and_balance.md` — headings `D1: Matter Return Source`, `D2: Gauge Return Source`, `D3: Static Self-Consistent Throat Balance`, and `D4: No New Numerator Magnitude`

## Open questions

The constitutive winner, no-leak behavior, complete spectrum, charge
normalization, cone locking, and carrier double-counting remain outside the
Stage-0 result. The return derivation does not solve the static throat, prove
existence or uniqueness, derive nonlinear wall functions, or close the local
gauge return. These limitations should be kept separate from the later frozen
branch verdict.

## Related memory

- `memory/sources/software/stage1-verdict.md`
