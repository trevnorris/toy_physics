# Directive pathA_C0 — Conditioning spike: make the throat-profile solver REACH a deep/realistic throat

**Status:** DRAFT (Claude-authored 2026-06-19; pending Codex design-review → confirm-pass → execution — user gate). First
step of **option C** (the throat-profile solve). Follows the Phase-1 solver reconnaissance (decision-13 §0/§13): the PDE
operators are MODEL-FAITHFUL and the solver + calibrate-predict harness are ~70–80% built + validated, but the solver has
NEVER converged on a deep/realistic throat — the blocker is NUMERICAL CONDITIONING (B2c stalls at `τ≈0.029` with "line
search failed" = a diagnosed conditioning floor, NOT a physical edge).

**Date:** 2026-06-19
**Owner:** Codex (codes + iterates until scripts exit 0). Claude reviews afterward.
**Trigger:** Phase-1 reconnaissance (decision-13 §13); user "dive in full force" + chose the conditioning spike. Chain:
…`pathA_21c` (force + sign DONE) → **this (`pathA_C0` conditioning spike)** → promote constitutive family to a calibrated
branch + wire multi-knob calibrate-predict → `pathA_22`.

## Why this step
Everything downstream (the force normalization, `G`, brane `c_γ`, `m̂0²·S_port`, the real B2c verdict, and the g−2/5PN/
multi-defect surplus) needs NUMBERS from a converged deep throat. The solver can produce them once it can reach the
regime. This is the cheapest attack on reachability BEFORE committing to a production linear-solver rebuild; its diagnostic
decides whether the rebuild is even needed.

## Scope & stance
NUMERICAL CONDITIONING ONLY. The hard rule: **do NOT alter the faithful PDE operators** (the stationary gauged-GNLS +
localized-Maxwell residual physics in `coupled_branch.py`/`operators.py` — confirmed model-faithful in Phase-1) and do NOT
touch frozen physics or the `physical_export_permitted` guard. Every regularization MUST be a convergence aid that VANISHES
at (or is provably consistent with) the physical solution — it may change the PATH to the solution, never the solution.
Codex codes, Claude reviews. Standing infra: CPU sparse-direct, GPU OFF; `timeout 600` per script (exit 124 = failure →
reformulate, never raise the cap); read-only `research/pde_audit/simulation/`; YAML/markdown human output, JSON only for
machine artifacts; additive — do not break the existing chunk-1a/1b/1c gates.

## Work items

### C0-1 — stiff-core safeguard (the `√ρ→0` Jacobian degeneracy)
Add a density floor and/or a `log ρ` (or `√ρ`) change of variable in the matter block so the matter-Jacobian nonlinear
diagonal (`~ρ³`) does not vanish and degenerate the operator to a near-singular Laplacian as the core empties. The floor
must sit BELOW any physical core density (so it is never active at the converged solution), OR the change of variable must
be exact. PROVE solution-invariance: show a tame already-converged case is unchanged to tolerance with the safeguard on.

### C0-2 — regularize the `k1 ∝ r⁴/R0⁵` wall-coupling blow-up
The matter↔wall return coefficient `k1 = 4 V_radial r⁴/R0⁵` (`coupled_branch.py`) diverges as `R0→0` (a deep throat),
wildly mis-scaling the wall/`R0` Jacobian rows. Regularize (e.g. a floored `R0`, `1/(R0⁵+ε⁵)`, or a reformulation) such
that the regularization is negligible/exact at the physical `R0`. PROVE it does not shift the converged `R0(w)`.

### C0-3 — Jacobi (row/column) scaling of the bordered Jacobian
Scale the matter, Maxwell, wall/`R0`, and the dense mass/`μ` border lanes so they share a magnitude scale before the
GMRES preconditioner sees them. This is a preconditioner/linear-algebra change — it must NOT change the residual or the
solution, only the conditioning of the linear solve. Report the condition-number / GMRES-iteration improvement.

### C0-4 — depth-homotopy continuation
Add a short continuation in throat DEPTH / aspect ratio (mirroring the existing `continuation_K_values` loop), carrying
warm starts step to step, to crawl from the current shallow regime into a genuinely deep throat. Adaptive/backtracking
step control if a step fails. This is the existing continuation pattern applied to a new homotopy parameter — not a new
solver.

### C0-5 — THE DIAGNOSTIC (the deliverable that decides the path)
On a single moderately-deep branch, with C0-1..4 active, report as the throat deepens: `D0` (the modal stiffness gap),
the smallest singular value of the bordered Jacobian, and GMRES iteration growth per Newton step. State a clear verdict:
**SPIKE_SUFFICIENT** (scaling + homotopy reach the regime; option C proceeds on this solver) OR
**PRODUCTION_SOLVER_REQUIRED** (the floor is fundamental → the decision-01 escape hatch: a production linear solver /
multigrid / PETSc, possibly mesh grading near the wall). Quantify how deep we got vs the prior `τ≈0.029` floor.

## Acceptance criteria (PASS/FAIL; exit-0 NECESSARY not SUFFICIENT)
1. A throat genuinely DEEPER than the prior conditioning floor is REACHED + CONVERGED (residual below the established
   tolerance), OR — if not — an HONEST `PRODUCTION_SOLVER_REQUIRED` verdict with the diagnostic evidence (a stalled solve
   reported as stalled, NEVER masked as converged).
2. Each regularization (C0-1, C0-2) PROVEN solution-invariant: a tame converged case is unchanged to tolerance with it on
   (before/after, or MMS). C0-3 proven residual/solution-neutral.
3. The C0-5 diagnostic table produced (`D0`, σ_min, GMRES growth vs depth) + the explicit SPIKE_SUFFICIENT /
   PRODUCTION_SOLVER_REQUIRED verdict + the depth reached vs `τ≈0.029`.
4. The faithful GNLS/Maxwell operators are UNTOUCHED (diff shows only conditioning/continuation/linear-algebra changes);
   frozen physics + `physical_export_permitted` untouched; chunk-1a/1b/1c gates still pass.

**Fail conditions (explicit):** altering the GNLS/Maxwell operator physics; a regularization that shifts the converged
solution; masking non-convergence as convergence; touching frozen physics or the export guard; raising the `timeout 600`
cap. VALID outcomes: a deeper converged throat (SPIKE_SUFFICIENT), OR an honest PRODUCTION_SOLVER_REQUIRED verdict — both
are wins because both tell us the path.

## Out of scope
The production linear-solver rebuild itself (only C0-5 decides if it's needed); promoting the placeholder constitutive
family to a calibrated branch + the multi-knob calibrate-predict harness (the NEXT step, once reachability is in hand);
the scale-map → `m̂0²·S_port` → B2c rerun (`pathA_22`); any change to the model operators or frozen physics.

## Review (orchestrator, after Codex)
Transliteration-fidelity on the new conditioning/continuation code (one clean agent); an adversarial pass: is the deeper
convergence REAL or masked (check the residual history, not just the exit code)? do the regularizations genuinely vanish
at the solution (verify the solution-invariance evidence, don't trust the claim)? is the Jacobi scaling solution-neutral?
is the diagnostic honest and the verdict supported by the σ_min/GMRES data? a diff-check that the faithful operators are
untouched. Claude reads only residuals. Then gate to the calibrated-branch + calibrate-predict step.
