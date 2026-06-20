# Directive pathA_C0 — Conditioning spike: make the throat-profile solver REACH a deep/realistic throat

**Status:** READY v2 (Claude-authored 2026-06-19; design-review SOUND-WITH-FIXES → all 11 fixes + 6 new-problem items
applied; Codex confirm-pass = SOUND-AS-IS 2026-06-19; PENDING USER GATE → execution). First step of **option C** (the
throat-profile solve). Follows the Phase-1 solver
reconnaissance (decision-13 §0/§13): the PDE operators are MODEL-FAITHFUL and the solver + calibrate-predict harness are
~70–80% built + validated, but the solver has NEVER converged on a deep/realistic throat — the blocker is NUMERICAL
CONDITIONING (B2c stalls at `τ≈0.029` with "line search failed" = a diagnosed conditioning floor, NOT a physical edge).

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

## THE SINGLE ARBITER PRINCIPLE (load-bearing — read first; governs every work item)
**The ORIGINAL, UNMODIFIED physical residual is the SOLE arbiter of both convergence AND solution-invariance.** Every
item below (density floor, `k1` clamp, Jacobi scaling, variable change, depth/auxiliary homotopy) is a PATH aid only. An
aid is admissible iff it satisfies BOTH:
  (i) **Residual-equivalence OR vanishing-at-solution.** Either it is a provably residual-equivalent reparametrization
      (`F_original(T(y)) == F_conditioned(y)` to tolerance), OR it is parameterized by a regularization scale `ε` that is
      driven to its physical limit (`ε→0` / clamp inactive / floor inactive) at the FINAL solve so the accepted state
      satisfies the ORIGINAL physical residual.
  (ii) **Original-residual gating.** Newton's merit function, the line-search acceptance test, and the convergence gate
      are computed on the ORIGINAL residual `‖F‖` — never on a scaled `‖S·F‖` or a floored/clamped surrogate (a step can
      reduce `‖S·F‖` while WORSENING `‖F‖`).
**Checkable faithful-operator boundary:** every accepted final state MUST be evaluated by the unmodified physical residual
entry point in `coupled_branch.py` and pass below the established B2c tolerance with ALL aids switched off. For ANY
variable transform, PROVE `F_original(T(y)) == F_conditioned(y)` to tolerance on THREE states: a tame case, a near-floor
case, and the deepest accepted final state. This residual-equality test (a diff plus a numeric check) IS the allowed-vs-
forbidden line — it is not left to judgment.

## Scope & stance
NUMERICAL CONDITIONING ONLY. The hard rule: **do NOT alter the faithful PDE operators** (the stationary gauged-GNLS +
localized-Maxwell residual physics in `coupled_branch.py`/`operators.py` — confirmed model-faithful in Phase-1) and do NOT
touch frozen physics or the `physical_export_permitted` guard. The frozen branch pins `a=1`, `L≈1.85`, the boundary class,
and the Hooke/constitutive family — these are NEVER continuation knobs. Codex codes, Claude reviews. Standing infra: CPU
sparse-direct, GPU OFF; `timeout 600` per script (exit 124 = failure → reformulate, never raise the cap); read-only
`research/pde_audit/simulation/`; YAML/markdown human output, JSON only for machine artifacts; additive — do not break the
existing chunk-1a/1b/1c gates.

## Work items

### C0-1 — stiff-core safeguard (the `√ρ→0` Jacobian degeneracy)
The matter-Jacobian nonlinear diagonal (`~ρ³`) vanishes and degenerates the operator to a near-singular Laplacian as the
core empties. Add a conditioning aid to keep the matter block well-posed there. TWO hard constraints:
  - **A density floor below "physical density" is UNKNOWABLE a priori** — a deep/empty-core throat may have physical
    density genuinely near zero, so a one-shot "floor sits below any physical core density" claim is NOT acceptable. Use
    instead an `ε`-regularization (e.g. a floored diagonal `ρ³+ε³`, or a soft floor) that is driven `ε→0` along a homotopy
    so the FINAL accepted state satisfies the original unfloored residual; OR a final original-residual polish step.
  - **A `log ρ` / `√ρ` change of variable is NOT automatically an exact reparametrization of THIS matter block** — the
    field is a COMPLEX `ψ` (or `(ρ, phase)` / gauge-coupled current), not a bare real `ρ`. A density-only rewrite that
    ignores the phase/current/gauge lanes is FORBIDDEN. A variable change is admissible ONLY if it preserves the
    phase/current/gauge lanes exactly and is proven residual-equivalent per the Single Arbiter Principle.
PROVE solution-invariance per C0-1's invariance test below (NOT merely a tame shallow case).

### C0-2 — regularize the `k1 ∝ r⁴/R0⁵` wall-coupling blow-up
The matter↔wall return coefficient `k1 = 4 V_radial r⁴/R0⁵` (`coupled_branch.py`) diverges as `R0→0` (a deep throat),
wildly mis-scaling the wall/`R0` Jacobian rows. A `k1` clamp / floored `R0` (`1/(R0⁵+ε⁵)`) is acceptable ONLY as a
continuation/preconditioner aid, UNLESS the FINAL state also satisfies the ORIGINAL `k1 = 4 V_radial r⁴/R0⁵` residual
below tolerance with the clamp inactive. Report clamp ACTIVATION (where/when it bit) and an `ε`-independence table.
**Invariance must compare downstream observables, not just `R0`:** a shifted `R0(w)` changes the BdG/Maxwell overlaps and
hence `D0` and the verdict quantities — so prove invariance of `{ψ, A, R0(w), μ, D0, B/Z/N overlaps}`, not `R0` alone.

### C0-3 — Jacobi (row/column) scaling of the bordered Jacobian
Scale the matter, Maxwell, wall/`R0`, and the dense mass/`μ` border lanes so they share a magnitude scale before the
GMRES preconditioner sees them. This is admissible ONLY as solving an EQUIVALENT linear system with invertible row/column
scalings inside the linear solve — the Newton merit function, line-search acceptance, and convergence gate MUST use the
ORIGINAL unscaled residual `‖F‖` (per the Single Arbiter Principle). Report the condition-number / GMRES-iteration
improvement SEPARATELY for the unscaled vs scaled Jacobian. NOTE: a scaled `σ_min` is NOT physical evidence of
conditioning — always report BOTH the original-bordered-Jacobian `σ_min`/condition number and the scaled values.

### C0-4 — depth-homotopy continuation
Crawl from the current shallow regime into a genuinely deep throat by continuation, carrying warm starts step to step
(mirroring the existing `continuation_K_values` loop), with adaptive/backtracking step control if a step fails. The depth
continuation is in the **frozen-family `τ` path** (the physical depth parameter), PLUS optionally solver-only AUXILIARY
homotopy parameters that RETURN to their physical/frozen values at the final solve. **FORBIDDEN as continuation knobs:**
`a`, `L`, `r_mouth`, `r_exit`, `w_max`, the boundary class, and the constitutive family — varying any of these is a FREEZE
VIOLATION, not a continuation ("aspect ratio" is NOT a legal depth knob). Require a DEPTH SEQUENCE that spans the prior
floor: at least one target AT or BELOW `τ≈0.029`, OR a documented backtracked failed attempt below it — a single
moderately-deep branch is NOT enough to support a trend verdict.

### C0-5 — THE DIAGNOSTIC (the deliverable that decides the path)
Across the depth sequence (C0-4), with C0-1..3 active, report a full diagnostic table: `D0` (the modal stiffness gap);
the smallest singular value `σ_min` of the bordered Jacobian (BOTH unscaled and scaled); GMRES residual curves +
iteration growth per Newton step; the Newton residual HISTORY (`‖F‖` per step); line-search `α` values and step norms;
condition number before/after C0-3; `min ρ`, `min R0`, and floor/clamp ACTIVATION; and the depth metric (`τ`, `R0_min`).
State a clear verdict:
  - **SPIKE_SUFFICIENT** — a throat genuinely deeper than `τ≈0.029` is REACHED and CONVERGED **by the ORIGINAL, unscaled
    physical residual at the existing B2c background gate/tolerance** (no relaxed tolerance; convergence NEVER judged
    solely by a transformed/scaled/floored residual). Option C proceeds on this solver.
  - **PRODUCTION_SOLVER_REQUIRED** — allowed ONLY after documented attempts with C0-1, C0-2, C0-3, AND C0-4 ALL active
    (including backtracked depth steps driven at or below the prior `τ≈0.029` floor), where the residual history, `σ_min`,
    and GMRES data show PERSISTENT failure with those aids. This is the decision-01 escape hatch: it may RECOMMEND a
    production linear solver / multigrid / PETSc / mesh grading near the wall, but MUST NOT implement any of them in this
    directive. Quantify how deep we got vs the prior `τ≈0.029` floor.

## Acceptance criteria (PASS/FAIL; exit-0 NECESSARY not SUFFICIENT)
1. EITHER a throat genuinely DEEPER than the prior `τ≈0.029` floor is REACHED + CONVERGED **by the original unscaled
   physical residual below the established B2c tolerance** (SPIKE_SUFFICIENT), OR an HONEST PRODUCTION_SOLVER_REQUIRED
   verdict gated as in C0-5 (C0-1..4 all active, persistent-failure evidence) — a stalled solve reported as stalled,
   NEVER masked as converged, and NEVER via a relaxed tolerance or a scaled/floored surrogate residual.
2. Each conditioning aid PROVEN admissible per the Single Arbiter Principle: (a) the faithful residual-equality test
   `F_original(T(y)) == F_conditioned(y)` to tolerance on a tame, a near-floor, AND the deepest accepted final state;
   (b) for `ε`-regularizations (C0-1, C0-2), an `ε`-independence table showing `{ψ, A, R0(w), μ, D0, overlaps}` stable as
   `ε→0` (invariance demonstrated AT/NEAR the deep regime, NOT only on a tame shallow case where the aid is dormant);
   (c) C0-3 scaling proven residual/solution-neutral with Newton/line-search/convergence on the ORIGINAL residual.
3. The C0-5 diagnostic table produced in full (`D0`, unscaled+scaled `σ_min`, GMRES growth, Newton residual history,
   line-search `α`/step norms, condition before/after C0-3, `min ρ`/`min R0`, floor/clamp activation, depth metric) +
   the explicit SPIKE_SUFFICIENT / PRODUCTION_SOLVER_REQUIRED verdict + the depth reached vs `τ≈0.029`.
4. The faithful GNLS/Maxwell operators are UNTOUCHED (diff shows only conditioning/continuation/linear-algebra changes);
   frozen physics (`a`, `L`, boundary class, constitutive family) + `physical_export_permitted` untouched; chunk-1a/1b/1c
   gates still pass.

**Fail conditions (explicit):** altering the GNLS/Maxwell operator physics; varying any frozen quantity (`a`, `L`,
`r_mouth`, `r_exit`, `w_max`, boundary class, constitutive family) as a "continuation" knob; a regularization/variable
change that shifts the converged solution or whose invariance is shown only on a dormant shallow case; gating Newton
convergence/line-search on a scaled or floored surrogate residual; declaring SPIKE_SUFFICIENT via a relaxed tolerance;
declaring PRODUCTION_SOLVER_REQUIRED before C0-1..4 are all active with persistent-failure evidence; masking
non-convergence as convergence; touching frozen physics or the export guard; raising the `timeout 600` cap; implementing
a production solver/multigrid/PETSc in this directive. VALID outcomes: a deeper converged throat (SPIKE_SUFFICIENT), OR a
properly-gated PRODUCTION_SOLVER_REQUIRED verdict — both are wins, but a negative verdict is valid ONLY AFTER the C0-5
effort gates are met (it must never be a face-saving early exit).

## Out of scope
The production linear-solver rebuild itself (only C0-5 decides if it's needed, and may only RECOMMEND it); promoting the
placeholder constitutive family to a calibrated branch + the multi-knob calibrate-predict harness (the NEXT step, once
reachability is in hand); the scale-map → `m̂0²·S_port` → B2c rerun (`pathA_22`); any change to the model operators or
frozen physics.

## Review (orchestrator, after Codex)
Transliteration-fidelity on the new conditioning/continuation code (one clean agent); an adversarial pass: is the deeper
convergence REAL or masked — check the residual HISTORY against the ORIGINAL unscaled residual, not just the exit code or
a scaled/floored surrogate? do the regularizations genuinely vanish at the solution — verify the `ε`-independence table
and the `F_original==F_conditioned` equality on the DEEP state, don't trust the claim? is the Jacobi scaling
solution-neutral with Newton/line-search on the original residual? is the depth continuation `τ`-only with no frozen
quantity varied (diff-check `a`/`L`/`r_mouth`/`r_exit`/`w_max`/boundary/constitutive)? is the verdict supported by the
`σ_min`/GMRES/residual-history data, and is any PRODUCTION_SOLVER_REQUIRED verdict properly gated (C0-1..4 all active)?
a diff-check that the faithful operators + export guard are untouched. Claude reads only residuals. Then gate to the
calibrated-branch + calibrate-predict step.
