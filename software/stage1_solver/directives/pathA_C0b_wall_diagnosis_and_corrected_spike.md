# Directive pathA_C0b — Wall diagnosis + corrected conditioning spike (supersedes pathA_C0 execution)

**Status:** DRAFT (Claude-authored 2026-06-19; pending Codex design-review → confirm-pass → execution. User has GATED the
"full corrected spike" path). Supersedes the `pathA_C0` v2 EXECUTION (the directive `pathA_C0_conditioning_spike.md` v2 was
SOUND; its execution was deficient — see below). Step of **option C** (the throat-profile solve), task #78.

**Date:** 2026-06-19
**Owner:** Codex (codes + iterates until scripts exit 0; APPLIES the fixes to the existing module). Claude reviews after.
**Trigger:** The `pathA_C0` v2 run reported `PRODUCTION_SOLVER_REQUIRED`; a 2-agent adversarial review found the verdict
NOT EARNED and the diagnosis likely wrong. This directive specifies the corrected run.

## What the v2 run got RIGHT (KEEP — do not rebuild)
`src/stage1_solver/patha_c0_conditioning_spike.py` already has a CORRECT Single-Arbiter core: Newton merit / line-search /
convergence gated on the ORIGINAL unscaled `coupled_branch.patha_closed_branch_residual`; ε confined to preconditioner
assembly; Jacobi scaling applied as `R·J·C dz = R·(−F)`, step `= C·dz`, with the line search on the unscaled residual;
faithful operators untouched (verified by `git diff`). REUSE all of this verbatim.

## What the v2 EXECUTION got WRONG (the corrected run MUST fix each — these are now FAIL conditions)
1. **No genuine depth crawl.** `prefer_existing_b2c_background_predictor=True` cold-loaded a pre-existing B2c background for
   7/8 τ; every "converged" row was an OLD B2c solution re-verified at **0 Newton iterations** (residuals matched to all
   digits). The spike solved ZERO new throats and never warm-started a crawl below the floor.
2. **Single-ε, ~single-step below-floor attack.** The aid loop `break`s on first failure (`patha_c0_…py:~1052`), so the ε
   schedule `[0.08→0.02→0]` never advanced past index 0; no depth backtracking; Newton returned after one failed line
   search. The directive's "C0-1..4 ALL active, persistent failure" bar was never met.
3. **Tautological admissibility.** `residual_equality_max_abs` is HARDCODED `0.0` (`:~1200`); `epsilon_independence=PASS`
   rests on a hardcoded `final_aids_inactive=True` and only inspects already-converged rows — the dormant-case FAIL pattern.
4. **Mis-diagnosis.** The bordered Jacobian is near-singular (cond≈1e20) even at the CONVERGED shallow τ=0.03, and
   conditioning IMPROVES (cond≈1e14) at the failed deeper τ — so the wall is NOT depth-driven conditioning. Recommending
   "multigrid / mesh-grading" targets a failure mode the data contradicts.
5. **Binding-constraint miss.** `min_R0` INCREASES (0.748→0.806) as τ↓, the k1 clamp goes INACTIVE at the deepest τ, and
   `min_ρ≈7e-6` is flat — so τ↓ moves AWAY from the empty-core (√ρ→0, R0→0) regime the program cares about. τ-depth ≠
   R0-depth.
6. **σ_min untrustworthy.** A 1-norm LU *estimate* (true SVD gated out by state size 1297 > 360) that INVERTS between
   healthy and broken states.

## THE SINGLE ARBITER PRINCIPLE (unchanged — still governs)
The ORIGINAL, unmodified physical residual (`patha_closed_branch_residual`) is the SOLE arbiter of convergence AND
solution-invariance. Every aid is PATH-only: residual-equivalent OR `ε→0`-vanishing at the final solve; Newton merit /
line-search / convergence on the ORIGINAL `‖F‖`, never a scaled/floored surrogate. (The v2 core already satisfies this —
preserve it.)

## Scope & stance
DIAGNOSIS + a GENUINE conditioning/continuation attack. Do NOT alter the faithful PDE operators in
`coupled_branch.py`/`operators.py`; do NOT touch frozen physics (`a=1`, `L≈1.85`, boundary class, constitutive family) or
`physical_export_permitted`. Diff must show only conditioning/continuation/linear-algebra/diagnostic changes. Additive —
chunk-1a/1b/1c gates must still pass. CPU sparse-direct, GPU OFF; `timeout 600` per script (exit 124 = failure →
reformulate, never raise the cap); standalone `python3` (no `$RT exec`); YAML/markdown human output, JSON only for machine
artifacts; NO commentary `python3 -c` scripts. **Runtime:** v2 already took ~465s; the true-σ curves + backtracking may
exceed 600s in one process — so SPLIT the run into separate ≤600s scripts (e.g. the crawl, then the linear/σ/fold
diagnostics reading the crawl's saved states + cached assembled matrices). A diagnostic that times out emits `NOT_MEASURED`
for its outputs (→ `DIAGNOSTIC_INCOMPLETE` per C0b-6), NEVER a verdict resting on missing evidence.

## Work items

### C0b-1 — GENUINE warm-started depth crawl (fix #1, #2, #5)
Below the prior `τ≈0.029` floor, set `prefer_existing_b2c_background_predictor=False` (or bypass it): each deeper τ MUST
warm-start from the PREVIOUS C0-converged state and take a real continuation step — NO cold-loading a pre-existing B2c
background as "the attempt", NO τ skipping. Add ADAPTIVE BACKTRACKING: on a failed τ step, halve the τ increment and retry
(to a documented minimum step), so the crawl genuinely probes how far it can get. At EACH attempted τ, run the FULL ε
schedule `[0.08→0.02→0]` AND the Jacobi scaling — do NOT break the aid loop on first failure; record every (ε, Newton-iter,
line-search-α, residual) so "persistent failure with C0-1..4 all active" is either produced or disproven. **Each ε probe at
a τ MUST start from that τ-attempt's SAME clean initial state (the warm start), NOT from a poisoned partial state left by a
prior failed ε** — unless a prior ε at that τ converged, in which case proceed from it. The final accepted state at each
converged τ must satisfy the ORIGINAL residual with ε inactive (already enforced by the v2 core).
**Checkable schema (required):** emit JSON `tau_attempts[]`, one per attempted τ, each with `target_tau`, `delta_tau`,
`backtrack_index`, `init.source` (must be `previous_c0_converged_state` below the floor — NOT an existing-B2c file),
`used_existing_b2c` (must be `false` below the floor), `start_from_tau`, and `epsilon_attempts[]` (one row per ε with its
Newton iterations, line-search α's, GMRES curve, and final original residual). A reviewer must be able to confirm from this
schema alone that the crawl warm-started and ran the full ε schedule + backtracking at each failed τ.

### C0b-2 — TRUE σ_min + null-vector identification (fix #4, #6)
Replace the 1-norm-LU σ_min estimate with a TRUE smallest singular triplet of the assembled bordered Jacobian (reuse the
existing `assemble_closed_coupled_colored_sparse_jacobian`). **Allowed methods (with fallback order), record which was
used:** `dense_svd` (the bordered system is ~1297-dim → a dense `numpy.linalg.svd` of the densified matrix is sub-second
and is the ROBUST DEFAULT for a near-singular matrix), else `svds_propack_SM`, else `eigsh_normal_shift_invert`. Record
`sigma_method`, the smallest singular value `σ_min`, the triplet residual `‖J v − σ u‖` (a quality check), and `status`
(`MEASURED` / `NOT_MEASURED`). **If no method yields a triplet with an acceptable residual, emit `NOT_MEASURED` and BLOCK
any σ-evidence-dependent verdict** (see C0b-6). Compute it at a CONVERGED shallow τ (e.g. 0.03) AND the deepest attempted τ.
**State layout (use the ACTUAL closed layout `5*cells + nw + 1`):** the RIGHT singular vector `v_min` decomposes over state
columns `field[0:5n]`, `R0[5n:5n+nw]`, `μ[-1]`; the LEFT singular vector `u_min` decomposes over residual ROWS `PDE rows`,
`wall rows`, `mass-constraint row` (mass is a RESIDUAL ROW, not a state lane). Report the energy fraction of `v_min` over
the three state groups and of `u_min` over the three row groups. For the gauge check, report `gauge_projection_status`
(`measured_global_phase` if a global-phase generator is available and projected, else `not_available`) and, if measured,
the gauge fraction. This decides whether the near-singularity is a BORDER/constraint artifact, a GAUGE/zero mode, or
genuine field-block ill-conditioning. Also report cond of the FIELD BLOCK ALONE (the square `field[0:5n]` submatrix) vs the
FULL bordered system.

### C0b-3 — FOLD / turning-point test (fix #4) — ROBUST, not binary
A fold and a persistent near-null space can COEXIST: a persistent gauge/border near-null mode can MASK a fold that lives in
the NEXT singular direction. So do not rely on a single σ_min or a raw det-sign. Track, along the τ sequence, the SMALLEST
FEW singular TRIPLETS (e.g. the smallest 3–5) of the state Jacobian at fixed τ, WITH vector continuity (overlap
`|⟨v_i(τ_k), v_j(τ_{k+1})⟩|` to follow each mode across τ). Identify any PERSISTENT modes (border/gauge per C0b-2) and
assess the fold AFTER projecting them out or accounting for them separately — i.e., look for a tracked singular value that
DECREASES toward zero / a tracked mode whose associated eigenvalue changes sign near the floor, in the complement of the
persistent null space. Distinguish:
  - **No complement singular value crosses zero (the small σ's are ~constant-tiny and tied to persistent border/gauge
    modes)** ⇒ PERSISTENT near-null space (structural; C0b-2 gives its source) — NOT a fold.
  - **A tracked complement singular value decreases and crosses ~0 near τ≈0.029** ⇒ a FOLD / TURNING POINT — Newton cannot
    march straight through it regardless of the linear solver; the correct tool is pseudo-arclength (Keller) continuation,
    NOT a production linear solver.
`det(J)` sign (via `slogdet` / LU diagonal-sign) may be reported as SUPPORTING evidence ONLY if it is STABLE under the
stated row/col scaling and pivoting (a near-singular nonsymmetric 1297×1297 determinant sign is otherwise numerical noise).
State which case holds (or `DIAGNOSTIC_INCOMPLETE` if the triplets can't be tracked), with the tracked-σ(τ) curves as the
primary evidence.

### C0b-4 — τ-depth vs R0-depth decoupling (fix #5)
Report, across the crawl, `R0_min` and `min_ρ` as functions of τ, and state plainly whether decreasing τ approaches the
empty-core regime (`R0→0`, `√ρ→0`) the program actually targets. If τ↓ moves AWAY from small R0 (as v2 showed), FLAG that
τ is the WRONG depth knob and identify EVIDENCE-SUPPORTED CANDIDATE knob(s) that could drive R0→0 (e.g. brane tension, the
wall/geometry depth parameter, or a direct R0 continuation) — or report `UNKNOWN_NOT_MEASURED`. Do NOT claim a driver
without a supporting scan or derivation, and do NOT change frozen physics; this is a RECOMMENDATION for the next step, not
an in-directive re-parameterization.

### C0b-5 — genuine admissibility checks (fix #3)
Replace the hardcoded admissibility with MEASURED checks: (a) `residual_equality_max_abs` = the actual max-abs difference
`|F_original(state) − F_conditioned(state)|` at ε=0 evaluated on THREE real states — a tame, a near-floor, and the deepest
attempted — each with a stated tolerance and `PASS/FAIL/NOT_MEASURED`; (b) ε-independence = an actual
`epsilon_independence_table[]` (one row per ε in the schedule) on a NEAR-FLOOR converged state, each row carrying the final
ORIGINAL residual at that ε, the per-state-group deltas `{ψ, A, R0(w), μ}` vs the ε=0 solution, the observable deltas
`{D0 (or its inputs), B/Z/N overlaps}` vs ε=0, the comparison tolerance, and `PASS/FAIL/NOT_MEASURED`. It must be a
near-floor converged state (NOT a dormant shallow one where the aid is inactive). If a check cannot be computed, report it
`NOT_MEASURED`, never a hardcoded PASS.

### C0b-6 — THE VERDICT (expanded; gated on the corrected evidence with FALSIFIABLE numeric thresholds)
Emit exactly ONE verdict, each with its `verdict_support` block citing the specific numeric evidence (σ values, fractions,
ratios, tracked-σ(τ), residuals). **If any evidence a verdict needs is `NOT_MEASURED` (σ_min unmeasurable, triplets
untrackable, a diagnostic timed out), the verdict is `DIAGNOSTIC_INCOMPLETE` — do NOT force one of the four substantive
verdicts.** Definitions: the **field block** = the square `field[0:5n]` submatrix; `dominant_fraction(v)` = the largest of
the three group energy fractions of the relevant singular vector; `cond_ratio = cond(field_block)/cond(full_bordered)`.
  - **SPIKE_SUFFICIENT** — a throat genuinely DEEPER than `τ≈0.029` reached + CONVERGED on the ORIGINAL unscaled residual
    at the B2c tolerance, via the GENUINE crawl (`init.source = previous_c0_converged_state`, `used_existing_b2c = false`;
    NOT a re-verified pre-existing background).
  - **FOLD_TURNING_POINT** — C0b-3 shows a TRACKED complement singular value decreasing to ~0 / a tracked-mode eigenvalue
    sign change near the floor (after persistent border/gauge modes are accounted for) ⇒ pivot to pseudo-arclength
    continuation; a production linear solver is NOT the fix. (`verdict_support`: the tracked-σ(τ) crossing + the persistent
    modes it is distinct from.)
  - **NEAR_NULL_SPACE_STRUCTURAL** — C0b-2 shows a PERSISTENT near-null vector (σ_min ≈ constant-tiny across τ, NOT
    crossing) whose energy is dominantly localized to the border (`R0`/`μ`) state group or the mass-constraint row, or to a
    measured gauge mode (`gauge_projection_status = measured_global_phase`) — EXACT GATE: `dominant_fraction(v_min) ≥ 0.7`
    on a border/gauge group AND `cond_ratio < 0.1` (the field block is ≥10× better conditioned than the full system ⇒ the
    ill-conditioning is the border/gauge, not the physics) ⇒ targeted fix (border rescaling / Schur complement / gauge
    constraint / null-space deflation), NOT necessarily a full rebuild. `NEAR_NULL_SPACE_STRUCTURAL(gauge)` is FORBIDDEN
    unless the gauge fraction was actually measured.
  - **PRODUCTION_SOLVER_REQUIRED** — allowed ONLY if (a) the genuine crawl with the FULL aid stack + backtracking
    persistently fails, AND (b) a TRUE σ_min (status `MEASURED`) shows the near-singularity is in the FIELD BLOCK (its
    energy dominantly in `field[0:5n]` with `dominant_fraction ≥ 0.7`, and `cond_ratio ≥ 0.1` — the field block carries
    within an order of magnitude of the full system's ill-conditioning, so the field block is itself ill-conditioned), AND
    (c) it is NOT a fold (no tracked crossing) and NOT a border/gauge artifact. Specify the
    evidence-supported solver capability (e.g. bordered-aware Schur/deflation; generic multigrid ONLY if the field block is
    shown mesh-frequency ill-conditioned).

## Acceptance criteria (PASS/FAIL; exit-0 NECESSARY not SUFFICIENT)
1. The genuine crawl ran below the floor, PROVEN by the `tau_attempts[]` schema (C0b-1): `init.source =
   previous_c0_converged_state` and `used_existing_b2c = false` below the floor; the full ε schedule + adaptive
   backtracking ran at each attempted τ (multi-ε `epsilon_attempts[]`, multi-step — NOT a single-ε/0-iter probe of a
   cold-loaded background); each ε probe started from the clean τ-attempt state.
2. TRUE σ_min (`sigma_method` = dense_svd / propack / normal-shift-invert, NOT 1-norm-LU; with triplet residual + status)
   at a shallow converged τ AND the deepest τ; + the corrected decomposition — `v_min` over `field[0:5n]`/`R0[5n:5n+nw]`/
   `μ[-1]`, `u_min` over PDE/wall/mass-constraint rows; + `gauge_projection_status`; + field-block-vs-bordered cond.
3. The fold test produced as TRACKED smallest-few singular triplets across τ (with vector-continuity overlaps), with
   persistent border/gauge modes accounted for, and the explicit structural call (FOLD vs PERSISTENT null space vs
   DIAGNOSTIC_INCOMPLETE); det-sign only as scaling-stable supporting evidence.
4. R0_min(τ), min_ρ(τ) reported + the explicit τ-vs-R0 depth statement (does τ↓ approach R0→0? if not, candidate knob(s) or
   UNKNOWN_NOT_MEASURED — no unsupported driver claim).
5. Admissibility checks MEASURED (not hardcoded) on tame + near-floor + deepest states, with the `epsilon_independence_table[]`
   (per-ε residual + state-group/observable deltas vs ε=0 + tolerance + PASS/FAIL); uncomputable ones marked NOT_MEASURED.
6. Exactly one C0b-6 verdict (one of the four substantive verdicts OR `DIAGNOSTIC_INCOMPLETE` when required evidence is
   NOT_MEASURED), each with its `verdict_support` numeric block; faithful operators + frozen physics + export guard
   untouched (diff); chunk-1a/1b/1c gates still pass; report + machine JSON produced; tests updated (no tautological
   self-checks — the v2 hardcoded-PASS pattern is itself a FAIL).

**Fail conditions (explicit):** any v2 short-circuit recurring (cold-load as "attempt" below floor; break-on-first-ε;
hardcoded admissibility; 1-norm-LU σ_min as the verdict basis); altering faithful operators / frozen physics / export
guard; gating Newton on a scaled/floored surrogate; declaring PRODUCTION_SOLVER_REQUIRED without a TRUE σ_min ruling out
fold + border/gauge; masking non-convergence as convergence; raising the timeout cap; implementing a production solver here.

## Deliverables
Updates to `src/stage1_solver/patha_c0_conditioning_spike.py` (the corrected crawl + true σ_min + fold test + genuine
admissibility) + its tests (real checks); `reports/pathA_C0b_wall_diagnosis.md` (the crawl log, σ_min(τ) + v_min
decomposition, det-sign/fold call, R0-vs-τ, admissibility table, the single verdict + evidence + the recommended next
step); the machine JSON.

## Review (orchestrator, after Codex)
Re-run the 2-agent adversarial check: did the crawl GENUINELY run (warm-started from the previous C0-converged state, full ε
schedule, backtracking, no prefer-existing below floor — read the `tau_attempts[]` JSON: `init.source`, `used_existing_b2c`,
`epsilon_attempts[]`)? is σ_min a TRUE MEASURED value (`sigma_method` ∈ {dense_svd, svds_propack_SM,
eigsh_normal_shift_invert}, NOT 1-norm-LU, with triplet residual + `MEASURED` status — else `NOT_MEASURED`)? is the
fold/null-space call supported by the TRACKED smallest-few singular triplets across τ (vector continuity, persistent
border/gauge modes accounted for), with det-sign only as scaling-stable support? are the admissibility checks MEASURED (not
hardcoded; `epsilon_independence_table[]` on a near-floor state)? is the verdict gated on its falsifiable numeric
`verdict_support` (or honestly `DIAGNOSTIC_INCOMPLETE`)? transliteration-fidelity on the changed code; diff-check faithful operators untouched. Claude
reads only residuals/diagnostics. Then gate the next step (pseudo-arclength / targeted fix / production solver — per the
verdict).
