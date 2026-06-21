# Directive pathA_C0g-build — gauge-fix (path-only) → re-confirm the fold → gauge-fixed pseudo-arclength to round the throat-radius turning point

**Status:** **B-1/B-2/B-4 + deepcrawl DONE & COMMITTED (`0c3fc98f`); B-3 RE-SCOPED + READY — awaiting the user EXECUTION
gate (2026-06-21).** See the "⭐ B-3 EXECUTION RE-SCOPE (AMENDMENT, 2026-06-21)" section below + the staged execution prompt
`_scratch/pathA_C0g_B3_execution_prompt.md`. **B-3 amendment design-review loop COMPLETE:** Codex design-review =
SOUND-WITH-FIXES (8 MUST-FIX + 3 NICE) → all applied → confirm-pass = STILL-NEEDS (2 residuals: stale `:1539` anchor in
Acceptance B-3(ii); a B-3.1 cold-load loophole) → both FIXED → Codex **re-confirm = SOUND-AS-IS** (Claude, 2026-06-21).
Logs `_scratch/codex_c0g_B3_{design_review,confirmpass,reconfirm}.log`.
**(Prior BUILD-directive loop, 2026-06-20/21):** Codex design-review = SOUND-WITH-FIXES → all 8
fixes APPLIED → confirm-pass = STILL-NEEDS (3 residuals) → 3 residuals FIXED + the `min_R0` anchor self-verified
→ Codex **re-confirm = SOUND-AS-IS** (Claude, 2026-06-20/21). Logs `_scratch/codex_c0g_build_directive_{review,confirmpass,reconfirm}.log`. The 8 design-review fixes: (1) B-2 relative/trend gates, not the refuted absolute
`cond>1e10`; (2) B-2 numeric precedence DISSOLVED→CONFIRMED→INCONCLUSIVE with explicit criteria; (3) dropped the leading
"gauge-fix lets Newton converge nearer the fold" (closer sampling is ATTEMPTED, not required); (4) forbid hardening the
`(1/ξ)` penalty in the faithful operator + require an `r0`-mode-preserved check; (5) gauge-aligned PATH-ONLY proof (not a
raw norm); (6) gauge projection inside EVERY pseudo-arclength predictor/corrector solve; (7) pre-pick `min_R0` + margin as
the depth metric; (8) B-4 concrete match-gate tolerances + coverage. Core confirmed SOUND (freeze discipline, Single
Arbiter, B-3 gated on B-2, LM/PTC exclusion, anchors). Review log `_scratch/codex_c0g_build_directive_review.log`. This
draft completes the standard loop: design-review → fixes → confirm-pass → user execution gate. Follows
`pathA_C0g_diag_fold_vs_conditioning_battery.md` (the diagnostic), which returned machine
`MIXED/INCONCLUSIVE` ADJUDICATED **FOLD-LEANING** (simple fold in the throat-radius `r0` continuation at τ_fold≈0.0291233;
conditioning refuted; sonic not testable). Step of option C, tasks #78/#79.

**Date:** 2026-06-20
**Owner:** Codex (DESIGNS + codes + iterates until exit 0; chooses the gauge-fixing method and the pseudo-arclength
predictor/corrector; this directive states REQUIREMENTS + ACCEPTANCE, not the implementation route). Claude reviews after
(fidelity + adversarial clean agents).

**Trigger.** C0g-diag established (dual-reviewed HONEST/faithful/no-gaming): the Newton crawl stalls just past τ=0.029125
because the throat-radius (`r0`) continuation hits a **simple fold / turning point** at τ_fold≈0.0291233 — five agreeing
lines (premise f_nn≈1.0 on the `r0` near-null mode; cosθ 0.51–0.64 transversality; σ_min² linear R²=0.9994; bordered-cond
TREND cond(Jb) flat≈570 while cond(J·Q_perp)→2.3e5, ratio→404; conditioning affirmatively refuted). SEPARATELY, the raw
solver carries a real **511-dim gauge near-null subspace** (U(1) phase + A-sector ∇λ) that is *not* the stall driver
(Newton-step gauge fraction 5e-27) but wrecks the raw Jacobian conditioning and must be removed regardless. Both reviews
converge on the build below; both reject an LM/PTC conditioning remedy (it would damp but never round a fold).

## Scope (this step BUILDS the fix — it WILL change solver/continuation logic, NEVER the physics)
C0g-build adds gauge-fixing and pseudo-arclength continuation to the **solver/continuation layer**. It **MUST NOT** change
the physical residual `patha_closed_branch_residual` (`coupled_branch.py:512`), the faithful PDE operators, the frozen
physics, or `physical_export_permitted`. **Single Arbiter Principle (binding):** the ORIGINAL, unmodified
`patha_closed_branch_residual` remains the SOLE arbiter of convergence and solution-identity; the gauge-fix is **PATH-ONLY**
(a constraint/projection acting in solver/preconditioner coordinates that is residual-equivalent at the solution — the
accepted physical state is identical, to tolerance, to one obtained without it); the pseudo-arclength constraint borders
the system for continuation but the physical block of the residual it drives to zero is the unchanged original residual.

**Work is GATED into four items; B-3 runs only if the B-2 re-confirm gate passes.**

### B-1 — Gauge-fix (PATH-ONLY), unconditional
Remove the 511-dim gauge near-null subspace from the solve: the U(1) global-phase exact-null mode (C0c-confirmed) and the
A-sector ∇λ gradient modes. **Codex designs the method** (candidates from the synthesis: a phase-pin / reference-point
constraint for U(1); tree-cotree or a Coulomb-gauge (`div A = 0`) constraint for the A-sector — Codex chooses + justifies).
It MUST act in solver/preconditioner coordinates or as bordering constraints, **NEVER inside the frozen residual**. **Do
NOT change `ξ` or harden the soft `(1/ξ)`grad-div penalty in the operator** — that penalty lives in the FAITHFUL operator
(`operators.py:429`) and editing it violates the operator/residual freeze; any grad-div-style regularization must be a
solver/preconditioner-coordinate aid only. **REQUIRED physical-mode-preservation check:** verify the `r0` fold near-null
mode is NOT removed by the gauge-fix — it must remain in the gauge-COMPLEMENT (use the existing gauge-complement machinery,
`patha_c0_conditioning_spike.py:4417`); a gauge-fix that lifts or absorbs the `r0` mode is a FAIL (it would mask the fold).
- **Acceptance B-1:** (i) the gauge-fixed Jacobian's gauge near-null subspace is gone — σ_min attributable to the gauge
  sector is lifted by orders (report the gauge-sector σ_min before/after); (ii) **PATH-ONLY proof (gauge-aligned, not a raw norm):** take a state solved WITH the gauge-fix and one
  WITHOUT it (at a τ where both converge); FIRST align them on the gauge orbit (explicit global-phase alignment for ψ + a
  gauge-orbit least-squares alignment for A), THEN require agreement on the GAUGE-INVARIANT quantities (|ψ|/ρ, the
  gauge-invariant field strengths) and on the `r0`, `μ` lanes to ≤ solve tol, with the ORIGINAL residual equal at both —
  i.e. the gauge-fix moved ONLY along gauge orbits, not physics (a raw `‖Δ‖` is NOT acceptable, since gauge-equivalent
  representatives legitimately differ);
  (iii) operators/frozen/export untouched (git diff).

### B-2 — Re-confirm the fold on the gauge-fixed crawl (the GATE to B-3)
Re-run the CHEAP C0g-diag Steps 1–3 + 6 (reuse the `_c0g_*` machinery) on the gauge-fixed crawl. ATTEMPT adaptive sampling
**closer to τ_fold** than the current last-convergeable τ=0.029125; if the crawl cannot get closer, do **NOT** treat that
as a failure — classify by the predefined trend/failure evidence below. (Gauge is NOT the stall driver, so closer
convergence is not guaranteed and MUST NOT be a hard requirement that falsely blocks B-3.) **All gates are
RELATIVE/trend-based — do NOT reinstate the self-defeating absolute `cond>1e10` bar the diagnostic refuted** (see the
report's "Post-hoc adjudication"). Decide, with numeric support, exactly one, evaluated in THIS precedence order:
- **FOLD_DISSOLVED** (checked FIRST) — with the gauge fixed, plain Newton produces **≥2 branch-continuous converged states
  BELOW τ_fold by a stated margin** (each ‖original residual‖ ≤ B2c tol) with NO `r0`-turning signature (σ_min not
  collapsing, cosθ not approaching transversality). ⇒ the "fold" was entangled with the gauge; STOP, do NOT build
  pseudo-arclength, report + re-gate with the user (a cheaper win than expected).
  **FOLD_DISSOLVED MUST be honestly contestable (no sandbagging — the C0f lesson):** the below-fold attempts MUST use the
  FULL Newton budget + τ-backtracking and be routed through the existing `crawl_persistent_failure` guard
  (`max_newton_iters_override=None`, `max_tau_backtracks>0`, `attempted_backtracking=True`). A 1-iteration / 0-backtrack
  below-fold attempt is a CRIPPLED-CONFIG artifact, NOT fold evidence. The dissolve branch must be genuinely REACHABLE in
  code (no source-string / precedence wiring that precludes it). Diagnostic: a τ that is ABOVE τ_fold and fails ONLY under
  a crippled budget but converges under the full budget REFUTES the dissolve-side "failure" — in that case FOLD_CONFIRMED
  is NOT supported on the dissolve side and the gate must re-sample (full budget) until it either crosses τ_fold
  (FOLD_DISSOLVED) or genuinely stalls below it under full budget (FOLD_CONFIRMED).
- **FOLD_CONFIRMED** — only if NO clean below-fold crawl exists (FOLD_DISSOLVED fails) AND the RELATIVE fold-trend gates
  pass: Step-2 cosθ stable and above the fold threshold; σ_min²(τ) linear-monotone with good fit (R² ≥ 0.95); cond(Jb)
  FLAT within a stated band across the sampled τ; and cond(J·Q_perp)/cond(Jb) (or σ_min(Jb)/σ_min(J·Q_perp)) INCREASING by
  at least a stated factor (e.g. ≥10×) toward τ_fold. ⇒ proceed to B-3.
- **STILL_INCONCLUSIVE** — neither set of criteria is met (trend flat/ambiguous, or fit poor); report + re-gate (do NOT
  force B-3).
**B-2 is a HARD GATE:** B-3 is implemented ONLY on FOLD_CONFIRMED.

### B-3 — Gauge-fixed pseudo-arclength continuation (only if B-2 = FOLD_CONFIRMED)
Implement Keller pseudo-arclength continuation on the gauge-fixed solve to round the τ_fold turning point and continue to
deeper throats: the continuation parameter becomes arclength `s`; τ becomes a solved unknown (free to decrease through the
fold and increase again); the bordered system = [original physical residual ; arclength constraint], **with the B-1 gauge
projection/constraints applied INSIDE every pseudo-arclength predictor AND corrector solve** (else the 511-dim gauge null
space returns and the bordered solve degenerates). The physical block driven to zero remains the ORIGINAL residual
(`coupled_branch.py:512`). **Codex designs** the predictor (tangent), the corrector (bordered Newton), step-size
adaptation, and the turning-point handling.
- **Acceptance B-3:** (i) the continuation **rounds τ_fold** — it produces converged states on BOTH sides of the fold
  (τ decreasing then increasing along `s`), each at ‖original residual‖ ≤ the B2c tolerance (Single Arbiter; NOT the
  bordered/least-squares norm); (ii) it reaches a throat **genuinely deeper** than the pre-fold branch, measured by the PRE-CHOSEN primary
  depth metric **`min_R0`** (recorded by the crawl as `min_R0_during_final` — `~:1875`/`~:8118`/`~:8219`; **locate by symbol
  with `rg`**) — report `min_R0` past the fold. **(NOTE: SUPERSEDED by the amendment's depth-metric correction — min_R0 is
  non-monotone/nudges UP, so "deeper" = SOLVED-τ continuation progress with min_R0 > 0, NOT a min_R0 decrease.)**
  Secondary metrics min ρ / throat depth reported for context; (iii) it
  traces ONE physical branch (no branch-jump; no frozen-param change — only τ and the state move along the natural
  continuation); (iv) operators/frozen/export untouched.

### B-4 — Analytic/sparse Jacobian assembly (unconditional, performance/conditioning)
After a color audit, provide analytic/sparse Jacobian assembly to replace (or augment) the colored-JVP dense path,
handling the dense μ/mass/wall rows (253 colors = deterministic radius-3 coloring, NOT a bug). Residual-equivalent with a CONCRETE
match gate: the assembled J must match the JVP of the ORIGINAL residual to concrete tolerances FIXED BEFORE execution
(default: relative ≤ 1e-8 AND absolute ≤ 1e-10 per block; tighten if the float64 JVP supports it), with
coverage = multiple random `Jv` probes PLUS the dense wall/mass/μ rows checked explicitly (the rows the colored path
handles specially); report the per-block max abs/rel mismatch. This is a perf/conditioning aid, NOT a physics change.

---

## ⭐ B-3 EXECUTION RE-SCOPE (AMENDMENT, 2026-06-21) — target the DEEPER near-singularity; honest 3-outcome adjudication

**Why this amendment exists.** The B-2 re-confirm did NOT return FOLD_CONFIRMED. With the B-1 gauge-fix + FULL Newton
budget + B-4 fast assembly, the crawl **CROSSED** the putative shallow fold τ_fold≈0.0291233 — 8 admissible
branch-continuous states down to **τ=0.0291132** (each ‖original residual‖ ≤ 1e-6, positive density, min_R0≈0.79). So the
machine result was **FOLD_DISSOLVED (shallow) → STILL_GOING (deeper)** = the directive's **STILL_INCONCLUSIVE** branch,
whose instruction is "report + re-gate with the user (do NOT force B-3)." That report was made; **the user RE-GATED B-3 on
2026-06-21** as the disambiguating tool for a NEW, genuinely-real near-singularity found ~1e-5 deeper. **B-3 is therefore
earned NOT as "round a confirmed fold" but as the ONLY tool that can distinguish a genuine fold from a hard near-singular
wall** at a state where plain Newton provably cannot proceed. The literal "B-3 only on FOLD_CONFIRMED" gate (B-2 §) is
superseded for this run by the user re-gate; every other constraint in this directive (Single Arbiter, gauge-projection-
inside-every-solve, freeze discipline, B-4 residual-equivalence, no LM/PTC, timeout policy) remains BINDING.

**The anti-sandbag contest is already satisfied (NOT a crippled-config artifact).** The deeper stall was produced under the
FULL machinery, recorded in `runs/pathA_C0g_deepcrawl/pathA_C0g_deepcrawl.json` (`persistent_failure_guard`):
`full_newton_budget=True`, `max_tau_backtracks=6`, `attempted_backtracking=True`, 21 failed attempts, residual stuck
~1.04e-6 just below τ=0.0291132. This is the contest-wall discipline DISCHARGED — the C0f/C0f2 trap does not apply here.

**The new target (re-confirm at execution from the deepcrawl JSON, do not hardcode).**
- Deepest CONVERGED state: **τ=0.0291132**, min_R0≈0.7979, min_ρ≈7.24e-6, μ≈2.0345, original-residual Linf ≤ 1e-6, admissible PASS.
- Near-singularity signature there: **σ_min(J·Q_perp) collapses 3.95e-4 → 7.03e-6 (factor 56)**; σ_min²(τ) fit
  `LINEAR_MONOTONE_FOLD_SUPPORT` r²=0.957, zero-crossing **τ≈0.0291139**; **branch_tangent = NO_TANGENT_REVERSAL** (this
  no-reversal-yet is exactly the ambiguity B-3 resolves). Below τ=0.0291132 the full-budget crawl persistently MISSES tol.
- **Warm-start** the pseudo-arclength continuation from the deepcrawl converged states under
  `runs/pathA_C0g_deepcrawl/gauge_fixed/gauge_fixed_crawl/states/` (deepest converged = `attempt_025_tau_0p0291132.npz`),
  GENUINE continuation, **NO cold-load** (see the no-cold-load rule under B-3.0).

**Depth-metric correction (binding for the B-3 acceptance below).** The original B-3 acceptance (ii) said "deeper = min_R0
DECREASING below the τ=0.029125 state." That is REFUTED by the data: **min_R0 is non-monotone and nudges UP with depth**
(0.7938→0.7986 over the converged ladder; `approached_small_core_regime: False`). The throat STAYS OPEN and slightly
widens — which SATISFIES the user's hard model constraint that the throat radius can never reach 0. Therefore:
- **"Deeper / progress" is measured by CONTINUATION PROGRESS: converged states (‖ORIGINAL residual‖ ≤ 1e-6) traced BEYOND
  the τ=0.0291132 stall that plain Newton cannot reach** — NOT by min_R0 decreasing.
- **min_R0 must be REPORTED at every B-3 state and must stay > 0** (a clear margin). A min_R0 trending toward 0 is a
  model-FORBIDDEN pathology to FLAG and HALT on, NOT a success signal.
- The PRE-FIXED progress margin (set BEFORE execution, not post-hoc) is on **SOLVED-τ progress past the stall** — use the
  corrector-SOLVED τ, NEVER the nominal/target τ. All counted states must be DISTINCT accepted continuation states (each
  from a genuine predictor→corrector step with recorded provenance; pairwise state separation above a small declared floor —
  no microscopic below-threshold duplicates). Then, by outcome:
  - **CONTINUES_NO_TURNING** requires **≥ 3 distinct accepted states with solved τ strictly below 0.0291132**, the deepest at
    solved τ ≤ 0.0291131 (each original-residual ≤ 1e-6, min_R0 > 0), one branch, no frozen-param change.
  - **ROUNDS_FOLD** does NOT require solved τ < 0.0291132 — at a genuine fold τ stops decreasing at ≈τ_fold≈0.0291139 and
    rises again. It requires reaching the turning point (minimum solved τ WITH a `dτ/ds` sign reversal) PLUS **≥ 3 distinct
    accepted converged states on the FAR side of the turning point** (arclength `s` beyond the reversal), each
    original-residual ≤ 1e-6, min_R0 > 0, branch-continuous, no frozen-param change.

**Staged execution (mirror the B-1/B-2/B-4 cadence — each stage reports, dual-reviewed before the next):**
- **B-3.0 — VALIDATE the machinery on the known-good branch FIRST (anti-bug / anti-tautology gate).** Implement the
  gauge-fixed Keller pseudo-arclength predictor (tangent) + corrector (bordered Newton), B-4 fast assembly for the
  direction, B-1 gauge projection INSIDE every predictor AND corrector solve. Run it on the ALREADY-CONVERGED region (e.g.
  from τ≈0.02914 down to τ=0.0291132) and require it to RE-TRACE THE SAME PHYSICAL STATES the plain Newton crawl already
  found — converged on the ORIGINAL residual ≤ 1e-6, gauge-invariant agreement with the recorded states to solve tol, no
  branch jump. **If B-3.0 does not reproduce the existing branch, STOP — the implementation is wrong; do not attempt the
  hard region.** This stage proves the continuation is correct before it is trusted to adjudicate the singularity.
  - **Method spec (MUST).** Reduced unknowns = (`Q_perp` coefficients, τ); the tangent is solved in the gauge complement and
    normalized in a DECLARED scaled metric; orientation by sign-continuity (positive dot with the previous accepted tangent;
    a secant fallback for the FIRST tangent); the arclength constraint uses the accepted previous state + tangent; the
    bordered-Newton corrector drives the UNCHANGED original residual (physical block) to ≤ 1e-6; adaptive `ds` with
    accept/reject as `σ_min(J·Q_perp)→0`; record the corrected `dτ/ds` sign at every accepted state.
  - **No cold-load (anti-tautology).** ONLY the single initial SEED state may be loaded to start the continuation; the
    recorded deepcrawl states may be read ONLY AFTER the run, for comparison — NEVER to initialize a predictor/corrector.
    Emit per accepted state: predecessor state id, predictor step, corrector iters, solved τ, and an explicit
    `not_initialized_from_comparison_artifact: true`.
  - Note: `_c0g_bordered_conditioning_for_state` (~`:5954`) is SVD-only diagnostic evidence (it builds a bordered matrix and
    takes its SVD — no predictor/corrector loop). B-3 must build a GENUINE continuation predictor/corrector, not reuse it for stepping.
- **B-3.1 — ATTEMPT the near-singularity at τ≈0.0291132.** Continue along arclength `s` through the σ_min collapse. Per
  step record: solved τ (now an unknown), the predictor sign `dτ/ds` (tangent direction), ‖original residual‖ Linf,
  min_R0, min_ρ, μ, σ_min(J·Q_perp), bordered-system conditioning, step size, wall seconds.
- **B-3.2 — ADJUDICATE honestly. Exactly ONE SCIENTIFIC verdict of three — but ONLY after measurement COMPLETES.** A
  per-script timeout, an interrupted/unfinished chunk, or too-few accepted states yields **NOT_MEASURED / HALT** (chunk +
  resume), NEVER one of the three scientific verdicts. NONE of the three is a can't-fail pass (the JOB is to disambiguate):
  - **ROUNDS_FOLD** — the tangent `dτ/ds` genuinely REVERSES sign (τ decreases then increases — a real turning point),
    yielding converged states (original-residual ≤ 1e-6) on BOTH sides of the turning point, satisfying the τ-progress
    margin above, min_R0 > 0 throughout, one branch, no frozen-param change. ⇒ genuine fold, ROUNDED; the throat solve
    continues past it. (The original B-3 acceptance (i)/(iii)/(iv) apply; acceptance (ii) is replaced by the τ-progress
    metric above.)
  - **CONTINUES_NO_TURNING** — pseudo-arclength continues SMOOTHLY past τ=0.0291132 with NO tangent reversal (just powered
    through a conditioning dip plain Newton couldn't, exactly as the shallow "fold" turned out), reaching the τ-progress
    margin. ⇒ the deeper near-singularity was ALSO conditioning, not a fold; report how deep it reached + re-gate whether
    to keep crawling.
  - **GENUINE_ENDPOINT / HARD_WALL** — declared ONLY after ALL of: B-3.0 PASS; B-4 reconfirm PASS; AND the FULL
    pseudo-arclength machinery (adaptive `ds` reduction + tangent reseeding + secant fallback, within the 600s-chunk policy)
    STILL cannot drive the corrected ORIGINAL residual ≤ 1e-6 — i.e. the bordered system itself goes singular, or the
    corrector genuinely stalls with τ free along arclength. ⇒ NOT a simple fold — a genuine branch endpoint or
    higher-codimension singularity; report it as a NUMERICAL/branch feature, noting min_R0≈0.798 ≫ 0 (BEFORE the physical
    empty-core/opening, so it is a solver/branch feature, NOT the model's throat closing). It is a VALID outcome ONLY when
    those exhaustion criteria are met; ABSENT them, classify as an implementation failure or NOT_MEASURED — endpoint must
    NOT become a face-saving exit.
- **B-3 acceptance is "the adjudication is honest and earned," NOT "ROUNDS_FOLD."** Claiming ROUNDS_FOLD requires an
  ACTUAL `dτ/ds` sign reversal + original-residual convergence on BOTH sides — never the bordered/arclength/least-squares
  norm, never "got a bit deeper" misread as rounding. The three outcomes are mutually exclusive and all reportable.

**B-3 anti-gaming checklist (the dual-review will enforce):** (a) B-3.0 reproduces the EXISTING branch before any hard-region
claim; (b) convergence is judged ONLY on the original `patha_closed_branch_residual` ≤ 1e-6 — the bordered/arclength norm is
NEVER the arbiter; (c) ROUNDS_FOLD shows a real tangent reversal, not a fabricated turning point or a branch jump
(state continuity checked); (d) the τ-progress margin and min_R0>0 are PRE-FIXED, checked, not asserted; (e) B-4 assembly
re-confirmed residual-equivalent before reuse; (f) gauge projection is inside EVERY predictor + corrector solve (else the
511-dim gauge null space returns and the bordered solve degenerates); (g) genuine warm-start, NO cold-load; (h) a per-script
timeout ⇒ NOT_MEASURED (chunk + resume), NEVER raise the 600s cap — if a single bordered Newton/linear solve cannot fit
600s even after chunking + B-4 + the modal/port form, HALT and bring it to the user.

---

## ⭐⭐ B-3 FOLLOW-UP (AMENDMENT, 2026-06-21) — decisive A-vs-B diagnostic → conditional gate-refine + fallback-wire → re-run

**Why.** The first B-3 run (commit pending; report `reports/pathA_C0g_B3_pseudoarclength.md`) BUILT a faithful gauge-fixed
Keller pseudo-arclength continuation (fidelity audit FAITHFUL_WITH_NOTES; Single Arbiter intact; gauge `Q_perp` inside
every solve; B-4 direction-only; no cold-load; freeze intact) and stopped HONESTLY at the B-3.0 validation gate
(`reproduction = FAIL` → `NOT_MEASURED`; B-3.1 not attempted). The dual review converged: the continuation reproduces the
recorded branch to ~1e-6 wherever the recorded states are reliable (3/4 comparison points PASS), and the lone FAIL is
**~85% a GATE-CALIBRATION ARTIFACT, not a continuation bug** — at the deepest point (τ=0.02911625, σ_min≈7.9e-5) it is
compared against an UNDER-CONVERGED recorded plain-Newton state (20 iters = max, residual 2.2e-7, step_norm still 1.08e-6),
while the continuation's OWN residual there is ~1e-10 (≈1000× tighter); the 6.82e-4 mismatch lives almost entirely in the
soft `r0`/`μ` fold (near-null) mode (ρ, curl A agree to 1e-6/1e-13). **The ONE missing decisive check (not run under the
no-commentary-script rule): project that mismatch onto the σ_min near-null subspace.** Separately: the continuation stalled
just ABOVE the recorded deepest τ and the directive-required fallback machinery (adaptive `ds` reduction / tangent reseed /
secant fallback) was NEVER exercised (`*_attempted: false`) — B-3.1 needs it wired in regardless.

**GUARD AGAINST GOALPOST-MOVING (binding):** refining a validation gate AFTER it failed can mask a real bug. Therefore the
DECISIVE diagnostic (Phase D) GATES the gate-refinement, with a PRE-COMMITTED, FALSIFIABLE threshold stated here BEFORE
execution. If Phase D refutes the artifact hypothesis, the gate is NOT touched and the continuation is debugged instead.

### Phase D — DECISIVE A-vs-B diagnostic (run FIRST; read-only on existing states + a fresh Jacobian)
On the existing artifacts ONLY (no new solve): the B-3.0 FAIL-point continuation state (`B3_0_accepted_009`, τ≈0.02911625,
under `runs/pathA_C0g_B3_pseudoarclength/B3_0/`) and the recorded deepcrawl state at the SAME τ
(`runs/pathA_C0g_deepcrawl/.../states/attempt_008_tau_0p02911625.npz`). **Codex designs** the projection; requirements:
gauge-ALIGN the two states first (global-phase for ψ + gauge-orbit LS for A — reuse the B-1 path-only alignment), form the
physical difference `Δ = x_cont − x_rec`, restrict to the gauge complement `Q_perp` (project out residual gauge), normalize;
assemble the scaled `J·Q_perp` near-null right-singular vectors at the FAIL point (top-k, k≈5, smallest σ — reuse the
existing dense-SVD machinery), and report **(i)** the overlap fraction `‖P_nullspace Δ̂‖` of the unit mismatch with the
top-k near-null subspace, **(ii)** the single-mode overlap `|⟨Δ̂, v_min⟩|`, **(iii)** corroborating facts: the recorded
state's own original-residual / Newton-iters / step_norm (is it under-converged?) and the continuation state's
original-residual; the stiff-lane (ρ, curl A) agreement.
- **PRE-COMMITTED verdict:** `ARTIFACT_CONFIRMED` iff top-k near-null overlap ≥ **0.90** AND the recorded state is
  demonstrably under-converged (residual > 10× the continuation's, or iters at max with step_norm ≳ 1e-6) AND the stiff
  lanes agree to tol. `ARTIFACT_REFUTED` iff overlap < 0.90 (the mismatch is NOT along the soft mode) ⇒ STOP, do NOT refine
  the gate, report a likely continuation defect to debug. `DIAGNOSTIC_INCOMPLETE` otherwise (report + re-gate).

### Phase G — Conditional B-3.0 gate refinement (ONLY if Phase D = ARTIFACT_CONFIRMED) — criteria PRE-STATED here
Replace the equality-to-recorded-state test (which treats unreliable references as ground truth) with these PRINCIPLED
criteria, FIXED BEFORE the re-run (not tuned to pass): B-3.0 PASSES iff ALL hold across the validation crawl —
(a) **reference-reliability filter:** the gauge-invariant 5e-5 equality test (incl. the soft `r0`/`μ` lanes) is applied ONLY
at recorded comparison states that are themselves WELL-CONVERGED (original-residual ≤ a stated factor ≪ 1e-6, e.g. ≤ 1e-8,
AND Newton iters not at max AND step_norm small) — i.e. trustworthy references; (b) at comparison states that FAIL the
reliability filter (under-converged near the singularity), require only: the continuation's own original-residual ≤
min(recorded residual, 1e-6), branch-continuity (no L2 jump), AND stiff-lane (ρ, curl A) agreement to tol — the soft
`r0`/`μ` equality is DROPPED there because the reference is provably unreliable along that mode (established by Phase D);
(c) branch-continuity and continuation-residual ≤ 1e-6 hold at EVERY accepted state. **The σ_min floor / residual factor /
iters-max / step_norm thresholds are stated in the execution prompt BEFORE running.** Single Arbiter unchanged.

### Phase F — Wire in the full fallback machinery (unconditional for the re-run)
Make the corrector/continuation actually EXERCISE adaptive `ds` reduction + tangent reseed + secant fallback before any
TARGET_MISSED/stall is recorded (the directive already requires these for any GENUINE_ENDPOINT call; the first run never
fired them). A stall is only honest after these are genuinely attempted (record `*_attempted: true` with counts).

### Phase R — Re-run B-3.0 → B-3.1 → B-3.2 (the actual disambiguation)
With the refined gate (Phase G) + full machinery (Phase F): B-3.0 must PASS honestly (no cold-load; provenance intact);
then B-3.1 attempts the hard region; then B-3.2 returns ONE honest verdict (ROUNDS_FOLD / CONTINUES_NO_TURNING /
GENUINE_ENDPOINT / NOT_MEASURED) per the amendment's solved-τ progress + exhaustion criteria. All prior Single-Arbiter /
gauge / freeze / timeout / no-cold-load constraints REMAIN BINDING.

**Acceptance (follow-up):** Phase D returns a pre-committed verdict with the numeric overlap; the gate is refined ONLY on
ARTIFACT_CONFIRMED with the pre-stated criteria; Phase F fallbacks genuinely fire; Phase R returns an honest B-3.2 verdict
on the ORIGINAL residual. **Fail conditions (follow-up):** refining the gate without ARTIFACT_CONFIRMED; choosing the
σ_min/residual thresholds post-hoc to pass; the projection computed without gauge-alignment; declaring ARTIFACT_CONFIRMED
on overlap < 0.90; any prior Fail condition.

---

## Key facts (anchors — Codex re-confirm at execution)
**⚠️ Line numbers drift — LOCATE BY SYMBOL with `rg`, do not trust exact lines.** Current (2026-06-21) spot-checks:
`BACKGROUND_RESIDUAL_TOL` at `patha_c0_conditioning_spike.py:71` (was `:46`); `min_R0_during_final` recorded at `~:1875`/
`~:8118`/`~:8219` (was `:1539`). Re-confirm every anchor below by symbol at execution.
- Original residual / arbiter: `patha_closed_branch_residual` (`coupled_branch.py:512`); state pack
  `[psi_real,psi_imag,a0,ar,aw,r0,mu]` (`:137`, unpack `:166-174`).
- Gauge generators to reuse for B-1 / B-2 (analytic, NOT from SVD modes): `_c0c_generators_for_state` (phase),
  `_c0d_scalar_gradient_matrix`/`_c0d_build_gauge_subspace` (A-sector ∇λ), `_c0e_coupled_gauge_matrix`, plus the C0g
  gauge-complement wrapper (`_c0g_*`) — all in `patha_c0_conditioning_spike.py`.
- Explicit sparse J: `assemble_closed_coupled_colored_sparse_jacobian` (`preconditioners.py:476`); grad closures
  `operators.py:254`/`:270`; curl `F_rw` `:431` (`:429`=div). Soft gauge penalty `(1/ξ)grad(div A)`, ξ=1.
- C0g-diag verdict + numeric support: `reports/pathA_C0g_diag_fold_vs_conditioning.md`; converged + stalled states under
  `runs/pathA_C0f2_timing_rerun/states/`. τ_fold≈0.0291233 (Step-3 zero-crossing).
- Continuation crawl + Single-Arbiter + warm-start machinery: `patha_c0_conditioning_spike.py` (genuine warm-start,
  `prefer_existing_b2c_background_predictor=False`, NO cold-load).

## Acceptance criteria (PASS/FAIL; exit-0 NECESSARY not SUFFICIENT)
**⭐ 2026-06-21 B-3 RE-SCOPE OVERRIDE (governs #2/#3 + the Fail conditions for the B-3 run):** B-2 returned
`FOLD_DISSOLVED (shallow) → STILL_GOING (deeper)`, NOT `FOLD_CONFIRMED`. Therefore the criterion "B-3 implemented ONLY on
FOLD_CONFIRMED" (and the matching Fail condition) is SUPERSEDED for this run by the user re-gate (see the amendment). For
the B-3 run, acceptance is the **honest 3-outcome adjudication of the amendment** (ROUNDS_FOLD / CONTINUES_NO_TURNING /
GENUINE_ENDPOINT, or NOT_MEASURED/HALT) — NOT "rounds τ_fold." A ROUNDS_FOLD claim still requires an actual `dτ/ds`
reversal + original-residual convergence on both sides; the target is the DEEPER near-singularity (~τ=0.0291132), not the
refuted shallow τ_fold≈0.0291233.
1. B-1 gauge-fix removes the gauge near-null subspace AND is PROVEN path-only (same physical state with/without, lane-wise
   on the ORIGINAL residual); operators/frozen/export untouched (diff).
2. B-2 re-confirm explicitly returns one of FOLD_DISSOLVED / FOLD_CONFIRMED / STILL_INCONCLUSIVE (in that precedence) with
   numeric RELATIVE-TREND support (σ_min² fit + the bordered-cond TREND/ratio toward τ_fold — NOT the refuted absolute
   >1e10 bar). **[DONE: FOLD_DISSOLVED→STILL_GOING; the "B-3 only on FOLD_CONFIRMED" clause is superseded per the override above.]**
3. **[Re-scoped per the amendment]** If B-3 runs: it executes the staged B-3.0→B-3.1→B-3.2 and returns ONE honest verdict
   (ROUNDS_FOLD / CONTINUES_NO_TURNING / GENUINE_ENDPOINT / NOT_MEASURED) with numeric support, judged on the ORIGINAL
   residual ≤ tol, one branch, no frozen-param change. (The old "rounds τ_fold + quantified deeper throat" wording applies
   only to the ROUNDS_FOLD branch and is replaced by the amendment's solved-τ progress metric.)
4. B-4 analytic assembly is residual-equivalent (matches the original-residual JVP to tol).
5. Single Arbiter throughout (no scaled/bordered/least-squares surrogate used as the convergence/identity arbiter); NO
   LM/PTC/Sobolev conditioning remedy built (out of scope by the verdict); depth continuation moves only τ + state along
   the natural branch.
6. CPU; `timeout 600` per script (chunked + resumable; a timeout ⇒ NOT_MEASURED, NEVER raise the cap — if a single
   bordered Newton/linear solve at the resolution we need cannot fit 600s EVEN AFTER chunking + B-4 assembly + the
   modal/port form, HALT and bring it to the user, per the standing timeout flag); no `python3 -c`; YAML/markdown human
   output + JSON machine artifacts; chunk-1a/1b/1c gates pass; report + machine JSON emitted.

## Fail conditions
Touching the frozen residual / faithful operators / frozen physics / `physical_export_permitted`; a gauge-fix that changes
the physical solution (NOT path-only — fails the with/without equality); ~~building B-3 when B-2 ≠ FOLD_CONFIRMED~~
**(SUPERSEDED for the 2026-06-21 B-3 run by the user re-gate — see the amendment + the Acceptance override; declaring
ROUNDS_FOLD without an actual `dτ/ds` reversal + both-sides original-residual convergence IS still a fail)**;
pseudo-arclength that branch-jumps or changes a frozen parameter; judging convergence on the bordered/arclength/least-squares
norm instead of the ORIGINAL residual; **initializing any corrector from a recorded comparison state (cold-load disguised as
warm-start)**; building an LM/PTC/PTC-style/Sobolev remedy (scope creep the verdict disfavors); cold-loading states instead
of genuine warm-start/continuation; **declaring GENUINE_ENDPOINT without meeting the full exhaustion criteria (B-3.0 PASS +
B-4 reconfirm + adaptive-ds/reseed/secant tried)**; masking NOT_MEASURED; raising the timeout cap unilaterally.

## Out of scope (later, gated on this unblocking the crawl)
- The high-resolution production solve (RunPod A100/H100 GPU lever — POST-cure only; useless at 16×16).
- Then: promote the constitutive family → calibrated branch → multi-knob calibrate-predict (R0/J/W → anchor → SURPLUS) →
  `pathA_22`.

## Review (orchestrator, after Codex)
Fidelity agent (code-vs-spec, term-by-term): is the gauge-fix PATH-ONLY (in solver/preconditioner/bordering coords, never
the frozen residual) and is the with/without physical-equality proof real (not a tautology); is the ORIGINAL residual the
sole convergence/identity arbiter everywhere (not the bordered/least-squares norm); is the pseudo-arclength bordered system
correct (tangent predictor, bordered-Newton corrector, original residual in the physical block); is B-4 assembly genuinely
residual-equivalent (matches the original-residual JVP); are operators/frozen/export untouched (diff)? Adversarial agent:
is B-2's gate honestly applied (FOLD_CONFIRMED genuinely earned by the relative fold-trend gates — not asserted, and
closer-to-fold sampling attempted, not required; FOLD_DISSOLVED
not hidden to justify building B-3); does B-3 ACTUALLY round the fold and reach a quantified deeper throat (converged on
the original residual both sides), or is "rounded" over-read; any branch-jump / frozen-param drift / can't-fail gate /
hardcoded depth claim / surrogate-norm convergence? Then gate the next phase (constitutive family → calibrate-predict).
