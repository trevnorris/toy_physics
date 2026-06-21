# Directive pathA_C0g-build — gauge-fix (path-only) → re-confirm the fold → gauge-fixed pseudo-arclength to round the throat-radius turning point

**Status:** **B-3 FIRST RUN DONE & COMMITTED (`46a3ea21`) → NOT_MEASURED at the B-3.0 gate (honest). GLM-REFRAMED:
the "decisive" overlap test was near-tautological → the CURRENT plan is the "⭐⭐ B-3 FOLLOW-UP v2" CHARACTERIZATION BATTERY.
EXECUTE via `_scratch/pathA_C0g_B3_followup_execution_prompt.md` once the v2 design-review loop closes SOUND-AS-IS, then the
user EXECUTION gate (2026-06-21).**
- **B-3 first run:** faithful gauge-fixed Keller pseudo-arclength (fidelity FAITHFUL_WITH_NOTES); reproduces the recorded
  branch to ~1e-6 at the reliable points but "FAILs" the B-3.0 gate at the deepest point (compared vs an under-converged
  recorded Newton state) → `NOT_MEASURED` (honest; dual-reviewed). Report `reports/pathA_C0g_B3_pseudoarclength.md`.
- **⚠️ The earlier "Phase D overlap → gate-refine → re-run" follow-up is SCRAPPED** (GLM + Claude + Codex agree: the overlap
  test is a NECESSARY consequence of two near-roots, so it can't separate gate-artifact from continuation-bug). **CURRENT =
  the CHARACTERIZATION BATTERY** (see "⭐⭐ B-3 FOLLOW-UP v2"): C0 read-offs → C1 re-converge-the-FAIL-state-TIGHT-from-both →
  C2 mode-characterization (+ bordered fold-transversality) → C3 line-scan → C4 resolution check; a 6-outcome tool-selection
  table (Case 0 INCONCLUSIVE / 1 under-converged-ref / 2 fold / 3 bifurcation / 4 wall / 5 gauge-or-discretization); **do NOT
  pre-commit to pseudo-arclength** — only Case 2 (with a user gate) lets the built pseudo-arclength continue. GLM consult DONE
  (`_scratch/pathA_C0g_B3_followup_for_GLM.md`); v2 design-review loop: Codex AGREES with GLM, SOUND-WITH-FIXES (12) → fixes
  applied → confirm-pass [pending]. Logs `_scratch/codex_c0g_B3_followup_v2_design_review.log`. Timeout: 600s cap LIFTED for
  solver runs → forward-progress monitoring (see that policy section).
- **(History: prior B-3 RE-SCOPE amendment loop + the scrapped Phase-D follow-up loop, 2026-06-21)** — design-reviewed
  SOUND-AS-IS at the time; the Phase-D one was later superseded by the GLM reframe.
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

## ⭐⭐ B-3 FOLLOW-UP v2 (AMENDMENT, 2026-06-21, GLM-REFRAMED) — CHARACTERIZE the singularity FIRST, then pick the tool

**This supersedes the earlier "decisive A-vs-B diagnostic (Phase D) → gate-refine → re-run" follow-up, which a GLM tertiary
review (and Claude+Codex on reflection) judged FLAWED.** Why the prior Phase D was scrapped: the proposed "decisive" test —
project `Δ = x_cont − x_rec` onto the near-null subspace of `J·Q_perp` and require overlap ≥ 0.90 — is **near-tautological**.
If `x_cont` and `x_rec` are BOTH near-roots of the same residual at the same τ, then `J·Δ ≈ R(x_cont) − R(x_rec) ≈ 0` by the
mean-value theorem, so `Δ` lies in the near-null space REGARDLESS of whether the continuation tracked the branch correctly or
drifted there via a bug (checked: `‖J·Δ‖/‖Δ‖ ≈ 2.2e-7/6.8e-4 ≈ 3e-4` — high overlap was guaranteed). The "false-confirmation
guard" just restated "both are roots." So Phase D is NECESSARY-not-SUFFICIENT and cannot separate A (gate artifact) from B
(continuation bug). **Lesson recorded:** a decisive test must not be a necessary consequence of the thing it claims to confirm.

**The reframe (GLM): characterize the singularity, don't optimize the gate.** The B-3.0 gate failure is a SYMPTOM; the disease
is an UNCHARACTERIZED near-singularity at τ≈0.0291132. The framing was also incomplete — there are FIVE possibilities, not two:
(1) gate artifact / under-converged reference (the recorded state hit the 20-iter cap, residual 2.2e-7, step_norm 1.08e-6 =
slow-not-stalled at cond~1e5); (2) physical simple FOLD (σ_min→0 with a tangent reversal reachable past the stall);
(3) physical BIFURCATION / branch point (the non-monotone r0/μ + r0 nudging UP + "valid roots differing from Newton's" are
hallmarks; pseudo-arclength is the WRONG tool here — it silently branch-switches); (4) conditioning WALL (σ_min bottoms out
>0, never reaches 0 — read-off: the deepest CONVERGED σ_min is a FINITE 7e-6 at τ=0.0291132, already BELOW the σ_min² fit's
zero-crossing 0.0291139, so the fit is unreliable and a wall is live); (5) incomplete GAUGE fix (a residual r0/μ gauge
freedom the 511-dim fix missed) or a 16×16 DISCRETIZATION artifact. **Plus an operational (0) INCONCLUSIVE / bug** — an
implementation, provenance, scaling, or gate defect is a real EXECUTION outcome that must STOP for debug, NOT be forced into
(1)–(5). **Pseudo-arclength is right only for (2); each other case needs a different tool (or a stop). So we MUST characterize
before committing the expensive deep solve.**

**Two read-offs already done (Claude, from existing logs; Codex re-confirms exact in C0):** raw σ_min(τ) — finite 7e-6 at the
deepest converged τ, fit zero-crossing unreliable (wall vs fold UNRESOLVED); arclength-metric balance — near the stall the
metric is x-DOMINATED (`‖dx‖²≈3.6e-6` vs `(1000·dτ)²≈6e-7`), so pseudo-arclength was NOT degenerating into τ-continuation.

### THE CHARACTERIZATION BATTERY (run cheap→expensive, GATED; Codex designs each; each step reports before the next)
- **C0 — confirm the two read-offs exactly** (raw σ_min(τ) array shape: bottoming-out vs trending-to-0; exact `‖dx‖²` vs
  `(1000·dτ)²` per B-3.0 step from the actual scaled-state vectors). Cheap, on existing artifacts.
- **C1 — RE-CONVERGE the FAIL-point state (THE decisive A-vs-B test; cheap).** At τ=0.02911625 run plain (gauge-fixed) Newton
  warm-started SEPARATELY from (i) the recorded deepcrawl state `attempt_008_tau_0p02911625.npz` and (ii) the B-3.0
  continuation state `B3_0_accepted_009`. **TIGHT convergence required (MF1):** original residual ≤ **1e-11** (NOT 1e-6 — the
  recorded state already sits at 2.2e-7, so a 1e-6 bar lets the solver "stop immediately" on the stale state) AND
  Newton step/update norm ≤ **1e-9**; the iteration budget (≥100) is a BUDGET, not success. Emit the full iteration trace
  (residual + step_norm per iter). Report both end-states + their gauge-invariant (ρ, curlA, r0, μ) comparison. **Pre-committed
  reading (MF2 — asymmetric/partial outcomes are EXPLICIT, not forced into 1/2/3):**
  - BOTH tight-converge to the SAME state (gauge-inv ≤ 5e-5) ⇒ **Case 1** (under-converged-reference ARTIFACT at THIS τ —
    "same root here," NOT "continuation globally fine"; C2/C4 still decide the deeper fold-vs-wall).
  - BOTH tight-converge to gauge-inv DIFFERENT states ⇒ **Case 3** (genuine multiple roots / bifurcation).
  - NEITHER reaches tight tol at the budget ⇒ a real deep singularity (**Case 2 or 4**; C2 distinguishes).
  - **ONE tight-converges and the other does NOT, OR both improve but miss tight tol, OR one start drifts down the soft mode**
    ⇒ **Case 0 / INCONCLUSIVE (conditioning evidence)** — do NOT force a physical case.
  (NICE: re-converge not just the deepest FAIL point but EVERY unreliable comparison state, for the gate refinement.)
- **C2 — CHARACTERIZE the near-null mode(s) at the stall (informs the TOOL — SUPPORTING evidence, not sole arbiter).** From the
  dense SVD of `row_scale·J·col_scale·Q_perp` at the deepest tight-converged state: (i) report the smallest ~10 σ + a
  predefined near-null CLUSTER criterion (by σ gap/ratio, e.g. a >10× gap to the next σ) — one isolated near-null σ is
  fold-LIKE, a cluster is bifurcation-LIKE, but **count is SUPPORTING only (MF3)**; (ii) the dominant right-singular vector's
  SECTOR decomposition (r0/μ/ψ/A energy fractions); (iii) **a bordered FOLD-TRANSVERSALITY check (MF3):** form the left null
  vector `w` and test `wᵀ F_τ ≠ 0` (simple fold) vs ≈ 0 (bifurcation/degenerate) via the bordered system / bordered
  conditioning — this, not the σ-count, is the fold-vs-bifurcation discriminant; (iv) **GAUGE test (MF4 — avoid
  metric-tautology):** the mode comes from `J·Q_perp`, so projecting it back onto the SAME analytic gauge basis used to build
  `Q_perp` can be forced small — require BOTH the scaled-metric projection onto the REMOVED gauge basis AND an INDEPENDENT
  expanded gauge-candidate test (a broader ∇λ/phase family) PLUS a curl / residual-equivalence check; high physical-metric
  overlap ALONE does not prove "missed gauge"; (v) localization (throat-concentrated ⇒ fold-like vs extended ⇒
  bifurcation-like). Read-only on the tight C1 states.
- **C3 — LINE-SEGMENT residual scan (SUPPORTING physical check, not proof; MF5).** GAUGE-ALIGN x_rec to x_cont first; evaluate
  `‖R(x(t),τ)‖` for `x(t)=(1−t)x_rec + t·x_cont`, t∈[0,1], 20–50 points, at τ=0.02911625, normalized against the TIGHT C1
  endpoint residuals (not the stale ones). A clear SPIKE ⇒ bifurcation / barrier / bug / INCONCLUSIVE; FLAT ⇒ CONSISTENT with
  one near-root manifold but NOT proof (a straight line is not a solution-manifold path; coordinate curvature + tiny σ_min can
  flatten it). SUPPORTING evidence only.
- **C4 — RESOLUTION check (is it even real?; MF6).** RE-LOCATE the comparable stall/feature at 24×24 (and 32×32 if affordable)
  — same gauge/scaling/tolerances, find the feature's τ per grid (it may SHIFT), then compare NORMALIZED σ_min trends + mode
  localization across resolutions. σ_min trend GROWS with resolution ⇒ DISCRETIZATION-artifact-SUSPECT → re-gate (not "the rest
  is moot"); PERSISTS (normalized) ⇒ physical. (The more expensive step; gate it on C1–C3 not already settling the question.)

### TOOL SELECTION (the battery's OUTPUT decides; do NOT pre-commit to pseudo-arclength)
| Characterization result | Right tool / next step |
|---|---|
| **Case 0 — INCONCLUSIVE / implementation-provenance-scaling-or-gate bug** (mixed C1/C2/C3 evidence; asymmetric C1; provenance mismatch) | **STOP — debug / re-gate with the user. Do NOT select pseudo-arclength, deflation, or LM.** |
| Case 1 — under-converged reference (C1: both tight→same state) | Re-converge the recorded deep states tightly; refine the B-3.0 gate (criteria below); the deeper fold-vs-wall question STAYS open → C2 transversality + C4 decide. |
| Case 2 — physical simple fold (C2: isolated near-null mode, `wᵀF_τ≠0`, localized; C4: persists) | Gauge-fixed pseudo-arclength (already built) — with the Phase-F fallbacks ACTUALLY wired (below). |
| Case 3 — bifurcation (C1: distinct tight roots; C2: cluster / `wᵀF_τ≈0`; C3: barrier) | Deflation / branch-switching — NOT plain pseudo-arclength. New build, re-gate with the user. |
| Case 4 — conditioning wall (C1: neither tight-converges; C2: isolated mode but σ_min bottoms out) | Levenberg–Marquardt / trust-region damping — re-gate with the user (LM is allowed ONLY as this gated Case-4 tool). |
| Case 5 — gauge mode (C2: independent-test gauge overlap high) OR discretization (C4: normalized σ_min grows) | Extend the gauge fix, or refine the grid — the "singularity" is not physical. |

**Only Case 2 lets the already-built pseudo-arclength continue (and only with a user gate); every other case STOPS for a
re-gate.** **If (and only if) Case 1 or Case 2**, the principled B-3.0 gate refinement applies (reference-reliability filter:
a recorded comparison state is RELIABLE iff original-residual ≤ 1e-8 AND Newton-iters < 20 AND step_norm ≤ 1e-7; missing
metadata ⇒ NOT reliable; at unreliable references require only continuation-residual ≤ min(recorded,1e-6) + numeric
branch-continuity + stiff-lane (ρ,curlA) agreement ≤ 5e-5; drop soft r0/μ equality ONLY at points C1 tight-re-convergence +
C3 show are genuine same-manifold) + the Phase-F fallback wiring (adaptive ds / tangent reseed / secant fallback with
ACTUAL-USE evidence, not booleans — note the existing bug where `secant_fallback_attempted` can be true while the predictor
still gets `None`), THEN the pseudo-arclength re-run with an honest B-3.2 verdict (ROUNDS_FOLD / CONTINUES_NO_TURNING /
GENUINE_ENDPOINT / NOT_MEASURED). GENUINE_ENDPOINT requires fallbacks demonstrably FIRED-and-FAILED (evidence, not flags).

### MONITORING & TIMEOUT POLICY (user-granted 2026-06-21 — supersedes the 600s cap FOR ALL THROAT-SOLVER NUMERICAL RUNS)
**This policy SUPERSEDES every earlier `timeout 600` / "never raise the cap" mention in THIS directive for solver runs** (the
prior B-3 amendment's anti-gaming item (h), Acceptance #6, the standing-flag references) — those are legacy / `.wl`-SymPy-only.
The hard cap was a forcing-function for the `.wl`/SymPy DERIVATION scripts (which KEEP it); the user LIFTED it for the
throat-solver NUMERICAL runs. Replace it with FORWARD-PROGRESS monitoring:
- NO hard wall-clock cap; a multi-hour run is FINE **as long as it is making forward progress** toward results we need.
- Every solver script MUST emit INCREMENTAL progress (append per accepted state / per Newton iter: τ, original-residual,
  min_R0, σ_min, iters, wall-time) to its run JSON/log so progress is observable WHILE running, and MUST checkpoint resumably.
- Each script MUST carry an INTERNAL no-forward-progress guard with an EXACT pre-stated criterion: STOP + report STALLED when
  **N = 30** consecutive Newton/continuation iterations show neither (a) original-residual decreasing by ≥ a relative
  `ε_res = 1%` over the window, NOR (b) the continuation parameter advancing (Δτ or Δarclength ≥ `ε_adv = 1e-9`). For the C1
  FIXED-τ re-converge (where τ does not advance), the guard is criterion (a) ALONE (residual must keep decreasing); reaching
  the iteration budget without tight tol ⇒ report NOT-tight-converged (a C1 outcome), not silently "converged." (Tune N/ε in
  the run only UP, with the value logged; never to mask a stall.)
- The orchestrator MONITORS the emitted progress periodically and TaskStops a run that is genuinely stuck (no progress), but
  lets a progressing run continue however long it needs.
- This applies to the Path-A throat-solver numerical work ONLY; the `.wl`/SymPy derivation/audit scripts KEEP the 600s cap.
  ([[feedback-script-timeout-policy]] updated; the decisions/13 standing timeout flag updated.)

**Acceptance (follow-up v2):** C0 read-offs confirmed; C1 returns a pre-committed case reading with the two re-converged
states compared gauge-invariantly; C2 characterizes the mode(s) (count, sector, gauge-overlap, localization); C3/C4 as gated;
the TOOL is chosen by the table from the EVIDENCE (not pre-committed); any gate refinement happens ONLY in Case 1/2 with the
principled criteria; Single Arbiter / gauge-inside-every-solve / freeze / no-cold-load all hold; forward-progress monitoring
replaces the 600s cap. **Fail conditions (follow-up v2):** committing to pseudo-arclength before the battery characterizes the
singularity; refining the gate outside Case 1/2; treating the near-tautological overlap as decisive; post-hoc threshold
tuning; any prior Fail condition (freeze/operators/export/Single-Arbiter/cold-load).

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

**⭐⭐ 2026-06-21 B-3 FOLLOW-UP v2 ACCEPTANCE OVERRIDE (governs the CURRENT run — supersedes the criteria below for it):**
The current run is the CHARACTERIZATION BATTERY, NOT a pseudo-arclength fold-round. Its acceptance is the "Acceptance
(follow-up v2)" block in the v2 section (C0 read-offs confirmed; C1 TIGHT re-converge-from-both with a pre-committed
6-outcome reading incl. Case 0/INCONCLUSIVE; C2 mode-characterization + bordered fold-transversality; C3/C4 as SUPPORTING
evidence; the TOOL chosen from the table by EVIDENCE, never pre-committed; gate refinement only in Case 1/2). Criteria #2/#3
below (old B-1/B-2/B-3 pseudo-arclength) are HISTORY for this run; #1 (B-1 path-only) and #4 (B-4 residual-equivalent) still
hold where reused; #5 holds with the LM reconciliation noted; #6's 600s cap is REPLACED by the forward-progress monitoring
policy. Pseudo-arclength's ROUNDS_FOLD/etc. acceptance applies ONLY if the battery returns Case 2 AND the user gates it.
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
5. Single Arbiter throughout (no scaled/bordered/least-squares surrogate used as the convergence/identity arbiter). **LM/PTC
   reconciliation (MF12):** NO LM/PTC/Sobolev remedy may be BUILT during THIS characterization run or without a user re-gate;
   LM/trust-region is permitted ONLY later as the explicitly-gated Case-4-SELECTED tool (per the tool-selection table). Depth
   continuation moves only τ + state along the natural branch.
6. CPU. **[For the B-3 FOLLOW-UP v2 solver runs the 600s cap is SUPERSEDED by the user-granted forward-progress monitoring
   policy — see that section; no hard wall-clock cap, multi-hour OK if progressing, internal no-progress guard + orchestrator
   monitoring instead.]** (Legacy, still binding for `.wl`/SymPy derivation scripts:) `timeout 600` per script (chunked +
   resumable; a timeout ⇒ NOT_MEASURED, NEVER raise the cap — if a single
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
**No-cold-load scope (MF10):** the C0–C3 DIAGNOSTICS legitimately LOAD recorded/continuation artifacts (C1 warm-starts from
BOTH; C3 reads both end-states) — that is allowed and expected; the prohibition is that any CONTINUATION / tool CLAIM (and the
B-3.0 reproduction) may NOT be INITIALIZED from a comparison artifact except the single declared seed. Report artifact
provenance for every loaded state (which file, used as seed vs comparison).
**Follow-up v2 Fail conditions:** committing to pseudo-arclength before the battery characterizes the singularity; forcing a
Case-0/asymmetric/inconclusive C1 outcome into a physical case; relaxing the C1 tight-tol (≤1e-11) to let the stale state
"converge"; treating C2 σ-count or C3 flatness as PROOF rather than supporting evidence; refining the gate outside Case 1/2.

## Out of scope (later, gated on this unblocking the crawl)
- The high-resolution production solve (RunPod A100/H100 GPU lever — POST-cure only; useless at 16×16).
- Then: promote the constitutive family → calibrated branch → multi-knob calibrate-predict (R0/J/W → anchor → SURPLUS) →
  `pathA_22`.

## Review (orchestrator, after Codex)
**For the B-3 FOLLOW-UP v2 CHARACTERIZATION BATTERY (current run):** Fidelity agent (code-vs-spec): is C1 a GENUINE tight
re-converge (≤1e-11 + step-norm ≤1e-9, full trace) from BOTH seeds — NOT a stop-on-the-stale-state; is the gauge-invariant
comparison correct; is C2's bordered fold-transversality (`wᵀF_τ`) and the INDEPENDENT-basis gauge test implemented (not the
metric-tautological self-projection); is C3 gauge-aligned + normalized to tight endpoints; is C4 re-locating the feature
per-grid (not same-nominal-τ); is the no-progress guard exact (N=30, ε_res=1%, ε_adv=1e-9; C1 residual-only); Single
Arbiter + freeze + provenance intact (diff)? Adversarial agent: is the Case reading HONEST (asymmetric/partial C1 → Case 0,
not forced; σ-count + line-scan flatness treated as SUPPORTING not proof; gauge overlap not metric-tautological); is the TOOL
chosen from EVIDENCE with NO pre-commitment to pseudo-arclength (only Case 2 + user gate proceeds); any forced case /
can't-fail / surrogate-norm / hidden cold-load? **Then STOP and bring the Case + chosen tool to the user (re-gate)** — only a
clean Case 2 + user gate runs the built pseudo-arclength.
**(Legacy review notes, for the pseudo-arclength run IF Case 2 gates it):** is the gauge-fix PATH-ONLY; the bordered system
correct (tangent predictor, bordered-Newton corrector, original residual in the physical block); B-4 residual-equivalent;
does it ACTUALLY round the fold (converged both sides on the original residual) or is "rounded" over-read; any branch-jump /
frozen-param drift / surrogate-norm convergence? Then gate the next phase (constitutive family → calibrate-predict).
