# Directive pathA_C0g-build — gauge-fix (path-only) → re-confirm the fold → gauge-fixed pseudo-arclength to round the throat-radius turning point

**Status:** **READY — awaiting the user execution gate.** Loop complete: Codex design-review = SOUND-WITH-FIXES → all 8
fixes APPLIED → confirm-pass = STILL-NEEDS (3 residuals) → 3 residuals FIXED + the `min_R0` anchor self-verified (`:1539`)
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
  depth metric **`min_R0`** (already recorded by the crawl as `min_R0_during_final`, `patha_c0_conditioning_spike.py:1539`) — report `min_R0` past
  the fold and require a stated numeric improvement margin below the τ=0.029125 state's `min_R0`, FIXED BEFORE execution
  (not selected post-hoc; secondary metrics min ρ / throat depth reported for context); (iii) it
  traces ONE physical branch (no branch-jump; no frozen-param change — only τ and the state move along the natural
  continuation); (iv) operators/frozen/export untouched.

### B-4 — Analytic/sparse Jacobian assembly (unconditional, performance/conditioning)
After a color audit, provide analytic/sparse Jacobian assembly to replace (or augment) the colored-JVP dense path,
handling the dense μ/mass/wall rows (253 colors = deterministic radius-3 coloring, NOT a bug). Residual-equivalent with a CONCRETE
match gate: the assembled J must match the JVP of the ORIGINAL residual to concrete tolerances FIXED BEFORE execution
(default: relative ≤ 1e-8 AND absolute ≤ 1e-10 per block; tighten if the float64 JVP supports it), with
coverage = multiple random `Jv` probes PLUS the dense wall/mass/μ rows checked explicitly (the rows the colored path
handles specially); report the per-block max abs/rel mismatch. This is a perf/conditioning aid, NOT a physics change.

## Key facts (anchors — Codex re-confirm at execution)
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
1. B-1 gauge-fix removes the gauge near-null subspace AND is PROVEN path-only (same physical state with/without, lane-wise
   on the ORIGINAL residual); operators/frozen/export untouched (diff).
2. B-2 re-confirm explicitly returns one of FOLD_DISSOLVED / FOLD_CONFIRMED / STILL_INCONCLUSIVE (in that precedence) with
   numeric RELATIVE-TREND support (σ_min² fit + the bordered-cond TREND/ratio toward τ_fold — NOT the refuted absolute
   >1e10 bar); B-3 is implemented ONLY on FOLD_CONFIRMED.
3. If B-3 runs: pseudo-arclength rounds τ_fold (converged states both sides, original-residual ≤ tol), reaches a quantified
   deeper throat, traces one branch with no frozen-param change.
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
the physical solution (NOT path-only — fails the with/without equality); declaring FOLD_CONFIRMED without the
relative fold-trend numeric support; building B-3 when B-2 ≠ FOLD_CONFIRMED; pseudo-arclength that branch-jumps or changes a
frozen parameter; judging convergence on the bordered/least-squares norm instead of the ORIGINAL residual; building an
LM/PTC/Sobolev remedy (scope creep the verdict disfavors); cold-loading states instead of genuine warm-start/continuation;
masking NOT_MEASURED; raising the timeout cap unilaterally.

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
