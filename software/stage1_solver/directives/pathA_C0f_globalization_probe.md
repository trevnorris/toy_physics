# Directive pathA_C0f — Is the τ≈0.029 "wall" just a crippled-solver config? (globalization-first probe + fold detector)

**Status:** READY (Claude-authored 2026-06-20; design-review = SOUND-WITH-FIXES → all 9 fixes applied [WALL_WAS_CONFIG
requires an ACCEPTED default crawl to τ=0.028 not a sweep; numeric fold rule + FOLD_RISK; "passing 0.029 ≠ sufficient";
reporting/adapter code only; report both solver_converged + b2c_linf_pass; smoothness wording + positive-R0 caveat;
ignore dead C0e curl labels; full default-config table incl. aids shown path-only; 1−P_G range]). Then the user EXECUTION
gate. Converged plan: `software/stage1_solver/_scratch/pathA_C0f_agreed_plan_for_GLM.md` (Claude+Codex §5–§7) + GLM 5.2
review (run-defaults-first; cheap dψ/dτ fold detector; smoothness; specify the depth target). Follows `pathA_C0e` whose
`GAUGE_FRAMING_REFUTED` verdict was REJECTED (k-biased-metric false negative). Step of option C, task #78.

**Date:** 2026-06-20
**Owner:** Codex (codes + iterates until exit 0; drives the EXISTING closed Newton crawl + adds bounded diagnostic/
reporting code; extends the C0e machinery for the gauge re-confirm). Claude reviews after (fidelity + adversarial agents).
**Trigger (the reframe, verified):** the C0b run that DEFINED the τ≈0.029 "wall" used `max_newton_iters_override=2`
(C0 default = 20, `patha_c0_conditioning_spike.py:282`; closed-Newton default = 18, `patha_closed_newton.py:81`),
`max_tau_backtracks=0` (default 5, `:153`), and `depth_sequence=[0.03, 0.02899]` (default = 8-element fine-halving ending
at τ=0.028, `:133-142`); the spike's own `PRODUCTION_SOLVER_REQUIRED` logic (`:2180-2182`) requires full budget + real
backtracking + below-floor failures, NONE of which C0b had. Separately, C0e-0 showed the near-null subspace is NOT what
blocks the Newton step (full step ~entirely in the physical complement, fraction 1.8e-21; removing the gauge part changes
‖F‖ not at all; clean solve 7.9e-13). So the "wall" is most likely a crippled-config / step-length artifact, NOT physics.
Cheapest-first test of that.

## Scope (this step RUNS the solver — boundary is Single-Arbiter + reporting-only, not "no-touch")
C0f RUNS the EXISTING closed Newton crawl (which ALREADY has Armijo backtracking line search, `newton.py:189` /
`patha_c0_conditioning_spike.py:562`, `max_line_search_iters=20`) with the CODE DEFAULTS, plus bounded diagnostic +
REPORTING/ADAPTER code only. It MUST NOT modify the solver logic, the physical residual `patha_closed_branch_residual`,
the faithful PDE operators, frozen physics, or `physical_export_permitted`; it MUST NOT change any FROZEN PHYSICS
parameter (`a`, `L`, `r_mouth`, `r_exit`, `w_max`, boundary kinds, constitutive shapes) — depth continuation is
**τ-only**. The **Single Arbiter Principle** holds: the ORIGINAL unmodified `patha_closed_branch_residual` is the SOLE
arbiter of convergence; merit / line-search / acceptance are judged on the original ‖F‖, never a scaled/floored surrogate;
the default conditioning aids (ε preconditioner schedule, Jacobi scaling, wall predictor, final-only EOS continuation)
are PATH-ONLY and the FINAL accepted solve is at zero ε (must be reported + shown path-only). C0f adds NO new numerics:
**trust-region/dogleg/LM is NEW numerics, DEFERRED to the gated C0g**; pseudo-arclength is DEFERRED (only if a fold is
detected). CPU; `timeout 600` per script (split if needed; at 16×16 the dense solve is sub-second/iter so the default
8-step crawl fits; timeout → NOT_MEASURED); standalone `python3`; no commentary `python3 -c`; YAML/markdown human output,
JSON only for machine artifacts; chunk-1a/1b/1c gates must still pass.

## Key facts (confirmed in the consult + GLM review + Codex design-review; Codex re-verify at execution)
- **The residual is SMOOTH on the monitored domain** (GLM/Codex, code-checked): the ψ-density nonlinear term is polynomial
  — `quintic_enthalpy=(5K/4)·ρ⁴`, `ρ=|ψ|²` (`physics.py:46-47`) ⇒ `(5K/4)·|ψ|⁸·ψ`, C^∞ in `(ψ_r,ψ_i)`; NO `√ρ`, no `|ψ|`,
  no density-division in `coupled_pde_residual` (`coupled_branch.py:381-388`). The full closed residual additionally has
  smooth positive-`R0` rational terms — smooth on the positive-`R0` domain the crawl monitors (min_R0≈0.75 in C0b). With
  `Jδ=−F` consistent (clean solve), the L2 merit `½‖F‖²` has slope `−‖F‖²` at α=0, so SOME small α MUST reduce ‖F‖ unless
  the Jacobian/residual pair is stale/inconsistent. ⇒ the α=1 overshoot (C0e-0: 5.08×) is pure large-step nonlinearity;
  the predicted-vs-actual residual should show a clean O(α²) gap. A DEVIATION would be the real surprise.
- **The fold is the one genuine residual risk.** A turning point near τ≈0.029 cannot be passed by any globalization (the
  branch turns back). Pseudo-arclength (the right tool) does NOT exist in the codebase (string literal only, `:2135`).
  C0f-0a is a cheap, NUMERIC dψ/dτ detector run before any sweep.
- **The near-null cluster is gauge-dominated and NOT the proximate blocker** (C0c phase mode confirmed; C0e-0). Mode 2
  (~17% transverse) is a real but DEFERRABLE candidate (not blocking the Newton step). Gauge deflation is cheap INSURANCE
  for later/higher-resolution linear solves, not part of C0f.

## Work items (cheapest-first; later steps run ONLY if the earlier one stalls)

### C0f-0 — run the DEFAULT config crawl (the cheapest, most decisive test)
From the converged τ=0.03 state, run the EXISTING closed Newton τ-continuation crawl with the CODE DEFAULTS (full Newton
budget = 20, `max_tau_backtracks = 5`, the default 8-element fine-halving `depth_sequence` ending at τ=0.028, the existing
Armijo line search) — do NOT apply any C0b override. **Report the FULL default config used** (incl. the ε preconditioner
schedule, Jacobi scaling, wall predictor, final-only EOS continuation), sourced from code, and confirm each aid is
PATH-ONLY (Newton residual/JVP/merit/convergence on the ORIGINAL `patha_closed_branch_residual`; FINAL accepted solve at
zero ε). C0f may add REPORTING/ADAPTER code ONLY (recompute final L2 via the original residual; derive smallest accepted α
from the zero-ε attempt history) — solver logic UNCHANGED. Report, per τ in the sequence: `solver_converged` (the zero-ε
Newton solve internally converged) AND `b2c_linf_pass` (Linf ≤ 1e-6) — report BOTH; a Linf≤1e-6 row that did NOT
internally converge is INCOMPLETE/suspect, not success. Also per τ: final Linf and L2 (original residual), Newton iters
used, backtracks used, smallest accepted α. **Depth target:** passing τ=0.029 is encouraging but NOT sufficient — the
verdict requires `solver_converged` AND `b2c_linf_pass` (Linf ≤ 1e-6) THROUGH τ=0.028. If it reaches 0.028, note how much
deeper it goes before any new stall.

### C0f-0a — fold detector (cheap, NUMERIC; run before any sweep)
From the converged states the C0f-0 crawl produces (and/or recomputed states at τ=0.0295 and 0.029), compute the
normalized solution sensitivity `‖Δx/Δτ‖ / ‖x‖` over CONSECUTIVE converged τ-intervals, with a per-lane breakdown
(ψ, R0, μ). FOLD_DETECTED requires BOTH: (a) `‖Δx/Δτ‖` grows monotonically by ≥ a stated factor (default ≥10×) across the
last 2–3 intervals approaching the stall, AND (b) Newton fails at the next-smaller τ for ALL backtracks/step sizes.
Growth without (b), or only one interval available, → `FOLD_RISK` (report the numbers, inconclusive). Insufficient
converged intervals → `DIAGNOSTIC_INCOMPLETE` for the fold test. Report the actual factor and the lane breakdown.

### C0f-1 — merit sweep with predicted-vs-actual (DIAGNOSTIC ONLY; run only if C0f-0 still stalls at some τ)
The default crawl ALREADY backtracks, so this step does NOT "find an α the solver missed" — it DIAGNOSES why a stall
remains. At the first τ where C0f-0 stalls, solve `Jδ=−F` on the reassembled UNSCALED original-residual Jacobian at the
exact stalled state (report linear rel-resid) and sweep `α ∈ {1, ½, ¼, …, 2⁻²⁰}`: report, on the ORIGINAL residual, BOTH
`‖F(x+αδ)‖` (L2 and Linf) AND the predicted linear residual `‖F+αJδ‖`. For a smooth residual + good direction, actual
tracks predicted with a clean O(α²) gap and some α reduces ‖F‖ — that is LOCALIZED GLOBALIZATION EVIDENCE (the stall is
step-control / line-search budget), NOT a declared fix. A DEVIATION (no α reduces ‖F‖, or actual ≫ predicted even for tiny
α) ⇒ stale/inconsistent Jacobian or predictor/branch geometry — report it as the localized blocker.

### C0f-2 — unbiased gauge re-confirmation (cheap; settle the C0e mislabel; run in parallel)
Replace the k-biased curl classifier with the dimensionless transverse-energy fraction `1−P_G` (and `1−P_cpl`) for the 4
Maxwell-lane modes + the phase mode, on the saved stall matrix (reuse the C0e machinery — but IGNORE the old C0e
`outcome`/`GAUGE_FRAMING_REFUTED` labels; recompute C0f classifications from `1−P_G`/`1−P_cpl`/controls). Report `1−P_G`,
`1−P_cpl`, the A-normalized transverse fraction, and the boundary-vs-interior curl SPLIT (to confirm the ~53% boundary
curl is one-sided stencil-closure leakage, not interior field). Add a MIXED control: an exact discrete gradient + a known
ε·transverse perturbation at several k, to show such a mode lands at `1−P_G≈ε²` (NOT a "transverse" band) — demonstrating
the curl gate's k-bias directly. Keep raw curl ONLY as a reported stiffness observable, never the classifier. Expected:
modes 1/3/4 gauge (`1−P_G` ≈ 0.02–0.7%), mode 2 the one transverse candidate (~17%).

### C0f-3 — VERDICT
Exactly one:
  - **WALL_WAS_CONFIG / GLOBALIZATION_FIXES_IT** — the DEFAULT-config crawl (C0f-0), with its EXISTING accepted Armijo
    line search, ACTUALLY reaches `solver_converged` AND Linf ≤ 1e-6 THROUGH τ=0.028 (not merely a diagnostic α found in
    C0f-1), with `‖Δx/Δτ‖` bounded (no fold). ⇒ the crawl is UNBLOCKED; report depth reached + recommend continuing the
    crawl toward the physics-target depth, with minimal gauge deflation as optional insurance.
  - **FOLD_DETECTED** — C0f-0a's numeric rule (a)+(b) both hold near the stall. ⇒ recommend C0g = pseudo-arclength
    continuation (NOT yet in the codebase); NOT a production solver. (`FOLD_RISK` if only (a): report, do not over-claim.)
  - **GLOBALIZATION_INSUFFICIENT** — the default crawl stalls, the C0f-1 sweep finds NO α reducing the true residual
    (and/or actual ≫ predicted), smaller-τ increments still plateau ≫ 1e-6 with CLEAN linear solves and bounded
    `‖Δx/Δτ‖`. ⇒ localize via the predicted-vs-actual gap (stale Jacobian / predictor geometry / assembly bug);
    recommend the specific follow-up. Still NOT a production linear solver.
  - **DIAGNOSTIC_INCOMPLETE** — required evidence NOT_MEASURED (timeout, too few converged intervals for the fold test).
Report falsifiable numeric `verdict_support`: the per-τ crawl table (solver_converged/b2c_linf_pass/Linf/L2/iters/
backtracks/α + the full default config shown path-only), the `‖Δx/Δτ‖` fold numbers + lane breakdown, the merit-sweep
predicted-vs-actual table (if run), and the `1−P_G` gauge re-confirmation.

## Acceptance criteria (PASS/FAIL; exit-0 NECESSARY not SUFFICIENT)
1. C0f-0 runs the DEFAULT config (no C0b overrides; state the actual `max_newton_iters`, `max_tau_backtracks`,
   `depth_sequence`, line-search, ε schedule, Jacobi scaling, wall predictor, EOS continuation — sourced from code — and
   confirm each aid path-only with the FINAL accepted solve at zero ε). Per-τ it reports BOTH `solver_converged` AND
   `b2c_linf_pass` on the ORIGINAL residual, plus Linf/L2/iters/backtracks/smallest-α. Convergence judged on the
   unmodified residual (Single-Arbiter). Solver logic UNCHANGED (reporting/adapter code only).
2. C0f-0a reports the NUMERIC `‖Δx/Δτ‖/‖x‖` over consecutive converged intervals + ψ/R0/μ lane breakdown, and applies the
   (a)+(b) rule → FOLD_DETECTED / FOLD_RISK / DIAGNOSTIC_INCOMPLETE (no qualitative "diverging" call).
3. C0f-1 (only if a stall remains) reports the α-sweep with BOTH actual `‖F(x+αδ)‖` and predicted `‖F+αJδ‖` (from the
   reassembled UNSCALED original-residual Jacobian) + the linear rel-resid, and states direction-good (tracks, some α
   reduces) vs deviates — explicitly as DIAGNOSTIC evidence, NOT a declared fix.
4. C0f-2 reports `1−P_G`/`1−P_cpl` + the boundary/interior curl split + the mixed-control demonstration of the curl gate's
   k-bias; recomputed labels (NOT the dead C0e curl labels); raw curl NOT used as a classifier.
5. Exactly one C0f-3 verdict with falsifiable numeric `verdict_support`; WALL_WAS_CONFIG only if the ACCEPTED crawl
   reaches solver_converged + Linf≤1e-6 through τ=0.028 with bounded `‖Δx/Δτ‖`.
6. Faithful operators + frozen physics + export guard + SOLVER LOGIC UNTOUCHED (diff); depth continuation τ-only; NO
   trust-region/dogleg/LM, NO pseudo-arclength, NO gauge deflation implemented; chunk gates pass; report + machine JSON.

**Fail conditions:** keeping the C0b 2-iter / 0-backtrack overrides in C0f-0; judging convergence on a scaled/floored
surrogate; declaring WALL_WAS_CONFIG from a shallow τ gain, from passing 0.029 only, or from a C0f-1 diagnostic α rather
than an ACCEPTED crawl through τ=0.028; a qualitative (non-numeric) fold call; declaring WALL_WAS_CONFIG while `‖Δx/Δτ‖`
satisfies the fold rule; conflating b2c_linf_pass with solver_converged; re-using the dead C0e curl labels or the curl as
a classifier; implementing trust-region / pseudo-arclength / gauge deflation / any solver-logic change; changing a frozen
physics parameter; altering operators/frozen/export; masking NOT_MEASURED; raising the timeout cap.

## Out of scope (gated follow-ons)
C0g = trust-region/dogleg/LM (only if C0f returns GLOBALIZATION_INSUFFICIENT with good direction) OR pseudo-arclength
continuation (only if FOLD_DETECTED); minimal gauge deflation as linear-solve insurance; the production-solver question;
the high-resolution profile solve; the calibrated branch + calibrate-predict + `pathA_22` (all only after the crawl is
unblocked toward the physics target).

## Review (orchestrator, after Codex)
Fidelity agent: did C0f-0 ACTUALLY run with code defaults (not the 2-iter override) — show the config; is convergence
judged on the unmodified `patha_closed_branch_residual`; are the aids shown path-only with zero-ε final solve; is the
dψ/dτ a real numeric state finite-difference with the lane breakdown; is the predicted-vs-actual sweep genuine (real F
evals + real `F+αJδ` from the unscaled Jacobian); is `1−P_G` the dimensionless energy fraction (not the curl ratio); are
operators/frozen/export/SOLVER-LOGIC untouched (diff); no trust-region/pseudo-arclength/deflation implemented? Adversarial
agent: is WALL_WAS_CONFIG EARNED (an ACCEPTED crawl genuinely reaching solver_converged + Linf≤1e-6 through τ=0.028, no
fold masked by the numeric rule, no surrogate residual, b2c vs solver_converged not conflated), or over-optimistic
(declared off a shallow gain / a diagnostic α / a deeper wall or fold lurking)? Is FOLD_DETECTED / FOLD_RISK /
GLOBALIZATION_INSUFFICIENT, if returned, supported by the numeric dψ/dτ + predicted-vs-actual evidence? Then gate the next
step (continue the crawl to the physics target, or C0g).
