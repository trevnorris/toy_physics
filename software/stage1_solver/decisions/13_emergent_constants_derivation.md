# Decision 13 — B2c verdict is UNDETERMINED; the real next chunk is the EMERGENT-CONSTANTS derivation

**Date:** 2026-06-19
**Status:** DECIDED direction (user-driven, 2026-06-19). This is the **resume-here record after /compact** for the
Path-A build. Supersedes the "rigorous MISS" reading of B2c (decision-12 B2c STATUS block — now flagged superseded).
**Mechanism:** B2c 3-round build/review → two audits (`pathA_17` validity, `pathA_18` dimensional) → two independent
verification agents → user methodology call (derive the emergent constants before `m̂0²·S_port`).

---

## 0. STATUS / NEXT ACTION (resume here — 2026-06-19, pathA_21c EXECUTING)
**Where we are:** `pathA_19`/`pathA_20`/`pathA_20b` EXECUTED + REVIEWED + COMMITTED (§8/§9/§10). `pathA_21` EXECUTED +
reviewed (§11) — negatives HONEST, but P1 force was a RESTATEMENT + the P5 spec wasn't computable → spawned `pathA_21b`.
`pathA_21b` EXECUTED + 4-agent reviewed (full ledger recorded with the pathA_21c review) — **big wins:** G1 stationary BVP genuinely CLOSED + codeable; drain
velocity genuinely Gauss-solved (`r`-power emerges from the surface measure); sign honest; residuals G2/G4/G5/G6 honest;
branch-realization to-do = a finite **4-item list** (`R0` selector, `J`-value, kernel shape, brane zero-mode). **ONE
surviving overclaim (adversarial-caught):** the inter-defect force COEFFICIENT was still a heuristic product
`F=m·N·Q2·v1` dressed as a Π_cross integral → spawned `pathA_21c`.
**`pathA_21c` (force + SIGN from the GNLS momentum-BALANCE / Noether tensor) — directive READY + NOW EXECUTING.**
Design-reviewed SOUND-WITH-FIXES → all 9 fixes APPLIED (momentum-BALANCE law `∂_t g_i+∂_j Π_ij=f_i^body` since
`V_conf`/`Z(w)`/`J_ext` break translation invariance; pinned sign convention; stress-improvement ambiguity controlled;
δρ/quantum/V_conf/Maxwell derive-or-residualize) **+ the CALIBRATE-PREDICT reframe** (force STRUCTURE/`r`-power/SIGN =
target-blind PREDICTIONS to derive; overall NORMALIZATION = a labeled CALIBRATION KNOB, NOT "derived") → confirm-pass
**SOUND-AS-IS**. Run: `_scratch/pathA_21c_execute.log` (session 019ee287). Execute prompt: `_scratch/pathA_21c_execute_prompt.md`.
**CALIBRATE-PREDICT REFRAME (user reaffirmed 2026-06-19, now load-bearing for the whole frontier):** the throat-profile
branch data (`R0`/`J`/kernel) and the force normalization are **CALIBRATION KNOBS**, not things to derive ab initio. Cross
the "everything needs the profile solve" wall by **pick a branch family → calibrate to a trusted anchor (Newtonian `G` /
GR quadrupole `54/5`) → predict the held-out SURPLUS** (g−2, 5PN, multi-defect). The DERIVED-FORM GATE forbids only
calling a calibrated thing "derived," NOT calibration itself. ([[feedback-calibrate-predict-methodology]], updated.)
**`pathA_21c` DONE + reviewed — FORCE MILESTONE.** The inter-defect force is FINALLY a genuine derivation (the
pathA_21/21b shortcuts are superseded): Noether stress `Π_ij=m_GNLS ρ v_iv_j+δ_ij P+σ_Q,ij`; balance law
`∂_t g_i+∂_jΠ_ij=f_i^body` VERIFIED against the parent Euler identity (the sub-identities `dP=ρdh`, `∂_jσ_Q,ij=ρ∂_iQ`
machine-verified — adversarial-confirmed REAL); force = a REAL surface integral [convective −4/3 + Bernoulli +1/3 = −1] →
`F_12=−(m_GNLS N_∞,3 Q1Q2/4π r²) r̂` (inverse-square, bulk R⁻³, structure+power EMERGE from the Gauss substitution);
**SIGN = FORCE_ATTRACTIVE_DERIVED for like drains** (target-blind, no positivity smuggle; full sign =
`SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE`). Calibrate-predict ledger honest (3 predictions / 1 knob = normalization).
Minor carried harness-tightening flag: the angular `4/3`/`1/3` are asserted (hardcoded `(1+1/d)`,`1/d`), not ∮n_in_j dΩ
integrated — genuineness rests on the hand-verified prose; tighten one engine later ([[feedback-transliteration-fidelity-audit]]).
**So: the analog produces an inverse-square ATTRACTIVE force for like drains, derived target-blind from the superfluid
stress tensor.** NOT yet committed.
**PHASE-1 SOLVER RECONNAISSANCE (3 clean agents) — THE WHALE REFRAMED:** the throat-profile solve (option C) was audited.
The PDE OPERATORS are MODEL-FAITHFUL (real quintic EOS, covariant Laplacian, gauge coupling, localized Maxwell H=Z,
ψ-sourced current — coupled_branch.py/operators.py); the "engineering-smoke" label = parameter values + constitutive
closure (V_conf/Z(w)/R0 shapes), NOT the operators. The solver + calibrate-predict harness are **~70-80% BUILT + VALIDATED**
(Newton-Krylov; chunk 1a/1b/1c DONE; B2c prototype: knob→real-solve, extraction map, anchor+root-finder, held-out-surplus
skeleton w/ target-blind firewall). **The blocker is NUMERICAL CONDITIONING, not physics and not a from-scratch build:**
no quantum-core safeguard (√ρ→0 degenerates the matter Jacobian); the `k1∝r⁴/R0⁵` wall coupling blows up unclamped as
R0→0; uniform grid no grading; non-production preconditioner. LIVE EVIDENCE: B2c STALLS at `tau≈0.029` = a diagnosed
conditioning floor, NOT a physical edge. The R0 constitutive family is a posited placeholder (a calibration knob not yet
promoted to a frozen branch).
**NEXT ACTION (resume here) — the CONDITIONING SPIKE (the next attack on the whale; task #78):** (1) density floor /
log-ρ in the matter block + clamp the `k1∝1/R0⁵` blow-up (convergence aids that vanish at the solution); (2) Jacobi
row/col scaling of the bordered Jacobian; (3) short depth-homotopy continuation; + a DIAGNOSTIC on one moderately-deep
branch (D0, smallest Jacobian singular value, GMRES growth vs depth) that DECIDES whether scaling+homotopy suffices or a
production linear-solver swap (PETSc/multigrid) is needed. Directive READY v2 at
`directives/pathA_C0_conditioning_spike.md`: **design-review = SOUND-WITH-FIXES → all 11 fixes + 6 new-problem items
applied; Codex confirm-pass = SOUND-AS-IS (both logs under `_scratch/pathA_C0_directive_{review,confirmpass}.log`,
gitignored).** The directive now turns on ONE load-bearing rule — the **Single Arbiter Principle**: the ORIGINAL,
unmodified physical residual is the SOLE arbiter of both convergence AND solution-invariance; every aid (floor/`k1`
clamp/Jacobi scaling/variable change/depth homotopy) is PATH-only and must be residual-equivalent OR vanish (`ε→0`) at the
final solve. Key teeth added: deep-regime invariance (not a dormant shallow case); `log ρ`/`√ρ` forbidden unless it
preserves the COMPLEX matter block's phase/current/gauge lanes; depth continuation is `τ`-only (no `a`/`L`/`r_mouth`/
`r_exit`/`w_max`/boundary/constitutive — freeze violation); `SPIKE_SUFFICIENT` only on the ORIGINAL unscaled residual at
the B2c tolerance; `PRODUCTION_SOLVER_REQUIRED` only after C0-1..4 all active with persistent-failure evidence.
**pathA_C0 EXECUTION RAN 2026-06-19 (Codex) → verdict reported `PRODUCTION_SOLVER_REQUIRED` → ADVERSARIAL REVIEW (2
clean agents) = VERDICT NOT EARNED / SHORT-CIRCUITED. DO NOT ACCEPT IT. Outputs are untracked WIP (NOT committed):
`src/stage1_solver/patha_c0_conditioning_spike.py`, `tests/test_patha_c0_conditioning_spike.py`,
`reports/pathA_C0_conditioning_spike.md`, `runs/pathA_C0_conditioning_spike/…json`.**
*What is SOUND (reusable):* the Single-Arbiter machinery — original-residual gating, Jacobi neutrality (`R·J·C`, line
search on unscaled), faithful-operator safety (no operator/freeze/export edits; verified by `git diff`). *What is BROKEN
(why the verdict fails the gate):* (1) **no genuine depth crawl** — `prefer_existing_b2c_background_predictor=True`
cold-loads pre-existing B2c backgrounds for 7/8 τ; every "converged" row is an OLD B2c solution re-verified at **0 Newton
iters** (residuals match to all digits) → the spike solved ZERO new throats; (2) **below-floor attack is single-ε /
~single-step** — the aid loop `break`s on first failure (`patha_c0_…py:1052`) so the ε schedule `[0.08→0.02→0]` never
advances and there is no backtracking → the gate's "C0-1..4 all active, persistent failure" bar is UNMET (the face-saving
exit the directive forbids); (3) **tautological admissibility** — `residual_equality_max_abs=0.0` is HARDCODED (`:1200`),
`epsilon_independence=PASS` rests on hardcoded `final_aids_inactive=True` + only checks already-converged rows (the
dormant-case FAIL pattern); (4) **MISDIAGNOSIS** — the bordered Jacobian is near-singular (cond≈1e20) even at the
CONVERGED shallow τ=0.03 and conditioning IMPROVES (cond≈1e14) at the failed deeper τ → the wall is an intrinsic
near-null-space, most likely the unscaled dense mass/μ **border lanes**, a **gauge/zero mode**, or a **τ-FOLD/turning
point** (which would make a production solver the WRONG fix — pseudo-arclength continuation is the right tool, NOT
multigrid/mesh-grading); (5) **binding-constraint MISS** — `min_R0` INCREASES (0.748→0.806) as τ↓ and the k1 clamp goes
INACTIVE at the deepest τ, and `min_ρ≈7e-6` is ~flat → **τ-depth ≠ R0-depth**; the empty-core (√ρ→0, R0→0) regime the
program cares about was NEVER approached; (6) **σ_min untrustworthy** — a 1-norm LU *estimate* (true SVD gated out by
state size 1297 > 360), and it INVERTS (healthy looks more singular than broken).
**pathA_C0b EXECUTED 2026-06-20 (corrected) → verdict `DIAGNOSTIC_INCOMPLETE` (HONEST — the anti-fake gate refused a
substantive verdict off a bounded crawl) → FIDELITY-AGENT VERIFIED clean (no v2 sins; faithful ops untouched; the
decomposition is trustworthy — slicing direction proven by a planted-null-mode unit test). VERIFIED FINDING (this is real
progress):** the τ≈0.029 "wall" is **NOT the mass/μ border** (REFUTED my earlier hypothesis) and **NOT a fold** (REFUTED —
tracked triplets, overlaps ~0.997–0.9999, no complement crossing). It is a **PERSISTENT near-null SUBSPACE in the PHYSICS
FIELD BLOCK**: true dense-SVD gives ~5 tracked modes (σ_min≈1.4e-12 deep / 2.1e-15 shallow, then 4.7e-11, 2.6e-10, 5.9e-10,
2.8e-8), ALL with `v_min` energy = 1.0 in `field[0:5n]` and `u_min` = 1.0 in the PDE rows, present at BOTH τ; `cond_ratio`
(field-block/full) rises 0.097 (shallow) → 1.01 (deep). **Most likely identity (UNCONFIRMED — the gauge/symmetry projection
was `not_available`/unimplemented):** an unfixed U(1) global-phase (gauge) zero mode (± other symmetry modes) of the
stationary GPE+Maxwell BVP. **If so, the fix is gauge-fixing / null-space DEFLATION (a CONSTRAINT) — CHEAP, NOT a production
solver and NOT pseudo-arclength** → potentially unblocks the throat solve cheaply. C0b outputs (verified honest):
`src/.../patha_c0_conditioning_spike.py`, `tests/test_patha_c0_conditioning_spike.py`, `scripts/pathA_C0b_{crawl,diagnostics}.py`,
`reports/pathA_C0b_wall_diagnosis.md`, `runs/pathA_C0b_wall_diagnosis/…json`.
**pathA_C0c EXECUTED 2026-06-20 → verdict `MIXED` → FIDELITY-AGENT VERIFIED TRUSTWORTHY (all 8 checks real: dual
independent annihilation paths assembled-vs-JVP, recomputed SVD from saved matrices, real planted-generator test, honest
verdict logic, dense-σ MATCH retiring the stencil caveat, faithful ops untouched). VERIFIED FINDING (big clarification):**
the worst near-null mode (mode 0, σ≈7.4e-14 ≈ exact zero) **IS the global U(1) gauge PHASE zero mode** — CONFIRMED by BOTH
gates (annihilation `‖J·g_phase‖/(σ_max‖g‖)=1.49e-12` ≤1e-8 via assembled AND autodiff-JVP; overlap=1.0; energy 99.99% in
the `psi_imag` lane; equivariance `J·g_phase≈G·F(x)`=2.7e-17 at the non-converged τ). The OTHER 4 modes (σ 4.7e-11..2.8e-8)
are **MAXWELL A-FIELD lane (ar/aw, ~99%) near-null modes**, labeled `UNEXPLAINED_STIFFNESS` ONLY because the single crude
`maxwell_residual_gauge` probe (one fixed `λ`) doesn't span them (overlaps 0.11–0.77). **KEY CAVEAT (fidelity agent): those
4 are almost certainly NOT genuine stiffness — they're very likely Maxwell-sector GAUGE modes a proper `∇λ` basis would
span; do NOT call them stiffness / conclude "production solver" until tested against a FULL Maxwell gauge basis.** So the
whole τ≈0.029 wall may be ENTIRELY GAUGE (U(1) phase + Maxwell A-sector) ⇒ a CHEAP combined gauge-fix/deflation could
dissolve it. C0c outputs (verified, committed): `src/.../patha_c0_conditioning_spike.py` (+C0c additions),
`scripts/pathA_C0c_nullmode.py`, `tests/…`, `reports/pathA_C0c_nullmode_identification.md`, `runs/pathA_C0c_…json`.
**pathA_C0d EXECUTED 2026-06-20 → verdict `MIXED_GAUGE_PLUS_RESIDUAL` (3/4 "genuine stiffness") → FIDELITY-AGENT VERIFIED
the CODE faithful (term-by-term ops; the 4 modes freshly re-SVD'd from the saved C0b matrix, NO cold-load; ξ=1.0 verified
so the dropped 1/ξ factor is inert; negative control circular = weak but not wrong) BUT ADVERSARIAL + Codex-consult +
GLM-review judged the VERDICT OVERSTATED.** Why: C0d gated gauge-vs-stiffness on `P_G` + a hand-picked
weighted-divergence-RESIDUAL threshold (0.1), which is a gauge-FIXING-operator response, NOT the gauge-INVARIANT field
strength. Modes 1/3/4 are ~99.9% pure gradient (P_G≈0.999 = curl-free = GAUGE); only mode 2 (P_G=0.83, ar-dominated, ~17%
non-gradient) is a partial exception. The 0.1 gate is the SOLE discriminator (mode 3 clears it by only 0.022); mode 1
sits 99.7% in the penalty's OWN blind subspace (G_harm) yet is called the worst stiffness (the gate reads `grad` of a
near-zero, mesh-oscillatory divergence). C0d outputs committed: `src/.../patha_c0_conditioning_spike.py` (+C0d additions),
`scripts/pathA_C0d_maxwell_gauge.py`, `tests/…`, `reports/pathA_C0d_maxwell_gauge_identification.md`, `runs/pathA_C0d_…json`.
**3-REVIEWER CONSULT (Claude+Codex+GLM) → CONVERGED PLAN** (`_scratch/pathA_C0e_agreed_plan_for_GLM.md`: §6 Claude+Codex,
§8 GLM + Claude's C5): the decisive test is the per-mode gauge-INVARIANT **curl fraction `‖Z·F_rw‖/‖A‖`**, with **mode 2 as
the swing GATE** (curl≈0 ⇒ all 5 modes gauge ⇒ the wall is entirely cheap; curl O(1) ⇒ ONE transverse mode, NOT a
production solver). CONFIRMED vs code: the COUPLED gauge generator `(δψ=i(q/ħ)λψ, δa0=0, δA=∇λ)` is complete (no temporal
δa0 for the stationary ansatz; a0 elliptic via Gauss → no near-null modes); BCs do NOT restrict the gauge subspace
(deflation SAFE); the discrete de Rham complex is INEXACT (so `P_G` can't bound the curl — measure it directly); no second
continuous symmetry. **Claude C5 (carried; Codex vets at design-review):** the σ≈1e-11 vs k²~O(1–10) **9-order gap** ⇒
likely ODD-EVEN/CHECKERBOARD decoupling from the mismatched div/grad stencils (or a small Maxwell-row scaling), NOT smooth
modes ⇒ if so, a consistent/staggered STENCIL fix is a STRUCTURAL alternative to deflation that generalizes across τ.
Added by GLM: a Step-0 full-Newton-step ‖F‖ check (rules globalization in/out); adaptive re-SVD at each new stall; the
no-deflation reference exists ONLY at τ≈0.029.
**pathA_C0e EXECUTED 2026-06-20 → verdict `GAUGE_FRAMING_REFUTED` → FIDELITY-AGENT VERIFIED the CODE faithful BUT the
VERDICT was REJECTED (k-biased-metric FALSE NEGATIVE) by adversarial + Codex + Claude.** Why rejected: the curl gate
`‖Z·F_rw[v]‖/‖A[v]‖` is structurally `‖∂A‖/‖A‖` (a derivative-over-value ratio ≈ a WAVENUMBER, NOT a content fraction),
so it AMPLIFIES the ~0.1% high-k remnant of a 99.9%-gradient near-null mode → it mislabeled gauge modes 1/3/4 TRANSVERSE.
Three convergent tells: (1) C0e's own Jv budget shows every mode's near-null-ness is a curl-term≈penalty-term CANCELLATION
with λ=SMOOTH / `SMOOTH_K2` (= the smooth-gradient gauge story, per the operator identity curl-curl+(1/ξ)grad-div on ∇λ →
(1/ξ)∇(Lap λ), small for smooth λ); (2) curl/energy INVERSION (mode 4 has 5× more transverse energy than mode 3 yet 5×
less curl — impossible if curl tracked physics); (3) ~53% of each mode's curl norm sits on the 1-cell boundary ring
(closure leakage). **Honest gauge picture (back to the 3-reviewer view): modes 1/3/4 are GAUGE; only mode 2 (P_G=0.83,
17% transverse) is a genuine candidate. Unbiased discriminator = dimensionless `1−P_G`. Drop the raw curl gate as a
classifier.** C0e outputs (faithful code; verdict rejected — annotate on commit): `src/.../patha_c0_conditioning_spike.py`
(+C0e), `scripts/pathA_C0e_curl_gate.py`, `tests/…`, `reports/pathA_C0e_gauge_invariant_curl_gate.md`, `runs/pathA_C0e_…json`.
**⭐⭐ THE REFRAME (two findings that change everything; Claude+2 agents+Codex converged):**
**(A) The τ≈0.029 "wall" is NOT a proven wall — it's a CRIPPLED-SOLVER ARTIFACT.** The C0b run that DEFINED the stall
used `max_newton_iters_override=2` (code default = 18, `patha_closed_newton.py:81`), `max_tau_backtracks=0`, and
`depth_sequence=[0.03, 0.02899]` (ONE τ step of 0.001). Independently verified via `runs/pathA_C0b_wall_diagnosis/pathA_C0b_diagnostic.json`.
So the "wall" = from a converged τ=0.03, take one τ step, allow 2 Newton iters, never backtrack → "fail". B2c gate is
Linf≤1e-6; stall sits at Linf=5.32e-5 (~53× above) = "not converged in 2 iters", NOT "unreachable". FIVE rounds (C0–C0e)
analyzed a 2-iteration budget as if it were physics.
**(B) The near-null subspace is NOT what blocks the Newton step (C0e-0, fidelity-verified genuine):** the full step
overshoots ‖F‖ 5.08× while the near-null subspace contributes ~nothing (fraction 1.8e-21); removing it changes nothing;
linear solve clean (7.9e-13). ⇒ the proximate blocker is GLOBALIZATION/step-length, not gauge-conditioning (Codex caveat:
α=1 overshoot alone isn't fully diagnostic; if `Jδ=−F` is consistent the L2 merit slope is −‖F‖² at α=0 so SOME small α
must reduce ‖F‖ — C0b even found accepted α=1/128 — so the real question is whether the trust radius is large enough to
CONVERGE vs plateau). **Production solver OFF the table for unblocking the crawl; gauge deflation drops to cheap INSURANCE.**
**GLM REVIEWED the C0f plan → ENDORSED + 3 refinements folded in:** run-the-DEFAULT-config FIRST (the cheapest test;
C0b overrode every default); a cheap NUMERIC `dψ/dτ` FOLD DETECTOR (the one genuine residual risk — a turning point near
τ=0.029 that no globalization passes; pseudo-arclange doesn't exist in the codebase yet); and the residual is SMOOTH (EOS
`(5K/4)ρ⁴`, ρ=|ψ|², NO √ρ in `coupled_pde_residual` — `physics.py:46-47`) so the α=1 overshoot is pure large-step
nonlinearity and SOME small α MUST reduce ‖F‖ (⇒ Q2 nonsmoothness worry eliminated). GLM also retracted its own earlier
"use curl fraction as a gate" rec. Plan: `_scratch/pathA_C0f_agreed_plan_for_GLM.md`; consult `_scratch/pathA_C0f_consult.log`.
**pathA_C0f EXECUTED 2026-06-20 → verdict `DIAGNOSTIC_INCOMPLETE` (HONEST), key data is ENCOURAGING + RE-VERIFIED:**
(1) **GAUGE QUESTION SETTLED** — the UNBIASED `1−P_G` re-confirm: modes 1/3/4 are GAUGE (`1−P_G_A`≈7e-4/2e-4/1e-3), mode 2
the LONE transverse candidate (0.17); the MIXED control reproduced `ε²/(1+ε²)=9.9e-3` EXACTLY (k-independent) while raw
curl grew with k ⇒ DEFINITIVE proof the C0e curl gate was k-biased and its REFUTED verdict WRONG. (2) **The crawl is
GENUINE** (Claude-verified: `prefer_existing_b2c_background_predictor=false`; each τ warm-starts from the prior converged
state — NO C0-v2 cold-load) and with the DEFAULT config CONVERGES CLEANLY (0 backtracks, Linf~1e-12) through τ=0.029125 —
PAST where the crippled C0b run "stalled" ⇒ the "wall" is dissolving under a non-crippled solver. (3) **WRINKLE:** a SINGLE
τ-step (0.0290625) hit the **600s PER-ATTEMPT timeout** (SIGALRM; not the cumulative crawl) before finishing → no stalled
state saved → merit sweep NOT_MEASURED. (4) **YELLOW FLAG:** `‖dψ/dτ‖` growing 3.2× & accelerating, **R0-dominated**, as
τ→0.029 (below the 10× gate ⇒ FOLD_RISK not FOLD_DETECTED). The step that timed out was the SMALLEST one — ambiguous
between "slow/expensive solve" and "emerging fold." C0f outputs: `src/.../patha_c0f_globalization_probe.py`,
`scripts/pathA_C0f_globalization_probe.py`, `tests/test_patha_c0f_globalization_probe.py`, `reports/pathA_C0f_globalization_probe.md`,
`runs/pathA_C0f_globalization_probe/…json`.
**pathA_C0f2 EXECUTED 2026-06-20 (chunked/timed, user-authorized >600s) → fidelity + adversarial reviewed (faithful;
numbers trustworthy).** VERIFIED: per-τ cost = **96–97% Jacobian assembly** (colored JVP build → ~100× sparse-assembly
target); the crawl converges CLEANLY to **τ=0.029125** (GENUINE warm-start, NO cold-load) then STALLS at 0.0290625 —
**plain Newton+Armijo is effectively STUCK** (merit sweep: best step 0.05% reduction; α=1 OVERSHOOTS 13 orders, 4× worse
than start). The J⁻¹-amplified near-null direction is **NON-gauge (mode 2 transverse)** ⇒ **gauge-only fix is DEAD**;
**production solver OFF** (linear solve clean to 1e-12 — a GLOBALIZATION problem, not linear-solve). **FOLD-ONSET LIVE
~30%** (the backtrack-exhaustion fold-gate FIRED; growth 3.2× R0-dominated; reachable-Linf "improvement" is turning-point
creep, |Δτ|→0). Timing extrapolation (~1.67h to τ=0.02) is MOOT until the stall is cured.
**⭐⭐ CONSOLIDATED FINDINGS + RESEARCH BRIEF (read FIRST on resume):
`reports/pathA_throat_solver_findings_and_research_brief.md`** — self-contained: the PRECISE numerical problem
(gauged-GPE + Maxwell continuation BVP; near-singular J cond≈1e13; near-null = phase gauge + A gradient modes + mode-2
transverse), the full C0→C0f2 chain, the VERIFIED diagnosis, the OPEN fold-vs-ill-conditioning question + the cheap
det(J)/σ_min(τ)/tangent probe that settles it, the C0g candidates, and the literature search terms. C0f2 outputs:
`src/.../patha_c0f2_timing_rerun.py`, `scripts/pathA_C0f2_timing_rerun.py`, `tests/…`, `reports/pathA_C0f2_timing_rerun.md`,
`runs/pathA_C0f2_timing_rerun/…json`.
**⭐ POST-/COMPACT PLAN — RESEARCH + DESIGN + C0g-diag EXECUTION ALL DONE (2026-06-20). Resume at the (4) C0g BUILD block below.**
(1) ✅ **ONLINE LITERATURE RESEARCH DONE** (5 parallel agents) → `reports/pathA_throat_solver_literature_synthesis.md`
(the convergent strategy: gauge-fix path-only FIRST → branch fold-vs-conditioning; PTC/LM + Sobolev preconditioner for
conditioning; pseudo-arclength for a fold; the NEW sonic/horizon-critical-point hypothesis). (2) ✅ **DESIGN LOCKED**
through the Claude+Codex consult (§5 of the synthesis — code-grounded corrections: gauge handling PATH-ONLY, scipy LM a
capped probe, no det sign-flip observable from one-sided states, 253 colors = deterministic radius-3 coloring not a bug)
+ the **GLM 5.2 review (§6 — the AUTHORITATIVE diagnostic-battery spec):** a premise-validating **Step 0 framing check at
the ACTUAL stall τ=0.0290625** (the mode-2-drives-the-stall premise was only verified at the OLD τ=0.02899 — MUST
re-check; if near-null fraction ~0, the whole battery premise collapses), **wᵀF_τ** as the cheap decisive fold test
(bordered cond(Jb) demoted to confirmatory), **σ_min²(τ) linear fit** (fold√ vs bifurcation; predicts τ_fold from
one-sided data), scipy-LM **+ branch-overlap gate** (both-sides signal), and the closed sonic question (ψ-form operator
is ALWAYS elliptic ⇒ NO principal-symbol type change; Mach = physical CONTEXT, the Jacobian test = the DISCRIMINANT).
**(3) ✅ C0g-diag EXECUTED + DUAL-REVIEWED (2026-06-20, Codex two gated passes).** Steps 0–3 (premise + primary fold
tests) → fidelity-review `FAITHFUL_WITH_MINOR_NOTES`; Steps 4–7 + Step-8 verdict → adversarial-review (both: HONEST,
faithful, scope-clean, NO gaming). Report `reports/pathA_C0g_diag_fold_vs_conditioning.md`; JSON
`runs/pathA_C0g_diag_fold_vs_conditioning/…`; scripts `pathA_C0g_{step0_premise,state_svd_ftau,steps4_6_7,step5_scipy_probe,aggregate_*}.py`;
C0g analysis code in `patha_c0_conditioning_spike.py` (+1589/−0; operators/physics/frozen/export UNTOUCHED).
**⭐⭐ VERDICT — the white-whale stall is a SIMPLE FOLD (turning point) in the throat-radius (`r0`) continuation at
τ_fold≈0.0291233; NOT conditioning, NOT sonic.** Machine verdict = `MIXED/INCONCLUSIVE` (literal contract) but ADJUDICATED
**FOLD-LEANING** (`primary_fold_support=True, primary_conditioning_support=False`, NO disagreements, NO gray items). FIVE
agreeing lines: Step-0 premise HOLDS (f_nn≈1.0 at the REAL stall τ=0.0290625 — GLM's re-validation; the step rides the
`r0` near-null mode, gauge fraction 5e-27 ⇒ **C0c/d/e GAUGE framing REFUTED as the stall DRIVER**); Step-2 cosθ 0.51–0.64
stable (wᵀF_τ≠0 ⇒ simple fold not bifurcation); Step-3 σ_min² linear R²=0.9994 → τ_fold; Step-6 bordered-cond TREND
(cond(Jb) flat@570, cond(J·Q_perp)→2.3e5, ratio→404 accelerating); conditioning refuted. Step-4 Mach: zero current
(ψ≈real, static empty-flow background) ⇒ sonic NOT TESTABLE (→ NON_SONIC). **Why not clean FOLD_CONFIRMED: ONLY the
Step-6 absolute >1e10 bar blocked it, and that bar is SELF-DEFEATING (adversarial-confirmed) — σ_min∝√(τ−τ_fold) so
cond>1e10 needs sampling ~1e-15 from τ_fold but Newton stalls ~1.6e-6 short BY THE DEFINITION of a fold; the directive
even said "no absolute O(1) thresholds" then set one. The bordered-cond TREND, not the absolute, is the discriminant.**
The 511-dim gauge near-null subspace is a SEPARATE, independently-real conditioning issue (NOT the stall driver) to fix
regardless.
**(4) NEXT ACTION (resume here) — the C0g BUILD (both reviews converge):** (a) **gauge-fix FIRST** — PATH-ONLY (far-field
ψ phase pin + A-sector gauge handling in solver/preconditioner coords, NEVER the frozen residual); kills the 511-dim gauge
null space AND lets the crawl sample closer to τ_fold. (b) re-run cheap Steps 1–3 + 6 on the gauge-fixed crawl nearer
τ_fold → if σ_min² stays linear-monotone and cond(J·Q_perp) crosses ~1e8–1e10 with cond(Jb) flat, promote MIXED →
FOLD_CONFIRMED. (c) **gauge-fixed PSEUDO-ARCLENGTH continuation** to round the fold. **SKIP the LM/PTC conditioning detour**
(evidence disfavors — LM damps but never rounds a fold). + analytic/sparse Jacobian assembly unconditionally (after the
color audit; 253 colors = deterministic radius-3, NOT a bug). → **C0g-BUILD directive READY**
(`directives/pathA_C0g_build_gaugefix_then_pseudoarclength.md`; design-review SOUND-WITH-FIXES → 8 fixes → confirm-pass
STILL-NEEDS(3) → 3 fixes → re-confirm SOUND-AS-IS; committed) — its B-1/B-2/B-3/B-4 structure: gauge-fix path-only →
re-confirm the fold via RELATIVE trend gates (the absolute cond>1e10 bar dropped as self-defeating) → gauge-fixed
pseudo-arclength (only if FOLD_CONFIRMED) → analytic/sparse assembly. **NEXT = the USER EXECUTION GATE → Codex executes**
([[feedback-directive-design-review]]) → then constitutive family → calibrated branch → calibrate-predict (R0/J/W → anchor
→ SURPLUS) → `pathA_22`.
**LESSON (candidate memory): a "gauge-invariant" metric built as derivative/value is dimensionally a wavenumber and
amplifies high-k remnants — use the dimensionless energy fraction; and ALWAYS check the solver CONFIG (iter budget,
backtracking) before diagnosing a "wall" as physics.**
**⏱ STANDING FLAG — `timeout 600` cap (RAISE WITH ME, DON'T DECIDE ALONE — user asked to be flagged 2026-06-19):** the cap
is currently a FORCING FUNCTION and is NOT binding (C0/C0b respect it by SPLITTING into ≤600s scripts; a timeout degrades
to `NOT_MEASURED`/`DIAGNOSTIC_INCOMPLETE`, never a fake). It WILL legitimately bind at the **real high-resolution profile
solve** (post-diagnosis, CPU-only/GPU-off → genuinely O(grid) cost). WHEN that step arrives, BEFORE touching the cap walk
the ladder IN ORDER: (1) coarse→fine continuation with RESUMABLE CHECKPOINT CHUNKS, each ≤600s (keep the cap per-chunk);
(2) the better solver (the production-solver step, IF the C0b verdict calls for one); (3) the smaller modal/port
formulation (fewer DOF/solve). Raising the per-chunk cap is the LAST rung and a USER-LEVEL call ([[feedback-script-timeout-policy]],
[[feedback-never-alter-calibrated-process]]). **TRIGGER to bring it to the user:** a SINGLE irreducible Newton/linear-solve
at the resolution we actually need cannot fit in 600s EVEN AFTER chunking + a better solver + the modal/port form. Until
that exact condition, keep 600s. (Note placed here because it surfaces at the high-res-solve juncture; do not delete on
compact.)
**Discipline reminder:** Codex derives/codes + applies fixes, Claude reviews; orchestrator owns directives/decisions.
The DERIVED-FORM GATE binds (no hand-inserted field/`r`-power, no convention sign, no `x==x` posing as a check, no
restatement to fake BVP closure). VALID expected outcomes: a derived far-field force with interior factors flagged, an
honest `ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL`, and G6/`α_J`/`ħ`/`2π` remaining residuals. Commit pathA_21 + its 21b
corrections together so the ledger lands honest (commit only when the user asks).

---

## 1. The headline: B2c does NOT give a validated hit-or-miss
B2c computed, on the self-consistent Path-A background, that the model's `P0` is **6.7–9.6 orders of magnitude below**
the GR quadrupole target `54/5 = 10.8` (measured at every converged τ: `P0=2.795e-9` at τ=1 … `1.22e-6` at
τ=0.029). The numbers are real. **But the comparison is `R_norm = m̂0²·S_port·P0 − 54Gc_s⁵/(5a⁵c⁵)`, and B2c pinned
`m̂0²·S_port = 1`.** Two audits + two independent verification agents established:

- **NOT a dimensional/units bug.** With units restored, the dimensional gaps the audit found are **order-neutral**
  (they equal `×1` under the natural-unit pins `a=c_s=ħ=m=1`); restoring them shifts NOTHING numerically. (This
  *corrected* `pathA_18`, which had overstated "8 inconsistencies" → really 2, and wrongly assumed a 4-spatial `G`.)
- **NOT a validated "missing physics" miss either.** `m̂0²·S_port` is **dimensionful** — it is the conversion factor
  `T_target/[P0]`, NOT a free dimensionless number (both verification agents, independently). Pinning it to `1` is a
  **scale choice** that *forces* the `D0→0` knife-edge. The pre-reg bakes this in: `T := 54Gc_s⁵/(5a⁵c⁵)/(m̂0²S_port)`,
  `D0 = N0/T` (`docs/pathA_preregistration.md` §253). So `m̂0²·S_port` directly sets whether calibration lands at a
  knife-edge (miss) or naturally (match).
- **The ~9-order "miss" is exactly the magnitude `m̂0²·S_port` would supply** (`10.8/P0 = 3.86e9` at τ=1). Whether the
  *derived* `m̂0²·S_port` is ≈1 (→ real miss) or carries large factors (→ match) is **plausible either way** (the
  healing length `a` is microscopic while the GR comparison is at the orbital/radiation scale, so a large
  `(scale-ratio)^n` is on the table) and **cannot be settled by dimensional analysis** — only by deriving it.

**Verdict: UNDETERMINED, hinging on the un-derived dimensionful normalization `m̂0²·S_port`.** B2c is NOT committed as
a miss. The B2c machinery (assembly/dual-engine/verdict-logic/warm-start/negative-controls) was fidelity-clean and is
retained; only the *interpretation* changed.

## 2. The foundational reframe (user-driven): the constants are EMERGENT, derive them first
`c`, `c_s`, `G` are NOT fundamental inputs in this model — they are **outputs of the 4D superfluid PDE**, and their
familiar "3+1" dimensions are the **brane (observed) dimensions**, while the **bulk PDE works with 4D-dimensioned
quantities**. The bulk→brane reduction (integrating out the transverse `w`-direction, width ~`a`) is where the
scale/dimension factors live — and `m̂0²·S_port` is almost certainly a piece of *that same bulk→brane normalization*.
So **deriving the emergent constants is logically PRIOR to `m̂0²·S_port`**; once the constants + reduction are in
closed form, `m̂0²·S_port` should largely fall out.

What the verification already established about the constants:
- **`c_s` (sound speed):** emergent from the EOS, `c_s² = 5Kρ⁴/m`; `[c_s]=L/T`; **value drifts with the background
  `ρ`** (not constant).
- **`G`:** emergent from defect back-reaction. Bulk 4-spatial Laplacian → `1/r²` (`[G_bulk]=L⁴T⁻²M⁻¹`); but the
  observed sector is **brane-projected to 3 spatial dims** → `1/r`, standard Burke-Thorne/Peters, `[G]=L³T⁻²M⁻¹`
  (confirmed from the parent action: defects live on the 3D brane, `w` transverse). **User hypothesis (likely right,
  must DERIVE): `G` reduces to the standard 3D form after the brane reduction, with the transverse-scale (`a`)
  factors being the content of the emergent `G`.**
- The GR target `54Gc_s⁵/(5a⁵c⁵)` is NOT dimensionless on its own (it carries `M⁻¹L⁻²T⁻²` with 3D `G`); it is used as
  the Peters **pure number `54/5`** with constants pinned to 1, and `m̂0²·S_port` silently absorbs the dimension.

## 3. The velocity structure (refined with the user, 2026-06-19): THREE scales + the standing-wave bridge
The model has **three distinct velocities, all `L·T⁻¹`** (so dimensional analysis alone cannot separate them — only the
dynamics can). The earlier "c vs c_s = wave-vs-terminal-velocity" framing is superseded by this richer picture:
- **`v_b` — background condensate flow velocity** `v_b=(ħ/m)∇θ` (drift of the medium). This is the **gravitational-sector
  variable**: gradients of `v_b` and `ρ` ARE the analog field; `v_b=c_s` is the acoustic horizon. Varies in space, not
  constant. **`ρ` is ALSO not constant and is coupled to `v_b`** via stationary continuity `∇·(ρv_b)=0` + quantum-
  Bernoulli `½mv_b²+μ(ρ)+V+Q=const` (flow speeds up → `ρ` drops → `c_s` drops). So `ρ, v_b, c_s` are ONE coupled
  profile, not independent dials; genuine constants are `K, m, ħ, a`, the Bernoulli head, the mass flux, and the
  asymptotic `ρ₀/c_{s,0}` (the `c_s=1` pin = `c_{s,0}`). (Its full `G` role → `pathA_20`.)
- **`c_s` — bulk sound speed** (density/phonon waves through the medium). EOS-set, `∝ρ²`.
- **`c_γ` — photon/gauge-wave speed** (massless gauge excitation on the brane); the brane light cone.

**The standing-wave bridge (the key insight).** A massive particle (throat) is a **standing wave of the photon/gauge
field** — two counter-propagating `c_γ`-waves. Consequences, to be derived: (a) rest internal oscillation = Compton
frequency `ω₀=m*c_γ²/ħ`, so `E_rest=m*c_γ²` is trapped-wave energy; (b) driving the envelope at `v` Doppler-shifts the
components and slows the internal clock as `ω₀/γ`, **freezing at `v→c_γ`** (= relativistic time dilation, the user's
"oscillations stop at light speed"); (c) the envelope cannot outrun its constituent waves, so the **terminal-velocity
ceiling is `c=c_γ`** — `c` is the terminal velocity AND the photon speed, unified because matter is standing light.
This makes the old `c≠c_s` hypothesis and the existing `c=c_s` result COMPATIBLE: `c=c_γ` is forced;

**The one surviving open number: `c_γ/c_s`.** Particle ceiling = photon speed is forced; whether the photon (brane
gauge) speed equals the bulk sound speed depends on whether gauge + density share the acoustic metric, or the
localization profile `Z(w)`/width `a` rescales the brane gauge cone (the EM sector already shows `μ₀^eff=μ₀/Z_int`).
**`c_γ/c_s` (closed form) IS the `(c/c_s)³` tail factor `R_tail`** (=1 iff they coincide). Let the math decide.

**Is DEFECT mass emergent? (user idea, 2026-06-19; sharpened by Codex directive-review into a falsifiable test.)** The
idea: a defect is a SINK pulling superfluid inward, and that inflow strength is the defect's mass; the same inflow
shared between two defects is their attraction, so **`m_defect` and `G` are two faces of one quantity** (→ `pathA_21`).
**Two masses must NOT be conflated:** `m_GNLS` (the constituent mass in the parent action's kinetic/EOS/current terms)
is part of the EXACT action and stays `[M]` unless an action-level rewrite keeps every term homogeneous; only the
DEFECT mass `m_defect` (a throat branch property) might be emergent. The ontology says defect mass ↔ drainage/VOLUME
deficit (volume flux), charge ↔ vorticity flux — so do not assume mass = number flux; frame-tag number flux
`J=T⁻¹`, volumetric `Q_vol=ρ⁻¹J`, mass flux. Honest relation: `m_defect=α_J ħ J/c_γ²` (`α_J` dimensionless/branch
data); `m_defect=J` only after `ħ=c=1` (→ using the not-yet-derived `c` to drop `M` is CIRCULAR; de Broglie route is a
note to `pathA_20`). **Falsifiable:** `{L,T}` (mass emergent) is allowed ONLY if the whole action rewrites homogeneously
AND a boundary/source/Noether/Hamiltonian derivation ties `m_defect` to the inflow; else RETAIN `{L,T,M}` and record
`J` as a conserved rate — a valid negative result. **`a`-pin:** `a` is a mouth-radius collective MOMENT (not a
fundamental coordinate), particle-dependent + deformation-fragile; the conserved invariant is the flux `J` (Gauss, in a
no-leakage region — but the throat bottom may be open, so no-net-accretion must be a derived BC or a logged gap). Assess
re-pinning `a→J`, but do NOT claim it changes `m̂0²·S_port` (may be dimensionally neutral) — that value is `pathA_22`.

## 4. NEXT CHUNK — the emergent-constants derivation (the agreed next step)
**Execution split (user call, 2026-06-19 — foundation first):** the derivation is run as FOUR gated sub-directives, not
one: **`pathA_19` (FOUNDATION** — base set / mass fork M-vs-`J` / flux-invariant re-pin of `a` / pin over-determination
/ dictionary / paper-prose reconciliation; this is split out first because it can change the base dimensional system
everything else sits on; ✅ DONE 2026-06-19, see §8) → **`pathA_20`** (`c_s` + `c`, the velocity structure; FINALIZED
against the pathA_19 base — now also carries **S2b** the throat as a transonic/choked drain → the flux law
`J_crit(ρ₀, geometry)` set by the sonic point = acoustic horizon, replacing the "constant `J`" assumption, and **S3**
the `ħ`-provenance fork [`ħ=m c_s0 a/√2` ⇒ `ħ` may be EMERGENT, not fundamental — decides whether `M` ever collapses])
→ **`pathA_21`** (`directives/pathA_21_emergent_G_mass_bridge.md`: `G` from defect back-reaction + 4D→3D reduction;
the **mass-bridge** `m_defect=α_J ħ J/c²` [= `E=mc²` with the inflow as the energy ≡ the standing-wave `ħω₀`] DERIVED +
the **M-collapse re-test** of the pathA_19 negative + the equivalence-principle check + the `m↔G` unification) →
**`pathA_22`** (natural-unit→physical scale map → derive `m̂0²·S_port` → re-run B2c → real verdict). Each gated by the
user. Conceptual basis for S2b/S3/the bridge co-developed with the user 2026-06-19 (see §3 + the directives).

A full **reasoned + dimensional derivation** of the model's constants, from the parent action outward (NOT assumed —
derived; much already exists scattered across `research/pde_ledger/paper/parts/part01_parent_geometry.tex`, the
sound-speed ledger, and `research/4d_2_5pn/paper/4d_2_5pn.tex` — so partly consolidate + make rigorous + dimension-
check, partly fill gaps):
1. **`c_s`** from the EOS (`c_s²=5Kρ⁴/m`); its dimension + state-dependence.
2. **`c`** — what it physically IS (test the `c≠c_s` / terminal-velocity hypothesis), its definition, and its
   relation to `c_s` (the `c/c_s` ratio).
3. **`G`** — DERIVE from the defect back-reaction + the 4D→3D brane reduction: closed-form in superfluid params
   (`K,a,ρ,ħ,m`), bulk vs brane dimensions, the transverse-scale factors; TEST the "reduces to 3D `L³T⁻²M⁻¹`"
   hypothesis explicitly.
4. **The natural-unit → physical-scale map:** what the pins `a=c_s=ħ=m=1` actually stand in for — which determines
   whether `m̂0²·S_port` is genuinely ≈1 or carries orders.
Run everything through the committed harness `software/stage1_solver/src/stage1_solver/dimensional_check.py`
(machine-verify); produce a consolidated "constants & dimensions" reference doc (4D-bulk vs 3D-brane frame tags +
reduction maps). Discipline: Codex derives/codes, Claude reviews; dimensional-check + fidelity + adversarial.

**THEN:** derive `m̂0²·S_port` (should largely follow from #3/#4) → re-do the B2c `R_norm(τ)=0` calibration with the
derived normalization → real hit-or-miss verdict.

## 5. GATE-A freeze implication (methodology decision, do NOT do unilaterally)
`m̂0²·S_port=1` is currently part of the **GATE-A freeze** (decision-11 §2/§3, hash `ed3585…`). Deriving it means
**un-pinning a frozen convention → a documented freeze amendment (new hash)**. This is a methodology call to settle
with the user when the constants derivation is in hand (it may show the pin was fine, or that it must be replaced by
the derived value).

## 6. Standing step locked in (user request)
**Dimensional-consistency check (units restored, SymPy harness) BEFORE trusting any number** is now a STANDING gate:
`docs/pathA_preregistration.md` "STANDING VERIFICATION STEPS"; memory `[[feedback-dimensional-consistency-check]]`;
harness `dimensional_check.py`. Pairs with `[[feedback-transliteration-fidelity-audit]]`.

## 7. Audit trail (distilled; raw outputs regenerable, in gitignored `_scratch/`)
- `pathA_17` (validity): the "miss" = the pinned `m̂0²·S_port`; `54/5`←Burke-Thorne is correctly derived; no dropped
  4π/Jacobian. Verdict: ARTIFACT for the full claim, REAL only inside the pinned convention.
- `pathA_18` (dimensional, units restored): built the reusable harness; flagged dimensional gaps — BUT overstated
  (8→2 distinct; 4-spatial-`G` wrong). Honest caveat: dimensional analysis fixes no decimal order (factors =1 in
  natural units).
- Verification agent A (adversarial): 8→2 real issues, both order-neutral natural-units conventions, not code bugs;
  "if anything STRENGTHENS the real-missing-physics reading" (i.e. dimensions don't rescue the miss).
- Verification agent B (independent re-derivation): `[G]=L³T⁻²M⁻¹` effective-3D (brane projection, from the parent
  action); `m̂0²·S_port` is the dimensionful conversion factor `T_target/[P0]`, not dimensionless; `D0=K−B0−Z0` is a
  legal subtraction (`K`=modal stiffness eigenvalue) within the frozen normalization conventions.

## 8. pathA_19 (FOUNDATION) execution + review ledger — 2026-06-19
**Run:** Codex gpt-5.5 @ xhigh (`_scratch/pathA_19_execute.log`, 162k tokens). Deliverables: extended harness group in
`src/stage1_solver/dimensional_check.py` (+509, side-by-side), `reports/pathA_19_dimensional_foundation.md` (F5 ref
doc), `tools/pathA_19_foundation_dimensional_crosscheck.wl` (dual-engine). Python 17/17 algebraic checks pass;
Mathematica 20 checks PASS; full suite 92 passed; `git diff --check` clean. Acceptance = **PASS_WITH_NAMED_RESIDUALS**
(exit-0 NOT treated as acceptance).

**F1 mass-fork verdict — RETAIN `{L,T,M}` (honest negative result):** `m_GNLS` stays an explicit parent-action mass
`[M]` (appears in kinetic operator, current, Madelung velocity, Euler eq, sound-speed law). `m_defect` emergence is
**NOT derived** — the ontology gives drainage/volume-deficit *scaling* (`brane_bulk_ontology.tex:1267-1297`) but no
boundary source / Noether charge / Hamiltonian energy theorem tying `m_defect` to the inflow rate. `ħJ/c_γ²=M` is a
*dimensional conversion only*, explicitly NOT a derivation. So mass-as-inflow is rejected **for this gate** and carried
to pathA_21 as `INFLOW_MASS_SOURCE_MISSING` (BLOCKS_MASS_EMERGENCE) with a concrete derivation target.

**F2 flux + a-pin:** `[J_number]=T⁻¹` in BOTH 4D-bulk (closed 3-surface) and brane (2-surface) frames; `Q_vol`=`L⁴T⁻¹`
(bulk)/`L³T⁻¹` (brane); `m_GNLS J`=`MT⁻¹`. Gauss shape-independence holds ONLY with no enclosed source/leakage — but
projection makes `S_leak` and the throat bottom is open/closed/connected-undecided → `NO_NET_ACCRETION_BC_UNDERIVED`
carried. `a` confirmed a mouth-radius collective MOMENT (`a0=R0(0)`, `a(t)`=mouth average), not fundamental →
`A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT`; use `J` as the invariant scale-map label.

**F3 pins / healing / dictionary:** 4 pins (`a=c_s=ħ=m=1`) on 3 base dims ⟹ one null relation
**`a=ħ/(m_GNLS·c_s0)`**. GNLS core balance ⟹ **`ξ_h=√2·ħ/(m_GNLS·c_s0)`** with `h0=(5K/4)ρ0⁴=m_GNLS c_s0²/4`; so if
`a`≡healing core then `a=ξ_h/√2` (convention/branch factor). Independent set = {`ħ`, `m_GNLS`, `K`, chosen `ρ0`};
derived = {`c_s0`, `ξ_h`, `a`(if core-identified), `m_defect`(blocked)}. 4D dictionary homogeneous; the 3D GR target is
a downstream conversion problem, not a base-system change.

**F4 paper-prose reconciliation:** `part01`/`pde.tex` 4D action/EOS/current/projection AGREE with the harness
dictionary; `em_fields.tex:1717-1786` flagged **WRONG-3D-CONVENTION** (ρ₀ as kg m⁻³, pressure/enthalpy-per-mass,
`V=πa²L` throat volume) — legacy 3D/SI prose, not the 4D number-density dictionary.

**Three R_norm dimensional flags (carried, NOT repaired here):** formal 4D R_norm target is NOT dimensionless (needs
`L T² M`); observed 3D GR target needs `L² T² M`; a TRUE `{L,T}` base FAILS the R_norm gate (needs `L² T²`). These
corroborate the B2c-undetermined verdict (§1) and belong to pathA_21/pathA_22 — pathA_18 behavior preserved.

**Review:** two clean transliteration-fidelity agents (Python module; Mathematica `.wl`) → both **FIDELITY-CLEAN**;
`.wl` confirmed INDEPENDENT (no `Import`/`Get` of Python results; own representation); no tautology/can't-fail gate (the
Maxwell `c²` factor and the `{L,T}`-gate rejection both demonstrably *can* fail); side-by-side scoping CLEAN. Two
prose-only (non-machine-checked) claims — the `√2` in `ξ_h` and `h0=(5K/4)ρ0⁴=m c_s0²/4` — hand-verified correct. One
harmless dead helper (`_dim_dict`).

## 9. pathA_20 (c_s + c velocity structure) execution + review ledger — 2026-06-19
**Run:** Codex gpt-5.5 @ xhigh. Deliverables: harness group in `dimensional_check.py` (+478, side-by-side,
`--patha20-velocity`), `reports/pathA_20_velocity_constants.md`, `tools/pathA_20_velocity_constants_crosscheck.wl`.
Python 21 dim + 5 alg checks pass; Mathematica PASS; full suite 92 passed. Acceptance = `PASS_WITH_NAMED_RESIDUALS`.

**Three `UNDETERMINED` verdicts — ALL VERIFIED HONEST WALLS (adversarial premature-punt audit), not punts:**
- `C_GAMMA_RATIO_UNDERDETERMINED`: the parent gauge action gives the photon the **bulk Minkowski cone** (`Z(w)` is an
  overall prefactor on `F_{MN}F^{MN}` → cancels from the principal symbol; it renormalizes the coupling `μ0_eff`/`q_eff`,
  NOT the cone), while `c_s` is the **emergent acoustic cone** — and the sources NEVER relate them (`em_fields.tex` only
  ASSERTS `c=c_s` weak-field by fiat). Carried: `λ_γ=c_γ/c_s`, tail `(c/c_s)³=λ_γ³`. [Caveat: the report's stated reason
  over-credits "`Z(w)` unsolved"; the real reason is structural (gauge-bulk-cone vs emergent-acoustic-cone) — sharpen
  this when pathA_21 carries it forward.]
- `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`: the model-faithful throat profile (with `Q`, `V_conf`, throat
  geometry `R0`, bottom-topology BC) is explicitly unsolved (`pde.tex:2845-2849`). Codex DID derive the conditional
  ideal-Euler nozzle (`c_*/c_s0=1/√3`, `ρ_*/ρ0=3^(-1/4)`, `J_crit/(ρ0 c_s0 A_*)=3^(-3/4)`, independently verified
  correct for the `Kρ⁵` EOS) and correctly REFUSED to promote it (dropping `Q` fails at the throat). [Caveat: report
  prose overstates "ODE cannot be closed" — the classical one WAS closed; the wall is the `Q`/`V_conf`/geometry/BC profile.]
- `HBAR_PROVENANCE_UNDETERMINED`: NO `ħ`-free substrate relation anywhere in the 4 cited papers (grep-confirmed);
  circulation `Γ∈ℤ` is an IMPOSED classical input label, not a derived substrate quantum. Refusing emergence = correct
  anti-tautology; `UNDETERMINED` over `FUNDAMENTAL` is the honest conservative call.

**Clean derivations:** `[c_s]=[c_γ]=L/T`, `c_s∝ρ²` (`d ln c_s/d ln ρ=2`), `[J]=T⁻¹` (bulk+brane). **Standing-wave
`c=c_γ` GENUINELY NON-CIRCULAR** (`ω₀=c_γ k_⊥` from the trapped-mode BC; `ω₀/γ` clock from a boosted wave operator;
NO `E=mc²`/Compton premise) — a sketch (asserts boost-covariance) but non-circular. Mass-bridge recorded candidate-only;
`M` not collapsed; scope CLEAN.

**Review:** Python + `.wl` transliteration both FIDELITY-CLEAN; `.wl` confirmed INDEPENDENT. Non-blocking: the
"group-velocity" check is dimensionally vacuous (doesn't differentiate the dispersion); the sonic `1/√3` is asserted as
a literal not derived-in-code (Bernoulli is prose-only → dual-engine agreement on it isn't meaningful) — both labeled
conditional, neither changes a verdict.

**STRATEGIC FINDING:** the emergent constants CANNOT be closed by dimensional analysis + the existing source equations —
the NUMBERS (`c_γ/c_s`, the flux law, and downstream `G`/`α_J`) all need (a) the SOLVED stationary defect profile and
(b) for `c_γ/c_s`, the coupled GNLS+Maxwell linearization (does the photon ride `g_μν(c_s)` or `η_MN`?); `ħ` needs NEW
substrate microphysics. This is the SAME bottleneck as the core Path-A solver → the emergent-constants chunk and the
main solver work converge here. FORK presented to the user (§0).

## 10. pathA_20b (c_γ vs c_s coupled linearization) execution + review ledger — 2026-06-19
**Run:** Codex gpt-5.5 @ xhigh. Deliverables: harness group in `dimensional_check.py` (+486, side-by-side,
`--patha20b-cgamma-cs`), `reports/pathA_20b_cgamma_cs_linearization.md`, `tools/pathA_20b_cgamma_cs_crosscheck.wl`.
Python 11 dim + 7 alg; Mathematica PASS; full suite green. Acceptance = `PASS_WITH_NAMED_RESIDUALS`.

**Verdicts (both UNDETERMINED, VERIFIED HONEST + sharper than pathA_20 — adversarial audit):**
- `bulk_verdict = C_GAMMA_BULK_UNDERDETERMINED` (`BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`). The coupled principal
  symbol FACTORIZES — off-diagonal GNLS↔gauge couplings are lower-order London/plasma terms (`δJ⁰=q★δρ`,
  `−(q★/m)ρ0 δA^i` algebraic in `δA`) below the `∂²` cone — giving `det P = ħ(ω²−c_s²k²)·(C_E(ω²−c_bulk²k²))²`. Photon
  cone `c_bulk²=C_B/C_E` (`Z/μ0` cancels). DERIVED, verified (not asserted).
- **KEY PHYSICS: `c_γ=c_s` is NOT forced.** The parent gauge metric is bulk Minkowski `η_MN` (no `K/ρ0/m`); `c_s` is an
  emergent Bogoliubov speed from a DIFFERENT sector; the matter sector is NON-RELATIVISTIC (parabolic, 1st-order in `t`)
  so it imprints NO coordinate light cone on the shared `t` → sharing `t` does NOT lock the cones. **`c_bulk/c_s` is a
  CALIBRATABLE normalization freedom (the `η` time-vs-space ratio), not a derived number and not a permanent gap.** If
  `c_bulk` fixed, `λ_γ=c_γ/c_s ∝ ρ0⁻²`.
- `brane_verdict = C_GAMMA_RATIO_STILL_UNDERDETERMINED` (`BRANE_ZERO_MODE_REDUCTION_UNDERIVED` +
  `BRANE_PHOTON_CONE_REQUIRES_PROFILE`) — the observed photon cone needs the zero-mode reduction / solved profile (parent
  paper declares the Maxwell/mixed spectrum OPEN: part01:1502,944; pde.tex:541-565). pathA_21 consumes this, `λ_γ` symbolic.

**Background legality:** legal only with a neutralizing external source `J_ext0⁰=−q★ρ0` (jellium; sourced
pde.tex:370-374) — handled, not punted.

**Review:** Python + `.wl` both FIDELITY-CLEAN; `.wl` INDEPENDENT (own CAS, computes the determinant natively); the
factorization + `c_bulk²=C_B/C_E` + the "unsourced `c_γ=c_s`" claim VERIFIED DERIVED; negative control BINDING (a forced
`c_γ=c_s` genuinely fails); coupled-symbol GENUINELY-COUPLED; scope CLEAN. **Non-blocking caveats (recorded):** (1) the
`lambdaLogSlope=−2` harness check is TAUTOLOGICAL in code (hardcodes `ρ0⁻²` rather than deriving it from `c_s∝ρ0²`) — the
CLAIM is hand-verified correct but the check doesn't back it; (2) the `forced-equals` check is cosmetic (binding teeth =
the residual-independence check) and the `coupled-det` check restates the block factorization → "7/7 algebraic" is NOT 7
independent confirmations; (3) the bulk residual is better framed as a CALIBRATABLE FREEDOM than a gap. [If we revisit
this harness, make `lambdaLogSlope` a real derivation.]

**STRATEGIC STATE:** pathA_19/20/20b now CONVERGE — the emergent-constants program has extracted everything derivable
from the symbolic action + dimensions. `c_γ/c_s` is a CALIBRATION KNOB; the flux law, brane cone, `G`, and `α_J` all need
the SOLVED THROAT PROFILE; `ħ`-emergence needs new substrate physics. Remaining moves: (A) `pathA_21` SYMBOLIC (complete
the unknowns/knobs map) or (C) the throat-profile solve (closes the numbers). User decision (§0).

## 11. pathA_21 (G + mass-bridge + profile-solve spec) execution + review ledger — 2026-06-19
**Run:** Codex gpt-5.5 @ xhigh (`_scratch/pathA_21_execute.log`, 223k tokens, exit 0). Deliverables: harness group in
`dimensional_check.py` (+887/-0, side-by-side, `--patha21-emergent-g`; 16 dim + 4 alg checks),
`reports/pathA_21_emergent_G_mass_bridge.md`, `tools/pathA_21_emergent_G_mass_bridge_crosscheck.wl`. Acceptance string =
`PASS_WITH_NAMED_RESIDUALS`. **NOT YET COMMITTED** (review found corrections needed → pathA_21b first).

**Verdicts as the run reported them:** P1 `G_FREE_PROFILE_FUNCTIONAL_DERIVED_CONDITIONAL_REDUCED_3D`
(`C_F,12=(m_GNLS·N_∞,3·Q1·Q2/4π)·I_F,12`, `r⁻²` reduced-3D / `r⁻³` bulk, attractive); P2 `MASS_BRIDGE_FORM_NOT_DERIVED`
+ `EP_NOT_DERIVED`; P3 `{L,T,M}` retained; P4 `NEWTON_G_FORM_NOT_DERIVED` (conditional `m↔G` with named `W_eff`); P5 a
35-row spec (16 profile-solve / 14 pathA_22 / 5 new-physics).

**Review = 5 clean agents (2 fidelity, 2 adversarial, 1 spec-concreteness). Outcome: TWO real shortfalls behind exit-0.**
- **NEGATIVES ALL CONFIRMED HONEST (trustworthy):** P2/EP/P3/P4 re-confirmed by independent adversarial agents. The only
  mass law in the sources is the volume-deficit scaling `M~ρ0πa²L` (`brane_bulk_ontology.tex:1294-1302`,
  `em_fields.tex:393-400`); the historical inflow was *tuned* to match Newton (`em_fields.tex:133`), not derived; the
  matter Noether charge is particle NUMBER not mass/energy; EP appears only as `m_G` reused on both sides of
  `em_fields.tex:95-98` (imposed by symbol reuse). `Γ∈ℤ` is an imposed classical label → no ħ-free relation → `{L,T,M}`
  correctly retained. No premature punts, no hidden restatements in P2/P3.
- **P1 IS A RESTATEMENT (overclaimed).** The `1/r²` drain field is HAND-INSERTED (`velocity_from_1=-q1/(4π r²)` at
  `dimensional_check.py:1965`), never solved from continuity (no dsolve/Green's fn anywhere); the "derivation" check
  reduces to `x==x` (`:1971-1976`); the attractive SIGN comes from `positive=True` declarations + choosing `Q_i>0`
  (`:1959-1964`), a convention the directive explicitly forbids. The `.wl` mirrors the same tautology. The report did NOT
  hard-PASS it (carries `PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED`), but the verdict label `..._DERIVED_CONDITIONAL` +
  prose overclaim. → corrected in pathA_21b.
- **P5 SPEC IS HONEST BUT NOT COMPUTABLE.** Schema-complete, anchors mostly real (𝔅 complete), but a residual ledger not
  a closed BVP: 6 gaps — G1 no stationary field-eq set (no time-indep GNLS, no stationary Maxwell, no `R0(w)`-selecting
  eq), G2 `R0` BC a placeholder, G3 `Π_cross` integrand never written, G4 no `J`-selecting BC, G5 reduction kernel
  `W(w)/χ_N` undefined, G6 brane `c_γ` needs the underived zero-mode. One overstated anchor: `α_H,ω` cites
  `pde.tex:318-391` for a "canonical Hamiltonian" that isn't there (only the action is). → addressed in pathA_21b.
- **Harness:** fidelity-clean, dual-engine genuinely INDEPENDENT (separate JSON, no cross-read), additive (887/0,
  pathA_18/19/20/20b untouched). Caveat (carried): the only check of `h=2πħ` is an `expr===expr` identity (no physics);
  several P1/P4 algebra checks are definitional cancellations.

**Net:** pathA_21 CONFIRMED the wall rather than breaching it (consistent with the pathA_20/20b convergence). Nothing new
rigorously DERIVED; the value is the re-confirmed honest negatives + the residual ledger. The "distrust-all-clean"
instinct caught the P1 overclaim that exit-0 + green checks masked.

## 12. STATUS / NEXT ACTION (resume here — 2026-06-19, post pathA_21 review)
**Where we are:** pathA_21 executed + reviewed (§11). Negatives solid; P1 overclaim + non-computable spec found.
`pathA_21b` directive DRAFTED (`directives/pathA_21b_force_closure_and_profile_bvp.md`) to (1) DERIVE the inter-defect
force honestly (solve continuity for the drain field — no hand-inserted `1/r²`; explicit `Π_cross` surface integral; sign
from the pressure/compressibility response or honest residual; fix the P1 verdict label) and (2) turn the P5 ledger into
a CODEABLE stationary-throat BVP, carrying P2/EP/P3/P4 + the pathA_20/20b residuals UNCHANGED.
**DESIGN-REVIEW = SOUND-WITH-FIXES → all 8 fixes applied.** KEY refinement (correct, consistent with pathA_20/21): the
draft's "close G1–G5" bar would PRESSURE A FAKE BVP — only **G1 (stationary GNLS+Maxwell) is genuinely
transcribable/closeable**; **G2 (`R0(w)` free-boundary selector), G4 (the `J`-value), G5 (kernel shape)** are
BRANCH-REALIZATION data the parent does NOT select → NAMED RESIDUALS, not closures. So each gap is now framed
"source-anchored closure OR named branch-realization residual" (faking a closure the parent doesn't provide = FAIL); the
sign/`positive=True` loophole, the P1-label resolution, the Π_cross↔P1 contradiction, and the "produces numbers"
overpromise were also fixed. So the HONEST expected pathA_21b result is: G1 closed + the inter-defect force solved in the
far field (sign from EOS or residual) + G2/G3/G4/G5/G6/`α_J`/`ħ` as named, enumerated residuals → option C is codeable
for the transcribable sector and has an explicit branch-realization to-do list.
**CONFIRM-PASS = SOUND-AS-IS** (`_scratch/pathA_21b_directive_confirmpass.log`; all 8 fixes APPLIED, new-problem scan
clean). Directive EXECUTION-READY; execute prompt STAGED at `_scratch/pathA_21b_execute_prompt.md`.
**NEXT ACTIONS (in order):** (1) **USER GATES execution of pathA_21b**; (2) FIRE
`codex exec --sandbox workspace-write -m gpt-5.5 -c model_reasoning_effort=xhigh - < software/stage1_solver/_scratch/pathA_21b_execute_prompt.md > software/stage1_solver/_scratch/pathA_21b_execute.log 2>&1`
(backgrounded, NEVER shell-timeout-wrapped) → review (transliteration + distrust-restated-target adversarial +
completeness-critic); (3) gate to **option C** (the throat-profile SOLVE, now specified by the BVP) → then `pathA_22`
(scale-map → `m̂0²·S_port` → re-run B2c). [Superseded by §0 + §13; pathA_21+21b+21c COMMITTED 2026-06-19.]

## 13. pathA_21b + 21c execution/review + Phase-1 solver reconnaissance ledger — 2026-06-19 (resume detail in §0)
**pathA_21b** (force closure + codeable BVP): EXECUTED + 4-agent reviewed. G1 stationary GNLS+Maxwell BVP genuinely CLOSED
+ codeable; drain velocity genuinely Gauss-solved (`r`-power from the surface measure); residuals G2 (`R0` selector) / G4
(`J`-value) / G5 (kernel) / G6 (brane cone) honest branch-realization residuals; branch-realization to-do = a finite
4-item list. Caught overclaim: the force COEFFICIENT was still a heuristic product → pathA_21c.
**pathA_21c** (force + sign from the Noether stress tensor): EXECUTED + 2-agent reviewed — **FORCE MILESTONE.**
`Π_ij=m_GNLS ρ v_iv_j+δ_ij P+σ_Q,ij`; balance law `∂_t g_i+∂_jΠ_ij=f_i^body` VERIFIED to reproduce the parent Euler
identity (sub-identities `dP=ρdh`, `∂_jσ_Q,ij=ρ∂_iQ` MACHINE-VERIFIED — adversarial-confirmed REAL); force = a REAL surface
integral [convective −4/3 + Bernoulli +1/3 = −1] → `F_12=−(m_GNLS N_∞,3 Q1Q2/4π r²) r̂` (inverse-square; bulk R⁻³;
structure+power EMERGE from the Gauss substitution, NOT the heuristic product); **SIGN = FORCE_ATTRACTIVE_DERIVED for like
drains** (target-blind, no positivity smuggle); full sign = `SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE`. Calibrate-predict
ledger honest (3 predictions: structure/power/sign — vs 1 knob: normalization `I_F,12^full`/`Θ_Q`). Carried negatives
verbatim. **Net: the analog yields an inverse-square ATTRACTIVE inter-defect force for like drains, derived target-blind
from the superfluid stress tensor — gravity's sign+power-law signature out of the PDE.** Minor carried harness-tightening:
the angular `4/3`/`1/3` are asserted (hardcoded `(1+1/d)`,`1/d`), not `∮n_in_j dΩ`-integrated — genuineness rests on the
hand-verified prose; tighten one engine later. The Noether tensor is posited-then-verified-by-divergence ("representative"),
not a from-scratch executed variation.
**PHASE-1 SOLVER RECONNAISSANCE (3 clean agents) — THE WHALE REFRAMED.** The throat-profile solve (option C) was audited:
(a) the PDE OPERATORS are MODEL-FAITHFUL (real quintic EOS, covariant Laplacian, gauge coupling, localized Maxwell H=Z,
ψ-sourced current — `coupled_branch.py`/`operators.py`); "engineering-smoke" = parameter values + constitutive shapes
(`V_conf`/`Z(w)`/`R0`), NOT the operators. (b) the solver + calibrate-predict harness are ~70–80% BUILT + VALIDATED
(Newton-Krylov; chunk 1a/1b/1c DONE incl. a self-caught overclaim remediation; B2c prototype: knob→real-solve, extraction
map, anchor root-finder, held-out-surplus skeleton w/ target-blind firewall). (c) **the blocker is NUMERICAL CONDITIONING,
not physics and not a from-scratch build:** no quantum-core safeguard (`√ρ→0` degenerates the matter Jacobian); the
`k1∝r⁴/R0⁵` wall coupling blows up unclamped as `R0→0`; uniform grid (no wall grading); non-production preconditioner.
LIVE EVIDENCE: B2c stalls at `τ≈0.029` = a diagnosed conditioning floor. (d) the R0 constitutive family is a posited
placeholder (a calibration knob, not yet a frozen branch). → THE ATTACK = the **conditioning spike**, directive DRAFTED at
`directives/pathA_C0_conditioning_spike.md` (density floor/log-ρ + clamp `k1` + Jacobi scaling + depth-homotopy + a
diagnostic that decides spike-sufficient vs production-solver-required). NEXT after /compact: design-review pathA_C0 →
confirm-pass → execute → review.
