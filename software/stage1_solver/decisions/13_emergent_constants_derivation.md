# Decision 13 — B2c verdict is UNDETERMINED; the real next chunk is the EMERGENT-CONSTANTS derivation

**Date:** 2026-06-19
**Status:** DECIDED direction (user-driven, 2026-06-19). This is the **resume-here record after /compact** for the
Path-A build. Supersedes the "rigorous MISS" reading of B2c (decision-12 B2c STATUS block — now flagged superseded).
**Mechanism:** B2c 3-round build/review → two audits (`pathA_17` validity, `pathA_18` dimensional) → two independent
verification agents → user methodology call (derive the emergent constants before `m̂0²·S_port`).

---

## 0. STATUS / NEXT ACTION

> **⭐ LIVE STATE (2026-06-25) — thin pointer; canonical "you are here" = repo `STATUS.md`.**
> - **Gravity sector is BUILT & GR-matched (calibrated) — DO NOT re-derive it.** The conservative PN ladder (1PN→4PN + 2.5PN
>   radiation-reaction) lives in `research/4d_*pn*`; the moving-drain no-aberration result is the "Static Limit Theorem" in
>   `research/1pn_orbital_dynamics`. Calibration is fine; first-principles is not required.
> - **pathA_28 (monopole/dipole gravitational radiation) DONE = `MONOPOLE_DIPOLE_RETURN_CONDITIONAL`** (dual-engine; arbiter PASS +
>   fidelity CLEAN; adversarial CONCERNS ⇒ it is a verified **constraint-spec, not a falsifiable test**). To avoid GR-forbidden
>   ℓ=0/ℓ=1 radiation the brane↔bulk return must deliver **`R0=−M0`** (kills the raw `O(ω¹)` monopole) and **`R1=−D1`** (kills the
>   raw `O(ω³)` dipole — net momentum-rate *including the carried odd wake*; global conservation alone is NOT enough). Artifacts:
>   `tools/pathA_28_monopole_{sympy.py,.wl}` + `reports/pathA_28_monopole*.{md,yaml}` + `reports/pathA_28_cancellation_condition.yaml`.
> - **pathA_29 (TRACK 3 gate-1: brane↔bulk return) DONE + VERIFIED 2026-06-25 = `RETURN_RESIDUAL_PREDICTION`** (tri-reviewed:
>   arbiter re-run reproduced + adversarial CLEAN + fidelity FAITHFUL). Given the **drain premise** (Z<0 = gravity is the inflow),
>   **1/r² Newtonian gravity SURVIVES the finite slab** (both DC-sink return completions solve a normalizable m=0 transverse zero
>   mode → p=2 via a counterfactual-guarded 3D-radial `dsolve`); the drain comes **bundled with an unavoidable bounded
>   monopole/dipole c_s-radiation residual ∝ ε0=1−𝒯₀(0)** = the **falsifiable departure from GR** (Birkhoff forbids it; the drain
>   breaks brane mass conservation). NOGO reachable via a derived delocalizing warp (p=3). **Sharpens but does NOT close**
>   `pde_ledger` open-item #9 (records the residual prediction; the 1/r² range item passes for the localizing flat-slab family).
>   Artifacts: `directives/pathA_29_brane_bulk_return.md` (v3) + `tools/pathA_29_brane_bulk_return_{sympy.py,.wl}` +
>   `reports/pathA_29_{brane_bulk_return.md,results.yaml}`. NEXT = the full **nonlinear** brane↔bulk return closure (downstream track-3).
> - **⭐ ACTIVE PUSH = the MOVING-THROAT PDE (a ~6-gate ladder; master checklist =
>   `research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md`; resume memory = `project-moving-throat-pde-push`).**
>   **Gate 1 `pathA_30` ✅ `DN_UNITTEST_BC_DEPENDENT`** (frozen-wall reduction reproduces the finite-throat DtN `−(ω/c_s)tan` + ladder;
>   committed `f460fc63`). **Gate 2 `pathA_31` ✅ `BREATHING_CALIBRATED`** (the distributed lift's breathing mode reproduces the legacy
>   `(a,L)` collective closure; operator-projected M/K, independent HF routes, structure gate computed+able-to-fail, dual-engine;
>   2-mode truncation clean for order-unity wall stiffness; v1 was tri-review-REJECTED for pass-by-construction → remediated; committed
>   `765db5f0`). Both land **CALIBRATED** (toy-model contract). **⭐ NEXT = Gate 3 (grouped-`P2` / ℓ=2)** — `η_2m` channels, grouped
>   response matrix, isotropy `a₂=b₂=0` (scaffold §11.3; the quadrupole sector first appears). Then Gate 4 (`54/5`) → Gate 5
>   (scalar/dipole + cross-ℓ unify pathA_28/29) → Gate 6 (nonlinear closure = the WALL); PN match-back = the decisive falsifier.
> - **NEXT = TRACK 3: the brane↔bulk return / brane parent action** — the keystone the calibrated PDE (`research/pde_ledger/`) is
>   gated on, and where the gravity sector's last falsifier actually lives (can an admissible return deliver `R0=−M0`, `R1=−D1`?).
> - **The model-mechanics corrections that keep getting lost** (nothing is static; three distinct speeds `c_s`/`v_r`/`c_γ`; gravity
>   = the FLOW not radiation; throat-soliton `J_w=0` is fine + AC→DC retired; photon ≠ gravity) live in memory
>   `[[project-model-mechanics-corrections]]` and `docs/conceptual_foundation.md` §3 (v4 block). The EM-re-founding / little-arrows /
>   GNLS polar-smectic material below is the **history** that LED here.

**⭐⭐⭐ 2026-06-23 — EM RE-FOUNDING NOW EXECUTING (RESUME HERE FIRST; front door = repo `STATUS.md`).** Chasing `λγ` exposed the
canonical EM sector DRIFTED from the single-medium concept (`reports/pathA_cgamma_of_rho_derivation.md`: EM = fundamental gauge
field on a flat metric, DECOUPLED, Type-4). Re-founded EM **medium-native**: **LIGHT = the brane's in-plane SHEAR waves** (our 3D
space = an elastic drumhead; shear ON THE BRANE not the bulk → bulk stays a pure fluid → magnetism preserved), via **MacCullagh
rotational elasticity** as the template to test (Cosserat/spinning-substructure mechanism). Accepted near-theorem: a single
**scalar** superfluid cannot carry transverse light. **Full physical picture + MacCullagh template (§11) + λγ subtlety (§13) +
honest conceptual costs (§14) + the Stage-1 leak finding (§7.1) = `decisions/15`.**

**Directive `pathA_23` v5 — THREE-WAY SOUND** (Codex design-review v3 NOT-SOUND→10 fixes → v4 Codex confirm SOUND → GLM tertiary
NOT-SOUND→3 req fixes+10 concerns → v5 Codex confirm SOUND-AS-IS). **NOW EXECUTING stage-by-stage**, each stage tri-reviewed
(orchestrator re-run + fidelity audit + adversarial review on clean agents) before its gate.
- **Stage 0 ✅ GATED:** `ACTION_SPECIFIED_CLASSIFIED`, `NEW_PARENT_ACTION`; constitutive form **POSTULATED** (honest — no
  independently-motivated microstructure available yet; declining the reverse-engineered gyrostat trap) ⇒ **conditional-verdict
  rule active**; coupling `SYMMETRY_ALLOWED_POSTULATED`; dual-engine 25/25. (`reports/pathA_23_stage0_action_and_contracts.md`.)
- **Stage 1 ✅ GATED (after a rework):** `LEAK_CONDITIONS_DEFERRED`. First attempt was TAUTOLOGICAL (hardcoded `∝k_a` sources)
  → reworked to DERIVE the bulk-stress projection. **KEY FINDING (concept stressor):** the interface traction
  `T_na = T_wa + (T_ww δ−T_ab)∂_b u_w` is **generically transverse**; no-leak holds ONLY if `P_T T_wa=0` AND isotropic `T_ab` at
  the brane — which are **NOT generic near a draining defect** (`T_wa=mρ v_w v_a` = fluid turning into the throat). So **the leak
  is EXPECTED**; survival rides on the **Stage-3 throat solve** (magnitude small OR non-fine-tuned projection cancellation);
  `v_n→bulk-vorticity` feedback = top Stage-3 priority. **DO NOT bank on no-leak.** (`reports/pathA_23_stage1_kinematic_leak_audit.md`;
  detail in `decisions/15` §7.1.)
- **Stage 2 ✅ GATED (after a rework) — THE CRUX:** `FAIL_UNSPECIFIED_SUBSTRUCTURE` (tri-reviewed: re-run 32/32 both engines +
  `FIDELITY_CLEAN` + adversarial `REWORK_VERDICT_HONEST`). First attempt returned `FAIL_NO_TRANSVERSE_STIFFNESS` but was caught
  TAUTOLOGICAL (it fixed `W=W(det F)=W(J)` at the input, forcing `μ_shear=0`) → reworked to keep the **deviatoric shear
  invariant present** (`U_∥=½K_br θ²+μ_br e⟨ab⟩e⟨ab⟩`) and route μ_br through a genuine able-to-fail substructure classifier
  (fluid-facts→0, network-facts→>0, actual-record→UNDETERMINED). **VERDICT: brane-shear EM is NOT derivable from the current
  medium specification** — μ_br is genuinely undetermined (record motivates cohesion/elasticity but never pins persistent
  neighbor-memory / a network free energy; §14 C9: substructure not fixed by GP/NLS). **TRILEMMA (verified algebraically):**
  μ_br=0→no light; μ_br>0 (Cauchy)→light BUT a stray longitudinal "second photon" (`FAIL_CAUCHY_STRAY_LONGITUDINAL`, Stage 4);
  MacCullagh curl-only (only clean-transverse form)→needs reverse-engineered gyrostats + the C5 gauge obstruction. ⇒ clean light
  needs an **extra postulated ingredient**. (`reports/pathA_23_stage2_constitutive_form.md`; `decisions/15` §15.)

**USER DECISION (2026-06-23) = PROCEED CONDITIONALLY (Path 1):** adopt the **rotational/MacCullagh form as a POSTULATE**
(`U_∥=½μ_R(∇×u)²`, CONDITIONAL throughout — the only form giving clean transverse light); gyrostat substructure = acknowledged GAP.
- **Stage 3 ✅ GATED:** `LEAK_BOUNDED_CONDITIONAL` (tri-reviewed: re-run 36/36 + `FIDELITY_CLEAN`; adversarial flagged it too soft).
  Postulated stress is antisymmetric (needs a postulated gyrostatic spin reservoir to close angular momentum); interface exchange is
  generically transverse; no-leak holds only under `ε_leak≪1` + an unmotivated cancellation/impedance price (else concept-fatal
  `FAIL_LEAK_BREAKS_MAGNUS`). (`reports/pathA_23_stage3_noleak_closure.md`.)
- **Stage 3b ✅ (verification):** `OVER_COUNT_CONFIRMED_CURVATURE_LOCALIZED` (`SIGMA_R_NOT_A_BULK_SOURCE`, `LIGHT_FREE_SLIPS_NO_LEAK`;
  tri-reviewed: 30/30·29/29 + `FIDELITY_CLEAN` + adversarial `RESTS_ON_MODELING_CHOICE_NEEDS_GLM`). Retired the "intrinsic-to-light
  fatal leak" reading: given the (separate-sector) action, `σ^R` is brane-internal not a bulk source; flat density-preserving light
  free-slips; leak is **curvature-localized** (∝|K|L_mix, far-field-vanishing). CAVEATS: model-contingent (rests on brane⊥bulk
  separation — single-medium vs membrane); defect/throat leak **relocated not retired** (still open, throat solve doesn't exist);
  report framing to downgrade. (`reports/pathA_23_stage3b_overcount_and_curvature.md`.)

**PIVOT (2026-06-23, user) → directive `pathA_24` = NOW EXECUTING (the "little-arrows" T1–T5 ladder; resume: live state in repo
`STATUS.md` "Next step").** The Stage-2/3 walls + "why is there a brane?" converged on the **little-arrows mechanism**: brane = a
**DOMAIN WALL** between two `±w` domains of a polar order parameter ("little arrows", carried by the medium); we = bound zero modes;
it dodges the single-well `U(ρ)`. A defect = a **PUNCTURE** → **CHARGE = puncture DIRECTION `±w`** (a binary orientation, **NOT a
winding / topological charge** — that's magnetism; this corrects the earlier "topological puncture" phrasing), **MASS** = trapped
geon standing wave, **THROAT SIZE** = tension-vs-holding-open balance; **charge⊥mass**. Light = **MacCullagh rotational-elastic**
brane shear (not Cauchy). THE PRIZE = a fundamental charged-massive-particle model uniting EM↔gravity at the defect.
**STATUS:** `pathA_24` reworked to **v2.2**, FULLY REVIEWED (Codex design-review ×5 → `SOUND-AS-IS` + GLM tertiary
`SOUND-WITH-CONCERNS` folded) & **committed `95ed2b86`**. Executing rung-by-rung (each: Codex codes+runs dual-engine → orchestrator
arbiter-re-run + transliteration-fidelity audit + adversarial → user gate; **AI prose never establishes a math fact**). **T0 DONE &
committed `f0c2745f`** (freeze `8fa41ac51e88` target-blind: O(4)-isotropic polar-OP, no easy axis, 0 independent params; dual-engine
dim-check PASS: bulk modes = 3 transverse Goldstones at `c_s` + gapped amplitude). **NEXT ACTION = T1** (stable-wall make-or-break:
profile + `σ_brane` + stability spectrum + `w`-emergent? + confinement). Honest prior: isotropic baseline likely → emergent-`w` but
UNSTABLE wall (sphere-of-vacua) = the three-way no-win leg; expect ≤2 of {light, stable wall, emergent-`w`}. Picture = `docs/conceptual_foundation.md`
(canonical) + `decisions/15` §17–§18. Memory `[[project-brane-existence-defect-structure]]`. pathA_23 Stages 4–6
(spectrum/`u_w`, Maxwell/C5, charge/cone-λγ) remain downstream and may be re-framed by pathA_24's wall result. (Also outstanding: a
GLM tertiary on the Stage-3b separate-sector question — likely subsumed by pathA_24's domain-wall derivation.) DEFERRED (parked,
`decisions/15` §9): mouth-inflow-vs-brane-leak. Falsifiable throughout; a clean break is a valid result. The dense 2026-06-21 block
below is PRIOR CONTEXT (gravity-sector verdict — still valid); `54/5` is ABSORBED (Gate 4 `GENUINE_BLOCKED`); the EM anchor `λγ` is
what this frontier pins.

---

**2026-06-21 status (prior):** (resume here — 2026-06-21 post-/compact; ⭐ PIVOT to pathA_22 m̂0²·S_port; pathA_22a DONE → `TUNABILITY_CHANNEL_PRESENT` [≥2 knobs χ_Q+c_γ/c_s], fix landed+closed (5h); pathA_22b GATED directive FINAL; **GATES 0+1 DONE + dual-reviewed + remediated** → `0a=MHAT_DIMENSIONFUL_CONFIRMED` (no pathA_22a flip); `0b=DOES_NOT_CANCEL (NOT_ESTABLISHED)` ⇒ W_eff stays on the Gate-4 critical path; Gate 1 `Z_int=BLOCKED` but CONTAINED; Gate 2 `λγ=STILL_TUNABLE` (C_B/C_E=1 computed ⇒ 1 dial) → **USER PIVOT to CALIBRATE-PREDICT (option B)** → **`decisions/14` BUILT** (value-provenance + calibration map, Codex+GLM reviewed: ~7 genuine inputs `{n,K,m,ρ0,μ0,q_*,λγ}`; everything else DERIVED or BRANCH-DETERMINED = compute-don't-calibrate) → **β resolved `BETA_GENUINE_GAP`** (λγ is the 7th genuine input; `c_s` derived, only `c_γ` free) + m̂ inconsistency HARDENED in pde.tex. ▶ **Gate 3 DONE (2026-06-22): χ_Q COMPUTED ≈ 0.712** [±0.0008 numerical + ~8.5% one-sided Z_w-reference definitional systematic; exact spherical-Hankel l=2 outgoing-DtN on the frozen finite-core branch; even-`nw` converged (jackknife-stable, `nr`-independent, flat tail to `nw=320`); FORM ω⁵/l=2/(1/27) derived exactly; **8-round adversarial saga** — caught seed-circularity → R⁵ placement-artifact → understated-error → premature-knob → false-3pt-convergence → spurious-W-lane-parity, then converged on even-`nw` once swept wide enough (user's "need more data" call flipped knob→computed)]. KNOWN ISSUE: shared W-lane stencil has an exact odd-`nw` null (Gate-3 used even-`nw`; other gates TBD) — see `decisions/14` §1c note. ▶ **Gate 4 DONE (2026-06-22): the gravity ratio `g_mhat²/g_G` is `GENUINE_BLOCKED`** (`BLOCKED_NEEDS_SOURCE_MAP_PROVENANCE`, dual-reviewed fidelity `FAITHFUL` + adversarial `GENUINE_BLOCKED`). `K_stress=χ_N·ρ_∞4` derives; `K_source` has NO target-blind kernel (all 22 `m̂` sites in `pde.tex` are target-facing); `α_J` underived + doesn't cancel. ⇒ the model does NOT derive its own gravity coupling (`m̂0`,`α_J` imported as target-matched normalizations) ⇒ `g_G` calibrated on Newton `G`, the GR-quadrupole `54/5` is **ABSORBED not predicted**, and the **EM-sector anchor (pins `λγ`) is now LOAD-BEARING**. Falsifiable payoff = held-out surplus (g−2, 5PN, ringdown, multi-defect) on the shared derived `χ_Q`+`P0/D0`+`c_s`. Full closeout = `decisions/14` §6 + pathA_22b Gate-4 report §. ▶ NEXT ACTION = **the EM-sector anchor** (to pin `λγ` and close the count) → then the held-out-surplus predictions. (Gate 5 verdict-assembly is now an honest closeout — quadrupole absorbed — not a prediction test. The only way to recover `54/5` as a prediction is ADDING a source-map/mass-bridge postulate to the action: a modeling choice, user-elected.) (full χ_Q detail = `decisions/14` χ_Q row + pathA_22b Gate-3 report §.) (5o) post-pivot/decision-14/β; (5n) Gate-2+pivot; (5m) Gate-1; (5l) Gate-0; (5k) directive; (5g) methodology; (5a)–(5f) throat history.)
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
**⚠️ THIS C0g-DIAG VERDICT WAS LATER REFUTED BY THE C0g BUILD (see (4) below). KEPT FOR HISTORY.** The diag called a SIMPLE
FOLD at τ_fold≈0.0291233; the gauge-fixed full-budget fast-assembly crawl then CROSSED it (8 admissible states deeper) ⇒
the "fold at 0.0291233" was a CONDITIONING DIP + crippled Newton budget, NOT a turning point. The σ_min²=0.9994 "linear
fold scaling" is RETRACTED as a near-null conditioning dip. (The C0f/C0f2 "crippled-solver artifact" lesson recurred —
distrust-clean + the honest dissolve-contest caught it.) Original (now-refuted) diag verdict text follows:
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
**(4) (HISTORY — superseded by (5) below) — the C0g BUILD: B-1/B-2/B-4 + deeper-crawl DONE & DUAL-REVIEWED; B-3
(pseudo-arclength) was THEN the earned next build (now: built+ran→NOT_MEASURED→GLM-reframed to the battery in (5)).** Directive `directives/pathA_C0g_build_gaugefix_then_pseudoarclength.md`
(design-review→8 fixes→confirm/re-confirm SOUND-AS-IS; B-2 AMENDED so FOLD_DISSOLVED uses FULL Newton budget + the
`crawl_persistent_failure` guard — the anti-sandbag fix). Executed STAGED, every step fidelity+adversarially reviewed
(all HONEST/faithful/scope-clean; freeze intact — operators/physics/`(1/ξ)`/residual/export UNTOUCHED throughout):
- **B-1 gauge-fix (PATH-ONLY, verified):** removes the 511-dim gauge null space (analytic generators → `Q_perp`; step
  direction only). PROVEN path-only — with/without gives the SAME physical state (gauge-invariant Δ 3.5e-18, equal original
  residual); `r0` mode preserved (stays in the complement). Single Arbiter = original `patha_closed_branch_residual` ≤1e-6.
- **B-4 analytic/sparse Jacobian (verified residual-equivalent):** `torch.func.jacfwd` of the UNCHANGED residual, matched
  to the colored-JVP at ~1e-14 (rel≤1e-8/abs≤1e-10), **~74–80× faster** assembly (was 96–97% of cost). Used only for the
  Newton search DIRECTION; convergence still on the original residual.
- **B-2 honest re-contest → `FOLD_DISSOLVED` (EARNED_WITH_CAVEATS), then deeper-crawl → `STILL_GOING` (HONEST_STILL_GOING).**
  With gauge-fix+full-budget+fast-assembly the crawl CROSSES the putative τ_fold≈0.0291233 — **8 admissible states down to
  τ=0.0291132 (~3× deeper)**, each original-residual ≤1e-6, branch-continuous, positive density, `min_R0`≈0.79 (throat STAYS
  OPEN — and nudges UP with depth, satisfying the user's open-throat constraint), NOT near the empty-core/opening regime yet.
  ⇒ the first "fold" is REFUTED. **BUT a NEW, genuinely-real near-singularity appears ~1e-5 deeper (~τ=0.0291132):**
  σ_min(J·Q_perp) collapses 3.95e-4→7.03e-6 (factor 56, on the AUTHORITATIVE Jacobian; σ_min² fit r²=0.96 extrapolates a
  zero-crossing ≈0.0291139), and full-budget warm-started attempts below it persistently MISS tol (~1.04e-6) — **but with NO
  branch-tangent reversal.** So σ_min says "near-singular," tangent says "no turning point yet" = the SAME fold-vs-conditioning
  ambiguity, now DEEPER and unsettled. Adversarial read: plain Newton (even gauge-fixed/full-budget/fast) CANNOT resolve
  σ_min→0 at any step size; **pseudo-arclength (B-3) is the ONLY tool that distinguishes "rounds a genuine fold" from "hard
  near-singular wall"** — so B-3 is now GENUINELY EARNED (not spent on the refuted shallow fold).
- **B-3 NOT built (held for the user gate).** **NEXT = author the B-3 execution (gauge-fixed Keller pseudo-arclength through
  τ≈0.0291132), STAGED + dual-reviewed**, per the directive's B-3 spec (gauge projection INSIDE every predictor/corrector;
  ORIGINAL residual the sole arbiter; `min_R0` depth metric + margin). Outcome resolves it: rounds the turning point →
  deeper throats (real further unblock); or cannot continue → genuine branch endpoint to understand (note: this is a
  NUMERICAL/branch feature BEFORE the physical empty-core/opening, since min_R0≈0.79 ≫ 0). Then → constitutive family →
  calibrated branch → calibrate-predict (R0/J/W → anchor → SURPLUS) → `pathA_22`.
- **UNCOMMITTED-AS-OF-RESUME?** No — all of B-1/B-2/B-4/deepcrawl committed 2026-06-21 (this commit). Files: C0g build/deep
  code in `patha_c0_conditioning_spike.py` + the B-4 autodiff assembly in `preconditioners.py`; scripts
  `pathA_C0g_build_B1B2.py`/`_finalize_timeout.py`/`build_B4_B2.py`/`deepcrawl.py`; reports `pathA_C0g_build_B1B2.md` +
  `pathA_C0g_deepcrawl.md`; JSON under `runs/pathA_C0g_build_B1B2/` + `runs/pathA_C0g_deepcrawl/` (gitignored).

**(5) ⭐⭐ HISTORY (2026-06-21) — the B-3 FOLLOW-UP v2 CHARACTERIZATION BATTERY RAN → Case 1 VERIFIED. [resume = block (5b) below].**
B-3 (gauge-fixed Keller pseudo-arclength) was BUILT + RAN (committed `46a3ea21`) → honest **`NOT_MEASURED`** at the B-3.0
validation gate (reproduces the branch to ~1e-6 except at the deepest point, where it's compared vs an UNDER-CONVERGED recorded
Newton state). **GLM tertiary consult REFRAMED it:** the proposed "decisive" Phase-D overlap test was NEAR-TAUTOLOGICAL (two
near-roots' Δ is automatically in J's near-null space by the mean-value theorem — can't separate gate-artifact from
continuation-bug; [[feedback-decisive-test-not-tautological]]). Claude+Codex AGREE with GLM. **The plan is now a CHARACTERIZATION
BATTERY** — C0 read-offs → **C1 re-converge the FAIL state TIGHT (≤1e-11) from BOTH seeds [the decisive cheap test]** → C2
mode-characterization (+ bordered fold-transversality `wᵀF_τ` + independent-basis gauge test) → C3 line-scan → C4 24×24
resolution check — feeding a **6-outcome tool-selection table** (Case 0 INCONCLUSIVE / 1 under-converged-ref / 2 simple-fold /
3 bifurcation / 4 conditioning-wall / 5 gauge-or-discretization). **DO NOT pre-commit to pseudo-arclength** — only Case 2 (with
a user gate) runs the already-built pseudo-arclength; every other case STOPS for a re-gate.
- **RESULT (RAN 2026-06-21, committed `50b9f459`): Case 1 at τ=0.02911625** — both seeds (recorded `attempt_008` + continuation
  `B3_0_accepted_009`) tight-converged (orig-residual ~2e-13, step ~2e-11) to the SAME gauge-invariant root (r0_linf ~5e-12) ⇒
  the B-3.0 gate "FAIL" was a **STALE-REFERENCE artifact, NOT a buggy continuation.** Dual-reviewed: **FIDELITY-CLEAN** (all
  C0–C4 ops faithful; Single Arbiter = original residual; frozen physics zero-diff) + **VERDICT-EARNED** (genuine independent
  two-start descent — different files AND different start-residuals 2.22e-7 vs 2.10e-9 — not a cold-load short-circuit; C2
  tension disclosed). Verdict `STOP_RECONVERGE_DEEP_STATES_AND_GATE_REFINE_THEN_C2_C4`; no tool pre-committed.
- **STILL OPEN:** the genuinely-DEEPER near-singularity (~τ=0.0291132, σ_min(J·Q_perp) still dropping steeply, 7.03e-6 at the
  deepest converged τ, 7e-7 above the σ² zero-crossing 0.0291139). C2 measured one rung ABOVE it leans FOLD (isolated near-null,
  `wᵀF_τ`≈0.92 fold-transversal, dominant **r0**) but localization dissents (EXTENDED) and **C4/resolution was never run** ⇒
  fold-vs-wall-vs-discretization UNSETTLED.
- **EXECUTE (history):** `_scratch/pathA_C0g_B3_followup_execution_prompt.md` (the Case-1 run, done).

**(5b) ⭐⭐ DEEP-POINT characterization (USER-GATED "characterize the deep point") — RAN 2026-06-21, DUAL-REVIEWED. [resume = (5c)]**
**RESULT (committed this batch): the τ≈0.0291159 wall is a GENUINE near-singularity — the contest-the-wall trap did NOT recur, and
it is NOT a discretization artifact.** Stage 1 = a clean MONOTONE tight→non-tight bracket: only τ=0.02911625 tight (a cold-loaded
prior anchor — see caveat), all deeper τ STALL (2.22e-8 at τ=0.02911594, worsening with depth). **Adversarial-VERIFIED the stalls
are FULL-BUDGET (`c1_max_iters=120`, 12-step halving line-search α→0.5¹¹≈4.9e-4 then TOTAL collapse, σ_min monotonically →3.7e-7)
— plain Newton genuinely cannot cross. The 3×-recurring crippled-solver trap is REFUTED here.** Stage 2 = isolated near-null,
r0-dominant (0.994), `wᵀF_τ`≈0.92 FOLD_TRANSVERSAL, non-gauge — a fold-LIKE mode at the bracket edge. Stage 3 (C4) = the 24×24
relocation was a GENUINE full attempt that converged 4 orders then hit the SAME wall at the SAME τ (throat open, min_R0≈0.798) ⇒
**the feature REPRODUCES under grid refinement = NOT a 16×16 discretization artifact** (the run UNDER-read this as merely
"insufficient"; adversarial promoted it to physical-supporting). Fidelity = FIDELITY-CLEAN (full budget, Single Arbiter, frozen
physics untouched). **TWO honest caveats (adversarial):** (a) the tight bracket-edge is a COLD-LOADED prior anchor — no FRESH deep
tight state was reproduced this run; (b) Stage-2 "FOLD_SIGNATURE_SHARPENS_OR_PERSISTS" is a VACUOUS self-comparison (anchor vs the
PRIOR run's C2 of the SAME anchor — byte-identical, can't-fail) → DROP it; only a deeper tight state (which pseudo-arclength can
produce) gives a real sharpens-vs-dissolves read ([[feedback-decisive-test-not-tautological]]).
**VERDICT REVISED (adversarial, Claude concurs): the run's `STOP_NO_REMEDY` was TOO CONSERVATIVE.** A genuine full-budget σ_min→0
stall that is NOT discretization and NOT gauge, with a fold-like transversal mode but NO demonstrated branch-tangent reversal, is
the EXACT precondition that EARNS the pseudo-arclength (B-3) tool — the only tool that distinguishes "rounds a genuine fold"
(tangent reversal) from "hard near-singular wall." σ_min→0 ALONE stays fold-vs-wall AMBIGUOUS ([[feedback-contest-wall-verdict]]),
and B-3 is the discriminator. So B-3 is NOW GENUINELY EARNED (not spent on a refuted shallow fold). Deep-point ORIGINAL spec was:
Re-converge `attempt_022..028` TIGHT (FULL budget) = contest-the-wall → C2 at the deepest TIGHT → C4 24×24 → 6-case map.

**(5c) ⭐⭐ RESUME HERE (2026-06-21) — RUN the pseudo-arclength (B-3), USER-GATE PENDING.** The deep-point battery EARNED it.
Bring to the user: gate to run the already-built gauge-fixed Keller pseudo-arclength (B-3) THROUGH the τ≈0.0291159 stall and read
**tangent-reversal (genuine fold → continue to deeper throats) vs hard-wall (corrector also fails → LM/trust-region or rethink).**
Secondary asks the next B-3 run must honor (adversarial): (i) DROP the tautological Stage-2 self-comparison — compare only against a
genuinely DEEPER tight state B-3 produces; (ii) reproduce ≥1 FRESH deep tight state (don't rest on the cold-loaded anchor);
(iii) credit the 24×24 reproduction as physical-supporting + optionally run 32×32 to confirm. NOTE: the FIRST B-3 run's B-3.0
"FAIL" was the Case-1 STALE-REFERENCE artifact (now refuted) — the pseudo-arclength continuation itself WORKS (reproduces the branch
to τ≈0.0291162); the task is to push it PAST the stall. **Build/run only after the user gates** (Case-2-tool rule).

**(5c) RESULT — B-3 pseudo-arclength RAN 2026-06-21 (user-gated), DUAL-REVIEWED. Outcome A GENUINE_FOLD — EARNED-BUT-OVERSTATED.**
FIDELITY-CLEAN (fold detector fail-able — planted-wall→Outcome B; secant-fallback bug FIXED+unit-tested; Single Arbiter held;
frozen physics 0-diff). Adversarial = **EARNED-BUT-OVERSTATED**. The numerical fold is real: from the anchor (τ=0.02911625, fresh
tight reconverge) the gauge-fixed pseudo-arclength produced 6 corrector-accepted points, `dτ/ds` through 0 at τ_fold=0.029115870
(2 down / 3 up), Single Arbiter held (max orig-residual 1.33e-7 ≤1e-6), σ_min(J·Q_perp) recovered 12.5× (1.0e-5→1.26e-4 monotone).
**BUT two corrections (override the report prose):** (1) **LIKELY TURNED EARLY** — it turned at σ_min≈1.0e-5, yet the deep-point
fixed-τ run reached σ_min≈4.8e-7 at τ=0.0291150 (DEEPER + ~20× more singular); a genuine fold turns AT the σ_min minimum, so
fold-vs-wobbled-shallow is UNRESOLVED (needs smaller-ds). (2) **NO DEEPER THROAT** — min_R0 RISES 0.7975→0.7992 (throat OPENS),
min_rho flat ~7.2e-6; the report's "Fresh Deep Tight State" is mislabeled (branch-accepted 6.7e-8, not tight; not deeper than
fixed-τ Newton). Report annotated with the correction (`reports/pathA_C0g_B3_run_pseudoarclength.md` header). Code/report committed.

**(5d) ⭐⭐ RESUME HERE (2026-06-21) — STRATEGIC PIVOT: τ does NOT control throat depth.** Across the deepcrawl + deep-point +
B-3 runs, **min_R0 is STUCK ~0.798 for the ENTIRE τ range** ([[project-pathA-build]]: "τ-depth ≠ R0-depth; empty-core NEVER
approached" — now CONCRETE). So pushing the τ-continuation — fold or wall — does NOT deliver the deep/empty-core throat the program
wants (R0→small). The τ≈0.0291159 obstruction is now CHARACTERIZED (a fold, modulo the turned-early check) and is a DEAD END for
throat-deepening. **The open question to bring to the user (DONE — see options): what actually controls throat depth?** Candidates:
the R0 CONSTITUTIVE FAMILY (the posited placeholder knob, never promoted to a frozen branch — the most likely real depth knob);
a different continuation parameter; or the calibrate-predict reframe (branch data = KNOBS, so maybe a deep throat is SELECTED via a
constitutive choice, not crawled to via τ). Cheap optional side-quest: smaller-ds B-3 re-run to settle turned-early-vs-genuine-deep-fold
(low priority if pivoting off τ). The report "Deep/Tight" mislabel + min_R0 caveat are corrected in decisions/13 + the report header;
a future Codex pass can clean the report body.

**(5e) ⭐⭐ RE-SCOPE (USER-GATED "do we even need a deep throat?") — CLAUDE(2 agents)+CODEX CONSENSUS: REFRAME-SOUND-WITH-CAVEATS;
GLM TERTIARY PENDING.** Two clean read-only investigations (calibrate-predict pipeline + per-observable physics) AND a Codex
adversarial refutation attempt ALL conclude: **NO downstream observable requires an EMPTY-CORE throat (R0→0, ρ_center→0); the
converged OPEN finite-core family (min_R0≈0.798) is SUFFICIENT to proceed** to constitutive-family → calibrate-predict → `pathA_22`.
Evidence: extraction calibrates over τ (root-find `R_norm(τ)=0`), the only geometric gate is OPEN-throat (`R0(L)>0`, `hard_cap`
REJECTED — the OPPOSITE of empty core); the 1c run extracted the full chain at finite core R0∈[0.99,1.0]; FORCE (pathA_21c)
"does NOT need the solved interior profile" (far-field + normalization knob); κ_PV positively REQUIRES a FINITE core (a0≈0.79;
empty core makes its functional singular); 5PN/g−2/moving-throat need finite-aperture overlaps + the `D0=K−B0−Z0→0` near-pole gate,
never ρ_center→0; canonical ontology (`pde_audit_full.md`, `defect_interaction_map.md`) REQUIRES a finite OPEN conduit (hard-cap
forbidden). **The "empty-core" target was introduced in the C0b solver directive as a re-interpretation of "deep/realistic" with NO
downstream derivation behind it ⇒ the τ-deep-throat "white whale" was largely chasing a numerical feature off the critical path.**
**THE ONE SERIOUS CAVEAT (Codex, → GLM):** the force derivation assumes `core scale ≪ aperture a ≪ r12`; the frozen family has
`a=1` (FROZEN_A) and min_R0≈0.798, so `0.798 ≪ 1` is FALSE — the current branch may NOT be in the controlled point-defect regime.
This needs SCALE SEPARATION (core≪a), NOT an empty core — satisfiable by a smaller core / larger aperture / explicit finite-aperture
(1/a) correction. Lesser risks: the one-knob B2c root isn't currently found (wants near-pole `D0≪K` — naturalness); emergent
constants need "a profile," not an empty one. **GLM prompt:** `_scratch/reframe_empty_core_GLM_consult_prompt.md` (user runs GLM;
focus = is core≪a a HARD prerequisite for pathA_21c/pathA_22, and can the finite-core calibrate-predict branch proceed now with
scale-separation as a parallel correction). After GLM: if 3-way consensus holds → PIVOT off the τ-deep-throat hunt → promote the
constitutive family + run calibrate-predict on the finite-core branch (the scale-separation fix is the new, cheaper sub-problem).

**(5f) ⭐⭐ GLM TERTIARY VERDICT (2026-06-21): AGREE — 3-WAY CONSENSUS. The REAL critical path is `m̂0²·S_port`, NOT the throat solve.**
GLM AGREES the empty core is OFF the critical path (independently: κ_PV is **structurally undefined** at R0→0 — the virial lock
`E_w+2E_f=3E_PV` breaks, not merely singular). GLM DOWNGRADES the scale-separation caveat from "hard prerequisite" to a **soft
residual on one calibration knob**: in the force, `v₂∝Q₂/a²` and the control-surface area `∝a²` CANCEL, so the force STRUCTURE
(`Q₁Q₂/r₁₂²`) + SIGN survive; the O(1) `core/a` correction shifts ONLY the prefactor, which is ALREADY a labeled calibration knob.
⇒ add `SCALE_SEPARATION_RESIDUAL` to the force ledger (alongside the 4 existing residuals), fix IN PARALLEL via a different wall law
`S_Σ` (R0 is set by `S_Σ`, NOT τ — do NOT resume the τ crawl; the fold is irrelevant). **THE REAL BLOCKER (GLM Q4#1 — CONFIRMED by the COMMITTED REAL
B2c code output `reports/patha_b2c_calibration_report.md` (Jun 19; dual-engine cross-check ~1e-16, error bars, negative-control
gates all failed-as-expected; verdict MISS): measured P0=2.795e-9 (τ=1) … 1.224e-6 (τ=0.029125), 6.7–9.6 orders below, R_norm≈−10.8
at EVERY τ ⇒ NO root with m̂0²·S_port=1. Not prose/GLM-reasoning — actual calibration output.): P0 sits 6.7–9.6 orders BELOW the
`54/5=10.8` target at EVERY converged τ (P0=2.795e-9 at τ=1 … 1.22e-6 at
τ=0.029); `R_norm=m̂0²·S_port·P0−54Gc_s⁵/(5a⁵c⁵)` with `m̂0²·S_port` PINNED=1; the ~9-order miss = exactly what the DIMENSIONFUL
`m̂0²·S_port` would supply (10.8/P0=3.86e9 at τ=1). So `R_norm(τ)=0` has NO root at any τ — the fold + scale-sep are SIDESHOWS.** The
critical path is DERIVE `m̂0²·S_port` (the bulk→brane natural-unit→physical scale map = `pathA_22`), logically PRECEDED by the
emergent-constants derivation (c_s/c/G — task #71; §13 already lays out: scale-map → `m̂0²·S_port` → re-run B2c → real verdict).
Whether derived `m̂0²·S_port`≈1 (→ real miss) or carries orders (→ match) is the decisive, FAIL-ABLE test (line 406). **⇒ PIVOT:
STOP throat-solver work; the min_R0≈0.798 finite-core family is the branch to extract from; resume the emergent-constants →
`m̂0²·S_port` (pathA_22) path.** [Open methodology Q for Claude+Codex at pathA_22 entry: is `m̂0²·S_port` a DERIVE-target (dimensionful
unit-map, decision-13's reading) or can the dimensionless surplus-RATIOS (P2/P0, P4/P0, R_pole) be predicted with it as a free overall
scale knob? — settle before claiming a verdict. GLM consult prompt + response archived under `_scratch/reframe_empty_core_*`.]
- **EXECUTE:** `codex exec --sandbox workspace-write` on **`_scratch/pathA_C0g_B3_deeppoint_execution_prompt.md`**, launched as a
  `Bash run_in_background:true` task ([[feedback-background-process-launch]]); then fidelity+adversarial dual-review, bring the
  Case + tool to the user. **Design-review:** Codex read-only on `_scratch/pathA_C0g_B3_deeppoint_design_review_prompt.md`
  (IN PROGRESS at last checkpoint — apply MUST-FIX before launch). Monitoring policy (no 600s cap) applies.
- **Contract:** the directive's "⭐⭐ B-3 FOLLOW-UP v2 (AMENDMENT, 2026-06-21, GLM-REFRAMED)" section
  (`directives/pathA_C0g_build_gaugefix_then_pseudoarclength.md`). **Loop CLOSED:** design-review (Codex AGREES w/ GLM) =
  SOUND-WITH-FIXES (12) → fixes applied → confirm-pass SOUND-AS-IS (1 trivial Case-0 enumeration fix applied after).
- **TIMEOUT:** the 600s cap is LIFTED for solver runs (user, 2026-06-21) → forward-progress monitoring (no wall-clock cap;
  emit incremental progress; internal no-progress guard N=30 / ε_res=1% / ε_adv=1e-9; orchestrator watches, stops only a
  genuinely-stuck run). `.wl`/SymPy keep the cap. (See the standing-flag update above + [[feedback-script-timeout-policy]].)
- Read-offs DONE (confirm exact in C0): raw σ_min(τ) finite 7e-6 at deepest (σ_min² fit zero-crossing UNRELIABLE — wall-vs-fold
  OPEN); arclength metric x-dominated near the stall (method not degenerating to τ-continuation).
- **COMMIT STATE:** the v2 directive checkpoint is committed (see the commit after `46a3ea21`); the execution prompt lives in
  `_scratch/` (gitignored). LESSON: a "decisive test" must be able to FAIL for the bad case ([[feedback-decisive-test-not-tautological]]).
**LESSON (now a memory): the C0f/C0f2 "crippled-solver artifact" trap RECURRED at the diag layer — a clean fold signature
(σ_min² linear R²=0.999 + cosθ transversality) was really a near-null CONDITIONING dip a budget-2/0-backtrack solver
couldn't push through. ALWAYS contest a "wall"/"fold" verdict with the FULL solver budget + faster assembly BEFORE building
the fix; the honest dissolve-contest (full budget, reachable DISSOLVED branch) is what caught it. Also: a "gauge-invariant"
metric built as derivative/value is dimensionally a wavenumber and amplifies high-k remnants — use the dimensionless energy
fraction.**
**⏱ STANDING FLAG — `timeout 600` cap — UPDATED 2026-06-21 (USER DECISION):** the user LIFTED the 600s cap **for the Path-A
throat-solver NUMERICAL runs** (it was a forcing-function meant for the `.wl`/SymPy DERIVATION scripts, which KEEP the cap).
For solver runs the cap is REPLACED by **forward-progress monitoring**: NO hard wall-clock cap; a multi-hour run is fine *as
long as it makes forward progress*; every solver script EMITS incremental progress + checkpoints resumably + carries an
INTERNAL no-forward-progress guard (stop + report STALLED if N consecutive steps make no measurable progress); the
orchestrator MONITORS the emitted progress and TaskStops a genuinely-stuck run but lets a progressing one run as long as it
needs. So the old "walk the ladder, raising the per-chunk cap is the last rung / user-level call" is now MOOT for solver runs
(the user has made that call). **Still binding for `.wl`/SymPy derivation/audit scripts:** keep `timeout 600`, a timeout ⇒
`NOT_MEASURED`, never raise. ([[feedback-script-timeout-policy]] updated.) (Do not delete on compact.)
**Discipline reminder:** Codex derives/codes + applies fixes, Claude reviews; orchestrator owns directives/decisions.
The DERIVED-FORM GATE binds (no hand-inserted field/`r`-power, no convention sign, no `x==x` posing as a check, no
restatement to fake BVP closure). VALID expected outcomes: a derived far-field force with interior factors flagged, an
honest `ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL`, and G6/`α_J`/`ħ`/`2π` remaining residuals. Commit pathA_21 + its 21b
corrections together so the ledger lands honest (commit only when the user asks).

**(5g) ⭐⭐ RESUME HERE (2026-06-21, post-/compact) — pathA_22 METHODOLOGY SETTLED (Claude+Codex) + DIMENSIONAL-SKELETON-FIRST (USER-GATED).**
The open methodology Q [(5f)] is RESOLVED by Claude+Codex (Codex consult `_scratch/codex_pathA22_methodology_consult.log`,
prompt `_scratch/pathA22_methodology_codex_consult_prompt.md`; grounded in real code + pre-reg §H, NOT prose):
- **`m̂0²·S_port` is a DERIVE-target, NOT a free scale knob.** Decisive: it's the dimensionful conversion factor `T_target/[P0]`
  (pre-reg lines 277–278 citing pathA_17/18; decision-13 §1 lines 421/423), and pinning it to 1 is a SCALE CHOICE that forces
  the `D0→0` knife-edge. The "predict scale-free RATIOS instead" escape hatch is DEAD: the only scale-free Stage-1 held-outs
  (P2=P4=R_pole=0) are THIN internal-consistency checks sharing the calibrated `D0` (pre-reg 258–264), not external benchmarks.
- **Codex CORRECTION (kept):** the two levers COEXIST, not compete. `D0` is the BRANCH lever (`P0=N0/D0`, set by the
  self-consistent background; frozen family can't reach the knife-edge → the pinned-convention MISS). `m̂0²·S_port` is a
  SEPARATE physical unit/source-map conversion multiplying the solved `P0` in `R_norm`. Verdict hinges on the DERIVED value:
  ≈1 → MISS stands; ~3.86e9 (=10.8/P0) → the SAME finite-core branch becomes a MATCH. **"Fail-able, not tunable"** — value from
  physics, NOT reverse-engineered from `10.8/P0`. Code lanes confirmed: scale enters ONLY `R_norm`/`observable_residuals`
  (`patha_extraction.py:477,544`); `gamma_eff` uses it but isn't held-out; secondary symbolic lane hardcodes `R_norm=P0−54/5`
  (scale-unaware, validates the pinned algebra only).
- **THE CRUX BOTH ENGINES ONLY ASSERTED (did NOT prove) → the gate of pathA_22:** is the final `m̂0²·S_port·P0` vs
  `54·G·c_s⁵/(5a⁵c⁵)` comparison genuinely DIMENSIONLESS (the ~9 orders fixed by model structure → a REAL fail-able test) or
  does a FREE overall scale survive (toy model, `a=c_s=ħ=m=1` are pins → the GR-quadrupole "match" would be TUNABLE, not a
  prediction)? Dimensional bookkeeping (units restored): target dim = `M⁻¹L⁻²T⁻²` (3D `G`); `m̂0²=(L⁻¹T⁻¹M⁻¹ᐟ²)²=M⁻¹L⁻²T⁻²`
  ALREADY carries the target dim ⇒ `S_port·P0` must be DIMENSIONLESS ⇒ `R_norm=0 ⟺ ξ·(S_port·P0)=54/5` with
  `ξ:=m̂0²/[G·c_s⁵/(a⁵c⁵)]` dimensionless. FAIL-ABLE iff `ξ` is a pure number fixed by the fundamental params (K,m,ħ,a,ρ₀ via
  the emergent c_s/c/G + m̂0 + S_port closed forms); TUNABLE if a free scale survives. This is [[feedback-dimensional-consistency-check]]
  + [[feedback-rescope-blocker-vs-downstream-need]] applied to the test itself ([[feedback-decisive-test-not-tautological]]: the
  pathA_22 test is only worth doing if it CAN fail).
- **SCOPE of pathA_22 (bigger than one directive; rests on pathA_21 honest-negatives):** emergent `G` (`NEWTON_G_FORM_NOT_DERIVED`,
  4D→3D brane reduction + `W_eff` kernel undefined = ledger gap G5 — THE big piece); mass-bridge (`MASS_BRIDGE_FORM_NOT_DERIVED`);
  source factor `m̂0` (`SCALE_MAP_SOURCE_FACTOR_UNDERIVED`); port convention `S_port` (frozen, not derived); the natural→physical
  scale map. The "branch packet" (`R0,N0,D0`,overlaps) input is supplied by the EXISTING converged finite-core branch — does NOT
  re-introduce the deep-solve white whale.
- **USER DECISION (2026-06-21):** **DIMENSIONAL-SKELETON FIRST** — do the cheap fail-able-vs-tunable check (units restored,
  SymPy/`dimensional_check.py`) BEFORE investing in the hard `G`/mass-bridge/`W_eff` derivations. If TUNABLE → GR-quadrupole
  anchor isn't a real external prediction (learn cheaply, re-think). If FAIL-ABLE → proceed to derive `G` + the scale map. This
  DEFERS the **GATE-A freeze amendment** (un-pinning `m̂0²·S_port=1`, hash `ed3585…`, decision-11 §2/§3 — §5 below: "settle with
  the user, do NOT do unilaterally") until we know the prize is real. GLM tertiary review = LATER (review the completed
  scale-map derivation), NOT now.
- **NEXT ACTION (resume here):** draft the dimensional-skeleton directive (`directives/pathA_22a_dimensional_skeleton_*.md`) →
  Codex design-review (read-only) → fixes → confirm-pass → execute (Codex workspace-write; SymPy/`dimensional_check.py`; KEEPS the
  `timeout 600` cap — derivation script, NOT a solver run) → dimensional-fidelity + adversarial review → bring FAIL-ABLE/TUNABLE
  verdict to the user. Discipline: Codex derives/codes, Claude reviews; DERIVED-FORM GATE binds (no `x==x` posing as a check).

**(5h) ⭐⭐ pathA_22a DIMENSIONAL-SKELETON EXECUTED + DUAL-REVIEWED (2026-06-21). VERDICT: `TUNABILITY_CHANNEL_PRESENT` (EARNED) —
the GR-quadrupole anchor is NOT a clean fail-able test as currently framed; ≥2 un-pinned normalization knobs survive.** Directive
`directives/pathA_22a_dimensional_skeleton.md` (design-review SOUND-WITH-FIXES→10 fixes→confirm-pass SOUND-AS-IS). Executed by
Codex (workspace-write): `reports/pathA_22a_dimensional_skeleton.md` + `dimensional_check.py` group `--patha22a-dimensional-skeleton`
+ tests + `tools/pathA_22a_dimensional_skeleton_crosscheck.wl`. Key reduction (verified): with units restored,
`m̂0²` ALONE carries the full target dim `M⁻¹L⁻²T⁻²`, and `R_norm=0 ⟺ P0·χ_Q·g_mhat²·λγ⁵/g_G = 54/5`
(`G=(a·c_s²/m_GNLS)·g_G`, `m̂0=(c_s/(a²√m_GNLS))·g_mhat`, `c=λγ·c_s`). **DUAL-REVIEW (2 clean agents):**
- **Adversarial (verdict): headline EARNED + STANDS.** The crux reconciliation: Stage 105 "EXACT FIXING of χ_Q=1" is
  BRANCH-CONDITIONAL — it pre-chose the canonical `σ_Q^can` + matched the free-space spherical-Hankel fingerprint (Stage 104),
  so `χ_Q=1` holds ONLY IF the Path-A outgoing branch IS that canonical mode (UNPROVEN; the notes carry `Δ_Q:=χ_Q−1` as the open
  deviation; cards say "not an unconditional actual-branch theorem"; project tag `OUTGOING_DTN_BRANCH_UNDERIVED`). So `χ_Q/S_port`
  is genuinely a class-(d) free knob for Path-A. **BUT the report's "SOLE class-(d) knob" is OVERSTATED:** `λγ=c_γ/c_s` enters as
  `λγ⁵` and is equally un-pinned (pathA_20b REJECTED forcing `c_γ=c_s`; decision-13 §3 "the one surviving open number") ⇒ **≥2
  tunability channels** (χ_Q AND λγ⁵), plus (c) residuals (`Θ_Q, α_J, W_eff/a, branch-kernel`) that go tunable under
  calibrate-predict. Not gamed (under-counts class-(d), the conservative direction).
- **Fidelity (dimensions): one FAITHFUL-BUT-WRONG dimension.** The harness ASSERTED `N0 = reduced_stiffness (=[K])`, but the code's
  actual `N0 = P²/Δ²` gives `[N0]=[K]T² = reduced_mass` (whole N-tower shifted T²; CROSS-VALIDATED: same building blocks give
  `[Z2]=[K]T²=[M]`, `[Z4]=[K]T⁴=[M]T²`, MATCHING the dictionary's Z/B contracts). ⇒ faithful `[P0]=T²` not 1 ⇒ the
  `HOMOGENEITY_PASS` is SPURIOUS at the N0/P0 link (passes only because P0=1 by fiat); the real `(c_s/a)²` factor was mislabeled a
  "convention" (masked because `m̂0`'s dim is itself a back-fit `=√(target/(P0·S_port))`). Negative controls genuinely fire
  (not hardcoded); ξ-algebra + χ_Q≡S_port mapping FAITHFUL; target dim FAITHFUL. Also: `pde.tex` internally inconsistent on `m̂`
  (`eq:outgoing-BT-target` dimensionful vs `eq:outgoing-natural-source-map` dimensionless).
- **FIX IN FLIGHT (Codex, `_scratch/pathA_22a_fix_prompt.md` → `codex_pathA_22a_fix.log`):** build N-tower dims from `P²/Δ²`
  (derive, not assert) + carry `(c_s/a)²` on P0 so the gate tests the real chain (report PASS-honestly OR a revealed gap);
  correct the report to ≥2 channels + the `m̂` inconsistency flag. Then a fidelity confirm-review.
- **⭐ STRATEGIC SCOPING PAYOFF (what this gate bought):** to make the GR-quadrupole anchor a REAL fail-able test, the critical
  path must DERIVE — not calibrate — `χ_Q` (the actual Path-A outgoing DtN branch coefficient = `OUTGOING_DTN_BRANCH_UNDERIVED`),
  `λγ=c_γ/c_s` (the brane zero-mode reduction), AND the `G`/source-map forms (`g_G`, `g_mhat`). NOT just `G`. These are the
  named sub-problems of pathA_22; they are also SHARED with the downstream observables (g−2, 5PN). The minimal combined target is
  `ξ = m̂0²·S_port/[G·c_s⁵/(a⁵c⁵)]`. Honest note (pre-reg §H): the GR-quadrupole is the ANCHOR, not held-out; its predictive
  value is limited — the real external surplus is DOWNSTREAM. **NEXT: bring verdict to user; gate the next chunk (derive
  χ_Q-outgoing-DtN + c_γ/c_s brane-reduction + G/source-map forms → minimal ξ → real B2c verdict; vs weigh leaning on downstream
  external surplus). GATE-A freeze amendment still deferred until a derived value is installed.**
- **pathA_22a FIX LANDED + CLOSED (2026-06-21):** N0-tower dims now DERIVED from `P²/Δ²` (not asserted); `P0` dimensionless only
  after explicit `(c_s/a)²` ⇒ `HOMOGENEITY_PASS` HONESTLY; report corrected to ≥2 channels + the `m̂` pde.tex inconsistency;
  tests 5+143 pass; `.wl` re-validated via `math -script` (18/18 dim incl. expected-negative, 1/1 alg). (`wolframscript` exits 255
  in this env — use `math -script`.) Task #81 done.

**(5i) ⭐⭐ RESUME HERE (2026-06-21) — pathA_22b plan: USER chose HOLISTIC ξ; Codex design-review = `UNSOUND` as a SINGLE directive →
DECOMPOSE + GLM direction-consult PENDING (user runs GLM).** Directive drafted `directives/pathA_22b_minimal_combination_xi.md`
(target-blind `P0·χ_Q·g_mhat²·λγ⁵/g_G ?=? 54/5`). Codex design-review (`_scratch/codex_pathA_22b_directive_design_review.log`):
- **VERDICT UNSOUND as one directive** — over-scoped + leans on an UNPROVEN cancellation; high risk of a pathA_21-style overclaim
  masked by green checks. **DECOMPOSE into gated sub-steps w/ HARD STOPS:** χ_Q → λγ → `g_mhat²/g_G` (prove the `W_eff/Z(w)`
  kernel cancellation FIRST, as its own deliverable) → only then assemble. Still targets ξ (user's choice stands); refines HOW.
- **Realistic outcome = honest negatives / conditional residuals, NOT a real MATCH/MISS**, unless the missing brane/profile/
  source-map equations already exist in sources unextracted. Per-factor: `χ_Q` partially tractable (canonical ω⁵ closed; ACTUAL
  Path-A branch equivalence UNPROVEN); `λγ` blocked on brane zero-mode reduction (reproduces pathA_20b negative); `g_mhat²/g_G`
  not tractable until the cancellation is proven.
- **Claude synthesis:** `χ_Q`-actual-branch + `λγ` + `g_G` all appear blocked on RELATED underived structure (brane zero-mode /
  `W_eff` reduction + actual outgoing-DtN) ⇒ a COMMON derivation-side wall (analogous to the throat solve we re-scoped off).
- **Codex+Claude AGREE on the refined plan (decompose); 4 DIRECTION-LEVEL questions → GLM** (prompt
  `_scratch/pathA_22b_direction_GLM_consult_prompt.md`, user runs GLM, paste response): (1) does `W_eff/Z(w)` genuinely cancel in
  `g_mhat²/g_G` (source-map vs back-reaction adjoint) or wishful? (2) can actual `χ_Q` come from the EXISTING finite-core branch
  via low-ω linear-response/DtN, or need a new dynamic solve? (3) is `λγ` derivable or an irreducible calibration choice? (4)
  STRATEGIC: attack the common brane-reduction blocker directly vs re-weigh downstream — and do the needed equations already exist
  in the parent action/pde.tex unextracted? **NEXT: user runs GLM → digest → revise the directive into the decomposed gated form
  per Codex+GLM → confirm-pass → execute first gate. Do NOT revise/execute before GLM (direction may change the decomposition).**

**(5j) ⭐⭐ RESUME HERE (2026-06-21) — GLM RESPONDED (constructive reframe): the "common wall" is TWO-LAYERED, one layer BREACHED →
PROCEED with a decomposed GATED plan (0→4). Codex SOURCE-VERIFICATION of GLM's claims IN FLIGHT.** GLM full response archived in
chat; consult prompt `_scratch/pathA_22b_direction_GLM_consult_prompt.md`; Codex verify `_scratch/codex_pathA_22b_glm_digest_verify.log`
(prompt `_scratch/pathA_22b_glm_digest_codex_verify_prompt.md`). GLM's reframe (claims to VERIFY vs sources — GLM ran NO code):
- **Q1 (cancellation) PARTIAL:** the `W_eff/Z(w)` brane kernel CANCELS in the ratio `g_mhat²/g_G` (same `Z(w)`, measure `√g_w dw`,
  volume powers) ⇒ **W_eff OFF the critical path for the RATIO**; the residual is an O(1) FIELD-CONTENT ratio `K_stress`(∂field)
  vs `K_source`(field value) — NOT adjoint, does NOT cancel — computable from the EXISTING profile. (Codex's "unproven" stands for
  ratio=O(1), but the KERNEL cancellation is the cheap win.)
- **Q2 (χ_Q) AGREE:** obtainable via a LINEAR-RESPONSE outgoing-DtN BVP around the EXISTING stationary branch (ω⁵ depends only on
  `P₀,χ_Q` per `pde.tex:2068`; canonical `σ_Q^can=9a⁵/(8Ω_Q⁵)` `pde.tex:1987`; framework `pde.tex:1941-2229`) — NOT a new nonlinear
  solve, does NOT re-introduce the off-path deep solve.
- **Q3 (λγ) AGREE — DERIVABLE, an EVALUATION not a derivation:** `μ₀^eff=μ₀/Z_int` already derived (`pde.tex:557`); `λγ=c_γ/c_s=
  √((C_B/C_E)·m/(5Kρ⁴))/√Z_int`; `Z_int=∫dw Z(w)` = one integral over the EXISTING exported `Z_w` (`m1c_background_export.py:169`);
  `C_B,C_E` parent gauge sector (`pde.tex:541-565`). pathA_20b's `BRANE_ZERO_MODE_REDUCTION_UNDERIVED` = bookkeeping gap (treated an
  existing formula as future work), not physics.
- **Q4 (strategic) AGREE:** equations ALREADY EXIST (parent action/pde.tex) — extraction/evaluation, NOT new physics. Realistic
  outcome converts TUNABLE knobs → DERIVED values or NAMED gaps (still decisive: unfalsifiable→falsifiable).
- **GLM's "ONE THING MISSED" (a Gate-0 PREREQUISITE):** the `pde.tex` `m̂` inconsistency (dimensionful `eq:outgoing-BT-target`
  `[m̂]=L⁻¹T⁻¹M⁻¹ᐟ²` vs dimensionless `m̂=1+O(a²/r²)` `eq:outgoing-natural-source-map`; code uses dimensionful,
  `dimensional_check.py:4280-4287`). If the DIMENSIONLESS reading is right, `m̂` is a FIELD VALUE not a normalization ⇒ the pathA_22a
  `TUNABILITY_CHANNEL_PRESENT` "tunability" becomes a PROFILE property, NOT a free knob — reconcile FIRST, in Gate 0.
- **GLM GATE ORDER:** Gate 0 (SymPy) reconcile `m̂` THEN prove `Z(w)` kernel cancellation (success⇒W_eff off path; failure⇒wall
  real) → Gate 1 (one integral) `λγ` via `Z_int` → Gate 2 (linear BVP) `χ_Q` → Gate 3 (real derivation) the O(1) field-content
  ratio `g_mhat²/g_G` → Gate 4 assemble `P0·χ_Q·g_mhat²·λγ⁵/g_G ?=? 54/5`.
- **NEXT:** Codex source-verify (in flight) → digest with Claude → rewrite `pathA_22b` into the gated form (Gates 0–4) → confirm-pass
  → execute Gate 0. Each gate is a hard stop; honest negatives still = progress (knob→derived/named-gap). [DONE → see (5k).]

**(5k) ⭐⭐⭐ RESUME HERE (2026-06-21, post-/compact) — pathA_22b GATED DIRECTIVE FINAL + design-review `SOUND-AS-IS`. NEXT ACTION =
EXECUTE GATE 0.** Codex SOURCE-VERIFICATION (`_scratch/codex_pathA_22b_glm_digest_verify.log`) TEMPERED GLM (net "CHANGE then
PROCEED gated"): GLM's reframe only PARTIALLY holds —
- **Q1 cancellation = a REAL PROOF OBLIGATION, NOT a given:** `Z(w)` (Maxwell weight `pde.tex:357-416`) ≠ `W(w)` (brane kernel
  `:277,496-505`); a ratio of WEIGHTED AVERAGES `∫Z·K_stress/∫Z·K_source` does NOT cancel `Z` unless BOTH reduce to the SAME
  factorizable scalar `I_Z=∫√g_w Z dw` × SEPARATE field-content kernels (sources do NOT establish this). Can FAIL ⇒ W_eff back on.
- **Q3 λγ NOT just an evaluation:** `μ₀^eff=μ₀/Z_int` is a COUPLING normalization (`pde.tex:541-558`), NOT a photon-SPEED law;
  `c_bulk²=C_B/C_E` with the bulk→brane SPEED normalization UNSPECIFIED (pathA_20b). `Z_int` computable; `λγ` still needs the speed
  derivation. **Q2 χ_Q** = a real NEW linear-response DtN solve (not an existing extraction); `σ_Q^can=4a⁵/(27c_s⁵)`, fingerprint
  `a⁵/(27c_s⁵)` (GLM mis-cited). **Q4/m̂** CONFIRMED — reconcile first (may flip the tunability reading).
- **FINAL gated directive `directives/pathA_22b_minimal_combination_xi.md` (v2)** — design-review `SOUND-AS-IS`
  (`_scratch/codex_pathA_22b_v2_design_review.log`). VERIFIED GATE ORDER (each a HARD STOP, separate execution, honest outcomes all
  reachable): **GATE 0** (SymPy, cheapest/fail-fast) = 0a reconcile `m̂` (`MHAT_DIMENSIONFUL_CONFIRMED`/`..._DIMENSIONLESS_REFRAME`)
  + 0b prove-or-fail the `Z`-kernel cancellation (`CANCELS`/`DOES_NOT_CANCEL`; gates Gate 4) → **GATE 1** evaluate `Z_int` from
  existing `Z_w` (coupling artifact ONLY) → **GATE 2** derive `λγ` speed-reduction (`DERIVED`/`STILL_TUNABLE`/`CONDITIONAL_ON C_B/C_E`)
  → **GATE 3** `χ_Q` via linear outgoing-DtN BVP around the FROZEN existing branch (`=1`/`Δ_Q≠0`/`NEEDS_DYNAMIC_SOLVE`) → **GATE 4**
  field-content ratio `g_mhat²/g_G` ONLY IF 0b=`CANCELS` (else `BLOCKED_NEEDS_W_EFF`) → **GATE 5** assemble `P0·χ_Q·g_mhat²·λγ⁵/g_G
  ?=? 54/5` (`REAL_MATCH`/`REAL_MISS`/`FAIL_ABLE_PENDING_X`). `P0` = TARGET-BLIND input from the EXISTING finite-core branch (no
  re-solve). Discipline: target-blind, DERIVED-FORM gate, dimensional-check, dual-engine via `math -script` (NOT `wolframscript`=255),
  `timeout 600`, additive, no commit. Per-gate review = dimensional-fidelity + adversarial; Gate 5 → GLM tertiary; then user
  (incl. GATE-A freeze amendment when a derived value is installed).
- **▶ EXECUTE-NEXT (post-/compact):** write `_scratch/pathA_22b_gate0_execute_prompt.md` pointing at the directive's GATE 0 only →
  `codex exec --sandbox workspace-write -m gpt-5.5 -c model_reasoning_effort=high` as a `Bash run_in_background:true` task (NOT
  timeout-wrapped) → on exit, dimensional-fidelity + adversarial review (esp. "is 0b a real proof, can it emit DOES_NOT_CANCEL?")
  → bring Gate 0 outcome to user before Gate 1. Codex invocation pattern + all context in this block. [DONE → see (5l).]

**(5l) ⭐⭐⭐ RESUME HERE (2026-06-21, post-/compact) — GATE 0 EXECUTED + dual-reviewed + remediated. OUTCOMES STAND; NEXT ACTION =
EXECUTE GATE 1 (`Z_int` evaluation), USER-GATE PENDING.** Codex ran Gate 0 (`_scratch/codex_pathA_22b_gate0_execute.log`,
exit 0) → 2 clean review agents (dimensional-fidelity + adversarial) → Codex remediation (`_scratch/codex_pathA_22b_gate0_fix.log`,
exit 0). Artifacts: `src/stage1_solver/patha22b_gate0.py`, `tools/pathA_22b_gate0_crosscheck.wl`,
`tests/test_patha22b_gate0.py` (7 pass), `reports/pathA_22b_minimal_combination_xi.md` (Gate 0 section).
- **0a = `MHAT_DIMENSIONFUL_CONFIRMED`** (SOUND, non-rigged). Dimensionful `m̂0` (`[m̂0]=L⁻¹T⁻¹M⁻¹ᐟ²`) is the scale carrier;
  the dimensionless `m̂=1+O(a²/r²)` (`pde.tex:2095-2099`) reconciles as the dimensionless profile/source-map factor `g_mhat`
  (`m̂=m̂0·g_mhat`). `m̂0²·Γ5=G/c⁵` balances exactly (`Γ5=T⁵` via the (c_s/a)²-normalized dimensionless `P0` per pathA_22a — **no
  `[P0]` conflict**, the highest-value fidelity check). The dimensionless reading FAILS the odd law (`T⁵≠L⁻²T³M⁻¹`) → real
  discrimination, not tautology. **Does NOT flip pathA_22a** (`m̂0²` still the carrier; `g_mhat` an un-pinned dimensionless residual).
- **0b = `DOES_NOT_CANCEL (NOT_ESTABLISHED)`** (SOUND, CONSERVATIVE direction = forces the HARDER W_eff path, not a lazy out).
  Verified honest source-read: sources establish NEITHER route (a) shared factored scalar `I_Z` NOR route (b) pointwise-proportional
  kernels for the real `g_G`/`g_mhat²`. **⇒ W_eff/full transverse profile STAYS on the Gate-4 critical path** (Gate 4 NOT eligible
  for the reduced-ratio shortcut unless it proves route (a) or (b) from the action-level kernels).
- **Review found the outcomes sound + 5 fixes (all applied; outcomes UNCHANGED):** (1) ⭐ added BOTH Z-independence routes — proportional
  kernels ALSO null the Z-dependence of the ratio, so Gate 4 must test route (a) AND route (b), not just factorization; (2) annotated
  `0b` as `NOT_ESTABLISHED` (absence-of-proof, not proof-of-absence); (3) renamed `z_kernel_cancellation_proof`→`..._source_assessment`
  (Gate 0 surveys source availability; action-level kernel derivation deferred to Gate 4); (4) `.wl` proportional-control parity (now
  asserts route (b) Z-independence via a two-point Z → ratio=2); (5) pinned the `Γ5=T⁵` attribution to the (c_s/a)²-normalized `P0`.
- **▶ EXECUTE-NEXT:** Gate 0 work is UNCOMMITTED (offer to commit). On user GO → GATE 1 = evaluate `I_Z=Z_int=∫Z(w)dw` from the
  EXISTING exported `Z_w` (`m1c_background_export.py:166-170`), report as a COUPLING-normalization artifact ONLY (`μ₀^eff=μ₀/Z_int`),
  do NOT promote to `λγ` (that's Gate 2's speed-reduction). Same pattern: write `_scratch/pathA_22b_gate1_execute_prompt.md` → Codex
  `--sandbox workspace-write` background (NOT timeout-wrapped) → dimensional-fidelity + adversarial review → user before Gate 2. [DONE → see (5m).]

**(5m) ⭐⭐⭐ RESUME HERE (2026-06-21, post-/compact) — GATE 1 EXECUTED + dual-reviewed + remediated → `BLOCKED` but CONTAINED.
NEXT ACTION = EXECUTE GATE 2 (`λγ` photon-cone speed reduction — first load-bearing factor), USER-GATE PENDING.** Codex ran Gate 1
(`_scratch/codex_pathA_22b_gate1_execute.log`, exit 0) → 2 clean review agents → Codex remediation
(`_scratch/codex_pathA_22b_gate1_fix.log`, exit 0). Artifacts: `src/stage1_solver/patha22b_gate1.py`,
`tools/pathA_22b_gate1_crosscheck.wl`, `tests/test_patha22b_gate1.py` (15 pass w/ gate0), report Gate-1 section.
- **OUTCOME = `BLOCKED_NEEDS_DECAYING_Z_PROFILE_OR_FLOOR_PROVENANCE`.** Quadrature is FAITHFUL (exact arithmetic, genuine independent
  `.wl`), measure = flat `∫Z dw` per `pde.tex:289-295` (the `√g_w` variant is NOT computable from the export — flagged). But the
  exported `Z_w` is a PLACEHOLDER: `Z=floor+amp·exp(Gaussian)` with `floor=0.8` an UNDOCUMENTED solver-config constant
  (`coupled_branch.py:188-192`, differs across presets) with NO source support — `pde.tex:277` defines `Z(w)` LOCALIZED over
  `(-∞,+∞)`, which the floor violates (makes `∫` divergent). **72.9% of the finite-box `Z_int=2.03 L` is `floor×domain` (1.48 L);
  only 27.1% (0.55 L) is the localized Gaussian** → `Z_int` is domain-dependent, NOT a pinned constant. The `±0.0011 L` was grid-
  resolution only (relabeled); true uncertainty unbounded.
- **CONTAINED — does NOT gate the ξ verdict.** Both review agents independently confirmed `Z_int` is NOT a factor in
  `P0·χ_Q·g_mhat²·λγ⁵/g_G`; it feeds ONLY the coupling normalization `μ₀^eff=μ₀/Z_int`, `q_eff=q_*/√Z_int` (off the ξ critical
  path). So Gate 1 BLOCKED = a contained side-quantity, carry SYMBOLICALLY (never the numeric 2.03), program proceeds. Gate 1 also
  decisively confirmed `Z_int` ≠ `λγ` (the GLM-reframe conflation is dead). UNBLOCK (only if ever needed): export a genuinely decaying
  `Z(w)` or get floor provenance + physical `w`-extent.
- **▶ EXECUTE-NEXT:** Gate 1 work UNCOMMITTED (offer to commit). On user GO → GATE 2 = derive `λγ=c_γ/c_s` from the bulk→brane
  PHOTON-CONE SPEED map (NOT the coupling `Z_int`); pathA_20b left `c_bulk²=C_B/C_E` with the speed normalization UNSPECIFIED — close
  it or carry `C_B/C_E` explicitly; do NOT force `c_γ=c_s` (pathA_20b's negative control rejects it). Outcomes:
  `LAMBDAGAMMA_DERIVED`/`STILL_TUNABLE_LAMBDAGAMMA`/`CONDITIONAL_ON_<C_B/C_E or speed-map>`. Same pattern: write
  `_scratch/pathA_22b_gate2_execute_prompt.md` → Codex background → dimensional-fidelity + adversarial review → user before Gate 3. [DONE → see (5n).]

**(5n) ⭐⭐⭐ RESUME HERE (2026-06-21, post-/compact) — GATE 2 DONE + dual-reviewed + remediated → `STILL_TUNABLE_LAMBDAGAMMA`;
USER DECISION = PIVOT TO CALIBRATE-PREDICT (option B). NEXT ACTION = build the KNOB INVENTORY + observable→knob dependency map +
minimal calibration set (NOT grind Gate 3/4/5 ab-initio).** Codex ran Gate 2 (`_scratch/codex_pathA_22b_gate2_execute.log`) → 2
reviewers DISAGREED on severity (fidelity=PASS-disclosed vs adversarial=3 BLOCKERS: tautological `κ_B/κ_E==κ_B/κ_E`, over-labeled
CONDITIONAL, computable-but-skipped C_B/C_E) → Codex remediation (`_scratch/codex_pathA_22b_gate2_fix.log`) did the REAL `F²`
expansion. Artifacts: `src/stage1_solver/patha22b_gate2.py`, `tools/pathA_22b_gate2_crosscheck.wl`,
`tests/test_patha22b_gate2.py` (21 pass w/ gate0/1), report Gate-2 section.
- **OUTCOME = `STILL_TUNABLE_LAMBDAGAMMA`.** Flat/unwarped EM metric (`η_MN=diag(-1,+1,+1,+1,+1)`, `pde.tex:242`) ⇒ explicit
  `F_MN F^MN=−2F_ti²+2F_ij²` ⇒ `C_E=C_B=Z_int/μ0` ⇒ **`C_B/C_E=1` COMPUTED** (the fake 2nd knob collapsed). `λγ` rides on ONE
  source-unpinned dial: `β_bulk_to_brane` = the bulk-metric→acoustic-speed identification (`λγ=β·√(m_GNLS/(5Kρ0⁴))=β/c_s`). Negative
  control downgraded to a GUARD (was hardcoded-False). Gate 5 protected: `λγ` carried as explicit open knob ⇒
  `FAIL_ABLE_PENDING_LAMBDAGAMMA_SPEED_MAP`, never folded into `REAL_MATCH`.
- **STRATEGIC PIVOT (user, 2026-06-21):** Gates 0–2 cemented that the load-bearing ξ constants are NOT pinned by current sources
  (0b cancellation not established; 1 Z_int placeholder; 2 λγ tunable). User chose **option B = calibrate-predict** (the program's
  spine): calibrate a MINIMAL set on a trusted anchor, then PREDICT the held-out surplus; do NOT shoot in the dark deriving every
  constant. **User REQUIREMENT: every knob must be a NAMED PHYSICAL quantity** — pin down what we're calibrating, no fudge factors.
- **THE KNOBS (physical inventory) in `ξ·χ_Q·P0 = P0·χ_Q·g_mhat²·λγ⁵/g_G`:** `P0` = static quadrupole response of the throat
  (COMPUTED from the finite-core background, NOT a knob); `χ_Q` = the defect's outgoing-radiation coupling / "antenna" normalization
  (=1 only on the canonical free-space branch; actual branch unproven — would be Gate 3 DtN); `λγ=c_γ/c_s` = speed-of-light /
  speed-of-sound of the medium (EM flat ⇒ c_bulk=1; the c_γ↔c_s identification is the open dial); `g_mhat²/g_G` = ratio of
  gravity-analog source coupling to brane-scale coupling (needs `W_eff` transverse profile = Gate 4). ~3 physical numbers, each a
  property of the model/medium/defect SHARED across observables (which is what makes calibrate-predict have predictive content).
- **▶ EXECUTE-NEXT:** (a) commit Gate 2; (b) build the calibrate-predict setup: a KNOB-INVENTORY + OBSERVABLE→KNOB DEPENDENCY MAP
  (which held-outs — g−2, 5PN, multi-defect, ringdown — depend on which knobs + in what combination) + the MINIMAL CALIBRATION SET
  (anchors: Newtonian `G`, GR quadrupole `54/5`) → confirms whether the knobs are genuinely SHARED (calibrate-predict sound) and what
  the held-out predictions are. Verify the cross-observable sharing against pre-reg + project_defect_regimes before claiming surplus. [DONE → see (5o).]

**(5o) ⭐⭐⭐ RESUME HERE (2026-06-21, post-/compact) — calibrate-predict pivot EXECUTED: `decisions/14` built + β resolved + traps
hardened. NEXT ACTION = COMPUTE `χ_Q` (outgoing-DtN), prompt STAGED.** This block supersedes (5n)'s "build the inventory" next-action.
- **`decisions/14` (value-provenance + calibration map) BUILT, Codex-reviewed (6 fixes) + GLM-reviewed (§5).** Classifies every
  load-bearing value as DERIVED / INPUT (genuine calibration-necessity) / BRANCH-DETERMINED (computable, NOT free — compute don't
  calibrate) / BENCHMARK. **~7 genuine free inputs: `{n=5, K, m_GNLS, ρ0, μ0, q_*, λγ}`** (+`ħ` fixed external const + conventions).
  Everything in the ξ verdict is DERIVED (`c_s, a-cond, P0, C_B/C_E=1`) or BRANCH-DETERMINED (`χ_Q, λγ, g_G, g_mhat`).
- **GLM tertiary (decision-14 §5):** ADOPT — `g_G`,`g_mhat` are NOT independent (one `W_eff` computation; verdict sees only the ratio
  `g_mhat²/g_G`) ⇒ NEVER calibrate them independently (hidden-freedom trap). REJECT GLM's "cheap ratio via `Z(w)` cancellation" — our
  Gate 0b already found `DOES_NOT_CANCEL`; the ratio stays blocked on `W_eff` (GLM's KK intuition = a thing to TEST at Gate 4).
  REFRAME (user's point, vindicated): the dimensionless quadrupole `54/5` is closer to a PREDICTION than an anchor — calibration is for
  DIMENSIONFUL scales; "tunability" was the artifact of leaving branch-determined factors UNCOMPUTED.
- **β RESOLVED = `BETA_GENUINE_GAP`** (β-status analysis + adversarial review `VERDICT_SOUND`): the parent action treats EM as an
  INDEPENDENT gauge field on the flat bulk metric, NOT tied to the acoustic speed ⇒ `λγ=c_γ/c_s` is a genuine CALIBRATION INPUT (not a
  prediction, not a removable convention; nothing to compute — filling it = ADDING a postulate). **`c_s` is DERIVED (only `c_γ`
  free)**; `[λγ]=1` computationally confirmed in Gate 2 (units-restored, dual-engine). Counter-construction (`c_γ=c_s`) is the
  SUPERSEDED legacy `em_fields.tex`; canonical Path-A re-postulates EM independent.
- **TRAPS HARDENED (committed 857cb1f4):** pde.tex got a "Dimensional note" reconciling the m̂ inconsistency (law's m̂ = dimensionful
  carrier `m̂0`; source-map's `m̂=1+O(a²/r²)` = dimensionless profile factor; do not conflate) — no equation/physics changed. Report
  hardened with the em_fields-superseded provenance line.
- **THE UNDER-DETERMINATION (decision-14 §2):** ξ verdict has branch-determined unknowns `χ_Q, g_G, g_mhat` (+input `λγ`); the 2
  gravity anchors (Newtonian `G`, quadrupole) leave it under-determined until `χ_Q` (DtN) and the `g_G`/`g_mhat` pair (`W_eff`) are
  COMPUTED. The honest task is COMPUTE-don't-calibrate the branch-determined factors.
- **▶ EXECUTE-NEXT (post-/compact):** launch `_scratch/pathA_22b_chi_Q_execute_prompt.md` (Gate 3 = `χ_Q` via outgoing-DtN around the
  FROZEN finite-core branch; SCOPE-AND-SOLVE with 2 guards: NON-VACUOUS (must be able to return `χ_Q≠1`) + FEASIBILITY-HONEST (report
  `NEEDS_DYNAMIC_SOLVE`/`BLOCKED` if the frozen branch can't support a clean outgoing-DtN, `pde.tex:2845-2849`)) → `codex exec
  --sandbox workspace-write -m gpt-5.5 -c model_reasoning_effort=high` as a `Bash run_in_background:true` task (NOT timeout-wrapped) →
  on exit, dimensional-fidelity + adversarial review → user. Then `g_G`/`g_mhat` via `W_eff` (Gate 4, the hard blocker) remains.

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
