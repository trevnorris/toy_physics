# STATUS — where the Path-A program is (single front door)

**This file is the canonical "you are here."** It is a thin pointer, not a copy — the detail lives in the linked docs.
Updated at every milestone (same moment `software/stage1_solver/decisions/13` §0 is updated). Last update: **2026-06-23.**

---

## Current state (one paragraph)

The toy model targets the GR-quadrupole verdict `P0·χ_Q·g_mhat²·λγ⁵/g_G = 54/5`. Under the **calibrate-predict** discipline
(every value DERIVED or a declared calibration gap — no silent knobs), the factors now stand: `P0` extracted, **`χ_Q ≈ 0.712`
COMPUTED** (Gate 3), `λγ` a genuine model gap (`BETA_GENUINE_GAP`), and — **Gate 4 (2026-06-22) → `GENUINE_BLOCKED`** — the
gravity ratio `g_mhat²/g_G` is **not derivable** from the current action (no target-blind source-map kernel; `α_J` doesn't
cancel; all 22 `m̂` sites in `pde.tex` are target-facing, dual-reviewed). So the model does not derive its own gravity coupling:
`g_G` is calibrated on Newton's `G`, and the **`54/5` quadrupole is an ABSORBED calibration anchor, not a prediction**. The
verdict closes only with the **EM-sector anchor** (which pins `λγ`) — now load-bearing. The falsifiable payoff is the **held-out
surplus** (g−2, 5PN, ringdown, multi-defect), riding the shared *derived* `χ_Q` + `P0/D0` bundle + `c_s`.

## ⭐ ACTIVE FRONTIER (2026-06-23): the EM re-founding — NOW EXECUTING (Stage 3 next; CONDITIONAL path)

Pinning `λγ` exposed that the EM sector had **drifted** from the single-medium concept (canonical EM = a fundamental gauge field
on a flat metric, *decoupled* from the medium; `reports/pathA_cgamma_of_rho_derivation.md`, Type-4). Re-founded EM **medium-native**:
- A single **scalar** superfluid **cannot carry transverse light** (scalar→spin-0; fluid→no shear). Accepted.
- Hypothesis (falsifiable): substructure → elasticity on surfaces; **our 3D space = an elastic brane**; **LIGHT = the brane's
  in-plane SHEAR waves**, shear **on the brane not the bulk** → bulk stays a pure fluid → magnetism (Magnus) preserved. Template:
  **MacCullagh rotational elasticity** (Cosserat/spinning-substructure mechanism). Three-sector: gravity=bulk flow,
  magnetism=bulk swirl, light=brane shear.

**Directive `pathA_23` v5 — THREE-WAY SOUND** (Codex design-review → GLM tertiary → Codex confirm). **NOW EXECUTING** stage-by-stage,
each tri-reviewed (re-run + fidelity + adversarial) before its gate:
- **Stage 0 ✅** — `NEW_PARENT_ACTION`; constitutive form **POSTULATED** (⇒ conditional-verdict rule active); dual-engine 25/25.
- **Stage 1 ✅** — `LEAK_CONDITIONS_DEFERRED` (after a rework that fixed a tautological first attempt). **KEY FINDING:** the brane↔bulk
  interface traction `T_na = T_wa + (T_ww δ−T_ab)∂u_w` is **generically transverse**; no-leak needs `T_wa=0` + isotropic `T_ab` at
  the brane, which are **NOT generic near a draining defect** → **the leak is EXPECTED; survival rides on the Stage-3 throat solve**
  (small magnitude OR non-fine-tuned projection cancellation). **Do not bank on no-leak.** (`decisions/15` §7.1.)
- **Stage 2 ✅ — THE CRUX** — `FAIL_UNSPECIFIED_SUBSTRUCTURE` (tri-reviewed; after a rework that fixed a tautological first
  attempt — the first try forced `μ_shear=0` by choosing a fluid EOS `W=W(J)` at the input). With the **deviatoric shear
  invariant present** and μ_br routed through a genuine able-to-fail classifier, **the medium's record does not determine μ_br**
  → **brane-shear EM is NOT derivable from the current single-medium spec.** Verified **trilemma:** μ_br=0→no light; μ_br>0
  (Cauchy)→light + a stray longitudinal "second photon" (`FAIL_CAUCHY_STRAY_LONGITUDINAL`, Stage 4); MacCullagh curl-only (only
  clean-transverse form)→reverse-engineered gyrostats + C5. ⇒ clean light needs an **extra postulated ingredient**. So `λγ` is
  **not derivable medium-natively** — it stays a genuine free input. (`decisions/15` §15.)

**Full physical picture + MacCullagh template + λγ subtlety + honest costs + Stage-1 leak (§7.1) + Stage-2 crux result (§15) = `decisions/15`.**

## Next step — USER CHOSE Path 1 (2026-06-23): proceed CONDITIONALLY, Stage 3 next

The crux is unmet as a *derivation*, so the route proceeds as a **conditional construction** (user-elected):
1. **Adopt the rotational/MacCullagh form as an explicit POSTULATE** (`U_∥=½μ_R(∇×u)²` — the only form that can give clean
   transverse light), **verdict CONDITIONAL flagged throughout**; the gyrostat substructure stays an acknowledged GAP, not a claim.
2. **Fire off Stage 3** (D1): constitutive no-leak closure — **decides the Stage-1 leak** — using the postulated stress/couple-stress
   + the Stage-1 bulk channels. Then tri-review → gate.
3. **Then:** 4 spectrum/`u_w` (the MacCullagh-vs-residual-longitudinal + `u_w`-fifth-force test) → 5 Maxwell/gauge (**the C5
   obstruction**)/Lorentz → 6a–d charge/energy/cone-`λγ`/leftovers. Able-to-fail at every stage (C5 / charge firewall / cone can
   each kill it). **Gate rule:** CONDITIONAL verdict; do NOT update `decisions/14`/papers without explicit user acceptance.

(Deferred, parked: why fluid flows *into* the mouth but leaks *back* into the brane — `decisions/15` §9.)

## Map — what you want → which doc holds it

| You want… | Look here |
|---|---|
| **The EM re-founding physical picture + MacCullagh template + Stage-1 leak finding (§7.1)** | `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md` |
| The EM re-founding execution plan (directive v5, EXECUTING) | `software/stage1_solver/directives/pathA_23_em_medium_native.md` |
| Stage 0 result (action + contracts) / Stage 1 result (leak audit) | `software/stage1_solver/reports/pathA_23_stage0_action_and_contracts.md` / `..._stage1_kinematic_leak_audit.md` |
| Full current state + resume-after-`/compact` pointer | `software/stage1_solver/decisions/13_emergent_constants_derivation.md` §0 |
| Every value classified DERIVED / INPUT / gap / benchmark | `software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md` (Gate-4 closeout = §6) |
| The defect-regime + held-out-surplus roadmap | `docs/defect_interaction_map.md` (CURRENT STATUS banner at top) |
| Per-gate derivation detail (Gates 0–4) | `software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md` |
| The directive that ran Gate 4 | `software/stage1_solver/directives/pathA_22b_gate4_ratio_or_blocked.md` |
| The calibrate-predict methodology | `software/stage1_solver/decisions/09_calibrate_predict_methodology.md`; `docs/pathA_preregistration.md` |

## The verdict equation (reference)

```
P0 · χ_Q · g_mhat² · λγ⁵ / g_G  =  54/5
 ✓     ✓     gap1     gap2  cal-on-G        (✓ = derived; gap1 g_mhat absorbs 54/5; gap2 λγ ← EM anchor)
G = (a·c_s²/m_GNLS)·g_G ,  m̂0 = (c_s/(a²·√m_GNLS))·g_mhat ,  c = λγ·c_s
```
