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

## ⭐ ACTIVE FRONTIER (2026-06-23): the EM re-founding — NOW EXECUTING (Stage 2 next)

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

**Full physical picture + MacCullagh template + λγ subtlety + honest costs + Stage-1 leak finding = `decisions/15`.**

## Next step (post-/compact) — USER CHOSE (A): Stage 2, the constitutive crux (TRACK B)

1. **Fire off Stage 2** (D2): the brane constitutive form. **Track B** — Codex MUST attempt an *independently-motivated*
   microstructure derivation (rotational/MacCullagh vs Cauchy vs Cosserat); postulate only as honest fallback ⇒ CONDITIONAL.
   (Fold 2 small Stage-1 report honesty notes into the exec prompt.) Then tri-review → gate.
2. **Then:** Stage 3 (constitutive leak closure — **decides the Stage-1 leak**) → 4 spectrum/`u_w` → 5 Maxwell/gauge (C5)/Lorentz
   → 6a–d charge/energy/cone-`λγ`/leftovers. **Gate rule:** a POSTULATED form ⇒ CONDITIONAL verdict; do NOT update
   `decisions/14`/papers without explicit user acceptance.

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
