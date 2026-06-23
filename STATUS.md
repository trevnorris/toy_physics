# STATUS — where the Path-A program is (single front door)

**This file is the canonical "you are here."** It is a thin pointer, not a copy — the detail lives in the linked docs.
Updated at every milestone (same moment `software/stage1_solver/decisions/13` §0 is updated). Last update: **2026-06-23.**

> **New to the model / need the physical picture? Read `docs/conceptual_foundation.md` FIRST.** It is the plain-language,
> native-terms statement of what the medium, the brane, the four sectors (gravity=drain, magnetism=swirl, electric
> charge=puncture direction ±w, light=brane shear), and the particle/defect (an extended throat — there are NO point particles)
> physically ARE. This `STATUS.md` is the *program state*; that doc is the *conceptual vision*.

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

## ⭐ ACTIVE FRONTIER (2026-06-23): EM re-founding → PIVOT to "why the brane exists" + the defect = a brane PUNCTURE (pathA_24)

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
  **not derivable medium-natively** — it stays a genuine free input. (`decisions/15` §15.) → **USER chose Path 1: postulate MacCullagh, CONDITIONAL.**
- **Stage 3 ✅** — `LEAK_BOUNDED_CONDITIONAL` (tri-reviewed; adversarial flagged too soft). No-leak holds only under an unmotivated
  `ε_leak≪1`/impedance price (else concept-fatal `FAIL_LEAK_BREAKS_MAGNUS`); postulated stress needs a postulated gyrostatic spin reservoir.
- **Stage 3b ✅ (verification)** — `OVER_COUNT_CONFIRMED_CURVATURE_LOCALIZED`: given the (separate-sector) action, `σ^R` is brane-internal
  not a bulk source; flat density-preserving light **free-slips** (inviscid bulk → free slip); leak is **curvature-localized** (∝|K|L_mix,
  far-field-vanishing). Retires the "intrinsic-to-light fatal leak" reading. CAVEATS (adversarial, fair): model-contingent (brane⊥bulk
  separation — single-medium vs membrane); defect/throat leak **relocated not retired** (throat solve still needed). (`decisions/15` §16.)

**Full physical picture (MacCullagh §11, λγ §13, costs §14, Stage-1 leak §7.1, Stage-2 crux §15, Stage-3/3b §16, brane-existence §17, defect=puncture §18) = `decisions/15`.**

## Next step — CONCEPT SETTLED (2026-06-23): test the "little-arrows" brane mechanism (the post-`/compact` plan)

The brane two-state mechanism is now an **agreed working hypothesis** (reasoned with the user): the medium's constituents are
**polar "little arrows"** that align (→ define `w`, → why space is 3D), flip across the wall (= our brane), **lie flat in-plane →
shear = LIGHT**, point `±w` in the bulk → shear-free → **magnetism**, and whose two `±w` domains give the **two charge signs**
(puncture direction, no winding). It dodges the single-well `U(ρ)` (the wall lives in the orientation field). Bonus hypothesis:
**dark energy** = a bulk↔brane cycle (matter drains brane→bulk; tension throttles a bulk→brane *areal* leak → accelerating
expansion; matter-vs-area scaling gives the decel→accel crossover). **Full picture + the test plan: `docs/conceptual_foundation.md`
(read first; §2 mechanism, §5 dark energy, §8 the T1–T5 ladder).**

**The plan (resume here after `/compact`):**
1. **Read `docs/conceptual_foundation.md`** top-to-bottom for grounding (esp. §2, §5, §8).
2. **Rework directive `pathA_24`** to encode the little-arrows mechanism + the **T1–T5 ladder** (T1 = does a polar field form a
   *stable* wall — the brane make-or-break, gates everything; T2 = does the wall carry 2-polarization in-plane shear = light;
   T3 = bulk stays shear-free; T4 = puncture→two charge signs + finite throat; T5 = the dark-energy cycle). The current
   Codex-SOUND `pathA_24` v1.1 is generic-domain-wall and must be **reworked to this concrete mechanism**, then re-run through the
   pipeline: Codex design-review (`xhigh`) → GLM tertiary → execute stage-by-stage, tri-reviewed (re-run + fidelity + adversarial)
   → user gate. **T1 FIRST** — if a polar field can't make a stable wall, the rest is moot and that is itself the result.
3. **Gate rule:** any postulated ingredient ⇒ CONDITIONAL; falsification is welcome; do NOT update `decisions/14`/papers without
   explicit user acceptance.

**Downstream / outstanding:** pathA_23 Stages 4–6 (spectrum/`u_w`, Maxwell/C5, charge/cone-`λγ`) — may be re-framed by pathA_24's wall
result; a GLM tertiary on the Stage-3b separate-sector question (likely subsumed by pathA_24's domain-wall derivation).
(Deferred, parked: why fluid flows *into* the mouth but leaks *back* into the brane — `decisions/15` §9.)

## Map — what you want → which doc holds it

| You want… | Look here |
|---|---|
| **The conceptual vision — what the medium / brane / 4 sectors / defect physically ARE (read first)** | `docs/conceptual_foundation.md` |
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
