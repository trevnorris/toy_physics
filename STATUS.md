# STATUS — where the Path-A program is (single front door)

**This file is the canonical "you are here."** It is a thin pointer, not a copy — the detail lives in the linked docs.
Updated at every milestone (same moment `software/stage1_solver/decisions/13` §0 is updated). Last update: **2026-06-23.**

> **New to the model / need the physical picture? Read `docs/conceptual_foundation.md` FIRST.** It is the plain-language,
> native-terms statement of what the medium, the brane, the four sectors (gravity=drain, magnetism=swirl, electric
> charge=puncture direction ±w, light=brane shear), and the particle/defect (an extended throat — there are NO point particles)
> physically ARE. This `STATUS.md` is the *program state*; that doc is the *conceptual vision*.

---

## ⭐⭐ LATEST STATE (2026-06-24) — READ THIS + `conceptual_foundation.md` §4/§5/§6.1 FIRST

> **POST-/COMPACT WARNING TO FUTURE-ME:** the active frame is the **throat-soliton / 4D-light unification synthesis** worked out in
> conversation on 2026-06-24. It is *easy to forget after compaction* and the user assumes you still hold it — that mismatch has
> derailed us before. **Before any model work, READ:** `docs/conceptual_foundation.md` **§4** (throat/light), **§5** (multi-brane
> cosmology), **§6.1** (our own prior work) + the memory **`project-light-is-4d-throat-hypothesis`**.

**B4 RESULT:** `FAIL_NOT_CODIM1` (final, dual-engine verified, adversarial-SOUND). The GNLS **density-smectic brane is a candidate
NO-GO** — the medium's cubic invariant `U'''=15Kρ0²>0` makes a multi-`k` crystal beat the codim-1 lamella.

**R/C RESULT (2026-06-24, track 2b — DONE, tri-reviewed, dual-engine, earned): `RC_DENSITY_SMECTIC_LIGHT_NOGO`.** The density-smectic
brane is now a **confirmed NO-GO across ALL admitted drivers:** R→`FAIL_NOT_CODIM1` (bilinear/isotropic, cubic survives); Cdiv→
`FAIL_ADMISSION` (Goldstone integrate-out gives O(k⁰) directional EOS-stiffness softening for any λ_Cdiv≠0); Cpin→`FAIL_LIGHT_STARVED`
(χ_Cpin<0 DOES open a genuine codim-1 density-smectic window — the cubic IS escapable — but the same coupling pins P ⊥-plane, `P_∥=0`
at the static Ω_u=0 minimizer, starving in-plane light). **So density CAN form a brane but only by killing light** → strong earned
evidence for the synthesis: the brane must be a **light-confining shear surface, not a density modulation.** A first-class
falsification: it removes a rival and points the same direction as the throat-soliton picture. **pathA_25 density route CLOSED;** the
only remaining throat-soliton track is **(2a) the drain-sector derivation — PAUSED (user, 2026-06-24).** Artifacts:
`reports/pathA_25_gateRC_cubic.md` (+`_results.yaml`), `tools/pathA_25_gateRC_cubic_{sympy.py,.wl}`; Family-C `k→0` tautology RESOLVED.

**THE SYNTHESIS (the thing not to forget):** particles = **throat-solitons**. Light = an **intrinsically (3+1)D brane shear wave**
(2 polarizations; confined because the bulk is **shear-free** and light *is* shear) whose 4D-ness is the brane's **extrinsic
curvature** at throats — NOT a bulk excursion. A throat traps light as a **sub-cutoff standing wave** (mode larger than the neck →
evanescent → genuinely bound) = **MASS**. The trapped wave holds the throat **open** against **brane tension + superfluid
backpressure** (force-balance soliton, radius `R*(E)`); the puncture lets superfluid **drain into the bulk = GRAVITY** (gravity is the
**FLOW**, not the curvature); **±w puncture direction = CHARGE**. A **stack of branes** partitions the bulk into RS-like **finite
slabs** → the neighbor provides the drain's **return** (steady circulation) + sets gravity's range. **DEATHS to respect (don't
over-claim):** the trapped-wave lepton mass *tower* is FALSIFIED (1:9:25 ≠ 1:207:3477) → one soliton mass, **no** lepton spectrum;
B4 density-smectic FAILED.

**ALREADY BUILT — don't rediscover (`conceptual_foundation.md` §6.1):** Paper-7 `notes/inner_throat/` (force-balance + 4D GNLS solver,
STALLED on `L/a` self-selection); `notes/lepton_mass_notes.md` (trapped-wave functional, tower falsified); `software/stage1_solver/`
throat (a *different* object — frozen-geometry gravity-closure profile solve, abandoned/re-scoped). **Recurring wall: nothing in either
self-selects the throat size.**

**⭐ pathA_26 DERRICK CHECK — DONE (2026-06-24): `THROAT_DRAIN_DESTABILIZED`** (full gauntlet: Codex r1→r6 GREEN, GLM SOUND, Codex r7
confirm; executed dual-engine; arbiter-reran deterministic; adversarial `SOUND` by 10³–10⁴×; fidelity `ISSUES`→remediated `c_w²`-fix).
Artifacts: `software/stage1_solver/reports/pathA_26_derrick.md` (+`_results.yaml`), `tools/pathA_26_derrick_sympy.py`+`.wl`; prompt
`_scratch/pathA_26_derrick_execute_prompt.md` (v7). Two halves:
- **Phase A = `A_CANDIDATE_MIN` (the SOLID, real result):** the conservative throat-soliton EXISTS **generically** (interior PD min at
  `a*≈1.89,L*≈1.81`, virial residuals 0, both Hessian eigs >0, fluid-only collapses; adversarial confirmed 75/75 coupling combos). The
  trapped wave's outward pressure (`ω∝1/a`: squeeze the throat → freq↑ → energy↑ → pushes back) balances tension+vacuum inward. **NOT
  Derrick-forbidden.** The throat-soliton energy balance closes.
- **Phase C = `C_GENERICALLY_UNSTABLE` → top-line `THROAT_DRAIN_DESTABILIZED`. ⚠ DO NOT READ AS "PARTICLE KILLED" (user+orchestrator
  interpretation 2026-06-24, post-investigation):** the destabilization is a static divergence ONLY above a critical drain
  `gcrit≈0.006`; **for weak drain the throat is STABLE** (`gcrit>0`). The `DESTABILIZED` verdict came from a directive choice I made —
  demanding robustness over an O(1) drain box (`g` up to ~89) — that does **NOT reflect that gravity is the weakest force**. Prior repo
  work (investigation 2026-06-24) shows the drain/gravity sector is the runt: `E_w:E_f:E_PV=11:2:5` (feed ≈11% of rest E); `g_G`
  **calibrated not derived**, bare `P0~10⁻⁹` (~9 orders below GR target); analog horizon only at the **transonic** limit `v_b→c_s`;
  on the EXACT trapped mode the trans-brane current `J_w=0` (drain formally ~0). **So a particle's physical drain is tiny → it sits deep
  in the stable corner; the instability plausibly marks the BLACK-HOLE/transonic regime, where the drain SHOULD dominate (a consistency
  feature).** Your intuition ("you can hold it open; the tiny flow doesn't break it") is vindicated structurally.
- **THE PRECISE GAP (what was never done):** nobody mapped a particle's PHYSICAL drain strength onto pathA_26's dimensionless
  `g`-axis to prove `g_phys ≪ gcrit`. The pathA_26 drain law is **parameterized, not derived** (the `g≈89` is a synthetic worst-case,
  not a particle). Closing this rigorously = the drain-sector derivation. **Deeper open question (separate from stability):** with
  `J_w=0` on the exact mode, whether the "gravity = drain" mechanism produces a **nonzero, right-sized** drain at all is itself unbuilt
  (needs the dynamic AC→DC rectification the notes specced-not-built).

**▶ NEXT ACTION — ONE REMAINING TRACK (2b is DONE; see R/C RESULT above):**
1. **(2a) Drain-sector derivation — the REAL throat-soliton frontier (big, open) — ⏸ PAUSED (user, 2026-06-24); resume here next.**
   Derive the actual drain law (steady flow + density depletion + backpressure built INTO the equilibrium, not an external perturbation;
   the user's picture = wave-out vs tension+backpressure-in from the start), tied to the physical gravity coupling → (i) pin `g_phys` vs
   `gcrit≈0.006` rigorously, (ii) resolve whether the drain is even nonzero/right-sized (the `J_w=0` rectification question). Reuses
   `stage1_solver` infra + `reports/pathA_throat_solver_literature_synthesis.md`. Needs a SCOPING/design phase first (not a single execute).
2. ~~**(2b) R/C cubic verdict**~~ — ✅ **DONE: `RC_DENSITY_SMECTIC_LIGHT_NOGO`** (see R/C RESULT above). Density-smectic brane route
   CLOSED across all admitted drivers; Family-C `k→0` tautology resolved. Pre-registered prediction verified (R=no; C escapes the cubic
   but starves light). The brane is therefore a light-confining shear surface, not a density modulation — which is exactly 2a's object.

**COMMIT STATE (2026-06-24):** B4 artifacts committed `581c0b2f`; pathA_26 Derrick committed `8e0205c4`; **Gate R/C (2b) committed this
milestone** (`reports/pathA_25_gateRC_cubic.md` + `_results.yaml`, `tools/pathA_25_gateRC_cubic_{sympy.py,.wl}`, the §B.0 G0 tautology
fix, + these STATUS/ledger doc edits; `_scratch/` prompts+logs gitignored). **Commit only when the user asks; stage explicit paths.**

**Process discipline (unchanged):** Codex codes / Claude reviews; **dual-engine** (Mathematica: Codex needs `--sandbox
danger-full-access` — workspace-write CAN'T run it; OR the orchestrator runs `math` directly as arbiter); **review ordering** = iterate
Codex to GREEN → one GLM pass → fold → Codex to green; crux execution prompts get the full gauntlet; reports-only; `codex exec … -c
model_reasoning_effort=xhigh` backgrounded, never wrap the session in `timeout`.

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

> **PIVOT (2026-06-23): pathA_24 T1 FALSIFIED the little-arrows domain-wall brane → now the GNLS polar-smectic superfluid candidate
> + a consistency-gate program. See the "Next step" section below for the live state.** The pathA_23 / little-arrows history below
> is the context that led here. Conceptual home: `docs/conceptual_foundation.md` (v3) + `docs/medium_requirements_and_prior_art.md`.

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

## Next step — PIVOT to the GNLS polar-smectic superfluid + consistency-gate program (resume here after `/compact`)

**Read first:** `docs/conceptual_foundation.md` (v3 — §0.6 analog reframe, §2 v3 brane update, §8 new plan) + the full writeup in
`docs/medium_requirements_and_prior_art.md` (requirements list A/B/C, prior-art survey, candidate structure, consistency gates).

**What happened (2026-06-23):**
- **T0 ✅ (`f0c2745f`)** froze the GNLS + polar-OP action (`reports/pathA_24_T0_freeze.md`, hash `8fa41ac51e88`).
- **T1 ❌ FALSIFIED (`2fa91886`) — `T1_FAIL_NO_STABLE_WALL`, tri-reviewed GENUINE** (arbiter re-run both engines identical;
  `FIDELITY_CLEAN`; adversarial `T1_FAIL_GENUINE`). A *static* polar-vector domain wall has a **connected `S³` vacuum manifold
  (π₀=0)** → it spreads to infinite width (`σ_L→0`) and unwinds with zero barrier: no stable wall, no flat core, no confinement.
  The three-way no-win (emergent-`w` / stable-wall / light-capable-core) is now **demonstrated**.
- **Prior-art survey (5 agents)** + two user reframes → the pivot. Survey verdicts: kinematics-without-dynamics is a **universal
  wall** (analog gravity + Volovik ³He both stall there — same as our `g_G`-calibrated finding); the emergent-axis obstruction is
  **structural** (Davies–George–Volkas + ³He independently confirm T1); the **smectic mechanism** is the escape; light has two
  rival routes (continuum MacCullagh — ours, must beat a negative-energy instability — vs Wen lattice string-nets); charge-as-
  extended-puncture gets independent corroboration (Wen string-ends).

**The new candidate + plan (the GNLS polar-smectic superfluid):** KEEP the GNLS (gravity/magnetism/sound/`c_s`); ADD the polar
orientation field (light + charge) + a **non-local/roton layering driver** that gives a **smectic** (1D density-modulated) brane as
the spontaneous **ground state** (fixes T1's no-stable-wall + emergent-axis). **Density now modulates** (honest change). Test by
**hunting no-gos among the consistency gates** (analog, not derivation): **Gate L** (light on the layer — bounded-below MacCullagh,
no-longitudinal, leak-free — THE CRUX), Gate S (magnetism), Gate B (brane↔gravity compat), Gate Q (two charge signs), Gate K
(cone-lock `c_γ≈c_s`, likely a calibration gap), Gate T (throat/mass). **Inherited walls (concede):** dynamics/`G`/`α`
(calibrate-predict), emergent-axis/why-3D, Lorentz/preferred-frame (toy-analog scope).

**LIVE STATE (2026-06-23):** the directive is BUILT — `pathA_25` (v4, review-gauntlet SOUND: Codex design-review → GLM tertiary → 2
Codex confirm-passes), committed `6cdaa821`. **G0 (structure freeze) DONE** = `SECOND_MEDIUM_DRIFT_AT_FREEZE` (5 independent inputs:
`c_L1,c_L2` smectic driver + `μ_br,J_Pu,κ_Pu` light sector; tri-reviewed genuine; committed `77fd0e72`). Per §14 the drift is a
record-and-proceed finding (the 5 inputs = the calibration budget). **NEXT = Gate B4** (does the baseline Family-L driver make a stable
codim-1 emergent-axis smectic? — the T1-replacement): draft its execution prompt → design-review with **BOTH Codex (xhigh) AND GLM**
(crux gate, §18) → execute → tri-review → user gate. **Live ledger + resume block: `software/stage1_solver/reports/pathA_25_STATUS.md`.**
Methodology: "specify the FULL structure (postulated) + test CONSISTENCY / hunt a no-go," not "freeze minimal + test derivation."

**Methodology (locked):** Codex codes + runs + dual-engine; **AI prose never establishes a math fact** — orchestrator arbiter
re-run + transliteration-fidelity audit + adversarial review on clean agents; user gate per gate.

**Downstream / outstanding:** pathA_23 Stages 4–6 — likely re-framed by the smectic-brane result. (Deferred, parked: why fluid
flows *into* the mouth but leaks *back* into the brane — `decisions/15` §9.)

## Map — what you want → which doc holds it

| You want… | Look here |
|---|---|
| **The conceptual vision — what the medium / brane / 4 sectors / defect physically ARE (read first)** | `docs/conceptual_foundation.md` |
| **Requirements list + prior-art survey + the GNLS polar-smectic candidate + consistency gates (the live frontier)** | `docs/medium_requirements_and_prior_art.md` |
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
