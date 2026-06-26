# STATUS — where the Path-A program is (single front door)

**This file is the canonical "you are here."** It is a thin pointer, not a copy — the detail lives in the linked docs.
Updated at every milestone (same moment `software/stage1_solver/decisions/13` §0 is updated). Last update: **2026-06-26.**

> **New to the model / need the physical picture? Read `docs/conceptual_foundation.md` FIRST.** It is the plain-language,
> native-terms statement of what the medium, the brane, the four sectors (gravity=drain, magnetism=swirl, electric
> charge=puncture direction ±w, light=brane shear), and the particle/defect (an extended throat — there are NO point particles)
> physically ARE. This `STATUS.md` is the *program state*; that doc is the *conceptual vision*.

---

## ⭐⭐ LATEST STATE (2026-06-26) — READ THIS FIRST

**DYNAMICAL-GRAVITY SECTOR = BUILT & GR-MATCHED; speed-of-gravity / aberration worry RESOLVED.** The conservative PN two-body ladder
**1PN→4PN + 2.5PN radiation is already derived, audited, and GR-matched** (calibrated / controlled-reduction) in `research/4d_*pn*`.
**DO NOT re-derive it.** (memory `project-pn-gravity-ladder`.)

**END GOAL = a fully CALIBRATED PDE delivering GR + EM.** Calibration is fine; first-principles is NOT required; **existing-in-any-shape
is the win.** The central spec is `research/pde_ledger/` (the 253-stage audited ledger of the target moving-throat PDE); every
calibration is a constraint that feeds it. (memory `project-calibrated-pde-goal`.)

**pathA_28 (monopole/dipole radiation) — DONE = `MONOPOLE_DIPOLE_RETURN_CONDITIONAL`.** A VERIFIED CONSTRAINT-SPEC (dual-engine;
arbiter PASS + fidelity CLEAN; adversarial CONCERNS = it is a constraint-spec, not a falsifiable test). **Handoff:** to avoid
GR-forbidden monopole/dipole gravitational radiation, the brane↔bulk return must deliver `R0 = −M0` (net mass-rate `M0 = ∫S_leak`)
and `R1 = −D1` (net dipole/momentum-rate `D1 = ∫x_i S_leak + ∫j_i`, including the carried odd wake). Artifacts under
`software/stage1_solver/` (tools / reports `pathA_28_monopole*`).

**⭐ TRACK 3 GATE-1 (brane↔bulk return, `pathA_29`) — DONE + VERIFIED (2026-06-25) = `RETURN_RESIDUAL_PREDICTION`.** The keystone
pathA_28 handed off — earned on the **4th execution (v3b)** after a full tri-review (orchestrator arbiter re-run reproduced +
adversarial CLEAN + fidelity FAITHFUL). **Given the drain premise (`Z<0` = gravity IS the inflow): 1/r² Newtonian gravity SURVIVES
the finite slab** — both admissible DC-sink return completions (de-structuring/absorbing + Bloch stack) genuinely solve a normalizable
`m=0` transverse zero mode → `p=2` via a real 3D-radial `dsolve` (counterfactual-guarded: a wrong `1/r⁵` → nonzero residual,
rejected). **And the drain comes bundled with an UNAVOIDABLE bounded monopole/dipole `c_s`-radiation residual `∝ ε0 = 1−𝒯₀(0)`**,
tied to the gravity strength — **the falsifiable departure from GR** (GR forbids monopole/dipole GW via Birkhoff/mass-conservation;
the drain breaks brane mass conservation). NOGO is genuinely reachable via a derived delocalizing warp (`p=3`). **Sharpens but does
NOT close `pde_ledger` open-item #9** (records the residual-radiation prediction); the **gravity-range (1/r²) item passes** for the
localizing flat-slab family. Artifacts: `software/stage1_solver/{directives/pathA_29_brane_bulk_return.md (v3), tools/pathA_29_*,
reports/pathA_29_*}`. **⭐ THE ACTIVE PUSH (user, 2026-06-25: "push until a wall"): complete the full nonlinear moving-throat PDE / brane↔bulk return
closure.** Central spec = `research/pde_ledger/` (253 stages, ALL algebra done; single OPEN item = **"actual branch realization"** — the
solved nonlinear PDE must RETURN the grouped-`P2` / quadrupole-normalization data `m̂₀²P₀ = 54Gc_s⁵/(5a⁵c⁵)`). **This is a ~6-gate
ladder, NOT one step** — master checklist = `research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md`.

- **Step 1 (reconcile the two open-item framings) ✅ DONE** = `moving_throat_pde_open_item_reconciliation.md`: pathA_29's "#9 /
  `R0=−M0`" (ℓ=0/1 forbidden channels) and the ledger's "branch realization" (ℓ=2 quadrupole) are **complementary multipole sectors of
  one return law**, not the same and not independent; a new cross-ℓ consistency gate falls out (Gate 5).
- **Gate 1 (`pathA_30`, frozen-wall D/N unit test) ✅ DONE = `DN_UNITTEST_BC_DEPENDENT`** (tri-review CLEAN: arbiter re-run both engines
  + fidelity + adversarial). The geometry lift's frozen reduction genuinely gives the const-coeff Helmholtz resonator (DtN
  `−(ω/c_s)tan`, half-shifted ladder, Robin counterfactual, dual-engine agreed; Bogoliubov `k⁴` deferred under `kξ≪1`). **D/N is a
  banked calibration input** (`bc_derivation_emitted=False`); earning it → `PASS` is an optional later upgrade. 2 non-blocking NITs
  logged in the ladder. Artifacts: `software/stage1_solver/{directives,tools,reports}/pathA_30_*`.
- **Gate 2 (`pathA_31`, scalar breathing `η₀₀`) ✅ DONE = `BREATHING_CALIBRATED`** (tri-review CLEAN after a remediation). The
  distributed ℓ=0 reduction genuinely reproduces the legacy `(a,L)` collective closure: derived `L₀`-harmonic profiles, `M_AB/K_AB`
  from operator projection (not `∂²E_geom`), two genuinely-independent HF routes match, structure gate computed + able-to-fail,
  dual-engine agreed. The 2-mode truncation is clean (`V₂`-overlap `o_1=0.993` at the physical `β_L0=1.85`) across the order-unity
  wall-stiffness window (`K_η/T_w≲2.6`), failing only for sharp walls. CALIBRATED ⇐ `μ_η,T_w,K_η` are calibration inputs. **Caveat:**
  the overlap doesn't guard profile-correctness (HF mismatch does); clean truncation shown *for assumed order-unity wall stiffness*.
  v1 was tri-review-REJECTED (HF `x−x`, hardcoded counterfactual, gamed threshold) → remediated. Artifacts: `pathA_31_*`.
- **Gate 3 (`pathA_32`, grouped-`P2` / ℓ=2 sector) ✅ DONE = `ISOTROPY_CALIBRATED`** (full tri-review CLEAN after a remediation). The
  distributed lift's ℓ=2 grouped-`P2` sector satisfies the isotropy / lane-degeneracy theorem at the linearized isotropic reference:
  the three grouped lanes `{20,21,22}` collapse to common conservative coefficients (raw-`D` lane defects = 0, the **PRIMARY** gate),
  the reduction is SO(3)-covariant (angular Gram = `I₅`; computed `−Δ_S²` eigenvalue `λ_m=6` per harmonic; the angular stiffness
  `K₂=∫[T_w β₂'² + (K_η+6T_Ω)β₂²]` — the term that dropped at ℓ=0 — is now alive), and the gate is genuinely able-to-fail (8
  counterfactual probes, each flips under ablation). **The quadrupole sector first appears here.** CALIBRATED ⇐ the wall constants
  `μ_η,T_w,K_η,T_Ω,β₂` + the symbolic radial scalars are calibration inputs. **Two-tier gate:** raw-`D` lane collapse is the
  verdict-bearing PRIMARY test; the published `a₂=b₂=a₄=b₄=0` (normalized `u`-defects) is a necessary-but-not-sufficient cross-check
  (a per-lane order-independent prefactor cancels in the normalized response). v1 was tri-review-REJECTED — the
  adversarial-with-ablation leg caught two pass-by-construction defects fidelity missed (dual-engine gaming = vacuous `x−x` on three
  byte-identical Mathematica copies; 5/8 probes typed their FAIL booleans + a tautological able-to-fail aggregate) → remediated
  (honest per-harmonic `.wl` + per-lane `D` cross-engine; each probe computed-from-mutated-input with a self-ablation) → re-verified
  (arbiter byte-for-byte; fidelity + adversarial re-ablation EARNED). Artifacts: `pathA_32_*`.
- **Gate 4 (`pathA_33`, quadrupole `54/5` normalization) ✅ DONE = `QUAD_CALIBRATED`** (full tri-review CLEAN after a remediation). On
  Gate-3's isotropic outgoing branch, the ℓ=2 sector delivers the **EARNED predictive surplus cleanly SEPARATED from the CALIBRATED
  magnitude**: the outgoing fingerprint `1/9, 4/81, 1/27` is **DERIVED** from the DtN `Ŷ₂^out=−3/Λ₂^out` Hankel series (not hardcoded),
  the prefactor algebra follows from `P(ω)=D₀N(ω)/D^cons(ω)²` (the `−2D₂N₀` term), `P0_target_scaling=a⁻⁵`, and the derived `χ_Q=1`
  closes `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵ ⟺ γ_quad^eff=2G/5c⁵`. **The headline result — the earned/calibrated split:** `54/5 = 2·27/5` — the
  `27` is **EARNED** (`derived_in_gate`, from the fingerprint), the `2/5` + `G` + assembled magnitude are **CALIBRATED**
  (`external_bridge_input`), classified by a 4-way **PROVENANCE** partition (NOT `G→λG` invariance, since `54/5` is G-invariant yet
  calibrated). `G` is `GENUINE_BLOCKED` — the PDE delivers the FORM/branch, **not** Newton's `G`. **Dimensional milestone:** a genuine,
  `μ̂₀`-FREE, able-to-fail dim check (`[P₀^phys]=(c_s/a)²·N₀/D₀` must be dimensionless from sourced `[N₀]=L⁻¹M`, `[D₀]=L⁻¹T⁻²M`) caught
  that the handoff's `P₀=N₀/D₀` silently drops a `(c_s/a)²` factor (natural-units trap). v1 was tri-review-REJECTED — the
  adversarial-with-ablation leg caught two pass-by-construction sub-controls fidelity + arbiter both passed: (a) the dim gate was STILL
  tautological (it **back-solved the FREE carrier `μ̂₀`** to force homogeneity); (b) the per-probe `self_ablation` was a **constant, not a
  re-run** → remediated (`μ̂₀`-free dim gate + corrupt-`[N₀]` probe §3d′ + real two-verdict `self_ablation`) → re-verified EARNED
  (corrupting `[N₀]` now fires `FAIL_DIMENSIONAL`; `[G]` doesn't). Artifacts: `pathA_33_*`.
- **⭐ NEXT = Gate 5 (scalar/dipole side-conditions + cross-ℓ unification):** from the *same* PDE, reconcile the ℓ=0/1 brane↔bulk return
  (pathA_28/29 — `R0=−M0`, `R1=−D1`, bounded residual `∝ε₀`) with the ℓ=2 quadrupole channel (the new cross-ℓ consistency gate). Then
  Gate 6 (full nonlinear closure = the WALL — where the symbolic branch data + the literal `54/5` magnitude must be RETURNED by the
  solved PDE); PN match-back threaded through = the decisive falsifier. Don't re-derive the audited PN ladder (`research/4d_*pn*`).
- **⚠️ Tracked debt (Gate-1–3 dimensional checks):** the `pathA_30/31/32` dim checks were found **VACUOUS** (typed-tuple ledgers
  disconnected from the real expressions) — a retrofit debt to redo after Gate 4 using the real `pathA_18` harness (the same harness
  the Gate-4 `μ̂₀`-free check now uses). Physics OOM-risk judged low.

Process lesson banked (memory `project-model-mechanics-corrections` + `feedback-offload-review-gauntlet`): every compute gate runs the
full Codex→GLM→Codex directive gauntlet + dual-engine + tri-review, with the **review rounds offloaded to a gauntlet-runner agent** so
the orchestrator context survives. pathA_30 took Codex(8)→reconfirm→GLM(Bogoliubov `k⁴`+5)→Codex post-GLM(2)→close-out before a clean
execution.

**Model-mechanics reminders** (memory `project-model-mechanics-corrections`): nothing is static; **three speeds** — `c_s` = speed
gravitational changes propagate (∝ρ²); `v_r` = field *strength*, not a speed; `c_γ` = light. Gravity = the **flow between drains**;
changes propagate at `c_s`, and uniform motion tracks the **current** position → **no aberration**. Throat-soliton has **no sloshing**
(`J_w=0` expected; AC→DC retired); gravity is a separate background de-structuring drain.

**COMMIT STATE (2026-06-26):** pathA_28 (`8cf6f1f1`), pathA_29 (`145c8426`), pathA_30 Gate-1 (`f460fc63`), pathA_31 Gate-2
(`765db5f0`), the prior docs-sync (`810e01f7`), and **pathA_32 Gate-3 (`6711167a`)** committed. The **pathA_33 Gate-4 milestone**
— directive, tools, dual engines, reports + this STATUS / ladder / decisions-13 / conceptual-foundation sync — is **committed as
this HEAD commit** (hash pinned at the next sync, per the `810e01f7` convention). **Commit only when the user asks; stage explicit paths.**

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
