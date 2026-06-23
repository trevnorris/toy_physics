# Directive pathA_24 — WHY the brane exists + the defect as a brane PUNCTURE (charge/mass/throat unification)

**Status:** DRAFT v0 (roadmap, written 2026-06-23 before `/compact`). **NOT yet executable** — must go through the standing process:
Codex design-review (`-c model_reasoning_effort=xhigh`) → apply fixes → GLM tertiary → apply fixes → Codex confirm → execute
stage-by-stage, each tri-reviewed (orchestrator re-run + fidelity audit + adversarial review on clean agents) before its user gate.
**Owner:** orchestrator (Claude). **Engine:** Mathematica leads (domain-wall profiles, zero-mode spectra, surface stress tensors,
puncture/throat energetics, Coulomb-field correspondence), SymPy cross-checks; ≤2 concurrent `math -script`.
**Resume/context:** `decisions/15` §17 (brane=domain wall) + §18 (defect=puncture: charge/mass/throat) — the physical picture this
directive tests; program state `decisions/13` §0; front door `STATUS.md`. Memory `[[project-brane-existence-defect-structure]]`,
`[[project-geon-throat-hypothesis]]`, `[[project-em-brane-shear-picture]]`.

---

## Decision + why now
pathA_23 reached two converging walls: **Stage 2** (the medium does not *derive* the brane shear — `FAIL_UNSPECIFIED_SUBSTRUCTURE`;
proceeding only on a *postulated* MacCullagh law ⇒ CONDITIONAL) and **Stage 3/3b** (no-leak holds only for a *separate-sector* model;
is that single-medium-legitimate?). Both reduce to one unanswered foundational question the user has intentionally left open until now:
**why is there a 3D brane in the 4D bulk, and what keeps us on it?** Answering it (user's surface-tension instinct → a **domain wall /
phase interface**) plausibly **derives** the brane, its **internal structure** (→ whether it carries shear = light, retiring the
Stage-2 postulate), the **separate-sector legitimacy** (the wall is a genuine emergent object), AND — the prize — a **fundamental model
of the charged massive particle**: a **puncture** through the tense brane, with **charge = the topological puncture**, **mass = the
trapped geon standing wave**, and **throat size = the tension-closing vs holding-open balance** (`decisions/15` §17–§18).

## Falsification stance (load-bearing — [[feedback-falsification-is-the-goal]])
A genuine test. Plausible deaths: (DA1) the medium cannot support a stable codim-1 wall without ad-hoc structure; (DA2) the wall is
structureless ⇒ fluid membrane ⇒ tension+bending but **no shear** (Stage-2 problem persists — light still not derived); (DA3) the wall
binds no zero modes ⇒ no confinement; (DB1) the puncture's tension-energy does **not** reproduce the Coulomb field / a sensible charge;
(DB2) the tension-vs-holding-open balance gives **no finite** throat size (self-energy still divergent or unstable); (DB3) charge comes
out **mass-dependent** (breaks universality) or the charge⊥mass decoupling fails. Any terminal FAIL is first-class and reported plainly,
never rescued. A suspiciously clean "it all derives" is re-checked HARDER. Honest expectation: this may show the single-scalar GNLS needs
**more structure** — that is itself a valuable, decision-relevant result about what the single-medium concept actually requires.

## Honest classification
Likely a **`NEW_PARENT_ACTION` refinement** (adding the structure that lets the medium self-localize a wall) — do NOT call the brane
"derived" until a mechanism is shown without reverse-engineering. Same provenance discipline as pathA_23: any added structure labeled
derived(independently-motivated) / postulated / ad hoc; a postulated ingredient ⇒ CONDITIONAL verdict, no `decisions/14`/paper updates
without explicit user acceptance. Carry the §14 conceptual costs (classical scope; substructure beneath the mean-field).

---

## PHASE A — Why the brane exists (the domain-wall / confinement derivation)
**A1. Can the current medium self-localize a brane?** Confirm the gap: `U(ρ)=(K/4)ρ⁵` is single-well ⇒ no domain wall from the scalar
alone (the brane is currently imposed via `V_conf`/`Z`/`W`/`B_ℓ`/`k_w`). Pin precisely what is imposed vs. derivable. Outcomes:
`BRANE_IMPOSED_NOT_DERIVED` (expected baseline) / `BRANE_SELF_LOCALIZES_AS_IS` (surprising — re-check hard).
**A2. Minimal structure for a stable wall + tension.** Find the least-added structure giving a stable codim-1 wall with a finite
**surface tension** `σ_brane`: degenerate vacua (double-well in an added/!modified field) vs. a 2nd order-parameter component
vs. self-trapping by the collective drain network. Classify each by parsimony + single-medium fidelity (avoid dualistic drift,
[[project-single-medium-concept-vs-math-drift]]). Derive the wall profile + `σ_brane`. Outcomes:
`WALL_DERIVED(mechanism, tension)` / `FAIL_NO_STABLE_WALL` / `WALL_NEEDS_ADHOC_STRUCTURE`.
**A3. Confinement (why we don't flow into the bulk).** Show the wall is a potential well binding **zero modes** at `w=0` (us); the cost
to leave = the gap. Connect to the existing `V_conf`/`Z`/`W`/`B_ℓ`/`k_w` as the *effective* description the wall would derive. Outcomes:
`CONFINEMENT_FROM_WALL` / `FAIL_NO_BOUND_ZERO_MODES`.
**A4. Does the wall carry IN-PLANE SHEAR? (the Stage-2 link — make-or-break for deriving light).** Compute the wall's surface
constitutive content: structureless wall = tension + bending only (`FAIL_FLUID_MEMBRANE_NO_SHEAR`, = Stage-2 problem persists);
structured wall (internal director/broken-symmetry/texture) = test whether in-plane **shear modulus `μ_brane`** emerges. If shear is
**derived** from independently-motivated wall structure ⇒ this **retires the Stage-2 postulate** (CONDITIONAL → derivation). Outcomes:
`WALL_DERIVES_SHEAR(μ_brane)` / `FAIL_FLUID_MEMBRANE_NO_SHEAR` / `SHEAR_STILL_POSTULATED(→CONDITIONAL)`.

## PHASE B — The defect as a brane puncture (charge / mass / throat)
**B1. Puncture energetics + throat size from the balance.** Model a puncture (throat) through the tense wall; compute the
surface-tension **restoring (closing) stress/energy** vs. the **holding-open** agent (geon standing-wave pressure and/or drain flow).
Solve the **balance** → the **finite throat radius `a`** and throat profile. Outcomes: `THROAT_SIZE_DERIVED(a)` /
`FAIL_NO_FINITE_BALANCE` (collapses or runs away).
**B2. Charge = the topological puncture (universality + firewall).** Show charge is carried by the puncture topology (that it punctures
+ winding/orientation), **mass-independent** → universality; reconcile with the existing `η_Q e_*` ontology + charge firewall
(`pde.tex:279-312`). Outcomes: `CHARGE_IS_PUNCTURE_TOPOLOGY` / `FAIL_CHARGE_MASS_DEPENDENT` / `FAIL_CHARGE_FIREWALL`.
**B3. Electric field energy = brane-tension deformation energy.** Test whether the tension/deformation field around the puncture
reproduces the **Coulomb `1/r²` field** and a **finite self-energy `~q²/a`** (using `a` from B1) → candidate **resolution of the
classical self-energy divergence**. Outcomes: `TENSION_REPRODUCES_COULOMB` / `FAIL_NOT_COULOMB` / `SELF_ENERGY_STILL_DIVERGENT`.
**B4. Charge ⊥ mass decoupling.** Test the user's assumption that the geon (mass) and the puncture-tension (charge) sectors do not
directly interact *at this level*. Outcomes: `DECOUPLING_HOLDS` / `FAIL_CHARGE_MASS_COUPLED` (then characterize the coupling — possibly
a feature, not a death).
**B5. Payoff back to the program.** If A4 derives `μ_brane` and B sets the throat: feed `c_γ²=μ_brane/ρ_brane` + brane tension into
**`λγ`** and the GR-quadrupole verdict bookkeeping; update the value map ONLY after a full pass + explicit user acceptance of any
remaining CONDITIONAL status. Connect the throat structure to the held-out surplus (g−2, particle-mass relations).

**B6. CALIBRATE-PREDICT: brane tension ↔ electric charge (the runway — user, 2026-06-23).** The high-leverage play: the EM coupling is
currently un-derived (Gate 4: gravity `g_G` absorbs `G` with zero surplus; `λγ` is a gap). If the puncture's tension-energy ties to
charge, **calibrate the single, universal brane tension `σ` on the measured `e`/`α`**, then predict the held-out surplus.
- **FIRST GATE — surplus accounting (the Gate-4 lesson, non-negotiable).** Derive the tension→charge map and **count free constants vs.
  independent downstream predictions.** If it carries a free O(1) constant per prediction ⇒ calibrating `σ` on `e` merely ABSORBS `e`
  (zero surplus = the `g_G`-absorbs-`G` trap) ⇒ report `EM_COUPLING_ABSORPTIVE_NO_SURPLUS` plainly. Constant-free / fewer constants than
  predictions ⇒ `EM_COUPLING_PREDICTIVE(surplus=N)` ⇒ runway. This ratio IS the test — settle it before any calibration number is quoted.
- **Leverage prerequisite:** the cascade (one `σ` → many predictions) requires Phase A's wall structure to **relate the moduli**
  (`σ`, `μ_brane`→`c_γ`, `ρ_brane`); three independent knobs ⇒ no cascade, no runway.
- **Charge quantization/universality** (same `e`, mass-independent) must **emerge** from the puncture topology — a prediction, not a fit
  (`FAIL_CHARGE_MASS_DEPENDENT` if it doesn't, cf. B2).
- **Held-out checklist (target-blind; check AFTER calibrating `σ` on `e`):** classical electron radius `r_e≈2.8 fm` (does throat size
  `a` land there? *which* scale it picks — `r_e` vs Compton `λ_C`, differing by `α` — is itself diagnostic); `λγ=c_γ/c_s≈1`
  (GW170817 cone-lock → **closes the verdict gap**); lepton mass ratios (throat-size↔mass via the balance); `g−2`; **stretch** = `α≈1/137`
  as a brane/bulk moduli ratio (moonshot — only meaningful if the map is constant-free; do NOT bank).
- **Units restored throughout** (the brane tension is a 3D energy density — a wall in 4 spatial dims is a 3-surface; the `σ ↔ ε₀/e²` map
  needs Jacobian/volume factors exact, cf. the `m̂0²·S_port` saga). Requires the **charge⊥mass decoupling** (B4) to calibrate cleanly.
  Outcomes feed B5; CONDITIONAL gate rule applies.

---

## Discipline (binding, same as pathA_23)
Target-blind / able-to-fail at every step; never soften/rescue a FAIL; a clean derivation re-checked HARDER. No tautology / no
reverse-engineered structure (a wall structure chosen *only* to give shear/charge is circular). Honest provenance + conditional-verdict
rule. Dimensional consistency with units restored. Dual-engine symbolic checks + transliteration-fidelity audit after every
computational step. Codex designs routes + writes/runs all scripts; Claude reviews only. No paper / `decisions/*` edits during execution
(reports → `reports/`; orchestrator owns canonical docs). Scripts ≤10 min with `timeout 600` on the scripts (never the Codex session).

## Deliverables
Per stage: `reports/pathA_24_<phase><n>_*.md` (verdict-first) + Mathematica `.wl` (+ SymPy cross-check). Final (gated, only after a full
pass + conditional-status acceptance if applicable): the brane-existence + defect-structure result with honest provenance; updates to
`decisions/15`, `decisions/14`, `decisions/13` §0, `STATUS.md`; limitations ledger.

## Review plan
1. Codex design-review of this DRAFT (is it able-to-fail, non-tautological, non-circular, correctly ordered; are DA/DB deaths reachable).
2. GLM tertiary (foundational directive — the brane-existence + charge/mass unification is concept-level). 3. Execute stage-by-stage
with per-stage tri-review + user gate.
