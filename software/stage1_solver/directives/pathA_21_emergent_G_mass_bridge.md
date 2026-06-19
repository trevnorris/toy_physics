# Directive pathA_21 — Emergent constants, Step 3: `G`, the mass-bridge, and the `m↔G` unification

**Status:** ⚠️ DRAFT — to be FINALIZED after `pathA_20` (`c_s`+`c`) lands + is reviewed. It CONSUMES `pathA_20` outputs:
the derived `c` (= `c_γ`) and the `c_γ/c_s` result (closed-form OR `C_GAMMA_RATIO_UNDERDETERMINED`); the `ħ`-provenance
verdict (`HBAR_EMERGENT` / `HBAR_FUNDAMENTAL` / `HBAR_PROVENANCE_UNDETERMINED`) + the `h`/`2π` decomposition assessment;
and the **`flux_law_verdict`** (`TRANSONIC_CHOKED` with `J_crit` / `NONTRANSONIC_NO_CHOKED_FLUX` with the alternate law /
`STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`) — NOT an unconditional `J_crit`; each case routes differently below.
Do NOT fire until `pathA_20` is gated. This is where the `pathA_19` honest negative (`INFLOW_MASS_SOURCE_MISSING`, base
retained `{L,T,M}`) gets RE-TESTED with the now-available ingredients. If `pathA_20` returned the underdetermined
verdicts, this step inherits those as blocking residuals (it cannot manufacture the missing profile/operator data).
**Date:** 2026-06-19
**Owner:** Codex (DERIVES + codes; iterates until scripts exit 0). Claude reviews afterward.
**Trigger:** decision-13 §4 item 3, Step 3 of 4. Chain: `pathA_19` (foundation) → `pathA_20` (`c_s`+`c`) →
**this (`G` + mass-bridge + `m↔G`)** → `pathA_22` (scale-map → `m̂0²·S_port` → re-run B2c).

## Why this step (context)
Two things converge here. (1) `pathA_19` could NOT collapse `M` because doing so required `c` (circular) and a verdict
on `ħ`; `pathA_20` supplies both. So the mass-as-inflow hypothesis becomes TESTABLE for the first time. (2) The user's
ontology: a defect's mass is the attractive pull it exerts via condensate inflow — so **the inter-defect force, the
defect mass, and `G` are three faces of the same inflow `J`**. Derive them as such, or report the blocking gap.

## Scope & stance
DERIVATION. No model-formula changes, no freeze touch, no `m̂0²·S_port` un-pin (that is `pathA_22`). Extend
`dimensional_check.py` (new group; leave pathA_18 + pathA_19 + pathA_20 groups intact). Same infra constraints as
`pathA_19`/`pathA_20` (read-only sim dir; no `$RT exec`; `timeout 600`; ≤2 MMA seats; YAML/md). **Dual-engine required**
where MMA can verify the force law, the `G` reduction, and the dimensional algebra. Use the `pathA_19` base + `pathA_20`
outputs (`c`, the `c_γ/c_s` result, the `ħ`-verdict + `h`/`2π` assessment, the `flux_law_verdict`).

## Work items

### P1 — the inter-defect force from inflow (mass-as-attraction, the user's ontology)
- From the stationary Bernoulli/pressure field of a single drain (continuity + quantum-Bernoulli + EOS, the
  `pathA_20` profile), derive the field a SECOND defect feels, and the resulting force between two drains of fluxes
  `J₁, J₂` at separation `r`. Test whether it is attractive and its `r`-scaling (expect an inverse-power law from the
  drain pressure field). Dimension-check the force.
- This is the physical content of "mass is measured by the attractive force between two objects, directly from the
  inflow." The force law is the bridge between `J` and a gravitational coupling.

### P2 — derive `m_defect` from inflow (re-test the `pathA_19` negative result)
- Produce the relation `pathA_19` demanded and could not find: an action-level / boundary-source / Noether-charge /
  Hamiltonian-energy derivation that ties the throat rest (gravitational) mass `m_defect` to the inflow rate `J`.
  Target FORM (from `pathA_20`): **`m_defect = α_J ħ J / c²`** (= `E=mc²` with `ħ·(α_J J)` the rest energy = the
  standing-wave `ħω₀`). DERIVE `α_J`; do not assert it. **`2π`/`h` resolution:** if `pathA_20` S3 found `J` is a
  CYCLE-count rate (`ν`-like), the form is `m_defect = α_J h J_ν/c² = 2π α_J ħ J_ν/c²`; DERIVE whether the `2π` is
  physical or absorbed into `α_J` (this factor is calibration-critical — a stray `2π` shifts the GR-anchor match).
- **Equivalence-principle check (a derivation, if it lands):** show the INERTIAL mass (resistance to accelerating the
  throat — added mass / momentum of the entrained flow) equals the SAME inflow quantity as the gravitational/source
  mass. If both equal `α_J ħ J/c²`, the equivalence principle is DERIVED here, not assumed.
- A negative result is still valid: if no such relation exists, keep `m_defect` as a separate datum and report the
  precise obstruction (sharper than `pathA_19`'s, now that `c`/`ħ` are known).

### P3 — the M-collapse test (does the base set reduce to `{L,T}`?)
- With the derived `c` (`pathA_20`) and the `ħ`-verdict (`pathA_20` S3), test whether `M` is now ELIMINABLE:
  substitute the bridge + any `ħ`-emergence relation and check whether every mass reduces to combinations of
  `{L,T}`-quantities + already-derived constants, OR whether an irreducible dimensionful constant remains.
  - If `ħ` is emergent AND the bridge derives: `M` collapses → base `{L,T}`, mass-as-inflow PROVEN. Update a NEW
    harness representation (side-by-side; do not mutate the `pathA_19` `{L,T,M}` dictionary).
  - If `ħ` stays fundamental OR the bridge fails: retain `{L,T,M}`; record exactly which input blocks collapse.
- Either outcome is a valid PASS with a named verdict. This directly resolves the `pathA_19`
  `INFLOW_MASS_SOURCE_MISSING` residual one way or the other.

### P4 — emergent `G` and the `m↔G` unification
- Derive the effective Newton constant `G` from defect back-reaction on the condensate + the 4D→3D (bulk→brane)
  reduction (transverse width ~`a`). Honest prior from the verification work: effective-3D `[G]=L³T⁻²M⁻¹` (brane
  projection). Dimension-check in the `pathA_19` base.
- Show `m_defect` and `G` are **two faces of one inflow quantity** (the force in P1 `∝ G m₁ m₂/r²` with both `m` and
  the coupling sourced by `J`). State the closed-form `G(K, ρ₀, m, ħ, J, a, c_γ/c_s, …)` to the extent derivable;
  flag any piece that must wait for `pathA_22`'s scale-map.

## Acceptance criteria
1. P1 inter-defect force derived from the inflow field + dimension-checked; attractive/`r`-scaling reported.
2. P2 `m_defect`↔`J` relation derived (with `α_J`) OR a sharpened obstruction reported; equivalence-principle check
   resolved (derived / not).
3. P3 explicit M-collapse verdict: base `{L,T}` (mass-as-inflow proven) OR `{L,T,M}` retained with the named blocker.
4. P4 emergent `G` derived + dimension-checked; the `m↔G` (both = inflow) relation stated; `pathA_22` hand-offs flagged.
5. Dual-engine `.wl` agreement; new harness group passes; scripts exit 0 within `timeout 600`; pathA_18/19/20 groups
   untouched.

**Acceptance is PASS/FAIL with NAMED RESIDUALS: exit-0 is NECESSARY, NOT SUFFICIENT.** Every rejected hypothesis leaves
a named residual + source + downstream consequence. A retained `{L,T,M}` / unproven bridge is a VALID negative result.

## Out of scope
The scale-map → `m̂0²·S_port` → B2c rerun (`pathA_22`); any freeze amendment (methodology call with the user). Do not
touch `m̂0²·S_port` here.

## Review (orchestrator, after Codex)
Transliteration-fidelity (one clean agent per new script); independent re-derivation of the force law, the `m_defect↔J`
bridge, and `[G]`; adversarial pass (distrust-all-clean) targeting the M-collapse verdict and the equivalence-principle
claim (is it derived or smuggled?). Claude reads only residuals. Then gate to `pathA_22`.
