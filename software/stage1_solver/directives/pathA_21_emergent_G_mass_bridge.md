# Directive pathA_21 — Emergent constants, Step 3: `G`, the mass-bridge, the `m↔G` unification (SYMBOLIC forms + profile-solve spec)

**Status:** Design-reviewed (Codex `gpt-5.5` → SOUND-WITH-FIXES; all fixes applied 2026-06-19; pending confirm-pass before
execution — gated by the user). This is the SYMBOLIC spec-completion step: pathA_19/20/20b extracted everything derivable
from the symbolic action + dimensions; the remaining NUMBERS all bottleneck on the SOLVED THROAT PROFILE. So pathA_21
DERIVES the FORMS (`G`, `m_defect`, the inter-defect force as profile-integral expressions) carrying
`λ_γ`/`α_J`/profile-integrals SYMBOLIC, re-tests the `M`-collapse, and PRODUCES THE PROFILE-SOLVE SPECIFICATION (the exact
list of profile quantities the future solve must compute). Concrete consumed inputs:
- base **`{L,T,M}`**, `m_GNLS` (constituent) ≠ `m_defect` (throat) — pathA_19;
- `[J]=T⁻¹` is the conserved invariant LABEL, but its VALUE = `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA` — pathA_20;
- `ħ`-provenance = `HBAR_PROVENANCE_UNDETERMINED` (NOT shown emergent) + `H_2PI_RATE_CLASSIFICATION_UNDETERMINED` (is `J`
  cycle-rate `J_ν` or angular-rate `J_ω`, hence where the `2π` sits) — pathA_20;
- `c=c_γ` is the standing-wave ceiling. In the strict bulk transverse zero-mode, `c_γ` is represented by `c_bulk` with
  `c_bulk²=C_B/C_E`; the OBSERVED brane `c_γ/c_s` remains the `brane_verdict` residual and is carried as `λ_γ` until the
  zero-mode/profile reduction closes it (`c_γ/c_s` is a CALIBRATION KNOB) — pathA_20b;
- healing: `a=ħ/(m_GNLS c_s0)`, `ξ_h=√2 ħ/(m_GNLS c_s0)`, `h0=(5K/4)ρ0⁴=m_GNLS c_s0²/4`; `c_s²=5Kρ0⁴/m_GNLS` — pathA_19.
**Date:** 2026-06-19
**Owner:** Codex (DERIVES + codes; iterates until scripts exit 0). Claude reviews afterward.
**Trigger:** user decision A (decision-13 §0/§10, 2026-06-19). Chain: `pathA_19` → `pathA_20` → `pathA_20b` →
**this (`G` + mass-bridge + spec)** → `pathA_22` (scale-map → `m̂0²·S_port` → re-run B2c).

## Why this step
(1) Re-test the `pathA_19` mass-as-inflow negative now that `c` exists and the `ħ` verdict is in. (2) The user's ontology:
a defect's mass = the attractive pull via condensate inflow → the inter-defect force, `m_defect`, and `G` are three faces
of the inflow `J`. (3) Everything bottlenecks on the throat profile, so the highest-value output is a CONCRETE
specification of what the profile solve must compute.

## Scope & stance
SYMBOLIC DERIVATION (forms, not numbers). Carry `λ_γ`, `α_J`, `J`-value, and profile integrals SYMBOLIC; do NOT
manufacture missing profile/operator data (inherited as blocking residuals). No model-formula changes, no freeze touch,
no `m̂0²·S_port` un-pin (`pathA_22`). Mass-symbol discipline: `m_GNLS` (EOS/healing/Madelung/action) vs `m_defect`
(throat rest/gravitational). Extend `dimensional_check.py` side-by-side (new group; leave pathA_18/19/20/20b groups +
`dimensional_dictionary()`/D1-D3 intact). Same infra constraints (read-only sim dir; never touch
`physical_export_permitted`; no `$RT exec`; `timeout 600`; ≤2 MMA seats; YAML/md). Do NOT collapse `M` unless an
INDEPENDENT `ħ`-free relation actually appears (pathA_20b left `ħ` UNDETERMINED → expected honest outcome is `{L,T,M}`).

**DERIVED-FORM GATE (anti-restatement — the core acceptance teeth).** A positive P1/P2/P4 form is accepted ONLY when the
report gives a SOURCE-EQUATION CHAIN from the parent equations to a named profile integral/relation. It is a FAIL to
introduce the target formula as an assumption, to define `α_J := m_defect c²/(ħ J)` (rearrangement), or to define `G`
only by rearranging `F = G m₁ m₂/r²`. A missing chain → a NAMED RESIDUAL (a VALID outcome), never a PASS. Dual-engine
`.wl`/SymPy agreement validates DIMENSIONS/ALGEBRA ONLY; it cannot upgrade a restatement into a derivation — every
non-algebraic P1/P2/P4 derivation needs a human-readable proof trace with source equations + residuals.

## Work items

### P1 — the inter-defect force from inflow (FORM of the force coefficient `C_F`, WITHOUT `G`)
- From the stationary drain pressure/Bernoulli field of a single defect, derive the FORM of the field a second defect
  feels and the force between two drains of fluxes `J₁,J₂` at separation `r`. Report the force coefficient
  `C_F(profile; J₁,J₂,…)` and the `r`-power — do NOT introduce `G`. Determine attractiveness FROM the pressure/density
  response (a sign from convention or from choosing `α_J>0` is NOT acceptable — else carry the sign as a residual).
- Keep the FULL stationary-profile dependencies (`Q`, `V_conf`, geometry, leakage/topology BC). Any ideal-Euler/nozzle
  reduction (which pathA_20 refused to promote) must be labelled CONDITIONAL and NON-CLOSING. Dimension-check the force;
  FLAG which profile integrals are unresolved (→ P5).

### P2 — `m_defect` from inflow (re-test the pathA_19 negative; FORM + `α_J`)
- **Candidate FORM to test:** `m_defect = α_J ħ J / c²` (= `E=mc²` with the inflow as the energy = standing-wave `ħω₀`).
  DERIVE OR REJECT it. A positive result REQUIRES an independent action-level / boundary-source / Noether / Hamiltonian
  relation defining `α_J` as a profile functional that does NOT contain `m_defect`, `G`, or the target bridge. Restatement
  or `α_J` by rearrangement → `MASS_BRIDGE_FORM_NOT_DERIVED` (a valid residual, NOT a PASS).
- **`2π`/`h`:** since `J` cycle-vs-angular is UNDETERMINED, use SEPARATE symbols `J_ν` (cycle) and `J_ω` (angular), give
  both forms (`α_J ħ J_ω` and `α_J h J_ν = 2π α_J ħ J_ν`), and show exactly where the `2π` sits + which residual blocks
  choosing it. `α_J` must NOT silently absorb the `2π`.
- **Equivalence-principle check (no smuggling):** compute `m_inertial` (accelerated-throat kinetic response /
  entrained-flow momentum) and `m_source` (far-field drain/force coefficient from P1) SEPARATELY, from separate sources.
  EP is DERIVED only if both reduce to the SAME profile integral with the SAME normalization. Otherwise `EP_NOT_DERIVED`
  (valid outcome).

### P3 — the `M`-collapse re-test (does the base reduce to `{L,T}`?)
- With `c` derived and `ħ` = UNDETERMINED, test whether `M` is eliminable. The anti-tautology gate binds: `M` collapses
  ONLY with a genuine INDEPENDENT `ħ`-free relation (absent in pathA_20b). **Honest expected outcome: retain `{L,T,M}`.**
  If no independent mass bridge is derived (P2), the correct resolution is to RETAIN or RENAME `INFLOW_MASS_SOURCE_MISSING`
  with a sharper blocker — do NOT force closure; `{L,T,M}` retained is a PASS. Only with a real `ħ`-free relation, collapse
  to `{L,T}` in a NEW side-by-side representation (never mutate the pathA_19 `{L,T,M}` dictionary).

### P4 — emergent `G` (FORM, derived AFTER P1/P2) + the `m↔G` unification
- Do NOT define `G` by rearranging the Newton form. Extract a universal `G_eff` ONLY after: (a) P1 gives `C_F` without
  `G`, and (b) P2 independently defines the masses — AND only if the force is attractive, INVERSE-SQUARE, factorizes into
  the same source quantities, and is universal for the branch class (not pair/instance-specific except via allowed
  universal branch data). Otherwise report `NEWTON_G_FORM_NOT_DERIVED` or `FORCE_NOT_NEWTONIAN` (record the observed power
  law + consequence).
- If extracted: express `G` as `G(K, ρ0, m_GNLS, ħ, λ_γ, α_J, W_eff; {profile integrals})` and dimension-check
  (`[G]=L³T⁻²M⁻¹` effective-3D). The 4D→3D reduction width is a NAMED profile/reduction quantity `W_eff` (a
  source-anchored branch relation) — do NOT set it `=a` or `=ξ_h/√2` without a source relation (`a` is a branch moment,
  not an invariant width). Show `m_defect` and `G` are two faces of one inflow quantity (FORM); flag `pathA_22` hand-offs.

### P5 — the PROFILE-SOLVE SPECIFICATION (the deliverable that makes this step worth doing)
Consolidate, from P1–P4, the EXACT machine-checkable table of every quantity the future stationary-throat-profile solve
must compute. **Required columns per row (no unnamed placeholders like "the profile integral for α_J"):**
`symbol` ; `definition` (domain + measure + integrand + profile fields + BC/source equation) ; `dimension` ;
`frame` (4D-bulk / brane / reduced-3D) ; `source anchor` (file:line or equation label) ; `closes which output`
(`C_F`/`α_J`/`G`/`brane c_γ`/`J`-value/…) ; `status` (`known` / `profile-solve` / `pathA_22` / `new-physics`) ;
`residual if absent` ; `downstream consumer`. Include the branch-data-set `𝔅` items from `pde.tex:2515-2566` where
relevant. DISTINGUISH quantities needed by pathA_21's forms from those deferred to pathA_22's scale map (do not drift
into scale-map work).

## Acceptance criteria (PASS/FAIL with NAMED RESIDUALS; exit-0 NECESSARY not SUFFICIENT)
1. P1 force-coefficient `C_F` FORM derived from the inflow field WITHOUT `G`, with a source-equation chain;
   attractiveness from the pressure/density response (or residual); `r`-power reported; full-profile deps kept (ideal
   reductions labelled conditional); profile integrals FLAGGED.
2. P2 mass-bridge DERIVED OR REJECTED per the derived-form gate (`α_J` an independent profile functional, or
   `MASS_BRIDGE_FORM_NOT_DERIVED`); `J_ν`/`J_ω` + `2π` placement explicit; EP resolved (derived with matching integrals,
   or `EP_NOT_DERIVED`).
3. P3 explicit `M`-collapse verdict: `{L,T,M}` retained (named blocker, residual retained/renamed) OR `{L,T}` (only with
   a real `ħ`-free relation, new side-by-side representation).
4. P4 `G` derived ONLY via P1+P2 (not by rearranging Newton) with the inverse-square/attractive/universal conditions, OR
   `NEWTON_G_FORM_NOT_DERIVED`/`FORCE_NOT_NEWTONIAN`; if derived, `[G]=L³T⁻²M⁻¹` checked + `W_eff` source-anchored;
   `m↔G` (both = inflow) FORM stated; `pathA_22` hand-offs flagged.
5. P5 the profile-solve specification: the full-schema table (above), specific enough to drive the option-C solve, with
   the `𝔅` items and the pathA_21-vs-pathA_22 split.
6. Dual-engine `.wl` agreement for ALGEBRAIC/DIMENSIONAL claims only; non-algebraic derivations carry a human-readable
   proof trace (`.wl` agreement is NOT proof of a derivation); new harness group passes; scripts exit 0 within
   `timeout 600`; pathA_18/19/20/20b groups untouched.

**Fail conditions (explicit):** restating a target form as derived; `α_J := m_defect c²/(ħ J)`; `G` by rearranging
`F=Gm₁m₂/r²`; EP by naming the same formula twice; a P5 row with an unnamed placeholder. A retained `{L,T,M}`, a
`MASS_BRIDGE_FORM_NOT_DERIVED`, an `EP_NOT_DERIVED`, or a profile-dependent (symbolic) `α_J`/`G` are VALID expected
outcomes — the win is the rigorously-derived FORMS + the P5 spec, not numbers.

## Out of scope
The scale-map → `m̂0²·S_port` → B2c rerun (`pathA_22`); the throat-profile SOLVE itself (later option C — pathA_21 only
SPECIFIES it); any freeze amendment. Do not touch `m̂0²·S_port`.

## Review (orchestrator, after Codex)
Transliteration-fidelity (one clean agent per new script); independent re-derivation of `C_F`, the `m_defect↔J` bridge +
`α_J`'s defining relation, and `[G]`; adversarial pass (distrust-all-clean AND distrust-restated-target — for each of
P1/P2/P4, is there a real source-equation chain, or a restatement? is `G` extracted only after a `G`-free force law? is
EP derived from two separate masses?); plus a check that the P5 spec is concrete + complete (would it drive the profile
solve?). Claude reads only residuals. Then gate to `pathA_22`.
