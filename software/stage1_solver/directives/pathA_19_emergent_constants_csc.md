# Directive pathA_19 — Emergent constants, Step 1 (FOUNDATION): re-found the dimensional system

**Date:** 2026-06-19 (rev. after Codex design-review `_scratch/pathA_19_directive_review.md`: SOUND-WITH-FIXES, 11
fixes applied — symbol split `m_GNLS`/`m_defect`, test-don't-declare mass emergence, PASS/FAIL residual gates, no
`m̂0²·S_port` here, dual-engine scoped to algebra).
**Owner:** Codex (DERIVES from the action + extends the dimensional harness + writes the foundation reference doc;
iterates until scripts exit 0). Claude reviews afterward (fidelity + dimensional + adversarial).
**Trigger:** decision-13 §4 (emergent-constants chunk), Step 1 of 4. Split out as its OWN gated step (user call,
2026-06-19): it may change the base dimensional system everything downstream sits on. Chain: **this (foundation)** →
`pathA_20` (`c_s`+`c`) → `pathA_21` (`G`) → `pathA_22` (scale-map → `m̂0²·S_port` → re-run B2c). Each gated.

## Why this step (context — do not re-litigate)
The B2c verdict is UNDETERMINED because the GR comparison pins the **dimensionful** factor `m̂0²·S_port=1` (decision-13
§1). We must derive the model's "constants" (PDE outputs, not fundamental inputs); the dimensional FOUNDATION comes
first. Three foundational questions (raised with the user, 2026-06-19) must be ANSWERED FROM THE ACTION as falsifiable
tests — **declaring the expected answer is NOT acceptance** (see PASS/FAIL gates below):
1. Is mass primitive, or is the **DEFECT** mass `m_defect` emergent from the defect inflow? (NOT the constituent
   GNLS mass `m_GNLS`, which is part of the exact action — see F1's symbol split.)
2. Is `a=1` a valid pin, or is `a` a particle-dependent collective moment, with the conserved invariant the flux `J`?
3. Are the natural-unit pins `a=c_s=ħ=m=1` over-determined (one hidden relation), and what is the EXACT relation?

Two facts the pre-work surfaced: (a) the committed harness `dimensional_check.py` already encodes a 4D-bulk base
(`ρ=L⁻⁴`, `K_eos=M·L¹⁸T⁻²`, `M` a base dimension) — the HARNESS is the authority; the research-paper PROSE is what uses
wrong 3-D conventions (e.g. `em_fields.tex` `ρ₀` as `kg·m⁻³`); (b) the exact parent GNLS action contains `m_GNLS`
explicitly in the kinetic operator, EOS sound speed, number current, Madelung velocity, and Euler/vorticity identities
(`research/pde_ledger/paper/parts/part01_parent_geometry.tex` ~174–203, 270–291; `research/pde/paper/pde.tex` ~326–352)
— that `m` is action content, not bookkeeping.

## Scope & stance
DERIVATION + harness extension. No model-formula changes, no GATE-A freeze touch, **no mention/derivation/change of
`m̂0²·S_port` (that is `pathA_22`)**. **Extend `dimensional_check.py` by ADDING a new check group** — do NOT mutate
`dimensional_dictionary()`, the global `D`, the D1–D3 checks, or the pathA_18 output behavior (committed standing
check). If you experiment with an `{L,T}` representation, do it SIDE-BY-SIDE (a separate representation/projection map),
never by changing the pathA_18 verdict. Infra: read-only on `research/pde_audit/simulation/`; never touch
`physical_export_permitted`; no `$RT exec-*` (MANIFEST race) — standalone `python3` + `math -script`; each script
`timeout 600` (exit 124 = reformulate); ≤2 concurrent `math -script` seats; YAML/markdown for human/LLM output.
**Dual-engine scope:** a `.wl` cross-check is REQUIRED for the ALGEBRAIC/DIMENSIONAL claims only (flux dimension,
pin-counting linear algebra, dictionary homogeneity, action-term homogeneity). It is NOT a verifier for prose
classification, boundary-condition choices, or whether the action contains a defect-mass source — those require
human-readable residuals.

## Work items

### F1 — The mass question: split the symbols, then TEST (do not declare)
- **Symbol split (mandatory).** Use `m_GNLS` for the constituent/inertial mass in the parent action's kinetic/EOS/
  current terms; use `m_defect` (or `M_defect`) for a throat's rest/gravitational mass. State plainly: the current
  action keeps `m_GNLS = [M]` UNLESS an explicit action-level reparameterization is derived that keeps every term
  homogeneous (see the homogeneity gate).
- **The hypothesis to TEST (not assume):** is `m_defect` emergent from the defect inflow? The existing ontology says
  defect MASS is associated with the drainage/VOLUME deficit (volume flux into the throat), and CHARGE with vorticity
  flux (`research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` ~297–302, 1967–1975) — so do NOT assume mass = the
  number flux. Frame-tag and distinguish the candidates (see F2): number flux `J=T⁻¹`, volumetric flux `Q_vol=ρ⁻¹J`,
  mass flux.
- **The honest relation (non-natural-unit form):** if `m_defect` is the rest-frequency/inflow rate, write
  `m_defect = α_J · ħ · J / c_γ²` with `α_J` dimensionless / branch-normalization data. Natural-unit equality
  `m_defect = J` is permitted ONLY after `ħ=1` and `c_γ=1` are established — which is `pathA_20`'s job, so using the
  de Broglie / standing-wave `c` argument to eliminate `M` HERE is circular. Downgrade the de Broglie/standing-wave use
  to a **conditional consistency note handed to `pathA_20`**, not a premise.
- **Action-homogeneity gate for `{L,T}`.** An `{L,T}` base set is permitted ONLY if `[ħ]`, `[K]`, the kinetic term, the
  gauge coupling, the Maxwell sector, the wall action, and `R_norm` ALL rewrite consistently (machine-checked). If they
  don't, `{L,T}` is merely a natural-unit *representation*, not a new base-dimensional system — say so.
- **Falsification rule (REQUIRED).** If no boundary/source/Noether-charge/Hamiltonian-energy derivation ties
  `m_defect` to the inflow, the foundation verdict MUST retain `{L,T,M}` and record `J` as a conserved RATE label only
  — a clean negative result, not a failure. Record the named residual that blocks emergence + its source + the
  downstream consequence.

### F2 — The flux `J` (frame-tagged) and the `a`-pin
- **Frame-tagged definition.** From continuity `∂_tρ+∇·(ρv_b)=0`, current `j=ρv_b`. Define and dimension-pin, distinctly:
  the 4D-bulk closed-3-surface number flux, the brane-projected 2-sphere number flux, the volumetric flux
  `Q_vol=ρ⁻¹J`, and the mass flux (× mass-per-particle/conversion). Confirm `[J_number]=T⁻¹` (robust: `L⁻⁴·L³` or
  projected `L⁻³·L²`).
- **Shape-independence WITH its limit.** Gauss gives surface-independent flux ONLY in a region with no enclosed
  sources/leakage. The exact projection yields an open-system leakage source (`part01_parent_geometry.tex` ~315–330;
  `pde.tex` ~512–539), and the throat bottom may be open/closed/connected (`brane_bulk_ontology.tex` ~1998–2039). So
  the "standing-wave recirculation / no net accretion" reconciliation must be DERIVED as a boundary condition OR
  recorded as a carried-forward gap — it may not be accepted as a note.
- **Re-examine the `a=1` pin.** `a` is a mouth-radius collective MOMENT, not a fundamental throat coordinate
  (`a0=R0(0)`, `a(t)` = mouth average of `R`; `part01_parent_geometry.tex` ~447–455, `pde.tex` ~633–648, 1076–1118),
  and it is particle-dependent + ill-defined under deformation. Assess whether the conserved flux `J` is the better
  invariant to hold fixed, with `a` as derived geometry. **Do NOT claim this changes `m̂0²·S_port`** (it may be
  dimensionally neutral); instead state only the dimensional dependencies and which scale-map inputs `pathA_22` must
  consume. Do not derive or change the normalization in this step.

### F3 — Base set, pins, over-determination, dictionary
- Resolve the pin over-determination: with `{L,T,M}` the four pins `a=c_s=ħ=m=1` over-fix 3 base units → one physical
  relation; with `{L,T}` the counting must be REDONE, not inherited.
- **The hidden relation is a DERIVATION TARGET, not an expected answer.** `a≈ħ/(m c_{s,0})` is the expected
  healing-length scale up to convention-dependent factors; derive the EXACT relation/factor from the GNLS core balance
  (`part01_parent_geometry.tex` ~174–203), and PERMIT a corrected factor or a different core-radius relation.
- Independent-vs-derived constants table, **conditional on the F1 symbol split**: `ħ, m_GNLS, K`, chosen state `ρ₀` are
  independent in the current action/harness; `c_{s,0}, a` derived; `m_defect` emergent only if F1's derivation +
  `ħ,c` conversion succeed.
- Confirm the harness dictionary against the action-derived dimensions; flag (do not edit pathA_18) any correction in
  the NEW group.

### F4 — Paper-prose reconciliation
Table of every place the papers state a quantity's dimension, marked AGREES / WRONG-3D-CONVENTION / AMBIGUOUS, with
file:line (prevents the 3D-attractor error from propagating into `pathA_20+`).

### F5 — Foundation reference doc
One markdown "dimensional foundation" doc under `software/stage1_solver/` (`reports/` or `notes/`): the base-set
verdict (F1, with the `m_GNLS`/`m_defect` split and the emergence verdict incl. any negative result), the frame-tagged
flux `J` + the `a`-pin assessment (F2), the pin over-determination + the DERIVED healing-length relation +
independent-vs-derived table (F3), the action-derived dictionary (frame-tagged 4D-bulk vs 3D-brane), the paper-prose
reconciliation (F4), and an explicit "gaps / carried-forward" section (EOS-from-GNLS; the no-accretion BC if unproven;
`m↔G` unification → `pathA_21`; what `pathA_20` must use).

## Acceptance criteria (PASS/FAIL with residuals — exiting 0 is necessary, NOT sufficient)
Every rejected hypothesis must leave a NAMED residual + source citation + downstream consequence. A negative result
(retain `{L,T,M}`) is a valid PASS if properly derived/recorded.
1. F1: symbols split (`m_GNLS` vs `m_defect`); the mass-emergence hypothesis either DERIVED (boundary/source/Noether/
   Hamiltonian, with the `m_defect=α_J ħ J/c_γ²` relation and the action-homogeneity check for any `{L,T}` claim) or
   explicitly REJECTED with the blocking residual recorded and `{L,T,M}` retained. No de Broglie/`c`-based elimination
   of `M` (circular) — that is a note to `pathA_20`.
2. F2: `J` defined frame-tagged (number/volumetric/mass) + dimension-pinned; Gauss shape-independence stated WITH the
   no-leakage condition; the no-net-accretion claim derived as a BC or logged as a gap; `a`-pin assessed (J as
   invariant, `a` as moment) WITHOUT any `m̂0²·S_port` claim — only dimensional dependencies + `pathA_22` inputs.
3. F3: over-determination resolved; healing-length relation DERIVED (exact factor, correction allowed) not assumed;
   independent-vs-derived table conditional on F1; harness dictionary confirmed or corrections flagged in the NEW group.
4. F4: paper-prose reconciliation table (AGREES/WRONG-3D/AMBIGUOUS + file:line).
5. Dual-engine: `.wl` reproduces the ALGEBRAIC/dimensional claims (flux dim, pin-counting, dictionary + action-term
   homogeneity); non-algebraic physics judgments carry human-readable residuals, not `.wl` "verification".
6. New harness group added + passing; `dimensional_dictionary()`, D1–D3, pathA_18 output all UNCHANGED; scripts exit 0
   within `timeout 600`.
7. F5 foundation reference doc committed-ready, frame-tagged, with the honest gaps/carried-forward section.
8. NO change to any model formula, the GATE-A freeze, or `m̂0²·S_port`.

## Out of scope (later gated sub-steps)
`c_s`+`c` (`pathA_20`, draft staged); emergent `G` + 4D→3D reduction + the `m↔G` unification (`pathA_21`); the
scale-map → `m̂0²·S_port` (`pathA_22`); B2c rerun; any GATE-A freeze amendment.

## Review (orchestrator, after Codex)
One clean agent runs the **transliteration-fidelity audit** on every new derivation/harness script; a
**dimensional-consistency** review independently re-derives `[J]`, the pin over-determination, and the action-homogeneity
of any `{L,T}` claim; an **adversarial** pass (distrust-all-clean) targets specifically (a) whether F1 actually
DERIVED `m_defect`-emergence or just asserted it (and whether a negative result was honestly taken), and (b) the
`a`-pin assessment. Claude reads only residuals. Only after this passes do we gate to `pathA_20`.
