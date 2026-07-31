# DEFECT REGISTER

**One row per known defect. A list, not machinery.** ⛔ Do not grow tooling around this file — the moment
it needs a schema, a checker, or a manifest, it has become the thing that killed the census.

**Status vocabulary:** `OPEN` · `RETIRED` (fixed, with the commit) · `FALSIFIED` (a closed negative
result — ⭐ keep, these are the expensive ones) · `QUARANTINED` (known-wrong, deliberately left in v2
reference material) · `HYGIENE` (citation/provenance, not physics).

Every row carries a locus. ⛔ If you cannot cite one, it is not a defect yet — it is a suspicion.

---

## A. Pin-shaped identifications

⭐⭐ **The governing pattern of this project: two quantities that share a dimension get silently
equated.** **Thirteen** instances found — and ⛔ **A11 was created by me while writing this very
register**, which is the strongest evidence available that this is a systematic pull rather than a set
of past mistakes. This is not a bug that recurs — it is *the* systematic failure mode, and
it is why `[[feedback-matching-number-is-not-evidence]]` exists. A new one should be assumed present in
any sector that has not been checked.

| id | the identification | evidence | status |
|---|---|---|---|
| **A1** | `a = ħ/(m c_s0)` — the units pin read as a physical throat radius | `stage004:122-137` (rewritten) | ✅ **RETIRED** `407eed94` |
| **A2** | `ℓ = √(2κ_χ/λ_χ) = √(κ_χ/2a_B) = δ` — closure-card length identified with the ledger's kink width | `ledger_stage044_parent_action_reconciliation.md:104` | **OPEN** |
| **A3** | `L_W := L` — auxiliary mixed-tube length identified with the physical throat length; the spec itself calls it a *"Frozen identification"* and *"not independently confirmed"* | `two_throat_simulation_handoff_spec.md:464` | **OPEN** |
| **A4** | implicit `β·a = 1` — `L0/a = 37/20` and `β·L0 = 37/20` together equate the wall-response inverse length with the inverse mouth radius, while the register says geometry does **not** derive `β` | `parameter_register.md:101-110,:178` vs `notes/stage013_pathA31_breathing_source_map.md:124-133` | **OPEN** |
| **A5** | `c_w` ↔ `c_s` — two artifacts treat support waves in the **same cavity**, and nothing relates their wave speeds. ⭐ **Walked back 2026-07-31 — the sources support less than the previous sharpening asserted:** what the two sides demonstrably share is a **geometric cavity**, the throat interior — a cylinder of mouth radius `a` and `w`-extent `L` (*"a cylindrical cavity of radius $a$ and length $L$"*; *"the region $0 < w < L$ … is approximately a straight 4D cylinder"* whose interior *"is the domain in which the cavity modes of the electromagnetic sector live"*) — ⛔ **not** an identity of the modes living in it. `pathA_26`'s wave is a **declared scalar surrogate** (*"Minimal wave field: scalar surrogate in 4D"*), and its own Scope section limits it to a stand-in: *"The scalar wave stands in only for the two-derivative scaling of the trapped support sector."* The lepton notes instead write a **medium-acoustic** support mode, `ω_supp = c_s χ_w/a`; and the physical support sector is **brane shear**. ⚠ The engine's `c_w = 1.0` is a **natural-unit numerical slice** — *"these values only select the able-to-fail numerical slice"* — ⛔ **not** evidence that `c_w = c_s`. ⇒ **No bridge** is established among the scalar surrogate, the acoustic mode, and the brane-shear support sector — and that missing bridge is exactly what makes this an open cross-artifact identification question | `pathA_26_derrick.md:21`, `:95-97` (Scope) vs `notes/lepton_mass_notes.md:154-167`; cavity = throat interior: `research/em_fields/paper/em_fields.tex:333`, `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex:643-646`; surrogate: `notes/inner_throat/inner_throat_hard_mode.md:459-475`; natural-unit slice: `software/stage1_solver/tools/pathA_26_derrick_sympy.py:37-39,:45` | **OPEN** |
| **A6** | the length chain `a = 20ℓ`, `ℓ = δ` — collapses every particle-sector length onto one scale | `parameter_register.md:200` + A2 | **OPEN** |
| **A7** | `ξ = ħ/(m c_s)` (stage011) vs `ξ_h = √2 ħ/(m c_s0)` (stage004) — two "healing lengths" a factor `√2` apart | `parameter_register.md:294` (R27), `:94`; `ledger_stage011_..._sympy_audit.py:430,:521,:658` | ✅ **NOT A DEFECT** — investigated 2026-07-31, see below |
| **A8** | `r_e` called *"the throat-**body** size"* while the balance `ke²/a = m_e c²` is written with the **mouth** radius `a` | `conceptual_foundation.md:449-455` | **OPEN** |
| **A9** | `W_slab` (brane slab width) merged with the `L/a` (sleeve aspect ratio) self-selection debt, without an equation | `two_throat_simulation_handoff_spec.md:160`; `parameter_register.md:163,:286` | **OPEN** |
| **A10** | sonic horizon `r ~ (J/ρc_s)^{1/3}` identified with the mouth radius | *(orchestrator error, 2026-07-31, retracted same session)* | ✅ **RETIRED** — listed so it does not return |
| **A11** | ⭐⭐ **CREATED BY ME, 2026-07-31, while writing the document that warns about this.** S16/S22 treat the `a` in the worldtube expansion `O(a²/r²)` as **the throat mouth radius**. The 1PN source defines `a` only as a **profile width / characteristic support scale** and supplies **no bridge** to the mouth radius; the bridge report explicitly warns that model `a` *"is not an invariant reduction width"* | `4d_1pn_full.tex:734-765`; `pathA_21_emergent_G_mass_bridge.md:87` | **OPEN** |
| **A12** | **A third healing-length convention.** The register claimed exactly two; the lepton notes use a third, `ξ_h = ħ/(√2 m c_s)` — differing by `√2` and `2` from the other two | `notes/lepton_mass_notes.md:1310-1317,:1501-1508` | **OPEN** |
| **A13** | **`χ_B` ontology contradiction.** The plan (and the sim spec) make `χ_B = r_B e^{iθ_B}` **complex**; the staged ontology defines it as a **real** scalar in `[0,1]` with ordered density `n_B = χ_B n`. The complex form makes that density complex and silently adds a phase DOF. ⛔ **Must be resolved BEFORE S1/S5/S12, not merely listed** — S1 and S12 use `χ_B` as a real material fraction so `χ_B·n` is real, while S5 writes it complex. ⭐ Clean option: use `r_B ∈ [0,1]` in the ordered-density and conversion balances, reserving `χ_B = r_B e^{iθ_B}` for the complex field. ⚠ Then pick **one** branch and count it: the inertial `{Z_χ, κ_χ, λ_χ}` action, **or** the staged real field with its dissipative mobility. ⛔ `{a_B, κ_B}` cannot parameterize an action containing `Z_χ` — they are not equivalent | `ledger_stage006_two_phase_chiB_ontology.md:60-100` vs `two_throat_simulation_handoff_spec.md:46` | **OPEN** |

### ⭐⭐ A7 investigated — it is not a defect, and it explains why A1 survived so long

**Verdict: a legitimate second convention, correctly firewalled.** Both lengths are real and each is
natural in its own context:

- `ξ_h = √2 ħ/(m c_s0)` — from **core balance**, `√(ħ²/(2 m h0))` with `h0 = m c_s0²/4` (`stage004`).
- `ξ = ħ/(m c_s)` — from the **BdG dispersion**. The deferred `k⁴` term's ratio is
  `ħ²k²/(4m²c_S²) = (kξ/2)²`, and that identity *forces* `ξ = ħ/(m c_S)` with no `√2`. It is
  load-bearing — it defines stage011's validity window `kξ ≪ 1` — not decorative.

`parameter_register.md:294` already flags it (*"source's no-√2 convention"*) and R27 is a **firewall**
doing its job (`ξ ≠ ℓ_c`, never substituted). ⇒ No action. Row kept so the question is not re-opened.

⭐⭐ **But this is the origin story for A1, and it is worth more than the row.** The pin's value
`ħ/(m c_s0)` was **never arbitrary — it *is* a healing length, in a standard convention.** That is
precisely why identifying `a` with "the GNLS healing core" felt natural and survived eleven months.

⇒ **The pin's defect was never that the number was wrong. It was a category error:** a **medium**
length used as a **defect** radius. ⛔ A future check that asks *"is this value right?"* would have
passed it. The question that catches it is ***"is this quantity the same KIND of thing?"*** — one number
for the whole medium, or one per particle.

⚠ **Apply this test to every open row A2–A13.** Each is two lengths that agree numerically or dimensionally. The
question is not whether the numbers match; it is whether both sides are indexed by the same thing.

### ⭐⭐ **`A-CAND`** — a candidate pin, flagged BEFORE it is made (2026-07-31)

⭐ **The register's first pin written down in advance rather than caught afterwards.** Every row above was
found *after* the identification had already been made and had propagated. ⚠ This one is ⛔ **NOT a
finding** and ⛔ **NOT yet an identification** — it is a cross-sector pattern that *looks* like one shape,
recorded before it is made because making it prematurely is precisely this section's governing failure
mode. ⛔ It is deliberately not a numbered instance: nothing may cite it as one.

- **EM side:** the far-field **FORM** is target-blind EARNED — *"The `1/R²` falloff and the `s₁s₂` product
  of the far-field shell are target-blind EARNED"* — while the **SIGN** is `R1_REQUIRED(bc_selection)`,
  *"the sign is NEITHER earned NOR calibrated"*, `outcome_not_invariant` across the four BC classes
  `{V,M,J,MIXED}`, and the **MAGNITUDE** `Q_E` is `R1_REQUIRED(magnitude)`.
- **Gravity side:** the `1/r²` law **and** the attractive sign are ⛔ CONDITIONAL pending the S14a drain
  bridge.
- **⛔ The shared shape is a QUESTION to be tested, not a stated shape:** *do the two sectors share one
  shape, or only one sentence?* ⛔ It may **not** be written as *"far-field FORM earned; SIGN and MAGNITUDE
  routed through the throat interior"* — that is **false on the gravity side**, whose magnitude also carries
  an **external/calibrated** debt: `CHARTER.md#clean-close-not-validation` says a clean v3 close shows at most that the *"far
  field is computable given supplied mass, multipoles, compactness and a calibrated normalization"*.
  ⇒ ⛔ Keep the two sectors' classes **as written above and separate** until the test is run: EM's SIGN and
  MAGNITUDE are `R1_REQUIRED` on the interior solve; gravity's `1/r²` and sign are CONDITIONAL on the S14a
  bridge **and** its normalization is calibrated. ⛔ **And there is no presently earned common form.** EM
  earns a candidate `1/R²` form only ***within Q1's postulated G0 closure*** — the plan's own words, ⛔
  *"not from primitives"* (`V3_STEP_PLAN.md#q2-earned-within-g0`) — while gravity's analogous `1/r²` form is the very
  thing called CONDITIONAL two bullets up. ⇒ **The shared shape is itself unearned**, and that is precisely
  what the knit must test.

⛔ **What may make it wrong — and this carries at least the weight of the pattern itself:** gravity's
`1/r²` is conditional because a **derivation chain was severed** (*"Correcting the drain to the dynamical
`Γ_B` law severed the chain to these"*), whereas EM's sign is open because a **boundary-condition class is
unselected**. A severed chain and an unselected BC class may be different **KINDS** of open, and welding
them into one shape would be a pin of exactly the type this section catalogues.

⭐ **The test to apply at the knit is A7's — *"is this the same KIND of thing?"* (`:58`) — ⛔ not *"do these
look alike?"***

**Evidence:** `notes/stages/ledger_stage031_puncture_deflection_field_identity_source.md:25`;
`parameter_register.md:329`, `:143`; `pde_ledger_v3/CHARTER.md#conditional-s14a`.

**Status: OPEN — CANDIDATE, to be tested at the knit (PHASE 5).**

---

## B. Falsified — closed negative results

⛔ **These are the most expensive things in the corpus and the easiest to lose in a restart.** They do
not look like progress. Import them first.

| id | result | evidence | status |
|---|---|---|---|
| **B1** | **The lepton mass tower is falsified.** The support-mode family gives `F_j ∝ (2j+1)²` = `1:9:25`; observed is `206.77`, `3477.37`. *"decisively ruled out"* | `notes/lepton_mass_notes.md:424-473`; `conceptual_foundation.md:487-488` | **FALSIFIED** |
| **B2** | **Gate L returned a no-go.** `FAIL_COUPLE_STRESS_NOGO` — the gate the medium survey had flagged as *"Highest risk; the most likely no-go"*. The GNLS polar-smectic program is *"SUPERSEDED at the brane-existence level"* | `notes/stage030_pathA35_gateL_source_map.md:16`; `conceptual_history.md:340`; `medium_requirements_and_prior_art.md:172-183` | **FALSIFIED** | {#B2}

⭐ **B1's consequence is the largest open question in the project.** The falsified picture — leptons as
support-mode excitations `j` of *one* throat — was what made charge universality automatic. Killing it
leaves the model with **no family label**: nothing explains how mass differs by 207× while charge is
exactly identical. ⇒ See **D1**.

---

## C. Open negatives and gaps

| id | the gap | evidence | status |
|---|---|---|---|
| **C1** | **`THROAT_DRAIN_DESTABILIZED`** — the reduced throat is unstable under its own drain. ⚠ **Read `:101-123` before citing the top line**: it says explicitly *"must NOT be read as 'the gravitating particle is killed'"*. Conservative existence is generic (75/75 scans); instability is above `gcrit≈0.006` while the tested box demanded `g≈89`. **The real gap: `g_phys` was never mapped.** | `software/stage1_solver/reports/pathA_26_derrick.md:58-72`, `:101-123` | **OPEN** |
| **C2** | **`L/a` contradiction.** Solved `L*/a* = 0.9575`; frozen `37/20 = 1.85` is what everything downstream uses. ⭐ Unsupported from four directions: admitted ansatz, adjudicated `free_choice`, contradicted by the solve, and in tension with **C5** | `pathA_26_derrick.md:43` vs `em_u1_body_definition.md:16`, `spec:73` | **OPEN** |
| **C3** | **`m_defect` is a GAP** — `INFLOW_MASS_SOURCE_MISSING` / `BLOCKS_MASS_EMERGENCE`. Only a dimensional bridge `α_J ħJ/c_γ²` exists; pathA_21's verdict is `MASS_BRIDGE_FORM_NOT_DERIVED`, *"No action-level, boundary-source, Noether, or Hamiltonian equation"* | `stage004:181`; `parameter_register.md:152-153`; `pathA_21_emergent_G_mass_bridge.md:98` | **OPEN** |
| **C4** | **The geon has no equation.** *"its profile is a **declared OPEN input** [POSTULATE: mass mechanism; profile OPEN]"* — the thing that carries the mass | `two_throat_simulation_handoff_spec.md:219` | **OPEN** |
| **C5** | **Charge universality is required but not earned.** The conceptual story says charge is pure direction with universal magnitude; the computed Coulomb magnitude is geometric (`3Q_E²ℓtanh²(b/ℓ)/(8πRb²)`) with `Q_E` `[CALIBRATED]`. Manifest: `q_h/Q_E universality is required but not earned from the current b/ell ledger` → `FALSIFIABLE_TENSION` | `conceptual_foundation.md:440-450,:590`; `spec:258,:260,:478` | **OPEN** |
| **C6** | **No closed parent action.** *"Closure status — no closed parent action exists."* The BVP is a *"template, not yet a closed BVP"*; well-posedness `[OPEN]`; status `UNRESOLVED(closure)` | `spec:324`, `:429`, `:439`, `:441` | **OPEN** |
| **C7** | **`W_slab` is `FREE-UNREDUCED`** — *"double-well selects NO width"*. The kink gives one interface; the brane is a finite slab. R19 `PENDING` | `parameter_register.md:163,:286` | **OPEN** |
| **C8** | **The wall rests on two postulated constants.** `δ` and `σ_wall` are genuinely derived — but from `{a_B, κ_B}`, which are postulated, not medium primitives | `ledger_stage006_two_phase_chiB_ontology.md:136-141`; `parameter_register.md:156,:287` | **OPEN** |
| **C10** | ⭐⭐ **The Spin Problem — a structural tension INSIDE the gravity sector, and it is in v3's scope.** The **inertial** sector (correct 1PN precession) constrains the throat to a **compact, "stubby"** geometry with a small aspect ratio. The **spin** sector needs the gravitomagnetic potential to fall off as a **dipole**, and a compact source gives the wrong scaling — a simple vortex yields a *"gravitomagnetic **monopole**"*, called *"physically inadmissible"*. Verdict in the note: *"You cannot get frame dragging from a compact 4D bubble; you **need** the tail."* ⇒ **Inertia wants compact; spin wants extended.** ⚠ The proposed fix — a composite *"Ion-Vortex Complex"* (stubby head + infinite vortex-filament tail) — is **uncited, unaudited, and referenced by nothing else in the corpus**; it also re-introduces **quantized circulation**, which the charge sector explicitly disclaims (*"not an additive winding"*). ⛔ Meanwhile `conceptual_foundation.md:589` still lists spin as *"**not yet placed in the picture**"* — so the corpus holds "unsolved" and "solved" simultaneously | `research/4d_1pn_bridge/notes/tadpole.md:1-19,:150-174`; `docs/conceptual_foundation.md:589` | **OPEN** |
| **C9** | **Live contradiction between tracked documents** on `K0c`/`K_eta`/`T_Omega`: `FREE-UNREDUCED`, PENDING, *"do NOT assert DERIVED"* vs *"likely DERIVED manifestations"*. A tier-1-vs-tier-3 disagreement about the same three symbols | `parameter_register.md:170` vs `notes/stage023_pathA34_nullspace_underdetermination_source_map.md:250-253` | **OPEN** |

---

## D. The question with no row of its own

| id | | status |
|---|---|---|
| **D1** | ⭐⭐ **What makes a muon a muon?** Mass differs 207× while charge is *exactly* identical. **B1** killed the *support-only* family label. ⛔ **Corrected 2026-07-31 — the stronger claim was wrong:** it is **not** true that no slope survives. `lepton_mass_notes.md:863` says outright that *"the old support-only `1:9:25` falsifier was **too simple**, because once throughput and geometry are allowed to respond dynamically, the family ladder changes"*, and §6/§13 develop later routes (`φ_j = R_j^{3/2}/√ν_j`, the low-harmonic benchmark `Wν = R^{9/5}`). ⇒ The honest statement: **no unique, target-blind, empirically successful mass–radius law survives; several mutually incompatible conditional slopes remain** — and `:857` notes the turbine route bridges the gap only *"in a regime very close to the critical point `s=3/4`"* and *"does not naturally produce"* the ratios | `notes/lepton_mass_notes.md:736,:857,:863,:2195` | **OPEN** |

⚠ An earlier `a ∝ m^{1/3}` claim (orchestrator, 2026-07-31) was **wrong**: it dropped `B` from
`F = A/a + B/a² + Ca³` while keeping the `18/11` factor that requires `B ≠ 0`. Three independent
reviewers, SymPy-verified. ⇒ Recorded so it is not revived.

---

## E. Quarantined — known-wrong, deliberately left in v2

⭐ The scope rule these sit under: **fix what computes, quarantine what only narrates.** ⛔ Their presence
here is not permission to cite them.

| id | | evidence | status |
|---|---|---|---|
| **E1** | 13 loci tag a **physical** radius `CONV` so it is not counted — several under headline claims *"ZERO new counted knobs"*, *"Part-II counted CALIB set UNCHANGED"*. ⇒ **v2's irreducible count is understated** | `parameter_register.md:186,:200,:306,:308,:310,:312,:314,:629,:669,:710,:731,:745` | **QUARANTINED** |
| **E2** | `stage004`/`stage005` engines still **assert** the retired pin relation via `expect_zero` | `..._stage004_sympy_audit.py:302,:306` | **QUARANTINED** — the stage note now records them as stale |
| **E3** | `midway_knob_audit_codimension_sympy.py` still hardcodes **ten** scalars including `a`; the registry now carries **nine** ⇒ they disagree | `scripts/midway_knob_audit_codimension_sympy.py:91-97` | **OPEN** (created by `407eed94`) |
| **E4** | pin conflations in the stage notes | `stage016:184`; `stage023:407,:481` | **QUARANTINED** |
| **E5** | `acceptance_check.py`'s header still says the payload was *"copied verbatim"*; it was recomputed | `reduction/acceptance_check.py:11-12` | **OPEN** (cosmetic) |

---

## F. Corpus hygiene

| id | | evidence | status |
|---|---|---|---|
| **F1** | ⭐ **Zero cross-artifact citations resolve to a locus** — measured on both pilot stages, both engines | — | **HYGIENE** |
| **F2** | False provenance: stage016's engines assert `M̃`/`K̃`/`T̃_Ω` are `CONSUMED-from-011/012/013`; those stages contain **none** of those symbols | — | **HYGIENE** |
| **F3** | Wrong locus in four tracked files: stage016's dimension literals are at `:314-325`, not `:355-366` (one spelling uses an en-dash, so a plain grep misses it) | `parameter_register.md:182-184`; `rewrite_reference_table.md:205`; `measure_register_sufficiency.md:100` | **HYGIENE** |
| **F4** | ⛔ **`_validate_loci` checks only that a cited range fits inside the file** — never that the cited lines still say what is cited. A line shift passes every gate **silently** | `reduction/registry_read.py:440-455` | **OPEN** — structural, and the reason F1–F3 went unnoticed |

⭐ **F4 is the one to fix if any.** It is the mechanism by which F1–F3 became possible, and it is the only
hygiene row that is a *live* hazard rather than a historical error.
