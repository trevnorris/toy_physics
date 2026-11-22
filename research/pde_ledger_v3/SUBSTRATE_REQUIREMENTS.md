# Substrate requirements — what the force sectors oblige the unbuilt steps to deliver

**Status: the light sector is closed for now (user decision, 2026-10-19). Passes 1 (S9, S10), 2 (S11,
S11b-A/B, S11c PARTIAL; accepted `22be0137`), O2 (accepted `72866fcf`) and S9b (accepted `ff3bac63`) are
complete. The light-closing pass (2026-10-19), including the oracle-comparison record and throat energy
split, is accepted at `e070c7f3`. The S12 record and register pass (2026-10-20) are accepted at `ac909aaf`;
their constructions and comparator are accepted at `a3da43d9` and `fe485b73`.**
The S14a part 1 record and register pass (2026-10-21) are accepted at `4a4dfaf5`;
their constructions and comparator are accepted at `94fddf2f` and `de81b9fc`.
The S14a part 2 record and register pass (2026-10-21) are accepted at `79511272` (round 2 clear).
By the user's "Record as requirement" decision, no build ran: the slab's own dynamics must select
the slow limit and far behaviour (`R-S14a-01`); S14/S15/S20 remain CONDITIONAL on that requirement.
Twenty-five entries, all OPEN. S16’s record and register pass were accepted at `4f907e9d`.
S16’s constructions and comparator are accepted at `1d990ee6` and `f08db34a`.
The S15 pass adds the conditional gravity observable to `R-S21-01`.
Pass 2 was populated 2026-10-05 from the kept results and their stated
conditions; it does not close the substrate or upgrade S11c's unresolved work. Later sectors remain to
be read as they close. The two 2026-08-07 prior-art entries retain their provenance; the second route is
recorded under Population passes.
O2 adds sources to six existing entries and one separate supplier/budget entry with no owner named;
its four adopted premises retain their
conditional-input labels in the O2 pass below. Entry counts for that pass are retrieved in
`steps/_measurements/O2_record_measurements.md`, M11. S9b adds one S21 requirement and sources
six existing entries; its before/after census is `steps/_measurements/S9b_record_measurements.md`, M10.
The light-closing pass adds one S22 entry and oracle-comparison evidence to `R-S21-01`; its lookups are
`directives/_measurements/light_sector_closeout_lookups.md`.
S12 merges its printed regional/face compatibility objects into `R-S12-02` and adds scoped source or
decision provenance to four existing entries; it adds no entry and closes none. Its retrievals and
before/after census are `steps/_measurements/S12_record_measurements.md`, M12/M15.
S14a part 1 merges projected/shell compatibility into `R-S12-02` and adds scoped
source/provenance to four existing entries, with no new entry or status change. Its census and
pass are retrieved in `steps/_measurements/S14a_1_record_measurements.md`, M12/M16.
S14a part 2 adds `R-S14a-01`, the slow-limit/far-behaviour obligation: 18 → 19 OPEN entries.
Selection is owed jointly by S12, S8/S22 and the face/interior laws with no owner named
(user's V2 decision, `12aa65d5`). No existing entry gains a source or changes target/status;
`R-S12-02` retains its scope at `978c359c`. This pass's retrievals and exact
ancillary diff are `steps/_measurements/S14a_2_record_measurements.md`, M12/M15.
S16 adds six entries (19 → 25 OPEN) and scoped source/provenance to `R-S1-01`, `R-S1-02`,
`R-S12-02` and `R-S21-01`. The user selected the passive-factor and inertial filings;
the list author filed compactness, face momentum, density response and S22's reduction-debt account.
The user's 2026-10-24 "Require all four (Recommended)" decision now supplies admission for those
four list-author filings, settling the record-list/register disagreement; their original filers
and owners remain. The inherited coordinate/support conditions merge into the existing entries;
no physical response is derived. Both the Y1–Y4 repair and this W1 repair are 25 → 25.
Retrievals, both-route applications and the exact ancillary diff are
`steps/_measurements/S16_record_measurements.md`, M01–M19.

## Why this exists

v3 is requirements-first: each sector states what it needs, and the knit happens last. The light sector
(S9, S10, S11, S11b) is built on **S1–S8, none of which have been run** — there are no step records for
them. That is the design, not an oversight: the substrate gets built once every sector has said what it
requires, so it can be checked against all of them at once rather than rebuilt per sector.

The bet only pays if the requirements are **captured**. Measured 2026-08-05: exactly two obligations are
recorded, both from S11, both about one quantity, both as inline anchors in a 1200-line plan. Everything
else a built sector assumes about the substrate is implicit prose spread across four step records.

`DEFECT_REGISTER.md` records what is *wrong*. `ANSATZ_LEDGER.md` records what is *postulated and not
derived*. Neither records **what a later step must deliver for an earlier result to stand.** That is this
file.

## The failure this prevents

At Phase 5 the sectors are knitted and the substrate is built against their combined requirements. If a
requirement was never written down, one of two things happens: the substrate is built without it and the
sector's result silently loses its footing, or it is rediscovered late and the substrate is rebuilt. Both
were survivable with one sector. With five they are not.

The S7 entry below is the shape of the problem. Standing at S7, the honest answer is "the slab width is
not selected" — a real result, bankable, and **insufficient**. Only S11 knows that the question must be
asked again with the width made position- and time-dependent. Nobody at S7 would think to ask.

## Entry schema

| field | meaning |
|---|---|
| `id` | stable identifier, cited from both ends |
| `source` | the step that needs it |
| `target` | the step that must deliver it |
| `requirement` | **the object required**, named — not a derivation path for it |
| `on failure` | what breaks in the source step if the target cannot deliver |
| `status` | OPEN · DELIVERED · RETIRED (with the commit that changed it) |

Name the object; do not prescribe how the target should obtain it. A requirement that specifies a
derivation route manufactures arguments about the route — see `CLAUDE.md` rule 3.

Distinguish carefully:
- a **requirement** is an obligation on a *future* step;
- a **postulate** is a value or form assumed *here* and not derived — that belongs in `ANSATZ_LEDGER.md`;
- a **defect** is something already wrong — that belongs in `DEFECT_REGISTER.md`.

A postulate with a named retirement condition generates a requirement. `B_comp` below is exactly that.

---

## Entries

⭐ **By target step** — what each unbuilt step owes the sectors that are already banked:

| target | id | the object owed | source |
|---|---|---|---|
| **S1** | `R-S1-02` | the substructure's shear response as a function of `χ_B`, both phases | S9, O2 conditional premises, S9b brane optical shear only, S12, S14a part 1 and S16 support identification only |
| **S1** *(register inference)* | `R-S1-03` | whether the substructure's microdynamics is time-reversible | S11b-B / unified S11b |
| **S1.5** | `R-S1.5-01` | the mode content of the linearised GNLS, including the bulk acoustic branch | S9, S11, S11b-A/B, S9b Part C bulk acoustic input only |
| **S6** | `R-S6-01` | the brane's compression modulus `B_comp` | S11 |
| **S6** | `R-S6-02` | the first variation of the substrate action at `v₀ = 0` | S9, S10, S11, S11b-A/B, S11c uniform |
| **S6** | `R-S1-01` | the brane's spatial dimension `D_brane` ⚠ *id kept; target corrected to S6* | S10, S11, S9b, S12, S14a part 1 and S16 coordinate dimension only |
| **S7** | `R-S7-01` | the slab-width flat direction under a non-uniform width | S11, unified S11b |
| **S8** | `R-S8-01` | the form of the brane's quadratic stiffness functional | S9, S10, S11, S11b-B, unified S11b, S11c uniform, S9b supplied transverse input |
| **S8** | `R-S8-02` | the quadratic operator on `u` **and** `h` together | S10, S11 |
| **S8** | `R-S8-03` | the **sign** of the physical transverse stiffness | S9, S10, S11, unified S11b, S11c uniform, S9b local optical sign/ratio gate |
| **S8** | `R-S8-04` | what carries the brane's **internal angular momentum** | S9, S10, S11, O2 carried obligation |
| **S8** | `R-S8-05` | the **frame** the brane's rotational stiffness is measured against | S9, S10, S11, O2 carried obligation |
| **S8** *(register inference for original entry)* | `R-S8-06` | material displacement and original quadratic slab inertia; O2 adds the material identification only | S11, S11b-B, O2 and S9b identification |
| **S12** *(register inference; no live owner named)* | `R-S12-01` | a reservoir and power budget if a non-passive interface response is adopted | S11b-B / unified S11b, O2 carried condition |
| **S12** | `R-S12-02` | native drain/return, printed regional/face and projected/shell compatibility on the real-fraction branch (A13 OPEN), and separate boundary data | S11b-A/B, S11c uniform, O2, S12, S14a part 1; S16 supplied projection/Keep S12 provenance only |
| **S21** | `R-S21-01` | conditional deflection and radar-slope profile requirements, with Part C rewrites and implied exchange; comparison of force sectors' flow requirements | S9b A–C; flow-only oracle comparison (evidence only); S12, S14a part 1 and S16 framing decision only; S15 conditional signed shell-crossing test |
| **S22**, with Q2 and S12 | `R-S22-01` | how the energy a throat releases splits among light, the extra longitudinal wave and the bulk | S11 departure; photon survey E32–E34, E38 |
| **no owner named** | `R-O2-01` | the physical net-power supplier and its stated budget, if a conditional closure requires net supply | O2 |
| **jointly S12** (source/flux forms), **S8/S22** (material law), **no owner named** (face/interior laws) | `R-S14a-01` | the slow limit and far behaviour selected by the slab's own dynamics | S14a part 2; requirement for S14/S15/S20; user V2 decision, `12aa65d5` |
| **no owner named** | `R-S16-01` | surroundings field, separation from the body's own field, and normalization for a passive factor | User's selected "Record what it needs" option; S16 list `:49–53` |
| **no owner named**; S22 receives the debt in its table | `R-S16-02` | compactness relative to both separation and surroundings-field variation scale | Orchestrator's S16 list filing `:158–161`, E3/F3; admission: user's 2026-10-24 "Require all four (Recommended)" |
| **no owner named** | `R-S16-03` | momentum exchanged at each slab face, split by face motion | Orchestrator's S16 list filing `:162–163`, E2; admission: user's 2026-10-24 "Require all four (Recommended)" |
| **S8** | `R-S16-04` | momentum–material relation needed for an inertial factor | User's selected "Count material; mass open" option; S16 list `:54–56`, O2 `:133` |
| **S8/S22** | `R-S16-05` | slab density response | Orchestrator's S16 list filing `:164–166`, O2 `:138`; admission: user's 2026-10-24 "Require all four (Recommended)" |
| **S22** (debt-table account only; physical suppliers not assigned by this entry) | `R-S22-02` | missing response-reduction inputs and width bridge, with compactness cross-referenced to `R-S16-02` | Orchestrator's S16 list filing `:28–30`; plan A11/debt rows; admission: user's 2026-10-24 "Require all four (Recommended)" |

⚠ Entries below are in the order they were found, ⛔ not in step order. All **twenty-five** are **OPEN**.
⭐ `R-S8-01`, `-03`, `-04`, `-05` are one family: the stiffness functional's **form**, its **sign**, its
**mechanical admissibility**, and its **reference frame**. ⛔ Delivering the form does not deliver the
other three.
⚠ `R-S1-01`'s id encodes its original target; the id is stable and the **target** field governs.
⚠ `R-S14a-01`'s id names its source step, not its deliverer; its **target** field governs.

### R-S6-01 — `B_comp` must be retired or re-affirmed

- **source** S11 (`steps/S11_stray_longitudinal.md`, moves 1–5) · **target** S6 · **status** OPEN
- **requirement** — the brane's compression modulus `B_comp` (`Q.brane.B_comp`), as a derived quantity or
  an explicitly re-affirmed postulate.
- **on failure** — S11's knob count stays an upper bound and the longitudinal mode keeps a postulated
  constant underneath it. Entered as a postulated knob on the user's explicit call (2026-08-02) precisely
  so the retirement would be visible when it happens.
- **note** — recorded at `V3_STEP_PLAN.md` `{#s6-b-comp-callback}`, written into S6 on purpose: a note
  living only in S11's record is a note nobody reads on arrival at S6.
- **pass-2 scope** — No S11b record rests on `B_comp`: B's enlarged energy basis is not a derivation of
  `B_comp` from the substrate. Its scorecard limits the
  static series-compliance claim to impermeable faces; with active transfer the static row is `μ_θ = 0`.
  This entry does not identify `B_comp` with `B_eff` or impose a series law on that permeable case.

### R-S7-01 — the slab-width flat direction, under a non-uniform width

- **source** S11 (`steps/S11_stray_longitudinal.md:90–92`); unified S11b
  (`steps/S11b_interface_coupling_law.md:156–157`) · **target** S7 · **status** OPEN
- **requirement** — whether the slab-width flat direction survives when the width is made **position- and
  time-dependent**. A uniform-width statement does not answer it.
- **on failure** — S11 says, *"Two channels accommodating additive shares of one strain, i.e.
  **springs in series**, so compliances add and `B_comp` is **softer than either channel alone**"*, then
  *"If that flatness survived, `B_wall = 0`, hence `B_comp = 0`, and S10's zero would never lift"*
  (`steps/S11_stray_longitudinal.md:84–88`). S11b-B's whole P3 row instead says *"compliances add at
  `ω → 0`"* (prediction) and *"⚠ partial — impermeable only; with active transfer the static row is
  `μ_θ = 0`"* (outcome) (`steps/S11bB_interface_assembly.md:172`). The records disagree
  on the unqualified series-law failure statement; this is unresolved. S11's own answer is that gradients
  lift the flat direction — a wave
  modulates the width, tilting and stretching the interfaces at cost proportional to `σ_wall|∇W|²`, so
  flat at `k=0` and stiff as `k²`. S7 must confirm or refute that.
- **named freeze** — S11 says its cone is *"computed with the wall width FROZEN"*
  (`steps/S11_stray_longitudinal.md:90–92`); the unified record calls the same imported-width identity
  *"a recorded **freeze**, ⛔ not a fix"* (`steps/S11b_interface_coupling_law.md:156–157`).
- **note** — `V3_STEP_PLAN.md` `{#s7-b-comp-callback}`. The charge anchor rests on this.

### R-S1-01 — the brane's spatial dimension

- **source** S16 (coordinate dimension only; `steps/S16_moving_region.md`, §§1/6; S16 M01/M04/M14/M17); S14a part 1 (coordinate dimension only; `steps/S14a_1_projection.md`, §§1–4/8); S12 (coordinate dimension only; `steps/S12_order_source.md`, §§1–2/7); S9b (coordinate dimension only; `steps/S9b_light_bending.md`, §§2–4/8); S10 (`steps/S10_two_transverse_photons.md`); S11
  (`steps/S11_stray_longitudinal.md`, moves 2–3 / finite census) · **target** S6 · **status** OPEN
- **requirement** — the brane's spatial dimension `D_brane`, as a derived quantity or an explicitly
  re-affirmed postulate.
- **on failure** — S10's headline reads *"light having exactly two polarisations is a statement that our
  space is three-dimensional."* Read backwards it says the opposite: `D_brane = 3` went in and `D−1 = 2`
  came out. Without a delivered `D_brane` the sentence is an assumption restated, ⛔ not a result.
- **note** — ⚠ **The target is S6, not S1.** S1 owns `D = 4` and the two-phase split; `D_brane` is a
  property of the **wall**, and S5–S6 are what construct it — S6 (the kink) is the first step at which a
  codimension-1 surface exists to have a dimension.
  ⛔⛔ **An earlier draft of this entry ruled out the obvious route on bad reasoning**, and a leg caught
  it. It argued that because S10's *mode-count computation* contains no codimension, `D_bulk − 1` is not
  the route. ⚠ **Non sequitur.** That the count does not *use* a codimension says nothing about whether
  `D_brane` can be *derived* as one elsewhere — a domain wall in a `D = 4` bulk is exactly codimension 1,
  and that is the natural derivation. ⇒ ⭐ **this entry does not prescribe a route**, per the schema; it
  asks for the object.
- **pass-2 consumer** — S11's selected homogeneous three-mode census and transverse/longitudinal
  separation use `D_brane = 3`. Its proper-rotation invariant count explicitly warns that the sector
  separation is not dimension-independent. This is the same dimension obligation, not a new knob.

- **S9b scope** — The kept optical/mass construction uses three far-field spatial coordinates,
  coordinate `d³x`, spherical `r² sin θ` measure and the `4πr²` flux. It consumes `D_brane=3`
  without deriving it. If that dimension is not delivered or re-affirmed, these coordinate
  functionals and the flux/velocity relation do not transfer unchanged to the substrate.
  This adds no dimensional derivation or selected wall thickness (S9b §§2/5; M1/M19).

- **S12 scope** — The supplied balances and printed regional/face objects have three in-plane
  coordinates and one bulk coordinate. Density's `L⁻⁴` assignment is shown by the SymPy
  construction (`scripts/S12_order_source.py:219–222`) and its printed, unjoined
  `primitive_function_dimensions` table, and by the Wolfram source's `unit` switch
  (`mathematica/S12_order_source.wl:29–37`), not a Wolfram primitive-unit transcript table. This adds a
  coordinate-domain consumer only. If that dimensional setting is not delivered or reaffirmed,
  those objects do not transfer unchanged. No dimension or brane construction is derived
  (`steps/S12_order_source.md`, §§1–2; its retrievals M01/M04/M06/M10/M17).

- **S14a part 1 scope** — The projection uses three in-plane coordinates and one bulk coordinate;
  its labelled spherical shell uses the printed `4πr²` area. This adds a coordinate-domain consumer
  of the existing dimension obligation, not a dimensional derivation or selected thickness.
  Failure to deliver/reaffirm that setting prevents transferring those projected/shell objects
  unchanged (`steps/S14a_1_projection.md`, §§1–4; its retrievals M01/M03/M04/M07).


- **S16 scope** — The supplied balances use part 1's projected domain, "in-plane coordinates
  x^i (three) and time t" (`directives/S16_decision_list.md:63–64`), and the shell geometry
  and measures use that same coordinate dimension. If this setting is not delivered or
  reaffirmed, the printed regional objects cannot transfer unchanged. This is a coordinate-domain
  consumer of the existing S6 obligation, not a dimension derivation; control shape, retained
  moment order and independent dimensional-audit units are separate choices (S16 §§1/6;
  S16 M01/M04/M06–M07/M14/M17).

### R-S1-02 — the substructure's shear response, as a function of the order parameter

- **source** S16 (inherited support identification only; `steps/S16_moving_region.md`, §§1/6; S16 M01/M04/M14/M17); S14a part 1 (support identification only; `steps/S14a_1_projection.md`, §§1/8); S12 (support identification only; `steps/S12_order_source.md`, §§1–3/7); S9b (brane optical shear only; `steps/S9b_light_bending.md`, §§2–4/8); O2 (`steps/O2_steady_brane_balance.md`, register handoff); S9 (`steps/S9_light_requires_shear.md`) · **target** S1 · **status** OPEN
- **requirement** — the shear response of the substructure as a function of `χ_B`, across **both** phases.
- **on failure** — the two halves fail differently, and both are fatal:
  - **ordered phase carries no shear** ⇒ `μ_R = 0`, the brane has no transverse sector, and **photons do
    not exist**;
  - **disordered phase carries shear** ⇒ light is not confined to the brane, it leaks into the bulk — an
    energy sink with no observational room — **and** the throat's trapped brane-shear standing wave
    radiates away, so the outward pressure vanishes and **the geon closes**.
- **note** — ⭐ **Three independent consumers** (photon propagation, photon confinement, geon support), and
  it is **one object**, not three. S9 records it as its own LIVE falsifier and as a knit question: can one
  substructure be ordered-and-shear-bearing in one phase and unstructured-and-shear-free in the other?
  ⚠ Stressed hardest at the throat, where the brane is bent into `±w`. Currently the only machine check
  on the bulk half anywhere in the corpus is **dimensional**.
- **O2 scope** — Premises 1 and 4 adopt optical-regime elastic shear in the brane and the postulated
  shear-free scalar bulk with face-normal mechanical loading. O2 derives neither phase response.
  The existing shear-response retirement question gains this conditional consumer, not a new
  entry for an adopted value. Failure to realize those responses prevents promoting O2's conditional
  material/load object to that substrate. No flowing pressure law, tangential bulk stress or optical
  compatibility result follows (`directives/O2_SHARED_PHYSICS.md`, §§2, 3.1, 3.3, 7).

- **S9b scope** — The supplied light branch is the transverse displacement of shear-carrying
  brane material, with the real positive optical ratio `μ⊥/ρbr`. A kept A–C optical wave is
  conditional on that brane shear response; without it the reported propagating optical
  interpretation fails its own gate. S9b therefore consumes the brane optical half of this
  existing phase-response requirement. It derives no `χ_B` response, bulk shear-free result
  or confinement result; O2 P4 is only carried as an adopted premise, unused by A–C
  (S9b §§2/7; M1/M19).

- **S12 scope** — U1 and its flagged support identification are adopted inputs: order vanishes
  outside the independently bounded shear-holding region, and remains live inside. This adds
  only the conditional material interpretation as a consumer of the existing phase-response
  obligation. Failure of that identification removes the claimed shear-region interpretation;
  the algebra does not establish shear response, confinement or a constitutive relation to
  `χ_B`. U1's re-ordering clause enters no computed object and remains OPEN
  (`steps/S12_order_source.md`, §§1–3; retrievals M01/M03/M04).

- **S14a part 1 scope** — U1's shear-holding-region/support identification is inherited as an
  adopted input, not derived from the projection. It adds only this conditional material
  interpretation to the existing phase-response obligation. If the identification fails, the
  shear-region interpretation fails; no shear modulus, confinement result or profile follows.
  U1's re-ordering clause enters no computed object as an enforced clause
  (`steps/S14a_1_projection.md`, §1; its retrievals M01/M03/M04).


- **S16 scope** — S16 supplies part 1's projected material balance on part 1's domain.
  Its interpretation as material in the shear-holding region inherits the flagged, adopted
  U1 support identification (part-1 record §1). If that identification fails, the interpretation
  fails while the conditional balance algebra remains. This adds a consumer of the existing
  phase-response obligation for support identification only: no shear response, confinement,
  order profile or material law is derived, and U1's re-ordering clause is not enforced
  (S16 §§1/6; S16 M01/M04/M14/M17).

### R-S1.5-01 — the GNLS's own mode content

- **source** S9b (Part C bulk acoustic input only; `steps/S9b_light_bending.md`, §§2–4/8); S9; S11 (`steps/S11_stray_longitudinal.md`, move 6 / `KW_ZERO_LOCUS`);
  S11b-A (`steps/S11bA_interface_response.md`, computed bulk response) and S11b-B
  (`steps/S11bB_interface_assembly.md`, breathing slice) · **target** S1.5 · **status** OPEN
- **requirement** — the mode content of the linearised GNLS about uniform `ρ₀`: which excitations it
  carries and their polarisation.
- **on failure** — S9's opening move is *"a single-component scalar superfluid cannot carry transverse
  light"*, and the entire requirements-first framing of the light sector rests on it. ⛔ It is cited from
  **one external review document** (`decisions/15`) and **no script in this repo executes it**; the S9
  rebuild did not add one. If the GNLS does carry a transverse excitation, light needs no substructure and
  S9 has no bill to present.
- **note** — a weaker in-repo form exists (`v = (ħ/m)∇θ ⇒ ∇×v ≡ 0`) but is never wired into a no-shear
  proof. ⚠ This is a **premise awaiting execution**, not a postulate: it is the kind of claim a CAS can
  settle, and nothing has been asked to.
- **pass-2 consumer** — the same mode-content question now includes the rest-bulk acoustic branch used
  in `q² = ω²/c_s0² − k²`, its polarisation and the regime in which that acoustic approximation applies.
  S11's threshold, A's radiation resistance/added mass, and B's radiative breathing load rest on that
  supplied branch. Failure leaves those results conditional on a bulk operator the GNLS has not yet
  delivered. This requests neither an EOS exponent nor a numerical light/sound speed ratio.

- **S9b scope** — Part C uses the supplied bulk acoustic quantity
  `c_s²=nKρ^(n−1)/m`, including its asymptotic value, and its fixed-ratio response uses the
  exact nonnegative sound-speed roots and their domains. If the bulk does not deliver that
  acoustic branch on the stated domain, the fixed-ratio optical rewrite cannot be read as
  a light/sound response of that substrate. This is a consumer of the existing acoustic-mode
  requirement only: no GNLS mode derivation, EOS exponent, light/sound identification or
  numerical ratio is supplied by S9b (S9b §§3–4; M1/M18).


### R-S6-02 — the reference state must be an equilibrium

- **source** S9, S10; S11 (`steps/S11_stray_longitudinal.md`, homogeneous action);
  S11b-A/B (`steps/S11b_interface_coupling_law.md`, standing limits); S11c
  (`steps/S11c_PARTIAL_CLOSEOUT.md`, conditional uniform result;
  `steps/S11c_d_profile_conditioned_scattering.md`, uniform check) · **target** S6 · **status** OPEN
- **requirement** — the **first variation of the substrate's action**, evaluated at the state S9 and S10
  linearise about: the unstrained brane at rest, `v₀ = 0`.
- **on failure** — every result in the light sector is a linear response, and a linear response about a
  state that is not an equilibrium is not a linear response at all: the expansion has a term linear in the
  displacement that nothing cancels. `ω² = (μ_R/ρ_br)k²`, the mode count, and the dimensions all inherit
  the defect. ⛔ Nothing in the sector tests it, because the substrate is absent from S9's and S10's
  actions.
- **note** — ⚠ **The standing warning applies here**: ask *whose* equilibrium, and check that the
  reference state is not being held in place by an assumption rather than by the dynamics. S6's kink is
  the candidate stationary solution; S6 has no record yet.
- **pass-2 scope** — the kept homogeneous and selected uniform results use the rest reference state.
  S11b-B states verbatim: *"⛔⛔ **Passivity and reciprocity are properties of fluctuations about an
  EQUILIBRIUM. This model's reference state is not one.** It carries a steady background normal flow
  `v₀` — the dark-energy drain — ⇒ a **driven steady state with a reservoir attached.**"*
  (`steps/S11bB_interface_assembly.md:46–48`). In that quote S11b's `v₀` is the bulk's normal drain,
  standard name `v_bulk_normal_0`; it is distinct from the brane rest velocity `v₀ = 0` in this requirement
  (`directives/S11b_SHARED_PHYSICS.md:101–104`; `directives/S11b_unified_decisions.md:33–40`). For S11b, the
  rest-frame linearisation discards a relative correction of order `O(v₀|q_n|/ω)`, where `v₀` again denotes
  the steady background normal flow (`steps/S11bB_interface_assembly.md:150–154`), and the unified record
  calls that correction uncarried
  (`steps/S11b_interface_coupling_law.md:158–164`). The kept S11c results state that the bulk is at rest
  (`steps/S11c_PARTIAL_CLOSEOUT.md:7`; `steps/S11c_d_profile_conditioned_scattering.md:9,23`). This is an
  unresolved tension between the driven reference state and the recorded rest-state calculations. The
  rest-state obligation remains here; live drain/return data are `R-S12-02`, which also records the
  uncarried correction and its record-named owner, and the power supply is `R-S12-01`. No equilibrium result
  is promoted to a result about the driven state.

### R-S8-01 — the form of the brane's quadratic stiffness functional

- **source** S9b (`steps/S9b_light_bending.md`, §§2/8, supplied transverse input only); S9, S10; S11 (`steps/S11_stray_longitudinal.md`, moves 2–3 / FORM control);
  S11b-B (`steps/S11bB_interface_assembly.md`, energy basis / breathing stability); unified S11b
  (`steps/S11b_interface_coupling_law.md`, energy quotient); S11c
  (`steps/S11c_PARTIAL_CLOSEOUT.md`, conditional uniform result) · **target** S8 · **status** OPEN
- **requirement** — the quadratic brane Lagrangian's **stiffness functional**, as delivered by the
  substructure rather than chosen.
- **on failure** — this is the one input the light sector's central claim is **measurably** sensitive to,
  and S9's own controls prove it:
  - the gradient-elastic form `−½ μ_R Σ(∂_i u_j)²` — an ordinary elastic solid — carries **the same two
    transverse modes at the same `c² = μ_R/ρ_br`**, and *also* propagates the longitudinal;
  - the divergence-only form makes the roles **swap**: the transverse pair drops to `ω² = 0` and the
    longitudinal propagates.
  ⇒ ⭐ **What curl-only buys is the ABSENCE of a propagating longitudinal mode, ⛔ not the presence of the
  transverse ones.** If S8 delivers any other form, S9's and S10's transverse results may survive while
  the no-longitudinal claim — Maxwell's third demand, and the whole reason the sector exists — does not.
- **note** — S9 classifies the curl-only form as **postulated (structural)**, justified by *"forced by
  Maxwell's no longitudinal mode"*. ⚠ That is a justification from the target, ⛔ not from the substrate.
  Defect `B2` closes the route from a polar substructure `P` to `μ_R` — ⚠ **the whole quantity, not only
  its magnitude**, as `DEFECT_REGISTER.md#B2` scopes it. ⛔ But it closes **one route**, and it says
  nothing about the **form**, so this requirement is live.
- **pass-2 consumers** — S11's unchanged transverse root rests on adding the trace invariant. As a
  first-route source only, its record says: *"The form control's `4/3` independently reproduces the
  corpus's rejected-Cauchy-branch coefficient"* (`steps/S11_stray_longitudinal.md:126–127`); the record
  names no prior-art objection to that branch.
  S11b-B needs the complete quadratic stored energy on `u, θ, e_W` under its whole stated symmetry group
  (`directives/S11bB_SHARED_PHYSICS.md:353–362`). The unified record's ten-dimensional quotient is likewise
  defined under its whole stated symmetry group (`directives/S11b_SHARED_PHYSICS.md:280–288`) and has different
  valid representatives; no individual representative coefficient is the required object. Its
  constrained breathing stiffness is `K₀ = B_ρ⁽³⁾ − 2CW₀ + k_W W₀²`. The stability claim is only the
  `k = 0`, impermeable, zero-reciprocal-traction slice. On that slice B hands forward `K₀ > 0` as
  *"an **empirical** constraint — the brane is here"* (`steps/S11bB_interface_assembly.md:192`). Failure
  to supply this energy leaves the
  decoupling and slice stability conditional; no sign or magnitude of an individual scalar coefficient
  is inferred. S11c's uniform result inherits that qualification, not a nonuniform confinement law.

- **S9b scope** — Parts A–C assume an isotropic transverse spatial stiffness with local
  `cγ²=μ⊥/ρbr`, on their supplied leading-eikonal optical domain. S9b therefore gains a source on
  the existing requirement for a substrate-delivered transverse stiffness input. It neither derives
  that input nor widens this entry into a completed live stress/relaxation or nonuniform operator law.
  The latter stays the unfinished S11c handoff; other unselected material response forms are not
  newly entered. No curl-only, no-longitudinal, mechanical-admissibility or physical-frame result
  follows from the optical comparison (`steps/S9b_light_bending.md`, §§2/6–8; measurements M1/M8).

### R-S8-02 — the in-plane and out-of-plane sectors must decouple at quadratic order

- **source** S10; S11 (`steps/S11_stray_longitudinal.md:278–280`, in-plane, frozen-wall-width spectrum)
  · **target** S8 · **status** OPEN
- **requirement** — the quadratic operator on the brane's **full** displacement, in-plane `u` **and**
  out-of-plane `h` together, and whether it is block-diagonal in that split.
- **on failure** — ⛔ **S10's headline number changes** — but ⚠ **not by the mechanism an earlier draft of
  this entry named, and a leg corrected it.** Mixing is **not** the failure mode: under the rotation group
  that acts on the brane, `h` and `u_L` are both **scalars** and `u_T` is a **vector**, so `h` can mix with
  `u_L` freely and the transverse count stays `D − 1` regardless. ⭐ **The condition that gives 3 is `h`
  being DEGENERATE WITH THE TRANSVERSE PAIR** — the brane's elasticity being isotropic in `D+1` rather
  than in `D`. That is what *"belonging to the same elastic sector"* has to mean, and it is what S8 must
  settle.
- **note** — ⭐ **This is the most load-bearing open item in the sector**, and it was flagged inside S10 at
  the moment the identification was made rather than found later: `h ≠ u_L`, stated in three places in the
  corpus, **user-confirmed as the picture** — but ⛔ confirmed as a picture, not computed. `V3_STEP_PLAN`
  puts *"transverse and longitudinal sectors; the reduced `h`/`u_L` operator"* in S8, so the object is
  already scheduled.
  ⚠⚠ **And v3 is not starting from nothing — a leg found the computation already exists in v2 and this
  entry did not cite it.** `research/pde_ledger_v2/paper/stages/stage_030.tex:99-125` builds the coupled
  `(u_L, h)` scalar block with stiffness `K = [[B_eff, C_hu], [C_hu, K_h]]` — explicitly **not**
  block-diagonal, with `C_hu` a registered free-unreduced parameter — and `R79` records that the mixed
  poles are cone-coincident only when `C_hu = 0`. ⇒ ⛔ *"nothing computes that"* is true of **v3** and
  misleading about the **corpus**; S8 should start from `stage_030`, ⛔ not from scratch.
  ⭐ Note this is the `(u_L, h)` block — the scalar sector — so it bears on the **longitudinal** slot and
  the charge anchor, ⛔ and not directly on the transverse count, which is what the corrected on-failure
  above turns on.
- **pass-2 source** — S11's result here is the one its record names as the in-plane, frozen-wall-width
  spectrum with `hBranon` excluded: *"⚠ Scope caveats, none touching the decoupling: (i) this is the
  **in-plane, frozen-wall-width** spectrum (`WALL_WIDTH_FIELDS={}`, `hBranon` excluded,
  `INTERFACE_EQUATIONS_SUPPLIED={}`) — it does not decide inhomogeneous mode conversion,"*
  (`steps/S11_stray_longitudinal.md:278–280`). S11's decoupling statement, with its own limit, reads:
  *"Thus the homogeneous \(D=3\) quadratic transverse and longitudinal eigenbranches have zero linear
  cross-block under the selected action. This does not establish decoupling at nonlinear order, on a
  nonuniform slab, at an interface or defect, or after additional allowed fields are introduced."*
  (`steps/S11_stray_longitudinal.md:74–76`). This entry's object is the joint `u`/`h` operator that the
  spectrum leaves out.
- **excluded field** — in that caveat S11's record freezes the wall width and excludes `hBranon`
  (`steps/S11_stray_longitudinal.md:279–280`); its Mathematica audit binds the token as
  `"OUT_OF_PLANE_FIELD_EXCLUDED" -> hBranon`
  (`mathematica/S11_stray_longitudinal_mathematica_audit.wl:1082`). If S8 does not deliver the joint
  operator, S11's spectrum is not extended beyond that stated scope. The record does not say what its
  roots or census would be with `hBranon` included, and this register does not infer it.

### R-S8-03 — the SIGN of the physical transverse stiffness

- **source** S9b (local optical sign/ratio gate; `steps/S9b_light_bending.md`, §§2–4/8); S9, S10; S11 (`steps/S11_stray_longitudinal.md`, unchanged transverse branch);
  unified S11b (`steps/S11b_interface_coupling_law.md`, transverse mode); S11c
  (`steps/S11c_PARTIAL_CLOSEOUT.md`, conditional uniform result;
  `steps/S11c_d_profile_conditioned_scattering.md`, uniform check) · **target** S8 · **status** OPEN
- **requirement** — the **sign** of the physical transverse stiffness in the quadratic brane
  Lagrangian, as delivered by the substructure: `μ_R` in S9/S10's selected action, `μ_⊥` in the enlarged
  S11b action. `R-S8-01` asks for the stiffness **functional**; this asks for its transverse sign.
- **on failure** — ⛔ **the transverse sector is two exponentially growing modes rather than two waves,
  and every mode count in S9 and S10 is unchanged.** S10's `XFORM_SIGNFLIP` control measures exactly this:
  with the sign flipped, `ω² = −μ_R k²/ρ_br` and **every nullity is identical to the baseline's** — the
  count cannot tell a wave from an instability. What distinguishes them is a single emitted object,
  `ROOT2_Q3_SIGN`.
- **note** — ⭐ **Found by a review leg reading the register against S10's own controls**, ⛔ not by the
  pass that wrote the register: `R-S8-01` is entirely about form and never mentions sign or positivity,
  so the pass consolidated one and dropped the other.
  ⚠ **It is not closed by defect `B2`**, which shuts the route for deriving `μ_R` — its sign as well as
  its magnitude — from a polar substructure `P`. Other substrate routes remain open: a step that
  delivers the stiffness functional delivers its sign with it, so the retirement condition is exactly
  as live as `R-S8-01`'s.
- **pass-2 object** — in the enlarged S11b action this is the physical transverse stiffness `μ_⊥`:
  WL-`μ_R` and SymPy-`(μ_R + μ_S/2)` are its two representatives under the unified record's invertible
  coefficient map. The obligation concerns the sign of `μ_⊥`, not a basis-dependent sign of `μ_S` or
  either engine's bare `μ_R`. Stable real-frequency uniform modes require `μ_⊥ ≥ 0`; decoupling alone
  does not exclude growth. S11c's selected finite-current check keeps its tested-input and review scope.

- **S9b scope** — Every kept A–C observable is gated on a real positive optical ratio
  `μ⊥/ρbr`, with the nonnegative speed chart and its engine-specific ray conditions.
  The conditional wave/ray interpretation fails outside those gates; supplying the stiffness
  form alone does not supply this sign condition. S9b is therefore a source on this existing
  physical transverse-sign requirement. Its local strict-positive ratio domain is kept
  separate from the earlier uniform `μ⊥≥0` statement. PY declares positive `ρbr`; WL keeps
  real nonzero `ρbr`, so this pass does not infer unconditional `μ⊥>0` for both engines.
  No sign is derived from the substrate and no existing status closes
  (S9b §2; M1/M19, complete restriction domains M18).


### R-S8-04 — what carries the brane's internal angular momentum

- **source** O2 (`steps/O2_steady_brane_balance.md`, register handoff); S9, S10, S11 · **target** S8 · **status** OPEN
- **requirement** — the object in the substructure that carries **internal angular momentum** on the brane,
  or the couple-stress it supports. ⛔ Not a mechanism, ⛔ not a model: the object, or a statement that there
  is none.
- **on failure** — ⛔⛔ **the curl-only stiffness functional is not an admissible continuum mechanics.** An
  energy in `(∇×u)²` alone has an **antisymmetric** Cauchy stress, and balance of angular momentum forces
  the Cauchy stress to be **symmetric** unless the medium carries distributed couples or internal spin. If
  the substructure supplies neither, the light sector's central form is inadmissible **regardless of its
  mode content** — S9's and S10's mode counts, dimensions and speeds would all be computed from a
  functional no medium can have.
- **note** — ⚠ **This is the objection that sank MacCullagh's aether**, and it is the one part of that
  theory the 19th century never answered: Stokes pressed it, and Kelvin's gyrostatic models were attempts to
  supply exactly this object. ⇒ ⭐ **prior art is the oracle here** — it tells us the obligation is real and
  that answers exist, ⛔ it tells us nothing about whether **ours** delivers one ⇒ `CLAUDE.md` rule 16.
  ⭐ Known families a delivered answer might fall into, ⛔ **none assumed and none prescribed**: continua
  with an independent microrotation degree of freedom (Cosserat/micropolar), and media with stored internal
  angular momentum.
  ⛔ **One family is the wrong one and should not be reached for:** modern *odd elasticity* buys an
  antisymmetric modulus tensor by making the solid **active and non-conservative**. MacCullagh's medium is
  **conservative** — it has a genuine energy functional — so a non-conservative realisation would be
  answering a different question.
  ⚠ `R-S8-01` asks for the **form** and is silent on admissibility; a substructure could deliver the
  curl-only form and still owe this.
- **O2 scope** — OPEN `𝒜_rot^live` carries this existing S8 obligation, including whether internal
  angular momentum/couple stress is present. O2 selects no stress symmetry or carrier and does not
  establish a curl-only live stress. The original conservative-MacCullagh objection keeps its
  original domain; it is not newly proved for an unspecified viscoelastic stress. The O2 addition
  carries the explicitly required carrier/admissibility question within its named material and
  corresponding power accounting; it adds no flowing constitutive response form. The accounting
  remains conditional with that question OPEN
  (`directives/O2_SHARED_PHYSICS.md`, §§3.2, 4, 6).

### R-S8-05 — the frame the brane's rotational stiffness is measured against

- **source** O2 (`steps/O2_steady_brane_balance.md`, register handoff); S9, S10, S11 · **target** S8 · **status** OPEN
- **requirement** — the frame with respect to which the brane's rotational stiffness is defined: **what
  `∇×u` is measured against.**
- **on failure** — for an infinitesimal rigid rotation `u = ω × r` the curl is `2ω ≠ 0`, so a curl-only
  energy is **nonzero when the whole medium is turned**. ⇒ if the reference is external, the medium knows
  an absolute orientation, and every result in the sector inherits a preferred-orientation signature that
  something later must hide. ⛔ That is a falsifiable consequence, ⛔ not a philosophical discomfort.
- **note** — ⭐ **A brane may escape this where a bulk aether cannot**, and the difference is the whole
  point of the entry. For a bulk medium *"rotation relative to what"* has no local answer — the 1839
  objection. For a **domain wall**, the wall supplies its own local frame (its normal, its induced
  geometry), rotation relative to the wall is meaningful and local, and rotating everything rotates the wall
  too ⇒ the stiffness would be **relative**, and the objection would not arise.
  ⚠⚠ **That is a HYPOTHESIS, ⛔ not a result.** Nothing in the corpus computes it. The slab-in-bulk
  calculation is where it is settled ⇒ **S11b-A** (interface response) is the existing partial artifact,
  built under the old pattern and on the rebuild list.
  ⚠ **A flowing medium does NOT answer this.** Flow bears on the **velocity** frame — whether a background
  drift is detectable — and the analog-gravity route hides it in an effective metric. ⛔ Orientation is a
  **separate** objection: being carried is not being turned. ⇒ ⭐ do not let one argument discharge both.
- **O2 scope** — OPEN `𝒜_rot^live` retains the physical rotational frame as an existing S8
  obligation explicitly carried by the spec. Neither radial flow, the graph normal nor LAB_HELD
  chooses that material frame. O2's named force/power accounting retains the frame question OPEN;
  this addition supplies no frame or new rotational-response form. O2 does not extend the original
  curl-only objection to an unspecified live stress
  (`directives/O2_SHARED_PHYSICS.md`, §§1, 3.2, 6, 7).

### R-S1-03 — the substructure's microscopic time-reversibility

- **source** S11b-B (`steps/S11bB_interface_assembly.md`, limits of the passive region), unified S11b
  (`steps/S11b_interface_coupling_law.md`, conditional Onsager–Casimir test) · **target** S1
  (**register inference; no owner named by the records**) · **status** OPEN
- **requirement** — whether the substructure's microdynamics is time-reversible, the premise of the
  conditional Onsager–Casimir relation `Λ_X(ω) = −Λ_V(ω)`.
- **on failure** — that relation cannot be imposed on the physical interface as an unconditional law.
  The conditional calculation and the independently computed passivity region still stand; neither
  selects a reciprocal medium. The records explicitly say microscopic reversibility is not postulated.
- **note** — this is a first-route obligation from B's hand-forward: `Λ_X(ω) = −Λ_V(ω)` is
  *"**conditional** on microscopic time-reversibility, which the model does not postulate"*
  (`steps/S11bB_interface_assembly.md:190`). The records name no step that owns the missing premise, so S1
  is a register inference. Onsager–Casimir is a second-route candidate only, not a reproduced prior result.
  The equilibrium issue is `R-S6-02`; the driven state and power budget are separately `R-S12-02` and
  `R-S12-01`.

### R-S8-06 — material displacement and the slab's quadratic inertia

- **source** S9b (material identification only; `steps/S9b_light_bending.md`, §§2/8); O2 (material identification only; `steps/O2_steady_brane_balance.md`, register handoff); S11 (`steps/S11_stray_longitudinal.md`, move 1 / finite census); S11b-B
  (`steps/S11bB_interface_assembly.md`, breathing quadratic) · **target** S8
  (**register inference for the original S11/S11b records; O2 names S8 for the material identification**) · **status** OPEN
- **requirement** — `u` as the material displacement of the stuff whose density is `ρ_br`, and the
  quadratic kinetic form of that material and the thickness degree of freedom (B's `μ_W`) on the
  original homogeneous/slab domain. O2 is a consumer of the material identification only.
- **on failure** — if `u` is a director rather than that displacement, S11's continuity identification
  `δρ_br = −ρ_br ∇·u` and its compression argument do not follow. Its finite census requires nonzero
  `ρ_br`; B's breathing-root interpretation also rests on the stated inertial model. A stiffness
  functional alone (`R-S8-01`) does not supply this identification or the kinetic form.
- **note** — this asks for the field identity and inertia, not numerical benchmark values of `ρ_br`
  or `μ_W`, and does not identify the thickness mode with S10's out-of-plane displacement. The cited
  records name no future owner for this object; S8 is a register inference
  (`steps/S11_stray_longitudinal.md:32–38`; `steps/S11bB_interface_assembly.md:80–85`).
- **O2 scope** — O2 accounts for the material whose displacement is `u`, density is `ρ_br` and live
  background velocity is `V`. This recorded/supplied identification fixes which material the named
  conditional balance describes; without it the stated material object and its force/velocity
  pairings are not identified. It gains O2 as a consumer on that identification alone
  (`directives/O2_input_contract.md:55–71`; `directives/O2_SHARED_PHYSICS.md:107–109`).
  The original homogeneous quadratic/slab and thickness inertia obligations keep their domains.
  O2 adds no flowing/embedded inertial-response requirement: `ℐ_br^live` and the momentum map remain
  unselected OPEN closure inputs, treated by the same criterion as O1 and O7 in the pass below.
  No live identification of momentum density with `ρ_br V` follows from the homogeneous anchor.
  S8's original assignment remains a register inference; O2's input contract names S8 for material
  identity while supplying no new live kinetic law.

- **S9b scope** — The supplied advection and coordinate mass law identify `u`, `ρbr` and `V`
  with the same material; this is a condition on the reported local optical object. Failure of the
  identification removes that object's material interpretation. S9b adds this source only: P5
  momentum density/carriage is scoped to held Part D, and Parts A–C derive no live kinetic/inertia
  law. The original homogeneous/slab inertia domain is unchanged (S9b record §§2/7; measurements M1).

### R-S12-01 — the reservoir and its power budget

- **source** O2 (`steps/O2_steady_brane_balance.md`, register handoff); S11b-B (`steps/S11bB_interface_assembly.md`, standing rule); unified S11b
  (`steps/S11b_interface_coupling_law.md`, passivity region) · **target** S12
  (**register inference**) · **status** OPEN
- **requirement** — a named reservoir and a stated power budget for any adopted non-passive interface
  coupling.
  The unified record names the candidate: *"The model supplies a candidate reservoir — the
  background drain `v₀`."* (`steps/S11b_interface_coupling_law.md:62–63`). That `v₀` is the bulk's normal
  drain, standard name `v_bulk_normal_0` (`directives/S11b_SHARED_PHYSICS.md:101–104`;
  `directives/S11b_unified_decisions.md:33–40`). B says the reference state carrying that flow is *"a
  **driven steady state with a reservoir attached.**"* (`steps/S11bB_interface_assembly.md:46–48`).
- **on failure** — the finite-memory velocity channel outside the passive region cannot be inherited
  as a physically supplied response. The computed region remains a classification, not a prohibition,
  and naming `v_bulk_normal_0` alone does not supply the missing power.
- **note** — this is the records' explicit condition for retaining that channel, not a requirement
  that the model choose it. A numerical observational bound also needs matter-to-compression coupling;
  B says that is unbuilt and reports only a structural test, so no numerical-bound requirement is added.
- **owner/target** — the records do not name S12; that target is a register inference. They hand the
  reservoir-and-budget condition to *"S11b-C"* (`steps/S11bB_interface_assembly.md:195–197`).
- **O2 scope** — O2 carries this same conditional non-passive-interface obligation: if a later
  closure adopts that response, name its reservoir and state its power budget. No such response is
  selected. Its live owner is **no owner named**; the historical S11b-C → S11c handoff has no recorded
  live successor, and S12 remains only the original register inference. This object is separate
  from S12's additional non-variational source partners, premise 1's relaxation power, and the
  net-supplier/budget obligation `R-O2-01`
  (`directives/O2_SHARED_PHYSICS.md:312–317`; `directives/O2_input_contract.md:459–465`).
- **S12 pass** — Carried separately in the S12 inventory, without a new computed power consumer.
  U3's endless material supply in the bulk IR is not a named reservoir or a power budget. The
  source compatibility objects deliver neither. The live owner remains **no owner named**;
  this does not assign the whole obligation to S12 (`steps/S12_order_source.md`, §§3/5/7).

### R-S12-02 — background drain/return and separate boundary data

- **source** S16 (supplied projection and "Keep S12" provenance only; `steps/S16_moving_region.md`, §§1/6; S16 M01/M04/M14); S14a part 1 (`steps/S14a_1_projection.md`, §§2–5/7–8; its retrievals M03/M05–M10); S12 (`steps/S12_order_source.md`, §§1–2/5–7; retrievals M01/M03/M07/M11); O2 (`steps/O2_steady_brane_balance.md`, register handoff); S11b-A/B (`steps/S11b_interface_coupling_law.md:158–164`, background-flow limit); S11c
  (`steps/S11c_PARTIAL_CLOSEOUT.md:7`, uniform result;
  `steps/S11c_d_profile_conditioned_scattering.md:21–25`, selected uniform check)
  · **target** S12 · **status** OPEN
- **requirement** — the native drain/return functions and their separate boundary data; the
  printed regional order-source objects and normalized face-source objects, with the separately
  printed total-material face balances, as compatibility requirements on conversion, face and
  order-flux forms. The objects are `general_order/region/{1,2,3}/source`,
  `general_order/face/{1,2}/{trace_coefficient,source_contribution}` and
  `material/face/{1,2}/requirement`. The distinct labelled conditional specialization is
  `conditional_fast_ordering/region/{1,2,3}/source` and its two face objects. Their full values,
  domain and cross-engine limits are S12 §2/4 and retrievals M03/M07–M09.
  **Branch condition:** these printed compatibility requirements rest on the **real-fraction
  branch; A13 stays OPEN**. The decision list says **"The committed law is written for a real
  ordered fraction; a phase-carrying branch would add content not computed here"**. The plan
  requires **"the conversion balance must use the branch chosen at S1"**. The printed objects
  remain on the real-fraction branch until S1's gate decides; the owners remain S1, S5 and S12
  (S12 record §6; retrievals M18/M19). This qualification travels with the future-form handoff.
- **on failure** — the kept selected uniform check cannot be transferred to a flowing background: it is
  a strict-rest-bulk calculation, and its finite samples do not prove flowing-background behaviour
  (`steps/S11c_d_profile_conditioned_scattering.md:21–25`). S11b records why the rest-frame limitation
  matters: its linearisation discards a relative correction of order
  `O(v₀|q_n|/ω)` (`steps/S11bB_interface_assembly.md:150–154`), where that record's `v₀` is the steady
  background normal flow, standard name `v_bulk_normal_0` (`directives/S11b_SHARED_PHYSICS.md:101–104`;
  `directives/S11b_unified_decisions.md:33–40`). The unified record calls the correction uncarried and
  unbounded (`steps/S11b_interface_coupling_law.md:158–164`) and carries the limit
  to *"C (and to the nonlinear program for the DC/harmonic/sideband radiation audit)"*; S11c-a inherited
  it as a standing limit but *"computes geometry only and does not address it"*
  (`steps/S11c_a_interface_shape_derivatives.md:257–260`).
- **owner** — S12 owns the drain/return data (`steps/S11c_PARTIAL_CLOSEOUT.md:21`) and is next for
  dynamical bulk-to-brane order conversion with separate boundary data
  (`steps/S11c_PARTIAL_CLOSEOUT.md:33`). Under the user's 2026-10-20 **"Keep S12"** decision,
  S12 stays OPEN for these forms: this record closes only its requirements part and delivers
  none of the forms. They come once the brane is defined, after the force sectors (S12 record,
  §§5–7; retrieval M18). The earlier assignment/framing disagreement is quoted there.
- **O2 scope** — Premise 2 makes dynamical order conversion the drive, with no separate body force.
  O2's local source/controller and boundary/domain inventories remain distinct, as do additional
  S12 momentum/energy partners and their reaction/supply systems. O2 adds a consumer for the stated
  drain identity and distinct source/boundary accounting; without those roles its adopted drive is
  unidentified or conflated with another load/partner. The original native drain/return functions
  retain their sourced domains, and O2 demands no new functional forms for them. Local in-plane
  carried `j_n V^i` does not deliver
  these data, the live `𝒥_map`, bulk-direction carry, energy partners or a `GM` coupling
  (`directives/O2_SHARED_PHYSICS.md`, §§3.3, 5, 6, 10). Those remain OPEN handoffs, not extra entries
  merely because a later profile solve would need them.
- **S12 scope and rest-on condition** — The conditional result rests on these printed compatibility
  objects, not on any chosen means of meeting them. Below and above, the order source is the
  regional order-flux divergence; inside it includes live material transport of order; each face
  carries the printed source coefficient and a distinct material requirement. If candidate data
  violate one of those objects on the stated domain, they cannot realize that conditional result.
  Source, face and flux forms are not selected, and separate boundary data are not supplied. The
  native-data obligation remains OPEN under the user's 2026-10-19 requirements-first framing.
  The printed constraints are not its delivery. The fast-ordering object is conditional, not an
  earned limit; concentrated input fields are outside the current object. The six unrestricted
  coefficient differences are labelled representation diagnostics confined off their respective
  delta supports, not open differences in the S12 distribution (contract item 3; retrieval M17).
  Both raw expressions and every unpaired/ledger leaf limit remain in S12 §4 and its retrievals.
  This merges source/face/flux compatibility into the existing native-conversion obligation,
  without duplicating an entry or prescribing a law, material, premise or candidate.

- **S14a part 1 projected/shell compatibility** — The required objects also include the printed
  projected material and order balances, their densities/in-plane fluxes, both interior and
  closed-projection face terms, concentrated source parts and projected order source; in the
  labelled spherical specialization, the printed face/core channels, `J`, signed order lost,
  storage, shell balance, face-channel-minus-order-loss difference, outer flux/radius dependence
  and the conditioned projected-velocity/volume-flux relations. These are the object requirements
  listed by native tag/binding in the S14a part 1 pass below. Their full values/domain are
  `steps/S14a_1_projection.md`, §§2–4, with all cross-engine differences/limits in §5 and retrievals
  M03/M05–M10. Future data violating them cannot realize the reported conditional projection.
  The requirements retain the real-fraction/A13 condition above, the live two-face function class,
  full concentrated parts, spherical restriction where used and nonzero-density velocity domain.
  No means of meeting them is named. S12's forms and separate boundary data remain OPEN under
  "Keep S12"; these constraints deliver neither. Unresolved exact differences and unpaired
  decompositions are not silently promoted to exact-equivalence requirements.

### R-O2-01 — the net-power supplier and its stated budget

- **source** O2 (`steps/O2_steady_brane_balance.md`, energy accounting and register handoff)
  · **target** no owner named · **status** OPEN
- **requirement** — the physical supplier of net power and its stated power budget, represented by
  the OPEN `𝒮_E,net` / `𝒫_E,supply` inputs, accompanying any conditional closure requiring net supply.
- **on failure** — such a closure cannot be claimed physically supplied. Naming the drain alone
  supplies neither available energy nor a budget; O2 measures no net supply and imposes no passive
  sign or relaxation law.
- **owner** — the spec calls this an **O2 requirement; source/holder/supplier forms retain their
  owners** (`directives/O2_SHARED_PHYSICS.md:282`). It names no future owner of the net-supplier/budget
  object. S12 source forms and Q2/S22 holder forms keep their separate assignments; none is assigned
  this entire obligation.
- **scope** — this is an accounting requirement on conditional relations requiring net supply,
  fixed by adopted premise 1, not a derived supplier, budget value or additive energy-channel law.
  It is separate from the conditional non-passive-interface reservoir obligation `R-S12-01` and
  from explicit reference/relaxation power with its unassigned reference owner
  (`directives/O2_SHARED_PHYSICS.md:270–282, 307–317`; `directives/O2_input_contract.md`, §8).
- **S12 pass** — The energy-partner inventory carries this obligation separately from
  `R-S12-01` and from conversion partners. No S12 object computes net supply or consumes a
  conditional net-supply closure, so it adds no computed source and changes no owner or status
  (`steps/S12_order_source.md`, §§3/5/7).

---


### R-S21-01 — conditional deflection and radar-slope profile requirements

- **source** S9b (`steps/S9b_light_bending.md`, §§3–6; measurements M3–M7/M11–M13/M17–M19);
  S16 framing decision only (`steps/S16_moving_region.md`, §§1/6; user 2026-10-19,
  carried in its decision list 2026-10-21; S16 M01/M14);
  S12 framing decision only (`steps/S12_order_source.md`, §§6–7; user 2026-10-19, retrieved M01);
  S14a part 1 framing decision only (`steps/S14a_1_projection.md`, §§1/7–8; user 2026-10-19,
  carried in its decision list 2026-10-21, retrieved M01)
  · S15 conditional signed shell crossing (`steps/S15_two_drain_region.md`, §§2–5)
  · **target** S21 (integration/knit; S9b record §7, O2 record §7 and the record directive's
  2026-10-08/09 user decisions) · **status** OPEN
- **requirement** — On S9b's isolated spherical, steady far-zone, leading-eikonal, LAB_HELD
  domain and selected G100/G010/G020/G001 grades, the printed Part B profile condition for
  **each** observable, deflection and radar logarithmic slope, relative to independent orbital
  `GM`, for every far-zone impact parameter. PY prints for each
  `GM=r[2c₀²rδ′+c₀²ξw′²−2rVV′+V²]/4` as a radial identity. WL prints its `GM` anchor rule
  plus the all-radius profile equation (S9b Appendix A), with derivatives, convergence, Abel
  invertibility and constructed-ray/branch/traversal domains. Their complete cross-engine
  equivalence is not established; neither form is discarded.
  Part C carries the separately printed constant-speed condition
  `GM=r[c₀²ξw′²−2rVV′+V²]/4`, fixed-ratio condition
  `GM=r[c₀²(n−1)rf′+c₀²ξw′²−2rVV′+V²]/4`, and power condition
  `GM=r[2c₀²srf′+c₀²ξw′²−2rVV′+V²]/4`, with WL's corresponding anchor/profile equations.
  Fixed-ratio and power objects retain both `ρ>0` and `ρ≤0` labels and exact response domains;
  `ρ₀>0` remains supplied. The first-order-f expressions supply no finite-density or f–ε extension.
  Each condition includes its printed implied `jn`, with `ρbr` symbolic and differentiated:
  the coordinate mass expression is `jn=−ρbr V′−Vρbr′−2ρbr V/r`, with both signed-velocity
  components, zero-locus differentiability/flux joins, excluded nonreal branches and printed
  integration-constant domains retained. The complete implied-exchange values remain uncompared.
  The requirements keep δ, ξw, ρbr, f, n, s, GM and all remaining printed inputs free;
  they select no profile, exponent, flux value/sign or GM normalization. M13 retrieves all 14
  within-engine deflection/radar condition pairs: every PY paired payload is byte-identical,
  while every WL paired payload differs (M13).
  WL `GM`, `ProfileEquation`, `ProfileEquationRadius`, outer `Domain`/`BranchDomain` and case
  `Domain` match literally; `InverseWeight` and `ReferenceInverse` differ in every pair.
  `VelocityReduction/SquaredVelocityEquation` additionally differs for general B, R_FLOW_ONLY B
  and combined C_CONSTANT/C_FIXED_RATIO/C_POWER; `VelocityReduction/Cases/0/Domain` differs
  for R_SPEED_ONLY B, R_TILT_ONLY B, both R_BULK_ONLY C pairs and all conditional forward B/C
  pairs. These paths are relative to `Cases/0`, within each density part's `Object` and the
  outer forward `Object` where present. M13 retains every identical subtree and differing
  payload; density labels/label conditions, bulk premise and forward premise match within
  their pairs. Whether the differing fields state the same condition is OPEN, no owner named,
  separately from the item-5 counting decision below. The matching GM/profile fields alone
  do not settle full-condition equivalence. (M13.)
  Also owed is the **comparison of force sectors' flow requirements**, as the user's
  2026-10-19 framing states (`directives/S12_decision_list.md`, "Filing and framing"; S12
  record M01/M18): each force states the flow it needs, S21 compares them, and a mismatch
  reopens S9b. This is an OPEN S21 deliverable, not a comparison performed by S12 or a
  cross-sector compatibility result. It matches this entry's by-target row.
- **on failure** — The conditional reference-matching statement for the affected observable
  cannot be realized on that domain. The computed conditional functionals and the reported
  cross-engine coverage do not become a GR-matching result.
- **S15 gravity input to the comparison.** If inverse-square signed momentum crossing into the
  shell around drain 2 is sought, each named component must satisfy `𝒟Q + 2Q = 0` along the
  separation interval, including the field-response departure. Domain `d > R > r_in > 0`,
  moving spherical boundaries, adopted projected P5/P6 and the sufficient derivative conditions
  of the [record](steps/S15_two_drain_region.md) apply. Both engines retain live general-field
  integrals for all three net components on the power-law family. Physical field responses,
  nonzero amplitude and attractive direction are open. The eventual Newtonian observable is the
  companion-induced free acceleration. The transfer from shell to core and its inertia remain
  open. This filing supplies the gravity-side conditional to S21. It selects no shared profile.
- **Flow comparison decision, user 2026-10-19** — Each force states the flow it needs; S21
  compares those requirements, and a mismatch reopens S9b. S12 states the source compatibility
  requirements and leaves forms open; it fits no flow and adds no optical calculation or
  profile to this entry. This decision provenance adds no claim of cross-sector compatibility.
- **evidence and limits** — `B_CONDITION_*`, `B_IMPLIED_JN_*`, `C_*_CONDITION_*`,
  `C_*_IMPLIED_JN_*`, and the `R_FLOW_ONLY_`, `R_SPEED_ONLY_`, `R_TILT_ONLY_`,
  `R_BULK_ONLY_` restrictions (S9b §§3–4, measurements M3/M4/M11). Flow-only constrains
  `r(V²−2rVV′)/4`; speed-only `c₀²r²δ′/2`; tilt-only `c₀²rξw′²/4`; bulk-only
  fixed ratio `c₀²(n−1)r²f′/4` and power `c₀²sr²f′/2`, each relative to GM.
  Speed-, tilt- and bulk-only exchange prints zero; flow-only retains the signed live-density
  exchange. Conditions/exchange/γ values have no formed cross-engine comparison; the domains
  and every printed difference retain record §6's limits. This entry asserts no realized match.
- **oracle comparison (evidence only, 2026-10-19)** — the accepted flow-only check (`8d91f448`;
  record `steps/S9b_flow_oracle_comparison.md`) compares the flow-only restriction with Hamilton–Lisle
  (2008) and Visser (1998).
  - Two identifications are needed: the orbital `GM` with each oracle's mass parameter, and S9b's clock
    and speed units with the oracle's. The check supplies neither, and prints both as `UNMATCHED`. Under
    both, the kept flow-only squared speed equals the river model's and Visser §7's. The signed velocity
    matches for one radial orientation per component, and the frame identification that would choose it is
    also `UNMATCHED`.
  - On the nonzero-velocity components, in both engines, the exchange `jn` vanishes only on
    `ρbr = C₁/r^(3/2)`, the density law Visser §7 quotes. On the zero-velocity components the exchange
    vanishes with the density left free.

  The comparison selects no profile and no exchange value, and it delivers nothing of this entry. The
  requirement, its free inputs and its failure consequence are unchanged.


**Counting the two observables under item 5.** The rule is “It counts the deflection
and the radar slope as separate tests only where the printed objects show their conditions differ.”
M17 retrieves that rule and both engines' complete γ objects, including their attached gates,
for general B, R_FLOW_ONLY, R_SPEED_ONLY, R_TILT_ONLY and conditional F_FLOW/F_LIVE.
On each object's on-domain `GM≠0` stratum, PY's γ-difference solution set is `{0}` and WL's
relation contains `gammaDifference==0 && GM!=0`. Thus the two observables are **not counted
as separate tests at these selected grades on that stratum**, in either engine. This is a
count of the printed reference-comparison constraints, not a claim of cross-engine equivalence
or observational redundancy outside this model point. (M17; paired conditions M13.)

At `GM=0`, PY separately prints the full real interval when both observable compatibility
functionals vanish, and `EmptySet` on the other three printed vanishing/nonvanishing combinations.
WL prints `GM==0` with both functionals zero, without a constraint on `gammaDifference`, and
no further compatible alternative. Each engine's individual γ objects retain the corresponding
unconstrained/empty singular strata. The record does not replace these by γ-difference zero:
γ is unidentifiable on the compatible singular stratum, and the incompatible strata supply no
joint γ solution. These strata establish no additional separate γ test under item 5; two
independent tests are not counted from the unrestricted difference. Outside the respective gates,
the objects are `NOT_ESTABLISHED`. (M17.)

Whether the differing WL `InverseWeight`, `ReferenceInverse` and velocity-reduction fields in
M13 state the same condition remains **a separate OPEN question, no owner named**. That question
does not put the above counting decision on hold. Nor are the differing strings treated as a
shown physical difference of conditions: the printed γ strata are part of the counting evidence.
Part C/bulk-only condition rewrites retain M13's field distinctions and earn no additional
separate test here; no unprinted Part C γ object is invented. (M13/M17.)

**Engine-specific restrictions, constants, exchange and domains.** M18 retrieves both
observables' complete `CONDITION` and `IMPLIED_JN` objects for all five restrictions below;
no compact domain omission is used as evidence. In the following readable transcriptions,
`r=s9b_r` for PY, `ρbr=rho_br_profile[r]` for PY and `rhoBr[r]` for WL. Engine-qualified
constants are not identified across engines.

**PY flow-only (`R_FLOW_ONLY_B_*`).** Its printed solution dictionary has
`squared = 2*GM/r`, `integrating_factor = 1/r`,
`general_solution: s9b_velocity_squared(r) = r*(2*GM/r² + C_PY)`, and `residual = 0`,
where `C_PY` names the printed real `s9b_velocity_integration_constant`.
The `constants` field retains `C_PY=0`; its `C_PY≠0` branch is labelled
`EXCLUDED_BY_SUPPLIED_V_OVER_C0_O_EPSILON_HALF`, with the printed limit of `C_PY` retained.
The `nonreal_branch` is `2*GM/r<0`, labelled `EXCLUDED_BY_REAL_MATERIAL_VELOCITY`.
The two positive-component entries print, respectively,

| PY velocity equality | PY implied exchange equality |
|---|---|
| `V = −√2 √GM / √r` | `jn = √2 √GM ρbr′/√r + (3/2)√2 √GM ρbr/r^(3/2)` |
| `V = +√2 √GM / √r` | `jn = −√2 √GM ρbr′/√r − (3/2)√2 √GM ρbr/r^(3/2)` |

Their attached domain keeps the solved ray/applicability gate, `2*GM/r≥0` throughout
`(s9b_b_min,∞)`, and `2*GM/r>0` on each named positive component. The condition itself keeps
its exterior interval and the δ=ξ=0-substituted gate, and prints `GM=r(V²−2rVV′)/4`.
Four further entries retain sign pairs `(-1,-1),(-1,1),(1,-1),(1,1)` at zero radius `z`:
`z>s9b_b_min`, `2*GM/z=0`, matching real one-sided limits of
`σ √2 √GM/(√q (q−z))`, `V(z)=0`, `V′(z)` equal to the right-hand limit, and
`jn(z)=−ρbr(z)` times that right-hand limit. These are the printed differentiable-zero-join
objects, including zero intervals; their limits are not evaluated in this record. The density
function retains PY's positive declaration. (M18, tuple paths `0`, `1/0..5`, `2`.)

**WL flow-only (`R_FLOW_ONLY_B_*`).** Write `I(r)` only as an abbreviation for the printed
`Inactive[Integrate][(-4*GM)/tailRadius^3,{tailRadius,r,Infinity}]`, without evaluating it.
WL prints `GeneralSquaredVelocity = r*C_WL − r*I(r)` with
`GeneralConstantDomain: Element[C_WL,Reals]`, and `SquaredVelocity = −r*I(r)`;
`C_WL` names `velocityIntegrationConstant`. The retained `ConstantDomain` is `C_WL==0`;
real nonzero `C_WL` branches are excluded by the supplied `V/c0=O(r^(-1/2))` counting.
Each positive-component case lists `IntegrationConstants -> {C_WL}`, domain
`C_WL==0 && −r*I(r)>0`, and its own substituted `AlongRayDomain`. It prints

| WL velocity | WL implied exchange |
|---|---|
| `−√(−r I(r))` | `[-ρbr(4GM+5r²I(r))−2r³I(r)ρbr′]/[2r²√(−r I(r))]` |
| `+√(−r I(r))` | `[ρbr(4GM+5r²I(r))+2r³I(r)ρbr′]/[2r²√(−r I(r))]` |

The zero-interval case has `C_WL==0 && −r*I(r)==0`, `V=V′=jn=0` and the same listed
constant/domain. `NoRealVelocityDomain` retains the negative-squared branch. `ZeroSetJoining`
and each `AtComponentBoundary` retain differentiable real velocity/flux and finite matching
one-sided limits. The enclosing condition retains `r>bFar`, `rAnchor>bFar`, functional
regularity/convergence/Abel-invertibility requirements, its `BranchDomain`, the anchor GM
rule and all-radius profile equation; each velocity case carries `AlongRayDomain` with its
own substitution. WL keeps real nonzero density and its optical-ratio gate. No evaluation of
`I(r)`, identification with PY's `2GM/r`, or identification of constants is claimed.
These flow-condition/exchange value and domain differences remain OPEN across engines,
no owner named. (M18; M6/M7.)

**Zero-velocity restrictions.** The following table reports each engine's own condition and
exchange-domain fields for both observables. PY's full attached gate at condition path `3`
is obtained by that restriction's substitutions; its radial interval is `(s9b_b_min,∞)`.
For bulk-only rows the condition is inside `OBJECT_ON_BOTH_PARTS`, with both density-label
conditions and the exact response domain. WL's enclosing `Domain` and `BranchDomain`, and
each case's `AlongRayDomain`, remain attached in every row; its `Domain` includes exterior
radius/anchor and functional requirements. The residual equalities below are WL's additional
`Cases/0/Domain` entries in `IMPLIED_JN` (also inside the velocity reduction of `CONDITION`),
not replacements for those gates. M18 preserves every domain verbatim.

| Restriction | PY printed condition and exchange | WL printed condition and exchange-case domain |
|---|---|---|
| Speed-only, V=ξ=0 | `GM=c₀²r²δ′/2`; `VELOCITY_FIXED_BY_RESTRICTION`, `V=0`, attached speed-only gate, `jn=0` | GM anchor rule plus `2δ′+rδ″=0`; deflection domain `−2GM/(c₀²r)+rδ′==0`, radar domain `−2GM/(c₀³r)+rδ′/c₀==0`; `V=V′=jn=0` |
| Tilt-only, δ=V=0 | `GM=c₀²rξ′²/4`; `VELOCITY_FIXED_BY_RESTRICTION`, `V=0`, attached tilt-only gate, `jn=0` | GM anchor rule plus `ξ′(ξ′+2rξ″)=0`; deflection domain `−2GM/(c₀²r)+ξ′²/2==0`, radar domain `−2GM/(c₀³r)+ξ′²/(2c₀)==0`; `V=V′=jn=0` |
| Bulk-only fixed ratio, V=ξ=0 | `GM=c₀²(n−1)r²f′/4`; fixed-velocity tuple, `V=jn=0`, exact fixed-ratio/substituted-ray domain on both labels | anchor rule plus `(n−1)(2f′+rf″)=0`; deflection domain `−2GM/(c₀²r)+(n−1)rf′/2==0`, radar domain `−2GM/(c₀³r)+(n−1)rf′/(2c₀)==0`, conjoined with the exact fixed-ratio domain on each label; `V=V′=jn=0` |
| Bulk-only power, V=ξ=0 | `GM=c₀²sr²f′/2`; fixed-velocity tuple, `V=jn=0`, exact power/substituted-ray domain on both labels | anchor rule plus `2sf′+rsf″=0`; deflection domain `−2GM/(c₀²r)+rsf′==0`, radar domain `−2GM/(c₀³r)+rsf′/c₀==0`, conjoined with the exact power domain on each label; `V=V′=jn=0` |

In **each** of these four restrictions PY emits no integration-constant field or solved-velocity
dictionary: velocity is prescribed zero. WL explicitly prints `IntegrationConstants -> {}` and
status “prescribed by the labelled substitution”, with no `ConstantDomain` field. Neither engine
selects δ, ξ or f, or prints a solved-profile integration constant for those live profiles here.
The fixed-ratio domain retains real local/asymptotic sound-speed squares, local square ≥0,
asymptotic square >0 and the exact nonnegative-root ratio real and >0. The power domain retains
`c₀ densityPower[1+f,s]` real and >0; `ρ₀>0` and the exact ray substitutions travel with both.
Both keep `ρ>0` and `ρ≤0` labels; no germ-equivalent response replaces an exact domain.

For **each** flow/speed/tilt optical restriction, deflection/radar/residual value samples form
but retain printed numeric differences and unresolved domains; γ rows compare markers only,
and full conditions/exchange have no formed value comparison. For **each** bulk-only response,
only label/domain comparisons form, not condition/exchange values. The WL/PY constant,
squared-velocity, signed-exchange, zero-join and domain representations above are retained as
OPEN differences, no owner named; no equivalence is supplied by the zero exchange entries.
(M18; M5–M7, with the counting rule separately applied in M17.)

### R-S22-01 — how the energy a throat releases splits among light, the extra wave and the bulk

- **source**
  - S11 (`steps/S11_stray_longitudinal.md`, departure table and the paragraphs after it). The table gives the
    question "whether matter directly observes the second cone" the owner "derived matter/interface
    coupling". The record says "A second cone is only a departure if matter COUPLES to it", and describes
    matter as "a throat held open by a trapped brane-shear standing wave".
  - The photon survey, E32–E34 and E38 (`PHOTON_POLARIZATION_SURVEY.md` §4). Each is "In apparent conflict",
    under conditions that differ by entry. In the survey's words, the classification rests on "E32–E34 on
    matter coupling and, for the blackbody comparison, ordinary thermal-radiation participation; E38
    additionally on the extra branch's cosmological relativistic population, energy/perturbation mapping and
    decoupling history". (C6, C7.)
- **target** S22 (the nonlinear throat solve), together with Q2 (the mouth source and its held datum) and
  S12 (the throat's order-conversion source and its boundary data). These owners were chosen by the user on
  2026-10-19; the S11 record itself names only "derived matter/interface coupling". (C8.)
  The user's 2026-10-20 **"Keep S12"** decision retains every owner. S12 stays OPEN for its
  forms, delivered only after the brane is defined following the force sectors; this requirements
  record delivers none of those forms (S12 record §§5–7, retrieval M18).
- **status** OPEN
- **requirement** — Consider a throat held open by its trapped shear wave. The trapped-wave state can change,
  the throat can move, and its mouth can change shape. When any of these happens, the object is the energy
  released, with its frequency content and rate, and how it splits among three channels:
  - the free transverse waves (light);
  - the extra longitudinal wave (S11's second cone);
  - the bulk.

  The entry also requires the reverse: how much of the energy each channel carries in is taken up by a
  throat, at each frequency.
- **on failure** — S11's departure cannot be classified. The record makes the second cone a departure only
  if matter couples to it, and without these objects nothing decides whether matter observes it. E32–E34
  and E38 stay "in apparent conflict", with their condition undecided.
- **note: what a delivered split would and would not decide** —
  - **Each measurement has its own conditions** (survey §4):
    - **E33**, the laboratory blackbody, needs the extra wave to be an additional ordinary thermal mode,
      with the usual occupation, energy and dispersion, and a calibrated emission and absorption response.
    - **E34**, the CMB spectrum, needs its power to enter the calibrated CMB radiance rather than cancel in
      the sky/reference comparison.
    - **E38**, `N_eff`, needs it to be extra gravitating radiation, populated and relativistic in the early
      universe. The survey adds that a weak population or earlier decoupling is not generically excluded.
  - **The split alone decides none of them.** Whether the extra wave is populated depends on several things:
    - its emission and uptake rates, compared with each measurement's timescale;
    - its escape into the bulk;
    - its initial population and the expansion history.

    What a populated wave then carries depends on its own speed. In S11 that speed is set by `B_comp/ρ_br`
    (`R-S6-01`), not by light's `μ_R/ρ_br`. So a share comparable to light's is not by itself a conflict,
    and neither is equilibrium with matter. Nor does a small share by itself rule one out.
  - **The falsifier** is each measurement, applied under its own conditions to the population and power
    that these objects, the extra wave's speed and its escape into the bulk predict. That application is
    later work. This entry states no bound; compare the note in `R-S12-01`.
  - **The quantum part.** The survey's conditions also assume the usual quantum occupation. The plan
    assigns that to a separate project (`V3_STEP_PLAN.md:1198–1204`; survey [P8]), so a classical split
    does not settle that part. (C7.)
- **relation** —
  - S22's linear test 1 asks how much light passing a supplied defect profile converts into the
    longitudinal mode. S11c asked how much light converts where the brane changes (`V3_STEP_PLAN.md:508`).
    Both concern incoming light; neither is a throat's own emission.
  - S22's linear test 2 asks whether the longitudinal mode radiates into the bulk, and needs the interface
    law (`V3_STEP_PLAN.md:1181`). That escape enters the extra wave's population. So applying any
    measurement under this entry also needs test 2's result.
  - The S11 record adds that at homogeneous linear order "the mode's bulk leakage and direct observability
    require the brane–bulk interface law". S11b supplies that law conditionally, for a uniform brane only.
  - The S11 record states that `c_L` is not a gravitational wave, so this entry is not about the carrier of
    gravity. (C6, C8.)

### R-S14a-01 — the dynamically selected slow limit and far behaviour

- **source** S14a part 2 (`steps/S14a_2_slow_limit.md`, §§1–5/7; retrievals M01–M13/M16/M17),
  for the conditional consumers S14, S15 and S20.
- **target** jointly S12 (source and flux forms), S8/S22 (the slab's material law), and
  **no owner named** for the face and interior laws, by the user's V2 decision at `12aa65d5`.
- **status** OPEN
- **requirement** — the slow limit and far behaviour selected by the slab's own dynamics,
  for S14a's still-underived controlled return and profile-dependent `J→Q` obligations.
  The selected steady-throughput and exchange-live decisions keep spreading, exchange and
  storage live. No average, solution family, density gate, order of limits, law, material,
  profile, premise or candidate is selected here.
- **on failure** — S14, S15 and S20 remain **CONDITIONAL, not carried results**.
- **decision and owners** — The user chose **"Record as requirement"** (2026-10-21), then
  **"All the slab's laws (Recommended)"** for V2. The latter is recorded in
  `directives/_measurements/S14a_2_record_r0_review_disposition.md`, "The user's decision on V2",
  at `12aa65d5` (M16). It replaces item 4's S12-only target sentence. The earlier
  **"(the open S12 forms)"** option text and accepted S12's different ownership statements are
  quoted in record §§4/6. No missing owner is invented. S12's source/flux forms and separate
  boundary data remain OPEN under **"Keep S12"**; the complete return-law/map obligations
  are not repartitioned. Recording this requirement delivers none of the forms, the slow
  return law or the source map.
- **evidence and limits** — The step-list pass-1 Q1–Q3 and pass-2 R1–R3 dispositions and the
  separately scoped Codex/Grok instruments in `directives/_measurements/S14a_2_review_witnesses/`
  are retrieved in M03–M10. Each is a review-leg instrument, not an accepted build; their
  constructions/checks differ as record §§2–3 states. The material-face zeros made from
  identical operands or typed zeros are displays, not checks. The pass-2 adjudication's
  dynamics reading remains a reading. No build ran by the user's decision. The inherited
  real-fraction/A13, live-face and part-1 domain/comparison qualifications remain (M11–M13).
- **entry rather than merge** — The register keys an entry by its deliverer and merges only
  the same requirement with the same deliverers. None of the 18 entries at `978c359c` covers
  this selection object with V2's joint deliverers. `R-S12-02` covers S12's native forms/data
  and compatibility objects, not the jointly owed selection; `R-S22-01` concerns the throat
  energy split, and `R-O2-01` concerns a net-power supplier/budget. Their objects and targets
  remain as recorded. This is a separate entry under the schema and method, not a second
  obligation added to `R-S12-02`'s target (record §5; M12/M17).

### R-S16-01 — the surroundings object needed for a passive factor

- **source** user's selected 2026-10-21 "Record what it needs" option: "What one needs
  (a named field of the surroundings, a rule separating it from the body's own field, and its
  normalization) goes to the register as a requirement" (`directives/S16_decision_list.md:49–53`;
  pass-1 disposition, "The user's decisions"). `steps/S16_moving_region.md`, §§1/6
  (S16 M01–M02/M04/M07/M14/M17).
- **target** no owner named
- **status** OPEN
- **requirement** — a named field of the surroundings, its separation from the body's own field,
  and its normalization, for the open passive-factor/gradient-form interpretation.
- **on failure** — received momentum cannot be promoted to that passive-factor interpretation.
  The supplied regional balances and received-momentum object remain conditional objects.
- **standing** — explicit future requirement selected by the user; no passive factor or field is
  assumed or derived. G7 is an oracle under its own conditions, not this object's delivery (S16 §5).
  This is distinct from source normalization in S22's debt account and from S21's comparison of
  force-sector flow requirements. No owner is inferred from those neighbouring interfaces.
- **admission ground** — the selected option expressly requires this future object for a passive
  interpretation; this is not an unselected closure form inferred from a later owner.

### R-S16-02 — the body's compactness in both comparisons

- **source** orchestrator's S16 list filing `directives/S16_decision_list.md:158–161`,
  "The body's compactness, separately", from pass-2 E3 and pass-3 F3; filing/framing `:24–30`;
  source assumption `research/4d_1pn_full/paper/4d_1pn_full.tex:697–700`;
  `steps/S16_moving_region.md`, §§4/6 (S16 M01/M04/M14).
- **target** no owner named for physical delivery; S22 is the named recipient of this debt in its table
- **status** OPEN
- **requirement** — the body's size small relative to inter-body separation and to the scale of
  variation of a field of the surroundings, kept as two comparisons.
- **on failure** — the compact-body response approximation cannot be claimed. S16's conditional
  region balances and count moments do not determine either comparison.
- **standing** — the original filing is the list author's separate small-body requirement from
  E3/F3, not an imposed compactness value or a derived property of the regional fields.
- **admission ground** — user's decision of **2026-10-24**, **"Require all four (Recommended)"**:
  **"You require all four on the register now, and each entry cites today's decision as its reason."**
  `directives/_measurements/S16_record_r1_review_disposition.md`, "The user's decision" (S16 M19).
  This explicit user requirement admits the entry; record-list item 7 is not its admission ground.
  The decision adds no owner and selects no law, form or material. The list-author source above
  and this entry's target remain unchanged.
  S22's debt-table recipient is not evidence of a physical compactness deliverer. The C10/tail
  question remains S19's. Count moments, control radii, profile width and mouth radius are not identified.

### R-S16-03 — momentum exchange at each slab face

- **source** orchestrator's S16 list filing `directives/S16_decision_list.md:162–163`,
  "For the faces", from pass-2 E2's direction to keep "the face-resolved momentum accounting
  OPEN"; supplied momentum input `:84–85`. `steps/S16_moving_region.md`, §§1–2/6
  (S16 M01–M02/M04/M07/M14/M17).
- **target** no owner named
- **status** OPEN
- **requirement** — the momentum exchanged at each slab face and its split by the motion of
  that face, separately identifiable within the projected total exchange.
- **on failure** — S16's undivided force-density input supplies no face-resolved momentum account.
  Its moving control-surface pieces do not deliver a slab-face split.
- **standing** — the list author's original filing retains the missing face-resolved accounting
  object identified by E2; the undivided supplied balance is not a defect or a face-resolution law.
- **admission ground** — user's decision of **2026-10-24**, **"Require all four (Recommended)"**:
  **"You require all four on the register now, and each entry cites today's decision as its reason."**
  `directives/_measurements/S16_record_r1_review_disposition.md`, "The user's decision" (S16 M19).
  This explicit user requirement admits the entry; record-list item 7 is not its admission ground.
  The decision adds no owner and selects no law, form or material. The list-author source above
  and this entry's target remain unchanged.
  The faces' law remains unowned. `R-S12-02` retains S12's source/flux/data and material-face scope;
  it does not assign this momentum object or the faces' law to S12.

### R-S16-04 — the momentum–material relation needed for inertia

- **source** user's selected 2026-10-21 "Count material; mass open" option: "The law it needs
  goes to the register, as O2 keeps it (S8)" (`directives/S16_decision_list.md:54–56`;
  pass-1 disposition, "The user's decisions"); list filing `:164–166`;
  `directives/O2_SHARED_PHYSICS.md:133`; `steps/S16_moving_region.md`, §§1/6
  (S16 M01/M04/M14).
- **target** S8, with S1.5 antecedents when available as the O2 row states
- **status** OPEN
- **requirement** — the relation of projected momentum density to the material fields and motion
  needed to define an inertial factor.
- **on failure** — the separately printed momentum and material-count centre relations do not
  supply an inertial factor or a relation between momentum and count-centre acceleration.
- **admission ground** — the selected option expressly requires the missing momentum–material
  relation for that future inertial interpretation; no kinetic form is selected (S16 M01–M02/M17).
- **entry rather than merge** — `R-S8-06` requires material displacement identification and the
  original homogeneous quadratic/slab inertia. Its O2 scope expressly excludes a live momentum-map
  widening. The present explicitly user-required live object has a different domain and requirement;
  no homogeneous kinetic form or constant inertia is prescribed here.

### R-S16-05 — the slab density response

- **source** orchestrator's S16 list filing `directives/S16_decision_list.md:164–166`:
  "The slab's density response is also OPEN there", citing
  `directives/O2_SHARED_PHYSICS.md:138`; `steps/S16_moving_region.md`, §§1/6
  (S16 M01/M04/M14).
- **target** S8 (substrate reduction) and S22 (completion), as the O2 row states
- **status** OPEN
- **requirement** — the slab density response with live bulk, flow, embedding,
  thickness/projection dependence and gradients, needed for the open material/mass accounting.
- **on failure** — that density-response interpretation remains unavailable. S16's supplied
  projected count balance and its live-density derivatives still supply no such response.
- **standing** — the original source and filer are the list author's O2 `:138` density-response
  handoff. Live projected count and its derivatives do not deliver that response. O2's earlier
  "No new density-response entry" remains its historical application of the rest-on test;
  the S16 admission below is a later explicit user requirement, not a new dependence of O2's object.
- **admission ground** — user's decision of **2026-10-24**, **"Require all four (Recommended)"**:
  **"You require all four on the register now, and each entry cites today's decision as its reason."**
  `directives/_measurements/S16_record_r1_review_disposition.md`, "The user's decision" (S16 M19).
  This explicit user requirement admits the entry; record-list item 7 is not its admission ground.
  The decision adds no owner and selects no law, form or material. The list-author source above
  and this entry's target remain unchanged. No bulk equation of state, factorization or response
  form is selected here (S16 M01/M04/M14/M19).

### R-S22-02 — the missing response-reduction inputs in S22's debt table

- **source** orchestrator's S16 filing `directives/S16_decision_list.md:28–30`: "What stays a
  debt for S22's table is what this step does not deliver";
  `V3_STEP_PLAN.md:732–741,1303–1305`; `steps/S16_moving_region.md`, §§4/6
  (S16 M01/M04/M14).
- **target** S22 for the debt-table account; no physical supplier is assigned here where the
  sources name none
- **status** OPEN
- **requirement** — the outstanding reduction-data inventory: `M`, `a` (`a_WT` in the plan),
  the multipoles of `M`, source normalization and any suppression of boundary fluxes, with the
  `a_WT`/`a_mouth` bridge separately identified under A11. Compactness is filed once as
  `R-S16-02` and cross-referenced in this table account.
- **admission ground** — user's decision of **2026-10-24**, **"Require all four (Recommended)"**:
  **"You require all four on the register now, and each entry cites today's decision as its reason."**
  `directives/_measurements/S16_record_r1_review_disposition.md`, "The user's decision" (S16 M19).
  This explicit user requirement admits the entry; record-list item 7 is not its admission ground.
  The decision adds no owner and selects no law, form or material. The list-author source above
  and this entry's target remain unchanged; no physical supplier or new regional-balance premise
  is inferred from the decision.
- **on failure** — the conditional response-side monopole reduction cannot be treated as
  delivered or interior independent. Count, count moments and live boundary fluxes are printed
  objects; none is substituted for a missing mass object or for flux suppression.
- **scope** — this target files the debt account; it neither assigns S22 every physical supplier
  nor directs S22 to suppress a flux. The core interior/count/momentum remains an S22 OPEN handoff.
  C3's drain–mass gap remains a defect-register gap, not a delivered mass relation. No law,
  material, premise or candidate for meeting any listed debt is selected.

## Population passes

**Pass 1 — S9 and S10. ⭐ DONE 2026-08-06**, alongside S10's step record and from the same reading.
Seven entries added, against the two that existed. The calibration it produced, for pass 2:

- ⭐ **Consolidate rather than duplicate.** S9 states its shear requirement twice — once for the ordered
  phase and once for the disordered — and a third consumer (the geon's trapped mode) appears elsewhere in
  the record. It is **one object**, `R-S1-02`, with three consumers, ⛔ not three entries.
- ⚠ **Distinguish a requirement from an ansatz by asking whether a retirement condition is LIVE.** `μ_R`'s
  **magnitude** is postulated and defect `B2` closed the polar-`P` route to deriving it ⇒ ansatz, ⛔ no
  entry. Its **form** has no such closure ⇒ `R-S8-01`, a live requirement. The two sit in the same
  sentence of S9's record and go to different registers.
- ⛔ **A requirement is not everything a later step would find useful.** A simulation needs `ρ_br` and
  `μ_R` separately rather than only their ratio; that is a *consumer's* need pointing **forward**, ⛔ not
  an obligation on which a banked result rests. It got no entry.
- ⭐ **The load-bearing ones are the identifications nobody computed.** `R-S8-02` — that the in-plane and
  out-of-plane sectors decouple — is the single entry on which S10's headline number depends, and it was
  flagged inside S10 **at the moment the identification was made**. ⇒ read a record for its *flagged
  identifications* first; they are where the requirements are.


### Pass 2 — S11, S11b-A/B and S11c PARTIAL · populated 2026-10-05; accepted `22be0137`

Read `steps/S11_stray_longitudinal.md`, all three S11b records (the unified record and the historical A/B
records), `steps/S11c_PARTIAL_CLOSEOUT.md`, and the a/b/c1/c2/d records in full. Applied the schema,
rest-on test and both routes below: **11 → 15 entries**, four new and seven existing entries gaining
sources (`R-S7-01`, `R-S1-01`, `R-S1.5-01`, `R-S6-02`, `R-S8-01`, `R-S8-02`, `R-S8-03`). New entries:
S1 `R-S1-03`; S8 `R-S8-06`; S12 `R-S12-01`/`02`.
All use schema status OPEN; CONDITIONAL is a qualification of the source result, not a status value.

S11c contributes only conditions of its kept uniform/selected-uniform results, with the closeout and d
record's review limits retained. The clean-condition packet sources no entry: the closeout marks only
uniform decoupling CONDITIONAL (`steps/S11c_PARTIAL_CLOSEOUT.md:7`), while its no-leak reading is OPEN
(`V3_STEP_PLAN.md:532–534`; `steps/S11c_d_profile_conditioned_scattering.md:64`). Its UNRESOLVED direct
mixed term, numerical loss, and OPEN nonuniform confinement/material, observable and real-throat questions
source no entry.
**Supplied objects.** The unified spec's names are the §1b branch prescription and the §4 face closure.
It heads §1b *"COMPLEX FREQUENCY — SUPPLIED"* (`directives/S11b_SHARED_PHYSICS.md:113`) and §4 *"THE FACE
CLOSURE — SUPPLIED, and the affinity is the canonical driver"* (`directives/S11b_SHARED_PHYSICS.md:208`).
Its supplied column reads *"§2 bulk acoustics; the §1b branch prescription; the §4 closure, affinity `𝒜`
and its normalization; the §5 term list (⛔ but its **closedness** is tested); the §6 balance-law route and
virtual-displacement rule"* (`directives/S11b_SHARED_PHYSICS.md:995`). B's record lists five supplied
objects: *"The affinity, the branch prescription, the mass balance, the virtual-displacement rule and the
row structure are **supplied**"* (`steps/S11bB_interface_assembly.md:136–137`). It then says *"⭐ Their
verification comes from the independent reviewer derivations above, ⛔ **not from the engines
agreeing.**"* (`steps/S11bB_interface_assembly.md:138–139`). What the reviewers re-derived is a different
list: *"the affinity `𝒜 = μ_s − δp/ρ_m` (**three times**), the complex-frequency continuation (**three
times**), the interfacial mass balance, the `Z_perm` reduction, the virtual-displacement rule, the
ten-element basis, the dimension table, and `K₀`"* (`steps/S11bB_interface_assembly.md:126–129`).
Separately, B says *"**Radiated-energy direction is inherited** from the supplied continuation"*
(`steps/S11bB_interface_assembly.md:158–159`). This register draws no postulate inference from
"supplied", and leaves open whether any supplied object belongs in the ansatz ledger; v3 has no ansatz
ledger yet.
Reading a/b/c1/c2 did not promote their pre-repair closure language or deferred comparisons to accepted
nonuniform physics. The clean-condition packet does not establish a support law; the
exploratory throat/EM documents supply candidates only. No new requirement is inferred just because
a future scattering calculation would need an input. S11's deferred interface/spectral and nonlinear
questions are not recast as dependencies of its closed homogeneous census.

**Bulk EOS exponent (`n`, rule 6).** No kept S9, S11, S11b or S11c result read here rests on selecting
`n`: S9's supplied bulk-sound formula carries the inheritance, while its light result uses the brane
stiffness/inertia ratio; S11 and S11b use `c_s0` in the bulk acoustic branch, and the kept S11c uniform
checks fix or vary effective sound speed without deriving the EOS. **No S2 entry is added.**
`n_eos = 5` is inherited. Its original light-side justification took light's phase speed to be the bulk
sound speed (`research/1pn_optics/paper/1pn_optics.tex:720–737`, `:975–978`); that paper's competing
constructions and reasons for preferring `n = 5` are at `:2165–2343`. Whether, and where, `n` enters
light bending and delay under v3's light picture remains open — **no owner named** in the records read.
`V3_STEP_PLAN.md:208–250` retains the provisional medium EOS and its half-two parent-action debt;
**O-02 stays open**. This pass does not decide whether `K` and the exponent count as one entry or two.

**Second route, by step.** S11 has **no identified prior-art route**: its record names the
rejected-Cauchy-branch coefficient but no objection to that branch, so the Cauchy branch is a candidate
only. The FORM outcome stays solely in `R-S8-01` as a first-route source, quoted from the S11 record.
S11b-A, S11b-B/unified S11b, and S11c-a/b/c1/c2/d each have **no identified prior-art route** meeting the
directive's record-and-objection test. The S11c closeout, which is a record, does name MacCullagh, in an
item it marks OPEN: *"**OPEN — material audit/Q2/S22:** uniform confinement is derived within the
supplied linear model; nonuniform confinement and the physical admissibility of rotational stiffness are
not established. That is the useful but limited MacCullagh distinction."*
(`steps/S11c_PARTIAL_CLOSEOUT.md:27`). That item is listed under *"Claims future work must not inherit"*
(`steps/S11c_PARTIAL_CLOSEOUT.md:23`), and an OPEN item sources no entry under the rest-on test, so it
sources none here. Onsager–Casimir is a candidate only; `R-S1-03` instead rests on B's
first-route hand-forward, *"conditional on microscopic time-reversibility"*. A's memoryless reduction is
by construction, and internal cross-engine reproductions are not prior art. The two pre-existing
prior-art entries are not re-certified by this pass. The obsolete rebuild schedule is replaced by this
completed population pass. Record qualifications/disagreements are reported at STOP, not adjudicated
here.

**Directive/register disagreement on the second route — unresolved.** The second-route results above
apply the directive's test: *"Run it only on a prior result that a record in the reading list identifies
as reproduced by the sector, with the objection taken from that record or from a source the record
cites."* (`directives/S11_requirements_pass2_directive.md:56–57`). This register's own second route, under
"Method for population passes" below, found `R-S8-04` and `R-S8-05` outside the records — *"⛔ Neither is
anywhere in S9's or S10's records, because the records only capture what their authors thought to
doubt"* — and says *"⇒ ⭐ **Run this route for every sector with identified prior art**, ⛔ not only the
light sector."* The directive also says *"Where this directive and the register disagree, the register
wins, and you report the disagreement."* (`directives/S11_requirements_pass2_directive.md:21–22`). Both
sides are quoted here; this pass does not resolve the disagreement.

#### Inputs a future nonuniform calculation must define

These are input assignments used by the S11c calculations, not physical values earned by a kept
result and not requirements on a substrate step. Fixed-input uniform checks retain their stated
domain; they do not select those inputs for a real defect. No values are copied here. The table has
**11 grouped rows**; every owner below is from the closeout, with **no owner named** where it assigns
none to that input. The structural objects above are distinct from these numerical assignments.

Input-map provenance **P**: `steps/S11c_d_profile_conditioned_scattering.md` →
`_measurements/S11c_d_sympy_builder_report.md`, “Retained user-approved solver/export contract,” items 2–3 →
`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_d_channel_preflight_input.json`
(`unit_frame`, `profiles`, `parameters`). P records the original development case; the d record's
“omega=3 benchmark” and “uniform check” sections record later frequency/speed choices. None is a
calibration of the draining model.

| Input family (values not reproduced) | Where recorded | Owner named by the closeout |
|---|---|---|
| Supplied thickness/modulus profiles `w₁`, `m₁`, their end data, contrast `eta_bg` and width scale `L_W` (with the physical `sigma_W` binding) | P, `profiles` / `parameters`; builder report, retained contract 2–3; d, benchmark | Q2/Q3/S22 for a physical holder/charge profile |
| Reference thickness `W_0` | P, `parameters`; a, background ansatz | no owner named |
| Density assignments `rho_br`, `rho_m` and the selected density/anchoring case | P, `parameters`; d, uniform check / benchmark | no owner named |
| Physical transverse stiffness, represented in the input map by `mu_R`, `mu_S` | P, `parameters`; unified S11b, transverse mode / representative split | no owner named for these absolute assignments |
| Thickness inertia `mu_W` | P, `parameters` | no owner named |
| Scalar stored-energy inputs `B_rho_3`, `C`, `k_W` | P, `parameters`; S11b-B, breathing stability | no owner named |
| Gradient/mixed stored-energy inputs `kappa_W`, `kappa_theta`, `kappa_theta_W`, `G_W_u`, `G_theta_u` and any remaining constitutive bindings | P, `parameters`; b, energy basis; builder report, retained complete-parameter-map contract | no owner named |
| Face-response amplitudes `Lambda_A_0`, `Lambda_V_0`, `Lambda_X_0` and memory times `tau_A`, `tau_V`, `tau_X` | P, `parameters`; S11b-A/B, interface response | no owner named for these numerical assignments |
| Effective bulk sound speed / actual transverse phase-speed ratio, distinguished from the bare S9 ratio | P, `c_s0` and brane coefficients; d, uniform check / benchmark; `_measurements/S11c_d_near_unity_uniform_continue_result.md`, coverage | S20a (calibration); S22/R10 (derivation) |
| Unit frame, incident frequency `omega` and the two tangential momenta | P, `unit_frame` / `parameters`; d, uniform check / benchmark | no owner named |
| Numerical regulator, quadrature/basis/cutoff and source/domain intervals | d, benchmark and controls; `_measurements/S11c_d_numerical_radiating_balance_report.md`, controls and their limits | no owner named |

### O2 pass — conditional live steady balance · populated 2026-10-07; accepted with the O2 record, `72866fcf`

Source: `steps/O2_steady_brane_balance.md`, §§1–9; governing inputs and accepted constructions/comparator
are identified there. Commands and literal retrieval output: `steps/_measurements/O2_record_measurements.md`,
M1, M4, M6–M12; generator `scripts/O2_record_measurements.py`. The pass follows this register's schema,
rest-on test and both routes. **15 → 16 entries**, all OPEN. **New entries by target: no owner named,
`R-O2-01`.** Six existing
entries gain the O2 source: `R-S1-02`, `R-S8-04`, `R-S8-05`, `R-S8-06`, `R-S12-01`, `R-S12-02`.
Original target inferences and review qualifications are retained; no DELIVERED or RETIRED status is assigned.
Repair 4 against preserved r3 `7ce97e8a` is **16 → 16**: no new entries by target and no newly sourced
existing entry in this repair. The O2 source on `R-S8-06` remains limited to material identification.

**First route / one rest-on criterion.** Enter an O2 source only for a sourced condition on which the
presently reported conditional object rests: its stated material/phase identifications and its
explicitly carried accounting/admissibility obligations. An unspecified response form needed only
to select or solve a later closed material state is an OPEN handoff, not an additional delivery
condition of this unsolved named object. A future closure owner alone does not satisfy this test.
This applies alike to every S8 deferral in the record; the table below states each application.

The material identity, adopted shear-bearing/shear-free setting, expressly carried angular-momentum/
frame questions, conditional non-passive-interface condition and distinct native source/boundary
accounting merge into the six existing objects. `R-S8-06` gains only O2's supplied material
identification; its original quadratic/slab inertia domain is retained, and the r1 flowing-inertia
widening is withdrawn. `R-S12-01` carries only the conditional non-passive
interface reservoir obligation, with no live owner named. The distinct O2 net-supplier/budget
requirement is `R-O2-01`, with no owner named for its delivery; spec §6 expressly requires that object
to accompany any relation requiring net supply. It is entered on that rest-on condition, not simply
because a future solve might use it. Relaxation power remains separately explicit and its reference
owner unassigned. None of these sources proves the substrate can deliver the object.

| S8-related dependence named by the O2 record | Application of the same criterion / O2 entry effect |
|---|---|
| Material identity: `u`, `ρ_br`, `V` (contract §2; spec §3.1) | The stated conditional object assumes this identification. O2 sources `R-S8-06` for this identification only, with no new kinetic form. |
| Angular-momentum/couple carrier and physical rotational frame, `𝒜_rot^live` (contract §3; spec §§3.2, 6) | The sources expressly carry the existing S8 admissibility/frame obligations into the material and paired power accounting. O2 sources `R-S8-04`/`05` for those questions only; no new rotational-response form or historical-objection extension. |
| Live inertia `ℐ_br^live` and momentum map (spec §3.2) | Unselected closure forms; no O2 inertial widening of `R-S8-06`. The material identification above supplies neither `ρ_br V` nor a live kinetic law. |
| O1 `ℳ_⊥` stiffness response (spec §3.2) | Unselected closure form, just as live inertia. No O2 addition to `R-S8-01`/`03`: the supplied optical ratio and premise 1 do not select a general response functional or establish its physical sign. |
| O7 `ℛ_br` density response (spec §3.2) | Unselected closure form, just as live inertia and O1. No new density-response entry; retaining/differentiating live `ρ_br` in the coordinate mass input uses no density-response law. |
| `𝒫_br^cons`, `𝒯_br^cons`; energy-reference/improvement convention including `C_ref` (spec §§3.2, 8.6) | OPEN antecedents/conventions, not selected live response forms used by the reported result. No O2 addition for their S1.5/S8 completion. |
| Full live stress `𝒯_br^live` and normal response `𝒩_br^live` (spec §3.2) | Their forms remain OPEN; premise 1 supplies character only. No O2 response-form entry or widening of the original quadratic/stiffness entries. S5–S8 ingredients and S8/S22 completion remain handoffs. |
| O4 `ℰ_h^live` and its S8/S22 completion (spec §3.2) | OPEN live equation/identification, not an assumed equality with O2's normal content. No O2 equation or coupled-operator requirement is added. |
| Material energy storage/transport `ℰ_br^live`, `𝒥_E^live` (spec §6) | OPEN forms. O2 carries their force/power compatibility and accounting conditions without selecting a stored-energy/transport law or entering an S8 energy-form requirement. |
| Reference/strain evolution `ℛ_ref/strain^live`, with inferred S8/S22 links (spec §3.2) | OPEN form under adopted premise 1; owner remains unassigned. No new S8 response-law entry. Explicit relaxation power and conditional net-supplier accounting remain separate obligations, not delivery of this form. |

The four user-selected choices have the exact provenance label **adopted substrate input to a
conditional model (2026-10-06)** (`directives/O2_premise_decision_list.md`, premises 1–4). They are
postulates for this conditional result, not newly derived requirements. The existing live retirement
condition for the two-phase shear response is retained in `R-S1-02`; the new premise choices supply
no additional named retirement conditions requiring duplicate entries. Their provenance is kept here
for S21; this pass edits neither an ansatz ledger nor the defect register.

| Adopted premise | Content retained for S21 | Unresolved form / owner boundary |
|---|---|---|
| 1 — material reference | Optical shear response elastic; relaxation under steady load; explicit relaxation power, sign OPEN, net supply with named supplier/budget | `ℛ_ref/strain^live` and full live response OPEN; relaxation/reference owner unassigned. Optical consequences deferred; S9's no-dissipation/frequency-independent limits are revisited. |
| 2 — drive | Dynamical order-conversion drain; no separate external body force | O2 chooses the absent separate-body-force representation. Source/return and partners S12; S14a/S14 source bridge and S16 response matching remain OPEN, with S16 assignment an inference. |
| 3 — exchange | Converted material carries local brane material velocity `V`; closes transported momentum only | In-plane carried `j_n V^i` in the outward-loss convention; bulk-direction reduction and additional reaction/energy partners remain OPEN through `𝒥_map`, normal response and S12. Momentum density remains OPEN. |
| 4 — bulk | Postulated shear-free scalar bulk; native face-normal mechanical load; no independent tangential bulk stress | Signed normal-load amplitude, geometry/projection and complete support partition OPEN; tilted loading can project in-plane. External supports and exchange remain distinct. |

O1/O3–O7, A13, live stress/energy/reference forms, missing grades/scales, core selection, source/response
`GM` matching and S21's integration/sort remain the record's OPEN handoff. They are not entered merely
because a later closed profile calculation would use them. The banked claim is not a solution that
already rests on those deliveries. In particular, the emitted mass residual is unsolved, no induced-
measure mass law is derived, and no steady equilibrium, power sign, passive response, optical outcome
or physical holder is earned. Existing rest-reference `R-S6-02` is not extended to a flowing equilibrium;
quadratic form/sign, mode-count, compression, width and EOS entries gain no new O2 source on that basis.
All cross-engine OPEN-content differences remain open; they are not recorded as substrate defects.
The pass inherits the record's whole difference ledger and its comparison bounds (§§3–5), carried
with each applicable accounting/admissibility condition in its Part D handoff (§8). In particular,
the record's literal counts agree only at the joint stored `(role, orientation, open_free)`
classification tuple count level (M7). The comparator prints empty head/role/orientation differences
at the unioned balance-role inventory level (M6), and zeros at the balance closed-part level (M8).
Its single `hold_normal`/`OPEN_free` comparison does not compare the three stored occurrences
per engine separately. These bounds travel with every condition below:
the power-accounting conditions retain the eight one-sided declared-role actions, the four paired
`energy_balance` action-orientation differences `[["argument",1]]` (stream line 207), and the stored
PY density/flux nesting in `OPEN_JointPowerAccounting` in `energy_power` and `energy_balance`
(lines 202/206; M6/M12). None settles compatibility or judges either engine to duplicate power.
The comparator acceptance assigns the WL-only `xi_w''` interpretation to the O2 record, sub-step 7;
the directive allows retrieval only. It remains an **undischarged sub-step-7 obligation returned
to the orchestrator under M1**, as quoted and bounded in record §§4/9 (M1). It is not unowned,
not computed here, and not a new substrate-delivery entry. The handoff carries that status alongside
every admissible-dependence and force/power condition it bears on.

**Second route.** Read the O2 record against scoping §5 and its source qualifications. That section
identifies surface-mechanics, diffuse-interface/evolving-reference and river/acoustic-flow **candidate
oracles**, with explicit dimension, model and access limits. O2 reproduces no identified prior-art
result; no oracle was applied by this record. There is therefore no O2 reproduced result whose historic
objection can generate a new entry under this register's second route. The earlier MacCullagh
admissibility/frame entries remain OPEN and gain only O2's explicit carried S8 obligations; they are
not re-certified. No new search, oracle check or prior-art model was introduced.

**Directive/register relation.** No operative disagreement is found for this O2 pass: directive item 5
says to preserve the distinction between a derived requirement and an adopted substrate premise;
this register says “a **postulate** is a value or form assumed *here* and not derived” and “A postulate
with a named retirement condition generates a requirement.” O2's premise provenance is retained as
conditional adopted input; the same rest-on criterion above excludes every unselected S8 closure form,
including live inertia, O1 and O7, while retaining the sourced identifications and explicit carried
obligations. No closure-form delivery is inferred merely from its future owner.
The earlier pass-2 second-route disagreement, quoted above, remains unresolved and is not adjudicated
or erased by this pass.


### S9b pass — conditional light bending and radar slope · populated 2026-10-19; accepted with the S9b record, `ff3bac63`

Source: `steps/S9b_light_bending.md`, Parts A–C, §§1–8 and appendices. Retrieval commands and
literal output: `steps/_measurements/S9b_record_measurements.md`, M1–M22; generator
`scripts/S9b_record_measurements.py`. **16 (`caa0710d`, before S9b) → 17 entries**, all OPEN;
round 2 kept 17 → 17 from `d1511b20`; round 3 kept 17 → 17 from `4760e6e2`; round 4 keeps 17 → 17 from `d877e857`. **New entry by target:
S21, `R-S21-01`.** Existing entries gaining a source over `caa0710d`: **`R-S1-01`, `R-S1-02`, `R-S1.5-01`,
`R-S8-01`, `R-S8-03`, `R-S8-06`**. The first three and `R-S8-03` gained their source in
round 2; `R-S8-01/06` already carried it in `d1511b20`. Rounds 3–4 add no source gains. M10 retrieves
the pinned censuses and source lists.
No status is changed to DELIVERED or RETIRED. This pass was accepted with the S9b record at `ff3bac63`.

**First route / rest-on test.** The reported optical functionals rest on the supplied local
transverse stiffness and on the identification of displacement, density and advected velocity
with the same material. Those consumers remain in R-S8-01/06. The kept positive optical-ratio
gate additionally rests on R-S8-03; coordinate d³x/spherical flux on R-S1-01; the material's
optical shear on the brane half of R-S1-02; and the supplied Part C acoustic branch on the
acoustic part of R-S1.5-01. Failure of each breaks the scoped interpretation stated in its
entry. These are six existing requirements with S9b consumers, without widening them into
unselected full material laws. The separately reported conditional reference-matching statement rests on the
printed Part B profile condition for each observable and Part C rewrites, with their printed
exchange/domain restrictions. That condition is an obligation for S21's integration/knit and
passes the test; it does not establish that the substrate satisfies it. `R-S21-01` specifies
only the required objects and failure consequence, no means of meeting them, no balance or
candidate material, and no effect of the adopted density link. The deflection/radar
conditions retain M13's explicit distinction, detailed in R-S21-01: PY paired payloads are
byte-identical; WL GM/profile fields match but inverse and velocity-reduction fields differ.
Whether the differing WL fields state the same condition remains a separate OPEN question
(M13). Item 5, applied to each engine's γ-difference zero on GM≠0, does not count the two
observables as separate tests there; the printed GM=0 unconstrained/empty strata earn no
additional separate γ test. This is not a count awaiting WL field-equivalence work (M17).
The source's full difference and unresolved-coverage ledger travels with the requirement.

| Candidate/source dependence | Application of the same rest-on test |
|---|---|
| Supplied transverse stiffness input, material identity, positive optical ratio, dimension, brane shear and Part C acoustic input | Selected inputs: R-S8-01/06, R-S8-03, R-S1-01, R-S1-02 and R-S1.5-01, with the scoped failure consequences above. No live inertia/full-stress or bulk-confinement widening. |
| Part B/Part C matching conditions and single-mechanism restrictions | Derived conditional requirements relative to independent GM: R-S21-01, with complete values still uncompared across engines. No realized matching outcome. |
| No-far-zone-loss forward premise | Adopted conditional input (user 2026-10-06/07/08), re-tested at knit; not a derived general exchange requirement. Its Φ–GM relation remains open. No separate entry. |
| w-parity, P5, P6, density link | Adopted conditional inputs from D2, with dates/scope retained in record §7 and re-tested at knit. Scoped to held Part D, used by no A–C computed object. No derived requirement or effect of the link inferred. |
| O2 P1–P4 | Adopted substrate inputs to a conditional model (2026-10-06). S9b records them without deriving a mechanical/relaxation/loading law; its optical dispersion is separately supplied. O2's existing entries/status are not re-certified. |
| S8/O1/O7 constitutive responses; S11c nonuniform operator; O3–O6 maps, supports, reactions | Unselected closures or excluded work. Owners remain in record §7; no new closure-form entry merely because a later physical realization needs them. |
| Positive-real optical ratio/constructed ray; restricted response domains | Travel with R-S21-01 and add the scoped sign consumer to R-S8-03. No substrate sign result is derived; PY/WL brane-density assumptions and all nine same-named branch-row differences remain OPEN (M19). |
| Historical limits, normal-drain correction, direction/polarization/coupled/subprincipal exclusions | Remain explicitly scoped handoffs. S12 drain/return ownership and the existing R-S12-02 qualification persist; S9b does not deliver that correction or add a flowing-equilibrium source to R-S6-02. |
| O2 derivative-key interpretation and other engine differences | Undischarged interpretation owned by the orchestrator; an unresolved difference alone establishes no substrate delivery requirement. No entry. |

The source dates and complete adopted contents (including live α, ccomp, η and ζ in their held
Part-D scopes) remain in S9b §7 for re-testing at the knit. Adopted values/forms are postulates,
not derived requirements. No ansatz or defect register is changed. The source's statement that
premises are re-tested at the knit supplies no new A–C constitutive law; held Part-D choices
have no present A–C rest-on consumer to duplicate as a new entry.


**Rest-on census, every entry.** This applies the same failure test to the kept A–C object,
not to a future completed substrate or a proposed Part D solve. M1/M18/M19 supply the kept
inputs and gates; M10 retrieves the schema, entries, source lists and this pass.

| Entry | S9b consumer / disposition |
|---|---|
| R-S6-01 | No `B_comp` input or compression-mode result is used; no source gain. |
| R-S7-01 | ℓ is a supplied reduction scale, not a wall-width solution; no nonuniform-width flat-direction result is consumed; no source gain. |
| R-S1-01 (target S6) | Coordinate d³x, spherical geometry and `4πr²` flux consume D_brane=3; source added. |
| R-S1-02 | Propagating brane optical shear is consumed; source added on that half only, without a bulk-confinement claim. |
| R-S1.5-01 | Part C's exact fixed-ratio light/sound response consumes the supplied bulk acoustic branch; source added on that branch only. |
| R-S6-02 | The kept local LAB_HELD construction does not assert a substrate equilibrium/first variation at the unstrained rest state. Its historical uniform anchor is a supplied input, not a newly banked equilibrium result; no source gain. |
| R-S8-01 | Supplied local transverse stiffness form; existing S9b source retained. |
| R-S8-02 | No coupled fluctuating `(u,h)` spectrum, block diagonality or h/transverse degeneracy is claimed. ξ is a supplied background embedding here; the full coupled operator is outside the kept object, so no source gain. |
| R-S8-03 | Every A–C observable's real positive optical-ratio gate; source added, engine-qualified. |
| R-S8-04 | No curl-only live energy or internal angular-momentum/couple-stress result is consumed by the supplied local symbol; original admissibility question stays OPEN, no new source. |
| R-S8-05 | No physical rotational-stiffness frame is selected by LAB_HELD or radial flow; original frame question stays OPEN, no new source. |
| R-S1-03 | No Onsager–Casimir/microscopic time-reversal result is claimed or used; no source gain. |
| R-S8-06 | Displacement/density/advection material identity; existing S9b source retained without a live kinetic-law claim. |
| R-S12-01 | No non-passive interface response is used in A–C; held P1/P6 supply no A–C power law, so no reservoir/budget consumer added. |
| R-S12-02 | Native normal drain/return functions and their boundary data are not inputs to these conditional optical functionals. `jn` is a live coordinate exchange or an implied output; its identification with `v_dr` remains OPEN. No native drain-data source gain is inferred from a future profile solve. |
| R-O2-01 | No closure requiring computed net supply is kept in A–C; no net-power supplier/budget consumer added. |
| R-S21-01 | New in the original S9b pass: engine-specific conditional optical conditions/restrictions and their exchange/domain requirements. Retained, all qualifications above apply; no means of satisfying them selected. |

**Second route.** the repository identifies more than GR/PPN references. The inventory below
names the prior-art objects relevant to the optical, flow, embedding and carried mechanical
questions. `survey` means repository-root `docs/medium_requirements_and_prior_art.md` at
`36f59ff1`; `scoping` means `directives/O2_steady_brane_balance_scoping.md` §5 at that revision.
M15 retrieves both sources and `directives/prior_art_4d_bulk_search_prompt.md`; M10 retrieves
the pre-S9b register's O2 second route. The access qualifications below are the sources' reports,
not fresh literature verification. **No oracle is applied in the S9b record.**

| Identified object and repository source | Use in this pass | Open handoff / owner |
|---|---|---|
| GR/PPN deflection and logarithmic radar references, A–C spec (M1) | Supplied comparison inputs in the accepted constructions; no new oracle check or GR solution reproduced here. | OPEN: applicability and second-route questions for each named oracle; no owner named. |
| Cahill (2003) inflow; Hamilton–Lisle (2008), *The river model of black holes*; Sbitnev inflow; Unruh acoustic sinks — survey line 246, cautions 263–266; 4D search prompt | All identified, none applied. Survey scopes them to 3D with no off-sheet bulk and flags Cahill's status. No inflow profile or GM relation is imported. | OPEN: applicability and second-route questions for each named oracle. The river/acoustic-flow check, assigned by the written 2026-10-19 decision below (M21), ran after record acceptance (Codex-written, G1-reviewed, accepted `8d91f448`) and compared Hamilton–Lisle and Visser only. Its outcome and open questions are in the light-closing pass below. |
| Hamilton–Lisle (2008), AJP 76, 519–532; Visser (1998), *Acoustic black holes*; Unruh (1981), *Experimental Black-Hole Evaporation?* — scoping §5 | None applied. Scoping reports manuscripts read for Hamilton–Lisle/Visser, only abstract/metadata for Unruh. GR river and barotropic/inviscid/irrotational scalar acoustics require matching time, velocity, dimension, EOS/exchange and independent mass calibration; they supply no shear-sheet/off-sheet law. | OPEN: applicability and second-route questions for each named oracle. The river/acoustic-flow check, assigned by the written 2026-10-19 decision below (M21), ran after record acceptance (Codex-written, G1-reviewed, accepted `8d91f448`) and compared Hamilton–Lisle and Visser only. Its outcome and open questions are in the light-closing pass below. |
| Gurtin–Murdoch (1975), *A continuum theory of elastic material surfaces*; Scriven (1960), *Dynamics of a fluid interface*; Helfrich (1973), *Elastic Properties of Lipid Bilayers* — scoping §5; O2 register second route | None applied. Surface stress/traction, Newtonian surface-fluid and curvature-sector candidates; Gurtin–Murdoch abstract/bibliography only, Scriven/Helfrich full texts unverified. Ordinary 2D surfaces in 3D bulk do not establish this 3D-brane/4D-bulk law. | OPEN: applicability and second-route questions for each named oracle; no owner named. |
| Abels–Garcke–Grün (2012), *Thermodynamically Consistent, Frame Indifferent Diffuse Interface Models…*; Rodriguez–Hoger–McCulloch (1994), *Stress-dependent finite growth in soft elastic tissues* — scoping §5; O2 register second route | None applied. Diffusive exchange-momentum and evolving-reference/growth candidates; manuscript read for the former, metadata only/full text unverified for the latter. No two-fluid, mobility, growth or reference law is adopted. | OPEN: applicability and second-route questions for each named oracle; no owner named. |
| MacCullagh (1839) curl-only light / Whittaker Ch. V; Cauchy–Navier elasticity — survey attribution table; register R-S8-04/05 (M10) | None applied as an oracle. Existing light-carrier, longitudinal/compression, angular-momentum and objectivity obligations remain OPEN; the supplied optical ratio does not re-certify them. | OPEN: applicability and second-route questions for each named oracle; no owner named. |
| Structural acoustics (fluid-loaded plate, radiation impedance, added mass); Fruchart–Scheibner–Vitelli (2023) passivity/Onsager–Casimir; Cembranos–Dobado–Maroto (2003) branon fluctuation — survey attribution table | None applied. Adjacent interface, driven-response and embedding context; no interface law, drive or curvature dynamics is computed in A–C. | OPEN: applicability and second-route questions for each named oracle; no owner named. |
| Analog/emergent gravity (Unruh, BLV, BEC); Volovik/³He; Sakharov/Volovik induced gravity — survey prior-art table/unified findings | None applied. The survey's kinematics-versus-dynamics and undelivered-G qualifications remain; optical conditional requirements are not a derivation of gravity or G. | OPEN: applicability and second-route questions for each named oracle; no owner named. |
| Rubakov–Shaposhnikov (1983), Akama (1982), RS2, DGP (2000); material-medium superfluid-vacuum class (Volovik, Sbitnev, Consoli, Huang, Fedi) — survey 4D-bulk table | None applied. Adjacent embedding/material-bulk context only; no localization or material model is supplied to A–C by these names. | OPEN: applicability and second-route questions for each named oracle; no owner named. |

**Open prior-art handoffs.** Every named oracle in the table whose applicability check is
unperformed remains an individual OPEN item, including each member of a grouped row. The
river/acoustic-flow check assigned by the written decision below (M21) has since run for Hamilton–Lisle
and Visser (light-closing pass); every other check retains **no owner named**. Each carries both questions: does the oracle apply to this model point,
with its live profiles, domains, retained order and physical identifications; and, under the
register's second route, which of its failure modes or historical objections transfer here, and
have those objections been answered? Both questions remain open for each oracle. The source
access qualifications and applicability cautions in the table are not answers to either question.
The supplied GR/PPN comparison references do not constitute an applied oracle check.

**User decision, 2026-10-19:** the written source is
`directives/_measurements/S9b_record_r2_review_disposition.md` at `4760e6e2`, lines 19–21,
retrieved verbatim in M21. The user chose **"Queue it after record"** for the check of the
flow-only condition against the river/acoustic-flow prior art. The check is new computation,
done **after the record is accepted and outside it: Codex writes the check and it is reviewed
under G1**. This queues the check without selecting an oracle's applicability, importing a flow
profile, or answering its second-route question. The identified river/inflow/acoustic candidates
in the table remain open with that handoff; no oracle is applied here.

No oracle is applied and no oracle result is stated in this record. No new requirement entry is
entered on an unapplied oracle result; the second route remains open with the handoffs above.
Existing MacCullagh obligations remain OPEN; the conditional optical functionals and M13's bounded
condition comparison do not discharge them. Applicability checks remain unperformed, not negative
results. The register's earlier pass-2 disagreement stays unresolved (M10).

**Directive/register relation.** No new operative disagreement for S9b: the directive expressly
conditions the optical entry on this register's rest-on test and keeps postulates distinct.
The older pass-2 second-route disagreement remains quoted above, unresolved. O2's induced-measure
“relative O(ε)” statement and S9b D4's “No source supplies that order” are retained as a source
tension in S9b §8; this pass claims only coordinate-d³x exchange. Part D is held by the later user
scope despite its earlier planned description. These qualifications do not close any entry.

### Light-closing pass — the oracle comparison and the extra wave · 2026-10-19; accepted `e070c7f3`

**Sources:**
- `steps/S9b_flow_oracle_comparison.md`, on the accepted check `8d91f448`;
- the S11 departure (`steps/S11_stray_longitudinal.md`);
- photon survey §4, E32–E34 and E38.

Lookups: `directives/_measurements/light_sector_closeout_lookups.md`, C1–C11.

**Census:** 17 entries at `8d91f448` → 18, all OPEN (C10). One entry is new, `R-S22-01`. `R-S21-01`
gains an evidence line only; its source, requirement and status are unchanged. No other entry changes, and
no status changes.

**The user's decision.** On 2026-10-19 the user chose to close the light sector for now, after a
reassessment against the charter's goal ("ONE medium supports the LINEAR part of every force"). The user's
words: "follow your recommendations to wrap up the light sector before we move on to S12". This pass is that
wrap-up. Everything still open in the light sector travels in the entries above. The parked work is:
- S11c, PARTIAL;
- the heavy backlog runs, in `DEFERRED_HEAVY_RUNS.md`;
- the applicability checks of every other named oracle, which still have no owner named.

**First route / rest-on test.**
- **The oracle comparison** asks nothing new of the substrate. Its conditional correspondence rests on an
  identification the check does not supply (S9b's orbital `GM` equals each oracle's mass parameter) and on
  R-S21-01's own objects. So it adds evidence to `R-S21-01` and creates no new entry.
- **The S11 departure and the thermal entries** need an object that no record carried as a requirement:
  the throat's energy split. The phrases "thermal", "blackbody", "Stefan" and "N_eff" had no hits in this
  register, the plan or the defect register before this pass (C9). That object is now `R-S22-01`.

**Second route.** The comparison names four objections to the analog: what drives the flow, where the
material goes, scalar sound versus transverse light, and analog kinematics versus dynamics. None is answered
here. The comparison record lists them with their owners:
- the flow-drive question separates S12's drain source/boundary realization (`R-S12-02`)
  from the force balance in O2's carried handoffs;
- the material-destination question separates which density/exchange alternative the drain
  produces (S12, `R-S12-02`) from whether the knit can meet it
  (`R-S21-01`, S21), as the oracle record states;
- the third stays a scoped S9b exclusion;
- the fourth stays in the S9b pass's analog-gravity row.

These are handoffs, not entries. No owner is added to any other oracle check.

**Prior art.** The positive light results reproduce known prior art: MacCullagh's curl-only medium, and the
squared flow speed of the river model and of Visser's §7 acoustic flow. The flow match holds only under the
unsupplied mass, clock and speed-unit identifications. The register's
MacCullagh obligations (`R-S8-04`, `R-S8-05`) remain OPEN.


### S12 pass — local order-source compatibility and separate inventories · 2026-10-20; accepted `ac909aaf`

Source: `steps/S12_order_source.md`, §§1–7. Its constructions are accepted at `a3da43d9`, and the
comparator with its engine-row note at `fe485b73`. The record and register pass were accepted in
round 3 under "Accept, note both" at `ac909aaf` (`directives/_measurements/S12_record_r3_review_disposition.md`); J1/J2 are repaired here.
Commands and literal stdout: `steps/_measurements/S12_record_measurements.md`, generated by
`scripts/S12_record_measurements.py`. M01/M03/M07–M11 retrieve the conditional object and inventories;
M12 retrieves the schema, rest-on method and entry census; M13 the second-route sources.

**Census:** 18 entries before (`560c6c51`) → 18 after, all OPEN. **New entries by target: none.**
Four existing entries gain S12 source/provenance: `R-S1-01`, `R-S1-02`, `R-S12-02` and `R-S21-01`.
The last gains the user's framing decision only, not new optical evidence. No status changes.
The round-2 repair of `4c37473c` was also 18 → 18, with no additional source gain. It restored the sourced
S12 owners, including S12 in `R-S22-01`'s user-chosen target. `R-S12-01` retains its original
S12 register-inference/no-live-owner qualification. `R-S21-01`'s requirement field now explicitly
owes the comparison of force sectors' flow requirements, matching its by-target row and citing
the 2026-10-19 user framing; no comparison is performed here (M12/M15/M18).
The present repair of `fa5d1396` is also 18 → 18, with no additional source gain. It adds the
real-fraction/A13 qualification to the S12 handoffs and reports the accepted-pass edits below;
it changes no owner, computed object or comparison claim (M12 for the count; M19/M20 for the edits).

**Edits to already-accepted population passes.** Relative to the accepted register at `e070c7f3`,
the net edits to pre-existing Population passes are the following two edits to the light-closing
pass. Pass 1, Pass 2, O2 and S9b have no net text edit; their accepted blocks each have one literal
occurrence in the current register (M20). Existing-entry additions and sources are accounted for
in the census and rest-on table, separately from these amendments to accepted pass text.

| Edit | Accepted text at `e070c7f3` | Current text | Source that settles it / retrieval |
|---|---|---|---|
| Light-closing acceptance heading | **"### Light-closing pass — the oracle comparison and the extra wave · 2026-10-19; review pending"** | **"### Light-closing pass — the oracle comparison and the extra wave · 2026-10-19; accepted `e070c7f3`"** | `light_sector_closeout_review_r2_disposition.md`, Outcome: **"The light-sector closeout packet is accepted as revision 2"**; filed acceptance commit `e070c7f3`. Retrieve M20. |
| Light-closing second-route handoffs | **"the first two belong to S12 (`R-S12-02`) and to O2's carried handoffs;"** | **"the flow-drive question separates S12's drain source/boundary realization (`R-S12-02`) from the force balance in O2's carried handoffs;"** and **"the material-destination question separates which density/exchange alternative the drain produces (S12, `R-S12-02`) from whether the knit can meet it (`R-S21-01`, S21), as the oracle record states;"** | Accepted `steps/S9b_flow_oracle_comparison.md`, second-route items 1–2, lines 156–166: drain source/boundary data with S12 and force balance with O2's carried handoffs; what the drain produces **"is again S12's object (R-S12-02), and whether the knit can meet it belongs to R-S21-01"**. Retrieve M20. |

Both accepted and current blocks, with their literal line breaks, are retrieved in M20. The handoff
wording is a substantive amendment to an accepted pass, not only an acceptance-status update.
It preserves the oracle record's split and changes no computed result or owner.

**Two additional status-only edits retained from the S12 deliverable (J1).** The S9b oracle record's
status changed from "its own review is pending" to accepted at `e070c7f3`; the register header's
light-closing sentence likewise records acceptance at `e070c7f3`. These are separate from the two
light-closing Population-pass edits above. Source: `directives/_measurements/S12_record_r3_review_disposition.md`, J1.

**User decision, 2026-10-20: "Keep S12".** **"S12 stays open for the forms. This record closes only
S12's requirements part, and the forms come once the brane is defined, after the force sectors.
Every owner the sources and you gave stays as written, and the record says it delivers none of
the forms yet."** This replaces round 0's C1 owner obligation. S12 retains the native drain/return
functions, their boundary data, the throat's order-conversion source and its other sourced
assignments, all OPEN. None of those forms is delivered by this record (M18).

**First route / one rest-on criterion.** The reported object is conditional on supplied balances,
coordinate dimension, the adopted U1 support identification and the stated regional function class.
It does not rest on a selected interior-order, face, flux, drain/return, momentum-partner or energy
law. The printed compatibility objects do support a conditional requirement: any future realization
on this domain must meet those objects, or it cannot realize the reported source. They merge into
the native-conversion requirement `R-S12-02`, with no choice of how to meet them. Its status stays
OPEN: the source constraints supply neither native source functions nor boundary data.
S12 retains the forms/data under "Keep S12"; the following requirements constrain that future
delivery without selecting a form or claiming it has been delivered.
**Branch condition for every requirement below:** the printed compatibility objects rest on
the **real-fraction branch; A13 stays OPEN**. The decision list says **"The committed law is
written for a real ordered fraction; a phase-carrying branch would add content not computed here"**;
the plan requires **"the conversion balance must use the branch chosen at S1"**. These are the
real-fraction objects until S1's gate decides, with the sourced S1/S5/S12 owners retained. This
condition also qualifies `R-S12-02` and its by-target summary (M18/M19).

| Printed object required | Existing entry / target / status | Source and failure meaning |
|---|---|---|
| `general_order/region/1/source` below and `region/3/source` above | `R-S12-02` / S12 / OPEN | M03/M07: the respective live order-flux divergences. Data violating either cannot realize the reported exterior source. |
| `general_order/region/2/source` inside | `R-S12-02` / S12 / OPEN | M03/M07: live order transport plus order-flux divergence. Data violating it cannot realize the interior object. |
| `general_order/face/{1,2}/trace_coefficient` and `/source_contribution` | `R-S12-02` / S12 / OPEN | M03/M07/M08: normalized source parts at each live face, retaining the supplied-domain and diagnostic qualifications. No face-concentrated source is set to zero. |
| `material/face/{1,2}/requirement` | `R-S12-02` / S12 / OPEN | M03/M07: distinct total-material conditions at both faces; printed, not solved. Violation prevents locally conserved material in this function class. |
| `conditional_fast_ordering/region/{1,2,3}/source` and both face trace/source objects | `R-S12-02` / S12 / OPEN | M03/M07: only the labelled specialization, general in two offsets. Violation prevents that conditional specialization; it is not a universal fast-ordering limit. |

These requirement objects specify no law, material, premise or candidate as a means of satisfying
them. Their standing is distinct from the supplied balances, adopted inputs and flagged support
identification that qualify their domain. U1–U4 remain **"an adopted substrate input to a conditional
model"** each, as the source says; they are not newly derived substrate requirements. U1's unencoded
re-ordering clause, face-concentrated input fields and the generality of the live functions remain
OPEN handoffs. No named retirement condition makes an uncomputed constitutive form a new entry.
The two existing dimension/shear obligations have live retirement questions already, so they gain
only the limited conditional consumers stated in their own entries.

**Rest-on census, every entry:**

| Existing entry | S12 dependence / effect |
|---|---|
| R-S6-01 | No compression-modulus input or result; no source gain. |
| R-S7-01 | Live faces impose no selected width or nonuniform-width flatness result; no source gain. |
| R-S1-01 (target S6) | Three in-plane coordinates and four-dimensional material balance; coordinate-only source added, no dimensional derivation. |
| R-S1-02 | Adopted U1 shear-region/support identification; conditional interpretation source added, no computed shear response or confinement. |
| R-S1.5-01 | No GNLS mode/acoustic input; no source gain. |
| R-S6-02 | No stationary/equilibrium state or substrate first variation; no source gain. |
| R-S8-01 | No stiffness functional supplied or derived; no source gain. |
| R-S8-02 | No coupled `(u,h)` operator or block-diagonality input; no source gain. |
| R-S8-03 | No optical-stiffness sign calculation or propagating optical gate; no source gain. |
| R-S8-04 | No angular-momentum/couple-stress result used; inherited inventory only, no source gain. |
| R-S8-05 | No rotational reference frame selected; inherited inventory only, no source gain. |
| R-S1-03 | No microscopic time-reversal or Onsager–Casimir result used; no source gain. |
| R-S8-06 | No displacement identification or quadratic slab inertia used; material conservation alone does not supply either; no source gain. |
| R-S12-01 | Separate conditional non-passive-interface reservoir/budget retained. No computed S12 power consumer, no source gain, no live owner added. |
| R-S12-02 | Printed regional/face conversion compatibility and distinct material face conditions merge here; S12 source added. Native functions/data remain OPEN. |
| R-O2-01 | Separate net-supplier/budget condition retained; no S12 net-power calculation or consumer, no source gain or owner change. |
| R-S21-01 | User framing provenance added: each force states its flow requirement, S21 compares, a mismatch reopens S9b. No new optical evidence or fitted flow. |
| R-S22-01 | No throat energy split, matter coupling or thermal/cosmological transfer result computed or used; carried inventory only, no source gain. |

**Second route.** Apply the register's broader rule, including prior-art candidates outside S12's
record. Read the O2 scoping §5 and the S9b second-route inventory (M13). Gurtin–Murdoch,
Abels–Garcke–Grün, Scriven, Helfrich and Rodriguez–Hoger–McCulloch identify surface/exchange/reference
candidates; Hamilton–Lisle, Visser, Unruh and the inherited Cahill/Sbitnev class identify flow
candidates. Their model/dimension/exchange and source-access restrictions travel with each row.
For **each named candidate individually**, applicability, historical failure modes/objections that
transfer and whether they are answered here remain OPEN, **no owner named for a new S12 check**.
S12 reproduces no identified candidate result: its balances are supplied and its support model is
adopted. Thus no derived S12 requirement is entered on an oracle match, and none of those historical
obligations is claimed discharged. The remaining S9b candidates each keep their original OPEN
handoffs; the positive light-sector oracle check does not become an S12 oracle check.

The inherited river/acoustic objections are still carried, not answered: flow drive separates
S12's drain/source-boundary realization (`R-S12-02`) from O2's carried force-balance handoffs;
material destination separates S12's drain-produced density/exchange choice (`R-S12-02`)
from S21's knit compatibility (`R-S21-01`). Scalar acoustics versus transverse light is a scoped
S9b exclusion; kinematics versus dynamics stays in the analog-gravity handoff. Local material/order
identities name compatibility but select no flow, closure or response dynamics. Existing MacCullagh
admissibility/frame entries stay OPEN. No fresh oracle computation or literature check is done here.

**Directive/register relation and source tensions.** No new operative disagreement with the S12
record contract: its item 7 expressly defers to this rest-on test. Its 2026-10-19 framing says S12
**"states what any conversion law must satisfy and leaves the forms open"**, while the older
`R-S12-02` requires **"the native drain/return functions and their separate boundary data"**.
The older obligation remains OPEN; the current constraints merge into it without claiming its
delivery. The earlier owner line **"S12 owns the drain/return data"** and the requirements-only
framing are kept side by side. The user's 2026-10-20 **"Keep S12"** resolves their ownership
reading: **"This record closes only S12's requirements part"**, and the forms **"come once the
brane is defined, after the force sectors"**. S12 stays open; the record delivers none of its
forms. The reviewed record's replacement **"OPEN, no live owner named"** disagreed with the
oracle's **"The drain's source and boundary data belong to S12 (R-S12-02)"**; its removal of
S12 from the energy-input target disagreed with the user-chosen **"S22 (the nonlinear throat
solve), together with Q2 (the mouth source and its held datum) and S12 (the throat's
order-conversion source and its boundary data)"**. The record quotes both sides of these
disagreements and the A13 owner disagreement; all sourced owners are restored (M18).
Stage 045's controllers close globally **"despite the local sink"**, whereas the current
law supplies **"total material"** conservation and rules out **"The frozen-wall mass sink with a
remote return"**. Whether the controllers retain a role after re-hosting is OPEN in the
S12 source/controller inventory. Both sides and the split oracle handoffs are retrieved in
M01/M11/M17/M18.
The earlier pass-2 second-route disagreement stays unresolved: the old directive says **"Run it only
on a prior result that a record in the reading list identifies as reproduced by the sector"**, and
this register says **"Run this route for every sector with identified prior art"** (M13). The
register's broader rule governs this pass. All comparator differences/coverage limits stay in
S12 §4 and M06–M09; the accepted engine-row comparison carries no information.

### S14a part 1 pass — closed projection and two-channel shell compatibility · 2026-10-21; accepted `4a4dfaf5`

Source: `steps/S14a_1_projection.md`, §§1–8. Constructions accepted `94fddf2f`; comparator accepted
`de81b9fc`. This record/pass is accepted at `4a4dfaf5` (record review round 1). Commands and literal stdout:
`steps/_measurements/S14a_1_record_measurements.md`, generated by
`scripts/S14a_1_record_measurements.py`. M01/M03/M05–M10 retrieve the conditional printed objects,
their domains and comparison limits; M12 retrieves the schema/rest-on method and census;
M13 the second-route sources; M16 this pass and the exact ancillary diff against `f5d283c1`.

**Census:** 18 entries before (`f5d283c1`) → 18 after, all OPEN. **New entries by target: none.**
Four existing entries gain scoped source/provenance: `R-S1-01`, `R-S1-02`, `R-S12-02`, `R-S21-01`.
The last gains user-framing provenance only, not new optical evidence or a selected flow.
No target or status changes. The S12 pass itself changes only acceptance status/provenance and
known J1/J2: its account now names the S9b oracle status edit and register-header light-closing
status edit; its round-3 `18 → 18` sentence cites S12 M12. S12 acceptance is `ac909aaf` (M15/M16).

**First route / rest-on test.** The printed objects rest on supplied material/order balances,
the coordinate dimension, adopted U1 support identification and the live regional/two-face function
class. Spherical symmetry is a labelled restriction, not an earned limit. Each is kept at that
standing, not entered as a new derived substrate property. No chosen interior, face, source,
flux or profile form is required by these identities. Their printed compatibility objects do
constrain future data that would realize the conditional projection; these merge into the
existing `R-S12-02` compatibility requirement. No method of satisfying them is prescribed.

**Common domain of every handoff below:** the real-fraction branch, **A13 OPEN** at S1's gate,
propagated through S5/S12; two live graph faces and the inherited regional field class; full
endpoint-concentrated parts. Shell objects additionally keep the labelled spherical specialization;
velocity objects keep nonzero projected density. Every comparison qualification in the record §5
travels with its object. In particular the seven exact integral residuals, shell-domain difference,
external dimension-context differences, unjoined decompositions and unpaired components are not
resolved by this pass. Source/flux forms and distinct boundary data remain S12's, OPEN under
the 2026-10-20 **"Keep S12"** decision. These requirements deliver none of the forms.

| Required printed object; native tag family / bound references | Entry / target / status | What fails if future data violate the object |
|---|---|---|
| `General/{material,order}/projected density`, `/projected in-plane flux`, `/closed projected balance/operands`; corresponding Wolfram `general projection/...` objects (b51/b52/b59/b73/b74/b79) | `R-S12-02` / S12 / OPEN | The reported closed-slab material/order balance is not realized on this domain. |
| Each `General/{material,order}/face {0,1}/{interior leaving term,closed-projection leaving term}`, `/concentrated face parts`; material face-equation coefficients; `General/order/projected source` (b49/b53–b58/b71/b75–b78/b82) | `R-S12-02` / S12 / OPEN | The reported relative face transport, full concentrated parts or projected source are not realized. Material face requirements remain distinct; none is solved or imposed by this pass. |
| The separately constructed `SphericalSpecialization` projected counterparts (b127/b129–b137/b147/b149–b155/b158) | `R-S12-02` / S12 / OPEN | The labelled radial projection is not realized. No transfer outside its restriction is required. |
| Spherical `face channel`, `core channel`, `J`, `order lost`, `outward outer flux`, `storage`, `shell balance`, `face channel minus order lost` (b166–b171/b173–b175; Wolfram `spherical specialization / shell/...`) | `R-S12-02` / S12 / OPEN | The reported two-channel definition and shell/source/storage relations are not realized. Storage and shell-balance exact comparison remain unresolved as recorded, not replaced with an exact-match claim. |
| Spherical `outer flux radius derivative`, `outer projected velocity`, `outer volume flux`, with the printed predicates (b172/b176/b177) | `R-S12-02` / S12 / OPEN | The reported live-radius flux/conditioned-velocity relation is not realized. No zero-density velocity value is owed. |

The difference `face channel minus order lost` is required **as the printed difference**, not set
to zero. Agreement of bound dimension values constrains their interpretation only at their
printed domains; primitive dimension assignments and unpaired audit trees are not new laws.
Closed differential-operand integrals retain b50/b72/b128/b148's unresolved exact comparison
and printed sample evidence: of each unresolved part's 54 stored samples, only the 18 quadratic-form
samples can move when a derivative is taken at another point; the 18 affine and 18 cubic samples
do not make that check. The same limitation qualifies b171, b173's unresolved lhs and b174;
b173 has 36 quadratic samples across both sides, of which only 18 belong to its unresolved lhs.
No new exact equivalence is required by assertion (M03/M07–M10/M17).
The pending controlled return law, mass link, core solution, source-map completion and profile
knit enter no completed part-1 result as supplied forms; they retain their step assignments and
do not generate duplicate entries merely because they remain OPEN (record §7, M01/M11).

**Rest-on census, every existing entry:**

| Entry | S14a part 1 dependence/effect |
|---|---|
| `R-S6-01` | No compression modulus used; no source gain. |
| `R-S7-01` | Live faces select no width/flatness result; no source gain. |
| `R-S1-01` | Three in-plane/one bulk coordinates and spherical area; coordinate-domain consumer only; source added, no dimension derived. |
| `R-S1-02` | Adopted U1 shear-region/support identification; conditional interpretation source only. No shear response/confinement derived. |
| `R-S1.5-01` | No GNLS linear-mode/acoustic input; no source gain. |
| `R-S6-02` | No equilibrium/first-variation input; no source gain. |
| `R-S8-01` | No stiffness functional used; no source gain. |
| `R-S8-02` | No coupled displacement/width quadratic operator used; no source gain. |
| `R-S8-03` | No optical-stiffness sign/gate result; no source gain. |
| `R-S8-04` | No new angular-momentum/couple-stress consumer; inherited second-route obligation remains. |
| `R-S8-05` | No rotational frame selected or used; inherited second-route obligation remains. |
| `R-S1-03` | No microscopic time-reversal/Onsager input; no source gain. |
| `R-S8-06` | No displacement/quadratic-inertia identification used; number conservation does not supply it. O2 comparison remains OPEN; no source gain. |
| `R-S12-01` | No non-passive response or computed power consumer; distinct carried condition, no live owner newly assigned. |
| `R-S12-02` | Projected/face/source and labelled shell compatibility objects merge here; source added, native forms/data remain OPEN. |
| `R-O2-01` | No net-power/supplier result or adopted net-supply closure; no source gain/owner change. |
| `R-S21-01` | Framing provenance only: each force states flow needs, S21 compares, mismatch reopens S9b. No new optical profile evidence or comparison result. |
| `R-S22-01` | No energy split or conversion-work/thermal-transfer result used; no source gain. Core exclusion keeps its existing owner, not a new computed energy requirement. |

**Second route.** Apply the register's broader rule for every sector with identified prior art.
Read O2 scoping §5 and the O2/S9b/light-closing source rows (M13), rather than treating this
record's own questions as an exhaustive prior-art inventory. Each named candidate remains
individual: Gurtin–Murdoch, Abels–Garcke–Grün, Scriven, Helfrich, Rodriguez–Hoger–McCulloch,
Hamilton–Lisle, Visser, Unruh and the inherited Cahill/Sbitnev class. For each, applicability,
transferable historical failures/objections and whether they are answered here remain OPEN,
**no owner named for a new S14a part 1 check**. The original access, dimension, material and
exchange limitations remain as the source rows state them. No part-1 result is claimed to
reproduce an applied oracle, so no new requirement is founded on such a match.

Every other candidate in the inherited S9b inventory retains its source row and original OPEN
handoff (M13). The river/acoustic objections are carried, not answered: S12 drain/source-boundary
realization versus O2's force balance; S12's produced density/exchange versus S21 knit compatibility;
scalar sound versus transverse light as S9b's scoped exclusion; analog kinematics versus dynamics
as its existing handoff. Neither global nor local conservation chooses a drive or response law.
MacCullagh's `R-S8-04/-05` remain OPEN. No literature search, new oracle check or candidate
selection is performed by this pass.

**Directive/register relation and source tensions.** The current contract expressly defers to
this rest-on test and both-route method; no operative disagreement is found. The historical
second-route directive says **"Run it only on a prior result that a record in the reading list
identifies as reproduced by the sector"**, while the register says **"Run this route for every
sector with identified prior art"** (M13). The register wins. The material-conservation and
S13-filing disagreements are quoted with both sides in record §6 (M01/M11/M15), not resolved by
new physics or an owner reassignment. The plan's four bridge objects are not all closed:
S14/S15/S20 remain **"CONDITIONAL, not carried results"**. No frozen DC completion or profile
is adopted to close them.

### S14a part 2 pass — dynamically selected slow limit and far behaviour · 2026-10-21; accepted `79511272` (round 2 clear)

Source: `steps/S14a_2_slow_limit.md`, §§1–7. Commands and literal stdout:
`steps/_measurements/S14a_2_record_measurements.md`, generated by
`scripts/S14a_2_record_measurements.py`. M01/M02 retrieve the contract and both step-list
revisions; M03 the review dispositions, lookup existence and assembly-audit excerpt; M04–M10 each
witness's complete source and saved stdout; M11 the accepted part-1/S12 sources and owners; M12 the schema/method/census;
M13/M14 the plan/comparison inventories and second-route handoffs; M15 this pass and the
exact ancillary diff against `978c359c`; M16 the record round-0 dispositions, V2 user decision
and committed-generator provenance; M17 the record round-1 dispositions. No build or witness rerun is performed here.

**Census:** 18 entries before (`978c359c`) → 19 after, all OPEN. **One new entry by target:**
`R-S14a-01`, jointly S12 (source/flux forms), S8/S22 (material law), and **no owner named** for
face/interior laws, as V2 assigns. **Existing entries gaining a source: none.** No existing
entry changes target/status; `R-S12-02` retains its block and by-target row at `978c359c`.
The part-1 pass changes only its acceptance heading/provenance to `4a4dfaf5`; its other
content is unchanged.

**First route / rest-on test.** The user chose **"Record as requirement"** after the two
step-list reviews. The object required is **the slow limit and far behaviour selected by
the slab's own dynamics**, for S14a's still-underived controlled return and profile-dependent
`J→Q` obligations. Its consumers are **S14, S15 and S20**: without that delivery they remain
**CONDITIONAL, not carried results**. This is a named future delivery condition for those
claims, not an assumed form or a defect in an accepted part-2 result. No such result exists.
The register method makes an assertion about a step with no record a requirement; this one
is explicitly required by the user, not entered merely because a future solve could use it.

The register's merge test is applied to the selection object and its joint deliverers,
not just S12's forms. No existing entry covers that same requirement with those deliverers.
`R-S12-02` covers S12's native source/flux forms, separate data and compatibility objects;
it does not cover selection owed jointly with S8/S22 and unowned face/interior laws.
`R-S22-01`'s joint energy-split obligation and `R-O2-01`'s unowned supplier/budget obligation
are different objects. The method keys each entry by its deliverer, so selection has its own
entry, **`R-S14a-01`**; `R-S12-02` is unchanged from `978c359c` (M12/M17).

The user's V2 choice **"All the slab's laws (Recommended)"** (2026-10-21) replaces item 4's S12-only target
sentence. Its source is `directives/_measurements/S14a_2_record_r0_review_disposition.md`,
"The user's decision on V2", committed at `12aa65d5` (M16). The selected orchestrator option says:

> "The selection is owed jointly by S12's source and flux forms, the slab's material law
> (S8/S22), and the face and interior laws, which have no owner yet. The register says so and
> states the missing owner as missing; it invents none. The record reports the disagreement
> with both sides quoted."

The new entry's target field and by-target row state those joint deliverers for selection. Existing S12
source/flux/boundary obligations remain with S12 under **"Keep S12"**. No missing owner is
invented, and the selection requirement does not repartition the complete return-law/map
obligations. Recording it delivers none of the forms, the slow return law or the source map.
The inherited real-fraction/A13, live-face and part-1
comparison qualifications remain. No new averaging operation, solution family, density gate,
limit order, drain/return law, material, profile, premise or candidate is prescribed.

The evidence is pass 1 Q1–Q3 and pass 2 R1–R3 as filed, not an accepted build or a universal
theorem: three pass-1 Codex instruments, one pass-1 Grok instrument, two pass-2 Codex
instruments and one pass-2 Grok instrument. Their constructions/checks differ, as record
§§2–3 states. In particular, pass-2 Grok's gate example checks one projected material balance,
not both native balances. The pass-2 orchestrator reading is quoted as a reading. The user
decision, not an extrapolation from these probes, supplies the requirement's standing.

All 18 existing entries retain their sources/targets/statuses. No new computed consumer of
their dimension, shear, stiffness, inertia, equilibrium, mode, width, power or energy-split
objects is supplied by this retrieval-only record. The repeated flow-framing decision keeps
S21's existing comparison obligation; no new S21 source/comparison result is added.

**Second route.** Apply the register's broader rule for every sector with identified prior
art. M14 retrieves O2 scoping §5 and the O2/S9b/light-closing/part-1 source inventories. For
each candidate individually—Gurtin–Murdoch, Abels–Garcke–Grün, Scriven, Helfrich,
Rodriguez–Hoger–McCulloch, Hamilton–Lisle, Visser, Unruh, Cahill and Sbitnev—applicability,
transferable historical objections/failures and whether those are answered here remain OPEN,
no owner named for a new part-2 check. All other individually named S9b candidates keep their
original OPEN handoffs and filed access/model qualifications. No candidate is applied or
selected, no literature verification is performed, and no part-2 prior result is reproduced.

The four inherited river/acoustic objections remain distinct: S12 source/boundary realization
versus O2's carried force balance; S12 produced density/exchange versus S21 knit compatibility;
the scoped scalar-sound/transverse-light exclusion; and analog kinematics versus dynamics.
None is answered here. MacCullagh's `R-S8-04/-05` stay OPEN. No new requirement rests on an
unapplied oracle; both-route questions are preserved, not treated as negative findings.

**Directive/register relation and disagreements.** The current contract defers to this schema,
rest-on test and both routes. Its original S12-only target sentence and the earlier option's
**"(the open S12 forms)"** gloss disagree with accepted S12's material/face/interior ownership
statements; both sides are quoted in record §6 (M01/M11/M16). The V2 user decision quoted above
supersedes that target sentence with joint delivery. The historical population
directive says **"Run it only on a prior result that a record in the reading list identifies
as reproduced by the sector"** (direct retrieval M11); this register says **"Run this route
for every sector with identified prior art"** (M12). The register governs. The pointwise-rate,
DC/slow-limit, historical material-conservation and S13-filing tensions remain quoted with
both sides in record §6 (M02/M03/M11/M13). The contract requests **"a committed generator script"**
and says **"No commit and no push"**; the list's own no-commit instruction governs this repair.
The round-1 reviewed generator and measurements are preserved at `85904c4e`; the repaired record,
generator, measurements and register pass were subsequently committed and accepted at `79511272`
(round 2 clear; `directives/_measurements/S14a_2_record_r2_review_disposition.md`, retrieved in S16 M15).
No physics disagreement is closed by this pass.

### S16 pass — live moving-region balances and separate response requirements · 2026-10-24; accepted at `4f907e9d`

Source: `steps/S16_moving_region.md`, §§1–8. The governing list is revision 4 at `0464b80a`;
the constructions are accepted at `1d990ee6`, the comparator at `f08db34a` with its citation
containment claim. Commands and literal stdout: `steps/_measurements/S16_record_measurements.md`,
generated by `scripts/S16_record_measurements.py` (M01–M19). Round 0 is preserved at `4ae7e6cf`
and the reviewed Y1–Y4 repair at `fa1fa413` (not accepted). The W1 user decision is recorded
at `4e5a94bc`; the W1 repair and S16 record/pass were accepted at `4f907e9d` (M17/M19).

**Census:** 19 entries at `4d4d22cb` → 25, all OPEN. **Six new entries by target:** three with
no physical owner named (`R-S16-01/-02/-03`; compactness is also owed to S22's debt table),
one S8 (`R-S16-04`, retaining the O2 antecedent qualification), one S8/S22 (`R-S16-05`),
and one S22 debt-table account (`R-S22-02`, without assigning missing physical suppliers).
**Four existing entries gaining source/provenance:** `R-S1-01` (coordinate dimension only),
`R-S1-02` (inherited support identification only), `R-S12-02` (supplied projected material object
and "Keep S12" only) and `R-S21-01` (the 2026-10-19 framing only). None changes requirement,
target or status. The Y1–Y4 repair is **25 → 25**, with no new entries and two additional
existing-entry S16 sources. This W1 repair is also **25 → 25**: no new entries or additional
S16 consumer sources; four entries gain the user's 2026-10-24 admission decision. No entry
closes; `R-S14a-01` remains OPEN with its existing joint target.

**First route / rest-on test.** As the schema says, a requirement is a future obligation;
an adopted value/form belongs to the ansatz ledger, and an already-wrong object to the defect
register. The selected user options expressly file the passive needs (`R-S16-01`) and inertial
relation (`R-S16-04`). The orchestrator's list separately files compactness from E3/F3
(`R-S16-02`), face momentum from E2 (`R-S16-03`), density response from O2 `:138`
(`R-S16-05`) and the S22 debt account (`R-S22-02`). Those four retain their original sources,
filers and targets, but their admission now rests on the user's **2026-10-24** selection of
**"Require all four (Recommended)"**. The option text, written by the orchestrator and selected
by the user, is:

> You require all four on the register now, and each entry cites today's decision as its reason. The register is
> where we collect what the brane must support before designing it, and these are laws or properties the slab and
> bodies must have. The record also notes that my decision list and the register disagreed and that you settled it.

Each of `R-S16-02`, `-03`, `-05` and `R-S22-02` quotes the dated decision in its admission
field (M19). It adds an admission ground, not a physical owner or a selected law, form or material.
No substrate response is derived and no new dependency is imposed on the regional balances.

**Record-list/register disagreement, settled by the user.** Record-list item 7 says
"The record and the pass file each one as that section states"; item 9 says "Where this list
and the register disagree, the register wins, and the disagreement is reported" (M01).
The register's generic method says "Each becomes an entry keyed by the step that must deliver",
but its earlier applications say "A requirement is not everything a later step would find useful",
"an OPEN item sources no entry under the rest-on test", and "An unspecified response form needed
only to select or solve a later closed material state is an OPEN handoff, not an additional delivery
condition of this unsolved named object. A future closure owner alone does not satisfy this test."
O2's density row explicitly says "No new density-response entry". S14a part 2's forward object
instead "is explicitly required by the user, not entered merely because a future solve could use it"
(M14/M19). The prior S16 admission on item 7 followed one side of that disagreement and did not
establish admission under those earlier applications. The user's 2026-10-24 decision settles
admission of these four as explicit requirements. Item 7 is not their admission ground; the
earlier passes are not retrospectively changed and no general ambiguity in the method is decided.

The same rest-on test admits two dependencies of the present conditional object, inherited
through the supplied part-1 balance. Its three in-plane coordinates and shell measures consume
the existing dimension obligation `R-S1-01` at S6; failure to deliver/reaffirm that domain
prevents unchanged transfer. Its shear-holding-region interpretation consumes the flagged,
adopted U1 support identification in `R-S1-02` at S1; failure removes that interpretation,
not the conditional algebra. They merge as scoped sources, as the earlier passes did. No new
dimension, shear law, confinement result or enforced U1 re-ordering clause follows (M01/M04/M14/M17).

No further entry has a printed object as its requirement. The supplied material/momentum balances,
their integrals, count-centre definitions and moments are conditional identities on the list's
domain, not a selected material state. No printed nonzero boundary flux is required to vanish.
The count projection and its S12 provenance already fall within `R-S12-02`; source and face forms
stay OPEN. General `G`, `Pi` and `f` forms have no owner assigned here, while S12 gives slab material
completion to S8/S22 and names no face-law owner. Neither the unowned momentum balance nor the
core's excluded content is assigned to S12. Control shape, moment truncation and independent
dimension units are labelled supplied/construction choices, distinct from the inherited coordinate
dimension and support identification tested above. No O2 premise is silently adopted.

**Merge test.** `R-S16-01` is distinct from source normalization and from S21's knit;
`R-S16-02` is distinct from count moments and the width bridge; `R-S16-03` is distinct from S12's
material-face compatibility. `R-S16-04` does not widen the original homogeneous/slab object in
`R-S8-06`. `R-S16-05` has no existing density-response entry. S22's existing `R-S22-01` concerns
energy partition, not the reduction-debt account. Compactness is cross-referenced there once.
The record's requirements table names each target and OPEN status; no means of meeting them is filed.

**Second route.** Read every gravity prior-art row G1–G13 with its stated conditions and access
qualifications (record §5, M05). G1/G2 attribute the standard transport algebra; their abstract
any-entity statement is compared only by retrieval with S16's integral/measure and rest/motion
objects. No external result is imported as a premise and no new residual is computed. G3/G4 and
G5/G6/G8 have other objects or conditions; G7 is the prior result for the open passive requirement,
not an identification; G9/G11 give no S16 force-law result; G10 is an unexecuted architectural
sketch; G12 gives posits, not a mass derivation. The quoted condition table in record §5 keeps
each comparison open where computation or a condition change would be needed.

G13's objections are each preserved: moving-drain damping (Wang), unexplained energy disappearance,
drag and matter creation (the last three from a secondary source), and universal pulsation phase
(secondary source). S16 supplies no energy balance, mass-creation relation, sink-force law or
pulsation model and answers none of them. Drag remains a question beside the passive/inertial
requirements without an assigned applicability-check owner; energy destination remains beside
S12's exchange and `R-S22-01`'s split without a new energy result; mass creation remains beside
C3 and S22's debt account; synchrony is an unapplied-oracle question with no owner named. None is
a defect found here or a new assumed property of an unbuilt substrate. No new entry follows merely
from these objections before applicability/reproduction under matching conditions is established.
This applies the broader "every sector with identified prior art" rule, not the historical
reproduced-result-only restriction. G7's no-external-body-force assumption is retained; O2's
different drive premise is not imported into S16 to make it match.

**Limits and disagreements.** The pass inherits the complete record §3 ledger: exact/as-computed
versus resolved residuals, open domain differences, native text/route differences, every unjoined
tag and unpaired component, below-count bindings and normalized-domain exclusions. It does not
claim whole-stream identity or a motion-free sum of `localEvolution` pieces. The plan/O2 interface
describes a compact-body calibrated-potential reduction; the later user list supplies no potential,
mass or compactness and computes only regional objects. Both sides are quoted in record §7. The
plan's S19 trust gate and charter-feedback handoffs remain OPEN beside that delivery; record §4
and `STATUS.md` retain each with no owner named for its disposition. Neither is decided here.
The list heading "Requirements the user's decisions file" over-attributes its additions:
the selected options file passive-factor needs ("goes to the register as a requirement") and
the inertial relation ("The law it needs goes to the register, as O2 keeps it (S8)"). Compactness,
face momentum and density response are the list author's original filings at `:158–166`, not
2026-10-21 selected user options. The separate 2026-10-24 decision now admits all four list-author
filings without changing those filers. Record §7 quotes both sides; the individual entry
sources/admission grounds govern (M01–M02/M04/M14/M17/M19). The
governing list's "committed generator script"/"No commit and no push" conflict is reported there;
the latter governs this authoring STOP. The list's revision heading still says "Revision 3" while
its body and acceptance say revision 4. Accepted S12's historical mass-link owner is "S13";
the current S16 list retains "C3, a gap", and the later user filing folds S13 into S14a.
Record §7 quotes those distinct source assignments without inferring a new mass-link owner.
No source is repaired by this pass. The record-list item-7/item-9 disagreement with the earlier
rest-on applications is reported above and settled for these four entries by the user's dated
decision (M19). Original filers and owners remain distinct from that admission ground.

### S15 pass — conditional flow around two moving drains

The first route adds the gravity-side test to the existing S21 comparison in `R-S21-01` on the
user’s “Keep it general” decision. The three net family crossings retain unspecified fields.
The new provenance covers the signed shell observable, its domain and the required total
departure if that observable is to scale as the inverse square. The body-force identification,
free acceleration and all-separation profile classification remain open.

The second route checked G7, G9 and G13. Their body-force and damping results have different
observables or assumptions. The existing material and supplier budgets and S19 retain the
relevant open questions. The S14a/S16 requirements and S22 source and normalization debts keep
their recorded owners. Census 25 → 25, all OPEN. One existing entry gains a source, `R-S21-01`
(target S21). No new entry by any target. [Record, prior art and handoffs](steps/S15_two_drain_region.md)
· [Native evidence](steps/_measurements/S15_record_measurements.md)

### Method for population passes

Method for a pass, per step record: read it for every statement of the form *this assumes*, *this rests
on*, *this is postulated at*, *this is deferred to*, *X enters at S**n***, and for every named retirement
condition. Each becomes an entry keyed by the step that must deliver. Where the record asserts something
about a step that has no record yet, that is a requirement by definition.

Two things to watch for, both seen already:
- A requirement can point **sideways**, not only backwards — S11's questions reduce to the S11b interface
  law, which is a sibling step, not a substrate one.
- The same object can be required by several sectors. Charge and magnetism ride the same brane–bulk
  coupling as light, so an entry may gain sources rather than being duplicated.

⭐⭐ **A SECOND ROUTE, added 2026-08-07 — ⛔ the method above would never have found `R-S8-04` or `R-S8-05`.**
Both came from asking, of a **known prior result our sector reproduces**, *why was it rejected in its own
time, and has that objection been answered here?* ⛔ Neither is anywhere in S9's or S10's records, because
the records only capture what their authors thought to doubt — and the sector reproduces MacCullagh's
algebra so cleanly that the objection to it never came up.

⇒ ⭐ **Run this route for every sector with identified prior art**, ⛔ not only the light sector.
⚠ It is the sharpest use of `CLAUDE.md` rule 16 available: the prior work's **failure modes** transfer as
obligations even where its **results** transfer as corroboration — and a result that matches prior art
tells you nothing about whether you inherited its problems.
⚠⚠ **The tell that this was overdue:** the closer a sector's agreement with prior art, the *less* anyone
thinks to ask what that prior art could not do. ⇒ ⛔ a clean reproduction is exactly when to run it.

**Honest sizing, now measured rather than guessed.** Pass 1 read two records and produced **seven**
entries, bringing the file from two entries to **nine total**, from a third of the sector. ⚠ That rate
is the argument for
running the pass **as each step closes**: reconstructing it at Phase 5, across five sectors, is precisely
where an implicit assumption goes missing.

⛔ **And pass 1 found no way to make this mechanical.** Every one of the seven came from reading prose for a
*flagged identification*, a *postulate with a live retirement condition*, or a *consequence stated about a
step that does not exist yet*. ⇒ ⛔ there is no grep for it, and a pass that produces entries quickly is a
pass that missed some.

## When this is done

Phase 5 knits the sectors and builds the substrate. At that point this file is the checklist the
substrate is tested against — the thing that answers *"is the brane–bulk we built sufficient to carry all
five forces"* by enumeration rather than by recollection.

## What two inward drains must reproduce in the Newtonian limit

**0PN scope, clarified 2026-11-14.** Take both mouths to be open and draining inward. As clarified on 2026-11-16, the isolated inward mass rate at fixed ambient conditions is a single function of mouth radius, $D(a)$. A known rate in the 0PN comparison means the rate for that radius under this common law. A controlled change of rate at unchanged radius tests sensitivity to a prescribed boundary condition. Compare the target mouth alone and with a companion at separation $d$, under the same distant brane and replenishment conditions. The gravity calculation ends at the opening. Its observable is the companion-induced acceleration of the free mouth. Any interior contribution at leading order must appear in the [response measured at the mouth](#what-a-later-throat-model-must-reproduce-at-the-mouth). The internal derivation is deferred. The mechanism that maintains the opening belongs to Q2/S22. [Source excerpts](_measurements/newtonian_drain_requirements_sources.txt) · [V3 step plan](V3_STEP_PLAN.md)

**Current deliverable.** Gather the physical inputs and independent checks for a future free two-mouth simulation. The orbit run is deferred until a larger machine is available. Each mouth radius must select its drain rate, and its interaction with the brane must determine its motion. The Newtonian orbit supplies the eventual comparison.

| What the 0PN gravity test needs | How to determine it |
| --- | --- |
| **A field produced by the companion.** | Use the companion's radius and the intake selected for that radius to calculate the brane flow incident on the target. Relate the isolated mass rate $D(a)$ to the exterior volume rate $Q$ through the flux balance. The reduced exterior with prescribed strengths gives an attractive inverse-square matter-current contribution. [Slow-limit record](steps/S14a_2_slow_limit.md) · [S15 scope](directives/S15_decision_list.md) |
| **Force at the mouth.** | Over the actual opening, integrate the companion-directed momentum carried by the entering brane material and the stress transmitted across it. Use the local crossing speed relative to the mouth. Subtract the same reading for the target alone. The flat-mouth potential-flow example gives an attractive inverse-square load for prescribed radial intake and scalar pressure. [Finite-mouth calculation](../../docs/drain_driven_initial_acceleration.md#paired-minus-alone-load-across-the-3d-mouth-at-0pn) · [Mouth integral](../../docs/drain_driven_initial_acceleration.md#sideways-momentum-entering-the-mouth-at-unchanged-total-intake) |
| **Free motion of the mouth.** | Relate the change in load to the translational response of the same radius and inflow state. Compare the resulting acceleration with the Newtonian dependence on source mass and separation. [Physical free-fall comparison](../../docs/drain_driven_initial_acceleration.md#physical-free-fall-comparison-at-0pn) |

**Criterion for the eventual orbit test.** In the weak, slow, far-field regime, let $a_\parallel$ be the target’s acceleration component toward the companion after subtracting its isolated baseline. Let $M_{\rm source}$ be the companion’s independently determined inertial mass. For Newtonian gravity, $d^2 a_\parallel/M_{\rm source}$ must approach one positive value of $G$ for either choice of target across the intended family of throats, with transverse acceleration and drag subleading at Newtonian order. At fixed masses, doubling $d$ should give one quarter of the leading acceleration. This tests the load obtained by integrating flow and stress at the mouth. Source strength, passive response and inertia must be identified independently before comparing their ratios. The Newtonian coefficient is a comparison with the output. S21 takes the resulting flow requirements into its comparison with light. [V3 step plan, S15–S16 and S21](V3_STEP_PLAN.md) · [Gravity mass ratios](../../docs/cross_sector_research_burdens_and_compatibility_gates.md)

The fixed-radius exterior accepts different weak inward currents at the same distant density and phase. Two trial laws that coincide over the entire reference mouth predict different intake ratios when its radius doubles. The physical mouth still owes the normal-current response that selects $D(a)$ for this ambient state. [One-mouth calculation and limits](../../docs/drain_driven_initial_acceleration.md#does-one-isolated-mouth-select-its-intake-from-its-radius)

### Does the drain field reach the companion?

Bulk material settles and re-orders into the brane, providing a small, approximately uniform replenishment per unit brane area. On that background, a source drain must produce a disturbance that reaches a moving companion. The inverse-square requirement applies to the companion’s **net free acceleration**. Its material-speed profile remains open. The eventual momentum balance includes what the companion receives, what leaves through its sleeve and slab faces, and what accumulates in the surrounding material. [Documented replenishment](../../docs/oriented_throat_return_source_correction.md) · [Force comparison](../../docs/drain_driven_initial_acceleration.md#passing-the-measured-field-to-the-release-test)

Whether a charged throat locally thins the slab and thereby enhances replenishment is an open hypothesis for the brane-calculation stage, after the force requirements have been collected. The finite-slab and prescribed-intake examples illustrate how to distinguish a background flow, a change in total intake and a change in its angular distribution. [Finite-slab comparison](../../docs/drain_driven_initial_acceleration.md#measuring-the-flow-due-to-the-source-throat) · [Prescribed-intake exterior](../../docs/drain_driven_initial_acceleration.md#how-a-change-at-the-mouth-appears-in-the-exterior)

The earlier flat-slab calculation gives a benchmark at distances much greater than the slab thickness. It prescribed the geometry and assumed an admissible drain. Both its absorbing-cell and periodic-stack cases gave an inverse-square radial field. In its control where the field could spread into an unbounded fourth spatial direction, the power changed to inverse cube. These are results for the specified geometries and boundary conditions. The old calculation accounts for the uniform areal supply separately from its source-localized slab-boundary response. The native rates of re-ordering, transport and stress near a draining throat remain open. [Slab assumptions and range](../../software/stage1_solver/directives/pathA_29_brane_bulk_return.md) · [Flat-slab results and control](../../software/stage1_solver/reports/pathA_29_brane_bulk_return.md)

Specify the range of separations and orbital rates over which the moving solution has the required leading response. In S14a's separate spreading-slab probe, averaging at fixed radius before taking the far limit gave a different outer flux from taking the far limit first. Its slow-response probes also distinguish a long-time average from a low-frequency limit. For the orbit, the relevant question is what field reaches the companion **as it moves**. A steady profile can serve if the dynamics justify it on that range. Any cancellation, delay or change of radial power there must be carried through to the predicted acceleration. The selection requirement draws on S12's drain and return laws, S8/S22's material law, and the laws for the slab faces and interior whose supplier is still open. [S14a probes and open inputs](steps/S14a_2_slow_limit.md)

### What the sleeve and each slab face must report

The brane-side mouth integral gives the incoming load. A downstream calculation can supply the other terms in the following balance. An independent total force measurement at the held mouth can supply their combined effect.

Compare the target in the same ambient medium, with matched initial motion and drain and return laws. Allow the companion to change the actual throughput, surface shape and local flow. In each state, use a volume enclosing the target and excluding the companion. Label nonoverlapping portions of its boundary as the lateral surface, upper face, lower face and sleeve outlet. Use the same rule for drawing them in both states. Let $\Delta$ mean the paired value minus the isolated value, projected along the same axis toward the companion. [S15 field and region definitions](directives/S15_decision_list.md) · [Moving-region balance](steps/S16_moving_region.md)

- **Upper and lower faces.** The material balance already gives the signed number crossing each moving graph per projected area. Its relative crossing speed uses $u^w-\partial_t w_s-u^i\partial_i w_s$, with opposite outward orientations on the two faces. To obtain carried *momentum*, identify the transferred mass, its in-plane velocity and the appropriate area measure in each state. Where the isolated O2 premise applies, the carried component is the signed outward mass current times the local brane material velocity. A companion can change this product even at fixed total throughput. Mechanical traction contributes separately, including the in-plane projection of a normal stress on a sloping face. [Material face traces](steps/S14a_1_projection.md) · [O2 exchange convention and domain](directives/O2_SHARED_PHYSICS.md)
- **Sleeve outlet.** Supply the companion-induced change in the integrated axial momentum carried across this cut and in the stress transmitted across it. Identify where that momentum goes and which system takes any conversion reaction. The core and outlet data remain open. A scalar order-loss rate describes conversion in the same material and cannot determine the departing material’s axial momentum. [S12 source and domain](steps/S12_order_source.md) · [S15 core and face handoffs](directives/S15_decision_list.md)

**Consistency check.** Choose the positive axis toward the companion. With inward lateral current $I$, signed outward total currents $L_s$ on the other cuts, independently accounted volume force $B$, and momentum $P_{\rm env}$ assigned to the enclosed surroundings, the axial rate assigned to the throat must agree with

$$
\Delta I-\sum_{s\in\{+, -,\,\mathrm{sleeve}\}}\Delta L_s
+\Delta B-\Delta\dot P_{\rm env}.
$$

Project the momentum derivatives after taking them in the ambient frame. Each $L_s$ uses current relative to the moving cut and includes its stress contribution. A face exchange already represented in the projected force density must be counted once in this balance. Moving the lateral surface through a region where the laws apply should give the same throat result once storage and exchange in the intervening region are included. This gives a check independent of the Newtonian target. [Moving-region balance](steps/S16_moving_region.md) · [Projected exchange convention](directives/S15_decision_list.md)

**How to interpret the result.** Exchange can reinforce, reduce or cancel the companion-directed contribution from the exterior. Under a demonstrated positive scalar translational inertia in this limit, exact cancellation of the total at inverse-square order removes the leading Newtonian attraction. A net outward coefficient gives repulsion. An inward coefficient still has to give the same mass-normalized acceleration when source and target are exchanged. For a coupled translational response, apply that test to the free acceleration obtained from the full response. An unspecified current or inertia leaves the test open. [Passive and inertial requirements](steps/S16_moving_region.md) · [Source excerpts](_measurements/newtonian_drain_requirements_sources.txt)

### How much gravity a throat sources per unit of inertial mass

Two equal throats can reproduce an inverse-square dependence while leaving the mass dependence undetermined. In a physical target comparison the mouth radius, inflow profile and effective mass belong to the same draining state. Independent identification means reading the intake, source field and translational response of that state by separate measurements. Assigning extra bare mass to an unchanged drain tests a mechanically loaded object. An isolated throat needs three independently identifiable quantities in the same ambient medium:

- **Drain rate.** In the spherical S14a accounting, $J=J_{\rm face}+J_{\rm core}$ counts material leaving through the faces and inward through the core boundary. The record separately gives the signed loss of brane order, storage, and outer flux at a finite radius. Recovering the far-field volume intake $Q$ from these data requires the profile and the slow return dynamics. [S14a flux definitions](steps/S14a_1_projection.md) · [Slow-limit requirement](steps/S14a_2_slow_limit.md)
- **Active source strength.** Specify the surrounding gravity-sector field and its normalization, common to all throats. Read its leading far-field amplitude from the isolated source. If drain velocity supplies the whole leading field, its monopole amplitude can be read from $Q$. If stress and return also source that field, include their contributions. Calling an observed pair acceleration a source mass would fold the target’s response into the source measurement. [S15 source and normalization](V3_STEP_PLAN.md) · [Mass and force normalization](../../software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md)
- **Inertial mass.** Supply the translational inertia of the same draining state, including its bound fields and momentum exchange, independently of a gravitational pair fit. The dimensional candidate $\alpha_J\hbar J/c_\gamma^2$ leaves the relation open. Fitting $\alpha_J$ for each throat would use the mass the relation is meant to predict. C3 records this gap. [Mass-bridge verdict](../../software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md) · [C3](DEFECT_REGISTER.md)

**Two comparisons distinguish source scaling from free fall.** In the weak-field regime, keep one target fixed *in type and preparation* and compare its free acceleration toward two alternative sources at the same separation. The ratio of accelerations must equal the ratio of source inertial masses. If one source has twice the independently measured inertial mass, it must produce twice the leading acceleration. Then use one source and compare two freely moving target types at the same separation and initial velocity. They must have the same leading acceleration. In both comparisons the companion supplies the field, the drain and return laws stay the same, and the comparison uses the gravity-sector response. These tests use ratios, so a common choice of units or a single calibration of $G$ cannot fix a failed ratio. A successful ratio test still leaves the absolute value of $G$ to an independent prediction or an explicitly identified calibration. [Active, passive and inertial response](../../docs/cross_sector_research_burdens_and_compatibility_gates.md) · [Source excerpts](_measurements/newtonian_drain_requirements_sources.txt)

#### How drain rate and retained momentum must scale for free fall

Take two open, slowly moving mouths far apart in the same ambient brane. Use the weak, radial gravity response after subtracting the isolated target at the same initial velocity. Let $D_i>0$ be the **mass crossing mouth $i$ per unit time** and $M_i>0$ the translational inertia of that same radius and inflow profile, including the brane material that moves with it. Steady throughput can coexist with constant $M_i$. Let $Q_s$ be the positive *volume rate inferred from the source's exterior field*. For a steady spherical monopole in three brane dimensions at ambient mass density $\rho_\infty$, the same inward mass flux crosses every sphere in a source-free annulus. Dividing by its area and by $\rho_\infty$ gives $U_s=Q_s/(4\pi r^2)$ toward the source at leading order. Flux conservation alone fixes the angular mean. The local speed takes this form when the monopole dominates. [Exterior derivation](../../docs/held_drain_pair.md) · [Flux definitions and scope](steps/S14a_2_slow_limit.md)

At fixed target preparation, let $e_t$ be the signed net force driving translation per $D_t U$ in a weak incident flow $U$. The mouth balance counts inward carried momentum and traction, then subtracts momentum passed on or stored outside the translating mouth and its bound material. For a stationary preparation, this net force can be read from the change in the external force needed to hold the center, with the sign reversed. Independently, the response of that same draining mouth to a known weak force gives its inertia. Acceleration-dependent exchange enters that response once. Assuming a scalar linear response in the slow regime, negligible coupling between the two accelerations at this order, and a velocity monopole supplying the whole leading drive, the radial acceleration is

$$
 a_{t\leftarrow s}
   =\frac{e_t D_t}{M_t}\frac{Q_s}{4\pi r^2}+o(r^{-2}),\qquad
 \lim_{r\to\infty}\frac{r^2 a_{t\leftarrow s}}{M_s}
   =\frac{1}{4\pi}\frac{e_tD_t}{M_t}\frac{Q_s}{M_s}.
$$

The limit refers to the far-field regime of that exterior. The right side comes from flux, momentum transfer and inertia. Comparing it with a single positive $G$ gives two independent requirements. Holding the source fixed, every target must have the same $e_tD_t/M_t$. Holding the target fixed, every source must have the same $Q_s/M_s$. Together these conditions suffice for the leading radial Newtonian law under the stated assumptions. A leading contribution that depends on both profiles beyond the incident speed requires checking the full pair coefficient. [Momentum balance](../../docs/drain_driven_initial_acceleration.md#one-isolated-mouth-at-the-instant-a-sideways-force-is-applied) · [Mass and profile dependencies](../../software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md)

**When all mouth intake is supplied by the exterior monopole**, $Q_i=D_i/\rho_\infty$. For a family of throats that can each act as source and target, the source requirement then fixes a common $\lambda=D_i/M_i$. The target requirement in turn fixes a common positive $e_i=e$, giving $G=e\lambda^2/(4\pi\rho_\infty)$. For example, doubling the intake would require twice the translational inertia and the same retained force per unit incident momentum rate. Comparing targets at only one source tests the product $e_iD_i/M_i$ and permits compensating changes in the two factors. [Symbolic ratio checks](_measurements/free_fall_mass_scaling.txt)

**What is established.** At eight mouth radii a prescribed, stationary radius-one target gave successive load ratios of $1.99999196$ and $1.99995211$ as its intake doubled twice with the source fixed. This measures the carried momentum and stress on the brane side. The missing physical relation is how the **net** force per incident speed and the independently determined translational inertia scale together across actual draining states. The existing solver prescribes the centers. The normalized inertia and the net force that moves a free center remain open, as does source strength per unit inertia. A finer grid for the same held problem cannot supply that relation. [Held target-rate measurements](../../docs/drain_driven_initial_acceleration.md#changing-the-targets-rate-at-fixed-companion-preparation) · [Free-mouth response](../../docs/drain_driven_initial_acceleration.md#what-would-turn-the-measured-held-force-into-a-release-acceleration)

The archived radius comparison takes its mean speed from a radius-one mouth whose rate was set to $0.25$. Under that prescribed speed profile, the measured brane reaction per unit intake increases as the radius grows. Applying the measured ratio to physical mouths requires an independently justified radius-to-intake law and compatible mouth profiles. If the full inertia per unit intake is common, the remainder of the translational response per unit intake must decrease by the same amount. The one-body force law identifies the unmeasured total response to a known sideways force on that same opening. [Mouth inertia comparison and its assumptions](../../docs/drain_driven_initial_acceleration.md#can-the-measured-brane-reaction-supply-the-inertia-of-a-draining-mouth)

For unequal mouths, the ideal steady exterior delivers equal and opposite momentum currents including pressure. In the physical balance, internal transfers cancel when the two regions are combined. The signed interaction-induced outward momentum flux and rate of storage outside the translating mouths, less any independent applied force, must have zero inverse-square coefficient for the two mouths to exert equal and opposite forces at that order. A nonzero inward drive on each is a separate requirement. [Unequal-mouth calculation and signed balance](../../docs/held_drain_pair.md#do-unequal-mouths-exert-equal-and-opposite-0pn-forces)

The archived radius-1 and radius-1.5 benchmark now has signed brane-side loads for both mouths with each isolated at the same position. Both point inward. The archives lack the physical support forces and the receiving-side momentum transfer needed to infer them. [Matched unequal-mouth readout](../../docs/held_drain_pair.md#what-can-be-read-for-the-two-held-unequal-mouths)

**Open 0PN requirement for a held mouth.** At each radius, keep the same ambient state and physical inward-drain law. Let $H_x(a,U)$ be the total external force required to hold the geometric center still in weak sideways incident flow $U$. The missing response is $e_tD_t=-\partial_U H_x(a,0)$ when the companion enters at leading order only through its local incident velocity. An independently derived force law at the opening or signed total holding-force readings at $\pm U$ and $\pm U/2$ would test it. In the periodic benchmark the incoming momentum rate is fixed by the feed. For a single neutral carrier of the same mass on both sides, a common differentiable phase trace gives equal tangential velocities. With equal pointwise normal number current, the carried sideways momentum passes across the opening. The projected parent confinement term gives a distributed wall load, and the wall tension can transmit a separate endpoint traction whose response remains undetermined. [Incoming momentum](../../docs/held_drain_pair.md#what-does-the-signed-one-mouth-flow-test-measure) · [Transmission at the opening](../../docs/held_drain_pair.md#what-passes-on-if-the-carrier-phase-is-continuous-at-the-mouth) · [Endpoint traction](../../docs/held_drain_pair.md#does-the-parent-confinement-term-give-a-sideways-force-at-the-mouth) · [Two-density calibration and missing receiving-side measurement](../../docs/_measurements/mouth_stress_response/calibration.md)

#### What a later throat model must reproduce at the mouth

**Open requirement, 2026-11-17.** For a mouth of radius $a$ draining at its physical rate $D(a)$, the brane calculation can use two independently determined responses. Let $C_a=-H'_a(0)$ be the net pull per unit weak incident brane-material speed $U$, and let $\chi_a$ be the acceleration per unit weak applied force on the same mouth in resting brane. For smooth release from rest, the predicted acceleration toward the source is $\chi_a C_a U$ to first order when $U$ supplies the leading drive and coupling between the two accelerations is subleading. If a companion also drives an independent bulk flow at leading order, measure its force response separately. If the response has memory on the time scale of interest, retain the measured force and displacement histories. [Separate brane and bulk drives](../../software/stage1_solver/reports/mouth_tow_and_companion.md) · [Mouth balance and release identity](../../docs/drain_driven_initial_acceleration.md#one-mouth-translation-under-an-applied-force-or-an-established-brane-flow)

**What this requires of material entering the throat.** In a steady comparison with other applied loads unchanged, write $C_{\rm in}$ for the incoming sideways momentum per unit time and incident speed, including stress. The signed net rate passed on per unit speed, including stresses on the remaining cuts, must then be $C_{\rm in}-C_a$. Under the adopted conversion premise the material carries its local in-plane brane velocity across conversion. Any change in that momentum must have a corresponding reaction. If carriage matches across the opening, a nonzero pull must come from transmitted stress or another identified momentum exchange. A later interior theory must reproduce the net transfer and push response at the opening. The 0PN balance allows different internal distributions with those same observables. [Local conversion rule](directives/O2_SHARED_PHYSICS.md#L219) · [Transfer requirement and its algebraic check](../../docs/_measurements/paired_gravity_transfer.txt)

**Test without an interior model.** Hold only the center. Apply opposite weak brane flows at full and half strength with the distant bulk at rest and the local bulk free to respond. Measure the total holder force. In resting brane, free the same mouth and apply opposite known pushes at both strengths to determine $\chi_a$. Use those two measurements to predict a separate release in established flow. Repeat for a different radius with its own $D(a)$ and the same source. Newtonian free fall requires the same positive $\chi_a C_a$ for both targets. On the velocity-monopole branch, the independently measured $Q_s/M_s$ must also be common across sources to give a single $G$. The physical holder force and push response remain open. The prescribed uniform receiver inlet fixes the sideways carrier velocity at the opening and cannot supply this calibration. [Force probes and release test](../../docs/drain_driven_initial_acceleration.md#what-would-turn-the-measured-held-force-into-a-release-acceleration) · [Inlet comparison](../../software/stage1_solver/reports/mouth_inlet_trace.md)

**Conversion and the ideal drag cancellation.** Keeping the local in-plane material velocity across conversion is compatible with the previously calculated ideal exterior balance at the level of integrated steady momentum. Carriage and inter-medium stress cancel across the internal cut. Equal intake/discharge mass rate $D$ then gives total reaction $D(U_{\rm br}-U_{\rm bulk})$ only when the mouth takes all remaining mechanical reactions. Signed momentum transferred to separately supported surroundings must be subtracted. Normal pressure on a curved face can transmit an in-plane force even without bulk shear; its amplitude and transmission to the geometric opening remain open. The physical boundary condition must determine those responses, retain the pressure/shape endpoint contribution, and identify every reaction system. The symbolic check supplies no matched flow or inertia. [Balances, live pressure moment and conditional radius scaling](../../software/stage1_solver/reports/mouth_conversion_drag.md)

**Energy condition for that ideal branch.** Let $H_f=e_f+p_{\infty,f}/\rho_f$ and $H_b=e_b+p_{\infty,b}/\rho_b$ denote the specific internal energy plus pressure work on the two sides. In steady translation, the energy arriving through the ideal exterior cuts minus the mechanical power $D(U_{\rm br}-U_{\rm bulk})V$ is $D[H_f-H_b+((U_{\rm br}-V)^2-(U_{\rm bulk}-V)^2)/2]$. The physical interface must account for this remainder through storage and outgoing energy, with any additional supply naming its source and power budget. The conversion momentum premise fixes neither $H_f-H_b$ nor that energy relation. Order energy already included in these heads must be counted once. Time-dependent exterior pressure also carries the energy needed to change the surrounding flow; dropping it violates the finite-shell energy balance. The ideal calculation specifies a requirement, with physical energy supply and the isolated-pair budget still open. [Energy calculation and checks](../../software/stage1_solver/reports/mouth_energy_balance.md) · [Adopted energy obligations](directives/O2_SHARED_PHYSICS.md#6-energy-input-and-forcepower-pairing)

**Pair energy in the ideal exterior.** With fixed symbolic rates $D_i=D(a_i)$, the matched common-domain interaction kinetic energy is $D_1D_2/(4\pi\rho_3d)+D_1D_2/(4\pi^2\rho_4d^2)$. It increases during approach at the same rate as the ideal attractive load does mechanical work. Including unsteady pressure and the transport of both moving cuts, the two inner regions together must therefore make available an interaction increment of twice that mechanical work. This is a conditional requirement on storage, conversion power and other energy outputs relative to isolated controls, not a measured power source or twice the baseline drain power. The physical boundary law must account for it if this fixed-rate monopole branch applies. Unknown interior energy and response dipoles were not assigned; comparisons using different cut domains require their own self-energy corrections. [Two-drain energy test](../../software/stage1_solver/reports/two_drain_energy.md)

**What the order equations supply for that power.** The historical bare potential $a_B\chi^2(1-\chi)^2$ has equal-depth homogeneous phases. A negative signed conversion rate can either raise or lower its energy along the transition; it establishes no available mechanical power. The shear gate and interface-gradient term contain possible stored energy, but their live response and energy destination remain open, as do the core/mixing terms. The recorded $\int\mu_\chi D_t\chi\,d^4X$ is partial order work, not the whole moving-region balance. At fixed $D_i=D(a_i)$, an effective net energy yield $\ell_i$ at the mouth would have to satisfy $D_i\Delta\ell_i=P_{\rm mech}$ for each mouth on the ideal branch above. A fixed yield at unchanged rate supplies no interaction increment; a change in yield, stored energy or another output must account for it. An independently determined boundary energy-response law could supply this datum without resolving the interior. No such law or physical energy source has been selected. [Order-conversion energy assessment](../../software/stage1_solver/reports/order_conversion_energy.md)

**What the known carrier action adds.** Its neutral-sector energy current is $\mathbf J_E=-(\hbar^2/m)\operatorname{Re}(\partial_t\psi^*\nabla\psi)$, with $\partial_t\mathcal E+\nabla\cdot\mathbf J_E=n\partial_tV_{\rm conf}$ and $\mathcal E=\hbar^2|\nabla\psi|^2/(2m)+V_{\rm conf}n+U(n)$. A moving cut uses $\mathbf J_E-\mathcal E\mathbf W$. This determines a carrier energy readout from complete boundary fields; it does not select the physical receiving trace or the live brane's additional energy partners. Local stress work uses its application velocity and cannot in general be inferred from the resultant center force alone. A translated confinement profile contributes reciprocal profile work, $\int n\partial_tV_{\rm conf}=-\mathbf F_{\rm profile}\cdot\mathbf V$, which cannot also be counted as independent conversion power. The held feed's prescribed phase frequency gives its carrier energy input, not the missing moving-mouth energy law. The two-drain power condition remains unverified for physical mouths. [Energy current and boundary-data assessment](../../software/stage1_solver/reports/mouth_stress_energy.md)

**Stationary-mouth test for a future solver.** Hold mouth 1 still. In the brane and receiving bulk, use disjoint spheres centered on the respective mouths. Keep those centers fixed relative to each mouth when repeating the test at a second radius. For nested stationary surfaces, their outward flux difference includes the change in energy stored between them and any intervening sources. In the local check, displacing the measuring sphere changed the flux slope at both radii. [Independent local review](../../docs/_measurements/stationary_mouth_energy_review/local.md) · [Measured controls](../../docs/_measurements/stationary_mouth_energy_review/local_evidence.json)

Measure the interaction energy leaving the stationary mouth's region during approach and recession. Divide the difference by twice the signed approach speed. In the slow-motion limit, the ideal prediction for this slope is the attraction magnitude measured independently from the total static holding reaction at the same separation. Compare the full- and half-speed slopes. If they disagree beyond uncertainty, reduce the speed or bound the finite-speed error before judging the limiting slope. The physical comparison requires a matched receiving-bulk energy current and the total force on the holder. [Test definition and isolated controls](../../software/stage1_solver/reports/stationary_mouth_energy_test.md) · [Local review and user disposition](../../docs/_measurements/stationary_mouth_energy_review/round2/status.md)

**Open 0PN requirement if a trapped shear wave participates in free fall.** For the same radius-selected inward branch \(D(a)\), determine whether the localized mode follows the geometric center under a weak push, for how long it remains confined, and how much it contributes to the total translational response at the gravity timescale. If it also transmits the companion's pull, determine the signed change in the cycle-averaged traction doing work on a translation. Identify the exterior momentum that supplies that load and compare with the independent total center-holding reaction. Radial pressure alone has zero sideways resultant for an isolated isotropic mouth. In a conditional quadratic-pressure example, interference of the standing strain with an in-phase, flow-induced dipole gives a load linear in incident speed. A dipole in quadrature gives zero at that order. The physical mode mixing, net traction, inertia, radius scaling and power source remain open. The mode could supply inertia even if another stress supplies the pull. This is a test for the later throat model, with radial support and the interior profile still deferred. [Conditional derivation and signed tests](../../software/stage1_solver/reports/trapped_wave_mouth_response.md) · [Candidate radial support](V3_STEP_PLAN.md#L917)

For the push measurement, save independent time histories of the geometric center and the localized shear pattern in the lab frame. A finite restoring response can give an early lag followed by coincident centers. One coincident snapshot or exact half-strength scaling cannot establish that the wave moves with the mouth on the gravity timescale. A reduced calculation gives a displacement-only test for the assumed restoring response and states when the two contributions to inertia can be distinguished. [Four-push calculation and measurement](../../software/stage1_solver/reports/trapped_wave_push.md)

For a quantitative tracking bound, specify the companion drive band and the allowed mismatch. In the conservative rim-driven model with relative stiffness \(k\) and trapped-wave inertia \(M_w\), tracking each prepared harmonic acceleration through angular frequency \(\nu_g\) to fractional accuracy \(\epsilon\) requires \(k/M_w\ge(1+1/\epsilon)\nu_g^2\). Centered radial support permits arbitrarily small positive relative stiffness. The native bound remains open. The total force and the geometric rim motion in one weak-push experiment could constrain this ratio under the stated response law. Abrupt release also requires control of the excited relative oscillation. [Tracking and release bounds, including the mouth-only test](../../software/stage1_solver/reports/wave_follow_bound.md)

The three-reading fit predicts the reserved fourth reading exactly on the archived single-mode example. The current draining-mouth archive supplies zero matched total-force and independently moving-center readings. The physical rate bound and its fourth check await that one-mouth response. [Fourth-reading result and inputs needed](../../software/stage1_solver/reports/wave_fourth_reading.md)

An effective boundary law may supply the missing sideways reaction without resolving the throat. It must determine the net force transmitted at the cut as the rim moves on the prepared inward branch. In the passive wall reduction, inertia beyond the cut contributes to the slow force per acceleration with weight equal to the squared static displacement of the cut per unit rim translation. The present equations leave that contribution open even for a reflecting end. On the receiving side, even a specified pointwise normal current and common far response can admit different integrated sideways carriage and stress under the same slow acceleration. [Wall force law](../../software/stage1_solver/reports/rim_end_load.md) · [Receiving-cut calculation and missing traces](../../software/stage1_solver/reports/cut_acceleration_flux.md)

For opposite weak rim accelerations, record the total downstream sideways momentum current and the change in lab momentum between that cut and the rim on the same draining state. Include the wall and any trapped wave in both readings where present, and close the impulse balance with the feeding face, side boundary and external supports. The signed receiver calculation resolves the neutral carrier contribution. The wall stress transmitted across the downstream cut and the wall and wave momentum inside that volume still require their own local response. [Signed cut and storage measurement](../../software/stage1_solver/reports/rim_cut_inventory.md)

The projected wall equation reads the cut load from the wall tension and displacement slope. Its stored momentum uses the wall velocity throughout the short region. For a single transverse standing packet on an ordered patch, the translation momentum uses the mode-weighted envelope velocity, while the cut current depends on the local gradient of the sideways shift. Identify these quantities on the prepared drain or supply an equivalent calibrated response. The axial shape and the rim acceleration alone do not fix the physical split. [Wall and trapped-wave calculation](../../software/stage1_solver/reports/wall_wave_cut.md)

The archived radius-one test at inward rate 0.25 has both acceleration signs, both half-strengths and the undriven control. It records the commanded rim position and feeding-side momentum rates. Its grid ends at the rim. To measure the split, use a receiving region on this same draining state with an independently responding wall, a normalized trapped shear field, and the continuation's load at the outer boundary. Record the wall geometry and phase-resolved shear field through the near region and on both sides of the cut, in a common lab frame. Compare the half-difference of the opposite drives at matched times and wave phase or over matched cycles. The half-strength pair tests linear scaling. Integrate the momentum density inside and the traction across the cut separately, using the calibrated inertia and stress laws. Check their impulse balance with the other measured boundary and applied loads. [Signed runs and saved fields](../../docs/_measurements/signed_mouth_wall_wave_coverage.txt) · [Local observables](../../software/stage1_solver/reports/wall_wave_cut.md)

The short receiving-segment solver can carry the signed rim drives forward with independently moving wall and wave coordinates and measure momentum inside and across an interior cut. The numerical benchmark closes this projected balance. Applying it to the physical drain requires the locally normalized standing mode, wall and shear coefficients, wall–wave coupling, incremental carrier loading and a tested downstream force response for this same state. Use a sampling window long enough for a disturbance to reach the cut at the measured speeds. [Implementation and signed benchmark](../../software/stage1_solver/reports/signed_mouth_collar.md)

The amplitude-restoring surface energy supplies an absolute force on the prescribed intake profile. For the radius-one, intake-0.25 state it can be evaluated during the signed drives, including the density response. Its use as a receiving-wall source requires an independent transfer relation. If the rim motion alone is prescribed, a common rescaling of all projected inertias, tensions, coupling and terminal load preserves the unloaded segment's motion while rescaling the momentum budget. A receiving-side traction or equivalent independent normalization is needed to resolve that ambiguity. [Same-mouth force and calibration limit](../../software/stage1_solver/reports/mouth_feed_profile_force.md)

For the slow free-fall response, the unknown normalization can survive in the acceleration itself. In a passive lateral field attached to the mouth, scaling its inertia and tension together preserves the static held load and the field's motion for a prescribed rim trajectory. Its momentum for that trajectory scales with the normalization. The physical mouth must select the ratio of total pull per incident speed to co-moving inertia, including any pull transmitted through the trapped wave. Compare that ratio across physical radii on their inward branches. [0PN cancellation test and its assumptions](../../software/stage1_solver/reports/free_fall_response_ratio.md)

#### Can a mouth respond to the companion's flow without drag from its own motion?

**Working 0PN requirement.** Take the same open mouth of radius $a$ on its inward-draining branch $D(a)$. Choose one brane axis. In the frame of the distant ambient material, let $U$ be the companion's incident brane velocity at the target location, evaluated with the target absent, and $V$ the velocity of the target's geometric center. Let $F_x$ be the total drive on that center, positive along the chosen axis. At constant speed and steady translating inertia it is the negative of the external force required to maintain the prescribed motion. Both tests concern translation of the opening.

| Independent response | Preparation and reading |
| --- | --- |
| $A(a)=\partial_U F_x(0,0)$ | Hold the center still. Place the weak, distant companion on either side and subtract the target alone at the same position. Use the component of incident flow along the fixed axis. Interpreting the result as a derivative in $U$ assumes the companion enters at leading order through its local velocity. A uniformly driven feed requires independent matching to that preparation. |
| $B(a)=\partial_V F_x(0,0)$ | Translate the isolated mouth at small constant speeds of opposite sign through the same resting ambient state. Read the signed total external force at the same position after preparation transients, and reverse its sign. The far boundary, radius and physical intake law must match. |

Use the odd difference divided by twice the speed in each test and repeat at half the incident or translation speed. Material drains inward throughout. For smooth, isotropic response in the slow limit the two slopes give $F_x(U,V)-F_x(0,0)=A U+B V+o(|U|+|V|)$ when other leading drives are absent. If the mouth senses only $U-V$, differentiation gives $B=-A$. Then $A>0$ would give an isolated moving target a force opposite its velocity. For a constant-mass body this force removes mechanical energy at rate $A|\mathbf V|^2$ to leading order in speed. A drag of this size surviving at Newtonian order would have to be accounted for in the orbital energy and angular momentum. The velocity-driven Newtonian mechanism needs positive radial response and subleading **net** tangential acceleration. Subtracting an isolated control in the radial comparison cannot remove a physical torque from the pair.

**What the equations decide.** A change of observer transforms the mouth, incident material and distant reservoir together. It preserves both velocities measured relative to that reservoir, so it supplies no identity between $A$ and $B$. The held periodic feed sets a phase gradient. The isolated drift calculation translates the intake profile with the outer bath fixed. Writing the incoming carried momentum per unit time in the mouth frame subtracts $D(a)V$ from the ambient-frame reading at linear order. That shift also applies to material transmitted through the opening and cancels between the two carried contributions for equal crossing mass rates. It supplies no value for the total drive. The actual force also involves the stress difference and any independent rim force. The physical surface law and the tangential exchange involving the receiving material are open in the body equations. The archived isolated exterior load changed sign on successive lateral refinements. The present boundary data determine neither physical slope nor the relation $A+B=0$. [Source commands and readings](_measurements/mouth_flow_and_motion_sources.txt) · [Transmission and rim forces](../../docs/held_drain_pair.md#what-passes-on-if-the-carrier-phase-is-continuous-at-the-mouth) · [Isolated drift balance](../../docs/drain_driven_initial_acceleration.md#drift-under-lateral-refinement-at-fixed-thickness-resolution)

**Input for a future simulation.** Supply one mouth-level momentum-transfer law for the net force on the geometric center, or independent total-force readings in both preparations to calibrate its linear response. Record the inertia of that same draining state separately. The existing radius-to-intake, source-strength and mass-scaling requirements above then connect the force to the predicted free acceleration. The later trajectory can test work and torque with those inputs held fixed.

#### Does the existing tangential slip law supply the mouth force?

The one-body model contains a candidate dissipative force at the boundary of the throat wall (option E5). It requires the normal velocities of material and wall to agree. Applying that condition over the **open aperture** would give zero relative normal current there, contrary to the required $D(a)>0$. It can describe a wall beside an opening that continues to drain. The model also lists a permeable option (E3), which leaves the velocity of the crossing material free relative to the moving pattern. Neither option specifies the sideways force retained by the draining aperture. The boundary study left the physical choice unresolved. [Boundary definitions](../../software/em_charge_attribute/u1_body_dynamics_inputs.yaml#L194) · [Boundary study](../../software/em_charge_attribute/reports/u2_boundary_adjudication_artifacts/stage_1_production_v12/production_summary.md)

For one tangential component on the wall, let $v$ and $V_\parallel$ be the local material and wall velocities. The slip coefficient $\gamma$ and the exterior traction coefficient $K$ are positive open inputs in the one-body calculation's normalization. With the scalar exterior equation written $K v+\gamma(v-V_\parallel)=T$, the archived translation problem has $T=0$. Here $T$ denotes a possible *applied exterior traction in that equation*. A companion-driven calculation would have to supply it. The Rayleigh density $\gamma(v-V_\parallel)^2/2$ gives the wall traction

$$
 f_\parallel=\gamma(v-V_\parallel)
             =\frac{\gamma}{K+\gamma}(T-K V_\parallel).
$$

On the archived branch the wall therefore feels a traction opposite its imposed motion. The combined mechanical power on wall and material from this slip term is $-\gamma(v-V_\parallel)^2$. For a held wall, the traction per unit projected incident speed is instead $\gamma\,\partial_{U_\parallel}T/(K+\gamma)$ at rest. The archived boundary equation supplies no mapping from the distant companion's flow to that local traction. A driven family that retains $T(0,V_\parallel)=0$ for all isolated translations has $\partial_{V_\parallel} T=0$ there. More generally, if the omitted forcing responds to both velocities, the two local slopes are equal and opposite precisely when $\partial_{U_\parallel}T+\partial_{V_\parallel}T=K$ at fixed coefficients. These are component tractions on the wall. Inferring a force on the geometric center would additionally require how this traction is transmitted to the rim and the remaining momentum balance at the opening. [Exact response, inputs and checks](_measurements/mouth_slip_traction.txt) · [Force transmitted to the rim](../../docs/held_drain_pair.md#does-the-parent-confinement-term-give-a-sideways-force-at-the-mouth)

**Requirement at the opening.** For the same radius-selected inward flow, a mouth law must give the net sideways reaction to independently varying the incident brane flow and the center velocity. At the opening take the difference of the two **total** momentum currents, each including the material crossing relative to the moving cut and its transmitted stress. Add independent rim forces and account once for momentum stored in or exchanged with material outside the translating body. Equal mass rates can still carry different tangential momenta. Their carried difference vanishes on the stationary, equal-current, single-carrier branch with a common differentiable tangential phase. An independent total mechanical reaction in each preparation would measure the combined effect. A wall-slip coefficient alone leaves both whole-mouth responses open.

#### What do two specified radii predict at release from rest?

Prepare two open, inward-draining mouths at rest a distance $d$ apart in the same uniform distant brane. Allow the held pair to settle. Remove both supports without an impulse. Preserve positions, velocities, independent field initial data, preparation history, and the intake and reservoir laws. Solve anew for the accelerations, field time derivatives and constraint stresses at release. Unsteady pressure can change immediately. Use one signed laboratory axis pointing from mouth 1 to mouth 2. The inputs are the radii $a_1,a_2$, separation and ambient state. Write $\alpha_i$ for the signed initial acceleration of geometric center $i$ in the gravity-sector comparison.

| Link in the prediction | Equation or independent reading | Status for physical mouths |
| --- | --- | --- |
| Radius to inward rate | Apply the same $D(a)$ to each isolated drain. The inward material current and its physical boundary law must agree. | Open. The archived boundary condition admits different prescribed rates at fixed radius. [Radius and intake](../../docs/drain_driven_initial_acceleration.md#does-one-isolated-mouth-select-its-intake-from-its-radius) |
| Source to incident brane flow | Infer the positive inward volume rate $Q_i$ from the source-alone current through spheres in a source-free far zone. For a steady dominant spherical monopole at uniform density, $U_{2\to1}=+Q_2/(4\pi d^2)$ and $U_{1\to2}=-Q_1/(4\pi d^2)$ to leading order. | Conditional exterior prediction. The relation between physical mouth intake and exterior strength still needs the intervening material budget. [Source field](../../docs/drain_driven_initial_acceleration.md#measuring-the-flow-due-to-the-source-throat) |
| Incident flow to force retained | In a steady held pair read the two signed total support forces $H_i$, or integrate the net transfer to each translating opening, including carried momentum, stresses and rim reaction. Use the total pair supports for the absolute release. Subtract matched isolated controls separately to identify the interaction. If other forces are present, their contributions must be identified to attribute the interaction to gravity. | Open. The prescribed-center archives give brane-side loads on both unequal mouths. They supply no total support reaction. [Unequal-mouth readings](../../docs/held_drain_pair.md#what-can-be-read-for-the-two-held-unequal-mouths) |
| Force to acceleration | At the same prepared state determine $\chi_{ij}=\partial\alpha_i/\partial f_j$ for independent weak forces along the lab axis on each center, with both free. In linear response to removing the supports, $\alpha_i(0^+)=-\sum_j\chi_{ij}H_j$ to first order in the released support forces. This assumes a smooth instantaneous force-to-acceleration response and zero acceleration with forces $H_j$ applied in the same state. Measure the onset response, including any acceleration-dependent brane load once. A steady drift or slow acceleration measurement alone cannot fix this response. | Open. The physical inertia and cross-response are unmeasured. The separate finite-companion benchmark demonstrates this readout for specified boundaries. [Pair holding and probes](../../docs/drain_driven_initial_acceleration.md#holding-force-and-small-force-release-test-with-a-finite-companion) · [Physical probe status](../../docs/drain_driven_initial_acceleration.md#what-would-turn-the-measured-held-force-into-a-release-acceleration) |

In the far zone, **if** the isolated gravitational drive vanishes by symmetry, the local incident velocity supplies the whole leading interaction, $\chi_{ii}$ approaches $1/M_i$, and cross-acceleration is subleading, the signed predictions reduce to $\alpha_1=+A(a_1)Q_2/(4\pi M_1 d^2)$ and $\alpha_2=-A(a_2)Q_1/(4\pi M_2 d^2)$ at leading order. The drift contribution $B V$ vanishes at the prepared instant. A release from rest cannot settle the later drag question. Only after the physical inputs are fixed can $d^2\alpha_1/M_2$ and $-d^2\alpha_2/M_1$ be compared across radii and separation to one positive $G$. With the present inputs, neither initial acceleration has a physical numerical prediction. [Source and target scaling](#how-drain-rate-and-retained-momentum-must-scale-for-free-fall)

### What must persist once both throats are orbiting

The eventual comparison uses two freely moving drains in a uniform ambient medium. Give them initial positions and velocities with nonzero orbital angular momentum. Their centers must then move according to the drain, brane and bulk equations. Use the independently determined inertial masses $M_1,M_2$ and the coupling from the mass-scaling comparison above. The Newtonian relations below are tests of the resulting trajectories. The earlier orbital solver supplies $G$ in its Poisson equation and serves as a reference for such trajectories. [Earlier solver](../1pn_orbital_dynamics/scripts/verify_newtonian.py) · [S15 moving-drain scope](directives/S15_decision_list.md)

**Motion of the pair.** Write $\mathbf r=\mathbf X_2-\mathbf X_1$ and $r=|\mathbf r|$. In the isolated, constant-mass Newtonian limit, the relative acceleration must approach $-G(M_1+M_2)\mathbf r/r^3$, and the mass-weighted center must move uniformly. Track both in the original ambient frame. A matched isolated moving throat can identify drag already present in that frame. Its effect counts in the orbital trajectory even if it cancels in a pair-minus-isolated diagnostic. Net momentum carried into the bulk must be included when interpreting any recoil. [One-body momentum and radiation channels](../../docs/em_u1_body_definition.md) · [Moving-region momentum accounting](steps/S16_moving_region.md)

**Period and orbital losses.** For a bound pair, measure the period $T$ and the semimajor axis $a$ of the *relative* trajectory. The Newtonian target is $T^2/a^3=4\pi^2/[G(M_1+M_2)]$. A period alone can conceal a slowly shrinking orbit. Track the specific orbital angular momentum $\mathbf h=\mathbf r\times\dot{\mathbf r}$ and the specific energy $\varepsilon=|\dot{\mathbf r}|^2/2-G(M_1+M_2)/r$ over successive orbits. Both must be constant at leading Newtonian order. These diagnostics use the same masses and coupling throughout. A sustained leading loss of energy or angular momentum would fail this comparison even with the right inward acceleration at release. [Keplerian reference dynamics](../1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex)

**Open input.** The moving-drain model owes the change in source strength and inertia over an orbital period and the net work and torque exchanged with the surrounding medium. Steady throughput may coexist with a steady body mass. If a computed trajectory drifts, compare the drift with storage and flux through the sleeve and faces, the return flow and any radiation, and check convergence before assigning it to the medium. The relevant intake, return and radiation responses are open. Whether a smaller departure is post-Newtonian belongs to the later velocity-dependent work. [One-body response status](../../software/em_charge_attribute/reports/u1_body_dynamics.md) · [S15 handoff to S19](directives/S15_decision_list.md) · [Source excerpts](_measurements/newtonian_drain_requirements_sources.txt)

**Attribution in a simulation with charge.** The model also permits an orientation-dependent static interaction with the same inverse-square falloff under Q1's postulated closure. A trajectory for a single pair of charged throats therefore measures their combined interaction. To attribute an acceleration to gravity, track the contributing momentum currents and stresses. Where the same ambient medium supports matched throats pointing to either side of the brane, compare orientations at matched masses and drain strengths, with material flowing inward in every case. An orientation average retains *all* orientation-independent contributions. The existing one-sided slab supports only the $+w$ orientation. [Static puncture result and scope](V3_STEP_PLAN.md#q2-earned-within-g0) · [Body and ambient symmetry](../../docs/em_u1_body_definition.md) · [Two-body force decomposition](../../docs/em_phaseC_force_decomposition.md)

### Can gravity and light use the same radial flow?

| In the isolated spherical exterior | Radial speed in the recorded calculation |
| --- | --- |
| Compact three-dimensional drain benchmark, with steady prescribed volume intake $Q>0$, constant ambient density and no exterior leakage | $\lvert v_r\rvert=Q/(4\pi r^2)$. [Reduced drain profile](../../software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md) |
| S9b's retained nonzero flow-only optical solution, with fixed local shear-wave speed, zero tilt and its stated velocity counting | $\lvert V\rvert=\sqrt{2GM/r}$ for an independently determined orbital $GM$, in the specified lab clock and speed units. [Light-bending record](steps/S9b_light_bending.md) · [Profile comparison](steps/S9b_flow_oracle_comparison.md) |

If both speeds represent the **same material velocity** around the same source in the same radial coordinate, clock and far-zone interval, their different radial powers preclude a match over that interval. Adjusting an amplitude can match a single radius. The Newtonian orbit criterion above constrains the acceleration of a throat and leaves the material-speed profile open. S9b's more general optical conditions also admit changes in local shear-wave speed and brane tilt. [Drain-force scope](directives/S15_decision_list.md) · [Optical profile comparison](steps/S9b_flow_oracle_comparison.md)

For S21, identify how the finite-radius projected material velocity relates to the brane velocity that transports light, including the density, face exchange and surface motion used in that identification. Carry the resulting profiles and independent orbital $GM$ to **both** light deflection and the logarithmic slope of radar delay over the common far zone, on S9b's recorded domains. The compatibility question stays open under `R-S21-01`. [Optical requirements](#r-s21-01--conditional-deflection-and-radar-slope-profile-requirements) · [S21 remit](V3_STEP_PLAN.md) · [Source excerpts](_measurements/newtonian_drain_requirements_sources.txt)

<a id="joint-requirements-for-0pn-gravity-and-static-electric-charge"></a>
<a id="joint-requirements-for-0pn-gravity-static-electric-charge-and-light"></a>

## Joint requirements for one brane and bulk medium

**Material class selected by the user, 22 November 2026.** Require a supersolid feeding brane with ordered transverse shear and frictionless coherent flow. On conversion into the fourth dimension require a disordered coherent superfluid with zero static shear modulus and zero shear viscosity. The conversion obeys mass and momentum balance and transmits pressure. Future microscopic laws for the bulk must realize this phase. Record 5 supplies a historical demonstration of shear relaxation in classical thermostatted particles and is retired as a physical receiver. [Model premises and 0PN tests](../../docs/zero_pn_free_fall_roadmap.md#material-assumed-throughout)

**Working synthesis for light, gravity, charge, magnetism and gravitomagnetism.** The ideal draining exterior delivers an inverse-square momentum load toward a companion. The localized Maxwell sector gives a long-range electric field once supplied with a charge. The assumed shear action on a uniform three-dimensional brane supports two transverse polarizations. Moving charge and moving mass test those same mouth and material laws in motion. [0PN calculation](../../docs/zero_pn_gravity.md) · [Electric charge checkpoint](../../docs/em_charge_attribute_requirements.md#november-checkpoint--what-a-charge-law-must-explain) · [Light-sector closeout](#light-closing-pass--the-oracle-comparison-and-the-extra-wave--2026-10-19-accepted-e070c7f3)

### What the common-medium search has established

We want one set of local material laws to predict what a distant probe sees and how a free mouth moves. The two-radius test gave an early constraint. On the drain-driven gravity branch, source strength must scale with independently measured mass. Elementary charge must keep the same magnitude as mouth radius and stored energy vary. Making intake, area or stored energy directly proportional to both strengths fails for unequal masses. A single medium could still provide a continuous inward current and a separate signed state. [Two-radius calculation](#first-comparison-across-two-mouth-radii)

- **Light.** An ordered three-dimensional medium with the assumed rotational shear stiffness and inertia carries two transverse polarizations. A compressible model also has a longitudinal wave. Which modes ordinary matter emits and detects, and how much energy escapes into the shear-free bulk, depend on the material couplings. In particular, a wave that reverses under exchange of the slab faces has zero linear coupling to ordinary face-symmetric shear in the symmetric resting slab. An instrument made sensitive to opposite face rotations can detect that wave. Identifying it with ordinary light requires an independent material readout. [Light requirements](steps/S9_light_requires_shear.md) · [Probes and power](../../software/stage1_solver/reports/optical_matter_and_bulk.md)
- **Gravity.** Continuity, carried momentum and pressure give an attractive inverse-square load in the ideal two-drain exterior. At a free mouth, some of the incoming momentum can pass into the receiving bulk or be stored in bound material. The same held load can give different free accelerations. We need the intake selected at each physical radius, the pull delivered to the translating mouth, and its independently determined inertia. In the flow-driven limit, source flow per source mass and pull per incident speed divided by target inertia must each be common across mouths. There is also a cross-check with light. With constant light speed and zero tilt, the flow-only optical solution requires a radial speed proportional to \(r^{-1/2}\), while a steady spherical constant-density drain with conserved exterior flux gives \(r^{-2}\). A common material theory must explain the optical profiles as well as free fall. [0PN result](../../docs/zero_pn_gravity.md) · [Optical comparison](#can-gravity-and-light-use-the-same-radial-flow)
- **Electric charge.** The relaxed scalar bend gives attraction between like orientations in the tested slab. A geometric wrapping count can unwind, and the existing conversion equations do not determine where an electric unit goes. Circulating a trapped wave can provide angular momentum, but pair formation permits either handed assignment. Two lower-level candidates made more progress. A radius-scaled responding layer with unequal gains for the two circular polarizations selected one supplied flux spacing. Models with angular contacts on the two material faces sustained a signed turn and carried transverse power. Prescribed motion of oppositely prepared mouths drove a distant detector sensitive to those turns. With strong rim attachment, weak coupling to incoming material let the turn reach the exterior. Strong coupling to the feed screened it. The open questions are what selects the same sign convention across independent mouths, what sets the electric unit, and what keeps or transfers it as material drains. [Charge tests](../../docs/em_charge_attribute_requirements.md#november-checkpoint--what-a-charge-law-must-explain) · [Flux selection](../../software/stage1_solver/reports/mouth_free_flux.md) · [Moving-source and detector test](../../software/stage1_solver/reports/moving_mirror_mouths.md) · [Transmission and screening](../../software/stage1_solver/reports/mouth_collar_charge_transmission.md)
- **Why we tested the trapped wave.** A localized rim excitation can let one passing parcel affect the next in a trial material. The tested conversion and standing-wave budgets fell short of paying for a flux unit at the sampled settings. We asked whether a probe could distinguish the energy to add the first excitation from the energy to add the second. Resolving those lines would test a property of storage. The matched scalar preparations give the same addition energies on either mirror. Wall motion and radiation broaden the lines. At the same intake, two wall profiles with the same speed contrast and end-to-end travel time gave opposite peak counts when the slow region shifted along the bend. **The prescribed-wall calculation is parked.** Reopen it with a known force and independent wall displacement versus depth and time on that drain, or a material equation that predicts that response. The fixed-end inertia leaves the 0PN inertia of a translating mouth open. [Conversion work](../../software/stage1_solver/reports/mouth_conversion_twist_work.md) · [Rim storage](../../software/stage1_solver/reports/mouth_rim_retained_excitation.md) · [Wall ripple test](../../software/stage1_solver/reports/mouth_wall_ripple.md)
- **Magnetism and gravitomagnetism.** A moving-mouth shear calculation recovered the magnetic force's directional and velocity dependence with an independently adjustable coupling. The proposed curl of displacement also has different full time-reversal parity from the Maxwell magnetic field. The static charge must fix the moving current and force with the same normalization. For frame dragging, the target is a gyroscope response when the centers of many inward-draining mouths carry a known total angular momentum. Its sign and strength must follow from that angular momentum and the same gravitational coupling that accounts for free fall. [Magnetic force test](../../software/em_charge_attribute/magnetism_moving_throat_result.md) · [Gravitomagnetic requirement](V3_STEP_PLAN.md#s19--gravitomagnetism--not-a-cite-and-record-step--it-has-a-known-open-problem)

**Next material decision.** Specify what can remain attached to an open mouth as constituent material passes through, what stores its orientation, energy and sideways momentum, and what is exported to the bulk. For the angular candidate, independently determine how strongly material at the rim and in the incoming flow resists a small imposed relative turn. Use those responses to predict whether a mirror pair gives opposite distant readings or both are screened. For gravity, the independent calibration is the held pull in a known flow and the acceleration under a known push of the same mouth in resting brane. Then reverse mouth orientation at matched intake and change physical radius with the coefficients held fixed. The electric source should reverse and retain its unit magnitude. At fixed separation, gravitational acceleration should persist on reflection, be proportional to source mass, and be independent of target mass. The angular field and any retained layer can add calculable inertia in a chosen model. Count these contributions once in that push calibration. [Angular-field inertia](../../software/stage1_solver/reports/mouth_orientation_and_motion.md#what-does-this-supply-for-0pn) · [Standing method](../../memory/topics/foundations.md#how-we-use-the-force-requirements) · [Comparison plan](#plan-for-testing-one-medium)

### Requirements by sector

| What the common model must supply | How to tell whether it does |
| --- | --- |
| **Gravity source strength.** At fixed ambient conditions an open mouth of radius \(R\) selects an inward mass rate \(D(R)\). Its far-field strength must scale with its independently measured inertial mass. | Measure intake at the mouth and flux through surrounding spheres at several radii. Define the exterior volume rate \(Q_s\) from those spheres. Account for exchange between them before identifying \(Q_s=D_s/\rho_\infty\). Compare \(Q_s/M_s\) across physical source mouths. [Source and mass requirements](#how-drain-rate-and-retained-momentum-must-scale-for-free-fall) |
| **Gravity response.** Momentum delivered by the flow must accelerate an unheld mouth toward the companion at the same rate for different target masses in the same source field. | Obtain the pull per unit weak incident material speed \(C_t\) from the opposite holding reaction, then determine the slow inertia \(M_t\) of the same mouth under a calibrated push in resting brane. In the scalar, flow-driven limit, \(a_t=(C_t/M_t)Q_s/(4\pi d^2)\). A common positive \(C_t/M_t\) and common \(Q_s/M_s\) give one \(G\) when source and target are exchanged. Keep any leading bulk drive or coupled acceleration in the response if present. [Held and free response](../../software/stage1_solver/reports/free_fall_response_ratio.md) |
| **Sustained motion and energy.** A steady drain must supply its interaction work and preserve Newtonian motion to leading order. | Balance energy, momentum and torque in the mouth, brane and bulk. The prescribed-approach benchmark with one stationary mouth requires excess inner power \(2Fv\), for attractive load magnitude \(F\) and approach speed \(v\). For freely moving mouths, test for leading drag, mass drift or orbital energy loss using the same parameters. [Energy benchmark](../../docs/_measurements/stationary_mouth_energy_review/round2/status.md) · [Orbital criteria](#what-must-persist-once-both-throats-are-orbiting) |
| **Electric source and magnitude.** Each elementary charged mouth must select one signed, additive unit coupled to the long-range field. Mirroring the mouth reverses the sign. Its magnitude persists as radius and trapped energy vary on the same charged branch. | Hold the charge assigned to the flowing material fixed. Subtract ordinary charge in transit and read the excess flux with one distant probe and one field normalization. Test the mirror and multiple radii. The existing reflected carrier gives equal probe readings. Assigning opposite charges as inputs bypasses this test. [Mirror probe](../../software/em_charge_attribute/reports/mirrored_mouth_probe.md) · [Charge targets](../../docs/toy_model_ontology_summary.md#113-additive-charge-and-conservation) |
| **Electric stability and conservation.** The charge must stay localized on the lifetime of the state. Any discharge must carry the lost unit into an identifiable recipient. | Compare a relaxing mouth with a held mouth at matched intake. Sum the change of excess charge bound to the mouth and layer and the current through the surrounding boundary. Establish what protects the charged state or what it costs to discharge. The same unwinding geometry admits conserved trial currents giving one unit or zero at the bulk edge. [Conversion test](../../software/em_charge_attribute/reports/conversion_bridge.md) · [Charge spreading](../../software/em_charge_attribute/reports/bulk_edge_charge.md) |
| **Electric force.** The electric contribution must repel like signs and attract opposite signs with a common inverse-square strength for elementary units. | Use the signed source determined above to predict the total electric reaction on a second physical mouth, including stresses and support work. Test both signs and unequal radii with the same parameters. A passive midpoint bend gave the opposite sign pattern in the relaxed benchmark. [Two-mouth results](../../docs/em_charge_attribute_requirements.md#november-checkpoint--what-a-charge-law-must-explain) |
| **Light propagation and polarization.** The ordered brane must support two stable transverse light polarizations with the same speed in its homogeneous, isotropic rest state. | Derive the material inertia and positive transverse stiffness about an equilibrium slab. The supplied curl-only action gives \(c_\gamma^2=\mu_R/\rho_{\rm br}\) and two transverse directions for three brane dimensions. Check the full spectrum including out-of-plane motion, and identify the angular momentum carrier and physical reference frame for rotational stiffness. [Shear requirement](steps/S9_light_requires_shear.md) · [Polarization count and its assumptions](steps/S10_two_transverse_photons.md#open-requirements-on-interpreting-the-count-as-light) |
| **Light confinement.** The same material must carry optical shear in its ordered phase and lose that shear response in the disordered bulk. A curved or thinning brane can still convert light into other motion. | Derive the shear response through the order transition. At a nonuniform mouth, account for incident energy in reflected and transmitted light, longitudinal motion and bulk motion, including any work by the drain or its support. The uniform linear model decouples transverse light from the bulk. The nonuniform calculation yielded no accepted loss value or upper bound. [Two-phase requirement](steps/S9_light_requires_shear.md) · [S11c partial closeout](steps/S11c_PARTIAL_CLOSEOUT.md) |
| **The extra longitudinal wave.** When the slab resists compression, the homogeneous model also supports a longitudinal material wave. Its visibility depends on how matter emits and absorbs it. | For a throat that moves, changes shape or changes its trapped-wave state, determine the energy released per frequency and time into free transverse light, longitudinal motion and the bulk. Determine uptake from each channel as well. The resulting rates, escape and population enter the comparison with thermal radiation. [Longitudinal spectrum](steps/S11_stray_longitudinal.md) · [Emission and absorption requirement](#r-s22-01--how-the-energy-a-throat-releases-splits-among-light-the-extra-wave-and-the-bulk) |
| **Light in a gravitational field.** One set of material profiles around a source must reproduce both weak light deflection and the logarithmic slope of radar delay, using an independently determined orbital \(GM\). | Use the same far-zone material velocity, local shear-wave speed and brane tilt for both observables across impact parameters. At constant wave speed and zero tilt, the retained nonzero optical solution requires a speed proportional to \(r^{-1/2}\). The steady, constant-density drain with no exterior leakage gives \(r^{-2}\). If they refer to the same material velocity, the profiles disagree over a radial interval. The general optical calculation leaves speed and tilt responses open. [Optical requirements and domains](#r-s21-01--conditional-deflection-and-radar-slope-profile-requirements) · [Common-flow comparison](#can-gravity-and-light-use-the-same-radial-flow) |

**Moving charge and moving mass**

| What the same mouth must supply | How to tell whether it does |
| --- | --- |
| **Magnetism.** The measured electric source of a stationary mouth must determine the current when that mouth translates and the force on a second moving charge. | Use the same static electric strength and optical normalization for parallel and antiparallel translations, with inward intake throughout. The direct shear calculation recovered the reference angular dependence, inverse-square force and order in the two velocities. Its coefficient ratio \(q_T^2/(\rho_{\rm br}A_E)\) depends on the moving-mouth coupling \(q_T\) and the static pair coefficient \(A_E\). Determine both independently and include active momentum flux in the total force. Reversing translation on the draining branch tests the motion-dependent signal. Under full time reversal the proposed magnetic field, the curl of material displacement, is even. Maxwell's magnetic field is odd. Carry that measurable distinction into the comparison. [Moving-throat calculation and limits](../../software/em_charge_attribute/magnetism_moving_throat_result.md) |
| **Gravitomagnetism.** The same inertial masses that source the Newtonian comparison must source frame dragging when many mouths carry angular momentum in a rotating body. | Give the mouth centers orbital motion. Keep the intake of each mouth inward. Reverse the body's angular momentum at matched mass distribution and compare the far-zone precession of a test gyroscope. Look for a dipole with sign following angular momentum and strength fixed by the independently measured angular momentum and gravitational normalization. Trace the angular momentum in the material, trapped waves and exported flux. Test whether its distribution is compatible with the compact inertial source. An extended vortex tail is one proposed resolution. Its necessity and its energy cost are open. [S19 question](V3_STEP_PLAN.md#s19--gravitomagnetism--not-a-cite-and-record-step--it-has-a-known-open-problem) · [Compactness and spin tension, C10](DEFECT_REGISTER.md) |

**Compatibility test.** Use one set of brane, throat and conversion laws. First change only one mouth's orientation, holding masses, radii, inward intakes and ambient conditions matched. The gravitational pull should persist and the electric force should reverse. Measure all four orientation pairings and trace the corresponding momentum budgets, since an orientation-independent force can have other sources. Then change the physical radius on the same elementary-charge branch. Gravity must track the independently measured source inertia while the magnitude of the electric unit stays fixed. Compare absolute strengths using the same parameters in both sectors. [Four-orientation comparison](../../docs/em_phaseC_force_decomposition.md#3-the-four-channel-decomposition-the-honest-bookkeeping) · [Gravity attribution](#what-must-persist-once-both-throats-are-orbiting)

**Charge and light on that same throat.** Once a stationary mouth's charge has been determined, use the same coupling to predict how an accelerating mouth excites the far-field transverse waves. If an electric theory introduces a gauge field, identify how its propagating polarizations relate to the two brane shear polarizations. A proposed trapped shear wave must also remain bound, or live long enough to supply the assumed throat support, with losses entering the same energy budget. [Free and trapped light](../../docs/toy_model_ontology_summary.md#9-light) · [Shear and electric modes](../../docs/em_charge_attribute_requirements.md#2-requirements-inversion-what-the-charge-attribute-must-provide)

Count stresses, stored momentum and propagating modes once. For example, order-front tension that moves with the slab already contributes to its bending stiffness. Reopen bend-mediated gravity if an isolated free mouth shows a surviving \(1/r\) midpoint displacement and an independently accounted opposite reaction. Reopen charge selection when a specified coupling predicts a mirror-odd, conserved electric source whose localized state can be tested. An effective law at the mouth could meet these tests before the full interior is resolved. [Order-front stress](../../software/stage1_solver/reports/order_face_normal_traction.md)

**Normal response of a proposed material side wall.** To use the geometric solver's reflecting side wall for an attachment-force prediction, establish its normal current response with the physical material on the far side. Measure the material sites and the effective carrier boundary in the same coordinate system. A cold sheet that reflects individual constituents can still transmit a steady carrier current when both sides share one phase. If there is appreciable crossing on the gravity timescale, include its momentum and the reaction on the material in the mouth budget. Axial intake and the trapping of transverse light have their own tests. [Neutral surface and normal-current comparison](../../software/stage1_solver/reports/mouth_surface_registration.md)

At the same dense ordered-to-converted interface, measure the mass and tangential momentum carried by each direction of crossing, the traction on the ordered material, and local momentum storage. A finite normal pressure or zero net mass crossing alone leaves the sideways reaction open. The cold continuous carrier, dilute gas and warm dense sample currently constrain different parts of this requirement. [Neutral interface tests and stopping decision](../../software/stage1_solver/reports/mouth_converted_interface_gate.md)

### First comparison across two mouth radii

Take two physical mouth preparations in the same ambient brane, with independently determined positive inertias \(M_1,M_2\), exterior volume flow rates \(Q_1,Q_2\), and the same elementary electric magnitude \(e>0\). On the flow-driven gravity branch the source target is \(Q_i=k M_i\) with one positive constant \(k\). Test the proposal that a single mouth amplitude \(x_i\) sets both strengths through fixed proportionalities \(Q_i=a x_i\) and \(|q_i|=b x_i\). Calibrating the electric unit on mouth 1 then predicts \(|q_2|=eM_2/M_1\). For unequal inertias this conflicts with one common electric unit. Equal inertias leave the test inconclusive. SymPy compared the two source rows. Mathematica eliminated the amplitudes. These are symbolic preparations. No physical mouths of different mass were measured. [Commands and literal results](_measurements/charge_gravity_source_scaling.json)

This rules out using intake, mouth area or stored energy as a common *proportional* source for both forces on that branch. One material field can still have a size-dependent flow and a normalized signed state. For example, a conserved signed count would have the right size dependence if the mouth bound it to an electric current. The tested wrapping can unwind, and the existing action leaves its electric coupling and destination undecided. [Source scaling](#how-drain-rate-and-retained-momentum-must-scale-for-free-fall) · [Conversion-layer currents](../../software/em_charge_attribute/reports/conversion_bridge.md)

### Plan for testing one medium

1. **Test the resting medium first.** Fix one candidate's fields, equations, constraints, face conditions, equilibrium slab and any reservoir. Check that the slab is a stationary solution and mark each independently calibrated parameter. Linearize that candidate about the flat slab. Compute the static response to a weak charge probe. In the same linear equations, compare the waves excited by a time-dependent electric current with those excited by a transverse mechanical drive. Determine whether both excite the same light modes. Identify the long-range static response, which propagating modes each probe excites, their energies and speeds, and which can escape into the bulk. An existing scalar realization has an electrically visible extra wave under its stated inertia and source assumptions. A later candidate must account for its visibility. Compare the gravity-signal and shear-light speeds separately. The existing benchmark calibrates their equality. [Earlier coupled-spectrum test](../../software/stage1_solver/reports/pathA_39_scalar_admixture_screen.md) · [Speed comparison](../../software/stage1_solver/reports/pathA_40_cone_lock.md) · [Candidate equations and spectrum plan](../../docs/cross_sector_research_burdens_and_compatibility_gates.md#cheapest-decisive-calculation)
2. **Identify what the mouth adds.** For one stationary inward-draining mouth and its mirror, derive the electric source and the mechanical traction from the same proposed interface laws. Two candidate additions have different jobs. One trial coupling pairs the electric potential with a conserved defect current. It must explain which unit the mouth binds, the energy keeping it localized and the current at the bulk edge during unwinding. A law relating material motion and stress at the mouth must give the transmitted force, stored momentum and power that conversion supplies. A trapped wave can contribute to this inertia. Equal steady pull already permits different slow accelerations. Record any required coefficient, constitutive function or boundary datum and an independent way to determine it. [Charge-carrier test](../../software/em_charge_attribute/reports/conversion_bridge.md) · [Force and inertia counterexample](../../software/stage1_solver/reports/free_fall_response_ratio.md)
3. **Reserve a second preparation for prediction.** With parameters fixed, determine the first mouth's total holding reaction in weak flow and its acceleration under a known push in resting brane. Combine the readings to predict free fall. Repeat at another physical radius and with the mirror. Test \(Q_s/M_s\), \(C_t/M_t\), the signed electric unit and its lifetime. Then use the measured charge and mass currents in the motion tests above. Check the brane flow, local light speed and tilt against both optical observables using the same independently determined orbital \(GM\). The radial powers already disagree for the simple drain and the flow-only optics model. A proposed remedy through variable light speed needs a material law for how shear stiffness and inertia respond to the local state. Tilt needs its own mechanical balance. [Gravity response and optical profiles](../../docs/zero_pn_gravity.md) · [Common-flow comparison](#can-gravity-and-light-use-the-same-radial-flow)
4. **Judge each proposed addition by what else it predicts.** Record which observations fixed its parameters and which stayed available as tests. A common set of stable equations passing the tests establishes compatibility on the tested regimes. A freely adjustable coefficient leaves that observable undetermined. A contradiction at fixed parameters rules out the particular candidate on that domain. Before proposing another force, identify the missing source, stress, symmetry or boundary relation, and show which independent observable a change would affect. A claim about every one-medium theory requires a contradiction across the allowed class. [Charter and falsification standard](CHARTER.md) · [Integration remit](V3_STEP_PLAN.md#s21--what-light-and-gravity-share)

### Can electric stress travel in the material wave that carries light?

In a working homogeneous-brane test, an independent Maxwell field coupled to material polarization can preserve the Coulomb tail and change the light speed while elastic shear travels on a different branch. A second candidate identifies electric displacement with the curl of an axial material coordinate. Its two transverse waves and its static electric stress share an energy density. A smooth, single-valued coordinate gives zero electric flux through a sphere surrounding a mouth. A circle-valued transition between two local descriptions allows integer flux in units set by one common calibration. Giving that phase uniform positive stiffness throughout the brane opens a gap in both transverse waves. These calculations supply tests for a proposed material law. [Resting-brane comparison](../../software/stage1_solver/reports/resting_brane_electric_shear.md)

For the next mouth law, let the mouth relax without assigning its flux. Compare its mirror and a neutral control, then change radius with the same coefficients. A candidate must select the unit and its sign, account for the current when the state changes, and leave the distant brane able to carry gapless transverse waves. A measurement of net electric stress on an enclosing surface still needs a physical rim response to predict acceleration. The 0PN mass and intake calibrations and the later moving-source tests for magnetism and gravitomagnetism remain independent requirements. [One-mouth electric currents](../../software/em_charge_attribute/reports/conversion_bridge.md) · [0PN gravity requirements](../../docs/zero_pn_gravity.md)

One specific candidate uses the geometric phase of a circularly polarized trapped wave on the saved closed mouth sections. If the polarization follows the projected face normal, either helicity supplies a signed geometric count of magnitude two while that normal covers its sphere once. The full polarization pair has zero net count. A proposed material coupling confined to an ordered mouth section gives the exterior flux in units of its elementary spacing as the polarization count plus the total signed winding of vortices in the relative phase. This gives a concrete energy test: at fixed coefficients and intake, can the mouth select a helicity and a stable flux state when vortices and loss of coherence are allowed? The earlier unwinding path passes through a point where the polarization direction becomes undefined. The new law must specify how charge traverses that event. [Conditional flux law and controls](../../software/stage1_solver/reports/resting_brane_electric_shear.md#4-can-the-trapped-transverse-wave-supply-that-phase)

### Energy selection of one electric flux spacing

A fixed-radius, round-mouth calculation now lets both helicity amplitudes, their vortices, and their coherence vary across integer electric flux sectors. A thin surface with fixed material coefficients chose zero, one or two spacings at different radii. A second candidate has a three-dimensional responding region whose width is proportional to mouth radius and whose ordering energy depends on the curvature of nested material sections. With one common set of coefficients, it favors one spacing for all radii of that mathematical family. Equal helicity gains leave both charge signs equally favorable on each mouth. Opposite signs on mirrored mouths require a reproducible handed preference with the same sign in the local-normal convention. No such preference has been measured. Holding the region's physical width fixed allows a sufficiently large mouth to prefer zero flux. [Energy comparison and controls](../../software/stage1_solver/reports/mouth_coherence_energy.md)

The original comparison fixed the angular flux distribution. Charge transport during a slip, the stress delivered to a free mouth, and the inertia associated with the new binding energy are open. [Original guarded results](./_measurements/mouth_coherence_energy.json)

With the electric field free to redistribute, the neutral sector can lower its energy by developing patches of flux of opposite signs. An all-field energy comparison nevertheless selects exactly one of the assumed flux spacings at every fixed radius in an example with a thicker responding region and an 18-to-1 ratio of circular gains. The example establishes sufficiency. What fixes the spacing and the actual material gains remains open. [Energy bounds, controls and domain](../../software/stage1_solver/reports/mouth_free_flux.md#a-comparison-over-all-admissible-fields)

The relation between mouth orientation and electric sign must hold for mouths formed independently. Pairwise neutral creation and local angular momentum balance allow either relative charge on two upward mouths from different pairs. A simple passive alignment field strong enough to favor like charges on same-side mouths also attracts them on that branch. A candidate shared handed state must explain how distant preparations inherit the same sign, predict an independently measurable circular response near the mouth, and preserve the propagation of both transverse light polarizations in the resting brane. Its stress also enters the 0PN force and inertia budget. [Formation and alignment tests](../../software/stage1_solver/reports/mouth_free_flux.md#what-selects-the-sign) · [Guarded results](./_measurements/mouth_free_flux.json)

### Can the two material faces give an electric sign without selecting a trapped-wave polarization?

Partner swaps can restore every directed bond at a draining cut. In a separate model a two-component internal state compares the material normals on neighboring patches. Its signed loop sum for the two physical faces is \((s,-s)\) on a wrapping mouth of orientation \(s\). An equally coupled loop probe measures zero. A locally bound pair can carry the sum of the constituent phases. If only this pair has bulk phase stiffness in the ordered brane, one opposite-face connection remains gapless and the computed mouth count is its smallest nonzero integer. The two exterior stresses add in the force between matched sources. The sign, transparent propagation and physical force depend on the explicit constitutive assumptions in the working calculation. [Material connection tests](../../software/stage1_solver/reports/mouth_material_connections.md) · [Guarded results](./_measurements/mouth_material_connections.json)

Determine the microscopic angular response and attachment energy on a neutral patch. At finite pair density, test stiffness under common and opposite phase twists, register the face orientation across independent sites, and compare mechanical shear and flux-sensitive probes of the same propagating mode. Determine the cost and recipient of a compensating unit during conversion. A signed material turn still requires an electric normalization and a coupling to the force on a moving mouth. [Single-mouth and light requirements](#can-electric-stress-travel-in-the-material-wave-that-carries-light)

### Can the paired material carry both the signed field and ordinary light?

The required common and opposite phase stiffnesses occur in explicit three-dimensional models at finite density. Local binding alone leaves the bulk outcome open. In a reaction-equilibrated ideal mixture, heating across the calculated pair-condensation temperature erased both responses while most constituents remained bound. An interacting constituent mean-field calculation gives the same distinction between the ordered and normal state under its own parameter choices. [Finite-density calculations and controls](../../software/stage1_solver/reports/mouth_pair_bulk_stiffness.md)

Physical face reflection gives a sharper light test. If a single compact field changes the sign of the mouth flux under reflection, its linear transfer from ordinary reflection-even transverse displacement vanishes in a symmetric resting state. A face-sensitive mechanical signal can still follow the signed curl. For the tested odd thickness oscillator its ordinary share of in-plane wave power vanishes quadratically with wave number. A lower law identifying an odd internal motion as light must predict its restoring stress and energy transport. A proposal breaking face reflection must independently recover both the mouth sign and the far-brane optical spectrum. [Mirror test and energy budget](../../software/stage1_solver/reports/mouth_pair_bulk_stiffness.md#can-ordinary-transverse-shear-be-that-same-electric-wave)

### Can internal rotations carry the signed field and be read as light?

An explicit array of two counterrotating material faces transmits long-wave power through local spring torques. In its linear regime two transverse waves carry this energy, and there are also observable angular motions with zero restoring force. With separate full-turn contacts and finite attachment to a prescribed radial mouth geometry, the zero-temperature energy favors a single signed material turn throughout an ordered exterior. Sensing only the relative angle allows a primitive unit to escape that readout. When its site rolls equilibrate, a weak tilt of a prepared mouth patch transmits a circulating torque of opposite sign on the mirror and couples to the two-wave subspace at finite wavelength. The calculated overlap decreases toward zero with wave number for this localized tilt. [Hypotheses, bounds and probe tests](../../software/stage1_solver/reports/mouth_torsional_light.md) · [Measurements](./_measurements/mouth_torsional_light.json)

The material must supply the angular registrations, the mechanism transmitting loop torques, their recoil and the ordered-to-disordered change in stiffness. Determine whether known optical matter responds to the face-odd rotation and which ordinary displacement transports its power. The existing face-even drive gives zero linear response in this coordinate. An independent physical electric unit, the complete free-mouth force and the 0PN inertia remain open. [Optical displacement requirement](./steps/S9b_light_bending.md) · [Mirror test](../../software/stage1_solver/reports/mouth_pair_bulk_stiffness.md#can-ordinary-transverse-shear-be-that-same-electric-wave)

### Which material can detect a signed wave and confine its power?

A neutral, strongly bound two-face pair on a closed loop provides a signed material readout in the proposed angular medium. A weak loop twist induces opposite circulating constituent currents. A common real strain produces zero signed current. Resonant absorption from the symmetric ground state into the low paired band is forbidden by face exchange for the tested preparation. Higher separated states permit absorption and also supply the virtual response below threshold. A two-site open link requires an independent time reference to separate a physical drive from a change of local phase convention. A two-face instrument with explicit absorbers verifies finite-distance power transfer between matching drives and detectors. The earlier ordinary optical displacement remains decoupled from the signed branch in a symmetric, resting slab. [Material readouts and controls](../../software/stage1_solver/reports/optical_matter_and_bulk.md)

Specify the dynamic bulk torque as well as its relaxed static stiffness. At a flat face, whether power crosses depends on frequency, conserved in-plane wave number and the available bulk receivers. One sliding-contact model at zero static bulk stiffness had a guided mode with no outward power while its tails held 64.02% of total mechanical energy. Drag on those same degrees of freedom attenuated the wave. Disengaging the contacts gives a different boundary law. The ordered material must meet the intended localization and attenuation requirements over its optical band, including any normal acoustic route. A finite-wavelength signal from a driven mouth still requires the relation between moving signed flux, compensating current and wave power. [Guided mode and energy budget](../../software/stage1_solver/reports/optical_matter_and_bulk.md#loss-of-order-and-optical-power)

### A moving source for the signed angular wave

A local constitutive law must specify the time registration of material turns and the angular traction delivered by a moving mouth. For a prescribed compact texture, the signed current crossing lattice planes accounts for the change in dipole. Finite exterior cuts can retain the same count for different interior crossing histories. Two finite spherical contacts in fixed opposite sectors, with prescribed small relative translation and equilibrated angular slips, launch outgoing waves through the exterior curl stress. A symmetric constituent loop detects the resulting signed flux with a current amplitude approaching inverse distance. Common translation greatly suppresses the total long-wave power and has an exact directional zero for the equatorial detector. [Source, wave power and material readout](../../software/stage1_solver/reports/moving_mirror_mouths.md)

Fix the torque-to-slip amplitude and phase independently over the relevant frequencies and angular patterns. A positive internal reference inertia and drag can change both the radiation and the drive work even when the static two-mouth fields agree. A held-mouth scattering calculation measures the torque-to-slip response. A passive translation-sensitive inertia leaves held angular scattering identical but changes outgoing pair power. For the tested isotropic extension, a second calibration using prescribed motion of one mouth fixes the missing response and predicts the pair without fitting its radiation. Physical charge calibration, the optical displacement, confinement, channel-dependent loss and the 0PN translational momentum balance remain open. [Measured outcomes and review](./_measurements/moving_mirror_mouths.json)

### Does the orientation-bearing material persist at an open mouth?

Track orientation-bearing sites separately from the material throughput. For an isolated round mouth moving along an axis, a conserved layer with relative tangential displacement \(uP_n\hat z\) has a fractional site-density dipole \(2u\cos\theta/R\). Measure this dipole, its local sources and the residence times of labeled sites during weak motion. An equilibrium site mass per area and area modulus, together with the held-patch orientation spectrum and bond stiffness, predict the isolated response under local retention and common bond endpoints. A net replacement dipole can evade this constraint even at zero total replacement imbalance. Two otherwise identical imposed registrations gave opposite wave phase at one tested frequency. A matched held-mouth force difference can test that registration independently. [Material and force tests](../../software/stage1_solver/reports/mouth_orientation_and_motion.md)

In the regular slow limit the isolated radial transverse angular field contributes \(Cq^2/(6\pi c^2R)\) to the inertial response. A retained material layer contributes \(4\pi R^2\Sigma\). Count these terms once in the push calibration on the same draining state. Establish whether the angular material survives there and whether it actually carries optical power. The companion-directed load and the rest of the inertia must still be determined before forming a free-fall acceleration. [Two derivations and scope](../../software/stage1_solver/reports/mouth_orientation_and_motion.md#what-does-this-supply-for-0pn) · [Guarded results](./_measurements/mouth_orientation_and_motion.json)

### What orientation reaches the outside as material is replaced?

Measure three quantities independently on one drain. Carrier IDs and front/back entry and exit readings test whether the orientation sites travel with the inflow. The two phase components of a weak tilt on neighboring sites test whether angular motion transfers between carriers. The decay of an excited state relative to the local preferred axis, together with the distribution of contact histories, tests whether initially uniform inflow can become wrapped. In the uniform six-contact vacuum receiver example the sign appears on the inner cut once excited-state retention falls below one half at the opposed pole. The hypothesis consumes more than a quarter of a transition quantum per participating site on average to reach that branch. Confirm a material sink for that energy at the local transition frequency. [One-mouth tests](../../software/stage1_solver/reports/mouth_turnover_and_orientation.md)

An inner preferred-orientation winding still needs a physical exterior flux. The supplied uniform outer orientation stays neutral and the feeding layer carries a compensating index. Determine how the physical angular link stiffness and preferred phase behave when the reduced material polarization vanishes, then measure the equilibrated link flux on both sides. A probe reading the full spin–receiver state or a single receiver can give zero even when the preferred spin winds. Count any net material reaction and energy transfer once in a future draining-mouth momentum balance. [Measurements and scope](./_measurements/mouth_turnover_and_orientation.json)

### Which part of the mouth couples its orientation to the angular field?

A free angular field can carry the inner signed unit past the feeding layer. In a model that couples to the dominant local orientation, weak attachment to incoming material and strong rim attachment transmit the unit. Strong attachment to the feed screens it. Ordinary orientation-independent transfer with the same preparation gives neutral ground states throughout the tested receiver family. An excitation that aligns with the local material prefers opposite phases on mirror loops. One excitation shared by all rim sites has insufficient energy to pay for a unit at the tested alignment ratios and hopping scales. Determine the held-patch phase response, the splitting and occupation of rim modes, and how their available energy varies with mouth radius at matched material parameters. The source-free exterior and identification of a physical electric probe require independent tests. [Free-link comparison and bounds](../../software/stage1_solver/reports/mouth_collar_charge_transmission.md)

A single delay before renewing contacts can exert torque can give strong attachment at the rim and weak attachment in the feed when their residence times differ. The signed winding survives replacement in the tested mean-energy ground state if a persistent local reference is supplied. Measure the torque per surviving contact against its age, its residence distribution and the energy spent on registration. A separate measurement of exit momentum and total reaction under steady lateral motion tests whether those contacts give their acquired sideways momentum to the bulk or transfer it back before leaving. Direct export produces a linear drag in the tested law and constrains its use for 0PN free fall. [Renewing contact and push test](../../software/stage1_solver/reports/mouth_registering_rim.md)

An unloading contact can give its sideways momentum back before exit. At the same angular preparation and intake, the measured exit momentum can agree with several different rim forces if another part of the medium takes some of the reaction. Fix the motion and ownership of that support. With one incident stream supplying only sideways kinetic power, a stationary receiving support and no draw on stored or inward-flow energy, the dimensionless held-pull coefficient \(c\) and towing coefficient \(g\) must satisfy \(c^2\le 2g\). A finite \(c\) at vanishing \(g\) requires an additional power channel in this steady uniform-flow comparison. The part of the rim actuator power odd in both incident speed and rim speed isolates the work to trace through pressure, conversion and outflow at the same intake. [Release reactions and energy balance](../../software/stage1_solver/reports/mouth_release_reaction.md)

Pressure work and kinetic transport cancel in that mixed part on the ideal exterior cuts. The translated neutral-carrier feed separately supplies only the uniform stream kinetic increment. Test whether the steady inward energy drop pays for the center work by changing the ordinary loss in conversion. If the loss per transferred mass is the same function \(L(D)\) as for a resting mouth, fixing intake with a resting bulk requires a mixed static head change \(+UV\). Fixing the two reservoir heads instead predicts a leading mixed intake change \(-UV/L'(D_0)\). Determine the rest slope independently, resolve heat, storage and other powered supports, and retain a separate force readout for the moving center. [Mixed power and pressure–intake test](../../software/stage1_solver/reports/mouth_stream_power.md)

A passive order front attached to the mouth can select the normal throughput from the pressure difference applied directly across it. Test its independently measured mobility and tension against the intake slope. For a thin front spanning the opening, measure the normal depth slope around the rim and the free decay rate relative to the material. A center-depth reading alone depends strongly on the spatial discharge profile. Establish how much of the far-field head reaches the conversion surface. A sideways contribution from normal traction needs its load variation across the opening and its mechanical reaction on the center. [Pressure-selected branch and edge test](../../software/stage1_solver/reports/mouth_front_intake.md)

For stationary passive conversion, identify the net energy supplied per converted constituent at the physical brane and bulk reservoirs. The saved intake of 0.25 would consume 0.05123 per constituent on its finest stored preparation. A stationary wave of fixed action can exert a load while the continuing power must come from another energy flux. Also determine whether the mouth and intact faces belong to the same connected regular order level. With one coherent velocity and reciprocal order force, their relative normal currents on such a level have weights proportional to mobility times density squared divided by the order gradient magnitude. A full theory must explain how the intake stays concentrated at the mouth under the same distant driving conditions. [Coupled energy and current requirements](../../software/stage1_solver/reports/mouth_coupled_conversion.md)

Distinguish conversion relative to a material face from coherent current through a fixed plane. To localize drainage at a stationary mouth, match the normal current, phase access and traction under the same reservoir conditions. A local curvature threshold alone also opens unpunctured regions of the same curvature. At matched curvature the neutral carrier and surviving order surface respond identically beside a prepared puncture and on the intact bend. A viable normal transport barrier must allow ordered feed to cross the downstream conversion layer, even where the same intermediate order occurs on intact material. Independently measure what selects normal phase stiffness and how momentum and pressure cross the intact moving interface. [Bent puncture and downstream transfer](../../software/stage1_solver/reports/mouth_bent_puncture_transport.md) [Matched pressure, phase and mobility tests](../../software/stage1_solver/reports/mouth_face_conversion_lock.md)

**Working boundary for the 0PN force tests.** Intact slab faces pass no steady normal carrier current. A localized downstream melting exit remains open. In a coherent crystalline field at matched curvature and shear order, imposing a dislocation opens a channel at either patch. The free crystal selects the same sealed state at both. The material origin of the distinction and the forces supporting a stationary flare remain open. [Crystalline connectivity test](../../software/stage1_solver/reports/mouth_crystal_connectivity.md)

**Energy condition for a steadily draining physical mouth.** At any radius, positive flow from ordered material through a stationary, finite-mobility passive conversion layer costs positive power. For a smooth unit-change order marker, the throughput squared bounds the product of the conversion power and the integrated mobility weighted by density and the marker derivative. A constant spatial phase offset supports coherent current with zero net port power. Identify the independent specific energy supplied by the prepared phases or an alternative reversible conversion law before selecting an equilibrium depth and flare. Release the radial constraint when testing whether one or several radii are stable at a given preparation. [Two-radius energy and free-shape test](../../software/stage1_solver/reports/mouth_free_shape_selection.md)

**Primary conversion hypothesis for 0PN.** The ordered brane and disordered receiving bulk share a coherent carrier. Use the reversible Hamiltonian order momentum of Route A in the next tests. Front inertia, registration with the geometric wall, and trapped-wave occupation require independent material readings. The fueled finite-mobility Route B remains parked. Zero viscosity in the receiving bulk by itself does not choose a conversion law. [Reversible order and coherent flow](../../software/stage1_solver/reports/mouth_reversible_front.md)

**Support at a narrowing waist.** At fixed absolute wave action, two conditional sealed ducts lowered their energy by forming an interior constriction. The coherent boundary's neutral curvature law and pressure from the registered carrier leave a residual shortening force. At tight held shapes, relaxing the wave-off order field in the whole interior depleted the crystalline order even on the axis. The same preparation must specify the change in stored shear-wave energy, carrier stress and crystal strain as the order responds and opposite sides approach one another. Measure the signed total shape force there before assigning an inertia or extrapolating across radii. [Training radii and implied material test](../../software/stage1_solver/reports/mouth_lip_free_fall_gate.md)

The fixed-action, sealed-mouth test now measures both signed force directions while the scalar order and coherent fields respond. At two entrances the shared width near 0.45 still gives net inward and shortening forces; the occupied mode is concentrated downstream. A lower-level source of support must predict the missing *total* shape force in both directions. An independent test of a crystalline lip is its wave-off bending moment at two curvatures, including through-thickness stress beyond the neutral order front. The physical mouth and blind radius remain open. [Coupled waist measurement](../../software/stage1_solver/reports/mouth_coupled_waist_balance.md)

**Shape support after conversion is powered or reversible.** The conversion law must admit positive current at the physical supply. A free mouth also needs restoring work when the waist expands and the flare contracts at fixed enclosed volume and fixed length. Uniform hydrostatic pressure and a stored energy depending only on those two quantities leave the tested turning flares with negative stiffness along that displacement. Measure the spatially dependent wave stress at independently determined occupation, the normal momentum fluxes on both sides and reactions at the material ends. An added specific energy per constituent can supply conversion without supplying an isotropic pressure jump. [Reversible and downhill comparison](../../software/stage1_solver/reports/mouth_conversion_routes.md)

At fixed radius and intake, the companion-directed first moment of the flat feeding pressure is positive in the paired-minus-isolated archive. For a weakly bowed, constant-tension surface with a fixed spherical edge, the dipole of the change in inward rim slope measures the first moment of the actual local pressure jump divided by tension and projected opening volume. This permits a test for bulk cancellation from the rim alone. Use the full angular slope, independently measure the tension and baseline slope, and account for all other momentum transfers when relating this pressure contribution to free fall. [First-moment calculation, controls and limits](../../software/stage1_solver/reports/mouth_front_pressure_dipole.md)

At the curved lip impose the flux relative to the local surface. Extending the saved flat current there gives outward crossing on a positively bowed front. A geometric deformation of the current preserves inward flux but violates the neutral carrier’s circulation condition in the tested paired fields. Determine the surface phase together with the upstream and downstream velocities, check the pressure jump on that same surface, and verify local one-way crossing separately from integrated intake. These are mouth-level constraints on the 0PN boundary law. [Responding rim and local phase condition](../../software/stage1_solver/reports/mouth_responding_rim.md)

Allowing both reservoirs to adjust on one rounded mouth supplies a single energetically stable phase under fixed bath values when the longitudinal carrier passes coherently through the conversion region. In the radius-one trial at intake 0.25 the flow crosses inward over the sampled interior and turns toward the bulk at the tangent lip. The homogeneous carrier’s pressure is continuous. For a conversion surface obeying both \(-K\Delta h=q\) and \(-\sigma\Delta h=p_f-p_b\), determine independently what generates the normal pressure difference \((\sigma/K)q\). Its integral is 5/24 for the trial coefficients. This local load must be accounted for before using a responding front to infer a two-mouth force. [Isolated curved-mouth result and checks](../../software/stage1_solver/reports/mouth_curved_phase.md)

At the same intake, moving the flare farther along the fourth spatial direction changes the phase offset needed to carry the current. Establish where the minimum radius lies and how the order surface attaches. For a single smooth phase subject only to barotropic, quantum and scalar potential forces plus reversible order force, the ratio of order-energy derivative to constituent density must be constant along each connected regular order level. A level with grazing flow at the edge and inward crossing inside violates that requirement under passive material relaxation. Determine the order transport, additional stress or drive that provides the missing reaction and power if the actual mouth has this attachment. [Stress, phase and two flared geometries](../../software/stage1_solver/reports/mouth_order_and_flare.md)

Treat the observer’s projected opening, the counting cut, the material’s loss of shear order and the standing wave’s escape as separate observables of one throat. With zero wall preference and passive order relaxation, the same imposed phase flow has both a near-waist and a deeper half-order surface. On each fixed carrier both are stable against the tested small axisymmetric order variations. A shear probe between them can distinguish the preparations. Compare transverse inertia, standing-mode retention and leakage under common boundary conditions and identify the momentum reaction that allows the longitudinal carrier to remain coherent during conversion. [Two depths on one throat](../../software/stage1_solver/reports/mouth_free_order_depth.md)

For a wave to help hold a throat open, test where its energy oscillates, how much returns through the opening, and how much travels into the receiving material. A finite-time release on one held-open flared throat distinguishes the two order depths at the same intake. The deeper state accepts much more strain energy for the same imposed displacement. Much of that energy enters a weakly ordered layer where the trial wave speed becomes small. Independently determine the transverse inertia through conversion and whether the draining material transports transverse momentum. A resonance that survives these tests can then be evaluated for its load on the wall and its contribution to 0PN inertia. [Wave retention and escape](../../software/stage1_solver/reports/mouth_order_wave_escape.md)

Require a calibration that distinguishes material inertia from stress relaxation through conversion. In a homogeneous quarter-order patch two passive laws with identical complex wavenumber predict different forces for the same displacement. Measure complex local traction and strain at matched frequencies and intake, or track the separately identified material velocities and momentum exchange. Determine whether the measurement resolves a locally uniform patch before the order changes appreciably. A choice of kinetic law fixed this way can enter a test of the standing wave and its mechanical reaction. [Mixed-order shear test](../../software/stage1_solver/reports/mouth_mixed_order_ripple.md)

In the converting layer, discriminate between loss of existing shear-bearing contacts and loss accompanied by formation of new, relaxed contacts. At matched local order, measure the surviving fraction of entry contacts among those still present and the traction after a held material-frame strain step. For equal-stiffness contacts that follow the imposed shear rigidly, with unbiased removal and stress-free birth, these two normalized readings must agree. A shared formation rate per existing contact that depends only on order also bounds the change between two preparations using the time each parcel spends at each order level. Check that bound before carrying a reference-renewal law into the standing-wave calculation. [Single-mouth test and bounds](../../software/stage1_solver/reports/mouth_contact_age.md)

Test whether load-bearing features follow the independently resolved carrier deformation and keep a common shear reference during uninterrupted overlap. For equal initial row spacing and pitch, this mechanism predicts high retention on the near and deep order paths at intake 0.25. An independently sliding structure or reference renewal at fixed occupancy can change the predicted force. Measure the relative motion in units of the feature pitch and test a reference-sensitive material signal or local incremental traction on the same parcel. The native scalar order perturbation alone cannot supply a contact lifetime. [Four-dimensional carrier and material-reference tests](../../software/stage1_solver/reports/mouth_contact_observation.md)

Measure the short-time and relaxed shear stiffness at fixed local order, and track the attachment displacement across a break. In a passive material with mobile attachment points, load can relax within uninterrupted contacts or reappear in a freshly formed link as its support returns to rest. At intake 0.25, matched-order trial materials with identical unforced crossing histories and different attachment mobility transmit different held shear. The physical substrate’s relaxation law and absolute stiffness are independent inputs for the throat load. [Unforced crossing and opposite-face force](../../software/stage1_solver/reports/mouth_contact_crossing_force.md)

Distinguish relaxation of attached features from displacement recovery during a break. With unchanged support stiffness and drag, the detached time equals the intact time divided by the relaxed force fraction. Check a second detached duration with the same loading history or account for differing times spent attached. For a fixed-reference, order-dependent quadratic shear energy, three matched onset order rates under zero and opposite local shear increments can reveal the signed strain. A resolved curvature and an independent material response are essential. Absolute force also needs calibrated stiffness and the relation of the sensed strain to transmitted load. [Gap-duration experiment and signed order readout](../../software/stage1_solver/reports/mouth_detached_recovery.md)

For a weak displacement-driven contact, measure the intact reaction over the whole duration of a prospective break and reload. If the same linear material behind the contact persists during detachment, a zero-force boundary during the gap predicts its reconnection displacement and later reaction. Test this at two gap durations against an independent force or a calibrated order-rate readout. A change in neighbor drag restricted to the gap is invisible in every intact reading. The coupled patch provides one material realization and a failure control at intake 0.25. [Reconnection force generated by neighboring material](../../software/stage1_solver/reports/mouth_neighbor_patch.md)

At matched radius and intake, determine the conversion depth, tangential momentum passed to the rim, and resistance of the surrounding material to sideways motion. A contact that breaks and reforms can transmit a companion-directed pull. In the deeper preparation a freely translating support also gives a positive force requirement during slow acceleration. Every nonnegative constant support resistance in this trial leaves steady towing drag. Test the independent receiving-side force under opposite tows and held incident flows, accounting for momentum exported and forces on other supports. Measure the total low-frequency inertia independently before predicting a release or testing radius independence. [Two separations and slow acceleration](../../software/stage1_solver/reports/mouth_pair_neighbor_impulse.md)

A coherent longitudinal carrier joining the three-dimensional feed and four-dimensional receiver supplies independent held and accelerating reactions on a wall whose entrance and flare share one center. With freely measured phase-driven intake, the weak two-mouth energy gives a leading reciprocal attraction proportional to the two intakes. For universal free fall, establish from the material which wall portions move together, how much mass and stress they bring, and why the same ambient state selects the required intake-to-total-inertia ratio across radii. A common axial depth and phase difference give approximate scaling in the tested range and a measurable departure at larger radius. To test the mechanical attachment, independently measure entrance and flare displacement and local restoring force in the same preparation. A lower law for shear loss, trapped waves or charge must preserve that measured transmission and its momentum budget. [Conditional 0PN mechanism](../../software/stage1_solver/reports/mouth_carrier_0pn.md) · [Radius and attachment tests](../../software/stage1_solver/reports/mouth_carrier_0pn_methods.md#how-much-load-must-travel-through-the-attachment)

At common weak incident bias, determine the axial force conjugate to lengthening an open throat and the physical restoring load. For the tested rounded mouths, constant positive excess energy per confined volume and wall tension cannot, in any nonnegative mixture made stationary at radius 1, also hold the resolved depths required for equal weak intake per resting carrier inertia at radii 2 and 4. In the slender limit an excess volume cost has the required leading scaling. Independently fix the axial force, occupation and radius dependence of any trapped mode or other material channel before attributing the missing load to it. [Axial work, equilibrium and inverse test](../../software/stage1_solver/reports/mouth_carrier_depth_selection.md)

A trapped wave can transmit a sideways restoring force between entrance and flare and add translational inertia through the same material coupling. A conservative scalar example does both at radius 1. Its bulk confinement contrast is too small for a bound mode at radius 0.5. Equal wave action and equal prepared energy per throat volume both miss the carrier-only axial force requirement across the tested sizes. Determine the actual transverse mode’s confinement, loss and stored action across radii independently of gravity, then find the occupied shape and its total moving inertia at that same intake. [Wave loading, sideways response and size limits](../../software/stage1_solver/reports/mouth_wave_flare_work.md)

A standing wave that contributes to 0PN inertia must persist against radiation onto the feeding brane as well as transmission into material losing shear order. At radius 1, a local exponential continuation fully reflects a plane wave at 12.2 if transverse inertia falls with shear order. A six-cycle torque deeper in the same curved preparation still returned 91.37–96.61% of its input work to the feeding brane over times 44–48 across three downstream terminations or positions. With constant transverse inertia, the weak layer held 80.30% over times 14–18. Fix the material kinetic law, identify the interior pattern and its coupling to outgoing brane channels, and establish a long-lived occupation or a source of compensating power before assigning it a gravity load. [Native shear release and radius tests](../../software/stage1_solver/reports/mouth_native_shear_ringdown.md)

At the same intake 0.25, a degree-four torsional pattern with inertia proportional to shear order kept 98.33–98.37% of its initial energy in ordered throat material through the 194–198 window. After a localized external torque, its early frequency predicted the motion at a separate interior patch there to 0.224% relative error (fit at 60–100, prediction at 170–198). The feeding brane carried little net energy; the downstream loss depended strongly on the continuation. For the degree-one low-frequency torque, changing the inertia remaining at zero order from zero to one-sixteenth reduced ordered-throat retention from 97.70% to 1.68% at times 44–48. Fix the local shear speed through conversion independently. Establish the surviving wave’s incremental sideways force per stored energy and what selects its action in a forming or draining mouth. A scalar pressure and stationary symmetric intake do no direct work on this angular motion. [Wave selection and material test](../../software/stage1_solver/reports/mouth_shear_pattern_release.md)

At radius 1 and intake 0.25, an occupied native torsional mode with its axis parallel to a relative sideways shift can assist or oppose the mechanical alignment of entrance and flare. In the two tested displacement shapes, degrees two through four give negative directional wave stiffness throughout. With degree one, a prescribed bend confined to the deeper region restores alignment. Letting the same total offset redistribute into the collar gives a negative-stiffness direction. An independent wave-empty wall with area energy \(\sigma A\) makes the two-shape directional work positive for all directions when its tension per stored wave energy exceeds the respective marginal values (0.02126 for degree four and 0.06415 for degree one, rounded, inverse wall-area units). Test the actual incremental stiffness of the unoccupied draining material and the occupied mode’s energy without fitting them to common acceleration. Distinguish the wall’s baseline stress from a tension already attributed to the wave. [Bend, both modeled reactions and material threshold](../../software/stage1_solver/reports/mouth_native_wave_bend.md)

At that same radius and intake, the relaxed draining carrier amplifies every small combination of the two bend shapes even with zero wave energy. Opposite finite displacements, the work including the reservoirs, and the direct pressure on the wall give the same sign. If uniform geometric wall tension alone opposes this carrier, it must exceed 0.010237878678124177 in model energy per unit wall area. With an occupied wave, add the carrier and wave stiffness matrices before testing positivity against wall tension. The conversion-layer stress depends on the selected depth and its absolute tension. Its back-reaction on the flow and the wall’s own strain response require independent inputs before this becomes a complete stiffness or a free-fall law. [Fixed-intake force measurement](../../software/stage1_solver/reports/mouth_empty_bend_stiffness.md)

For the same two bends at radius 1 and intake 0.25, the order-weighted four-dimensional curl elasticity can supply static alignment under a specified tangential attachment. The later-opening flare with deep conversion restores both bends at unit ordered modulus even with a finite ordered-wall spring of unit strength. The original deep flare at the same modulus retains an amplifying bend even with full attachment. Both shallow preparations fail at this modulus, and normal-only contact admits zero shear work. Measure sustained material traction per curl strain and wall traction per tangential slip independently, and establish which order depth forms. No new radiative mode or geometric wall tension is required by this particular static mechanism. The physical attachment, reciprocal conversion stress, wave loading and freely translating response remain open. [Relaxed strain and attachment comparison](../../software/stage1_solver/reports/mouth_bend_static_shear.md)

The static alignment on the later-opening deep flare can persist under replacement if the contact reference moves downstream. In a quasi-static held-bend test at intake 0.25, a marker carried along the wall reached the end of the deeper bend in 26.0591 model time units. At a renewal time equal to that transit, fixed sites lost a restoring direction while transported contacts continued to oppose both bends. The steady directional-work crossing occurred near renewal times of 0.08403 and 0.08720 transits for straight-reference and zero-traction entry respectively, conditional on carrier-speed transport, the common unit modulus and attachment, and the tested geometry. Measure reference transport speed and relaxation independently, and identify the energy and momentum partner before using a sustained force for free translation. [Held force under material replacement](../../software/stage1_solver/reports/mouth_bend_renewal.md)

At the straight intake of 0.25, test the source of work for material contacts during a sustained bend. A +0.01 deep displacement requires 5.711 × 10⁻⁸ model power under the straight-entry traveling-reference law. Its mean spring energy entering across the waist accounts for at most 14.2%. With traction-free entry this contribution is zero. At the original equal-frequency drive the scalar carrier passes sideways momentum to the curved wall and supplies zero net power from the reservoirs. At matched intake, measure the bent-minus-straight change in energy carried into and out of the conversion region, the change in pre-existing dissipation, and the force transmitted to material along its path. Identify the physical recipient of the opposite momentum before assigning the restoring force to a freely falling mouth. [Held-bend energy and momentum](../../software/stage1_solver/reports/mouth_bend_flow_reaction.md)

At fixed intake 0.25, bending the deeply ordered later-opening mouth lowers the loss in a passive conversion layer. To pay the full contact renewal in the formal weak-bend calculation requires transferring about 16% of the saving for straight-prepared contacts or 15% with traction-free entry. The current carrier supplies zero net reservoir power. Measure the independent work maintaining the straight conversion, the fraction transferred to contacts under bending, and the opposite force on material. Track stored energy and heat separately. The order-dependent shear and wall forces change the conversion loss at the same quadratic order. Enforce physical order bounds downstream before interpreting a finite-bend result with this feedback. [Conversion budget and transfer test](../../software/stage1_solver/reports/mouth_bend_conversion_power.md)

Distinguish work entering contacts from reduced loss of their stored energy. Supplying an independent trial energy per constituent at fixed intake moved the conversion layer and reduced passive contact losses, with signed work still leaving the contacts. A conversion-driven change of relaxed length can transfer net work in a partially converted patch when its reaction is charged back to the conversion supply. Calibrate that change from marked contacts with fixed slip, opposite prestrains and measured transit time. Count both end stress-work rates, ordinary relaxation, and energy carried into and out of the test. Verify how much work reaches ordered material or moves the mouth before shear resistance vanishes. The sampled curved-mouth preparations fell short of paying for full renewal. [Energy-transfer measurement](../../software/stage1_solver/reports/mouth_conversion_contact_transfer.md)

Resolve contact identity and both endpoints during conversion. An extended junction with two fixed relaxed offsets can shift its combined force-free position as its two attachments soften at different rates. This can emit previously stored elastic energy into a connected ordered medium. To count it as a sustained power source, identify how each new junction acquires its initial strain. Measure force versus slip on matched copies to locate the unloaded position, record breakup, and integrate force times ordered-end velocity on the unperturbed path. [Carrier traversal, passive comparison and work bound](../../software/stage1_solver/reports/mouth_contact_tracking.md)

In a single tagged junction the fraction of entry strain energy delivered to an ordered receiver changed from 4.58643% to 25.08964% under the same two breakup events when the stiffness of each live attachment was held constant. Calibrate the force-versus-slip slope on surviving contacts at several order fractions and the capture width of the same features. A fraction of engaged contacts and a stiffness per surviving contact affect the macroscopic modulus independently. Track energy released at rupture separately from work done on ordered material before rupture. [Tagged 0PN measurement](../../software/stage1_solver/reports/mouth_tagged_contact.md)

Record the instantaneous transverse force change for both signs of a weak step on a surviving site at known registration. Under common affine displacement a patch stiffness proportional to order and a live fraction proportional to order fix only the mean response of surviving sites. A central binding well gives a softer off-center attachment as its capture edge approaches while preserving that same patch scaling. Compare the transverse force slope of centered and offset contacts at matched local order, and test the off-center slope against the registry-restoring force divided by the signed offset with a minus sign. [Surviving-attachment comparison](../../software/stage1_solver/reports/mouth_surviving_attachment.md)

Calibrate attachment forces with two displacement components present. At half shear order a site held one eighth of a pitch off center exerted 75% of the centered sideways force in a central binding well and 100% under separate forces in the two directions. Both laws reproduce the same uniaxial readings. Measure the actual offset during conversion, the transverse restoring force, and any change in the force holding the offset under opposite transverse steps. In a reversible well the latter change integrated along the offset predicts the difference in sideways holding work. [Matched forces and mobile-site controls](../../software/stage1_solver/reports/mouth_matched_attachment.md)

A single-site force calibration requires independent positions of both connected features and a reaction transducer at controlled actual offset. Compare the mean of opposite sideways steps to a zero-step run at equal elapsed time and retain the odd component separately. An initially prestrained continuum patch with evolving shear order has its own even stress-rate change under these steps. Account for that and other surrounding stresses in the transducer calibration. Two single-axis force curves and the mean patch stiffness leave the tagged force response under both displacements undetermined. [Contact-free brane and identifiability checks](../../software/stage1_solver/reports/mouth_offset_brane_readout.md)

With both pad centers and gauges resolved, a registered pair at one-eighth pitch offset gives opposite changes in its holding force for two different crossed-strand layouts under the same sideways steps. Matched brane-only readings isolate the pair in this frozen-order trial. A substrate model should fix the positions of its anchors and its axial force law independently. For these layouts, the zero-step force already distinguishes the geometries, allowing the much smaller sideways change to serve as a prediction. [Two-site force calibration](../../software/stage1_solver/reports/mouth_two_site_clamp.md)

In the half-head radius-one preparation at intake 0.25, weak angular disturbances of the native scalar order decay with the carrier held fixed. This excludes a linear instability of that field as the source of attachment sites in this state. Resolve the neutral material’s pitch, persistent connected sites and force scale independently before using a zero-step holding force to select a layout for the mouth. [Independent native anchor test](../../software/stage1_solver/reports/mouth_anchor_order_spectrum.md)

A neutral two-layer substrate with short-range binding can form persistent site pairs and predict the one-eighth-pitch holding force from an independently measured neutral-pair tension curve once the angular dependence is specified. Two interactions with identical normal tension differed in the sideways reaction, and strict pair survival varied by preparation. Calibrate a weak tangential displacement at a known normal preload, separate the direct contact from surrounding material, and test its survival and force scale as order is lost before transferring this input to 0PN. [Neutral patch formation and clamp](../../software/stage1_solver/reports/neutral_patch_sites.md)

Carrying the same neutral pair identities along the 0PN flow distinguishes continued registration from force transmission. With cross-layer bond strength proportional to local order, the two tagged pairs stayed registered to 1% order in four thermal histories even as their binding fell below the thermal scale. A trial connection to the held waist transmitted 4.88% or 4.05% of the imposed sideways impulse during transit. Determine the conversion dependence of the binding, the grip strength and shear travel time independently. To connect the measurement to free fall, obtain the lateral loading from the companion flow at matched intake and resolve the isolated-mouth background. [Pair survival and transmitted rim force](../../software/stage1_solver/reports/neutral_patch_conversion.md)

The measured feeding-field difference at fixed intake supplies a weak lateral drive for the converting patch. At eight radii both trial bond laws transmit a companion-directed force to a held waist. At sixteen radii each transmits 24.93% of its eight-radius increment. Matching the phase across the feeding face and narrower waist and exposing one or both material layers are constitutive assumptions. Resolve each layer’s relative velocity and drag, the local normal current and traction at the join, and the number and angular coverage of active contacts. These quantities set the whole-mouth strength of this channel. [Matched companion and isolated rim readings](../../software/stage1_solver/reports/neutral_patch_companion_rim.md)

An independent test of flow penetration should read the signed traction on the wall-facing sites and their mean tangential motion relative to the local carrier during conversion. At fixed grip and intake both tested bond laws predict a backward change in bath force on that layer when only the neighboring layer is exposed, and a forward change when both are exposed. At eight radii the respective lags relative to the same companion flow were 73.83–76.37% and 9.81–10.22%. In the trial the force follows independently from the outgoing grip load minus the interlayer force on those sites. A physical calibration requires material motion and momentum exchange resolved independently of the imposed exposure. [Matched layer slip and signed drag](../../software/stage1_solver/reports/neutral_patch_layer_slip.md)

A two-way local shear model with identical friction on all sites measures the force on the wall-facing layer from incoming shear minus local fluid momentum storage. At the same intake and grip, the force changed sign with the material behind the layer and with how much lateral motion preceded the transit. The incremental force reaching the held rim remained companion-directed in all twelve fine-grid comparisons. Determine the time-dependent ratio of traction to velocity behind the sites and the incident shear history independently of the paired force. A one-mouth weak pulse, with separate fluid-storage and boundary readings, would constrain both its propagation and return. [Resolved site–carrier momentum exchange](../../software/stage1_solver/reports/neutral_patch_reciprocal_flow.md)

Use the decay of the shear pulse as well as its early rise when calibrating the material behind the sites. Two interior velocity and traction readings distinguished finite resistance, free slip, a fixed backing and extended viscous diffusion in one-mouth trials. Four further drive shapes and an initial-flow release predicted the signed held-rim and wall-layer forces across 32 companion comparisons at the same intake within these material models. An independent stress law and density profile across the loss of order are still needed to choose a physical response. A whole-mouth momentum and inertia budget remains necessary for free fall. [Blind one-mouth calibration and prospective force tests](../../software/stage1_solver/reports/neutral_patch_shear_pulse.md)

At the same intake the reversible sideways stiffness at fixed material reference follows shear order. It fell by a factor of 95 on the saved converting parcel. Two open elastic receivers gave late signed traction errors of 21.58% and 24.81% even for the closest prior result, the deep viscous preparation. On the free-backed pulse the interval retained positive strain while its reconstructed mean stress and both recorded face tractions reversed. Measure the force retained after a small strain step at several orders, along with the moving mass or pulse speed and the momentum carried by converting material. These inputs would connect the measured companion flow to a physical receiving-side force. [One-mouth order and shear response](../../software/stage1_solver/reports/neutral_patch_order_pulse.md)

At fixed shear order, compare the internal force during a small held strain with the response to a separately driven pulse. The two candidate bond laws gave different force histories on the same four gate levels. For the layered material at 1% order, 64 held-strain measurements retained 56.2% of their initial mean force late in the hold. They predicted the independent pulse with 16.01% full-record error and too much impulse after release. In the rested patch, restricting its overall orientation suppressed most of the central material’s force loss. Determine how strongly surrounding ordered material restricts that rotation, and how the force after release depends on the time since arrival. These are independent material inputs to the companion-force calculation. [Fixed-order hold and independent pulse test](../../software/stage1_solver/reports/neutral_patch_strain_hold.md)

Measure the angle and torque of a converting patch relative to the immediately adjoining material, distinguishing within-layer turning from relative motion across layers. In direct thermal torque drives, the relative turn after one transit grew as order fell in both 24-site models. The composite turn and within-layer turns have different cold restoring stiffnesses. Under a bonded planar-sheet hypothesis, surrounding material can support extra interlayer force over a transit. A static curl-free planar turn and a fully surrounded three-dimensional turn have different torques at infinity. Determine which faces stay shear-connected, the extent of tangential slip and how torque returns to the ordered brane before assigning a whole-mouth response. [Torque from neighboring material](../../software/stage1_solver/reports/neutral_patch_neighbor_turn.md)



### How much energy can an occupied standing mode release?

Finite trapping can admit an angular dipole at the reference mouth size. Determine the wave-only radial and axial load and resolve the phase of orthogonal dipole patterns. Two force readings distinguish two mode-family actions when their load coefficients are independent. A signed local spectrum also requires a coupling to the angular field. In a model with positive rim springs calibrated to the stated test values, every fixed frequency-ranked mixture of bound modes falls short of one exterior unit on both finite lattices at six tested radii when the scaled collar cost is included. Independent inputs still needed are the patch and neighboring-mode splittings, the wave-only loads in the driven drain, and any extra work of conversion or attachment during a change in flux. [Standing-wave budget and scope](../../software/stage1_solver/reports/mouth_standing_wave_budget.md)

### Can conversion supply signed work to a stationary mouth's angular field?

Measure the internal excitation and translational energy entering and leaving a patch at matched intake. Independently measure torque as the angular field moves. A steady discharge can give a mirror-opposite change in outgoing excitation even when the net quasistatic work on every full turn vanishes. In the equal-site-energy model tested here, the greatest available conversion and standing-wave work together leave an energy deficit for nonzero exterior flux at each of six sampled radii on two finite lattices. [Conversion energy and bounds](../../software/stage1_solver/reports/mouth_conversion_twist_work.md)

A proposed shift of the mean site energy with deposited excitation permits net work over a driven turn. Compare the transmission phase at an isolated site before and after a controlled change in deposition using the same transit time, a coherent reference and both internal channels. The existing resonant contact changes the interference amplitude but gives zero common phase at weak dwell even for a retained receiver. An added shift gives a signed phase in that same experiment. Establish the actual magnitude, residence time and relaxation of deposited energy, its sustaining energy source, and a mechanism to form and retain one unit on release. [Guarded evidence and assumptions](./_measurements/mouth_conversion_twist_work.json)

### Which rim material can affect a later parcel?

Identify a degree of freedom that stays attached as matter crosses the mouth. Determine the first two excitation energies, the exchange amplitudes with a passing parcel, and the spectral weight transferred into each accessible bulk channel at that frequency. Compare the time over which a localized preparation survives with the interval between visits to the same patch. A finite-band example supports a persistent mode. A second channel at its energy removes permanent binding at the tested couplings. The spectrum and arrival rate together predict what remains locally when the next parcel arrives. [Storage and leakage tests](../../software/stage1_solver/reports/mouth_rim_retained_excitation.md)

Compare a weak probe of a loaded patch with the same patch empty using the same bypass reference. A localized excitation with finite capacity and detuned exchange can shift the mean of both internal phases even when the deposited parcel has left. The leading mean energy shift cancels for a harmonic receiver with oscillator exchange matrix elements. The tested scalar memory and its phase are the same on the matched mirror. Predict opposite exterior fluxes from the spatially arranged material and count loading, storage and angular work separately. Connecting this candidate to the archived standing wave requires identifying the material observable that gives these measured transitions. [Equations, energy budgets and working evidence](./_measurements/mouth_rim_retained_excitation.json)

### How can a trapped mode acquire distinguishable addition energies?

With a separate monopole maintaining the opening, the scalar well supports an initially empty dipole at the archived radius. Its energy is linear in action at fixed shape. Relaxing the depth lowers the second incremental cost on both the reflecting-end comparison and the smooth well connected to the scalar exterior. For the smooth well the predicted angular-frequency separation is 0.00406283876392 at the trial action quantum 0.03. Establish the wave action scale and the effective wall constraint on the time scale of a material probe. A global dipole and its mirror have identical addition energies. Coupling to a local carrier and to a mirror-opposite exterior flux require separate evidence. [Trapped-mode spectrum](../../software/stage1_solver/reports/mouth_trapped_wave_spectrum.md)

At a fixed radius and matched intake, measure the radial displacement with a dipole angular pattern under a calibrated force at the difference and sum frequencies of the support and stored modes. The displacement per force, including phase, enters the energy transfer. Also measure axial motion at twice the support frequency and outgoing scalar energy in the upper sideband. Determine the intrinsic line widths and phase noise, including bulk longitudinal channels. Two passive bend responses give decay rates on opposite sides of the same static spacing. Storing energy in this mode adds terms to the force and momentum budget of a moving throat, but does not select how much of the external pull accelerates the mouth. [Evidence and limits](./_measurements/mouth_trapped_wave_spectrum.json)

For a coherently supported lossless dipole coupled to one bend, the two exact combination-frequency wall readings carry one independent complex response. Applying both drives in phase prepares a wave with zero wall displacement for a subsequent ringdown. At the same stationary current, two passive bend laws give two or one peaks under the same specified readout. To determine the material linewidth, measure the force from the draining carrier on a prescribed angular displacement, include work by the support, and check intrinsic wave loss and dephasing. A stationary flux calibration alone fixes none of these. [Driven bend and ringdown](../../software/stage1_solver/reports/mouth_driven_bend.md)

A moving reflecting carrier wall supplies a calculable contribution to bend damping. Under prescribed dipolar motion, integrate the fluid traction against wall velocity and check the outgoing perturbation-energy flux at the upstream and downstream cuts. An angular perturbation preserves the integrated intake to first order. At weak steady intake, most of the tested wall radiation also occurs at rest. Inferring a linewidth requires the radiative response of the actual brane-side inlet and receiving bulk, the bend’s static stiffness and inertia, and other wave-loss and phase-noise channels. The finite-width one-arm test retained two peaks under the earlier selective readout. [Carrier radiation and support power](../../software/stage1_solver/reports/mouth_carrier_bend_work.md)

For one receiving arm, outgoing carrier waves can be matched across both apertures to a radially extended brane and a four-dimensional bulk half-space. Reflection depends on the mixture of transverse profiles at the opening; locally decaying fields also transmit power through interference. In the prescribed bend calculation, the brane and bulk usually pass most of the incident traveling power at either drive. A pure incident radial mode can undergo much stronger reflection. The emitted energy depends on slab thickness, carrier law and face conditions. With those specified, two peaks remain compatible with longitudinal radiation. Establish the inertia and stiffness of the same bend and its remaining loss and readout channels independently to predict the physical linewidth. [Open-reservoir test](../../software/stage1_solver/reports/mouth_open_reservoirs.md)

The radial breathing inertia constrains a different integral of the material mass than the sideways bend. At the same intake, two positive mass profiles with identical radial normalization, side tension, scalar support and carrier radiation produce opposite answers to the two-excitation resolution test. For the stated selective readout and test stiffness the changeover is near projected wall mass 13. A wall law must determine the depth-weighted inertia and restoring force independently of the spectrum. Measure the same bend’s static force and low-frequency force slope with the known field responses removed. Its fixed endpoints preclude identifying this inertia with the mass of a freely translating mouth. [Bend projection and linewidths](../../software/stage1_solver/reports/mouth_bend_coefficients.md)

The static and low-frequency field-force subtractions for this bend are now available. To extract the wall inertia independently at the same intake, the wall-bearing medium must also supply its total actuator reaction or the momentum entrained along the bend. Under a local elastic and inertial wall law, doubling the weak drive frequency quadruples the in-phase residual force drop at fixed displacement. The two candidate wall profiles pass this test with different total forces and opposite peak-count verdicts. [Held-bend force calculation](../../software/stage1_solver/reports/mouth_bend_holding_force.md)

At fixed intake the resolved longitudinal carrier has a calculable transverse velocity profile and contributes to the already included carrier inertia. The additional wall inertia depends on which brane material follows the normal bend. With the previous radial normalization and tension, any local attached-mass density that decreases monotonically into the arm falls short of the measured high-inertia two-peak threshold. Attachment proportional to the computed trapped-wave pressure gives one peak. Independent depth-resolved wall momentum or a material law for that momentum is needed to identify the realized profile. [Depth-resolved bend test](../../software/stage1_solver/reports/mouth_bend_depth_participation.md)

For this fixed-end bend, two equal-width momentum readings, one at the entrance and one at mid-depth, leave enough freedom in the wall-mass distribution to give either peak count at the same intake and radial mass. A sufficient independent input for the stipulated local, normally moving wall is its velocity-phase normal momentum, projected onto the imposed motion and integrated over both halves. The solver supplies the surrounding carrier’s contribution separately. This test leaves the inertia for free translation at 0PN open. [Wall momentum and peak-count comparison](../../software/stage1_solver/reports/mouth_bend_band_momentum.md)

For a thin reversible order front that supplies the side tension, an independent measurement of its intrinsic ripple speed gives a route to normal mass per wall area, equal to tension divided by that speed squared. Uniform properties with the archived radial inventory give one peak at the reference intake. At the same stiffness the fine-grid two-peak branch requires a largest-to-smallest mass-per-depth ratio above 5.45779 for distinct peaks, even at optimal placement. The live order dynamics, identification with the wall, and transmission of material inertia into the bend all require independent support. [Full-depth wall material test](../../software/stage1_solver/reports/mouth_bend_full_material.md)

A measured end-to-end wall travel time and speed contrast can both agree across preparations with opposite peak counts. To use ripples for the bend inertia, retain spatially resolved wall motion and conjugate force, and remove the complex fluid and supporting-wave loads on the same motion and frequencies. A distant quarter of the wall can carry a force even when clamped. Under the local order-front hypothesis the required material moment weights the inverse squared intrinsic speed by the squared bend shape. [Local ripple test](../../software/stage1_solver/reports/mouth_wall_ripple.md)
