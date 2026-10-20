# Substrate requirements — what the force sectors oblige the unbuilt steps to deliver

**Status: the light sector is closed for now (user decision, 2026-10-19). Passes 1 (S9, S10), 2 (S11,
S11b-A/B, S11c PARTIAL; accepted `22be0137`), O2 (accepted `72866fcf`) and S9b (accepted `ff3bac63`) are
complete. The light-closing pass (2026-10-19) adds the oracle comparison and the throat energy split; its
review is pending.**
Eighteen entries, all OPEN. Pass 2 was populated 2026-10-05 from the kept results and their stated
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
| **S1** | `R-S1-02` | the substructure's shear response as a function of `χ_B`, both phases | S9, O2 conditional premises, S9b brane optical shear only |
| **S1** *(register inference)* | `R-S1-03` | whether the substructure's microdynamics is time-reversible | S11b-B / unified S11b |
| **S1.5** | `R-S1.5-01` | the mode content of the linearised GNLS, including the bulk acoustic branch | S9, S11, S11b-A/B, S9b Part C bulk acoustic input only |
| **S6** | `R-S6-01` | the brane's compression modulus `B_comp` | S11 |
| **S6** | `R-S6-02` | the first variation of the substrate action at `v₀ = 0` | S9, S10, S11, S11b-A/B, S11c uniform |
| **S6** | `R-S1-01` | the brane's spatial dimension `D_brane` ⚠ *id kept; target corrected to S6* | S10, S11, S9b coordinate dimension only |
| **S7** | `R-S7-01` | the slab-width flat direction under a non-uniform width | S11, unified S11b |
| **S8** | `R-S8-01` | the form of the brane's quadratic stiffness functional | S9, S10, S11, S11b-B, unified S11b, S11c uniform, S9b supplied transverse input |
| **S8** | `R-S8-02` | the quadratic operator on `u` **and** `h` together | S10, S11 |
| **S8** | `R-S8-03` | the **sign** of the physical transverse stiffness | S9, S10, S11, unified S11b, S11c uniform, S9b local optical sign/ratio gate |
| **S8** | `R-S8-04` | what carries the brane's **internal angular momentum** | S9, S10, S11, O2 carried obligation |
| **S8** | `R-S8-05` | the **frame** the brane's rotational stiffness is measured against | S9, S10, S11, O2 carried obligation |
| **S8** *(register inference for original entry)* | `R-S8-06` | material displacement and original quadratic slab inertia; O2 adds the material identification only | S11, S11b-B, O2 and S9b identification |
| **S12** *(register inference; no live owner named)* | `R-S12-01` | a reservoir and power budget if a non-passive interface response is adopted | S11b-B / unified S11b, O2 carried condition |
| **S12** | `R-S12-02` | background drain/return and separate boundary data | S11b-A/B, S11c uniform, O2 |
| **S21** | `R-S21-01` | conditional deflection and radar-slope profile requirements, with Part C rewrites and implied exchange | S9b A–C; flow-only oracle comparison (evidence only) |
| **S22**, with Q2 and S12 | `R-S22-01` | how the energy a throat releases splits among light, the extra longitudinal wave and the bulk | S11 departure; photon survey E32–E34, E38 |
| **no owner named** | `R-O2-01` | the physical net-power supplier and its stated budget, if a conditional closure requires net supply | O2 |

⚠ Entries below are in the order they were found, ⛔ not in step order. All **eighteen** are **OPEN**.
⭐ `R-S8-01`, `-03`, `-04`, `-05` are one family: the stiffness functional's **form**, its **sign**, its
**mechanical admissibility**, and its **reference frame**. ⛔ Delivering the form does not deliver the
other three.
⚠ `R-S1-01`'s id encodes its original target; the id is stable and the **target** field governs.

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

- **source** S9b (coordinate dimension only; `steps/S9b_light_bending.md`, §§2–4/8); S10 (`steps/S10_two_transverse_photons.md`); S11
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


### R-S1-02 — the substructure's shear response, as a function of the order parameter

- **source** S9b (brane optical shear only; `steps/S9b_light_bending.md`, §§2–4/8); O2 (`steps/O2_steady_brane_balance.md`, register handoff); S9 (`steps/S9_light_requires_shear.md`) · **target** S1 · **status** OPEN
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

### R-S12-02 — background drain/return and separate boundary data

- **source** O2 (`steps/O2_steady_brane_balance.md`, register handoff); S11b-A/B (`steps/S11b_interface_coupling_law.md:158–164`, background-flow limit); S11c
  (`steps/S11c_PARTIAL_CLOSEOUT.md:7`, uniform result;
  `steps/S11c_d_profile_conditioned_scattering.md:21–25`, selected uniform check)
  · **target** S12 · **status** OPEN
- **requirement** — the native drain/return functions and their separate boundary data.
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
  (`steps/S11c_PARTIAL_CLOSEOUT.md:33`).
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

---


### R-S21-01 — conditional deflection and radar-slope profile requirements

- **source** S9b (`steps/S9b_light_bending.md`, §§3–6; measurements M3–M7/M11–M13/M17–M19)
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
- **on failure** — The conditional reference-matching statement for the affected observable
  cannot be realized on that domain. The computed conditional functionals and the reported
  cross-engine coverage do not become a GR-matching result.
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
  - The photon survey, E32–E34 and E38 (`PHOTON_POLARIZATION_SURVEY.md` §4). Each is "In apparent conflict"
    under the condition that the extra longitudinal branch couples to matter and acts as an additional
    radiating mode. (C6, C7.)
- **target** S22 (the nonlinear throat solve), together with Q2 (the mouth source and its held datum) and
  S12 (the throat's order-conversion source and its boundary data). These owners were chosen by the user on
  2026-10-19; the S11 record itself names only "derived matter/interface coupling". (C8.)
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

### Light-closing pass — the oracle comparison and the extra wave · 2026-10-19; review pending

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
- the first two belong to S12 (`R-S12-02`) and to O2's carried handoffs;
- the third stays a scoped S9b exclusion;
- the fourth stays in the S9b pass's analog-gravity row.

These are handoffs, not entries. No owner is added to any other oracle check.

**Prior art.** The positive light results reproduce known prior art: MacCullagh's curl-only medium, and the
squared flow speed of the river model and of Visser's §7 acoustic flow. The flow match holds only under the
unsupplied mass, clock and speed-unit identifications. The register's
MacCullagh obligations (`R-S8-04`, `R-S8-05`) remain OPEN.


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
