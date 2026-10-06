# Substrate requirements — what the force sectors oblige the unbuilt steps to deliver

**Status: passes 1 (S9, S10) and 2 (S11, S11b-A/B, S11c PARTIAL) complete; pass-2 review pending.**
Fifteen entries, all OPEN. Pass 2 was populated 2026-10-05 from the kept results and their stated
conditions; it does not close the substrate or upgrade S11c's unresolved work. Later sectors remain to
be read as they close. The two 2026-08-07 prior-art entries retain their provenance; the second route is
recorded under Population passes.

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
| **S1** | `R-S1-02` | the substructure's shear response as a function of `χ_B`, both phases | S9 |
| **S1** *(register inference)* | `R-S1-03` | whether the substructure's microdynamics is time-reversible | S11b-B / unified S11b |
| **S1.5** | `R-S1.5-01` | the mode content of the linearised GNLS, including the bulk acoustic branch | S9, S11, S11b-A/B |
| **S6** | `R-S6-01` | the brane's compression modulus `B_comp` | S11 |
| **S6** | `R-S6-02` | the first variation of the substrate action at `v₀ = 0` | S9, S10, S11, S11b-A/B, S11c uniform |
| **S6** | `R-S1-01` | the brane's spatial dimension `D_brane` ⚠ *id kept; target corrected to S6* | S10, S11 |
| **S7** | `R-S7-01` | the slab-width flat direction under a non-uniform width | S11, unified S11b |
| **S8** | `R-S8-01` | the form of the brane's quadratic stiffness functional | S9, S10, S11, S11b-B, S11c uniform |
| **S8** | `R-S8-02` | the quadratic operator on `u` **and** `h` together | S10, S11 |
| **S8** | `R-S8-03` | the **sign** of the physical transverse stiffness | S9, S10, S11, unified S11b, S11c uniform |
| **S8** | `R-S8-04` | what carries the brane's **internal angular momentum** | S9, S10, S11 |
| **S8** | `R-S8-05` | the **frame** the brane's rotational stiffness is measured against | S9, S10, S11 |
| **S8** *(register inference)* | `R-S8-06` | the material displacement and quadratic inertia of the slab | S11, S11b-B |
| **S12** *(register inference)* | `R-S12-01` | a named reservoir and its power budget | S11b-B / unified S11b |
| **S12** | `R-S12-02` | background drain/return and separate boundary data | S11b-A/B, S11c uniform |

⚠ Entries below are in the order they were found, ⛔ not in step order. All **fifteen** are **OPEN**.
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

- **source** S10 (`steps/S10_two_transverse_photons.md`); S11
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

### R-S1-02 — the substructure's shear response, as a function of the order parameter

- **source** S9 (`steps/S9_light_requires_shear.md`) · **target** S1 · **status** OPEN
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

### R-S1.5-01 — the GNLS's own mode content

- **source** S9; S11 (`steps/S11_stray_longitudinal.md`, move 6 / `KW_ZERO_LOCUS`);
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

- **source** S9, S10; S11 (`steps/S11_stray_longitudinal.md`, moves 2–3 / FORM control);
  S11b-B (`steps/S11bB_interface_assembly.md`, energy basis / breathing stability;
  `steps/S11b_interface_coupling_law.md`, energy quotient); S11c
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

### R-S8-02 — the in-plane and out-of-plane sectors must decouple at quadratic order

- **source** S10, S11 · **target** S8 · **status** OPEN
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
- **pass-2 source** — S11's selected-action statement says its homogeneous-sector decoupling does not
  establish decoupling *"after additional allowed fields are introduced"*
  (`steps/S11_stray_longitudinal.md:75–76`), and its scope caveat explicitly excludes `hBranon` from the
  frozen-wall-width in-plane spectrum (`steps/S11_stray_longitudinal.md:279–280`). It therefore sources the
  same open obligation for the joint `u`/`h` operator.

### R-S8-03 — the SIGN of the physical transverse stiffness

- **source** S9, S10; S11 (`steps/S11_stray_longitudinal.md`, unchanged transverse branch);
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

### R-S8-04 — what carries the brane's internal angular momentum

- **source** S9, S10, S11 · **target** S8 · **status** OPEN
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

### R-S8-05 — the frame the brane's rotational stiffness is measured against

- **source** S9, S10, S11 · **target** S8 · **status** OPEN
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

- **source** S11 (`steps/S11_stray_longitudinal.md`, move 1 / finite census); S11b-B
  (`steps/S11bB_interface_assembly.md`, breathing quadratic) · **target** S8
  (**register inference; no owner named by the records**) · **status** OPEN
- **requirement** — `u` as the material displacement of the stuff whose density is `ρ_br`, and the
  quadratic kinetic form of that material and the thickness degree of freedom (B's `μ_W`).
- **on failure** — if `u` is a director rather than that displacement, S11's continuity identification
  `δρ_br = −ρ_br ∇·u` and its compression argument do not follow. Its finite census requires nonzero
  `ρ_br`; B's breathing-root interpretation also rests on the stated inertial model. A stiffness
  functional alone (`R-S8-01`) does not supply this identification or the kinetic form.
- **note** — this asks for the field identity and inertia, not numerical benchmark values of `ρ_br`
  or `μ_W`, and does not identify the thickness mode with S10's out-of-plane displacement. The cited
  records name no future owner for this object; S8 is a register inference
  (`steps/S11_stray_longitudinal.md:32–38`; `steps/S11bB_interface_assembly.md:80–85`).

### R-S12-01 — the reservoir and its power budget

- **source** S11b-B (`steps/S11bB_interface_assembly.md`, standing rule); unified S11b
  (`steps/S11b_interface_coupling_law.md`, passivity region) · **target** S12
  (**register inference**) · **status** OPEN
- **requirement** — a named reservoir and a stated power budget for any adopted non-passive interface
  coupling. The candidate reservoir is the background drain `v₀` (standard name `v_bulk_normal_0`,
  `directives/S11b_SHARED_PHYSICS.md:101–104`; `directives/S11b_unified_decisions.md:33–40`), as recorded in
  `steps/S11bB_interface_assembly.md:46–48`.
- **on failure** — the finite-memory velocity channel outside the passive region cannot be inherited
  as a physically supplied response. The computed region remains a classification, not a prohibition,
  and naming `v_bulk_normal_0` alone does not supply the missing power.
- **note** — this is the records' explicit condition for retaining that channel, not a requirement
  that the model choose it. A numerical observational bound also needs matter-to-compression coupling;
  B says that is unbuilt and reports only a structural test, so no numerical-bound requirement is added.
- **owner/target** — the records do not name S12; that target is a register inference. They hand the
  reservoir-and-budget condition to *"S11b-C"* (`steps/S11bB_interface_assembly.md:195–197`).

### R-S12-02 — background drain/return and separate boundary data

- **source** S11b-A/B (`steps/S11b_interface_coupling_law.md:158–164`, background-flow limit); S11c
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

---

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

### Pass 2 — S11, S11b-A/B and S11c PARTIAL · populated 2026-10-05; review pending

Read `steps/S11_stray_longitudinal.md`, all three S11b records (the unified record and the historical A/B
records), `steps/S11c_PARTIAL_CLOSEOUT.md`, and the a/b/c1/c2/d records in full. Applied the schema,
rest-on test and both routes below: **11 → 15 entries**, four new and six existing entries gaining
sources (`R-S7-01`, `R-S1-01`, `R-S1.5-01`, `R-S6-02`, `R-S8-01`, `R-S8-03`). New entries:
S1 `R-S1-03`; S8 `R-S8-06`; S12 `R-S12-01`/`02`.
All use schema status OPEN; CONDITIONAL is a qualification of the source result, not a status value.

S11c contributes only conditions of its kept uniform/selected-uniform results, with the closeout and d
record's review limits retained. The clean-condition packet sources no entry: the closeout marks only
uniform decoupling CONDITIONAL (`steps/S11c_PARTIAL_CLOSEOUT.md:7`), while its no-leak reading is OPEN
(`V3_STEP_PLAN.md:532–534`; `steps/S11c_d_profile_conditioned_scattering.md:64`). Its UNRESOLVED direct
mixed term, numerical loss, and OPEN nonuniform confinement/material, observable and real-throat questions
source no entry.
The unified spec marks the §1b branch prescription supplied and lists both it and the §4 closure among the
supplied objects (`directives/S11b_SHARED_PHYSICS.md:113`; `directives/S11b_SHARED_PHYSICS.md:995`); B's
record likewise includes the branch
prescription among the supplied objects (`steps/S11bB_interface_assembly.md:136–137`) and separately says
that radiated-energy direction is inherited from the supplied continuation
(`steps/S11bB_interface_assembly.md:158–159`). From that supplied status, this register infers that the
interface closure and outgoing-branch prescription have no named retirement condition and routes them to
the ansatz ledger as candidate postulates; v3 has no ansatz ledger yet.
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
directive's record-and-objection test. Onsager–Casimir is a candidate only; `R-S1-03` instead rests on B's
first-route hand-forward, *"conditional on microscopic time-reversibility"*. A's memoryless reduction is
by construction, and internal cross-engine reproductions are not prior art. The two pre-existing
prior-art entries are not re-certified by this pass. No directive/register method disagreement was found;
the obsolete rebuild schedule is replaced by this completed population pass. Record
qualifications/disagreements are reported at STOP, not adjudicated here.

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
