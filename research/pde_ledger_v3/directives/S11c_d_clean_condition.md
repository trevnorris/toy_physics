# Light leakage — the clean condition (S11c-d re-scope proposal)

**Status:** PROPOSED · **v5** (2026-10-01) · Codex-revised (v4–v5) from the orchestrator-written v1–v3 · ⛔ not
governing until review-cleared. The committed v4 baseline is `5693e861`.
- Rounds 1–4 and every accepted disposition are recorded in
  `directives/_measurements/S11c_d_clean_condition_review_disposition.md`.
- Repo citations below are reproduced verbatim, with their commands, in
  `directives/_measurements/S11c_d_clean_condition.md`.

**Plain summary.** We stop asking *"does light leak at a non-uniform slab?"* and ask instead *"under what clean
condition does it provably not leak?"* There is one candidate: **in an achiral medium, light whose motion is a
pure twist about a symmetric non-uniformity cannot enter a reflection-even face/bulk channel at linear order.** This is a standard
symmetry argument, and it gives an exact zero rather than a small number. It
contains no light/bulk speed ratio. It is the natural candidate for **light trapped in a throat** (particle
stability), which is the case that needs an exact zero. It does **not** cover all passing light: at an oblique
hit, one polarization still converts. Three things are open and become the requirement's tests:
1. whether *our* operator obeys the symmetry (computed by
   `directives/S11c_d_zinvariant_operator_blocks_directive.md`);
2. whether it still does with the drain flow on;
3. whether a real throat has the symmetry and supports a trapped twist mode.

---

## 0 · Why this exists

The user's directive, 2026-09-30, verbatim: *"It's a toy model. We start with the answer and work backwards. Is
there a clean condition in which we can show that light doesn't leak? If so, then we use that."*

Light is not observed to leak out of our three-dimensional space, and particles are observed to be stable.
Under requirements-first (`research/pde_ledger_v3/CHARTER.md:14–17`), that observation is treated as a
**requirement**. The program's job is to:

1. find a condition under which the model meets the requirement;
2. adopt that condition as a requirement with falsifiers;
3. test whether it clashes with any other requirement.

A clash is a falsification, so this is ⛔ not a rescue.

S11c-d asked the forward question on a **development input**: a strict rest-bulk, `LAB_HELD` question that "does
not yet answer leakage in the calibrated, draining medium" (flow-calibration assessment `:3`). Its bare input
ratio is `c_γ/c_s = 0.1`, and its selected modal ratio is `≈ 0.122` (`:11–13`). The cone lock `λγ = 1` is a
**calibrated, uncommitted** equality whose observational target is known to about 1 part in 10¹⁵
(`V3_STEP_PLAN.md:1107–1111`, `:1116–1126`).

## 1 · The condition — hypothesis H (a symmetry argument)

**H-planar.** Take a background that depends on a **single** in-plane direction (call it 1). The constitutive law
and face laws must transform covariantly, and every background **value** (tensor data, support, boundary data,
domain, response laws) must be **invariant** under the reflection of a second in-plane direction (3 → −3).
Covariance of a law is not enough: the chosen background values must be fixed by the symmetry. The S11c background
under the directive's restriction R1 is an example. For perturbations independent of direction 3:
- the in-plane displacement component `u_3` is **odd** under that reflection;
- every scalar (`θ`, `ζ_±`, `δp_s`, `μ_s`, `𝒜_s`, `J_s`, `V_s`, bulk `φ`) and the components `u_1`, `u_2` are
  **even**;
- every polar vector the face laws use (`n̂_s`, `v_face`, `v_bulk`, `t_s`) has an even part and an odd part, and
  under the reflection its odd part is its direction-3 component. On this background, the face normals have no
  direction-3 component.

A linear operator that commutes with the reflection cannot connect odd to even, so `u_3` evolves on its own at
linear order.

**Every incidence direction needs one more premise.** Suppose every background value is also invariant under the
full `O(2)` of rotations and reflections about direction 1. R1 states exactly this symmetry; it does not erase
otherwise allowed transverse components merely because their component labels contain 2 or 3. Then **each
tangential Fourier component** of an oblique
perturbation can be rotated into the form above, and a superposition is classified component by component. ⇒ At
such a planar interface, **for each incidence direction, the polarization perpendicular to the plane of
incidence** (TE-like) is decoupled. The polarization **in** the plane of incidence (TM-like) is not decoupled at
oblique incidence. Without full `O(2)` invariance, only incidence planes that contain an actual mirror of the
background are protected.

**H-round.** Take a background invariant under all rotations and reflections of the three in-plane coordinates
about a point (`O(3)`), *if* such a throat background exists. That means an achiral medium, a round throat, a parity-even
constitutive law, and a background flow with no azimuthal (swirl) component. At each angular order `(ℓ, m)` with
`ℓ ≥ 1`:
- **twist-type (toroidal)** displacements, tangent to the spheres `r = const` and divergence-free, have
  inversion parity `(−1)^{ℓ+1}`;
- every scalar field and the non-twist (spheroidal) displacements have parity `(−1)^ℓ`.

Rotation invariance conserves `(ℓ, m)`, and inversion conserves parity. So twist-type displacements evolve
independently at linear order.

**What reaches the bulk is a further step.** S11c-b "performs no curved-bulk response solve" (spec `:95–97`);
its supplied laws couple the slab to bulk trace operands through the face quantities (spec `:145–148`; S11c-a
`:343–354`, `:365–366`). The proposed closure criterion is **parity block-diagonality**, not the absence of every
twist-dependent face quantity. Odd-to-odd kinematics are allowed: for example, the direction-3 components of
`v_face,s` and `δ_v x_s` belong to the odd sector with `u_3`. What the symmetry forbids is an odd↔even block in
the constructed weak operator, its formal coordinate map, or the separately typed virtual/test kinematic map.
On class `P`, the supplied rest-frame potential-flow equations restrict the bulk traces, their normal jets,
density and current. An equivariant closure on that restricted domain cannot connect opposite parities, but Part
1 does not supply or test the missing curved-bulk closure. The v5 bulk-pullback equations are computed from the
supplied acoustic inputs in
`directives/_measurements/S11c_d_clean_condition_v5_author_scripts/bulk_pullback_pairing_audit.py`, with literal
stdout beside it.

**Premises the argument needs.** They do not all have the same status:

| # | premise | where it stands |
|---|---|---|
| P1 | The medium is achiral: the constitutive law has no parity-odd term. | **Supplied, not testable by Part 1.** The accepted basis is constructed with "in-plane `O(3)` isotropy and parity" (spec `:115–121`) and is recorded as the `O(3)`-Kronecker family (`:35–37`). A chiral extension is outside that supplied basis. |
| P2 | Every operand the face laws use transforms covariantly under the reflection (or `O(3)`), **and** every background value and support/boundary datum is invariant under it. This includes scalar, polar-vector and axial-vector operands, their normal jets, and every support field. No further field is present (e.g. a microrotation or director field, among the listed spin-carrier candidates at `native_light…:112–120`). | The face laws carry scalar and vector operands (S11c-a `:343–354`; S11c-b `:145–148`). Their transformation and the background-value invariance remain applicability conditions; Part 1 checks the represented operator on R1. |
| P3 | The background flow respects the symmetry: normal drain; radial in-plane flow at a throat; no swirl. | ⚠ Untested. `v_bulk_normal_0` "appears in no derived operator" (spec `:90–91`). The drain flow is **absent** from the S11c-b operator. |
| P4 | The truncations, the constraint fold (pin B), both anchorings, and the sign conventions in our operator do not break the reflection. | ⚠ Not yet computed in the operator. Round-1 leg scripts report pin B and both anchorings reflection-even on the planar class (stdout excerpts in the grounding file). Four upstream sign/coordinate repairs from S11c-d are unreviewed, so a convention error that breaks a reflection is exactly what Part 1 can catch. |

## 2 · What the condition covers

| case | covered by H? | consequence |
|---|---|---|
| light trapped at a throat (particle stability) | **linear order only, conditionally** — in an achiral medium, on an `O(3)`-symmetric, swirl-free throat background that admits a normalizable pure twist-type bound mode. No represented physical throat exists yet: the `h_±` graphs "are … not a complete nonlinear throat topology" (ontology `:315`). The support mode must be "a spectrally normalizable bound state or acceptably long-lived resonance of the complete variable-coefficient transverse operator" (ontology `:957`). | exact zero at linear order, if the conditions hold. ⛔ Not particle stability by itself: see the nonlinear gates of R-LEAK-1 |
| passing light, TE-like component at a planar one-direction interface; the twist part of passing light at a round scatterer | **yes**, subject to P1–P4 | no linear leak |
| passing light, TM-like at oblique incidence; the non-twist part at a round scatterer | **no** | converts → **R-LEAK-2** (a bound) |
| second order (nonlinear) | **no** — the square of a twist-type field is even and can source scalar deformation | half-two inventory |

**Transfer to the calibrated model.** H contains no speed coefficient. If the completed calibrated, draining
operator has the same symmetry and field content, H carries over to it unchanged. ⛔ That antecedent is not yet
shown: the drain is absent from every derived operator (P3), and the pilot "does not yet answer leakage in the
calibrated, draining medium" (assessment `:3`).

## 3 · Proposed requirements (to file in the register after review)

**R-LEAK-1 — trapped light (particle stability).** The trapped transverse brane-shear standing mode that "helps hold
each throat open" (ontology summary `:26`; the related statement at `:362` says it "helps hold the aperture open")
loses energy to the bulk at a rate below observational limits.
H supplies its adopted **linear** clean condition: the medium is achiral and the mode is pure twist-type about a
throat that is `O(3)`-symmetric in the brane coordinates and carries no swirl. ⛔ H alone does not deliver particle stability: a twist-type field's square
is even and can source scalar motion at second order. So R-LEAK-1 also needs a nonlinear zero or bound.

**Operator falsifier of the conditional linear selection rule:**
- **F1.** On a symmetric background with the drain flow carried live and the bulk closed, our linear operator
  contains a reflection-odd↔reflection-even block between the twist sector and the closed slab/face/bulk system.

**Applicability/failure tests for the proposed R-LEAK-1 realization:**
- **F2.** The model's spin carrier violates the adopted achiral/`O(3)` condition at linear order. Possibilities
  include a swirl, a chiral constitutive term, a non-round throat, or a net mixed `a`–`w` circulation: the last is
  an in-plane vector that selects a direction, so a nonzero value is not fixed by all rotations. On angular
  momentum: "A real linearly polarized standing wave
  can have zero time-averaged angular momentum. Two degenerate modes with a relative phase can form a circularly
  polarized bound pattern that carries angular momentum" (`native_light…:2377–2381`). Two degenerate twist-type
  modes in quadrature are such a pattern, and they leave the background symmetric at linear order.
- **F2b.** The spin carrier is an **additional field** (microrotation/director; `native_light…:113–116`). It
  enlarges the field content without breaking `O(3)`, so its parity and couplings must be classified before H
  can be applied.
- **F2c.** The required trapped mode is not pure twist-type. In particular, a trapped chiral-shear realization
  that mixes toroidal and poloidal sectors mixes opposite inversion parities and is not protected by H. The
  zero-helicity pure-twist representative and the nonzero-helicity mixed representative are computed in
  `directives/_measurements/S11c_d_clean_condition_v5_author_scripts/round_spin_carrier_audit.py`, with literal
  stdout beside it.
- **F5.** No normalizable twist-type bound mode exists on the throat background.
- **F6.** The full nonlinear throat, represented by the parent fields rather than `h_±` graphs (ontology `:315`),
  is not `O(3)`-symmetric.

**Nonlinear gates of R-LEAK-1** (half two; each must give a zero or a bound below the observational limit):
- **N1.** A twist-type standing mode must be able to hold a throat open.
- **N2.** Second-order leakage of the twist mode, including that from the throat deformation and mean flow it
  induces. The v5 round-sector script also computes nonzero `ℓ=0` and `ℓ=2` scalar projections for an explicit
  Coriolis image; that example identifies channels to retain, not a universal coefficient.

`native_light…:1343` already records that "the de-structured bulk carries no comparable shear channel … is not
sufficient by itself." H is the proposed **sufficient condition at linear order**. N1 and N2 remain.

**R-LEAK-2 — passing light.** For each kind of non-uniformity (throats; drain-flow gradients near masses), the
leak per encounter of the unprotected component stays below observational limits. H supplies **no generic
symmetry-protected zero** here; special incidences, coefficients or other symmetries might. It needs:
1. the applicable observational bounds, assembled. These are not done, and the bound that applies to gradual
   leakage, as opposed to a specific decay channel, has not been worked out;
2. an estimate of the per-encounter conversion at the calibrated point with the drain flow live.

## 4 · Options considered and not adopted

- **Bulk sound much faster than light.** This gives rapid kinematic **suppression** (exponential for suitable
  analytic profiles), not a zero: a localized profile still supplies some momentum transfer. It is also closed: `λγ = 1` is
  constrained by GW170817, given that `c_s` is the gravity-change speed (`V3_STEP_PLAN.md:1116–1126`). It reopens
  only if C13 (what a gravitational wave is) moves gravity signals off `c_s`.
- **A flow horizon**: bulk inflow toward the brane faster than `c_s`. By the standard acoustic-horizon argument
  it would block all outward propagation, but no operator in this program carries the drain, so it is
  unverified here. Not adopted:
  - the model's drain at throats carries material **out** of the brane into the bulk, with distributed return
    inward (ontology summary `:100`, `:1366`);
  - it does not protect a particle: no mechanism has been shown to return leaked mode energy coherently to the
    same particle.
- **Smallness only**: throat ≪ wavelength, weak fluid loading. Not clean; usable for R-LEAK-2 only.

## 5 · The computation

`directives/S11c_d_zinvariant_operator_blocks_directive.md` has two parts:

- **Part 1** uses the existing S11c-b slab/face operator on the planar one-direction background (restriction R1).
  - It constructs the complete weak operator payload and derives its coordinate domain mechanically from every
    free perturbation coordinate actually present. The manifest types every coordinate or records why it is held;
    no authored input list defines the map.
  - It emits the full formal coordinate map, its supplied potential-flow pullback (pressure, velocity, density and
    current, at the flat reference faces), and the virtual/test kinematic map. Every labelled entry is printed.
  - It runs pinned FORM controls K1 and K2 at stored energy and G3 as a fixed transverse background datum introduced
    after R1, so the face/trace/kinematic construction is also exercised. These names identify controls, not results.
  - It runs both engines, plus a raw-first comparator with separately labelled coefficient- and convention-mapped
    diagnostics. Completion under the guard is not presumed.
  - It inherits S11c-b's scope: **no drain flow and no bulk solve, both declared freezes**.
  - It also yields the TM-like blocks that R-LEAK-2 will need.

  What it can establish is the block structure of the rest-frame slab/face operator. R-LEAK-1 additionally needs
  the drain-on, round, closed-bulk case.
- **Part 2** is a **scope note** for that case: governing equations present/missing, especially how the drain
  enters; whether the three-direction jets can represent an `O(3)` profile; measured cost under the guard. Then
  it stops for the user's go/no-go. How the drain enters the operator is an **orchestrator spec question** (P3)
  before any build.

## 6 · What happens to S11c-d

- The `ω = 3` central benchmark (PID 4097233) was resumed on 2026-09-30 after its containment identity was checked
  (equal-speed feasibility `:47`). It remains a **labelled development benchmark**: rest bulk,
  `c_γ/c_s ≈ 0.12`, `LAB_HELD`, and supports no claim about the model. Its record forbids another job alongside it
  (assessment `:5`), so Part 1 waits for its completion. This packet makes no claim about the PID's present OS state.
- The passing-light question becomes R-LEAK-2, a bound check. It is no longer a matter of more precision at the
  development point.

## 7 · Prior art — an oracle, ⛔ not a premise (M3)

Existence and abstract were checked 2026-09-30 (measurements file §P). No full text was read, and no quoted
number or formula was verified.

| reference | relevance (from the abstract) |
|---|---|
| Molz & Beamish, JASA 99, 1894 (1996) | SH0 and L1 attenuation in an alumina membrane jumps when the surrounding helium freezes, attributed to radiation into the solid; SH0 generates only shear waves in the solid |
| Demma, Cawley & Lowe, JASA 113, 1880 (2003) | SH0 at thickness steps/notches, free plate |
| Kubrusly, von der Weid & Dixon, NDT&E Int. 108 (2019) | symmetric discontinuities create only modes sharing the incident mode's symmetry |
| Peyton et al., NDT&E Int. 158 (2026), doi 10.1016/j.ndteint.2025.103534 | SH0 → S0 Lamb conversion at a finite defect, free plate |
| Gu & Fuller, JASA 90, 2020 (1991) | sound radiation from subsonic wave scattering at discontinuities on fluid-loaded plates |
| "…fundamental torsional mode from axi-symmetric defects … in pipes", JASA 127(6), 3440 | title only |
| Friedland & Giannotti, PRL 100, 031602 (2008) | braneworld photon escape and stellar-cooling bounds; ⛔ its mapping to this model is **not** established (withdrawn as a falsifier, 2026-09-30) |

None of these has permeable faces, a drain flow, or a curl-only (MacCullagh) constitutive law. They may be used
only to check our computed result.
