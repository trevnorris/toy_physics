# Light leakage — the clean condition (S11c-d re-scope proposal)

**Status:** PROPOSED · **v2** (2026-09-30) · orchestrator-written · ⛔ not governing until review-cleared.
- Round 1: Codex + Grok; v1 preserved at `7b38e9dc`; all twelve findings accepted
  (`directives/_measurements/S11c_d_clean_condition_review_disposition.md`).
- Repo citations below are reproduced verbatim, with their commands, in
  `directives/_measurements/S11c_d_clean_condition.md`.

**Plain summary.** We stop asking *"does light leak at a non-uniform slab?"* and ask instead *"under what clean
condition does it provably not leak?"* There is one candidate: **light whose motion is a pure twist about a
symmetric non-uniformity never moves anything that the bulk can feel, so it has no linear channel into the
bulk.** This is a standard symmetry argument, and it gives an exact zero rather than a small number. It
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
Under requirements-first (`CHARTER.md`), that observation is a **requirement**. The program's job is to:

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

**H-planar.** Take a background that depends on a **single** in-plane direction (call it 1), with the
constitutive law, face laws and all background data invariant under the reflection of a second in-plane direction
(3 → −3). The S11c background with every profile jet along directions 2 and 3 set to zero is an example. For
perturbations independent of direction 3:
- the in-plane displacement component `u_3` is **odd** under that reflection;
- every scalar (`θ`, `ζ_±`, `δp_s`, `μ_s`, `𝒜_s`, `J_s`, `V_s`, bulk `φ`) and the components `u_1`, `u_2` are
  **even**;
- every polar vector the face laws use (`n̂_s`, `v_face`, `v_bulk`, `t_s`) has an even part and an odd part, and
  under the reflection its odd part is its direction-3 component. On this background, the face normals have no
  direction-3 component.

A linear operator that commutes with the reflection cannot connect odd to even, so `u_3` evolves on its own at
linear order.

That background is also invariant under rotations about direction 1. So **each tangential Fourier component** of
an oblique perturbation can be rotated into this form, and a superposition is classified component by component.
⇒ At a planar one-direction interface, **for each incidence direction, the polarization perpendicular to the
plane of incidence** (TE-like) is decoupled. The polarization **in** the plane of incidence (TM-like) is not
decoupled at oblique incidence.

**H-round.** Take a background invariant under all rotations and reflections of the three in-plane coordinates
about a point (`O(3)`), *if* such a throat background exists. That means a round throat, a parity-even
constitutive law, and a background flow with no azimuthal (swirl) component. At each angular order `(ℓ, m)`:
- **twist-type (toroidal)** displacements, tangent to the spheres `r = const` and divergence-free, have
  inversion parity `(−1)^{ℓ+1}`;
- every scalar field and the non-twist (spheroidal) displacements have parity `(−1)^ℓ`.

Rotation invariance conserves `(ℓ, m)`, and inversion conserves parity. So twist-type displacements evolve
independently at linear order.

**What reached the bulk in each case is a further step.** In S11c-b the bulk enters the slab only through face
quantities (`δp_s`, `n̂_s·v_bulk,s`, `J_s`, `V_s`), and S11c-b "performs no curved-bulk response solve" (spec `:145–148`, `:95–97`). The
proposed interpretation, for review: suppose **no face-facing quantity depends on the twist-sector field, and the
twist-sector rows contain no bulk operand**. Then no bulk closure can connect that sector to the bulk, because the
bulk acts on the slab only through those quantities. Part 1 prints exactly those dependences.

**Premises the argument needs.** Each is a thing to check, ⛔ not to assume:

| # | premise | where it stands |
|---|---|---|
| P1 | The constitutive law has no parity-odd term. | The S11c-b energy basis is "the O(3)-Kronecker field-bilinear invariant family" (record `:35`). |
| P2 | Every operand the face laws use, and every support/boundary datum, transforms covariantly under the reflection (or `O(3)`). This covers scalars (`δp_s`, `μ_s`, `𝒜_s`, `J_s`, `V_s`), polar vectors (`n̂_s`, `v_face,s`, `v_bulk,s`, `t_s`) and any axial vector. No further field is present (e.g. a microrotation or director field, a listed spin-carrier candidate, `native_light…:113–116`). | The face laws carry vector operands (S11c-a `:350–354`; S11c-b `:145–148`). Both round-1 legs found these covariant on the planar background, as leg evidence about the term structure (`…review_scripts/`). |
| P3 | The background flow respects the symmetry: normal drain; radial in-plane flow at a throat; no swirl. | ⚠ Untested. `v_bulk_normal_0` "appears in no derived operator" (spec `:90–91`). The drain flow is **absent** from the S11c-b operator. |
| P4 | The truncations, the constraint fold (pin B), both anchorings, and the sign conventions in our operator do not break the reflection. | ⚠ Not yet computed in the operator. Round-1 legs found pin B and both anchorings reflection-even on the planar class. Four upstream sign/coordinate repairs from S11c-d are unreviewed, so a convention error that breaks a reflection is exactly what Part 1 can catch. |

## 2 · What the condition covers

| case | covered by H? | consequence |
|---|---|---|
| light trapped at a throat (particle stability) | **conditionally** — on an `O(3)`-symmetric, swirl-free throat background that admits a normalizable twist-type bound mode. No represented physical throat exists yet: the `h_±` graphs "are … not a complete nonlinear throat topology" (ontology `:315`). The support mode must be "a spectrally normalizable bound state or acceptably long-lived resonance of the complete variable-coefficient transverse operator" (ontology `:957`). | exact zero at linear order, if the conditions hold → **R-LEAK-1** |
| passing light, TE-like component at a planar one-direction interface; the twist part of passing light at a round scatterer | **yes**, subject to P1–P4 | no linear leak |
| passing light, TM-like at oblique incidence; the non-twist part at a round scatterer | **no** | converts → **R-LEAK-2** (a bound) |
| second order (nonlinear) | **no** — the square of a twist-type field is even and can source scalar deformation | half-two inventory |

**Transfer to the calibrated model.** H contains no speed coefficient. If the completed calibrated, draining
operator has the same symmetry and field content, H carries over to it unchanged. ⛔ That antecedent is not yet
shown: the drain is absent from every derived operator (P3), and the pilot "does not yet answer leakage in the
calibrated, draining medium" (assessment `:3`).

## 3 · Proposed requirements (to file in the register after review)

**R-LEAK-1 — trapped light.** The trapped transverse brane-shear standing mode that "helps hold each throat open"
(ontology summary `:26`, `:362`) is twist-type about a throat that is `O(3)`-symmetric in the brane coordinates
and carries no swirl.

Falsifiers:
- **F1.** On a symmetric background with the drain flow carried live and the bulk closed, our linear operator
  couples twist-type displacements to any bulk-facing quantity.
- **F2.** The model's spin carrier forces the background to break `O(3)` at linear order. Possibilities: a swirl,
  a chiral constitutive term, or a non-round throat. Note that only a **circulating** (complex, travelling)
  twist mode carries angular momentum. A real standing twist mode carries none (`native_light…:2377`, "A trapped
  standing wave is not automatically spinning"). A circulating twist mode leaves the background symmetric at
  linear order.
- **F2b.** The spin carrier is an **additional field** (microrotation/director; `native_light…:113–116`). It
  enlarges the field content without breaking `O(3)`, so its parity and couplings must be classified before H
  can be applied.
- **F3.** A twist-type standing mode cannot hold a throat open. This is nonlinear; half two.
- **F4.** Second-order leakage of the twist mode, including that from the throat deformation it induces. This is
  nonlinear; half two.
- **F5.** No normalizable twist-type bound mode exists on the throat background.
- **F6.** The full nonlinear throat, represented by the parent fields rather than `h_±` graphs (ontology `:315`),
  is not `O(3)`-symmetric.

`native_light…:1343` already records that "the de-structured bulk carries no comparable shear channel … is not
sufficient by itself." R-LEAK-1 is the proposed **sufficient condition at linear order**.

**R-LEAK-2 — passing light.** For each kind of non-uniformity (throats; drain-flow gradients near masses), the
leak per encounter of the unprotected component stays below observational limits. H supplies **no generic
symmetry-protected zero** here; special incidences, coefficients or other symmetries might. It needs:
1. the applicable observational bounds, assembled. These are not done, and the bound that applies to gradual
   leakage, as opposed to a specific decay channel, has not been worked out;
2. an estimate of the per-encounter conversion at the calibrated point with the drain flow live.

## 4 · Options considered and not adopted

- **Bulk sound much faster than light.** This gives kinematic **suppression**, exponential for smooth profiles,
  not a zero: a localized profile still supplies some momentum transfer. It is also closed: `λγ = 1` is
  constrained by GW170817, given that `c_s` is the gravity-change speed (`V3_STEP_PLAN.md:1116–1126`). It reopens
  only if C13 (what a gravitational wave is) moves gravity signals off `c_s`.
- **A flow horizon**: bulk inflow toward the brane faster than `c_s`. By the standard acoustic-horizon argument
  it would block all outward propagation, but no operator in this program carries the drain, so it is
  unverified here. Not adopted:
  - the model's drain at throats carries material **out** of the brane into the bulk, with distributed return
    inward (ontology summary `:100`, `:1366`);
  - it would also trap gravity changes;
  - it does not protect a particle, because leaked energy returns to the brane, not to the particle.
- **Smallness only**: throat ≪ wavelength, weak fluid loading. Not clean; usable for R-LEAK-2 only.

## 5 · The computation

`directives/S11c_d_zinvariant_operator_blocks_directive.md` has two parts:

- **Part 1** uses the existing S11c-b slab/face operator on the planar one-direction background (restriction R1).
  - It builds the complete linear operator on perturbations independent of direction 3, with every bulk-facing
    quantity printed as an explicit row/column.
  - It runs pinned FORM controls (K1: a fixed-axis term; K2: a single-Levi-Civita term).
  - It runs both engines, plus a blockwise comparator. All jobs run under the guarded runner.
  - It inherits S11c-b's scope: **no drain flow and no bulk solve, both declared freezes**.
  - It also yields the TM-like blocks that R-LEAK-2 will need.

  What it can establish is the block structure of the rest-frame slab/face operator. R-LEAK-1 additionally needs
  the drain-on, round, closed-bulk case.
- **Part 2** is a **scope note** for that case: governing equations present/missing, especially how the drain
  enters; whether the three-direction jets can represent an `O(3)` profile; measured cost under the guard. Then
  it stops for the user's go/no-go. How the drain enters the operator is an **orchestrator spec question** (P3)
  before any build.

## 6 · What happens to S11c-d

- The suspended `ω = 3` central benchmark (PID 4097233) stays a **labelled development benchmark**: rest bulk,
  `c_γ/c_s ≈ 0.12`, `LAB_HELD`. It supports no claim about the model. Whether to finish it is the user's call.
  ⚠ Its own containment says not to launch another job alongside it (assessment `:5`). Part 1's engines wait on
  that decision.
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
