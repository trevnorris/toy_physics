# Light leakage — the clean condition (S11c-d re-scope proposal)

**Status:** PROPOSED · orchestrator-written 2026-09-30 · under two-leg review (Codex + Grok) · ⛔ not governing
until review-cleared. Every repo citation below is reproduced verbatim, with its command, in
`directives/_measurements/S11c_d_clean_condition_lookups.md`.

**Plain summary.** We stop asking *"does light leak at a non-uniform slab?"* and ask instead *"under what clean
condition does it provably not leak?"* There is one: **light whose motion is a pure twist about a symmetric
non-uniformity never moves the slab's faces in or out, so it has no way to push the bulk.** This is a standard
symmetry argument. It is exact rather than merely small. It does not depend on the light/bulk speed ratio or on
the drain flow's strength. It covers **light trapped in a round throat**, which is the particle-stability case
and the one that needs an exact zero. It does **not** cover all passing light: at an oblique hit, one
polarization still converts. Passing light therefore becomes a bound check. Whether *our* operator actually
obeys the symmetry is not assumed. It is computed by the companion directive
`directives/S11c_d_zinvariant_operator_blocks_directive.md`.

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

S11c-d asked the forward question at a **development model point**, which does not transfer to the model:
- bare input ratio `c_γ/c_s = 0.1`; selected modal ratio `≈ 0.122` (flow-calibration assessment `:11–13`);
- strict `v_bulk_normal_0 = 0` (amendment `:375–379`);
- `LAB_HELD` anchoring.

The calibrated model sets `λγ = 1`, and observation constrains that to about 1 part in 10¹⁵
(`V3_STEP_PLAN.md:1116–1126`).

## 1 · The condition — hypothesis H (a symmetry argument, ⛔ not yet computed in our operator)

**H-planar.** Take a background invariant under translation along an in-plane direction `z` and under the
reflection `z → −z`. The S11c interface background `W_bg(y)`, `μ_R,bg(y)` is an example. For perturbations with
no `z`-dependence:
- the `z`-polarized in-plane displacement `u_z` is **odd** under the reflection;
- every scalar field (`θ`, `ζ_±`, the face-response variables, bulk `φ`/`δp`) and the other displacement
  components `u_x`, `u_y` are **even**.

A linear operator that commutes with the reflection cannot connect odd to even. So `u_z` evolves on its own
at linear order.

The interface background is also invariant under rotations about its normal (`y`). Any oblique perturbation can
therefore be rotated into this form. ⇒ At a planar interface, **for every incidence direction, the polarization
perpendicular to the plane of incidence** (TE-like) is decoupled. The polarization **in** the plane of incidence
(TM-like) is not decoupled at oblique incidence: it has a component along `∇W_bg`.

**H-round.** Take a background invariant under all rotations and reflections of the three in-plane coordinates
about a point (`O(3)`). That means a round throat, a constitutive law with no parity-odd term, and a background
flow with no azimuthal (swirl) component. At each angular order `(ℓ, m)`:
- **twist-type (toroidal)** displacements, tangent to the spheres `r = const` and divergence-free, have
  inversion parity `(−1)^{ℓ+1}`;
- every scalar field and the non-twist (spheroidal) displacements have parity `(−1)^ℓ`.

They therefore evolve independently at linear order. Twist-type trapped light has **no linear channel** into
anything that pushes the bulk.

**Premises the argument needs.** Each is a thing to check, ⛔ not to assume:

| # | premise | where it stands |
|---|---|---|
| P1 | The constitutive law has no parity-odd term. | The S11c-b energy basis is "the O(3)-Kronecker field-bilinear invariant family" (record `:33–35`). It is built from Kronecker contractions. |
| P2 | Every brane–bulk coupling enters through scalar face/bulk quantities (`ζ_±`, face flux/permeability responses, `θ`, `φ`/`δp`). | S11c-b field content: `u` has "three in-plane components, no w-component" (spec `:59–66`). |
| P3 | The background flow respects the symmetry: normal drain; radial in-plane flow at a throat; no swirl. | ⚠ Untested. `v_bulk_normal_0` "appears in no derived operator" (spec `:90–91`). The drain flow is **absent** from the S11c-b operator. |
| P4 | The truncations, the constraint fold (pin B), the anchorings and the sign conventions in our operator do not break the reflection. | ⚠ Untested. S11c-d exposed four upstream sign/coordinate repairs that are still unreviewed. A convention error that breaks a reflection is exactly what the computation can catch. |

## 2 · What the condition covers

| case | covered by H? | consequence |
|---|---|---|
| light trapped at a throat (particle stability) | **yes**, if the trapped mode is twist-type and the throat is `O(3)`-symmetric and swirl-free | exact zero at linear order → **R-LEAK-1** |
| passing light, TE-like polarization at a planar interface; the twist part of passing light at a round scatterer | **yes** | no linear leak |
| passing light, TM-like polarization at oblique incidence; the non-twist part at a round scatterer | **no** | converts → **R-LEAK-2** (a bound, not a zero) |
| second order (nonlinear) | **no** | half-two inventory |

H holds for any speed ratio, any profile shape and any frequency. For the round case it also holds for any
symmetric drain-flow magnitude, though P3 is still untested. ⇒ unlike the S11c-d benchmark, a confirmed H
**transfers to the calibrated model** (`λγ = 1`, drain flow on), wherever the symmetry holds.

## 3 · Proposed requirements (to file in the register after review)

**R-LEAK-1 — trapped light.** The trapped transverse brane-shear standing mode that "helps hold each throat
open" (ontology summary `:26`, `:362`) is twist-type about a throat that is `O(3)`-symmetric in the brane
coordinates and carries no swirl.

Falsifiers:
- **F1.** On a symmetric background with the drain flow carried live, our linear operator couples twist-type
  displacements to any field that pushes the bulk (the computation, §5).
- **F2.** The model's spin carrier forces the background to break `O(3)` at linear order. Possibilities: a swirl,
  a chiral constitutive term, or a non-round throat. The spin carrier is open (`native_light…:451`; its
  candidates `:113–116` include "trapped chiral shear"). A circulating twist mode (`m ≠ 0`) carries angular
  momentum while leaving the background symmetric at linear order. A symmetry-breaking carrier makes the leak
  nonzero, and it must then be bounded.
- **F3.** A twist-type standing mode cannot hold a throat open. This is nonlinear; half two.
- **F4.** Second-order leakage of the twist mode, including that from the throat deformation it induces. This is
  nonlinear; half two.

`native_light…:1343` already records that "the de-structured bulk carries no comparable shear channel … is not
sufficient by itself." R-LEAK-1 is the proposed **sufficient condition at linear order**.

**R-LEAK-2 — passing light.** For each kind of non-uniformity (throats; drain-flow gradients near masses), the
leak per encounter of the unprotected polarization stays below observational limits. It needs:
1. the applicable observational bounds, assembled. These are not done, and the bound that applies to gradual
   leakage, as opposed to a specific decay channel, has not been worked out;
2. an estimate of the per-encounter conversion at `λγ = 1` with the drain flow live.

No clean zero exists for this case.

## 4 · Options considered and not adopted

- **Bulk sound much faster than light**, which would make smooth non-uniformities kinematically safe. Closed:
  `λγ = 1` is constrained by GW170817, given that `c_s` is the gravity-change speed (`V3_STEP_PLAN.md:1116–1126`).
  It reopens only if C13 (what a gravitational wave is) moves gravity signals off `c_s`.
- **A flow horizon**: bulk inflow toward the brane faster than `c_s`. Exact for all light, but not adopted:
  - the model's drain at throats carries material **out** of the brane into the bulk, with distributed return
    inward (ontology summary `:100`, `:1366`);
  - it would also trap gravity changes;
  - it does not protect a particle, because leaked energy returns to the brane, not to the particle.
- **Smallness only**: throat ≪ wavelength, weak fluid loading. Not clean; usable for R-LEAK-2 only.

## 5 · The computation

`directives/S11c_d_zinvariant_operator_blocks_directive.md` has two parts:

- **Part 1** uses the existing S11c-b operator at the planar interface. It builds the complete linear operator on
  `z`-independent perturbations and prints every block, with FORM controls that add a symmetry-breaking term to
  the action. It inherits S11c-b's scope: **the drain flow is absent, a declared freeze**. It also yields the
  TM-like blocks that R-LEAK-2 will need.
- **Part 2** covers the round background with the drain flow live, but only as a **scope note** (governing
  equations present/missing, measured cost). Then it stops for the user's go/no-go. How the drain flow enters the
  slab/face/bulk equations is currently unspecified (P3). That is a spec question for the orchestrator before
  any build.

## 6 · What happens to S11c-d

- The suspended `ω = 3` central benchmark (PID 4097233) stays a **labelled development benchmark**: rest bulk,
  `c_γ/c_s ≈ 0.12`, `LAB_HELD`. It supports no claim about the model. Whether to finish it is the user's call.
  ⚠ Its own containment says not to launch another job alongside it (assessment `:5`). Part 1 waits on that
  decision.
- The passing-light question becomes R-LEAK-2, a bound check. It is no longer a matter of more precision at the
  development point.

## 7 · Prior art — an oracle, ⛔ not a premise (M3)

Existence and abstract were checked 2026-09-30 (lookups §P). No full text was read, and no quoted number or
formula was verified.

| reference | relevance |
|---|---|
| Molz & Beamish, JASA 99, 1894 (1996) | a uniform SH0 plate mode does not radiate into fluid helium; it radiates once the helium is solid (acquires shear) |
| Demma, Cawley & Lowe, JASA 113, 1880 (2003) | SH0 at thickness steps/notches, free plate |
| Kubrusly, von der Weid & Dixon, NDT&E Int. 108 (2019) | symmetric discontinuities convert only to modes of the same symmetry |
| Peyton et al., NDT&E Int. 158 (2026), doi 10.1016/j.ndteint.2025.103534 | SH0 → S0 Lamb conversion at a finite defect, free plate |
| Gu & Fuller, JASA 90, 2020 (1991) | a subsonic plate wave radiates once it scatters off a discontinuity |
| "…fundamental torsional mode from axi-symmetric defects … in pipes", JASA 127(6), 3440 | title only: torsional waves at axisymmetric defects |
| Friedland & Giannotti, PRL 100, 031602 (2008) | braneworld photon escape and stellar-cooling bounds; ⛔ its mapping to this model is **not** established (withdrawn as a falsifier, 2026-09-30) |

None of these has permeable faces, a drain flow, or a curl-only (MacCullagh) constitutive law. They may be used
only to check our computed result.
