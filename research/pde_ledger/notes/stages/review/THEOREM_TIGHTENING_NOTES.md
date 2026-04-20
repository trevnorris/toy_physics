# Moving-Throat PDE — Theorem Tightening Notes

## Purpose

This file turns the first hardening-tranche review into an edit list for the
notes themselves.

The issue here is no longer symbolic correctness. The issue is statement scope:
some results are exact only

- on a selected branch,
- inside an explicit closure family,
- to first order,
- on a numerical bracket,
- or within the reduced/coherent model rather than the full PDE.

Those distinctions should stay visible in the note prose.

## High-Priority Tightening Targets

## 1. Reduced-sector and outgoing foundations

### Stage `003`

- Keep the perturbative pole-shift claim explicitly away from resonance.
- Prefer wording like:
  `on the weak-coupling separated-pole branch varpi^2 > Omega_eta^2`.
- Avoid wording that sounds globally exact for the resonant case.

### Stage `004`

- Keep the outgoing dressing explicitly first-order in `Pi_out`.
- Keep the scalar rescue statement tied to the derivative-coupling assumption.
- Avoid wording that sounds like a full multimode Maxwell completion theorem.

## 2. Continuum/coherent support chain

### Stages `020-027`

- Keep “exact” attached to the reduced continuum operator and selected-branch
  formulas, not to the completed PDE.
- Make sure “physical branch” means:
  the root selected inside the reduced continuum model by the stated continuity
  condition.
- Avoid wording that suggests the PDE has already been proved to realize that
  reduced branch.

### Stage `028`

- Keep the main theorem explicitly tied to the coherent local-kernel hypothesis.
- Prefer:
  `for the first coherent local D/N kernel, the reduced branch lands exactly on
  tracking`.
- Avoid:
  `the moving-throat PDE lands exactly on tracking`, unless the PDE hypothesis
  is restated right there.

### Stage `029`

- Keep the ordering theorem on the constructive split-`U` tracking branch.
- Avoid global comparisons that suppress the domain
  `0 < R < 1`, `0 < xi < 1`, `delta > 0`.

### Stage `030`

- Keep `R_target` independence of `zeta` tied to the coherent tracking branch.
- Avoid wording that makes this independence sound generic for the unrestricted
  rank-2 closure.

### Stage `031`

- Keep the compensation theorem as a reduced-level existence statement.
- Preferred wording:
  `there is no reduced-level support no-go`.
- Keep visible that the unresolved PDE question is whether the realized physical
  `zeta` actually reaches `zeta_req`.

## 3. Coupled source/support fixed point

### Stage `041`

- Phrase the main theorem as existence plus interval bounds, with a small-`Xi`
  asymptotic.
- Avoid wording that implies a closed-form or global uniqueness theorem for
  `Pe_*`.

## 4. Canonical outgoing closure

### Stage `089`

- Keep the GR-target statement conditional on the canonical outgoing DtN branch
  and the strict point-particle limit.
- Avoid prose that sounds like unconditional recovery of the GR coefficients by
  the full PDE.

## 5. Positive-source and explicit Family-1 mouth closure

### Stage `108`

- Keep the positive-source theorem tied to a positive normalized localized mouth
  source on the first D/N interval.
- Avoid reading the exclusion of `g_+^F1` as a theorem for arbitrary sign-changing
  or nonlocalized mouth data.

### Stages `125-127`

- Keep “unique regular branch” explicitly inside the explicit Family-1
  positive-mouth closure.
- Prefer:
  `inside the explicit exponential positive-mouth family`.
- Avoid shorthand that sounds like the full PDE has a unique mouth branch.

### Stage `128`

- Keep the stage marked as consolidation/status, not as a new theorem.

## 6. Mouth deformations and rigidity

### Stages `129-130`

- Keep all deformation formulas marked first-order or perturbative.
- Avoid the word `exact` for `delta Pi`, `delta S_q`, or the rigidity kernel
  unless explicitly qualified as `exact at first order`.

### Stage `131`

- Keep the two-family evaluation as evidence/example, not as exhaustive
  classification of positive mouth deformations.

### Stage `132`

- Keep the stage marked as status/consolidation.
- Avoid implying that mouth rigidity has been proved for all positive profiles.

### Stages `133-135`

- Keep Stage `133` exact only for the explicit canonical branch data.
- Keep Stages `134-135` marked as first-order and mouth-only.
- Keep the Stage-`135` Picard update labeled as a consistency check or one-step
  nonlinear probe, not a convergence theorem.

### Stage `136`

- Keep the stage marked as status.

## 7. Co-evolving core-mouth fixed point

### Stage `137`

- Keep the co-evolving map exact, but distinguish that the defect transport law
  is first-order near the canonical branch.
- Avoid wording that sounds like the compensated branch has already been solved
  again at this stage.

### Stage `138`

- Keep the fixed-point statement tied to the frozen-traction slice
  `Sigma0 = Sigma0^*`.
- If uniqueness is stated, qualify it as:
  `the solved positive fixed point on the analyzed branch/window is unique`.
- Avoid unqualified global uniqueness claims unless separately proved.

### Stage `139`

- Keep the condition `g_fp(Sigma0) = g_*` exact, but keep the realized branch
  values explicitly numerical.
- Preferred wording:
  `the exact renormalized compensation condition has a numerically located
  branch at ...`.
- Avoid language implying a closed-form formula for the renormalized point.

### Stage `140`

- Keep the stage marked as status/consolidation.
- Make sure quoted summary numerics match `138-139` exactly.

## 8. Linear grouped anisotropy and microscopic orbit structure

### Stage `153`

- Keep the result as a linear grouped outlet map on the compensated isotropic
  branch.
- Avoid wording that sounds like the actual physical anisotropies are already
  solved.

### Stage `168`

- Keep the monomial statement explicitly linearized/reference-branch.
- Preferred wording:
  `zero defect is equivalent, at first grouped weak-axisymmetric order, to
  preservation of the three direct microscopic monomials`.
- Avoid wording that sounds like finite nonlinear monomial preservation has
  already been proved.

### Stage `169`

- Keep the theorem explicitly tangent-space / infinitesimal.
- Avoid prose that makes the finite orbit statement sound already complete at
  Stage `169`; that belongs to Stage `170`.

### Stage `170`

- Keep the quotient theorem explicitly inside the positive coherent microscopic
  sector.
- Keep visible that the remaining open question is dynamical branch selection:
  whether the actual PDE branch preserves those invariants and stays on a single
  similarity orbit.
- Avoid treating Stage `170` as the final PDE branch theorem. It is the final
  invariant-structure theorem currently proved.

## Phrase Replacements

- Replace `exact` with `exact on the coherent branch` when the branch matters.
- Replace `unique branch` with `unique regular branch inside the explicit
  Family-1 closure` when appropriate.
- Replace `the PDE does` with `the reduced coherent model does` when the
  completion theorem is still open.
- Replace `solves exactly` with `is characterized exactly and then solved
  numerically` when the final values come from numerical root-finding.
- Replace `proof complete` with `reduced/invariant structure complete` when the
  dynamical branch theorem is still missing.

## Suggested Edit Order

1. Tighten branch scope in `028-031`, `089`, `108`, `125-127`, `168-170`.
2. Tighten perturbative wording in `129-135`.
3. Tighten numerical fixed-point wording in `138-140`.
4. Then do a final read-through of the status stages `128`, `132`, `136`, `140`
   so the summaries match the tightened theorem scope exactly.
