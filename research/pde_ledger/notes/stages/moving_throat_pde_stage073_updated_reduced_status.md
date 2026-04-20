# Moving-Throat PDE — Stage 73: Updated Reduced Theorem Status After the Loading-Ratio Extraction

## Purpose

Stages 71–72 complete the explicit Family-1 support/source side under the natural minimal isotropic passive/outgoing quadrupole branch.

This note records the clean status update.

Script-backed status:
- `scripts/moving_throat_pde_stage073_updated_reduced_status_sympy_audit.py`
  rechecks the reduced-status verdict using the carried minimal-module and
  Family-1 threshold data without introducing new free constants.

The key new fact is that the explicit Family-1 branch no longer merely has a large admissible ceiling.
Under the natural contact-plus-pole realization of the minimal isotropic conservative quadrupole precursor, it selects

`rho_alpha = 4/3,`
`zeta_req = 1/3,`
`Pe_req = 0,`

and therefore succeeds with large margin.

So the remaining reduced theorem gap is no longer on the explicit support/source side.

---

## 1. What is now settled inside the reduced Family-1 branch

The following are now fixed:

1. the explicit Family-1 support/source branch has a hard constructive ceiling;
2. the outgoing normalization factors cancel from the explicit support theorem;
3. the support theorem depends only on `rho_alpha = alpha_req/alpha_mix`;
4. the natural contact-plus-pole interpretation of the minimal isotropic conservative quadrupole precursor gives

   `rho_alpha = 4/3`;

5. that value lies strictly inside the exact Family-1 success region;
6. and on the explicit Family-1 transport map it already succeeds at zero transport bias.

So the explicit Family-1 support/source branch is no longer the place where the reduced program can fail first.

---

## 2. What is still genuinely open

The remaining reduced question is now narrower and more structural:

> does the actual moving-throat grouped-`P2` / geometry branch really realize the minimal isotropic conservative quadrupole precursor in the natural contact-plus-pole way?

Equivalently, the remaining reduced task is not to strengthen the explicit support/source side again.
It is to derive the conservative quadrupole module itself from the real grouped-`P2` / geometry branch of the moving throat.

That is exactly the point where the reduced program now reconnects to the already-frozen 3PN and 2.5PN theorem ledgers.

---

## 3. Best current expert verdict after Stage 73

For the explicit Family-1 branch, the reduced moving-throat PDE program has advanced one step further than Stage 70:

- the explicit support/source side is finished,
- the outgoing branch no longer leaves a free loading-ratio datum under the natural minimal isotropic contact-plus-pole identification,
- and the remaining reduced bottleneck is now the derivation of that minimal isotropic conservative quadrupole module from the actual grouped-`P2` / geometry throat branch.

That is the right place for the next derivation phase.
