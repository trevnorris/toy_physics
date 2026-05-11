# Moving-Throat PDE — Stage 096: Actual Branch Verdict for the Geometry-Lane Check

## Purpose

Stages 77–78 were the direct answer to the open question left at Stage 76:

> on the actual moving-throat branch, does the geometry lane stay dynamically inert through `O(omega^4)`?

This note records the verdict.

Script-backed status:
- `scripts/moving_throat_pde_stage147_geometry_lane_check_verdict_sympy_audit.py`
  rechecks the isotropic `l=0 <-> l=2` decoupling verdict and collapses the
  obstruction formula to the clean `3/4 + 1/4` conservative module.

---

## 1. Actual reduced-branch answer

Inside the present reduced hierarchy, the actual natural branch is the isotropic branch selected by:

- the isotropic quadratic wall operator,
- the grouped real `P2` quadrupole carrier,
- and the 3PN result that the leftover scalar contribution is a **static** geometry completion rather than a second dynamic quadrupole pole.

On that branch, Stage 77 proves exactly that the scalar/geometry `l=0` lane and the grouped real `l=2` bundle are block diagonal in the quadratic wall theory.
So the geometry lane cannot contribute dynamic `omega^2` or `omega^4` moments to the isotropic quadrupole module.

Therefore

`eps_2 = eps_4 = 0`

on the actual isotropic branch in the present hierarchy.

---

## 2. Consequence for the conservative quadrupole module

With the contamination numbers zero, the Stage-75 obstruction formula collapses to

`c_pole = 1/4`,

`c_geom = 3/4`.

So the actual branch realizes

`Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

And therefore the already-carried support/source conclusions follow directly:

`rho_alpha = 4/3`,

`zeta_req = 1/3`.

So the geometry-lane check comes out **clean** on the actual isotropic branch.

---

## 3. What is still open after the check

This does **not** solve the full moving-throat PDE.
It removes one remaining reduced ambiguity.

What remains open is no longer whether the geometry lane contaminates the minimal isotropic quadrupole module. It does not, on the current actual branch.

The remaining serious open question is deeper and more physical:

> does the completed moving-throat PDE really realize the natural isotropic grouped-`P2` one-pole branch itself, with the passive/outgoing normalization required by the 2.5PN bridge?

So the geometry-lane check is now finished at reduced level.
The live gap is again the same narrow one isolated by the 2.5PN program: the final passive/outgoing quadrupole normalization on the true moving-throat branch.
