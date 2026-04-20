# Moving-Throat PDE — Stage 72: Explicit Family-1 Verdict for the Minimal Isotropic Passive/Outgoing Quadrupole Branch

## Purpose

Stage 71 extracted the exact loading ratio selected by the natural contact-plus-pole realization of the minimal isotropic quadrupole precursor:

`rho_alpha = 4/3,`
`zeta_req = 1/3,`
`Pi_tr = (4/3) C_mix.`

The next honest step is to compare that branch directly against the explicit Family-1 support/source windows already frozen in Stages 62–70.

This stage shows that the match is not merely possible.
It is **strongly inside** the explicit Family-1 success region.

The main result is:

- the minimal isotropic passive/outgoing branch lies well below the explicit Family-1 loading-ratio ceiling,
- it lies in the exact symmetric-lowest-twin regime,
- and on the explicit Family-1 transport map it already succeeds at **zero** transport bias.

So under the natural contact-plus-pole identification, the explicit Family-1 support/source side is not the remaining reduced bottleneck anymore.

---

## 1. Exact comparison to the Family-1 loading-ratio window

From Stage 69, at the natural shell-weighted normalization `lambda_mu = 1` and in the unblocked limit,

`rho_suff^(chi) ≈ 3.46622291347846,`

`rho_fail^(chi) ≈ 3.46752913273870,`

`rho_max^(F1)   ≈ 3.46752922945601.`

From Stage 71, the minimal isotropic passive/outgoing quadrupole branch selects

`rho_alpha^(min) = 4/3 ≈ 1.33333333333333.`

So the exact margins are

`Delta_suff := rho_suff^(chi) - 4/3`
`           ≈ 2.13288958014513,`

`Delta_fail := rho_fail^(chi) - 4/3`
`           ≈ 2.13419579940537,`

`Delta_max  := rho_max^(F1) - 4/3`
`           ≈ 2.13419589612268.`

These are not small margins. The minimal isotropic branch sits far below even the guaranteed-success threshold.

---

## 2. Exact support-ratio and regime comparison

In the unblocked limit,

`zeta_req = rho_alpha - 1,`

so the same branch gives

`zeta_req^(min) = 1/3.`

The explicit Family-1 support ceiling from Stage 63/64 is

`zeta_max^(F1) ≈ 2.46752922945601.`

So the support-ratio margin is

`zeta_max^(F1) - zeta_req^(min)`
`≈ 2.13419589612268.`

And since

`0 < zeta_req^(min) = 1/3 < 1,`

the branch lies in the exact Stage-35 symmetric-lowest-twin window:

`C_mix < Pi_tr < 2 C_mix.`

So the minimal isotropic passive/outgoing branch does not need non-twin asymmetry at all.

---

## 3. Zero-transport-bias result on the explicit Family-1 branch

Stage 62 proved that the explicit Family-1 transport map obeys

`zeta_F1(Pe) = A_F1 Omega_Pe^2,`

with

`A_F1 ≈ 1.00005192880220,`

and that:

- if `zeta_req < A_F1`, the demand is already met at zero transport bias,
- if `A_F1 <= zeta_req <= zeta_max^(F1)`, a unique constructive `Pe_req` is needed.

For the minimal isotropic passive/outgoing branch,

`zeta_req^(min) = 1/3 < A_F1 ≈ 1.00005192880220.`

So the required transport bias is exactly

`Pe_req = 0`

on the explicit Family-1 branch.

Equivalently, the branch succeeds before any additional transport-driven overlap boost is turned on.

That is a very strong statement: the explicit Family-1 support/source branch already supports the minimal isotropic outgoing quadrupole demand in its zero-bias state.

---

## 4. What this changes

Before this step, the explicit Family-1 reduction had one unresolved outgoing-branch datum left: the loading ratio `rho_alpha`.

Stage 71 fixed that datum on the natural minimal isotropic passive/outgoing branch.
Stage 72 shows what that means numerically:

> **the explicit Family-1 support/source branch passes the minimal isotropic passive/outgoing quadrupole test with large margin.**

So for the explicit Family-1 branch, the remaining reduced theorem gap is no longer support sufficiency.
It is now the deeper question of whether the actual moving-throat grouped-`P2` / geometry branch really realizes the minimal isotropic contact-plus-pole conservative quadrupole module.

---

## 5. Best current theorem statement after Stage 72

Under the natural unblocked contact-plus-pole realization of the minimal isotropic passive/outgoing quadrupole branch,

`Y_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2),`

the explicit Family-1 support/source branch satisfies:

`rho_alpha = 4/3,`
`zeta_req = 1/3,`
`Pe_req = 0,`

and therefore lies safely inside the exact Family-1 success region.

So the explicit branch-level reduced theorem is now:

> **if the actual passive/outgoing moving-throat quadrupole branch realizes the minimal isotropic conservative precursor in the natural contact-plus-pole way, then the explicit Family-1 support/source side already succeeds without requiring transport bias or non-twin asymmetry.**
