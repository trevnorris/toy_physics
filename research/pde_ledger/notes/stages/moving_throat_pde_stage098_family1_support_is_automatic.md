# Moving-Throat PDE — Stage 098: The Explicit Family-1 Support Test Is Automatic on the Actual Isotropic Branch

## Purpose

Stage 80 reduced the actual isotropic passive/outgoing quadrupole branch to one scalar normalization defect `N_Q`.

But the explicit branch still has a support/source theorem sitting beside that radiative theorem.
The next honest question is therefore:

> does the explicit Family-1 support/source side still need to be checked separately once the actual branch loading ratio is fixed?

This stage shows that the answer is effectively **no**.

On the actual isotropic branch,

`rho_alpha = 4/3`,

so even with finite blocking the exact support demand is so small that any explicit branch with support ceiling `zeta_max > 1` already passes it throughout the admissible blocked regime.

The Family-1 branch has

`zeta_max^(F1) ≈ 2.46752922945601 > 1`,

so its support/source side is automatic.

---

## 1. Exact blocked demand on the actual isotropic branch

Stage 68 gave the exact support-ratio demand

`zeta_req = (rho_alpha - 1) / (1 - eps_blk (2 - rho_alpha)).`

Inserting the actual isotropic branch value

`rho_alpha = 4/3`

gives

`zeta_req^(act)(eps_blk) = (1/3)/(1 - (2/3) eps_blk) = 1/(3 - 2 eps_blk).`

So the entire blocked support demand is one monotone increasing function of the blocking parameter.

---

## 2. Automatic support theorem for any branch with `zeta_max > 1`

Assume a branch has support ceiling `zeta_max > 1`, and use the same admissible blocking window already frozen earlier,

`0 <= eps_blk < 1/zeta_max`.

Because `zeta_req^(act)` increases with `eps_blk`, its worst value on the admissible window is bounded by

`zeta_req^(act) < 1 / (3 - 2/zeta_max).`

Now compare that to `zeta_max` itself. The inequality

`1 / (3 - 2/zeta_max) < zeta_max`

is equivalent to

`3 zeta_max (zeta_max - 1)/(3 zeta_max - 2) > 0,`

which is true for every `zeta_max > 1`.

Therefore:

> for the actual isotropic branch with `rho_alpha = 4/3`, any explicit support/source family whose constructive ceiling satisfies `zeta_max > 1` already passes the support test throughout the full admissible blocked regime.

So the support/source side is automatic once the actual isotropic branch is fixed.

---

## 3. Family-1 specialization

For the explicit Family-1 branch,

`zeta_max^(F1) ≈ 2.46752922945601 > 1`.

Therefore the automatic theorem above applies immediately.

At zero blocking,

`zeta_req^(act)(0) = 1/3,`

and even at the edge of the admissible Family-1 blocked window,

`eps_blk -> 1 / zeta_max^(F1)`,

the demand remains bounded by

`zeta_req^(act) < zeta_max^(F1)`.

Numerically,

`zeta_req^(act) < 0.456730991107963 < 2.46752922945601`.

So the explicit Family-1 support/source branch no longer carries an independent reduced theorem burden.
It is already safe once the actual isotropic outgoing branch is accepted.

---

## 4. What remains after this step

After Stages 80–81, the reduced theorem split is now fully sharp:

- **support/source theorem:** automatic on the actual isotropic branch,
- **radiative theorem:** controlled by the single normalization defect `N_Q`.

So the explicit Family-1 branch has now dropped out of the active uncertainty ledger.
The only remaining reduced theorem question is whether the completed moving-throat PDE gives

`N_Q = 1`.
