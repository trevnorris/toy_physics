# Moving-Throat PDE — Stage 085: Exact Cancellation of the Outgoing-Normalization Factors in the Selected Quadrupole-Demand Product

## Purpose

Stage 084 left the reduced moving-throat PDE in the form of one master residual

`R_quad = zeta_req(Pi_tr,C_mix,eps_blk) - zeta_phys(Pe_*),`

but the selected branch still looked as if it depended on several separate outgoing-normalization objects at once:

- the target quadrupole normalization `N_Q^(target)`,
- the selected static prefactor `P_{0,-}`,
- the transfer factor `beta_0`,
- the selected overlap `s_-`,
- and the selected conservative eigenvalue `lambda_-`.

The next honest step is to ask whether those are really independent in the final support test.

This stage shows that they are not.

The main result is the exact pair of identities

`Pi_tr = (N_Q^(target)/beta_0) alpha_req,`

`C_mix = (N_Q^(target)/beta_0) alpha_mix,`

and therefore

`Pi_tr / C_mix = alpha_req / alpha_mix.`

Using the selected-mode normalization relation from Stage 030,

`N_Q^(target) = mhat_-^2 beta_0 s_- / lambda_-,`

the same identities become

`Pi_tr = mhat_-^2 (s_- / lambda_-) alpha_req,`

`C_mix = mhat_-^2 (s_- / lambda_-) alpha_mix.`

So once the outgoing quadrupole branch is required to hit the `2G/(5 c^5)` normalization target, the explicit support theorem no longer depends separately on `(mhat_-, beta_0, s_-, lambda_-)`.
It depends only on the **loading ratio** `alpha_req / alpha_mix`.

That is the cleanest bridge yet between the Stage-030 selected quadrupole normalization stack and the explicit Family-1 support/source branch.

---

## 1. Exact product identities

From Stages 035–036, the selected-branch target ratio is

`R_target = N_Q^(target) A / ( beta_0 kappa_0^2 ),`

with the exact D/N constant

`kappa_0^2 = 8 / pi^2.`

The exact dimensionless total and mixed loadings are

`G_tr = 8 alpha_req / (pi^2 A),`

`M_mix = 8 alpha_mix / (pi^2 A).`

By construction,

`Pi_tr = R_target G_tr,`

`C_mix = R_target M_mix.`

Inserting the formulas above and using `kappa_0^2 = 8/pi^2` gives the exact cancellations

`Pi_tr = (N_Q^(target)/beta_0) alpha_req,`

`C_mix = (N_Q^(target)/beta_0) alpha_mix.`

So the selected demand ratio is simply

`Pi_tr / C_mix = alpha_req / alpha_mix.`

This identity is exact.

---

## 2. Spectral form using the selected-mode normalization stack

Stage 030 gave the selected-mode outgoing normalization relation

`mhat_-^2 P_{0,-} = N_Q^(target),`

with

`P_{0,-} = beta_0 s_- / lambda_-`.

Therefore

`N_Q^(target) / beta_0 = mhat_-^2 s_- / lambda_-`.

Substituting this into the exact product identities yields the spectral forms

`Pi_tr = mhat_-^2 (s_- / lambda_-) alpha_req,`

`C_mix = mhat_-^2 (s_- / lambda_-) alpha_mix.`

So the same common spectral factor multiplies both the demanded product and the mixed baseline.

That is why the explicit support test loses all separate dependence on the outgoing normalization amplitudes once the branch is normalized.

---

## 3. Exact selected-demand law in pure loading-ratio form

Stage 052 gave the exact selected support demand

`zeta_req = (Pi_tr - C_mix) / [ C_mix - eps_blk (2 C_mix - Pi_tr) ].`

Using the cancellation identities above, the entire right-hand side reduces to

`zeta_req`
`= (alpha_req - alpha_mix)`
`  / [ alpha_mix - eps_blk (2 alpha_mix - alpha_req) ].`

Equivalently,

`zeta_req`
`= [ alpha_req/alpha_mix - 1 ]`
`  / [ 1 - eps_blk ( 2 - alpha_req/alpha_mix ) ].`

So the explicit support demand is now a function only of the loading ratio

`rho_alpha := alpha_req / alpha_mix`.

In the unblocked limit `eps_blk = 0`, this collapses to the especially simple law

`zeta_req = rho_alpha - 1.`

So on the unblocked branch the support-ratio demand is literally just the selected total loading divided by the mixed baseline, minus one.

---

## 4. What this changes conceptually

This stage removes one layer of apparent complexity from the remaining theorem gap.

Before this step, the outgoing quadrupole branch seemed to require us to know separately

- the selected conservative stiffness `lambda_-`,
- the overlap factor `s_-`,
- the transfer factor `beta_0`,
- the source map `mhat_-`,
- the demanded product `Pi_tr`,
- and the mixed baseline `C_mix`.

After the exact cancellations, the explicit support theorem depends on those objects only through the single ratio

`rho_alpha = alpha_req / alpha_mix`.

So the moving-throat support/source side is not being asked to solve the full outgoing normalization problem independently. Once the normalized quadrupole branch is imposed, it is being asked only one thing:

> **how much larger is the selected total directional loading than the mixed baseline?**

That is the right direct bridge between the selected quadrupole branch and the explicit Family-1 support/source program.

---

## 5. Best current theorem statement after Stage 085

Once the outgoing quadrupole branch satisfies the selected-mode normalization target,

`mhat_-^2 Gamma_{5,-} = 2 G / (5 c^5),`

the explicit support-demand problem is equivalent to the single ratio test

`rho_alpha = alpha_req / alpha_mix`.

More precisely,

`zeta_req = zeta_req(rho_alpha, eps_blk)`

with

`zeta_req(rho_alpha, eps_blk)`
`= (rho_alpha - 1) / [ 1 - eps_blk (2 - rho_alpha) ].`

So the remaining support theorem gap is no longer “determine several separate outgoing coefficients.”
It is:

> **determine the normalized selected-branch loading ratio `rho_alpha` of the actual passive/outgoing moving-throat quadrupole branch.**
