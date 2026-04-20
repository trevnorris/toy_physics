# Moving-Throat PDE — Stage 66: Direct Operator-Selected Family-1 Window for the Surviving Quadrupole Branch

## Purpose

Stage 65 compressed the whole remaining reduced moving-throat PDE into one master residual

`R_quad = zeta_req - zeta_phys(Pe_*),`

with `Pe_*` selected by the exact support/source fixed-point equation.

The next honest step is to evaluate that master residual on the first explicit branch itself, rather than continuing to carry the Family-1 window through intermediate handoffs.

This stage does that.

The main results are:

1. the operator-selected transport bias is monotone in the wall/source strength `Xi`,
2. therefore the explicit Family-1 support ratio and the corresponding branch-product thresholds are also monotone in `Xi`,
3. inserting the natural Family-1 shell-weighted and Jensen wall data reproduces the Stage-61/63/64 windows directly from the coupled operator,
4. and the natural branch already lies extremely close to the hard Family-1 ceiling.

So the explicit Family-1 branch is now fully closed at the operator-selected level.

---

## 1. Exact monotonicity of the operator-selected branch

The stationary branch is selected by

`Pe_* = Xi Delta(Pe_*;kappa,eta).`

Differentiating implicitly gives

`dPe_*/dXi = Delta(Pe_*;kappa,eta) / [ 1 - Xi partial_Pe Delta(Pe_*;kappa,eta) ].`

On the constructive stable branch,

`Delta > 0,`

`partial_Pe Delta > 0,`

and the fixed-point branch remains stable only while

`1 - Xi partial_Pe Delta > 0.`

Therefore

`dPe_*/dXi > 0.`

Because the explicit lowest-lane support ratio

`zeta_phys(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)`

is strictly increasing in `Pe`, the operator-selected support ratio is strictly increasing in `Xi`.

So stronger wall/source gain always pushes the explicit Family-1 branch upward toward its hard ceiling.

---

## 2. Exact Family-1 operator data

For the explicit Family-1 branch,

`kappa_F1 = 12321/5,`

`eta_F1 = 37,`

and the exact support/source strength reduces to

`Xi_F1 = W_wall = 1369 Upsilon_w = 136900 Theta_w.`

Using the exact Stage-41 support/source brackets,

`Pe_-^(F1) = Xi_F1 Delta_0(kappa_F1,eta_F1),`

`Pe_+^(F1) = Xi_F1 Delta_inf(kappa_F1,eta_F1),`

with

`Delta_0(kappa_F1,eta_F1) ≈ 1.73302079021525e-4,`

`Delta_inf(kappa_F1,eta_F1) ≈ 2.01447565540522e-2.`

So the explicit operator-selected Family-1 support ratio lies between

`zeta_-^(F1)(Xi_F1) = zeta_phys(Pe_-^(F1),37;12321/5),`

and

`zeta_+^(F1)(Xi_F1) = zeta_phys(Pe_+^(F1),37;12321/5).`

---

## 3. Natural shell-weighted and conservative floor windows

Stage 60 gave the exact wall-depth data

`Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2,`

`Theta_w^(J)   ≈ 0.927552032539308 lambda_mu^2.`

Therefore the operator-selected Family-1 strengths are

`Xi_F1^(chi) ≈ 556995.768726174 lambda_mu^2,`

`Xi_F1^(J)   ≈ 126981.873254631 lambda_mu^2.`

At `lambda_mu = 1`, the exact bracketed transport windows become

`Pe_-^(chi) ≈ 96.5285247264385,`

`Pe_+^(chi) ≈ 11220.5441626259,`

`Pe_-^(J)   ≈ 22.0062226330754,`

`Pe_+^(J)   ≈ 2558.01892349205.`

These are exactly the Stage-61 transport windows, now derived directly from the master support/source operator.

---

## 4. Direct operator-selected support-ratio windows

Inserting those operator-selected transport bounds into the exact Family-1 support map gives

`zeta_-^(chi)(1) ≈ 2.46622291347846,`

`zeta_+^(chi)(1) ≈ 2.46752913273870,`

`zeta_-^(J)(1)   ≈ 2.44257571477179,`

`zeta_+^(J)(1)   ≈ 2.46752736855058.`

Compare these with the hard Family-1 ceiling

`zeta_max^(F1) ≈ 2.46752922945601.`

So on the natural shell-weighted branch the guaranteed-success lower window already sits less than about `1.31e-3` below the hard ceiling, while the guaranteed-failure upper window differs from the ceiling only in the seventh decimal place.

This is the operator-selected version of the earlier explicit-branch verdict:

> the Family-1 branch is not support-starved; it already drives the explicit support ratio essentially to its maximal constructive value.

---

## 5. Direct operator-selected branch-product windows

Using the exact inverse map

`Pi = C_mix Q(zeta;eps_blk),`

`Q(zeta;eps_blk) = [1 + (1 - 2 eps_blk) zeta] / [1 - eps_blk zeta],`

these operator-selected support windows become direct selected-branch product windows

`Pi_suff^(F1) = C_mix Q( zeta_-^(F1); eps_blk ),`

`Pi_fail^(F1) = C_mix Q( zeta_+^(F1); eps_blk ).`

At vanishing blocking `eps_blk = 0`, these reduce simply to

`Pi_suff^(F1) / C_mix = 1 + zeta_-^(F1),`

`Pi_fail^(F1) / C_mix = 1 + zeta_+^(F1).`

So for `lambda_mu = 1`,

`Pi_suff^(chi)/C_mix ≈ 3.46622291347846,`

`Pi_fail^(chi)/C_mix ≈ 3.46752913273870,`

`Pi_suff^(J)/C_mix   ≈ 3.44257571477179,`

`Pi_fail^(J)/C_mix   ≈ 3.46752736855058,`

while the hard ceiling remains

`Pi_max^(F1)/C_mix ≈ 3.46752922945601.`

So the explicit operator-selected Family-1 window in the actual quadrupole branch variable is now fully closed.

---

## 6. What Stage 66 changes

This stage removes the last remaining ambiguity in the explicit Family-1 branch phase.

The earlier stages had already shown that support depth and source asymmetry could in principle reach the required window. What was still missing was a direct operator-selected statement.

Now we have it.

The coupled support/source operator itself selects a branch window that is already essentially saturated at the Family-1 ceiling once the natural shell-weighted wall datum is inserted.

So the remaining reduced theorem gap is no longer on the explicit support/source side at all.
It is entirely on the outgoing quadrupole-normalization side:

> what value of `Pi_tr` does the actual passive/outgoing moving-throat quadrupole branch produce, and does it stay below the explicit Family-1 ceiling? 
