# Moving-Throat PDE — Stage 081: Final Explicit Family-1 Quadrupole-Demand Window in the Branch-Product Variable `Pi_tr`

## Purpose

Stage 63 translated the explicit Family-1 wall-depth verdict into the support-ratio demand variable `zeta_req`.

The last honest step in this explicit-branch phase is to put the result back into the exact quadrupole branch-product language used by Stages 34–35.

This stage does that.

The main result is that, for the explicit Family-1 branch, the final selected quadrupole branch is constrained by exact thresholds of the form

`Pi_tr <= Pi_suff^(chi)(lambda_mu;eps_blk)`  -> guaranteed success,

`Pi_tr >= Pi_fail^(chi)(lambda_mu;eps_blk)`  -> guaranteed failure,

with a hard explicit Family-1 ceiling

`Pi_tr < Pi_max^(F1)(eps_blk)`

required for the branch to be reachable at all.

At vanishing blocking (`eps_blk = 0`) and with the natural wall normalization `lambda_mu = 1`, the explicit branch already gives

`Pi_suff^(chi) / C_mix ≈ 3.46622291347846,`

while the hard ceiling is only

`Pi_max^(F1) / C_mix ≈ 3.46752922945601.`

So the entire remaining explicit Family-1 theorem gap is now compressed to a very narrow `Pi_tr / C_mix` window.

---

## 1. Exact inversion of the Stage-35 support-demand law

Stage 35 gave

`zeta_req = (Pi_tr - C_mix) / [ C_mix - eps_blk (2 C_mix - Pi_tr) ],`

with `eps_blk` the blocking ratio carried from the selected quadrupole branch.

Solving this exactly for `Pi_tr` gives

`Pi_tr = C_mix Q(zeta_req;eps_blk),`

where

`Q(zeta;eps_blk)`
`:= [ 1 + (1 - 2 eps_blk) zeta ] / [ 1 - eps_blk zeta ].`

This map is exact.

Its anchor values are

`Q(0;eps_blk) = 1,`

`Q(1;eps_blk) = 2,`

and its derivative is

`dQ/dzeta = (1 - eps_blk) / (1 - eps_blk zeta)^2 > 0`

whenever `eps_blk zeta < 1`.

So the product demand grows strictly with the support-ratio demand.

---

## 2. Explicit Family-1 success and failure thresholds in `Pi_tr`

Insert the Stage-63 Family-1 thresholds.

Define

`Pi_suff^(chi)(lambda_mu;eps_blk)`
`:= C_mix Q( zeta_suff^(chi)(lambda_mu); eps_blk ),`

`Pi_fail^(chi)(lambda_mu;eps_blk)`
`:= C_mix Q( zeta_fail^(chi)(lambda_mu); eps_blk ),`

and similarly for the conservative floor,

`Pi_suff^(J) := C_mix Q( zeta_suff^(J); eps_blk ),`

`Pi_fail^(J) := C_mix Q( zeta_fail^(J); eps_blk ).`

Then the explicit Family-1 branch theorem is

`Pi_tr <= Pi_suff^(chi)`  -> guaranteed success,

`Pi_tr >= Pi_fail^(chi)`  -> guaranteed failure,

with the conservative version obtained by replacing `(chi)` with `(J)`.

So the Family-1 explicit-branch result is now written directly in the same branch-product variable that the quadrupole normalization program actually uses.

---

## 3. Hard explicit Family-1 ceiling in product form

Stage 62 gave the hard support-ratio ceiling

`zeta_req <= zeta_max^(F1) ≈ 2.46752922945601.`

Therefore the explicit Family-1 branch can only be reached if

`Pi_tr < Pi_max^(F1)(eps_blk)`

with

`Pi_max^(F1)(eps_blk)`
`:= C_mix Q( zeta_max^(F1); eps_blk )`
` = C_mix [ 1 + (1 - 2 eps_blk) zeta_max^(F1) ] / [ 1 - eps_blk zeta_max^(F1) ].`

This requires the denominator to remain positive, i.e.

`eps_blk < 1 / zeta_max^(F1) ≈ 0.405263689711371.`

So even before the final PDE-side quadrupole normalization is solved, the explicit Family-1 branch has an exact finite product ceiling.

---

## 4. Numerical illustration at `eps_blk = 0`

If the final quadrupole branch stays close to the unblocked limit `eps_blk = 0`, then

`Q(zeta;0) = 1 + zeta.`

So at `lambda_mu = 1` the explicit thresholds become

`Pi_suff^(chi) / C_mix = 1 + zeta_suff^(chi)(1)`
`                      ≈ 3.46622291347846,`

`Pi_fail^(chi) / C_mix = 1 + zeta_fail^(chi)(1)`
`                      ≈ 3.46752913273870,`

while the hard Family-1 ceiling is

`Pi_max^(F1) / C_mix = 1 + zeta_max^(F1)`
`                    ≈ 3.46752922945601.`

For the conservative lower envelope,

`Pi_suff^(J) / C_mix ≈ 3.44257571477179,`

`Pi_fail^(J) / C_mix ≈ 3.46752736855058.`

So the natural shell-weighted Family-1 branch at `lambda_mu = 1` already pushes the explicit product threshold essentially all the way to the hard ceiling.

---

## 5. Final explicit verdict for this phase

The explicit Family-1 branch now has a completely closed reduced theorem statement.

- Stage 61 showed that wall-depth supply is not the dominant open issue.
- Stage 62 converted the remaining bottleneck into a hard support-ratio ceiling `zeta_max^(F1)`.
- Stage 63 translated the Stage-61 wall verdict into the demanded support-ratio variable `zeta_req`.
- This stage finally expresses the result directly in the selected-branch product variable `Pi_tr`.

So the remaining explicit theorem gap is now as narrow as it can be without solving the final outgoing quadrupole normalization branch itself:

> does the completed moving-throat quadrupole branch produce a selected product `Pi_tr` that stays below the explicit Family-1 ceiling `Pi_max^(F1)(eps_blk)` and, more sharply, below the natural success threshold `Pi_suff^(chi)(lambda_mu;eps_blk)`?

That is the clean finish line for the present explicit-branch phase.
