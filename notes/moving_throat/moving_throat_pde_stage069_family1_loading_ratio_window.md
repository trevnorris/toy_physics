# Moving-Throat PDE — Stage 69: Exact Family-1 Success/Failure Window in the Pure Loading-Ratio Variable

## Purpose

Stage 68 showed that, once the selected outgoing quadrupole branch is normalized, the explicit support-demand problem depends only on the loading ratio

`rho_alpha := alpha_req / alpha_mix`.

The next honest step is to rewrite the entire explicit Family-1 theorem directly in that variable.

This stage does that.

The main result is that the explicit Family-1 branch no longer needs the product variables `(Pi_tr, C_mix)` at all.
It is governed exactly by the ratio window

`rho_alpha <= rho_suff^(chi)(lambda_mu;eps_blk)`  -> guaranteed success,

`rho_alpha >= rho_fail^(chi)(lambda_mu;eps_blk)`  -> guaranteed failure,

with a hard constructive ceiling

`rho_alpha < rho_max^(F1)(eps_blk).`

At `eps_blk = 0` and `lambda_mu = 1`, the natural shell-weighted Family-1 window is simply

`rho_alpha <= 3.46622291347846`  -> guaranteed success,

`rho_alpha >= 3.46752913273870`  -> guaranteed failure,

with the hard ceiling

`rho_alpha < 3.46752922945601.`

So the whole explicit Family-1 support theorem has collapsed to a very narrow pure loading-ratio window.

---

## 1. Exact ratio form of the Stage-64 product map

Stage 64 gave the exact product inversion

`Pi_tr = C_mix Q(zeta;eps_blk),`

with

`Q(zeta;eps_blk)`
`= [ 1 + (1 - 2 eps_blk) zeta ] / [ 1 - eps_blk zeta ].`

By Stage 68,

`Pi_tr / C_mix = rho_alpha = alpha_req / alpha_mix.`

Therefore the explicit Family-1 thresholds are immediately

`rho_suff^(chi)(lambda_mu;eps_blk)`
`:= Q( zeta_suff^(chi)(lambda_mu) ; eps_blk ),`

`rho_fail^(chi)(lambda_mu;eps_blk)`
`:= Q( zeta_fail^(chi)(lambda_mu) ; eps_blk ),`

and similarly for the conservative lower envelope,

`rho_suff^(J) := Q( zeta_suff^(J) ; eps_blk ),`

`rho_fail^(J) := Q( zeta_fail^(J) ; eps_blk ).`

The hard Family-1 constructive ceiling is

`rho_max^(F1)(eps_blk)`
`:= Q( zeta_max^(F1) ; eps_blk )`
` = [ 1 + (1 - 2 eps_blk) zeta_max^(F1) ] / [ 1 - eps_blk zeta_max^(F1) ].`

So the explicit support theorem is now entirely a statement about `rho_alpha`.

---

## 2. Exact unblocked Family-1 window

In the unblocked limit,

`eps_blk = 0,`

so

`Q(zeta;0) = 1 + zeta.`

Therefore the explicit Family-1 branch reduces to the exact unblocked loading-ratio thresholds

`rho_suff^(chi) = 1 + zeta_suff^(chi),`

`rho_fail^(chi) = 1 + zeta_fail^(chi),`

`rho_max^(F1) = 1 + zeta_max^(F1).`

Using the Stage-63/64 values at `lambda_mu = 1`,

`zeta_suff^(chi)(1) ≈ 2.46622291347846,`

`zeta_fail^(chi)(1) ≈ 2.46752913273870,`

`zeta_suff^(J)(1)   ≈ 2.44257571477179,`

`zeta_max^(F1)      ≈ 2.46752922945601,`

gives

`rho_suff^(chi)(1;0) ≈ 3.46622291347846,`

`rho_fail^(chi)(1;0) ≈ 3.46752913273870,`

`rho_suff^(J)(1;0)   ≈ 3.44257571477179,`

`rho_max^(F1)(0)     ≈ 3.46752922945601.`

So on the natural shell-weighted branch the guaranteed-success threshold lies only about

`1.30631597755e-3`

below the hard constructive ceiling, and the guaranteed-failure threshold differs from the ceiling only in the seventh decimal place.

---

## 3. Exact blocking condition in ratio form

The same denominator condition from Stage 64 remains necessary:

`eps_blk < 1 / zeta_max^(F1) ≈ 0.405263689711371.`

So the explicit Family-1 loading-ratio ceiling exists only while the blocked branch stays below that exact limit.

On any admissible branch,

`d rho_max^(F1) / d eps_blk`
`= zeta_max^(F1) [ zeta_max^(F1) - 1 ] / [ 1 - eps_blk zeta_max^(F1) ]^2 > 0,`

so blocking always raises the required loading ratio for a given support-ratio demand.

In other words, blocking hurts the support theorem exactly the way the product formulation suggested, but the statement is now written directly in the true reduced variable `rho_alpha`.

---

## 4. Final explicit Family-1 theorem in its cleanest form

The explicit Family-1 moving-throat support/source branch succeeds on the normalized outgoing quadrupole branch iff the selected loading ratio satisfies

`rho_alpha < rho_max^(F1)(eps_blk),`

and, more sharply, it is guaranteed to succeed once

`rho_alpha <= rho_suff^(chi)(lambda_mu;eps_blk)`.

At the natural shell-weighted normalization `lambda_mu = 1` and in the unblocked limit, that means

`rho_alpha <= 3.46622291347846`  -> guaranteed success,

`rho_alpha >= 3.46752913273870`  -> guaranteed failure,

with the absolute constructive ceiling

`rho_alpha < 3.46752922945601.`

So the explicit Family-1 theorem gap is no longer about the outgoing normalization factors separately and no longer about the product variables separately.
It is just this:

> **does the actual passive/outgoing moving-throat quadrupole branch produce a normalized loading ratio `rho_alpha = alpha_req/alpha_mix` below about `3.4675`?**
