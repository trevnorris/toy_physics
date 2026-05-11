# Moving-Throat PDE — Stage 087: Final Reduced Finish-Line for the Explicit Family-1 Outgoing Quadrupole Branch

## Purpose

Stages 65–69 compressed the full reduced moving-throat PDE to one surviving quadrupole residual, then to one explicit product window, and finally to one pure loading-ratio criterion.

This stage states the cleanest final verdict of that reduction chain.

The main conclusion is:

> **for the explicit Family-1 support/source branch, the remaining reduced theorem gap is exactly the normalized loading ratio**
>
> `rho_alpha = alpha_req / alpha_mix`
>
> **of the actual passive/outgoing moving-throat quadrupole branch.**

Everything else on the explicit support/source side is now fixed.

---

## 1. What has completely dropped out

Once the selected outgoing quadrupole branch is required to satisfy

`mhat_-^2 Gamma_{5,-} = 2 G / (5 c^5),`

Stages 68–69 show that the explicit Family-1 support theorem no longer depends separately on

- the selected conservative overlap `s_-`,
- the selected conservative stiffness `lambda_-`,
- the outgoing transfer factor `beta_0`,
- the source-map normalization `mhat_-`,
- the product variables `Pi_tr` and `C_mix`,
- or the intermediate Peclet demand `Pe_req`.

All of those collapse to the single ratio

`rho_alpha = alpha_req / alpha_mix`.

So the support/source side of the reduced moving-throat PDE is no longer a multivariable closure problem.
It is a one-number test.

---

## 2. Exact final Family-1 criterion

For the explicit Family-1 branch, the exact normalized success/failure theorem is

`rho_alpha <= rho_suff^(chi)(lambda_mu;eps_blk)`  -> guaranteed success,

`rho_alpha >= rho_fail^(chi)(lambda_mu;eps_blk)`  -> guaranteed failure,

with the hard constructive ceiling

`rho_alpha < rho_max^(F1)(eps_blk)`.

At the natural shell-weighted normalization `lambda_mu = 1` and in the unblocked limit,

`rho_alpha <= 3.46622291347846`  -> guaranteed success,

`rho_alpha >= 3.46752913273870`  -> guaranteed failure,

`rho_alpha < 3.46752922945601`   -> absolute constructive ceiling.

So the explicit support/source branch is already fully solved in reduced form.

---

## 3. What is still genuinely open

The only remaining reduced question is not on the support/source side anymore.
It is on the outgoing quadrupole side:

> **what value of `rho_alpha = alpha_req/alpha_mix` is actually selected by the passive/outgoing moving-throat quadrupole branch?**

That value still requires the actual branch solve.

So the program is now in the cleanest state it has reached so far:

- the explicit support/source side is finished,
- the explicit Family-1 ceiling is finished,
- the wall-depth side is finished,
- the outgoing normalization factors have been shown to cancel out of the explicit support theorem,
- and the surviving reduced theorem gap is exactly one normalized loading ratio on the physical quadrupole branch.

---

## 4. Best current expert verdict

At this point, the full reduced moving-throat PDE is no longer missing a large qualitative ingredient.
It is missing one sharp quantitative datum.

The strongest honest statement is therefore:

> **within the explicit Family-1 support/source branch, the reduced moving-throat PDE program is complete up to one final passive/outgoing quadrupole loading-ratio calculation.**

That is as close as the reduced program can get to a full explicit write-up without solving the last outgoing branch itself.
