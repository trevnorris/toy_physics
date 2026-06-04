# Moving-Throat PDE — Stage 063: Parent-Overlap Threshold Theorem and Exact Microscopic Success/No-Go Tests

## Purpose

Stage 062 expressed the Stage-061 microscopic gain directly in parent 4D variables:

`G_micro = rho_* g_phi^2 O_(sigma phi)^2 / [ m c_(s*)^2 K_X N_(sigma sigma) ]`

or, equivalently,

`G_micro = [ rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) ] C_(sigma phi)^2.`

The next honest step is to combine that exact parent formula with the Stage-061 operator phase diagram and turn the support/source theorem gap into explicit thresholds on:

- the parent confinement-loading amplitude `g_phi`,
- and the source/support coherence factor `C_(sigma phi)^2`.

This stage does that.

The main results are:

1. inserting the parent gain into the Stage-061 thresholds gives the exact parent fail/succeed conditions

   `g_phi^2 <= g_(phi,fail)^2`  -> fail,

   `g_phi^2 >= g_(phi,suff)^2`  -> succeed,

   with

   `g_(phi,fail)^2`
   `= m c_(s*)^2 K_X N_(sigma sigma) G_fail / [ rho_* O_(sigma phi)^2 ],`

   `g_(phi,suff)^2`
   `= m c_(s*)^2 K_X N_(sigma sigma) G_suff / [ rho_* O_(sigma phi)^2 ];`
2. equivalently, in coherence-factor form,

   `g_(phi,fail)^2`
   `= m c_(s*)^2 K_X G_fail / [ rho_* N_(phi phi) C_(sigma phi)^2 ],`

   `g_(phi,suff)^2`
   `= m c_(s*)^2 K_X G_suff / [ rho_* N_(phi phi) C_(sigma phi)^2 ];`
3. solving instead for the required source/support alignment gives the exact coherence thresholds

   `C_fail^2 = m c_(s*)^2 K_X G_fail / [ rho_* g_phi^2 N_(phi phi) ],`

   `C_suff^2 = m c_(s*)^2 K_X G_suff / [ rho_* g_phi^2 N_(phi phi) ],`

   with

   `C_fail^2 <= C_suff^2;`
4. because `0 <= C_(sigma phi)^2 <= 1`, there is an immediate Cauchy-based no-go theorem:

   if

   `rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) < G_fail,`

   then the branch fails **for every possible** source profile in the chosen support channel;
5. the exact “best-case” reachable gain at fixed `g_phi` is therefore

   `G_max(g_phi) = rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X),`

   attained only at perfect alignment `C_(sigma phi)^2 = 1`;
6. inserting `kappa = K_X L^2 / T_X` and the Stage-061 formulas

   `G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

   `G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ],`

   the parent amplitude thresholds become

   `g_(phi,fail)^2`
   `= m c_(s*)^2 T_X N_(sigma sigma) Pe_req / [ rho_* L^2 O_(sigma phi)^2 Delta_inf(kappa,eta) ],`

   `g_(phi,suff)^2`
   `= m c_(s*)^2 T_X N_(sigma sigma) Pe_req / [ rho_* L^2 O_(sigma phi)^2 Delta_0(kappa,eta) ];`
7. in coherence-factor form this is

   `g_(phi,fail)^2`
   `= m c_(s*)^2 T_X Pe_req / [ rho_* L^2 N_(phi phi) C_(sigma phi)^2 Delta_inf ],`

   `g_(phi,suff)^2`
   `= m c_(s*)^2 T_X Pe_req / [ rho_* L^2 N_(phi phi) C_(sigma phi)^2 Delta_0 ],`

   so the prefactor `K_X` cancels exactly, surviving only through the shape function `Delta_(0,inf)(kappa,eta)`.

So the remaining PDE-side theorem gap is no longer “find some support/source gain.” It is now the exact comparison of one parent confinement-loading amplitude and one source/support coherence factor against the explicit fail/succeed thresholds above.

---

## 1. Exact parent thresholds on the confinement-loading amplitude `g_phi`

Stage 061 proved the microscopic operator phase diagram

`G_micro <= G_fail(kappa,eta)`  -> fail,

`G_micro >= G_suff(kappa,eta)`  -> succeed.

Insert the parent formula from Stage 062,

`G_micro = rho_* g_phi^2 O_(sigma phi)^2 / [ m c_(s*)^2 K_X N_(sigma sigma) ].`

Then the fail/succeed thresholds on the parent loading amplitude are exact:

`g_(phi,fail)^2`
`= m c_(s*)^2 K_X N_(sigma sigma) G_fail / [ rho_* O_(sigma phi)^2 ],`

`g_(phi,suff)^2`
`= m c_(s*)^2 K_X N_(sigma sigma) G_suff / [ rho_* O_(sigma phi)^2 ].`

Therefore

- if `g_phi^2 <= g_(phi,fail)^2`, the chosen parent branch cannot reach the target inside this operator family,
- if `g_phi^2 >= g_(phi,suff)^2`, the branch is guaranteed to reach it,
- only the interval `g_(phi,fail)^2 < g_phi^2 < g_(phi,suff)^2` still requires the full fixed-point solve.

---

## 2. Exact threshold on the profile coherence factor

Using

`O_(sigma phi)^2 = N_(sigma sigma) N_(phi phi) C_(sigma phi)^2,`

the gain can be rewritten as

`G_micro = [ rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) ] C_(sigma phi)^2.`

So the exact required coherence floors are

`C_fail^2 = m c_(s*)^2 K_X G_fail / [ rho_* g_phi^2 N_(phi phi) ],`

`C_suff^2 = m c_(s*)^2 K_X G_suff / [ rho_* g_phi^2 N_(phi phi) ].`

These control the same three-zone structure:

- if `C_(sigma phi)^2 <= C_fail^2`, fail;
- if `C_(sigma phi)^2 >= C_suff^2`, succeed;
- only the narrow intermediate interval still needs the full root solve.

Because `G_fail <= G_suff`, one has

`C_fail^2 <= C_suff^2.`

So the unresolved PDE data have been split very cleanly: one parent strength scale and one profile-alignment factor.

---

## 3. Exact Cauchy no-go theorem

The coherence factor satisfies

`0 <= C_(sigma phi)^2 <= 1`.

So at fixed parent loading amplitude `g_phi`, the **largest possible** gain in the chosen support channel is

`G_max(g_phi) = rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X),`

obtained only at perfect alignment `C_(sigma phi)^2 = 1`.

Therefore there is an exact no-go theorem:

if

`G_max(g_phi) < G_fail(kappa,eta),`

then the branch fails for **every possible** source profile in that support channel.

Equivalently,

if

`rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) < G_fail(kappa,eta),`

then no profile engineering of `chi_sigma` can rescue the branch.

This is the first exact parent-overlap no-go theorem in the program.

---

## 4. Exact amplitude thresholds in terms of `Pe_req`, `Delta_0`, and `Delta_inf`

Stage 061 already gave

`G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

`G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ],`

with

`kappa = K_X L^2 / T_X.`

Insert these into the parent amplitude thresholds.

For failure:

`g_(phi,fail)^2`
`= [ m c_(s*)^2 K_X N_(sigma sigma) / (rho_* O_(sigma phi)^2) ]`
`  [ Pe_req / (kappa Delta_inf) ]`
`= m c_(s*)^2 T_X N_(sigma sigma) Pe_req`
`  / [ rho_* L^2 O_(sigma phi)^2 Delta_inf(kappa,eta) ].`

For sufficiency:

`g_(phi,suff)^2`
`= m c_(s*)^2 T_X N_(sigma sigma) Pe_req`
`  / [ rho_* L^2 O_(sigma phi)^2 Delta_0(kappa,eta) ].`

So the explicit prefactor `K_X` cancels exactly, surviving only inside the geometry shape functions `Delta_inf` and `Delta_0` through `kappa`.

This is a very clean structural result: the amplitude threshold is tension/length/compressibility controlled, while the baseline support stiffness enters only through the support-shape response.

---

## 5. Matched-profile and normalized-lane special case

If the source and support channels coincide and are normalized,

`chi_sigma = chi_phi,`

`N_(sigma sigma)=N_(phi phi)=1,`

`O_(sigma phi)=1,`

`C_(sigma phi)^2 = 1,`

then the parent thresholds simplify to

`g_(phi,fail)^2 = m c_(s*)^2 K_X G_fail / rho_*,`

`g_(phi,suff)^2 = m c_(s*)^2 K_X G_suff / rho_*.`

Equivalently,

`g_(phi,fail)^2 = m c_(s*)^2 T_X Pe_req / [ rho_* L^2 Delta_inf(kappa,eta) ],`

`g_(phi,suff)^2 = m c_(s*)^2 T_X Pe_req / [ rho_* L^2 Delta_0(kappa,eta) ].`

So in the best possible aligned lowest-lane closure, the theorem gap reduces to a single confinement-loading amplitude compared against two explicit geometry-controlled thresholds.

---

## 6. What Stage 063 changes

After Stage 063, the remaining theorem gap is no longer phrased in terms of the abstract microscopic gain `G_micro`.

It is now:

1. compute the parent confinement-loading amplitude `g_phi`,
2. compute the profile coherence `C_(sigma phi)^2` on the true moving-throat branch,
3. compare them against the exact thresholds above.

That is stronger than the Stage-061 statement because it pushes the support/source theorem gap all the way back to parent-action overlap data.

In other words, the unresolved PDE question is now no longer “is the microscopic gain big enough?” It is

> does the completed moving-throat branch generate sufficient **parent confinement loading** and sufficient **source/support coherence** to cross the explicit fail/succeed surfaces above?
