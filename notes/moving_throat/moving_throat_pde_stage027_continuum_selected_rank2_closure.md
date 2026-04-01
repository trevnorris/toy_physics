# Moving-Throat PDE — Stage 27: Continuum-Selected Rank-2 Closure and the Exact Quadratic Branch Equation

## Purpose

Stage 26 extracted the actual support/BdG loading data from the continuum operator:

- the support direction factor `R_phi`,
- and the physical support baseline `M_supp`.

That means the Stage-24/25 rank-2 problem is no longer a family of possible closures. The continuum kernel now supplies the data needed to write the **actual physical selected-branch equations**.

This stage does that.

The main result is that the continuum-selected wall branch is pinned by one exact quadratic equation for the softening depth `xi`, and the 2.5PN normalization test then becomes one exact residual equation evaluated on that root.

So the selected-branch problem is now no longer “support-tied or tracking?” in the abstract. It is an explicit continuum-selected closure with two exact special surfaces:

- the **minimal-kernel source-tied surface**, and
- the **interference-matched tracking surface**.

Generically, the physical kernel lands in between them.

---

## 1. Exact continuum-selected branch data

Collect the actual continuum outputs from Stages 22 and 26.

### Mixed lane

`M_mix`
`= 8 Z_W (1 + rho_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_W^(split)) ],`

`R_U = [ 1 + rho_0/(1+delta_U) ] / (1 + rho_0).`

### Support/BdG lane

`M_supp`
`= 8 Z_phi (1 + sigma_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_phi^(split)) ],`

`R_phi = [ 1 + sigma_0/(1+delta_U) ] / (1 + sigma_0).`

### Wall/source geometry

The source vector remains the original D/N direction

`s = s_0 (1,t)^T,`

with

`t = kappa_1 / kappa_0,`

`t^2 = lambda_0 = 2/9.`

So the actual mixed and support direction ratios entering the selected-branch problem are

`q = t R_U,`

`r = t R_phi.`

The anisotropy ratio is still the Stage-22 shifted wall value `delta = delta_split`.

---

## 2. Exact physical selected-branch equation

Stage 24 gave the exact support-loading theorem

`n_req(xi,delta;m,q,r)`
`= [ xi(delta + xi) - m( delta + (1 + q^2) xi ) ]`
`  / [ delta + (1 + r^2) xi - m (q-r)^2 ].`

The actual physical branch is obtained by setting

`m = M_mix,`

`n_req = M_supp,`

`q = t R_U,`

`r = t R_phi.`

So the exact continuum-selected branch equation is

`M_supp`
`= [ xi(delta + xi) - M_mix( delta + (1 + lambda_0 R_U^2) xi ) ]`
`  / [ delta + (1 + lambda_0 R_phi^2) xi - M_mix lambda_0 (R_U - R_phi)^2 ].`

This is already a complete exact statement of the selected wall branch.

---

## 3. Exact quadratic theorem for the softening depth

The key structural simplification is that the continuum-selected branch equation is exactly quadratic in `xi`.

Rearranging gives

`xi^2 + B_cont xi + C_cont = 0,`

with

`B_cont`
`= delta`
`  - M_mix (1 + lambda_0 R_U^2)`
`  - M_supp (1 + lambda_0 R_phi^2),`

`C_cont`
`= - delta (M_mix + M_supp)`
`  + lambda_0 M_mix M_supp (R_U - R_phi)^2.`

So the physical softening depth is

`xi_phys`
`= [ - B_cont + sqrt( B_cont^2 - 4 C_cont ) ] / 2,`

where the `+` branch is selected because it reduces continuously to

`xi_phys = 0`

when `M_mix = M_supp = 0`.

This is the exact quadratic selected-branch theorem.

So the physical branch is now no longer a search over a continuous family. It is an explicit algebraic root of the continuum kernel data.

---

## 4. Exact continuum-selected normalization test

Stage 25 gave the exact rank-2 normalization law

`F_(q,r,t)(xi,delta;m)`
`= [ delta + (1 + q r) xi ]^2`
`  [ delta + (1 + r t) xi - m(q-r)(q-t) ]^2`
`  / [ (1 - xi) D_(q,r)^2 ],`

with

`D_(q,r)`
`= [ delta + xi - m q(q-r) ]^2 + [ m(q-r) + r xi ]^2.`

Substituting the continuum-selected data gives

`F_cont(xi)`
`= [ delta + (1 + lambda_0 R_U R_phi) xi ]^2`
`  [ delta + (1 + lambda_0 R_phi) xi`
`    - lambda_0 M_mix (R_U - R_phi)(R_U - 1) ]^2`
`  / [ (1 - xi) D_cont(xi)^2 ],`

with

`D_cont(xi)`
`= [ delta + xi - lambda_0 M_mix R_U (R_U - R_phi) ]^2`
`  + lambda_0 [ M_mix (R_U - R_phi) + R_phi xi ]^2.`

The actual normalization theorem gate is now simply

`R_target = F_cont( xi_phys ).`

So the whole selected-branch problem has collapsed to two exact scalar equations:

1. the quadratic branch-selection equation for `xi_phys`,
2. the normalization residual equation `R_target - F_cont(xi_phys) = 0`.

---

## 5. Exact special surfaces

### 5.1 Minimal-kernel source-tied surface

If the support-interference ratio vanishes,

`sigma_0 = 0,`

then

`R_phi = 1,`

and Stage 26 shows that the support direction is exactly source-tied.

So the **minimal** continuum kernel lands on the source-tied Stage-24/25 closure.

### 5.2 Interference-matched tracking surface

If

`g_B g_R = g_W g_S,`

equivalently

`sigma_0 = rho_0,`

then

`R_phi = R_U,`

so the support and mixed directions coincide exactly.

In that case the continuum-selected branch equation collapses to

`M_mix + M_supp = G_q(xi,delta),`

with

`q = t R_U,`

`G_q(xi,delta) = xi(delta+xi) / [ delta + (1 + q^2) xi ].`

So on the tracking surface the two rank-1 loadings merge into a single one-direction branch with exact total baseline

`M_tot = M_mix + M_supp.`

---

## 6. Exact mismatch penalty and why the generic branch is intermediate

The direction mismatch is

`Delta_R := R_U - R_phi`
`= delta_U (sigma_0 - rho_0)`
`  / [ (1+delta_U)(1+rho_0)(1+sigma_0) ].`

The quadratic branch equation shows that the genuine rank-2 penalty enters through

`lambda_0 M_mix M_supp Delta_R^2`

in `C_cont`.

So the mismatch penalty is exact, positive, and quadratic in the support/mixed direction split.

This yields the structural conclusion:

- source-tied is one exact special surface,
- tracking is another exact special surface,
- and the generic extended continuum kernel sits between them with a genuine positive rank-2 mismatch penalty.

So the physical kernel does not generically choose one of the two Stage-24 extremes. It defines a continuum-selected intermediate closure.

---

## 7. Best current theorem statement after Stage 27

The support-direction bottleneck is now resolved at the reduced-theorem level.

### What is now exact

- the actual continuum-selected support direction `R_phi`,
- the actual continuum-selected support baseline `M_supp`,
- the exact quadratic physical branch equation for `xi_phys`,
- the exact continuum-selected normalization function `F_cont`,
- the minimal-kernel source-tied surface,
- the interference-matched tracking surface,
- and the exact mismatch penalty `lambda_0 M_mix M_supp (R_U - R_phi)^2`.

### The sharp current conclusion

> The minimal Stage-20 continuum kernel lands on the source-tied closure.
> The first symmetry-allowed extended continuum kernel lands generically on an exact intermediate closure.
> Exact tracking is a special interference-match surface, not the generic outcome.

So the next theorem step is no longer to decide which abstract closure to use.
It is to evaluate the continuum-selected residual

`R_target - F_cont(xi_phys)`

on concrete kernel data from the moving-throat operator.
