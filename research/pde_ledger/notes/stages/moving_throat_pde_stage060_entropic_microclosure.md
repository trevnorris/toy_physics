# Moving-Throat PDE — Stage 060: Entropic Source Microclosure and the Microscopic Support/Source Gain

## Purpose

Stage 058 introduced the coupled support/source branch equation

`Pe = Xi Delta(Pe;kappa,eta),`

but the coupling `Xi` was still phenomenological:

`Xi = mu_sigma Lambda_phi^2 L^2 / (D_sigma T_X).`

The next honest step is to derive that quantity from the first explicit microscopic closure that is still compatible with the parent 4D ontology:

- exact projected continuity for the source density,
- a positive scalar source density carried along the throat axis,
- an explicit support/source free-energy coupling,
- and an Onsager/Fokker–Planck transport law that preserves positivity.

This stage does that.

The main results are:

1. the first explicit microscopic source/support free energy can be written as

   `F[sigma,phi]`
   `= int_0^L ds [ Theta_sigma sigma (log(sigma/sigma_*) - 1) - Lambda_phi sigma phi`
   `               + (T_X/2) phi_s^2 + (K_X/2) phi^2 ] + (K_m/2) phi(0)^2;`
2. its Euler–Lagrange variations are

   `mu_sigma^(chem) = delta F / delta sigma = Theta_sigma log(sigma/sigma_*) - Lambda_phi phi,`

   `-T_X phi_ss + K_X phi = Lambda_phi sigma,`

   with Robin/Neumann support boundary conditions

   `T_X phi_s(0)=K_m phi(0),` `phi_s(L)=0;`
3. the minimal positive-density Onsager current is

   `J = -M_sigma sigma partial_s mu_sigma^(chem)`
   `  = -D_sigma partial_s sigma + M_sigma Lambda_phi sigma partial_s phi,`

   where the exact Einstein relation is

   `D_sigma = M_sigma Theta_sigma;`
4. under the same affine-drop reduction already implicit in the Stage-058 average-drift closure,

   `phi(s) ~= phi(0) + [Delta phi] s/L,`

   the stationary zero-flux branch is exactly the Stage-056 exponential family,

   `sigma(s) = C exp[(Lambda_phi Delta phi)/(Theta_sigma L) s],`

   or, in normalized coordinates `x=s/L`,

   `Sigma_Pe(x) = Pe exp(Pe x)/(exp(Pe)-1);`
5. the transport bias is therefore no longer phenomenological:

   `Pe = (Lambda_phi/Theta_sigma) Delta phi;`
6. using the Stage-058 support normalization

   `Phi = T_X phi/(Lambda_phi L^2),`
   `Delta phi = (Lambda_phi L^2 / T_X) Delta(Pe;kappa,eta),`

   the branch equation becomes

   `Pe = Xi_micro Delta(Pe;kappa,eta),`

   with the exact microscopic coupling

   `Xi_micro = Lambda_phi^2 L^2 / (Theta_sigma T_X)`
   `        = chi_sigma Lambda_phi^2 L^2 / T_X,`

   where `chi_sigma := 1/Theta_sigma;`
7. equivalently, the phenomenological Stage-058 coupling is now explained as

   `Xi = mu_sigma Lambda_phi^2 L^2/(D_sigma T_X)`
   `  = Lambda_phi^2 L^2/(Theta_sigma T_X)`

   by the Einstein relation `D_sigma = mu_sigma Theta_sigma;`
8. the source dynamics is passive: if `phi` is held fixed or slaved quasi-statically by the support Euler–Lagrange equation, then

   `dF/dt = - int_0^L ds J^2/(M_sigma sigma) <= 0`

   under no-flux boundaries.

So the operator strength is no longer an arbitrary ratio of mobility and diffusion. On the first explicit microscopic closure, it collapses to one entropic-support gain

`Xi_micro = chi_sigma Lambda_phi^2 L^2/T_X.`

---

## 1. Microscopic source/support free energy

Take the throat-axis source density `sigma(s,t) > 0` on `s in [0,L]` and the support field `phi(s,t)`.

The first explicit free-energy functional that preserves positivity of `sigma` and reproduces the Stage-058 support operator is

`F[sigma,phi]`
`= int_0^L ds [ Theta_sigma sigma (log(sigma/sigma_*) - 1) - Lambda_phi sigma phi`
`               + (T_X/2) phi_s^2 + (K_X/2) phi^2 ] + (K_m/2) phi(0)^2.`

Here:

- `Theta_sigma` is the source entropic/chemical scale,
- `sigma_*` is a reference density,
- `Lambda_phi` is the support/source coupling,
- `T_X` is the axial support tension,
- `K_X` is the baseline support stiffness,
- `K_m` is the Robin mouth spring.

Variation with respect to `sigma` gives the exact chemical potential

`mu_sigma^(chem) = Theta_sigma log(sigma/sigma_*) - Lambda_phi phi.`

Variation with respect to `phi` gives the support equation

`-T_X phi_ss + K_X phi = Lambda_phi sigma,`

with boundary conditions

`T_X phi_s(0)=K_m phi(0),`

`phi_s(L)=0.`

So the Stage-058 support law is now embedded directly into one explicit free-energy closure.

---

## 2. Positive-density Onsager transport law

Use the minimal gradient-flow current

`J = -M_sigma sigma partial_s mu_sigma^(chem),`

with mobility `M_sigma > 0`.

Expanding the chemical potential gradient gives

`J = -M_sigma sigma [ Theta_sigma partial_s log sigma - Lambda_phi partial_s phi ]`
`  = -M_sigma Theta_sigma partial_s sigma + M_sigma Lambda_phi sigma partial_s phi.`

So the transport law takes the drift-diffusion form

`J = -D_sigma partial_s sigma + M_sigma Lambda_phi sigma partial_s phi,`

with exact Einstein relation

`D_sigma = M_sigma Theta_sigma.`

This is the microscopic origin of the Stage-058 phenomenological pair `(D_sigma,mu_sigma)`.

---

## 3. Recovery of the Stage-056 exponential family

Stage 058 already used an average-drift closure in which the support field enters only through the end-to-end drop.

Within that same lowest-lane closure, replace the support profile by its affine-drop reduction

`phi(s) ~= phi(0) + [Delta phi] s/L,`

so that `partial_s phi ~= Delta phi / L` is constant.

On the stationary zero-flux branch `J=0`, the exact transport law gives

`partial_s sigma = (Lambda_phi Delta phi)/(Theta_sigma L) sigma.`

Therefore

`sigma(s) = C exp[(Lambda_phi Delta phi)/(Theta_sigma L) s].`

Normalize on `[0,L]` and pass to `x=s/L in [0,1]`. Then

`Sigma_Pe(x) = L sigma(Lx)`
`           = Pe exp(Pe x)/(exp(Pe)-1),`

with the exact Peclet number

`Pe = (Lambda_phi/Theta_sigma) Delta phi.`

So the Stage-056 exponential branch is not just a convenient guess. It is the exact stationary branch of the first positive-density Onsager closure after the same affine-drop reduction already used implicitly in Stage 058.

---

## 4. Exact microscopic coupling `Xi_micro`

Stage 058 defined the dimensionless support field by

`Phi = T_X phi/(Lambda_phi L^2),`

so the end-to-end physical support drop is

`Delta phi = (Lambda_phi L^2/T_X) Delta(Pe;kappa,eta).`

Insert this into the exact Peclet law above:

`Pe = (Lambda_phi/Theta_sigma) Delta phi`
`   = (Lambda_phi^2 L^2/(Theta_sigma T_X)) Delta(Pe;kappa,eta).`

Therefore the Stage-058 branch equation is recovered with the exact microscopic coupling

`Xi_micro = Lambda_phi^2 L^2/(Theta_sigma T_X).`

Equivalently, with the source susceptibility

`chi_sigma := 1/Theta_sigma,`

one has

`Xi_micro = chi_sigma Lambda_phi^2 L^2/T_X.`

And because `D_sigma = M_sigma Theta_sigma`, this also reproduces the Stage-058 phenomenological form:

`Xi_micro = mu_sigma Lambda_phi^2 L^2/(D_sigma T_X).`

So the apparent mobility/diffusion ambiguity is gone. The coupling is one entropic-support gain.

---

## 5. Exact passivity / Lyapunov identity

The same free-energy closure supplies an exact dissipation law.

At fixed `phi` one has

`dF/dt = int_0^L ds mu_sigma^(chem) partial_t sigma.`

Using continuity `partial_t sigma = -partial_s J` gives

`dF/dt = - int_0^L ds mu_sigma^(chem) partial_s J`
`     = - [ mu_sigma^(chem) J ]_0^L + int_0^L ds (partial_s mu_sigma^(chem)) J.`

Now insert the Onsager current

`J = -M_sigma sigma partial_s mu_sigma^(chem).`

Then

`(partial_s mu_sigma^(chem)) J = - J^2/(M_sigma sigma).`

So the exact free-energy identity is

`dF/dt = - [ mu_sigma^(chem) J ]_0^L - int_0^L ds J^2/(M_sigma sigma).`

Under no-flux boundaries `J(0)=J(L)=0`,

`dF/dt = - int_0^L ds J^2/(M_sigma sigma) <= 0.`

If `phi` is slaved quasi-statically by the support Euler–Lagrange equation, the same identity applies to the full coupled free energy because the `phi`-variation term vanishes on shell.

So this microscopic closure is automatically passive.

---

## 6. What Stage 060 changes

The operator problem has now advanced in a concrete way.

Before Stage 060, the support/source strength was the phenomenological ratio

`Xi = mu_sigma Lambda_phi^2 L^2/(D_sigma T_X).`

After Stage 060, the same quantity is an explicit microscopic gain:

`Xi_micro = chi_sigma Lambda_phi^2 L^2/T_X.`

That is a much tighter result because:

- it is independent of the separate values of mobility and diffusion,
- it is tied to one explicit free-energy closure,
- it preserves positivity of the source density,
- it reproduces the Stage-056 exponential family exactly in the same lowest-lane reduction already used earlier,
- and it is automatically passive.

So the remaining theorem gap is now sharper again:

> compute the actual moving-throat values of `chi_sigma`, `Lambda_phi`, `T_X`, and `L`, form `Xi_micro`, and compare it directly to the exact support thresholds from Stage 059.
