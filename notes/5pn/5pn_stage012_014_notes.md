# 5PN continuation notes — Stage 12–13 normalized monomial bridge

These two stages connect the earlier explicit Stage-5 primitive prototype

a) directly to the Stage-10/11 microscopic similarity-orbit package, and
b) into formulas that can be used on the actual moving-throat branch.

## Stage 12 — normalized prototype already contains the Stage-11 quotient coordinates

Using the coherent-kernel dictionary

- `K = K_eta^(eff)`
- `M = mu_eta`
- `G_U = c_etaU / sqrt(mu_eta mu_U)`
- `G_W = lambda_W / sqrt(mu_eta mu_W)`
- `R = gamma lambda_W / sqrt(mu_U mu_W)`
- `Omega_U^2 = K_U / mu_U`
- `Omega_W^2 = K_W^(eff) / mu_W`
- `delta_U = pi^2 T_U / (L^2 K_U)`

one gets the exact normalized coherent ratios

\[
\chi_0 = \frac{R G_U}{\Omega_U^2 G_W},
\qquad
\epsilon_\eta = \frac{M G_U^2}{K\Omega_U^2},
\qquad
\epsilon_W = \frac{R^2\sigma}{\Omega_U^2\Omega_W^2},
\qquad
Z_W = \frac{M G_W^2}{K\Omega_W^2}.
\]

The direct Stage-168/169 monomials become

\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{R G_U}{\Omega_U^2 G_W}\right)^{1+\delta_{U,*}}
\delta_U^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\frac{M G_W^2}{K\Omega_W^4}
\left(\frac{R^2\sigma}{\Omega_U^2\Omega_W^2}\right)^{E_*}
\delta_U^{-F_*},
\]

\[
\epsilon_\eta = \frac{M G_U^2}{K\Omega_U^2}.
\]

So the Stage-5-style normalized prototype already contains the exact Stage-11 quotient coordinates once the split-`U` variable `delta_U` is admitted.

The Stage-10 slippages collapse to the normalized drift formulas

\[
\Sigma_Z = d\ln M + 2 d\ln G_W - d\ln K - 4 d\ln\Omega_W,
\]
\[
\Sigma_\chi = d\ln R + d\ln G_U - d\ln G_W - 2 d\ln\Omega_U,
\]
\[
\Sigma_\eta = d\ln M + 2 d\ln G_U - d\ln K - 2 d\ln\Omega_U,
\]
\[
\Sigma_\epsilon = 2(d\ln R - d\ln\Omega_U - d\ln\Omega_W),
\qquad
\Sigma_\delta = d\ln\delta_U.
\]

So the extra raw mass/stiffness bookkeeping mostly cancels out of the actual defect coordinates.

## Stage 13 — exact zero-defect kernel in normalized Stage-5 variables

In the normalized variables

\[
(d\ln G_W,
 d\ln G_U,
 d\ln R,
 d\ln K,
 d\ln M,
 d\ln\Omega_U,
 d\ln\Omega_W,
 d\ln\delta_U),
\]

the direct monomial-drift matrix is

\[
M_{\rm norm}=
\begin{pmatrix}
-(1+\delta_U) & 1+\delta_U & 1+\delta_U & 0 & 0 & -2(1+\delta_U) & 0 & 1+\chi_0 \\
2 & 0 & 2E_* & -1 & 1 & -2E_* & -(4+2E_*) & -F_* \\
0 & 2 & 0 & -1 & 1 & -2 & 0 & 0
\end{pmatrix}.
\]

Its rank is exactly `3`, so the normalized zero-defect tangent space has dimension `5`.

The exact zero-defect equations are

\[
(1+\delta_U)(d\ln R + d\ln G_U - d\ln G_W - 2 d\ln\Omega_U)
+
(1+\chi_0)d\ln\delta_U = 0,
\]

\[
d\ln M + 2 d\ln G_U - d\ln K - 2 d\ln\Omega_U = 0,
\]

\[
d\ln M + 2 d\ln G_W - d\ln K - 4 d\ln\Omega_W
+
2E_*(d\ln R - d\ln\Omega_U - d\ln\Omega_W)
-
F_* d\ln\delta_U = 0.
\]

These solve triangularly:

### Tracking
\[
d\ln\delta_U
=
-\frac{1+\delta_U}{1+\chi_0}
\bigl(d\ln R + d\ln G_U - d\ln G_W - 2 d\ln\Omega_U\bigr).
\]

### Dressing
\[
d\ln M
=
 d\ln K - 2 d\ln G_U + 2 d\ln\Omega_U.
\]

### Nontracking
\[
d\ln\Omega_W
=
\frac{d\ln G_W - d\ln G_U + (1-E_*)d\ln\Omega_U + E_* d\ln R - \tfrac{F_*}{2}d\ln\delta_U}{E_*+2}.
\]

So once a candidate moving-throat branch gives the drifts of

- `K`
- `G_U`
- `G_W`
- `R`
- `Omega_U`

it automatically fixes the drifts required in

- `delta_U`
- `M`
- `Omega_W`

to stay tangent to the exact similarity orbit.

## Practical Stage-5 absolute-slope form

Writing the Stage-5 primitive slopes as

- `dK`
- `dM`
- `d(lambda_U)`
- `d(lambda_W)`
- `d(lambda_R)`
- `d(Omega_U)`
- `d(Omega_W)`
- `d(delta_U)`

one gets the exact compatibility formulas

\[
d(\delta_U)
=
-\delta_U\frac{1+\delta_U}{1+\chi_0}
\left[
\frac{d\lambda_R}{\lambda_R}
+
\frac{d\lambda_U}{\lambda_U}
-
\frac{d\lambda_W}{\lambda_W}
-
2\frac{d\Omega_U}{\Omega_U}
\right],
\]

\[
dM
=
M\left[
\frac{dK}{K} - 2\frac{d\lambda_U}{\lambda_U} + 2\frac{d\Omega_U}{\Omega_U}
\right],
\]

\[
\frac{d\Omega_W}{\Omega_W}
=
\frac{1}{E_*+2}
\left[
\frac{d\lambda_W}{\lambda_W}
-
\frac{d\lambda_U}{\lambda_U}
+
(1-E_*)\frac{d\Omega_U}{\Omega_U}
+
E_*\frac{d\lambda_R}{\lambda_R}
-
\frac{F_*}{2}\frac{d(\delta_U)}{\delta_U}
\right].
\]

So the Stage-10/11 similarity-orbit theorem is now directly usable in the Stage-5 primitive deformation language.

## Immediate consequence for the 5PN program

The next honest theorem gate is now smaller than before:

1. extract the actual branch drifts of
   `K, lambda_U, lambda_W, lambda_R, Omega_U`
   from the moving-throat PDE,
2. use the formulas above to predict the required co-drifts of
   `delta_U, M, Omega_W`
   if the branch is tangent to the exact zero-defect similarity orbit,
3. then compare those predictions with the actual reduced PDE branch.

If they agree, the weak-axisymmetric first-order obstruction is pure similarity-gauge and the calibration survives. If they fail, the moving-throat branch leaves the exact monomial-preserving orbit and the calibration breaks for a concrete microscopic reason.

## Stage 14 — the BdG primitive drifts are exactly support-blind for the Stage-11 monomials

If the primitive drift space is extended back to the full Stage-5 list

- `lambda_B`
- `varpi`
- `lambda_U`
- `lambda_W`
- `lambda_R`
- `K`
- `M`
- `Omega_U`
- `Omega_W`
- `delta_U`

then the direct Stage-11 monomial-drift matrix acquires **two exact zero columns** in the `dln lambda_B` and `dln varpi` directions.

So the weak-axisymmetric direct monomials are exactly support-blind:

\[
\partial_{\ln \lambda_B}
(\delta\ln \mathfrak C_{{\rm tr},*},
 \delta\ln \mathfrak C_{{\rm nt},*},
 \delta\ln \epsilon_\eta)
=0,
\]
\[
\partial_{\ln \varpi}
(\delta\ln \mathfrak C_{{\rm tr},*},
 \delta\ln \mathfrak C_{{\rm nt},*},
 \delta\ln \epsilon_\eta)
=0.
\]

That means the extended primitive zero-defect tangent space has dimension

\[
2 + 5 = 7,
\]

namely

1. two BdG-support-blind directions,
2. plus the five normalized similarity directions from Stage 13.

This is important conceptually.

The Stage-10/11 similarity-orbit theorem constrains only the mixed/wall/U normalization problem. It does **not** constrain the explicit BdG support drifts. So monomial-rigidity alone cannot finish the full conservative 5PN problem.

Those BdG directions must still be controlled, if at all, by the separate conservative front-end conditions:

- `K1 = 0`,
- the hidden-even consistency slot,
- and the `O(omega^4)` single-pole / grouped-response test.

So the Stage-5/6 conservative extraction theorem and the Stage-10/11 similarity-orbit theorem are complementary rather than redundant.
