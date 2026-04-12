# 5PN continuation notes — Stages 315–318

These stages finally turn the Stage 312–314 microscopic similarity-orbit theorem into a **usable actual-branch tester** in the reduced continuum variables that the moving-throat branch is expected to hand us.

The key reduction is that the actual coherent branch does **not** need the full eight-variable microscopic kernel state to test weak-axisymmetric orbit lock. It only needs the five-state packet

\[
(\chi_0,\ \delta_U,\ \widehat Z_W,\ \epsilon_W,\ \epsilon_\eta),
\]

where

\[
\widehat Z_W := \frac{\lambda_W^2 \mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
= \frac{Z_W}{\Omega_W^2}.
\]

So the actual-branch weak-axisymmetric theorem gate is now formulated directly in the language of the coherent continuum branch rather than in the older full microscopic kernel ledger.

## Stage 315 — exact continuum monomial drift map

The direct Stage-312 monomials become, in the reduced continuum variables,

\[
\mathfrak C_{{\rm tr},*}
=
\chi_0^{1+\delta_{U,*}}\,\delta_U^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\widehat Z_W\,\epsilon_W^{E_*}\,\delta_U^{-F_*},
\]

\[
\epsilon_\eta = \epsilon_\eta.
\]

Therefore the finite quotient packet is exactly

\[
q_{\rm tr}
=
(1+\delta_{U,*})\ln\!\frac{\chi_0}{\chi_{0,\rm ref}}
+
(1+\chi_{0,*})\ln\!\frac{\delta_U}{\delta_{U,\rm ref}},
\]

\[
q_{\rm nt}
=
\ln\!\frac{\widehat Z_W}{\widehat Z_{W,\rm ref}}
+
E_*\ln\!\frac{\epsilon_W}{\epsilon_{W,\rm ref}}
-
F_*\ln\!\frac{\delta_U}{\delta_{U,\rm ref}},
\]

\[
q_\eta
=
\ln\!\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
\]

At infinitesimal level, the branch-adapted monomial drifts are

\[
\Sigma_{\rm tr}
=
(1+\delta_{U,*})\,d\ln\chi_0
+
(1+\chi_{0,*})\,d\ln\delta_U,
\]

\[
\Sigma_{\rm nt}
=
 d\ln\widehat Z_W + E_*\,d\ln\epsilon_W - F_*\,d\ln\delta_U,
\]

\[
\Sigma_\eta = d\ln\epsilon_\eta.
\]

So the actual orbit-lock tester now acts on the reduced five-drift packet

\[
(d\ln\chi_0,\ d\ln\delta_U,\ d\ln\widehat Z_W,\ d\ln\epsilon_W,\ d\ln\epsilon_\eta).
\]

## Stage 316 — exact support-blindness of the continuum orbit tester

The coherent support-enhancement branch still has

\[
M_{\rm tr}=M_{\rm mix}\,S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\]

and the explicit twin-family support tower still determines the physical support ratio

\[
\zeta_n^{(\rm twin)}
=
\frac{1}{(2n+1)^2\bigl(1+x\,n(n+1)\bigr)}.
\]

But the actual continuum quotient packet and the actual defect packet are exactly blind to

- the support-enhancement ratio \(\zeta\),
- the total baseline \(M_{\rm tr}\),
- and the explicit twin-family stiffness \(x\).

Symbolically,

\[
\partial_\zeta(q_{\rm tr},q_{\rm nt},q_\eta)=0,
\qquad
\partial_{M_{\rm tr}}(q_{\rm tr},q_{\rm nt},q_\eta)=0,
\qquad
\partial_x(q_{\rm tr},q_{\rm nt},q_\eta)=0,
\]

and likewise for the physical defect triplet

\[
(\Theta_1,\Xi_1,\mathcal R_1).
\]

So support compensation and support-harmonic selection belong entirely to the **isotropic normalization** side of the 5PN endgame. They do **not** move the coherent weak-axisymmetric branch on or off the exact similarity orbit.

## Stage 317 — exact reduced actual-branch orbit tester

The five-drift packet maps linearly to the monomial-drift triple through the exact reduced projector

\[
M_{\rm red}:
(d\ln\chi_0,\ d\ln\delta_U,\ d\ln\widehat Z_W,\ d\ln\epsilon_W,\ d\ln\epsilon_\eta)
\mapsto
(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta).
\]

Then the physical defect triplet is

\[
\Theta_1 = -C_{\rm tr}\,\Sigma_{\rm tr},
\]

\[
\Xi_1 = A_{\rm tr}\,\Sigma_{\rm tr}+\Sigma_{\rm nt},
\]

\[
\mathcal R_1 = -\Xi_1 - \frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\Sigma_\eta,
\]

with

\[
C_{\rm tr}=
\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})},
\qquad
A_{\rm tr}=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}.
\]

The reduced tester has exact rank `3`, so its zero-defect kernel has dimension `2` inside the five-drift continuum space.

A convenient exact basis is

\[
v_\chi
=
\left(
1,
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}},
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}},
0,
0
\right),
\]

\[
v_{\epsilon_W}
=
(0,0,-E_*,1,0).
\]

Equivalently, exact orbit lock on the actual coherent branch is the codrift system

\[
d\ln\delta_U
=
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}\,d\ln\chi_0,
\]

\[
d\ln\widehat Z_W
=
- E_*\,d\ln\epsilon_W + F_*\,d\ln\delta_U,
\]

\[
d\ln\epsilon_\eta = 0.
\]

So the actual coherent branch lies on the exact similarity orbit iff its reduced drift vector lies in the 2-plane spanned by \(v_\chi\) and \(v_{\epsilon_W}\).

## Stage 318 — packaged actual-branch tester API

The new tester is now ready in two equivalent forms.

### Finite tester
Input:

\[
(\chi_0,\delta_U,\widehat Z_W,\epsilon_W,\epsilon_\eta)
\]

plus a reference branch state and the frozen reference exponents.

Output:

\[
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

### Infinitesimal tester
Input:

\[
(d\ln\chi_0,\ d\ln\delta_U,\ d\ln\widehat Z_W,\ d\ln\epsilon_W,\ d\ln\epsilon_\eta).
\]

Output:

either

\[
(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)
\]

or

\[
(\Theta_1,\Xi_1,\mathcal R_1).
\]

The exact zero-defect packet is the 2-plane above, so the physical branch no longer needs to be checked against a larger microscopic state space.

## Net result after Stages 315–318

These stages close the exact step we had been missing after Stage 312–314:

1. the microscopic similarity-orbit theorem is now expressed directly in the reduced continuum variables of the actual coherent branch;
2. support enhancement and twin/non-twin support selection are proven to be orthogonal to the weak-axisymmetric orbit-lock problem;
3. the actual branch tester is now a concrete five-variable input / three-variable output object;
4. and the exact zero-defect branch inside that five-variable continuum state space is a 2-plane.

So the next theorem gate is no longer “how do we test the real branch?”
It is narrower:

> extract the actual moving-throat branch drifts of
> \(\chi_0,\delta_U,\widehat Z_W,\epsilon_W,\epsilon_\eta\),
> feed them into the new tester,
> and check whether the branch lies in the exact reduced orbit plane while simultaneously satisfying the isotropic normalization branch conditions.
