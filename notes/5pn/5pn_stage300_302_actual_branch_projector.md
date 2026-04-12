# 5PN continuation notes — Stages 300–302

These stages continue directly from Stages 296–299.

The earlier reduction already did two key things:

1. Stage 296–298 showed that the physical first-order grouped weak-axisymmetric defect is carried by the three-scalar normal form
   
   - `Theta_1`,
   - `Xi_1 = P1/P0`,
   - `R_1`,
   
   with exact microscopic coordinates `(Sigma_tr, Sigma_nt, Sigma_eta)`.

2. Stage 299 showed that these are exactly the first logarithmic drifts of the three direct microscopic monomials
   
   - `C_tr,*`,
   - `C_nt,*`,
   - `epsilon_eta`,
   
   and that the zero-defect branch is the tangent space of the exact five-parameter monomial-preserving similarity orbit.

What was still missing was the **actual branch tester**:

> given an actual moving-throat microscopic branch state or branch tangent, how do we decide whether it is on-orbit or off-orbit without re-solving the whole drift system every time?

Stages 300–302 supply that missing tester.

---

## Stage 300 — exact finite quotient coordinates for the actual branch

Instead of working only infinitesimally, introduce the exact finite quotient coordinates between an actual microscopic branch state and a reference branch state:

\[
Q_{\rm tr}=\ln\!\frac{\mathfrak C_{{\rm tr},*}}{\mathfrak C_{{\rm tr},*,\rm ref}},
\qquad
Q_{\rm nt}=\ln\!\frac{\mathfrak C_{{\rm nt},*}}{\mathfrak C_{{\rm nt},*,\rm ref}},
\qquad
Q_{\eta}=\ln\!\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
\]

With the direct microscopic monomials from Stage 299,

\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2 T_U}{L^2 K_U}\right)^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{\rm (eff)}(K_W^{\rm (eff)})^2}
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{\rm (eff)}}\right)^{E_*}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{-F_*},
\]

\[
\epsilon_\eta=
\frac{c_{\eta U}^2}{K_U K_\eta^{\rm (eff)}},
\]

these quotient coordinates are exact on the positive microscopic state space.

The key Stage-300 theorem is that their first linearization reproduces the Stage-299 monomial-drift vector exactly:

\[
\delta Q = M_*\,\delta x,
\]

where `M_*` is the exact rank-3 monomial-drift matrix from Stage 299.

So the finite quotient language and the infinitesimal defect language are now stitched together exactly.

---

## Stage 301 — exact affine decomposition of an actual branch tangent

Let the microscopic drift vector be ordered as

\[
\delta x =
(d\ln\lambda_W,
 d\ln c_{\eta U},
 d\ln\gamma,
 d\ln K_U,
 d\ln K_\eta^{\rm (eff)},
 d\ln K_W^{\rm (eff)},
 d\ln\mu_W,
 d\ln T_U).
\]

The Stage-299 compatibility equations can be solved exactly for the three dependent drifts

- `d ln K_eta^(eff)`,
- `d ln T_U`,
- `d ln mu_W`,

in terms of the five free drifts

- `d ln lambda_W`,
- `d ln c_etaU`,
- `d ln gamma`,
- `d ln K_U`,
- `d ln K_W^(eff)`,

plus the three quotient residuals `(q_tr, q_nt, q_eta)`.

The exact affine solve is

\[
d\ln K_\eta^{\rm (eff)} = 2\,d\ln c_{\eta U} - d\ln K_U - q_\eta,
\]

\[
d\ln T_U = d\ln K_U
+\frac{q_{\rm tr}-(1+\delta_{U,*})(d\ln\gamma+d\ln c_{\eta U}-d\ln K_U)}{1+\chi_{0,*}},
\]

\[
d\ln\mu_W
=
q_{\rm nt}-q_\eta - d\ln K_U + 2d\ln K_W^{\rm (eff)} + 2d\ln c_{\eta U} - 2d\ln\lambda_W
\]
\[
\qquad
+E_*\bigl(d\ln K_U + d\ln K_W^{\rm (eff)} - 2d\ln\gamma - 2d\ln\lambda_W\bigr)
\]
\[
\qquad
+\frac{F_*\bigl((1+\delta_{U,*})(d\ln K_U-d\ln c_{\eta U}-d\ln\gamma)+q_{\rm tr}\bigr)}{1+\chi_{0,*}}.
\]

This produces an exact decomposition

\[
\delta x_{\rm actual} = \delta x_{\rm orbit} + N\,q,
\qquad
q=(q_{\rm tr},q_{\rm nt},q_\eta)^T,
\]

where

- `delta x_orbit` lies in the five-dimensional similarity-orbit tangent space,
- `N` is a convenient exact `8 x 3` normal-basis matrix satisfying
  \[
  M_* N = I_3.
  \]

A particularly simple exact normal basis is

\[
n_{\rm tr} = \left(0,0,0,0,0,0,\frac{F_*}{1+\chi_{0,*}},\frac{1}{1+\chi_{0,*}}\right)^T,
\]

\[
n_{\rm nt} = (0,0,0,0,0,0,1,0)^T,
\]

\[
n_{\eta} = (0,0,0,0,-1,0,-1,0)^T.
\]

So any actual branch tangent now has a unique exact split:

- five tangent directions that preserve the three monomials,
- plus three normal coordinates `(q_tr, q_nt, q_eta)` measuring the off-orbit failure.

---

## Stage 302 — exact first-order defect compiler from quotient residuals

The three quotient residuals are not just convenient coordinates. They map exactly into the physical first-order defect triplet.

Using the reference-branch coefficients

\[
C_{\rm tr}
=
\frac{\chi_{0,*}\delta_{U,*}}{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})},
\qquad
A_{\rm tr}
=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})},
\]

and `epsilon_eta,*`, the exact first-order map is

\[
\Theta_1 = - C_{\rm tr}\, q_{\rm tr},
\]

\[
\Xi_1 = A_{\rm tr}\, q_{\rm tr} + q_{\rm nt},
\qquad
\Xi_1 = \frac{P_1}{P_0},
\]

\[
\mathcal R_1 = -\Xi_1 - \frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\, q_\eta.
\]

So in matrix form,

\[
\begin{pmatrix}
\Theta_1\\[3pt]
\Xi_1\\[3pt]
\mathcal R_1
\end{pmatrix}
=
\begin{pmatrix}
-C_{\rm tr} & 0 & 0\\
A_{\rm tr} & 1 & 0\\
-A_{\rm tr} & -1 & -\dfrac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}
\end{pmatrix}
\begin{pmatrix}
q_{\rm tr}\\[3pt]
q_{\rm nt}\\[3pt]
q_\eta
\end{pmatrix}.
\]

The determinant of this defect matrix is nonzero on the physical reference branch,
so the map is exactly invertible. The inverse is

\[
q_{\rm tr} = -\frac{\Theta_1}{C_{\rm tr}},
\qquad
q_{\rm nt} = \Xi_1 + \frac{A_{\rm tr}}{C_{\rm tr}}\Theta_1,
\qquad
q_\eta = -\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}(\mathcal R_1+\Xi_1).
\]

So the actual moving-throat branch can now be tested in either of two exactly equivalent first-order languages:

1. finite quotient residuals `(q_tr, q_nt, q_eta)`,
2. physical defect triplet `(Theta_1, Xi_1, R_1)`.

And the zero-defect theorem is now the cleanest it has been:

\[
\Theta_1 = \Xi_1 = \mathcal R_1 = 0
\iff
q_{\rm tr}=q_{\rm nt}=q_\eta=0.
\]

---

## Net result after Stage 302

The continuation point is now smaller again.

We no longer need to test an actual moving-throat weak-axisymmetric branch against the raw eight-dimensional drift space.

The branch only needs to provide either:

- the three finite quotient residuals
  \[
  (Q_{\rm tr},Q_{\rm nt},Q_\eta),
  \]
  whose first linearization is `(q_tr,q_nt,q_eta)`,

or equivalently

- the three physical first-order defects
  \[
  (\Theta_1,\Xi_1,\mathcal R_1).
  \]

So the next honest theorem gate is no longer algebraic compression. It is finally a real branch test:

> compute the actual moving-throat branch state (or branch tangent), form its three exact quotient residuals against the reference coherent branch, and see whether those residuals vanish.

That is the cleanest actual-branch projector we have so far.
