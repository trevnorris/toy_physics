# 5PN Stages 187–192 — Finite Similarity-Orbit Action, Reference-Transport Laws, and Exact Residual Diagnostics

This block continues directly from the finite orbit interface at Stage 186.
The earlier stages had already shown two things:

1. the weak-axisymmetric zero-defect branch is exactly the finite similarity orbit \
   \(\mathcal G_*\), and
2. the full finite branch-selection problem can be written as a test on the dependent
   triple \((T_U,K_\eta,\mu_W)\).

What was still missing was the **finite transport law** itself: how the five free microscopic
coordinates move a branch along an exact orbit, how to compare an actual candidate branch to a
reference orbit point, and how to localize any failure into exact multiplicative residuals.

## Stage 187 — exact finite similarity-orbit action

The five-parameter multiplicative similarity orbit \(\mathcal G_*\) is now written explicitly as a
finite action on the full microscopic state. If
\[
(\lambda_W,c_{\eta U},\gamma,K_U,K_W)
\to
(e^{\Lambda}\lambda_W,e^C c_{\eta U},e^{\Gamma}\gamma,e^U K_U,e^W K_W),
\]
then the dependent triple coevolves by
\[
K_\eta' = e^{2C-U}K_\eta,
\]
\[
T_U' = e^{U-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)}T_U,
\]
\[
\mu_W'
=
\exp\!\Bigl[
2C-U+2W-2\Lambda
-E_*\bigl(2\Gamma+2\Lambda-U-W\bigr)
-F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)
\Bigr]\mu_W.
\]
The three quotient monomials
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)
\]
are preserved **exactly**, not only infinitesimally.

## Stage 188 — exact group law, inverse, and parameter recovery

The finite orbit action is an exact abelian five-parameter group:
composition adds the five generator exponents, and inversion negates them. If two states are on
one orbit, the orbit parameters are recovered exactly from the free-coordinate ratios,
\[
\Lambda=\ln\frac{\lambda_W'}{\lambda_W},
\qquad
C=\ln\frac{c_{\eta U}'}{c_{\eta U}},
\qquad
\Gamma=\ln\frac{\gamma'}{\gamma},
\qquad
U=\ln\frac{K_U'}{K_U},
\qquad
W=\ln\frac{K_W'}{K_W}.
\]
So once those five free-coordinate ratios are known, the dependent triple on the same orbit is
predicted uniquely.

## Stage 189 — exact reference-orbit transport laws

Relative to a reference orbit point, the exact dependent-coordinate transport laws are
\[
R_{K_\eta}^{(\mathrm{orbit})}=\frac{R_c^2}{R_U},
\]
\[
R_{T_U}^{(\mathrm{orbit})}=R_U\left(\frac{R_U}{R_\gamma R_c}\right)^{\frac{1+\delta_{U,*}}{1+\chi_{0,*}}},
\]
\[
R_{\mu_W}^{(\mathrm{orbit})}
=
\frac{R_{K_\eta}^{(\mathrm{orbit})}R_W^2}{R_\lambda^2}
\left(\frac{R_\gamma^2R_\lambda^2}{R_UR_W}\right)^{-E_*}
\left(\frac{R_{T_U}^{(\mathrm{orbit})}}{R_U}\right)^{F_*}.
\]
This is the exact finite coevolution law of the dependent triple along a fixed orbit.

## Stage 190 — exact dependent residual coordinates

A general candidate branch with the same five free-coordinate ratios can be factored as
\[
R_{T_U}^{(\mathrm{actual})}=R_{T_U}^{(\mathrm{orbit})}m_T,
\qquad
R_{K_\eta}^{(\mathrm{actual})}=R_{K_\eta}^{(\mathrm{orbit})}m_K,
\qquad
R_{\mu_W}^{(\mathrm{actual})}=R_{\mu_W}^{(\mathrm{orbit})}m_\mu,
\]
where \((m_T,m_K,m_\mu)\) is the **dependent residual mismatch triple**.
The invariant ratios then collapse exactly to
\[
\frac{\mathfrak C_{{\rm tr},*}^{\mathrm{actual}}}{\mathfrak C_{{\rm tr},*}^{\mathrm{ref}}}=m_T^{1+\chi_{0,*}},
\qquad
\frac{\epsilon_\eta^{\mathrm{actual}}}{\epsilon_{\eta}^{\mathrm{ref}}}=\frac{1}{m_K},
\qquad
\frac{\mathfrak C_{{\rm nt},*}^{\mathrm{actual}}}{\mathfrak C_{{\rm nt},*}^{\mathrm{ref}}}=\frac{m_\mu}{m_Km_T^{F_*}}.
\]
So the free-coordinate transport along an orbit drops out completely; the quotient coordinates are
nothing but the logarithmic chart of the dependent residual triple.

## Stage 191 — factorized actual-branch interface

The actual candidate branch now admits an exact factorization
\[
\text{actual branch}
=
(\text{reference orbit point})
\times
(\text{free-coordinate orbit transport})
\times
(\text{dependent residual mismatch}).
\]
Restoration to the same orbit at fixed free-coordinate ratios is achieved by dividing only by the
residual mismatch ratios:
\[
K_\eta^{(\mathrm{restore})}=\frac{K_\eta^{(\mathrm{actual})}}{m_K},
\qquad
T_U^{(\mathrm{restore})}=\frac{T_U^{(\mathrm{actual})}}{m_T},
\qquad
\mu_W^{(\mathrm{restore})}=\frac{\mu_W^{(\mathrm{actual})}}{m_\mu}.
\]
So orbit lock is exactly the statement
\[
m_T=m_K=m_\mu=1.
\]

## Stage 192 — diagnostic signatures of each failure channel

The three quotient coordinates now have a direct physical interpretation:
\[
q_{\rm tr}=(1+\chi_{0,*})\ln m_T,
\qquad
q_\eta=-\ln m_K,
\qquad
q_{\rm nt}=\ln m_\mu-\ln m_K-F_*\ln m_T.
\]
That gives three clean signatures:

- pure \(T_U\) residual mismatch turns on \(q_{\rm tr}\) and, via \(F_*\), also \(q_{\rm nt}\),
- pure \(K_\eta\) residual mismatch turns on \(q_\eta\) and also \(q_{\rm nt}\),
- pure \(\mu_W\) residual mismatch turns on \(q_{\rm nt}\) only.

So once the PDE delivers an actual candidate branch, the pattern of
\[
(q_{\rm tr},q_\eta,q_{\rm nt})
\]
identifies exactly which dependent coevolution law failed, if any.

## Bottom line after Stage 192

The branch-selection problem is no longer merely
“compute the actual branch and see if the quotient coordinates move.”
It is now a **factorized finite comparison problem**:

1. choose a reference point on the exact similarity orbit,
2. use the five free-coordinate ratios to predict the orbit-transported dependent triple,
3. compare the actual dependent triple to that prediction,
4. read off the residual mismatch ratios \((m_T,m_K,m_\mu)\),
5. and infer the quotient coordinates and restoration map immediately.

So the next PDE theorem gate is now even sharper:

> compute the actual branch values of the eight microscopic coordinates, form the free-coordinate
> ratios and the dependent residual mismatch triple, and test whether
> \(m_T=m_K=m_\mu=1\).

If yes, the branch stays on a single exact \(\mathcal G_*\) orbit. If not, the failure is localized
exactly into the \(T_U\), \(K_\eta\), and/or \(\mu_W\) transport laws.
