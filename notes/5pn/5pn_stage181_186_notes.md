# Stage 181–186 notes

This block continues the post-170 / post-180 branch-selection program in the most natural direction: it upgrades the **first-order** orbit-lock language to an **exact finite** law for the dependent microscopic triple.

## Main new result

The exact Stage-168 monomials

\[
\mathfrak C_{{\rm tr},*},\qquad
\mathfrak C_{{\rm nt},*},\qquad
\epsilon_\eta
\]

can be solved **exactly** for the dependent microscopic coordinates

\[
(T_U, K_\eta, \mu_W)
\]

once the five free microscopic coordinates

\[
(\lambda_W, c_{\eta U}, \gamma, K_U, K_W)
\]

and the invariant triple are fixed.

That gives the exact finite single-orbit law:

\[
K_\eta^{(\mathrm{orbit})} = \frac{c_{\eta U}^2}{K_U\,\epsilon_{\eta,*}},
\]

\[
T_U^{(\mathrm{orbit})}
= \frac{L^2 K_U}{\pi^2}
\left[
\frac{\mathfrak C_{{\rm tr},*}}
{(\gamma c_{\eta U}/K_U)^{1+\delta_{U,*}}}
\right]^{\!1/(1+\chi_{0,*})},
\]

\[
\mu_W^{(\mathrm{orbit})}
=
\frac{\mathfrak C_{{\rm nt},*} K_\eta^{(\mathrm{orbit})} K_W^2}{\lambda_W^2}
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W}\right)^{-E_*}
\left(\frac{\pi^2 T_U^{(\mathrm{orbit})}}{L^2 K_U}\right)^{F_*}.
\]

So the finite similarity orbit is no longer abstract: the dependent triple is an exact algebraic function of the free microscopic point and the invariant triple.

## Exact finite mismatch coordinates

For any candidate branch with the **same five free microscopic coordinates** as the orbit point, define

\[
T_U = m_T\,T_U^{(\mathrm{orbit})},\qquad
K_\eta = m_K\,K_\eta^{(\mathrm{orbit})},\qquad
\mu_W = m_\mu\,\mu_W^{(\mathrm{orbit})}.
\]

Then the exact invariant ratios are

\[
\frac{\mathfrak C_{{\rm tr},*}}{\mathfrak C_{{\rm tr},*}^{(\mathrm{orbit})}} = m_T^{1+\chi_{0,*}},
\qquad
\frac{\epsilon_\eta}{\epsilon_{\eta,*}} = \frac{1}{m_K},
\qquad
\frac{\mathfrak C_{{\rm nt},*}}{\mathfrak C_{{\rm nt},*}^{(\mathrm{orbit})}} = \frac{m_\mu}{m_K m_T^{F_*}}.
\]

So the finite branch-selection problem is exactly three-dimensional.

## Exact logarithmic chart

If we write

\[
\tau := \ln m_T,
\qquad
\kappa := \ln m_K,
\qquad
\mu := \ln m_\mu,
\]

then the quotient coordinates are **exactly**

\[
q_{\rm tr} = (1+\chi_{0,*})\tau,
\qquad
q_\eta = -\kappa,
\qquad
q_{\rm nt} = \mu - \kappa - F_*\tau.
\]

So the Stage-179 first-order formulas are not merely infinitesimal approximations; they are the exact logarithmic chart of the finite mismatch ratios.

## Exact restoration map

Given the finite quotient coordinates, the exact restoration to the same orbit is achieved by changing only the dependent triple:

\[
T_U^{(\mathrm{restore})} = T_U\,e^{-q_{\rm tr}/(1+\chi_{0,*})},
\]

\[
K_\eta^{(\mathrm{restore})} = K_\eta\,e^{q_\eta},
\]

\[
\mu_W^{(\mathrm{restore})}
= \mu_W\,e^{-q_{\rm nt}+q_\eta-F_* q_{\rm tr}/(1+\chi_{0,*})}.
\]

This returns the candidate branch to the exact orbit with the same free microscopic coordinates and the same invariant triple.

## Finite adiabatic-elastic orbit-lock criterion

Under the adiabatic-elastic boundary rule, the exact finite branch-selection problem is:

\[
\text{orbit lock}
\iff
m_T = m_K = m_\mu = 1
\iff
q_{\rm tr}=q_{\rm nt}=q_\eta=0.
\]

So after this block the remaining theorem gap is completely concrete: once the PDE gives the actual microscopic branch values, one can test orbit lock directly by comparing the candidate dependent triple to the exact orbit values above.
