
# Stages 193–198 — Pairwise orbit-lock, cocycle laws, and minimal-data verdict

This block continues the finite orbit-lock chain after Stages 181–192 by making the
orbit criterion fully **reference-independent** and **PDE-ready**.

## Stage 193 — exact pairwise orbit criterion

Given two positive microscopic states `x` and `y`, with shared branch constants
`(chi0_*, deltaU_*, E_*, F_*)`, the five free-coordinate ratios
\[
R_\lambda,\ R_c,\ R_\gamma,\ R_U,\ R_W
\]
determine the exact orbit-predicted dependent-coordinate ratios
\[
R_{K_\eta}^{(\mathrm{orbit})}=\frac{R_c^2}{R_U},
\]
\[
R_{T_U}^{(\mathrm{orbit})}
=
R_U\left(\frac{R_U}{R_\gamma R_c}\right)^{\frac{1+\delta_{U,*}}{1+\chi_{0,*}}},
\]
\[
R_{\mu_W}^{(\mathrm{orbit})}
=
\frac{R_{K_\eta}^{(\mathrm{orbit})}R_W^2}{R_\lambda^2}
\left(\frac{R_\gamma^2R_\lambda^2}{R_UR_W}\right)^{-E_*}
\left(\frac{R_{T_U}^{(\mathrm{orbit})}}{R_U}\right)^{F_*}.
\]

Comparing with the actual dependent-coordinate ratios defines the pairwise residual triple
\[
m_T^{(x\to y)}=\frac{R_{T_U}^{(\mathrm{act})}}{R_{T_U}^{(\mathrm{orbit})}},
\qquad
m_K^{(x\to y)}=\frac{R_{K_\eta}^{(\mathrm{act})}}{R_{K_\eta}^{(\mathrm{orbit})}},
\qquad
m_\mu^{(x\to y)}=\frac{R_{\mu_W}^{(\mathrm{act})}}{R_{\mu_W}^{(\mathrm{orbit})}}.
\]

The invariant ratios between `y` and `x` are then exactly
\[
\frac{\mathfrak C_{{\rm tr},*}(y)}{\mathfrak C_{{\rm tr},*}(x)}
=
\bigl(m_T^{(x\to y)}\bigr)^{1+\chi_{0,*}},
\]
\[
\frac{\epsilon_\eta(y)}{\epsilon_\eta(x)}
=
\frac{1}{m_K^{(x\to y)}},
\]
\[
\frac{\mathfrak C_{{\rm nt},*}(y)}{\mathfrak C_{{\rm nt},*}(x)}
=
\frac{m_\mu^{(x\to y)}}{m_K^{(x\to y)}\bigl(m_T^{(x\to y)}\bigr)^{F_*}}.
\]

So `x` and `y` lie on the same exact similarity orbit iff
\[
m_T^{(x\to y)}=m_K^{(x\to y)}=m_\mu^{(x\to y)}=1,
\]
equivalently iff the three invariant ratios are all unity.

## Stage 194 — multiplicative cocycle and additive quotient law

For three states `x,y,z`, the residual ratios compose exactly:
\[
m_T^{(x\to z)}=m_T^{(x\to y)}m_T^{(y\to z)},
\qquad
m_K^{(x\to z)}=m_K^{(x\to y)}m_K^{(y\to z)},
\qquad
m_\mu^{(x\to z)}=m_\mu^{(x\to y)}m_\mu^{(y\to z)}.
\]

The logarithmic quotient coordinates therefore add:
\[
q_{\rm tr}^{(x\to z)}=q_{\rm tr}^{(x\to y)}+q_{\rm tr}^{(y\to z)},
\]
\[
q_{\rm nt}^{(x\to z)}=q_{\rm nt}^{(x\to y)}+q_{\rm nt}^{(y\to z)},
\]
\[
q_\eta^{(x\to z)}=q_\eta^{(x\to y)}+q_\eta^{(y\to z)},
\]
with inverse laws under reversal,
\[
q^{(y\to x)}=-q^{(x\to y)}.
\]

So a sequence of PDE branch snapshots can be tracked either multiplicatively in the
residual ratios or additively in the quotient coordinates.

## Stage 195 — pairwise quotient-to-observable compiler

From the residual ratios,
\[
q_{\rm tr}=(1+\chi_{0,*})\ln m_T,\qquad
q_\eta=-\ln m_K,\qquad
q_{\rm nt}=\ln m_\mu-\ln m_K-F_*\ln m_T.
\]

Composing with the Stage-170 linear observable map gives the small pairwise observable
signature
\[
\Theta_1^{(\mathrm{lin})}
=
-\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
\,q_{\rm tr},
\]
\[
\Xi_1^{(\mathrm{lin})}
=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}\,q_{\rm tr}+q_{\rm nt},
\]
\[
(\mathcal R_1+\Xi_1)^{(\mathrm{lin})}
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,q_\eta.
\]

The pure-channel signatures become:

- pure `T_U` mismatch:
  \[
  q_{\rm tr}=(1+\chi_{0,*})\ln z_T,\qquad
  q_{\rm nt}=-F_*\ln z_T,\qquad
  q_\eta=0;
  \]

- pure `K_\eta` mismatch:
  \[
  q_{\rm tr}=0,\qquad
  q_{\rm nt}=-\ln z_K,\qquad
  q_\eta=-\ln z_K;
  \]

- pure `\mu_W` mismatch:
  \[
  q_{\rm tr}=0,\qquad
  q_{\rm nt}=\ln z_M,\qquad
  q_\eta=0.
  \]

This is the finite pairwise uplift of the earlier Stage-192 channel diagnostics.

## Stage 196 — exact inversion from invariant ratios

The invariant-ratio packet
\[
(R_{\rm tr},R_{\rm nt},R_\eta)
:=
\left(
\frac{\mathfrak C_{{\rm tr},*}(y)}{\mathfrak C_{{\rm tr},*}(x)},
\frac{\mathfrak C_{{\rm nt},*}(y)}{\mathfrak C_{{\rm nt},*}(x)},
\frac{\epsilon_\eta(y)}{\epsilon_\eta(x)}
\right)
\]
already contains the full orbit-lock verdict.

The exact inversion is
\[
m_T = R_{\rm tr}^{1/(1+\chi_{0,*})},
\qquad
m_K = \frac{1}{R_\eta},
\qquad
m_\mu = R_{\rm nt}\,m_K\,m_T^{F_*}.
\]

So once the invariant ratios are known, the five free-coordinate ratios are not needed
for the verdict or for the dependent-coordinate restoration.

## Stage 197 — canonical orbit-distance quadratic form

Write the logarithmic residuals as
\[
t=\ln m_T,\qquad k=\ln m_K,\qquad \mu=\ln m_\mu.
\]

Then
\[
\begin{pmatrix}
q_{\rm tr}\\ q_{\rm nt}\\ q_\eta
\end{pmatrix}
=
A
\begin{pmatrix}
t\\k\\\mu
\end{pmatrix},
\qquad
A=
\begin{pmatrix}
1+\chi_{0,*} & 0 & 0\\
-F_* & -1 & 1\\
0 & -1 & 0
\end{pmatrix}.
\]

The canonical scalar orbit-distance is
\[
D^2=q_{\rm tr}^2+q_{\rm nt}^2+q_\eta^2
=
\begin{pmatrix}
t&k&\mu
\end{pmatrix}
Q
\begin{pmatrix}
t\\k\\\mu
\end{pmatrix},
\qquad
Q=A^TA,
\]
i.e.
\[
Q=
\begin{pmatrix}
(1+\chi_{0,*})^2+F_*^2 & F_* & -F_*\\
F_* & 2 & -1\\
-F_* & -1 & 1
\end{pmatrix}.
\]

Its principal minors are
\[
(1+\chi_{0,*})^2+F_*^2>0,
\]
\[
2(1+\chi_{0,*})^2+F_*^2>0,
\]
\[
(1+\chi_{0,*})^2>0,
\]
so `Q` is positive definite on the constructive branch.

Therefore
\[
D^2=0
\iff
m_T=m_K=m_\mu=1.
\]

So the entire pairwise orbit-lock failure can be summarized by one exact reference-free
positive scalar.

## Stage 198 — minimal-data orbit verdict

The full orbit-lock verdict can be reached from **any one** of three exact packets:

1. residual mismatch ratios:
   \[
   (m_T,m_K,m_\mu),
   \]

2. invariant ratios:
   \[
   (R_{\rm tr},R_{\rm nt},R_\eta),
   \]

3. quotient coordinates:
   \[
   (q_{\rm tr},q_{\rm nt},q_\eta).
   \]

They are exactly interconvertible:
\[
R_{\rm tr}=m_T^{1+\chi_{0,*}},\qquad
R_{\rm nt}=\frac{m_\mu}{m_K m_T^{F_*}},\qquad
R_\eta=\frac{1}{m_K},
\]
\[
q_{\rm tr}=\ln R_{\rm tr},\qquad
q_{\rm nt}=\ln R_{\rm nt},\qquad
q_\eta=\ln R_\eta,
\]
\[
m_T=e^{q_{\rm tr}/(1+\chi_{0,*})},\qquad
m_K=e^{-q_\eta},\qquad
m_\mu=e^{q_{\rm nt}-q_\eta+F_*q_{\rm tr}/(1+\chi_{0,*})}.
\]

So future PDE numerics only need to provide whichever packet is cleanest. From that
packet one can reconstruct the dependent-coordinate restoration map and the scalar
distance `D^2`.
