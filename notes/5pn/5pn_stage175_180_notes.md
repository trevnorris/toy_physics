# 5PN continuation — Stages 175–180

This block continues directly from the adiabatic-wall / adiabatic-elastic orbit theorem.
The earlier result had already reduced the live first-order branch-selection problem to the exact
Stage-170 quotient test
\[
M_*\,\delta\mathbf x = 0,
\]
with quotient coordinates
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
\]

The missing step was to turn that theorem into an **exact branch-selection compiler**:
for any candidate microscopic drift, isolate the pure-similarity part, isolate the true quotient
failure, and write the unique microscopic correction that would restore single-orbit lock.

## Stage 175 — exact orbit/quotient projectors

Using the exact Stage-170 pivot block in the dependent coordinates
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu),
\]
with
\[
P_{(T,K_\eta,\mu)} = M_*[:,(T,K_\eta,\mu)],
qquad
\det P = 1+\chi_{0,*}>0,
\]
one gets exact complementary projectors
\[
Q = E P^{-1} M_*,
\qquad
O = I-Q.
\]
They satisfy
\[
Q^2=Q,
\qquad
O^2=O,
\qquad
QO=OQ=0,
\qquad
M_*O=0,
\qquad
M_*Q=M_*.
\]
So every microscopic drift splits uniquely as
\[
\Delta\mathbf x = \Delta\mathbf x_{\rm orbit}+\Delta\mathbf x_{\rm fail},
\]
with
\[
\Delta\mathbf x_{\rm orbit}\in\ker M_*,
\qquad
\Delta\mathbf x_{\rm fail}=Q\Delta\mathbf x.
\]
A particularly sharp result is that the entire quotient-failure piece has support only in the
**dependent triple**
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu),
\]
not in the five free similarity directions.

## Stage 176 — observable-to-microscopic correction compiler

The Stage-166/167 observable inversion gives
\[
\delta\ln\mathfrak C_{{\rm tr},*}
= -\frac{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}{\chi_{0,*}\delta_{U,*}}\,\Theta_1,
\]
\[
\delta\ln\mathfrak C_{{\rm nt},*}
= \Xi_1+\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}\,\Theta_1,
\]
\[
\delta\ln\epsilon_\eta
= -\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\,(\mathcal R_1+\Xi_1).
\]
Composing this with the Stage-175 projector gives the exact microscopic quotient correction
supported only on the dependent triple:
\[
\Delta_T^{(q)}=
\frac{\delta\ln\mathfrak C_{{\rm tr},*}}{1+\chi_{0,*}},
\]
\[
\Delta_{K_\eta}^{(q)}=-\delta\ln\epsilon_\eta,
\]
\[
\Delta_\mu^{(q)}
=
\delta\ln\mathfrak C_{{\rm nt},*}
-\delta\ln\epsilon_\eta
+\frac{F_*}{1+\chi_{0,*}}\,\delta\ln\mathfrak C_{{\rm tr},*}.
\]
So once the observable defect triple is known, the exact microscopic dependent-coordinate
correction needed to represent it is already fixed.

## Stage 177 — exact finite restoration operator

Because the Stage-170 fibre equations are linear in the finite log-ratios, the same projector
logic gives an exact **finite** restoration operator.
For any candidate finite drift \(\Delta\mathbf x\), define
\[
\mathbf q = M_*\Delta\mathbf x,
\qquad
\Delta\mathbf x_{\rm fail}=E P^{-1}\mathbf q,
\qquad
\Delta\mathbf x_{\rm restore}=-\Delta\mathbf x_{\rm fail}.
\]
Then
\[
M_*(\Delta\mathbf x+\Delta\mathbf x_{\rm restore})=0.
\]
So any candidate branch can be returned to a single exact \(\mathcal G_*\)-orbit by adjusting
only the dependent triple \((T_U,K_\eta,\mu_W)\).

## Stage 178 — adiabatic-elastic branch decomposition

Under the strengthened boundary rule
\[
\delta\ln\Theta_w=0,
\qquad
\varepsilon_L=\varepsilon_v=\varepsilon_T=0,
\]
the thermal wall channel and the scalar off-bundle escape are both removed. The remaining
first-order branch-selection problem is therefore purely microscopic and eight-dimensional.
For any adiabatic-elastic candidate branch drift \(\Delta\mathbf x_{AE}\), the exact split is
\[
\Delta\mathbf x_{AE} = \Delta\mathbf x_{\rm orbit}+\Delta\mathbf x_{\rm fail},
\]
with
\[
M_*\Delta\mathbf x_{\rm orbit}=0,
\qquad
M_*\Delta\mathbf x_{\rm fail}=M_*\Delta\mathbf x_{AE}.
\]
And the weak-axisymmetric observables depend only on the quotient piece:
\[
\Theta_1
= -\frac{\chi_{0,*}\delta_{U,*}}{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
\,q_1,
\]
\[
\Xi_1
= \frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}\,q_1 + q_2,
\]
\[
\mathcal R_1+\Xi_1
= -\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,q_3,
\]
where \(\mathbf q=M_*\Delta\mathbf x_{AE}\).
So the adiabatic-elastic branch is orbit-locked iff
\[
\Delta\mathbf x_{\rm fail}=0
\iff
M_*\Delta\mathbf x_{AE}=0.
\]

## Stage 179 — exact dependent-coordinate mismatch formulas

Comparing an arbitrary candidate branch to the exact Stage-170 single-orbit law gives the three
microscopic mismatch coordinates
\[
m_T := \Delta_T-\Delta_T^{\rm orbit},
\qquad
m_{K_\eta}:=\Delta_{K_\eta}-\Delta_{K_\eta}^{\rm orbit},
\qquad
m_\mu:=\Delta_\mu-\Delta_\mu^{\rm orbit}.
\]
These are not arbitrary. They are exactly the quotient drifts:
\[
\delta\ln\mathfrak C_{{\rm tr},*}=(1+\chi_{0,*})m_T,
\]
\[
\delta\ln\epsilon_\eta=-m_{K_\eta},
\]
\[
\delta\ln\mathfrak C_{{\rm nt},*}=m_\mu-F_*m_T-m_{K_\eta}.
\]
So the remaining dynamical theorem gap has been localized completely:
**the PDE only has to prove that the dependent microscopic coordinates follow the exact
Stage-170 orbit law generated by the five free similarity directions.**

## Stage 180 — consolidated adiabatic-elastic orbit-lock verdict

The strengthened boundary rule has now been compiled all the way down to a single exact
microscopic finish line.

- adiabatic wall freezes the thermal wall datum,
- elastic/no-fraying removes the scalar off-bundle escape,
- and the remaining first-order defect is nothing but the mismatch of the dependent triple
  \((\Delta_T,\Delta_{K_\eta},\Delta_\mu)\) from the exact single-orbit law.

Equivalently, the adiabatic-elastic moving-throat branch is locked to one exact
\(\mathcal G_*\)-orbit iff
\[
\Theta_1=\Xi_1=\mathcal R_1+\Xi_1=0,
\]
or microscopically iff
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu)
\]
follow the Stage-170 orbit law.

## Bottom line after Stage 180

The theorem gap is now narrower than “solve the PDE.” It is:

> show that the completed moving-throat dynamics makes the dependent microscopic triple
> \((T_U,K_\eta,\mu_W)\) follow the exact Stage-170 orbit law on the adiabatic-elastic branch.

If that happens, the first-order weak-axisymmetric defect vanishes automatically and the branch
stays on a single exact \(\mathcal G_*\)-orbit.
