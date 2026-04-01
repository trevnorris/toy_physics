# Moving-Throat PDE — Stage 163: Effective Transfer-Shape Collapse and the Actual Continuum-Kernel Slope Law

## Purpose

Stage 162 reduced the remaining weak-axisymmetric grouped `2.5`PN defect to the outgoing-weighted average of the wall-normalized transfer-shape slopes,
\[
\Xi_1 = 2\sum_r \rho_r^{(N)}\,\tau_r.
\]
That was already very sharp, but it still left the theorem gate phrased as a weighted sum over individual ports.

At the same time, the earlier continuum-kernel stages already gave an explicit one-port isotropic moving-throat branch in terms of the dimensionless ratios
\[
\epsilon_\eta,\qquad \epsilon_W,\qquad \rho,\qquad Z_W,
\]
together with
\[
K_0=\frac{K_\eta^{(\mathrm{eff})}}{\mu_\eta},
\qquad
\beta_0
=
\frac{\mu_W}{\mu_\eta}\frac{K_\eta^{(\mathrm{eff})}}{K_W^{(\mathrm{eff})}}
\frac{Z_W(1+\rho)^2}{(1-\epsilon_W)^2}.
\]

So the next honest step is now clear:

1. collapse the weighted sum of port slopes to **one** effective wall-normalized transfer shape,
2. then evaluate that shape explicitly on the actual minimal continuum branch.

That is what this stage does.

The main result is that the whole remaining weak-axisymmetric grouped defect is the logarithmic slope of a single effective transfer shape,
\[
\boxed{
\Xi_1 = \frac{\delta\ln \mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A},
\qquad
\mathcal T_{\mathrm{eff},A}^2:=\sum_r \mathcal T_{A,r}^2=\frac{N_{A,0}}{K_A},
}
\]
and, on the actual minimal one-port isotropic continuum branch,
\[
\boxed{
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,
\frac{1-\epsilon_{\eta,A}}{R_{\mathrm{target},A}}.
}
\]

So the remaining theorem gap is no longer “compute many raw outgoing-port slopes.”
It is:

> determine whether the actual weak-axisymmetric moving-throat branch keeps this single effective transfer shape rigid.

---

## 1. Exact collapse from many ports to one effective transfer shape

Stage 162 already gave the exact portwise factorization
\[
N_{A,0}^{(r)}=K_A\,\mathcal T_{A,r}^2.
\]
Summing over active outgoing ports gives
\[
N_{A,0}
=
\sum_r N_{A,0}^{(r)}
=
K_A\sum_r \mathcal T_{A,r}^2.
\]
So it is natural to define the effective wall-normalized transfer shape
\[
\boxed{
\mathcal T_{\mathrm{eff},A}^2
:=
\sum_r \mathcal T_{A,r}^2
=
\frac{N_{A,0}}{K_A}.
}
\]

Now perturb weak-axisymmetrically:
\[
\delta\ln\mathcal T_{A,r}
=
\epsilon\lambda_A\,\tau_r.
\]
Then
\[
\frac{\delta\ln \mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A}
=
\frac{2\sum_r \mathcal T_r^2 \tau_r}{\sum_s \mathcal T_s^2}
=
2\sum_r \rho_r^{(N)}\tau_r,
\]
because the Stage-162 weights are exactly
\[
\rho_r^{(N)}=\frac{N_{0}^{(r)}}{N_0}=\frac{\mathcal T_r^2}{\sum_s \mathcal T_s^2}.
\]
Therefore
\[
\boxed{
\Xi_1
=
2\sum_r \rho_r^{(N)}\tau_r
=
\frac{\delta\ln \mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A}.
}
\]
Equivalently, if
\[
\delta\ln \mathcal T_{\mathrm{eff},A}
=
\epsilon\lambda_A\,\tau_{\mathrm{eff}},
\]
then
\[
\boxed{
\tau_{\mathrm{eff}}=\sum_r \rho_r^{(N)}\tau_r,
\qquad
\Xi_1=2\tau_{\mathrm{eff}}.
}
\]

So the whole remaining grouped defect is the slope of one effective transfer shape.

---

## 2. Exact one-port effective transfer shape on the actual minimal continuum branch

On the minimal isotropic continuum branch there is only one active outgoing port, so
\[
\mathcal T_{\mathrm{eff},A}=\mathcal T_A,
\qquad
N_{A,0}=\beta_{0,A},
\qquad
K_A=K_{0,A}.
\]
Using the Stage-21 continuum formulas,
\[
K_{0,A}=\frac{K_{\eta,A}^{(\mathrm{eff})}}{\mu_{\eta,A}},
\]
\[
\beta_{0,A}
=
\frac{\mu_{W,A}}{\mu_{\eta,A}}
\frac{K_{\eta,A}^{(\mathrm{eff})}}{K_{W,A}^{(\mathrm{eff})}}
\frac{Z_{W,A}(1+\rho_A)^2}{(1-\epsilon_{W,A})^2},
\]
we get immediately
\[
\boxed{
\mathcal T_A^2
=
\frac{\beta_{0,A}}{K_{0,A}}
=
\frac{\mu_{W,A}}{K_{W,A}^{(\mathrm{eff})}}
\frac{Z_{W,A}(1+\rho_A)^2}{(1-\epsilon_{W,A})^2}.
}
\]
Since
\[
\Omega_{W,A}^2=\frac{K_{W,A}^{(\mathrm{eff})}}{\mu_{W,A}},
\]
this is equivalently
\[
\boxed{
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2}.
}
\]

This is the first exact actual-port formula for the Stage-162 transfer shape.

A useful physical consequence is immediate:

the wall–U dressing ratio \(\epsilon_\eta\) drops out of the direct continuum transfer-shape formula. At the direct port level, the defect only sees

- the wall-to-mixed overlap ratio \(Z_W\),
- the interference ratio \(\rho\),
- the mixed-sector blocking ratio \(\epsilon_W\),
- and the mixed frequency scale \(\Omega_W\).

So the wall–U dressing lane affects the port transfer only indirectly, not directly.

---

## 3. Direct weak-axisymmetric slope law in continuum variables

Now perturb the one-port continuum expression weak-axisymmetrically.

Write
\[
\delta\ln Z_{W,A}=\epsilon\lambda_A\,\zeta_W,
\]
\[
\delta\ln \Omega_{W,A}^2=\epsilon\lambda_A\,\omega_W,
\]
\[
\delta \rho_A=\epsilon\lambda_A\,\rho_1,
\]
\[
\delta\epsilon_{W,A}=\epsilon\lambda_A\,\varepsilon_W.
\]
Then from
\[
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2},
\]
we obtain
\[
\boxed{
\Xi_1
=
\frac{\delta\ln\mathcal T_A^2}{\epsilon\lambda_A}
=
\zeta_W-\omega_W+\frac{2\rho_1}{1+\rho}+\frac{2\varepsilon_W}{1-\epsilon_W}.
}
\]
Equivalently,
\[
\boxed{
\tau
=
\frac12\,\Xi_1
=
\frac12\left(
\zeta_W-\omega_W+\frac{2\rho_1}{1+\rho}+\frac{2\varepsilon_W}{1-\epsilon_W}
\right).
}
\]

So on the actual minimal continuum port, the remaining weak-axisymmetric grouped defect is controlled by four microscopic drift channels:

1. wall-to-mixed overlap drift \( \zeta_W \),
2. mixed-frequency drift \( \omega_W \),
3. interference-ratio drift \( \rho_1 \),
4. mixed blocking drift \( \varepsilon_W \).

This is the first exact actual-port answer to the Stage-162 “compute \(\tau_r\)” question.

---

## 4. Selected-branch reformulation

The same transfer shape can be rewritten in terms of the selected-branch continuum placement variables.

Stage 21 gave
\[
R_{\mathrm{target},A}
=
\Lambda_A
\frac{(1-\epsilon_{\eta,A})(1-\epsilon_{W,A})^2}{Z_{W,A}(1+\rho_A)^2},
\]
with
\[
\Lambda_A
=
\frac{27\pi^2 G c_s^5 K_{W,A}^{(\mathrm{eff})}}{20 a^5 c^5 \mu_{W,A}}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\Omega_{W,A}^2.
\]
So the one-port transfer shape becomes
\[
\boxed{
\mathcal T_A^2
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,
\frac{1-\epsilon_{\eta,A}}{R_{\mathrm{target},A}}.
}
\]

This is the sharpest bridge between the direct port language and the selected-branch placement language.

At linear grouped \(P_2\) order, scalar observables such as \(a\) and \(c_s\) do not shift on the isotropic branch, so the front factor is inert at this order. Therefore, if we define
\[
\delta\epsilon_{\eta,A}=\epsilon\lambda_A\,\eta_1,
\qquad
\delta\ln R_{\mathrm{target},A}=\epsilon\lambda_A\,\mathcal R_1,
\]
then
\[
\boxed{
\Xi_1
=
-\frac{\eta_1}{1-\epsilon_\eta}
-
\mathcal R_1.
}
\]

So the same weak-axisymmetric defect has a second exact interpretation:

it is the mismatch between

- the wall–U dressing drift, and
- the selected-branch demand-ratio drift.

This is the point where the direct port language and the selected-branch placement language become exactly equivalent.

---

## 5. Zero-defect theorem on the actual one-port continuum branch

The direct zero-defect condition is now explicit:
\[
\boxed{
\delta\ln\!\left[
\frac{Z_W(1+\rho)^2}{\Omega_W^2(1-\epsilon_W)^2}
\right]=0.
}
\]
Equivalently, in selected-branch variables,
\[
\boxed{
\delta\ln R_{\mathrm{target}}
=
-\frac{\delta\epsilon_\eta}{1-\epsilon_\eta}.
}
\]

This gives two useful corollaries.

### 5.1 Target-rigid branch

If the selected-branch demand ratio is weak-axisymmetrically rigid,
\[
\delta\ln R_{\mathrm{target}}=0,
\]
then zero grouped defect requires
\[
\boxed{
\delta\epsilon_\eta=0.
}
\]
So on a target-rigid branch, the wall–U dressing ratio itself must be lane-rigid.

### 5.2 Wall–U-rigid branch

If the wall–U dressing ratio is weak-axisymmetrically rigid,
\[
\delta\epsilon_\eta=0,
\]
then zero grouped defect requires
\[
\boxed{
\delta\ln R_{\mathrm{target}}=0.
}
\]
So on an \(\epsilon_\eta\)-rigid branch, the selected demand ratio must be lane-rigid.

Thus the one-port continuum theorem is now genuinely bidirectional.

---

## 6. Consequence for the grouped weak-axisymmetric quadrupole defect

Stage 156 already fixed the grouped lane pattern of the remaining weak-axisymmetric normalization defect:
\[
\Delta_Q^{(20)}=\epsilon\,\Xi_1,
\qquad
\Delta_Q^{(21)}=\frac{\epsilon}{2}\,\Xi_1,
\qquad
\Delta_Q^{(22)}=-\epsilon\,\Xi_1.
\]
Because \(\Xi_1=2\tau\) on the one-port branch, this becomes
\[
\boxed{
\Delta_Q^{(20)}=2\epsilon\,\tau,
\qquad
\Delta_Q^{(21)}=\epsilon\,\tau,
\qquad
\Delta_Q^{(22)}=-2\epsilon\,\tau.
}
\]
So once the actual one-port transfer-shape slope is known, the full weak-axisymmetric grouped normalization pattern is fixed immediately.

This is the cleanest explicit continuation of the Stage-162 theorem chain.

---

## 7. What Stage 163 changes

Stage 162 already reduced the last grouped defect to the outgoing-weighted average transfer-shape slope.
Stage 163 sharpens that result in two decisive ways.

First, it shows that the many-port weighted average is itself exactly the slope of a single effective transfer shape,
\[
\mathcal T_{\mathrm{eff}}^2=\frac{N_0}{K}.
\]

Second, on the actual minimal one-port continuum branch, that transfer shape is no longer abstract. It is explicitly
\[
\mathcal T^2
=
\frac{Z_W(1+\rho)^2}{\Omega_W^2(1-\epsilon_W)^2}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,
\frac{1-\epsilon_\eta}{R_{\mathrm{target}}}.
\]

So the remaining theorem gate is no longer:

> “compute all raw outgoing-port slopes.”

It is now:

> determine whether the actual weak-axisymmetric moving-throat branch keeps this single continuum transfer shape rigid.

That is the direct continuation point.
