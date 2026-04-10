# Step 9 — Effective transfer shape collapse and the one-port continuum branch

## Goal

Step 8 reduced the remaining quartic anomaly layer to the outgoing-weighted transfer-shape slope
```math
\Xi_1 = 2\sum_r \rho_r^{(N)}\tau_r,
\qquad
N_{A,0}^{(r)} = K_A\,\mathcal T_{A,r}^2.
```
That was already sharp, but it still carried the theorem gate as a weighted many-port average.

The next honest question is therefore:

> **does that weighted average collapse to the slope of one single effective transfer shape, and if so, what is that shape on the actual minimal continuum branch?**

This step answers that.

---

## Inputs carried forward

### From Step 8

For each active outgoing port,
```math
N_{A,0}^{(r)} = K_A\,\mathcal T_{A,r}^2,
```
with weak-axisymmetric slope
```math
\delta\ln \mathcal T_{A,r} = \epsilon\lambda_A\,\tau_r.
```
So the remaining grouped defect is already
```math
\Xi_1 = 2\sum_r \rho_r^{(N)}\tau_r,
\qquad
\rho_r^{(N)} = \frac{N_{A,0}^{(r)}}{\sum_s N_{A,0}^{(s)}}.
```

### From the moving-throat continuum branch

On the minimal isotropic continuum branch there is only one active outgoing port and the actual static transfer coefficient is written as
```math
N_{A,0}=\beta_{0,A},
\qquad
K_A=K_{0,A}.
```
The continuum formulas carried forward are
```math
K_{0,A} = \frac{K_{\eta,A}^{(\mathrm{eff})}}{\mu_{\eta,A}},
```
```math
\beta_{0,A}
=
\frac{\mu_{W,A}}{\mu_{\eta,A}}
\frac{K_{\eta,A}^{(\mathrm{eff})}}{K_{W,A}^{(\mathrm{eff})}}
\frac{Z_{W,A}(1+\rho_A)^2}{(1-\epsilon_{W,A})^2}.
```
These are the only branch inputs needed here.

---

## Step 9A — Exact collapse from many ports to one effective transfer shape

Summing the portwise factorization gives
```math
N_{A,0}
=
\sum_r N_{A,0}^{(r)}
=
K_A\sum_r \mathcal T_{A,r}^2.
```
So define the **effective wall-normalized transfer shape** by
```math
\boxed{
\mathcal T_{\mathrm{eff},A}^2
:=
\sum_r \mathcal T_{A,r}^2
=
\frac{N_{A,0}}{K_A}.
}
```
Now perturb weak-axisymmetrically:
```math
\delta\ln \mathcal T_{A,r} = \epsilon\lambda_A\,\tau_r.
```
Then
```math
\frac{\delta\ln \mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A}
=
\frac{2\sum_r \mathcal T_r^2\tau_r}{\sum_s \mathcal T_s^2}
=
2\sum_r \rho_r^{(N)}\tau_r,
```
because
```math
\rho_r^{(N)}
=
\frac{N_0^{(r)}}{N_0}
=
\frac{\mathcal T_r^2}{\sum_s \mathcal T_s^2}.
```
Therefore
```math
\boxed{
\Xi_1
=
\frac{\delta\ln \mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A}.
}
```
Equivalently, if
```math
\delta\ln \mathcal T_{\mathrm{eff},A} = \epsilon\lambda_A\,\tau_{\mathrm{eff}},
```
then
```math
\boxed{
\tau_{\mathrm{eff}}=\sum_r \rho_r^{(N)}\tau_r,
\qquad
\Xi_1 = 2\tau_{\mathrm{eff}}.
}
```

So the whole remaining grouped defect is already the slope of **one** effective transfer shape.

---

## Step 9B — Exact one-port continuum transfer shape

On the actual minimal isotropic continuum branch there is only one active outgoing port, so
```math
\mathcal T_{\mathrm{eff},A}=\mathcal T_A,
\qquad
N_{A,0}=\beta_{0,A},
\qquad
K_A=K_{0,A}.
```
Hence
```math
\mathcal T_A^2 = \frac{\beta_{0,A}}{K_{0,A}}.
```
Substituting the continuum formulas gives immediately
```math
\boxed{
\mathcal T_A^2
=
\frac{\mu_{W,A}}{K_{W,A}^{(\mathrm{eff})}}
\frac{Z_{W,A}(1+\rho_A)^2}{(1-\epsilon_{W,A})^2}.
}
```
Using
```math
\Omega_{W,A}^2 = \frac{K_{W,A}^{(\mathrm{eff})}}{\mu_{W,A}},
```
this becomes
```math
\boxed{
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2}.
}
```

This is the first exact actual-port formula for the effective transfer shape on the real continuum branch.

A useful physical consequence is immediate:

**the wall–U dressing ratio `\epsilon_\eta` drops out of the direct continuum transfer-shape formula.**

At the direct port level the defect only sees

- the wall-to-mixed overlap ratio `Z_W`,
- the interference ratio `\rho`,
- the mixed-sector blocking ratio `\epsilon_W`,
- the mixed frequency scale `\Omega_W`.

So the wall–U dressing lane affects the port transfer only indirectly.

---

## Step 9C — Direct weak-axisymmetric slope law in continuum variables

Perturb the one-port continuum formula weak-axisymmetrically by
```math
\delta\ln Z_{W,A} = \epsilon\lambda_A\,\zeta_W,
```
```math
\delta\ln \Omega_{W,A}^2 = \epsilon\lambda_A\,\omega_W,
```
```math
\delta \rho_A = \epsilon\lambda_A\,\rho_1,
```
```math
\delta\epsilon_{W,A} = \epsilon\lambda_A\,\varepsilon_W.
```
Then from
```math
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2},
```
we obtain the exact one-port weak-axisymmetric defect law
```math
\boxed{
\Xi_1
=
\frac{\delta\ln \mathcal T_A^2}{\epsilon\lambda_A}
=
\zeta_W - \omega_W + \frac{2\rho_1}{1+\rho} + \frac{2\varepsilon_W}{1-\epsilon_W}.
}
```
Equivalently,
```math
\boxed{
\tau
=
\frac12\,\Xi_1
=
\frac12\left(
\zeta_W - \omega_W + \frac{2\rho_1}{1+\rho} + \frac{2\varepsilon_W}{1-\epsilon_W}
\right).
}
```

So on the actual minimal continuum port, the remaining grouped defect is controlled by **four** direct microscopic drift channels:

1. wall-to-mixed overlap drift `\zeta_W`,
2. mixed-frequency drift `\omega_W`,
3. interference-ratio drift `\rho_1`,
4. mixed blocking drift `\varepsilon_W`.

This is the first exact actual-port answer to the Step-8 “compute `\tau_r`” question.

---

## Step 9D — Selected-branch reformulation

The same one-port transfer shape can be rewritten in terms of the selected-branch placement variables.

The carried demand-ratio formula is
```math
R_{\mathrm{target},A}
=
\Lambda_A
\frac{(1-\epsilon_{\eta,A})(1-\epsilon_{W,A})^2}{Z_{W,A}(1+\rho_A)^2},
```
with
```math
\Lambda_A
=
\frac{27\pi^2 G c_s^5 K_{W,A}^{(\mathrm{eff})}}{20 a^5 c^5 \mu_{W,A}}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\Omega_{W,A}^2.
```
So the transfer shape becomes
```math
\boxed{
\mathcal T_A^2
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}
\frac{1-\epsilon_{\eta,A}}{R_{\mathrm{target},A}}.
}
```

At linear grouped weak-axisymmetric order, scalar observables such as `a` and `c_s` do not shift on the isotropic branch, so the front factor is inert at this order. If we define
```math
\delta\epsilon_{\eta,A} = \epsilon\lambda_A\,\eta_1,
\qquad
\delta\ln R_{\mathrm{target},A} = \epsilon\lambda_A\,\mathcal R_1,
```
then the same defect becomes
```math
\boxed{
\Xi_1
=
-\frac{\eta_1}{1-\epsilon_\eta} - \mathcal R_1.
}
```

So the one-port defect has a second exact interpretation:

it is the mismatch between

- the wall–U dressing drift, and
- the selected-branch demand-ratio drift.

---

## Step 9E — Exact one-port zero-defect theorem

The direct zero-defect condition is now explicit:
```math
\boxed{
\delta\ln\!\left[
\frac{Z_W(1+\rho)^2}{\Omega_W^2(1-\epsilon_W)^2}
\right] = 0.
}
```
Equivalently, in selected-branch variables,
```math
\boxed{
\delta\ln R_{\mathrm{target}}
=
-\frac{\delta\epsilon_\eta}{1-\epsilon_\eta}.
}
```
So the one-port continuum theorem is genuinely bidirectional.

Two useful corollaries follow immediately.

### Target-rigid branch

If the selected-branch demand ratio is weak-axisymmetrically rigid,
```math
\delta\ln R_{\mathrm{target}}=0,
```
then zero grouped defect requires
```math
\boxed{\delta\epsilon_\eta = 0.}
```
So on a target-rigid branch, the wall–U dressing ratio itself must be lane-rigid.

### Wall–U-rigid branch

If the wall–U dressing ratio is weak-axisymmetrically rigid,
```math
\delta\epsilon_\eta = 0,
```
then zero grouped defect requires
```math
\boxed{\delta\ln R_{\mathrm{target}}=0.}
```
So on an `\epsilon_\eta`-rigid branch, the selected demand ratio must be lane-rigid.

---

## Step 9F — Direct quartic anomaly gate on the actual one-port branch

The carried quartic anomaly target is still
```math
\Xi_1 = \Lambda_1,
\qquad
\Lambda_1 \approx 0.279605891931464.
```
So the exact one-port direct-continuum anomaly condition is
```math
\boxed{
\zeta_W - \omega_W + \frac{2\rho_1}{1+\rho} + \frac{2\varepsilon_W}{1-\epsilon_W}
=
\Lambda_1.
}
```
Equivalently,
```math
\boxed{
\tau = \frac{\Lambda_1}{2} \approx 0.139802945965732.
}
```
The selected-branch form is just as sharp:
```math
\boxed{
-\frac{\eta_1}{1-\epsilon_\eta} - \mathcal R_1 = \Lambda_1.
}
```

A useful **reference balanced split** in the direct continuum variables is obtained by sharing the defect equally among the four rescaled drift channels. That gives
```math
\boxed{
\zeta_W = \frac{\Lambda_1}{4},
\qquad
\omega_W = -\frac{\Lambda_1}{4},
\qquad
\rho_1 = \frac{1+\rho}{8}\,\Lambda_1,
\qquad
\varepsilon_W = \frac{1-\epsilon_W}{8}\,\Lambda_1.
}
```
Numerically,
```math
\frac{\Lambda_1}{4} \approx 0.0699014729828660.
```
A corresponding balanced split in the selected-branch variables is
```math
\boxed{
\eta_1 = -\frac{1-\epsilon_\eta}{2}\,\Lambda_1,
\qquad
\mathcal R_1 = -\frac{\Lambda_1}{2}.
}
```
Numerically,
```math
-\frac{\Lambda_1}{2} \approx -0.139802945965732.
```
These are not unique solutions; they are just the cleanest reference branch targets.

---

## Step 9G — Consequence for the grouped weak-axisymmetric quadrupole defect

The grouped weak-axisymmetric pattern is still
```math
\Delta_Q^{(20)} = \epsilon\,\Xi_1,
\qquad
\Delta_Q^{(21)} = \frac{\epsilon}{2}\,\Xi_1,
\qquad
\Delta_Q^{(22)} = -\epsilon\,\Xi_1.
```
Because `\Xi_1 = 2\tau` on the one-port branch, this becomes
```math
\boxed{
\Delta_Q^{(20)} = 2\epsilon\,\tau,
\qquad
\Delta_Q^{(21)} = \epsilon\,\tau,
\qquad
\Delta_Q^{(22)} = -2\epsilon\,\tau.
}
```
So once the actual one-port transfer-shape slope is known, the full weak-axisymmetric grouped normalization pattern is fixed immediately.

---

## Reduced verdict

Step 8 said the remaining quartic anomaly layer is the outgoing-weighted transfer-shape slope.
Step 9 sharpens that result in two decisive ways.

First, the many-port weighted average is itself exactly the slope of a **single effective transfer shape**,
```math
\mathcal T_{\mathrm{eff}}^2 = \frac{N_0}{K}.
```

Second, on the actual minimal one-port continuum branch, that transfer shape is no longer abstract. It is explicitly
```math
\mathcal T^2
=
\frac{Z_W(1+\rho)^2}{\Omega_W^2(1-\epsilon_W)^2}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}
\frac{1-\epsilon_\eta}{R_{\mathrm{target}}}.
```

So the remaining theorem gate is no longer
> “compute all raw outgoing-port slopes.”

It is now
> **determine whether the actual weak-axisymmetric moving-throat branch keeps this single continuum transfer shape rigid, and if not, which of its four direct microscopic drift channels carries the required `\Lambda_1`-sized defect.**

That is the direct continuation point.
