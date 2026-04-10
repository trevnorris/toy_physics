# Step 10 — Coherent tracking-branch substitution and the support-blindness theorem

## Goal

Step 9 reduced the remaining quartic anomaly layer to the weak-axisymmetric slope of a **single** one-port continuum transfer shape,
```math
\mathcal T^2
=
\frac{Z_W(1+\rho)^2}{\Omega_W^2(1-\epsilon_W)^2}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\frac{1-\epsilon_\eta}{R_{\mathrm{target}}},
\qquad
\Xi_1 = \delta\ln \mathcal T^2.
```
The next honest move is therefore not to invent a new transport layer by hand. It is to substitute the **actual coherent local D/N tracking branch** into that one-port law and see what survives.

That is what this step does.

The decisive outcome is:

> **coherent support can raise the steady normalization baseline, but it drops out identically from the first grouped weak-axisymmetric defect.**

So the next live g-2 calculation is no longer “how much coherent support is present?” It is “how the mixed/outgoing placement variables drift across the grouped weak-axisymmetric branch.”

---

## Inputs carried forward

### From Step 9

The one-port continuum transfer-shape law is
```math
\mathcal T^2
=
\frac{Z_W(1+\rho)^2}{\Omega_W^2(1-\epsilon_W)^2}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\frac{1-\epsilon_\eta}{R_{\mathrm{target}}}.
```
Its grouped defect is
```math
\Xi_1 = \delta\ln \mathcal T^2,
```
and the quartic anomaly target remains
```math
\Xi_1 = \Lambda_1,
\qquad
\Lambda_1 \approx 0.279605891931464,
\qquad
\tau = \frac{\Lambda_1}{2} \approx 0.139802945965732.
```

### From the coherent local D/N tracking branch

The exact coherent branch data are
```math
R_{\mathrm{tr}}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
```
```math
\epsilon
=
\epsilon_W^{(\mathrm{split})}
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right],
```
```math
M_{\mathrm{tr}} = M_{\mathrm{mix}}\,S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
```
```math
R_{\mathrm{target}}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\qquad
\Lambda = \frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\Omega_W^2.
```

---

## Step 10A — Exact coherent tracking-branch substitution

Substituting the coherent-branch demand ratio into the selected-branch transfer-shape law gives immediately
```math
\mathcal T^2
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\frac{1-\epsilon_\eta}{R_{\mathrm{target}}}
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}.
```
So on the actual coherent tracking branch, the direct transfer shape depends only on

- the wall-to-mixed overlap ratio `Z_W`,
- the mixed frequency scale `\Omega_W`,
- the common interference ratio `\chi_0`,
- the split mixed blocking ratio `\epsilon`.

The support-enhancement factor `S(\zeta;\epsilon)` does **not** appear.

---

## Step 10B — Exact support-blindness theorem

The coherent support lane changes the total baseline through
```math
M_{\mathrm{tr}} = M_{\mathrm{mix}}\,S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
```
with
```math
\frac{\partial S}{\partial \zeta} = \frac{1-\epsilon}{(1-\zeta\epsilon)^2} > 0
```
on the physical branch. So support is a strict baseline enhancer.

But the direct transfer shape itself is
```math
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
```
and the selected-branch demand ratio is
```math
R_{\mathrm{target}}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
```
Neither contains `\zeta`. Therefore
```math
\boxed{
\frac{\partial \ln \mathcal T^2}{\partial \zeta}=0,
\qquad
\frac{\partial \ln R_{\mathrm{target}}}{\partial \zeta}=0,
\qquad
\frac{\partial \Xi_1}{\partial \zeta}=0.
}
```
This is the exact **support-blindness theorem**.

A useful contrast is that the steady normalization product still feels support:
```math
R_{\mathrm{target}} M_{\mathrm{tr}}
=
\frac{8\Lambda(1-\epsilon)}{\pi^2}
S(\zeta;\epsilon).
```
So support can help the steady normalization test, but it cannot repair or spoil the first grouped weak-axisymmetric defect.

---

## Step 10C — Exact weak-axisymmetric defect law in physical branch variables

To avoid confusion with the split blocking ratio `\epsilon`, write the grouped perturbation parameter as `s\lambda_A`.
Take
```math
\delta\ln Z_{W,A}=s\lambda_A\,\zeta_Z,
\qquad
\delta\ln \Omega_{W,A}^2=s\lambda_A\,\omega_W,
```
```math
\delta\chi_{0,A}=s\lambda_A\,\chi_1,
\qquad
\delta\epsilon_{W,A}=s\lambda_A\,\varepsilon_W,
\qquad
\delta\delta_{U,A}=s\lambda_A\,\delta_{U,1},
```
```math
\delta\epsilon_{\eta,A}=s\lambda_A\,\eta_1.
```
Because
```math
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right],
```
its weak-axisymmetric drift is
```math
\boxed{
\epsilon_1
=
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}.
}
```
Then
```math
\ln \mathcal T^2
=
\ln Z_W + 2\ln(1+\chi_0) - \ln \Omega_W^2 - 2\ln(1-\epsilon),
```
so the grouped defect becomes
```math
\boxed{
\Xi_1
=
\zeta_Z - \omega_W + \frac{2\chi_1}{1+\chi_0} + \frac{2\epsilon_1}{1-\epsilon}.
}
```
Expanding `\epsilon_1` gives the full physical-branch law
```math
\boxed{
\Xi_1
=
\zeta_Z - \omega_W + \frac{2\chi_1}{1+\chi_0}
+
\frac{2}{1-\epsilon}
\left[
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}
\right].
}
```
So the first grouped defect is carried only by the mixed/outgoing placement drifts

1. `\zeta_Z` — wall-to-mixed overlap drift,
2. `\omega_W` — mixed-frequency drift,
3. `\chi_1` — common interference-ratio drift,
4. `\varepsilon_W` — bare mixed-blocking drift,
5. `\delta_{U,1}` — split-`U` axial drift.

The coherent support ratio `\zeta` is absent.

---

## Step 10D — Selected-branch reformulation

The Step-9 selected-branch identity remains exact on the coherent branch:
```math
\Xi_1 = -\frac{\eta_1}{1-\epsilon_\eta} - \mathcal R_1,
\qquad
\mathcal R_1 := \delta\ln R_{\mathrm{target},A}.
```
Since
```math
R_{\mathrm{target}}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
```
its weak-axisymmetric drift is
```math
\boxed{
\mathcal R_1
=
\omega_W - \frac{\eta_1}{1-\epsilon_\eta}
- \zeta_Z - \frac{2\chi_1}{1+\chi_0} - \frac{2\epsilon_1}{1-\epsilon}.
}
```
So the direct-port and selected-branch descriptions still agree exactly.

---

## Step 10E — Tracking-factor drift is not sufficient

The coherent tracking factor is
```math
R_{\mathrm{tr}}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)}.
```
Its weak-axisymmetric logarithmic drift is
```math
\boxed{
\Theta_1
=
-\frac{\chi_0(1+\chi_0)\,\delta_{U,1} + \delta_U(1+\delta_U)\,\chi_1}
{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}.
}
```
So `\Theta_1` depends only on `(\chi_1,\delta_{U,1})`, whereas `\Xi_1` still depends on `(\zeta_Z,\omega_W,\varepsilon_W)` as well.

A simple explicit slice is
```math
\chi_1 = 0,
\qquad
\delta_{U,1}=0.
```
Then
```math
\Theta_1 = 0,
```
but
```math
\Xi_1
=
\zeta_Z - \omega_W
+
\frac{2}{1-\epsilon}
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W,
```
which need not vanish.

So exact tracking-factor rigidity is **not** sufficient to kill the grouped defect.

---

## Step 10F — Coherent-branch quartic anomaly gate

The carried quartic target is still
```math
\Xi_1 = \Lambda_1,
\qquad
\Lambda_1 \approx 0.279605891931464,
\qquad
\tau = \frac{\Lambda_1}{2} \approx 0.139802945965732.
```
So the coherent tracking branch must satisfy
```math
\boxed{
\zeta_Z - \omega_W + \frac{2\chi_1}{1+\chi_0}
+
\frac{2}{1-\epsilon}
\left[
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}
\right]
=
\Lambda_1.
}
```
This is the first exact coherent-branch quartic anomaly gate written entirely in physical branch variables.

---

## Reduced verdict

Step 9 said the remaining quartic anomaly layer lives in one actual one-port transfer shape.
Step 10 sharpens that result in the most useful possible way.

First, on the actual coherent local D/N tracking branch,
```math
\mathcal T^2 = \frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
```
with `\epsilon = \epsilon_W[1 - (2/11)\delta_U/(1+\delta_U)]`.

Second, the coherent support lane drops out identically from
```math
\mathcal T^2,
\qquad
R_{\mathrm{target}},
\qquad
\Xi_1.
```
So support can raise the steady baseline but cannot carry the first grouped defect.

Third, the full coherent weak-axisymmetric defect is now carried only by the mixed/outgoing placement drifts
```math
(\zeta_Z,\omega_W,\chi_1,\varepsilon_W,\delta_{U,1}),
```
and exact tracking rigidity by itself is not enough to kill it.

So the next clean move is no longer “analyze coherent support.”
It is:

> **push the coherent branch one layer deeper, down to the microscopic coherent-kernel slippages, and separate the true tracking and nontracking defect coordinates.**

That is the direct continuation point.
