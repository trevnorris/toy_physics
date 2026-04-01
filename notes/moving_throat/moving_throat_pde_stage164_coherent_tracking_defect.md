
# Moving-Throat PDE — Stage 164: Coherent Tracking-Branch Weak-Axisymmetric Defect Law and the Support-Blindness Theorem

## Purpose

Stage 163 reduced the remaining grouped weak-axisymmetric `2.5`PN defect to the logarithmic slope of a single effective transfer shape,
\[
\Xi_1=\frac{\delta\ln \mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A},
\qquad
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2}
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}\,
\frac{1-\epsilon_{\eta,A}}{R_{\mathrm{target},A}}.
\]

Earlier coherent-kernel stages already gave the **actual** one-port tracking branch of the first local D/N throat kernel:
\[
R_{\mathrm{tr}}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\]
\[
\epsilon
=
\epsilon_W^{(\mathrm{split})}
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right],
\]
\[
M_{\mathrm{tr}}
=
M_{\mathrm{mix}}\,S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\]
\[
R_{\mathrm{target}}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
\]

So the next honest step is now very sharp:

1. substitute the actual coherent tracking-branch map into the Stage-163 transfer-shape formula,
2. compute the exact weak-axisymmetric defect law in the **physical** branch variables,
3. and determine whether coherent support can help cancel the grouped defect.

The main result is that it cannot.

On the coherent local D/N tracking branch,
\[
\boxed{
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}\,
\frac{1-\epsilon_\eta}{R_{\mathrm{target}}},
}
\]
and therefore the weak-axisymmetric grouped defect is
\[
\boxed{
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}+\frac{2\epsilon_1}{1-\epsilon},
}
\]
with the exact split-blocking drift
\[
\boxed{
\epsilon_1
=
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}.
}
\]

Most importantly, the support-enhancement ratio \(\zeta\) drops out **identically**:
\[
\boxed{
\partial_\zeta \ln \mathcal T^2 = 0,
\qquad
\partial_\zeta \ln R_{\mathrm{target}} = 0,
\qquad
\partial_\zeta \Xi_1 = 0.
}
\]

So coherent support can enhance the baseline and help the steady normalization test, but it cannot repair or spoil the weak-axisymmetric grouped defect at first order.
That defect is carried only by the mixed/outgoing placement variables.

---

## 1. Exact coherent tracking-branch substitution

On the coherent local D/N branch the common interference ratio is
\[
\rho_0=\sigma_0=\chi_0,
\]
the mixed blocking ratio is the split quantity
\[
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right],
\]
and the selected-branch demand ratio is
\[
R_{\mathrm{target}}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
\]

Substituting this directly into the Stage-163 selected-branch transfer-shape law gives
\[
\mathcal T^2
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}\,
\frac{1-\epsilon_\eta}{R_{\mathrm{target}}}
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\]
because
\[
\Lambda
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}\,\Omega_W^2.
\]

So on the actual coherent branch the transfer shape depends only on

- the wall-to-mixed overlap \(Z_W\),
- the mixed frequency scale \(\Omega_W\),
- the common interference ratio \(\chi_0\),
- and the split mixed blocking ratio \(\epsilon\).

The support-enhancement factor \(S(\zeta;\epsilon)\) and the total baseline \(M_{\mathrm{tr}}\) do **not** enter.

---

## 2. Exact support-blindness theorem

The coherent branch support enhancement is
\[
M_{\mathrm{tr}} = M_{\mathrm{mix}}\,S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon}.
\]
But the transfer shape itself is
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}.
\]

Therefore
\[
\boxed{
\frac{\partial \mathcal T^2}{\partial \zeta}=0,
\qquad
\frac{\partial \ln \mathcal T^2}{\partial \zeta}=0.
}
\]
Using
\[
R_{\mathrm{target}}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]
we also have
\[
\boxed{
\frac{\partial R_{\mathrm{target}}}{\partial \zeta}=0,
\qquad
\frac{\partial \ln R_{\mathrm{target}}}{\partial \zeta}=0.
}
\]

Hence the Stage-163 defect scalar is exactly support-blind:
\[
\boxed{
\partial_\zeta \Xi_1 = 0.
}
\]

### Physical meaning

This is a very sharp separation of roles:

- the coherent support lane changes the **baseline** \(M_{\mathrm{tr}}\),
- but it does not change the effective transfer shape,
- so it can help satisfy the steady normalization demand,
- yet it cannot help cancel the weak-axisymmetric grouped defect.

That is the strongest theorem-level narrowing obtained in this stage.

---

## 3. Exact weak-axisymmetric defect law in physical branch variables

Now perturb the coherent tracking branch weak-axisymmetrically:
\[
\delta\ln Z_{W,A}=\epsilon\lambda_A\,\zeta_Z,
\qquad
\delta\ln \Omega_{W,A}^2=\epsilon\lambda_A\,\omega_W,
\]
\[
\delta\chi_{0,A}=\epsilon\lambda_A\,\chi_1,
\qquad
\delta\epsilon_{W,A}=\epsilon\lambda_A\,\varepsilon_W,
\qquad
\delta\delta_{U,A}=\epsilon\lambda_A\,\delta_{U,1},
\qquad
\delta\epsilon_{\eta,A}=\epsilon\lambda_A\,\eta_1.
\]

Then
\[
\Xi_1
=
\frac{\delta\ln\mathcal T_A^2}{\epsilon\lambda_A}
=
\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}+\frac{2\epsilon_1}{1-\epsilon},
\]
where \(\epsilon_1\) is the weak-axisymmetric drift of the split blocking ratio
\[
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right].
\]
Differentiating gives the exact transport law
\[
\boxed{
\epsilon_1
=
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}.
}
\]
So the final exact coherent-branch defect law is
\[
\boxed{
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}
+\frac{2}{1-\epsilon}
\left[
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}
\right].
}
\]

This formula is the direct physical continuation of Stage 163.

### Interpretation

The first-order grouped defect is carried only by five microscopic drift channels:

1. wall-to-mixed overlap drift \( \zeta_Z \),
2. mixed-frequency drift \( \omega_W \),
3. common interference-ratio drift \( \chi_1 \),
4. bare mixed-blocking drift \( \varepsilon_W \),
5. split-\(U\) axial drift \( \delta_{U,1} \).

The coherent support amplitude ratio \(\zeta\) is absent.

---

## 4. Selected-branch reformulation

The Stage-163 selected-branch identity remains exact:
\[
\Xi_1
=
-\frac{\eta_1}{1-\epsilon_\eta}
-
\mathcal R_1,
\qquad
\mathcal R_1:=\frac{\delta\ln R_{\mathrm{target},A}}{\epsilon\lambda_A}.
\]
Using the coherent-branch expression for \(R_{\mathrm{target}}\),
\[
R_{\mathrm{target}}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]
we obtain
\[
\boxed{
\mathcal R_1
=
\omega_W
-\frac{\eta_1}{1-\epsilon_\eta}
-\zeta_Z
-\frac{2\chi_1}{1+\chi_0}
-\frac{2\epsilon_1}{1-\epsilon}.
}
\]
So the direct-port and selected-branch descriptions continue to agree exactly on the coherent tracking branch.

---

## 5. Tracking-factor drift is not the defect carrier

The coherent tracking factor is
\[
R_{\mathrm{tr}}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{\chi_0+\delta_U+1}{(1+\chi_0)(1+\delta_U)}.
\]
Its weak-axisymmetric logarithmic drift is
\[
\delta\ln R_{\mathrm{tr},A}
=
\epsilon\lambda_A\,\Theta_1,
\]
with the exact coefficient
\[
\boxed{
\Theta_1
=
-\frac{\chi_0(1+\chi_0)\,\delta_{U,1}
+\delta_U(1+\delta_U)\,\chi_1}
{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}.
}
\]

This is useful because it shows something nontrivial:

- the tracking-factor drift depends only on \((\chi_1,\delta_{U,1})\),
- but the grouped defect \(\Xi_1\) also depends on \((\zeta_Z,\omega_W,\varepsilon_W)\).

So exact tracking-factor rigidity,
\[
\Theta_1=0,
\]
is **not** sufficient to kill the grouped defect.

For example, if
\[
\chi_1=\delta_{U,1}=0,
\]
then
\[
\Theta_1=0,
\]
but
\[
\Xi_1
=
\zeta_Z-\omega_W
+
\frac{2}{1-\epsilon}
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W,
\]
which need not vanish.

Thus the weak-axisymmetric grouped defect is not controlled solely by the tracking factor.
It is controlled by the transfer-shape variables themselves.

---

## 6. Best current theorem statement after Stage 164

On the actual coherent local D/N tracking branch:

1. the effective transfer shape is exactly
   \[
   \mathcal T^2
   =
   \frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
   \]
   with
   \[
   \epsilon
   =
   \epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right],
   \]
2. the coherent support-enhancement ratio \(\zeta\) drops out identically,
   \[
   \partial_\zeta \ln\mathcal T^2
   =
   \partial_\zeta \ln R_{\mathrm{target}}
   =
   \partial_\zeta \Xi_1
   =0,
   \]
3. the full weak-axisymmetric grouped defect is
   \[
   \Xi_1
   =
   \zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}
   +
   \frac{2}{1-\epsilon}
   \left[
   \left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
   -
   \frac{2\epsilon_W}{11(1+\delta_U)^2}\delta_{U,1}
   \right],
   \]
4. and exact tracking-factor rigidity is not enough to kill it.

So the remaining weak-axisymmetric normalization problem is now narrower than “compute the support branch.”

It is:

> compute the actual weak-axisymmetric drifts of the mixed/outgoing placement variables
> \((Z_W,\Omega_W^2,\chi_0,\epsilon_W,\delta_U)\)
> on the coherent moving-throat branch.

That is the direct continuation point.
