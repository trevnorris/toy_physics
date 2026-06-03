# Moving-Throat PDE — Stage 182: Microscopic Coherent-Kernel Slippage Decomposition and the Exact Tracking/Nontracking Split

## Purpose

Stage 181 reduced the weak-axisymmetric grouped defect on the actual coherent local D/N tracking branch to the physical placement variables
\[
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}
+\frac{2\epsilon_1}{1-\epsilon},
\qquad
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right],
\]
with
\[
\epsilon_1
=
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}.
\]
It also showed that the coherent support ratio \(\zeta\) drops out identically, and that exact tracking-factor rigidity is not enough to kill \(\Xi_1\).

So the next honest step is now very small:

> push the Stage 181 physical variables one layer deeper, down to the actual microscopic coherent-kernel couplings, and identify which microscopic slippages really carry the grouped defect.

That is what this stage does.

The main result is that the coherent weak-axisymmetric defect depends on only four microscopic slippage combinations:
\[
\Sigma_Z,
\qquad
\Sigma_\chi,
\qquad
\Sigma_\epsilon,
\qquad
\Sigma_\delta,
\]
while the selected-branch form introduces one additional dressing slippage \(\Sigma_\eta\).
Explicitly,
\[
\boxed{
\Xi_1
=
\Sigma_Z
+\frac{2\chi_0}{1+\chi_0}\,\Sigma_\chi
+\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\,\Sigma_\delta
\right].
}
\]
Moreover, the tracking-factor drift depends on only one particular combination
\[
\boxed{
\Sigma_{\rm tr}:=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi,
}
\]
through
\[
\boxed{
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr}.
}
\]
So the grouped defect splits exactly into a tracking piece and a genuinely nontracking mixed/outgoing piece.

The strongest sharpened conclusion is therefore:

- coherent support remains absent,
- exact tracking rigidity means \(\Sigma_{\rm tr}=0\),
- but the grouped defect still survives unless the nontracking microscopic slippages also vanish.

---

## 1. Microscopic coherent-kernel variables

Stage 30 already expressed the coherent local D/N branch in terms of the microscopic couplings
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
\qquad
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}},
\]
\[
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\qquad
\delta_U=\frac{\pi^2 T_U}{L^2 K_U},
\qquad
\Omega_W^2=\frac{K_W^{(\mathrm{eff})}}{\mu_W}.
\]
The support lane enters through
\[
\zeta=\frac{\lambda_\phi^2 K_W^{(\mathrm{eff})}}{\lambda_W^2 K_\phi^{(\mathrm{eff})}},
\]
but Stage 181 already showed that \(\zeta\) drops out of the grouped defect.

So the natural microscopic weak-axisymmetric logarithmic drifts are
\[
\delta\ln \lambda_{W,A}=\epsilon\lambda_A\,\lambda_1,
\qquad
\delta\ln c_{\eta U,A}=\epsilon\lambda_A\,c_1,
\qquad
\delta\ln \gamma_A=\epsilon\lambda_A\,\gamma_1,
\]
\[
\delta\ln K_{U,A}=\epsilon\lambda_A\,\kappa_U,
\qquad
\delta\ln K_{\eta,A}^{(\mathrm{eff})}=\epsilon\lambda_A\,\kappa_\eta,
\qquad
\delta\ln K_{W,A}^{(\mathrm{eff})}=\epsilon\lambda_A\,\kappa_W,
\]
\[
\delta\ln \mu_{W,A}=\epsilon\lambda_A\,\mu_1,
\qquad
\delta\ln T_{U,A}=\epsilon\lambda_A\,\tau_1.
\]
Because the grouped weak-axisymmetric branch keeps the scalar size variables lane-invariant at first order, the axial scale \(L\) is held fixed here. If one later wants an explicit effective axial-size drift, it can be restored by replacing
\[
\Sigma_\delta=\tau_1-\kappa_U
\]
with
\[
\Sigma_\delta=\tau_1-\kappa_U-2\ell_{L,1}.
\]
At the present level we keep \(\ell_{L,1}=0\).

---

## 2. Exact physical-branch drifts from the microscopic variables

Differentiating the Stage-30 coherent ratios gives the physical Stage 181 drift variables directly.

### 2.1 Direct port variables

For the wall-to-mixed overlap ratio,
\[
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\]
we get
\[
\boxed{
\zeta_Z=2\lambda_1-\kappa_\eta-\kappa_W.
}
\]

For the mixed frequency scale,
\[
\Omega_W^2=\frac{K_W^{(\mathrm{eff})}}{\mu_W},
\]
we get
\[
\boxed{
\omega_W=\kappa_W-\mu_1.
}
\]

For the interference ratio,
\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\]
we get
\[
\boxed{
\chi_1=\chi_0\,(\gamma_1+c_1-\kappa_U).
}
\]

For the wall/internal dressing ratio,
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
\]
we get
\[
\boxed{
\eta_1=\epsilon_\eta\,(2c_1-\kappa_U-\kappa_\eta).
}
\]

For the unsplit mixed blocking ratio,
\[
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}},
\]
we get
\[
\boxed{
\varepsilon_W=\epsilon_W\,(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W).
}
\]

For the split-\(U\) axial ratio,
\[
\delta_U=\frac{\pi^2 T_U}{L^2 K_U},
\]
we get
\[
\boxed{
\delta_{U,1}=\delta_U\,(\tau_1-\kappa_U).
}
\]

### 2.2 Support variables are absent

The microscopic support quantities
\[
\lambda_\phi,
\qquad
K_\phi^{(\mathrm{eff})}
\]
enter the coherent branch only through \(\zeta\), and Stage 181 already proved that \(\zeta\) drops out of
\[
\mathcal T^2,
\qquad
R_{\mathrm{target}},
\qquad
\Xi_1.
\]
So the underlying support drifts
\[
\delta\ln\lambda_{\phi,A},
\qquad
\delta\ln K_{\phi,A}^{(\mathrm{eff})}
\]
do not enter the defect at all.

This is the microscopic version of the support-blindness theorem.

---

## 3. Exact microscopic slippage variables

The Stage 181 physical drifts naturally reorganize into the following microscopic slippages:
\[
\boxed{
\Sigma_Z:=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W=\zeta_Z-\omega_W,
}
\]
\[
\boxed{
\Sigma_\chi:=\gamma_1+c_1-\kappa_U=\frac{\chi_1}{\chi_0},
}
\]
\[
\boxed{
\Sigma_\eta:=2c_1-\kappa_U-\kappa_\eta=\frac{\eta_1}{\epsilon_\eta},
}
\]
\[
\boxed{
\Sigma_\epsilon:=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W=\frac{\varepsilon_W}{\epsilon_W},
}
\]
\[
\boxed{
\Sigma_\delta:=\tau_1-\kappa_U=\frac{\delta_{U,1}}{\delta_U}.
}
\]

These are the only microscopic slippages that survive into the grouped defect.

---

## 4. Exact microscopic form of the grouped defect

Stage 181 gave
\[
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}
+\frac{2}{1-\epsilon}
\left[
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-\frac{2\epsilon_W}{11(1+\delta_U)^2}\delta_{U,1}
\right].
\]
Using the microscopic slippages above, this becomes
\[
\boxed{
\Xi_1
=
\Sigma_Z
+\frac{2\chi_0}{1+\chi_0}\,\Sigma_\chi
+\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\,\Sigma_\delta
\right].
}
\]
This is the first exact microscopic slippage law for the coherent branch.

### Interpretation

The grouped weak-axisymmetric defect is carried by exactly four microscopic combinations:

1. wall-normalized overlap/frequency slippage \(\Sigma_Z\),
2. interference/co-loading slippage \(\Sigma_\chi\),
3. mixed-blocking similarity slippage \(\Sigma_\epsilon\),
4. split-\(U\) axial-ratio slippage \(\Sigma_\delta\).

No support-lane microscopic variable appears.

---

## 5. Exact selected-branch demand slippage law

The Stage 180/181 selected-branch identity is
\[
\Xi_1=-\frac{\eta_1}{1-\epsilon_\eta}-\mathcal R_1.
\]
So in microscopic slippage form,
\[
\boxed{
\mathcal R_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta
-\Xi_1.
}
\]
Substituting the microscopic defect law gives
\[
\boxed{
\mathcal R_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta
-\Sigma_Z
-\frac{2\chi_0}{1+\chi_0}\,\Sigma_\chi
-\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\,\Sigma_\delta
\right].
}
\]
So the selected-branch demand ratio is controlled by one additional dressing slippage \(\Sigma_\eta\) beyond the four direct transfer-shape slippages.

---

## 6. Exact tracking-factor slippage and the tracking/nontracking split

Stage 181 already showed that the tracking factor is
\[
R_{\mathrm{tr}}=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}

ahead of the actual grouped defect.
Its logarithmic drift was
\[
\Theta_1
=
-\frac{\chi_0(1+\chi_0)\,\delta_{U,1}+\delta_U(1+\delta_U)\,\chi_1}
{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}.
\]
Substituting the microscopic slippages gives
\[
\boxed{
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}
\Big[(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi\Big].
}
\]
So the natural microscopic tracking slippage is
\[
\boxed{
\Sigma_{\rm tr}:=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi.
}
\]
Then
\[
\boxed{
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr}.
}
\]
Therefore exact tracking rigidity is simply
\[
\boxed{
\Theta_1=0
\iff
\Sigma_{\rm tr}=0
}
\]
on the physical branch \(\chi_0>0\), \(\delta_U>0\).

Now rewrite the grouped defect using
\[
\Sigma_\chi=
\frac{\Sigma_{\rm tr}-(1+\chi_0)\Sigma_\delta}{1+\delta_U}.
\]
Then
\[
\boxed{
\Xi_1
=
\Sigma_Z
+\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)}\,\Sigma_{\rm tr}
+\frac{2\epsilon_W}{1-\epsilon}\,rac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+
\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
}
\]
This is the exact tracking/nontracking split.

### Consequence

If \(\Sigma_{\rm tr}=0\), the grouped defect reduces to
\[
\boxed{
\Xi_1^{(\Theta=0)}
=
\Sigma_Z
+\frac{2\epsilon_W}{1-\epsilon}\,rac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+
\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
}
\]
So exact tracking rigidity is necessary but still not sufficient. The grouped defect survives unless the nontracking microscopic slippages also vanish or cancel.

This is the sharpest exact version of the Stage 181 conclusion.

---

## 7. Zero-defect corollaries

The microscopic zero-defect condition is now completely explicit.

### 7.1 Direct transfer-shape rigidity

A sufficient condition for
\[
\Xi_1=0
\]
is
\[
\boxed{
\Sigma_Z=\Sigma_\chi=\Sigma_\epsilon=\Sigma_\delta=0.
}
\]
That means:

- the wall-normalized overlap/frequency ratio is rigid,
- the interference ratio co-loads with the wall/internal stiffness,
- the mixed blocking similarity is rigid,
- and the split-\(U\) axial ratio is rigid.

### 7.2 Tracking-rigid branch

If the branch is exactly tracking-rigid,
\[
\Sigma_{\rm tr}=0,
\]
then a sufficient condition for zero grouped defect is the smaller set
\[
\boxed{
\Sigma_Z=\Sigma_\epsilon=\Sigma_\delta=0.
}
\]
because \(\Sigma_{\rm tr}=0\) plus \(\Sigma_\delta=0\) already forces \(\Sigma_\chi=0\).

### 7.3 Full selected-branch rigidity

A sufficient condition for both
\[
\Xi_1=0
\qquad\text{and}\qquad
\mathcal R_1=0
\]
is
\[
\boxed{
\Sigma_Z=\Sigma_\chi=\Sigma_\epsilon=\Sigma_\delta=\Sigma_\eta=0.
}
\]
So the selected-branch rigidity problem is now reduced to a finite microscopic slippage ledger.

---

## 8. Best current theorem statement after Stage 182

On the actual coherent local D/N tracking branch:

1. the weak-axisymmetric grouped defect is completely independent of the support microscopic variables \((\lambda_\phi,K_\phi^{(\mathrm{eff})})\);
2. it depends only on the four microscopic mixed/outgoing slippages
   \((\Sigma_Z,\Sigma_\chi,\Sigma_\epsilon,\Sigma_\delta)\);
3. the selected-branch demand ratio introduces one additional dressing slippage \(\Sigma_\eta\);
4. the tracking-factor drift is carried by exactly one combination
   \(\Sigma_{\rm tr}\);
5. and exact tracking rigidity still leaves a nontracking defect unless
   \(\Sigma_Z,\Sigma_\epsilon,\Sigma_\delta\)
   also vanish or cancel.

So the remaining theorem gap is now smaller than Stage 181 left it.

It is no longer merely

> compute the weak-axisymmetric drifts of \((Z_W,\Omega_W^2,\chi_0,\epsilon_W,\delta_U)\).

It is now

> compute the microscopic coherent-kernel slippages
> \((\Sigma_Z,\Sigma_\chi,\Sigma_\epsilon,\Sigma_\delta,\Sigma_\eta)\)
> on the actual moving-throat branch, and test whether the nontracking set vanishes.

That is the direct continuation point.
