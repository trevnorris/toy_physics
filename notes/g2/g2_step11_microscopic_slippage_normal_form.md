# Step 11 — Microscopic slippage decomposition and exact triangular normal form

## Goal

Step 10 proved that the coherent local D/N tracking branch is **support-blind** at the level of the first grouped weak-axisymmetric defect:
```math
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\qquad
\Xi_1=\delta\ln \mathcal T^2,
```
with
```math
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right].
```
So the next honest move is not to revisit coherent support. It is to push the defect one layer deeper, down to the actual **microscopic coherent-kernel slippages**.

That is what this step does.

The main result is that the whole coherent weak-axisymmetric problem collapses to three branch-adapted scalars,
```math
\Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta,
```
rather than a broad raw-placement drift ledger.

---

## Inputs carried forward

From the coherent branch,
```math
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U=\frac{\pi^2 T_U}{L^2K_U},
\qquad
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
```
```math
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}},
\qquad
\frac{Z_W}{\Omega_W^2}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}.
```
The grouped defect from Step 10 was
```math
\Xi_1
=
\Sigma_Z
+\frac{2\chi_0}{1+\chi_0}\,\Sigma_\chi
+\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\,\Sigma_\delta
\right],
```
and the tracking-factor drift came from
```math
R_{\rm tr}=\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}.
```

---

## Step 11A — Direct microscopic log coordinates

Perturb the microscopic coherent-kernel couplings multiplicatively:
```math
\gamma_A=\gamma e^{s\lambda_A\gamma_1},
\quad
c_{\eta U,A}=c_{\eta U}e^{s\lambda_A c_1},
\quad
\lambda_{W,A}=\lambda_W e^{s\lambda_A\lambda_1},
```
```math
K_{U,A}=K_U e^{s\lambda_A\kappa_U},
\quad
K_{\eta,A}=K_\eta e^{s\lambda_A\kappa_\eta},
\quad
K_{W,A}=K_W e^{s\lambda_A\kappa_W},
```
```math
\mu_{W,A}=\mu_W e^{s\lambda_A\mu_1},
\qquad
T_{U,A}=T_U e^{s\lambda_A\tau_1}.
```
Then the direct logarithmic slippages are
```math
\boxed{\Sigma_\chi=\delta\ln\chi_0=\gamma_1+c_1-\kappa_U,}
```
```math
\boxed{\Sigma_\delta=\delta\ln\delta_U=\tau_1-\kappa_U,}
```
```math
\boxed{\Sigma_\eta=\delta\ln\epsilon_\eta=2c_1-\kappa_U-\kappa_\eta,}
```
```math
\boxed{\Sigma_\epsilon=\delta\ln\epsilon_W=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,}
```
```math
\boxed{\Sigma_Z=\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right)=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W.}
```
So the first grouped defect already depends only on five direct microscopic slippages.

---

## Step 11B — Exact tracking combination

The tracking-factor drift depends only on the single combination
```math
\boxed{
\Sigma_{\rm tr}:=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi.
}
```
A direct logarithmic expansion of
```math
R_{\rm tr}=\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
```
gives
```math
\boxed{
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr}.
}
```
So exact tracking rigidity is simply
```math
\Theta_1=0 \iff \Sigma_{\rm tr}=0
```
on the physical branch.

---

## Step 11C — Genuine nontracking transfer-shape slippage

The grouped defect still contains one genuinely nontracking scalar. Define
```math
\boxed{
\Sigma_{\rm nt}
:=
\Sigma_Z
+\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
}
```
Then the coherent grouped defect separates exactly into tracking and nontracking pieces:
```math
\boxed{
\Xi_1
=
A_{\rm tr}\,\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
A_{\rm tr}:=\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)}.
}
```
So once the universal tracking feed-through is subtracted off, the grouped defect is controlled by `\Sigma_{\rm nt}` alone.

---

## Step 11D — Dressing slippage and exact triangular normal form

The selected-branch relation
```math
R_{\rm target}\,\mathcal T^2=\Lambda_0(1-\epsilon_\eta)
```
with scalar front factor `\Lambda_0` inert at grouped weak-axisymmetric order gives
```math
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta,
\qquad
\Sigma_\eta=\delta\ln\epsilon_\eta.
}
```
Putting everything together, the reduced observable drifts take the exact triangular form
```math
\boxed{
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
}
```
```math
\boxed{
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
}
```
```math
\boxed{
\mathcal R_1+\Xi_1=-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta,
}
```
with
```math
C_{\rm tr}
=
\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}.
```
This is the strongest exact compression so far.

The whole coherent weak-axisymmetric problem is no longer a five-slippage bookkeeping problem. It is a **three-scalar normal form**.

---

## Step 11E — Exact inverse reconstruction formulas

Because the normal form is triangular, it can be inverted exactly.

### Tracking slippage from `\Theta_1`
```math
\boxed{
\Sigma_{\rm tr}
=
-\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U}\,\Theta_1.
}
```

### Nontracking slippage from `(\Theta_1,\Xi_1)`
Using `A_{\rm tr}/C_{\rm tr}=2(1+\chi_0+\delta_U)/\delta_U`,
```math
\boxed{
\Sigma_{\rm nt}
=
\Xi_1+\frac{2(1+\chi_0+\delta_U)}{\delta_U}\,\Theta_1.
}
```
So `\Sigma_{\rm nt}` is the grouped defect with the universal tracking feed-through removed.

### Dressing slippage from `(\mathcal R_1,\Xi_1)`
```math
\boxed{
\Sigma_\eta
=
-\frac{1-\epsilon_\eta}{\epsilon_\eta}\,(\mathcal R_1+\Xi_1).
}
```
So the selected-branch residual directly measures the dressing mismatch after the direct transfer-shape defect is added back.

---

## Step 11F — Quartic anomaly gate in the branch-adapted coordinates

The carried quartic anomaly target remains
```math
\Lambda_1\approx 0.279605891931464.
```
So the exact coherent-branch quartic gate is now
```math
\boxed{
A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt}=\Lambda_1.
}
```
This has two immediate specializations.

### If exact tracking rigidity holds
```math
\Sigma_{\rm tr}=0
\quad\Longrightarrow\quad
\boxed{\Sigma_{\rm nt}=\Lambda_1.}
```
So on a tracking-rigid branch the missing `O(f^4)` anomaly layer is purely the **nontracking transfer-shape slippage**.

### If tracking rigidity and selected-branch rigidity both hold
If in addition `\mathcal R_1=0`, then
```math
\boxed{
\Sigma_\eta
=
-\frac{1-\epsilon_\eta}{\epsilon_\eta}\,\Lambda_1.
}
```
So a target-rigid selected branch would require a compensating dressing slippage of precisely that size.

---

## What this changes

This is a real simplification.

Before this step, the remaining g-2 problem still looked like a loose collection of microscopic placement drifts.
After this step, it is exact and minimal:

1. `\Sigma_{\rm tr}` — tracking slippage,
2. `\Sigma_{\rm nt}` — nontracking transfer-shape slippage,
3. `\Sigma_\eta` — dressing slippage.

That is the correct reduced coordinate system for the next moving-throat test.

---

## Continuation point

The next clean move is now very sharp:

> determine whether the actual moving-throat branch preserves or drives the three branch-adapted scalars `(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)`, and then push them one step further down to exact branch composites / microscopic monomials.
