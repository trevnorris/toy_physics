# 5PN continuation notes — Stages 309–311

These stages continue directly from the Stage-306–308 physical branch observable compiler.
The key question there was whether the actual physical coherent branch can satisfy the nontrivial orbit-lock condition
\[
 d\ln \Omega_W^2-d\ln Z_W-
 \frac{2\chi_0}{1+\chi_0}d\ln\chi_0-
 \frac{2\epsilon}{1-\epsilon}d\epsilon = 0
\]
on the same branch where the support baseline is strong enough.

The answer is sharper than a generic “maybe.”

1. The first-order grouped weak-axisymmetric defect is carried only by the **mixed/outgoing placement variables**, not by the coherent support-enhancement baseline.
2. Exact tracking rigidity is **necessary but not sufficient**.
3. The full coherent first-order problem collapses to an exact **triangular three-scalar normal form**.
4. Equivalently, first-order zero defect is just invariance of three exact branch composites.

This is exactly the direction the later moving-throat notes had already pointed toward: the remaining weak-axisymmetric theorem gap is not support amplitude, but the drift of the physical mixed/outgoing placement variables and the corresponding exact branch composites. fileciteturn80file6turn80file13

## Stage 309 — physical tracking-branch transfer-shape defect law

On the coherent local D/N tracking branch,
\[
R_{\rm tr}=\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\qquad
\epsilon
=\epsilon_W\!\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
\[
M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\]
\[
R_{\rm target}
=\Lambda_0\,\Omega_W^2\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\qquad
\Lambda_0=\frac{27\pi^2Gc_s^5}{20a^5c^5}.
\]
The exact transfer shape is therefore
\[
\mathcal T^2
=\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}
=\Lambda_0\,\frac{1-\epsilon_\eta}{R_{\rm target}}.
\]

The first exact theorem is support-blindness:
\[
\partial_\zeta \ln \mathcal T^2=0,
\qquad
\partial_\zeta \ln R_{\rm target}=0,
\qquad
\partial_\zeta \Xi_1=0.
\]
So coherent support can raise the baseline through \(M_{\rm tr}\), but it cannot repair or spoil the first-order grouped defect.

Writing the weak-axisymmetric drifts as
\[
\delta\ln Z_W=\epsilon\lambda_A\zeta_Z,
\qquad
\delta\ln \Omega_W^2=\epsilon\lambda_A\omega_W,
\]
\[
\delta\chi_0=\epsilon\lambda_A\chi_1,
\qquad
\delta\epsilon=\epsilon\lambda_A\epsilon_1,
\qquad
\delta\epsilon_\eta=\epsilon\lambda_A\eta_1,
\]
one gets the exact physical defect laws
\[
\Xi_1
=\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}+\frac{2\epsilon_1}{1-\epsilon},
\]
\[
\mathcal R_1
=\omega_W-\zeta_Z-\frac{2\chi_1}{1+\chi_0}-\frac{2\epsilon_1}{1-\epsilon}-\frac{\eta_1}{1-\epsilon_\eta},
\]
so that
\[
\mathcal R_1+\Xi_1=-\frac{\eta_1}{1-\epsilon_\eta}.
\]

The exact split-blocking drift is
\[
\epsilon_1
=
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\epsilon_{W,1}
-\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}.
\]

## Stage 310 — microscopic coherent-kernel slippages

The Stage-309 physical drift variables are exact functions of the microscopic coherent-kernel drifts
\[
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1).
\]
The exact port-side translations are
\[
\zeta_Z=2\lambda_1-\kappa_\eta-\kappa_W,
\qquad
\omega_W=\kappa_W-\mu_1,
\]
\[
\chi_1=\chi_0(\gamma_1+c_1-\kappa_U),
\qquad
\eta_1=\epsilon_\eta(2c_1-\kappa_U-\kappa_\eta),
\]
\[
\epsilon_{W,1}=\epsilon_W(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W),
\qquad
\delta_{U,1}=\delta_U(\tau_1-\kappa_U).
\]

This motivates the exact microscopic slippage variables
\[
\Sigma_Z:=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W=\zeta_Z-\omega_W,
\]
\[
\Sigma_\chi:=\gamma_1+c_1-\kappa_U=\frac{\chi_1}{\chi_0},
\qquad
\Sigma_\eta:=2c_1-\kappa_U-\kappa_\eta=\frac{\eta_1}{\epsilon_\eta},
\]
\[
\Sigma_\epsilon:=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W=\frac{\epsilon_{W,1}}{\epsilon_W},
\qquad
\Sigma_\delta:=\tau_1-\kappa_U=\frac{\delta_{U,1}}{\delta_U}.
\]
Then the grouped weak-axisymmetric defect becomes
\[
\Xi_1
=\Sigma_Z
+\frac{2\chi_0}{1+\chi_0}\Sigma_\chi
+\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\Sigma_\delta
\right].
\]

The exact tracking combination is already isolated:
\[
\Sigma_{\rm tr}:=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi,
\]
with
\[
\Theta_1
=-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr}.
\]
So exact tracking rigidity kills \(\Theta_1\), but it does **not** kill \(\Xi_1\) unless the remaining nontracking slippages also cooperate. That is the explicit microscopic version of the note-side warning that exact tracking-factor rigidity is not sufficient by itself. fileciteturn80file6

## Stage 311 — exact triangular normal form and direct branch composites

The exact branch-adapted nontracking slippage is
\[
\Sigma_{\rm nt}:=
\Sigma_Z
+\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon
-
\left[
\frac{2\chi_0}{1+\delta_U}
+\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
\]
Then the three observables collapse to the exact triangular system
\[
\Theta_1=-C_{\rm tr}\,\Sigma_{\rm tr},
\qquad
C_{\rm tr}:=\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)},
\]
\[
\Xi_1=A_{\rm tr}\,\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
A_{\rm tr}:=\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)},
\]
\[
\mathcal R_1+\Xi_1=-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta.
\]
This inverts exactly to
\[
\Sigma_{\rm tr}
=-\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U}\,\Theta_1,
\]
\[
\Sigma_{\rm nt}
=\Xi_1+\frac{2(1+\chi_0+\delta_U)}{\delta_U}\,\Theta_1,
\]
\[
\Sigma_\eta
=-\frac{1-\epsilon_\eta}{\epsilon_\eta}(\mathcal R_1+\Xi_1).
\]

So the coherent first-order problem is no longer a five-slippage bookkeeping problem. It is an exact three-scalar normal form.

There is also an exact direct branch-composite theorem. Define
\[
C_*:=\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U},
\qquad
B_*:=\frac{2(1+\chi_0+\delta_U)}{\delta_U},
\]
\[
\mathfrak T_*:=R_{\rm tr}^{-C_*},
\qquad
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*}.
\]
Then their logarithmic drifts are exactly
\[
\delta\ln \mathfrak T_* = \Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak N_* = \Sigma_{\rm nt},
\qquad
\delta\ln \epsilon_\eta = \Sigma_\eta.
\]
Therefore
\[
\Theta_1=\Xi_1=\mathcal R_1+\Xi_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0
\iff
\delta\ln R_{\rm tr}=0,
\ \delta\ln(\mathcal T^2 R_{\rm tr}^{B_*})=0,
\ \delta\ln\epsilon_\eta=0.
\]
So the actual first-order coherent theorem gate is now exactly the invariance of three branch composites, not the adjustment of the support baseline. That matches the compact continuation point already recorded in the moving-throat ledger. fileciteturn80file13turn80file15turn80file11

## Bottom line after Stages 309–311

The program did **not** hit a dead end here.
It hit a sharper separation of roles.

- The coherent support lane can still rescue the **steady normalization baseline**.
- But it is exactly absent from the **first-order grouped weak-axisymmetric defect**.
- The first-order defect is carried only by the mixed/outgoing placement variables and their exact branch composites.

So the next honest theorem gate is now very small:

> compute the actual first grouped weak-axisymmetric drifts of
> \[
> R_{\rm tr},
> \qquad
> \mathcal T^2 R_{\rm tr}^{B_*},
> \qquad
> \epsilon_\eta
> \]
> on the physical moving-throat branch.

If all three are invariant, the coherent first-order grouped defect vanishes automatically. If not, the failure is localized immediately into tracking, nontracking transfer shape, or dressing. fileciteturn80file13turn80file15
