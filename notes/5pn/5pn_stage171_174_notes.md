# 5PN Stages 171–174 — Adiabatic Wall, Elastic Branch Selection, and Orbit Locking

These stages implement the requested bridge from the `$g-2$` track back into the moving-throat / 5PN branch-selection problem.

## Stage 171 — Adiabatic wall bundle transport

Impose the adiabatic wall constraint
\[
\delta\ln\Theta_w=0.
\]
Using the exact Stage-149/150 inversion and bundle transport laws, the isotropic drift family collapses to
\[
\delta\ln\rho_w=\delta\ln c_{s,w}=\delta\ln\ell=0,
\]
\[
\delta\ln a=\delta\ln L_W=\frac12\,\delta\ln K_s,
\]
\[
\delta\ln c_s=\frac12\,\delta\ln K_s+\frac15\,\delta\ln P_0,
\]
\[
\delta\ln\mathcal Z_q=\delta\ln K_q-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln v_{w0}=-\frac34\,\delta\ln K_s+\frac12\,\delta\ln K_q,
\]
\[
\delta\ln\mathcal T_m=-\frac54\,\delta\ln K_s+\frac12\,\delta\ln K_q-\frac25\,\delta\ln P_0.
\]
The parent invariants remain tangent-compensated:
\[
\delta\ln r_c=\delta\ln\mathfrak r=\delta\ln\mathfrak g=0.
\]
So the adiabatic wall removes the wall-depth / thermal-fraying isotropic drift channel, but leaves a 3-parameter isotropic family labelled by
\[
(\delta\ln K_s,\ \delta\ln K_q,\ \delta\ln P_0).
\]

## Stage 172 — Adiabatic-elastic slippage collapse

The scalar off-bundle slippages are
\[
\varepsilon_L=\delta\ln L_W-\frac12\,\delta\ln K_s,
\]
\[
\varepsilon_v=\delta\ln v_{w0}+\frac34\,\delta\ln K_s-\frac12\,\delta\ln K_q,
\]
\[
\varepsilon_T=\delta\ln\mathcal T_m+\frac54\,\delta\ln K_s-\frac12\,\delta\ln K_q+\frac25\,\delta\ln P_0.
\]
On the adiabatic wall branch these need not vanish automatically, but if we impose the stronger elastic/no-fraying rule
\[
\varepsilon_L=\varepsilon_v=\varepsilon_T=0,
\]
then the scalar normal defect vanishes identically:
\[
\delta_\perp=-\varepsilon_\perp=0.
\]
So the adiabatic-elastic boundary condition kills the first scalar off-bundle obstruction completely.

## Stage 173 — Why adiabatic wall alone is not enough for orbit locking

Stage-169/170 says the reduced weak-axisymmetric defect vanishes iff the microscopic grouped drift is tangent to the exact similarity orbit $\mathcal G_*$, equivalently iff the three quotient coordinates are preserved:
\[
\delta\ln\mathfrak C_{{\rm tr},*}=0,
\qquad
\delta\ln\mathfrak C_{{\rm nt},*}=0,
\qquad
\delta\ln\epsilon_\eta=0.
\]
These are encoded by the rank-3 monomial-drift map
\[
M_*\,\delta\mathbf x.
\]
Stage 173 shows explicitly that `\delta\ln\Theta_w=0` by itself does **not** imply
\[
M_*\delta\mathbf x=0.
\]
There are microscopic drift directions, such as a pure `\Delta\lambda_W` or pure `\Delta c_{\eta U}` perturbation, that leave the wall-depth condition untouched but still move the quotient coordinates. So the adiabatic-wall condition removes one failure channel, but does not by itself prove orbit locking.

## Stage 174 — Adiabatic-elastic orbit theorem

Combining the Stage-172 scalar result with the Stage-166 and Stage-169/170 quotient-closure theorem gives the clean unified statement:

If we impose
\[
\delta\ln\Theta_w=0,
\qquad
\varepsilon_L=\varepsilon_v=\varepsilon_T=0,
\]
then the scalar off-bundle source vanishes, and the remaining first-order defect is zero **iff** the branch preserves the three quotient coordinates:
\[
\delta\ln\mathfrak C_{{\rm tr},*}=0,
\qquad
\delta\ln\mathfrak C_{{\rm nt},*}=0,
\qquad
\delta\ln\epsilon_\eta=0.
\]
Equivalently,
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
M_*\delta\mathbf x=0
\iff
\delta\mathbf x\in T_{\mathrm{id}}\mathcal G_*
\iff
\text{the branch stays on a single exact }\mathcal G_*\text{ orbit.}
\]

So the requested “unified loop” result is now explicit:

- the adiabatic wall condition freezes the thermal wall channel,
- the elastic/no-fraying condition removes the scalar off-bundle obstruction,
- and the remaining branch-selection test is **exactly** whether the moving-throat branch preserves
  \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\), i.e. stays on one \(\mathcal G_*\) orbit.

## Direct continuation point

The next clean theorem gate is no longer to manipulate isotropic bundle transport. That part is closed. The next gate is to compute the actual physical-branch drift vector and test whether its projection under `M_*` vanishes. In other words:
\[
M_*\delta\mathbf x\stackrel{?}=0.
\]
If yes, the adiabatic-elastic moving-throat branch is orbit-locked. If not, the failure is localized immediately into the tracking, nontracking-transfer, or dressing quotient directions.
