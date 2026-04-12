# 5PN continuation notes — Stages 296–299

This continuation takes the Stage-295 collapse
\[
u_2^{(1)},\qquad \frac{P_1}{P_0},
\]
and pushes it all the way down to the actual coherent moving-throat branch variables and then to the exact microscopic monomial / similarity-orbit structure.

A deliberate firewall was kept throughout:

- **do not modify the parent GNLS bulk PDE,**
- continue only on the moving-throat / boundary-response side,
- and keep the live 5PN / 2.5PN / 4PN obstruction on the actual grouped branch rather than inventing a new bulk medium.

## Stage 296 — actual-branch static slope compiler

Files:
- `5pn_stage296_actual_branch_static_slope_compiler.py`
- `5pn_stage296_actual_branch_static_slope_compiler_output.txt`

Main result:

The physical weak-axisymmetric continuation point is now compiled directly in the sharpest static grouped form:
\[
D_{01}=K_1-B_0^{(1)}-Z_0^{(1)},
\qquad
D_{21}=-(M_1+B_2^{(1)}+Z_2^{(1)}),
\]
\[
u_2^{(1)}=-\frac{D_{21}+u_2D_{01}}{D_0},
\qquad
\frac{P_1}{P_0}=\frac{N_1}{N_0}-\frac{D_{01}}{D_0}.
\]

On the canonical even-preserving branch,
\[
D_{21}=-\frac{D_{01}}{9},
\qquad
u_2=\frac19,
\]
so the conservative slope vanishes exactly,
\[
u_2^{(1)}=0,
\]
and the whole remaining grouped defect becomes
\[
\Xi_1=\frac{P_1}{P_0}=\bar\nu_N-\kappa_1,
\]
with
\[
\bar\nu_N=\sum_r w_r \nu_r,
\qquad
\nu_r=\delta\ln N_{A,0}^{(r)}/(\epsilon\lambda_A).
\]

So the actual static continuation point is no longer “all grouped perturbations.”
It is

1. one conservative static slope `D01/D0`,
2. one outgoing-transfer static slope `N1/N0`,
3. and, on the even-preserving branch, only their difference `Xi_1`.

## Stage 297 — coherent tracking-branch transfer shape and support-blindness

Files:
- `5pn_stage297_coherent_tracking_branch_transfer_shape.py`
- `5pn_stage297_coherent_tracking_branch_transfer_shape_output.txt`

Main result:

On the actual coherent local D/N tracking branch, the exact transfer shape is
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\qquad
\epsilon
=
\epsilon_W\!\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
with selected-branch identity
\[
R_{\mathrm{target}}\mathcal T^2=\Lambda_0(1-\epsilon_\eta).
\]

The exact weak-axisymmetric grouped defect is
\[
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}+\frac{2\epsilon_1}{1-\epsilon},
\]
where
\[
\epsilon_1
=
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}.
\]

The strongest theorem from this stage is the exact support-blindness result:
\[
\partial_\zeta \ln \mathcal T^2
=
\partial_\zeta \ln R_{\mathrm{target}}
=
\partial_\zeta \Xi_1
=
0.
\]

So coherent support can move the steady normalization baseline, but it cannot repair or spoil the first weak-axisymmetric grouped defect at all. That defect is carried only by the mixed/outgoing placement variables.

## Stage 298 — microscopic coherent-kernel slippage normal form

Files:
- `5pn_stage298_microscopic_coherent_slippage_normal_form.py`
- `5pn_stage298_microscopic_coherent_slippage_normal_form_output.txt`

Main result:

The coherent weak-axisymmetric defect depends only on the microscopic slippages
\[
\Sigma_Z,\qquad
\Sigma_\chi,\qquad
\Sigma_\eta,\qquad
\Sigma_\epsilon,\qquad
\Sigma_\delta,
\]
with
\[
\Sigma_Z=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W,
\]
\[
\Sigma_\chi=\gamma_1+c_1-\kappa_U,
\qquad
\Sigma_\eta=2c_1-\kappa_U-\kappa_\eta,
\]
\[
\Sigma_\epsilon=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,
\qquad
\Sigma_\delta=\tau_1-\kappa_U.
\]

Then the exact branch-adapted coordinates are
\[
\Sigma_{\rm tr}=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi,
\]
\[
\Sigma_{\rm nt}
=
\Sigma_Z
+\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-
\left[
\frac{2\chi_0}{1+\delta_U}
+
\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
\]

The observable drift ledger becomes the exact triangular normal form
\[
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta.
\]

So the full first weak-axisymmetric grouped normalization problem is no longer a five-slippage bookkeeping problem. It is exactly the vanishing of three branch-adapted scalars:
\[
\Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta.
\]

The inverse reconstruction is exact:
\[
\Sigma_{\rm tr}
=
-\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U}\,\Theta_1,
\]
\[
\Sigma_{\rm nt}
=
\Xi_1+\frac{2(1+\chi_0+\delta_U)}{\delta_U}\,\Theta_1,
\]
\[
\Sigma_\eta
=
-\frac{1-\epsilon_\eta}{\epsilon_\eta}(\mathcal R_1+\Xi_1).
\]

So the sharp zero-defect statement is now
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0.
\]

## Stage 299 — direct microscopic monomials and similarity orbit

Files:
- `5pn_stage299_microscopic_monomial_orbit_bridge.py`
- `5pn_stage299_microscopic_monomial_orbit_bridge_output.txt`
- `fivepn_stage296_299_common.py`

Main result:

The three branch-adapted defect coordinates are exactly the first logarithmic drifts of three direct microscopic monomials:
\[
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln \epsilon_\eta=\Sigma_\eta.
\]

The exact direct monomials are
\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2 T_U}{L^2 K_U}\right)^{1+\chi_{0,*}},
\]
\[
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*},
\]
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}.
\]

The exact monomial-drift map on the microscopic grouped drift vector
\[
\delta\mathbf x
=
(\delta\ln\lambda_W,\delta\ln c_{\eta U},\delta\ln\gamma,\delta\ln K_U,
\delta\ln K_\eta^{(\mathrm{eff})},\delta\ln K_W^{(\mathrm{eff})},
\delta\ln\mu_W,\delta\ln T_U)
\]
is the rank-3 matrix
\[
M_*
\delta\mathbf x
=
\begin{pmatrix}
\delta\ln \mathfrak C_{{\rm tr},*}\\
\delta\ln \mathfrak C_{{\rm nt},*}\\
\delta\ln \epsilon_\eta
\end{pmatrix},
\]
with exact minor
\[
\det M_*^{(\delta\ln T_U,\delta\ln K_\eta^{(\mathrm{eff})},\delta\ln\mu_W)}
=
1+\chi_{0,*}>0.
\]

So
\[
\mathrm{rank}(M_*)=3,
\qquad
\dim\ker M_* = 5.
\]

This is the exact tangent space of a five-parameter microscopic similarity orbit preserving
\[
\mathfrak C_{{\rm tr},*},\qquad
\mathfrak C_{{\rm nt},*},\qquad
\epsilon_\eta
\]
exactly.

So the final reduced weak-axisymmetric closure theorem is now:
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x \in \ker M_*.
\]

Equivalently, on the coherent local D/N tracking branch, the actual grouped weak-axisymmetric moving-throat tangent must lie inside the exact five-dimensional monomial-preserving similarity orbit.

## What changes after Stage 299

The continuation point is now much smaller than when we started this pass.

It is no longer:

> compute all grouped-lane perturbations of
> \(
> K_A,\ M_A,\ c_{A\alpha},\ \varpi_{A\alpha},\ \Omega_{U/W,A,r},\ R_{A,r},\ g_{U/W,A,r}
> \)
> and then somehow infer \(u_2^{(1)}\) and \(P_1/P_0\).

It is now:

1. compile the actual branch static slopes into
   \[
   \kappa_1=\frac{D_{01}}{D_0},
   \qquad
   \bar\nu_N=\frac{N_1}{N_0},
   \qquad
   \Xi_1=\bar\nu_N-\kappa_1,
   \]
   on the even-preserving branch;
2. compute the three direct weak-axisymmetric monomial drifts
   \[
   \delta\ln \mathfrak C_{{\rm tr},*},\qquad
   \delta\ln \mathfrak C_{{\rm nt},*},\qquad
   \delta\ln \epsilon_\eta,
   \]
   on the actual coherent moving-throat branch;
3. test whether the physical branch tangent lies in the exact five-dimensional similarity kernel.

So the next real theorem gate is not more algebra. It is an **actual branch-tangent extraction problem** on the moving-throat PDE side.
