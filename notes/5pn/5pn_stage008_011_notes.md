
# 5PN continuation notes — stages 8 through 11

These stages continue the weak-axisymmetric grouped-`P2` program from the Stage-7 scalar amplitude
\[
\Xi_1=\frac{P_1}{P_0}.
\]

The overall effect is that the remaining first-order grouped normalization problem is no longer a generic “compute lots of microscopic drifts” problem. It has collapsed to a small sequence of exact equivalent formulations.

## Stage 8 — direct outgoing-port co-loading

The remaining grouped weak-axisymmetric defect can be written directly in terms of the actual outgoing-port slopes:
\[
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1)
=
\bar\nu_N-\kappa_1,
\]
where
\[
\nu_r = 2(\mathfrak p_r-\mathfrak d_r)
\]
is the logarithmic slope of
\[
N_{0}^{(r)}=\frac{P_r^2}{\Delta_r^2}.
\]

Equivalently, if
\[
N_{A,0}^{(r)} = K_A \mathcal T_{A,r}^2,
\qquad
\delta\ln \mathcal T_{A,r}=\epsilon\lambda_A\,\tau_r,
\]
then
\[
\nu_r=\kappa_1+2\tau_r,
\qquad
\Xi_1 = 2\sum_r \rho_r^{(N)}\tau_r.
\]

So the exact zero-defect condition is now:
\[
\sum_r \rho_r^{(N)}\tau_r=0.
\]

A stronger sufficient condition is per-port co-loading:
\[
\tau_r=0
\qquad
\text{for every active outgoing port }r.
\]

Under upper-leg and coupling rigidity, the transfer-shape slope collapses to the raw mixed-leg slope and recovers the old square-root mixed-leg law.

## Stage 9 — one effective transfer shape and the actual continuum branch

The many-port weighted sum collapses exactly to one effective transfer shape:
\[
\mathcal T_{\mathrm{eff},A}^2
=
\sum_r \mathcal T_{A,r}^2
=
\frac{N_{A,0}}{K_A},
\qquad
\Xi_1
=
\frac{\delta\ln\mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A}.
\]

On the actual one-port continuum branch,
\[
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2},
\]
so
\[
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\rho_1}{1+\rho}+\frac{2\varepsilon_W}{1-\epsilon_W}.
\]

In selected-branch form,
\[
\mathcal T_A^2
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}\,
\frac{1-\epsilon_{\eta,A}}{R_{\mathrm{target},A}},
\]
which yields
\[
\Xi_1
=
-\frac{\eta_1}{1-\epsilon_\eta}
-\mathcal R_1.
\]

On the coherent local D/N branch,
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\qquad
\epsilon
=
\epsilon_W\!\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
and the defect becomes
\[
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}
+\frac{2\epsilon_1}{1-\epsilon}.
\]

The support ratio drops out identically:
\[
\partial_\zeta \ln \mathcal T^2 = 0.
\]

So the coherent defect is support-blind, and exact tracking rigidity by itself is not enough to kill it.

## Stage 10 — microscopic coherent-kernel slippages and exact triangular normal form

The coherent branch depends only on the microscopic slippages
\[
\Sigma_Z,\quad
\Sigma_\chi,\quad
\Sigma_\epsilon,\quad
\Sigma_\delta,
\]
with one additional dressing slippage
\[
\Sigma_\eta
\]
entering the selected-branch form.

The exact microscopic grouped-defect law is
\[
\Xi_1
=
\Sigma_Z
+\frac{2\chi_0}{1+\chi_0}\Sigma_\chi
+\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\Sigma_\delta
\right].
\]

The exact tracking combination is
\[
\Sigma_{\rm tr}
=
(1+\chi_0)\Sigma_\delta
+
(1+\delta_U)\Sigma_\chi,
\]
with
\[
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr}.
\]

Defining the genuine nontracking slippage
\[
\Sigma_{\rm nt},
\]
the coherent problem takes the exact triangular normal form
\[
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta.
\]

So on the constructive coherent branch,
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0.
\]

This is the first exact three-scalar normal form of the full coherent weak-axisymmetric problem.

## Stage 11 — direct microscopic monomials, similarity orbit, and quotient closure

The three branch-adapted coordinates can be written as logarithmic drifts of three exact microscopic monomials:
\[
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln \epsilon_\eta=\Sigma_\eta,
\]
with
\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
\]
\[
\mathfrak C_{{\rm nt},*}
=
\left(\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}\right)
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}\right)^{E_*}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{-F_*},
\]
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}.
\]

The exact monomial-drift map is the rank-3 matrix
\[
M_*\,\delta\mathbf x
=
\begin{pmatrix}
\delta\ln \mathfrak C_{{\rm tr},*}\\
\delta\ln \mathfrak C_{{\rm nt},*}\\
\delta\ln \epsilon_\eta
\end{pmatrix},
\]
with
\[
\det M_*^{(\tau,\kappa_\eta,\mu_1)}=1+\chi_{0,*}>0,
\]
so
\[
\dim\ker M_* = 5.
\]

There is an exact five-parameter multiplicative similarity orbit \(\mathcal G_*\) preserving the three monomials exactly, and the scripts show
\[
\ker M_* = T_{\mathrm{id}}\mathcal G_*.
\]

More strongly, the exact finite invariant-fibre equations are
\[
M_*\,\Delta \mathbf x = 0,
\]
and solving them reproduces the exact orbit exponents. Therefore the finite level sets of
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)
\]
are precisely the similarity orbits \(\mathcal G_*\).

So the coherent weak-axisymmetric zero-defect theorem can now be written as
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x \in T_{\mathrm{id}}\mathcal G_*,
\]
and, at finite level,
\[
\mathcal M_+/\mathcal G_*
\cong
(\mathbb R_{>0})^3
\]
with quotient coordinates
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
\]

## What this means for the 5PN program

At this point the weak-axisymmetric normalization problem is no longer an algebraic bottleneck.

There are now three equivalent ways to state the next theorem gate:

1. compute the direct grouped scalar
   \[
   \Xi_1=\frac{P_1}{P_0},
   \]
2. compute the branch-adapted defect coordinates
   \[
   (\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta),
   \]
3. or test whether the actual moving-throat weak-axisymmetric branch is tangent to the exact monomial-preserving similarity orbit \(\mathcal G_*\).

That last form is the sharpest current continuation point.
