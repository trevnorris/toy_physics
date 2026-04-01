# Moving-Throat Translation Dictionary
## From 4D fields to reduced macroscopic variables, microscopic variables, and the three monomial invariants

## 0. Why this document exists

This document is the **ledger** handoff. It is the compact dictionary that turns the raw 4D/moving-throat PDE variables into the reduced response variables actually used in the later derivations.

It is meant to let a fresh session answer questions of the form:

- what exactly is a macroscopic variable here?
- what are the microscopic kernel variables?
- how do the grouped `P2` response coefficients arise?
- what are the direct coherent-branch defect coordinates?
- and what are the **three exact monomial invariants** that now carry the final closure theorem?

This document is structured from coarse to fine:

1. core 4D variables,
2. moving-throat and grouped-`P2` macroscopic/reduced variables,
3. coherent-branch microscopic variables,
4. microscopic slippages,
5. branch-adapted defect coordinates,
6. the three monomial invariants,
7. the similarity orbit / quotient theorem.

The intended use is simple: if a new session can read **this ledger + the PDE engine document**, it should be able to reconstruct the current theorem status without needing the conversational history.

---

## 1. Notation firewall and reading rules

### 1.1 Exact vs reduced vs open

- **Exact**: follows directly from the declared action, exact definitions, or exact algebra.
- **Controlled reduction**: follows only after a stated ansatz or reduction.
- **Protocol closure**: fixed only inside the declared response hierarchy.
- **Open**: still depends on the completed moving-throat PDE.

### 1.2 Notation firewall

The following notational separations are non-negotiable.

1. Electric charge is carried by
   \[
   \eta_Q,\ q_\star,\ q_{\rm eff},
   \]
   not by circulation.
2. The historical gravity-side bare `q=1` is really
   \[
   \kappa_\rho=1,
   \]
   not electric charge.
3. Grouped labels `20/21/22` refer to grouped real `P2` lanes, not spacetime indices.
4. The mixed channels
   \[
   A_w,\ J^w,\ F_{\mu w},\ E_w,\ C_a
   \]
   are suppressed only in the strict far-field brane reduction. They remain microscopic degrees of freedom.

---

## 2. Core 4D variables

### 2.1 Coordinates and fields

\[
x^M=(t,x,y,z,w),\qquad
\mathbf X=(x,y,z,w),\qquad
\mathbf x=(x,y,z).
\]

\[
\psi(\mathbf X,t),\qquad \rho=|\psi|^2,\qquad
A_M=(A_0,A_i),\qquad F_{MN}=\partial_MA_N-\partial_NA_M.
\]

### 2.2 Charge variables

\[
\eta_Q=\pm 1,\qquad
q_\star=\eta_Q e_\star,\qquad
e_\star>0.
\]

After zero-mode canonical normalization,
\[
q_{\rm eff}=\frac{q_\star}{\sqrt{Z_{\rm int}}},
\qquad
e_{\rm eff}=\frac{e_\star}{\sqrt{Z_{\rm int}}},
\qquad
Z_{\rm int}=\int Z(w)\,dw.
\]

### 2.3 Old collective throat variables

\[
a(t),\qquad L(t).
\]

These are **collective moments** of the moving-throat shape field, not the fundamental geometry variables.

---

## 3. Moving-throat geometry variables

### 3.1 Level-set and shape-field representation

\[
\Sigma(\mathbf X,t)=r-R(\Omega,w,t),
\qquad
r=\sqrt{x^2+y^2+z^2},
\qquad
\Omega=\mathbf x/r\in S^2.
\]

The throat surface is \(\Sigma=0\).

The stationary reference throat is
\[
\Sigma_0(\mathbf X)=r-R_0(w).
\]

### 3.2 Wall displacement

\[
R(\Omega,w,t)=R_0(w)+\eta(\Omega,w,t).
\]

### 3.3 Harmonic decomposition

\[
\eta(\Omega,w,t)
=
\eta_0(w,t)Y_{00}(\Omega)
+\sum_{m\in P_2({\rm real})}q_{2m}(w,t)\,Y_{2m}^{\rm real}(\Omega)
+\eta_{\ge 3}.
\]

The grouped real `P2` set is
\[
\{20,\ 21c,\ 21s,\ 22c,\ 22s\}.
\]

### 3.4 Monopole normalization bridge

With
\[
Y_{00}=\frac{1}{2\sqrt\pi},
\]
the physical mouth-average shift \(\delta a\) and the normalized monopole coefficient satisfy
\[
\boxed{
q_{00}(0,t)=2\sqrt\pi\,\delta a(t).
}
\]

---

## 4. Brane/macroscopic variables

### 4.1 Projected brane observables

\[
\rho_{\rm brane}=\int W(w)\rho(\mathbf x,w,t)\,dw,
\qquad
\mathbf j_{\rm brane}=\int W(w)\mathbf j_{xyz}(\mathbf x,w,t)\,dw,
\]
\[
\mathbf v_{\rm brane}=\mathbf j_{\rm brane}/\rho_{\rm brane}.
\]

### 4.2 Leakage

\[
S_{\rm leak}
=
-\left[Wj^w\right]_{-\infty}^{+\infty}
+\int W'(w)j^w\,dw.
\]

### 4.3 Brane velocity potential and Poisson hook

\[
\mathbf v_{\rm brane}=\nabla_3\varphi+\mathbf v_T,\qquad \nabla_3\cdot\mathbf v_T=0.
\]

In the quasi-static regime this gives the Poisson hook for \(\varphi\).

---

## 5. Family-1 / core–mouth macroscopic variables

These variables belong to the earlier core–mouth / compensated-branch part of the moving-throat program. They are not the final monomial invariants, but they remain part of the current dictionary.

### 5.1 Support/mouth source and stiffness variables

\[
K_s,\qquad K_q,\qquad \lambda,\qquad g_s,\qquad g_q.
\]

On the explicit throat-core branch,
\[
K_s=4\pi a^2\!\left(\frac{H_w\ell}{3}+\frac{\hbar^2}{15m_\psi\rho_w\,\ell}\right),
\]
\[
K_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi^2 c_s^2}{4L_W^2},
\]
\[
\lambda=-q_\star v_{w0}\mathcal I_{sq},
\]
\[
g_s=\mathcal T_m\frac{4\pi a^2\ell}{3},
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi}{\sqrt2\,L_W^{3/2}}.
\]

### 5.2 Dimensionless Family-1 ratios

\[
\boxed{
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}.
}
\]

### 5.3 Traction and source moments

The canonical mouth-traction variable is
\[
\boxed{
\Sigma_0=\frac{20}{9}\,\widehat T_m^2.
}
\]

The mixed loading ratio is
\[
\boxed{
\mathcal R=\frac{(\mathfrak g-\mathfrak r)^2}{1+\mathfrak r^2}.
}
\]

The support-shape functional is
\[
\mathcal S[\Sigma]=\int_0^1 K_q(x)\,\Sigma(x)\,dx,
\]
with the reference quadrupole kernel used in the fixed-point audits
\[
K_q(x)=\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}.
\]

The mouth slope variable is
\[
\boxed{
\Pi=\Sigma_0\,[1-\mathcal R\,\mathcal S].
}
\]

The mouth source moments are
\[
\boxed{
M_s=\Sigma_0,
\qquad
M_q=-\Sigma_0\,\mathcal R.
}
\]

### 5.4 Useful lower-branch exact identities

The parent compensation condition is
\[
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2.
\]

On the lower compensated branch,
\[
\mathfrak g_-(\mathfrak r)=\mathfrak r-\frac12\sqrt{1+\mathfrak r^2}.
\]

The first off-family normal coordinate is
\[
\delta_\perp
=
\delta\mathfrak g
-
\mathfrak g'_-(\mathfrak r_\ast)\,\delta\mathfrak r.
\]

This later collapses to explicit logarithmic microscopic imbalance channels.

---

## 6. Microscopic variables behind the core–mouth branch

### 6.1 Throat-core microscopic variables

\[
a,\qquad L_W,\qquad \rho_w,\qquad c_{s,w},\qquad c_s,\qquad \mathcal Z_q,\qquad \mathcal T_m,\qquad v_{w0}.
\]

### 6.2 Healing-lock shell variables

\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\qquad
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\,\ell}
\]
on the carried healing-locked shell branch.

### 6.3 Exact lower-branch drift laws

On the exact lower compensated branch,
\[
\delta\ln L_W=\delta\ln a,
\]
\[
\delta\ln v_{w0}
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
+\delta\ln c_s
-\frac52\,\delta\ln a,
\]
\[
\delta\ln \mathcal T_m
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
-\delta\ln c_s
-\frac32\,\delta\ln a.
\]

These are exact reduced lower-branch transport identities.

---

## 7. Grouped real `P2` macroscopic response variables

For each grouped lane \(A\in\{20,21,22\}\), define the conservative operator
\[
D_A^{\rm(cons)}(\omega)=D_{A,0}+D_{A,2}\omega^2+D_{A,4}\omega^4+O(\omega^6).
\]

### 7.1 Normalized grouped response moments

\[
Y_A(\omega)=\frac{D_{A,0}}{D_A^{\rm(cons)}(\omega)}
=
1+u_2^{(A)}\omega^2+u_4^{(A)}\omega^4+O(\omega^6).
\]

Then
\[
u_2^{(A)}=-\frac{D_{A,2}}{D_{A,0}},
\qquad
u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2}.
\]

### 7.2 Grouped weighted trace/anomaly variables

For any grouped triple \(x_A\),
\[
x_{\rm bar}=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\]
\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}{2}.
\]

Applied to the grouped response moments,
\[
u_{\rm bar,2},\ a_2,\ b_2,\qquad
u_{\rm bar,4},\ a_4,\ b_4
\]
are the isotropic trace and the two anisotropy defects.

The grouped isotropy gate is
\[
a_2=b_2=a_4=b_4=0.
\]

### 7.3 Outgoing prefactor data

On the isotropic branch,
\[
P_0=\frac{N_0}{D_0},
\]
\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]
\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

The universal 2.5PN/4PN normalization target uses only \(P_0\) at leading order.

---

## 8. Coherent local-kernel branch variables

These are the central macroscopic variables of the later coherent-branch and invariant program.

### 8.1 Effective stiffnesses

\[
K_{U1}=K_U(1+\delta_U),
\]
\[
K_\eta^{(\mathrm{eff})}=K_\eta+6T_\Omega,
\]
\[
K_W^{(\mathrm{eff})}=K_W+\frac{\pi^2 T_W}{4L^2},
\]
\[
K_\phi^{(\mathrm{eff})}=K_\phi+\frac{\pi^2 T_\phi}{4L^2}.
\]

### 8.2 Dimensionless coherent ratios

\[
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
}
\]
\[
\boxed{
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}},
}
\]
\[
\boxed{
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
}
\]
\[
\boxed{
\delta_U=\frac{\pi^2 T_U}{L^2 K_U},
}
\]
\[
\boxed{
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
}
\]
\[
\boxed{
\zeta=\frac{\lambda_\phi^2 K_W^{(\mathrm{eff})}}{\lambda_W^2 K_\phi^{(\mathrm{eff})}},
}
\]
\[
\boxed{
\Lambda=\frac{27\pi^2 G c_s^5 K_W^{(\mathrm{eff})}}{20 a^5 c^5 \mu_W}.
}
\]

Useful coherent identities:
\[
\rho_0=\sigma_0=\chi_0,
\qquad
\epsilon_\phi=\zeta\epsilon_W,
\qquad
Z_\phi=\zeta Z_W.
\]

### 8.3 Split blocking ratio and tracking factor

The split mixed blocking ratio is
\[
\boxed{
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right].
}
\]

The exact tracking factor is
\[
\boxed{
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)}.
}
\]

On the constructive coherent branch,
\[
\frac{1}{1+\delta_U}<R_{\rm tr}<1.
\]

### 8.4 Mixed, support, and total baselines

\[
\boxed{
M_{\rm mix}
=
\frac{8 Z_W (1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)}.
}
\]

Using \(Z_\phi=\zeta Z_W\) and \(\epsilon_\phi^{\rm(split)}=\zeta\epsilon\),
\[
\boxed{
M_{\rm supp}
=
\frac{8\zeta Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\zeta\epsilon)}.
}
\]

The total baseline is
\[
\boxed{
M_{\rm tr}=M_{\rm mix}+M_{\rm supp}=M_{\rm mix}S(\zeta;\epsilon),
}
\]
with support-enhancement factor
\[
\boxed{
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon}.
}
\]

### 8.5 Transfer shape and normalization demand ratio

The coherent transfer shape is
\[
\boxed{
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}.
}
\]

The exact selected-branch demand ratio is
\[
\boxed{
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
}
\]

Equivalent identity:
\[
\boxed{
\mathcal T^2
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,
\frac{1-\epsilon_\eta}{R_{\rm target}}.
}
\]

A key coherent-branch fact is that \(R_{\rm target}\) is independent of the support ratio \(\zeta\).

---

## 9. Microscopic variables and their grouped weak-axisymmetric drifts

The key positive microscopic state vector is
\[
x=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U).
\]

Its grouped weak-axisymmetric log-drift vector is
\[
\delta\mathbf x
=
\begin{pmatrix}
\lambda_1\\
c_1\\
\gamma_1\\
\kappa_U\\
\kappa_\eta\\
\kappa_W\\
\mu_1\\
\tau_1
\end{pmatrix}
=
\begin{pmatrix}
\delta\ln\lambda_W\\
\delta\ln c_{\eta U}\\
\delta\ln\gamma\\
\delta\ln K_U\\
\delta\ln K_\eta^{(\mathrm{eff})}\\
\delta\ln K_W^{(\mathrm{eff})}\\
\delta\ln\mu_W\\
\delta\ln T_U
\end{pmatrix}_{\rm grp}.
\]

---

## 10. Microscopic slippage variables

The exact coherent-kernel slippages are

\[
\boxed{
\Sigma_Z
=
2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W
=
\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right),
}
\]
\[
\boxed{
\Sigma_\chi
=
\gamma_1+c_1-\kappa_U
=
\delta\ln\chi_0,
}
\]
\[
\boxed{
\Sigma_\eta
=
2c_1-\kappa_U-\kappa_\eta
=
\delta\ln\epsilon_\eta,
}
\]
\[
\boxed{
\Sigma_\epsilon
=
2\gamma_1+2\lambda_1-\kappa_U-\kappa_W
=
\delta\ln\epsilon_W,
}
\]
\[
\boxed{
\Sigma_\delta
=
\tau_1-\kappa_U
=
\delta\ln\delta_U.
}
\]

The tracking combination is
\[
\boxed{
\Sigma_{\rm tr}
=
(1+\chi_0)\Sigma_\delta + (1+\delta_U)\Sigma_\chi.
}
\]

The genuine nontracking transfer-shape slippage is
\[
\boxed{
\Sigma_{\rm nt}
=
\Sigma_Z
+\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
}
\]

---

## 11. Observable defect variables

The first grouped weak-axisymmetric observables are

- \(\Theta_1\): tracking-factor drift,
- \(\Xi_1\): grouped transfer-shape drift,
- \(\mathcal R_1\): selected-branch demand-ratio drift.

The exact triangular normal form is

\[
\boxed{
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr},
}
\]
\[
\boxed{
\Xi_1
=
\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)}\,\Sigma_{\rm tr}
+\Sigma_{\rm nt},
}
\]
\[
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta.
}
\]

So:

- \(\Sigma_{\rm tr}\) carries tracking failure,
- \(\Sigma_{\rm nt}\) carries nontracking transfer-shape failure,
- \(\Sigma_\eta\) carries dressing failure.

---

## 12. Intermediate exact branch composites

Stage 167 packaged the same three directions into exact branch composites.

Define
\[
B_*=\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}},
\qquad
C_*=\frac{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}{\chi_{0,*}\delta_{U,*}}.
\]

Then
\[
\mathfrak T_* = R_{\rm tr}^{-C_*},
\qquad
\mathfrak N_* = \mathcal T^2 R_{\rm tr}^{B_*},
\qquad
\mathfrak D = \epsilon_\eta,
\]
with
\[
\delta\ln\mathfrak T_*=\Sigma_{\rm tr},
\qquad
\delta\ln\mathfrak N_*=\Sigma_{\rm nt},
\qquad
\delta\ln\mathfrak D=\Sigma_\eta.
\]

These are useful, but the final closure is even sharper in the direct microscopic monomials below.

---

## 13. The three final monomial invariants

These are the final reduced invariants that matter.

### 13.1 Tracking monomial

\[
\boxed{
\mathfrak C_{{\rm tr},*}
=
\chi_0^{\,1+\delta_{U,*}}\,
\delta_U^{\,1+\chi_{0,*}}.
}
\]

It satisfies
\[
\boxed{
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr}.
}
\]

### 13.2 Nontracking monomial

Define
\[
\boxed{
E_*
=
\frac{2\epsilon_{W,*}}{1-\epsilon_*}\,
\frac{11+9\delta_{U,*}}{11(1+\delta_{U,*})},
}
\]
\[
\boxed{
F_*
=
\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+
\frac{4\epsilon_{W,*}\delta_{U,*}}{11(1-\epsilon_*)(1+\delta_{U,*})^2}.
}
\]

Then the nontracking monomial is
\[
\boxed{
\mathfrak C_{{\rm nt},*}
=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*}.
}
\]

It satisfies
\[
\boxed{
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt}.
}
\]

### 13.3 Dressing invariant

The third invariant is simply
\[
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
}
\]
with
\[
\boxed{
\delta\ln\epsilon_\eta=\Sigma_\eta.
}
\]

---

## 14. Similarity orbit and exact quotient theorem

### 14.1 Invariant map

Define the invariant map
\[
\boxed{
\mathcal I(x)=\bigl(\mathfrak C_{{\rm tr},*}(x),\ \mathfrak C_{{\rm nt},*}(x),\ \epsilon_\eta(x)\bigr).
}
\]

### 14.2 Exact monomial-drift matrix

The exact finite/infinite log-drift map is
\[
\begin{pmatrix}
\delta\ln\mathfrak C_{{\rm tr},*}\\
\delta\ln\mathfrak C_{{\rm nt},*}\\
\delta\ln\epsilon_\eta
\end{pmatrix}
=
M_*\,
\delta\mathbf x,
\]
with
\[
\boxed{
M_*=
\begin{pmatrix}
0 & 1+\delta_{U,*} & 1+\delta_{U,*} & -(2+\chi_{0,*}+\delta_{U,*}) & 0 & 0 & 0 & 1+\chi_{0,*}\\
2(1+E_*) & 0 & 2E_* & F_*-E_* & -1 & -(2+E_*) & 1 & -F_*\\
0 & 2 & 0 & -1 & -1 & 0 & 0 & 0
\end{pmatrix}.
}
\]

A useful minor gives
\[
\det M_*^{(\tau,\kappa_\eta,\mu_1)}=1+\chi_{0,*}>0,
\]
so the map has rank \(3\) and kernel dimension \(5\).

### 14.3 Exact five-parameter similarity orbit

Choose free co-scalings for
\[
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})})
\]
and determine the remaining three by monomial preservation:
\[
K_\eta^{(\mathrm{eff})}\mapsto e^{\,2C-U}K_\eta^{(\mathrm{eff})},
\]
\[
T_U\mapsto
e^{\,U-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)}T_U,
\]
\[
\mu_W\mapsto
e^{\,M(\Lambda,C,\Gamma,U,W)}\mu_W,
\]
where
\[
M(\Lambda,C,\Gamma,U,W)
=
2C-U+2W-2\Lambda
-
E_*(2\Gamma+2\Lambda-U-W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U).
\]

This defines the exact five-parameter similarity family \(\mathcal G_*\).

### 14.4 Exact finite quotient theorem

On the positive-coupling state space
\[
\mathcal M_+
=
\{(\lambda_W,c_{\eta U},\gamma,K_U,K_\eta^{(\mathrm{eff})},K_W^{(\mathrm{eff})},\mu_W,T_U)>0\},
\]
the fibres of \(\mathcal I\) are exactly the \(\mathcal G_*\)-orbits:
\[
\boxed{
\mathcal I(\widetilde x)=\mathcal I(x)
\iff
\widetilde x\in \mathcal G_*\cdot x.
}
\]

So
\[
\boxed{
\mathcal M_+/\mathcal G_*
\cong (\mathbb R_{>0})^3
\quad\text{with quotient coordinates}\quad
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
}
\]

This is the cleanest final reduced closure.

---

## 15. Final theorem ledger

### 15.1 Infinitesimal weak-axisymmetric closure

\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln\mathfrak C_{{\rm tr},*}
=
\delta\ln\mathfrak C_{{\rm nt},*}
=
\delta\ln\epsilon_\eta
=0.
}
\]

Equivalently,
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x\in T_{\rm id}\mathcal G_*.
}
\]

### 15.2 Finite closure

\[
\boxed{
\text{The coherent weak-axisymmetric defect vanishes exactly when the actual microscopic branch stays on one }
\mathcal G_*\text{-orbit,}
}
\]
that is, exactly when the three quotient coordinates
\[
(\mathfrak C_{{\rm tr},*},\ \mathfrak C_{{\rm nt},*},\ \epsilon_\eta)
\]
are preserved.

### 15.3 What is still open

The remaining theorem gap is now purely **branch-selective**:

> Does the actual completed moving-throat branch preserve the three exact quotient coordinates?

All algebraic compression is finished. The only remaining unknown is the true branch dynamics of the full PDE.

---

## 16. Minimal source anchors

This dictionary was distilled from the final moving-throat stages anchored by

- `moving_throat_pde_full.md`
- `moving_throat_pde_stage147_microscopic_log_channels.md`
- `moving_throat_pde_stage148_exact_branch_drifts.md`
- `moving_throat_pde_stage165_microscopic_coherent_slippage.md`
- `moving_throat_pde_stage166_triangular_normal_form.md`
- `moving_throat_pde_stage167_branch_invariant_coordinates.md`
- `moving_throat_pde_stage168_microscopic_monomials.md`
- `moving_throat_pde_stage169_similarity_orbit_closure.md`
- `moving_throat_pde_stage170_orbit_quotient_closure.md`

The intended rule is:

- use the **PDE engine** document for the equations and branch construction,
- use **this dictionary** for the reduced variable ledger and the final invariant theorem.
