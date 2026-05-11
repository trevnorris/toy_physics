# Moving-Throat PDE Engine
## Raw 4D GNLS + localized Maxwell + moving-throat geometry, with the grouped real `P2` bridge

## 0. Why this document exists

This document is the **engine** handoff. It is meant to be sufficient for a fresh session to reconstruct the actual moving-throat PDE program without needing the conversational history.

It does four things:

1. states the exact **4D parent action** and the exact bulk equations already fixed by the program,
2. defines the **topological defect / throat** as a moving brane–bulk surface,
3. records the **boundary data and mode structure** needed for the finite-throat internal branch,
4. and states the reduced grouped real `P2` output that the completed PDE must eventually supply.

The document is intentionally explicit about **claim status**.

- **Exact** means it follows directly from the declared action, exact definitions, or exact algebra.
- **Controlled reduction** means it follows only after a stated ansatz or low-frequency/small-body reduction.
- **Protocol closure** means it is fixed only within the declared response hierarchy, not yet from a fully solved moving-throat PDE.
- **Open** means the full PDE still has to decide it.

The strongest current reading is:

- the parent 4D field theory is exact,
- the moving-throat geometry lift and linearized coupled bulk/interface skeleton are explicit,
- the grouped real `P2` conservative/output bridge is explicit,
- and the last serious unresolved theorem gate is still the **actual branch selection and normalization** of the outgoing quadrupole channel on the true moving-throat solution.

---

## 1. Core ontology and kinematic stage

### 1.1 Bulk spacetime and indices

The fundamental arena is a `(4+1)`-dimensional spacetime with coordinates
\[
x^M=(t,x,y,z,w),
\qquad
M,N\in\{0,1,2,3,4\}.
\]

The bulk spatial coordinates are
\[
\mathbf X=(x,y,z,w)\in\mathbb R^4,
\]
and the brane spatial coordinates are
\[
\mathbf x=(x,y,z)\in\mathbb R^3.
\]

Brane indices are
\[
\mu,\nu\in\{0,1,2,3\},
\]
bulk spatial indices are
\[
i,j\in\{x,y,z,w\},
\]
and brane spatial indices are
\[
a,b\in\{x,y,z\}.
\]

The flat bulk metric is
\[
\eta_{MN}=\mathrm{diag}(-1,+1,+1,+1,+1),
\]
with d’Alembertian
\[
\Box_5=-\partial_t^2+\nabla_3^2+\partial_w^2.
\]

### 1.2 Fields

The core dynamical variables are

\[
\psi(\mathbf X,t),\qquad A_M(x),\qquad \Sigma(\mathbf X,t),
\]
or equivalently the shape field
\[
R(\Omega,w,t)\quad\text{with}\quad \Sigma=r-R(\Omega,w,t),
\]
where
\[
r=\sqrt{x^2+y^2+z^2},\qquad \Omega=\mathbf x/r\in S^2.
\]

The bulk density is
\[
\rho=|\psi|^2.
\]

The localized gauge field is
\[
A_M=(A_0,A_i),\qquad
F_{MN}=\partial_M A_N-\partial_N A_M.
\]

The old finite-dimensional throat variables survive only as collective moments:
\[
a(t),\qquad L(t).
\]

### 1.3 Corrected charge ontology

The corrected electric-charge bookkeeping is
\[
\eta_Q=\pm 1,\qquad
q_\star=\eta_Q e_\star,\qquad
e_\star>0.
\]

After canonical zero-mode brane normalization,
\[
q_{\rm eff}=\frac{q_\star}{\sqrt{Z_{\rm int}}},
\qquad
e_{\rm eff}=\frac{e_\star}{\sqrt{Z_{\rm int}}},
\qquad
Z_{\rm int}=\int_{-\infty}^{+\infty} Z(w)\,dw.
\]

Important firewall:

- electric-charge sign is carried by \(\eta_Q\),
- circulation belongs to the magnetic/vortical sector,
- the historical gravity-side `q=1` is the mass-dressing coefficient \(\kappa_\rho=1\), not electric charge.

---

## 2. Exact 4D parent action

### 2.1 Matter sector: gauged 4D GNLS

The exact matter Lagrangian density is
\[
\boxed{
\mathcal L_\psi
=
\frac{i\hbar}{2}\left(\psi^*D_t\psi-\psi D_t\psi^*\right)
-\frac{\hbar^2}{2m}(D_i\psi)^*(D_i\psi)
- V_{\rm conf}(\mathbf X;\Sigma)\,\rho
- U(\rho).
}
\]

The covariant derivatives are
\[
D_t\psi=\partial_t\psi+\frac{i q_\star}{\hbar}A_0\psi,
\qquad
D_i\psi=\partial_i\psi-\frac{i q_\star}{\hbar}A_i\psi.
\]

The frozen stiff-polytropic equation of state is
\[
P(\rho)=K\rho^5,
\qquad
U(\rho)=\frac{K}{4}\rho^5,
\qquad
h(\rho)=\frac{dU}{d\rho}=\frac{5K}{4}\rho^4,
\]
so the bulk sound speed is
\[
c_s^2(\rho)=\frac{1}{m}\frac{dP}{d\rho}=\frac{5K}{m}\rho^4.
\]

### 2.2 Localized Maxwell sector

The exact localized Maxwell Lagrangian density is
\[
\boxed{
\mathcal L_{\rm EM}
=
-\frac{Z(w)}{4\mu_0}F_{MN}F^{MN}
-\frac{1}{2\xi\mu_0}(\partial\!\cdot\!A)^2
- A_M J_{\rm ext}^M.
}
\]

The total source entering Maxwell’s equation is
\[
J_{\rm tot}^M = J_\psi^M + J_{\rm ext}^M.
\]

The key bookkeeping rule is that varying the covariant matter action already generates the dynamical matter current \(J_\psi^M\), so one must not double-count it in the explicit EM source term.

### 2.3 Legacy finite-dimensional geometry closure

The original reduced geometry sector was
\[
\mathcal L_{\rm geom}
=
\frac12 M_a \dot a^2 + \frac12 M_L \dot L^2
- E_{\rm geom}(a,L) - E_{\rm ratio}(a,L),
\]
with optional penalty
\[
E_{\rm ratio}(a,L)=\frac{\kappa}{2}(L-\alpha a)^2,
\]
and damped laws
\[
M_a\ddot a+\Gamma_a\dot a = -\frac{\partial H_{\rm tot}}{\partial a},
\qquad
M_L\ddot L+\Gamma_L\dot L = -\frac{\partial H_{\rm tot}}{\partial L}.
\]

This survives only as the **lowest collective truncation** of the moving-throat field.

---

## 3. Exact bulk equations already fixed by the parent theory

### 3.1 Gauged 4D GNLS equation

Variation with respect to \(\psi^\ast\) gives
\[
\boxed{
i\hbar D_t\psi
=
\left[
-\frac{\hbar^2}{2m}D_iD_i
+V_{\rm conf}(\mathbf X;\Sigma)
+h(\rho)
\right]\psi.
}
\]

### 3.2 Exact current and continuity

The bulk number current is
\[
\boxed{
j^i=\frac{\hbar}{m}\,\Im(\psi^\ast D_i\psi),
}
\]
and exact continuity is
\[
\boxed{
\partial_t\rho+\partial_i j^i=0.
}
\]

Where \(\rho>0\), define the bulk velocity by
\[
j^i=\rho v^i.
\]

### 3.3 Exact localized Maxwell equation

Variation with respect to \(A_N\) gives
\[
\boxed{
\partial_M\!\left(Z(w)F^{MN}\right)
+\frac{1}{\xi}\partial^N(\partial\!\cdot\!A)
=
\mu_0 J_{\rm tot}^N.
}
\]

The Bianchi identities are
\[
\partial_{[L}F_{MN]}=0.
\]

Gauge invariance requires current conservation:
\[
\partial_M J_{\rm tot}^M = 0.
\]

### 3.4 Madelung rewrite and Euler-like form

Write
\[
\psi=\sqrt{\rho}\,e^{i\theta}.
\]

The gauge-invariant velocity is
\[
\boxed{
v_i=\frac{\hbar}{m}\partial_i\theta-\frac{q_\star}{m}A_i.
}
\]

The quantum potential is
\[
\boxed{
Q(\rho)=-\frac{\hbar^2}{2m}\frac{\nabla_4^2\sqrt{\rho}}{\sqrt{\rho}}.
}
\]

The exact Euler-like bulk equation is
\[
\boxed{
m(\partial_t+v_j\partial_j)v_i
=
q_\star(E_i+v_j B_{ij})
-\partial_i\!\left(V_{\rm conf}+h(\rho)+Q(\rho)\right),
}
\]
with
\[
E_i=-\partial_tA_i-\partial_iA_0,\qquad
B_{ij}=F_{ij}.
\]

### 3.5 Exact vorticity–gauge identity

Define the bulk vorticity 2-form
\[
\Omega_{ij}=\partial_i v_j-\partial_j v_i.
\]

Away from phase singularities,
\[
\boxed{
\Omega_{ij}=-\frac{q_\star}{m}F_{ij}.
}
\]

This is why circulation belongs to the magnetic/vortical sector rather than the electric-charge dictionary.

### 3.6 Exact mixed-sector gauge invariants

The mixed fields
\[
E_w = F_{w0} = -\partial_t A_w - \partial_w A_0,
\]
\[
C_a = F_{aw} = \partial_a A_w - \partial_w A_a,
\]
are exact gauge invariants under
\[
A_0\to A_0-\partial_t\chi,\qquad
A_a\to A_a+\partial_a\chi,\qquad
A_w\to A_w+\partial_w\chi.
\]

These mixed channels are **suppressed only in the strict far-field zero-mode brane reduction**. They remain part of the microscopic ontology and are essential for the honest outgoing bridge.

---

## 4. Projection versus reduction

### 4.1 Projection: exact brane observables

For a normalized brane weight \(W(w)\),
\[
\int W(w)\,dw=1,
\]
the projected brane observables are
\[
\rho_{\rm brane}(\mathbf x,t)=\int W(w)\rho(\mathbf x,w,t)\,dw,
\]
\[
\mathbf j_{\rm brane}(\mathbf x,t)=\int W(w)\mathbf j_{xyz}(\mathbf x,w,t)\,dw,
\]
\[
\mathbf v_{\rm brane}=\mathbf j_{\rm brane}/\rho_{\rm brane}.
\]

Projected continuity is exact:
\[
\boxed{
\partial_t\rho_{\rm brane}+\nabla_3\cdot \mathbf j_{\rm brane}=S_{\rm leak},
}
\]
with
\[
\boxed{
S_{\rm leak}
=
-\left[Wj^w\right]_{-\infty}^{+\infty}
+\int W'(w)j^w\,dw.
}
\]

### 4.2 Poisson hook

With the Helmholtz decomposition
\[
\mathbf v_{\rm brane}=\nabla_3\varphi+\mathbf v_T,\qquad \nabla_3\cdot \mathbf v_T=0,
\]
one obtains the exact longitudinal identity
\[
\boxed{
\rho_{\rm brane}\,\nabla_3^2\varphi
=
S_{\rm leak}
-\partial_t\rho_{\rm brane}
-(\nabla_3\rho_{\rm brane})\cdot(\nabla_3\varphi+\mathbf v_T).
}
\]

In the quasi-static longitudinal regime this reduces to a Poisson equation on the brane.

### 4.3 Zero-mode Maxwell reduction

Under the controlled far-field zero-mode assumptions
\[
A_w\approx 0,\qquad
\partial_w A_\mu\approx 0,\qquad
J^w\approx 0,\qquad
F_{\mu w}\approx 0,
\]
integration over \(w\) gives an effective brane Maxwell sector
\[
\partial_\mu F^{\mu\nu}=\mu_0^{\rm eff}J_{\rm eff}^\nu,
\qquad
\mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}}.
\]

This is a **controlled reduction**, not a denial of the mixed-core structure.

---

## 5. Topological defect / throat definition

### 5.1 Moving throat as a level-set surface

The smallest honest moving-throat lift is
\[
\boxed{
\Sigma(\mathbf X,t)=r-R(\Omega,w,t).
}
\]

The finite throat surface is
\[
\Sigma(\mathbf X,t)=0.
\]

The sign convention is

- exterior: \(\Sigma>0\),
- interior/support region: \(\Sigma<0\).

The reference stationary throat is
\[
\Sigma_0(\mathbf X)=r-R_0(w).
\]

### 5.2 Reference stationary throat

The reference branch is a finite throat with

- mouth at \(w=0\),
- finite depth \(0\le w\le L_0\),
- mouth radius
  \[
  a_0=R_0(0),
  \]
- bottom closure through either
  \[
  R_0(L_0)=0
  \]
  or an equivalent regular bottom condition.

The old collective variables are recovered as averages:

\[
a(t)=\frac{1}{4\pi}\int_{S^2}R(\Omega,0,t)\,d\Omega,
\]
and if \(W_b(\Omega,t)\) is the bottom defined by \(R(\Omega,W_b(\Omega,t),t)=0\),
\[
L(t)=\frac{1}{4\pi}\int_{S^2}W_b(\Omega,t)\,d\Omega.
\]

### 5.3 Moving-wall confinement coupling

Promote the confinement potential to the level-set form
\[
\boxed{
V_{\rm conf}(\mathbf X;\Sigma)=V_{\rm wall}\!\left(\frac{\Sigma(\mathbf X,t)}{\ell_c}\right).
}
\]

Linearizing around the reference throat,
\[
R(\Omega,w,t)=R_0(w)+\eta(\Omega,w,t),
\]
gives
\[
\delta\Sigma=-\eta,
\]
hence
\[
\boxed{
\delta V_{\rm conf}
=
-\frac{V'_{\rm wall}(\Sigma_0/\ell_c)}{\ell_c}\,\eta.
}
\]

This is the direct wall–bulk coupling that drives the linearized matter and gauge sectors.

---

## 6. Throat mode decomposition and harmonic content

### 6.1 Real spherical-harmonic decomposition on the mouth sphere

Write the wall displacement as
\[
\eta(\Omega,w,t)
=
\eta_0(w,t)Y_{00}(\Omega)
+
\sum_{m\in P_2({\rm real})} q_{2m}(w,t)\,Y_{2m}^{\rm real}(\Omega)
+
\eta_{\ge 3}(\Omega,w,t).
\]

The grouped real `P2` set is
\[
\{20,\ 21c,\ 21s,\ 22c,\ 22s\}.
\]

### 6.2 Monopole normalization bridge

With
\[
Y_{00}=\frac{1}{2\sqrt\pi},
\]
the physical mouth-average shift \(\delta a\) is related to the normalized monopole coefficient by
\[
\boxed{
q_{00}(0,t)=2\sqrt\pi\,\delta a(t).
}
\]

A useful axisymmetric split is
\[
\eta_0(w,t)=\alpha_a(w)\,\delta a(t)+\alpha_L(w)\,\delta L(t)+g(w,t),
\]
where \(g(w,t)\) is the remaining axisymmetric geometry lane orthogonal to the collective \(a\) and \(L\) moments.

### 6.3 Why the grouped real `P2` bundle matters

The grouped real `P2` sector is the first nontrivial harmonic family beyond the monopole. It is the literal moving-throat realization of the conservative grouped real quadrupole payload already exposed by the 2PN/3PN/2.5PN hierarchy.

---

## 7. Finite-throat support branch and boundary data

### 7.1 Exact D/N finite-throat support problem

On the finite throat interval
\[
z\in[0,L],
\]
the minimal internal support equation is
\[
\psi''+k^2\psi=0.
\]

The finite-throat D/N branch imposes

- **Dirichlet** at the mouth:
  \[
  \boxed{\psi(0)=0,}
  \]
- **Neumann** at the bottom:
  \[
  \boxed{\psi'(L)=0.}
  \]

The eigenmodes are
\[
\psi_j(z)=A_j\sin(k_j z),
\qquad
k_j=\frac{\pi}{L}\left(j+\frac12\right),
\qquad
j=0,1,2,\dots
\]
with frequencies
\[
\boxed{
\omega_j=c_s k_j
=
\frac{\pi c_s}{L}\left(j+\frac12\right).
}
\]

### 7.2 Mouth DtN operator

For prescribed mouth datum \(\psi_m\), the finite-throat D/N branch gives the exact mouth derivative
\[
\boxed{
Z_{00}(\omega)
=
-\frac{\omega}{c_s}\tan\!\left(\frac{\omega L}{c_s}\right).
}
\]

Its poles are exactly the D/N ladder above.

### 7.3 Trapped-branch round-trip closure

The scalar round-trip factor is
\[
R_{\rm rt}=r_0 r_L e^{2ikL}.
\]

For the D/N branch,
\[
r_D=-1,\qquad r_N=+1,
\]
so on the exact D/N ladder
\[
R_{\rm rt}=1,
\qquad
\phi_0\equiv 0 \pmod{2\pi}.
\]

This is the exact trapped-support closure used downstream.

---

## 8. Distributed wall action and reduced wall sectors

### 8.1 Minimal distributed wall action

The smallest passive distributed geometry action used in the moving-throat program is
\[
\boxed{
S_\eta^{(2)}
=
\frac12\int dt\,dw\,d\Omega\,\sqrt{\gamma_0}
\left[
\mu_\eta(w)(\partial_t\eta)^2
-
T_w(w)(\partial_w\eta)^2
-
T_\Omega(w)\,\eta(-\Delta_{S^2})\eta
-
K_\eta(w)\eta^2
\right].
}
\]

This yields the modal operator
\[
\mu_\eta q_{lm,tt}
-
\partial_w(T_w \partial_w q_{lm})
+
\bigl[K_\eta+l(l+1)T_\Omega\bigr]q_{lm}
=
S_{lm}^{(\psi,A)}+f_{lm}^{\rm ext}.
\]

### 8.2 Axisymmetric reduction back to \((a,L)\)

Using the two-mode axisymmetric truncation
\[
\eta_0(w,t)=2\sqrt\pi\,[\alpha_a(w)\delta a(t)+\alpha_L(w)\delta L(t)],
\]
one recovers the old reduced matrix system
\[
L_{\rm red}^{(0)}
=
\frac12 M_{AB}\dot Q^A\dot Q^B
-
\frac12 K_{AB}Q^A Q^B,
\qquad Q^A=(\delta a,\delta L),
\]
with
\[
M_{AB}=4\pi\int dw\,\mu_\eta\,\alpha_A\alpha_B,
\]
\[
K_{AB}=4\pi\int dw\,[T_w \alpha_A'\alpha_B' + K_0 \alpha_A\alpha_B].
\]

### 8.3 Grouped real `P2` reduction

For one grouped real quadrupole component,
\[
\eta_{2m}(\Omega,w,t)=\beta_2(w) q_{2m}(t) Y_{2m}^{\rm real}(\Omega),
\]
one obtains
\[
L_{2m}=\frac12 M_2 \dot q_{2m}^2 - \frac12 K_2 q_{2m}^2,
\]
with
\[
M_2=\int dw\,\mu_\eta \beta_2^2,
\]
\[
K_2=\int dw\,[T_w(\beta_2')^2 + (K_\eta+6T_\Omega)\beta_2^2].
\]

On an isotropic reference throat the grouped real `P2` channels are degenerate before additional matter/gauge anisotropy is turned on.

---

## 9. Linearized coupled bulk/interface skeleton

### 9.1 Stationary background

Take a stationary reference solution
\[
\psi(\mathbf X,t)=e^{-i\mu_0 t/\hbar}\psi_0(\mathbf X),
\qquad
A_M(\mathbf X,t)=A_{M0}(\mathbf X),
\qquad
R(\Omega,w,t)=R_0(w).
\]

### 9.2 Linearized BdG-type matter sector

Write
\[
\psi=e^{-i\mu_0 t/\hbar}(\psi_0+\delta\psi),
\qquad
\delta\rho=\psi_0^\ast\delta\psi+\psi_0\delta\psi^\ast.
\]

Schematically, the linearized matter sector is
\[
i\hbar \partial_t
\begin{bmatrix}
\delta\psi\\ \delta\psi^\ast
\end{bmatrix}
=
L_{\rm BdG}
\begin{bmatrix}
\delta\psi\\ \delta\psi^\ast
\end{bmatrix}
+
C_A[\delta A_M]
+
C_\eta[\eta],
\]
with wall source entering through
\[
\delta V_{\rm conf}
=
-\frac{V'_{\rm wall}(\Sigma_0/\ell_c)}{\ell_c}\,\eta.
\]

### 9.3 Linearized Maxwell sector

Write
\[
A_M=A_{M0}+\delta A_M.
\]

Then
\[
\boxed{
\partial_M\!\left(Z(w)\delta F^{MN}\right)
+\frac{1}{\xi}\partial^N(\partial\!\cdot\!\delta A)
=
\mu_0\,\delta J^N.
}
\]

The important ontology point is that the linearized mixed channels
\[
\delta A_w,\qquad \delta J^w,\qquad \delta F_{\mu w}
\]
all remain active.

### 9.4 Linearized geometry PDE

Varying the distributed wall action gives
\[
\boxed{
\mu_\eta \partial_t^2\eta
-
\partial_w(T_w\partial_w\eta)
-
T_\Omega \Delta_{S^2}\eta
+
K_\eta \eta
=
S_\eta^{(\psi)}[\delta\psi,\delta\psi^\ast]
+
S_\eta^{(A)}[\delta A_M]
+
f_{\rm ext}.
}
\]

Projecting onto harmonics gives separate modal equations for `l=0`, grouped real `l=2`, and higher multipoles.

---

## 10. Reduced conservative kernels produced by the PDE

The moving-throat PDE is expected to reduce, lane by lane, to the grouped real `P2` reduced bundle.

For each grouped lane \(A\in\{20,21,22\}\), define:

- wall/worldtube amplitude \(q_A\),
- stable BdG support modes \(X_{A\alpha}\) with frequencies \(\varpi_{A\alpha}\),
- localized brane-like gauge coordinates \(U_{A,r}\) with frequencies \(\Omega_{U,A,r}\),
- mixed `A_w/F_{\mu w}/J^w` coordinates \(W_{A,r}\) with frequencies \(\Omega_{W,A,r}\),
- internal mixed-sector couplings \(R_{A,r}\).

### 10.1 Conservative BdG moments

\[
B_{A,0}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^2},
\qquad
B_{A,2}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^4},
\qquad
B_{A,4}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^6}.
\]

### 10.2 Conservative Maxwell/mixed moments

For each port \(r\),
\[
\Delta_{A,r}=\Omega_{U,A,r}^2\Omega_{W,A,r}^2-R_{A,r}^2,
\]
\[
S_{A,r}=\Omega_{U,A,r}^2+\Omega_{W,A,r}^2,
\]
\[
Q_{A,r}=g_{U,A,r}^2\Omega_{W,A,r}^2+2g_{U,A,r}g_{W,A,r}R_{A,r}+g_{W,A,r}^2\Omega_{U,A,r}^2,
\]
\[
G_{A,r}=g_{U,A,r}^2+g_{W,A,r}^2.
\]

Then
\[
Z_{A,0}^{(r)}=\frac{Q_{A,r}}{\Delta_{A,r}},
\]
\[
Z_{A,2}^{(r)}=\frac{Q_{A,r}S_{A,r}-G_{A,r}\Delta_{A,r}}{\Delta_{A,r}^2},
\]
\[
Z_{A,4}^{(r)}=\frac{Q_{A,r}(S_{A,r}^2-\Delta_{A,r})-S_{A,r}G_{A,r}\Delta_{A,r}}{\Delta_{A,r}^3}.
\]

Summing over ports gives \(Z_{A,n}\).

### 10.3 Outgoing-transfer moments

Define
\[
P_{A,r}=\Omega_{U,A,r}^2 g_{W,A,r}+R_{A,r}g_{U,A,r}.
\]

Then the outgoing-transfer moments begin with
\[
N_{A,0}^{(r)}=\frac{P_{A,r}^2}{\Delta_{A,r}^2},
\]
\[
N_{A,2}^{(r)}=\frac{2P_{A,r}\big(P_{A,r}S_{A,r}-\Delta_{A,r}g_{W,A,r}\big)}{\Delta_{A,r}^3},
\]
\[
N_{A,4}^{(r)}
=
\frac{
\Delta_{A,r}^2 g_{W,A,r}^2
-2\Delta_{A,r}P_{A,r}^2
-4\Delta_{A,r}P_{A,r}S_{A,r}g_{W,A,r}
+3P_{A,r}^2 S_{A,r}^2
}{\Delta_{A,r}^4}.
\]

Summing over ports gives \(N_{A,n}\).

### 10.4 Conservative grouped-lane operator

The full conservative grouped-lane operator is
\[
D_A^{\rm(cons)}(\omega)=D_{A,0}+D_{A,2}\omega^2+D_{A,4}\omega^4+O(\omega^6),
\]
with
\[
\boxed{
D_{A,0}=K_A-B_{A,0}-Z_{A,0},
}
\]
\[
\boxed{
D_{A,2}=-\big(M_A+B_{A,2}+Z_{A,2}\big),
}
\]
\[
\boxed{
D_{A,4}=-\big(B_{A,4}+Z_{A,4}\big).
}
\]

---

## 11. Grouped real `P2` normalized response and isotropy

Define the normalized grouped response
\[
Y_A(\omega)=\frac{D_{A,0}}{D_A^{\rm(cons)}(\omega)}
=
1+u_2^{(A)}\omega^2+u_4^{(A)}\omega^4+O(\omega^6).
\]

Then
\[
\boxed{
u_2^{(A)}=-\frac{D_{A,2}}{D_{A,0}},
}
\qquad
\boxed{
u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2}.
}
\]

For any grouped triple \(x_A\), define the weighted trace/anomaly variables
\[
x_{\rm bar}=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\]
\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}{2}.
\]

Applied to \(u_2^{(A)}\) and \(u_4^{(A)}\), the isotropy gate is
\[
\boxed{
a_2=b_2=a_4=b_4=0.
}
\]

On the isotropic branch the three grouped lanes collapse to common coefficients
\[
D_{20,n}=D_{21,n}=D_{22,n}=D_n,
\qquad
N_{20,n}=N_{21,n}=N_{22,n}=N_n.
\]

---

## 12. Outgoing `l=2` bridge and the remaining theorem gap

### 12.1 Compact outgoing `l=2` fingerprint

The normalized outgoing compact `l=2` branch has the low-frequency fingerprint
\[
\boxed{
\widehat Y_2^{\rm(out)}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
}
\]

### 12.2 Outgoing prefactor

The grouped-lane outgoing prefactor is
\[
P_0=\frac{N_0}{D_0},
\]
and more generally
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

On the isotropic branch, the leading odd coefficient is controlled only by the static outgoing prefactor \(P_0\).

### 12.3 Universal quadrupole-normalization condition

The invariant normalization condition isolated by the 2.5PN bridge is
\[
\boxed{
m_{\hat 0}^{\,2} P_0
=
\frac{54\,G\,c_s^5}{5\,a^5 c^5}.
}
\]

Equivalently, at leading point-particle order on the natural source-map branch \(m_{\hat 0}=1+O(a^2/r^2)\),
\[
P_0=\frac{54\,G\,c_s^5}{5\,a^5 c^5}.
\]

This is the **single sharp normalization target** still left for the completed moving-throat PDE.

---

## 13. What is already solved, and what is still open

### 13.1 Solved at the reduced-theory level

- exact 4D GNLS + localized Maxwell parent equations,
- exact projection identities,
- exact level-set throat lift,
- exact finite-throat D/N support branch,
- exact distributed wall reduction,
- exact reduced grouped-lane conservative bundle formulas,
- exact grouped-response conversion formulas,
- exact statement of the outgoing quadrupole normalization target.

### 13.2 Still open

- the fully solved nonlinear moving-throat PDE,
- the true branch data \((K_A,M_A,B_{A,n},Z_{A,n},N_{A,n})\) of the completed solution,
- the exact selection of the isotropic/passive/outgoing quadrupole branch on that solution,
- and the final hit on
  \[
  m_{\hat 0}^{\,2}P_0=\frac{54Gc_s^5}{5a^5c^5}.
  \]

---

## 14. Minimal source anchors

This engine document was distilled from the exact/reduced stack anchored by

- `4d_summary.md`
- `4d_em_fields_summary.md`
- `4d_plasma_summary.md`
- `4d_1pn_bridge_summary.md`
- `moving_throat_pde_full.md`
- `moving_throat_pde_stage164_microscopic_log_channels.md`
- `moving_throat_pde_stage165_exact_branch_drifts.md`
- `moving_throat_pde_stage185_microscopic_monomials.md`
- `moving_throat_pde_stage186_similarity_orbit_closure.md`
- `moving_throat_pde_stage187_orbit_quotient_closure.md`

The intended use is: a fresh session should be able to start from **this document alone**, then only drill into the stage files if it wants the longer derivation path.
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

Stage 184 packaged the same three directions into exact branch composites.

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
- `moving_throat_pde_stage164_microscopic_log_channels.md`
- `moving_throat_pde_stage165_exact_branch_drifts.md`
- `moving_throat_pde_stage182_microscopic_coherent_slippage.md`
- `moving_throat_pde_stage183_triangular_normal_form.md`
- `moving_throat_pde_stage184_branch_invariant_coordinates.md`
- `moving_throat_pde_stage185_microscopic_monomials.md`
- `moving_throat_pde_stage186_similarity_orbit_closure.md`
- `moving_throat_pde_stage187_orbit_quotient_closure.md`

The intended rule is:

- use the **PDE engine** document for the equations and branch construction,
- use **this dictionary** for the reduced variable ledger and the final invariant theorem.
