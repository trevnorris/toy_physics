# Moving-Throat PDE Program Compact

## 0. Purpose

This document is the primary compact working master for the moving-throat PDE
program. It is meant to support new sessions and new derivation work without
requiring a full walk through the entire stage tree.

It is now a fully drafted compact program ledger: parent theory, moving-throat
geometry, reduced wall/support/outgoing stack, Family-1 mouth/core closure,
co-evolving branch, grouped real `P2` outgoing bridge, coherent local-kernel
reduction, invariant structure, theorem ledger, and the explicit realization gap.

This is **not** a paper draft.
It is a compact program ledger.

---

## 1. Reading Rules and Status Firewall

### 1.1 Claim-status tags

Every major statement in this document should be read under one of these tags.

- **Exact**: follows directly from the declared action, exact definitions, or exact algebra.
- **Exact Within Closure**: exact inside a stated reduced closure family, branch, or hierarchy.
- **Reduced / Controlled Reduction**: follows only after a stated ansatz, low-frequency reduction, projection, or branch restriction.
- **Effective Closure**: a physically motivated closure choice that is not yet a unique theorem of the completed PDE.
- **Numerically Located**: defined exactly but realized at values currently obtained by numerical solve rather than closed form.
- **Open**: still depends on the actual completed moving-throat PDE branch.

### 1.2 Non-negotiable notation firewall

The following separations are structural.

1. Electric charge is carried by
   \[
   \eta_Q,\qquad q_\star,\qquad q_{\rm eff},
   \]
   not by circulation.
2. The historical gravity-side bare `q=1` is the mass-dressing coefficient
   \[
   \kappa_\rho=1,
   \]
   not electric charge.
3. Grouped labels `20/21/22` denote grouped real `P2` lanes, not spacetime indices.
4. The mixed channels
   \[
   A_w,\qquad J^w,\qquad F_{\mu w},\qquad E_w,\qquad C_a
   \]
   are suppressed only in the strict far-field brane reduction. They remain part
   of the microscopic ontology and are required for the honest outgoing bridge.

### 1.3 Present theorem-status summary

Current best reading of the program:

- the parent `4+1` field-theory block is fixed at the exact declared-action level,
- the moving-throat geometry lift and reduced wall/support/outgoing program are explicit but not all fully realized by the completed PDE,
- the Family-1 mouth/core and coherent invariant structures are strong reduced closures,
- the main remaining theorem gap is full PDE branch realization, not algebraic compression.

---

## 2. Parent Theory and Exact Bulk Equations

### 2.1 Arena, coordinates, and indices

**Status:** `Exact`

The fundamental arena is a `(4+1)`-dimensional spacetime with coordinates
\[
x^M=(t,x,y,z,w),
\qquad
M,N\in\{0,1,2,3,4\}.
\]

Bulk spatial coordinates are
\[
\mathbf X=(x,y,z,w)\in\mathbb R^4,
\]
while brane spatial coordinates are
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

The bulk metric is
\[
\eta_{MN}=\mathrm{diag}(-1,+1,+1,+1,+1),
\]
with d'Alembertian
\[
\Box_5=-\partial_t^2+\nabla_3^2+\partial_w^2.
\]

### 2.2 Core fields

**Status:** `Exact`

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
r=\sqrt{x^2+y^2+z^2},
\qquad
\Omega=\mathbf x/r\in S^2.
\]

The bulk density is
\[
\rho=|\psi|^2.
\]

The localized gauge field is
\[
A_M=(A_0,A_i),
\qquad
F_{MN}=\partial_MA_N-\partial_NA_M.
\]

The old finite-dimensional throat variables
\[
a(t),\qquad L(t)
\]
survive only as collective moments of the moving-throat field, not as the
fundamental geometry variables.

### 2.3 Charge ontology

**Status:** `Exact`

The corrected electric-charge bookkeeping is
\[
\eta_Q=\pm 1,
\qquad
q_\star=\eta_Q e_\star,
\qquad
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

The firewall is:

- electric-charge sign is carried by \(\eta_Q\),
- circulation belongs to the magnetic/vortical sector,
- the historical gravity-side `q=1` is really \(\kappa_\rho=1\).

### 2.4 Exact parent action

**Status:** `Exact`

The exact parent theory presently fixed by the program is
\[
S=\int dtd^4X\,(\mathcal L_\psi+\mathcal L_{\rm EM}),
\]
with the geometry sector encoded through the confinement coupling
\[
V_{\rm conf}(\mathbf X;\Sigma).
\]

#### 2.4.1 Matter sector: gauged GNLS

\[
\boxed{
\mathcal L_\psi
=
\frac{i\hbar}{2}\left(\psi^*D_t\psi-\psi D_t\psi^*\right)
-\frac{\hbar^2}{2m}(D_i\psi)^*(D_i\psi)
-V_{\rm conf}(\mathbf X;\Sigma)\,\rho
-U(\rho).
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
h(\rho)=\frac{dU}{d\rho}=\frac{5K}{4}\rho^4.
\]

The corresponding bulk sound speed is
\[
c_s^2(\rho)=\frac{1}{m}\frac{dP}{d\rho}=\frac{5K}{m}\rho^4.
\]

#### 2.4.2 Localized Maxwell sector

\[
\boxed{
\mathcal L_{\rm EM}
=
-\frac{Z(w)}{4\mu_0}F_{MN}F^{MN}
-\frac{1}{2\xi\mu_0}(\partial\!\cdot\!A)^2
-A_MJ_{\rm ext}^M.
}
\]

The full source entering Maxwell's equation is
\[
J_{\rm tot}^M=J_\psi^M+J_{\rm ext}^M.
\]

Important bookkeeping rule:
varying the covariant matter action already generates the dynamical matter
current \(J_\psi^M\), so it must not be double-counted in the explicit Maxwell
source term.

### 2.5 Exact field equations already fixed by the parent theory

#### 2.5.1 Gauged GNLS equation

**Status:** `Exact`

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

#### 2.5.2 Exact current and continuity

**Status:** `Exact`

The bulk number current is
\[
\boxed{
j^i=\frac{\hbar}{m}\,\Im(\psi^\ast D_i\psi).
}
\]

Exact continuity is
\[
\boxed{
\partial_t\rho+\partial_i j^i=0.
}
\]

Where \(\rho>0\), define the bulk velocity by
\[
j^i=\rho v^i.
\]

#### 2.5.3 Exact localized Maxwell equation

**Status:** `Exact`

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
\partial_MJ_{\rm tot}^M=0.
\]

#### 2.5.4 Madelung rewrite and Euler-like form

**Status:** `Exact`

Write
\[
\psi=\sqrt{\rho}\,e^{i\theta}.
\]

The gauge-invariant bulk velocity is
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

The exact Euler-like equation is
\[
\boxed{
m(\partial_t+v_j\partial_j)v_i
=
q_\star(E_i+v_jB_{ij})
-\partial_i\!\left(V_{\rm conf}+h(\rho)+Q(\rho)\right).
}
\]

Here
\[
E_i=-\partial_tA_i-\partial_iA_0,
\qquad
B_{ij}=F_{ij}.
\]

#### 2.5.5 Exact vorticity–gauge identity

**Status:** `Exact`

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

So circulation belongs to the magnetic/vortical sector rather than the electric-charge dictionary.

#### 2.5.6 Exact mixed-sector gauge invariants

**Status:** `Exact`

The mixed fields
\[
E_w=F_{w0}=-\partial_tA_w-\partial_wA_0,
\]
\[
C_a=F_{aw}=\partial_aA_w-\partial_wA_a,
\]
are exact gauge invariants under
\[
A_0\to A_0-\partial_t\chi,
\qquad
A_a\to A_a+\partial_a\chi,
\qquad
A_w\to A_w+\partial_w\chi.
\]

These mixed channels are suppressed only in the strict far-field zero-mode brane
reduction. They remain part of the microscopic ontology and are essential to the
honest outgoing bridge.

#### 2.5.7 Cold-start projection and zero-mode reduction hooks

**Status:** `Reduced / Controlled Reduction`

For a normalized brane weight \(W(w)\),
\[
\int W(w)\,dw=1,
\]
the exact projected brane observables are
\[
\rho_{\rm brane}(\mathbf x,t)=\int W(w)\rho(\mathbf x,w,t)\,dw,
\]
\[
\mathbf j_{\rm brane}(\mathbf x,t)=\int W(w)\mathbf j_{xyz}(\mathbf x,w,t)\,dw,
\qquad
\mathbf v_{\rm brane}=\mathbf j_{\rm brane}/\rho_{\rm brane}.
\]

Projected continuity is exact:
\[
\boxed{
\partial_t\rho_{\rm brane}+\nabla_3\cdot \mathbf j_{\rm brane}=S_{\rm leak},
}
\]
with leakage term
\[
\boxed{
S_{\rm leak}
=
-\left[Wj^w\right]_{-\infty}^{+\infty}
+\int W'(w)j^w\,dw.
}
\]

Under the Helmholtz split
\[
\mathbf v_{\rm brane}=\nabla_3\varphi+\mathbf v_T,
\qquad
\nabla_3\cdot\mathbf v_T=0,
\]
the exact longitudinal identity is
\[
\boxed{
\rho_{\rm brane}\,\nabla_3^2\varphi
=
S_{\rm leak}
-\partial_t\rho_{\rm brane}
-(\nabla_3\rho_{\rm brane})\cdot(\nabla_3\varphi+\mathbf v_T).
}
\]

In the quasi-static longitudinal regime this becomes the brane Poisson hook for
\(\varphi\).

Under the controlled far-field zero-mode assumptions
\[
A_w\approx 0,
\qquad
\partial_w A_\mu\approx 0,
\qquad
J^w\approx 0,
\qquad
F_{\mu w}\approx 0,
\]
integration over \(w\) gives the effective brane Maxwell sector
\[
\boxed{
\partial_\mu F^{\mu\nu}=\mu_0^{\rm eff}J_{\rm eff}^\nu,
\qquad
\mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}}.
}
\]

This is a controlled reduction, not a denial of the mixed-core structure. It is
the short cold-start bridge between the exact parent theory and the reduced brane
language used later in the document.

### 2.6 What is fixed and what is not at this level

**Fixed at this level**

- the parent field content,
- the exact gauged GNLS plus localized Maxwell action,
- the exact bulk equations and exact mixed-sector observables,
- the corrected charge ontology.

**Not fixed at this level**

- the detailed moving-throat branch geometry,
- the reduced wall/support/outgoing hierarchy,
- the actual outgoing quadrupole normalization branch,
- the full PDE realization of the reduced coherent closures.

---

## 3. Moving-Throat Geometry and Throat-Mode Decomposition

### 3.1 Moving throat as a level-set / shape-field lift

**Status:** `Effective Closure`

The smallest moving-throat lift currently used by the program is the hybrid
level-set / shape-field representation
\[
\boxed{
\Sigma(\mathbf X,t)=r-R(\Omega,w,t),
}
\]
with
\[
r=\sqrt{x^2+y^2+z^2},
\qquad
\Omega=\mathbf x/r\in S^2.
\]

The finite throat surface is
\[
\Sigma(\mathbf X,t)=0.
\]

Sign convention:

- exterior region: \(\Sigma>0\),
- interior/support region: \(\Sigma<0\).

This lift is not yet claimed as a unique theorem of the parent PDE.
It is the smallest effective geometry choice that:

- keeps the `4+1` ontology intact,
- makes the mouth \(S^2\) multipoles explicit,
- and lets the old collective variables emerge as moments rather than as
  fundamental geometry degrees of freedom.

### 3.2 Reference stationary throat and recovery of \((a,L)\)

**Status:** `Exact Within Closure`

The reference stationary throat is
\[
\Sigma_0(\mathbf X)=r-R_0(w),
\]
with:

- mouth at \(w=0\),
- finite depth \(0\le w\le L_0\),
- mouth radius
  \[
  a_0=R_0(0),
  \]
- bottom closure through
  \[
  R_0(L_0)=0
  \]
  or an equivalent regular bottom condition.

The old collective variables are recovered as geometry moments:
\[
a(t)=\frac{1}{4\pi}\int_{S^2}R(\Omega,0,t)\,d\Omega.
\]

If \(W_b(\Omega,t)\) is defined by \(R(\Omega,W_b(\Omega,t),t)=0\), then
\[
L(t)=\frac{1}{4\pi}\int_{S^2}W_b(\Omega,t)\,d\Omega.
\]

So the old \((a,L)\) closure survives as the lowest geometry moments of the
distributed throat field.

### 3.3 Promoted confinement coupling

**Status:** `Exact Within Closure`

The old confinement potential is promoted to a moving-surface coupling by
\[
\boxed{
V_{\rm conf}(\mathbf X;\Sigma)=V_{\rm wall}\!\left(\frac{\Sigma(\mathbf X,t)}{\ell_c}\right),
}
\]
with \(\ell_c\) a wall-thickness scale.

Linearize around the reference throat by writing
\[
R(\Omega,w,t)=R_0(w)+\eta(\Omega,w,t).
\]
Then
\[
\delta\Sigma=-\eta,
\]
so the direct wall-bulk coupling is
\[
\boxed{
\delta V_{\rm conf}
=
-\frac{V_{\rm wall}'(\Sigma_0/\ell_c)}{\ell_c}\,\eta.
}
\]

This is the basic linear source through which the moving wall drives the matter
and gauge sectors in the reduced hierarchy.

### 3.4 Harmonic decomposition on the mouth sphere

**Status:** `Exact Within Closure`

Expand the wall displacement in real spherical harmonics:
\[
\eta(\Omega,w,t)
=
\eta_0(w,t)Y_{00}(\Omega)
+
\sum_{m\in P_2({\rm real})} q_{2m}(w,t)\,Y_{2m}^{\rm real}(\Omega)
+
\eta_{\ge 3}(\Omega,w,t).
}
\]

The grouped real `P2` set is
\[
\{20,\ 21c,\ 21s,\ 22c,\ 22s\}.
\]

With
\[
Y_{00}=\frac{1}{2\sqrt\pi},
\]
the physical mouth-average shift \(\delta a\) and the normalized monopole
coefficient are related by
\[
\boxed{
q_{00}(0,t)=2\sqrt\pi\,\delta a(t).
}
\]

A useful axisymmetric split is
\[
\eta_0(w,t)=\alpha_a(w)\,\delta a(t)+\alpha_L(w)\,\delta L(t)+g(w,t),
\]
where \(g(w,t)\) is the residual axisymmetric geometry lane orthogonal to the
collective \((a,L)\) directions.

So the geometry sector separates into:

- collective \(l=0\) throat motion \((\delta a,\delta L)\),
- residual axisymmetric geometry,
- grouped real `l=2` lanes,
- higher \(l\ge 3\) lanes.

### 3.5 Why the grouped real `P2` bundle matters

**Status:** `Reduced / Controlled Reduction`

The grouped real `P2` sector is the first nontrivial harmonic family beyond the
monopole. In the moving-throat program it is the literal geometry/support
realization of the conservative grouped quadrupole payload already exposed by
the earlier reduced PN hierarchy.

It is therefore not decorative harmonic bookkeeping.
It is the first throat-localized multipole bundle that the full PDE must supply
honestly.

### 3.6 Geometry-sector status at this point

**Fixed within the present geometry lift**

- level-set / shape-field throat definition,
- recovery of \((a,L)\) as collective moments,
- promoted confinement coupling,
- mouth-sphere harmonic decomposition.

**Still open**

- the fully realized moving-throat branch geometry of the completed PDE,
- whether the true branch preserves the same reduced harmonic separation cleanly
  beyond the present closure.

---

## 4. Reduced Wall and Finite-Throat Support Engine

### 4.1 Finite-throat D/N support branch

**Status:** `Exact Within Closure`

On the finite throat interval
\[
z\in[0,L],
\]
the minimal internal support equation is
\[
\psi''+k^2\psi=0.
\]

The selected finite-throat support branch imposes:

- Dirichlet at the mouth,
  \[
  \boxed{\psi(0)=0,}
  \]
- Neumann at the bottom,
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

This is the exact trapped-support ladder used downstream inside the chosen
finite-throat D/N closure.

### 4.2 Mouth DtN operator and trapped round-trip closure

**Status:** `Exact Within Closure`

For prescribed mouth datum \(\psi_m\), the finite-throat D/N branch gives the
exact mouth derivative
\[
\boxed{
Z_{00}(\omega)
=
-\frac{\omega}{c_s}\tan\!\left(\frac{\omega L}{c_s}\right).
}
\]

Its poles are exactly the D/N ladder above.

The scalar round-trip factor is
\[
R_{\rm rt}=r_0r_Le^{2ikL}.
\]

For the D/N branch,
\[
r_D=-1,
\qquad
r_N=+1,
\]
so on the exact D/N ladder
\[
R_{\rm rt}=1,
\qquad
\phi_0\equiv 0 \pmod{2\pi}.
\]

This is the trapped-support closure currently carried into the reduced program.

### 4.3 Minimal distributed wall action

**Status:** `Effective Closure`

The minimal passive distributed wall action presently used is
\[
\boxed{
S_\eta^{(2)}
=
\frac12\int dt\,dw\,d\Omega\,\sqrt{\gamma_0}
\left[
\mu_\eta(w)(\partial_t\eta)^2
-T_w(w)(\partial_w\eta)^2
-T_\Omega(w)\,\eta(-\Delta_{S^2})\eta
-K_\eta(w)\eta^2
\right].
}
\]

Here:

- \(\mu_\eta(w)\) is the effective wall inertia density,
- \(T_w(w)\) is axial wall stiffness,
- \(T_\Omega(w)\) is angular stiffness on the mouth sphere,
- \(K_\eta(w)\) is a local restoring potential.

These are fixed effective constitutive functions of the chosen reference throat.
They are not to be refit stage by stage to rescue downstream normalization.

From this point onward the wall amplitudes are written in a densitized
one-dimensional convention: after integrating over the reference sphere, the
remaining surface weight is absorbed into the effective axial coefficients and
modal amplitudes.

The resulting modal operator is
\[
\mu_\eta q_{lm,tt}
-\partial_w(T_w\partial_w q_{lm})
+\bigl[K_\eta+l(l+1)T_\Omega\bigr]q_{lm}
=
S_{lm}^{(\psi,A)}+f_{lm}^{\rm ext}.
\]

So the scalar lane \(l=0\) and grouped real `P2` lane \(l=2\) are already split
before any additional matter/gauge closure is imposed.

### 4.4 Axisymmetric reduction back to \((a,L)\)

**Status:** `Reduced / Controlled Reduction`

Using the two-mode axisymmetric truncation
\[
\eta_0(w,t)=2\sqrt\pi\,[\alpha_a(w)\delta a(t)+\alpha_L(w)\delta L(t)],
\]
the distributed wall theory reduces to
\[
L_{\rm red}^{(0)}
=
\frac12 M_{AB}\dot Q^A\dot Q^B
-\frac12 K_{AB}Q^AQ^B,
\qquad
Q^A=(\delta a,\delta L),
}
\]
with
\[
M_{AB}=4\pi\int dw\,\mu_\eta\,\alpha_A\alpha_B,
\]
\[
K_{AB}=4\pi\int dw\,[T_w\alpha_A'\alpha_B'+K_0\alpha_A\alpha_B].
\]

This is the conservative reduction back to the old \((a,L)\) matrix system.
So the distributed wall is a lift of the old closure, not a replacement that
forgets it.

### 4.5 Grouped real `P2` reduction and isotropic degeneracy

**Status:** `Reduced / Controlled Reduction`

For one grouped real quadrupole component,
\[
\eta_{2m}(\Omega,w,t)=\beta_2(w)q_{2m}(t)Y_{2m}^{\rm real}(\Omega),
\]
the reduced one-mode Lagrangian is
\[
L_{2m}=\frac12 M_2\dot q_{2m}^2-\frac12 K_2 q_{2m}^2,
\]
with
\[
M_2=\int dw\,\mu_\eta\beta_2^2,
\]
\[
K_2=\int dw\,[T_w(\beta_2')^2+(K_\eta+6T_\Omega)\beta_2^2].
\]

Before symmetry breaking or matter/gauge coupling is turned on, the grouped real
`P2` channels are degenerate on the isotropic reference throat.

That is the microscopic reason the grouped quadrupole block appears as a
degenerate bundle before additional anisotropy is introduced.

### 4.6 Wall/support-engine status at this point

**Fixed within the present reduced engine**

- finite-throat D/N support ladder,
- exact D/N mouth DtN operator,
- trapped round-trip closure on that branch,
- minimal distributed wall action,
- conservative reduction back to \((a,L)\),
- grouped real `P2` reduction and isotropic degeneracy.

**Still open**

- the full coupled matter/gauge renormalization of these reduced lanes,
- the outgoing odd response and quadrupole-normalization branch,
- the actual branch realized by the completed moving-throat PDE.

### 4.7 Linearized coupled bulk/interface skeleton

**Status:** `Reduced / Controlled Reduction`

Take a stationary reference solution
\[
\psi(\mathbf X,t)=e^{-i\mu_0 t/\hbar}\psi_0(\mathbf X),
\qquad
A_M(\mathbf X,t)=A_{M0}(\mathbf X),
\qquad
R(\Omega,w,t)=R_0(w).
\]

At the reduced linearized level:

- the matter sector is BdG-like,
- the geometry sector enters through \(\delta V_{\rm conf}\),
- and the mixed channels
  \[
  \delta A_w,\qquad \delta J^w,\qquad \delta F_{\mu w}
  \]
  remain active.

The linearized matter sector is schematically
\[
i\hbar \partial_t
\begin{bmatrix}
\delta\psi\\
\delta\psi^\ast
\end{bmatrix}
=
L_{\rm BdG}
\begin{bmatrix}
\delta\psi\\
\delta\psi^\ast
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
-\frac{V_{\rm wall}'(\Sigma_0/\ell_c)}{\ell_c}\,\eta.
\]

The linearized Maxwell sector is
\[
\boxed{
\partial_M\!\left(Z(w)\delta F^{MN}\right)
+\frac{1}{\xi}\partial^N(\partial\!\cdot\!\delta A)
=
\mu_0\,\delta J^N.
}
\]

The linearized geometry sector is
\[
\boxed{
\mu_\eta\partial_t^2\eta
-\partial_w(T_w\partial_w\eta)
-T_\Omega\Delta_{S^2}\eta
+K_\eta\eta
=
S_\eta^{(\psi)}+S_\eta^{(A)}+f_{\rm ext}.
}
\]

This is still a reduced-sector scaffold, not the fully solved coupled PDE.
But it is the point at which the wall, support, and mixed gauge sectors enter one
common linearized hierarchy.

### 4.8 Reduced conservative grouped-lane bundle

**Status:** `Exact Within Closure`

The moving-throat PDE is expected to reduce, lane by lane, to the grouped real
`P2` bundle.

For each grouped lane \(A\in\{20,21,22\}\), define:

- wall/worldtube amplitude \(q_A\),
- stable BdG support modes \(X_{A\alpha}\) with frequencies \(\varpi_{A\alpha}\),
- localized brane-like gauge coordinates \(U_{A,r}\) with frequencies \(\Omega_{U,A,r}\),
- mixed coordinates \(W_{A,r}\) with frequencies \(\Omega_{W,A,r}\),
- internal mixed-sector couplings \(R_{A,r}\).

On the stable separated-pole branch, eliminating the stable BdG modes gives the
exact conservative moments
\[
B_{A,0}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^2},
\qquad
B_{A,2}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^4},
\qquad
B_{A,4}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^6}.
\]

For each Maxwell/mixed port \(r\), define
\[
\Delta_{A,r}=\Omega_{U,A,r}^2\Omega_{W,A,r}^2-R_{A,r}^2,
\qquad
S_{A,r}=\Omega_{U,A,r}^2+\Omega_{W,A,r}^2,
\]
\[
Q_{A,r}=g_{U,A,r}^2\Omega_{W,A,r}^2+2g_{U,A,r}g_{W,A,r}R_{A,r}+g_{W,A,r}^2\Omega_{U,A,r}^2,
\qquad
G_{A,r}=g_{U,A,r}^2+g_{W,A,r}^2.
\]

Then
\[
Z_{A,0}^{(r)}=\frac{Q_{A,r}}{\Delta_{A,r}},
\qquad
Z_{A,2}^{(r)}=\frac{Q_{A,r}S_{A,r}-G_{A,r}\Delta_{A,r}}{\Delta_{A,r}^2},
\]
\[
Z_{A,4}^{(r)}=
\frac{Q_{A,r}(S_{A,r}^2-\Delta_{A,r})-S_{A,r}G_{A,r}\Delta_{A,r}}{\Delta_{A,r}^3}.
\]

Summing over ports gives \(Z_{A,n}\).

The full conservative grouped-lane operator is
\[
D_A^{\rm(cons)}(\omega)=D_{A,0}+D_{A,2}\omega^2+D_{A,4}\omega^4+O(\omega^6),
\]
with
\[
\boxed{
D_{A,0}=K_A-B_{A,0}-Z_{A,0},
}
\qquad
\boxed{
D_{A,2}=-\big(M_A+B_{A,2}+Z_{A,2}\big),
}
\]
\[
\boxed{
D_{A,4}=-\big(B_{A,4}+Z_{A,4}\big).
}
\]

So at this level the completed PDE is expected to supply a grouped conservative
bundle of static, quadratic, and quartic low-frequency coefficients lane by lane.

### 4.9 Conservative-lane status at this point

**Fixed within the present reduced hierarchy**

- the linearized coupled bulk/interface scaffold,
- the exact conservative BdG and Maxwell/mixed moment formulas on the selected
  stable reduced branch,
- the grouped conservative operator coefficients \((D_{A,0},D_{A,2},D_{A,4})\).

**Still open**

- the actual values of the branch data on the completed moving-throat solution,
- whether the physical branch lands on the isotropic/passive/outgoing route,
- the outgoing odd normalization data.

---

## 5. Family-1 Core/Mouth Closure Stack

### 5.1 Core/mouth variables and normalized Family-1 ratios

**Status:** `Exact Within Closure`

The core/mouth branch is parameterized by the microscopic throat-core quantities
\[
K_s,\qquad K_q,\qquad \lambda,\qquad g_s,\qquad g_q.
\]

The normalized Family-1 ratios are
\[
\boxed{
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g_c=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}.
}
\]

The parent compensation condition is
\[
1+\mathfrak r^2=4(\mathfrak g_c-\mathfrak r)^2,
\]
so on the compensated family
\[
\mathfrak g_c
=
\mathfrak r\pm\frac12\sqrt{1+\mathfrak r^2}.
\]

The lower compensated branch is
\[
\boxed{
\mathfrak g_-(\mathfrak r)=\mathfrak r-\frac12\sqrt{1+\mathfrak r^2}.
}
\]

### 5.2 Positive localized mouth-source theorem

**Status:** `Exact Within Closure`

Inside the first one-lane positive localized-source closure on the first D/N
interval, take a nonnegative normalized axial source profile
\[
\sigma(z)\ge 0,
\qquad
\int_0^L \sigma(z)\,dz=1.
\]

The normalized mouth-bias factor is
\[
\boxed{
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)\cos\!\left(\frac{\pi z}{2L}\right)\,dz.
}
\]

Because
\[
0\le \cos\!\left(\frac{\pi z}{2L}\right)\le 1,
\]
every positive normalized source law satisfies
\[
\boxed{
0\le \mathfrak g[\sigma]\le 1.
}
\]

So within this one-lane positive localized-source setup:

- the upper compensated Family-1 branch is impossible,
- the lower compensated branch is the only admissible compensated candidate.

This theorem does **not** cover arbitrary sign-changing, multimode, or
nonlocalized mouth data.

### 5.3 Explicit mouth boundary-layer law

**Status:** `Effective Closure`

The explicit mouth boundary layer is modeled by the positive source-density free energy
\[
\boxed{
F_{\rm mouth}[\sigma]
=
\int_0^L dz\,
\Big[
\Theta_\sigma\,\sigma\!\left(\ln\frac{\sigma}{\sigma_*}-1\right)
+
V_m(z)\sigma
\Big].
}
\]

Near the mouth, the effective potential is linearized as
\[
\boxed{
V_m(z)\approx V_1 z,
\qquad
V_1=\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-q_\star\left.\partial_zA_0\right|_{\rm m}.
}
\]

That linear potential is a reduced mouth-layer closure on the active interval
\([0,L]\), not a claim that the full throat potential is globally linear.

Using the positive Onsager current
\[
J_\sigma=-M_\sigma\,\sigma\,\partial_z\mu_\sigma^{\rm chem},
\]
the stationary zero-flux branch is exactly
\[
\boxed{
\sigma_\Pi(z)
=
\frac{\Pi e^{-\Pi z/L}}{L(1-e^{-\Pi})},
\qquad
\Pi=\frac{V_1L}{\Theta_\sigma}>0.
}
\]

So the previously ad hoc truncated exponential source family becomes the exact
zero-flux equilibrium branch of the reduced mouth-layer model.

### 5.4 Explicit core-to-mouth gain map

**Status:** `Exact Within Closure`

The explicit throat-core ansatz gives the actual coupled mouth gains
\[
\boxed{
M_s=\frac{Lg_s^2}{K_s\Theta_\sigma},
\qquad
M_q=-\frac{L(K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
}
\]

In normalized variables, define
\[
\boxed{
\Sigma_0:=\frac{Lg_s^2}{K_s\Theta_\sigma}.
}
\]
Then
\[
\boxed{
M_s=\Sigma_0,
\qquad
M_q=-\Sigma_0\,R_q,
\qquad
R_q=\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
}
\]

So the coupled Family-1 mouth law is
\[
\boxed{
\Pi=\Sigma_0\Big[1-R_q\,\mathcal S_q(\Pi)\Big].
}
\]

On the exact compensated branch,
\[
\boxed{
R_q=\frac14.
}
\]

So the outlet-consistent one-parameter mouth closure is derived rather than assumed.

### 5.5 Self-matched mouth susceptibility closure

**Status:** `Effective Closure`

The first self-matched mouth susceptibility closure identifies the source channel
with the same active shell layer as the shell/compliance mode:
\[
\boxed{
\Theta_\sigma=H_wJ_s.
}
\]

This is a same-layer closure.
It removes an otherwise free susceptibility scale rather than introducing a new
fit parameter. If the actual mouth source lives on a materially different layer,
the resulting prefactor need not stay unchanged.

Under this closure,
\[
\boxed{
\Sigma_0=\frac{20}{9}\,\widehat T_m^2,
\qquad
\widehat T_m:=\frac{\rho_w\ell\sqrt L\,\mathcal T_m}{\hbar c_{s,w}}.
}
\]

So the shell gain is explicit in terms of a normalized mouth-traction amplitude.

### 5.6 Actual Family-1 gains on the canonical point

**Status:** `Numerically Located`

For the concrete Family-1 branch,
\[
\mathfrak r_{F1}\approx 1.77799353547498.
\]

At the canonical compensation point \(\Pi_*\),
\[
\Pi_*\approx 1.50882951349316,
\qquad
\mathcal S_q(\Pi_*)\approx 0.658075937605429.
\]

The natural equal-normalized branch \(\mathfrak g_c=1\) gives
\[
M_s^{\rm nat,*}\approx 1.66854252965624,
\qquad
M_q^{\rm nat,*}\approx -0.242696939724365.
\]

The exact compensated branch gives
\[
M_s^{\rm comp,*}\approx 1.80594111095636,
\qquad
M_q^{\rm comp,*}\approx -0.451485277739090.
\]

So:

- the shell gain changes by about `8.23%`,
- the mixed-gain magnitude increases by about a factor `1.86`,
- while under the self-matched susceptibility closure the normalized mouth
  traction differs by only about `4.04%`.

This means the natural branch is on the correct sign side and not far from the
canonical compensated branch, but it does not reproduce it exactly.

### 5.7 Self-consistent explicit mouth branch

**Status:** `Exact Within Closure`

On the explicit positive mouth family, identify
\[
\mathfrak g_c=\mathfrak g_\Pi,
\]
where
\[
\boxed{
\mathfrak g_\Pi
=
\frac{2\Pi(2\Pi e^\Pi+\pi)}{(4\Pi^2+\pi^2)(e^\Pi-1)}.
}
\]

Then the mixed-to-shell ratio becomes
\[
\boxed{
R_q(\Pi)=\frac{(\mathfrak g_\Pi-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
}
\]
and the Family-1 mouth law closes to
\[
\boxed{
\Pi=\Sigma_0\Big[1-R_q(\Pi)\mathcal S_q(\Pi)\Big].
}
\]

Equivalently,
\[
\boxed{
\Sigma_0(\Pi)=\frac{\Pi}{1-R_q(\Pi)\mathcal S_q(\Pi)}.
}
\]

Under the self-matched closure,
\[
\boxed{
\widehat T_m(\Pi)=
\sqrt{\frac{9\Pi}{20\left[1-R_q(\Pi)\mathcal S_q(\Pi)\right]}}.
}
\]

So the explicit mouth branch is governed by a single scalar bias parameter
\(\Pi\), not a free gain pair.

### 5.8 Singular equal-normalized limit and unique regular canonical branch

**Status:** `Exact Within Closure`

Inside the explicit positive exponential mouth family:

1. the equal-normalized branch \(\mathfrak g_c=1\) is **not** a finite-bias branch;
2. it is reached only in the singular point-source limit \(\Pi\to\infty\);
3. and in the same limit the normalized mouth traction diverges.

So the equal-normalized branch is not a regular finite branch of the explicit
mouth-layer dynamics.

The lower compensated Family-1 value is
\[
\mathfrak g_-^{F1}\approx 0.758035078944663.
\]

Since
\[
\frac{2}{\pi}<\mathfrak g_-^{F1}<1,
\]
there exists a unique positive finite \(\Pi_*\) such that
\[
\mathfrak g_{\Pi_*}=\mathfrak g_-^{F1}.
\]

Numerically,
\[
\boxed{
\Pi_*\approx 1.50882951349316,
\qquad
\Sigma_0(\Pi_*)\approx 1.80594111095636,
\qquad
\widehat T_m(\Pi_*)\approx 0.901484054174205.
}
\]

So inside the explicit Family-1 positive exponential mouth-layer closure:

- the upper compensated branch is impossible,
- the equal-normalized branch is singular,
- the lower compensated branch is the unique regular finite-bias / finite-traction branch.

### 5.9 Family-1 mouth-stack status at this point

**Fixed within the explicit mouth closure stack**

- one-lane positive-source branch-selection theorem,
- explicit exponential mouth boundary layer,
- explicit core-to-mouth gain map,
- normalized gain ratio \(R_q\),
- self-matched susceptibility closure,
- self-consistent explicit mouth branch,
- unique regular canonical branch inside the explicit positive family.

**Still open**

- whether the actual moving-throat mouth layer is realized closely enough by this
  one-lane positive exponential closure,
- whether the same-layer susceptibility closure is the physically realized one,
- and how the full co-evolving mouth/core branch shifts the canonical point once
  backreaction is turned back on.

---

## 6. Co-Evolving Mouth/Core Branch

### 6.1 Exact canonical full-profile residual

**Status:** `Exact Within Closure`

On the explicit canonical Family-1 mouth branch, the exponential source
\[
\Sigma_*(x)=\frac{\Pi_*e^{-\Pi_*x}}{1-e^{-\Pi_*}},
\qquad x\in[0,1],
\]
generates exact shell and mixed profiles
\[
T_s(x;\Pi_*)
=
\frac{1-e^{-\Pi_* x}}{\Pi_*(1-e^{-\Pi_*})}
-
\frac{x e^{-\Pi_*}}{1-e^{-\Pi_*}},
\]
\[
T_q(x;\Pi_*)
=
A_q\sinh\!\left(\frac{\pi x}{2}\right)
-
C_q\cosh\!\left(\frac{\pi x}{2}\right)
+
C_q e^{-\Pi_* x},
\]
with
\[
C_q=
\frac{\Pi_*}{(1-e^{-\Pi_*})(\pi^2/4-\Pi_*^2)},
\]
\[
A_q=
\frac{C_q\left(\frac{\pi}{2}\sinh(\pi/2)+\Pi_*e^{-\Pi_*}\right)}
{\frac{\pi}{2}\cosh(\pi/2)}.
\]

Using the canonical compensated outlet data
\[
M_s^*=4\Sigma_m^*,
\qquad
M_q^*=-\Sigma_m^*,
\qquad
\Sigma_m^*\approx 0.451485277739090,
\]
the full mouth potential is
\[
\Phi_*(x)=4\Sigma_m^*\,T_s(x;\Pi_*)-\Sigma_m^*\,T_q(x;\Pi_*).
\]

Relative to the tangent exponential potential \(\Pi_*x\), define
\[
\boxed{
R_*(x):=\Phi_*(x)-\Pi_*x.
}
\]

Then the canonical branch satisfies
\[
\boxed{
R_*(0)=0,
\qquad
R_*'(0)=0,
\qquad
R_*''(0)
=
-\,3\Sigma_m^*\frac{\Pi_*}{1-e^{-\Pi_*}}<0.
}
\]

So the explicit exponential source is tangent-matched at the mouth, but the full
coupled potential is locally sublinear there. Within the mouth-only closure this
means the true self-consistent source broadens relative to \(\Sigma_*\).

### 6.2 First-order mouth-only correction and one-step nonlinear check

**Status:** `Reduced / Controlled Reduction`

The exact full positive source on the mouth-only branch is
\[
\Sigma_{\rm full}(x)\propto e^{-\Phi_*(x)}=e^{-\Pi_*x-R_*(x)}.
\]

Expanding to first order around \(\Sigma_*\) gives
\[
\boxed{
\delta\Sigma_{\rm act}(x)
=
-\Sigma_*(x)\,\widetilde R_*(x),
\qquad
\widetilde R_*(x)=R_*(x)-\langle R_*\rangle_*.
}
\]

Only two projected moment shifts matter:
\[
\delta\mathfrak g_{\rm act}
=
-\operatorname{Cov}_*(c,R_*),
\qquad
\delta\mathcal S_{\rm act}
=
-\operatorname{Cov}_*(K_q,R_*),
\]
where
\[
c(x)=\cos\!\left(\frac{\pi x}{2}\right),
\qquad
K_q(x)=\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}.
\]

On the canonical explicit branch,
\[
\boxed{
\delta\mathfrak g_{\rm act}\approx -0.0648069687666328,
\qquad
\delta\mathcal S_{\rm act}\approx -0.0388718368650403.
}
\]

Using the previously frozen canonical rigidity data,
\[
\mathfrak g_*' \approx 0.0714453558083195,
\qquad
A_T\approx -4.27263956256927,
\qquad
B_T\approx 0.134875005736706,
\]
the first-order retuning is
\[
\boxed{
\delta\Pi_{\rm act}\approx 0.907084414842908,
\qquad
\delta\widehat T_{m,\rm act}\approx 0.271653979462338.
}
\]

So the mouth-only corrected point is
\[
\boxed{
\Pi_{\rm corr}\approx 2.41591392833607,
\qquad
\widehat T_{m,\rm corr}\approx 1.17313803363654.
}
\]

This is a controlled first-order mouth-only correction, not an exact finite
nonlinear law. A one-step nonlinear Picard probe gives
\[
\Pi_1\approx 2.53914847609768,
\qquad
\widehat T_{m,1}\approx 1.21036942084359,
\]
which moves in the same direction and on the same scale, but does not by itself
prove full convergence of the mouth-only iteration.

### 6.3 Exact co-evolving Family-1 map

**Status:** `Exact Within Closure`

The next reduction lets the positive mouth source and the core loading ratio
co-evolve. For any normalized positive profile
\[
\Sigma(x)\ge 0,
\qquad
\int_0^1 \Sigma(x)\,dx=1,
\]
define the moments
\[
\mathfrak g[\Sigma]
=
\int_0^1 \Sigma(x)\cos\!\left(\frac{\pi x}{2}\right)\,dx,
\]
\[
\mathcal S[\Sigma]
=
\int_0^1 \Sigma(x)
\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}\,dx.
\]

On the explicit Family-1 core branch, the shell/mixed ratio is no longer free:
\[
\boxed{
\mathcal R[\Sigma]
=
\frac{(\mathfrak g[\Sigma]-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
\qquad
\mathfrak r_{F1}\approx 1.77799353547498.
}
\]

The shell and mixed kernels are
\[
\mathcal T_s[\Sigma](x)=\int_0^1 \min(x,y)\,\Sigma(y)\,dy,
\]
\[
\mathcal T_q[\Sigma](x)
=
\int_0^1
\frac{\sinh\!\left(\frac{\pi}{2}\min(x,y)\right)
\cosh\!\left(\frac{\pi}{2}(1-\max(x,y))\right)}{(\pi/2)\cosh(\pi/2)}
\Sigma(y)\,dy.
\]

So the exact co-evolving mouth potential is
\[
\boxed{
\Phi_{\Sigma_0}[\Sigma](x)
=
\Sigma_0\Big[\mathcal T_s[\Sigma](x)-\mathcal R[\Sigma]\mathcal T_q[\Sigma](x)\Big],
}
\]
and the reduced fixed-point equation is
\[
\boxed{
\Sigma(x)=
\frac{e^{-\Phi_{\Sigma_0}[\Sigma](x)}}
{\int_0^1 e^{-\Phi_{\Sigma_0}[\Sigma](y)}\,dy}.
}
\]

The lower compensated branch is still characterized by
\[
\mathfrak g[\Sigma]=\mathfrak g_*,
\qquad
\mathfrak g_*\approx 0.758035078944663,
\]
which is exactly equivalent to
\[
\boxed{
\mathcal R[\Sigma]=\frac14.
}
\]

Near the canonical branch,
\[
\boxed{
\delta\mathcal R
=
-\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+O(\delta\mathfrak g^2),
}
\]
so source broadening \((\delta\mathfrak g<0)\) drives the mixed loading ratio
above \(1/4\).

For any self-consistent profile on this reduced branch,
\[
\boxed{
\Pi[\Sigma]
=
\Phi'_{\Sigma_0}[\Sigma](0)
=
\Sigma_0\Bigl[1-\mathcal R[\Sigma]\mathcal S[\Sigma]\Bigr].
}
\]

Under the self-matched susceptibility closure,
\[
\boxed{
\Sigma_0=\frac{20}{9}\widehat T_m^2.
}
\]

### 6.4 Frozen-traction fixed point

**Status:** `Numerically Located`

Holding the old canonical traction fixed at
\[
\Sigma_0^*\approx 1.80594111095636,
\qquad
\widehat T_{m,*}\approx 0.901484054174204,
\]
the exact co-evolving map on the analyzed positive branch window converges to a
positive fixed point with
\[
\boxed{
\mathfrak g_{\rm fp}\approx 0.693352419668063,
\qquad
\mathcal S_{\rm fp}\approx 0.6216013167514007,
\qquad
\mathcal R_{\rm fp}\approx 0.2827139049082381,
\qquad
\Pi_{\rm fp}\approx 1.4885734438300713.
}
\]

So the branch survives and the mouth bias stays close to \(\Pi_*\), but the fixed
point no longer lands exactly on the compensated value \(\mathcal R=1/4\).

This is strong anti-tuning evidence inside the reduced closure: once backreaction
is allowed, the old canonical point does not protect itself.

### 6.5 Renormalized canonical branch under full co-evolution

**Status:** `Numerically Located`

The exact restoration condition is
\[
\mathfrak g_{\rm fp}(\Sigma_0)=\mathfrak g_*.
\]
On the analyzed positive branch window, this condition has a unique numerically
located root
\[
\boxed{
\Sigma_0^{\rm can}\approx 4.651033550168876,
\qquad
\widehat T_{m,\rm can}\approx 1.446708366456762.
}
\]

At that renormalized traction, the co-evolving reduced fixed point satisfies
\[
\boxed{
\mathfrak g_{\rm can}=\mathfrak g_*,
\qquad
\mathcal R_{\rm can}=\frac14,
\qquad
\mathcal S_{\rm can}\approx 0.6703621156734617,
\qquad
\Pi_{\rm can}\approx 3.8715643774790087.
}
\]

Relative to the original canonical point, exact reduced compensation costs roughly
\[
\boxed{
\frac{\Sigma_0^{\rm can}}{\Sigma_0^*}-1\approx 1.5754070949223031,
}
\]
\[
\boxed{
\frac{\Pi_{\rm can}}{\Pi_*}-1\approx 1.5659389234213572,
}
\]
\[
\boxed{
\frac{\widehat T_{m,\rm can}}{\widehat T_{m,*}}-1\approx 0.6048074946616844.
}
\]

So the compensated Family-1 branch survives the analyzed co-evolving closure, but
only as a renormalized finite-bias, finite-traction branch rather than at the old
canonical point.

### 6.6 Co-evolving mouth/core status at this point

**Fixed within the analyzed reduced closure**

- the exact co-evolving Family-1 fixed-point map,
- the first-order sign of the compensation defect transport,
- failure of exact compensation at the old canonical traction,
- the numerically located renormalized compensated branch on the analyzed
  positive window.

**Still open**

- whether the full moving-throat PDE dynamically realizes this reduced
  co-evolving branch closely enough for referee-level purposes,
- whether other microscopic branches exist outside the analyzed positive reduced
  window,
- and how the co-evolving reduced branch feeds the final outgoing grouped-`P2`
  normalization on the completed PDE solution.

---

## 7. Grouped Real `P2` Response and Outgoing-Normalization Bridge

### 7.1 Normalized grouped response and isotropy gate

**Status:** `Exact Within Closure`

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

### 7.2 Compact outgoing `l=2` fingerprint

**Status:** `Reduced / Controlled Reduction`

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

Within the reduced outgoing bridge, this is the compact passive/outgoing
quadrupole fingerprint that the moving-throat operator is trying to reproduce.

### 7.3 Outgoing prefactor and normalization target

**Status:** `Exact Within Closure`

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

On the isotropic branch, the leading odd coefficient is controlled by the static
outgoing prefactor \(P_0\).

The invariant normalization target isolated by the reduced 2.5PN bridge is
\[
\boxed{
m_{\hat 0}^{\,2}P_0
=
\frac{54Gc_s^5}{5a^5c^5}.
}
\]

On the natural point-particle source-map branch,
\[
m_{\hat 0}=1+O(a^2/r^2).
\]

### 7.4 Strongest honest reduced closure at this stage

**Status:** `Reduced / Controlled Reduction`

On the canonical outgoing DtN branch of the explicit compact passive/outgoing
grouped-`P2` model, one has
\[
\chi_Q=1.
\]

So in the strict point-particle limit, the reduced normalization stack closes on
the canonical branch:
\[
\hat m_0^{\,2}\chi_Q N_Q=1
\quad\Rightarrow\quad
N_Q=1
\]
on that branch.

This is the strongest honest reduced statement currently available:

- it is closed on the canonical outgoing DtN branch,
- it is conditional on that branch realization,
- and it is conditional on the strict point-particle limit.

It is **not** yet a theorem that the full moving-throat PDE must realize that
branch.

### 7.5 Outgoing-bridge status at this point

**Fixed at the reduced level**

- normalized grouped-response formulas,
- isotropy gate for the grouped conservative branch,
- compact outgoing `l=2` fingerprint,
- outgoing prefactor formulas,
- the sharp quadrupole-normalization target,
- the strongest reduced closure on the canonical outgoing DtN branch.

**Still open**

- whether the completed moving-throat PDE actually realizes the canonical
  passive/outgoing branch rather than a deformed one with \(\chi_Q\neq 1\),
- the actual branch selection and normalization on the true solution.

---

## 8. Coherent Local-Kernel Reduction and Microscopic Invariant Structure

### 8.1 Coherent-kernel effective ratios and tracking branch

**Status:** `Exact Within Closure`

Inside the first coherent local D/N support closure, the mixed and support lanes
couple through the same local wall/\(U\) density, so the coherent branch is
controlled by the effective stiffnesses
\[
K_{U1}=K_U(1+\delta_U),
\qquad
K_\eta^{(\mathrm{eff})}=K_\eta+6T_\Omega,
\]
\[
K_W^{(\mathrm{eff})}=K_W+\frac{\pi^2T_W}{4L^2},
\qquad
K_\phi^{(\mathrm{eff})}=K_\phi+\frac{\pi^2T_\phi}{4L^2},
\]
and the dimensionless ratios
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
\qquad
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}},
\]
\[
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U},
\]
\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\zeta=\frac{\lambda_\phi^2K_W^{(\mathrm{eff})}}{\lambda_W^2K_\phi^{(\mathrm{eff})}},
\]
\[
\Lambda=\frac{27\pi^2Gc_s^5K_W^{(\mathrm{eff})}}{20a^5c^5\mu_W}.
\]

The coherent local-kernel hypothesis forces
\[
\rho_0=\sigma_0=\chi_0,
\qquad
\epsilon_\phi=\zeta\epsilon_W,
\qquad
Z_\phi=\zeta Z_W,
\]
so the reduced Stage-27 branch lands exactly on the tracking surface
\[
\boxed{
R_{\rm tr}=R_U=R_\phi
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)}.
}
\]

On the constructive coherent branch \(\chi_0>0\), \(\delta_U>0\),
\[
\boxed{
\frac{1}{1+\delta_U}<R_{\rm tr}<1.
}
\]

### 8.2 Coherent placement map and support-compensation theorem

**Status:** `Exact Within Closure`

Define the split mixed blocking ratio
\[
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right].
\]
Then the coherent branch obeys the exact reduced placement map
\[
\boxed{
M_{\rm mix}
=
\frac{8Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
}
\]
\[
\boxed{
M_{\rm supp}
=
\frac{8\zeta Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\zeta\epsilon)},
}
\]
\[
\boxed{
M_{\rm tr}=M_{\rm mix}+M_{\rm supp}=M_{\rm mix}S(\zeta;\epsilon),
}
\]
with support-enhancement factor
\[
\boxed{
S(\zeta;\epsilon)
=
1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon}.
}
\]

The outgoing demand ratio is
\[
\boxed{
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
}
\]
so \(R_{\rm target}\) is independent of \(\zeta\). The support lane changes only
the available baseline, not the target.

The coherent tracking branch then obeys
\[
M_{\rm tr}=G_{\rm tr}(\xi,\delta;R_{\rm tr}),
\qquad
R_{\rm target}=F_{\rm tr}(\xi,\delta;R_{\rm tr}),
\]
where
\[
G_{\rm tr}(\xi,\delta;R)
=
\frac{9\xi(\xi+\delta)}{9\delta+(9+2R^2)\xi},
\]
\[
F_{\rm tr}(\xi,\delta;R)
=
\frac{[9\delta+(9+2R^2)\xi]^2[9\delta+(9+2R)\xi]^2}
{81(1-\xi)\,[9\delta^2+18\delta\xi+(9+2R^2)\xi^2]^2}.
\]

At fixed \((\xi,\delta)\), the coherent tracking branch is ordered by \(R\):
\[
\frac{dG_{\rm tr}}{dR}<0,
\qquad
\frac{dF_{\rm tr}}{dR}>0.
\]
So because \(R_{\rm tr}<1\) on the constructive split-\(U\) branch, the coherent
local kernel requires more total loading and gives less normalized response than
the old flat branch.

That deficit is still not a reduced-level no-go. On the stable-side domain
\[
0<\epsilon<1,
\qquad
0\le \zeta<1/\epsilon,
\qquad
0<\xi<1,
\]
the support-enhancement factor is strictly increasing and invertible:
\[
\boxed{
\frac{dS}{d\zeta}=\frac{1-\epsilon}{(1-\zeta\epsilon)^2}>0.
}
\]

If the mixed-only baseline lies below the tracking critical load
\[
M_{\rm crit}(\delta,R)
=
G_{\rm tr}(1,\delta;R)
=
\frac{9(1+\delta)}{9\delta+9+2R^2},
\]
and the mixed-only branch has not already reached the target, then the unique
required support ratio is
\[
\boxed{
\zeta_{\rm req}
=
\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)},
\qquad
S_{\rm req}:=\frac{M_{\rm req}}{M_{\rm mix}},
}
\]
with \(\zeta_{\rm req}<\zeta_{\rm crit}<1/\epsilon\). So there is no
reduced-level support no-go on the coherent tracking branch.

This theorem remains a reduced existence statement. The open PDE question is
whether the realized physical branch actually supplies \(\zeta\ge\zeta_{\rm req}\).

### 8.3 Microscopic slippage ledger and triangular defect normal form

**Status:** `Exact Within Closure`

The positive microscopic coherent state is
\[
x=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U),
\]
with grouped weak-axisymmetric log-drift vector
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
\end{pmatrix}_{\!\rm grp}.
\]

The five direct microscopic slippages are
\[
\Sigma_Z=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W
=
\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right),
\]
\[
\Sigma_\chi=\gamma_1+c_1-\kappa_U=\delta\ln\chi_0,
\qquad
\Sigma_\eta=2c_1-\kappa_U-\kappa_\eta=\delta\ln\epsilon_\eta,
\]
\[
\Sigma_\epsilon=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W=\delta\ln\epsilon_W,
\qquad
\Sigma_\delta=\tau_1-\kappa_U=\delta\ln\delta_U.
\]

The branch-adapted coordinates compress these to
\[
\boxed{
\Sigma_{\rm tr}=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi,
}
\]
\[
\boxed{
\Sigma_{\rm nt}
=
\Sigma_Z
+
\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+
\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
}
\]

Then the three reduced observable drifts take the exact triangular form
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

So the coherent grouped defect decomposes into:

- a tracking sector \(\Sigma_{\rm tr}\),
- a genuine nontracking transfer-shape sector \(\Sigma_{\rm nt}\),
- and a dressing sector \(\Sigma_\eta\).

### 8.4 Final direct monomial invariants

**Status:** `Reduced / Controlled Reduction`

At first grouped weak-axisymmetric/reference-branch order, the three defect
coordinates are the first logarithmic drifts of three direct microscopic
monomials.

First define the reference-branch exponents
\[
\boxed{
E_*=
\frac{2\epsilon_{W,*}}{1-\epsilon_*}\,
\frac{11+9\delta_{U,*}}{11(1+\delta_{U,*})},
}
\]
\[
\boxed{
F_*=
\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+
\frac{4\epsilon_{W,*}\delta_{U,*}}
{11(1-\epsilon_*)(1+\delta_{U,*})^2}.
}
\]

The final direct monomials are
\[
\boxed{
\mathfrak C_{{\rm tr},*}
:=
\chi_0^{\,1+\delta_{U,*}}\,
\delta_U^{\,1+\chi_{0,*}},
}
\]
\[
\boxed{
\mathfrak C_{{\rm nt},*}
:=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*},
}
\]
\[
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}}.
}
\]

Their first log drifts are exactly
\[
\boxed{
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln\epsilon_\eta=\Sigma_\eta.
}
\]

So at this reference-branch order the coherent zero-defect ledger is equivalent to
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln \mathfrak C_{{\rm tr},*}=0,
\quad
\delta\ln \mathfrak C_{{\rm nt},*}=0,
\quad
\delta\ln\epsilon_\eta=0.
}
\]

This is the sharpest honest linearized invariant statement: zero defect is
equivalent, at first grouped weak-axisymmetric/reference-branch order, to
preservation of the three direct microscopic monomials above.

### 8.5 Similarity orbit and finite quotient closure

**Status:** `Exact Within Closure`

Inside the positive coherent microscopic sector
\[
\mathcal M_+
=
\bigl\{
(\lambda_W,c_{\eta U},\gamma,K_U,K_\eta^{(\mathrm{eff})},K_W^{(\mathrm{eff})},\mu_W,T_U)>0
\bigr\}.
\]
define the invariant map
\[
\boxed{
\mathcal I(x)=\bigl(\mathfrak C_{{\rm tr},*}(x),\mathfrak C_{{\rm nt},*}(x),\epsilon_\eta(x)\bigr).
}
\]

The monomial-drift map is
\[
\begin{pmatrix}
\delta\ln\mathfrak C_{{\rm tr},*}\\[2pt]
\delta\ln\mathfrak C_{{\rm nt},*}\\[2pt]
\delta\ln\epsilon_\eta
\end{pmatrix}
=
M_*\,\delta\mathbf x,
\]
with
\[
\boxed{
M_*=
\begin{pmatrix}
0 & 1+\delta_{U,*} & 1+\delta_{U,*} & -(2+\chi_{0,*}+\delta_{U,*}) & 0 & 0 & 0 & 1+\chi_{0,*}\\[4pt]
2(1+E_*) & 0 & 2E_* & F_*-E_* & -1 & -(2+E_*) & 1 & -F_*\\[4pt]
0 & 2 & 0 & -1 & -1 & 0 & 0 & 0
\end{pmatrix}.
}
\]

On the constructive coherent branch,
\[
\det M_*^{(\tau,\kappa_\eta,\mu_1)}=1+\chi_{0,*}>0,
\]
so \(M_*\) has rank \(3\) and kernel dimension \(5\).

Choose free finite co-scalings
\[
(\Lambda,\ C,\ \Gamma,\ U,\ W)
\]
for
\[
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}),
\]
and determine the remaining three by exact monomial preservation:
\[
K_\eta^{(\mathrm{eff})}\mapsto e^{\,2C-U}K_\eta^{(\mathrm{eff})},
\]
\[
T_U\mapsto
e^{\,U-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)}T_U,
\]
\[
\mu_W\mapsto e^{\,M(\Lambda,C,\Gamma,U,W)}\mu_W,
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

This defines the exact five-parameter similarity orbit \(\mathcal G_*\). The
finite quotient theorem is:
\[
\boxed{
\mathcal I(\widetilde x)=\mathcal I(x)
\iff
\widetilde x\in \mathcal G_*\!\cdot x.
}
\]

Equivalently,
\[
\boxed{
\mathcal M_+/\mathcal G_*\cong(\mathbb R_{>0})^3
\quad\text{with quotient coordinates}\quad
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
}
\]

So within the positive coherent microscopic sector, finite weak-axisymmetric
defect motion is exactly motion in this three-dimensional quotient. This is an
invariant-structure theorem for the reduced coherent hierarchy, not yet a theorem
that the full PDE dynamically remains on one similarity orbit.

### 8.6 Coherent/invariant status at this point

**Fixed within the coherent reduced hierarchy**

- the coherent local-kernel tracking reduction,
- the exact coherent placement map and support-enhancement factor,
- the stable-side reduced support-compensation theorem,
- the microscopic slippage ledger and triangular defect normal form,
- the final direct monomial invariant coordinates,
- and the exact similarity-orbit / quotient closure in the positive coherent sector.

**Still open**

- whether the actual moving-throat PDE realizes the coherent local-kernel
  hypothesis closely enough on the physical branch,
- whether the realized branch reaches the required support ratio
  \(\zeta_{\rm req}\) before softening,
- and whether the full PDE branch preserves the three reduced coherent invariants
  closely enough that the quotient description is the right physical branch
  language rather than only the right reduced one.

---

## 9. Final Theorem Ledger

### 9.1 Exact parent-theory statements

**Status:** `Exact`

- The bulk arena, fields, charge ontology, parent action, GNLS equation,
  continuity law, Maxwell equation, Madelung rewrite, and mixed gauge-invariant
  identities in Sections `2-4` are the exact starting theory of this program.
- Nothing in the later reduced hierarchy changes those parent equations; all
  later reductions are branch-specific or closure-specific consequences.

### 9.2 Exact results within explicit reduced closures

**Status:** `Exact Within Closure`

- The finite-throat D/N support problem, the minimal distributed wall action, and
  the grouped conservative lane coefficients are exact inside the reduced
  moving-throat wall/support hierarchy.
- The one-lane positive localized-source theorem, the explicit exponential
  mouth-layer law, the exact Family-1 gain map, the self-consistent explicit
  mouth branch, and the singular exclusion of the equal-normalized branch are
  exact inside the stated positive-mouth closure.
- The co-evolving Family-1 fixed-point map, the equivalence
  \(\mathfrak g=\mathfrak g_*\iff \mathcal R=1/4\), and the mouth-slope identity
  are exact inside the analyzed co-evolving reduced closure.
- The coherent local-kernel tracking reduction, the coherent placement map, the
  support-enhancement factor, the microscopic slippage normal form, and the
  finite similarity-orbit / quotient theorem are exact inside the coherent
  reduced hierarchy and the positive coherent microscopic sector where stated.
- The grouped linear anisotropy bottleneck collapses exactly to the hidden-even
  and hidden-odd outlet combinations on the compensated isotropic branch.

### 9.3 Numerically located branch data

**Status:** `Numerically Located`

- The explicit Family-1 canonical mouth point is
  \[
  \Pi_*\approx 1.50882951349316,
  \qquad
  \widehat T_{m,*}\approx 0.901484054174205.
  \]
- The mouth-only first-order/full-profile correction gives the shifted comparison
  point
  \[
  \Pi_{\rm corr}\approx 2.41591392833607,
  \qquad
  \widehat T_{m,\rm corr}\approx 1.17313803363654,
  \]
  with a one-step nonlinear probe at
  \[
  \Pi_1\approx 2.53914847609768,
  \qquad
  \widehat T_{m,1}\approx 1.21036942084359.
  \]
- At the old canonical traction, the co-evolving reduced fixed point survives but
  shifts to
  \[
  \mathfrak g_{\rm fp}\approx 0.693352419668063,
  \qquad
  \mathcal R_{\rm fp}\approx 0.2827139049082381,
  \qquad
  \Pi_{\rm fp}\approx 1.4885734438300713.
  \]
- Exact reduced compensation is restored only on the renormalized co-evolving
  branch
  \[
  \Sigma_0^{\rm can}\approx 4.651033550168876,
  \qquad
  \widehat T_{m,\rm can}\approx 1.446708366456762,
  \qquad
  \Pi_{\rm can}\approx 3.8715643774790087.
  \]

### 9.4 Reduced or conditional endpoint statements

**Status:** `Reduced / Controlled Reduction`

- The compact outgoing quadrupole target matches the canonical GR-style
  fingerprint only on the canonical outgoing DtN branch and in the strict
  point-particle limit.
- The direct monomial compatibility law
  \[
  \Theta_1=\Xi_1=\mathcal R_1=0
  \iff
  \delta\ln\mathfrak C_{{\rm tr},*}
  =
  \delta\ln\mathfrak C_{{\rm nt},*}
  =
  \delta\ln\epsilon_\eta
  =
  0
  \]
  is the honest linearized zero-defect theorem at first grouped
  weak-axisymmetric/reference-branch order.
- The finite quotient theorem is exact inside the positive coherent microscopic
  sector, but it is still a reduced invariant-structure theorem rather than a
  completed full-PDE dynamical branch theorem.

### 9.5 Practical global verdict

**Status:** `Open`

- There is currently no known algebraic blocker, no known hidden tuning
  mechanism, and no known contradiction in the derivation chain through
  Stage `170`.
- The strongest remaining theorem gap is realization:
  whether the actual moving-throat PDE branch realizes the reduced outgoing,
  mouth/core, coherent-support, and invariant structures strongly enough for
  referee-level claims.

---

## 10. Open Realization Gap

The remaining gap is now concentrated and explicit.

### 10.1 What is no longer the main risk

- not dropped factors, sign mistakes, or stale symbolic algebra,
- not SymPy-vs-Mathematica disagreement,
- not hidden branch rescue inside the explicit Family-1 closure,
- not a physically nonsensical use of positivity, Onsager drift-diffusion, or
  grouped linear response.

### 10.2 What is still genuinely open

- Whether the completed moving-throat PDE actually realizes the canonical
  passive/outgoing DtN branch rather than a nearby deformed one.
- Whether the real mouth layer is close enough to the one-lane positive
  exponential closure, and whether the same-layer susceptibility closure is the
  physically realized one.
- Whether the full PDE branch dynamically realizes the renormalized co-evolving
  Family-1 fixed point rather than only the reduced fixed-point law.
- Whether the physical coherent-support ratio \(\zeta\) reaches the exact
  reduced requirement \(\zeta_{\rm req}\) before the branch softens out.
- Whether the true grouped weak-axisymmetric branch preserves the three coherent
  monomial invariants closely enough that the similarity-orbit quotient is the
  correct physical branch language rather than only the correct reduced one.

### 10.3 How to read the current program honestly

- The document now supports the claim that there is a mathematically coherent
  reduced derivation chain from the parent PDE setup to the present endpoint.
- It does not yet justify collapsing every reduced closure into an unconditional
  full-PDE realization theorem.
- So the right referee-facing language is:
  exact parent theory, exact statements inside explicit closures, numerically
  located reduced branches, and an explicitly isolated realization gap.

---

## 11. Quick Translation Dictionary

### 11.1 Core fields and brane reduction hooks

The exact parent fields are
\[
\psi(\mathbf X,t),
\qquad
\rho=|\psi|^2,
\qquad
A_M=(A_0,A_i),
\qquad
F_{MN}=\partial_MA_N-\partial_NA_M.
\]

The exact brane projection hooks are
\[
\rho_{\rm brane}=\int W(w)\rho(\mathbf x,w,t)\,dw,
\qquad
\mathbf j_{\rm brane}=\int W(w)\mathbf j_{xyz}(\mathbf x,w,t)\,dw,
\]
\[
\mathbf v_{\rm brane}=\mathbf j_{\rm brane}/\rho_{\rm brane},
\qquad
\partial_t\rho_{\rm brane}+\nabla_3\cdot\mathbf j_{\rm brane}=S_{\rm leak}.
\]

In the quasi-static longitudinal regime,
\[
\mathbf v_{\rm brane}=\nabla_3\varphi+\mathbf v_T,
\qquad
\nabla_3\cdot\mathbf v_T=0,
\]
and \(\varphi\) is the brane velocity-potential / Poisson-hook variable.

### 11.2 Moving-throat geometry variables

The moving throat is represented by the level-set / shape field
\[
\Sigma(\mathbf X,t)=r-R(\Omega,w,t),
\qquad
r=\sqrt{x^2+y^2+z^2},
\qquad
\Omega=\mathbf x/r.
\]

The throat surface is \(\Sigma=0\), and the wall displacement is
\[
R(\Omega,w,t)=R_0(w)+\eta(\Omega,w,t).
\]

The old collective variables
\[
a(t),\qquad L(t)
\]
are collective moments of the moving-throat field, not the fundamental geometry
variables.

### 11.3 Harmonic and grouped real `P2` variables

The mouth-sphere harmonic decomposition is
\[
\eta(\Omega,w,t)
=
\eta_0(w,t)Y_{00}(\Omega)
+\sum_{m\in P_2({\rm real})}q_{2m}(w,t)Y_{2m}^{\rm real}(\Omega)+\eta_{\ge 3}.
\]

The grouped real `P2` lanes are the weighted triplet
\[
A\in\{20,21,22\},
\]
with the full five-mode real bundle understood underneath.

The grouped conservative response variables are
\[
D_A^{\rm(cons)}(\omega)=D_{A,0}+D_{A,2}\omega^2+D_{A,4}\omega^4+O(\omega^6),
\]
\[
Y_A(\omega)=\frac{D_{A,0}}{D_A^{\rm(cons)}(\omega)}
=
1+u_2^{(A)}\omega^2+u_4^{(A)}\omega^4+O(\omega^6).
\]

The outgoing prefactor data are
\[
P_0=\frac{N_0}{D_0},
\qquad
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\qquad
P_4=\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
\]

### 11.4 Family-1 mouth/core variables

The throat-core closure variables are
\[
K_s,\qquad K_q,\qquad \lambda,\qquad g_s,\qquad g_q.
\]

Their normalized Family-1 ratios are
\[
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}.
\]

The mouth/core control variables used later are
\[
\Sigma_0=\frac{20}{9}\widehat T_m^2,
\qquad
\mathcal R=\frac{(\mathfrak g-\mathfrak r)^2}{1+\mathfrak r^2},
\]
\[
\mathcal S[\Sigma]=\int_0^1 K_q(x)\Sigma(x)\,dx,
\qquad
\Pi=\Sigma_0\,[1-\mathcal R\,\mathcal S].
\]

So:

- \(\widehat T_m\) is the normalized mouth traction,
- \(\Sigma_0\) is its exact reduced quadratic proxy,
- \(\mathcal R\) is the mixed-to-shell loading ratio,
- \(\mathcal S\) is the mouth-profile quadrupole overlap,
- \(\Pi\) is the mouth slope / bias variable.

### 11.5 Coherent local-kernel variables

The coherent reduced hierarchy is organized by
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
\qquad
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}},
\]
\[
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U},
\]
\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\zeta=\frac{\lambda_\phi^2K_W^{(\mathrm{eff})}}{\lambda_W^2K_\phi^{(\mathrm{eff})}}.
\]

The key coherent tracking factor is
\[
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\]
and the coherent support-enhancement factor is
\[
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\qquad
\epsilon=\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right].
\]

### 11.6 Microscopic grouped defect variables

The positive coherent microscopic state is
\[
x=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U).
\]

Its grouped weak-axisymmetric log-drift vector is
\[
\delta\mathbf x
=
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1)^T.
\]

The direct microscopic slippages are
\[
\Sigma_Z=\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right),
\qquad
\Sigma_\chi=\delta\ln\chi_0,
\qquad
\Sigma_\eta=\delta\ln\epsilon_\eta,
\]
\[
\Sigma_\epsilon=\delta\ln\epsilon_W,
\qquad
\Sigma_\delta=\delta\ln\delta_U.
\]

The branch-adapted coordinates are
\[
\Sigma_{\rm tr}=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi,
\]
\[
\Sigma_{\rm nt}
=
\Sigma_Z
+
\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+
\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
\]

The observable defect variables are

- \(\Theta_1\): tracking-factor drift
- \(\Xi_1\): grouped transfer-shape drift
- \(\mathcal R_1\): selected-branch demand-ratio drift

and they satisfy the exact triangular normal form
\[
\Theta_1 \leftrightarrow \Sigma_{\rm tr},
\qquad
\Xi_1 \leftrightarrow \Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1 \leftrightarrow \Sigma_\eta.
\]

### 11.7 Final invariant variables

The final reduced coherent invariant coordinates are
\[
\mathfrak C_{{\rm tr},*}
=
\chi_0^{\,1+\delta_{U,*}}\,
\delta_U^{\,1+\chi_{0,*}},
\]
\[
\mathfrak C_{{\rm nt},*}
=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*},
\]
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}}.
\]

Inside the reduced coherent hierarchy, these are the three quotient coordinates
carried by the final invariant theorem.
