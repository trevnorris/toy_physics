# Moving-Throat PDE — Phase 1 Linearized Finite-Throat Scaffold

## Purpose
This note is the first actual derivation step after the roadmap and Phase-0 theorem target.
It does **not** claim the full moving-throat PDE is solved.
Its job is narrower:

1. make the first honest geometry lift beyond the old collective variables \((a,L)\),
2. write the minimal coupled linearized bulk/interface problem that can generate the missing response data,
3. define the first unit-test branch and the first extraction formulas.

The design rule is simple: every new object introduced here must earn its keep by being needed to output one of the already-open observables:
- dynamic pole scales,
- grouped real \(P_2\) conservative data,
- geometry completion data,
- passive/outgoing quadrupole normalization.

---

## 1. Frozen input from the present stack

The current stack already fixes the following backbone.

### 1.1 Exact parent bulk sectors
- gauged 4D GNLS matter sector,
- localized 4+1 Maxwell sector,
- projection/open-system framework,
- geometry carried only through effective collective coordinates \((a(t),L(t))\).

### 1.2 What is still open
- the fully dynamical moving-throat PDE,
- the dynamic pole data
  \[
  \Omega_{1\perp},\ \Omega_{10},\ \Omega_0,\ \Omega_{20},\ \Omega_{21},\ \Omega_{22},\ \Omega_g,
  \]
- the microscopic grouped real \(P_2\) constitutive law,
- the microscopic geometry completion,
- the passive/outgoing STF quadrupole normalization.

### 1.3 Structural restriction
Do **not** impose the far-field brane-Maxwell reduction while building the microscopic PDE.
In particular, do not set to zero at the PDE level:
\[
A_w,\qquad \partial_w A_\mu,\qquad J^w,\qquad F_{\mu w}.
\]
Those are controlled brane-effective suppressions, not part of the microscopic ontology.

---

## 2. First geometry lift: from \((a,L)\) to a distributed throat field

The old parent files describe geometry only through \(a(t)\) and \(L(t)\).
That is enough for the reduced conservative hierarchy, but not enough for a genuine moving-throat PDE.

The minimal honest lift is to introduce a **level-set wall field**
\[
\Sigma(X,t)=0,
\qquad X=(x,y,z,w),
\]
with the following interpretation:
- \(\Sigma<0\): interior / throat region,
- \(\Sigma>0\): exterior medium,
- \(\Sigma=0\): moving throat interface.

A convenient linearized parameterization around a stationary finite throat is
\[
\Sigma(X,t)=\Sigma_0(X)-\eta(\Omega,s,t),
\]
where
- \(s\in[0,L_0]\) is the coordinate along the reference throat centerline,
- \(\Omega\in S^2\) labels directions on the mouth/cross-section,
- \(\eta\) is the normal wall displacement.

Here the old collective variables are recovered as the lowest moments:
\[
a(t)\sim a_0 + \frac{1}{4\pi L_0}\int_0^{L_0}\!ds\int_{S^2}\!\eta\,d\Omega,
\]
\[
L(t)\sim L_0 + \bigl[\text{cap/end displacement extracted from }\eta\bigr].
\]

So \((a,L)\) are not discarded; they reappear as coarse moments of the new field.

---

## 3. Mode content of the first geometry lift

To make contact with the already-solved PN hierarchy, the wall displacement should be decomposed into the channels the hierarchy already says are physically relevant.

Use the real spherical-harmonic expansion
\[
\eta(\Omega,s,t)
=
\eta_{00}(s,t)Y_{00}(\Omega)
+\sum_{m\in\{20,21c,21s,22c,22s\}}\eta_{2m}(s,t)Y_{2m}(\Omega)
+\eta_g(s,t)+\cdots.
\]

Interpretation:
- \(\eta_{00}\): breathing / scalar mouth response,
- \(\eta_{2m}\): grouped real \(P_2\) support channels,
- \(\eta_g\): the geometry lane not reducible to the grouped real mouth harmonics alone,
- higher multipoles: deferred unless the first extraction fails.

This is the first point where the 3PN grouped real \(P_2\) lane and the separate geometry lane become explicit wall variables rather than symbolic residuals.

---

## 4. Promote the confinement potential to depend on the distributed interface

The existing parent theory uses
\[
V_{\rm conf}(X;a,L).
\]

The first PDE lift is to replace this by
\[
V_{\rm conf}(X;\Sigma).
\]

At linear order around the stationary throat,
\[
V_{\rm conf}(X;\Sigma)
=
V_{\rm conf}(X;\Sigma_0)
+
\left(\frac{\partial V_{\rm conf}}{\partial \Sigma}\right)_0
(-\eta)
+O(\eta^2).
\]

So the geometry perturbation enters the bulk matter equation through the explicit source
\[
\delta V_{\rm conf}(X,t)= - \left(\frac{\partial V_{\rm conf}}{\partial \Sigma}\right)_0 \eta.
\]

This is the cleanest linear coupling because it is the distributed analogue of the already-frozen Hellmann–Feynman force formulas in the \((a,L)\) closure.

---

## 5. Linearized bulk/interface system

Choose a stationary finite-throat background
\[
\psi_0(X),\qquad A_{M0}(X),\qquad \Sigma_0(X),
\]
and perturb by
\[
\psi=\psi_0+\delta\psi,
\qquad
A_M=A_{M0}+\delta A_M,
\qquad
\Sigma=\Sigma_0-\eta.
\]

### 5.1 Matter sector
For actual solving it is best to use the Nambu/Bogoliubov form
\[
\delta\Psi
\equiv
\begin{pmatrix}
\delta\psi\\
\delta\psi^*
\end{pmatrix},
\]
so the linearized GNLS sector becomes schematically
\[
i\hbar\,\partial_t\delta\Psi
=
\mathcal L_{\psi,0}\,\delta\Psi
+\mathcal C_A[\delta A]
+\mathcal C_\Sigma[\eta].
\]
Here:
- \(\mathcal L_{\psi,0}\) is the stationary BdG-type operator built from \(\psi_0,A_{M0},\Sigma_0\),
- \(\mathcal C_A[\delta A]\) is the linear gauge coupling from \(\delta A_0,\delta A_i\),
- \(\mathcal C_\Sigma[\eta]\) is the wall-drive term generated by \(\delta V_{\rm conf}\).

At the level of bookkeeping,
\[
\mathcal C_\Sigma[\eta]
\sim
\begin{pmatrix}
\delta V_{\rm conf}\,\psi_0\\
-\delta V_{\rm conf}\,\psi_0^*
\end{pmatrix}.
\]

### 5.2 Gauge sector
The linearized Maxwell equation is
\[
\partial_M\!\bigl(Z(w)\,\delta F^{MN}\bigr)
+\frac1\xi\,\partial^N(\partial\!\cdot\!\delta A)
=
\mu_0\,\delta J_\psi^N,
\]
with
\[
\delta F_{MN}=\partial_M\delta A_N-\partial_N\delta A_M.
\]
The source \(\delta J_\psi^N\) is obtained from the linearized matter current.

### 5.3 Geometry sector
The distributed wall equation should be written as the field generalization of the old collective-coordinate force law. The clean linear form is
\[
\mathcal G_\Sigma\,\eta
=
f_\Sigma^{(\psi)} + f_\Sigma^{(A)} + f_{\rm ext},
\]
where
\[
f_\Sigma^{(\psi)}(\Omega,s,t)
=
-\left.\frac{\delta H_\psi}{\delta \Sigma(\Omega,s,t)}\right|_0,
\qquad
f_\Sigma^{(A)}(\Omega,s,t)
=
-\left.\frac{\delta H_A}{\delta \Sigma(\Omega,s,t)}\right|_0.
\]

This is the distributed Hellmann–Feynman generalization of
\[
F_a=-\int \rho\,\partial_a V_{\rm conf},
\qquad
F_L=-\int \rho\,\partial_L V_{\rm conf}.
\]

Important design choice:
- do **not** put phenomenological dissipation into \(\mathcal G_\Sigma\) at the first microscopic stage;
- let passive/outgoing odd parts emerge from the coupled bulk/interface problem and the outgoing branch conditions.

---

## 6. First branch conditions

The first branch to study should be the **compact, passive, outgoing** branch already isolated by the 2.5PN program.

At the linearized PDE level this means:

1. **compact finite throat:** the reference geometry has finite length \(L_0\) and finite radius \(a_0\);
2. **outgoing exterior behavior:** radiative channels satisfy outgoing-wave conditions rather than standing-wave growth or incoming-drive contamination;
3. **passivity:** the retarded response has nonnegative absorptive part in the physical channels;
4. **small-body regime:** extraction is done in the regime already used by the PN hierarchy;
5. **rotational invariance on the isotropic branch:** the grouped \(P_2\) block must collapse to a scalar operator when the background is isotropic.

---

## 7. First operational response operator

Use the same mouth/worldtube drive-measure idea already suggested by the bridge paper, but now on the linearized finite-throat PDE.

Choose a port basis
\[
\{P_A\}
=
\{P_{00},P_{20},P_{21c},P_{21s},P_{22c},P_{22s},P_g\}.
\]

Recommended effort variables:
- scalar enthalpy / breathing drive,
- grouped real \(P_2\) wall or support drive,
- geometry drive.

Recommended measured variables:
- normal matter flux,
- wall traction / generalized force,
- optional \(w\)-leakage flux,
- grouped response amplitudes.

Define projected amplitudes
\[
u_A(\omega)=\int_\Gamma \overline{P_A}\,u\,d\mu,
\qquad
j_A(\omega)=\int_\Gamma \overline{P_A}\,j\,d\mu,
\]
and the effective response matrix
\[
j_A(\omega)=\sum_B Z^{\rm eff}_{AB}(\omega)\,u_B(\omega).
\]

The low-frequency expansion is
\[
Z^{\rm eff}_{AB}(\omega)
=
Z^{(0)}_{AB} + i\omega Z^{(1)}_{AB} + \omega^2 Z^{(2)}_{AB} + i\omega^3 Z^{(3)}_{AB} + \cdots.
\]

The already-solved hierarchy tells us what to look for:
- scalar danger: \(i\omega\),
- dipole danger: \(i\omega^3\),
- universal quadrupole channel: \(i\omega^5\).

---

## 8. First extraction formulas

### 8.1 Dynamic pole data
The first dynamic observables are obtained from the pole structure of the relevant response denominators:
\[
\Omega_{1\perp},\ \Omega_{10},\ \Omega_0,\ \Omega_{20},\ \Omega_{21},\ \Omega_{22},\ \Omega_g.
\]

### 8.2 Grouped real \(P_2\) low-frequency data
From the even low-frequency expansion of the grouped channels,
\[
Y_{2m}(\omega)=1+u_2^{(2m)}\omega^2+u_4^{(2m)}\omega^4+\cdots,
\]
extract
\[
u_2^{(20)},\qquad u_2^{(21)},\qquad u_2^{(22)}.
\]

On the isotropic branch define
\[
\bar u_2=\frac{u_2^{(20)}+2u_2^{(21)}+2u_2^{(22)}}{5},
\]
\[
a_2=\frac{2u_2^{(20)}-u_2^{(21)}-u_2^{(22)}}{10},
\qquad
b_2=\frac{u_2^{(21)}-u_2^{(22)}}{2},
\]
and require
\[
a_2=b_2=0.
\]

### 8.3 Quadrupole normalization data
On the isotropic passive/outgoing quadrupole branch extract the canonical pair
\[
(\overline K_0,\overline K_2),
\]
then compute the effective odd quadrupole normalization
\[
\gamma_{\rm quad}^{\rm eff}.
\]
The target remains
\[
\gamma_{\rm quad}^{\rm eff}=\frac{2G}{5c^5}.
\]

---

## 9. First solvable benchmark branch

Before solving the full coupled wall-field problem, there is one benchmark that should be solved exactly as a unit test.

### 9.1 Geometry-frozen D/N throat benchmark
Take a straight finite throat of length \(L_0\), freeze wall motion, and solve the longitudinal support channel on the interval
\[
s\in[0,L_0]
\]
with
- Dirichlet condition at the mouth,
- Neumann condition at the cap.

Then the exact reference mouth operator is
\[
Z_{00}^{\rm DN}(\omega)
=
-\frac{\omega}{c_s}\tan\!\left(\frac{\omega L_0}{c_s}\right),
\]
with pole ladder
\[
\omega_n^{\rm pole}=
\frac{\pi c_s}{L_0}\left(n+\frac12\right).
\]

This benchmark matters because it gives an exact finite-throat reference for:
- the half-shifted pole ladder,
- the compact branch logic,
- and the later passive/outgoing low-frequency extraction.

### 9.2 Why this benchmark comes first
If the full linearized PDE cannot reduce to this exact unit-test branch in the frozen-wall limit, then the geometry lift or interface bookkeeping is wrong before any grouped \(P_2\) work begins.

---

## 10. Acceptance tests for the scaffold itself

The scaffold is acceptable only if it can be used to pose the following concrete theorem checks.

### 10.1 2PN gate
The linearized response problem must reproduce the distinction between:
- fixed zero-frequency conservative data,
- and genuinely open dynamic pole data.

### 10.2 3PN gate
The geometry-plus-grouped-\(P_2\) wall variables must make it possible, in principle, to derive
- the richer grouped real \(P_2\) constitutive law,
- the unique geometry completion.

### 10.3 2.5PN gate
The same PDE must be able to realize the compact/passive/outgoing branch and produce the orbital/worldtube STF quadrupole normalization route.

### 10.4 4PN gate
The extracted quadrupole normalization must be compatible with the hereditary bridge, which introduces no new independent normalization datum once the 2.5PN target is closed.

---

## 11. The immediate next calculations

The next useful calculations are now sharply ordered.

1. **Unit-test branch:** solve the geometry-frozen finite-throat D/N benchmark and rederive the low-frequency expansion of the mouth operator in the present notation.
2. **First geometry-on sector:** switch on only the scalar breathing field \(\eta_{00}\) and verify how the distributed Hellmann–Feynman force reproduces the old \((a,L)\) closure in the lowest-mode truncation.
3. **Grouped-\(P_2\) sector:** switch on \(\eta_{2m}\) one channel at a time, define the grouped response matrix, and test isotropy \((a_2=b_2=0)\) on the reference branch.
4. **Full coupled extraction:** compute \((\overline K_0,\overline K_2)\) on the isotropic passive/outgoing branch and test the quadrupole normalization target.

That is the first honest path from the present hierarchy toward a real moving-throat PDE theorem.
