# Paper 7 (Hard-Mode 4D) — Definition Freeze Sheet v0.1

Use this as the **single source of truth** for *what the model is*. Once you freeze these items, anything downstream (derivations, numerics, response extraction, scans) is well-defined and doesn’t drift.

---

## A) Coordinates, constants, and EOS

**Coordinates (4D space + time)**
[
\mathbf X=(x,y,z,w),\quad t,\quad \nabla_4=(\partial_x,\partial_y,\partial_z,\partial_w)
]

**Constants (simulation units may nondimensionalize later)**
[
\hbar,; m,; K,; c
]

**EOS (non-negotiable)**
[
P(\rho)=K\rho^5,\quad
c_s^2(\rho)=\frac{dP}{d\rho}=5K\rho^4,\quad
h(\rho)=\frac{5K}{4}\rho^4,\quad
U(\rho)=\frac{K}{4}\rho^5
]

---

## B) Fields and derived quantities

### Fluid field (always present)

* Complex order parameter: (\psi(\mathbf X,t))
* Density: (\rho = |\psi|^2)
* Current:
  [
  \mathbf J=\frac{\hbar}{2im}\left(\psi^*\nabla_4\psi-\psi\nabla_4\psi^*\right)
  ]
* Continuity (must hold): (\partial_t\rho+\nabla_4\cdot \mathbf J=0)

### Wave-support field (optional, separate experiment)

* Real scalar surrogate: (A(\mathbf X,t))

---

## C) Geometry DOFs and smooth-step conventions

**Geometry DOFs (baseline):**

* ([x]) (a) (radius)  ([x]) (L) (length)
* (optional later) ([,]) (a(\cdot)) shape field

**Smooth step** (frozen):
[
S(u)=\tfrac12\left(1+\tanh u\right)
]
**SmoothAbs for axial/endcap gate (frozen):**
[
|w|\ \to\ {\rm SmoothAbs}(w)=\sqrt{w^2+\varepsilon_w^2}.
]
**Localization regime assumption (frozen):**
[
\varepsilon_w \ll \delta_\parallel \ll L.
]
**Wall thickness / smoothing parameters (frozen symbols):** (\delta_r,\delta_\parallel,\varepsilon_w) and steepness exponent (p).

---

## D) Confinement potential (V_{\rm conf}) (THIS MUST BE EXPLICIT)

**One 4D model everywhere:**
[
V_{\rm conf}(\mathbf X;a,L)=V_{\rm brane}(w;\mathbf X)+V_{\rm wall}(\mathbf X;a,L)+V_{\rm cap}(\mathbf X;a,L)
]

### D1) Choose baseline family (freeze one as “primary”; keep the other as comparative)

* ([x]) **Family 1**: modulated brane trap (local weakening in throat) (primary)
* ([,]) **Family 2**: true 4D transverse tube (via (R_\gamma)) (comparative)

### D2) Gate definition (freeze)

Frozen baseline choice:
[
s(\mathbf X)=w,\quad r(\mathbf X)=R_3=\sqrt{x^2+y^2+z^2}.
]

**Gate (tube + endcaps):**
[
G(\mathbf X;a,L)=G_r(R_3;a)\,G_w(w;L)
]
[
G_r=1-S!\left(\frac{R_3-a}{\delta_r}\right),\qquad
G_w=1-S!\left(\frac{{\rm SmoothAbs}(w)-L/2}{\delta_\parallel}\right).
]
So (G\approx 1) deep inside, (G\approx 0) outside.

### D3) Family 1 explicit form (modulated (w)-trap)

**Brane trap (harmonic in (w), stiffness varies in (\mathbf X)):**
[
V_{\rm brane}(w;\mathbf X)=\tfrac12 m,\Omega_w^2(\mathbf X),w^2
]
[
\Omega_w^2(\mathbf X)=\Omega_{\rm out}^2-\big(\Omega_{\rm out}^2-\Omega_{\rm in}^2\big),G(\mathbf X;a,L)
]
**Throat wall + endcaps (soft barriers):**
[
V_{\rm conf}=V_{\rm brane}+V_{\rm wall}+V_{\rm cap},
]
[
V_{\rm wall}=V_0,S!\left(\frac{R_3-a}{\delta_r}\right)^p,\qquad
V_{\rm cap}=V_0,S!\left(\frac{{\rm SmoothAbs}(w)-L/2}{\delta_\parallel}\right)^p
]
The explicit endcap barrier (V_{\rm cap}) is included in the baseline; this mitigates drain/leak.

### D4) Family 2 explicit form (true 4D transverse radius)

Freeze anisotropic transverse radius:
[
R_\gamma=\sqrt{r_\perp^2+\left(\frac{w}{\gamma}\right)^2}
]
Then set (r \equiv R_\gamma) in the wall term and (optionally) keep a separate brane trap outside the throat.

---

## E) Wave confinement (\mu_A^2(\mathbf X;a,L)) (for wave-supported track)

Wave PDE:
[
\partial_t^2 A - c^2\nabla_4^2 A + \mu_A^2(\mathbf X;a,L)\,A = S_{\rm drive}
]

**Freeze a confinement form** (simple and consistent with the same gate):
[
\mu_A^2(\mathbf X;a,L) = \mu_{\rm in}^2 + (\mu_{\rm out}^2-\mu_{\rm in}^2),[1-G(\mathbf X;a,L)]
]
(So waves propagate inside and are suppressed outside, creating standing-mode capability.)

---

## F) Field equations (frozen)

### Fluid (4D GNLS)

[
i\hbar,\partial_t\psi
=====================

\left[
-\frac{\hbar^2}{2m}\nabla_4^2
+V_{\rm conf}(\mathbf X;a,L)
+\frac{5K}{4}|\psi|^8
\right]\psi
]

**Constrained stationary equation (fixed (N), baseline scans):**
[
-\frac{\hbar^2}{2m}\nabla_4^2\psi
+V_{\rm conf}\psi
+\frac{5K}{4}|\psi|^8\psi
=\mu\psi,
\qquad
N=\int |\psi|^2\,d^4X\ \text{fixed.}
]
**Notation note:** use a single symbol consistently for (\psi) (avoid script variable names like `psiStat` vs `psi0`).

### Wave (optional)

[
\partial_t^2 A - c^2\nabla_4^2 A + \mu_A^2(\mathbf X;a,L)\,A = S_{\rm drive}
]

---

## G) Energies and wall law (closure)

### G1) Energies (frozen)

Fluid Hamiltonian:
[
H_{\rm fluid}=\int\left[\frac{\hbar^2}{2m}|\nabla_4\psi|^2+V_{\rm conf}|\psi|^2+\frac{K}{4}|\psi|^{10}\right]d^4X
]

Wave Hamiltonian:
[
H_{\rm wave}=\int \tfrac12\left[(\partial_t A)^2+c^2|\nabla_4 A|^2+\mu_A^2A^2\right]d^4X
]

Geometry energy (must be explicit):
[
E_{\rm geom}(a,L)=P_{\rm vac},V(a,L)+\sigma,A(a,L)+\kappa_b,B(a,L)+\cdots
]

### G2) Freeze geometry measures (V(a,L),A(a,L))

Pick a baseline “4D tube” convention (recommended starting point):

**4D cylinder = (B^3(a)\times[0,L])**

* 4D volume:
  [
  V(a,L)=\mathrm{Vol}(B^3(a)),L=\left(\frac{4\pi}{3}a^3\right)L
  ]
* Boundary 3D “area” (hypersurface measure):
  [
  A(a,L)=\underbrace{\mathrm{Area}(S^2(a)),L}*{4\pi a^2 L}
  +\underbrace{2,\mathrm{Vol}(B^3(a))}*{2\cdot\frac{4\pi}{3}a^3}
  ]
  (If you later change the shape family, you must update these consistently.)
  [
  A(a,L)=(4\pi a^2)L+2\left(\frac{4\pi}{3}a^3\right)
  ]
  [
  E_{\rm geom}(a,L)=P_{\rm vac}V(a,L)+\sigma A(a,L)
  ]
  Optional derivatives:
  [
  \partial_a E_{\rm geom}
  =P_{\rm vac}(4\pi a^2 L)+\sigma(8\pi a L+8\pi a^2),
  ]
  [
  \partial_L E_{\rm geom}
  =P_{\rm vac}\left(\frac{4\pi}{3}a^3\right)+\sigma(4\pi a^2).
  ]

### G3) Two experiments (keep separate; both frozen)

* **Experiment F (fluid-only):**
  [
  H_{\rm tot}^{(F)}=H_{\rm fluid}+E_{\rm geom}
  ]
* **Experiment W (wave-supported):**
  [
  H_{\rm tot}^{(W)}=H_{\rm fluid}+H_{\rm wave}+E_{\rm geom}
  ]

### G4) Closure choice (freeze)

* ([,]) **Equilibrium**: (\partial_a H_{\rm tot}=0,;\partial_L H_{\rm tot}=0)
* ([,]) **Dynamics**:
  [
  M_a\ddot a+C_a\dot a=-\partial_a H_{\rm tot},\quad
  M_L\ddot L+C_L\dot L=-\partial_L H_{\rm tot}
  ]

### G5) Constraint choice for equilibrium (freeze)

* ([x]) Baseline: fixed norm (N=\int |\psi|^2\,d^4X)
* ([,]) Alternate/future: fixed effective mass target (m_{\rm eff}=H_{\rm tot}/c^2)
* ([,]) Other: ____________________

---

## H) Brane observable map (how 3D “appears”)

**Projection weight (freeze):**

* ([x]) (W(w)=|\chi_0(w)|^2) (preferred)
* ([,]) Other: ____________________

If (V_{\rm brane}) is harmonic far away with (\Omega_{\rm out}), ground-state width:
[
\ell_w=\sqrt{\frac{\hbar}{m\Omega_{\rm out}}},\quad
\chi_0(w)=\left(\frac{1}{\pi \ell_w^2}\right)^{1/4}\exp!\left(-\frac{w^2}{2\ell_w^2}\right)
]

**Brane-projected density:**
[
\rho_{\rm brane}(x,y,z,t)=\int W(w),|\psi(x,y,z,w,t)|^2,dw
]

**Numerical projection window (frozen):**
[
|w|\le W_{\rm proj},\quad {\cal N}_W=\int_{-W_{\rm proj}}^{W_{\rm proj}} W(w),dw,\quad
\widetilde W(w)=W(w)/{\cal N}_W\quad (|w|\le W_{\rm proj}),
]
Tail mass (projection error): (1-{\cal N}_W).

**Effort variable (frozen):**
[
u=\delta h(\rho_{\rm brane})\approx 5K\rho_{\rm brane,0}^3,\delta\rho_{\rm brane}
]

---

## I) Measurement region (\Gamma), ports (P_i), flux definition

### I1) Measurement region (\Gamma) (freeze one)

* ([x]) 2-sphere on brane: (R_3=r_{\rm port},; w=0)  (**recommended for clean multipoles**)
* ([,]) 2D disk on brane plane: (z=0,; r_\perp\le r_{\rm port},; w=0)
* ([,]) 3D slab near brane: (|w|\le w_{\rm port}) + region in ((x,y,z))

**Surface measure + normal (frozen):**
[
d\mu=r_{\rm port}^2\sin\theta\,d\theta\,d\phi,\quad \hat n=\hat R_3\ \text{(outward radial at } w=0).
]

**Quadrature (frozen):**
[
\theta: \text{Gauss--Legendre (}n_\theta=32\text{)},\quad
\phi: \text{uniform (}n_\phi=64\text{)}.
]

### I2) Port basis (freeze)

If using brane 2-sphere:
[
P_{lm}(\theta,\phi)=Y_{lm}(\theta,\phi),\quad l=0..l_{\max}
]
Minimal: (l_{\max}=2) (monopole + dipole + quadrupole).
Inner product uses conjugation for complex (Y_{lm}): (\int \overline{P_i} P_j d\mu).

### I3) Flux definition (j) (freeze)

Define a brane-projected current:
[
\mathbf J_{\rm brane}(x,y,z,t)=\int W(w),\mathbf J(x,y,z,w,t),dw
]
Then choose:

* ([,]) **Mouth/through-flux (Gauss form)**:
  [
  j(\mathbf s,t)=\mathbf J_{\rm brane}\cdot \hat n\quad \text{on } \Gamma=S^2
  ]
* ([x]) **Leakage flux into bulk** (monitor (J_w)):
  [
  J_w=\frac{\hbar}{2im}\left(\psi^*\partial_w\psi-\psi\,\partial_w\psi^*\right)
  ]
  [
  j_w(t)=\int_{\Omega_{xyz}} J_w\,d^3x
  \quad \text{evaluated at } w=\pm W_{\rm cutoff}\ (\text{baseline: } R_3<R_{\rm measure})
  ]
  Use outward-oriented cap fluxes: (j^{\rm out}_+(t)=\int J_w|_{w=+W_{\rm cut}} d^3x), (j^{\rm out}_-(t)=\int[-J_w|_{w=-W_{\rm cut}}] d^3x), (j_{\rm leak}=j^{\rm out}_+ + j^{\rm out}_-).
* ([,]) Multi-output: ((j_{\rm mouth},j_w))

### I4) Port amplitudes (frozen)

[
u_i(t)=\int_\Gamma \overline{P_i(\mathbf s)}\,u(\mathbf s,t)\,d\mu,\quad
j_i(t)=\int_\Gamma \overline{P_i(\mathbf s)}\,j(\mathbf s,t)\,d\mu
]

---

## J) Drive protocol (freeze the definition of “input”)

Primary definition (recommended):

* **Potential modulation drive near (\Gamma)**:
  [
  V_{\rm conf}(\mathbf X;a,L)\to V_{\rm conf}(\mathbf X;a,L)+\epsilon\cos(\omega t)\,f(\mathbf X)\,Y_{\ell m}(\theta(\mathbf X),\phi(\mathbf X))
  ]
  with
  [
  f(\mathbf X)=\exp\!\left(-\frac{w^2}{w_{\rm drive}^2}\right)\exp\!\left(-\frac{(R_3-r_{\rm port})^2}{r_{\rm drive}^2}\right),
  \quad R_3=\sqrt{x^2+y^2+z^2}.
  ]
  Angle map used for ports (matches `definitions.wl`):
  [
  \theta=\arccos\!\left(\frac{z}{\sqrt{x^2+y^2+z^2+\varepsilon_{\rm ang}^2}}\right),\quad
  \phi=\mathrm{atan2}(y,x).
  ]

Wave-supported track:

* [
  S_{\rm drive}=\epsilon\cos(\omega t),s(\mathbf X),P_k(\mathbf s(\mathbf X))
  ]

---

## K) Effective operator definition (frozen)

Frequency-domain response (port-to-port):
[
\mathbf j(\omega)=Z^{\rm eff}(\omega),\mathbf u(\omega)
]
Here (j_i) uses the mouth flux on (\Gamma): (j(\mathbf s,t)=\mathbf J_{\rm brane}\cdot\hat n), and leakage is a separate diagnostic.
Frozen extraction protocol:
* discard (N_d) drive periods: (N_d = nDiscardPeriodsDefault = 5)
* measure (N_m) drive periods: (N_m = nMeasurePeriodsDefault = 10)
* lock-in estimator: (\widehat s(\omega)=\frac{1}{T}\int_{t_0}^{t_0+T} s(t)e^{-i\omega t}dt), (T=N_m 2\pi/\omega)
  (note: for (A\cos\omega t), this returns (A/2); cancels in (Z^{\rm eff}))
* solve policy: direct solve when well-conditioned; else pseudoinverse / ridge (see (EstimateZeffRobust))

---

## L) Regression limit (freeze the “known-answer” setting)

Define the parameter limit where the model should reduce to the earlier cavity/selector-like behavior:

* Family 1: (\Omega_{\rm in}\to \Omega_{\rm out}) large (suppresses (w)-excitations)
* Family 2: (\gamma\to\infty) (suppresses (w) in transverse radius)
* Walls steep: (V_0) large, (\delta) small, (p) moderate/large
* Wave drive off (unless specifically testing wave modes)

This is a regression test setting, not a claim about reality.

---

## M) Verification suite (must pass before results)

Symbolic harness: `/mnt/data/master_checks.wl`
Required checks to pass:

* GNLS EOM, continuity, thermo identities
* dynamic continuity with explicitly time-dependent real (V_{\rm conf}(\mathbf X,t)) from (a(t),L(t))
* fluid force identity for (a) **and (L)**
* wave EOM, wave force identity (if wave track active)
* moving-wall work/energy balance with explicit (a(t),L(t))
* wave moving-wall work/energy balance with (\mu_A^2(\mathbf X;a(t),L(t))) (if wave track active)
* total-force additivity guard: (-\partial_{a,L}H_{\rm tot}) equals sum of sector forces (if wave track active)
* localization unit tests: (\partial_a V_{\rm conf}) wall example and (\partial_L V_{\rm conf}) endcap test with tolerance-based central leakage criterion due to SmoothAbs + tanh endcaps:
  [
  \frac{|(\partial_L V_{\rm conf})(w=0)|}{\max_{w\in[-L,L]} |(\partial_L V_{\rm conf})(w)|} < \epsilon
  \quad\text{under}\quad \varepsilon_w \ll \delta_\parallel \ll L.
  ]

---

## “Freeze Signatures” (fill these in once you commit)

* Primary (V_{\rm conf}) family: [x] F1 [ ] F2
* Axis definition (s(\mathbf X)): (s=w)
* Transverse radius (r(\mathbf X)): (r=\sqrt{x^2+y^2+z^2})
* (E_{\rm geom}) terms used: [x] (P_{\rm vac}) [x] (\sigma) [ ] (\kappa_b)
* Constraint: [x] fixed (N) [ ] fixed (m_{\rm eff})
* Projection weight (W(w)): [x] (|\chi_0|^2) [ ] other
* Measurement region (\Gamma): [x] sphere [ ] disk [ ] slab
* Surface measure/normal: [x] (d\mu=r_{\rm port}^2\sin\theta\,d\theta\,d\phi; \hat n=\hat R_3)
* Quadrature: [x] \theta GL (n_\theta=32), [x] \phi uniform (n_\phi=64)
* Flux definition (j): [x] mouth [x] leakage [x] both
* Port basis: [x] (Y_{lm}, l<=l_{\max}, ordering=PortsList) [ ] other
* Drive protocol: [x] Vconf -> Vconf + eps cos(ωt) exp(-w^2/wDrive^2) exp(-(R3-rPort)^2/rDrive^2) Y_lm(θ(x,y,z),φ(x,y))

---

If you want, I can also turn this into a **literal one-page LaTeX “Definition Freeze” appendix** that you can drop into Paper 7 verbatim (same content, but formatted as a compact checklist + equations).
