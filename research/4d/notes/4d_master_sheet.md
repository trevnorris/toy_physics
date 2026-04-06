## 0) Coordinates, indices, operators

We work in **4D space + time** with
[
\mathbf X=(x,y,z,w),\qquad t\in\mathbb R,
\qquad
\nabla_4=(\partial_x,\partial_y,\partial_z,\partial_w),
\qquad
\nabla_4^2=\partial_x^2+\partial_y^2+\partial_z^2+\partial_w^2.
]
Spatial indices (i,j\in{x,y,z,w}). Spacetime indices (M,N\in{0,x,y,z,w}) with (0\equiv t).

The **brane** is operationally the neighborhood of (w=0), and “brane observables” are obtained by a frozen projection map (below).

---

## 1) Fixed EOS and thermodynamics (non-negotiable)

You froze the “stiff” polytrope
[
P(\rho)=K\rho^5.
]
This implies the required thermodynamic closures
[
U(\rho)=\frac{K}{4}\rho^5,
\qquad
h(\rho)=U'(\rho)=\frac{5K}{4}\rho^4,
\qquad
P(\rho)=\rho h(\rho)-U(\rho)=K\rho^5,
\qquad
c_s^2(\rho)=\frac{dP}{d\rho}=5K\rho^4.
]
These identities were explicitly verified by the symbolic harness as consistent with the GNLS nonlinearity.

---

## 2) Fundamental dynamical variables (frozen set)

### 2.1 Matter/superfluid field

Complex order parameter:
[
\psi(\mathbf X,t)\in\mathbb C,\qquad \rho(\mathbf X,t)=|\psi|^2.
]
Canonical 4D mass/number current:
[
\mathbf J(\mathbf X,t)=\frac{\hbar}{2im}\left(\psi^*\nabla_4\psi-\psi\nabla_4\psi^*\right).
]
Exact continuity identity:
[
\partial_t\rho+\nabla_4\cdot\mathbf J=0.
]
All three were verified.

### 2.2 Geometry degrees of freedom

Reduced throat geometry:
[
a(t)>0,\qquad L(t)>0,
]
with shape-field upgrade deferred (frozen off for now).

### 2.3 True 4+1D gauge sector (required for “real EM”)

The unified program commits to a genuine 4+1D gauge field (A_M(\mathbf X,t)) (not “EM-from-fluid” as the defining model once C1 is adopted).

Field strength:
[
F_{MN}=\partial_M A_N-\partial_N A_M.
]


---

## 3) Frozen coupling: C1 minimal coupling + neutrality + localization

### 3.1 Gauge-covariant derivatives (C1)

Minimal coupling means:
[
D_t\psi=\left(\partial_t+\frac{i q}{\hbar}A_0\right)\psi,\qquad
D_i\psi=\left(\partial_i-\frac{i q}{\hbar}A_i\right)\psi,\quad i\in{x,y,z,w}.
]

### 3.2 Charge neutrality (your freeze: jellium)

You froze the background-subtracted charge density:
[
\boxed{
J^0(\mathbf X)=q\left(|\psi(\mathbf X,t)|^2-\rho_0\right).
}
]
Spatial charge current is (single-field version)
[
\boxed{
J^i(\mathbf X,t)=q,J_{\rm mass}^i(\mathbf X,t),
\qquad
J_{\rm mass}^i=\frac{\hbar}{m}\Im\left(\psi^* D^i\psi\right).
}
]
This is the current used as the source in Maxwell.

### 3.3 Gauge localization (your freeze: dielectric confinement (Z(w)))

You froze the kinetic prefactor localization strategy (“LZ”):
[
\boxed{
\mathcal L_{\rm EM}=
-\frac{1}{4\mu_0},Z(w),F_{MN}F^{MN},
\qquad
Z(w)=\exp!\left(-\frac{w^2}{\lambda_{\rm conf}^2}\right).
}
]
This is exactly one of the variationally well-posed localization mechanisms in the unified plan.

**Important note (frozen requirement, not yet “proven true”):** even with (Z(w)), the brane must be checked to actually see the intended (1/r^2) scaling in projected observables.

---

## 4) Frozen confinement potential family (F1) and exact shapes

### 4.1 F1 “modulated brane trap” ontology

Geometry is encoded **as a smooth potential** (V_{\rm conf}(\mathbf X;a,L)), not as hard boundary conditions, so that (\partial_a V_{\rm conf}) and (\partial_L V_{\rm conf}) exist and define forces.

### 4.2 Your specific shape freeze

You froze:

* **Bulk/brane trap:** harmonic in (w),
  [
  V_w(w)=\frac12 m\omega_z^2 w^2.
  ]
* **Radial throat envelope:** super-Gaussian (order 4),
  [
  E(R_3;a)=\exp!\left(-\frac{R_3^4}{a^4}\right),
  \qquad R_3=\sqrt{x^2+y^2+z^2}.
  ]

### 4.3 A concrete explicit (V_{\rm conf}) consistent with both the freeze sheet + your shape choice

To remain in the “Family-1 gate” framework (so (a,L) really exist as geometry DOFs), define a smooth interior gate
[
G(\mathbf X;a,L)=G_r(R_3;a),G_w(w;L),
]
with (G_r) taken from your super-Gaussian choice
[
G_r(R_3;a)=E(R_3;a)=\exp!\left(-\frac{R_3^4}{a^4}\right),
]
and (G_w) a smooth endcap gate (the hard-mode baseline uses tanh + SmoothAbs so (\partial_L) is well-behaved):
[
S(u)=\tfrac12(1+\tanh u),\qquad
{\rm SmoothAbs}(w)=\sqrt{w^2+\varepsilon_w^2},
]
[
G_w(w;L)=1-S!\left(\frac{{\rm SmoothAbs}(w)-L/2}{\delta_\parallel}\right),
\qquad
(\varepsilon_w\ll\delta_\parallel\ll L).
]

Then a clean “modulated trap” is:
[
\boxed{
V_{\rm conf}(\mathbf X;a,L)=\frac12 m\Big(\Omega_{\rm out}^2-\big(\Omega_{\rm out}^2-\Omega_{\rm in}^2\big),G(\mathbf X;a,L)\Big),w^2.
}
]

* “Stiff far field”: (\Omega_{\rm out}) large.
* “Locally weakened throat”: (\Omega_{\rm in}\ll \Omega_{\rm out}) (can be (0) if you want fully open corridor).

This is the direct, differentiable encoding of your physical statement “steep far-field confinement but locally weakened inside the throat.”

### 4.4 Sponge boundary condition (your freeze)

Matter-sector absorbing layer (CAP) is frozen as:
[
\boxed{
V_{\rm sim}(\mathbf X,t)=V_{\rm conf}(\mathbf X;a,L)-i,\Gamma_{\rm sponge}(\mathbf X).
}
]

---

## 5) Geometry energy and closure law (your freeze)

### 5.1 Geometry energy: surface term

You froze the **contractile closing force** as a surface-tension-like penalty:
[
\boxed{
E_{\rm geom}(a,L)=\sigma,\mathcal A(a,L)
}
]
(optionally with a small curvature term only as numerical regularization).

A minimal explicit “tube with endcaps” 4D hyperarea model (useful for code and paper bookkeeping) is:
[
\mathcal A(a,L)\approx 4\pi a^2 L+\frac{8\pi}{3}a^3,
]
so
[
\partial_a E_{\rm geom}=\sigma\left(8\pi a L+8\pi a^2\right),\qquad
\partial_L E_{\rm geom}=\sigma\left(4\pi a^2\right).
]
(If you later choose a different (\mathcal A(a,L)), that’s a single swap — the rest of the ledger stays intact.)

### 5.2 Closure law: dynamical “breathing” wall

You froze **dynamic closure** (not static balance). Geometry evolves via ODEs driven by the total generalized forces:
[
\boxed{
M_a\ddot a + C_a\dot a = F_a(t),
\qquad
M_L\ddot L + C_L\dot L = F_L(t),
}
]
with
[
\boxed{
F_a(t)=-\frac{\partial H_{\rm tot}}{\partial a},
\qquad
F_L(t)=-\frac{\partial H_{\rm tot}}{\partial L}.
}
]

---

## 6) Core coupled PDE system (the “full 4D model”)

### 6.1 Matter PDE: 4D stiff GNLS (verified backbone)

The symbolic harness verified that the GNLS Lagrangian produces:
[
\boxed{
i\hbar,\partial_t\psi=
\left[
-\frac{\hbar^2}{2m}\nabla_4^2
+V_{\rm conf}(\mathbf X;a,L)
+\frac{5K}{4}|\psi|^8
\right]\psi.
}
]

### 6.2 Gauged matter PDE (required once EM is “real EM”)

With C1 minimal coupling, the natural upgrade is:
[
\boxed{
i\hbar,D_t\psi=
\left[
-\frac{\hbar^2}{2m},D_iD_i
+V_{\rm sim}(\mathbf X;a,L)
+\frac{5K}{4}|\psi|^8
\right]\psi
;+;\eta_{\rm drive}(\mathbf X,t).
}
]
Where:

* (V_{\rm sim}=V_{\rm conf}-i\Gamma_{\rm sponge}) (frozen),
* (\eta_{\rm drive}) is the **wave-drive injection slot** (your drive protocol choice; its exact form is frozen operationally, not derived).

The existence of an explicit “drive slot” is part of the unified plan’s open-system stance.

### 6.3 Maxwell PDE with localization (required)

With your frozen (Z(w)) localization, the Maxwell equation bundle is:
[
\boxed{
\partial_M!\left(Z(w),F^{MN}\right)
+\frac{1}{\xi},\partial^N(\partial!\cdot!A)
= \mu_0,(J^N+J^N_{\rm ext}).
}
]

* (J^N) is the matter current defined above (with jellium in (J^0)).
* (J_{\rm ext}^N) is the **EM drive slot** (if you choose to pump cavity modes directly).
* (\xi) is the gauge-fixing parameter (Lorenz gauge family; (\xi=1) is the simplest baseline). The *value* of (\xi) is one of the few small items still worth explicitly freezing (see “Remaining items” at the end).

---

## 7) Madelung / hydrodynamic form (for derivations and “what falls out”)

This is the interpretation layer your scripts were exploring (and where Poisson/extra sectors can emerge *as limits*, not as imposed laws).

### 7.1 Madelung transform

[
\psi=\sqrt{\rho},e^{i\theta}.
]

Gauge-invariant 4D superfluid velocity:
[
\boxed{
v_i=\frac{\hbar}{m}\partial_i\theta-\frac{q}{m}A_i.
}
]

### 7.2 Continuity (bulk)

[
\boxed{
\partial_t\rho+\partial_i(\rho v_i)=0,
\qquad i\in{x,y,z,w}.
}
]
(This is the same identity as (\partial_t\rho+\nabla_4\cdot\mathbf J=0), rewritten.)

### 7.3 Euler / Hamilton–Jacobi form (bulk)

A standard gauged-Madelung Euler form is:
[
\boxed{
m\left(\partial_t v_i+v_j\partial_j v_i\right)
==============================================

-\partial_i!\Big(V_{\rm conf}+h(\rho)+Q\Big)
;+;q\Big(F_{i0}+v_jF_{ij}\Big),
}
]
with enthalpy (h(\rho)=\tfrac{5K}{4}\rho^4) (fixed) and quantum pressure
[
Q(\rho)= -\frac{\hbar^2}{2m}\frac{\nabla_4^2\sqrt{\rho}}{\sqrt{\rho}}.
]
Your script’s “Euler-like equation (component form, 4D space)” is exactly this structure expanded componentwise (including the stiff enthalpy gradient and the quantum-pressure third-derivative terms).

---

## 8) Brane projection and operational observables (frozen)

### 8.1 Projection operator (your freeze: ground-state weighting)

You froze weighted brane projection:
[
\boxed{
\mathcal P_W[f](x,y,z,t)=\int W(w),f(x,y,z,w,t),dw,
\qquad
W(w)=|\chi_0(w)|^2,
}
]
implemented with a finite window (|w|\le W_{\rm proj}) and renormalization:
[
{\cal N}*W=\int*{-W_{\rm proj}}^{W_{\rm proj}} W(w),dw,\qquad
\widetilde W(w)=W(w)/{\cal N}_W.
]
Tail mass (1-{\cal N}_W) is the projection-validity diagnostic.

### 8.2 Brane density and current

[
\boxed{
\rho_{\rm brane}(x,y,z,t)=\mathcal P_W[\rho],
\qquad
\mathbf J_{\rm brane}(x,y,z,t)=\mathcal P_W[\mathbf J_{xyz}],
}
]
(and the brane “velocity” is (\mathbf v_{\rm brane}=\mathbf J_{\rm brane}/\rho_{\rm brane}) when defined).

### 8.3 Brane-observed EM fields

From (F_{MN}), the brane-observed electric-type components are
[
\boxed{
E_i^{\rm brane}(x,y,z,t)=\mathcal P_W[F_{0i}](x,y,z,t),\qquad i\in{x,y,z}.
}
]
Magnetic-type components come from projected (F_{ij}) restricted to (i,j\in{x,y,z}) and mapped to the usual 3D pseudovector.

### 8.4 Measurement surface (your freeze: brane 2-sphere)

You froze:
[
\boxed{
\Gamma:\ R_3=r_{\rm port},\ w=0,
\qquad
d\mu=r_{\rm port}^2\sin\theta,d\theta,d\phi,
\qquad
\hat n=\hat R_3.
}
]

### 8.5 Ports and modal observables (frozen form)

With port basis (P_i(\theta,\phi)) (spherical harmonics (Y_{\ell m}) in the freeze sheet):
[
\boxed{
u_i(t)=\int_\Gamma \overline{P_i(\mathbf s)},u(\mathbf s,t),d\mu,
\qquad
j_i(t)=\int_\Gamma \overline{P_i(\mathbf s)},j(\mathbf s,t),d\mu.
}
]

Effort variable (frozen enthalpy perturbation proxy):
[
u=\delta h(\rho_{\rm brane})
\approx h'(\rho_{\rm brane,0}),\delta\rho_{\rm brane}
=====================================================

5K\rho_{\rm brane,0}^3,\delta\rho_{\rm brane}.
]

### 8.6 Flux definition (your freeze: multi-output mouth + leakage)

**Mouth flux:**
[
\boxed{
j_{\rm mouth}(\mathbf s,t)=\mathbf J_{\rm brane}\cdot\hat n\quad\text{on }\Gamma.
}
]

**Leakage flux into bulk (monitor (J_w)):**
[
\boxed{
J_w=\frac{\hbar}{2im}\left(\psi^*\partial_w\psi-\psi,\partial_w\psi^*\right),
}
]
[
\boxed{
j^{\rm out}*+(t)=\int*{\Omega_{xyz}} J_w\big|*{w=+W*{\rm cut}},d^3x,
\quad
j^{\rm out}*-(t)=\int*{\Omega_{xyz}} \big[-J_w\big|*{w=-W*{\rm cut}}\big],d^3x,
\quad
j_{\rm leak}=j^{\rm out}*+ + j^{\rm out}*-.
}
]

This is *exactly* the diagnostic structure you want if you’re investigating the “refill/leak” idea.

---

## 9) Force ledger identities (what the scripts “locked in”)

These are the pieces that make the geometry + PDE system **honest** (no forced stabilization).

### 9.1 Hellmann–Feynman geometry force identity (verified)

When geometry only enters via (V_{\rm conf}(\mathbf X;a,L)),
[
\boxed{
F^{(\psi)}_{a,L}
================

# -\frac{\partial H_\psi}{\partial(a,L)}

-\int d^4X\ \rho(\mathbf X,t),\partial_{a,L}V_{\rm conf}(\mathbf X;a,L).
}
]
Verified.

### 9.2 Moving-wall work/energy identity (verified)

For time-dependent confinement,
[
\boxed{
\partial_t\mathcal H+\nabla_4\cdot\mathbf S=\rho,\partial_t V_{\rm conf}.
}
]
This is the bookkeeping gate that makes “breathing walls” meaningful.

---

## 10) What “falls out” without forcing: Poisson sector + extra sectors

### 10.1 Poisson is a *limit identity* (not an imposed law)

From the earlier 3D emergence work: in the **linearized / low-Mach / predominantly potential-flow regime**, a velocity potential extracted from the longitudinal velocity satisfies a Poisson equation in the quasi-static limit,
[
\nabla^2\psi = \frac{S_\rho}{\rho_0}-\left\langle \frac{S_\rho}{\rho_0}\right\rangle
]
(with the mean subtraction required by periodic solvability). This was explicitly framed as “Poisson shows up as a diagnostic identity of the quasi-steady limit,” not as a forced evolution law.

In the *4D* program, the correct stance is:

* the full system is GNLS+Maxwell+geometry,
* Poisson-like behavior is something we test for in an appropriate projected/low-frequency regime,
* and it may come with additional terms (below).

### 10.2 The “extra” sectors that appear naturally in the projected theory

Once you define brane observables via projection (\mathcal P_W), you do **not** get a closed 3D perfect-fluid system automatically. Two big unavoidable “extra channels” are:

1. **Leakage source terms in the brane continuity**, because (w)-flux and weight gradients don’t vanish in general:
   [
   \partial_t\rho_{\rm brane}+\nabla_3\cdot\mathbf J_{\rm brane}
   =============================================================

   \text{(terms involving (J_w), (W'(w)), and cutoff fluxes)}.
   ]
   This is exactly the mathematical doorway for “refill/leak” physics — but it must be measured, not assumed.

2. **Reynolds/stress corrections in the brane momentum equation**, because projection does not commute with nonlinear products:
   [
   \mathcal P_W[\rho v_i v_j] \neq
   \rho_{\rm brane} v^{\rm brane}_i v^{\rm brane}*j,
   ]
   which yields an emergent stress (R*{ij}) source. This is one robust mechanism by which the brane can show **effective vorticity / transverse structure** even if the bulk flow is “mostly potential.”

Similarly, in the EM sector, brane observables see a projection of (F_{MN}), and **(w)-components** like (F_{iw}), (F_{0w}) can exist and act like “hidden channels” that only appear indirectly in 3D measurements.

This is exactly the kind of “extra structure” you said you’re curious about — and the model has a mathematically clean place for it to appear without being inserted by hand.

---

## 11) Remaining “small freezes” / items to explicitly note as not-yet-final

These are not alternate theories — they’re just a few remaining parameters/conventions that must be pinned down so a paper is airtight:

1. **Gauge-fixing choice**: pick (\xi) (start with (\xi=1) unless numerics force otherwise). 
2. **Exact implementation of the wave drive** (your chosen protocol): whether (\eta_{\rm drive}) drives (\psi) directly, or (J_{\rm ext}) drives the Maxwell sector, or both (but with one canonical “paper definition” frozen).
3. **Precise leakage integration region** (\Omega_{xyz}), and the cut locations (W_{\rm cut}), (R_{\rm measure}) (you already have a baseline form; it just needs the final parameter values).
4. **Numerical gauge-consistent discretization choice**: explicit covariant derivative vs link-variable method (the plan notes both; pick one for v1 and keep gauge-residual diagnostics). 
5. **Chemical potential implementation (fixed (\mu))**: how the “infinite ocean reservoir” is realized in the PDE (e.g., boundary reservoir, relaxation term, or initialization + driving). This is conceptually frozen by you; the exact term-level implementation should be frozen in code/paper once you pick the cleanest.

---

If you want, I can also convert the “Master sheet” above into a **compact LaTeX appendix** (one or two pages) in the style of the freeze sheet, so it can literally be dropped into the paper as “Model Definition (Frozen).”

