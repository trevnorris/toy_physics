## 0. Executive Summary

### 0.1 What changed

We are no longer building Paper 7 around a stitched “3D exterior + 4D interior” picture (or around an assumed inner impedance). The new direction is:

* **One unified 4D model** on (\mathbb{R}^4) with coordinates (\mathbf{X}=(x,y,z,w)).
* “Brane physics” is not a separate domain; it **emerges** because the field is strongly confined near (w=0) outside the throat region.
* The throat is not a boundary object; it is created by a **geometry-as-potential** construction (V_{\rm conf}(\mathbf{X};a,L,\dots)) that locally relaxes confinement and/or opens a 4D corridor.

This is the “correct-first, hard-mode” commitment: if the interesting behavior is genuinely 4D, we want it to appear naturally rather than be averaged away by design.

---

### 0.2 The core objective

Paper 7’s core job is to replace “assume a throat impedance” with a **governing interior+geometry model** that:

1. **Defines** the throat (geometry) in 4D via potentials and a geometry energy,
2. **Evolves or equilibrates** geometry DOFs (at minimum (a,L)) using a wall/closure law derived from the same total energy,
3. Produces an **emergent effective mouth response operator** (Z^{\rm eff}(\omega)) measurable on/near the brane (not imposed as a boundary condition),
4. Is algebraically and energetically consistent (conservation laws, work balance).

---

### 0.3 Two independent “support” mechanisms (separate experiments)

We are explicitly running **two separate physical experiments** for how a particle throat can stay open, and comparing results:

* **Experiment F: Fluid-stress supported throat**
  Opening is maintained primarily by superfluid pressure and momentum flux (through-flow / intake). No reliance on standing wave support.

* **Experiment W: Wave-stress supported throat**
  Opening is maintained by stored standing-wave energy (mode pressure / radiation stress), modeled as an **independent wave field** in 4D.

They can later be combined, but Paper 7 must first establish each in isolation so we can identify failure modes cleanly.

---

### 0.4 Scope note: particles now, 2PN/black-hole regime later

This paper is **particle-focused**. We are not trying to be 2PN-correct at this stage and we are not assuming black holes are wave-supported. The plan is:

* Build a correct, self-consistent 4D particle throat model first.
* Later, compare the resulting effective long-range behavior against PN/EFT expectations and adjust fundamentals if needed.

---

### 0.5 Deliverables (what Paper 7 must output)

Paper 7 is “done” when we can produce:

1. **A frozen master model**: explicit (V_{\rm conf}(\mathbf X;a,L)), explicit (E_{\rm geom}(a,L)), and the chosen field equations.
2. **A closure law**: equilibrium conditions (\partial_a H_{\rm tot}=0), (\partial_L H_{\rm tot}=0) (and/or dynamics (M_a\ddot a=-\partial_a H_{\rm tot})).
3. **An emergent effective operator** (Z^{\rm eff}_{ij}(\omega)) defined by a drive/measure protocol on the brane (ports).
4. **A failure map**: parameter regimes where equilibrium does not exist, is unstable, becomes resonance-dominated/nonlocal, or projection to brane observables breaks down.
5. **Mathematical validation**: symbolic checks confirming the PDEs, conservation laws, and force/work identities are consistent (already largely achieved via the Mathematica harness).

---

## 1. Hard-Mode 4D Philosophy and Non-Negotiables

### 1.1 Why full 4D

We are choosing hard-mode 4D because we explicitly want to allow:

* bulk interactions in (w) that are invisible to a purely 3D model,
* intermittent/“sporadic” 3D behavior that could be the shadow of 4D mode coupling,
* future extensions (e.g., reconnection/MHD-like phenomena) where extra-dimensional structure may create qualitatively new pathways.

Hard-mode is not “more complicated for its own sake”; it’s a commitment that we will not build in assumptions that prevent genuinely 4D effects from existing.

---

### 1.2 One PDE in 4D everywhere (no stitched interfaces)

We do **not** want a model where:

* “inside throat” is solved in 4D,
* “outside” is solved separately in 3D,
* and the two are matched by a boundary condition that implicitly assumes what should be derived.

Instead we solve a single 4D field problem and define “brane observables” as projections or measurements localized near (w=0).

---

### 1.3 The brane emerges via confinement in (w)

The “brane” is implemented as a strong confinement mechanism in the extra coordinate:

* Far from the particle/throat: the confining potential makes excitations in (w) expensive, so the dynamics effectively collapse to 3D behavior.
* Near the throat: the confinement is relaxed/opened in a controlled region so the field can explore 4D.

Operationally:
[
V_{\rm conf}(\mathbf X;a,L)=V_{\rm brane}(w;\text{trap})+V_{\rm throat}(\mathbf X;a,L;\text{shape}),
]
where the “mouth” is not a boundary—it is a **transition zone** in which the effective confinement changes.

---

### 1.4 Microphysics constraint: the (n=5) EOS is fixed

A non-negotiable continuity constraint from the established toy-universe calibration is the stiff polytrope:
[
P(\rho)=K\rho^5,\qquad
c_s^2(\rho)=\frac{dP}{d\rho}=5K\rho^4.
]
In GNLS form this corresponds to:
[
h(\rho)=\frac{5K}{4}\rho^4,\qquad U(\rho)=\frac{K}{4}\rho^5,
]
and therefore a nonlinear term (h(|\psi|^2)\propto |\psi|^8) in the 4D equation.

This is not a stylistic choice: it’s the anchor that keeps Paper 7 in the same “physics universe” as the hybrid/EOS foundations.

---

### 1.5 Wave support must be modeled independently (for comparison)

We will not bury wave-support inside fluid variables when we are trying to decide whether wave-support is a viable mechanism. For Paper 7’s decision-making, wave support is treated as an explicit second field:

* superfluid field (\psi(\mathbf X,t)) (GNLS or hydro),
* independent wave field (A(\mathbf X,t)) (scalar surrogate initially, Maxwell upgrade later if needed).

This lets us compare:

* equilibria that exist only with wave support,
* equilibria that exist with fluid-only support,
* and the differences in emergent (Z^{\rm eff}(\omega)) and stability.

---

### 1.6 Continuity requirement: preserve the EM/cavity selector as a unit-test limit

Even though we are going full 4D, the model must admit a controlled limit in which the earlier “EM/cavity” selector branch can be recovered (e.g., by enforcing strong confinement that suppresses (w)-structure and reproduces the cylindrical Bessel-root behavior).

This is not because that limit is “the truth” for all objects—it’s because Paper 7 must have at least one **known-answer regression test** that prevents branch drift and silent convention changes.

---

### 1.7 Correctness gates already in place (symbolic verification)

Before we do heavy derivations or numerics, we require internal consistency:

* EOM derived from the stated Lagrangian/Hamiltonian,
* continuity/current identity,
* breathing force identities from (-\partial_{a,L}H),
* moving-wall work/energy balance,
* thermodynamic consistency of (U\to h\to P\to c_s^2).

The current Mathematica harness already demonstrates these checks for the master GNLS + wave-field skeleton, which means the next work is about **freezing the explicit potentials and geometry measures**, not re-litigating algebra.

---

## 2. Geometry Without Boundaries: Potentials as Geometry

### 2.1 Coordinates and core notation

We work on full 4D space with coordinates
[
\mathbf{X}=(x,y,z,w), \qquad \nabla_4=(\partial_x,\partial_y,\partial_z,\partial_w).
]
Useful radii:

* brane 3D radius: (\displaystyle R_3=\sqrt{x^2+y^2+z^2})
* transverse 2D radius: (\displaystyle r_\perp=\sqrt{x^2+y^2}) (used only for specific “cylindrical branch” limits)
* optional anisotropic 4D transverse radius (Family 2): (\displaystyle R_\gamma=\sqrt{r_\perp^2+(w/\gamma)^2})

We represent “geometry” not as a hard boundary but as a smooth confinement potential with parameters:

* (a): throat “radius” (effective transverse size),
* (L): throat “length” (extent along a chosen axis in 4D),
* “shape” parameters: flare, rounding, taper, smoothness widths (\delta), steepness exponent (p), etc.

---

### 2.2 One governing idea: geometry = confinement potential

Instead of defining a domain (\mathcal{T}) with boundary conditions, we define a single smooth potential:
[
V_{\rm conf}(\mathbf{X};a,L,\text{shape})=V_{\rm brane}(w;\text{trap})+V_{\rm throat}(\mathbf{X};a,L,\text{shape}).
]

Interpretation:

* **Brane region:** (V_{\rm brane}(w)) strongly traps the system near (w=0), so excitations in (w) are suppressed and the physics looks effectively 3D.
* **Throat region:** (V_{\rm throat}) locally weakens that trap and/or creates a low-energy corridor into the bulk so the field explores 4D.

The “mouth” is the **transition zone** where (V_{\rm conf}) changes character (not a boundary surface).

---

### 2.3 Smooth gating functions (so forces are well-defined)

We use smooth step/gate functions so that derivatives like (\partial_a V_{\rm conf}) and (\partial_L V_{\rm conf}) exist and are numerically stable. A typical smooth step is:
[
S(u)=\tfrac12\big(1+\tanh u\big),
]
and a typical “inside-throat gate” is a product of radial and axial gates:
[
G(\mathbf X)=G_\perp(\mathbf X),G_\parallel(\mathbf X),
]
with (G\approx 1) inside the throat core and (G\approx 0) outside.

**Reason this matters:** the wall/breathing forces are obtained from
[
F_a=-\frac{\partial H}{\partial a},\qquad F_L=-\frac{\partial H}{\partial L},
]
and these are volume integrals involving (\partial_a V_{\rm conf}), (\partial_L V_{\rm conf}). Smooth gates prevent spurious singular forces and make the “wall law” unambiguous.

---

### 2.4 Potential Family 1: Modulated brane trap (4D access by relaxing (w)-confinement)

**Core idea:** keep a strong brane trap in (w) almost everywhere, but weaken it in a localized 4D region that defines the throat.

A canonical form:
[
V_{\rm brane}(w;\mathbf X)=\frac12 m,\Omega_w^2(\mathbf X),w^2,
]
with
[
\Omega_w^2(\mathbf X)=\Omega_{\rm out}^2-\big(\Omega_{\rm out}^2-\Omega_{\rm in}^2\big),G(\mathbf X),
\quad \Omega_{\rm out}\gg \Omega_{\rm in}\ge 0.
]

Then choose a throat gate (G(\mathbf X)) that turns on inside a tube-like region characterized by ((a,L)) (the exact tube orientation can be chosen later, e.g., along (z) or along (w), depending on the object class you’re modeling).

**What this buys us:**

* A clean dial for “how 4D is it”: (\Omega_{\rm in}/\Omega_{\rm out}).
* A controlled limit where the system becomes almost purely brane-like (take (\Omega_{\rm in}\to\Omega_{\rm out})).

---

### 2.5 Potential Family 2: True 4D transverse tube (4D cross-section physics)

**Core idea:** define the tube wall directly as a soft wall in an explicitly 4D transverse radius (or anisotropic variant).

Using (R_\gamma):
[
R_\gamma=\sqrt{r_\perp^2+\left(\frac{w}{\gamma}\right)^2}.
]
Then a tube wall can be represented by:
[
V_{\rm wall}(\mathbf X;a,L)=V_0,S!\Big(\frac{R_\gamma-a(\cdot)}{\delta}\Big)^p + \text{(endcap terms setting }L\text{)}.
]

A separate (V_{\rm brane}(w)) is still typically included to ensure the far field remains brane-trapped outside the mouth region.

**What this buys us:**

* Explicitly 4D transverse eigenstructure.
* A clean “EM/cavity limit knob”: (\gamma\to\infty) collapses (R_\gamma\to r_\perp), suppressing (w) as a transverse contributor.

---

### 2.6 Geometry parameters and “what is (L) in hard-mode 4D?”

In hard-mode 4D, “length” is not automatically tied to a single spatial coordinate unless we decide it is. We will interpret (L) as the characteristic extent of the low-potential corridor along a chosen axis (or along the centerline of a funnel). Operationally, (L) enters through the axial/endcap gating in (V_{\rm throat}).

This must be frozen before implementation because (\partial_L V_{\rm conf}) defines the generalized “endcap force.”

---

### 2.7 Required derivatives (the minimum geometry math we must lock)

For the wall law and breathing dynamics we require explicit expressions for:
[
\partial_a V_{\rm conf}(\mathbf X;a,L),\qquad \partial_L V_{\rm conf}(\mathbf X;a,L).
]

In practice:

* (\partial_a V_{\rm conf}) should be localized to the tube wall region (where the gate transitions),
* (\partial_L V_{\rm conf}) should be localized to the endcap/termination transition region, with any center leakage parametrically small.

These localization properties are sanity checks: if (\partial_a V) has support far away from the particle, the geometry definition is wrong.
With SmoothAbs + tanh endcaps, (\partial_L V_{\rm conf}) will not be exactly zero at (w=0); the unit test should therefore be tolerance-based, e.g., enforce that the center value is parametrically small relative to the endcap peak under (\varepsilon_w \ll d_w \ll L). One explicit criterion is:
[
\frac{|(\partial_L V_{\rm conf})(w=0)|}{\max_{w\in[-L,L]} |(\partial_L V_{\rm conf})(w)|} < \epsilon
\quad\text{under}\quad \varepsilon_w \ll d_w \ll L.
]

---

### 2.8 What is frozen vs what is still a choice (end of geometry section)

**Frozen (conceptual):**

* Geometry is encoded as (V_{\rm conf}=V_{\rm brane}+V_{\rm throat}), not boundary conditions.
* We require smooth gates so (\partial_{a,L}V_{\rm conf}) exist and define forces.

**Frozen baseline (Paper 7 fluid-only track):**

* **Baseline family:** Family 1 (Modulated Brane Trap + Soft Walls).
* **Axis choice:** (s=w) with (r=\sqrt{x^2+y^2+z^2}) as the transverse radius.
* **Gate choice:** tanh-based smooth gates with **SmoothAbs** for (|w|) in the axial/endcap gate.
* **Endcap barrier:** explicitly present (mitigates drain/leak risk at the mouth).

**Still to choose (remaining details):**

* The concrete functional forms of (G(\mathbf X)), (a(\cdot)) (flare), and the numerical stiffness parameters ((V_0,\delta,p)).

---

## 3. Superfluid Core Model in 4D (n=5 EOS)

### 3.1 Non-negotiable thermodynamics (the EOS anchor)

We adopt the stiff polytrope:
[
P(\rho)=K\rho^5,
]
which fixes:
[
c_s^2(\rho)=\frac{dP}{d\rho}=5K\rho^4,
\qquad
h(\rho)=\int^\rho \frac{dP}{\rho'}=\frac{5K}{4}\rho^4,
\qquad
U(\rho)=\int^\rho h(\rho'),d\rho'=\frac{K}{4}\rho^5.
]

These identities are structural: if we change them, we are no longer in the same calibrated toy-universe regime.

---

### 3.2 Primary PDE: 4D GNLS (recommended backbone)

Let (\psi(\mathbf X,t)) be the complex order parameter and define density (\rho=|\psi|^2). The 4D GNLS is:
[
i\hbar,\partial_t \psi
======================

\left[
-\frac{\hbar^2}{2m}\nabla_4^2
+V_{\rm conf}(\mathbf X;a,L)
+h(|\psi|^2)
\right]\psi,
]
with
[
h(|\psi|^2)=\frac{5K}{4}|\psi|^8.
]

So explicitly:
[
i\hbar,\partial_t \psi
======================

\left[
-\frac{\hbar^2}{2m}\nabla_4^2
+V_{\rm conf}(\mathbf X;a,L)
+\frac{5K}{4}|\psi|^8
\right]\psi.
]

**Why GNLS is the default:** it is conservative, has a clean energy functional, yields continuity + current automatically, and provides a built-in dispersive regularization (important near steep walls).

---

### 3.3 Energy functional (for forces, equilibrium, and bookkeeping)

The fluid energy (Hamiltonian) is:
[
H_{\rm fluid}[\psi;a,L]=
\int_{\mathbb R^4}
\left[
\frac{\hbar^2}{2m}|\nabla_4\psi|^2
+V_{\rm conf}|\psi|^2
+\frac{K}{4}|\psi|^{10}
\right],d^4X.
]

This is central because:

* the PDE can be stated as (i\hbar,\partial_t\psi=\delta H_{\rm fluid}/\delta\psi^*),
* breathing forces are (F_a=-\partial_a H_{\rm fluid}), (F_L=-\partial_L H_{\rm fluid}),
* energy/work balance with a moving wall comes from the same Hamiltonian structure.

---

### 3.4 Continuity and current (guaranteed structure)

Define density (\rho=\psi^*\psi). The GNLS implies:
[
\partial_t\rho+\nabla_4\cdot \mathbf J=0,
]
with current
[
\mathbf J=\frac{\hbar}{2im}\left(\psi^*\nabla_4\psi-\psi\nabla_4\psi^*\right).
]

This identity is not optional; it’s the conservation law that makes “flux through a mouth region” a meaningful observable later.

---

### 3.5 Hydrodynamic (Madelung) interpretation (what the terms mean)

Write (\psi=\sqrt{\rho},e^{i\theta}). Then define
[
\mathbf v=\frac{\hbar}{m}\nabla_4\theta.
]
GNLS decomposes into:

* continuity: (\partial_t\rho+\nabla_4\cdot(\rho\mathbf v)=0),
* an Euler-like relation with enthalpy (h(\rho)), confinement (V_{\rm conf}), and a dispersive “quantum pressure” term arising from (\nabla_4^2\sqrt{\rho}/\sqrt{\rho}).

We do not need to commit to a particular explicit “quantum stress tensor” form for wall forces, because we can compute generalized forces via (-\partial_{a,L}H_{\rm fluid}) (validated symbolically).

---

### 3.6 Baseline constraint choice (Paper 7 scans)

To define a “particle at rest” and the equilibrium family, we must decide which constraint we hold fixed during geometry variation. For Paper 7 baseline scans we fix (N), and keep (m_{\rm eff}) as an alternate scenario:

* **Baseline (fixed norm / particle content):**
  [
  N=\int |\psi|^2\,d^4X \quad \text{held fixed.}
  ]
* **Alternate (future extension): fixed effective rest mass target:**
  [
  m_{\rm eff}=\frac{H_{\rm tot}}{c^2} \quad \text{held fixed or matched to a scale.}
  ]

Both are valid as modeling choices, but they produce different variational problems; Paper 7 declares (N)-fixed as the baseline.

---

### 3.7 What is frozen vs what is still a choice (end of superfluid section)

**Frozen:**

* EOS (n=5) and its derived (h(\rho)), (U(\rho)), (c_s^2(\rho)).
* GNLS PDE structure and the fluid Hamiltonian (H_{\rm fluid}) in 4D.
* Continuity/current identity.

**Still to choose (must be frozen next):**

* The numerical/nondimensionalization conventions for ((\hbar,m,K)) in simulation units (so parameter scans are interpretable).

---

## 4. Independent Wave-Support Model (Separate Track)

### 4.1 Why we model wave support as an independent field

Wave support is being treated as an explicit, separate mechanism because we want to answer a clean question:

* **Does a stable particle-throat equilibrium exist without standing waves?**
* **Does a stable equilibrium exist only if standing-wave radiation stress contributes?**
* If both exist, **how do the equilibria and emergent response operators differ?**

If we buried “wave pressure” inside the superfluid model, we would lose the ability to attribute causality. So Paper 7 keeps wave support as its own field theory.

---

### 4.2 Minimal wave field: scalar surrogate in 4D

We introduce a real scalar field (A(\mathbf X,t)) intended as an “EM-like” stored-energy mode carrier (upgradeable later to Maxwell/Vector).

A clean baseline equation (no sources):
[
\partial_t^2 A - c^2\nabla_4^2 A + \mu_A^2(\mathbf X;a,L),A = 0,
]
or with a drive to excite a standing mode:
[
\partial_t^2 A - c^2\nabla_4^2 A + \mu_A^2(\mathbf X;a,L),A = S_{\rm drive}(\mathbf X,t).
]

Here:

* (c) is the model’s characteristic signal speed for the wave sector.
* (\mu_A^2(\mathbf X;a,L)) is a **confinement/mass** term that localizes the wave to the throat and creates reflective behavior without imposing hard boundaries.

---

### 4.3 Wave energy (what produces “pressure”)

Define the wave Hamiltonian:
[
H_{\rm wave}[A;a,L]
===================

\int_{\mathbb R^4}\frac12
\left[
(\partial_t A)^2
+
c^2|\nabla_4 A|^2
+
\mu_A^2(\mathbf X;a,L) A^2
\right]d^4X.
]

This is the entire “wave-stress support” mechanism: the geometry feels a generalized force derived from how the stored wave energy changes with geometry.

---

### 4.4 Geometry coupling and wave-induced breathing force

We couple waves to geometry by letting (\mu_A^2) depend on ((a,L)) (or by using a wave-wall potential shaped by the same gates as (V_{\rm conf})). Then the generalized forces are unambiguous:

[
F_a^{\rm wave}=-\frac{\partial H_{\rm wave}}{\partial a}
========================================================

-\int_{\mathbb R^4}\frac12,A^2,\frac{\partial \mu_A^2}{\partial a},d^4X,
]
[
F_L^{\rm wave}=-\frac{\partial H_{\rm wave}}{\partial L}
========================================================

-\int_{\mathbb R^4}\frac12,A^2,\frac{\partial \mu_A^2}{\partial L},d^4X.
]

This is the “optomechanics-style” identity we want: it is the cleanest way to represent radiation pressure / mode support without introducing extra stress-tensor conventions.

---

### 4.5 “Standing wave” as a controlled state, not a postulate

A standing-wave-supported throat requires a controlled mode population. Paper 7 will treat this explicitly:

* either via a drive (S_{\rm drive}) tuned near an eigenfrequency,
* or via initial data prepared in an eigenmode.

The key measurable is the stored wave energy (H_{\rm wave}) and its geometry derivatives; whether the wave is exactly “standing” is not assumed—it’s verified by the state (e.g., time-averaged energy density localization and oscillation pattern).

---

### 4.6 Separation from the fluid-only experiment

Paper 7 maintains strict separation:

* **Fluid-only track:** set (A\equiv 0), ignore (H_{\rm wave}).
* **Wave-supported track:** evolve/solve both (\psi) and (A), include (H_{\rm wave}) in the geometry closure.

Only after each track is understood independently do we allow “hybrid” runs.

---

### 4.7 Upgrade path (future work, not required for Paper 7 core)

Once the scalar surrogate is working:

* upgrade to a vector wave field (Maxwell-like),
* introduce polarization/mode families,
* couple wave dynamics more tightly to superfluid state (if desired),
* later extend to MHD/plasma variables.

Paper 7’s claim does not require this upgrade; it requires **a self-consistent wave-stress mechanism** that can be compared to fluid-only support.

---

## 5. Geometry Closure: The Wall Law (Most Important Section)

### 5.1 What the wall law must accomplish

The wall law is the closure that makes the throat a physical object rather than a fixed conduit. It must:

1. Provide an equilibrium ((a,L)) (or dynamics (a(t),L(t))) from **the same total energy** as the field equations,
2. Allow **both** opening mechanisms to be tested separately,
3. Make failure conditions explicit (no equilibrium, unstable equilibrium, nonlocal resonance domination, etc.).

The fundamental choice is whether we treat geometry as:

* reduced DOFs ((a,L)) (Paper 7 baseline), or
* a shape field (a(\cdot)) (future upgrade).

Paper 7 baseline is reduced DOFs because it is sufficient to derive (Z^{\rm eff}(\omega)) while keeping stability questions tractable.

---

### 5.2 Geometry energy (E_{\rm geom}(a,L))

We introduce a minimal geometry cost functional:
[
E_{\rm geom}(a,L)=P_{\rm vac},V(a,L)+\sigma,A(a,L)+\kappa_b,B(a,L)+\cdots
]

Where:

* (P_{\rm vac}V) represents the “vacuum work” cost of maintaining an open bulk-access region.
* (\sigma A) is a surface/hypersurface energy term (often stabilizing).
* (\kappa_b B) is an optional bending/curvature penalty (useful if you allow taper/flare as true shape DOFs).

**Non-negotiable requirement:** Paper 7 must define what it means by (V(a,L)) and (A(a,L)) in hard-mode 4D (and not silently change dimensional conventions mid-derivation).

---

### 5.3 Total energy for the two experiments

Define the fluid Hamiltonian (H_{\rm fluid}[\psi;a,L]) from Section 3 and wave Hamiltonian (H_{\rm wave}[A;a,L]) from Section 4. Then:

#### Experiment F (fluid-only supported throat)

[
H_{\rm tot}^{(F)}(\psi;a,L)=H_{\rm fluid}[\psi;a,L]+E_{\rm geom}(a,L).
]

#### Experiment W (wave-supported throat)

[
H_{\rm tot}^{(W)}(\psi,A;a,L)=H_{\rm fluid}[\psi;a,L]+H_{\rm wave}[A;a,L]+E_{\rm geom}(a,L).
]

This separation is deliberate and must remain explicit in the paper.

---

### 5.4 Equilibrium wall law (baseline closure)

At equilibrium, geometry minimizes total energy subject to a declared constraint (fixed (N) or fixed (m_{\rm eff}), etc.). The baseline equilibrium conditions are:
[
\frac{\partial H_{\rm tot}}{\partial a}=0,
\qquad
\frac{\partial H_{\rm tot}}{\partial L}=0,
]
with (\psi) (and (A) if present) taken as the stationary state for that geometry.

The practical computational interpretation is:

* for each ((a,L)), solve for the stationary (\psi) (and (A), if wave-supported),
* evaluate (H_{\rm tot}(a,L)),
* then locate minima/saddles and stability.

**Constrained stationary equation (fixed (N) baseline):**
[
-\frac{\hbar^2}{2m}\nabla_4^2\psi
+V_{\rm conf}(\mathbf X;a,L)\psi
+\frac{5K}{4}|\psi|^8\psi
=\mu\psi,
\qquad
N=\int |\psi|^2\,d^4X\ \text{fixed.}
]

**Notation note:** use a single symbol consistently in the derivation (e.g., (\psi) or (\psi_0(\mathbf X))) to avoid confusion with script variable names.

---

### 5.5 Dynamic wall law (optional but recommended for stability diagnostics)

If we want explicit stability and transient behavior, we add inertial and damping terms for the geometry DOFs:
[
M_a\ddot a+C_a\dot a = -\frac{\partial H_{\rm tot}}{\partial a},
\qquad
M_L\ddot L+C_L\dot L = -\frac{\partial H_{\rm tot}}{\partial L}.
]

This is not needed to define equilibrium, but it is a clean way to test:

* whether equilibrium is dynamically stable,
* whether perturbations grow (mode blow-up) or decay.

---

### 5.6 Force identities that remove ambiguity (key “correct-first” move)

Because geometry enters the model through potentials, the generalized forces can be written without choosing a particular stress tensor form:

For the fluid sector:
[
-\frac{\partial H_{\rm fluid}}{\partial a}
==========================================

-\int_{\mathbb R^4} |\psi|^2,\frac{\partial V_{\rm conf}}{\partial a},d^4X
\quad
(\text{when }a\text{ appears only through }V_{\rm conf}),
]
and similarly for (L).

For the wave sector:
[
-\frac{\partial H_{\rm wave}}{\partial a}
=========================================

-\int_{\mathbb R^4}\frac12,A^2,\frac{\partial \mu_A^2}{\partial a},d^4X
\quad
(\text{when }a\text{ appears only through }\mu_A^2),
]
and similarly for (L).

These identities have already been validated symbolically in the Mathematica harness, which is the strongest reason we adopted the “geometry as potential” strategy.

**Frozen geometry convention (Paper 7 baseline):**
We use the cylindrical approximation for the tube region,
[
B^3(a)\times\left[-\frac{L}{2},+\frac{L}{2}\right],
]
so:
[
V(a,L)=\left(\frac{4\pi}{3}a^3\right)L,
\qquad
A(a,L)=(4\pi a^2)L+2\left(\frac{4\pi}{3}a^3\right),
]
and
[
E_{\rm geom}(a,L)=P_{\rm vac}V(a,L)+\sigma A(a,L).
]
The corresponding derivatives are:
[
\partial_a E_{\rm geom}
=P_{\rm vac}(4\pi a^2 L)
+\sigma\left(8\pi a L+8\pi a^2\right),
]
[
\partial_L E_{\rm geom}
=P_{\rm vac}\left(\frac{4\pi}{3}a^3\right)
+\sigma(4\pi a^2).
]

---

### 5.7 Stability and “make it fail” criteria (built into the closure)

Paper 7 explicitly treats failure as a result, not an inconvenience. Typical failure modes include:

* **No equilibrium:** (H_{\rm tot}(a,L)) has no minimum in any physical parameter range.
* **Unstable equilibrium:** stationary point exists but has unstable directions (negative curvature) and grows under dynamics.
* **Collapse under through-flow:** fluid-only support leads to density depletion and loss of confinement (no sustained throat).
* **Wave support nonlocality:** equilibria exist only near resonant drive, yielding response dominated by sharp resonances and no low-(\omega) locality.
* **Parameter sensitivity:** equilibria depend strongly on numerical knobs (resolution/damping), indicating non-predictive closure.

These become explicit paper diagnostics rather than ad hoc interpretations.

---

### 5.8 What is frozen vs what is still a choice (end of wall-law section)

**Frozen:**

* Wall closure is derived from total energy via (\partial_{a,L}H_{\rm tot}), not imposed stress balance.
* Two experiments (fluid-only and wave-supported) remain distinct.

**Still to choose (must be frozen next):**

* Whether Paper 7 uses equilibrium-only or includes dynamic wall DOFs in the main text.
* The constraint definition for the variational problem in alternate scenarios (fixed (m_{\rm eff}) and how (m_{\rm eff}) is defined).

---

## 6. Verified Mathematics: What We’ve Proven Symbolically (Mathematica Harness)

### 6.1 Purpose of the verification suite

Before committing to any specific throat potential family or numerical implementation, we require that the *master equations* are internally consistent as a field theory:

* equations of motion truly follow from the stated variational principle,
* conservation laws follow from those equations (not from assumption),
* generalized geometry forces are consistent with the same Hamiltonian,
* energy/work bookkeeping remains correct when geometry changes in time.

This ensures Paper 7’s derivations are not “physics by notation.”

---

### 6.2 GNLS EOM derived from the Lagrangian (4D, (n=5))

We symbolically derived the GNLS equation by treating (\psi) and (\psi^*) as independent fields in the variational calculus and verifying that the Euler–Lagrange equations yield:

[
i\hbar,\partial_t\psi
=====================

\left[
-\frac{\hbar^2}{2m}\nabla_4^2
+V_{\rm conf}(\mathbf X;a,L)
+\frac{5K}{4}|\psi|^8
\right]\psi.
]

This fixes sign conventions and the critical factor (\frac{5K}{4}) arising from the (n=5) enthalpy term.

---

### 6.3 Continuity equation and current (exact identity)

Using the verified GNLS EOM, we verified the continuity equation:

[
\partial_t\rho+\nabla_4\cdot\mathbf J=0,\qquad
\rho=|\psi|^2,
]
with current
[
\mathbf J=\frac{\hbar}{2im}\left(\psi^*\nabla_4\psi-\psi\nabla_4\psi^*\right).
]

This is foundational for defining measurable fluxes later (ports, leakage, intake).
The harness also confirms that this continuity identity remains exact when (V_{\rm conf}) is explicitly time-dependent (e.g., (a(t)), (L(t))) so a real moving wall does not introduce source terms.

---

### 6.4 Breathing-force identity (Hellmann–Feynman form)

We verified that when geometry parameters ((a,L)) enter the fluid model only through (V_{\rm conf}), the generalized “breathing” forces are:

[
F_a^{\rm fluid}=-\frac{\partial H_{\rm fluid}}{\partial a}
==========================================================

-\int_{\mathbb R^4}\rho,\frac{\partial V_{\rm conf}}{\partial a},d^4X,
\qquad
F_L^{\rm fluid}=-\frac{\partial H_{\rm fluid}}{\partial L}
==========================================================

-\int_{\mathbb R^4}\rho,\frac{\partial V_{\rm conf}}{\partial L},d^4X.
]

This is the key mathematical reason we chose “geometry as potential”: it removes ambiguity about stress tensors at steep walls and makes the wall law computable.

We have now explicitly validated **both** generalized force identities (for (a) and for (L)) in the harness, and similarly for the wave sector via (\mu_A^2(\mathbf X;a,L)).

---

### 6.5 Wave field EOM and sign conventions (separate experiment)

For the independent wave-support field (A(\mathbf X,t)), we verified that the chosen Lagrangian produces the intended wave equation, with an explicit acknowledgment of Euler–Lagrange sign conventions for (L=T-V). Concretely, the harness checks that the derived EL equation matches the expected “standard wave operator” form up to the known overall sign convention.

This keeps the wave-support track algebraically consistent with the stated Hamiltonian.

---

### 6.6 Wave breathing-force identity

We verified the wave-analog of the breathing-force identity (when geometry dependence enters via (\mu_A^2)):

[
F_a^{\rm wave}=-\frac{\partial H_{\rm wave}}{\partial a}
========================================================

-\int_{\mathbb R^4}\frac12,A^2,\frac{\partial \mu_A^2}{\partial a},d^4X,
\qquad
F_L^{\rm wave}=-\frac{\partial H_{\rm wave}}{\partial L}
========================================================

-\int_{\mathbb R^4}\frac12,A^2,\frac{\partial \mu_A^2}{\partial L},d^4X.
]

This is what makes “wave-stress support” operational without invoking additional stress-tensor conventions.

---

### 6.6.1 Total-force additivity check (regression guard)

We additionally verify that differentiating the combined Hamiltonian reproduces the sum of sector forces:
[
-\partial_{a,L}H_{\rm tot}
=
\big(-\partial_{a,L}H_{\rm fluid}\big)+\big(-\partial_{a,L}H_{\rm wave}\big).
]
This “differentiate after summing” guard prevents bookkeeping errors when combining the fluid and wave sectors.

---

### 6.7 Thermodynamics and sound-speed consistency (n=5 closure)

We verified that the internal energy density (U(\rho)=\frac{K}{4}\rho^5) implies:

* enthalpy (h(\rho)=U'(\rho)=\frac{5K}{4}\rho^4),
* pressure (P(\rho)=\rho h(\rho)-U(\rho)=K\rho^5),
* sound-speed identity (dP/d\rho=\rho h'(\rho)) and thus (c_s^2=5K\rho^4).

This closes the loop between the GNLS nonlinear term and the hydrodynamic interpretation.

---

### 6.8 Moving wall work/energy balance

We verified an energy balance identity of the form:

[
\partial_t\mathcal H + \nabla_4\cdot \mathbf S = \rho,\partial_t V_{\rm conf},
]

where (\mathcal H) is the fluid Hamiltonian density and (\mathbf S) is the associated energy flux density (as constructed in the harness).

Interpretation:

* If (V_{\rm conf}) changes in time (because (a(t)), (L(t)), or other shape DOFs change), the energy change is exactly accounted for by the work term (\rho,\partial_t V_{\rm conf}).
* This is the bookkeeping guarantee we need once we allow breathing dynamics.

The harness explicitly substitutes (a\to a(t)), (L\to L(t)) and verifies the chain rule:
[
\partial_t V_{\rm conf}=(\partial_a V_{\rm conf})\dot a+(\partial_L V_{\rm conf})\dot L.
]

---

### 6.8.1 Wave moving-wall work balance (newly verified)

For the wave sector, the harness now verifies the moving-wall work identity:
[
\partial_t\mathcal H_A+\nabla_4\cdot\mathbf S_A
=
\frac12 A^2\,\partial_t\mu_A^2.
]
with canonical flux (\mathbf S_A=-c^2 A_t \nabla A) in the chosen normalization.

---

### 6.9 Madelung/hydrodynamic consistency check

We verified that the GNLS nonlinear potential term corresponds exactly to the enthalpy (h(\rho)) under the Madelung interpretation:
[
\frac{(5K/4)\rho^4\psi}{\psi}= \frac{5K}{4}\rho^4 = h(\rho).
]
This is the “microphysics consistency” check that the GNLS closure matches the intended (n=5) barotropic fluid.

---

### 6.10 Known harness issues and hardening notes

The current script set is strong enough for Paper 7’s “core derivation validity,” but we should still:

* fix mapping typos (so future checks don’t silently rely on a lucky cancellation),
* harden mixed-derivative substitution rules to cover both derivative orderings ((\partial_t\partial_x\psi) and (\partial_x\partial_t\psi)),
* add the (L)-localization analog of the wall localization test: verify (\partial_L V_{\rm conf}) is localized at endcaps,
* replace the generic tanh radial test with a frozen candidate of the actual throat gate (tube × endcaps) once Family 1/2 is frozen.

These are maintenance items, not conceptual gaps.

---

## 7. Observables: How 3D Physics Emerges From the 4D Solve

### 7.1 The key principle: “brane observables” are measurements, not separate equations

Since we solve one unified 4D model, a 3D observer’s “physics” must be defined as a measurement procedure. The brane is associated with the neighborhood of (w=0), but there is no hard boundary at (w=0).

Therefore, every “brane-side variable” must be defined by:

* selecting a measurement region near (w=0),
* projecting/weighting the 4D fields onto that region.

This is the correct hard-mode replacement for “match PDEs at an interface.”

---

### 7.2 Brane-projected density

Define a brane-weighted density:
[
\rho_{\rm brane}(x,y,z,t)=\int_{-\infty}^{\infty} W(w),|\psi(x,y,z,w,t)|^2,dw,
]
where (W(w)) is a fixed weight function concentrated near (w=0). Two canonical choices:

* (W(w)=\delta(w)) (idealized “exact brane slice,” used only formally),
* (W(w)=|\chi_0(w)|^2), where (\chi_0) is the ground state of the far-field brane trap (V_{\rm brane}(w)) (preferred for a physically meaningful projection).

This defines what “density on the brane” means without breaking the 4D model.

---

### 7.3 Brane effort variable: enthalpy perturbation

For small perturbations about a brane-background (\rho_{\rm brane,0}), define the canonical effort variable:
[
u \equiv \delta h(\rho_{\rm brane})
\approx \left.\frac{dh}{d\rho}\right|*{\rho*{\rm brane,0}}\delta\rho_{\rm brane}
=5K\rho_{\rm brane,0}^3,\delta\rho_{\rm brane}.
]

This choice is aligned with the (n=5) thermodynamic structure and is convenient for defining an impedance-like operator.

---

### 7.4 Brane flux variable: choose what is being “measured”

Flux is not unique in hard-mode 4D, and Paper 7 must explicitly state which flux is the conjugate of (u). Two important flux notions exist:

1. **Brane-mouth normal mass flux** (what a brane observer would interpret as “in/out of the mouth”):
   choose a 2D disk (\mathcal D) in ((x,y,z)) near the mouth and define
   [
   j_{\rm mouth}(t)=\int_{\mathcal D}\rho_{\rm brane},v_n,dS,
   ]
   where (v_n) is the normal component of the effective brane velocity.

2. **Bulk-leakage flux in the (w) direction** (how much 3D-projected density is leaving into the bulk):
   define the leakage current density component
   [
   J_w=\frac{\hbar}{2im}\left(\psi^*\partial_w\psi-\psi\,\partial_w\psi^*\right),
   ]
   and a “leakage monitor” by integrating across a (w)-surface,
   [
   j_w(t)=\int_{\Omega_{xyz}} J_w\,d^3x,
   ]
   possibly weighted near the mouth region. A concrete baseline is to integrate over a 3-ball (r<R_{\rm measure}) at (w=\pm W_{\rm cutoff}).

Which one is “the” flux depends on the experiment:

* in through-flow supported throats, both may matter (intake vs leakage),
* in wave-supported throats, mouth flux may be small even when wave energy is large.

---

### 7.5 Ports: turning fields into a finite response operator

To define an effective impedance/operator, we choose a small basis of “port patterns” localized on the brane near the mouth region. Let (\Gamma) be the chosen measurement region (typically a 2D disk or a thin 3D slab near (w=0)), and choose basis functions (P_i(\mathbf s)) on that region.

Define port amplitudes:
[
u_i(t)=\int_{\Gamma} P_i(\mathbf s),u(\mathbf s,t),d\mu,
\qquad
j_i(t)=\int_{\Gamma} P_i(\mathbf s),j(\mathbf s,t),d\mu,
]
where (d\mu) is the appropriate measure for (\Gamma).

Then define the emergent response operator by drive/measure:
[
j_i(\omega)=\sum_j Z^{\rm eff}_{ij}(\omega),u_j(\omega).
]

This is the hard-mode replacement for “inner DtN.”

---

### 7.6 Drive protocols (how we excite the system to measure (Z^{\rm eff}))

Because we are not imposing boundary conditions at a mouth, we define a drive as a localized source term or localized modulation in the near-brane region. Canonical examples:

* “enthalpy drive”: add a small localized forcing term to the GNLS chemical potential near the mouth region,
* “flux drive”: impose a phase bias or current-like source localized near the mouth region,
* “wave drive”: excite (A) with (S_{\rm drive}) to populate a standing mode (wave-supported track).

Each drive is then mapped into measured ((u_i,j_i)) to extract (Z^{\rm eff}).

---

### 7.7 What is frozen vs what is still a choice (end of observables section)

**Frozen (conceptual):**

* Brane observables are projections/measurements derived from the 4D solve.
* Effort variable uses the (n=5) enthalpy perturbation structure.
* Leakage flux is a baseline observable, defined by (J_w) and a fixed measurement surface (e.g., 3-ball at (w=\pm W_{\rm cutoff})).

**Still to choose (must be frozen next):**

* The explicit weight (W(w)) (delta-slice vs ground-state weighting).
* The measurement region (\Gamma) geometry and measure (d\mu).
* The precise flux definition used for the impedance operator (mouth flux vs leakage flux, or both in a multi-output operator).
* The port basis (P_i) and drive protocol conventions (including sign conventions).

---

## 8. Response Operator Extraction in Hard-Mode 4D

### 8.1 What replaces DtN in hard-mode 4D

In earlier “boundary-matching” approaches, an inner Dirichlet-to-Neumann (DtN) map is a boundary operator: prescribe a boundary value, read off a normal derivative, and you get an impedance-like object.

In hard-mode 4D, there is no privileged boundary surface. The replacement is an **operationally defined response operator**:

1. Choose a near-brane measurement region (\Gamma) (disk/slab near (w=0)).
2. Choose port basis functions (P_i) on (\Gamma).
3. Choose conjugate variables ((u,j)) (effort/flux) derived from 4D fields.
4. Drive the system with a small, localized forcing protocol.
5. Measure ((u_i(\omega),j_i(\omega))) and define
   [
   j_i(\omega)=\sum_j Z^{\rm eff}_{ij}(\omega),u_j(\omega).
   ]

So (Z^{\rm eff}) is a **measured effective operator**, not a theoretical boundary map.

---

### 8.2 Canonical effort/flux pair (and why this choice is stable)

We take the effort variable as enthalpy perturbation (consistent with the (n=5) EOS):
[
u(\mathbf s,t)=\delta h(\rho_{\rm brane})
\approx 5K\rho_{\rm brane,0}^3,\delta\rho_{\rm brane}.
]

Flux is the conjugate quantity chosen by the measurement definition. The two main candidates are:

* **Mouth flux** (brane-side normal mass flux through a chosen 2D disk (\mathcal D)):
  [
  j_{\rm mouth}(t)=\int_{\mathcal D}\rho_{\rm brane},v_n,dS.
  ]

* **Bulk leakage flux** (4D current through a (w)-oriented monitor surface):
  [
  j_w(t)=\int_{\Omega_{xyz}} J_w\,d^3x,
  ]
  optionally localized near the mouth region.

Paper 7 must declare which one defines (j) for the impedance operator. If both are important, treat this as a **multi-output operator** (same input ports, two measured outputs).

---

### 8.3 Port basis and projection

Let (\Gamma) be the chosen measurement region near the mouth (disk or thin slab). Choose a finite set of basis functions ({P_i(\mathbf s)}_{i=1}^N) with a clear physical meaning (monopole-like, dipole-like, quadrupole-like patterns, etc.).

Project to port amplitudes:
[
u_i(t)=\int_{\Gamma}P_i(\mathbf s),u(\mathbf s,t),d\mu,
\qquad
j_i(t)=\int_{\Gamma}P_i(\mathbf s),j(\mathbf s,t),d\mu.
]

The measure (d\mu) is determined by (\Gamma)’s dimensionality (2D disk vs 3D slab) and must be stated once and never changed silently.

---

### 8.4 Drive protocols (how we excite the operator)

We need drives that are:

* localized near the mouth region,
* small enough for linear response,
* defined without imposing boundary conditions.

Three standard drives:

**(D1) Enthalpy/chemical potential drive (fluid sector)**
Add a small forcing term to the GNLS potential near (\Gamma):
[
V_{\rm conf}\to V_{\rm conf}+\epsilon,f(\mathbf X)\cos(\omega t),
]
or equivalently add a small source term in the chemical potential sector (depending on implementation).

**(D2) Flux/phase-bias drive (through-flow style)**
Impose a controlled phase gradient or current-like source localized near (\Gamma), designed to create a net intake/through-flow.

**(D3) Wave drive (wave-supported track)**
Drive the wave field:
[
S_{\rm drive}(\mathbf X,t)=\epsilon,s(\mathbf X)\cos(\omega t),
]
and measure how the induced standing-wave energy modifies the effective operator.

Paper 7 should standardize one drive as the primary definition of (Z^{\rm eff}) and treat others as cross-checks.

---

### 8.5 Linear response extraction procedure

For each chosen frequency (\omega):

1. Find a baseline stationary state (\psi_0(\mathbf X)) (and (A_0) if wave-supported) for the chosen ((a,L)).
2. Apply a small drive (\epsilon) and evolve to a periodic steady state.
3. Compute port time series (u_i(t), j_i(t)) over multiple cycles.
4. Fourier extract the fundamental component at (\omega): (u_i(\omega), j_i(\omega)).
5. Assemble the matrix (Z^{\rm eff}(\omega)) by solving the linear system:
   [
   \mathbf j(\omega)=Z^{\rm eff}(\omega),\mathbf u(\omega).
   ]

If using multiple independent drives (e.g., exciting each port separately), the matrix is identified cleanly.

---

### 8.6 Passivity and consistency diagnostics (operator sanity)

To keep (Z^{\rm eff}) physically interpretable, we impose diagnostics:

* **Energy consistency:** the power input from the drive equals the measured energy accumulation + flux/leakage.
* **Passivity (when damping/leakage is included):** response should not produce net energy without input under the chosen conventions.
* **Robustness:** (Z^{\rm eff}) should not depend strongly on numerical knobs (resolution, small regularization terms) in a regime claimed to be predictive.

---

### 8.7 Low-frequency locality gate (EFT/PDE reduction test)

A key question is whether the measured operator admits a low-frequency expansion away from resonances:
[
Z^{\rm eff}*{ij}(\omega)\sim A*{ij}\omega^2+B_{ij}\omega^4+\cdots
\quad (\text{or the appropriate analytic structure for the chosen }u,j).
]

If (Z^{\rm eff}) is dominated by sharp resonances in the band of interest, then any local-in-time effective description will fail, and Paper 7 must say so explicitly.

---

### 8.8 What is frozen vs what is still a choice (end of response section)

**Frozen (conceptual):**

* (Z^{\rm eff}) is defined operationally by drive/measure in a near-brane region.
* Ports, effort/flux variables, and drive protocols define the operator.

**Still to choose (must be frozen next):**

* The measurement region (\Gamma) (disk vs slab) and its measure.
* Whether (j) is mouth flux, leakage flux, or a multi-output vector.
* The canonical port basis ({P_i}) and the canonical drive protocol used to define (Z^{\rm eff}).
* The damping/leakage conventions used for passivity testing.

---

## 9. Scenario Library (So We Don’t Mix Physics)

### 9.1 Why scenarios exist

Hard-mode 4D has many moving parts (geometry potentials, fluid PDE, wave PDE, closure law, observables). Without explicit scenario modes, it becomes too easy to accidentally mix assumptions and interpret results incorrectly.

Scenarios are the “physics states” under which we run the model. Each scenario declares:

* which fields are active,
* which drives are allowed,
* which geometry closure is used,
* which measurements define (Z^{\rm eff}).

---

### 9.2 Scenario S1: Particle baseline (fluid-only support)

**Purpose:** establish whether a particle-throat equilibrium exists supported by fluid pressure/momentum flux alone.

* Active fields: (\psi) only ((A\equiv 0)).
* Geometry closure: (H_{\rm tot}=H_{\rm fluid}+E_{\rm geom}).
* Drives: through-flow / phase bias permitted if this is the mechanism under test.
* Observables: brane-projected (u=\delta h), flux definition chosen (mouth, leakage, or both).

Key outputs:

* equilibrium ((a,L)) existence and stability,
* (Z^{\rm eff}(\omega)) for that equilibrium,
* failure conditions if equilibrium cannot be sustained.

---

### 9.3 Scenario S2: Wave-supported particle baseline

**Purpose:** test the “standing wave holds it open” mechanism cleanly.

* Active fields: (\psi) and independent (A).
* Geometry closure: (H_{\rm tot}=H_{\rm fluid}+H_{\rm wave}+E_{\rm geom}).
* Drives: wave drive (S_{\rm drive}) is primary; fluid drive optional but should be off in the baseline run.
* Observables: same (u) definition; measure how wave energy modifies (Z^{\rm eff}) and equilibrium.

Key outputs:

* whether equilibria appear that do not exist in fluid-only runs,
* how equilibrium geometry depends on stored wave energy,
* resonance dominance vs locality.

---

### 9.4 Scenario S3: Hybrid regime (both supports present)

**Purpose:** only after S1 and S2 are understood, explore combined effects.

* Active fields: (\psi) and (A).
* Geometry closure: full (H_{\rm fluid}+H_{\rm wave}+E_{\rm geom}).
* Drives: both fluid and wave drives allowed, but applied in controlled sequences.
* Observables: operator extraction with careful bookkeeping of energy channels (fluid vs wave vs geometry).

This scenario is likely rich but easy to misinterpret; it must not be the first thing we run.

---

### 9.5 Scenario S4: EM/cavity selector regression limit (unit test)

**Purpose:** ensure continuity with earlier selector behavior and prevent “branch drift.”

Operationally:

* choose potential parameters that strongly suppress 4D transverse structure (e.g., strong (w)-trap everywhere or a limit like (\gamma\to\infty) in Family 2),
* choose geometry/termination settings that act like a reflective cavity.

The goal is not to claim the universe is always in this limit; it’s to verify that the hard-mode framework contains the known-answer behavior as a controlled special case.

---

### 9.6 Scenario S5: Failure / breakdown regimes (paper-worthy)

Paper 7 treats failure outcomes as results. Scenario S5 is a structured way to generate those:

* No equilibrium for any reasonable (E_{\rm geom}) parameters.
* Equilibrium exists but is dynamically unstable.
* Fluid-only support collapses density and breaks confinement.
* Wave-supported equilibria exist only at resonance, destroying locality of (Z^{\rm eff}).
* Observables become ill-defined: brane projection weight (W(w)) fails to isolate “brane physics,” or leakage dominates.

This scenario’s output is a “phase diagram” of where the model is predictive vs not.

---

### 9.7 What is frozen vs what is still a choice (end of scenario section)

**Frozen:**

* Scenarios exist to prevent accidental mixing of mechanisms and to keep interpretation clean.
* S1 and S2 must be run and understood before S3.

**Still to choose (must be frozen next):**

* The exact definitions of drives for each scenario (especially phase-bias vs potential-modulation drives).
* The canonical measurement region/ports used consistently across scenarios.
* The “regression limit” parameter settings that define S4.

---

## 10. Implementation Plan (What We Build Next)

### 10.1 Freeze the explicit (V_{\rm conf}) first (this is the keystone)

Before any numerics, we must choose one explicit confinement family as the baseline and write it in final form:
[
V_{\rm conf}(\mathbf X;a,L)=V_{\rm brane}(w;\text{trap})+V_{\rm throat}(\mathbf X;a,L;\text{shape}).
]

**Immediate output of this step:**

* the explicit formula for (V_{\rm conf}),
* explicit formulas for (\partial_a V_{\rm conf}) and (\partial_L V_{\rm conf}),
* a qualitative “support map” showing where these derivatives are nonzero (wall region vs endcaps).

This is the point at which “hard-mode 4D geometry” becomes concrete.

**Status:** completed in `mathematica/inner_throat/throat_spec_derivation.wl` (Revised v2): explicit (V_{\rm conf}), fixed-(N) constrained stationary equation, and force-balance kernels.

---

### 10.2 Freeze (E_{\rm geom}(a,L)) and the hard-mode measures (V(a,L)), (A(a,L))

Write explicit expressions for:

* the 4D “volume-like” measure (V(a,L)),
* the 4D “area-like” measure (A(a,L)),

consistent with whatever geometric convention we adopt (tube vs funnel vs rounded corridor). Then freeze:
[
E_{\rm geom}(a,L)=P_{\rm vac}V(a,L)+\sigma A(a,L)+\kappa_b B(a,L)+\cdots
]

**Note:** this is where many hidden inconsistencies creep in (dimension changes, missing factors). We freeze it once and never silently redefine it later.

**Status:** frozen in the baseline cylindrical approximation (see Section 5.6) and implemented in `mathematica/inner_throat/throat_spec_derivation.wl`.

---

### 10.3 Extend the Mathematica harness to cover what we’re about to implement

The current harness already validates the core field theory. Before moving into numerics, add these as mandatory checks:

* **Mirror force check for (L)** (in addition to (a)), with a tolerance-based central-leakage criterion under (\varepsilon_w \ll d_w \ll L).
* **Explicit chain-rule moving wall check** with (a\to a(t)), (L\to L(t)), verifying:
  [
  \partial_t V_{\rm conf} = \partial_a V_{\rm conf},\dot a + \partial_L V_{\rm conf},\dot L.
  ]
* **Substitute a concrete (V_{\rm conf})** from the chosen family and verify (\partial_a V_{\rm conf}), (\partial_L V_{\rm conf}) have the expected localization.

This ensures we’re not implementing a potential whose derivatives create nonsense breathing forces.

---

### 10.4 Minimal numerical milestones (in hard-mode 4D order)

These are the smallest steps that still produce paper-meaningful results:

**M1 — Static ground state at fixed geometry**
For fixed ((a,L)), find a stationary (\psi_0(\mathbf X)) (and (A_0) if wave-supported). Practically this can be:

* imaginary-time relaxation (for GNLS),
* or constrained minimization of (H_{\rm fluid}) under fixed (N).

**M2 — Geometry equilibrium scan**
Compute (H_{\rm tot}(a,L)) on a coarse grid and locate minima/saddles for:

* fluid-only (H_{\rm fluid}+E_{\rm geom}),
* wave-supported (H_{\rm fluid}+H_{\rm wave}+E_{\rm geom}) (for a defined wave population).

**M3 — Linear response and operator extraction**
Choose a measurement region (\Gamma), ports (P_i), and a drive protocol. Extract:
[
j(\omega)=Z^{\rm eff}(\omega),u(\omega)
]
around one equilibrium.

**M4 — Stability diagnosis**
Perturb either:

* geometry DOFs ((a,L)),
* or small field perturbations,
  and measure growth/decay. If geometry dynamics are used, check whether the equilibrium is dynamically stable.

This milestone set is “enough for Paper 7”: equilibrium existence + effective operator + stability/failure story.

---

### 10.5 Parameter scan plan (keep it disciplined)

The scan axes should be chosen to diagnose physics, not overwhelm us. Suggested minimal scan dimensions:

* geometry energy parameters: (P_{\rm vac},\sigma) (and (\kappa_b) if needed),
* throat openness control: Family 1 (\Omega_{\rm in}/\Omega_{\rm out}) or Family 2 (\gamma),
* equilibrium geometry: (a,L),
* drive frequency (\omega) for (Z^{\rm eff}(\omega)).

Everything else (grid resolution, damping) is treated as a numerical convergence parameter and tested separately.

---

### 10.6 “Definition freeze” checklist (must be done before claiming results)

Before any plots or claims appear in Paper 7, we must have frozen:

1. (V_{\rm conf}(\mathbf X;a,L)) and its derivative support properties,
2. (E_{\rm geom}(a,L)) and the meaning of (V,A) in hard-mode 4D,
3. constraint type (fixed (N) vs fixed (m_{\rm eff})),
4. brane observable map (choice of weight (W(w)) and measurement region (\Gamma)),
5. port basis (P_i) and drive protocol that defines (Z^{\rm eff}).

As of `mathematica/inner_throat/throat_spec_derivation.wl` (Revised v2), items 1-3 are frozen for the baseline track. If any of the remaining items are ambiguous, the operator (Z^{\rm eff}) is not a well-defined object.

---

## 11. Connections to Earlier Papers (Only What Still Matters)

### 11.1 What we keep (hard constraints that anchor the universe)

Even though we’re discarding most “old scaffolding,” we retain:

* **EOS constraint:** (P(\rho)=K\rho^5) and its derived (h(\rho)), (U(\rho)), (c_s^2(\rho)).
  This keeps Paper 7 in the same microphysics regime as the hybrid model foundation.

* **Branch regression limit:** the hard-mode model must admit a controlled parameter regime that reproduces the earlier EM/cavity selector behavior (e.g., strong (w)-trap or (\gamma\to\infty) limit).
  This is a regression test, not a universal assumption.

* **Energy bookkeeping culture:** mass-equivalence (m_{\rm eff}=E_{\rm tot}/c^2) remains a defined diagnostic, even if not used as the primary constraint initially.

---

### 11.2 What we explicitly discard (to avoid contamination)

We are no longer relying on:

* a fixed inner cylinder DtN map as the “true interior,”
* a stitched interface where 4D “turns into 3D” at a boundary,
* impedance curves inserted by hand.

Those can remain as historical unit tests, but they are not the model.

---

### 11.3 How old DtN work is still useful (as regression scaffolding)

The existing inner/outer DtN machinery remains valuable in three ways:

1. **Regression/limit tests:** verify that in a parameter regime where the 4D model collapses to an effectively 3D waveguide/cavity, the response behaves similarly to earlier analytic DtN forms.

2. **Operator extraction methodology:** the philosophy of “define ports, drive, measure, build a matrix operator” carries over directly—only now it is extracted from a 4D PDE rather than a boundary solve.

3. **Failure comparisons:** if the new hard-mode model produces response features that old models cannot reproduce (or vice versa), that difference becomes a diagnostic of what physics was missing.

---

## 12. Future Extensions (After Paper 7 Is Solid)

### 12.1 Reintroducing 2PN/EFT comparisons (later, not now)

Once Paper 7 produces a stable particle-throat model and an emergent (Z^{\rm eff}(\omega)), we can compare:

* whether low-frequency locality holds (or fails) in a way compatible with EFT truncations,
* how far-field effective behavior deviates from PN expectations,
* what fundamental changes would be required to reduce mismatch.

This is deliberately postponed until the 4D model is working and defensible.

---

### 12.2 4D plasma/MHD reconnection path

If the long-term goal includes reconnection-like behavior, Paper 7’s architecture is the correct foundation because it already provides:

* a unified 4D geometry,
* a variationally consistent fluid sector,
* a clean place to add additional fields.

Future upgrades could include:

* vector gauge fields (Maxwell-like) instead of scalar (A),
* resistive/non-ideal terms with controlled gating (localized dissipation),
* helicity-like invariants and diagnostics adapted to the toy universe,
* multi-fluid or charged-fluid generalizations.

---

### 12.3 Multi-throat networks and junction manifolds

If composite objects require it, the hard-mode potential geometry approach naturally extends to:

* multiple nearby bulk corridors,
* branched throat manifolds (junctions),
* mode-coupling networks that could generate effective “fractional coupling weights” without altering global conservation structure.

This is explicitly “future work” until the single-throat case is nailed.

---

### 12.4 Observational signatures and interpretation (what to look for)

Once the model runs, we can ask what a brane observer would see:

* intermittent bursts in brane observables due to energy exchange with bulk modes,
* resonance-driven amplification windows,
* leakage signatures (energy/current leaving into (w)),
* regimes where the system is effectively 3D vs genuinely 4D.

These become the bridge from “hard-mode math” to “testable phenomenology” later.

---

## Appendix A. Canonical Notation and Conventions

### A.1 Coordinates, operators, and measures

We work in **4D space + time**:
[
\mathbf X=(x,y,z,w),\qquad t,\qquad \nabla_4=(\partial_x,\partial_y,\partial_z,\partial_w),\qquad \nabla_4^2=\partial_x^2+\partial_y^2+\partial_z^2+\partial_w^2.
]
The 4D volume measure is (d^4X=dx,dy,dz,dw).

Useful radii (used as needed by the geometry family):
[
R_3=\sqrt{x^2+y^2+z^2},\qquad r_\perp=\sqrt{x^2+y^2},\qquad R_\gamma=\sqrt{r_\perp^2+(w/\gamma)^2}.
]

### A.2 Fields and derived quantities

**Superfluid sector (always present).** Complex field (\psi(\mathbf X,t)) with density
[
\rho(\mathbf X,t)=|\psi|^2=\psi^*\psi.
]
The canonical 4D current is
[
\mathbf J=\frac{\hbar}{2im}\left(\psi^*\nabla_4\psi-\psi\nabla_4\psi^*\right),
]
and the continuity equation is
[
\partial_t\rho+\nabla_4\cdot\mathbf J=0.
]

**Wave-support sector (optional, separate experiment).** Real scalar field (A(\mathbf X,t)) used as a standing-wave energy carrier (surrogate for “EM-like” stored energy). Its confinement is encoded by a geometry-dependent (\mu_A^2(\mathbf X;a,L)).

### A.3 EOS and thermodynamic closures (fixed)

We impose the stiff polytrope
[
P(\rho)=K\rho^5,
]
which fixes
[
c_s^2(\rho)=\frac{dP}{d\rho}=5K\rho^4,\qquad
h(\rho)=\int^\rho\frac{dP}{\rho'}=\frac{5K}{4}\rho^4,\qquad
U(\rho)=\int^\rho h(\rho'),d\rho'=\frac{K}{4}\rho^5.
]

### A.4 Geometry as potential (no stitched interfaces)

Geometry enters through a smooth confining potential
[
V_{\rm conf}(\mathbf X;a,L)=V_{\rm brane}(w;\mathbf X)+V_{\rm throat}(\mathbf X;a,L),
]
with smooth gating to avoid singular wall forces.

We use a smooth step
[
S(u)=\tfrac12\left(1+\tanh u\right),
]
and a typical “inside-throat” gate (G(\mathbf X;a,L)\in[0,1]) such that (G\approx 1) inside the corridor and (G\approx 0) outside. The “mouth” is not a boundary surface; it is the transition region where the gates change.

### A.5 Governing equations (frozen forms)

**GNLS (4D, (n=5)).**
[
i\hbar,\partial_t\psi=
\left[-\frac{\hbar^2}{2m}\nabla_4^2+V_{\rm conf}(\mathbf X;a,L)+\frac{5K}{4}|\psi|^8\right]\psi.
]

**Wave surrogate (optional).**
[
\partial_t^2A-c^2\nabla_4^2A+\mu_A^2(\mathbf X;a,L),A=S_{\rm drive}.
]

### A.6 Energies and force conventions

Fluid Hamiltonian:
[
H_{\rm fluid}=\int_{\mathbb R^4}\left[
\frac{\hbar^2}{2m}|\nabla_4\psi|^2+V_{\rm conf}|\psi|^2+\frac{K}{4}|\psi|^{10}
\right]d^4X.
]
Wave Hamiltonian:
[
H_{\rm wave}=\int_{\mathbb R^4}\frac12\left[
(\partial_tA)^2+c^2|\nabla_4A|^2+\mu_A^2A^2
\right]d^4X.
]
Geometry energy (closure choice):
[
E_{\rm geom}(a,L)=P_{\rm vac}V(a,L)+\sigma A(a,L)+\kappa_b B(a,L)+\cdots.
]
Total energies:

* Fluid-only experiment: (H_{\rm tot}^{(F)}=H_{\rm fluid}+E_{\rm geom}).
* Wave-supported experiment: (H_{\rm tot}^{(W)}=H_{\rm fluid}+H_{\rm wave}+E_{\rm geom}).

**Generalized breathing forces** are defined by
[
F_a=-\partial_aH_{\rm tot},\qquad F_L=-\partial_LH_{\rm tot}.
]

### A.7 Brane observable map, ports, and operator sign conventions

Brane observables are defined by a fixed projection weight (W(w)) concentrated near (w=0):
[
\rho_{\rm brane}(x,y,z,t)=\int_{-\infty}^{\infty}W(w),|\psi(x,y,z,w,t)|^2,dw.
]
Effort variable (enthalpy perturbation):
[
u=\delta h(\rho_{\rm brane})\approx 5K\rho_{\rm brane,0}^3,\delta\rho_{\rm brane}.
]
Flux variable (j) is defined by the declared measurement convention (e.g., normal flux through a brane-region (\Gamma), and/or leakage monitor).

With port basis (P_i) on (\Gamma),
[
u_i(t)=\int_\Gamma P_i,u,d\mu,\qquad j_i(t)=\int_\Gamma P_i,j,d\mu.
]
In frequency space we define the effective response operator by
[
\mathbf j(\omega)=Z^{\rm eff}(\omega),\mathbf u(\omega),
]
with the convention that “positive real power input” corresponds to the chosen (u/j) sign pairing (checked operationally via work/energy balance).

---

## Appendix B. Symbolic Verification Suite: What Was Proven and Why It Matters

### B.1 Purpose

Paper 7’s “hard-mode 4D” architecture depends on the claim that:

1. the PDEs come from a single consistent variational structure,
2. conservation laws and energy bookkeeping follow from that same structure, and
3. geometry forces are computed unambiguously from energy derivatives (not from ad hoc stress formulas).

A Mathematica verification harness was used to validate these points symbolically.

### B.2 Verified results (GNLS sector)

The suite verifies that the stated GNLS Lagrangian produces the exact 4D equation of motion
[
i\hbar,\partial_t\psi=
\left[-\frac{\hbar^2}{2m}\nabla_4^2+V_{\rm conf}(\mathbf X;a,L)+\frac{5K}{4}|\psi|^8\right]\psi,
]
including the crucial (n=5) factor (\frac{5K}{4}).

It also verifies the exact continuity identity
[
\partial_t\rho+\nabla_4\cdot\mathbf J=0,
\qquad
\mathbf J=\frac{\hbar}{2im}\left(\psi^*\nabla_4\psi-\psi\nabla_4\psi^*\right).
]

Thermodynamic consistency checks confirm:
[
U(\rho)=\frac{K}{4}\rho^5,\quad h(\rho)=U'(\rho)=\frac{5K}{4}\rho^4,\quad
P(\rho)=\rho h(\rho)-U(\rho)=K\rho^5,\quad
\frac{dP}{d\rho}=\rho h'(\rho)=5K\rho^4.
]

Finally, the Madelung/hydrodynamic correspondence is verified at the level required here: the GNLS nonlinear potential term corresponds exactly to the enthalpy (h(\rho)).

### B.3 Verified results (geometry coupling and wall-work bookkeeping)

A central result is the **Hellmann–Feynman-style force identity**: when geometry parameters enter through (V_{\rm conf}), the generalized forces reduce to volume integrals involving (\partial_{a,L}V_{\rm conf}). Operationally this is the mathematical license to compute wall/breathing forces without introducing an ambiguous “stress tensor at the wall.”

The suite also verifies an energy/work balance identity for time-dependent confinement:
[
\partial_t\mathcal H+\nabla_4\cdot\mathbf S=\rho,\partial_tV_{\rm conf},
]
showing that changes in geometry (via time-dependent confinement) contribute exactly the expected work term. This is the consistency gate for any dynamic wall model (a(t),L(t)).

### B.4 Verified results (wave-support sector)

For the wave surrogate field, the suite verifies that the chosen (L=T-V) convention produces the expected wave equation structure (up to the standard overall Euler–Lagrange sign convention, handled explicitly in the check), and it verifies the corresponding geometry-force identity:
[
-\partial_aH_{\rm wave}=-\int\frac12A^2,\partial_a\mu_A^2,d^4X
]
(and similarly for (L)).

### B.5 What these checks do (and do not) guarantee

These symbolic checks guarantee:

* the master equations, conservation laws, and energy/force bookkeeping are mutually consistent;
* geometry forces are derivable from (-\partial_{a,L}H) without ambiguity.

They do not determine the *best* choice of (V_{\rm conf}), (E_{\rm geom}), projection (W(w)), ports, or drive protocols; those must be frozen by definition (and then tested numerically).

---

## Appendix C. “Make It Fail” Checklist and Diagnostics

This appendix defines the failure outcomes we treat as **paper-worthy results**, plus the diagnostics used to classify them.

### C.1 No equilibrium throat exists

**Failure statement:** For any physically reasonable geometry-energy parameters ((P_{\rm vac},\sigma,\kappa_b,\dots)), the total energy admits no stable equilibrium ((a,L)).

**Diagnostics:**

* compute (H_{\rm tot}(a,L)) on a coarse grid (fluid-only and wave-supported separately);
* verify no local minima exist in the interior of the scanned domain;
* confirm the absence of minima is not a numerical artifact (resolution/damping independence).

### C.2 Equilibrium exists but is dynamically unstable

**Failure statement:** Stationary points exist but are not stable under perturbations.

**Diagnostics:**

* compute the Hessian of the reduced energy landscape:
  [
  \mathcal H_{ij}=\partial_{q_i}\partial_{q_j}H_{\rm tot},\quad (q_1,q_2)=(a,L),
  ]
  and check eigenvalues (negative direction ⇒ instability);
* if using geometry dynamics, perturb (a,L) slightly and measure growth/decay.

### C.3 Fluid-only support collapses under through-flow / intake

**Failure statement:** In the fluid-only experiment, the same mechanisms that generate intake/through-flow drive density depletion or deconfinement, destroying the throat.

**Diagnostics:**

* monitor minimum/mean density in the throat region and brane region;
* monitor leakage into (w) via (J_w) and the evolution of (\rho_{\rm brane});
* verify that collapse persists under numerical refinement.

### C.4 Wave-supported equilibria exist only at resonance (nonlocal response)

**Failure statement:** Equilibria require finely tuned resonant wave excitation; the effective operator becomes resonance-dominated and fails any low-frequency locality gate.

**Diagnostics:**

* measure (Z^{\rm eff}(\omega)) over a band and identify whether response is smooth/expandable vs spiky/resonant;
* attempt a low-frequency fit (away from poles) and check whether a stable truncated expansion exists;
* quantify sensitivity of equilibrium ((a,L)) to small detuning of the drive frequency.

### C.5 Projection-to-brane observables becomes ill-defined

**Failure statement:** The chosen brane observable map does not isolate “3D physics” robustly; outcomes depend strongly on the projection weight (W(w)) or measurement-region definitions.

**Diagnostics:**

* repeat key measurements for at least two reasonable weights (W(w)) (e.g., ground-state weighting and a slightly broadened weighting);
* verify qualitative stability of (Z^{\rm eff}) and equilibrium classification under small changes in (\Gamma), port radius, and weighting width.

### C.6 Numerical non-identifiability (results depend on knobs)

**Failure statement:** Claimed physics changes materially with numerical resolution, regularization strength, artificial damping, or domain size.

**Diagnostics:**

* convergence tests in spatial resolution and timestep;
* check invariants (norm (N), energy drift in conservative runs, or consistent dissipation in damped runs);
* confirm operator extraction (Z^{\rm eff}) is stable under refinement.

### C.7 Global bookkeeping inconsistency

**Failure statement:** The energy/work accounting does not close, implying the model implementation is not faithful to the derived equations.

**Diagnostics:**

* evaluate the work balance residual (derived in Appendix B) in driven and/or breathing runs;
* reject results if the residual does not converge to zero (or to the expected controlled dissipation) with refinement.
