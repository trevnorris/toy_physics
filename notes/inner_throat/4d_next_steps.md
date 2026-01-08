## Table of Contents (additions)

- Force ledger for throat DOFs (non-negotiable bookkeeping)
- 4+1D Maxwell equations with localization (explicit PDE form)
- Neutrality / source specification (blocking decision)
- Wall ontology (explicit modeling stance)
- Flow support and back-pressure: required diagnostics
- 4+1D coupling constants and dimensional analysis (note)

## 0. Executive Summary

### 0.1 The physical hypothesis in one paragraph

We model a **time-dependent 4D throat/defect** that intermittently (or persistently) opens a corridor from the observable **3D brane** ((w!=!0)) into a **4D bulk** ((w) accessible). The ambient “vacuum” is treated as a **compressible superfluid substrate** whose density and flow can restructure near the defect. The throat is not assumed to be a passive, static object: it can be **actively maintained** by (i) a **4D electromagnetic (photon) cavity** whose standing modes exert stress that counteracts closure, and (ii) **superfluid intake/through-flow** whose momentum flux and pressure gradients contribute additional support. Because matter and fields can expand into the bulk on the far side of the throat, the **back-pressure and impedance** differ from purely 3D expectations. Observable EM and hydrodynamic signatures are taken to be **brane projections** of genuinely 4D interactions, so the model must specify an explicit projection/measurement map from bulk fields to brane observables.

---

### 0.2 What we must capture

This document exists to prevent us from “forgetting a piece” of the physical story. The full system must represent:

**Geometry and topology of the defect**

* A throat with time-dependent geometry (a(t)) (radius) and (L(t)) (length), and a consistent notion of “inside/outside” in 4D.
* A mechanism for the geometry to respond to forces (quasi-static constraints or explicit dynamics).

**Superfluid substrate (vacuum fluid)**

* A 4D dynamical field carrying density and flow, with a stiff EOS (the previously derived polytropic scaling).
* Intake/through-flow support: nonzero flux through/along the throat and the stresses it produces.

**Electromagnetic / photon sector in 4D**

* Actual **4+1D gauge fields** (A_M(\mathbf X,t)) (not just emergent potentials).
* A cavity/trapping mechanism so photon modes can localize in/near the throat and exert pressure/stress.
* A physically consistent source structure: charge neutrality of the ambient vacuum, with charge localized to excitations/defects.

**Bulk expansion and back-pressure**

* The receiving region in 4D must be represented so that outflow expands into greater dimensional volume, modifying density buildup and back-pressure.
* Boundary conditions/absorbers or reservoir modeling so the back-pressure is physical and not a numerical artifact.

**Brane observables**

* A precise map from 4D fields to 3D measurements (brane projection and/or restriction at (w!=!0)).
* Surface/port definitions and time-domain extraction rules for impedance-like response, with clear sign conventions.

**Energy bookkeeping and “E=mc² + density dependence”**

* A consistent energy budget including field energy, geometry energy, EM energy, and energy exchange.
* A way to represent energy shedding/redistribution when density changes (open-system or fixed-energy closure rather than a purely fixed-norm static state).

---

### 0.3 What success looks like

We will treat “success” as producing a *closed*, internally consistent equation stack with testable consequences, not merely a numerical equilibrium.

**Success criteria (theoretical)**

1. A unified variational formulation (S[\psi, A_M, a, L]) (or equivalent Hamiltonian) from which:

   * the gauged superfluid PDE,
   * the 4+1D Maxwell PDEs,
   * the geometry force laws,
   * and conservation/energy-balance identities
     all follow without contradictions.
2. A clearly stated coupling philosophy (we adopt **C1 minimal coupling** as the fundamental interaction, with earlier “hydrodynamic EM” reinterpreted as an effective reduction/limit).

**Success criteria (phenomenological / physical)**
3. A demonstrated stabilization mechanism for the throat that is genuinely dynamic:

* wave/cavity stress + flow stress can counteract closure,
* with explicit criteria for when the throat collapses or remains open.

4. A coherent explanation for why 4D interactions can appear as 3D EM/hydrodynamic signatures on the brane via projection.
5. A list of **derived effective parameters** (“new constants”) that emerge when the unified system is reduced to brane observables (analogous to how earlier work fixed specific exponents/coefficients).

**Success criteria (computational)**
6. Mathematica notebooks/scripts that:

* generate Euler–Lagrange equations,
* generate Noether currents and conservation laws,
* perform controlled reductions (mode decompositions in (w), brane effective theory),
* and produce check identities we can regression-test in Python solvers.

7. A minimal coupled simulation plan (even if crude) that can test:

   * cavity resonance and stress,
   * intake/back-pressure behavior,
   * stability of (a(t), L(t)),
   * and brane-observable response.

---

### 0.4 Deliverables and failure modes (explicitly)

**Primary deliverables**

* A master equation bundle (PDE/ODE system) with coupling terms and boundary/drive prescriptions.
* A Mathematica derivation pipeline and a parameter-extraction pipeline.
* A roadmap from the current “frozen operational definitions” to fully dynamic coupled runs.

**Explicit failure modes we want to detect early**

* The 4D Maxwell sector cannot reproduce brane-observed Coulomb behavior without added localization (then localization physics becomes mandatory and must be modeled explicitly).
* Charge neutrality cannot be implemented consistently with the chosen matter content (then we need a two-field matter model or defect-only charge model).
* The baseline geometry energy (E_{\rm geom}(a,L)=P_{\rm vac}V+\sigma A) cannot support the desired aspect ratio or dynamics (then shape/curvature terms or a different geometry model is required).
* The system only “works” with numerically fragile boundaries (then the back-pressure/outflow modeling is incorrect and must be redesigned).

---

## 1. Physical System Description

### 1.1 The arena: a 3D brane embedded in a 4D bulk

We assume the physical substrate exists in a **4D spatial bulk** with coordinates
[
\mathbf X=(x,y,z,w),\qquad t \text{ (time)}
]
and that ordinary observations correspond to a **3D brane slice** located near (w=0). In the far field (away from a throat/defect), the system is effectively brane-confined: dynamics in (w) are strongly suppressed so that physics appears 3D. Near a throat/defect, confinement in (w) is locally relaxed, allowing material and fields to explore the bulk coordinate (w).

Operationally, the distinction “3D vs 4D” is not philosophical; it is enforced by a **confinement structure** that determines whether the (w)-direction is energetically accessible.

---

### 1.2 The substrate: a compressible superfluid vacuum

The “vacuum” is modeled as a **compressible superfluid-like medium**. The medium has:

* a local density (amount of substrate available),
* a phase (which controls flow),
* a current/flux (which transports density and momentum),
* and an equation of state (EOS) that ties density to pressure/enthalpy.

The key physical point is that the throat is not supported by geometry alone. The **state of the substrate** (particularly density and flow) affects both:

1. how much energy is required to maintain an opening, and
2. how much stress/support the substrate can supply against closure forces.

This is where the “density matters” logic enters: if density drops, the available support changes and the system must respond (e.g., the defect relaxes, the geometry changes, energy is shed into waves/fields, etc.).

---

### 1.3 The throat/defect: a dynamically maintained corridor into the bulk

A throat is a localized region where the usual brane confinement is weakened, permitting motion and field structure in (w). Physically it behaves like a **corridor** or **puncture**:

* on the brane it appears as a localized mouth (a 2-sphere-like surface in 3D),
* in the bulk it has an extent along (w) (a “length”) and an effective transverse size (a “radius”),
* and it is treated as a *dynamical object* that can open, close, or oscillate.

We parameterize the throat geometry in the simplest reduced form by:
[
a(t)\ \text{(radius)},\qquad L(t)\ \text{(length)}.
]
This reduced model is not assumed to be ultimate truth; it is an operational starting point that will later be upgraded to a shape field (R(w,t)) if needed to capture flare/rounding physics and aspect-ratio selection.

---

### 1.4 Two stabilization mechanisms that must both be representable

We assume two physically distinct mechanisms can contribute to keeping the throat open. The model must allow each mechanism independently, and also allow them to coexist and interact.

#### 1.4.1 Resonant photon / cavity-field support (standing waves in 4D)

A central hypothesis is that there exists a **resonant cavity** associated with the throat—conceptually “photons” trapped in a standing-wave pattern. The key physical effects are:

* **Mode structure:** the cavity supports standing modes whose spectrum depends on geometry ((a,L)).
* **Field stress:** the stored field energy produces an effective stress (radiation pressure / Maxwell stress) that can push on the throat walls and endcaps, opposing closure.
* **4D nature:** the field lives in the 4D bulk; what we observe on the brane is a projection/slice.

This implies that any correct model must represent:

1. a genuine dynamical field in 4D (ultimately a gauge field),
2. a localization/trapping mechanism that makes a cavity possible, and
3. a coupling that lets field energy apply force to the geometry.

#### 1.4.2 Superfluid intake / through-flow support (non-equilibrium momentum flux)

In addition to wave/cavity support, the throat can be supported by **through-flow** of the substrate:

* When the throat opens, surrounding substrate can rush in (intake).
* This produces momentum flux and pressure gradients that can contribute to maintaining an opening (or destabilizing it, depending on regime).
* The system is inherently **non-equilibrium**: intake and outflow imply energy/mass exchange with a larger reservoir.

This implies the model must represent:

1. nonzero flux through the throat (especially along/through (w)),
2. a receiving region in the bulk so outflow can expand into 4D volume, and
3. a closure that permits exchange (open boundary conditions, reservoir coupling, or explicit port driving).

A purely static/closed equilibrium is therefore at best a limiting diagnostic; the full physical story is dynamic.

---

### 1.5 Bulk expansion and back-pressure: why 4D matters for impedance

A distinctive feature of the throat hypothesis is what happens *after* material passes through: it can expand into a **higher-dimensional receiving region**. In 4D, “available volume” increases more rapidly with distance than in 3D, and this changes the relationship between:

* outflow flux,
* density buildup in the receiving region,
* and the resulting back-pressure at the throat.

So the model must allow the following causal loop:

[
\text{throat opens} \Rightarrow \text{substrate flows into bulk}
\Rightarrow \text{bulk density distribution evolves}
\Rightarrow \text{back-pressure changes}
\Rightarrow \text{geometry dynamics respond}.
]

This is also where “impedance” language becomes meaningful: the effective (Z^{\rm eff}(\omega)) measured on the brane depends on how the bulk accepts and returns flux (dimensional expansion, trapping, dissipation, reflections).

---

### 1.6 Electromagnetic structure from the puncture: 4D interaction, 3D observation

We also hypothesize that the puncture/throat is associated with an **electromagnetic structure**—an effective “charge/field configuration” tied to the defect. The important point is not simply “there is an electric field,” but that:

* the interacting fields may be fully 4D, while
* the brane observer sees only a slice/projection of those interactions.

This requires:

1. a well-defined 4D gauge field sector,
2. a physically consistent sourcing mechanism (charge neutrality of the vacuum; charge localized to excitations/defects rather than the entire background), and
3. a projection/measurement rule that produces brane-observed (\mathbf E,\mathbf B) (or their operational substitutes) from the 4D fields.

This makes it possible for brane observables to reflect bulk behavior that is not obvious from 3D intuition.

---

### 1.7 Energy bookkeeping: (E=mc^2), density dependence, and energy shedding

A recurring physical constraint is the relationship between:

* the energy content attributed to a localized object (rest energy (mc^2)),
* the energy required to maintain an open throat,
* and the local substrate density that participates in supporting that throat.

The system cannot be treated as “rest-mass energy is fixed and geometry is fixed.” Instead, the hypothesis is:

* maintaining the throat consumes or stores energy in the throat region,
* the required energy depends on substrate density (availability and stiffness),
* if density decreases, the system’s capacity to sustain the throat decreases,
* the object may shed energy (into waves, EM radiation, bulk excitations, or geometry relaxation) to restore consistency.

This motivates a fully dynamic framework with explicit energy exchange, rather than a purely static closure.

---

### 1.8 Observables: what we measure on the brane

Finally, the physical system must connect to measurable quantities on the brane. Operationally we care about:

* **Brane density perturbations** (projected density changes near the mouth),
* **Flux through a measurement surface (\Gamma)** (mouth outflow/inflow),
* **Cavity signatures** (resonance peaks, phase lags),
* **Effective impedance / response operators** (Z^{\rm eff}(\omega)),
* **Field signatures** consistent with electric/magnetic behavior tied to the defect.

Because these are brane observables of a bulk system, the model must specify:

1. projection rules from 4D fields to 3D observables, and
2. a measurement protocol (ports, modes, time-domain extraction), so the theory makes concrete predictions.

---

### 1.9 Where this description points us

This physical description makes one conclusion unavoidable: the full model must be **fully dynamic** and must include an explicit 4D EM/gauge sector. Any intermediate “static equilibrium” calculations are only diagnostic subcases. The central mathematical task is therefore to assemble a unified coupled system:

* a charged superfluid PDE in 4D,
* a 4+1D gauge field PDE,
* and throat geometry dynamics,
  with a consistent energy budget, boundary/reservoir modeling, and brane projection map.

---

## 2. Modeling Choices and Assumptions

This section pins down the choices that convert the physical story (Section 1) into a mathematically closed system. The guiding principle is: **everything that couples must be derived from a single variational formulation** unless we explicitly label it as an effective/phenomenological add-on.

---

### 2.1 What is fundamental vs what is emergent

We separate two levels of description:

**Fundamental (this document’s target):**

* A **4D superfluid substrate** represented by a dynamical complex field (\psi(\mathbf X,t)).
* A genuine **4+1D gauge field** (A_M(\mathbf X,t)) (Maxwell sector in the bulk).
* A throat geometry described by **time-dependent DOFs** (a(t)), (L(t)) (and later, a shape field).

**Emergent/effective (prior work):**

* The earlier “EM-from-fluid” construction (Paper IV style) is treated as an **effective reduction** or approximation that may be recovered in a certain limit, but is **not** the defining EM model once we adopt C1. In the fundamental model, EM is not computed from fluid variables; it is solved as its own PDE and coupled by current and covariant derivatives.

This distinction is crucial: adopting C1 changes the “meaning” of the EM paper from “definition of EM” to “candidate effective mapping to be derived or falsified.”

---

### 2.2 Why we commit to fully dynamic evolution

We explicitly assume the throat is not a static object. The model must allow:

* time-dependent density (\rho(\mathbf X,t)),
* time-dependent flows (\mathbf J(\mathbf X,t)),
* time-dependent EM fields and stored cavity energy,
* and time-dependent geometry ((a(t),L(t))).

A static equilibrium may still exist as a special case, but the intended physical mechanisms (standing waves, intake support, bulk back-pressure) are fundamentally **non-equilibrium** and require dynamics.

Operational consequence: we must write evolution equations (PDEs/ODEs) and energy balance identities, not only stationary conditions.

---

### 2.3 Why we commit to 4D gauge fields (C1 minimal coupling)

We adopt **minimal coupling (C1)** as the primary interaction between the superfluid and the gauge sector. Concretely:

* (\psi) couples to the gauge potential (A_M) through covariant derivatives (D_M\psi).
* The gauge field is sourced by a conserved current (J^M) derived from (\psi).

This is chosen because:

* it is the most constraint-rich option (and therefore most likely to generate “unexpected” derived relations),
* it gives a clean route to stress/radiation pressure through the EM stress tensor,
* and it enables a coherent brane-projection picture: 4D interactions in (A_M) can yield 3D signatures.

---

### 2.4 Geometry representation choices

We adopt a staged geometry representation:

#### 2.4.1 Reduced DOFs ((a(t),L(t))) as a first dynamical geometry model

The initial geometry model uses only:

* radius (a(t)),
* length (L(t)),

and represents the throat region via a smooth confinement “gate” encoded in (V_{\rm conf}(\mathbf X;a,L)).

This is a pragmatic minimal model: it lets us couple geometry to field stresses without yet solving full moving-boundary PDEs.

#### 2.4.2 Planned upgrade: a shape field (R(w,t)) (flare/rounding physics)

We explicitly recognize that many geometric effects (preferred aspect ratio, flare, curvature costs) cannot be captured by a pure tube model with (V(a,L)) and (A(a,L)). A more physical geometry would be described by a function (R(w,t)) (or equivalent) with a geometric energy functional including curvature/rounding terms.

This upgrade is not optional if:

* the system requires an intrinsic preferred (L/a),
* or if endcap physics depends on non-tubular geometry.

For now, we keep ((a,L)) while keeping the shape-field extension visible.

---

### 2.5 Confinement strategy: geometry as a potential (Family-1 baseline)

We represent “where the throat is” and “what the walls are” using a **confinement potential** (V_{\rm conf}) rather than hard boundaries. The baseline choice is the existing **Family-1** pattern:

* a brane trap in (w) whose stiffness is modulated by a gate (G(\mathbf X;a,L)),
* a radial wall barrier (V_{\rm wall}(R_3;a)),
* and an endcap barrier (V_{\rm cap}(w;L)).

This approach is chosen because it:

* is differentiable in (a,L),
* yields clean force expressions (\int d^4X\,\rho\,\partial_{a,L}V_{\rm conf}),
* and is directly compatible with spectral PDE solvers and symbolic differentiation.

Alternative confinement families remain legitimate, but Family-1 is the frozen baseline for the combined system.

---

### Wall ontology (explicit modeling stance)

There are two distinct ontologies for what the “walls” are:

**(W1) Walls as geometry-as-potential (current hard-mode implementation):**
- The throat is defined by a smooth confinement potential \(V_{\rm conf}(\mathbf X;a,L)\).
- The energetic cost of sustaining a throat is represented separately by \(E_{\rm geom}(a,L)\).
- The wall is not a separate material interface; it is an energetic constraint.

**(W2) Walls as a material interface/defect in the fields (future upgrade):**
- The wall is a localized configuration (domain wall / soliton / topological defect) of the underlying field(s).
- \(a,L\) become collective coordinates (or are replaced by a shape field \(R(w,t)\)).
- Surface tension and curvature energies emerge from the field theory, rather than being prescribed.

For Paper 7 next steps we proceed with (W1) and treat (W2) as a future “shape-field” upgrade. Any claim about aspect-ratio selection must state whether it is expected to emerge from (W1) via EM/flow support or whether it requires (W2).

---

### 2.6 Geometry energy model: baseline and known limitation

We include a separate geometry energy (E_{\rm geom}) representing the energetic cost of sustaining a defect of size ((a,L)). The baseline form is the pressure+surface term model:
[
E_{\rm geom}(a,L)=P_{\rm vac}V(a,L)+\sigma A(a,L)\quad(+\ \text{optional higher terms}).
]

**Explicit limitation (kept in view):** with a tube-like (V) and (A), (\partial_L E_{\rm geom}) is typically independent of (L). Therefore, the baseline model does **not** by itself produce a preferred aspect ratio. If the physics requires aspect-ratio selection, it must come from:

* shape/curvature energy terms, and/or
* field contributions (EM cavity stress, intake stress), and/or
* constraints that reduce the independent DOFs.

---

### 2.7 Open vs closed system closures (and why this matters)

Because the system is dynamic and involves intake/outflow, a purely closed “fixed norm” closure is generally insufficient for the full model. We will consider closures explicitly rather than implicitly:

**Closed-system style (diagnostic-only):**

* fixed total norm (N=\int d^4X\,|\psi|^2) in an isolated box.

**Open-system styles (physically motivated):**

* fixed chemical potential (\mu) (reservoir coupling),
* fixed energy budget (H_{\rm tot}\approx mc^2) (rest-energy closure with density dependence),
* boundary/port driving plus dissipation (steady-state under pumping).

We will treat closure as part of the model specification, because it changes whether “density drops → energy sheds → geometry relaxes” can occur.

---

### 2.8 Gauge-field localization is not optional in 4D

A major consequence of true 4 spatial dimensions is that free Maxwell fields obey a different Gauss-law scaling than 3D intuition. If we want brane observers to see familiar 3D electromagnetic behavior, we will likely require **localization** of gauge modes near the brane or throat.

We therefore plan to include one of the following *variationally consistent* mechanisms:

* a position-dependent prefactor (Z(\mathbf X)) multiplying (F_{MN}F^{MN}) (effective medium),
* a Proca/Higgs-like mass (m_\gamma^2(\mathbf X)A_MA^M) that suppresses bulk propagation away from brane/throat,
* or explicit brane-localized terms that trap a zero mode.

This is a modeling choice we expect to become a derived constraint: the need for localization should be testable by matching brane-observed scaling laws.

---

### 2.9 Charge neutrality and sourcing strategy

If (\psi) is a vacuum condensate, naïvely setting (J^0=q|\psi|^2) would charge the entire vacuum, which is not acceptable. We therefore assume charge is **neutral in the far field** and **localized to defects/excitations**. We treat this as a required part of the combined system specification.

Acceptable strategies include:

* a neutralizing background (J^0=q(|\psi|^2-\rho_0)),
* a two-field model (neutral condensate + charged excitation field),
* defect-localized charge via gating (charge proportional to deviation from background within the throat region).

This choice strongly affects the derived effective EM sector and must be frozen explicitly before deriving final reduced equations.

---

### 2.10 What stays “separate modules” vs what is unified

We unify at the action level:

* (\psi), (A_M), and ((a,L)) belong to one coupled variational system.

We allow modularity only for:

* measurement/port extraction (a post-processing map from fields to observables),
* numerics (different discretizations for (\psi) vs (A_M)),
* and optional effective reductions (deriving brane-effective equations by integrating out (w)-modes).

In other words: the physics is unified; the workflow can be modular.

---

### 2.11 Summary of commitments

For the remainder of this note, we proceed with:

1. Fully dynamic evolution in 4D bulk.
2. True 4+1D gauge fields (A_M) with C1 minimal coupling to (\psi).
3. Dynamic throat geometry via ((a(t),L(t))) (with an explicit path to shape-field upgrades).
4. Explicit open-system closure (reservoir/drive/dissipation) as part of the model, not an afterthought.
5. A localization mechanism for gauge modes as a necessary ingredient to obtain correct brane-observed EM scaling.
6. Explicit charge neutrality strategy so the vacuum is not uniformly charged.

---

## 3. Master Variables, Domains, and Projection Map

This section defines the objects that exist in the unified model and the operational maps that turn 4D bulk dynamics into 3D brane observables.

---

### 3.1 Coordinates, domains, and index conventions

We work on a **4D spatial bulk** with coordinates
[
\mathbf X=(x,y,z,w),\qquad t\in\mathbb R.
]
We use:

* spatial indices (i,j\in{x,y,z,w}),
* spacetime indices (M,N\in{0,x,y,z,w}) where 0 denotes time,
* (\nabla_4 = (\partial_x,\partial_y,\partial_z,\partial_w)).

The **brane** is operationally associated with the neighborhood of (w=0). In computation and measurement, “brane observables” are either:

* evaluated at (w=0) (restriction), or
* projected with a weight (W(w)) that represents the brane-localized mode structure (projection).

---

### 3.2 Dynamical fields and degrees of freedom (the full state)

The fundamental state variables are:

#### (i) Superfluid/order parameter field

[
\psi(\mathbf X,t)\in\mathbb C.
]
Key derived quantities:
[
\rho(\mathbf X,t)=|\psi|^2,\qquad
\mathbf J_{\rm mass}(\mathbf X,t)=\frac{\hbar}{m}\Im!\left(\psi^* D\psi\right),\qquad
\mathbf J_{\rm ch}(\mathbf X,t)=q\,\mathbf J_{\rm mass}(\mathbf X,t)
]
where (D) is the gauge-covariant spatial derivative (defined in Section 4/5).
If we adopt the two-field neutrality option (N2), (\mathbf J_{\rm ch}) is instead computed from (\chi).

#### (ii) 4+1D gauge field

[
A_M(\mathbf X,t)\in\mathbb R,\quad M\in{0,x,y,z,w}
]
with field strength
[
F_{MN}=\partial_M A_N-\partial_N A_M.
]
Brane-observed EM-like quantities will be obtained by restriction or projection of components of (F_{MN}) (see §3.6).

#### (iii) Throat geometry degrees of freedom

We use the reduced geometry:
[
a(t)>0,\qquad L(t)>0
]
representing an effective radius and length of the throat corridor. These may be treated quasi-statically ((a,L) solve instantaneous force balance) or dynamically via a wall-law ODE; the choice is part of the coupling plan (Section 5).

---

### Force ledger for throat DOFs (non-negotiable bookkeeping)

We will track the generalized forces on the throat geometry DOFs explicitly as a **sum of channels**:

\[
F_a(t) \equiv -\frac{\partial H_{\rm tot}}{\partial a},\qquad
F_L(t) \equiv -\frac{\partial H_{\rm tot}}{\partial L}.
\]

With the unified model these decompose as:

\[
F_{a,L} = F^{(\psi)}_{a,L} + F^{(\rm EM)}_{a,L} + F^{(\rm flow)}_{a,L} + F^{(\rm geom)}_{a,L} + F^{(\rm loc)}_{a,L}.
\]

Definitions (each term must be either derived from the Hamiltonian or computed from a stress tensor consistently):

- Matter/confinement force:
  \[
  F^{(\psi)}_{a,L} = -\int d^4X\,\rho(\mathbf X,t)\,\partial_{a,L}V_{\rm conf}(\mathbf X;a,L).
  \]
- Geometry cost:
  \[
  F^{(\rm geom)}_{a,L} = -\partial_{a,L}E_{\rm geom}(a,L).
  \]
- EM cavity / gauge sector contribution (from the bare Maxwell term):
  \[
  F^{(\rm EM)}_{a,L} = -\partial_{a,L}H_{\rm EM}[A;\,a,L].
  \]
- Flow / momentum-flux support (only present in driven/open dynamics):
  \[
  F^{(\rm flow)}_{a,L} \text{ is extracted from the matter stress tensor or equivalent flux diagnostics (see below).}
  \]
- Localization forces (from \(Z(\mathbf X)\) and/or \(m_\gamma(\mathbf X)\) terms):
  \[
  F^{(\rm loc)}_{a,L} = -\partial_{a,L}H_{\rm loc}[A;\,a,L] \quad (\text{if localization is implemented by } Z(\mathbf X) \text{ or } m_\gamma(\mathbf X)).
  \]

**Rule:** any claim about “the throat is supported” must specify which subset of these terms balances which (e.g. EM stress balances \(\partial_L E_{\rm geom}\), intake balances radial collapse, etc.).

---

### 3.3 Geometry gate and throat “interior” definition

To define “inside the throat” in a differentiable way, we introduce a gate function
[
G(\mathbf X; a, L)\in[0,1]
]
that is near 1 in the interior region and near 0 outside. In the frozen baseline, (G) is separable into a 3D radial factor and a (w)-extent factor:
[
G(\mathbf X;a,L)=G_r(R_3;a)\,G_w(w;L),\qquad R_3=\sqrt{x^2+y^2+z^2}.
]

Operationally:

* (G_r) transitions near (R_3\approx a) with width (\delta_r),
* (G_w) transitions near (|w|\approx L/2) with width (\delta_\parallel),
* and smoothing functions ensure (\partial_a G), (\partial_L G) exist.

This gate is the common “shape function” that modulates:

* the confinement potential (V_{\rm conf}),
* (optionally) gauge localization terms,
* and any defect-localized source terms (charge neutrality strategy).

---

### 3.4 Confinement potential and the “brane vs bulk accessibility” map

We encode brane confinement vs bulk accessibility in a potential
[
V_{\rm conf}(\mathbf X; a,L)
]
whose role is:

* far from the defect: make motion into (w) energetically expensive (brane-like physics),
* in the defect region: relax confinement so that (w) becomes accessible (bulk corridor),
* and impose wall/endcap barriers that define radius and length.

In the baseline Family-1 interpretation, (V_{\rm conf}) includes:

* a (w)-trap term whose stiffness depends on (G),
* a radial wall barrier depending on (R_3-a),
* an endcap barrier depending on (|w|-L/2).

We treat (V_{\rm conf}) as the operational representation of “throat geometry in the bulk,” and it is the primary mechanism by which ((a,L)) couple into the field PDEs.

---

### 3.5 Brane projection operator (\mathcal P)

Because measurements are made on the brane, we define a map from bulk fields to brane observables.

There are two operationally valid projection styles:

#### (A) Restriction (“slice”)

Evaluate bulk fields at (w=0):
[
\mathcal P_0[f](x,y,z,t)=f(x,y,z,w=0,t).
]

This is appropriate when the brane is treated as a strict hypersurface.

#### (B) Weighted projection (“mode overlap”)

Project along (w) with a nonnegative weight (W(w)), normalized appropriately:
[
\mathcal P_W[f](x,y,z,t)=\int_{-\infty}^{\infty} W(w)\,f(x,y,z,w,t)\,dw
]
with (W(w)) chosen to represent the brane-localized far-field ground mode (e.g., the harmonic ground state of the far-field trap).

Operationally, the projection used in prior freezing includes:

* a finite (w)-window truncation,
* renormalization of the truncated weight,
* and a tail-mass diagnostic to ensure projection validity.

This projection map is essential to correctly represent “4D interactions that appear as 3D observables.”

---

### 3.6 Brane-observed EM and fluid observables (what we claim is measurable)

#### 3.6.1 Brane density and flow observables

From (\psi), define bulk density (\rho=|\psi|^2) and mass/number current (\mathbf J_{\rm mass}). The charge current is (\mathbf J_{\rm ch}=q\,\mathbf J_{\rm mass}) unless the two-field option (N2) is used. The brane-observed density is:
[
\rho_{\rm brane}(x,y,z,t)=\mathcal P[\rho](x,y,z,t)
]
using either restriction or weighted projection.

Brane-observed flux across a brane surface is derived from the brane-projected current and the outward normal of the measurement surface (see §3.7).

#### 3.6.2 Brane EM observables from 4D gauge fields

From the gauge field define (F_{MN}). Brane-observed EM-like components are obtained via (\mathcal P):

* electric-type components on the brane:
  [
  E_i^{\rm brane}(x,y,z,t)=\mathcal P[F_{0i}](x,y,z,t),\quad i\in{x,y,z},
  ]
* magnetic-type components on the brane (defined from spatial components (F_{ij}) restricted/projection to the brane, then mapped to a 3D pseudovector in the usual way).

Because the underlying field is 4D, additional components involving (w) (e.g. (F_{0w}), (F_{iw})) may be present. Whether these are directly observable or only indirectly inferred is part of the model-to-observable map. The important point is that the brane sees a projection of (F), not necessarily the full bulk structure.

---

### 3.7 Measurement surface (\Gamma), ports, and modal observables

To connect with the operational measurement framework, we define a brane surface (\Gamma) representing the throat mouth on the brane.

Baseline choice:

* (\Gamma) is a 2-sphere on the brane at radius (R_3=r_{\rm port}) and (w=0),
* with outward normal (\hat n=\hat R_3),
* and surface measure (d\mu=r_{\rm port}^2\sin\theta\,d\theta\,d\phi).

We then define a port basis ({P_i(\theta,\phi)}) (e.g., spherical harmonics) and modal observables by surface integration:

* “effort” modes (u_i(t)) built from a brane thermodynamic variable (u(\mathbf s,t)) (typically a linearized enthalpy/pressure proxy derived from (\rho_{\rm brane})),
  [
  u_i(t)=\int_\Gamma \overline{P_i(\mathbf s)}\,u(\mathbf s,t)\,d\mu
  ]
* “flux” modes (j_i(t)) built from normal flux through (\Gamma),
  [
  j_i(t)=\int_\Gamma \overline{P_i(\mathbf s)}\,\big(\mathbf J_{\rm brane}\cdot \hat n\big)(\mathbf s,t)\,d\mu
  ]

These are the time series used later for frequency-domain response/impedance extraction.

---

### 3.8 Drive, open-system coupling, and dissipation placeholders

Because the full system is dynamic and open, we reserve explicit slots in the model for:

* **drive/injection** (ports, external sources (S_{\rm drive}), boundary forcing),
* **dissipation** (localized damping layers, resistivity-like terms, radiation losses),
* **reservoir closure** (fixed (\mu) or fixed energy budget targets).

We do not fix these here, but we treat them as first-class parts of the system specification, because they determine whether intake, back-pressure, and steady-state cavity behavior can exist.

---

### 3.9 Summary: the full state and the observable map

The “master state” is:
[
\boxed{\ \big(\psi(\mathbf X,t),\ A_M(\mathbf X,t),\ a(t),\ L(t)\big)\ }
]
with geometry encoded into the fields through the gate and confinement potential.

The brane-observable map is:
[
\boxed{\ \text{bulk fields} \xrightarrow{\ \mathcal P\ \text{and}\ \Gamma,{P_i}\ }\ \rho_{\rm brane},\ u_i(t),\ j_i(t),\ E^{\rm brane},\ B^{\rm brane},\ldots\ }
]

---

## 4. Unified Variational Formulation

This section defines the single “master” variational object from which the coupled PDE/ODE system will be derived. The goal is that **all couplings and force laws are consequences** of one action (or equivalently one Hamiltonian), rather than an ad hoc assembly.

---

### 4.1 Master action: fields, geometry, and coupling in one place

We define a total action
[
S[\psi,A_M,a,L] = \int dt\int d^4X\,\Big(\mathcal L_\psi+\mathcal L_{\rm EM}+\mathcal L_{\rm loc}+\mathcal L_{\rm gf}+\mathcal L_{\rm drive/diss}\Big) - \int dt\,E_{\rm geom}(a,L)
]
where:

* (\psi(\mathbf X,t)) is the superfluid/order parameter field,
* (A_M(\mathbf X,t)) is the 4+1D gauge field,
* ((a(t),L(t))) are geometry DOFs entering through (V_{\rm conf}) and possibly through localization terms,
* (E_{\rm geom}(a,L)) is the geometry energy cost functional,
* (\mathcal L_{\rm drive/diss}) stands for explicit open-system coupling terms (kept modular but variationally consistent where possible).

We work in a nonrelativistic matter sector (GNLS-type) coupled to a relativistic-form Maxwell sector in 4+1D. This hybrid is standard in condensed-matter-like effective field theories and is sufficient for our goals (a fully relativistic matter sector can be pursued later if needed).

---

### 4.2 Minimal coupling (C1) and covariant derivatives

We adopt minimal coupling with charge parameter (q). Define
[
D_t \psi \equiv \left(\partial_t + \frac{i q}{\hbar}A_0\right)\psi,
\qquad
D_i \psi \equiv \left(\partial_i - \frac{i q}{\hbar}A_i\right)\psi,\quad i\in{x,y,z,w}.
]
The gauge field strength is
[
F_{MN}=\partial_M A_N-\partial_N A_M.
]

This choice ensures gauge invariance under:
[
\psi\to e^{-iq\lambda/\hbar}\psi,\qquad
A_M\to A_M+\partial_M\lambda.
]

---

### Neutrality / source specification (blocking decision)

Minimal coupling produces a conserved 4+1D current \(J^M\) from the matter field(s). However, we must avoid a uniformly charged vacuum. Therefore **\(J^0\) must be defined operationally** and frozen before claiming a physical EM model.

We will choose exactly one of the following and treat it as a spec-level decision:

**(N1) Neutralizing background**
\[
J^0 = q\left(|\psi|^2-\rho_0\right),\quad
J^i = \frac{q\hbar}{m}\Im(\psi^* D_i\psi).
\]
(\(\rho_0\) is the far-field background density.)

**(N2) Two-field model (recommended if we want a truly neutral condensate)**
- \(\psi\): neutral condensate controlling throat/geometry
- \(\chi\): charged excitation field sourcing EM
\[
J^0 = q|\chi|^2,\quad
J^i = \frac{q\hbar}{m_\chi}\Im(\chi^* D_i\chi).
\]

**(N3) Defect-gated charge (puncture-localized source)**
\[
J^0 = q\,G_{\rm src}(\mathbf X;a,L)\,\left(|\psi|^2-\rho_0\right)
\]
with \(G_{\rm src}\) a frozen gate localizing charge to the defect/throat region.

**Rule:** do not proceed to Mathematica “full coupled derivations” until one of (N1–N3) is selected and frozen.

---

### 4.3 Superfluid sector Lagrangian (gauged GNLS)

We take the matter Lagrangian density as
[
\mathcal L_\psi =
\frac{i\hbar}{2}\Big(\psi^* D_t\psi - \psi (D_t\psi)^*\Big)
-\frac{\hbar^2}{2m}|D\psi|^2
- V_{\rm conf}(\mathbf X;a,L)|\psi|^2
- U(|\psi|^2).
]

Here:

* (m) is the effective mass parameter,
* (V_{\rm conf}) is the confinement potential encoding the throat geometry (defined in §4.6),
* (U(\rho)) encodes the EOS (with (\rho=|\psi|^2)). Our stiff polytrope choice corresponds to a specific (U) such that (U'(\rho)\psi) reproduces the frozen GNLS nonlinearity (the “(n=5)” style scaling from earlier work).

This Lagrangian is chosen because it yields:

* the gauged GNLS equation for (\psi),
* a conserved gauge current as a Noether current,
* and a clear energy density and stress tensor for momentum/force derivations.

---

### 4.4 Electromagnetic sector in 4+1D (Maxwell + gauge fixing)

We take a Maxwell Lagrangian in the bulk:
[
\mathcal L_{\rm EM} = -\frac{1}{4\mu_0}F_{MN}F^{MN}.
]

To make the PDEs well-posed for symbolic/numerical work, we include a covariant gauge-fixing term (Lorenz gauge family):
[
\mathcal L_{\rm gf} = -\frac{1}{2\xi\mu_0}(\partial_M A^M)^2
]
where (\xi) is a gauge parameter ((\xi=1) is the simplest “Feynman gauge” choice).

This yields Maxwell’s equations in a form that can be directly derived and manipulated in Mathematica.

---

### 4+1D coupling constants and dimensional analysis (note)

In 4+1D the gauge coupling has different mass dimensions than in 3+1D. We therefore treat \(\mu_0\) (or equivalently a 5D coupling \(g_5\)) as an **effective parameter** to be matched after localization and brane reduction. All “observed 3D EM” statements must be made in terms of the localized brane effective coupling extracted from the zero-mode reduction, not by importing 3D dimensional intuition.

---

### 4.5 Gauge localization / trapping (required in 4D)

A key modeling requirement is that the brane observer should not generically see “4D Coulomb.” We therefore introduce a localization mechanism as a variational term. We keep it general in the master action so that different localization strategies can be compared.

Two standard, variationally consistent options are:

#### Option L1: Position-dependent EM prefactor (“dielectric”/warp factor)

[
\mathcal L_{\rm loc}^{(Z)} = -\frac{1}{4\mu_0}\big(Z(\mathbf X;a,L)-1\big)F_{MN}F^{MN}.
]
Here (Z(\mathbf X;a,L)) is large/small in regions that suppress bulk propagation and can be tied to the throat gate (G(\mathbf X;a,L)) and/or a brane-localizing profile. With this convention, (Z=1) gives no localization term and avoids double-counting the Maxwell kinetic energy.

#### Option L2: Proca/Higgs-like mass term (“bulk suppression”)

[
\mathcal L_{\rm loc}^{(m)} = +\frac{1}{2\mu_0}m_\gamma^2(\mathbf X;a,L)A_M A^M.
]
Here (m_\gamma(\mathbf X;a,L)) is chosen to be large away from the brane/throat, trapping low-energy gauge modes in the relevant region.

We will treat (Z) or (m_\gamma) as part of the model’s “derived-constant” pipeline: their form determines effective 3D EM behavior after reduction, and may be constrained by matching brane-observed scaling laws.

---

### 4.6 Confinement potential (V_{\rm conf}(\mathbf X;a,L)) and the gate

The throat geometry enters the action primarily through (V_{\rm conf}), which is assumed smooth in both (\mathbf X) and ((a,L)). The baseline Family-1 structure is:

[
V_{\rm conf}(\mathbf X;a,L)=V_{\rm brane}(\mathbf X;a,L)+V_{\rm wall}(R_3;a)+V_{\rm cap}(w;L),
\qquad R_3=\sqrt{x^2+y^2+z^2}.
]

* (V_{\rm brane}) enforces brane localization far from the throat and relaxes it inside the gated region (so the bulk becomes accessible).
* (V_{\rm wall}) imposes a radial barrier around the throat radius.
* (V_{\rm cap}) imposes endcaps near (|w|\approx L/2).

The gate (G(\mathbf X;a,L)) is a smooth interior indicator built from smooth-step functions, ensuring (\partial_a V_{\rm conf}) and (\partial_L V_{\rm conf}) exist. This is essential for clean geometry force derivations.

---

### 4.7 Geometry energy functional (E_{\rm geom}(a,L))

We include a geometry energy cost representing the energetic burden of sustaining a throat of size ((a,L)). The baseline model is:

[
E_{\rm geom}(a,L)=P_{\rm vac}V(a,L)+\sigma A(a,L)
\quad(+\text{optional corrections}).
]

We explicitly record the limitation: this tube-based (V,A) model typically does not select a preferred aspect ratio (L/a) by itself. If aspect ratio selection is required, the variational framework supports adding:

* curvature/rounding terms,
* shape-field energies (E_{\rm geom}[R(w)]),
* or coupling-dependent effective geometry terms from integrating out field modes.

We treat these as planned upgrades rather than hidden assumptions.

---

### 4.8 Charge neutrality strategy as part of the action

Minimal coupling implies a matter current that sources Maxwell. However, if the ambient vacuum condensate (\psi) is charged, the far field would be uniformly charged, which we do not want. Therefore, the action must encode a charge-neutral background.

Two compatible strategies (to be frozen later) are:

* **Neutralizing background:** replace the charge density by deviation from a reference density, e.g. (J^0 \propto (|\psi|^2-\rho_0)). This can be implemented by defining the conserved current with a subtraction or by coupling the EM field only to deviations.
* **Two-field matter sector:** keep (\psi) neutral (geometry/throat substrate) and introduce a separate charged excitation field (\chi) that couples to (A_M). This keeps the vacuum neutral by construction.

We keep the master action flexible enough to support either approach.

---

### 4.9 Drive and dissipation terms (open system)

To represent intake, cavity pumping, radiation losses, and steady states, we include (\mathcal L_{\rm drive/diss}) as a placeholder for terms such as:

* external currents (J^M_{\rm ext}) that drive (A_M),
* reservoir coupling/damping for (\psi),
* absorbing boundary layers for EM and matter waves.

Where possible we will encode these in a way consistent with conservation laws (e.g., dissipation represented by explicit non-Hamiltonian terms whose contribution to energy balance is tracked).

The key modeling commitment is that the system is not closed; drive/dissipation are part of the physical description and therefore part of the full equation stack.

---

### 4.10 What this formulation buys us

With this unified variational formulation:

1. The **gauged GNLS equation** for (\psi) follows from (\delta S/\delta\psi^*=0).
2. The **4+1D Maxwell equations** follow from (\delta S/\delta A_M=0), with sources fixed by Noether’s theorem.
3. The **geometry forces** follow from (-\partial_{a,L}H_{\rm tot}) (or from varying the action w.r.t. (a,L)), yielding clean terms like \(-\int d^4X\,\rho(\mathbf X,t)\,\partial_{a,L}V_{\rm conf}(\mathbf X;a,L)\) plus EM localization and geometry-energy contributions.
4. Conservation laws (charge continuity, energy balance, momentum flux) can be derived systematically, giving the compact “encapsulating” identities we expect to reveal hidden constraints and new derived parameters.

---

## 5. Governing Equations: Full Coupled PDE/ODE System

This section lists the explicit evolution equations implied by the unified action in Section 4, together with the key derived quantities (current, continuity, geometry forces) and the “open system” hooks (drive/dissipation) required to represent intake, resonance pumping, and energy shedding.

---

### 5.1 Superfluid matter PDE: gauged GNLS in 4D bulk

Let (\rho=|\psi|^2). With minimal coupling, define
[
D_t = \partial_t + \frac{i q}{\hbar}A_0,\qquad
D_i = \partial_i - \frac{i q}{\hbar}A_i,\quad i\in{x,y,z,w}.
]

The matter equation of motion is:
[
i\hbar D_t\psi =
\left[
-\frac{\hbar^2}{2m}D_iD_i
+V_{\rm conf}(\mathbf X;a,L)
+U'(\rho)
\right]\psi
+\mathcal N_{\rm drive/diss}[\psi]
]
where:

* (U'(\rho)) is the EOS/nonlinearity (frozen from the prior “stiff polytrope” choice),
* (V_{\rm conf}) is the throat confinement potential (Family-1 baseline),
* (\mathcal N_{\rm drive/diss}) represents optional non-Hamiltonian terms used to model reservoir coupling, damping, or injection.

**Remarks**

* In the purely conservative limit (\mathcal N_{\rm drive/diss}=0), this equation is Hamiltonian and admits standard energy and charge conservation identities.
* Intake/through-flow and energy shedding generally require (\mathcal N_{\rm drive/diss}\neq 0) and/or boundary forcing.

---

### 5.2 Conserved current and continuity (Noether)

Gauge invariance implies a conserved **charge** current. We keep mass/number and charge currents distinct:

Mass/number current from (\psi):
[
J_{\rm mass}^i = \frac{\hbar}{m}\Im\!\big(\psi^* D_i\psi\big),
\qquad i\in{x,y,z,w}.
]

Charge current sourcing Maxwell:
[
J_{\rm ch}^i = q\,J_{\rm mass}^i
]
or (J_{\rm ch}^i) derived from (\chi) if the two-field option (N2) is chosen.

The charge density is (J_{\rm ch}^0) (see neutrality strategy below). In the simplest charged-matter form it is (J_{\rm ch}^0=q|\psi|^2), but we do **not** generally accept that for a vacuum condensate.

The continuity equation is:
[
\partial_M J_{\rm ch}^M = 0
]
up to explicit source/sink terms if we introduce reservoir coupling that violates strict gauge-charge conservation (in which case the violation must be tracked and physically justified).

---

### 5.3 Charge neutrality and the physical form of (J^0)

To avoid a uniformly charged vacuum, we must choose one consistent implementation. We keep the PDE form general:

* **Neutralizing background option:**
  [
  J^0 = q\left(|\psi|^2-\rho_0\right),
  ]
  with (\rho_0) the far-field condensate density. This makes total charge vanish in the homogeneous vacuum.

* **Two-field option (recommended if needed):**
  Introduce a neutral condensate (\psi) (geometry/throat) and a charged excitation field (\chi) (EM sources):
  [
  J^M = J^M[\chi],
  ]
  while (\psi) couples to geometry and possibly to EM only indirectly (or through higher-order terms).

Which option we choose affects both Maxwell sourcing and what “charge” means in the throat. The system cannot be fully specified without freezing this choice.

---

### 5.4 Maxwell equations in 4+1D with localization and drive

Define (F_{MN}=\partial_MA_N-\partial_NA_M). The gauge field equations are:
[
\partial_M\!\left(Z(\mathbf X;a,L)\,F^{MN}\right)
+ m_\gamma^2(\mathbf X;a,L)\,A^N
+ \frac{1}{\xi}\,\partial^N(\partial_M A^M)
= \mu_0\,(J^N+J^N_{\rm ext}).
]

Here:

* (Z(\mathbf X;a,L)) is an optional localization prefactor (if used),
* (m_\gamma(\mathbf X;a,L)) is an optional Proca/Higgs-like localization mass (if used),
* (\xi) is the gauge-fixing parameter,
* (J^N) is the matter current (from (\psi) or from (\chi) if two-field),
* (J^N_{\rm ext}) represents external drive/current injection used to pump cavity modes or represent imposed EM forcing.
* If we freeze (LZ), set (m_\gamma=0); if we freeze (LM), set (Z=1).

### 4+1D Maxwell equations with localization (explicit PDE form)

We will treat the gauge field as **dynamical in 4+1D**:
\[
A_M(\mathbf X,t),\quad M\in\{0,x,y,z,w\},\qquad F_{MN}=\partial_M A_N-\partial_N A_M.
\]

Because naive 4D Gauss-law scaling is not acceptable for brane-observed EM, the gauge sector must be **localized**. We will implement localization in one of the following variationally well-posed forms (both are allowed; pick one and freeze it):

**(LZ) Position-dependent kinetic prefactor**
\[
\partial_M\!\left(Z(\mathbf X;a,L)\,F^{MN}\right)
+\frac{1}{\xi}\,\partial^N(\partial\!\cdot\!A)
= \mu_0\,(J^N+J^N_{\rm ext}).
\]

**(LM) Position-dependent Proca/Higgs mass**
\[
\partial_M F^{MN}
+ m_\gamma^2(\mathbf X;a,L)\,A^N
+\frac{1}{\xi}\,\partial^N(\partial\!\cdot\!A)
= \mu_0\,(J^N+J^N_{\rm ext}).
\]

Here \(\xi\) is the gauge-fixing parameter (Lorenz gauge if \(\xi=1\)), and the dependence of \(Z\) or \(m_\gamma\) on \((a,L)\) is what allows EM stress to push on the throat geometry through the force ledger.
We keep gauge-fixing unweighted by (Z) for simplicity; if numerical stiffness appears, consider the alternative (\mathcal L_{\rm gf}\propto -Z(\partial\!\cdot\!A)^2).

**Cavity interpretation:** the throat and/or brane localization profiles (Z) and/or (m_\gamma) are what allow standing modes (a cavity spectrum). Without localization, 4D Maxwell tends to leak into the bulk and generically produces different far-field scaling than a brane observer expects.

---

### 5.5 Geometry evolution: forces from total energy (dynamic throat)

Let the total Hamiltonian (energy) be:
[
H_{\rm tot} = H_\psi + H_{\rm EM} + H_{\rm loc} + E_{\rm geom}
]
where each term corresponds to Section 4.
Here (H_{\rm EM}) is the **bare Maxwell** energy from (\mathcal L_{\rm EM}), and (H_{\rm loc}) is the energy from localization terms ((Z) prefactor and/or (m_\gamma) mass term). This partition is used consistently in the force ledger; if localization is folded into (H_{\rm EM}), then set (H_{\rm loc}=0) and drop (F^{(\rm loc)}).

The generalized forces on geometry are:
[
F_a = -\frac{\partial H_{\rm tot}}{\partial a},\qquad
F_L = -\frac{\partial H_{\rm tot}}{\partial L}.
]
In the presence of explicit drive/dissipation, these are the **conservative** generalized forces; nonconservative contributions must be added separately through the chosen closure/drive model.

Expanding the matter contribution yields the operationally important identities:
[
F_a^{(\psi)} = -\int d^4X\,\rho(\mathbf X,t)\,\partial_a V_{\rm conf}(\mathbf X;a,L)
]
[
F_L^{(\psi)} = -\int d^4X\,\rho(\mathbf X,t)\,\partial_L V_{\rm conf}(\mathbf X;a,L)
]

EM/localization terms contribute additional geometry forces whenever (Z(\mathbf X;a,L)) and/or (m_\gamma(\mathbf X;a,L)) depend on ((a,L)). For example, a Proca term contributes:
[
F_a^{(m)} = -\int d^4X\,\frac{1}{2\mu_0}\,(\partial_a m_\gamma^2)\,A_MA^M
]
and similarly for (L).

The geometry energy contributes:
[
F_a^{(\rm geom)} = -\partial_a E_{\rm geom}(a,L),\qquad
F_L^{(\rm geom)} = -\partial_L E_{\rm geom}(a,L).
]

We then choose the geometry law:

* **Quasi-static geometry (instantaneous force balance):**
  [
  F_a=0,\qquad F_L=0.
  ]
* **Dynamic wall-law (recommended for full dynamics):**
  [
  M_a\ddot a + C_a\dot a = F_a,\qquad
  M_L\ddot L + C_L\dot L = F_L,
  ]
  where (M_{a,L}) and (C_{a,L}) are effective inertia and damping coefficients (possibly derived from or fit to microphysics).

This is where “non-static throat” is enforced directly.

---

### 5.6 Energy balance: how energy moves between sectors

In the conservative, closed limit (no external drive/dissipation), the unified system has:
[
\frac{d}{dt}H_{\rm tot}=0.
]

With drive and dissipation, the balance takes the schematic form:
[
\frac{d}{dt}H_{\rm tot} = P_{\rm drive} - P_{\rm diss} - \Phi_{\rm boundary}
]
where:

* (P_{\rm drive}) is power injected by (J^N_{\rm ext}) and/or reservoir coupling in the matter sector,
* (P_{\rm diss}) is explicitly modeled damping (e.g. absorbing layers, resistivity-like terms),
* (\Phi_{\rm boundary}) is net **outward** energy flux through the computational boundary.

**Physical interpretation:** this is the formal location where “density drops → energy sheds” lives. A drop in local density can correspond to energy leaving the throat region via EM radiation, matter waves, or geometry relaxation, depending on which flux terms dominate.

---

### 5.7 Intake / through-flow and back-pressure: required open-system structure

To represent intake support and 4D receiving-region back-pressure, the model must not be purely closed. At minimum, we require:

1. A means to generate sustained flow (nonzero (J^w) and/or mouth flux through (\Gamma)):

   * boundary/port forcing in the matter sector,
   * reservoir coupling at far-field boundaries,
   * or explicit source terms.

2. A receiving region in the bulk in which (\rho(\mathbf X,t)) can evolve under EOS pressure:

   * the confinement must allow expansion into (w) away from the throat,
   * boundary conditions must prevent unphysical reflections (absorbing layers or large domains).

3. A diagnostic definition of “back-pressure” at the throat:

   * operationally defined via pressure/enthalpy from (\rho),
   * and related to flux by a response function (impedance) once linearized or measured.

### Flow support and back-pressure: required diagnostics

Because intake/through-flow support is intrinsically dynamic, we must define at least one operationally computable diagnostic that captures it.

We will include both a **flux diagnostic** and a **back-pressure diagnostic**:

**(Dflow) Through-flow flux**
\[
\Phi_w(t) \equiv \int_{\Sigma_w} J_{{\rm mass},w}(\mathbf X,t)\,d^3S,
\]
where \(\Sigma_w\) is a frozen cross-section (e.g. a brane-centered mouth plane or a throat interior plane). Nonzero \(\Phi_w\) is the signature that the throat is in an intake/outflow regime.
If charge transport is needed, define \(\Phi^{(\rm ch)}_w\) analogously using \(J_{{\rm ch},w}\).

**(Dback) Back-pressure / enthalpy drop**
Given EOS \(P(\rho)\) and corresponding enthalpy proxy \(h(\rho)\), define:
\[
\Delta h(t) \equiv \langle h(\rho)\rangle_{\rm near\ mouth} - \langle h(\rho)\rangle_{\rm far\ bulk}.
\]
This directly encodes the “4D receiving volume changes back-pressure” effect: dilution in the bulk lowers \(\langle\rho\rangle_{\rm far}\) and thus changes \(\Delta h\).

**Rule:** if we claim “flow holds the throat open,” we must show correlation between \(F^{(\rm flow)}_{a,L}\) (from the force ledger) and either \(\Phi_w(t)\) or \(\Delta h(t)\) under controlled drive conditions.

These are not optional features; they are the mathematical embodiments of the physical statements in Section 1.5.

---

### 5.8 Brane observables and response extraction (interface to measurement)

Given (\psi) and (A_M), brane observables are constructed via the projection operator (\mathcal P) and measurement surface (\Gamma):

* (\rho_{\rm brane}=\mathcal P[|\psi|^2]),
* brane flux through (\Gamma): (j(\mathbf s,t)=\mathbf J_{\rm brane}\cdot\hat n),
* EM observables: (E_i^{\rm brane}=\mathcal P[F_{0i}]), and (B^{\rm brane}) from spatial components.

Port-mode time series:
[
u_i(t)=\int_\Gamma \overline{P_i}u\,d\mu,\qquad
j_i(t)=\int_\Gamma \overline{P_i}(\mathbf J_{\rm brane}\cdot\hat n)\,d\mu
]

Response extraction (lock-in, conditioning) is treated as a post-processing map:
[
\mathbf j(\omega)=Z^{\rm eff}(\omega)\,\mathbf u(\omega)
]
but now (Z^{\rm eff}) is a property of the full coupled dynamic system, not a property of the superfluid-only equilibrium branch.

---

### 5.9 Summary: the coupled “master equation bundle”

The full dynamic model is defined by the coupled system:

1. **Matter PDE (gauged GNLS in 4D):** (\psi(\mathbf X,t)).
2. **Gauge PDE (Maxwell + localization + gauge-fixing in 4+1D):** (A_M(\mathbf X,t)).
3. **Geometry ODEs (dynamic throat):** (a(t),L(t)) driven by total energy gradients.
4. **Open-system terms (drive/dissipation/boundaries):** required for intake, cavity pumping, and realistic back-pressure.
5. **Projection + measurement map:** bulk (\to) brane observables and response.

---

## 6. Key Physical Consistency Requirements

This section lists the non-negotiable constraints the unified system must satisfy to be physically meaningful. These are not “nice-to-have checks”; each one corresponds to a failure mode that can quietly invalidate conclusions if it is ignored.

---

### 6.1 Charge neutrality of the vacuum (no uniformly charged condensate)

**Requirement:** The far-field vacuum state must be electrically neutral on the brane and in the bulk (up to controlled, localized defect/excitation charge).

**Why:** If the condensate (\psi) directly sources charge density as (J^0=q|\psi|^2), then a homogeneous vacuum becomes a uniformly charged medium, producing unphysical large-scale fields.

**Operational enforcement options (must be frozen explicitly):**

* **Background subtraction:** (J^0=q(|\psi|^2-\rho_0)) (and associated consistent definition of spatial current).
* **Two-field matter content:** neutral (\psi) + charged excitation field (\chi).
* **Defect-localized charge:** source proportional to deviation from background *within* a gated region, e.g. (J^0 \propto G(\mathbf X)(|\psi|^2-\rho_0)), with careful conservation bookkeeping.

**Check:** In the far-field vacuum initial condition, the solved Maxwell sector must produce vanishing large-scale brane fields and a stable neutral state (modulo gauge).

---

### 6.2 Brane-observed EM scaling vs 4D Gauss law (localization is mandatory)

**Requirement:** The brane observer must see an effective EM behavior consistent with the intended phenomenology (often 3D Coulomb-like scaling), even though the underlying gauge field propagates in 4 spatial dimensions.

**Why:** In 4 spatial dimensions, free Maxwell fields generally imply different far-field scaling than 3D intuition. If we do nothing, the brane will not automatically reproduce familiar EM behavior.

**Operational enforcement:** Introduce a localization mechanism in the EM sector:

* position-dependent prefactor (Z(\mathbf X)) in (F^2), and/or
* Proca/Higgs-like mass (m_\gamma^2(\mathbf X)A^2), and/or
* brane-localized terms that trap a gauge zero-mode.

**Checks:**

* Compute effective brane field profiles from a localized source and verify intended scaling in the brane-projected observables.
* Verify gauge-mode confinement length in (w) is finite and stable against throat dynamics.

---

### 6.3 Back-pressure consistency in a bulk-expansion receiving region

**Requirement:** Once flux passes through the throat into the bulk, the receiving region must be modeled so that density buildup and back-pressure are physically determined by:

* EOS pressure response,
* 4D spreading geometry,
* and any confining potentials.

**Why:** If the bulk region is not represented (or is represented with reflecting boundaries), “back-pressure” becomes a numerical artifact and can fake stabilization/destabilization.

**Operational enforcement:**

* A bulk domain where confinement permits expansion (at least in the throat-connected region).
* Absorbing/outflow boundary conditions (or sufficiently large domains) for both matter and EM fields.
* Explicit diagnostic definition of back-pressure (e.g., brane-projected enthalpy/pressure near the mouth vs far-field reference).

**Checks:**

* Repeat runs with enlarged domains / stronger absorbers and confirm brane-observable response converges.
* Track mass/energy flux into bulk and verify conservation/balance identities.

---

### 6.4 Aspect ratio and endcap physics: why tube (E_{\rm geom}) alone cannot decide (L/a)

**Requirement:** If the physics predicts a preferred aspect ratio (e.g. a stable (L/a)), the model must contain at least one mechanism that can *select* it.

**Why:** The baseline geometry energy (E_{\rm geom}=P_{\rm vac}V+\sigma A) for a tube typically gives (\partial_L E_{\rm geom}) ~ constant, so geometry alone does not create an (L)-dependent restoring force or stable (L).

**How selection can occur (must be explicit):**

* Field-mediated selection: EM cavity modes exert (L)-dependent stress; flow/back-pressure gives (L)-dependent forces.
* Geometry upgrade: curvature/rounding/flare energy terms or shape-field (R(w,t)).
* DOF reduction: (L) is not independent but slaved to other physics (e.g., mode quantization constraint).

**Checks:**

* Demonstrate a stable equilibrium or bounded dynamics for (L(t)) under the chosen mechanism.
* Show that equilibrium persists under parameter variation and is not a boundary artifact.

---

### 6.5 Dimensional reduction sanity checks (recover correct 3D limits)

**Requirement:** In the “throat closed / strong confinement” regime, the system must reduce to a consistent effective 3D model on the brane.

**Why:** If the far-field does not reduce correctly, then brane observables are not meaningful and the throat interpretation is not anchored.

**Operational enforcement:**

* Strong (w)-trap far from the throat should suppress (w)-excitations.
* The brane projection (\mathcal P) should show small tail-mass outside the brane-localized mode.
* EM localization should produce brane-like behavior away from the defect.

**Checks:**

* Quantify (w)-mode leakage (tail diagnostics) in far field.
* Recover known 3D wave/field relations under brane confinement limits.

---

### 6.6 Symmetry and conservation checks (the “trust but verify” set)

**Requirement:** The unified system must satisfy the identities implied by its symmetries, and any symmetry-breaking terms must be accounted for explicitly.

**Key checks:**

* **Gauge symmetry / continuity:** verify (\partial_M J^M=0) (or track deviations if reservoir terms break it).
* **Energy balance:** verify (dH_{\rm tot}/dt) matches injected power minus dissipation minus boundary flux.
* **Momentum/stress balance:** verify consistent force transfer between fields and geometry via the derived stress tensors / energy gradients.
* **Parity and throat symmetry:** for symmetric configurations, verify left/right (w)-splits are consistent (e.g., equal and opposite contributions where expected).

These checks should be part of both the Mathematica derivation pipeline (symbolic identities) and the numerical regression suite.

---

### 6.7 Projection validity: “we only see a slice” must be controlled

**Requirement:** Brane observables must be stable with respect to the choice of projection operator, windowing, and resolution.

**Why:** If brane signals depend sensitively on projection truncation, then what we call “observable” is not robust.

**Operational enforcement:**

* Use normalized truncated weight (W(w)) and track tail mass.
* Compare restriction vs projection in regimes where both should agree.
* Ensure port integrals converge with quadrature resolution.

**Checks:**

* Convergence of (\rho_{\rm brane}), port signals, and extracted (Z^{\rm eff}) under projection-window changes.
* Tail-mass remains below threshold in the “brane-like” regime.

---

### 6.8 Numerical well-posedness requirements (avoid fake physics)

**Requirement:** The system must be posed and solved in a way that does not inject spurious stabilization.

**Practically:**

* Gauge fixing must be consistent and stable (especially for time evolution).
* Absorbers must be designed so they remove outgoing radiation/flow without reflection.
* Geometry integration must not be artificially damped/stiffened relative to physical expectations unless intentionally modeled.

**Checks:**

* Repeat key runs with alternative gauge parameters (\xi), solver resolutions, and absorber strengths; verify physical conclusions persist.
* Validate against limiting cases and analytic small-perturbation expectations where available.

---

### 6.9 Summary: the “must-pass” gate before interpreting results

A coupled dynamic 4D throat model is only interpretable if it passes these consistency requirements:

1. Neutral vacuum far field (no uniform charge).
2. EM localization sufficient for intended brane scaling.
3. Bulk receiving region and back-pressure modeled physically (not numerically).
4. Aspect ratio selection mechanism explicitly present (fields, shape energy, or DOF constraints).
5. Correct 3D limit recovery when throat is closed.
6. Conservation and symmetry identities verified (symbolic + numeric).
7. Projection/port observables robust.
8. Numerics well-posed and regression-tested.

---

## 7. Reduction Strategies and “What Pops Out”

This section describes how we plan to extract a compact, interpretable set of equations and derived parameters from the full coupled 4D system. The goal is not to “simplify for convenience,” but to (i) identify controlled limits where brane-effective dynamics emerge, and (ii) let the unified variational structure generate new constraints (the analogs of earlier “(\beta=3)” / “(n=5)” moments).

---

### 7.1 The reduction philosophy: derive, don’t assume

We treat reductions as **derived results** from the master action, not as prior assumptions. That means:

* start with the coupled system in 4D ((\psi, A_M, a(t), L(t))),
* identify relevant small parameters or separation of scales (strong brane confinement, slow geometry, localized gauge zero-mode, weak coupling, etc.),
* perform expansions or mode truncations in a way that preserves gauge invariance and conservation laws,
* and verify that the reduced equations reproduce brane observables produced by the full system in the same regime.

Reductions that fail these checks are not used.

---

### 7.2 Mode decomposition in (w): brane-localized zero modes and KK tower

The primary structured reduction is to separate brane-parallel and bulk directions via a mode expansion in (w).

#### 7.2.1 Matter field expansion (if brane-like regime exists)

In the strong-confinement regime, we expect (\psi) to be dominated by a brane-localized mode:
[
\psi(x,y,z,w,t) \approx \sum_{n} \psi_n(x,y,z,t)\,u_n(w; a,L)
]
where (u_n) are eigenmodes of the far-field (w)-confining operator (modified near the throat).

A “brane effective” truncation keeps (n=0) (and possibly the first few excited modes) and integrates over (w). This yields effective 3+1D dynamics with couplings given by overlap integrals like:
[
g_{\rm eff}\sim \int dw\,|u_0|^4,\qquad
V_{\rm eff}(x,y,z)\sim \int dw\,|u_0|^2 V_{\rm conf}.
]

#### 7.2.2 Gauge field expansion (localization determines everything)

Similarly, for the gauge field:
[
A_M(x,y,z,w,t) \approx \sum_{n} a_M^{(n)}(x,y,z,t)\,v_n(w; a,L)
]
where the (v_n) are determined by the localization mechanism ((Z(w)), (m_\gamma(w)), or brane terms).

The existence of a localized (n=0) gauge mode and its profile (v_0(w)) control:

* whether the brane sees “3D-like EM,”
* the effective coupling (q_{\rm eff}),
* and the leakage/radiative loss into the bulk.

---

### 7.3 Deriving effective constants (“new (\beta/n)-like outputs”)

Once mode profiles are determined, effective parameters on the brane are no longer arbitrary—they are integrals.

Examples of quantities that should *emerge*:

1. **Effective gauge coupling on the brane**
   [
  q_{\rm eff} \sim q \int dw\,W(w)\,|v_0(w)|\,|u_0(w)|^2
   ]
   (schematic; exact form depends on normalization conventions).

2. **Effective EM stiffness / permittivity-like constants**
   From the reduced (F^2) term:
   [
   \frac{1}{g_{\rm EM,eff}^2} \sim \int dw;Z(w),|v_0(w)|^2.
   ]

3. **Localization length(s)**
   The decay scale of (v_0(w)) and/or (u_0(w)) defines:
   [
   \ell_A,\ \ell_\psi
   ]
   which dictate coupling strength, leakage, and whether brane sees bulk effects.

4. **Geometry-coupled coefficients**
   Because (v_n(w)) and (u_n(w)) depend on ((a,L)), the reduction produces explicit functions:
   [
   q_{\rm eff}(a,L),\quad g_{\rm eff}(a,L),\quad \ell_A(a,L),
   ]
   which become inputs to the geometry force equations and can generate “unexpected” equilibrium or stability conditions.

These are exactly the sort of values we expect to “pop out” once the unified coupling is treated correctly.

---

### 7.4 Linearization about a dynamic background: stability and resonance conditions

To study cavity support and throat stability, we linearize around a background trajectory:
[
\psi=\psi_0+\delta\psi,\quad
A_M=A_{M0}+\delta A_M,\quad
a=a_0+\delta a,\quad
L=L_0+\delta L
]
where ((\psi_0,A_{M0},a_0,L_0)) may be:

* a steady state under drive/dissipation, or
* a slowly varying solution (adiabatic geometry).

Linearization yields:

* coupled wave equations for perturbations,
* a mode spectrum that depends on ((a_0,L_0)),
* and stability conditions of the form “net restoring vs net destabilizing” for (\delta a,\delta L).

This is where the “photon cavity holds throat open” claim becomes a calculable statement: the EM spectrum and stored energy determine stress, which determines geometry response.

---

### 7.5 Regimes and dominance: when each mechanism matters

We expect distinct regimes:

#### (i) Wave-dominated support (cavity pressure dominant)

* Strong localized EM modes,
* significant stored field energy,
* geometry stabilized primarily by Maxwell stress contributions.

#### (ii) Flow-dominated support (intake/back-pressure dominant)

* Strong matter flux through the throat,
* significant momentum/pressure support from the superfluid,
* geometry stabilized by hydrodynamic stress rather than EM.

#### (iii) Mixed regime (most realistic)

* Cavity and flow reinforce or compete,
* stability is set by the coupled energy balance and geometry forces.

Each regime should produce distinctive brane observables (phase relationships, impedance profiles, resonance shifts).

---

### 7.6 Extracting response functions and impedance from the full system

Once we have a steady or periodic operating point under drive, we compute brane observables (u_i(t), j_i(t)) from ports and define frequency-domain response:
[
\mathbf j(\omega)=Z^{\rm eff}(\omega)\,\mathbf u(\omega).
]

In the reduced description, (Z^{\rm eff}) is expected to depend on:

* bulk leakage (through localization length (\ell_A) and matter leakage),
* cavity resonance structure (EM mode spectrum),
* and back-pressure dynamics in the receiving region.

This provides a systematic way to connect the unified PDE system to testable brane signatures.

---

### 7.7 What we are explicitly looking to “pop out”

We expect the unified system to generate constraints and parameter relations that were not visible in separated-module reasoning. The primary targets are:

1. **A self-consistent aspect ratio selector** (L/a) emerging from EM mode structure + geometry energy + flow/back-pressure, rather than from geometry energy alone.
2. **A derived EM localization requirement** (e.g., a specific form or scale for (m_\gamma(w)) or (Z(w))) demanded by matching brane-observed scaling.
3. **A consistent energy closure** linking throat energy, density, and effective mass-energy (formalizing the “density drop → energy shed” intuition).
4. **Coupling renormalizations** (q_{\rm eff}(a,L)) and EOS-dependent stiffnesses that alter stability in non-obvious ways.
5. **New dimensionless control parameters** that partition the regime diagram (analogous to previous fixed-exponent results).

---

### 7.8 Summary: reduction as a discovery engine

The reduction pipeline is not an afterthought; it is the mechanism by which the unified 4D theory becomes:

* interpretable,
* testable on the brane,
* and capable of producing new derived relations.

---

## 8. Mathematica Program: How We Derive and Verify the System

This section lays out the Mathematica-first workflow for turning the unified action into (i) a verified equation bundle, (ii) compact conservation/force identities, and (iii) reduced “brane-effective” equations with derived constants. The intent is to use Mathematica not as a calculator, but as the **authoritative derivation engine** for the coupled system.

---

### 8.1 Symbolic representation of the action and automatic Euler–Lagrange derivations

**Goal:** Encode the full master action
[
S[\psi,A_M,a,L]
]
in Mathematica in a way that allows automated variation with respect to:

* (\psi^*) (matter PDE),
* (A_M) (Maxwell PDEs),
* (a(t),L(t)) (geometry force laws / wall dynamics).

**Implementation plan (WL structure):**

1. **Canonical definitions layer (frozen):**
   Reuse the existing operational definition stack (the definitions that already drive `master_checks.wl` and `throat_spec_derivation.wl`). This ensures:

   * the gate (G(\mathbf X;a,L)),
   * the smooth-step primitives,
   * the confinement potential (V_{\rm conf}(\mathbf X;a,L)),
   * and sign conventions
     are single-sourced.

2. **New “unified action” layer:**
   Create a dedicated WL file/notebook (e.g. `unified_action.wl`) that defines:

   * covariant derivatives (D_t, D_i),
   * (F_{MN}),
   * localization term choice ((Z) and/or (m_\gamma)),
   * matter EOS potential (U(\rho)),
   * gauge fixing term,
   * and (optionally) drive/dissipation placeholders as symbolic source terms.

3. **Variation engine:**
   Implement variation by:

   * explicit Euler–Lagrange operators for fields with derivatives (functional derivatives),
   * or a controlled “integration-by-parts then coefficient match” routine,
   * always under assumptions that keep expressions manageable (real parameters, smoothness, boundary decay or periodic BCs).

**Deliverable:** Mathematica outputs the coupled PDE set in a normalized “paper form” and in a machine-usable form (for code generation / checks).

---

### 8.2 Automated conservation law derivations

**Goal:** Produce the compact identities that “encapsulate” the model and often reveal hidden constraints—especially when geometry is dynamic and when EM localization is present.

We will derive, symbolically:

1. **Gauge current and continuity**
   From gauge symmetry / Noether:
   [
   \partial_M J^M = 0
   ]
   (or the exact modified statement if we include explicit source/sink terms). Mathematica should produce (J^M) in terms of (\psi) and (A_M) automatically from the action (not hand-written), so the current used in Maxwell is guaranteed consistent.

2. **Energy balance identity**
   Derive an explicit statement of the form:
   [
   \frac{d}{dt}H_{\rm tot}
   =
   P_{\rm drive} - P_{\rm diss} - \Phi_{\rm boundary}
   + \text{(work exchange terms between sectors)}.
   ]
This is where the “density change ↔ energy shedding ↔ geometry response” story becomes an enforceable identity.

3. **Momentum / stress transfer**
   Derive a stress-energy or momentum-flux statement (appropriate to the hybrid nonrelativistic matter + Maxwell sector) that cleanly shows:

* Maxwell stress contributions,
* matter stress contributions,
* and how geometry forces arise from total-energy gradients.

**Deliverable:** A set of symbolically verified conservation identities that become regression tests for numerics.

---

### 8.3 Consistency and limit tests (symbolic + numeric cross-checks)

**Goal:** Build a suite of Mathematica checks that fail loudly if we break the model.

Core tests:

1. **Closed-throat / strong confinement limit**
   Show that when the gate is “off” and brane confinement is strong, the system reduces to:

* brane-localized matter dynamics (no bulk leakage),
* localized gauge behavior consistent with intended brane EM.

2. **Neutral vacuum limit**
   For the chosen neutrality strategy (background subtraction or two-field):

* verify far-field (J^0) integrates to zero (or vanishes pointwise),
* verify Maxwell equations admit a trivial/no-field vacuum.

3. **Gauge-fixing consistency**
   Confirm the gauge-fixed Maxwell equations remain compatible with continuity (no hidden inconsistencies).

4. **Geometry-force sanity**
   Mathematica should reproduce the key force identities:
   [
   F_a^{(\psi)}=-\int d^4X\,\rho(\mathbf X,t)\,\partial_a V_{\rm conf}(\mathbf X;a,L),\quad
   F_L^{(\psi)}=-\int d^4X\,\rho(\mathbf X,t)\,\partial_L V_{\rm conf}(\mathbf X;a,L)
   ]
   and show the additional force contributions from localization terms when (Z) or (m_\gamma) depend on ((a,L)).

5. **Symmetry/parity checks**
   In symmetric configurations, verify expected parity relations (e.g. equal split of contributions for (w>0) and (w<0) where applicable).

**Deliverable:** A Mathematica “master checks” notebook for the unified system, analogous in spirit to the existing check files, but upgraded to the coupled gauge–matter–geometry model.

---

### 8.4 Parameter extraction and “new constants” pipeline

**Goal:** Use Mathematica to extract the reduced, brane-effective parameters that “pop out” when the unified system is reduced properly.

Planned pipeline:

1. **Mode decomposition in (w)**
   Symbolically define the (w)-direction Sturm–Liouville problems implied by:

* matter confinement / gate structure (for (u_n(w))),
* EM localization structure (for (v_n(w))).

2. **Overlap integrals → effective couplings**
   Compute integrals that produce:

* effective brane coupling (q_{\rm eff}),
* effective EM stiffness (or permittivity-like constants),
* localization length scales (\ell_A,\ell_\psi),
* geometry-dependent renormalizations (q_{\rm eff}(a,L)), etc.

3. **Geometry selection conditions**
   Derive conditions under which the combined stresses + geometry energy yield stable ((a,L)) behavior. This is where a preferred aspect ratio can emerge without being put in by hand.

4. **Dimensionless regime parameters**
   Construct a minimal set of dimensionless control parameters that govern:

* wave-dominated vs flow-dominated support,
* leakage-dominated vs trapped regimes,
* stable vs unstable geometry.

**Deliverable:** A Mathematica-generated “Derived Parameters Table” that becomes a new frozen reference (the analog of earlier exponent-fixing milestones).

---

### 8.5 Cross-validation with numerical solvers and regression against prior papers

**Goal:** Ensure the symbolic theory and the numerical solvers stay locked together.

1. **Equation export and code alignment**

* Export frozen parameter blocks (JSON) for Python runs (same philosophy as current config export).
* Export derived symbolic expressions (e.g. analytic derivatives, reduced coefficients) for direct numerical evaluation.

2. **Regression tests**

* Check that previously frozen operational definitions (ports, projection, surface measures, sign conventions) remain consistent with the new unified dynamics.
* Re-run known-limits tests (e.g., 1PN hybrid behaviors, brane-only behavior) as regression anchors.

3. **Identity checks in simulation**

* Numerics must report:

  * continuity residuals,
  * energy balance residuals,
  * gauge-condition residuals,
  * and geometry force residuals.
    The exact symbolic forms of these diagnostics come from Mathematica.

**Deliverable:** A stable “symbolic-to-numeric handshake” where Mathematica is the authority on the equations and Python/GPU code is the authority on time evolution, with automated checks ensuring neither silently drifts.

---

### 8.6 Summary: Mathematica as the discovery engine

Mathematica’s role is to:

* generate the coupled equations correctly,
* derive the conservation laws that govern energy and force exchange,
* expose what must be added (localization, neutrality, geometry upgrades) for consistency,
* and produce the reduced effective constants and regime parameters that should “pop out” only once the full system is combined properly.

---

## 9. Computational Plan

This section describes how we will move from the unified equation bundle to executable simulations and diagnostics, while keeping the system honest (no fake stabilization from numerics, no hidden inconsistency with the Mathematica-derived identities). The plan is staged: start with a minimal coupled prototype that is stable and checkable, then add complexity only when each module passes the consistency gates in Section 6.

---

### 9.1 Minimal coupled prototype (first runnable “full dynamics” system)

**Goal:** A simulation that evolves (\psi(\mathbf X,t)), (A_M(\mathbf X,t)), and (a(t),L(t)) together in 4D, even if the first version is coarse.

**Prototype feature set:**

1. **Matter PDE:** gauged GNLS in 4D with the frozen EOS and baseline (V_{\rm conf}^{(F1)}(\mathbf X;a,L)).
2. **Gauge PDE:** 4+1D Maxwell with gauge fixing (Lorenz gauge family) and **one** chosen localization mechanism (start with the simplest variational term).
3. **Geometry:** dynamic wall-law ODEs for (a(t),L(t)) with tunable inertia/damping, driven by energy gradients computed from fields.
4. **Open-system mechanism:** minimal drive/dissipation sufficient to sustain a cavity mode and allow intake/outflow without boundary artifacts.

**Why this is the correct “minimal” set:** It already contains all three coupled sectors and the core non-static physics, which is what the current equilibrium-only scans cannot capture.

---

### 9.2 Numerical scheme sketches (what we actually implement)

#### 9.2.1 Matter solver (gauged GNLS)

We maintain the philosophy used in prior work: a stable operator-split approach is acceptable, but it must respect the coupling.

Practical approach:

* Use a **split-step** or **Strang splitting** method:

  * kinetic step with covariant derivatives (requires careful handling of (\mathbf A)),
  * potential + nonlinearity step,
  * optional damping/reservoir step.
* Spatial discretization options:

  * pseudo-spectral FFT in all 4 dimensions (fast but periodic BCs),
  * or finite-difference in (w) with spectral in ((x,y,z)) (hybrid).

**Key implementation constraint:** Covariant derivatives must be implemented in a gauge-consistent way. Two reasonable routes:

* (a) compute ((\nabla - i q\mathbf A/\hbar)) explicitly on the grid (FD + staggered grid if needed),
* (b) use a link-variable (“Peierls phase”) discretization to preserve gauge invariance discretely.

We start with the simplest consistent approach and use gauge-residual diagnostics to decide if we need the link-variable upgrade.

#### 9.2.2 Gauge solver (4+1D Maxwell)

We need an explicit time evolution scheme for (A_M) or for field strengths.

Two viable approaches:

* **Potential-based wave equation form** in Lorenz gauge:
  [
  \Box A^N + \text{(localization terms)} = J^N + J_{\rm ext}^N
  ]
  discretized with a stable explicit method (leapfrog / RK) and absorbing layers.
* **Field-strength evolution** (FDTD-like), with gauge handled implicitly or separately.

We prefer the potential-based Lorenz-gauge wave-equation form initially because:

* it is derivable cleanly in Mathematica,
* it gives a direct PDE for each component (A_N),
* it is easier to add Proca/localization terms.

#### 9.2.3 Geometry ODE integrator

Evolve
[
M_a\ddot a + C_a\dot a = F_a,\qquad
M_L\ddot L + C_L\dot L = F_L
]
with forces computed from the fields (matter + EM + localization + geometry energy).

Use a robust ODE integrator (e.g., symplectic-ish if conservative, otherwise standard RK) and ensure time-stepping is synchronized with PDE steps.

---

### 9.3 Coupling strategy: how the sectors talk each timestep

At each timestep (t\to t+\Delta t):

1. Given ((\psi,A_M,a,L)), compute:

   * (V_{\rm conf}(\cdot;a,L)),
   * localization profiles (Z(\cdot;a,L)) or (m_\gamma(\cdot;a,L)),
   * matter current (J^N[\psi,A]) (with neutrality strategy enforced).

2. Update (A_M) using Maxwell + localization + gauge fixing + external drive.

3. Update (\psi) using gauged GNLS with the updated (A).

4. Compute geometry forces (F_a,F_L) from total energy gradients (or force identities).

5. Update (a,L) using the wall ODEs.

6. Emit diagnostics:

   * continuity residual,
   * gauge condition residual,
   * energy balance residual,
   * boundary fluxes,
   * brane observables and port signals.

This explicit “operator splitting across sectors” is acceptable if it passes the conserved-identity regression tests from Mathematica.

---

### 9.4 Boundary conditions and absorbers (preventing fake back-pressure)

Because the model depends strongly on outflow and bulk expansion, boundary handling is critical.

**Matter sector:**

* Absorbing layers in the outer region (complex absorbing potential or sponge damping).
* Large enough domain in (w) and in brane directions to avoid early reflection.

**Gauge sector:**

* Absorbing boundary layers for outgoing radiation (PML-like or sponge).
* Localization profiles should not unintentionally act as reflective walls unless that is the intended cavity.

**Checks:**

* Domain-size convergence tests: key observables must not change materially when domain is enlarged.
* Absorber-strength convergence tests: results must stabilize as absorbers are strengthened without distorting near-throat physics.

---

### 9.5 Diagnostics: what we measure to decide “physics vs artifact”

Diagnostics must be treated as first-class outputs, not afterthoughts.

**Core invariants and residuals**

* Charge/continuity residual: (|\partial_M J^M|).
* Gauge residual: (|\partial_M A^M|) (or chosen gauge condition).
* Energy balance residual: compare (dH_{\rm tot}/dt) to injected/dissipated/boundary flux.
* Geometry work terms: confirm (F_a\dot a + F_L\dot L) matches exchange terms.

**Physical diagnostics**

* Stored EM energy in cavity region vs time.
* Matter density/pressure profiles near throat and in receiving region.
* Net flux through mouth (port flux) and its phase relative to drive.
* “Back-pressure” proxy: brane-projected enthalpy/pressure near the mouth vs reference.
* Stability measures: boundedness of (a(t),L(t)), response to perturbations.

**Brane measurement outputs**

* Port time series (u_i(t)), (j_i(t)).
* Extracted (Z^{\rm eff}(\omega)) under periodic drive.
* Resonance tracking: shift of peak frequencies as geometry evolves.

---

### 9.6 Stability and falsification tests (we want ways to make it fail)

We explicitly plan tests designed to break the model:

1. **No-drive collapse test:** with cavity drive off and intake disabled, the throat should collapse unless geometry energy alone supports it (which we do not expect in general).
2. **Drive-only test:** with cavity pumping but no intake, verify whether EM stress can maintain a throat.
3. **Intake-only test:** with intake forcing but minimal cavity field, test whether flow stress alone can maintain a throat.
4. **Leakage test:** weaken localization and verify that brane-observed EM behavior and stabilization degrade in the predicted way.
5. **Neutrality test:** perturb background density and verify no spurious global charge field grows.
6. **Resolution / domain scaling:** ensure stability claims persist under refinement.

A mechanism that only works in a fragile numerical setting is not accepted.

---

### 9.7 Implementation milestones (what gets built first)

**Milestone A — Symbolic closure**

* Mathematica generates the exact PDE bundle and the diagnostic identities.
* A frozen parameter block is exported.

**Milestone B — Minimal coupled solver**

* 4D evolution of (\psi) coupled to evolved (A_M) (even if coarse).
* Basic absorbers and external drive.
* Diagnostic outputs pass basic sanity checks.

**Milestone C — Dynamic geometry**

* Add (a(t),L(t)) evolution with forces from fields.
* Verify energy/work consistency.

**Milestone D — Brane observables**

* Implement projection map (\mathcal P), measurement surface (\Gamma), ports, and response extraction.
* Produce first (Z^{\rm eff}(\omega)) in a coupled run.

**Milestone E — Derived constants**

* Run mode-reduction pipeline and extract effective parameters.
* Compare reduced predictions to full simulation outputs.

---

### 9.8 Summary: computational plan as a controlled ascent

We will build the unified simulation in a way that:

* keeps the variational derivation as the authority,
* treats open boundaries and absorbers as critical physics infrastructure,
* and emits diagnostics that prevent us from mistaking artifacts for stabilization.

---

## 10. Roadmap: Where We Are, Where We’re Going, and End-State Targets

This section anchors the document in the current project reality (what is already frozen and working), then lays out a practical path to the end-state: a unified, fully dynamic 4D throat model with gauge fields, geometry dynamics, and brane-observable predictions.

---

### 10.1 Current status (what exists today, what it proves, and what it doesn’t)

**Operational definitions are frozen and executable (Paper 7 Steps 0–5).**
We have a stable operational layer that defines:

* geometry DOFs ((a,L)) and sign conventions,
* a brane projection window and weight (W(w)) with tail-mass diagnostic,
* mouth surface (\Gamma), surface normal/measure, and port basis,
* drive protocol and response extraction pipeline (lock-in, conditioning),
* unit/consistency checks connecting WL definitions to executables.

**We have a working 4D GNLS solver on uniform 4D grids (GPU).**
This includes:

* split-step spectral evolution (for the matter-only model),
* imaginary-time relaxation to find ground states at fixed norm (diagnostic equilibrium finder),
* Family-1 confinement potential implementation,
* numerical derivatives (FD) and now analytic derivatives for (\partial_a V), (\partial_L V),
* equilibrium scan harness and diagnostic scripts (e.g., (I_a), (I_L), split checks).

**What the current equilibrium scans show.**

* Radial force balance can be close in the neighborhood (a\sim 0.8).
* Length/endcap balance is not close under the baseline geometry-energy + confinement choice: residuals remain dominated by (R_L) with the best scan point around ((a,L)\approx(0.8167,2.75)).
* Increasing grid size changes (I_L) and residuals but does not obviously converge to an equilibrium root in the scanned ranges.

**What the current equilibrium scans do *not* prove.**

* They do not test the core physical stabilization hypotheses (cavity stress + intake support), because:

  * there is no evolved gauge field sector,
  * there is no explicit open-system intake/back-pressure model,
  * ((a,L)) are not dynamically evolved under coupled stresses,
  * and “static equilibrium at fixed norm” is not the intended full physical closure.

So the present state is excellent for “operational plumbing” and “matter-only diagnostics,” but it is not yet the coupled physical system.

---

### 10.2 Immediate next milestones (the shortest path to “fully dynamic + gauge”)

**Milestone 1 — Freeze the neutrality strategy and localization strategy.**
Before we can write final coupled PDEs, we must explicitly choose:

* charge neutrality implementation (background subtraction vs two-field),
* gauge localization mechanism (e.g., (m_\gamma(w)) or (Z(w))).

This is the single most important “spec freeze” remaining, because it determines whether the gauge sector can reproduce intended brane-observed behavior.

**Milestone 2 — Mathematica derivation of the unified equation bundle.**
Implement the master action in WL and derive:

* gauged GNLS,
* 4+1D Maxwell + localization + gauge fixing,
* geometry force laws including EM/localization contributions,
* conservation laws (continuity, energy balance),
* and the diagnostic residual expressions used in simulation.

**Milestone 3 — Minimal coupled solver (coarse but honest).**
Build the first simulation that evolves:
[
(\psi(\mathbf X,t),\ A_M(\mathbf X,t),\ a(t),\ L(t))
]
with:

* absorbers/outflow boundaries,
* a basic cavity drive,
* and diagnostics reporting conservation/gauge residuals.

This solver is the first system that can answer the real question: “can a throat be maintained dynamically by wave + flow stresses?”

---

### 10.3 Medium-term milestones (where “new constants” and selection rules emerge)

**Milestone 4 — Dynamic stability mapping.**
Run parameter sweeps over:

* drive strength / frequency,
* intake forcing / reservoir coupling,
* localization length scales,
* geometry-energy parameters ((P_{\rm vac},\sigma), etc.),
  and map regimes:
* collapse,
* oscillatory breathing,
* metastable opening,
* sustained steady operation.

**Milestone 5 — Mode reduction and derived effective parameters.**
Perform (w)-mode decompositions for both matter and gauge sectors and extract:

* (q_{\rm eff}), (g_{\rm EM,eff}),
* localization lengths (\ell_A,\ell_\psi),
* geometry-dependent renormalizations (q_{\rm eff}(a,L)),
* and any intrinsic aspect ratio selection rule that emerges.

This is explicitly where we expect “beta/n-like” values to appear, but now as reduction outputs rather than assumptions.

**Milestone 6 — Brane-observable response and impedance predictions.**
Using the frozen port pipeline, compute:

* port time series (u_i(t),j_i(t)),
* resonance curves and resonance tracking,
* (Z^{\rm eff}(\omega)) under drive,
* and compare reduced-theory predictions to full 4D simulation.

---

### 10.4 End-state targets (what “done” looks like for this short paper’s purpose)

By the time this roadmap is executed, we want:

1. **A single unified coupled model**, written cleanly as:

* one master action (or equivalent),
* one PDE/ODE bundle,
* one brane projection/measurement map.

2. **A falsifiable stabilization claim**, expressed as:

* explicit conditions on drive, flow, and localization under which the throat persists,
* and explicit conditions under which it must collapse.

3. **A set of derived effective parameters**, with formulas and numerical values in representative regimes:

* effective brane coupling(s),
* localization lengths,
* geometry selection conditions (e.g., an emergent (L/a) relation),
* regime-control dimensionless numbers.

4. **A Mathematica-backed derivation pipeline**, producing:

* equations,
* conservation identities,
* reduction integrals,
* regression diagnostics.

5. **A minimal coupled simulation that passes the “must-pass” consistency gates** in Section 6:

* neutrality,
* localization,
* back-pressure realism,
* conservation laws,
* projection robustness.

---

### 10.5 Known open questions (kept visible, not buried)

These are the items we explicitly expect to resolve (or falsify) during the next phase:

* **Neutrality:** which matter content strategy is physically and mathematically cleanest?
* **Localization:** what minimal EM localization mechanism yields correct brane behavior without over-constraining the bulk?
* **Geometry selection:** what term(s) actually produce a stable (L) (and thus aspect ratio) in a way that survives refinement?
* **Open-system closure:** what is the correct “rest-energy + density” bookkeeping (how exactly does energy shed when density changes)?
* **Coupling strength scaling:** does the unified model generate realistic separations of scale between throat energetics and brane-observed EM effects?
* **Numerical robustness:** do stabilization results survive domain enlargement, absorber strengthening, and resolution refinement?

---

### 10.6 Summary: how this connects to where we want to end up

We are past the “definition chaos” phase: the operational layer is frozen and we can measure what we compute. The missing step is not more equilibrium scanning; it is assembling the **full dynamic coupled system** (superfluid + 4D gauge fields + dynamic geometry + open-system closure) and letting the unified structure tell us what must be true.

The end goal is a model that is:

* internally consistent by derivation,
* externally testable via brane observables,
* and capable of generating new derived constraints that either validate the throat hypothesis in a restricted regime or decisively rule it out.

---

## Appendix A. Notation, Units, and Index Conventions

### A.1 Coordinates and domains

* **Time:** (t\in\mathbb R)
* **Bulk spatial coordinates:**
  [
  \mathbf X=(x,y,z,w)\in\mathbb R^4
  ]
* **Brane:** operationally located near (w=0); brane observables are obtained by restriction (w=0) and/or projection with a weight (W(w)).

**3D radial coordinate on the brane:**
[
R_3 \equiv \sqrt{x^2+y^2+z^2}.
]

---

### A.2 Indices and summation conventions

We use Einstein summation unless otherwise stated.

* **Brane spatial indices:** (a,b,c \in {x,y,z})
* **Bulk spatial indices:** (i,j,k \in {x,y,z,w})
* **Spacetime indices:** (M,N \in {0,x,y,z,w}) where (0) denotes time.

Spatial derivatives:
[
\partial_i \equiv \frac{\partial}{\partial X^i},\qquad
\nabla_4 \equiv (\partial_x,\partial_y,\partial_z,\partial_w).
]

---

### A.3 Fields and derived quantities

**Matter / superfluid field**

* (\psi(\mathbf X,t)\in\mathbb C): order parameter / condensate field
* (\rho(\mathbf X,t) \equiv |\psi(\mathbf X,t)|^2): density

**Gauge field**

* (A_M(\mathbf X,t)\in\mathbb R): 4+1D gauge potential
* Field strength:
  [
  F_{MN}=\partial_M A_N-\partial_N A_M.
  ]

**Geometry degrees of freedom**

* (a(t)>0): effective throat radius
* (L(t)>0): effective throat length

---

### A.4 Covariant derivatives (C1 minimal coupling)

Charge parameter (q), effective mass parameter (m), and (\hbar) as usual.

Define:
[
D_t \equiv \partial_t + \frac{i q}{\hbar}A_0,\qquad
D_i \equiv \partial_i - \frac{i q}{\hbar}A_i,\quad i\in{x,y,z,w}.
]

---

### A.5 Currents and continuity

A standard gauge current derived from the matter sector (exact form depends on the neutrality strategy for (J^0)):

Mass/number current:
[
J_{\rm mass}^i = \frac{\hbar}{m}\Im\!\left(\psi^* D_i\psi\right)
\qquad i\in{x,y,z,w}.
]

Charge current sourcing Maxwell:
[
J_{\rm ch}^i = q\,J_{\rm mass}^i
]
or (J_{\rm ch}^i) derived from (\chi) if the two-field option (N2) is used.

Continuity (in the conservative/gauge-consistent limit):
[
\partial_M J_{\rm ch}^M = 0
]

**Neutrality note:** (J^0) is **not** automatically taken as (q|\psi|^2) for a vacuum condensate; see the neutrality strategy in the main text.

---

### A.6 Projection operators and brane observables

Two standard operational choices:

**Restriction (slice):**
[
\mathcal P_0[f](x,y,z,t)=f(x,y,z,w=0,t).
]

**Weighted projection (mode overlap):**
[
\mathcal P_W[f](x,y,z,t)=\int_{-\infty}^{\infty}W(w)\,f(x,y,z,w,t)\,dw
]
with (W(w)\ge 0) normalized (often with truncated-window renormalization in numerics).

Brane density:
[
\rho_{\rm brane}(x,y,z,t)=\mathcal P[\rho](x,y,z,t).
]

Brane electric-type components:
[
E^{\rm brane}_a(x,y,z,t)=\mathcal P[F_{0a}](x,y,z,t),\quad a\in{x,y,z}.
]

---

### A.7 Geometry gate and confinement potential

* Gate/interior indicator:
  [
  G(\mathbf X;a,L)\in[0,1]
  ]
  smooth in (\mathbf X) and ((a,L)).

* Confinement potential encoding throat geometry:
  [
  V_{\rm conf}(\mathbf X;a,L)
  ]
  (Family-1 baseline: brane trap + radial wall + endcaps, with smooth transitions).

Geometry energy:
[
E_{\rm geom}(a,L)=P_{\rm vac}V(a,L)+\sigma A(a,L)+\cdots
]
(with optional higher-order/shape corrections as needed).

---

### A.8 Gauge fixing and localization terms

Gauge fixing parameter (\xi) (Lorenz gauge family).

Localization options (either may appear):

* (Z(\mathbf X;a,L)) multiplying (F_{MN}F^{MN}), and/or
* (m_\gamma^2(\mathbf X;a,L),A_MA^M) (Proca/Higgs-like suppression away from brane/throat).

---

### A.9 Units and scaling conventions

We keep units flexible but consistent; parameters are defined in the frozen config.

Typical choices:

* **Dimensionless/natural units** in numerics: set (\hbar=1), (m=1), and scale lengths by a chosen (L_0).
* **Physical-unit presentation** in papers: restore (\hbar), (m), (q), and any EM constants ((\mu_0), etc.) as needed.

**Rule:** every derived identity should remain dimensionally consistent when units are restored, and the unit checks should be automated in the Mathematica pipeline.

---

### A.10 Sign conventions (operational)

* “Outward” normal on the mouth surface (\Gamma) is defined to point away from the throat mouth (brane radial outward).
* Flux sign conventions for port observables (j_i(t)) follow this normal.
* Any symmetry splits in (w) (e.g., (w>0) vs (w<0)) use the standard sign of (w).

(Exact operational sign checks are enforced in the frozen definitions layer and should be referenced verbatim when implementing diagnostics.)

---

## Appendix B. Explicit Forms for (V_{\rm conf}), Gates, and Smooth-Step Derivatives

This appendix records the explicit “building blocks” used in the baseline Family-1 confinement model and its analytic geometry derivatives. The goal is to have one place that precisely states the mathematical objects that appear repeatedly in the force laws and the solver.

---

### B.1 Smooth primitives (baseline)

#### B.1.1 Smooth-step

We use a smooth approximation to the Heaviside step:
[
S(u) \equiv \frac{1+\tanh(u)}{2}
]
so that:
[
S'(u)=\frac{1}{2}\operatorname{sech}^2(u)=\frac{1}{2\cosh^2(u)}.
]

#### B.1.2 Smooth absolute value

We use a smooth absolute value:
[
|w|*\epsilon \equiv \sqrt{w^2+\epsilon^2}
]
so that:
[
\frac{d}{dw}|w|*\epsilon = \frac{w}{\sqrt{w^2+\epsilon^2}}.
]

These are chosen to ensure differentiability of all gate/potential terms with respect to both coordinates and geometry parameters ((a,L)).

---

### B.2 Coordinates and shorthand

[
R_3 \equiv \sqrt{x^2+y^2+z^2},\qquad \mathbf X=(x,y,z,w).
]

Geometry parameters and smoothing widths:

* radius (a),
* length (L),
* radial transition width (d_r),
* endcap transition width (d_w),
* smoothing (\epsilon_W) in (|w|).

(Names correspond to frozen numeric parameters; e.g. (d_r\leftrightarrow P.dr), (d_w\leftrightarrow P.dw), etc.)

---

### B.3 Gate definition (Family-1 baseline)

Define radial and longitudinal transition arguments:
[
u_R \equiv \frac{R_3-a}{d_r},
\qquad
u_W \equiv \frac{|w|_{\epsilon_W}-\frac{L}{2}}{d_w}.
]

Define step functions:
[
S_R \equiv S(u_R),\qquad S_W \equiv S(u_W).
]

Define “inside factors” (1 inside, 0 outside):
[
\mathrm{GateR} \equiv 1-S_R,\qquad
\mathrm{GateW} \equiv 1-S_W.
]

Define the interior gate:
[
G(\mathbf X;a,L) \equiv \mathrm{GateR}\,\mathrm{GateW}.
]

Interpretation:

* (G\approx 1) for (R_3<a) and (|w|<L/2) (inside the throat corridor),
* (G\approx 0) outside either the radial wall or the endcaps.

---

### B.4 Brane trap modulation and confinement potential (V_{\rm conf})

We define a (w)-harmonic trap whose stiffness depends on whether we are inside the throat:

[
\Omega^2(\mathbf X;a,L)
\equiv
\Omega_{\rm out}^2-\big(\Omega_{\rm out}^2-\Omega_{\rm in}^2\big)\,G(\mathbf X;a,L)
]

Then the “brane confinement” piece is:
[
V_{\rm brane}(\mathbf X;a,L)=\frac{1}{2}m\,\Omega^2(\mathbf X;a,L)\,w^2
]

The radial wall barrier:
[
V_{\rm wall}(\mathbf X;a)=V_0\,S_R^{p}
]

The endcap barrier:
[
V_{\rm cap}(\mathbf X;L)=V_0\,S_W^{p}
]

So the baseline Family-1 total confinement potential is:
[
\boxed{
V_{\rm conf}(\mathbf X;a,L)=V_{\rm brane}(\mathbf X;a,L)+V_{\rm wall}(\mathbf X;a)+V_{\rm cap}(\mathbf X;L)
}
]

---

### B.5 Analytic derivatives (\partial_a V_{\rm conf}) and (\partial_L V_{\rm conf})

These are the load-bearing formulas used in geometry force computation.

#### B.5.1 Intermediate derivatives

First,
[
\partial_a u_R = -\frac{1}{d_r}
\qquad
\partial_L u_W = -\frac{1}{2d_w}
]

And
[
\partial_a S_R = S'(u_R)\,\partial_a u_R = -\frac{S'(u_R)}{d_r}
\qquad
\partial_L S_W = S'(u_W)\,\partial_L u_W = -\frac{S'(u_W)}{2d_w}
]

Since (\mathrm{GateR}=1-S_R) and (\mathrm{GateW}=1-S_W):
[
\partial_a \mathrm{GateR} = +\frac{S'(u_R)}{d_r}
\qquad
\partial_L \mathrm{GateW} = +\frac{S'(u_W)}{2d_w}
]

Therefore for (G=\mathrm{GateR}\mathrm{GateW}):
[
\partial_a G = (\partial_a \mathrm{GateR})\,\mathrm{GateW}
=\frac{S'(u_R)}{d_r}\,\mathrm{GateW}
]
[
\partial_L G = \mathrm{GateR}\,(\partial_L \mathrm{GateW})
=\mathrm{GateR}\,\frac{S'(u_W)}{2d_w}
]

#### B.5.2 Derivatives of (\Omega^2)

Let (\Delta\Omega^2=\Omega_{\rm out}^2-\Omega_{\rm in}^2). Then:
[
\partial_a \Omega^2 = -\Delta\Omega^2\,\partial_a G
\qquad
\partial_L \Omega^2 = -\Delta\Omega^2\,\partial_L G
]

#### B.5.3 Brane term derivatives

[
\partial_a V_{\rm brane} = \frac{1}{2}m\,w^2\,\partial_a\Omega^2
\qquad
\partial_L V_{\rm brane} = \frac{1}{2}m\,w^2\,\partial_L\Omega^2
]

#### B.5.4 Wall term derivative

[
\partial_a V_{\rm wall} = -\frac{V_0 p}{d_r}\,S_R^{p-1}\,S'(u_R)
]

#### B.5.5 Endcap term derivative

[
\partial_L V_{\rm cap} = -\frac{V_0 p}{2d_w}\,S_W^{p-1}\,S'(u_W)
]

#### B.5.6 Total derivatives

[
\boxed{
\partial_a V_{\rm conf}=\partial_a V_{\rm brane}+\partial_a V_{\rm wall}
}
]
(because only the brane and wall terms depend on (a) in this baseline),

[
\boxed{
\partial_L V_{\rm conf}=\partial_L V_{\rm brane}+\partial_L V_{\rm cap}
}
]
(because only the brane and endcap terms depend on (L) in this baseline).

---

### B.6 “Half-length factor” sanity diagnostic

Because (u_W) contains (L/2), a common implementation bug is a missing factor of (1/2) in (\partial_L u_W). A robust check is to compute an alternate derivative using (L) instead of (L/2) and compare.

Operationally:

* define a variant (u_W^{(\mathrm{alt})} = (|w|_\epsilon - L)/d_w),
* compute (\partial_L V_{\rm cap}^{(\mathrm{alt})}),
* then check that the integrals satisfy approximately:
  [
  I_L^{(\mathrm{alt})} \approx 2\,I_L^{(\mathrm{base})},
  ]
  for symmetric configurations.

This diagnostic catches factor-of-two errors in (L)-derivatives that otherwise masquerade as “physics.”

---

### B.7 Notes on integer exponent (p)

In implementation, (p) is typically an integer exponent. For analytic derivatives, consistency requires:

* (V_{\rm wall}\propto S_R^p), (V_{\rm cap}\propto S_W^p),
* (\partial_a V_{\rm wall}\propto p S_R^{p-1}\partial_a S_R),
* (\partial_L V_{\rm cap}\propto p S_W^{p-1}\partial_L S_W).

If code uses `p_int = int(p)`, the spec assumes (p) is already intended as an integer and frozen as such.

---

## Appendix C. Diagnostics, Residuals, and Regression Tests

This appendix defines the “must-log” diagnostics that keep the unified 4D system honest. The theme is simple: **every claim (stability, resonance, equilibrium, impedance) must come with residuals that prove it isn’t a numerical artifact.**

---

### C.1 Geometry force residuals (the core equilibrium/stability check)

For geometry DOFs ((a,L)), define the matter-sector force integrals
[
I_a(a,L) \equiv \int d^4X\, \rho(\mathbf X)\,\partial_a V_{\rm conf}(\mathbf X;a,L)
\qquad
I_L(a,L) \equiv \int d^4X\, \rho(\mathbf X)\,\partial_L V_{\rm conf}(\mathbf X;a,L)
]
with (\rho=|\psi|^2).

Then the “residual” form used operationally is:
[
R_a(a,L) \equiv \big(-I_a(a,L)\big) - \partial_a E_{\rm geom}(a,L)
]
[
R_L(a,L) \equiv \big(-I_L(a,L)\big) - \partial_L E_{\rm geom}(a,L)
]

A scalar summary is:
[
|R| \equiv \sqrt{R_a^2 + R_L^2}.
]

**Interpretation**

* (R_a\approx 0): radial wall support balances geometric “closure” tendency.
* (R_L\approx 0): endcap/length-direction support balances geometric “closure” tendency.
* In the fully dynamic model, (R_a) and (R_L) also diagnose whether the geometry ODEs are being driven toward equilibrium or away from it.

---

### C.2 Symmetry diagnostics in (w): split checks for endcap forces

To ensure (\partial_L V_{\rm conf}) and integration are not biased, split the (I_L) integral into (w>0) and (w<0) halves:

[
I_L^{(+)} \equiv \int_{w>0} d^4X\,\rho\,\partial_L V_{\rm conf},
\qquad
I_L^{(-)} \equiv \int_{w<0} d^4X\,\rho\,\partial_L V_{\rm conf}.
]

Define the relative split diagnostic:
[
\mathrm{IL_split_rel}\equiv
\frac{|I_L^{(+)}-I_L^{(-)}|}{\max(\epsilon, |I_L|)}.
]

**What it catches**

* Asymmetric grids, indexing mistakes, or masking/dtype issues.
* Bugs where one endcap is effectively stronger numerically than the other.

**Rule of thumb**

* For symmetric potentials and symmetric initial conditions, (\mathrm{IL_split_rel}) should be very small (often (\ll 10^{-3}) if the solver is behaving well).

---

### C.3 “Half-length factor” diagnostic for (\partial_L) implementation

Because the endcaps use (|w| - L/2), a classic failure mode is dropping the (1/2) factor in (\partial_L).

Define:

* baseline (u_W = (|w|_\epsilon - L/2)/d_w)  → gives baseline (\partial_L V),
* alternate (u_W^{(\mathrm{alt})} = (|w|_\epsilon - L)/d_w)  → gives (\partial_L^{(\mathrm{alt})} V).

Then compute:
[
I_L^{(\mathrm{alt})} \equiv \int d^4X\,\rho\,\partial_L^{(\mathrm{alt})} V_{\rm conf}.
]

Expected relationship (when everything is consistent):
[
I_L^{(\mathrm{alt})} \approx 2\,I_L.
]

A convenient scalar error measure is:
[
\mathrm{IL_Lhalf_err} \equiv \left|\frac{I_L^{(\mathrm{alt})}}{\max(\epsilon,2|I_L|)} - 1\right|.
]

**Interpretation**

* (\mathrm{IL_Lhalf_err}\approx 0): factor-of-two is correct.
* (\mathrm{IL_Lhalf_err}\approx 1) or (\approx 2): strong evidence of a missing factor or sign.

---

### C.4 Derivative-method convergence diagnostics ((\partial_a V), (\partial_L V))

When computing (\partial_{a,L} V_{\rm conf}), we want to detect cancellation, step-size issues, or precision issues.

**Methods**

* **Analytic derivatives** (preferred): compute (\partial_a V), (\partial_L V) directly from the smooth-step forms.
* **FD in float64** (fallback check):
  [
  \partial_L V \approx \frac{V(a,L+\Delta L)-V(a,L-\Delta L)}{2\Delta L},
  ]
  with (\Delta L \sim \text{(grid spacing in }w\text{)}) and (V) evaluated in float64.

**Recommended regression table**
For a fixed ((a,L)), log:

* method (analytic, fd64 with several (\Delta L)),
* (I_L), (R_L),
* (\mathrm{IL_split_rel}),
* and (if analytic) (\mathrm{IL_Lhalf_err}).

**What you’re looking for**

* fd64 values approaching analytic as (\Delta L) becomes “grid-scale” (not tiny).
* (I_L) should not change wildly with small method changes if the implementation is stable.

---

### C.5 Grid convergence tests (what “converged” should mean here)

For a chosen ((a,L)), repeat the same calculation at increasing grid size (N) (keeping physical box size fixed). Track:

* (I_a(N), I_L(N)),
* (R_a(N), R_L(N)),
* (\mathrm{IL_split_rel}(N)).

**Important caution**
For steep walls/endcaps (large (V_0), sharp (d_r,d_w)), convergence can be non-monotone because the wall resolves differently as (N) changes. That’s not automatically “bad,” but it means you should:

* also track convergence of *profiles* (e.g., brane density near the wall),
* and consider whether the physical model wants sharp barriers or a softer, physically motivated wall thickness.

---

### C.6 Brane projection validity diagnostics (slice vs weighted projection)

If you use a weighted brane projection (W(w)), you must log:

**Tail mass**
[
\mathrm{TailMass} \equiv 1 - \int_{|w|\le w_{\max}} W_{\rm trunc}(w)\,dw
]
(where (W_{\rm trunc}) is the renormalized truncated weight used numerically).

**Projection sensitivity**
Compare brane observables computed by:

* restriction (w=0),
* weighted projection (\mathcal P_W),

in regimes where confinement is strong and both should agree closely.

**Rule of thumb**
If tail mass is not small in the “brane-like” regime, the model is leaking into (w) even away from the throat, and your “brane observable” interpretation is unstable.

---

### C.7 Coupled-system residuals (when gauge fields and dynamics are turned on)

Once (A_M) and dynamic ((a(t),L(t))) are included, add:

#### C.7.1 Continuity / charge conservation residual

[
\mathcal R_{\rm cont} \equiv |\partial_M J^M|
]
(choose an (L^2) norm over the domain, and always log it normalized by a scale like (|J^0|/T) or (|J|/L)).

If the chosen open-system terms intentionally violate strict conservation, log the expected source term explicitly and compute:
[
\mathcal R_{\rm cont} \equiv |\partial_M J^M - S_{\rm cont}|.
]

#### C.7.2 Gauge-condition residual (if using Lorenz gauge)

[
\mathcal R_{\rm gauge} \equiv |\partial_M A^M|.
]

#### C.7.3 Energy balance residual

Track:

* total energy (H_{\rm tot}(t)),
* injected power (P_{\rm drive}(t)),
* dissipated power (P_{\rm diss}(t)),
* boundary flux (\Phi_{\rm boundary}(t)).

Define a running residual:
[
\mathcal R_E(t)\equiv
\left| \frac{dH_{\rm tot}}{dt} - P_{\rm drive} + P_{\rm diss} + \Phi_{\rm boundary} \right|.
]

This is the primary guardrail against “numerical stabilization.”

#### C.7.4 Geometry work consistency

If geometry evolves:
[
W_{\rm geom}(t)\equiv F_a\dot a + F_L\dot L
]
and it should appear consistently as exchange work between geometry and the fields in the energy accounting.

---

### C.8 Suggested “overnight run” logging format

For long runs, log at two cadences:

**Fast cadence (every few steps)**

* (a(t),L(t)), (\dot a,\dot L)
* (H_{\rm tot}(t))
* (\mathcal R_{\rm cont}(t)), (\mathcal R_{\rm gauge}(t)), (\mathcal R_E(t))

**Slow cadence (every many steps)**

* full profiles (slices in (w) and in (R_3))
* brane observables (port time series windows)
* resonance/response extraction outputs (if in driven mode)

And always write a compact JSON “summary row” at the end of each run with:

* parameter hashes (freeze SHA),
* grid/domain,
* solver dt,
* and the max/median of each residual.

---

### C.9 Minimal regression suite (what should run automatically)

1. **Analytic vs fd64 derivative comparison** at one representative ((a,L)).
2. **Half-length factor check** ((\mathrm{IL_Lhalf_err})).
3. **(w)-split symmetry check** ((\mathrm{IL_split_rel})).
4. **Grid convergence** at 2–3 (N) values you can afford.
5. **Projection tail-mass check** in a brane-like configuration.
6. **(When coupled)** continuity, gauge, and energy residual checks under:

   * no-drive,
   * drive-only,
   * intake-only,
   * mixed drive+intake.

If any of these fails, treat every “physics conclusion” downstream as suspect until fixed.

---

## Appendix D. Implementation Notes for Overnight “Intelligent” Searches

This appendix is a practical guide for turning the current scan/root-finding tooling into an overnight runner that (i) hunts efficiently for equilibria or operating points, (ii) avoids wasting GPU time, and (iii) produces logs that immediately tell us whether we found real physics or a numerical mirage.

---

### D.1 Why “intelligent search” is required (what brute scans miss)

A brute grid scan over ((a,L)) is expensive because each point requires:

* building (V_{\rm conf}(\cdot;a,L)),
* relaxing (\psi) (and later also solving for (A_M)),
* computing force integrals and residuals,
* repeating at multiple resolutions for convergence confidence.

Meanwhile, the residual landscape can have:

* narrow valleys,
* “flat” directions (especially in (L) if selection physics is weak),
* and strong dependence on closure assumptions (fixed norm vs reservoir).

So the overnight runner should:

* warm-start aggressively,
* adapt step sizes,
* and only refine where diagnostics indicate real progress.

---

### D.2 What an overnight runner should output (minimum artifacts)

Every run should write:

1. **A single JSON summary row per evaluated point**

   * ((a,L)), (R_a,R_L,|R|)
   * (I_a,I_L), (\partial_aE_{\rm geom},\partial_LE_{\rm geom})
   * symmetry checks: (\mathrm{IL_split_rel}), (\mathrm{IL_Lhalf_err})
   * solver stats: iterations, final error, time, GPU memory peak (if possible)
   * freeze hash + config path

2. **A CSV table** for quick sorting/plotting.

3. **Optional checkpoints**

   * (\psi) saved for best (k) points (to restart later),
   * in the coupled model, (A_M) checkpoints too.

4. **A “run ledger”**

   * what points were tried,
   * which were rejected early (and why),
   * what the best point was at each stage.

If this is present, you can stop and resume overnight runs without losing context.

---

### D.3 Smart evaluation: early exit and staged accuracy

Use a two-tier (or three-tier) evaluation strategy:

#### Tier 0: cheap precheck (optional but helpful)

Before any relaxation:

* compute geometric derivatives (\partial_aE,\partial_LE),
* check if ((a,L)) are inside admissible ranges,
* check if the confinement potential walls are resolvable on the current grid (e.g., (d_r) and (d_w) not below grid spacing).

Reject points that are obviously ill-posed.

#### Tier 1: probe solve (fast approximate)

* run imaginary-time relaxation for a small number of steps (or higher tolerance),
* compute approximate (R_a,R_L),
* if (|R|) is far worse than the current best (by a margin), discard without refinement.

#### Tier 2: full solve (expensive, only for promising points)

* run to strict tolerance,
* compute the full diagnostic suite,
* optionally run a short “stability poke” (small perturbations) to ensure the result isn’t a fragile minimizer artifact.

This staged approach is what makes overnight search feasible.

---

### D.4 Warm-start policy (the biggest speed win)

Warm-starting should be used systematically:

* **Within a local search:** when stepping ((a,L)) slightly, initialize (\psi) with the previous converged (\psi).
* **Across tiers:** promote the Tier-1 (\psi) as the initial guess for Tier-2.
* **Across restarts:** save the best (\psi) periodically and reload at the start of a new run.

If you implement only one “intelligent” feature, make it this.

---

### D.5 Adaptive ((a,L)) exploration strategies

There are multiple workable strategies; the runner should support at least two so we can cross-check.

#### Strategy S1: Trust-region Newton on residuals (local root find)

* compute (R(a,L)),
* estimate Jacobian (J=\partial(R_a,R_L)/\partial(a,L)) by finite differences *using probe solves*,
* compute Newton step (\Delta = -J^{-1}R),
* accept via line search.

**Key practical improvement:** Jacobian evaluations should use probe solves (few steps), then only after a candidate step is accepted should you do a full solve to validate the improvement.

#### Strategy S2: Coordinate descent with bracketed line search

Alternate:

* hold (L) fixed, search (a) to reduce (|R_a|),
* hold (a) fixed, search (L) to reduce (|R_L|).

This is robust when one direction is stiff and the other is flat.

#### Strategy S3: Bayesian / surrogate optimization (for expensive coupled runs)

For the coupled ((\psi,A)) problem, each evaluation will be much more expensive. Use:

* a Gaussian process or random forest surrogate for (|R|),
* propose new points via expected improvement,
* keep a “trust region” around the best known basin.

This is ideal for overnight on a big GPU cluster.

---

### D.6 Handling “flat (L)” and missing selection physics

If (R_L) remains stubbornly negative across wide (L) ranges (as current scans suggest), that may indicate the model lacks a strong (L)-selecting mechanism in the static closure.

The overnight runner should therefore include:

1. **A monotone diagnostic scan**

   * for fixed (a), scan (L) downwards and upwards,
   * record (I_L(L)), (dE/dL), and (R_L(L)),
   * detect whether (R_L) crosses zero or is monotone.

2. **Automatic “model deficiency flag”**

   * if (R_L(L)) never crosses zero and is nearly constant in sign/magnitude, label:

     * “no (L)-root in this closure”
   * and automatically recommend the next knob:

     * change closure (fixed (\mu) instead of fixed (N)),
     * add a shape/curvature term,
     * or turn on dynamic field support (EM cavity / intake).

This prevents burning compute on a root that cannot exist in the current model.

---

### D.7 Precision and cancellation controls (keep derivatives honest)

The runner should enforce:

* analytic (\partial_{a,L}V_{\rm conf}) as default,
* float64 accumulation for integrals when running on GPU,
* logging of (\mathrm{IL_Lhalf_err}) and (\mathrm{IL_split_rel}) at every evaluated point,
* automatic fallback comparisons to fd64 for a small fraction of points (spot checks).

If any derivative diagnostic spikes, the point is marked “suspect” and not treated as evidence of physics.

---

### D.8 Multi-resolution confirmation (“trust but verify”)

For any candidate “best” point (or any claimed equilibrium), run a **confirmation ladder**:

* N = 64 (or baseline) → full solve
* N = 80/96 → warm-start from interpolated (\psi) or from previous N result
* N = 128 (if possible) → final confirmation

Accept a claimed equilibrium only if:

* residuals do not change qualitatively,
* symmetry diagnostics remain small,
* and brane-observable outputs are stable.

This ladder is exactly what you did manually with the A100 run; the runner should automate it.

---

### D.9 A recommended overnight schedule (static-first, then coupled)

**Overnight run mode 1 (current state, matter-only):**

1. pick an (a)-band around the best current value,
2. scan (L) downward more aggressively (since best points are at lower L in the table),
3. run staged evaluation + trust-region search,
4. if no (L)-root emerges, stop early and label “needs selection physics.”

**Overnight run mode 2 (next phase, coupled):**

1. keep ((a,L)) dynamic with wall-law ODEs,
2. apply a simple cavity drive (J^N_{\rm ext}) in Maxwell,
3. allow intake/reservoir coupling,
4. search not for static roots but for bounded steady operation:

   * bounded (a(t),L(t)),
   * steady periodic response,
   * stable impedance signature.

This aligns with the physical system better than static equilibrium hunting.

---

### D.10 Minimal code organization recommendations

To make this maintainable:

* A single `RunnerConfig` that records:

  * domain, grid, solver tolerances,
  * search strategy,
  * warm-start/checkpoint policy,
  * logging paths.

* A single `evaluate_point(a,L, tier, psi0, A0)` function that:

  * runs solve,
  * returns a structured result dict including all diagnostics,
  * and optionally returns checkpoint fields.

* A search driver that chooses points and uses `evaluate_point`.

This modularity is what makes it safe to “leave it running overnight.”

---

## Appendix E. Open-System Closures, Drive Models, and Dissipation Templates

This appendix records concrete, “drop-in” mathematical options for making the unified throat system **fully dynamic and open**, so it can represent (i) cavity pumping, (ii) intake/through-flow, (iii) bulk receiving-region back-pressure, and (iv) energy shedding when density changes. These are the pieces that *must exist* once we move beyond static/fixed-norm equilibrium searches.

---

### E.1 Why closure matters (the core point)

A conservative GNLS + Maxwell system in a finite periodic box will generally:

* conserve total norm (or charge) and energy (modulo gauge work),
* reflect waves unless special boundaries exist,
* and often produce “equilibria” that depend on the box more than the physics.

To model the physical throat story, we need **explicit choices** for:

1. what is held fixed (norm? chemical potential? total energy budget?),
2. how energy and mass/charge enter and leave (drive + dissipation + boundary flux),
3. how the receiving region behaves (back-pressure emerges from EOS + geometry + leakage).

---

### E.2 Matter-sector closures (choices for (\psi))

#### E.2.1 Fixed norm (N) (diagnostic-only)

Hold
[
N \equiv \int d^4X\,|\psi|^2 = N_0
]
by renormalizing (\psi) each step (or in imaginary-time relaxation).
**Use:** finding ground states / diagnostic equilibria.
**Do not use:** modeling intake, shedding, or reservoir physics.

#### E.2.2 Fixed chemical potential (\mu) (reservoir coupling)

Evolve (\psi) with a dissipative relaxation toward a target chemical potential:
[
i\hbar D_t\psi = \cdots - i\hbar\,\gamma_\mu \big(\mu_{\rm loc}(\mathbf X,t)-\mu_0\big)\psi
]
where a practical definition of local chemical potential is:
[
\mu_{\rm loc} \equiv \frac{1}{\psi}\left[-\frac{\hbar^2}{2m}D_iD_i\psi + V_{\rm conf}\psi + U'(|\psi|^2)\psi\right]
]
(using a stabilized real-part extraction in numerics).

**Interpretation:** a bath tries to maintain (\mu\approx\mu_0) in chosen regions (often far from the throat).

**Implementation note:** apply this term only in a “reservoir zone” (S_{\rm res}(\mathbf X)) away from the throat:
[
\gamma_\mu(\mathbf X)=\gamma_{\mu0}\,S_{\rm res}(\mathbf X).
]

#### E.2.3 Fixed far-field density (\rho_0) (soft density pinning)

Add a soft restoring term that nudges density toward a target (\rho_0):
[
i\hbar D_t\psi = \cdots - i\hbar\,\gamma_\rho\,S_{\rm res}(\mathbf X)\big(|\psi|^2-\rho_0\big)\psi
]

This is easy and robust for sustaining an ambient “vacuum density” that can feed intake.

#### E.2.4 Energy-budget closure (formalizing “density drop → energy shed”)

If we want the rest-energy bookkeeping to directly constrain throat dynamics, we can enforce an energy budget at the *system* or *region* level.

Define a throat-region weight (W_{\rm throat}(\mathbf X;a,L)) (often built from the gate (G)). Define throat energy:
[
E_{\rm throat}(t)=\int d^4X\,W_{\rm throat}(\mathbf X;a,L)\,\mathcal H_\psi(\mathbf X,t)
]
and impose a rule such as:
[
\dot E_{\rm throat} = P_{\rm in} - P_{\rm out} - P_{\rm geom}
]
with each term computed from explicit fluxes and geometry work.

This is less “plug-and-play” than (\mu)-closure but matches your stated physical priority. Practically, it becomes a diagnostic constraint used to **tune** (\gamma_\rho,\gamma_\mu), and absorber strengths so that the observed energy flow matches the intended closure.

---

### E.3 Drive models (how you pump the cavity / excite ports)

#### E.3.1 EM-sector drive via external current (J^M_{\rm ext})

Drive the gauge field by adding:
[
\partial_M(\cdots F^{MN}) + \cdots = J^N + J^N_{\rm ext}.
]

A practical cavity-mode drive is a time-harmonic current localized near the mouth or along the throat corridor:
[
J^N_{\rm ext}(\mathbf X,t)=\Re\!\Big\{\hat J^N(\mathbf X)\,e^{-i\omega t}\Big\}
]
with spatial envelope (\hat J^N) chosen to match a desired mode symmetry (e.g., axial, dipole, quadrupole).

**Where to localize:** use gate-weighting for throat targeting:
[
\hat J^N(\mathbf X)=J_0\,G(\mathbf X;a,L)\,\hat p^N(\mathbf X)
]
with (\hat p^N) a polarization pattern.

#### E.3.2 Matter-sector drive via forcing term (\eta(\mathbf X,t))

Add a source to the GNLS:
[
i\hbar D_t\psi = \cdots + \eta(\mathbf X,t)
]
where (\eta) is localized at an intake port / reservoir boundary.

This can directly drive density waves or induce through-flow.

#### E.3.3 Port drive consistent with the measurement basis

If you already have port basis functions (P_i(\theta,\phi)) on (\Gamma), define a drive aligned to a port mode:

* choose a mouth shell near (\Gamma) in the brane coordinates and apply a forcing whose angular dependence matches (P_i).

This is the cleanest way to connect “drive amplitude” to “measured port effort/flux” and later to (Z^{\rm eff}(\omega)).

---

### E.4 Dissipation and absorbers (how energy leaves without reflecting)

#### E.4.1 Matter sponge / complex absorbing potential (CAP)

Add a damping term in an outer layer:
[
i\hbar\,\partial_t\psi = \cdots - i\hbar\,\gamma_\psi(\mathbf X)\psi
]
with (\gamma_\psi(\mathbf X)) nonzero only near boundaries.

This removes outgoing matter waves and prevents artificial back-pressure from reflections.

#### E.4.2 EM sponge / PML-lite

For Maxwell in potential form, a simple absorber is:
[
\partial_t^2 A^N + \gamma_A(\mathbf X)\,\partial_t A^N + \cdots = \text{sources}
]
with (\gamma_A) nonzero near boundaries.

This is not a full PML, but it is often sufficient for first prototypes if validated by domain enlargement tests.

#### E.4.3 Localized dissipation in the throat (to model conversion/loss)

If you want the throat to have intrinsic non-ideality (e.g., mode-to-heat conversion), you can use a throat-gated damping:
[
\gamma_{\rm throat}(\mathbf X)=\gamma_0\,G(\mathbf X;a,L).
]

Be explicit: throat-gated dissipation can stabilize (or destabilize) the system. It must be logged and treated as a physical knob, not a numerical patch.

---

### E.5 Receiving-region modeling (how 4D expansion creates back-pressure)

Back-pressure is not a boundary condition; it’s a **result** of:

* EOS pressure from (U(\rho)),
* geometric spreading in 4D,
* and any confining potentials in the receiving region.

To make this real, ensure:

1. the confinement beyond the throat actually permits bulk occupancy (the gate truly opens access),
2. the receiving region is large enough and/or has absorbers so that matter can spread or exit rather than reflect,
3. you measure back-pressure via a defined observable (e.g., brane-projected enthalpy near (\Gamma) relative to far-field baseline).

A useful operational diagnostic is:
[
\Delta h_{\rm back}(t) \equiv \overline{h}_{\rm near\,mouth}(t)-\overline{h}_{\rm far}(t)
]
where (h(\rho)) is the enthalpy implied by your EOS and the overlines denote averages in specified regions.

---

### E.6 Coupled energy accounting (what must be logged)

To support claims about shedding and stabilization, always log:

* (H_\psi(t)): matter energy (including confinement and nonlinearity)
* (H_{\rm EM}(t)): EM energy (plus localization contributions)
* (E_{\rm geom}(a(t),L(t)))
* injected powers:

  * (P_{\rm drive,\psi}), (P_{\rm drive,EM})
* dissipated powers:

  * (P_{\rm diss,\psi}), (P_{\rm diss,EM})
* boundary fluxes:

  * (\Phi_{\rm boundary,\psi}), (\Phi_{\rm boundary,EM})
* geometry work:
  [
  W_{\rm geom} = F_a\dot a + F_L\dot L.
  ]

Then the “sanity residual” is:
[
\mathcal R_E(t)=\left|\frac{d}{dt}(H_\psi+H_{\rm EM}+E_{\rm geom}) - P_{\rm drive} + P_{\rm diss} + \Phi_{\rm boundary}\right|.
]

If (\mathcal R_E) is not small relative to the power scale, don’t trust the run.

---

### E.7 Recommended “first coupled runs” recipe

A practical starting configuration that matches the physical intent:

1. **Neutrality:** background-subtracted charge or two-field (must be chosen).
2. **Localization:** choose *one* minimal localization mechanism (e.g., simple (m_\gamma(w)) profile).
3. **Drive:** EM current drive localized near the throat corridor (gate-weighted), single frequency (\omega).
4. **Matter reservoir:** far-field density pinning toward (\rho_0) in reservoir zones.
5. **Absorbers:** matter + EM sponge layers at domain edges.
6. **Geometry:** wall-law ODEs with moderate damping (so you see dynamics but avoid numerical blow-up).
7. **Outputs:** ports (u_i,j_i) + residuals (\mathcal R_{\rm cont}, \mathcal R_{\rm gauge}, \mathcal R_E).

This setup is specifically designed to let you see:

* whether EM pumping can produce stabilizing stress,
* whether intake/back-pressure emerges without reflections,
* and whether (a(t),L(t)) settle, oscillate, or diverge.

---

### E.8 What “success” looks like (in this open-system context)

A “successful” run is **not** “found a static root.” It is:

* (a(t),L(t)) remain bounded for long times (possibly periodic),
* brane observables show a stable resonance/response signature,
* energy accounting closes (small (\mathcal R_E)),
* and results persist under domain enlargement / absorber strengthening.

If those are satisfied, then we can safely proceed to reductions and “what pops out” parameter extraction.
