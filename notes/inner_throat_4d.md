Below is a **Paper 7 “PDE-first throat model” spec** you can reuse in future sessions. It’s written as a practical blueprint: what we assume, what we solve, what we measure/extract, and what the milestones/sanity checks are.

---

# Paper 7 Spec: 4D Throat PDE First, Then Brane Coupling

## Why this direction

We want a model that:

1. **Actually contains the 4D throat physics** (not just an imposed impedance curve),
2. Builds in the **two “keep it open” effects** you highlighted:

   * standing-wave (mode) pressure,
   * and the superfluid’s own flow/momentum flux pushing back,
3. Keeps **energy accounting explicit** (so we can define an effective mass via (E = m c^2)),
4. Automatically includes **density–flow self-coupling** (Bernoulli-like effects), and
5. Still produces the correct interface object to couple to the 3D brane: a **mouth response operator** (DtN / impedance / port matrix), extracted from the PDE solution.

The earlier cylinder DtN work becomes a **unit test** and “known-answer limit,” not the final model.

---

## 1) Geometry and degrees of freedom

### 1.1 Throat domain (4D “bulk” region)

Define a throat region (\mathcal{T}\subset\mathbb{R}^4) with coordinates, e.g.
[
\mathbf{X} = (x,y,z,w),
]
where ((x,y,z)) are brane-like spatial coordinates and (w) is the extra/bulk spatial coordinate.

Start with a **4D cylinder-like** throat for validation, then generalize to a **rounded funnel family**:

* “length” parameter (L),
* “radius” parameter (a),
* optional shape parameters for rounding/taper.

Operationally, implement geometry via a confining potential (V_{\text{wall}}(\mathbf{X};a,L,\text{shape})) that is ~0 inside (\mathcal{T}) and large outside.

### 1.2 Geometry DOF (what can move)

At minimum, allow the throat to “breathe”:

* (a(t)): effective radius
* (L(t)): effective length

Optional upgrade:

* (a(z,t)) (or (a(z,w,t))) for funnel deformation.

This is what lets the throat “stay open” by **force balance**, rather than by assumption.

---

## 2) Parent PDE: 4D Gross–Pitaevskii Equation (GPE)

### 2.1 Governing equation (bulk field)

Use a complex order parameter (\psi(\mathbf{X},t)) on (\mathcal{T}):
[
i\hbar,\partial_t \psi
======================

\left(
-\frac{\hbar^2}{2m}\nabla_4^2

* V_{\text{wall}}(\mathbf{X};a,L,\text{shape})
* g|\psi|^2
  \right)\psi
  ;+;S_{\text{drive}}.
  ]

- (\nabla_4^2) is the Laplacian in ((x,y,z,w)).
- (g>0) sets compressibility / interaction strength.
- (S_{\text{drive}}) is optional and represents forcing at the mouth or volume driving (used when we want a specific standing wave).

### 2.2 Hydrodynamic form (what it means physically)

Write (\psi=\sqrt{\rho},e^{i\theta}). Then:

* Velocity: (\mathbf{v}=(\hbar/m)\nabla_4\theta)
* Continuity:
  [
  \partial_t\rho+\nabla_4\cdot(\rho\mathbf{v})=0
  ]
* Bernoulli/Euler-like (including quantum pressure (Q)):
  [
  \hbar,\partial_t\theta + \frac{m}{2}v^2 + g\rho + V_{\text{wall}} - Q = 0,
  \quad
  Q=\frac{\hbar^2}{2m}\frac{\nabla_4^2\sqrt{\rho}}{\sqrt{\rho}}.
  ]

This is where your “recursive” effect lives:

* increasing flow speed (v) tends to reduce (\rho) (at fixed chemical potential),
* and the local sound speed (c_s) depends on (\rho) (in GPE: (c_s^2 \propto g\rho/m)).

So yes: **GPE already encodes the Bernoulli-like density changes and self-consistency**.

---

## 3) Energy, pressure, and “stays open” force balance

### 3.1 Fluid energy functional inside the throat

Define the GPE energy (for a fixed geometry):
[
E_{\text{fluid}}[\psi;a,L]
==========================

\int_{\mathbb{R}^4}
\left[
\frac{\hbar^2}{2m}|\nabla_4\psi|^2

* V_{\text{wall}}|\psi|^2
* \frac{g}{2}|\psi|^4
  \right]d^4X.
  ]

This energy automatically includes:

* standing-wave energy (gradient term),
* flow energy (phase gradients),
* compressibility/interaction energy (the (g|\psi|^4) term),
* density variations induced by flow.

### 3.2 Geometry “pushback” energy (vacuum/tension cost)

Add an explicit geometry cost (toy model analogue of vacuum work / throat surface cost):
[
E_{\text{geom}}(a,L)=P_{\text{vac}},V(a,L)+\sigma,A(a,L)+\kappa_b,(\text{curvature cost})+\cdots
]

* (P_{\text{vac}}V) is the simplest “energy cost to carve/maintain volume.”
* (\sigma A) is a surface-energy term (optional but often stabilizing).
* curvature cost helps for rounded funnel stability.

### 3.3 Total energy and equilibrium condition

Total:
[
E_{\text{total}}(a,L;\psi)=E_{\text{fluid}}[\psi;a,L]+E_{\text{geom}}(a,L).
]

The throat “stays open” when geometry and fluid equilibrate:
[
\frac{d}{da}\Big(E_{\text{total}}\Big)=0,\qquad
\frac{d}{dL}\Big(E_{\text{total}}\Big)=0,
]
with (\psi) chosen self-consistently for the given ((a,L)).

This captures your “two opening forces” automatically:

* standing wave contributes via (E_{\text{fluid}}) and wall stress,
* through-flow contributes via momentum/pressure terms embedded in the same functional (or via stress tensor evaluation).

### 3.4 Mass–energy equivalence for the toy model

Define an effective rest mass of the “particle/throat object”:
[
m_{\text{eff}}(a,L)\equiv \frac{E_{\text{total}}(a,L)}{c^2},
]
where (c) is the toy model’s characteristic signal speed (your “sound-speed-as-light-speed” role).

Density dependence naturally enters because (E_{\text{fluid}}) depends on (\rho); lower background density lowers available energy and changes equilibrium geometry.

---

## 4) Boundary conditions and physical scenarios

We’ll use distinct regimes as **unit tests / operating modes**:

### 4.1 Closed-cavity consistency mode (recover (L/a\approx 1.85))

Purpose: reproduce the earlier series’ variational selector.

* Geometry: cylinder-ish throat
* Mode: Dirichlet-like radial selection (x_{01}) and axial fundamental (k_z=\pi/L)
* No net through-flow (or negligible)
* Minimize (E_{\text{total}}) with respect to ((a,L)) under the same constraints used in the earlier enthalpy argument (often fixed “charge”/particle content (N=\int\rho,d^4X) or fixed effective mass scale).

**Expected result:** (L/a \rightarrow \sqrt{2}\pi/x_{01}\approx 1.8475) in the appropriate limit.

This is the key “don’t break the series” check.

### 4.2 Open-mouth / coupled mode (the actual Paper 7 goal)

Purpose: throat interacts with the 3D brane.

* Mouth region (\Gamma_{\text{mouth}}) couples to outer field
* Implement either:

  * driven boundary forcing (S_{\text{drive}}) localized near the mouth, or
  * a boundary condition that fixes chemical potential / phase at the mouth.

### 4.3 Through-flow mode (superfluid rushing through)

Purpose: include the second “opening” effect explicitly.

* Impose a flux or chemical potential difference across the throat:

  * Dirichlet in (\theta) (phase difference) between ends, or
  * fixed current at a boundary region.

This produces density changes and pressure gradients automatically.

---

## 5) The thing the 3D brane needs: the mouth response operator

Even with a full PDE, the coupling to the 3D brane is done through a compact interface object.

### 5.1 Define port variables at the mouth

Choose one or a few mouth “ports” (basis functions on the mouth cross-section). Minimal:

* (P_0): uniform (monopole) port
* (P_2): quadrupole-like port (zero mean)

Define port amplitudes (example):
[
u_i(t)=\int_{\Gamma_{\text{mouth}}} P_i(\mathbf{s}),\phi(\mathbf{s},t),dS,
]
and conjugate fluxes (e.g. normal current)
[
j_i(t)=\int_{\Gamma_{\text{mouth}}} P_i(\mathbf{s}),\mathbf{J}(\mathbf{s},t)\cdot \hat{n},dS,
\quad
\mathbf{J}=\frac{\hbar}{m}\Im(\psi^*\nabla\psi).
]

### 5.2 Frequency-domain response / impedance matrix

Drive port (i) sinusoidally at frequency (\omega) and measure response:
[
j_i(\omega)=\sum_j Z_{ij}(\omega),u_j(\omega).
]

This (Z(\omega)) is the **numerical DtN/impedance** extracted from the full PDE.

### 5.3 Outer (3D brane) operator

Outer side has its own DtN map (Z_{\text{out}}(\omega)) (multipole-by-multipole or in the same port basis). Matching at the interface is:
[
Z_{\text{in}}(\omega)\approx Z_{\text{out}}(\omega)
]
in the chosen basis, subject to coupling weights/symmetry.

This is where your previous Step 4 machinery plugs back in cleanly—except now (Z_{\text{in}}) is generated by the PDE solver rather than a closed-form waveguide formula.

---

## 6) 2PN/EFT compatibility gates (sanity checks)

We don’t need full 2PN GR matching now, but we do want to avoid building something that can’t even live in a PN-style EFT.

**Gate A: Low-frequency locality**
Check that (Z_{ij}(\omega)) admits a low-(\omega) expansion away from resonances:
[
Z_{ij}(\omega)=A_{ij}\omega^2 + B_{ij}\omega^4 + \cdots
]
(or the appropriate analytic structure depending on variables). This is what makes a local-in-time effective action possible.

**Gate B: Passivity**
With damping/leakage included, ensure energy flux is sign-definite under a consistent convention (no “free energy” artifacts).

**Gate C: Point-throat limit**
As (a\to 0), outer multipoles beyond monopole should be suppressed appropriately; no extra long-range “hair” should persist unless explicitly intended.

---

## 7) Numerical strategy and milestones

### Milestone 1: Cylinder validation (unit-test stage)

* Implement 4D GPE in cylinder geometry.
* Linearize / small-amplitude limit: extract (Z(\omega)) numerically.
* Compare against analytic cylinder DtN results (your existing Mathematica work) in the regime where the mapping should coincide.

### Milestone 2: Recover the series 1.85 selector (consistency stage)

* Run closed-cavity energy minimization with Dirichlet-like (x_{01}) branch.
* Verify the preferred (L/a) approaches (1.8475) in that limit.

### Milestone 3: Turn on through-flow (two-force stage)

* Impose flux/phase bias.
* Observe self-consistent density reductions (Bernoulli-like), changes in (c_s(\rho)), and how equilibrium ((a,L)) shifts.

### Milestone 4: Extract mouth operator in funnel geometry (Paper 7 core)

* Replace cylinder with rounded funnel family.
* Extract (Z(\omega)) (or small port matrix) from PDE.
* Identify regimes where (Z(\omega)) is PN-local vs resonance-dominated.

### Milestone 5: Couple to outer 3D brane solver (matching stage)

* Use outer DtN operator (static + radiating).
* Match in the port basis; compute mixing diagnostics and failure radii/regions.

---

## 8) “Make it fail” checklist (paper-worthy)

* No equilibrium ((a,L)) exists for any reasonable (P_{\text{vac}},\sigma).
* Equilibrium exists but is dynamically unstable (mode growth).
* Through-flow collapses density too much → loss of confinement / breakdown of assumptions.
* Resonances dominate in the frequency band of interest → no local EFT truncation.
* Extracted (Z(\omega)) depends strongly on numerical knobs (resolution/damping), meaning the model is not predictive.

---

## Next steps for the very next session

1. Decide the **minimal geometry** for the first full-PDE run: 4D cylinder vs “3D cylinder + extra coordinate as a parameter.” (Start with the simplest that still exercises (\nabla_4).)
2. Choose which “conservation constraint” defines a particle at rest:

   * fixed norm (N=\int|\psi|^2),
   * or fixed (E_{\text{total}}/c^2) target mass.
3. Write down the explicit (E_{\text{geom}}(a,L)) you want (start with (P_{\text{vac}}V+\sigma A)).
4. Implement Milestone 1 (cylinder) and verify the solver reproduces a known DtN response in the linear regime.
5. Implement Milestone 2 (closed-cavity selector) and confirm the 1.85 limit appears in the right regime.

That gives you a clean, defensible foundation for Paper 7: **a real 4D PDE throat model + an extracted mouth operator + energy/mass bookkeeping + the stability/identifiability story.**

===

Yep — this spec is *very* aligned with what `brane_bulk_ontology.tex` and `1pn_hybrid.tex` are implicitly asking for: **a real inner throat model** that (i) reproduces the “sphere outside / cylinder inside” ontology, and (ii) outputs the one thing the 3D brane needs: a **mouth response operator** (DtN / impedance) that can be plugged back into the 1PN + strong-field machinery.

Here’s what I think is strongest, and what I’d tighten so Paper 7 doesn’t accidentally drift away from your already-calibrated toy-universe constraints.

---

## What’s already solid

### The deliverable is exactly right: “inner PDE → mouth operator”

That’s the whole bridge between the ontology and the hybrid paper:

* Ontology: far field must look **monopolar** with corrections suppressed like ((a/r)^2), while the **interior** supports **cylindrical Bessel modes** and recovers the EM aspect ratio selection.
* Hybrid: strong-field deviations are controlled by **(\Lambda=L/a)** plus **impedance effects** (your words in the abstract are basically “give me (\Lambda) and (Z(\omega)), then I can talk about horizons/photon sphere corrections”).

So the spec’s structure — *PDE first, then extract (Z(\omega)), then match to outer brane* — is the correct “missing layer.”

### You explicitly included the two “keep it open” forces in the *right* place

Putting both standing-wave stress *and* through-flow momentum flux into the **same stress/energy accounting** is exactly how you avoid hand-waving. That part is conceptually clean.

---

## The main thing I would change (important): microphysics consistency with the (n=5) vacuum

Right now the spec chooses a **4D Gross–Pitaevskii equation (GPE)** as the parent PDE. That’s a great PDE *shape* (conservative, has a clean energy functional, supports waves + nonlinear density-flow coupling), but vanilla cubic GPE corresponds to an effective equation of state closer to a **polytrope with (n=2)** behavior, not the **(n=5)** stiff polytrope that the 1PN/optics calibration in `1pn_hybrid.tex` treats as non-negotiable.

So you have two good options:

1. **Keep the GPE form but generalize the nonlinearity** so the pressure law matches the calibrated polytrope.

   * i.e., replace (g|\psi|^2\psi) with a power-law term that yields (P(\rho)\propto \rho^5) (or equivalently a chemical potential (\mu(\rho)) consistent with your chosen EOS).
   * You still get a Hamiltonian PDE + stress tensor + clean (E[\psi]).

2. **Use the hydrodynamic variables directly** (continuity + Euler with (P\propto \rho^5)) and add a controlled dispersive/regularizing term only as needed (your “quantum pressure analog”).

   * This may actually fit your “Level A: local throat solve + matching” goal better, because it speaks the same language as the 1PN papers.

If Paper 7’s throat PDE uses microphysics inconsistent with the calibrated vacuum, you’ll end up with a “beautiful throat” that’s not the same universe as Papers I–VI. This is the single biggest coherence risk.

---

## Second key tightening: make the “recover (L/a\simeq 1.85)” check non-optional and *BC-accurate*

Your ontology is explicit that the EM aspect ratio comes from a **Dirichlet-type radial root** (x_{01}) (first zero of (J_0)) and an axial standing wave, giving
[
\frac{L}{a}=\frac{\sqrt{2}\pi}{x_{01}}\approx 1.85.
]
So Paper 7 needs a clearly labeled **EM-cavity limit** subsection where:

* the mouth behaves like a **high-impedance termination** (closed-ish cavity limit),
* the **radial BC matches the EM branch** (Dirichlet if that’s what produced the (x_{01}) selection),
* and the same minimization constraint is stated (fixed “charge” / fixed (N=\int\rho) / whatever your EM paper used).

Otherwise you’ll keep seeing the (correct but confusing) phenomenon your DtN experiments already showed: **changing Neumann↔Dirichlet changes the radial root and the “selected” (\Lambda)** — and then it looks like the throat model “disagrees” with the series when it’s really just a different branch.

---

## Third tightening: separate “impedance extraction” from “geometry evolution” in the workflow

You *can* let (a(t),L(t)) breathe in the full nonlinear run, but in practice you’ll get a cleaner Paper 7 if you stage it:

1. **Fixed-geometry inner PDE → linearized response** → extract (Z_{\text{in}}(\omega)).
2. **Then** add geometry evolution using wall stress / shape derivatives:

   * standing wave stress on the wall,
   * plus through-flow momentum flux contribution,
   * balanced against (E_{\text{geom}}(a,L)).

This keeps the impedance story (the thing the brane needs) from getting tangled with moving-boundary numerics on day 1.

---

## Fourth tightening: be explicit about the conjugate “port variables” so matching can’t drift

Right now the spec is a bit abstract about (u) vs (j). I’d lock this down early:

* If the brane couples via **pressure/enthalpy**, then (u) should literally be a boundary enthalpy/pressure amplitude and (j) the **normal mass flux**.
* If the coupling is via **phase/chemical potential** (superfluid language), then (u) should be (\delta\mu) or (\delta\theta), and (j) the conjugate current.

This matters because “passivity” and “PN-locality expansion” depend on choosing the correct conjugate pair (and consistent normals/signs). Your existing matching code already bakes in these conventions (and includes locality/passivity sanity checks) — it’s a good unit-test scaffold to keep around.

---

## Fifth tightening: bake in the symmetry/overlap selection rule (otherwise you’ll match the wrong channels)

Your Step 4 experience is a real physics constraint: certain inner modes won’t couple to certain outer multipoles under a given mouth projection. That needs to be promoted from “numerical gotcha” to “model statement.”

So in Paper 7, explicitly include:

* a port basis at the mouth,
* the overlap integrals that tell you which ((m,n)) couples to which (\ell),
* and a statement like: *“uniform monopole mouth forcing does not excite modes with zero disk-average.”*

That prevents false “no matching exists” conclusions that are really “wrong channel.”

---

## One meta-point to state plainly in the paper: matching alone won’t always *select* (\Lambda) without an anchor

Your own roadmap already anticipates this, and your cylinder matching experiments basically proved it:

* In conservative/static matching, you can get **monotone drift** (best at scan boundary) unless something pins the relevant frequency content (drive spectrum, termination, or a fixed “distance above cutoff” anchor).
* That’s not a failure — it’s identifiability. The outer world must decide which mode family matters.

I’d actually treat this as a “strength” of Paper 7: it tells you exactly what extra physics input is required for (\Lambda) to become an inference rather than a free knob.

(Your “near–far matching layer” roadmap in `inner_throat.md` already frames this correctly — especially the idea that the *interface operator* is the real physical output, and (\Lambda) selection is regime/anchor-dependent.)

---

## If you want a crisp “Paper 7 commit message”

If you had to summarize the paper in one sentence that is consistent with the ontology + hybrid:

**Paper 7 builds the first self-consistent inner throat PDE whose linearized mouth DtN map (Z(\omega)) and geometry-dependent energy bookkeeping reproduce the EM cavity (\Lambda\simeq 1.85) limit, quantify ((a/r)^2) far-field leakage, and provide the impedance object needed to propagate strong-field corrections through the 1PN hybrid framework.**

That’s the “spine.” Everything in your spec supports it — with the main caveat being: **make the throat microphysics match the calibrated (n=5) vacuum (or explicitly justify why the throat is a different effective medium).**

If you want, next I can rewrite the spec into a tighter “Paper 7 checklist” where each milestone is explicitly keyed to: (i) the ontology claims (multipole suppression + EM aspect ratio), and (ii) the hybrid paper needs (impedance + (\Lambda)-dependent strong-field corrections).

===

# Paper 7 “Full Picture” Master Equation Sheet

I’ll write this in a way that is **complete**, but also makes it obvious what parts are already implemented vs what’s still missing.

## 1) Inside the throat (bulk): pick one of two equivalent “parent” forms

### 1A) Hydrodynamic (compressible, barotropic) form — most direct

On the 4D bulk domain (\mathcal{T}\subset\mathbb{R}^4) with (\mathbf{X}=(x,y,z,w)):

**Continuity**
[
\partial_t\rho + \nabla_4\cdot(\rho \mathbf{v})=0.
]

**Euler (inviscid barotrope + confinement)**
[
\partial_t \mathbf{v} + (\mathbf{v}\cdot\nabla_4)\mathbf{v}
= -\nabla_4 h(\rho) ;-; \nabla_4 V_{\rm wall}(\mathbf{X};a,L,\text{shape}) ;+; (\text{optional regularization}).
]

**EOS constraint from `1pn_hybrid.tex`**
[
P(\rho)=K\rho^n,\qquad n=5.
]

Then
[
c_s^2(\rho)=\frac{dP}{d\rho}=nK\rho^{n-1}=5K\rho^4,
\qquad
h(\rho)=\int^\rho \frac{dP}{\rho'}=\frac{nK}{n-1}\rho^{n-1}=\frac{5K}{4}\rho^4.
]

✅ This is the cleanest “full picture” interior PDE.
⚠️ Missing detail: what “optional regularization” term you want (if any) when gradients steepen near the wall/mouth.

---

### 1B) Generalized NLS (GNLS) form — same physics + built-in dispersive regularization

Let (\rho=|\psi|^2). Use
[
i\hbar \partial_t \psi
======================

\left[
-\frac{\hbar^2}{2m}\nabla_4^2
+V_{\rm wall}(\mathbf{X};a,L,\text{shape})
+h(|\psi|^2)
\right]\psi,
]
with
[
h(\rho)=\frac{5K}{4}\rho^4
\quad\Rightarrow\quad
h(|\psi|^2)=\frac{5K}{4}|\psi|^8.
]

This ensures the Madelung form reproduces the barotropic Euler equation with (P\propto\rho^5) and includes a quantum-pressure-like dispersive term automatically.

✅ This fixes the EOS mismatch with standard cubic GPE.
⚠️ Missing detail: you still need to specify the wall model (below) and interface conventions.

---

## 2) Throat wall / geometry model (this is the critical missing piece)

There are *three levels*. Paper 7 should explicitly choose one.

### 2.1 Rigid wall (what current DtN/cylinder formulas assume)

Geometry fixed. Impenetrable boundary on (\Sigma_{\rm wall}):
[
\mathbf{v}\cdot \hat{n}=0
]
(or equivalent Neumann condition on your chosen perturbation variable).

✅ Already used implicitly in the cylinder DtN work.
❌ Not sufficient to claim “stays open” via force balance.

---

### 2.2 Soft wall via confinement potential (fixed mesh, but “breathes” by parameters)

Use a smooth (V_{\rm wall}(\mathbf{X};a(t),L(t),\dots)). No hard boundary; instead the “radius” lives inside the potential parameters.

Then “breathing” becomes an **energy/force balance problem**:
[
E_{\rm total}(a,L) = E_{\rm fluid}[\psi;a,L] + E_{\rm geom}(a,L),
]
with equilibrium defined by
[
\frac{\partial E_{\rm total}}{\partial a}=0,\qquad
\frac{\partial E_{\rm total}}{\partial L}=0.
]

✅ Numerically friendly (no moving mesh).
⚠️ Missing: choose a concrete (V_{\rm wall}) family consistent with the “rounded funnel” ontology.

---

### 2.3 Explicit wall law (stress balance) — the “real” throat wall dynamics

Introduce a geometric energy (minimal useful form):
[
E_{\rm geom} = P_{\rm vac},V(a,L) + \sigma,A(a,L) + \kappa_b(\text{curvature term})+\cdots
]

Then impose **normal traction balance** at the wall:
[
\hat{n}\cdot \mathbf{T}_{\rm fluid}\cdot \hat{n}
================================================

P_{\rm vac} + \sigma,H + (\text{curvature/bending terms}),
]
where (in the inviscid limit) the fluid momentum flux is
[
\mathbf{T}_{\rm fluid} = \rho,\mathbf{v}\mathbf{v} + P(\rho),\mathbf{I} + (\text{GNLS dispersive stress if using GNLS}).
]

✅ This is exactly where “standing wave pressure” + “through-flow momentum flux” become *the* opening force.
❌ Missing right now: you have not yet fixed (E_{\rm geom}) (terms + interpretation) nor written the final traction condition in the same variable conventions as your solver.

**This wall law is the #1 derivation/spec item before “full throat modeling” is coherent.**

---

## 3) Just outside the throat (brane side): background + perturbations

### 3.1 Background inflow profile (already in `1pn_hybrid.tex`)

You already have the steady, spherically symmetric inflow model:

Mass flux:
[
\dot{M}_{\rm flux} = 4\pi r^2\rho(r)u(r).
]

Euler + EOS gives the standard transonic structure (your paper writes it as):
[
\big(u^2-c_s^2\big)\frac{du}{dr}
================================

u\left[
\frac{2c_s^2}{r}
-\frac{d\Phi}{dr}
+\mathcal{F}_{\rm throat}(r;a)
\right],
]
and the sonic point conditions at (r=r_H).

✅ Present in the hybrid paper.
⚠️ Missing for Paper 7 PDE coupling: **an explicit analytic form** for the throat correction (\mathcal{F}_{\rm throat}), not just “it’s encoded in a solver script.”

**Best “full picture” way to define (\mathcal{F}_{\rm throat}): area transition**
Treat the near-mouth region as quasi-1D flow with an effective area (A_{\rm eff}(r)) that transitions from (4\pi r^2) (3D) to the throat intake area (whatever your throat geometry implies). Then the nozzle equation gives
[
\big(u^2-c_s^2\big)\frac{du}{dr}
================================

u\left[
c_s^2\frac{d}{dr}\ln A_{\rm eff}(r)
-\frac{d\Phi}{dr}
\right].
]
If (A_{\rm eff}=4\pi r^2), you recover the (2c_s^2/r) term. Deviations from spherical area become the explicit “throat correction.”

➡️ **This is a key missing derivation/spec**: pick (A_{\rm eff}(r)) consistent with the brane→bulk mouth geometry in the ontology.

---

### 3.2 Perturbations “just outside” (near-field PDE, variable coefficients)

This is what you were worried was missing — and yes, it’s needed if we want scans to represent strong-field near-mouth physics.

For perturbations on a background ((\rho_0(r),\mathbf{v}_0(r),c_s(r))), the correct outer-near PDE is the acoustic-on-flow form (3D):
[
(\partial_t+\mathbf{v}_0\cdot\nabla_3)\left(\frac{\rho_0}{c_s^2}(\partial_t+\mathbf{v}_0\cdot\nabla_3)\phi_1\right)
-\nabla_3\cdot(\rho_0\nabla_3\phi_1)=0.
]

✅ Equation is known and consistent with your framework.
❌ Missing: an **overlap/handoff rule** that tells us where we can replace this with the far-field operator below.

---

## 4) Far outside (scan-ready simplification)

In the far field (or in a simplified constant-medium matching stage), use Helmholtz on the brane:
[
(\nabla_3^2+k^2)U=0,\qquad k=\omega/c,
]
with spherical multipole DtN eigenvalues (static limit (-(\ell+1)/a), radiating Hankel ratio in the dynamic case).

✅ This is what the current Python outer matching code uses.

---

## 5) Interface at the mouth (must be pinned down once)

Pick a conjugate pair ((u,j)) **consistent with the hybrid thermodynamics**. The safest choice is:

* **Port variable**: enthalpy perturbation (u \equiv h_1) (or pressure (p_1), but then be consistent everywhere)
* **Flux**: normal mass flux (j \equiv (\rho_0 v_{n,1}))

Interface conditions at (\Gamma_{\rm mouth}):
[
h_{1,\rm in}=h_{1,\rm out},
\qquad
(\rho_0 v_{n,1})*{\rm in} = (\rho_0 v*{n,1})_{\rm out}.
]

Then define the extracted inner operator (matrix if using ports):
[
j_i(\omega)=\sum_j Z^{\rm in}_{ij}(\omega),u_j(\omega).
]

✅ The concept is already in your spec and code structure.
❌ Missing: the **canonical** variable choice + sign conventions + mouth basis (P_i) (so passivity/locality checks don’t drift).

---

# What we’re missing (derivation/spec checklist)

This is the “stop getting ahead of ourselves” list — in priority order.

## Must-have next (before “meaningful scans” beyond the toy operator testbed)

1. **Wall law selection + definition**

* Decide: rigid BC vs impedance wall vs deformable stress-balance wall.
* If deformable: write (E_{\rm geom}(a,L)) explicitly and derive the traction condition (or reduced (\partial E/\partial a=0) equilibrium condition).
* Define what parameters mean physically (vacuum pressure term, surface term, curvature term).

2. **Explicit near-mouth 3D→4D transition model**

* Specify (A_{\rm eff}(r)) (or an equivalent geometric transition function) so (\mathcal{F}_{\rm throat}) is not a black box.
* This is what converts “outside just near the mouth” from narrative to equations.

3. **Outer-near perturbation PDE handoff**

* Define an overlap criterion or (r_{\rm match}) where variable-coefficient outer-near can be replaced by far-field multipole DtN.

4. **Port variable conventions**

* Lock: (u=h_1) vs (u=\phi_1) vs (u=p_1).
* Lock flux definition and sign conventions.
* Lock mouth basis (P_i) and the selection rules (what couples to what).

## Next tier (once the above is pinned)

5. **Interior “parent PDE” commitment**

* Choose hydrodynamic vs GNLS as your baseline.
* If GNLS: write the stress tensor you will use in the wall traction balance (including dispersive stress).
* If hydrodynamic: choose the minimal regularization/dissipation needed for stable numerics.

6. **Clarify the interior transverse eigenproblem**

* The EM (J_0(x_{01}r/a)) branch implies a specific transverse geometry/BC.
* Right now the ontology language has a small ambiguity (“(|\mathbf{x}|\lesssim a)” but “area (\pi a^2)”). Paper 7 should state clearly what the transverse domain is that produces the EM selector.

---

## How this answers your “are we ahead of ourselves?” feeling

* **We have enough equations to do operator scans** in the “cylinder inner + far-field outer” surrogate sense.
* But to claim scans reflect the *full throat physics*, we first need: **wall law + explicit 3D→4D transition + near-field handoff + fixed port conventions**.

===

# Paper 7 Master Equations

## 1) Geometry via potentials (no stitched dimensional interfaces)

Work on full 4D space with coordinates
[
\mathbf{X}=(x,y,z,w),\qquad \nabla_4=(\partial_x,\partial_y,\partial_z,\partial_w).
]

### Brane confinement (makes the world effectively 3D far away)

Introduce a stiff confining potential
[
V_{\rm brane}(w);;\text{with minimum at }w=0,\ \text{large for }|w|\gg \ell_w.
]
So outside the throat, the field sits in the (w)-ground state and physics becomes effectively 3D.

### Throat opening (bulk channel)

Add a throat “hole” potential that locally *reduces* confinement and creates a 4D corridor:
[
V_{\rm throat}(\mathbf{X};a,L,\text{shape}).
]
Total confinement:
[
V_{\rm conf}(\mathbf{X})=V_{\rm brane}(w)+V_{\rm throat}(\mathbf{X};a,L,\text{shape}).
]

* Deep throat: bulk accessible (weak confinement in (w))
* Far field: brane-trapped (strong confinement in (w))

This is the hard-mode way to make “brane + bulk” literally one PDE.

---

## 2) Bulk superfluid dynamics with the **non-negotiable** (n=5) EOS

You can run this as hydrodynamics or GNLS. For “correct + stable numerics,” GNLS is usually better.

### EOS / enthalpy (fixed)

[
P(\rho)=K\rho^5,\quad c_s^2(\rho)=\frac{dP}{d\rho}=5K\rho^4,\quad
h(\rho)=\int^\rho \frac{dP}{\rho'}=\frac{5K}{4}\rho^4.
]

### 2A) 4D GNLS (recommended backbone)

Let (\rho=|\psi|^2). Use
[
i\hbar \partial_t\psi
=====================

\left[
-\frac{\hbar^2}{2m}\nabla_4^2
+V_{\rm conf}(\mathbf{X})
+h(|\psi|^2)
\right]\psi,
\qquad
h(|\psi|^2)=\frac{5K}{4}|\psi|^8.
]

Energy functional:
[
E_{\rm fluid}[\psi]=\int_{\mathbb{R}^4}
\left[
\frac{\hbar^2}{2m}|\nabla_4\psi|^2
+V_{\rm conf}|\psi|^2
+\frac{K}{4}|\psi|^{10}
\right]d^4X.
]

### 2B) 4D Hydrodynamics (parallel track)

[
\partial_t\rho+\nabla_4\cdot(\rho\mathbf v)=0,
]
[
\partial_t\mathbf v+(\mathbf v\cdot\nabla_4)\mathbf v
=-\nabla_4 h(\rho)-\nabla_4 V_{\rm conf}(\mathbf{X})+\mathbf f_{\rm reg}.
]

---

## 3) “Two opening forces,” modeled as **two separate experiments**

You were explicit: model them separately and see where they go. So we set up two *distinct* wall-balance closures.

### 3A) Fluid-stress support (through-flow / momentum flux)

From hydrodynamics, the inviscid barotropic stress tensor is
[
\mathbf T_{\rm fluid}=\rho,\mathbf v\mathbf v+P(\rho)\mathbf I.
]
Normal traction (opening pressure) on a wall with unit normal (\hat n):
[
\Pi_{\rm fluid}=\hat n\cdot \mathbf T_{\rm fluid}\cdot \hat n
= P(\rho)+\rho(\mathbf v\cdot\hat n)^2.
]

**In GNLS form**, instead of wrestling with explicit “quantum stress,” you can get the wall force cleanly by varying the confinement potential (this is the most robust hard-mode tactic):
[
F_a^{\rm fluid}=-\frac{\partial E_{\rm fluid}}{\partial a}
\approx -\int |\psi|^2,\frac{\partial V_{\rm throat}}{\partial a},d^4X,
]
(and similarly for (L)). This is the “no ambiguity” route.

### 3B) Wave-stress support (standing wave), **as a separate field**

You want this independent. Do it explicitly.

Introduce a 4D wave field (A(\mathbf X,t)) (scalar surrogate for “EM-like” support is fine to start):
[
\partial_t^2 A - c^2\nabla_4^2 A = S_{\rm drive},
]
with boundary/confinement behavior aligned to the same geometry (or its own waveguide confinement).

Energy:
[
U_{\rm wave}=\int \frac12\left(\frac{(\partial_t A)^2}{c^2}+|\nabla_4 A|^2\right)d^4X.
]

Radiation force on the breathing DOF:
[
F_a^{\rm wave}=-\frac{\partial U_{\rm wave}}{\partial a},
\qquad
F_L^{\rm wave}=-\frac{\partial U_{\rm wave}}{\partial L}.
]

This keeps “wave support” cleanly separable from the fluid PDE.

> If later you want “phonon-like” wave support, do the same thing but set (c\to c_s(\rho_0)) and define the wave region/drive accordingly. It’s still a separate track.

---

## 4) Geometry energy and wall closure (the one thing that must be frozen)

Define
[
E_{\rm geom}(a,L)=P_{\rm vac}V(a,L)+\sigma A(a,L)+\kappa_b(\text{curvature})+\cdots
]

Then you have two separate closure experiments:

### Experiment 1: Fluid-only particle throat

[
\partial_a(E_{\rm fluid}+E_{\rm geom})=0,
\qquad
\partial_L(E_{\rm fluid}+E_{\rm geom})=0.
]

### Experiment 2: Wave-supported particle throat

[
\partial_a(E_{\rm fluid}+U_{\rm wave}+E_{\rm geom})=0,
\qquad
\partial_L(E_{\rm fluid}+U_{\rm wave}+E_{\rm geom})=0.
]

Dynamic variant (if you want):
[
M_a\ddot a+C_a\dot a+\partial_a E_{\rm geom}=F_a^{\rm fluid}+F_a^{\rm wave},
]
and same for (L).

---

## 5) How the “brane observables” come out of 4D (no cheating)

Because we’re not stitching a 3D PDE, define **brane-projected observables** by projecting onto the brane’s transverse ground state (\chi_0(w)) induced by (V_{\rm brane}).

Define the effective 3D brane field:
[
\psi_{\rm brane}(x,y,z,t)=\int \chi_0(w),\psi(x,y,z,w,t),dw.
]

Then define mouth “ports” using integrals over the **brane region** (3D volume patch) or the **brane mouth surface** (2D disk) *depending on what you want to measure*. The key is: it’s derived from the 4D solution, not imposed.

For your impedance-style object, a common choice is:

* effort variable: (u=\delta h) evaluated on/near the brane ((w\approx 0)),
* flux: normal mass flux across a chosen brane-mouth 2D disk (in ((x,y,z))-space) **or** leakage flux into bulk across (w) (if you want “loss into 4D” explicitly).

You can then define
[
j_i(\omega)=\sum_j Z^{\rm eff}_{ij}(\omega),u_j(\omega)
]
where (Z^{\rm eff}) is now an **emergent effective operator** from the full 4D solve.

This is the correct hard-mode replacement for “inner DtN.”

---

## 6) Your side note about black holes vs particles (and what Paper 7 should do)

Completely consistent:

* **Paper 7 (particles):** standing-wave support is plausible, so we run both tracks (fluid-only and wave-supported) and see what equilibria exist and what effective operators look like.
* **Later (black holes / 2PN relevance):** you likely dial wave-support → negligible and push into the **through-flow dominated** regime; then compare “how off” we are from PN/EFT expectations and adjust fundamentals if needed.

Paper 7 doesn’t need to “be 2PN-correct” as long as it produces **internally consistent 4D dynamics + a clean 3D-effective observable map**.

---

# The immediate next step (hard-mode, zero ambiguity)

Freeze these three things, and everything else becomes derivable:

1. **Choose explicit (V_{\rm brane}(w))** (sets (\chi_0(w)) and defines what “brane observable” means).
2. **Choose explicit (V_{\rm throat}(\mathbf X;a,L))** (sets what “throat geometry” means in 4D).
3. **Choose (E_{\rm geom}(a,L))** and the constraint (fixed (N=\int|\psi|^2 d^4X) vs fixed (m_{\rm eff})).

Once those are fixed, you can run:

* fluid-only equilibria,
* wave-supported equilibria,
* and then compute the emergent effective mouth operator from the 4D solution.

If you want, I’ll propose two concrete families for (V_{\rm brane}(w)) and (V_{\rm throat}) that (i) preserve the EM/cavity branch as a *limit*, but (ii) are genuinely 4D and amenable to PDE numerics without hidden boundary hacks.

===

Here are **two concrete, numerically friendly 4D potential families** that (i) give you a clean “brane emerges as a trapped (w)-ground-state” far away, (ii) create a **true 4D bulk corridor** (the throat) without stitching separate PDEs, and (iii) preserve the **EM/cavity (\Lambda=L/a\simeq 1.85)** branch as a controllable limit.

I’ll write each as a **complete (V_{\rm conf}(\mathbf X;a,L,\text{shape}))** with a clear split:
[
V_{\rm conf}=V_{\rm brane}(w;\text{local trap strength});+;V_{\rm wall}(\text{tube shape}).
]

Throughout, let
[
r_\perp=\sqrt{x^2+y^2},\qquad z\ \text{axial},\qquad w\ \text{bulk coordinate}.
]
And use a smooth step (so forces are differentiable):
[
S(u)\equiv \tfrac12\big(1+\tanh u\big),\quad S(u)\approx 0\ (u\ll 0),\ S(u)\approx 1\ (u\gg 0).
]

---

## Family 1: Modulated (w)-trap + soft cylindrical wall (best for “EM branch survives”)

### 1A) Brane confinement (harmonic, but *locally weakened* inside the throat)

Take a harmonic trap in (w) whose stiffness depends on ((r_\perp,z)):
[
V_{\rm brane}(w;r_\perp,z)=\frac12,m,\Omega_w^2(r_\perp,z),w^2.
]

Define a smooth “inside-throat” gate (G(r_\perp,z)\in[0,1]) (1 inside, 0 outside):
[
G(r_\perp,z)=\Big[1-S!\Big(\frac{r_\perp-a(z)}{\delta_r}\Big)\Big];
\Big[1-S!\Big(\frac{|z|-L/2}{\delta_z}\Big)\Big].
]

Then interpolate trap frequency:
[
\Omega_w^2(r_\perp,z)=\Omega_{\rm out}^2-\big(\Omega_{\rm out}^2-\Omega_{\rm in}^2\big),G(r_\perp,z),
]
with (\Omega_{\rm out}\gg \Omega_{\rm in}\ge 0).

**Meaning:**

* Far from the throat: (\Omega_w\approx \Omega_{\rm out}) → strong brane trapping → effectively 3D.
* Deep throat: (\Omega_w\approx \Omega_{\rm in}) → bulk (w)-excursions become possible → genuinely 4D dynamics.

> **EM/cavity limit knob:** set (\Omega_{\rm in}\approx \Omega_{\rm out}) (strong (w)-confinement everywhere). Then the throat’s transverse structure is essentially in ((x,y)) only, and your (J_0(x_{01} r_\perp/a)) branch is cleanly recovered.

### 1B) Throat wall (soft “hard-wall” in (r_\perp), plus endcaps in (z))

Use a steep soft wall:
[
V_{\rm wall}(r_\perp,z;a,L)=V_0,S!\Big(\frac{r_\perp-a(z)}{\delta_r}\Big)^{p}
;+;
V_0,S!\Big(\frac{|z|-L/2}{\delta_z}\Big)^{p}.
]

A simple flared mouth profile:
[
a(z)=a_0\Big(1+\beta,e^{-z^2/z_m^2}\Big),
]
so the tube widens near (z=0) (the brane neighborhood), and approaches (a_0) deep inside.

### 1C) Total confinement

[
V_{\rm conf}(\mathbf X)=V_{\rm brane}(w;r_\perp,z)+V_{\rm wall}(r_\perp,z;a,L).
]

**Why this family is excellent**

* Keeps the **cylindrical EM selector** intact as a limit (because the wall depends on (r_\perp) only).
* Gives a controlled “how 4D is it?” dial: (\Omega_{\rm in}/\Omega_{\rm out}).
* Makes wall forces extremely clean via Hellmann–Feynman (you’ll use (\partial_a V_{\rm wall}), (\partial_L V_{\rm wall}), (\partial_a \Omega_w), etc., inside volume integrals).

---

## Family 2: True 4D tube in ((r_\perp,w)) with an “ellipticity” dial (best for “4D cross-section physics”)

This one bakes 4D geometry *directly* into the tube cross-section, while still letting you recover the EM branch by taking a parameter limit.

### 2A) Define a 4D transverse radius with an anisotropy parameter (\gamma)

[
R_\gamma \equiv \sqrt{r_\perp^2+\left(\frac{w}{\gamma}\right)^2 }.
]

* If (\gamma\to\infty): (R_\gamma\to r_\perp) → cross-section ignores (w) → **EM cylinder limit**.
* If (\gamma\sim 1): (w) contributes comparably → true 4D transverse structure.

### 2B) Throat wall as a soft wall in (R_\gamma) + axial endcaps

[
V_{\rm wall}(r_\perp,w,z;a,L)=
V_0,S!\Big(\frac{R_\gamma-a(z)}{\delta_\perp}\Big)^p
;+;
V_0,S!\Big(\frac{|z|-L/2}{\delta_z}\Big)^p.
]

Same flared mouth (a(z)) as above works well.

### 2C) Brane confinement as a *separate* stabilizer (recommended)

Even though the tube already confines in (w), you still want a brane-trap outside the mouth to ensure 3D behavior far away. Use a Pöschl–Teller (nice smooth tails, controllable ground state width):
[
V_{\rm brane}(w)=V_b,\tanh^2!\Big(\frac{w}{\ell_w}\Big)
\quad\text{(≈ quadratic near 0, saturates at }V_b\text{)}.
]

Then reduce this brane trap inside the throat with the same gate (G(r_\perp,z)):
[
V_{\rm brane}(w;r_\perp,z)=\big(1-\eta,G(r_\perp,z)\big),V_b,\tanh^2!\Big(\frac{w}{\ell_w}\Big),
\quad 0\le \eta\le 1.
]

### 2D) Total confinement

[
V_{\rm conf}(\mathbf X)=V_{\rm brane}(w;r_\perp,z)+V_{\rm wall}(r_\perp,w,z;a,L).
]

**Why this family is useful**

* It forces you to confront truly 4D transverse mode structure (good for your “sporadic 3D behavior from 4D interactions” goal).
* It has a clean “how much does (w) count as transverse geometry?” dial: (\gamma).
* It still has a clean EM-limit: (\gamma\to\infty) *and* strong brane trap ((\eta\to 0)).

---

## How these preserve (and test) the (\Lambda\simeq 1.85) cavity branch

To intentionally recover the EM/cavity selector, run a “cavity limit” where:

* (w) is strongly confined everywhere (Family 1: (\Omega_{\rm in}\approx\Omega_{\rm out}) large; Family 2: (\eta\approx 0) and/or (\gamma\to\infty)),
* walls are steep (large (V_0), small (\delta), moderate-to-large (p)),
* endcaps are reflective (your (z)-wall term does this smoothly).

Then the dominant transverse eigenfunction in the throat approaches the usual (J_0(x_{01}r_\perp/a)) form, and your (\Lambda) selector test is meaningful again.

---

## Quick guidance on which to start with

* **Start with Family 1** if you want: *“brane is 3D, throat becomes 4D only where needed”* and you care about retaining the EM selector as a crisp unit test.
* **Start with Family 2** if you want: *“the throat is geometrically 4D in its cross-section”* and you want to explore whether extra transverse structure creates new, possibly intermittent coupling channels.

If you want, next I can propose **one concrete “default parameter set”** (dimensionless) for each family that is numerically stable (smooth but stiff walls), and I’ll write down the **exact derivatives (\partial_a V_{\rm conf}), (\partial_L V_{\rm conf})** you’ll use to compute breathing forces without ambiguity.

===

Start by freezing the **one-page “closed system” of equations** that everything else is derived from. If we get that right, the rest (DtN/impedance extraction, stability, parameter scans, future MHD extensions) becomes *mechanical*.

Here’s the clean place to start, in hard-mode **fully 4D** and with **two support mechanisms as two separate experiments**.

---

## 0) First principle: make it a single 4D variational system

Write one total energy/Hamiltonian that depends on:

* superfluid field (\psi(\mathbf X,t)) in 4D,
* optional independent standing-wave field (A(\mathbf X,t)) in 4D,
* geometry DOFs ((a,L)) (or fields (a(\cdot)) later),
* a confining potential (V_{\rm conf}(\mathbf X;a,L)),
* a geometry energy (E_{\rm geom}(a,L)).

Then *everything* is derived by variation:
[
i\hbar,\partial_t\psi=\frac{\delta H}{\delta \psi^*},\qquad
\partial_t^2 A = -,\frac{\delta H}{\delta A},\qquad
M_a \ddot a + C_a \dot a = -\frac{\partial H}{\partial a}\quad(\text{or }\partial_a H=0).
]

This is the “hammer the math” foundation.

---

## 1) Freeze the geometry *as a potential*, not a stitched boundary

Pick one of your two families and **write it explicitly** as:
[
V_{\rm conf}(\mathbf X;a,L)=V_{\rm brane}(w;\text{trap}) ;+; V_{\rm throat}(\mathbf X;a,L;\text{shape}).
]

This avoids any 3D/4D seam. The “mouth” is just the region where the potential transitions.

**Deliverable here:** an explicit (V_{\rm conf}) plus its derivatives (\partial_a V_{\rm conf}), (\partial_L V_{\rm conf}) (because those are the breathing forces).

---

## 2) Freeze the superfluid PDE (n=5 is already fixed)

Use the hard constraint:
[
P(\rho)=K\rho^5,\quad
c_s^2(\rho)=5K\rho^4,\quad
h(\rho)=\frac{5K}{4}\rho^4.
]

### GNLS backbone (recommended)

Let (\rho=|\psi|^2). Define the fluid Hamiltonian:
[
H_{\rm fluid}[\psi;a,L]=\int_{\mathbb R^4}
\Big[
\frac{\hbar^2}{2m}|\nabla_4\psi|^2
+V_{\rm conf}(\mathbf X;a,L)|\psi|^2
+\frac{K}{4}|\psi|^{10}
\Big],d^4X.
]

Then the equation of motion is:
[
i\hbar,\partial_t\psi
=====================

\left[
-\frac{\hbar^2}{2m}\nabla_4^2
+V_{\rm conf}(\mathbf X;a,L)
+\frac{5K}{4}|\psi|^8
\right]\psi.
]

**This is your “fluid + Bernoulli + compressibility + regularizer” all in one.**

---

## 3) Freeze the *independent* wave-support PDE (the second experiment)

If you want wave stress as a separate mechanism, keep it as a separate field with its own stored energy.

A clean scalar surrogate:
[
\partial_t^2 A - c^2 \nabla_4^2 A + \mu_A^2(\mathbf X;a,L),A = S_{\rm drive}.
]

Wave energy:
[
H_{\rm wave}[A;a,L]=\int_{\mathbb R^4}\frac12\Big(\frac{(\partial_t A)^2}{c^2}+|\nabla_4 A|^2+\mu_A^2 A^2\Big),d^4X.
]

Where (\mu_A(\mathbf X;a,L)) is how you “confine” the wave to the throat (or make it reflective) using the same geometric parameters.

---

## 4) Freeze the wall / geometry energy (closure law)

Pick:
[
E_{\rm geom}(a,L)=P_{\rm vac},V(a,L)+\sigma,A(a,L)+\kappa_b(\text{bend})+\cdots
]
with explicit formulas for (V(a,L)) and (A(a,L)) **in your chosen geometric convention**.

Then the **two separate support experiments** are literally:

### Experiment F (fluid-only support)

[
H_{\rm tot}=H_{\rm fluid}+E_{\rm geom}
\quad\Rightarrow\quad
\partial_a H_{\rm tot}=0,;\partial_L H_{\rm tot}=0
\quad(\text{or dynamics}).
]

### Experiment W (wave + fluid support)

[
H_{\rm tot}=H_{\rm fluid}+H_{\rm wave}+E_{\rm geom}
\quad\Rightarrow\quad
\partial_a H_{\rm tot}=0,;\partial_L H_{\rm tot}=0.
]

**Key win:** you never need to hand-write stress tensors to get forces. You can compute
[
-\frac{\partial H_{\rm fluid}}{\partial a}
==========================================

-\int |\psi|^2,\frac{\partial V_{\rm conf}}{\partial a},d^4X,
\qquad
-\frac{\partial H_{\rm wave}}{\partial a}
=========================================

-\int \frac12,A^2,\frac{\partial \mu_A^2}{\partial a},d^4X
]
(and similarly for (L)), which is unambiguous and numerically stable.

---

## 5) Only after 0–4 are frozen: define the 4D→brane “observable map” and the mouth operator

Because we’re hard-mode, you don’t “set a 3D boundary condition.” You **define what a 3D observer measures** as a projection near (w=0). For example:

* brane density:
  [
  \rho_{\rm brane}(x,y,z,t)=\int W(w),|\psi(x,y,z,w,t)|^2,dw
  ]
  with a fixed weight (W(w)) concentrated near (w=0) (or the ground-state (|\chi_0|^2) if you want it).

* effort variable (enthalpy perturbation):
  [
  u=\delta h(\rho_{\rm brane})\approx 5K\rho_{\rm brane,0}^3,\delta\rho_{\rm brane}.
  ]

* flux variable: choose a specific brane-measurement surface (e.g., a 2D disk (\mathcal D) around the mouth in ((x,y,z))) and define
  [
  j=\int_{\mathcal D} \rho_{\rm brane},v_n,dS
  ]
  (or define a set of ports (P_i) and project).

Then define the emergent response:
[
j_i(\omega)=\sum_j Z^{\rm eff}_{ij}(\omega),u_j(\omega).
]

That’s your “mouth operator,” but now it’s correctly derived from a unified 4D model.

---

# So: where do we start?

Start by writing (in the paper) and freezing these **five blocks** in this order:

1. **Exact (V_{\rm conf}(\mathbf X;a,L))** (pick Family 1 or 2 as baseline, and write it explicitly).
2. **Exact (E_{\rm geom}(a,L))** and the convention for (V(a,L),A(a,L)).
3. **4D GNLS with (n=5)** (that part is already basically done).
4. **Independent wave field model** (A) (also explicit).
5. **Variational wall law**: equilibrium ((\partial_a H=0)) or dynamics ((M_a\ddot a=-\partial_a H)).

Once those are frozen, *then* we define the brane observable map and (Z^{\rm eff}(\omega)).

If you tell me which baseline you want as the “primary” (V_{\rm conf}) (Family 1 vs Family 2), I’ll write the **exact full set of equations as a single cohesive “Master Model” section** (with all definitions, the total Hamiltonian, and the derived EOMs), so you can drop it straight into Paper 7.

===
