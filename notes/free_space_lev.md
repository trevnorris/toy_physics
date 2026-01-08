## 0. Executive summary and scope

This note summarizes our current “free-space levitation” exploration inside the superfluid toy-physics framework developed across your project papers (EM mapping, brane/bulk ontology, and the 1PN/GR-matching program). The aim is not to claim a loophole in conservation laws, but to identify—using your model’s medium picture—what *would have to be true* for practical, rail-less hovering of macroscopic payloads (e.g., cargo transport) to be possible, and to derive the core scaling relations and numerical targets that decide whether the idea is dead-on-arrival.

A central distinction emerged early and stayed decisive: **support** versus **propulsion**. If a craft is supported by a nearby boundary (ground, plate, engineered layer), it can in principle hover using a static stress/traction distribution in the medium; this is a “ground-effect” mechanism, and mathematically it reduces to field-stress forces and image-method logic. However, in free space there is no boundary to push against. In that case, a sustained hover force (F=mg) must arise from exporting momentum into the surrounding medium—i.e., *propulsion by momentum flux*. In continuum language, the net force on a body equals the momentum flux of fields/excitations through a surface enclosing it, so a self-contained, static near-field configuration cannot provide net lift. Free-space hover therefore requires a “momentum sink”: either the craft radiates momentum away (waves/wakes) or ejects coherent momentum-carrying structures (the medium’s analog of a propeller pushing air).

This reframe makes a first, sharp viability test possible. For any momentum carrier with energy–momentum ratio (\varepsilon = E/P), the hover power is
[
\mathcal P = F,\varepsilon \simeq (mg),\varepsilon.
]
If the only available carriers behave like relativistic waves with (\varepsilon \sim c), then (\mathcal P \sim mgc), which is (\mathcal O(10^{12},\mathrm{W})) per ton—immediately impractical. The only path to engineering-feasible free-space hover is therefore to identify a momentum carrier with **(\varepsilon \ll c)**. In ordinary aerodynamics, the “carrier” is bulk air accelerated to modest velocities; the analog in a superfluid medium is a slow, momentum-rich excitation.

Within the toy model’s medium language, the most promising candidate for such a slow carrier is **vortex-ring–like structures** in an appropriate transverse sector: coherent, topological momentum packets whose energy-per-momentum is set by their translation speed (U), not by (c). Under standard vortex-ring scalings, (E/P \sim U), so hover power becomes (\mathcal P \sim (mg)U). This makes the engineering question concrete: for power budgets of (\sim 0.1)–(1) MW per ton, the effective “downwash” speed must be in the range (U\sim 10)–(100) m/s. The remaining unknowns shift from “can we cancel gravity?” to “can the medium support and couple to such momentum packets at the required rates without catastrophic losses or obvious electromagnetic signatures?”

In short, the scope of the current write-up is:

1. To document the key reframing from “gravity cancellation” to “vacuum aerodynamics / momentum flux,”
2. To record the core scaling laws and numerical targets that any free-space hover mechanism must satisfy, and
3. To map these targets into your model’s variables (circulation (\Gamma), core size (a), effective density (\rho_{\rm eff}), and defect-density scale (\rho_0)), so future sessions can test consistency against the ontology and GR-matching constraints.

We do not assume the existence of exotic energy sources or violations of conservation laws. We do assume, as a working hypothesis for this line of investigation, that the toy model’s vacuum medium may admit slow, coherent momentum carriers (vortex rings or analogous structures) and that engineered devices might couple to them—exactly the kind of “unnatural coupler” step that historically turned “impossible” into “engineering.”

---

## 1. Wright brothers framing: the “propeller” analogy for a superfluid vacuum

A useful historical analogy is not “people thought flying violated physics,” but rather: people lacked the *right actuation mechanism and the right mental model of the medium*. By the late 19th century much of the relevant physics (fluid pressure, lift, drag) existed in pieces, but practical flight required an engineering discovery: **a controllable way to exchange momentum with the air**. Birds provided the empirical clue—air can supply support forces if you shape flow and pressure—but the enabling technology was not “copy birds exactly.” It was the realization that a rotating blade can generate sustained circulation and accelerate a large mass of air downward at modest speed: the **propeller**. Propellers do not occur naturally in that form; they are a purpose-built coupler to a medium.

That is the right mindset for your superfluid-vacuum ontology. If the universe behaves like a medium with longitudinal (gravity-linked) and transverse (wake/vector-linked) degrees of freedom, then “hover” is not a mysterious negation of gravity. It is an engineering problem: **can we impose a controlled momentum flux in the medium** analogous to a rotor downwash in air?

This immediately clarifies what “new physics” would and would not mean in this program.

* The *math* of forces from stresses is not new. In any field/continuum theory, forces arise from gradients of energy density and from momentum fluxes in a stress tensor. If you can establish an asymmetric stress distribution, you can get traction and lift.

* The potential novelty is in the **coupling mechanism**: how a macroscopic device interacts with the vacuum medium’s available momentum carriers. If the only carrier is photon-like radiation, then the energy-per-momentum is too large ((\varepsilon \sim c)), and free-space hover becomes energetically hopeless. So the “Wright-brothers move” is to identify (or engineer access to) a different carrier—one that can transport momentum with (\varepsilon) set by a low speed (U), not by (c).

In your model language, this becomes a crisp statement: free-space hover requires that the vacuum medium admit a slow momentum carrier—e.g., vortex-ring–like packets in a transverse sector—and that matter can be engineered to couple to that sector strongly enough to generate and shed those packets efficiently. This is analogous to how a rotor couples to a compressible fluid: it does not “cancel gravity”; it produces sustained thrust by exporting momentum into the environment.

This framing also explains why earlier levitation pathways felt underwhelming. Boundary-supported levitation (ground-effect, image forces) can be powerful and efficient, but it is fundamentally “push on the boundary,” which makes it conceptually close to maglev. The unique value proposition of a superfluid vacuum, if it exists, is not rail-less ground levitation; it is **free-space momentum exchange** without needing air, rails, or reaction mass in the usual sense—i.e., a “vacuum propeller.”

So the Wright-brothers framing gives us the proper target for the rest of the analysis:

1. Identify the medium’s plausible momentum carriers and their (\varepsilon = E/P).
2. Determine whether any carrier has (\varepsilon) in the tens of m/s range (rather than (c)).
3. Translate that into concrete requirements in your toy-model parameters (circulation scales, core sizes, effective density, emission rates).
4. Only then ask whether the ontology and GR-matching constraints permit those requirements without contradicting obvious observations.

This provides a disciplined path: the “propeller” is not assumed to exist; it is parameterized, and the math tells us what properties it must have to make free-space hover plausible.

---

## 2. Physical objective in the toy model

### 2.1 Gravity as medium flow into sinks

Across the papers, the working picture is that what we call “gravity” is not a fundamental geometric axiom but an emergent dynamics of a medium (the superfluid) whose defects behave as sinks/sources and whose long-distance behavior can be matched to GR at 1PN/EIH order. In that language, a massive object is associated with a defect structure that sets up an inward-directed longitudinal flow pattern (or an equivalent scalar potential (\Phi)), and test bodies accelerate according to that effective potential:
[
\ddot{\mathbf x}\approx -\nabla\Phi \quad (\text{Newtonian limit}).
]
The “downward pull” near Earth is then a steady background acceleration field that, in medium terms, corresponds to a persistent flow/potential gradient sourced by Earth’s defect content.

### 2.2 What we are trying to achieve physically

The practical goal is not abstract “antigravity,” but a specific mechanical effect:

> Maintain a macroscopic payload at rest relative to Earth (hover) without physical contact and without relying on a supporting structure beneath it.

In ordinary engineering terms, hovering requires a steady upward force balancing weight:
[
F_\uparrow = mg.
]

In medium terms, this means we must create an *upward* momentum transfer to the craft from the medium (or an *equivalent* downward momentum transfer from the craft into the medium). There are only two broad physical ways to do that:

1. **Support by stress/traction**: Create a pressure/stress difference across the craft (e.g., lower pressure above, higher below), producing net upward traction. This is “lift by stress distribution.”

2. **Propulsion by momentum flux**: Continuously export momentum downward into the medium (waves, wakes, coherent structures), producing upward thrust as reaction.

Both are legitimate within a continuum framework. The crucial difference is the need for an external momentum partner (discussed in Section 3). A stress/traction method can work efficiently when there is a boundary or “reaction surface” nearby that participates in the stress field. A propulsion method can work in free space but must pay whatever energy cost is required to export momentum.

### 2.3 Why “countering the inflow” is deceptively expensive

A recurring intuition was: if gravity is inflow toward sinks, perhaps we can “impede” or “redirect” that inflow under a craft to reduce its effective weight. The model makes this idea easy to visualize, but the math reveals why it tends to be costly if interpreted as a *static cancellation* of the background gravitational field.

To reduce apparent weight in a static way, you must either:

* create an opposing field (an “upward (\nabla\Phi)”) comparable to Earth’s (g), or
* drastically alter the coupling between the craft’s defect content and (\Phi) (effectively changing gravitational charge).

In a Poisson-like theory, producing a local opposing field is equivalent to introducing huge effective sources/polarization, because Earth’s (\Phi/c^2) is tiny in ordinary weak-field conditions. This is why attempts to “cancel gravity directly” tend to resemble “mountain-equivalent” source strengths. In your GR-matching corner (1PN/EIH), the scalar sector is intentionally rigid, so these knobs are not freely adjustable without leaving the GR-matched regime.

This realization is what pushed us toward a Wright-brothers style reformulation: don’t try to statically cancel the Earth field; instead, identify an efficient medium coupling that yields thrust/lift via momentum exchange—analogous to how flight uses air as the momentum partner rather than trying to “cancel gravity.”

### 2.4 Candidate targets for practical utility

For cargo transport, the objective can be made quantitative with plausible engineering targets:

* Payload scale: (m \sim 10^3)–(10^5) kg (tons to tens of tons).
* Hover endurance: minutes to hours.
* Power budget: ideally (\lesssim 0.1)–(1) MW per ton (already very aggressive, but at least not astronomical).

These targets provide concrete benchmarks for deciding whether a proposed coupling mechanism has “legs.”

---

## 3. The key constraint: momentum conservation and the need for a momentum partner

### 3.1 Net force equals momentum flux through a control surface

In any continuum field theory (fluids, EM, elastic media), the net force on a body can be written as the surface integral of a stress tensor over a closed surface (S) surrounding the body:
[
F_i ;=; \oint_S T_{ij},n_j,dA,
]
where (T_{ij}) is the momentum-flux (stress) tensor, (n_j) is the outward normal, and repeated indices are summed.

This expression is the cleanest way to say: **you only get a net force if momentum flows across the boundary of your chosen control volume.**

### 3.2 Why self-contained static configurations don’t hover

If all fields/excitations are confined to the craft and decay symmetrically, the stress integral over a large enclosing surface cancels and the net external force is zero. Internal stresses can redistribute forces within the craft, but they cannot accelerate the craft’s center of mass without exchanging momentum with something outside the control volume.

This is not a special assumption about your toy model; it is a generic consequence of momentum conservation applied to a medium-plus-device system.

Practically, it means any proposal of the form “build a device that creates a static stress pattern around itself and thereby floats” will fail unless that stress pattern involves an external momentum partner.

### 3.3 The three allowable momentum partners

For sustained hover at rest relative to Earth, there are only three ways to supply the missing momentum partner:

**(A) A nearby boundary/structure (“ground effect”)**
If the craft’s field/stress configuration couples strongly to the ground (or an engineered layer), the boundary experiences the opposite force and supplies the reaction. This can be extremely efficient because it can be quasi-static. But it is inherently not free-space hover; it is a boundary-supported technology.

**(B) Momentum export into the medium (“propulsive hover”)**
If the craft emits waves/wakes/coherent packets downward, those carry momentum away and the craft experiences upward reaction thrust. This is free-space compatible, but it incurs a power cost set by the energy-per-momentum (\varepsilon = E/P) of the emitted carrier:
[
\mathcal P = F,\varepsilon \simeq (mg),\varepsilon.
]
If (\varepsilon\sim c), free-space hover is dead on arrival (terawatts per ton). The only viable path is (\varepsilon \ll c), i.e., a slow momentum carrier.

**(C) Long-range coupling to distant masses via a field**
Gravity itself is a long-range interaction: Earth and craft exchange momentum through the field. In principle, one could imagine engineering additional long-range fields that couple to Earth and supply lift. In practice, within the GR-matching weak-field regime, creating an additional long-range field comparable to Earth’s (g) looks equivalent to enormous effective source strengths unless one changes the coupling laws. This makes (C) the hardest path unless the toy model admits a new sector with unusually strong coupling at engineering scales.

### 3.4 Immediate consequence for our program

This section is the fulcrum of the whole levitation discussion:

* **If we require free-space hover**, we must pursue option (B): exporting momentum into the medium.
* The “Wright brothers move” is then to identify a momentum carrier with small (\varepsilon) (like air downwash), rather than radiation-like carriers with (\varepsilon\sim c).
* This is what motivated focusing on vortex-ring–like structures: in classic superfluid/hydrodynamics, vortex rings can carry large momentum for relatively modest energy, with (\varepsilon) set by their translation speed (U), not by (c).

Everything that follows in the technical sections is essentially an exploration of whether the toy model can support such carriers and whether their required parameters (speed, circulation, density, emission rate) land in an engineering-feasible regime.

---

## 4. What we learned from the EM and ontology papers

This section collects the “dictionary” pieces from your project that matter for levitation/propulsion, and explains why they point to a two-sector view (longitudinal gravity-like vs transverse wake-like) that can be exploited in different ways.

### 4.1 EM-paper dictionary: potentials as medium variables

In the EM mapping, the electromagnetic-like potentials are not taken as fundamental fields; they are **re-parameterizations of medium variables**. The key identifications we relied on are:

* A scalar potential built from fluid enthalpy (h) and kinetic energy density:
  [
  \phi_{\rm EM} \equiv \lambda\Big(h+\tfrac12 v^2\Big),
  ]
* A vector potential proportional to a flow-like field:
  [
  \mathbf A \equiv \lambda,\mathbf v.
  ]

With these, the “magnetic field” is automatically a vorticity-like quantity:
[
\mathbf B \equiv \nabla\times \mathbf A = \lambda,(\nabla\times\mathbf v)=\lambda,\boldsymbol\omega.
]

**Implication for levitation:** any EM-style force can be re-read as a **medium stress** force (pressure/traction from energy density gradients), and “magnetic structure” corresponds to structured circulation/vorticity in the medium variables.

### 4.2 Ontology-paper split: longitudinal bulk vs transverse brane wake

A crucial additional ingredient in the ontology paper is the explicit **decomposition of the medium response**:

* A bulk/longitudinal (potential) sector that is essentially irrotational at large scales and supports the gravity-like scalar potential (\Phi).
* A brane-confined **transverse wake sector** (\mathbf v_T), divergence-free on the brane, introduced precisely because the bulk being irrotational would otherwise kill vorticity and (naively) kill any “magnetic” structure.

Formally:
[
\mathbf v_{\rm 3D} = \nabla\Phi + \mathbf v_T,\qquad \nabla\cdot \mathbf v_T = 0,
]
and the ontology ties the vector potential to the transverse wake:
[
\mathbf A \equiv \kappa_A,\mathbf v_T,\qquad \mathbf B=\nabla\times\mathbf A.
]

The sourcing equation for (\mathbf A) takes a Poisson-like (magnetostatic-like) form:
[
\nabla^2 \mathbf A \sim -(\text{effective wake current}),
]
with a moving “charged” defect acting like a current element. We used this structure repeatedly because it means (at least in a quasi-static regime) the transverse sector behaves like a genuine stress/momentum carrier with its own Green’s function and boundary sensitivity.

**Implication for levitation:** if a device can excite and shape (\mathbf v_T) (or (\mathbf A)), it can generate stresses/tractions without necessarily needing to directly modify the gravity-like scalar (\Phi). That is the “vacuum aerodynamics” handle: use a momentum/stress-carrying transverse sector rather than trying to cancel the longitudinal gravity field.

### 4.3 Why “not just EM” is about coupling, not math

Once you adopt the stress-tensor viewpoint, the mathematics of traction is generic. What can be new is **how matter couples to the transverse wake sector**:

* If (\mathbf A) is literally the ordinary photon field, then your model is mostly a reinterpretation, and boundary reflection/levitation requires superconductors or active field control (maglev-like).
* If (\mathbf A\propto \mathbf v_T) is a distinct brane wake mode that *only looks Maxwell-like in vacuum*, then engineered materials might couple strongly to (\mathbf v_T) (pin it, screen it, reflect it) without being electrically superconducting. That would be a genuinely new constitutive possibility.

This is the lever we identified as the only plausible source of “Wright brothers–style novelty”: not a new conservation-violating force, but an engineered coupling to a medium mode that nature doesn’t hand us in a convenient form.

### 4.4 Where 1PN/hybrid matching constrains us

The 1PN and hybrid papers were written to reproduce GR/EIH behavior in the weak-field two-body regime. In that “GR-matching corner,” the scalar sector is intentionally rigid: masses, potentials, and optical mapping are tuned so that in ordinary weak fields you don’t have arbitrary freedom to dial couplings. That is good for consistency, but it also means:

* “reduce gravity by tweaking the scalar sector” is not a free knob near Earth, and
* any practical handle is more likely to live in **boundary conditions, constitutive response in matter, or non-conservative momentum exchange**—all of which are largely invisible to the original EIH matching program.

This motivates the path we took: explore boundary-supported stress first (as a sanity check), then pivot to free-space hover as a propulsion-by-momentum-flux problem.

---

## 5. Early levitation route and why we pivoted away

We started with two intuitive routes: (i) magnetic/diamagnetic levitation as a known existence proof, and (ii) “redirect gravitational inflow” as a medium picture. Both were useful, but both hit hard limits that forced a reframing.

### 5.1 Diamagnetic (frog) levitation: what it teaches and what it doesn’t

The frog-in-a-magnet experiment shows that levitation is possible in principle via electromagnetic fields: diamagnetic materials experience a force in a field gradient. The schematic scaling is:

* Diamagnetic force scales with susceptibility (\chi) and field-gradient strength, roughly like
  [
  F \sim \frac{\chi V}{\mu_0},\nabla!\left(\frac{B^2}{2}\right),
  ]
  where (\chi) is tiny for most materials.

**Lesson:** Levitation via field stresses exists.

**Problem:** Because (\chi) is so small, you need extremely large (B) and especially extremely large (\nabla(B^2)) to get macroscopic lift. That’s why frog levitation requires a very powerful magnet and is not an obvious path to general cargo transport.

This pushed us to ask a different question: not “can matter be nudged by EM gradients,” but “can we create a medium stress configuration that directly produces traction, like a pressure difference, without relying on tiny susceptibilities?”

That question led to the wake-stress / boundary-image approach (and later the free-space momentum-flux approach).

### 5.2 “Redirect inflow” / gravity-diode idea: why static gravity cancellation is expensive

In the medium picture, it’s tempting to think: gravity is flow into sinks; if we can impede that flow from below a craft, maybe it becomes weightless.

Mathematically, the analog is to treat the gravitational potential as a field satisfying a Poisson-like equation and to attempt to modify it by an “effective dielectric” (constitutive) response:
[
\nabla^2\Phi = 4\pi G\rho
\quad\longrightarrow\quad
\nabla\cdot\big(\varepsilon(\mathbf x)\nabla\Phi\big) = 4\pi G\rho,
]
or an anisotropic (\varepsilon_{ij}) for directional impedance.

**The harsh scaling:** To reduce Earth’s (g) locally by a significant fraction, you must create an *opposing field* comparable to (g), which (in any Poisson-like theory) corresponds to enormous effective source strengths.

Two simple back-of-envelope examples we used:

1. **Point-source equivalent:** to produce acceleration (g) at distance (r),
   [
   g \sim \frac{GM}{r^2}\quad\Rightarrow\quad M \sim \frac{g r^2}{G}.
   ]
   At (r=10) m, this is (M\sim 10^{13}) kg (mountain-scale).

2. **Infinite-sheet equivalent:** a sheet gives (g\sim 2\pi G\Sigma), so
   [
   \Sigma \sim \frac{g}{2\pi G}\sim 10^{10},\mathrm{kg/m^2},
   ]
   again implying absurd effective “mass charge” densities.

**Conclusion:** A static “gravity shield” that cancels Earth’s field is not merely hard—it is equivalent to engineering gigantic effective sources or polarizations. Within the GR-matching regime, there is no easy knob that avoids that conclusion.

### 5.3 The 1PN/hybrid rigidity effect: weak-field knobs are tiny

A related point is that in the GR-matching corner, Earth’s gravitational potential is extremely weak in relativistic units:
[
\frac{|\Phi_\oplus|}{c^2}\sim 10^{-9}.
]
So any correction that is genuinely “post-Newtonian small” (e.g., a 2PN-type effect) tends to be suppressed by additional powers of (|\Phi|/c^2), making it negligible for engineering unless you access a nonlinear/strong-field regime.

This reinforced the idea that “scalar-sector tweaks” are unlikely to give usable lift near Earth while staying consistent with your 1PN/EIH matching assumptions.

### 5.4 The pivot: from “cancel gravity” to “vacuum lift” and then to “vacuum propulsion”

These limits forced two pivots:

1. **Boundary-supported lift (wake cavity):**
   Instead of canceling (\nabla\Phi), create an **upward traction** from the transverse wake stress. That yields modest-looking pressure requirements for ton-scale loads, but it fundamentally relies on a reaction partner (ground/plate), making it “maglev-like.”

2. **Free-space hover as propulsion:**
   If we insist on free-space hover, we must export momentum into the medium. This is the Wright-brothers framing: look for a medium momentum carrier with low energy-per-momentum ((\varepsilon \ll c)). Radiation-like carriers fail ((\varepsilon\sim c)); vortex-ring–like carriers are attractive because they can have (\varepsilon\sim U) with (U) in the tens of m/s range—opening an engineering window.

That is the point where the project stops being “rebranded EM” and becomes a concrete mathematical question about whether your superfluid ontology admits slow, coherent momentum structures and whether they can be actuated efficiently and selectively.

---

## 6. Boundary-supported “wake cavity under a disk” (ground-effect levitation)

This section documents the full “wake cavity” route: support a payload by creating an **upward traction** (pressure difference) in the transverse/wake sector, with the **ground acting as the momentum partner** via boundary reflection (image forces). This route is valuable because it produces concrete numbers quickly and teaches us what is—and is not—possible without exporting momentum into free space.

### 6.1 Required support traction

For a craft of mass (m) supported over area (A=\pi R^2), the required pressure difference is
[
F=mg=\Delta p,A
\quad\Rightarrow\quad
\boxed{\Delta p=\frac{mg}{\pi R^2}}.
]

Numerical examples for (m=1000,\mathrm{kg}) ((mg\simeq 9.81\times 10^3,\mathrm{N})):

* (R=1,\mathrm{m}): (\Delta p \approx 3.1\times10^3,\mathrm{Pa}) (3 kPa)
* (R=2,\mathrm{m}): (\Delta p \approx 7.8\times10^2,\mathrm{Pa}) (0.8 kPa)
* (R=5,\mathrm{m}): (\Delta p \approx 1.25\times10^2,\mathrm{Pa}) (125 Pa)

**Key takeaway:** the *traction magnitude* needed for ton-scale support is modest if the supporting area is large.

### 6.2 Stress/energy interpretation and stored energy scale

In any field/continuum system, a subsystem with energy density (u) can exert stresses/pressures of order (u). We used the conservative scaling
[
\Delta p \sim u_{\rm wake}.
]

If the supporting stress occupies a working “cavity thickness” (\ell), the stored energy is
[
U_{\rm cav}\sim u_{\rm wake}(A\ell)\sim \Delta p,A,\ell=(mg),\ell.
]

So for 1 ton and (\ell=1,\mathrm{m}),
[
U_{\rm cav}\sim 10,\mathrm{kJ}.
]

This was an important intuition: hovering is not automatically a power monster if it is quasi-static and the medium can store stress with low dissipation; the hard part is *creating* and *maintaining* the stress configuration and ensuring it yields a net force on the craft.

### 6.3 EM-template amplitude check (“B-pressure”)

To sanity-check magnitudes, we temporarily used the electromagnetic pressure template
[
u_B=\frac{B^2}{2\mu_0}
\quad\Rightarrow\quad
\Delta p\simeq \frac{B^2}{2\mu_0}
\quad\Rightarrow\quad
\boxed{B\simeq \sqrt{2\mu_0\Delta p}}.
]

For (\Delta p\approx 125,\mathrm{Pa}),
[
B \approx 0.018,\mathrm{T}\quad (18~\mathrm{mT}).
]
This is far below “frog magnet” multi-tesla fields, highlighting that the frog experiment is hard mainly because diamagnetic coupling is tiny—not because stress-based support necessarily requires huge field energy density.

### 6.4 The image-force model: disk as a current loop above a boundary

To compute a concrete lift law, we modeled the craft as an effective loop current (I) of radius (R) at height (h) above a plane boundary (the ground). We parameterized the boundary’s reflectivity to the transverse wake mode by (0\le\eta\le 1):

* (\eta=1): perfect reflector (full image strength)
* (\eta=0): no reflection (transparent; no image force)

Then the ground’s effect is represented by an image loop at (-h) with current (-\eta I), producing repulsive lift. The exact expression for lift in terms of the mutual inductance (M(z)) of two coaxial loops separated by (z) is:
[
\boxed{F_\uparrow(h)=\eta, I^2\left(-\frac{dM}{dz}\right)_{z=2h}}.
]

This is a clean stress-tensor / Green’s-function result: lift is the derivative of interaction energy with respect to height.

### 6.5 Small-gap approximation and cm-hover scaling

For (h\ll R), the interaction has a simple asymptotic form leading to
[
F_\uparrow \approx \eta,I^2\frac{\mu_0 R}{2h},
]
so balancing weight gives
[
\boxed{I(h,\eta)\approx \sqrt{\frac{2h,mg}{\eta,\mu_0,R}}}.
]

**This was the practical “roads not rails” motivation:** if you only need a few centimeters of lift, (h) is small, and required source strength scales like (\sqrt{h}).

For a 1-ton payload, (R=5,\mathrm{m}), (h=2,\mathrm{cm}):

* (\eta=1): (I\sim 8,\mathrm{kA})
* (\eta=0.1): (I\sim 25,\mathrm{kA})

(Always remembering this is an *effective* “wake current,” not necessarily an electrical current.)

### 6.6 Why this does not give free-space hover

The image-force mechanism works because the boundary participates in the stress field: the craft pushes on the field, the field pushes on the ground, and the ground pushes back. In free space (no boundary), a static near-field stress configuration cannot yield a net force on an isolated craft (Section 3). Therefore:

* The wake cavity / image-force route is a **ground-effect support technology**.
* It is potentially valuable, but it does not satisfy the “free-space hover” requirement unless one adds a different momentum partner (i.e., exports momentum via propagation).

This conclusion is what motivated the pivot to a propulsion-by-momentum-flux approach.

---

## 7. “Not just EM” route: engineered coupling and wake penetration depth

Section 6 showed that boundary-supported levitation is mathematically straightforward but risks being “just maglev in different clothes.” The only way to make it meaningfully new (even as a boundary technology) is if the **transverse wake mode** has a materially distinct coupling to matter—so that engineered surfaces can reflect/pin it without being electrical superconductors or active-coil infrastructure.

This section formalizes that idea in one parameter: a wake penetration depth (\lambda_T).

### 7.1 Minimal linear response model inside matter

We assumed the transverse wake field (\mathbf A) (with (\mathbf A\propto\mathbf v_T)) induces an opposing response in a material layer:
[
\boxed{\mathbf J_{\rm ind} = -\chi_T,\mathbf A}.
]

Substituting into the quasi-static field equation gives a screened Helmholtz form inside the material:
[
(\nabla^2-\lambda_T^{-2})\mathbf A=0,
\qquad
\boxed{\lambda_T^{-2}=\mu_T,\chi_T},
]
where (\mu_T) is the transverse-sector analog of (\mu_0).

**Interpretation:** (\lambda_T) is the depth over which transverse wake fields penetrate the material. Small (\lambda_T) means strong screening/pinning and strong reflection—exactly what we need for large (\eta) in the image-force mechanism.

### 7.2 Half-space reflectivity (\eta(k)) from mode matching

For a planar interface, decompose into lateral Fourier modes with wavenumber (k). Above the interface (vacuum) modes decay with vertical rate (\alpha_1=k). Inside the screened medium, the vertical rate is
[
\alpha_2=\sqrt{k^2+\lambda_T^{-2}}.
]

Matching yields an effective reflection amplitude:
[
\boxed{\eta(k)=\frac{\alpha_2-\alpha_1}{\alpha_2+\alpha_1}
==========================================================

\frac{\sqrt{k^2+\lambda_T^{-2}}-k}{\sqrt{k^2+\lambda_T^{-2}}+k}}.
]

For a craft of characteristic radius (R), the dominant lateral scale is (k\sim 1/R). This converts (\lambda_T) into a single-number boundary reflectivity (\eta) relevant for lift.

A useful inversion (with (k\sim 1/R)) is:
[
\boxed{\lambda_T \approx \frac{1-\eta}{2\sqrt{\eta}},R}.
]

This showed that you do **not** need nanometer-scale penetration depths to get useful reflectivity at meter-scale craft sizes. Meter-scale (\lambda_T) can still give (\eta\sim 0.1)–0.7 depending on (R).

### 7.3 What makes this “not just EM”

If (\mathbf A) were literally the electromagnetic vector potential, then strong DC screening is essentially the Meissner effect (superconductivity), and engineered reflectors would require superconductors or active metasurfaces. That’s not fundamentally new.

In your ontology, however, (\mathbf A\equiv \kappa_A\mathbf v_T) is a **brane-confined transverse wake mode** introduced precisely to resolve the “bulk irrotational (\Rightarrow B=0)” obstruction. The *new* possibility is that certain microstructures could couple strongly to (\mathbf v_T) (pin/screen it), yielding small (\lambda_T) and large (\eta), without being electrical superconductors.

This is the Wright-brothers-like “unnatural coupler” idea: a metamaterial that does not exist in nature but forces the vacuum medium to satisfy a boundary condition.

### 7.4 Why this still doesn’t solve free-space hover

Even if a metamaterial pavement achieved (\eta\to 1), it would still be a boundary technology: it needs the ground to supply reaction force. It could be economically transformative (roads instead of rails), but it does not meet the requirement of free-space hovering. That requirement forces us beyond boundary reflection and into momentum export (Sections 8–11), which is where vortex-ring propulsion entered as the most promising “not-just-EM” path.

---

## 8. Pivot to true free-space hover: propulsion by exporting momentum into the vacuum

Once we accept the momentum-partner constraint (Section 3), “free-space levitation” becomes a propulsion problem: a craft must generate a sustained upward force by exporting momentum downward into the surrounding medium. In this section we formalize that as a general energy–momentum bookkeeping law and identify the fundamental “photon rocket” barrier that any viable mechanism must avoid.

### 8.1 Universal momentum and power relations

Let a device emit momentum-carrying packets (waves, wakes, coherent structures) downward into the medium. If each packet carries downward momentum (P_{\rm pkt}) and energy (E_{\rm pkt}), and the emission rate is (\dot N), then the upward thrust is
[
F = \dot N,P_{\rm pkt}.
]
The required hover condition is (F=mg).

The power spent on emission is
[
\mathcal P = \dot N,E_{\rm pkt}.
]

Eliminating (\dot N) yields the central relation:
[
\boxed{\mathcal P = F\left(\frac{E_{\rm pkt}}{P_{\rm pkt}}\right) \simeq (mg),\varepsilon},
]
where
[
\varepsilon \equiv \frac{E}{P}
]
is the energy-per-momentum of the chosen momentum carrier.

This is the cleanest “sanity check” for free-space hover: you don’t need detailed microphysics to know that if (\varepsilon) is too large, the power requirement is immediately fatal.

### 8.2 The photon/wave barrier: why radiation-like thrust is dead

For radiation-like carriers propagating at (c), the energy–momentum ratio is (\varepsilon \sim c). (For photons, (E=pc) exactly.)

Then
[
\mathcal P \sim (mg),c.
]

For (m=1000,\mathrm{kg}), (mg\approx 10^4,\mathrm{N}), this is
[
\mathcal P \sim 10^4 \times 3\times 10^8 \sim 3\times 10^{12},\mathrm{W\ per\ ton}.
]

That is the “photon rocket” bound. It does not depend on clever engineering; it is a consequence of momentum conservation plus the dispersion relation of the carrier. If the only way the vacuum medium can accept momentum is via fast wave modes with (\varepsilon\sim c), then practical free-space hover is dead on arrival.

### 8.3 The only route forward: find a slow momentum carrier

Therefore, any viable “vacuum hover” must use a carrier for which
[
\varepsilon \ll c.
]

In an ordinary fluid, the way to achieve (\varepsilon \ll c) is obvious: you accelerate bulk material (air) to modest speeds. For a superfluid vacuum, the analogous route is to export momentum in a carrier whose energy-per-momentum is set by a **low group speed**, not by (c). This is the structural reason we pivoted to vortex-like topological excitations: they can carry large momentum while moving at speeds that are not tied to the vacuum light speed.

This sets the target: identify a slow momentum carrier in the toy model and quantify whether the resulting power-per-ton can plausibly be in the (\sim 0.1)–(1) MW/ton range, rather than terawatts/ton.

---

## 9. Vacuum “propeller” candidate: vortex-ring thrust in a superfluid medium

Within superfluid and ideal-fluid theory, **vortex rings** are coherent, topological momentum packets. Their key property for our purposes is that their energy-per-momentum is set by their translation speed (U), which can be many orders of magnitude smaller than (c). This makes them the most promising “Wright-brothers propeller analog” inside a superfluid ontology.

### 9.1 Vortex ring as a momentum packet: impulse and energy scalings

For a thin-core vortex ring of radius (R), core radius (a), background density (\rho), and circulation (\Gamma), classic scalings are:

* Impulse (momentum):
  [
  \boxed{P_{\rm ring} \sim \rho,\Gamma,\pi R^2}
  ]
* Energy:
  [
  \boxed{E_{\rm ring} \sim \tfrac12,\rho,\Gamma^2,R,\ln!\left(\frac{8R}{a}\right)}
  ]

These are order-of-magnitude formulas; the logarithm encodes the “thin core” structure and changes slowly with (R/a).

The ring translation speed is
[
\boxed{U_{\rm ring} \sim \frac{\Gamma}{4\pi R}\left[\ln!\left(\frac{8R}{a}\right)-\mathcal O(1)\right]}.
]

### 9.2 The decisive ratio: (E/P\sim U)

Dividing energy by momentum yields
[
\frac{E_{\rm ring}}{P_{\rm ring}}
\sim
\frac{\Gamma}{2\pi R},\ln!\left(\frac{8R}{a}\right)
\sim U_{\rm ring}.
]

So for vortex rings,
[
\boxed{\varepsilon \equiv \frac{E}{P}\sim U_{\rm ring}}.
]

Plugging into the universal hover power law (Section 8) gives:
[
\boxed{\mathcal P \sim (mg),U_{\rm ring}}.
]

This is exactly the “propeller” structure: thrust times an effective downwash/exhaust speed. It is the mathematical reason vortex rings are the first candidate with genuine “legs” beyond EM.

### 9.3 Power-per-ton targets translate into ring speeds

For (m=1000,\mathrm{kg}), (mg\approx 10^4,\mathrm{N}), we get the useful conversion:
[
\mathcal P \approx 10^4,U_{\rm ring}\ \ (\mathrm{W\ per\ ton}).
]

Examples:

* (U_{\rm ring}=5,\mathrm{m/s}) → (\mathcal P\sim 50,\mathrm{kW/ton})
* (U_{\rm ring}=20,\mathrm{m/s}) → (\mathcal P\sim 200,\mathrm{kW/ton})
* (U_{\rm ring}=50,\mathrm{m/s}) → (\mathcal P\sim 500,\mathrm{kW/ton})
* (U_{\rm ring}=200,\mathrm{m/s}) → (\mathcal P\sim 2,\mathrm{MW/ton})

This gives an immediate engineering interpretation: if we want hover power budgets comparable to heavy industrial machinery, we should be looking for vortex-ring-like carriers with (U) in the **tens of m/s** range.

### 9.4 Mapping to your model variables: circulation is the lever

Your EM paper already defines circulation as the loop integral of the flow field:
[
\Gamma \equiv \oint \mathbf v\cdot d\boldsymbol\ell,
]
with units (\mathrm{m^2/s}). That is the same (\Gamma) that controls vortex ring dynamics.

From the speed formula,
[
U_{\rm ring}\approx \frac{\Gamma}{4\pi R},L,
\qquad
L\equiv \ln!\left(\frac{8R}{a}\right)-\frac12.
]
So to achieve a target (U), the required circulation is
[
\boxed{\Gamma \approx \frac{4\pi R}{L},U}.
]

For example, with (R=0.5,\mathrm{m}), (L\sim 20), and (U=20,\mathrm{m/s}), we found (\Gamma\sim 7,\mathrm{m^2/s}). This is not a tiny number; it encodes the requirement that the device impose a large coherent circulation on the vacuum medium—exactly the “vacuum propeller” actuation challenge.

### 9.5 Core-speed constraint and why core size matters

A vortex ring also has a characteristic swirl speed near its core of order
[
v_{\max}\sim \frac{\Gamma}{2\pi a}.
]

In your GR-matching picture the relevant characteristic wave speed in the medium is (c_s\sim c), so physically we require (v_{\max}\lesssim c) (or some fraction of (c) if we demand non-relativistic internal motion). This implies
[
\boxed{\Gamma \lesssim 2\pi a,c}.
]

Combining this with the target (U) relation implies that to achieve (U\sim 10)–(100) m/s with macroscopic ring radii (R), the core radius (a) cannot be arbitrarily large; it tends to land in atomic-to-nanometer scales for meter-scale rings unless one allows core swirl speeds close to (c). This is not yet a kill shot, but it is a clear requirement: “vacuum vortices” must be highly concentrated structures if they are to translate at useful speeds without requiring impossible circulation.

### 9.6 Why this is “not just EM”

This vortex-ring thrust route is qualitatively different from EM-based levitation:

* For EM waves/photons, (\varepsilon=E/P=c), making hover power astronomical.
* For vortex rings, (\varepsilon\sim U), which can be tens of m/s, making hover power potentially tractable.

So the “not just EM” novelty is precisely the existence of slow, momentum-rich topological excitations in the vacuum medium and the possibility of engineered couplers to generate them. This is the closest vacuum analog to the Wright brothers’ propeller insight: a purpose-built mechanism that exchanges momentum with the environment efficiently by choosing a carrier whose energy-per-momentum is small.

---

## 10. Ring emission rate and plume clearance (operational viability)

Sections 8–9 show that vortex rings can evade the photon-rocket power wall because (E/P \sim U), giving (\mathcal P \sim (mg)U). That is necessary but not sufficient: a practical system must also export momentum **without piling up** its own wake structures underneath it. This section derives the required ring emission cadence (\dot N) and introduces a simple “clearance” criterion.

### 10.1 Emission rate from momentum balance

Hover requires thrust (F=mg). If each ring carries downward impulse (P_{\rm ring}), and rings are emitted at rate (\dot N),
[
mg = \dot N,P_{\rm ring}
\quad\Rightarrow\quad
\boxed{\dot N = \frac{mg}{P_{\rm ring}}}.
]

Using the thin-core vortex-ring impulse scaling, mapped to your notation with an effective density (\rho_{\rm eff}) (the density relevant to the transverse momentum carrier),
[
\boxed{P_{\rm ring}\approx C_P,\rho_{\rm eff},\Gamma,\pi R^2},
]
with (C_P\sim\mathcal O(1)).

So
[
\boxed{\dot N \approx \frac{mg}{C_P,\rho_{\rm eff},\Gamma,\pi R^2}}.
]

This shows the basic operational lever: for fixed thrust, emission rate decreases with larger (\rho_{\rm eff}), larger circulation (\Gamma), and larger ring radius (R).

### 10.2 Eliminating (\Gamma) using the target downwash speed (U)

From Section 9, the ring translation speed is approximately
[
U \approx \frac{\Gamma}{4\pi R},L,\qquad
L\equiv \ln!\left(\frac{8R}{a}\right)-\frac12,
]
so
[
\Gamma \approx \frac{4\pi R}{L},U.
]

Insert this into (P_{\rm ring}):
[
P_{\rm ring}\approx C_P,\rho_{\rm eff},\pi R^2\left(\frac{4\pi R}{L}U\right)
============================================================================

\left(\frac{4\pi^2C_P}{L}\right)\rho_{\rm eff},U,R^3.
]

Therefore the emission rate becomes
[
\boxed{\dot N \approx \frac{mg,L}{4\pi^2C_P,\rho_{\rm eff},U,R^3}}.
]

**Critical scaling:** for fixed hover power target (fixed (U)),
[
\dot N \propto \frac{1}{\rho_{\rm eff} R^3}.
]
This is why ring size matters so much: doubling (R) cuts (\dot N) by (8\times).

### 10.3 A practical “cargo” baseline and numerical emission rates

We adopted a cargo-relevant hover target of roughly
[
U \sim 20,\mathrm{m/s}
\quad\Rightarrow\quad
\mathcal P \sim (mg)U \sim 200,\mathrm{kW/ton}.
]

Taking (L\sim 20), (C_P\sim 1), and (mg\approx 9810,\mathrm{N}), we get the useful approximation (for 1 ton):
[
\dot N \approx \frac{254}{\rho_{\rm eff}R^3},
]
with (\rho_{\rm eff}) in kg/m(^3) and (R) in meters.

Two illustrative ring sizes:

**(A) (R=0.05,\mathrm{m}) (5 cm rings):** (R^3=1.25\times 10^{-4})
[
\dot N \approx \frac{2.0\times 10^6}{\rho_{\rm eff}}\ \mathrm{s^{-1}}.
]
So:

* (\rho_{\rm eff}=1) → (2\times 10^6) rings/s (operationally hopeless)
* (\rho_{\rm eff}=10^3) → (2\times 10^3) rings/s (still extreme)
* (\rho_{\rm eff}=10^5) → (20) rings/s (sane)

**(B) (R=0.5,\mathrm{m}) (50 cm rings):** (R^3=0.125)
[
\dot N \approx \frac{2030}{\rho_{\rm eff}}\ \mathrm{s^{-1}}.
]
So:

* (\rho_{\rm eff}=1) → (2000) rings/s
* (\rho_{\rm eff}=100) → (20) rings/s
* (\rho_{\rm eff}=1000) → (2) rings/s

**Interpretation:** even with modest power targets, practical hover strongly prefers either (i) larger rings, or (ii) large effective density (\rho_{\rm eff}), or both.

### 10.4 Clearance / pile-up criterion

Even if thrust balances weight on average, the wake must “clear” so that momentum is actually transported away rather than accumulating underneath the craft.

A crude but useful diagnostic is the mean axial spacing of rings in the plume:
[
d \sim \frac{U}{\dot N}.
]

For (U=20) m/s:

* (\dot N=2000)/s → (d\sim 1) cm (strong interaction, turbulent plume)
* (\dot N=200)/s → (d\sim 10) cm (still strongly interacting)
* (\dot N=20)/s → (d\sim 1) m (manageable discrete packets)
* (\dot N=2)/s → (d\sim 10) m (very clean)

This suggests a practical regime for controlled ring propulsion is (\dot N\sim 1)–(50) rings/s per ton per nozzle, unless the intent is to operate as a continuous “jet” (in which case the system is no longer ring-dominated but turbulence-dominated).

### 10.5 Energy per ring and pulse scale

Given a fixed hover power (\mathcal P) and ring rate (\dot N),
[
E_{\rm ring}\approx \frac{\mathcal P}{\dot N}.
]

At (\mathcal P\sim 200) kW/ton:

* (\dot N=20)/s → (E_{\rm ring}\sim 10) kJ per ring
* (\dot N=2)/s → (E_{\rm ring}\sim 100) kJ per ring

These kJ–100 kJ pulses are within the scale of industrial pulsed power systems, reinforcing that the math is not immediately absurd if (\rho_{\rm eff}) and (R) fall in favorable ranges.

---

## 11. Anchoring (\rho_0) and core sizes using your defect mass relation

The emission-rate analysis depends on (\rho_{\rm eff}) (the density relevant to the momentum-carrying transverse sector). Your papers provide a natural density scale (\rho_0) through the defect mass formula; this section shows how we used that to generate plausible numerical anchors consistent with a “useful tech” objective.

### 11.1 Defect mass formula and the density–core-size relationship

Your model relates gravitational mass of a defect to the background density and the defect’s geometric volume:
[
m_G = \kappa_m,\rho_0,V
\approx \kappa_m,\rho_0,\pi a^2 L
= \kappa_m,\rho_0,\pi a^3\Lambda,
]
where (a) is a characteristic throat/core radius, (L) is an effective length scale, and (\Lambda\equiv L/a\approx 1.85) is the dimensionless ratio found in your project work.

Solving for (a) gives
[
\boxed{a(\rho_0) = \left(\frac{m_G}{\kappa_m\rho_0\pi\Lambda}\right)^{1/3}}.
]

This tells us: for a fixed “elementary mass” (m_G), larger (\rho_0) implies smaller defect cores.

### 11.2 Choosing an anchor appropriate to “known mass scales”

Since your goal is cargo tech rather than particle microphysics, we did not try to force exact identification with standard-model charges. But we still need a plausible anchor to estimate orders of magnitude. A conservative choice is:

* **Anchor:** one elementary defect corresponds to a proton-scale mass (m_G\sim m_p).

Then, taking (\kappa_m\sim 1) and (\Lambda\approx 1.85), we obtained:

* If (\rho_0=20) kg/m(^3) → (a\approx 0.24) nm (atomic scale)
* If (\rho_0=100) kg/m(^3) → (a\approx 0.14) nm
* If (\rho_0=1000) kg/m(^3) → (a\approx 0.066) nm

These are “atomic-ish” throat/core sizes, which is at least qualitatively consistent with the idea that the toy model’s defect microstructure lives at very small length scales, while macroscopic gravity emerges from many defects.

### 11.3 Linking (\rho_0) to hover viability through (\rho_{\rm eff})

The ring emission analysis is controlled by (\rho_{\rm eff}), not directly by (\rho_0). However, within your ontology it is natural (as a first approximation) to take (\rho_{\rm eff}\sim \rho_0) up to an order-unity sector factor—because both represent background “mass density” available to store momentum in the medium.

Under that identification, the same density range that yields atomic-scale defect cores ((\rho_0\sim 10)–(10^3) kg/m(^3)) also lands us in the operationally favorable emission-rate regime:

* For (R=0.5) m, (U=20) m/s:
  (\rho_{\rm eff}\sim 100) kg/m(^3) → (\dot N\sim 20) rings/s per ton
  (\rho_{\rm eff}\sim 1000) kg/m(^3) → (\dot N\sim 2) rings/s per ton

This is the main reason we concluded the idea “still has legs” mathematically: there is a self-consistent parameter window where (i) defect cores remain microscopic, (ii) hover power per ton is in the 0.1–1 MW range, and (iii) ring emission cadence is not absurd.

### 11.4 What remains unresolved (and what we must test next)

This anchoring step is not a proof; it is a plausibility check. The next-stage “kill tests” are:

1. whether (\rho_{\rm eff}) for the transverse momentum carrier can indeed be of order (\rho_0) (and not many orders smaller),
2. whether the required circulation (\Gamma\sim\mathcal O(1!-!10),\mathrm{m^2/s}) can be generated without producing catastrophic electromagnetic signatures (if the EM mapping ties (\Gamma) to ordinary charge), and
3. whether stable vortex-ring–like packets exist in the appropriate sector under your model’s constraints (bulk irrotationality vs brane transverse modes).

But crucially, **the bookkeeping and scaling analysis itself did not kill the concept**: it identified the narrow set of medium properties that must hold for free-space hover to remain feasible.

---

## 12. “Not just EM” consistency requirements and kill tests

Everything up to Section 11 establishes a *mathematically viable window* for free-space hover **if** the vacuum medium supports slow momentum carriers (vortex-ring–like packets) and devices can generate them efficiently. Section 12 is the reality gate: what must be true for this to be “not just EM,” and what checks could still kill the idea.

### 12.1 Why ordinary EM coupling is a potential show-stopper

In your EM mapping, circulation and vorticity are directly tied to EM-like fields, and in some parts of the dictionary, “charge” is expressed as proportional to (\rho_0 a^2 \Gamma) (up to constants). If the same circulation (\Gamma) we need for propulsion also implied large *ordinary* EM charge/current, then producing the required vortex structures would:

* generate enormous electric/magnetic signatures,
* couple strongly to ordinary matter and detectors,
* and be inconsistent with everyday observations (we would already see these effects).

Therefore, the “free-space hover has legs” window almost certainly requires one of the following:

1. **A distinct transverse momentum-carrying sector** (“dark wake” mode) that is *not* the photon field, even if it has Maxwell-like equations in vacuum; or
2. **A conditional coupling** where the effective mapping to ordinary EM charge is suppressed except in specific states/structures; or
3. **Mode separation**: the vortex rings exist in (\mathbf v_T) as a medium mode but do not map one-to-one onto electromagnetic (\mathbf B) in the way naive identification would suggest.

This is the core “not just EM” requirement: the momentum carrier must be slow and engineerable, but must not automatically behave like ordinary EM vorticity.

### 12.2 Why “we would have noticed” constraints are powerful

If a new vacuum momentum carrier couples strongly to generic matter in ordinary conditions, it tends to imply new forces, shielding effects, dissipation signatures, or radiation that would have shown up in laboratory physics.

To avoid that, any viable not-just-EM route must plausibly satisfy:

* **In ordinary matter/conditions:** coupling is weak (no anomalous forces).
* **In engineered structures/phases:** coupling becomes strong (enables actuation).

This is exactly the Wright-brothers logic: the medium is there all along, but you need an engineered coupling element (propeller) to access the useful regime.

### 12.3 Kill tests identified so far

These are the most direct ways the idea could fail even if the bookkeeping works:

**Kill test 1: Only (c)-speed carriers exist.**
If every momentum-carrying excitation accessible to matter has (\varepsilon=E/P \sim c), then (\mathcal P\sim mgc) and hover is dead.

**Kill test 2: (\rho_{\rm eff}) is too small.**
If the effective density participating in transverse momentum transport is orders of magnitude below the (\rho_0) scale, then practical ring sizes require absurd emission rates (thousands–millions per second per ton).

**Kill test 3: Required (\Gamma) implies unavoidable EM signatures.**
If generating (\Gamma\sim 1)–(10,\mathrm{m^2/s}) in the relevant sector necessarily produces large ordinary electromagnetic fields/charges/currents, then the mechanism is observationally excluded or becomes an EM rocket again.

**Kill test 4: No stable coherent vortex packets in the allowed sector.**
Even if a transverse mode exists, it might not admit stable ring solutions under the ontology constraints (bulk irrotationality, brane confinement, dissipation, etc.).

**Kill test 5: Actuation losses dominate.**
The power law (\mathcal P\sim mgU) is an ideal momentum-exchange limit. Real devices might dissipate far more power while trying to create rings, wiping out the advantage.

### 12.4 What would count as “evidence” inside the model

Within your framework, the strongest signs that we’re not just rebranding EM would be:

* A well-defined transverse sector with slow coherent structures that do not map to photon radiation,
* Material-dependent coupling rules that allow strong actuation only in engineered states,
* And predicted observables that differ from EM (e.g., momentum flux with minimal EM fields).

---

## 13. Practical architecture sketches (conceptual)

These are not implementation claims; they are conceptual “system layouts” that follow directly from the scaling laws. They help organize what must be built if the mechanism is real.

### 13.1 “Vacuum helicopter” layout: ring nozzles as thrust units

Treat each vortex-ring emitter as a thrust element that produces a downward momentum packet with speed (U). For a vehicle with total mass (m), the total thrust requirement is (mg), which can be distributed across (N_e) emitters.

* Per-emitter thrust:
  [
  F_e \approx \frac{mg}{N_e}.
  ]
* Per-emitter power:
  [
  \mathcal P_e \approx \frac{mg}{N_e},U.
  ]
* Per-emitter ring rate:
  [
  \dot N_e \approx \frac{mg}{N_e,P_{\rm ring}}.
  ]

This immediately suggests:

* many moderate emitters (redundancy, controllability),
* rather than one extreme central emitter.

### 13.2 Scale choice: ring radius vs plume control

From Section 10, (\dot N\propto 1/(\rho_{\rm eff}R^3)) at fixed (U). That strongly favors larger (R), but large rings may:

* be harder to generate cleanly,
* have stronger interactions with the craft geometry,
* and introduce stability/control issues.

A plausible architecture is multi-nozzle:

* each nozzle emits rings of (R\sim 0.1)–1 m,
* arranged in a symmetric array under the craft.

### 13.3 Stability and control (medium analog of rotorcraft)

Even in ideal momentum exchange, a hover platform must control:

* total thrust (altitude),
* differential thrust (pitch/roll),
* lateral thrust components (translation),
* and yaw torque balance.

A ring-emitter array provides these controls by modulating:

* emission rate (\dot N),
* ring speed (U) (via (\Gamma)),
* and emitter phasing/spatial distribution.

Conceptually this is closer to a multi-rotor drone than to a maglev train: free-space hovering demands active control.

### 13.4 Operational modes

Two broad operating modes emerge:

1. **Discrete-ring mode:** low (\dot N), coherent rings, meters of spacing in the plume; likely lower dissipation and more “structured” momentum transport.

2. **Continuous jet mode:** high (\dot N), strong interaction and turbulence, effectively a downward momentum jet. This may be easier to realize if coherence is hard, but may increase losses and noise (in the medium sense).

---

## 14. Summary of what we have so far and what’s next

### 14.1 What is established (within our assumptions)

1. **Free-space hover requires momentum export.**
   It cannot be achieved by a purely static self-contained field configuration without an external momentum partner.

2. **The key viability metric is (\varepsilon = E/P).**
   Hover power is (\mathcal P = (mg)\varepsilon). Radiation-like carriers with (\varepsilon\sim c) imply terawatts per ton and are dead.

3. **Vortex-ring–like carriers evade the photon wall.**
   For vortex rings, (E/P\sim U), yielding (\mathcal P\sim (mg)U). For cargo-relevant budgets (0.1–1 MW/ton), this corresponds to (U\sim 10)–100 m/s.

4. **Operational feasibility depends on emission rate and effective density.**
   Ring emission rate scales like
   [
   \dot N \approx \frac{mg,L}{4\pi^2 C_P,\rho_{\rm eff},U,R^3}.
   ]
   Reasonable ring rates (tens per second per ton) occur if (\rho_{\rm eff}) is tens–hundreds kg/m(^3) and ring radii are (\sim 0.1)–1 m.

5. **A plausible anchoring exists without immediate contradiction.**
   Using your defect mass formula with a proton-mass anchor, (\rho_0\sim 10)–(10^3) kg/m(^3) corresponds to atomic-scale defect cores, and that same density range is favorable for ring emission rates—so the math does not self-destruct at the level of dimensional analysis.

### 14.2 What remains unknown (and most likely to kill the concept)

* Whether the relevant transverse momentum carrier is **distinct from ordinary EM**, or whether creating (\Gamma) inevitably produces ordinary EM signatures.
* Whether (\rho_{\rm eff}) for the momentum carrier is indeed of order (\rho_0), or is effectively tiny.
* Whether stable coherent vortex-ring packets exist in the allowed sector (brane wake vs bulk irrotational constraints).
* Whether actuation can approach the ideal (\mathcal P\sim mgU) limit without enormous dissipative overhead.

### 14.3 Recommended next-session work plan

To move from “legs on paper” to “alive or dead,” the next tasks are:

1. **Derive (\rho_{\rm eff}) for the transverse momentum carrier** from your ontology and hybrid/1PN normalization (is it (\sim\rho_0), or suppressed?).

2. **Clarify sector identity:**
   Does the vortex-ring carrier live in the same field that maps to ordinary EM (making it likely excluded), or is there a “dark transverse wake” mode consistent with the ontology?

3. **Consistency constraints:**
   Work out whether large (\Gamma) excitations imply measurable EM fields/charges. If yes, quantify the signatures; if no, specify the coupling suppression mechanism.

4. **Package the results into a “viability worksheet”:**
   Inputs: (m), target (\mathcal P/m), chosen (R), chosen (a), assumed (\rho_{\rm eff}).
   Outputs: required (U), (\Gamma), (\dot N), energy per ring, plume spacing.
   This makes future exploration modular and prevents drift.

That completes the core narrative: we began with boundary-supported levitation and discovered it cannot deliver free-space hover, then pivoted to the only remaining consistent route—vacuum propulsion by momentum export—and found a mathematically plausible window if the toy model admits a slow, momentum-rich transverse excitation that does not collapse back into ordinary electromagnetism.

---

## Appendix contents to add

### Appendix A. Notation and dictionary (quick reference)

* Medium variables: (\rho_0), (h), (\Phi), (\mathbf v), (\mathbf v_T)
* EM mapping: (\phi_{\rm EM}=\lambda(h+\tfrac12 v^2)), (\mathbf A=\lambda\mathbf v), (\mathbf B=\nabla\times\mathbf A)
* Ontology mapping: (\mathbf A=\kappa_A\mathbf v_T), (\nabla^2\mathbf A\sim-) (wake current source)
* Defect geometry: core radius (a), length (L), ratio (\Lambda=L/a\simeq 1.85), packing constant (\kappa_m)
* Vortex-ring geometry: ring radius (R), core radius (a_{\rm ring}), log factor (L_{\log}=\ln(8R/a_{\rm ring})-\tfrac12)

---

### Appendix B. Momentum-flux derivation: why free-space hover requires momentum export

* Sketch derivation from local momentum conservation:
  [
  \partial_t(\text{momentum density})+\nabla\cdot T = 0
  ]
  leading to
  [
  \mathbf F = \oint T\cdot \mathbf n,dA
  ]
* Specialization to time-averaged steady hover: why internal/static configurations give zero net force.
* Explicit statement of the three momentum partners (boundary, radiation/excitation, long-range coupling).

---

### Appendix C. “Photon rocket” power bound

* Derive (P=Fc) for radiation thrust from (E=pc).
* Numerical example per ton:
  [
  P \sim mgc \approx 3\times10^{12}\ \text{W/ton}
  ]
* Generalization: any carrier with group speed (\approx c) has (\varepsilon=E/P\sim c).

---

### Appendix D. Boundary levitation math details (loop + image)

**D.1 Force from mutual inductance**

* Energy of two coupled loops:
  [
  U(z) = - I_1 I_2 M(z)
  ]
* Force:
  [
  F_z = -\frac{dU}{dz}= I_1I_2\frac{dM}{dz}
  ]
* Insert image current (I_2=-\eta I) at separation (z=2h):
  [
  F_\uparrow(h)=\eta I^2\left(-\frac{dM}{dz}\right)_{z=2h}
  ]

**D.2 Small-gap asymptotic**

* Show (-dM/dz \approx \mu_0 R/z) for (z\ll R)
* Derive:
  [
  F_\uparrow \approx \eta I^2\frac{\mu_0 R}{2h},
  \qquad
  I \approx \sqrt{\frac{2hmg}{\eta\mu_0R}}
  ]
* Include the numeric example (1 ton, (R=5) m, (h=2) cm → (I\sim 8) kA for (\eta=1)).

---

### Appendix E. Wake penetration depth and reflectivity (half-space match)

**E.1 Screening equation**

* Start from linear response:
  [
  \mathbf J_{\rm ind}=-\chi_T\mathbf A
  ]
* Derive Helmholtz form:
  [
  (\nabla^2-\lambda_T^{-2})\mathbf A=0,\quad \lambda_T^{-2}=\mu_T\chi_T
  ]

**E.2 Fourier mode matching**

* Above: decay rate (\alpha_1=k)
* Below: decay rate (\alpha_2=\sqrt{k^2+\lambda_T^{-2}})
* Reflection amplitude:
  [
  \eta(k)=\frac{\alpha_2-\alpha_1}{\alpha_2+\alpha_1}
  ]
* Dominant-mode estimate (k\sim 1/R)
* Inversion:
  [
  \lambda_T \approx \frac{1-\eta}{2\sqrt{\eta}},R
  ]
* Include a small table mapping (\eta\leftrightarrow\lambda_T) for a representative (R).

---

### Appendix F. Vortex-ring formulas and their mapping to toy-model parameters

**F.1 Vortex-ring impulse, energy, speed (thin core)**

* Impulse:
  [
  P_{\rm ring}\approx C_P\rho_{\rm eff}\Gamma\pi R^2
  ]
* Energy:
  [
  E_{\rm ring}\approx C_E\rho_{\rm eff}\Gamma^2 R,\ln!\left(\frac{8R}{a_{\rm ring}}\right)
  ]
* Speed:
  [
  U \approx \frac{\Gamma}{4\pi R}\left(\ln\frac{8R}{a_{\rm ring}}-\frac12\right)
  ]
* Show explicitly that (E/P \sim U) up to log/constant factors.

**F.2 Core swirl constraint**

* Near-core swirl speed:
  [
  v_{\max}\sim \frac{\Gamma}{2\pi a_{\rm ring}}
  ]
* Constraint (v_{\max}\lesssim \beta c) (choose (\beta=1) or conservative (\beta=0.1))
* Implication:
  [
  \Gamma \lesssim 2\pi a_{\rm ring}\beta c
  ]
* Combine with speed formula to get (U_{\max}(R,a_{\rm ring})).

---

### Appendix G. Ring emission rate derivation and scaling laws

**G.1 From thrust to ring rate**

* Hover: (mg=\dot N P_{\rm ring})
* Therefore:
  [
  \dot N = \frac{mg}{C_P\rho_{\rm eff}\Gamma\pi R^2}
  ]

**G.2 Eliminating (\Gamma) using target (U)**

* With (U \approx (\Gamma/(4\pi R))L_{\log}), obtain:
  [
  \dot N \approx \frac{mg,L_{\log}}{4\pi^2 C_P,\rho_{\rm eff},U,R^3}
  ]

**G.3 Clearance spacing**

* Define spacing:
  [
  d\sim \frac{U}{\dot N}
  ]
* Provide recommended “good” regimes (e.g., (\dot N\lesssim 50)/s per ton per nozzle).

**G.4 Energy-per-ring**

* With (\mathcal P\approx mgU):
  [
  E_{\rm ring}\approx \frac{\mathcal P}{\dot N}
  ]
* Include the kJ–100 kJ pulse scale examples.

---

### Appendix H. Anchoring (\rho_0) via defect mass and relating to (\rho_{\rm eff})

* Defect mass formula:
  [
  m_G=\kappa_m\rho_0\pi a^3\Lambda
  ]
* Solve:
  [
  a(\rho_0)=\left(\frac{m_G}{\kappa_m\rho_0\pi\Lambda}\right)^{1/3}
  ]
* Provide a short worked example with the “proton anchor” and (\rho_0=20,100,1000) kg/m³ giving atomic-scale (a).
* Discuss the assumption (\rho_{\rm eff}\sim \rho_0) and flag it as a key future kill-test.

---

### Appendix I. Compact “viability worksheet” (one page)

A plug-and-play checklist that future sessions can reuse:

**Inputs**

* payload (m)
* chosen hover power-per-mass (\mathcal P/m) (or choose (U))
* chosen ring radius (R)
* chosen core radius (a_{\rm ring}) (or core speed cap (\beta c))
* assumed (\rho_{\rm eff})
* constants (C_P\sim 1), (L_{\log}=\ln(8R/a_{\rm ring})-\tfrac12)

**Outputs**

* required (U=\mathcal P/(mg))
* required circulation (\Gamma=(4\pi R/L_{\log})U)
* check (v_{\max}=\Gamma/(2\pi a_{\rm ring})\le \beta c)
* ring impulse (P_{\rm ring}=C_P\rho_{\rm eff}\Gamma\pi R^2)
* ring rate (\dot N=mg/P_{\rm ring})
* spacing (d=U/\dot N)
* energy per ring (E_{\rm ring}=\mathcal P/\dot N)

---

### Appendix J. “Alive vs dead” decision checklist

A crisp list of conditions that must hold for free-space hover to remain plausible:

1. existence of a slow momentum carrier with (\varepsilon\sim U\ll c)
2. (\rho_{\rm eff}) not too small (quantify threshold for (\dot N\lesssim 50)/s per ton at chosen (R))
3. ability to generate (\Gamma\sim\mathcal O(1!-!10)) m²/s without catastrophic loss
4. coupling not identical to ordinary EM (or suppression mechanism) to avoid obvious signatures
5. stability of coherent packets and manageable plume interactions
