# Paper 2 Roadmap — Kinetic Plasma + Throat (4D) Extensions in the Brane–Bulk Toy Model

## 0. Goal and end-state

**Primary goal:** push the toy-model ontology as far as it can go toward **plasma + kinetic signatures**, while staying faithful to:

* **brane/bulk/throat ontology** (`brane_bulk_ontology.tex`),
* **field-theory EM mode** (`em_fields.tex`),
* **operator-correct localized non-ideality / gating** (`mhd.tex`),
* and the **hybrid/PN bookkeeping constraints** (`1pn_hybrid.tex`).

**End-state definition (what “done” looks like):**

1. A stable **brane EM field solver** (wave-equation completion) coupled to sources.
2. A kinetic “plasma” module (Vlasov or PIC) that reproduces canonical kinetic benchmarks:

   * Debye shielding, Langmuir oscillations,
   * Landau damping,
   * two-stream instability (and later Weibel / filamentation in ≥2D).
3. A controlled way to represent **4D throat/internal dynamics** that are partially hidden from 3D measurements, e.g. via a small set of throat-mode DOF per defect and/or an effective impedance/memory kernel.
4. Explicit, falsifiable **failure surfaces**: parameter regimes where the ontology or numerics stop being predictive.

---

## 1. Ontology recap (brane–bulk–throat)

### 1.1 3D effective decomposition of motion/wake

From `brane_bulk_ontology.tex`:

[
\mathbf{v}_{\mathrm{3D}} = \nabla\Phi + \mathbf{v}_T,
\qquad \nabla\cdot \mathbf{v}_T = 0.
\tag{BB-1}
]

Interpretation:

* **Longitudinal** part (\nabla\Phi): bulk-directed / compressive channel.
* **Transverse** part (\mathbf{v}_T): brane-parallel wake channel (the seed of the EM analogy).

### 1.2 Effective 3D continuity + Newtonian potential (gravity channel)

[
\partial_t \rho_{\mathrm{3D}} + \nabla\cdot(\rho_{\mathrm{3D}},\mathbf{v}*{\mathrm{3D}})=0,
\qquad
\nabla^2\Phi = 4\pi G*{\mathrm{eff}},\rho_{\mathrm{3D}}.
\tag{BB-2}
]

(We keep this present for bookkeeping and for later “grand-unification stress tests,” but Paper 2’s kinetic work is mostly EM-side.)

### 1.3 Brane vector potential as transverse wake (static / Poisson level)

From `brane_bulk_ontology.tex`:

[
\nabla\cdot\mathbf{A}=0,\qquad
\nabla^2\mathbf{A} = -\kappa_A,\mathcal{Q},\mathbf{u},\delta^{(3)}(\mathbf{r}).
\tag{BB-3}
]

and the identification
[
\mathbf{A} \equiv \kappa_A,\mathbf{v}_T.
\tag{BB-4}
]

with the Coulomb-like solution
[
\mathbf{A}(\mathbf{r}) = \frac{\kappa_A,\mathcal{Q}}{4\pi},\frac{\mathbf{u}}{r}.
\tag{BB-5}
]

### 1.4 Wake-mixing constraint from 1PN matching

From `brane_bulk_ontology.tex` and `1pn_hybrid.tex`:
[
\alpha^2 = \frac{3}{4}.
\tag{BB-6}
]

This is a structural constraint that tells us the toy ontology is already enforcing a particular longitudinal/transverse mixing that matches the EIH 1PN structure.

### 1.5 Finite-size / higher-order warning (2PN-ish)

From `1pn_hybrid.tex`:
[
\epsilon_{\text{fs}} \sim \left(\frac{a}{r_H}\right)^p.
\tag{H-1}
]

This is our reminder that once the dynamics probe scales (\sim a) (throat radius) or comparable, additional structure enters that is not captured by pointlike defect models.

---

## 2. EM equation library (dictionary mode vs field-theory mode)

Paper 2 should **commit** to field-theory mode for the EM sector, while treating dictionary mode as:

* a limiting case,
* a calibration/interpretation layer,
* and a sanity check for identities.

### 2.1 Dictionary mode (interpretation layer; not independent fields)

From `em_fields.tex`:

**Potentials from fluid variables**
[
\phi_{\mathrm{EM}} \equiv \lambda\left(h + \frac{1}{2}v^2\right),
\qquad
\mathbf{A} \equiv \lambda,\mathbf{v}.
\tag{EM-D1}
]

**Fields from potentials**
[
\mathbf{B} \equiv \nabla\times\mathbf{A},
\qquad
\mathbf{E} \equiv -\nabla\phi_{\mathrm{EM}} - \partial_t\mathbf{A}.
\tag{EM-D2}
]

**Homogeneous Maxwell identities (automatic from definitions)**
[
\nabla\cdot\mathbf{B}=0,
\qquad
\nabla\times\mathbf{E} = -\partial_t\mathbf{B}.
\tag{EM-D3}
]

**Key limitation:** if (\mathbf{A}\propto\mathbf{v}), the “EM” sector is not independent; you do not get genuine EM waves/topology/energy separate from the fluid.

### 2.2 Field-theory mode (independent gauge field)

From `em_fields.tex`:

**4-potential**
[
A^\mu=(\phi_{\mathrm{EM}}/c,,\mathbf{A}).
\tag{EM-F0}
]

**Flat d’Alembertian**
[
\Box ;\rightarrow; -\frac{1}{c^2}\partial_t^2 + \nabla^2.
\tag{EM-F1}
]

**Sourced wave equation**
[
\Box A^\mu = -\mu_0 J^\mu.
\tag{EM-F2}
]

**Lorenz gauge**
[
\partial_\mu A^\mu = 0
\quad\Leftrightarrow\quad
\frac{1}{c^2}\partial_t\phi_{\mathrm{EM}} + \nabla\cdot\mathbf{A}=0.
\tag{EM-F3}
]

**Component wave equations**
[
\Box\phi_{\mathrm{EM}} = -\mu_0 c^2\rho_e,
\qquad
\Box\mathbf{A} = -\mu_0\mathbf{J}.
\tag{EM-F4}
]

**Gauss law**
[
\nabla\cdot\mathbf{E} = \frac{\rho_e}{\epsilon_0},
\qquad \epsilon_0\equiv\frac{1}{\mu_0 c^2}.
\tag{EM-F5}
]

**Ampère–Maxwell**
[
\nabla\times\mathbf{B} - \frac{1}{c^2}\partial_t\mathbf{E} = \mu_0\mathbf{J}.
\tag{EM-F6}
]

### 2.3 Defect-based sources (particle/throat sources)

From `em_fields.tex`:

[
\rho_e(\mathbf{x},t)=q,\delta^3(\mathbf{x}-\mathbf{x}_d(t)),
\tag{EM-S1}
]

[
\mathbf{J}(\mathbf{x},t)=\rho_e(\mathbf{x},t),\mathbf{u}_d(t)
= q,\mathbf{u}_d(t),\delta^3(\mathbf{x}-\mathbf{x}_d(t)).
\tag{EM-S2}
]

**Effective charge from defect parameters**
[
q = \kappa_q,\rho_0\Gamma A = \kappa_q,\rho_0\pi a^2\Gamma.
\tag{EM-S3}
]

### 2.4 Lorentz/Magnus correspondence (force check)

From `em_fields.tex` (magnetic part written per unit length (L)):
[
\mathbf{f}_{L,\mathrm{mag}} = \frac{q}{L},\mathbf{u}\times\mathbf{B}.
\tag{EM-FOR1}
]

and its defect-expanded version
[
\mathbf{f}_{L,\mathrm{mag}} = \frac{\kappa_q\rho_0\pi a^2\Gamma B_0}{L},\mathbf{u}\times\hat{\mathbf{z}}.
\tag{EM-FOR2}
]

with the Magnus cross-term
[
\mathbf{f}_{M,u} = -\rho_0\Gamma,\mathbf{u}\times\hat{\mathbf{z}}.
\tag{EM-FOR3}
]

**Paper-2 usage:** we’ll keep the standard particle push equation (Lorentz force form) for kinetic simulations and use these as consistency/interpretation checks.

---

## 3. Operator-correct “gate” (localized non-ideality) — what carries over from `mhd.tex`

### 3.1 Scalar operator-correct diffusion lesson

When (\eta(\mathbf{x},t)) varies, the operator-correct diffusion form is divergence form:
[
\partial_t\psi = \nabla\cdot\big(\eta(\mathbf{x},t)\nabla\psi\big)
\quad\text{(preferred)}
\tag{G-1}
]
not (\eta\nabla^2\psi) (which introduces uncontrolled (\nabla\eta) artifacts).

### 3.2 Vector resistive operator lesson (MHD analogue)

The corresponding resistive operator is (schematically)
[
\text{(resistive term)};\propto; -\nabla\times\big(\eta(\mathbf{x},t),\nabla\times\mathbf{B}\big),
\tag{G-2}
]
which is the vector analogue that preserves clean dissipation structure.

### 3.3 How we use the “gate” in Paper 2

We will implement the gate primarily as one (or more) of:

* **Localized resistivity / non-ideal Ohm:** (\mathbf{E}+\mathbf{v}\times\mathbf{B}=\eta\mathbf{J}) inside a spacetime blob.
* **Localized kinetic collisions:** a collision operator (C[f]) active only in the gate region.
* **Localized phase-space diffusion:** a controlled proxy for turbulence/anomalous transport.

Paper 2 will keep the thesis: *gates are allowed, but only in operator-correct form, with explicit energy/entropy bookkeeping.*

---

## 4. Kinetic extension (how we get real kinetic plasma signatures)

### 4.1 Why kinetic requires velocity-space DOF

True kinetic signatures (Landau damping, phase-space trapping, two-stream, Weibel) require a distribution function (f_s(\mathbf{x},\mathbf{v},t)) or an equivalent particle ensemble.

### 4.2 Tier A (Rung A): 1D Vlasov–Poisson (fastest honest kinetic proof)

Use an electrostatic subset:
[
\partial_t f_s + v,\partial_x f_s + \frac{q_s}{m_s}E,\partial_v f_s = C[f_s],
\tag{K-A1}
]
[
E = -\partial_x\phi,\qquad
\partial_x^2\phi = -\frac{\rho_e}{\epsilon_0},\qquad
\rho_e=\sum_s q_s\int f_s,dv.
\tag{K-A2}
]

**Benchmarks (must reproduce):**

* Debye shielding,
* Langmuir oscillations,
* Landau damping rate (linear regime),
* two-stream growth and nonlinear saturation.

### 4.3 Tier B (Rung B): 2D/3D kinetic + EM fields (PIC or Eulerian Vlasov–Maxwell)

**Vlasov–Maxwell form:**
[
\partial_t f_s + \mathbf{v}\cdot\nabla_{\mathbf{x}} f_s

* \frac{q_s}{m_s}(\mathbf{E}+\mathbf{v}\times\mathbf{B})\cdot\nabla_{\mathbf{v}} f_s = C[f_s],
  \tag{K-B1}
  ]
  with moments
  [
  \rho_e = \sum_s q_s\int f_s,d^3v,\qquad
  \mathbf{J} = \sum_s q_s\int \mathbf{v}f_s,d^3v.
  \tag{K-B2}
  ]

**PIC equivalent:** particles ((\mathbf{x}_p,\mathbf{v}_p)) push with
[
\dot{\mathbf{x}}_p = \mathbf{v}_p,\qquad
\dot{\mathbf{v}}_p = \frac{q_s}{m_s}(\mathbf{E}+\mathbf{v}_p\times\mathbf{B}),
\tag{K-B3}
]
then deposit (\rho_e,\mathbf{J}) on the grid.

**Field equations (use Paper-2 EM field-theory mode):**

* either evolve (A^\mu) via (\Box A^\mu=-\mu_0 J^\mu) with gauge control,
* or evolve ((\mathbf{E},\mathbf{B})) with constraint control (e.g., Yee/CT).

**Benchmarks (additions):**

* magnetized particle motion (Larmor),
* Weibel/filamentation (anisotropic (f)),
* reconnection with kinetic effects (later).

### 4.4 Tier C (Rung C): optional bulk coupling (grand stress test)

Only after the brane kinetic-EM is validated:

* allow bulk channel to modulate effective inertia or background potential (slowly),
* track energy flow between brane EM, particle kinetic, and bulk/throat modes.

---

## 5. 4D throat/internal dynamics (hidden DOF that can make 3D look “unpredictable”)

### 5.1 Conceptual statement

In the ontology, EM happens in **4D in the throat** and is observed as a **3D projection/shadow** on the brane. That means a purely brane-local instantaneous closure (J^\mu=J^\mu[\text{brane vars}]) can be incomplete.

### 5.2 Minimal mathematical implementation (recommended)

Add a small set of internal mode coordinates (a_n(t)) per throat/defect (toy surrogate of throat cavity modes):
[
\ddot{a}_n + 2\gamma_n\dot{a}_n + \omega_n^2 a_n = g_n,\mathcal{S}[A^\mu,\text{mouth data}],
\tag{T-1}
]
where (\mathcal{S}) is a coupling functional at the throat mouth (e.g., depending on (\phi,\mathbf{A}), or local (\mathbf{E}\cdot\hat{n}), etc.).

Then allow polarization-like contributions to sources:
[
J^\mu = J^\mu_{\text{free particles}} + J^\mu_{\text{polar}}(a_n, A^\mu),
\tag{T-2}
]
which induces effective **frequency-dependent response** and **memory** when you integrate out the (a_n).

### 5.3 Diagnostics for “hidden channel”

* Energy budget closure improves when (a_n) are included.
* Resonant sweeps show phase lags/spikes in dissipation corresponding to (\omega_n).
* Brane-level EM response becomes dispersive in a controlled, model-predictive way.

---

## 6. What we still need to develop (math + numerics) and what can wait

### 6.1 Must-have for Paper 2 (non-negotiable)

1. **Brane EM solver in field-theory mode** (stable + constraint controlled).
2. **Charge conservation** at the discrete level (continuity check).
3. **Energy bookkeeping**: conserved in ideal runs; dissipates only through declared channels (gate, collisions, radiation).
4. A kinetic module (at least Tier A) that passes the canonical kinetic benchmarks.

### 6.2 Strongly recommended (to test “4D hidden interactions”)

5. Minimal throat internal-mode surrogate (few (a_n) oscillators) + coupling at the mouth.

### 6.3 Can wait (unless we push into these regimes)

* Full 2PN finite-size / intake dynamics: needed when (\lambda_D), skin depth, or interaction scales approach (a), or when near-field overlap makes pointlike approximations fail.
* Deep unification with bulk gravity channel: a later stress test.

---

## 7. Roadmap / milestones (start → finish)

### Phase 0 — Paper 2 framing and constants

**Deliverables:**

* A clear parameter dictionary: (\kappa_A,\kappa_q,\lambda, c) (or (c_{\rm EM})) and how they’re chosen.
* Explicit statement: dictionary mode vs field-theory mode.

### Phase 1 — EM core (field-theory mode)

**Implement:**

* Either (A^\mu) evolution: (\Box A^\mu=-\mu_0J^\mu) + Lorenz constraint control, or
* ((\mathbf{E},\mathbf{B})) evolution with (\nabla\cdot\mathbf{B}=0) control.

**Unit tests:**

* Vacuum wave propagation speed and dispersion.
* Gauge residual remains bounded (if using (A^\mu)).
* Gauss constraint satisfaction with known sources.

### Phase 2 — Particle/defect sourcing (non-kinetic baseline)

**Implement:**

* Point-defect sources (EM-S1–S3).
* Particle push in prescribed fields (sanity).

**Tests:**

* Field of a moving charge (qualitative Liénard–Wiechert-like behavior in the toy setting).
* Radiation/near-field behavior is numerically stable.

### Phase 3 — Kinetic Rung A (1D Vlasov–Poisson)

**Implement:**

* 1D Vlasov–Poisson (K-A1–K-A2) with optional gate collisions (C[f]) localized.

**Benchmarks:**

* Debye shielding.
* Langmuir oscillations.
* Landau damping (rate vs theory).
* Two-stream instability.

**Exit criteria:** the classic plots match and converge with resolution.

### Phase 4 — Kinetic Rung B (2D/3D PIC or Vlasov–Maxwell)

**Implement:**

* PIC push + deposit + EM fields, or Eulerian Vlasov–Maxwell.

**Benchmarks:**

* Weibel/filamentation (2D).
* Magnetized particle orbits (finite (\mathbf{B})).

### Phase 5 — Gate + topology tests (reconnection program)

**Implement:**

* Localized resistivity or kinetic collisions in a spacetime gate.
* Operator-correct forms (G-1, G-2).

**Benchmarks:**

* Reconnection: gate OFF vs ON.
* Connectivity change localized to gate.
* Energy dissipation localized to the declared channel.

### Phase 6 — 4D throat/internal modes (hidden channel)

**Implement:**

* Per-throat mode coordinates (a_n) (T-1) and source polarization (T-2).

**Benchmarks:**

* Resonant response sweeps.
* Energy closure improved and interpretable.
* Frequency-dependent effective response emerges.

### Phase 7 — Two-species plasma (first “full plasma” claim)

**Implement:**

* Two species (s\in{e,i}) in kinetic module.

**Benchmarks:**

* Charge separation effects.
* Skin depth / Hall-like behavior (as diagnostics, not necessarily as claims at first).

### Phase 8 — Optional: 2PN / finite-size corrections

Trigger this phase only if we intentionally push into (\lambda_D\sim a) or near-field overlap regimes.
Use `1pn_hybrid.tex` finite-size scaling (H-1) as the organizing parameter.

---

## 8. Diagnostics: what we must measure in every run

### 8.1 Conservation / correctness

* Discrete continuity residual: (\partial_t\rho_e+\nabla\cdot\mathbf{J}).
* (\nabla\cdot\mathbf{B}) residual.
* Gauge residual (if (A^\mu)).
* Total energy budget: field + particles + throat modes + declared dissipation.

### 8.2 Kinetic signature diagnostics

* (f(x,v,t)) snapshots (Tier A), phase-space holes/trapping (Tier A/B).
* Spectra: (\omega(k)) dispersion for Langmuir/EM waves.
* Growth/damping rates vs known linear theory.
* Debye length measurement.

### 8.3 Topology / reconnection

* Field-line connectivity diagnostics.
* Current sheet formation metrics.
* Dissipation localization to gate.

---

## 9. Expected failure surfaces (what to probe deliberately)

1. **Energy bookkeeping fails** unless throat modes are modeled → indicates real hidden bulk/throat channel.
2. **UV cutoff:** Debye length or skin depth (\lesssim a) → point-throat approximation breaks; 2PN/finite-size needed.
3. **Nonlocal / memory response:** brane-local closures can’t fit frequency response → need internal modes or retarded kernels.
4. **High defect density:** near-field overlap destroys linear superposition.
5. **Loss of universality:** constant re-tuning required across benchmarks → ontology not predictive in that regime.

---

## 10. Paper 2 proposed structure (write-up plan)

1. **Introduction + ontology** (brane/bulk/throat; why kinetic is the stress test).
2. **EM sector**: dictionary vs field-theory; adopt field-theory mode; relate to brane wake Poisson picture and its wave completion.
3. **Sources**: defects/particles; charge mapping (q) and calibration.
4. **Gate / non-ideality**: operator correctness; how (\eta(\mathbf{x},t)) or (C[f]) is implemented.
5. **Kinetic module**:

   * Tier A Vlasov–Poisson and benchmarks,
   * (optional) Tier B PIC/Vlasov–Maxwell.
6. **4D throat modes**: minimal oscillator surrogate; evidence for hidden channel effects.
7. **Results**: benchmark suite + reconnection + resonance sweeps.
8. **Failure surfaces**: explicit regimes of breakdown.
9. **Discussion**: what is derived vs assumed; roadmap to 2PN/finite-size and bulk coupling.

---

## 11. Practical build order (engineering)

**Module 1:** EM solver + constraint control + unit tests.

**Module 2:** Particle source + deposit + push + verification with analytic/synthetic tests.

**Module 3:** 1D Vlasov–Poisson (Tier A) and benchmark suite.

**Module 4:** 2D/3D PIC + EM (Tier B) and benchmark suite.

**Module 5:** Gate operators (resistivity/collisions) in operator-correct form + reconnection tests.

**Module 6:** Per-throat internal modes (a_n) + coupling + resonance sweeps.

**Module 7:** Two species + skin-depth/Hall diagnostics.

**Module 8 (optional):** finite-size / 2PN correction layer.

---

## 12. One-page “new session starter” checklist

* [ ] Decide EM evolution form: (A^\mu) vs ((\mathbf{E},\mathbf{B})).
* [ ] Fix mapping constants (\kappa_A,\kappa_q,\lambda, c_{\rm EM}).
* [ ] Implement EM core + vacuum wave + constraint tests.
* [ ] Implement particle sources (EM-S1–S3) + push.
* [ ] Implement Tier A Vlasov–Poisson and pass Landau + two-stream benchmarks.
* [ ] Add gate as localized (C[f]) or (\eta) with operator-correct form.
* [ ] Add internal throat modes (a_n) and test resonance/energy closure.
* [ ] Promote to 2D/3D PIC and re-run benchmarks.
* [ ] Add two species and run plasma benchmarks again.
