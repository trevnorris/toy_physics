# Two-Throat Simulation Handoff Specification
### A standalone, partially-closed problem statement for conditional electric and magnetic force-sign computations in a 4+1D superfluid analog

**Status: DRAFT v0 — 2026-07-19.** Authored as a standalone deliverable for two audiences: (a) a human computational physicist who wants to attempt the two-throat force computation from scratch, and (b) an AI-driven simulation ladder that will attempt it rung by rung. **Design contract:** the document consolidates every relation the committed sources presently fix and explicitly enumerates every action coefficient, coupling, boundary completion, and infrared choice that remains an OPEN input. The model does **not** yet have one fully closed interior action for the coupled `(χ, ρ, v, h, u)` throat system, so this is a standalone, **partially-closed handoff**, not a turn-key action. The upstream research audit trail (directives, sealed contracts, ledger stages) is *referenced by pointer* for provenance but is **never required reading** to identify either the committed equations or the missing closure data. Where a quantity is calibrated, postulated, or open, it is tagged as such inline.

Provenance tags used throughout: **[COMMITTED]** (a settled model commitment), **[DERIVED]** (computed and verified upstream), **[CALIBRATED]** (fixed by matching, not derived from first principles), **[POSTULATE]** (a declared new input), **[OPEN]** (undetermined by the committed model — a target of the simulation), **[CONVENTION]** (a units/frame choice).

---

## §0. Purpose, audience, and how to read this

### 0.1 What this document is for
This model is a toy **single-medium superfluid analog** in which gravity, light, electric force, and magnetism are meant to be four faces of one 4+1-dimensional compressible superfluid. Three of the four faces are built and calibrated (gravity, light) or structurally committed (the electric/magnetic *objects*). The one thing the committed model has been *proven not to determine* is the **sign of the electric and magnetic force between two charges** — because that sign runs through a boundary condition, `𝔅`, that the committed structure does not fix (see §0.3).

The purpose of this document is to hand a solver the **consolidated committed structure** for computing that sign by direct simulation of two throats — the model's word for two charged particles — after first solving the one-throat structure that the sign computation rests on. It is written so that the computation can be *attempted honestly*: the conditional `𝔅` branches are enumerated, and §2.6 names the additional OPEN closure data that a physicist must supply or a future simulation program must pin before a boundary-value problem can actually be instantiated. No conditional solve is allowed to masquerade as a derivation of its imposed inputs.

### 0.2 The two audiences and what each needs
- **A human computational physicist** needs: the committed governing sectors and the explicit OPEN closure register (§2), the conditional boundary-condition menu (§3), the frozen inputs and deferred quantities (§4), the measurement protocol that turns a solved field into a force sign (§5), and an honest account of what makes the full problem hard (§6). All of that is here without reference to the project's internal tooling. Supplying the OPEN closure card in §2.6 is a prerequisite to posing a concrete BVP.
- **An AI simulation ladder** additionally needs: a rung-by-rung build plan (§7) in which each rung becomes an independently-checkable numerical problem once its OPEN closure card is pinned, carrying the same anti-self-deception discipline that governed the upstream analytic work (manufactured-solution tests, convergence gates, dimensionless acceptance criteria, provenance guards, and fresh-agent recomputation). §7 supplies that conditional plan and forbids running past an unpinned coefficient.

### 0.3 The central open question, stated plainly
A particle is a **throat**: a puncture in the 3D brane whose wall structure extends into the 4th spatial dimension `w` on one side. "Which side" (`+w` or `−w`) is the sign of the charge. When a throat sits in the medium, the medium meets the throat's sleeve surface `Σ` under some boundary/constraint law `𝔅` — does the medium stick to the moving sleeve, slip past it, pass through it, lock to the brane's shear, or dissipate against it? **The committed model does not determine `𝔅`.** An exhaustive static adjudication (the "U2" analysis, referenced in §3) found that *every* route to the electric/magnetic sign passes through this undetermined `𝔅`, and that the committed structure leaves it open across all candidate classes. A simulation can therefore solve and test each imposed branch, and it can genuinely exclude a branch that admits no finite-energy, stable solution. It cannot infer merely from a conditional solution that nature enforces the imposed `𝔅_b`.

Accordingly this document does **not** pretend to a single answer. It organizes the two-throat problem as a small, finite family of **conditional branches** — one per candidate `𝔅` class — and asks the simulation to (i) test each instantiated branch for existence/stability and (ii) report the force sign *per surviving branch*. R1 may narrow the set by excluding inadmissible branches, but **unique selection of `𝔅` is not guaranteed by a one-body static solve**. If more than one branch survives, `𝔅` remains [OPEN] unless a downstream energetic or dynamical selection criterion is separately supplied and computed; that criterion itself requires the missing parent action/closure and is not smuggled into R1.

### 0.4 Scope and honest limits
- **In scope:** the static one-throat structure; the static two-throat force (→ electric sign); the moving one-throat response; the moving two-throat force (→ magnetic sign, and the test of whether a moving charge really carries a current `∝ V`).
- **Out of scope / deferred:** a full first-principles derivation of Newton's constant, the fine-structure constant, or the ~10⁴² gravity/EM hierarchy magnitude; these remain calibrated or open. The document targets the **form and sign** of the forces (target-blind, dimensionless), not their absolute magnitudes.
- **Genuine current limitation:** there is no single closed parent action for the coupled `(χ, ρ, v, h, u)` interior, sleeve, drain, and free boundary. Sections 2.6, 3.2, and 3.5 identify the missing coefficients/operators. A run is conditional on whatever values/forms it supplies for them; deriving those inputs is part of the unfinished program goal.
- **Falsification is a first-class outcome.** If a branch yields a null or wrong-sign result for the defect/antidefect force mechanism it is meant to realize, that is a genuine, welcome result that kills that mechanism on that branch — not a failure of the exercise (§5.4 names the kill).

### 0.5 Reading order
§1 (the model in one page) → §2 (equations) → §3 (the `𝔅` branch menu = the BC packet) → §4 (numbers) → §5 (how to measure a force sign) → §6 (why it's hard) → §7 (the simulation ladder). A physicist can stop after §6; the ladder consumer continues into §7.

### 0.6 Upstream provenance (reference only — not required)
The committed facts below trace to: `docs/conceptual_foundation.md` (the plain-language vision), `docs/em_u1_body_definition.md` v3.1 (the one-throat definition and the `𝔅` endpoint table), `docs/em_phaseC_force_decomposition.md` (the two-throat decomposition and rubric), and the `pathA_29/36/38/39` solver reports (gravity, light, static-electric, magnetism). The `𝔅`-undetermined verdict is the "U2" adjudication (directive + sealed production, wrap commit `5ceebb24`); the 16-branch conditional templates this BC packet is built from are the "v12 amendment" (wrap `6a6f317c`). None of these need to be opened to use this spec.

---

## §1. The model in one page

### 1.1 One medium
There is a single physical substance: a **compressible superfluid filling a 4+1-dimensional spacetime** — three ordinary space dimensions `x, y, z`, a fourth space dimension `w`, and time `t`. [COMMITTED] The bulk carrier is described by density `ρ(x,w,t)` and flow `v(x,w,t)` (equivalently a phase gradient) and is **inviscid and shear-free**. The wall/sleeve order `χ` and the reduced brane/throat coordinates `u` and `h` describe ordered configurations and collective modes of that same medium; they are not additional substances. There is no background force field added by hand: gravity is *flow*, light is a *ripple*, and electric/magnetic forces are the *static and moving structure of throats* in the medium.

### 1.2 The brane = our 3D space
Our observable 3-space is not the whole medium; it is a **brane**: a thin ordered-phase **domain wall** sitting at `w ≈ 0`, with bulk medium on both sides. [COMMITTED] The wall is where an order parameter (a "wall field" `χ`) interpolates between the two sides — a kink in `w`. Unlike the shear-free bulk, the **brane carries a rotational (MacCullagh) shear stiffness `μ_R`** [POSTULATE: stiffness value]: it resists in-plane twist. This single extra stiffness is what lets the brane carry transverse waves (light) that the bulk cannot.

### 1.3 A particle = a throat (a puncture)
A particle is a **throat**: a localized puncture of the brane whose wall structure does not simply close up but extends as a tube — the **sleeve** — into the bulk, on **one** side of the brane (either `+w` or `−w`). [COMMITTED] The sleeve is a tube of ordered-phase material of mouth radius `~a` and extent `~L` into `w`; the ratio `L/a = 37/20` [CALIBRATED: frozen branch ansatz, not self-selected]. Bundled with the sleeve is a trapped nonlinear core (a **geon**) that constitutes the particle's mass [POSTULATE: mass mechanism], and a bound flow field. Fundamentally the body is the field configuration `(χ,ρ,v)` at the puncture; the reduced variables `u` and `h` used below are collective projections of that same configuration, not extra body substances.

**Charge is a direction, not a winding.** [COMMITTED] The sign of the charge, `s = ±1`, is simply *which side* of the brane the sleeve extends into (`+w` vs `−w`). It is a topological/`Z₂` attribute (two disconnected choices), **not** an additive integer winding number. There is no native `U(1)` gauge charge and no native Gauss law in this model — that was computed, not assumed. The `s = −1` particle is *defined* as the `w`-mirror of the `s = +1` particle: `Φ₋(x, w) := R_w Φ₊(x, −w)`, where `R_w` is reflection `w → −w`.

### 1.4 The four forces, in the medium's own terms
All four are behaviors/collective sectors of the one medium; there is no separate "force field":

| Phenomenon | Native mechanism | Status |
|---|---|---|
| **Gravity** | The throat **drains**: a one-way radial inflow of medium toward the puncture. A second throat is carried along the inflow → attraction. The `1/r²` law is the far field of the drain through the finite-thickness brane slab. | [DERIVED + CALIBRATED] (`pathA_29`; GR-matched to 4PN by calibration) |
| **Light** | A **transverse in-plane shear ripple of the brane** (a MacCullagh rotational wave), propagating at `c_γ` with `c_γ² = μ_R / ρ_br`. Two transverse polarizations. The photon does not compress the medium. | [DERIVED on the postulated `μ_R`] (`pathA_36`) |
| **Electric force (static)** | The interaction of the two throats' **4D bodies beyond the mouths** — the sleeves/geon structure extending into `w`, mediated by a scalar (`h`) channel of the medium. Sign depends on `𝔅` (§0.3). | [COMMITTED object; sign OPEN] (`pathA_38`) |
| **Magnetism** | The **moving** 4D throat-body: the same sleeve structure in translation. Distinct from bulk vorticity / gravitomagnetism — it lives *in the throat*, not the bulk. Sign and current-law depend on `𝔅`. | [COMMITTED object; sign + current-law OPEN] (`pathA_39`) |

Mnemonic: **gravity = flow; light = brane shear; charge & magnetism = the throat (static vs moving).**

### 1.5 The three speeds (do not conflate)
The model has three distinct characteristic speeds, and confusing them is the classic error: **`c_s`** the bulk compression-wave speed (how fast a *change in the drain/gravity* propagates; `c_s ∝ ρ²`); **`c_γ`** the brane transverse-shear speed (the speed of light); and the field strength **`v_r`** (the local drain inflow speed — a field value, not a signal speed). Gravity is the *steady flow* between drains, not a compression ripple.

### 1.6 The open boundary fork that gates every downstream sign: `𝔅`
When the medium meets a throat's sleeve surface `Σ`, it obeys some boundary/constraint law `𝔅`. Five operationally-distinct **endpoint classes** span the committed possibilities (impermeable no-slip, impermeable free-slip, permeable, a nonholonomic "roll" lock to the brane shear, and a dissipative slip); mixtures are possible; and where `𝔅` physically lives (mouth collar vs bulk) is also open. `𝔅` is the decisive undetermined **boundary object** through which every downstream force sign runs; it is not the only OPEN coefficient in the unfinished coupled action (§2.6). §3 turns the candidate classes into a finite conditional branch menu. The frame throughout is the **medium rest frame** (asymptotic bulk at rest) [CONVENTION]; `V` is a throat's velocity relative to it — no Lorentz-covariant reinterpretation may be smuggled in by notation.

---

## §2. Consolidated brane+bulk equation set

This section assembles, in one place, the governing relations the committed sources actually fix: the bulk Madelung system, the quadratic brane sector, the reduced `h`/longitudinal operators, and the sector kernels. It is a **consolidation of frozen upstream results, not a re-derivation and not a closed coupled action**. Where the sources describe a relation in prose rather than writing the explicit PDE, the standard reduction is supplied and flagged as such (see §2.1). Section 2.6 then lists the action coefficients, cross-couplings, boundary completions, and infrared choices still [OPEN]. A solver may discretize the fixed sectors, but it cannot pose the coupled throat BVP until it supplies that closure card; this is a genuine limitation of the model, not an omitted implementation detail. The force sectors remain behaviours/collective coordinates of the **one medium**, not additional substances, and no Maxwell equations, `E`/`B` field pair, point source, or native current law is imported.

Symbols are collected in §2.6. The frame throughout is the **medium rest frame** (asymptotic bulk at rest) [CONVENTION] [from `docs/em_u1_body_definition.md:38`]; `V` is a throat's velocity relative to it.

---

### §2.1 The bulk medium (the one substance)

The single substance is a **compressible superfluid filling 4+1 dimensions** — three space dimensions `x, y, z`, a fourth space dimension `w`, and time `t` [COMMITTED] [from `docs/conceptual_foundation.md:102`; `docs/em_u1_body_definition.md:11`]. Its carrier is a complex order field

- **`ψ = √ρ · e^{iθ}`** (the GNLS carrier) [COMMITTED] [from `docs/medium_requirements_and_prior_art.md:141,145-146`],

with real number density `ρ(x, w, t)` and phase `θ(x, w, t)`. The flow is **potential flow set by the phase gradient**, `v = (ℏ/m) ∇θ` [COMMITTED] [from `docs/medium_requirements_and_prior_art.md:146` "flow / circulation / vortices"]. All four force phenomena are configurations and collective excitations of this one medium; there is no second substance.

**Governing equations (Madelung form of the GNLS carrier above).** The committed sources fix the carrier `ψ = √ρ e^{iθ}`, the quantum-pressure term, and the equation of state; the explicit continuity/Euler PDEs below are the **standard Madelung reduction of that stated carrier** (see caveat, end of §2.1), written out so the medium is discretizable:

1. **Continuity** (mass conservation of the one medium; divergence over all four space directions `x, y, z, w`):
   - `∂_t ρ + ∇·(ρ v) = 0`  [DERIVED — standard Madelung reduction of the committed carrier] [carrier from `docs/medium_requirements_and_prior_art.md:141,146`]
2. **Phase / Euler–Bernoulli (potential-flow momentum) equation**, carrying the internal energy from the equation of state plus the quantum-pressure (Madelung) term. With `μ(ρ) := U′(ρ) = (5K/4)ρ⁴` and `P(ρ) := ρμ − U = Kρ⁵`, the dimensionally consistent forms are
   - `ℏ ∂_t θ + (ℏ²/2m)|∇θ|² + μ(ρ) − (ℏ²/2m)(∇²√ρ)/√ρ = 0`,
   - equivalently, `∂_t v + (v·∇)v = −∇P/(mρ) + (ℏ²/2m²)∇[(∇²√ρ)/√ρ]`.
   These are [DERIVED — standard Madelung reduction] [carrier + quantum pressure from `docs/medium_requirements_and_prior_art.md:146`].
   - Gravity is the **steady flow** carried by this equation's Bernoulli pressure between draining throats, **not** a compression ripple [COMMITTED] [from `docs/conceptual_foundation.md:341-342`].
3. **Equation of state (single-well) and quantum pressure:**
   - internal-energy density `U(ρ) = (K/4) ρ⁵` — a **single well** (one stable bulk state) [COMMITTED] [from `docs/conceptual_foundation.md:165-166`; `docs/medium_requirements_and_prior_art.md:146`]
   - quantum-pressure (Madelung) energy `(ℏ²/8mρ)(∇ρ)²` [COMMITTED] [from `docs/medium_requirements_and_prior_art.md:146`]
4. **Compression-wave (bulk sound) speed** — the "speed of gravity" channel `c_s`:
   - `c_s² = 5 K ρ₀⁴ / m`, so **`c_s ∝ ρ²`** (density-dependent; slower in low density — the lensing/Shapiro mechanism) [COMMITTED] [from `docs/conceptual_foundation.md:122,334`]
5. **Shear-free property (load-bearing).** In the bulk the medium is essentially **inviscid and cannot sustain shear** — no sideways rigidity; it supports compression waves and potential flow only [COMMITTED] [from `docs/conceptual_foundation.md:117-119`]. This is exactly what protects magnetism (§2.4d) and confines light to the brane (§2.2).

**Three-speed caution (do not conflate)** [COMMITTED] [from `docs/conceptual_foundation.md:334-339`; see §1.5]: `c_s` = bulk compression-wave speed (`∝ ρ²`, how a *change* in the drain propagates); `c_γ` = brane transverse-shear speed = light (§2.2); `v_r` = local drain inflow speed (a *field strength*, not a signal speed). `v_r ≪ c_s` is weak gravity, not slow gravity.

> **Caveat (faithfulness).** Items 1 and 2 (the explicit continuity and Euler/Bernoulli PDEs) are **not quoted verbatim** from a committed source — the sources fix the carrier `ψ = √ρ e^{iθ}`, the quantum-pressure term, the EOS, and "potential flow / Bernoulli pressure" in prose, and the PDEs above are their canonical Madelung form. Items 3–5 (EOS, quantum pressure, `c_s`, shear-free) **are** quoted from committed sources.

---

### §2.2 The brane (our 3D space) and light

**The brane is a domain wall.** Our observable 3-space is a thin **ordered-phase domain wall** at `w ≈ 0` inside the bulk, with bulk medium on **both** sides [COMMITTED] [from `docs/conceptual_foundation.md:140-155`]. An order field (the **wall field `χ`**, a.k.a. `χ_B`) interpolates between the two degenerate bulk states — a **kink in `w`**; its energy per unit area is the surface tension [COMMITTED] [from `docs/conceptual_foundation.md:141-147`]. Because the single-well `U(ρ)` of §2.1 alone cannot form a wall, the two coexisting states come from the medium's added polar/smectic structure [COMMITTED] [from `docs/conceptual_foundation.md:165-168`; `docs/medium_requirements_and_prior_art.md:141-143`].

- **Slab caveat.** With the `χ_B ∈ [0,1]` double-well, "bulk on both sides" makes the brane a finite `χ_B = 1` **slab bounded by `χ_B = 0` bulk** — a **kink–antikink pair**, not a single wall; the slab width `W_slab` is an un-earned input (mapped onto the `L/a` self-selection open item) [COMMITTED + OPEN] [from `docs/conceptual_foundation.md:157-163`].

**MacCullagh rotational shear stiffness — the one extra property that makes light.** Unlike the shear-free bulk, the brane carries a **rotational (MacCullagh) shear stiffness `μ_R`** [POSTULATE: stiffness value] [from `docs/em_u1_body_definition.md:11`]. The stiffness is **curl-type**: energy `∝ (∇×u)²` in the in-plane displacement `u`, which is what yields exactly two transverse polarizations and no longitudinal photon (as opposed to Cauchy `(∂u)²` elasticity, which would give a spurious third longitudinal "photon") [COMMITTED] [from `docs/conceptual_foundation.md:365-369,385-389`].

**Primitive quadratic brane Lagrangian** (the frozen starting point of the light derivation; fields = in-plane displacement `u`, order-parameter phase `θ`, conjugate density perturbation `δρ_B`) [DERIVED starting point] [from `software/stage1_solver/reports/pathA_36_c5_phase_potential.md:29`]:

```
L = ½ ρ_br u̇²  −  ½ μ_R (curl u)²  −  ½ B (div u)²  +  J θ̇ δρ_B  +  ½ K_θ (grad θ)²
```

**Transverse sector — light [DERIVED, banked].** Each of the two transverse polarizations obeys [from `software/stage1_solver/reports/pathA_36_c5_phase_potential.md:108-114`]:

```
L_T = ½ ρ_br u̇_T²  −  ½ μ_R k² u_T²      →     ω² = (μ_R/ρ_br) k²,     c_γ² = μ_R / ρ_br
```

i.e. **two massless transverse photons at `c_γ² = μ_R/ρ_br`** — the model's light sector [DERIVED on the postulated `μ_R`] [from `pathA_36...:112-114`; `docs/conceptual_foundation.md:123,192`]. The photon does not compress the medium.

**Longitudinal sector — a characterized departure, NOT clean Maxwell [DERIVED verdict].** The longitudinal block does **not** close into a Maxwell Gauss chain on the committed (with-provenance) branch; the finite conjugate-density stiffness contributes a **Cauchy longitudinal modulus `B_eff = ρ_B0²/χ_c`**, leaving one physical longitudinal DOF (verdict `FAIL_CAUCHY_STRAY_LONGITUDINAL`) [DERIVED] [from `software/stage1_solver/reports/pathA_36_c5_phase_potential.md:10-23,33,18,64`]:

```
L_L = ½ ρ_br u̇_L²  −  C_J k u_L θ̇  +  ½ K_θ k² θ²  −  ½ B_eff k² u_L²,     with  C_J = −J ρ_B0,  B_eff = ρ_B0²/χ_c
```

This pathA36 block contains **only the longitudinal brane density mode `(u_L, θ)`**. It is neither the bulk GNLS sound mode `c_s` (an equality would require an [OPEN] coefficient bridge) nor the electric `h`-branon. Its extra longitudinal DOF is the computed departure from an isolated transverse-only light sector; it must not be tuned away to imitate Maxwell (that locus is `BY_TUNING`, not `WITH_PROVENANCE`) [from `software/stage1_solver/reports/pathA_36_c5_phase_potential.md:69-89`].

**Separate charge-coupled scalar `h` sector.** The `h`-branon is a distinct collective scalar of the same wall/throat medium, not `u_L` and not a second substance. At the reduced quadratic level, later source work assembles the scalar coordinates `Φ_s = (u_L,h)` as [DERIVED conditional reduced operator; coefficients partly OPEN] [from `software/stage1_solver/reports/pathA_39_stage4_field_classification.md:7-15,32-42,89-109`]:

```
L_h^(2)   = ½ M_h ḣ² − ½ K_h |grad h|²,                    K_h := M_h c_E²
L_mix^(2) = −C_hu grad(u_L)·grad(h)

Q_s(ω,k) = [ ρ_br ω² − B_eff k²      −C_hu k²       ]
           [ −C_hu k²                 M_h ω² − K_h k² ] .
```

Here `M_h > 0` is only a symbolic normalization, `c_E` is inserted from the pathA38 dynamic Green operator, and `C_hu` is sim-deferred. The stability condition for this reduced scalar block is `B_eff K_h − C_hu² > 0`. This block characterizes a possible reduced wave operator; it does **not** supply the missing parent `χ`/sleeve action or pin `M_h`, `K_h`, `c_E`, or `C_hu`. Those remain part of the closure card in §2.6.

---

### §2.3 The throat / sleeve / geon (a particle)

A particle is a **throat**: a localized puncture of the brane whose wall structure extends as a tube — the **sleeve** — into the bulk on **one** `w`-side [COMMITTED] [from `docs/em_u1_body_definition.md:12,22-25`; `docs/conceptual_foundation.md:427-429`]. The fundamental body is the field configuration `(χ,ρ,v)` at the puncture; `u` and `h` are reduced collective sectors of that configuration, not extra matter [COMMITTED body ontology] [from `docs/em_u1_body_definition.md:22-27`]:

- **Sleeve** — the wall-field `χ` wrapping into the bulk: a tube of ordered-phase material, mouth radius `~a`, extent `~L` in `w` on ONE side (`s = ±1`) [COMMITTED] [from `docs/em_u1_body_definition.md:23`].
- **Geon content** — the trapped nonlinear configuration constituting the mass; **its profile is a declared OPEN input** [POSTULATE: mass mechanism; profile OPEN] [from `docs/em_u1_body_definition.md:24`; `docs/conceptual_foundation.md:436`].
- **Bound flow field** — drain inflow + motion-induced flow [COMMITTED] [from `docs/em_u1_body_definition.md:25`].

**Puncture geometry.** Mouth radius `a`, `w`-extent `L`, with **`L/a = 37/20`** [CALIBRATED: a frozen branch ansatz, not self-selected] [from `docs/em_u1_body_definition.md:16`]. The throat radius cannot reach zero (a closed throat cannot drain) [COMMITTED] [from `docs/conceptual_foundation.md:427-429`]. As a finite-size matching estimate, the throat's nonzero size may be identified through the self-energy balance `k e²/a = m_e c²`, giving `a ∼ k e²/(m_e c²) = r_e` [CALIBRATED size comparison] [from `docs/conceptual_foundation.md:449`]. This is the throat's **finite-size resolution** of the divergent point-source idealization; it is not a point-source ontology and is never an interior-PDE or force-sign input.

**Charge = puncture direction, not a winding.** The sign `s = ±1` is simply **which side** (`+w` vs `−w`) the sleeve extends into — a topological/`Z₂` attribute, **not** an additive winding [COMMITTED] [from `docs/em_u1_body_definition.md:12`; `docs/conceptual_foundation.md:431-435,419-421`]. There is **no native `U(1)` gauge charge and no native Gauss law** — computed, not assumed (`NATIVE_P_NO_EMERGENT_GAUSS`) [COMMITTED — computed] [from `docs/em_u1_body_definition.md:15`].

**Body conjugation (definition).** The `s = −1` body ansatz is *defined* as the `w`-mirror of the `s = +1` ansatz:
- `Φ₋(x, w) := R_w Φ₊(x, −w)`,  with `R_w` = reflection `w → −w` [COMMITTED definition] [from `docs/em_u1_body_definition.md:29`].
- Whether the *ambient background* is itself `R_w`-symmetric is a **separate declared** [POSTULATE — U1 background]; the executable gravity slab (§2.4a) is one-sided, so ambient symmetry vs the one-sided slab is an [OPEN] seam [from `docs/em_u1_body_definition.md:31`].

---

### §2.4 The four sectors' governing operators (native terms)

Each sector is a distinct operator on the same medium. Gravity and light are built/calibrated; the electric and magnetic *objects* are committed and their localization/kernel forms are derived, but the **physical force sign is OPEN**, running through the sleeve-surface boundary/constraint operator `𝔅` (§0.3, §3).

#### §2.4a Gravity — the drain / brane-bulk return [DERIVED + CALIBRATED]

The throat **drains**: a one-way radial inflow of medium toward the puncture (the return is the separate global `S_leakage`, not a local throat loop) [COMMITTED] [from `docs/conceptual_foundation.md:203-204`]. The executable family and its far field [from `software/stage1_solver/reports/pathA_29_brane_bulk_return.md:7-18`]:

- **Geometry:** a flat finite slab, brane at `w = 0`, an adjacent return/absorber at `w = d`; the geometry is postulated, the return response **derived from projected continuity and the solved Helmholtz transport phase** [from `pathA_29...:7-8`].
- **Transport amplitude:** `T0(0) = 1/(ε0 + 1)`; the residual is bounded but lower order: `p_res(ℓ0) = 1`, `p_res(ℓ1) = 3` [from `pathA_29...:9`].
- **Far field = `1/r²` via the zero mode:** on the admissible DC-sink completions (`destructuring_absorbing`; `bloch_stack`), the solved compact-cell spectrum + 3D radial equation give **`p = 2`** (i.e. `1/r²` gravity) [DERIVED] [from `pathA_29...:11-13,27`].
- **Drain admissibility / signed source:** `Z = −M0·ε0/(ε0 + 1)`, with `Z < 0` the drain-admissibility premise [from `pathA_29...:18`].
- Anti-localizing warp `μ(w) = exp(2 k_warp w)` → non-normalizable zero mode, `p = 3`, `RETURN_NOGO` (the able-to-fail control) [from `pathA_29...:24-25`]. GR-matched to 4PN **by calibration** [CALIBRATED; see §1.4].

#### §2.4b Light — brane MacCullagh shear [DERIVED]

The transverse wave operator and `c_γ² = μ_R/ρ_br` are as in §2.2 (two transverse polarizations) [DERIVED] [from `software/stage1_solver/reports/pathA_36_c5_phase_potential.md:108-114`]. No separate operator is added here; light **is** the brane's transverse MacCullagh shear ripple.

#### §2.4c Static electric — the `±w` throat body [DERIVED object; physical sign OPEN]

The static electric interaction is the two throats' 4D bodies beyond the mouths, mediated by a scalar (`h`) channel of the medium, brane-localized so it comes out `1/r²` in 3D [COMMITTED object] [from `docs/conceptual_foundation.md:207,444-451`]. The frozen localization/kernel result [from `software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md:6-16`]:

- **Transverse wall-spectrum eigenproblem** (solve first; only then the radial factor): `O f = m² K_∥ f` [DERIVED] [from `pathA_38...:7`].
- **Zero-mode norm:** `N0 = 8/(3ℓ)`, from the integral density `2/(ℓ² cosh⁴(w/ℓ))` [DERIVED] [from `pathA_38...:8`].
- **Radial exponents:** `p_static = 2` and `p_dynamic = 2` → **`1/r²` Coulomb** [DERIVED] [from `pathA_38...:9`].
- **Compact source projections:** `q_h(+) = 2 Q_E tanh(b/ℓ)/b`, `q_h(−) = −2 Q_E tanh(b/ℓ)/b`; the pure-even gravity overlap and the orthogonalized no-monopole source are both `0` [DERIVED] [from `pathA_38...:11-12`].
- **Projected Green kernel:** `3 Q_E² ℓ tanh²(b/ℓ) / (8π R b²)` [DERIVED] [from `pathA_38...:13`].
- **`±w` sign matrix:** `U₊₊ = +3 Q_E² ℓ tanh²(b/ℓ)/(8π R b²)`, `U₊₋ = −3 Q_E² ℓ tanh²(b/ℓ)/(8π R b²)` [DERIVED relative-sign structure] [from `pathA_38...:16`].
- **Able-to-fail controls:** delocalized → `1/r³` (`p = 3`); ghost → instability; no-monopole; pinned-branon; Yukawa (all `FAIL_*`) [from `pathA_38...:19-24`]. `Q_E` is [CALIBRATED] and deferred [from `pathA_38...:30`].

`O` and `K_∥` in the pathA38 eigenproblem are symbolic operator/weight names in this handoff, not a supplied nonlinear `χ`/`h` parent PDE. The reported zero mode, projections, and kernel are fixed reduced results; constructing the interior throat operator that produces them requires the §2.6 wall/scalar/action closure [OPEN].

> The kernel above fixes a reduced relative-sign structure; the **absolute physical attract/repel sign is OPEN** — the U2 adjudication found that the total force sign runs through `𝔅` and the sleeve microstructure [OPEN; see §0.3 and §1.4].

#### §2.4d Magnetism — a conditional moving-source benchmark; sign + current-law OPEN

Magnetism is meant to be the **same** 4D throat-body under motion: an `O(V)` one-body source whose two-body exchange begins at `O(V₁V₂)`, charge-coupled and `±w`-sensitive; it is distinct from weak, mass-coupled gravitomagnetism [COMMITTED object; current law OPEN] [from `docs/conceptual_foundation.md:249-257`]. The frozen pathA39 Stage-2 kernel is **not a native derivation of that source**. It is [DERIVED conditional on the imported `j∝sV` ansatz — the ansatz R2/R4 must replace] [from `software/stage1_solver/reports/pathA_39_magnetic_force.md:87-103`]:

- **Imported source premise:** `j_T = q_A^T V`, `j_L = q_L V`, with `q_A^T = Nu·aT·s` and `q_L = Nu·aL·s`. Thus the calculation assumes both linearity in `V` and charge-oddness `∝s`; it does not compute them.

- **Static-projector operators and their Green inversions** [from `pathA_39...:11-14`]:
  - `O_T = k² μ_R`,  `G_T = 1/(k² μ_R)`
  - `O_L = k² ρ_B0²/χ_c`  (`= k² B_eff`),  `G_L = χ_c/(k² ρ_B0²)`
- **Exchange integral** (a genuine propagator inversion; sign and falloff fall out, not asserted) [from `pathA_39...:18-22`]:
  - `U_12 = −∫ d³k/(2π)³  j₁(−k) · G(k) · j₂(k)`,  with `F⁻¹[1/k²] = 1/(4π R)`.
- **Compact kernels** (`R = X₂ − X₁`, `n = R/R`, `D = V₁·V₂`, `A = (V₁·n)(V₂·n)`) [from `pathA_39...:27-32`]:
  - `K_T = (D + A)/(8π R)`,  `K_L = (D − A)/(8π R)`
  - `U_12 = −s₁ s₂ Nu²/(8π R) · [ aT² (D+A)/μ_R + aL² (D−A)/B_eff ]`,  `F_12 = −∇_R U_12`.
- **Measured falloff:** potential `R⁻¹`, body-center force `R⁻²` [DERIVED, measured not supplied] [from `pathA_39...:45`].
- **Conditional sign diagnostic:** `NO_CANCELLATION_BOTH_CHANNELS_ATTRACTIVE` — under the imported source premise, both stable channels are attractive for like sources; outward radial force = repulsion, inward = attraction [DERIVED conditional on imported `j∝sV`] [from `pathA_39...:47-55`].
- **Deferred amplitudes and frame:** `aT`, `aL` remain [OPEN]; the result is preferred-frame unless `c_E = c_γ` [from `docs/conceptual_foundation.md:296`]. The algebraic consistency residual `U_eom − U_integrate = 0` closes **within the conditional ansatz** [from `pathA_39...:43`].

> This kernel is a **post-source benchmark only**. It is barred from the ancestry of R2 and R4 until R2 independently computes the one-body moving source from the sleeve response. If that source is natively `∝sV`, the kernel becomes an eligible comparator; if the source is a departure or zero, the benchmark does not override it. The physical magnetic sign and current law remain `𝔅`-gated and [OPEN] (§0.3, §5.2, §7.2–§7.4).

---

### §2.5 The one-body force law (the frame the throat's dynamics live in)

The body is an **open system** (one-way intake; `𝔅` includes non-variational endpoints), so a single effective-action Hessian is **not** a well-posed definition. The one-body law is a **four-channel force balance over the collective coordinates `q_A ∈ {X, p}`** [COMMITTED] [from `docs/em_u1_body_definition.md:56-68`]:

```
F_A^var  +  F_A^flux  +  F_A^𝔅  +  F_A^rad  =  0,     with   F_A^𝔅 = F_A^constraint(𝔅) + F_A^Rayleigh(𝔅)
```

- **`F^var`** — variational channel, from a time-density Lagrangian `L_eff` on a declared control volume `Ω_c` plus surface actions on `∂Ω_c` and the sleeve surface `Σ`; supplies the inertia tensor `M_AB(𝔅)` and symplectic/Berry form `Ω_AB(𝔅)` (no `ρa³L` or `E/c²` pre-asserted) [from `docs/em_u1_body_definition.md:62`].
- **`F^flux`** — control-surface channel: momentum/mass/energy carried across `∂Ω_c` by drain and return; a `∝V̇` contribution is reported in `M_AB` **or** the intake coefficient `C_ṁ`, never both [from `docs/em_u1_body_definition.md:63`].
- **`F^𝔅` = `F^constraint` + `F^Rayleigh`** — the boundary/constraint channel: multiplier reactions for constraint endpoints (the E4 nonholonomic shear-lock reaction lives here; E1 no-slip lives here or in `F^var` depending on imposition), and the genuinely dissipative Rayleigh channel (E5) [from `docs/em_u1_body_definition.md:64-67`].
- **`F^rad`** — radiative channel (hereditary in general); the local form `F = M V̇ + Ω V` is explicitly the **sub-radiative slow-acceleration truncation**, not the definition [from `docs/em_u1_body_definition.md:68`].

**The tilt observable `p`** [POSTULATE] [from `docs/em_u1_body_definition.md:81-86`]:
- `p` = the in-brane projection of the sleeve's tangent = the **leading in-brane multipole of the sleeve `χ`-configuration**; `p` itself is the collective orientation coordinate (a signed `θ = s·p` is banned — it would insert charge-oddness by convention) [from `docs/em_u1_body_definition.md:81-83`].
- **Parity is computed, not assumed:** under body conjugation + ambient `R_w`, `p(−s, V) = p(+s, V)` — charge-**even**; if the response is analytic and the in-brane background isotropic, the leading nonzero `p(V)` is **linear in `V`** [from `docs/em_u1_body_definition.md:84`].
- **Tilt alone cannot establish a current** — current-likeness is decided only by the channel-by-channel coupling package; the honest outcomes are `O(V) ∝ sV` (the imported `j = sV` becomes computed), a characterized departure, or zero (the magnetism source dies natively) [from `docs/em_u1_body_definition.md:85-86`].

---

### §2.6 Dimensional conventions, OPEN action coefficients, and one-throat unknowns

**Unit / dimensional convention.** The upstream `pathA` reports carry all quantities **symbolically** (no numerical unit system is imposed); dimensional homogeneity is checked per-expression. The parameter register uses a **`(L, T, M)`** length–time–mass dimension triple to tag every quantity [CONVENTION] [from `research/pde_ledger_v2/notes/parameter_register.md:762`]. The frame is the medium rest frame (§2.1) [from `docs/em_u1_body_definition.md:38`].

**Parameter symbols used above** (collected for the calibration table, §4): medium `ℏ, m, K, ρ₀`; brane `ρ_br, μ_R, B, J, K_θ, C_J, ρ_B0, χ_c, B_eff (= ρ_B0²/χ_c)`; scalar sector `M_h, K_h, c_E, C_hu`; throat/sleeve `a, L (L/a = 37/20), ℓ, b, Q_E, s = ±1`; gravity slab `d, ε0, M0, T0`; conditional magnetism benchmark `Nu, aT, aL`; speeds `c_s, c_γ, v_r`. The **calibratable subset** of these (the frozen geometry, its derived dependents, and the calibrated physical bridges) is tabulated in §4; the [OPEN] closure-card coefficients introduced below have **no committed values** — where a few of them appear in §4 (e.g. `M_h`, `K_h`, `C_hu`, `g_nat`, `K_m`) they are shown *only* as OPEN placeholders, and the remainder (`Z_χ`, `κ_χ`, `V_χ`, `κ_bend`, `κ_anchor`, `g_χh`, `g_χdrain`, …) are not tabulated at all.

**Closure status — no closed parent action exists.** The committed sources fix the bulk Madelung equations (§2.1), the quadratic brane Lagrangian and its transverse/longitudinal sectors (§2.2), the reduced `h`/`u_L` operator (§2.2), and the sector kernels (§2.4). They do **not** fix one variational/non-variational system that couples those pieces to a free sleeve, geon, drain, and `𝔅`. Before any branch can be instantiated, the physicist or simulation campaign must provide a versioned **closure card** containing the following [OPEN] data. The symbols in parentheses are bookkeeping slots only; they do not assert a functional form or value.

| OPEN closure item | What must be supplied before posing the BVP |
|---|---|
| `χ` kinetic coefficient (`Z_χ`) | The time-kinetic normalization and any field dependence for the wall/sleeve order variable `χ`. |
| Wall stiffness / wall potential (`κ_χ`, `V_χ`) | The spatial gradient stiffness, double-well/polar-smectic potential, slab-width control, and their relation to the committed surface tension. `μ_R` is the brane's rotational shear stiffness and does not fill this slot. |
| Sleeve shape energy (`κ_bend`, `κ_anchor`) | Bending, anchoring/collar, surface-tension, and free-boundary coefficients or functionals that determine `Σ`; these are required for a free-geometry solve. |
| `χ↔h` coupling (`g_χh`) | The actual local/nonlocal coupling and parity by which sleeve orientation sources the distinct `h` scalar. The reduced source projection in §2.4c does not supply this parent coupling. |
| `χ↔drain` / geon coupling (`g_χdrain`) | The coupling of the wall/sleeve and geon to `(ρ,θ)` plus the throat sink `S_drain`, the global return `S_leakage`, and the associated mass/momentum/energy source terms. |
| Scalar-wave coefficients | `M_h`, `K_h` (or `c_E` with `K_h=M_hc_E²`), `C_hu`, and any nonlinear completion; require `M_h>0` and `B_effK_h−C_hu²>0` on a stable reduced branch. |
| EOS / wave-operator completion | The bulk `U(ρ)=(K/4)ρ⁵` and quantum pressure are fixed in §2.1, but the coupled wall/throat EOS, nonlinear `χ`, `h`, `u`, geon wave operators, interface/jump laws, and any coefficient bridges among `c_s`, `c_γ`, and `c_E` are [OPEN]. |
| Boundary-operator completion | The branch-specific traction, permeability/texture law, E4 constraint map, E5 Rayleigh law/coefficient, and mixture implementation listed precisely in §3.2. |
| Mouth-ensemble boundary data | Not merely "source vs defect" (§5.3(a)) but the **exact boundary functional** it implies: for the **fixed-source** ensemble, the imposed scalar-source flux datum on the mouth annulus (the Neumann/Robin source term and its magnitude); for the **fixed-defect** ensemble, the signed geometric/winding datum held fixed (the `±w` defect displacement profile and its topological charge) plus the natural (stress-free) condition on the conjugate variable. The two give opposite stiffness-dependence of the charge coupling; the committed model does not settle which applies, so the *specific* functional is a per-branch [OPEN] boundary-data input — a datum to supply, not a knob to tune the sign. |
| IR / return scheme | One- vs two-sided domain, return/absorber law, far-field truncation/extrapolation, treatment of zero modes/moduli, and (for moving/radiative rungs) the outgoing/absorbing prescription. |

This list is an **input register, not a license to fit the sign**. A chosen closure card makes a run conditional on that card; a future derivation may replace a card entry with a committed coefficient. If a required entry is absent, the correct rung status is `UNRESOLVED(OPEN closure item)`, not a numerical result.

**The eight OPEN one-throat field profiles (conditional BVP unknowns).** Once the closure card is supplied, a one-throat solve can produce the **tilted-sleeve field profiles**. They remain [OPEN] and unsolved at the present model state; the native tilted-sleeve BVP is their named discharge route, but that BVP is itself gated on both `𝔅` and the closure data above. Each enters the embedding as `δΦ(y) = p_i · profile_i(y; endpoint, ambient)` [OPEN — reduction debt; symmetry class computed in production, not assumed] [from `research/pde_ledger_v2/notes/parameter_register.md:749-777`]:

| # | Field profile (BVP unknown) | dims `(L, T, M)` | domain |
|---|---|---|---|
| 1 | `indexed_density_tilt_profile` | `(−4, 0, 0)` | scalar density tangent on `Ω_c` |
| 2 | `indexed_flow_tilt_response` | `(1, −1, 0)` | vector velocity tangent on `Ω_c` |
| 3 | `indexed_h_tilt_profile` | `(0, 0, 0)` | scalar `h` tangent on `Ω_c` |
| 4 | `indexed_phase_tilt_profile` | `(0, 0, 0)` | scalar phase tangent on `Ω_c` |
| 5 | `indexed_shear_tilt_profile` | `(1, 0, 0)` | brane displacement tangent on `Ω_c` |
| 6 | `indexed_sleeve_surface_normal_profile` | `(0, 0, 0)` | normal variation on `Σ` |
| 7 | `indexed_sleeve_tilt_profile` | `(0, 0, 0)` | `χ`-field tangent on `Ω_c` |
| 8 | `indexed_uw_tilt_profile` | `(1, 0, 0)` | normal displacement tangent on `Ω_c` |

[table from `research/pde_ledger_v2/notes/parameter_register.md:762-771`]

These eight are what a **closure-complete one-throat solve would compute** (§7 rung R1); until then they are counted reduction debt, and any output whose content descends entirely from them lands `UNRESOLVED` (the anti-padding rule) [from `research/pde_ledger_v2/notes/parameter_register.md:775-777`].

---

## §3. The `𝔅` branch menu = the boundary-condition packet

This section is the operational branch packet. Because the committed model does not determine `𝔅` (§0.3), the two-throat computation is not one problem but a **finite family of conditional templates**, one per candidate boundary-operator class. The packet enumerates the committed endpoint character and states exactly which mathematical pieces are still missing. It does not pretend that the endpoint prose is already a closed operator. The upstream static adjudication established that these are the candidate classes and excluded none; it did not prove that every class is physically admissible or uniquely selected.

### 3.1 What a branch is (and what it is not)
A **branch** is a hypothesis of the form *"suppose the medium meets the sleeve surface `Σ` under boundary operator `𝔅_b`, in ambient `A`."* Solving a branch tells you the physics **conditional on that hypothesis**. It does **not** assert that nature enforces `𝔅_b`: there is no "16 branches ⇒ 16 admissible boundary conditions" reading. A finite-energy/stability failure can legitimately exclude a branch and narrow the surviving set. Conversely, existence of a conditional solution is not evidence that nature chose its imposed operator. Unique selection may require a downstream energetic/dynamical criterion built from a future closed parent action and may remain [OPEN]. Throughout, branch labels are **non-evidential tags**: they name a hypothesis, carry a stable identifier for cross-referencing, and never count as evidence for selection. [COMMITTED discipline; from the U2 static/dynamical scope line.]

**Shared committed sectors, branch-specific boundary.** Every branch uses the same committed sector equations in §2 and the same selected §2.6 closure card; branches then differ in `𝔅_b` and ambient `A`. Because the coupled parent residual is [OPEN], the present model does not yet justify the stronger claim that a fully assembled interior discretization already exists and only a boundary swap remains.

### 3.2 The candidate axis: 9 operator classes
The boundary operator's committed possibilities are spanned by **five endpoint classes** `E1–E5`, plus **three mixture families**, plus a catch-all. The table gives the committed endpoint character and, in the last column, the exact closure still required:

| Class | Committed endpoint character on `Σ` | Native-root class | OPEN mathematical completion |
|---|---|---|---|
| **E1** impermeable + no-slip | `v·n̂ = V·n̂` and `v_∥ = V_∥` | holonomic kinematic BC | The free-surface geometry and its compatible `χ,h,u` interface/jump conditions from the parent action. |
| **E2** impermeable + free-slip | `v·n̂ = V·n̂`; tangential traction is zero | holonomic normal BC + natural traction | The stress/traction tensor derived from the closed action; “stress-free” is not computable before it is supplied. |
| **E3** permeable / phase-texture | medium passes through the translating `χ` texture; drain flux unimpeded | interface/texture branch | Permeability/phase-texture law, `χ` interface and jump conditions, and drain-flux/source matching. “No extra BC” is not a closure. |
| **E4** nonholonomic shear-lock ("roll") | throat velocity locks to brane transverse shear at the collar | multiplier / virtual-work reaction | The exact velocity-level map `g(V,u̇_T;χ,Σ)=0`, its collar support, multiplier pairing, and reaction-force functional. |
| **E5** dissipative slip | tangential slip with loss | non-variational Rayleigh channel | The Rayleigh functional, positive dissipation/slip coefficient or kernel, and traction–slip law. |
| **MIXTURE(F_E1_E4)** | E1 ∧ E4 | holonomic + nonholonomic | The completed E1 and E4 operators imposed as a positive conjunction, with compatibility checked. |
| **MIXTURE(F_E2_E4)** | E2 ∧ E4 | holonomic/natural + nonholonomic | The completed E2 and E4 operators imposed as a positive conjunction, with compatibility checked. |
| **MIXTURE(F_E4_E5)** | E4 ∧ E5 | nonholonomic + Rayleigh | The completed E4 and E5 operators imposed as a positive conjunction, with compatibility checked. |
| **OTHER** (catch-all) | residual complement | — | No defining operator; remains `UNRESOLVED(named operator definition/closure)`. |

Notes: [COMMITTED] E4 (the nonholonomic "roll" lock of the throat to the brane's shear) is the ontology's decisive "roll-vs-slip" pivot and is **distinct** from the fluid no-slip E1 — do not conflate them. The mixture families are the *only* admissible combinations the committed constraint structure generates (a positive conjunction of orthogonal boundary components); free-parameter blends (weighted `ζ`-mixtures) are **not** admissible. **OTHER has no operator to pose** and therefore **receives no BVP template** below; it is carried as the open record `UNRESOLVED(named operator definition/closure)` — a reminder that basis closure over `{E1–E5, mixtures}` is itself [OPEN].

### 3.3 The ambient axis: 2 closures
Independently of `𝔅`, the throat sits in one of two **ambient backgrounds**, which change the domain and the far-field matching (not the boundary operator):
- **`one_sided_pathA29`** [COMMITTED; gravity-sector executable]: the finite slab actually used to build gravity — brane at `w = 0`, a return/absorber at `w = d`, medium on the `0 < w < d` side **only**. It is asymmetric in `w`; its gravity transport result does not close the coupled throat action. **Because medium exists on only one side, this ambient can host only a `+w` sleeve** (a throat extending into `0 < w < d`); a `−w` throat (`s = −1`) has no medium to extend into. It therefore supports the single `+w` one-throat solve (R1, `s = +1`) and the like-charge `++` two-throat, but **cannot pose any `−w` configuration** — see the four-orientation constraint below (§3.4).
- **`two_sided_R_w_postulate`** [POSTULATE]: the conceptual `w`-symmetric background (medium on both sides of the brane; a `+w` and a `−w` throat see equivalent environments). This is a *postulated* ambient; a branch in this ambient is **doubly conditional** — on `𝔅_b` and on the `R_w` postulate — and can never, by itself, count as evidence.

### 3.4 The 16-branch packet
Crossing the **8 operator-bearing candidate classes** (E1–E5 + the three mixtures; OTHER excluded — no operator) with the **2 ambients** gives **16 conditional branch templates**. Each becomes a distinct tilted-sleeve BVP only after the common action/PDE closure card and the class-specific completion in §3.2 are supplied:

```
Branch identifier          𝔅 class            ambient
────────────────────────────────────────────────────────────────
U2:E1:one_sided_pathA29     E1                 one-sided slab
U2:E1:two_sided_R_w_postulate  E1              two-sided (postulate)
U2:E2:one_sided_pathA29     E2                 one-sided slab
U2:E2:two_sided_R_w_postulate  E2              two-sided (postulate)
U2:E3:one_sided_pathA29     E3                 one-sided slab
U2:E3:two_sided_R_w_postulate  E3              two-sided (postulate)
U2:E4:one_sided_pathA29     E4                 one-sided slab
U2:E4:two_sided_R_w_postulate  E4              two-sided (postulate)
U2:E5:one_sided_pathA29     E5                 one-sided slab
U2:E5:two_sided_R_w_postulate  E5              two-sided (postulate)
U2:MIXTURE(F_E1_E4):one_sided_pathA29      F_E1_E4    one-sided slab
U2:MIXTURE(F_E1_E4):two_sided_R_w_postulate F_E1_E4   two-sided (postulate)
U2:MIXTURE(F_E2_E4):one_sided_pathA29      F_E2_E4    one-sided slab
U2:MIXTURE(F_E2_E4):two_sided_R_w_postulate F_E2_E4   two-sided (postulate)
U2:MIXTURE(F_E4_E5):one_sided_pathA29      F_E4_E5    one-sided slab
U2:MIXTURE(F_E4_E5):two_sided_R_w_postulate F_E4_E5   two-sided (postulate)
```
(`OTHER` × 2 ambients = the 2 named-open records, no BVP.)

**Ambient ↔ orientation-census constraint (a genuine restriction, not a modeling choice).** Because the one-sided `pathA29` ambient hosts only `+w` sleeves (§3.3), its eight branches support the one-throat solve (R1, `s = +1`) and the like-charge `++` two-throat configuration only; the `{+−, −+, −−}` orientation pairs each require a `−w` throat and are **impossible** there. The full four-orientation census `{++, +−, −+, −−}` — and therefore the **electric sign**, which is intrinsically a like-vs-unlike *comparison* (§5.2) — is posable **only in the two-sided `R_w` ambient** (medium on both sides). Two consequences the consumer must carry: (i) the headline static-electric-sign result (rung R3) is **conditional on the `R_w` two-sided-ambient postulate** — the one-sided ambient can deliver R1 profiles and the `++` like-charge force magnitude, but not the sign comparison; (ii) any future work that argues the physical ambient really is one-sided would leave the electric *sign* undefined in this construction until the `−w` sector is given an admissible embedding. This does not change the 16-branch enumeration (each branch is still a well-defined template); it restricts which *rung* each ambient can serve.

### 3.5 The conditionally posed tilted-sleeve BVP template (per branch)
For each of the 16 branches, the committed structure supplies the following **template**, not yet a closed BVP:

- **Unknowns (the "tilt profiles").** The eight sleeve field profiles the solve produces — the throat's static deformed shape [OPEN, this is exactly what the solve computes]:
  `indexed_density_tilt_profile`, `indexed_flow_tilt_response`, `indexed_h_tilt_profile`, `indexed_phase_tilt_profile`, `indexed_shear_tilt_profile`, `indexed_sleeve_surface_normal_profile`, `indexed_sleeve_tilt_profile`, `indexed_uw_tilt_profile`. (These are the register's `GAP_OPEN_FIELD_PROFILE` leaves; they are **undetermined by any committed algebra** and are the reason a solve is needed.)
- **Committed interior pieces.** The static limit of the bulk Madelung system (§2.1), the brane quadratic sectors and distinct reduced `h`/`u_L` operator (§2.2), and the sector operators in §2.4, in the branch's ambient domain.
- **OPEN interior completion (required).** Supply every §2.6 closure-card entry: the `χ` kinetic/stiffness/potential terms; sleeve bending/anchoring/free-boundary energy; `χ↔h`, `χ↔drain`, geon, and source/return couplings; scalar/EOS/wave-operator completion; and all coefficients. Without these, there is no shared coupled Euler–Lagrange residual to discretize.
- **Boundary operator on `Σ`.** Supply the branch's exact mathematical completion from §3.2: E1's compatible interface conditions; E2's action-derived traction; E3's texture/permeability/jump law; E4's explicit `g` and multiplier reaction; E5's Rayleigh law/coefficient; or the completed conjunction for a mixture.
- **Ambient / asymptotic matching.** The far-field/return matching of the branch's ambient (§3.3): decay into the one-sided slab toward the `w = d` absorber, or `R_w`-symmetric decay on both sides for the two-sided postulate.
- **Zero-mode treatment.** The translational Goldstone (the sleeve's in-brane translation) is a degenerate direction of the interior Hessian; the closure card must choose and write an exact phase/modulus condition (for example, a declared orthogonality condition to the translational tangent) or the solve is singular. The need to quotient is [DERIVED]; the exact quotient is an IR-scheme input.
- **Well-posedness.** [OPEN] until the interior and boundary completions are supplied. Numerical convergence can test an instantiated BVP; it cannot establish well-posedness for an operator that has not been defined.

**Instantiation rule and honest limitation.** A concrete BVP record must reproduce the complete §2.6 closure card, the exact §3.2 endpoint completion, the domain/interface/jump/free-boundary laws, and the zero-mode/IR prescription before mesh generation. Until then its status is `UNRESOLVED(closure)`. Once instantiated, the eight profiles remain genuine solver unknowns; no profile, eigenvalue, response coefficient, branch admissibility, or branch selection is assumed. This explicit closure debt is the standalone handoff: a physicist can see precisely what must be supplied, but the committed model does not yet supply it.

### 3.6 How the packet feeds the force computation
After the closure gate, rung R1 (§7) attempts each branch's tilted-sleeve BVP and reports the eight profiles and admissibility **per branch**. Those profiles feed the one-body parity census (§5.1) and, in pairs, the two-body kernels and four-channel decomposition (§5.2). R1 may exclude branches with no converged finite-energy/stable solution, but it does not infer that a surviving imposed `𝔅_b` is nature's unique choice. Electric and magnetic force signs are therefore reported per surviving branch; unique physical selection remains [OPEN] unless a separate energetic/dynamical criterion is later defined and run.

---

## §4. Calibration and ansatz table — inputs, derived dependents, and deferred outputs

The two-throat computation is **target-blind and dimensionless**: it aims at the *form and sign* of the forces, never their absolute magnitude (§0.4). This section keeps three categories separate: (a) frozen calibration/ANSATZ inputs, which a run must declare and cannot claim to predict; (b) algebraic dependents of those inputs; and (c) quantities that only a fuller free-geometry, moving, or radiative solve could output. A wrong held-out *dimensionless* ratio is fatal; a wrong dimensionful absolute (`G`, `c`, `ℏ`) is calibration convention, not falsification.

### 4.1 Frozen branch ANSATZ inputs and rung ownership
`L/a = 37/20` and `ε_r = 1/20` are frozen clean-rational **branch ANSATZ inputs**, explicitly not self-selected [CALIBRATED placeholders]. A static R1 run with this geometry uses them; it does not reproduce them. Letting the sleeve geometry vary and dynamically self-select its aspect/healing ratios is a harder, sim-deferred free-geometry problem that additionally requires the OPEN sleeve action in §2.6. [Provenance: `research/pde_ledger/redteam_adversarial/ANSATZ_LEDGER.md` §1.]

| Symbol | Value | Meaning | Honest simulation role |
|---|---|---|---|
| `L/a` | `37/20 = 1.85` | throat axial extent / mouth radius (aspect ratio Λ*) | **R1 frozen input**; a future free-geometry solve may replace it with an output [ANSATZ §1 row 1]. |
| `r_F1` | `√(4107−100π²)/(10π) ≈ 1.77799` | normalized hybridization | Algebraic dependent of the input `L/a`, not an independent R1 prediction [row 1a]. |
| `r_c^F1` | `= r_F1² ≈ 3.16126` | Family-1 dependent | Algebraic dependent of `r_F1` [row 1b]. |
| `ε_r` | `1/20 = 0.05` | wall/healing width ratio | **R1 frozen input**; sets `ℓ = ε_r a`; a future free-geometry solve may replace it with an output [row 2]. |
| `ℓ` | `= ε_r·a` | support/healing scale | Algebraic dependent under the frozen R1 geometry [row 2a]. |
| `c0` | `3/4` (with `c1 = 1/4`) | conservative quadrupole ansatz | A moving/quadrupolar **R2/R4 ansatz input**, not a static R1 observable; a fuller moving solve could test or replace it [row 3]. |
| `χ_Q` | `1` | outgoing/retarded `ℓ=2` normalization | A radiative **R4 ansatz input** unless the outgoing solve computes it; never an R1 target [row 4]. |
| `L_W` | `:= L` | mixed-tube ↔ throat identification | Frozen identification that helps define `r_F1`; not independently confirmed by the same constructed geometry [row 7]. |
| `eps_2`, `eps_4` | `0` (posited) | geometry-contamination numbers | [OPEN]; report sensitivity or let a fuller free-geometry solve determine them [row 8]. |
| `g_nat`, `K_m` | branch ansätze (O(1)) | normalized mouth coupling; wall Robin closure | Entries of the OPEN closure card, not outputs of a BVP that already assumes them [rows 5, 6]. |

**What R1 can and cannot check.** With frozen `L/a` and `ε_r`, R1 can check residuals, boundary/interface conditions, finite energy, stability, parity, profile convergence, and sensitivity to explicitly declared closure inputs. It cannot count the imposed geometry as a reproduced result, cannot test a quadrupole/radiative normalization absent from its equations, and cannot self-select the throat dimensions. No entry in this table may be varied after seeing a force sign.

### 4.2 The calibrated physical bridges (dimensionful — convention, not prediction)
These fix the map from model units to physical units. Getting them "wrong" is not falsification; they are calibrated once and carried. [Provenance as cited.]

| Bridge | Relation | Status |
|---|---|---|
| Newton's constant `G` | fixed by matching the drain far-field to the PN/GR ladder | [CALIBRATED; `GENUINE_BLOCKED` from first principles — a from-throat derivation is sim-deferred] |
| Speed of light `c_γ` | `c_γ² = μ_R / ρ_br` (brane transverse-shear speed) | [DERIVED on the POSTULATED brane stiffness `μ_R`; `pathA_36:112-118`] |
| GR radiation normalization `γ_GR` | `2G / (5c⁵)` (quadrupole reaction) | held-out benchmark the 2.5PN sector matches [Peters 1964 / Maggiore; `benchmarks.yaml:20-42`] |
| Static Coulomb magnitude | projected Green kernel `∝ 3 Q_E² ℓ tanh²(b/ℓ)/(8π R b²)`, zero-mode norm `N0 = 8/(3ℓ)` | [CALIBRATED] magnitude/deferred `Q_E`; the `1/r²` form and `±w` relative-sign matrix are reduced targets (`pathA_38:8-16,30`). |
| Longitudinal brane stiffness | `B_eff = ρ_B0² / χ_c` | [DERIVED from the pathA36 finite-compressibility branch; `pathA_36:10-23,33`] |
| `h`/moving-sector coefficients | `M_h`, `K_h=M_hc_E²`, `C_hu`, `aT`, `aL` | [OPEN] sim-deferred parent-action and source coefficients; the pathA39 moving kernel is conditional on the imported current ansatz (§2.4d). |

### 4.3 Held-out benchmarks (the dimensionless judges)
The only first-class falsifiers are **held-out dimensionless ratios**. The two-throat exercise's own such ratio is the **electric/gravitational force ratio `F_e/F_g`** between two elementary throats — currently a *fit*, and the sharpest both-sectors test the simulation could eventually make target-blind. The proton/electron mass ratio `m_p/m_e = 1836.15267343` [CODATA 2018; `benchmarks.yaml:3-19`] is the standard external anchor. The ~10⁴² gravity/EM hierarchy magnitude is **[OPEN]** — explicitly *not* a target of this exercise, which fixes only signs and forms.

---

## §5. The two-throat measurement protocol

This section is the payoff: it turns a *solved field configuration* into a *force sign*. Given a converged one-throat solution (§7 rung R1) and a converged two-throat solution on a chosen `𝔅` branch (§3, §7 rungs R3/R4), the procedure here emits a single sentence — *"on this branch the force is attractive/repulsive, the charge-coupled sign is X, and the throat-current mechanism is confirmed/killed."* It is deliberately a **re-presentation of the committed force-decomposition rubric**, lifted out of the upstream process artifacts (`docs/em_phaseC_force_decomposition.md`, `docs/em_u1_body_definition.md`) into a form a computational physicist can run with a solver in hand. Nothing here presumes a sign; the whole point is that the sign is *measured off* the solved field, not assumed.

**Discipline carried throughout (do not violate).** [COMMITTED] The throat is an *open, driven* object: there is no closed point particle and, in general, no scalar pair-potential. The complete two-body force is the sum of four channels
```
F_pair = F_var + F_flux + F_𝔅 + F_rad
```
— variational (conservative), flux (momentum through the mouths), constraint/boundary reaction, and radiative — and a drop in one channel can be outweighed or reversed by another. The energy/potential language below (`U_ab`) captures *only* the conservative `F_var` piece and may be collapsed to a scalar potential only after an integrability check passes; some channels have no potential at all. [src: `docs/em_phaseC_force_decomposition.md:84-90, 127-142`] The frame is the **medium rest frame** (asymptotic bulk at rest) [CONVENTION]; `V` is a throat's velocity relative to it. No Lorentz-covariant reinterpretation and no imported Maxwell relation (in particular **no `j ∝ sV` inserted by hand**) may enter — every current-like structure is an *output* to be classified, never an input. [src: `docs/em_phaseC_force_decomposition.md:144-154`]

The protocol has four steps, run in order:
- **§5.1** — the one-body parity census (do this on a *single* solved throat first);
- **§5.2** — the four-channel decomposition and the two two-body kernels (the force-sign read-off, static → electric, moving → magnetism + the `j∝V` test);
- **§5.3** — the mouth-ensemble fork, the required ablations, and the predeclared gates (what proves a computed response is not vacuous);
- **§5.4** — the evaluation rubric and the **named KILL** (the predeclared landings, including the welcome null/wrong-sign outcome).

---

### §5.1 Protocol step A — the one-body parity census

Before any two-body work, characterize *one* solved throat. This is the "step A" one-body deliverable; it is a prerequisite, and a one-body coupling to the bulk is **not yet a force** — a two-body force needs a genuine cross term between two induced configurations (a shared mediator or a shared boundary/constraint reaction). [src: `docs/em_phaseC_force_decomposition.md:124-125, 242-244, 249-251`]

On the converged stationary one-throat field, compute:

1. **The post-mouth parity structure.** Split every field/profile into orientation-even and orientation-odd parts on a genuinely `R_w`-symmetric background:
   ```
   Φ_s = Φ_even + s·Φ_odd ,   v_s = v_even + s·v_odd ,   s ∈ {+1,−1}
   ```
   where `v_even` is the common radial **drain** (the gravity channel) and `v_odd` is the `R_w`-odd axial post-mouth component (the charge channel). Do this for each of: the axial drain, the sleeve/`χ` continuation, the wall displacement, and the surface flux across the control surface `∂Ω_c`. [src: `docs/em_phaseC_force_decomposition.md:100-110, 249-251`]

2. **Tag every parity claim by provenance.** Mark each `even`/`odd` result as **body-only**, **ambient-postulate-dependent** (relies on the postulated two-sided `R_w`-symmetric background), or **one-sided-asymmetry-map** (an artifact of the one-sided slab that the executable gravity solve `pathA_29` actually computes). This tag is load-bearing: a nonzero `s`-odd component that is really ambient asymmetry is *not* a charge effect. [src: `docs/em_phaseC_force_decomposition.md:138-143, 249-251, 331`]

3. **The tilt `p` and its velocity response `p(V)`.** [POSTULATE: the tilt observable] `p` is the in-brane projection of the sleeve's tangent (equivalently the leading in-brane multipole of the sleeve `χ`-configuration) — the collective orientation coordinate. Record the steady `p(V)` (distinct from the dynamic susceptibility `p(ω,V)`). **Parity is computed, not assumed:** under body conjugation plus ambient `R_w`, `p(−s,V) = p(+s,V)` — `p` is charge-**even**; and if the steady response is analytic and the in-brane background isotropic, the leading nonzero `p(V)` is **linear in `V`**. [DERIVED conditional on analyticity] [src: `docs/em_u1_body_definition.md:83-84`]

4. **The mediator couplings, channel by channel.** For each native collective sector — bulk GNLS compression `(ρ,θ)` at `c_s`; the distinct brane longitudinal mode `(u_L,θ)` with its own coefficient-dependent speed; brane transverse shear `u_T` at `c_γ`; the separate charge-coupled scalar `h`; wall/`χ` modes; and internal geon modes — record the one-body source strength and its `R_w` transformation. No bridge among those scalar speeds is assumed. **Crucial subtlety carried into step §5.2:** a charge-**even** tilt `p` can still produce a charge-**odd** mediator coupling once the mediators' own `R_w` transformations enter (even × odd = odd), and vice versa. The current-likeness verdict is therefore *never* read off `p` alone — it is decided channel-by-channel by the coupling package. [src: `docs/em_u1_body_definition.md:85`, `docs/em_phaseC_force_decomposition.md:239-240`]

**Gate for step A:** the census must be complete and every parity claim tagged before the two-body kernels are attempted. An untagged `s`-odd component silently contaminates the electric channel.

---

### §5.2 The two-body kernels and the four-channel decomposition (the force-sign read-off)

#### (a) The four force cells

Every two-body interaction organizes into a 2×2 grid of **{static / moving} × {charge-even / orientation-pair}**. This grid *is* the map from the solved field to the four force phenomena. [src: `docs/em_phaseC_force_decomposition.md:100-118`]

| field overlap | `s`-parity | velocity order | classification | status |
|---|---|---|---|---|
| radial / common drain | even | `V⁰` | **gravity** | [DERIVED + CALIBRATED] (`pathA_29`) |
| oriented post-mouth body/flow | `s1·s2` | `V⁰` | **static electric** candidate | [sign OPEN] |
| moving common drain | even | `V1·V2` | **gravitomagnetic** candidate | [sign OPEN] |
| moving oriented post-mouth body | `s1·s2` | `V1·V2` | **EM-magnetic** candidate | [sign OPEN] |

Naming discipline (enforce it): a *stationary* orientation-dependent attraction/repulsion is the **static electric** candidate — **not** magnetism. A term earns the name **magnetic** only if it is generated by *moving* the oriented throat-body, is bilinear in the two bodies' velocities at leading order, and **vanishes when either velocity is zero**. [src: `docs/em_phaseC_force_decomposition.md:119-125`]

#### (b) The static two-body kernel → the electric sign

Compute the two-throat interaction at zero velocity for **all four** orientation pairs `{s1,s2} ∈ {++, +−, −+, −−}`, using *one* fixed action/stress/flux prescription. **(This full census requires both `±w` throats, hence the two-sided `R_w` ambient; the one-sided `pathA29` ambient hosts only `+w` throats and can supply the `++` cell alone — §3.4. The electric-sign comparison is therefore `R_w`-postulate-conditional.)** Because the throat is open/driven, the primary object is the **channel-resolved force**, decomposed per channel `ch ∈ {var, flux, 𝔅, rad}` via the four-orientation **Hadamard decomposition**:
```
F^ch_{A,ab}(R) = (1/4) Σ_{s1,s2}  s1^a · s2^b · F^ch_A(s1,s2) ,   a,b ∈ {0,1}
  (a,b)=(0,0)   charge-even / mass channel   → gravity cell
  (1,0),(0,1)   should VANISH on a symmetric, exchange-symmetric background
  (1,1)         orientation-pair channel      → static-electric candidate
```
The `s1·s2` structure is what isolates the pieces: summing the like-pairs (`++`, `−−`, where `s1·s2 = +1`) against the unlike-pairs (`+−`, `−+`, where `s1·s2 = −1`) projects out the `(1,1)` orientation-pair channel and cancels the charge-even `(0,0)` mass/gravity channel. A nonzero `(1,0)` or `(0,1)` component is a **diagnostic of ambient asymmetry** (the seam between the postulated two-sided background and the one-sided slab), **not** a clean charge interaction — flag it, do not fold it into the electric channel. [src: `docs/em_phaseC_force_decomposition.md:127-143`]

**Reading off the electric sign.** Collapse the `(1,1)` channel to a scalar potential `U_11(R)` *only* after its conservative (`var`) contribution passes an integrability check (flux/constraint/radiative pieces may have no potential). [src: `docs/em_phaseC_force_decomposition.md:138-142`] Then classify the asymptotic orientation energy:
- **`U_11(R) ~ +A_q/R` with `A_q > 0`** — like-pairs (`s1·s2 = +1`) gain energy as they approach → **repel**; unlike-pairs lose energy → **attract**. This is the **Coulomb-sign landing**, sitting alongside the gravity cell `U_00(R) ~ −A_g/R`. It may be written down **only after** demonstrating (a) a gapless localized mode, (b) a nonzero monopole, and (c) the native sign from the full calculation. [src: `docs/em_phaseC_force_decomposition.md:71-79`]
- **`A_q < 0`** — the anti-Coulomb sign (like-pairs attract). This is a genuine, reportable outcome, **not** to be softened.
- **`A_q = 0`** (or no `1/R` tail — local healing at `R≈a` without a long-range tail) — no static electric force on this branch. A KILL-adjacent landing (see §5.4). [src: `docs/em_phaseC_force_decomposition.md:78-79`]

The magnitude and sign of `A_q` are [OPEN] under the committed model — they are exactly what the branch solve computes.

#### (c) The moving two-body kernel → the magnetic sign and the `j∝V` test

Translate the *complete* one-body solution at velocity and solve the induced two-body response; expand the interaction and read off the tensor structures **without inserting a current**:
```
U_12 = U_mass(R)                              # charge-even, V⁰  (gravity)
     + s1·s2·U_orientation(R)                 # orientation-pair, V⁰  (static electric)
     + V1_i·K_mass_ij(R)·V2_j                 # charge-even, V¹V¹  (gravitomagnetic)
     + s1·s2·V1_i·K_orientation_ij(R)·V2_j    # orientation-pair, V¹V¹  (EM-magnetic candidate)
     + other one-sided / parity-mixed structures + higher orders
```
[src: `docs/em_phaseC_force_decomposition.md:144-154`]

For each predeclared velocity geometry, repeat the solve for **all four orientation pairs** `{++,+−,−+,−−}` and compute **each force channel** `{var,flux,𝔅,rad}` before summing. Apply the same Hadamard projection in §5.2(b) channel by channel at fixed `(V₁,V₂)`; the moving `(1,1)` coefficient is the orientation-pair candidate. Parallel/antiparallel runs alone do not replace the four-orientation census. **(As in §5.2(b), the `s1·s2` census needs both `±w` throats → the two-sided `R_w` ambient; the magnetic sign is therefore `R_w`-postulate-conditional exactly as the electric sign is — §3.4. The one-sided ambient can supply the moving `++` cell only.)**

**Reading off the magnetic sign.** The magnetic candidate is the **`s1·s2·V1·V2` kernel `K_orientation`**. Its sign is the magnetic sign; it must satisfy the naming discipline — bilinear in the two velocities and **vanishing when either `V = 0`** (verified by the §5.3 ablation that removes translational velocity). [src: `docs/em_phaseC_force_decomposition.md:122-123, 262-263`]

**The `j∝V` test (does a moving charge carry a current?).** [OPEN] The desired `j ∝ s·V` behavior is an **output classification of `K_orientation`**, *not* an input — it may instead be zero, nonlocal, anisotropic, higher-order in `V`, or carry extra tensor structure. [src: `docs/em_phaseC_force_decomposition.md:152-154`] Execute the test by chaining the one-body census into the coupling package:
1. take the measured tilt response `p(V)` from §5.1 with its computed parity and velocity order; the conditional charge-even/linear expectation is not substituted for the result;
2. propagate it through each mediator coupling **channel by channel**, applying the mediator's own `R_w` transformation — a charge-even `p` yields a charge-**odd** coupling wherever the mediator is `R_w`-odd (even × odd = odd), and vice versa; [src: `docs/em_u1_body_definition.md:85`]
3. classify the resulting `O(V)` orientation-channel coupling as one of: **∝ `s·V`** (the imported current becomes a *computed* result), **characterized departure** (a different tensor structure), or **zero** (this magnetism source dies natively). All three are first-class. [src: `docs/em_u1_body_definition.md:86`]

Do **not** project onto `s·V` until the full response exists; classify `K_orientation` as exactly convection-like / convection-like-plus-departures / departure-only / null / unresolved-on-a-named-missing-closure. [src: `docs/em_phaseC_force_decomposition.md:255-258`]

---

### §5.3 The mouth-ensemble fork, required ablations, and predeclared gates

#### (a) The load-bearing fork — which mouth "ensemble" the throat imposes

This is the choice that most changes the answer, and it **interacts with `𝔅`** — it is keyed per branch, not chosen globally. [OPEN] The question is *what the throat actually holds fixed at its mouth*: [src: `docs/em_phaseC_force_decomposition.md:172-183`]
- **Fixed conventional source** (a normal, externally-imposed scalar source): the coupling response scales `~ J²/K`. A *stiffer* wall gives a **weaker** coupling.
- **Fixed displacement / winding / topological defect** (a fixed signed geometric boundary condition): the stored energy scales `~ K·q_top²`. A *stiffer* wall gives a **stronger** coupling.

The two ensembles give **opposite stiffness-dependence of the coupling character**: `J²/K` weakens with stiffness, while `Kq_top²` strengthens. This does **not** imply that the force sign must flip. Run fixed-source ↔ fixed-defect as a diagnostic exchange and report the response magnitude, sign, channel composition, and stiffness scaling in both cases. A sign or hierarchy that does flip is conditional on unresolved mouth physics; conversely, a response that is invariant even in the coupling character/stiffness dependence that the committed boundary structure says should change is evidence that the adjudication is vacuous. The committed model does not settle which ensemble applies (`pathA_38`'s reduced source projection does not decide it). [src: `docs/em_phaseC_force_decomposition.md:178-183, 268-273`]

#### (b) Required ablations (step D) — each able to fail at its own assert

A computed response is not credited until an ablation proves it is not vacuous — removing the structure that supposedly generates it must **flip or null** the response: [src: `docs/em_phaseC_force_decomposition.md:259-263`]
- ablate the common drain / mass source;
- ablate the `w`-odd sleeve / orientation structure;
- ablate the post-mouth axial flow;
- ablate the `h` coupling;
- ablate the brane transverse shear;
- ablate the longitudinal compression;
- ablate the E4 constraint (the nonholonomic roll-lock);
- ablate the outer control-surface flux / return.

Two ablations are mandatory and load-bearing: **the magnetic candidate must vanish when the orientation structure is removed**, and **removing translational velocity must eliminate every magnetic term while leaving the static orientation term intact**. [src: `docs/em_phaseC_force_decomposition.md:262-263`]

#### (c) Physical comparison (step E)

A Maxwell landing needs more than "`1/R²` + attraction." Check: vanishing at `V=0`; reversal under one *derived* current reversal; transverse/Darwin tensor structure; the **relative coefficient** to the static orientation force; propagation poles/speeds; longitudinal/scalar residue; preferred-frame pieces; and the dependence on the OPEN boundary/closure choices. [src: `docs/em_phaseC_force_decomposition.md:264-267`]

#### (d) Predeclared gates — each able to fail

- **Topology gate.** Is `±w` a genuine topological charge/winding, or only a `Z₂` orientation? Test whether the two sectors are disconnected, whether a finite-energy interpolation exists, and whether opposite throats can annihilate. Until shown, call them *conjugate geometric-defect candidates*, not established defect/antidefect topology. [src: `docs/em_phaseC_force_decomposition.md:306-309`]
- **Range gate.** Classify the potential explicitly — `1/R`, Yukawa, bulk `1/R²`, dipole, or null. Do not assume the desired `1/R` potential; local healing at `R≈a` is not a long-range tail. [src: `docs/em_phaseC_force_decomposition.md:310-311`]
- **Stability gate.** A correct sign obtained via a ghost, negative-norm mode, or unstable background is a **failure**, not a success. [src: `docs/em_phaseC_force_decomposition.md:312-314`]
- **Momentum-closure gate.** Verify total body-plus-medium momentum balance. In an open system the two body forces need **not** be equal-and-opposite on their own, but any imbalance must be carried by a **named** flux/return/radiation channel, never left implicit. [src: `docs/em_phaseC_force_decomposition.md:315-317`]

#### (e) The hierarchy read-off (step F) — from the *same* solved throat

If both channels land on the same `1/R` potential, report
```
F_electric / F_gravity = A_q / A_g
```
with the throughput `ṁ`, the signed boundary datum, the wall stiffness, the zero-mode normalizations, and the return closure all fixed by the *same* native throat solution — **no new parameter inserted after either force is read off.** [OPEN] `A_q/A_g` is **not** a universal constant; it carries the two bodies' effective masses and charge strengths, so the target is the *family* of ratios (electron pair `~4×10⁴²`, proton pair `~10³⁶`), not one number. If the two channels do **not** land on the same `1/R`, do not report a single ratio — compare the far radial-force coefficients directly and report the result as **distance-dependent / non-Coulombic.** [src: `docs/em_phaseC_force_decomposition.md:184-200, 268-273`]

---

### §5.4 The evaluation rubric and the named KILL

The protocol must accept every outcome, but evaluation must still be willing to conclude an outcome is *too structurally different to call electromagnetism.* The standing risk is the unfalsifiability / Lorentz-ether trap — letting "characterized departure" quietly mean *every* result is a success. Predeclare the questions and the landings, with numeric cut-lines wherever possible. [src: `docs/em_phaseC_force_decomposition.md:277-283, 319-320`]

**Predeclared questions.** [src: `docs/em_phaseC_force_decomposition.md:284-295`]
- Does a generic charge-odd `O(V)` source exist after the mass/drain channel is ablated?
- Is it `∝ s·V`, or does it carry additional tensor structure?
- Does the *same native force package* (same body fields, boundary operator `𝔅`, flux law, typed native roots) give the static sign **and** the moving sign?
- Does the transverse term have the Maxwell/Darwin tensor form **and relative coefficient** — not merely `1/R²`?
- How much residue lies in the propagating `h` / longitudinal-scalar modes?
- Are `c_E`, `c_γ`, and the coupling coefficients **tied**, or independently fitted? (the surplus-vs-knobs ledger)
- Do preferred-frame terms remain? (the Lorentz-violation falsifier)
- Is there any approximate additive conserved charge in the multi-throat theory?
- Which results rely on the postulated ambient `R_w` symmetry vs the one-sided slab actually computed?

**Predeclared landings.** [src: `docs/em_phaseC_force_decomposition.md:297-302`]
- **Maxwell-like low-energy EFT** — Darwin relative coefficient within a predeclared tolerance.
- **Scalar–vector characterized departure** — right sign/structure, wrong relative coefficient beyond tolerance.
- **Signed two-body force analog, but NOT EM** — orientation sign exists, but tensor structure / propagation fail.
- **The named KILL:** **a null or wrong-sign result means the throat-current mechanism FAILS.**

**On the KILL (first-class, welcome).** [COMMITTED methodology] If a branch yields no charge-coupled orientation channel, or a charge-coupled sign opposite to what the defect/antidefect wall-healing picture requires, that is a **genuine result that kills the mechanism on that branch** — it is not a defect of the exercise and must **never** be rescued, softened, or re-scoped away. A clean "it all works" across every branch would itself be suspicious. The whole protocol is built to be *able to fail per finding* (§0.4). [src: `docs/em_phaseC_force_decomposition.md:279-282,302`]

**Interpretation guards (do not trip these).** [src: `docs/em_phaseC_force_decomposition.md:326-337`]
- Do not call a `V⁰` orientation effect magnetic — it is the static electric candidate.
- Do not assume same-side flow ⇒ attraction / opposite-side ⇒ repulsion — it may be attraction-vs-zero or a return-closure artifact.
- Do not import viscous jet entrainment into the ideal superfluid — the native stress/action decides the sign.
- Do not infer `s·V` from geometry alone — a moving orientation can yield charge-even, higher-order, or no transverse source.
- Do not let the gravity drain masquerade as charge — keep the mass-only ablation, while still classifying genuinely `w`-odd drain-derived cross terms.
- Do not ignore the ambient seam — the executable background is one-sided; clean `s1·s2` parity uses the postulated two-sided background. Run both or give an explicit branch map.
- Do not treat one body's coupling to the bulk as a two-body force — it must descend from a genuine cross interaction or shared constraint/flux channel.
- Do not tune the post-mouth contribution independently after seeing the result — shared profiles/coefficients stay shared, or the hierarchy test loses meaning.
- Do not use a large throat-creation barrier as a force — only its separation-*dependent* pair deformation attracts; the isolated core is mass.
- Do not claim the hierarchy from a suggestive large ratio — power, normalization, and parameter relation must come from the same native throat solution.

**The output sentence.** After running §5.1 → §5.4 on a solved branch, emit exactly: *"On `𝔅`-branch B, the mouth ensemble is [source/defect], the static orientation coefficient is `A_q` [sign], giving [attraction/repulsion/none] between like charges; the moving orientation kernel `K_orientation` is [∝ s·V / characterized departure / null] with magnetic sign [X]; gates [topology/range/stability/momentum-closure] [pass/fail]; landing = [Maxwell-like / characterized departure / signed-but-not-EM / KILL]."* That sentence, one per branch, is the deliverable of the whole two-throat exercise.

---

## §6. Numerical-feasibility appendix — why it is hard, and what is reachable

The full object — a time-dependent, fully nonlinear, 4+1-dimensional solve of two moving throats radiating into an open medium — is **presently blocked first by the missing parent-action closure (§2.6) and, even after closure, by severe solver tractability** [COMMITTED scope decision]. This appendix states the numerical walls so a solver can see which reductions buy tractability and grades each ladder rung (§7) conditionally. The upstream program deliberately worked at the collective-coordinate / effective-action level for exactly these reasons [`docs/em_u1_body_definition.md:138`; the "in place vs missing" list at `docs/em_phaseC_force_decomposition.md:204-231` enumerates what the analytic pathA work does *not* compute: the cross-interaction, the nonlinear mouth boundary, and the derived sources].

### 6.1 The five tractability walls
1. **Dimensionality.** The medium lives in 4+1D `(x, y, z, w, t)`. A naive grid is 5-dimensional. Every reduction below is aimed at collapsing this.
2. **Scale separation (stiff geometry).** The throat carries at least four disparate lengths: mouth radius `a`, axial extent `L = 1.85 a`, healing/wall width `ℓ = ε_r a = a/20`, and — for two throats — the separation `R`. Resolving the `a/20` wall while reaching the far field over many `a` is a multi-scale meshing problem; the thin healing layer is where the boundary operator `𝔅` acts, so it cannot be coarsened away.
3. **Stiffness in time (disparate speeds).** The medium carries a bulk compression speed `c_s` (gravity-change propagation), the brane shear speed `c_γ` (light), and slow throat motion `V`. These differ by orders of magnitude, so an explicit time march is CFL-limited by the fastest wave while the physics of interest evolves on the slowest scale — classic stiffness. (This wall is absent from the *static* rungs, which is why the ladder front-loads them.)
4. **Open boundaries.** The throat is a one-way sink (the drain); the domain is unbounded and, under motion, radiates. Naive reflecting boundaries corrupt the answer. Absorbing treatments are needed — sponge/Rayleigh layers, complex absorbing potentials (CAP, `−iσ(x)`), matter-wave PML, or exterior complex scaling — but only where waves actually propagate. [`docs/boundary_and_noise_methods.md:70-107` catalogs the families.]
5. **The free/nonlinear mouth boundary.** The sleeve surface `Σ` where `𝔅` is imposed is itself part of the solution (a free boundary), and for the mixture/constraint branches `𝔅` is a velocity-level nonholonomic condition, not a simple Dirichlet/Neumann datum. This is the hardest structural feature.

### 6.2 The reductions that buy tractability
- **Static ⇒ no time-marching stiffness (wall 3 gone).** For the *static* one- and two-throat problems the disparate wave speeds do not create a CFL limit. Once a stable closure card supplies the missing action/operators, the intended reduction is an elliptic BVP; ellipticity must be checked for that instantiated operator. **On a stationary branch the open-boundary treatment also simplifies: there are no propagating waves to absorb, so a declared far-field/return matching scheme replaces CAP/PML, which is deferred to radiative rungs** [`docs/boundary_and_noise_methods.md:19-28` §0 decision of record].
- **Axisymmetry (wall 1 partially gone).** A single static throat is axisymmetric about its `w`-axis, collapsing `(x, y, z, w)` to `(r, w)` — an effectively **2D** elliptic BVP. This is the regime of rung R1.
- **Perturbative motion (wall 3 controlled).** A slowly-moving single throat is the static solution plus an `O(V)` linear response `p(V)` — a 2D static solve plus a linear correction, no full time march (rung R2).
- **Reduced two-body geometry.** Two static throats break axisymmetry but retain a symmetry plane; the problem is 3D but with exploitable structure (bispherical-type coordinates, or a fixed two-center mesh). This is the R3 regime — *plausibly* reachable, the first genuinely 3D rung.

### 6.3 Feasibility grade per rung (detail in §7)
| Rung | Problem | Reduced dimensionality | Grade |
|---|---|---|---|
| **R1** | static one-throat sleeve profiles, per branch | intended 2D `(r, w)` elliptic BVP | **likely after the §2.6 closure card is supplied** |
| **R2** | moving one-throat `p(V)` response | 2D static + linear `O(V)` correction | **likely after R1 closure** |
| **R3** | static two-throat → electric sign | intended 3D reduced (two-center) elliptic | **plausibly reachable after R1 closure** |
| **R4** | moving two-throat → magnetism sign + current test | 3D/4D + radiation (stiffness + open BC return) | **hard summit; may require HPC, may stay deferred** |

### 6.4 The honest guardrail
This ladder is designed for **partial ascent that still yields real physics once the relevant closure is supplied**. Even if R4 remains out of reach, R1–R3 can produce concrete, falsifiable results: per-branch sleeve profiles/admissibility, the tilt response, and the static electric `±w` force sign per surviving branch. The **full** nonlinear two-throat sim stays deferred; any question that genuinely requires it (e.g. the ~10⁴² magnitude, the full current law) stays [OPEN], and a null or wrong-sign result genuinely kills the mechanism it tests (§5.4). The ladder organizes the partially-closed program; it does not complete the missing theory.

---

## §7. The simulation build-plan (the ladder)

This section is the build plan for the AI-driven (or human-driven) simulation. It is a **ladder of four conditional rungs**: each becomes an independently-checkable numerical problem only after the required §2.6 closure data are pinned. The upstream analytic program was governed by an anti-self-deception discipline; §7.0 transfers that discipline to numerics, and **every rung must pass §7.0 before its result is banked**. The design goal is honest partial ascent: no rung is allowed to manufacture closure, unique `𝔅` selection, or a current law.

### 7.0 The numeric gauntlet (applies to EVERY rung — a banked result must clear all seven)
1. **Closure-card gate.** Before discretization, archive the exact action/PDE coefficient card required by §2.6: the `χ` kinetic coefficient; wall stiffness/potential; sleeve bending/anchoring/free-boundary functional; `χ↔h` and `χ↔drain`/geon/source couplings; `M_h,K_h,c_E,C_hu`; coupled EOS/wave operators and interface/jump laws; the branch-specific E1–E5 mathematical completion; and the IR/return/zero-mode/outgoing scheme. Identify every entry as [COMMITTED], [CALIBRATED], [POSTULATE], or [OPEN] supplied for this run. Missing entries stop the rung at `UNRESOLVED(closure)`; supplied OPEN entries make the result explicitly conditional and may not be fitted to the desired sign.
2. **Manufactured-solution (MMS) verification.** Before trusting any physical solve, verify the instantiated discretization on a **manufactured exact solution**: insert a known analytic field into the operator, carry the residual as a source term, solve, and confirm the numerical solution recovers the manufactured field at the scheme's **design order of accuracy**. An MMS "tooth" per operator (interior residual, each completed boundary-operator class `𝔅_b`, interface/jump/free-boundary residual, and zero-mode quotient) must reproduce its design order, or the discretization is wrong.
3. **Convergence gates (the ≥7-point rule).** No convergence claim — for a profile, a force, or a sign — may rest on fewer than seven grid resolutions. Fit the observable as `χ(n) = χ_∞ + c·n^(−p)` over **≥7 resolutions**, require the extrapolant `χ_∞` to be **jackknife-stable** (dropping any one resolution does not move it beyond its error bar), and require **independent convergence on each refinement axis** (spatial `h`, and for R2+ the perturbation amplitude and, for R4, the timestep) — a quantity that converges in `h` but drifts in the timestep is not converged. Only the gate-clean Richardson extrapolant is banked.
4. **Dimensionless, target-blind verdicts.** The banked output of every rung is a **held-out dimensionless quantity or a sign**, computed without tuning to the expected answer: a sign (`+`/`−`), a power-law exponent, a ratio (e.g. `F_e/F_g`), or a shape. Dimensionful absolutes (`G`, the Coulomb magnitude) are calibration and are **not** verdicts. The solver must be able to return either force sign.
5. **Provenance / ancestry guards (no smuggling).** Inputs are the committed native structure plus the explicit closure card: **no Maxwell equation, no `j = sV`, no imported like-repel prior, no point source, no `θ = s·p` convention**. In particular, the conditional pathA39 kernel (§2.4d) is barred from R2/R4 ancestry until R2 natively derives a source `∝sV`; it may then be used only as a comparator. `L/a` and `ε_r` are frozen §4.1 inputs, not outputs to “reproduce,” and no geometry or closure number may be adjusted after seeing a sign.
6. **Fresh-agent adversarial recompute.** Every banked verdict (each sign, each profile family, each ratio) must be **independently recomputed by a fresh solver instance / agent** that re-derives the decisive quantity from the discretization, and adversarially tries to break it (perturb the mesh, boundary implementation, zero-mode fix, and closure sensitivity). A verdict that only one implementation produces is not banked.
7. **Per-branch, non-evidential accounting.** Attempt all 16 branch templates (§3.4) for which a closure-complete operator is available; report results **per branch**. A branch with no finite-energy/stable solution may be explicitly excluded. A branch that lacks closure, fails convergence, or fails MMS is reported `UNRESOLVED`, not silently dropped. Survival does not select `𝔅`, and uniqueness is never inferred from conditional solves.

### 7.1 Rung R1 — the 16 static one-throat sleeve attempts (profiles + per-branch admissibility)
**Problem.** After the §7.0 closure-card gate, instantiate and solve each available conditional template in §3.5 over `(r,w)` (axisymmetric, 2D): the assembled static interior residual, completed `𝔅_b` on `Σ`, ambient matching, interface/free-boundary laws, and exact zero-mode quotient. **Output:** the eight sleeve profiles and the closure/residual/energy/stability record per branch.
**What R1 can decide.** R1 can legitimately **exclude** a branch with no nondegenerate finite-energy/stable solution and thereby narrow the candidate set. It reports every remaining solution as *admissible conditional on its imposed `𝔅_b`, ambient, and closure card*. Existence does not prove nature enforces that operator. Unique selection of `𝔅` is not guaranteed by this one-body static solve and remains [OPEN] if several branches survive; a later energetic/dynamical selection calculation would require a separately supplied criterion and closed action.
**Acceptance (beyond §7.0).** `L/a=37/20` and `ε_r=1/20` are declared frozen inputs, not validation targets. Require small interior/boundary/interface/free-surface residuals; finite energy and the stated scalar stability bound; ≥7-point convergence for all eight profiles and energy/stability observables; controlled domain/IR extrapolation; sensitivity maps for every supplied OPEN coefficient; and a fresh-agent recomputation. `c0` and `χ_Q` are absent from static R1 acceptance (§4.1).
**Grade:** likely achievable only after the missing closure is supplied (§6.3).

### 7.2 Rung R2 — the moving one-throat response `p(V)` (discharges the current-law question)
**Problem.** For each branch with an admissible R1 solution, compute the linear `O(V)` response of the sleeve to slow motion: the static R1 solution plus a linear-response solve for the **tilt** `p(V)` (the in-brane projection of the sleeve tangent, §2.5/§5.1) and the induced mediator couplings. No full time march is required for the quasistatic response (§6.2).
**What it decides.** (a) Is `p(V)` **linear in `V`** (the analyticity/isotropy expectation) or a departure? (b) Channel by channel, does the `O(V)` mediator coupling come out **`∝ sV` (a genuine charge current — the imported `j = sV` becomes *computed*)**, a **characterized departure**, or **zero (this magnetism source dies natively)**? All three are first-class outcomes (§5). This rung is where "does a moving charge carry a current?" is answered from the medium, replacing the imported current law.
**Acceptance:** §7.0; convergence in perturbation amplitude (the `O(V)` coefficient must be amplitude-independent as `V→0`); a channel-by-channel parity/source census (§5.1); and the zero-velocity ablation. The quadrupole coefficient `c0` belongs here or in R4 if that module is used: if prescribed from §4.1 it is an ansatz input, not a reproduced result. The pathA39 `j∝sV` kernel remains ancestry-barred unless this rung independently lands on that source. **Grade:** likely achievable after R1.

### 7.3 Rung R3 — the static two-throat force (⭐ the electric `±w` sign)
**Problem.** For each branch with an admissible R1 solution, place two static throats at separation `R` and solve **all four orientation cells** `{++,+−,−+,−−}` in the reduced two-center geometry (§6.2). **Because the four-orientation census needs both `±w` throats, R3's electric-sign comparison runs on the two-sided `R_w`-ambient branches only (§3.4); the one-sided branches yield the `++` like-charge magnitude but not the sign, so the electric sign is `R_w`-postulate-conditional.** For every cell compute the force separately in all four channels `ch∈{var,flux,𝔅,rad}` using the fixed closure-card prescription, then apply the four-orientation Hadamard decomposition of §5.2. A single momentum-flux surface is not the total force. **Output:** per branch and per separation, the channel-resolved `(00),(10),(01),(11)` force components, the total electric `(11)` sign, and its radial falloff.
**Why this is the headline.** This is the sign the committed model could not determine. R3 reports it **per surviving conditional branch**, not on an R1-selected physical branch, and adjudicates the mechanism against §5.4's rubric and named KILL.
**Acceptance:** §7.0; all four orientation cells and all four force channels converged separately; `(10)` and `(01)` vanish within tolerance on the symmetric ambient or are reported as asymmetry contamination; the conservative `var` channel passes the §5.2 integrability check before any `U_11` is quoted; total body-plus-medium momentum closes through named flux/constraint/radiative terms; orientation/drain/mediator ablations in §5.3 fire; and the fixed-source ↔ fixed-defect exchange is run as a diagnostic control. The exchange must show the appropriate change in coupling character/stiffness dependence; **no force-sign flip is prescribed**. **Grade:** plausibly reachable after closure (first genuinely 3D rung).

### 7.4 Rung R4 — the moving two-throat force (the magnetic sign + the current confirmation)
**Problem.** For each branch with an R1/R2 solution, run **all four orientation cells** `{++,+−,−+,−−}` for each predeclared velocity geometry (parallel, antiparallel, and the one-body-zero controls). **Like R3, the four-orientation census requires both `±w` throats, so R4's magnetic-sign comparison runs on the two-sided `R_w`-ambient branches only (§3.4); the magnetic sign is `R_w`-postulate-conditional.** For every cell compute `F_var`, `F_flux`, `F_𝔅`, and `F_rad` separately, apply the Hadamard decomposition at fixed `(V₁,V₂)`, and extract the bilinear orientation kernel of §5.2(c). This yields the magnetic sign per branch and tests the R2 source in a two-body setting. Radiation reintroduces stiffness and open-boundary absorption (§6.1 walls 3–4).
**Acceptance:** §7.0; all four orientation cells and channels converged; the magnetic `(11)` term vanishes when either velocity or the orientation structure is removed; the `var` piece passes integrability before potential language is used; total momentum closes including `F_rad`; timestep and absorbing-boundary/domain convergence pass; and all §5.3 gates/ablations are reported. The outgoing/retarded `χ_Q` belongs only to this radiative rung: if fixed to `1` it is an ansatz input, while a fuller outgoing solve may compute and replace it. The conditional pathA39 kernel may be compared only if R2 first derives `j∝sV`. **Grade:** the hard summit — may require HPC and may remain [OPEN]; a channel-complete quasi-static result is still real physics.

### 7.5 Ladder dependencies and the honest exit
```
R1 (static one-throat, 16 branches) ──narrows admissible set──▶ R2 (moving one-throat p(V), per branch)
      │                                                        │
      └────────────────────────┐                              │
                               ▼                              ▼
                    R3 (static two-throat: ELECTRIC SIGN per branch)   [current law j∝V?]
                               │
                               ▼
                    R4 (moving two-throat: MAGNETIC SIGN + current test per branch)
```
Each closure-complete rung is **independently bankable**: R1 can discharge the R49 profiles and exclude inadmissible branches; R3 can deliver the electric sign per surviving branch. It does not follow that R1 uniquely selects `𝔅`. Questions a rung cannot reach stay [OPEN], and a mechanism that a rung kills stays killed. The ladder's contract is to expose closure debt and gauntlet every instantiated numerical problem, not to call the unfinished parent action complete.

---

## §8. Provenance index (reference only)
The physics above is self-contained; this index lets a reader who wants the audit trail find it. Consolidated equation sources: `docs/conceptual_foundation.md` (§1–§4), `docs/em_u1_body_definition.md` v3.1, `software/stage1_solver/reports/pathA_{29,36,38,39}*`. The `𝔅` endpoint table: `em_u1_body_definition.md` §2. The two-throat decomposition + rubric: `docs/em_phaseC_force_decomposition.md` §3/§4/§6/§7. Calibration: `research/pde_ledger/redteam_adversarial/{ANSATZ_LEDGER.md,benchmarks.yaml}`. The tilt-profile leaves (R49): `research/pde_ledger_v2/notes/parameter_register.md`. The 16-branch conditional templates this BC packet is drawn from: the U2 v12 amendment (directive `directive_u2_boundary_adjudication.md` §4; wrap commit `6a6f317c`); the `𝔅`-undetermined verdict: the U2 adjudication (wrap `5ceebb24`). Open-boundary methods: `docs/boundary_and_noise_methods.md`.
