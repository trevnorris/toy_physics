# Conceptual Foundation — what each part of the model physically IS

**Read this FIRST.** This is the canonical plain-language statement of the model's physical vision — what the medium, the
brane, the four force-sectors, and the particle/defect actually *are*, in the model's **own native terms**. It exists so the
conceptual picture is never lost between work sessions and never has to be re-explained from scratch. It is a **living document**:
update it whenever the conceptual picture sharpens. It is deliberately separate from the math/derivation/verdict machinery
(those live in `software/stage1_solver/decisions/*`, the `pathA_*` directives, and `STATUS.md`).

**Status:** v3, 2026-06-23. (v2 added the little-arrows brane hypothesis.) **v3 records a pivot:** the little-arrows *domain-wall*
brane was **FALSIFIED** at pathA_24 T1 (a static polar-vector wall on the connected `S³` vacuum manifold spreads/unwinds —
tri-reviewed genuine, commit `2fa91886`). A prior-art survey + two user reframes led to the **GNLS polar-smectic superfluid**
candidate (the arrows now ride on a **smectic density-layering** brane, not a self-localizing wall) and a **methodology shift**:
we are building an **analog and testing whether one self-consistent structure satisfies all requirements** — NOT deriving the
universe. Full requirements list, prior-art survey, candidate structure, and the consistency-gate plan now live in
`docs/medium_requirements_and_prior_art.md` (read alongside this doc). Changes below: §0 reframe rule; §2 brane refinement (density
now modulates); §8 new plan.

---

## 0. The non-negotiable conceptual rules (how to think about this model)

These are load-bearing. Violating them is the recurring failure mode (importing standard-physics habits that don't fit).

1. **ONE medium, everything emergent.** There is a single substance — a compressible superfluid filling 4 spatial + 1 time
   dimensions. Gravity, magnetism, electricity, light, mass, charge, and even **cosmic expansion / dark energy** (§5) are all
   *behaviours of that one medium*. Nothing is layered on top as a separate fundamental field. When the math starts looking like
   "the medium PLUS a separate gauge field on its own metric," that is **drift**, and it is to be pruned, not accommodated.
2. **There are NO point particles.** A particle is an **extended** object — a finite-size throat/puncture in the brane. Nothing in
   this model is a zero-size point. The classical point-particle pathologies (infinite self-energy, etc.) are things this model
   *removes by construction*, so any argument that imports "point defect / point charge" machinery is using the wrong tool.
3. **Think in the medium's native mechanisms, not imported abstractions.** The real handles are: *flow*, *swirl*, *puncture
   direction*, *brane shear*, *surface tension*, *cohesion*. Before reaching for "topological winding number," "π₂ of a vacuum
   manifold," "gauge field," or "point particle," check whether a native mechanism already covers it. (In particular: **electric
   charge is a puncture *direction*, NOT a winding** — see §3. Conflating the two is a standard-physics reflex to avoid.)
4. **Falsification is the goal.** A result that BREAKS the concept is welcome and first-class. Never rescue or soften a failure to
   protect the picture. A clean "it all just works" is the *suspicious* outcome.
5. **Framing split.** In private/working context (like this doc) we engage the unification picture directly. In public
   papers, the framing is strict toy-analog only. This doc is the private conceptual vision.
6. **Analog, not derivation — find a consistent structure (added v3, load-bearing).** We are NOT trying to derive the medium /
   brane / constituents from first principles (they could be arbitrarily deep — e.g. 4D structures in their own 5D space;
   unknowable). We are building a mathematical analog and asking the finite question: **is there a single self-consistent
   superfluid structure that satisfies ALL our requirements at once?** We may **postulate the structure freely**; the test is
   internal **consistency**. This makes falsification *sharper*: the failure mode is a genuine **no-go — two requirements
   mutually incompatible in any structure** (the T1 three-way no-win was the first). The requirements list + the candidate
   structure live in `docs/medium_requirements_and_prior_art.md`.

---

## 1. The medium — the one substance

- A **compressible superfluid** filling **4 spatial dimensions + time** (coordinates `x, y, z, w`; `w` is the "extra"/bulk
  direction perpendicular to our space).
- It is **made of its own tinier constituents held together by their own cohesive forces.** It has *substructure* beneath the
  smooth mean-field description. This is not a decorative detail — it is the **source of the medium's material properties**:
  - **Cohesion → surface tension.** Cohesive forces between the constituents are exactly what let the medium have *surfaces/
    interfaces with tension* (the same reason water has a skin). This is the root of the brane (§2).
  - **Cohesion → transient/high-frequency elasticity.** At high frequency the medium can resist deformation elastically even
    though it flows freely at low frequency. This is the root of the brane being able to carry **shear** (light, §3).
- **Epistemic status — we *infer* the substructure from its effects; we do not observe it directly (keep this as a lens).** We
  only see the medium's behaviour at our scale and reason backward to what its constituents must be like. They could obey a far
  richer set of rules we know nothing about — we are *supposing* how they work. This is not optimal, but it offers a diagnostic
  **lens we should keep looking through**: when an interaction, force, or behaviour at our level doesn't make sense, the
  resolution may be that we haven't yet defined the superfluid's substructure richly enough. **Treat a puzzle in the emergent
  physics as a hint about the constituents, not automatically as a flaw in the picture** — and conversely, don't over-fit the
  substructure to one puzzle without checking it against the rest.
- Low-frequency, large-scale: it behaves as an (essentially inviscid) **fluid** — it flows and it can swirl, but in the bulk it
  **cannot sustain shear** (no sideways rigidity). *Keeping the bulk shear-free is load-bearing: it is what protects magnetism*
  (§3).
- Two wave speeds live in the medium and both depend on its density `ρ0`:
  - `c_s` = the **ripple speed** — the linear density/phonon wave speed of the medium. Gravitational/density ripples ride it; it
    plays the analog of "the speed of light" for the *gravitational* sector. (`c_s² = 5Kρ0⁴/m` in the current spec.)
  - `c_γ` = the **light speed** — the speed of the brane's shear wave (§3), `c_γ = √(μ_brane/ρ_brane)`.
  - Whether they lock (`λγ ≡ c_γ/c_s = 1`, as GW170817 requires of real light vs gravity) is an **open** derivation, not an
    assumption.

---

## 2. The brane — why our 3D space exists

> **⭐ v3 UPDATE (2026-06-23) — the brane mechanism changed; read this before the little-arrows prose below.**
> pathA_24 **T1 FALSIFIED** the little-arrows *domain-wall* brane: a static polar-vector order parameter has a **connected** vacuum
> manifold (`S³`, π₀=0), so a `+w`/`−w` wall **spreads to infinite width and unwinds with zero barrier** — no stable wall, no flat
> core, no confinement (tri-reviewed genuine; commit `2fa91886`). Independently corroborated by prior art (Davies–George–Volkas:
> isotropic 4-space gives no 3+1D zero modes without imposing the axis; ³He's A–B interface is externally maintained, not
> self-localizing). **The refinement:** the brane is now a **smectic (1D density-layering) structure** — a codim-1 surface that is
> the spontaneous **ground state** (not a metastable kink), with an emergent axis, immune to the spread-and-unwind failure. **The
> little arrows are KEPT** — they ride on the smectic layer and still do their two jobs (in-plane shear = light; `±w` flip across
> the layer = the two charge signs). **Honest change:** the **density now MODULATES** (the layering); the claim below that "the two
> states live in orientation, not density, so `U(ρ)` stays single-well and we never modulate density" is **superseded** — `U(ρ)`
> stays single-well *locally*, but a new non-local/roton interaction drives a finite-`k` density modulation. Full candidate
> structure (the **GNLS polar-smectic superfluid**), requirement→mechanism map, and the consistency gates:
> `docs/medium_requirements_and_prior_art.md`. The prose below is retained as the conceptual ancestor.

This is the most *imposed*, least *derived* part of the construction, and pinning it down is a current frontier.

**What it is.** Our 3D space is a **brane**: a 3-dimensional surface sitting at `w ≈ 0` inside the 4D bulk. The physical picture
(your water-surface instinct, made precise) is that the brane is a **domain wall / phase interface** in the one medium:

- The medium has (at least) **two degenerate stable states**. Wherever two of them meet, you get a codimension-1 surface — a
  wall — separating them.
- That wall carries an **energy per unit area = surface tension** (literally the water-surface intuition; the cohesion of §1 is
  where the tension comes from).
- It is **stable** because it interpolates between genuinely distinct states — it can't just relax away.
- **We are confined to it as bound zero-energy modes** sitting in the wall's potential well. Leaving the sheet (moving off into
  the bulk along `w`) costs energy, so the lowest-energy excitations — us, and everything we're made of — are stuck to the sheet.

**Bulk on BOTH sides (a real deduction, not a choice).** Electric charge has two signs, realized as a throat puncturing toward
`+w` or `−w` (§3, §4). For *both* directions to be available there must be bulk to drain into on **both sides** of the brane.
A one-sided "surface of a pond" (medium below, vacuum above) would allow throats to drain only one way → only one charge sign.
**So the brane is an interface in the *interior* of the medium, with bulk on both sides** — the ³He A–B interface flavour, not
the outer surface of a finite body.

**The obstacle (honest).** The current medium potential `U(ρ) = (K/4)ρ⁵` is a **single well** — one stable state — so as written
the medium *cannot* form a wall. That is exactly why the brane currently has to be *put in by hand* (a confinement potential
`V_conf`, the `w`-profiles `Z/W/B_ℓ`, the `k_w u_w²` restoring term — all *inputs*, not derived). To **derive** the brane, the one
medium needs structure that gives it **two coexisting states.**

**The working hypothesis (agreed 2026-06-23): the "little arrows".** The two states come from giving the medium's tiny
constituents the simplest non-trivial property beyond mass — a **polar orientation**: each constituent is a little *arrow*, with
a head and a tail (genuinely polar, **not** a headless axis). Then:

- **The arrows align** (lower energy aligned, like little magnets) along **one axis** in the 4D space — and *that axis is what we
  call `w`*. So `w` is not pre-existing geometry the arrows find; **`w` is *defined* by the alignment direction.**
- **Why our space is 3D:** align along a single axis in 4D and the surface perpendicular to it is 3D. The codimension-1-ness of
  our space falls straight out of the orientation being a *single* direction.
- **Two domains ⇒ a wall.** When the medium first settled it needn't have picked the same direction everywhere — domains of
  arrows-pointing-`+w` and arrows-pointing-`−w` form generically (the way any symmetry-breaking-on-cooling leaves domains with
  walls between them). **We live on a wall between a `+w` domain and a `−w` domain** — two genuine mirror-image states of the
  *same* medium (so the two sides, and the two charge signs, are symmetric).
- **The wall's own structure does triple duty.** To get from `+w` on one side to `−w` on the other, the arrows must *rotate
  through the wall*, and **at the wall centre they lie *flat, in the plane of the brane*.** Therefore:
  - in the **bulk**, arrows point ⊥ to our space (`±w`) → no in-plane structure → **bulk stays shear-free → magnetism preserved**;
  - on the **brane**, arrows lie *in* our space → the brane has genuine in-plane orientational structure → **it can carry in-plane
    shear → light**.
  The *same* arrow-field gives the wall's existence, the brane's shear (light), and the bulk's shear-freeness (magnetism) at once.
- **It dodges the single-well obstacle cleanly.** The two states live in the *orientation*, not the density, so `U(ρ)` can stay
  single-welled and monotonic — `ρ` ("how much medium is here") remains single-valued everywhere. We never fake a double-well.
- **Charge falls out.** A throat puncturing toward the `+w` domain vs the `−w` domain = the two charge signs — pure puncture
  *direction*, no winding (§3/§4) — and we now see *why there are exactly two*: there are exactly two arrow-domains to puncture
  into.
- **It is substructure, not drift.** Giving the constituents an orientation is acknowledging a property the *one* medium's
  particles already have (the §1 lens), not adding a second substance. Existence proof: superfluid ³He — its constituents carry
  an orientation and it is emphatically *one* superfluid with structured A/B interfaces. Firewall: the arrows must be *carried by*
  the medium's particles, not an independent decoupled field.

*This is the working hypothesis to test, not a settled result.* Make-or-break checks: the arrows must be **polar (`+w ≠ −w`)** or
there is no wall and no two charge signs; and the in-plane wall texture must actually yield **two shear polarizations at one
speed** (the photon signature) rather than the wrong count/dispersion. See §8 for the test plan.

**Alternatives considered and set aside** (kept for the record): **(a) two densities** — a liquid/vapour coexistence from the same
cohesion; clean, but "dense vs dilute" are not mirror images (asymmetric → asymmetric charge) and a pure density wall is a fluid
membrane (no in-plane shear → no light), so too weak to be fundamental. **(b) self-trapping by the collective drain network** —
"all the distant throats together = the gravitational field" as a self-consistent well; maximally native but a *persistence*
mechanism, not an *origin* (it's circular — the brane holds matter because matter deepens the brane), so likely a real
*secondary* effect that deepens the well, not the seed.

**The unification (why this question is the keystone).** Four problems are really *one* question — *what is the wall's internal
structure?*
- A **structureless** wall = a *fluid membrane*: tension + bending only, **no in-plane shear** ⇒ no light (this is the pathA_23
  "Stage-2" wall).
- A **structured** wall (carrying internal orientational order) **can** carry in-plane shear ⇒ it **derives light**, that internal
  order **is** the "substructure beneath the mean-field" the model keeps needing, and a genuine emergent wall is what makes the
  brane a legitimately separate sector from the bulk (so brane-shear can decouple from bulk-flow without it being an ad-hoc
  trick). One answer to "what is the wall's internal structure" settles brane-existence, light, the substructure question, and
  the no-leak/separate-sector question together.

---

## 3. The four sectors — keep them distinct

This is the heart of the model. The four phenomena are **four different mechanisms of the one medium.** Do not conflate them
(especially: electric ≠ magnetic; charge is a direction, not a swirl).

> **⭐ v4 UPDATE (2026-06-25) — model-mechanics corrections that keep getting LOST between sessions (read this; it is the recurring
> relapse into textbook defaults). Mirror of memory `[[project-model-mechanics-corrections]]`.**
>
> 1. **NOTHING in this model is static.** Every static-field simplification has FAILED because it deletes *the coupling that
>    selects the balance.* The acoustic horizon, the throat size, and the drain are all self-consistent **DYNAMICAL** balances —
>    not fixed surfaces or snapshots (throat self-selection "requires dynamics"; the horizon is just where rising `v_r` meets
>    falling `c_s`, both downstream of the *same* drain). If you catch yourself freezing a field, STOP.
> 2. **THREE distinct speeds — never conflate them:**
>    - `c_s` (**ripple**) = the speed gravitational **CHANGES** propagate ("the speed of gravity"). `c_s²=5Kρ⁴/m` ⇒ **`c_s ∝ ρ²`**
>      — density-dependent, *slower in low density.* **That density-dependence IS the lensing/Shapiro mechanism** (a mass depletes
>      ρ → lowers `c_s` → bends/delays signals).
>    - `v_r` (**inflow / drain**) = the **STRENGTH** of the (quasi-)static field (like the *value* of `g`), **NOT a propagation
>      speed.** `v_r ≪ c_s` is normal = **weak** gravity, NOT "*slow*" gravity.
>    - `c_γ` (**brane shear**) = **light** (`c_γ=√(μ_br/ρ_br)`), a separate speed. `λγ=c_γ/c_s=1` (the GW170817 cone-lock) is
>      still **open** (§7 #8).
> 3. **Gravity = the FLOW between draining defects** — carried by the flow + Bernoulli pressure, **NOT** by ripples/radiation.
>    But a *change* in it propagates at `c_s` (the flow can't re-aim at a distance without a pressure signal). A *uniformly
>    moving* drain drags its inflow with it → the field points to the drain's **CURRENT position** → **no aberration, no
>    Laplace/solar-system instability.** This is **already PROVEN** in `research/1pn_orbital_dynamics` (the "**Static Limit
>    Theorem**": elliptic/instantaneous near-zone; the 1PN precession comes from inertial dressing `m_eff(r)=m[1+σ(r)]`, not from
>    a potential lag) — that old tension is RESOLVED.
> 4. **Throat-soliton — NO sloshing.** The trapped light = **mass** is *statically bound* (mode larger than the neck → sub-cutoff
>    → evanescent) and carries **NO net trans-brane current** ⇒ `J_w=0` on the exact mode is **EXPECTED and FINE.** The photon's
>    only job is to hold the throat **OPEN.** The old **"AC→DC rectification" idea is RETIRED.** Gravity is a **SEPARATE** steady
>    *background* de-structuring drain through the open throat (ordered brane → unstructured bulk ground state, escaping through
>    the puncture). **Photon ≠ gravity.**
>
> **And on scope (so it never gets re-litigated): the conservative PN gravity ladder (1PN→4PN + the 2.5PN radiation-reaction term)
> is already BUILT & GR-matched (calibrated) in `research/4d_*pn*`** — don't re-derive it. The end goal is a fully **calibrated**
> PDE (`research/pde_ledger/`) delivering GR + EM. **Calibration is fine; first-principles is not required** (§0.6 — analog, not
> derivation).

| Sector | What it physically IS | Native quantity |
|---|---|---|
| **Gravity** | The medium's **inflow / drain** toward defects. Test bodies are carried inward by the flow. Largely unobservable as a flow. | flow velocity `v_r` |
| **Magnetism** | The medium's **swirl** — vorticity/rotation. Magnus force. Lives in the bulk and needs **no shear** (which is why a shear-free bulk is required). | vorticity / circulation (the *winding*) |
| **Electric charge** | The **puncture direction** — which way a throat punctures the brane into the bulk (`+w` vs `−w`). A **binary orientation**. No rotation, no swirl, no winding. | puncture orientation `±w` |
| **Light** | The brane's in-plane **rotational-elastic (MacCullagh) shear wave** — our 3D space resisting being *locally twisted*, the wave riding on the curl of an in-plane displacement (not the drumhead's up/down bowing, which is a separate scalar mode). Two transverse polarizations, no longitudinal mode; rigidity lives **on the brane, not the bulk**. | brane rotational-elastic shear |

**Why two polarizations and no third (the MacCullagh point):** if the brane's stiffness is **rotational/curl-type** — it costs
energy to *locally rotate* the medium, energy `∝ (∇×u)²` — then the wave equation has **exactly two transverse polarizations and
no longitudinal mode**: precisely electromagnetism. (Ordinary "Cauchy" elasticity, which resists *stretching/sliding strain*
`∝ (∂u)²`, would instead give a spurious *third*, longitudinal "photon" that real light does not have — so light needs the
*rotational* kind of stiffness specifically, not generic shear.) This is MacCullagh's 1839 rotational ether, the historical
ancestor of the EM field. It is the open "structured wall" question of §2 — see the deep-dive just below.

**The summary image (keep this):** *Our 3D space is a taut elastic sheet in the 4D bulk; gravity is how that sheet drapes and how
the bulk drains through it; light is the sheet shivering sideways within itself; and the bulk underneath stays a frictionless
fluid so magnetism never feels the stiffness.*

**What light actually is — and the two ways a wiggle can carry energy (read this; it's the part that's easy to get lost on).**
We were never going to "derive photons" out of the bare scalar superfluid — a scalar has only compression waves (the `c_s` ripple,
the gravity sector). We *postulate a mode* (the arrows) and ask whether it behaves like light. Picture the arrows lying flat in the
brane and wiggle them. There are **two physically different ways the wiggle can store energy, and only one is light:**
- **Frank / orientation wave (NOT light).** The arrows re-point while the medium stays put — like a stadium wave: the people
  (material) don't move, only their "pointing" propagates. The energy cost is for *neighbours to disagree* — an orientation
  gradient, `∝ (∇n)²` — and the restoring force is a **torque**. It can even have two transverse components and a fixed speed, so
  it can *look* like light on a dispersion plot — **but it is not**, because nothing material is displaced and it exerts no
  mechanical force across a surface.
- **MacCullagh rotational-elastic wave (THIS is light).** Here the medium resists being *locally rotated at all*; energy
  `∝ (∇×u)²`, the curl of a genuine displacement. This gives exactly two transverse polarizations, no longitudinal mode, Maxwell —
  *real* electromagnetism. It is special because ordinary matter resists squeezing and sliding but **not** rotating; to resist
  rotation you need substructure made of little orientable elements that have a direction to defend — **which is exactly what the
  arrows are.** The arrows are the rotation-resisting (gyrostat / Cosserat) elements MacCullagh needed in 1839 but couldn't name.
- **The decisive difference — torque vs traction.** If wiggling an arrow only pushes back as a *torque* on its orientation → Frank
  → not light. If it exerts a real *traction* (a mechanical force across a surface, because the arrow's rotation is rigidly tied to
  the medium's own motion) → MacCullagh → light. **This is why "the arrows are carried by the one medium, not a free decoupled
  field" is mechanically load-bearing, not just philosophy:** only if tilting an arrow *is* a rotation of the material does the
  orientation stiffness become rotational-elastic stiffness — Frank turns into MacCullagh. The honest hard part (the genuine
  coin-flip, and the same wall pathA_23 Stage-2 hit): getting it **curl-only** (no spurious longitudinal third "photon," which a
  generic *Cauchy* elastic wall would have) and closing the gauge structure. (Tested directly in `pathA_24` T2.)

**Why light comes in localized packets — an honest open frontier (no clean answer yet).** This is a *separate* and deeper question,
and worth keeping distinct from the above. Three things are easy to blur: (1) does the medium carry the light *mode* at all; (2)
does that mode have EM's *structure* (2 polarizations, no longitudinal, Maxwell); (3) why does light come in **localized, quantized
packets**. `pathA_24`/T2 only addresses (1)+(2). The packet question (3) — which is really *two* hard things, localization *and*
quantization (`E=ℏω`) — is **not solved in the model**, and there is no clean off-the-shelf analog. The honest landscape:
- A *linear* wave does **not** stay a packet — it disperses and spreads. So pure MacCullagh elasticity alone gives spread-out
  waves, not lumps; the packet has to come from somewhere else.
- **The closest real analog is the optical soliton:** in a *nonlinear* medium (an optical fibre) a pulse can hold itself together
  forever — the nonlinearity self-focuses it exactly enough to cancel spreading. Our medium *is* nonlinear (the GNLS potential), so
  a brane-shear packet self-trapping into a soliton is the most natural mechanical picture of "light that sticks together."
- **The honest tension:** real vacuum photons don't interact (light passes through light), whereas solitons rely on a nonlinearity
  that usually *does* make pulses interact. So either the nonlinearity is negligible for ordinary light (and "packets" are then a
  purely quantum field-quantization effect the *classical* model can't produce), or light is weakly soliton-like and
  self-interacting — which would be a **prediction/difference** from textbook EM (and notably, real vacuum *is* weakly nonlinear at
  extreme intensity via QED photon-photon scattering). A falsifiable fork, not a fudge.
- **An illuminating contrast (helps the intuition):** a *massive particle* is a **trapped** wave — the geon standing-wave caught in
  a throat (§4); localized because it's *caught*. A *photon* is a **free** wave — nothing traps it. So a photon's packet-ness can't
  be trapping-in-a-throat; it must be either quantum quantization or soliton self-trapping. Mass and "light-in-a-lump" are
  different mechanisms even though both are "wave energy." (This frontier is deliberately **out of `pathA_24`'s scope** — a later
  "light quantization / soliton" question; forcing it into the wall test would muddy a clean result.)

**On charge being a direction, not a winding (the recurring correction):** in textbook physics "topological charge" *is* a
winding number, so it is tempting to describe electric charge as a winding. **In this model it is not.** Winding/swirl is the
*magnetic* mechanism. Electric charge is purely *which way the throat pokes through the brane*. They are independent mechanisms.

---

## 4. The defect / particle — the throat

**A particle is a throat: a puncture through the brane into the bulk.** It is **extended**, with a finite size and internal
structure — never a point. (The throat radius can't reach 0; if it did, the throat would be *closed* and the medium couldn't flow
through it into the bulk.) The decomposition:

- **CHARGE = the puncture direction (`±w`).** That it punctures *at all*, and which way. This makes charge:
  - **quantized / universal** — every puncture is one puncture, so the magnitude is the same for every particle regardless of
    mass ("the fact that it punctures at all");
  - **two-signed** — `+w` vs `−w` (matter vs antimatter);
  - **mass-independent** — the direction doesn't care how big the throat is.
- **MASS = a trapped geon-like standing wave** — a standing pattern of (shear) waves held in the throat region between 3D and 4D.
  Trapped wave → quantized wavelength → rest energy → the particle's mass.
- **THROAT SIZE / structure = a balance** — the brane's surface tension trying to **close** the puncture vs. the standing wave /
  drain flow **holding it open**. This balance sets a **finite throat radius**, which is what gives a **finite self-energy**
  (the model's resolution of the classical point-charge divergence).
- **Charge ⊥ mass.** The charge sector (puncture direction) and the mass sector (the trapped wave) do not directly interact at
  this level — which is *why* the charge is the same regardless of mass.

**Different particles are different defect structures.** Working picture:
- **Electron = a simple drain** (a clean puncture).
- **Proton = a knot-like structure** (a more complex defect).

**Open hypothesis (new, 2026-06-23 — worth chasing):** the *structure* of a defect may, through its **interaction with the
bulk**, **prefer** one puncture direction over the other — so a simple drain tends to one `±w` (one charge sign) and a knot tends
to the other. That would tie defect *structure* to charge *sign* and could be a native handle on why the common matter particles
carry the signs they do.
- *Honest refinement:* antiparticles exist (a positron = the electron's structure punctured the *other* way; an antiproton = the
  knot punctured the other way). So structure cannot *rigidly lock* direction — the "preference" must be energetic/statistical
  (one direction lower-energy or more stable), with the opposite orientation being the antiparticle. The observed cosmic
  *abundance* asymmetry (why matter dominates) would then be a *separate* question — plausibly about a small asymmetry in the
  bulk or in the two brane-states — layered on top of the structure↔direction preference. Both halves are open and interesting.

**Sharpening (2026-06-24 — light's dimensionality, how the throat traps it, and the honest death). Working synthesis:**
- **Light is INTRINSICALLY a (3+1)D brane field** — so it has **2 transverse polarizations automatically**, and it **never leaves the
  brane**. *Why it can't escape:* the bulk is **shear-free** and light **is** shear → there is no medium for light off the brane. Its
  "4D-ness" is the brane's **extrinsic curvature/embedding** (the surface bends into `w` at a throat), **not** a bulk excursion or a 3rd
  polarization. (Like a 2D wave on a sheet of paper rolled into a 3D tube: intrinsically 2D, extrinsically 3D.) This resolves the old
  worry "how can a 3D object live in a 4D throat" — it's intrinsic-vs-extrinsic.
- **How the throat traps it (the localization mechanism, answering §7 #9):** the trapped mode expands into the 4D throat and becomes
  **larger than the throat's neck → below the neck's waveguide cutoff → evanescent in the neck → genuinely BOUND** (not a slowly-leaking
  Wheeler resonance); it sheds the non-resonant energy and locks in.
- **The balance, precisely:** equilibrium throat radius `R*` set by **outward** trapped-light pressure vs **inward** (brane TENSION from
  pushing into 4D + ground-state superfluid BACKPRESSURE) → a mass–radius relation `R*(E)`; a self-bound Q-ball/soliton.
- **Causal chain (correcting a slip — gravity is the FLOW, not the curvature):** standing wave = **mass** → holds the throat **open** →
  the puncture lets superfluid **drain into the bulk** → that inward **drain IS gravity**.
- **HONEST DEATH (do not over-claim):** the trapped-wave **mass TOWER is falsified** — it predicts mass ratios 1:9:25, reality is
  1:207:3477 — and the absolute mass scale (the wave amplitude) is undetermined. So this gives **one** soliton's mass, **not** a
  predictive lepton spectrum.

---

## 5. Cosmology — the bulk↔brane cycle and dark energy

A consequence of the little-arrows wall (hypothesis, 2026-06-23). **Medium cycles between bulk and brane in two opposing flows:**

- **Matter drains medium OUT (brane → bulk).** Gravity = medium draining into defects and down their throats into the bulk. So
  matter *removes* medium from our space. This flow scales with the **matter content**.
- **Tension leaks medium IN (bulk → brane).** The bulk is drawn toward the brane (toward the matter), but to enter, its arrows
  must rotate from `±w` to flat — and that **realignment costs energy**. The wall's surface tension is exactly the **throttle**:
  it stops the bulk from collapsing in all at once and admits it as a **slow, steady leak**. This flow happens across the **whole
  brane surface**, so it scales with **area**. *(Precision: the flat wall state is not lower-energy — the wall has positive
  tension, that's the throttle; the bulk comes in because the drains pull it, and the tension limits the rate.)*

**Expansion.** Our space *is* the wall. Medium leaking in at a roughly fixed areal density forces the wall to **grow its area** →
space expands. The medium is conserved — it just cycles. The **net** of the two flows sets the cosmic dynamics (expand / hold /
contract).

**Why it could be *accelerating*, and match cosmic history** (both fall out *if the signs work*):
- The inward leak scales with **area**, so influx ∝ area → exponential, de Sitter-like **accelerating** expansion — the
  dark-energy signature.
- The outward drainage scales with **matter**. So **early (matter-dense): drainage dominates → deceleration; later (matter dilutes
  with expansion): the steady areal leak wins → acceleration kicks in.** That crossover matches the observed history
  (deceleration → recent acceleration) from the *same* one-medium picture.

**Honest caveats (the math must confirm the signs and rates):** that the net leak is *inward* not outward; that influx-at-fixed-
density really means *area growth* (not just *densification*, which would instead drift `c_s` over time — a distinct, checkable
prediction); and that the acceleration crossover timing is even roughly sane. Any of these can break it. As a *hypothesis* it is
coherent and testable. **This makes cosmic expansion / dark energy another thing the one medium + little arrows would deliver.**

**Multi-brane / Randall–Sundrum refinement (2026-06-24).** "Where does the drained medium go?" was long hand-waved as an *infinite
bulk vacuum* actively pulling fluid in. Two reframings + one structural commitment: (i) it's **one medium** — draining = *de-structuring*
(ordered brane → unstructured bulk ground state), so conservation is global/automatic; (ii) a codim-1 brane draining into codim-0 bulk is
a 3-surface into a 4-volume → dimensionally vast reservoir; (iii) **the structural commitment we'd dodged: is `w` infinite or
bounded?** A **stack of parallel branes** (the smectic itself) **partitions `w` into finite slabs** (RS1-like: two branes bound a finite
bulk, generalized to a stack). Then a throat on one brane drains into the adjacent slab and the **neighboring brane provides the RETURN
(re-absorption)** → the drain closes into a **steady circulation** (resolving the "drain is an open/non-equilibrium system" worry). The
**`±w` puncture direction = which adjacent slab = charge sign**; the **inter-brane spacing sets gravity's coupling/range**. So the brane
stack does triple duty: light substrate + bulk partition + conservation closure. **Honest constraint (able-to-fail):** recovering
long-range `1/r²` gravity from a finite slab is the classic braneworld problem — it needs **warping / graviton-zero-mode localization**;
if that can't be arranged, the finite-slab picture breaks. **KEY:** the steady drain-flow field around a throat *is* the gravity field of
that mass — so solving conservation directly yields the gravitational profile of a particle.

> **⭐ MILESTONE (2026-06-25) — `pathA_29` (track-3 gate-1) DID exactly this, and it CHECKS OUT: `RETURN_RESIDUAL_PREDICTION`, tri-review-verified.**
> Modelling the brane↔bulk return as a finite slab with a **DC-sink completion** (de-structuring/absorbing, cross-checked against a
> periodic Bloch stack) and taking the drain as the **premise** (gravity = inflow, `Z<0`), the steady-conservation solve gives **two
> results at once:**
> 1. **Long-range `1/r²` Newtonian gravity SURVIVES the finite slab.** The transverse zero mode is normalizable for the localizing
>    (flat/RS) geometry → `p=2` (solved from a real 3D-radial equation, counterfactual-guarded), so the braneworld-localization worry
>    above is — *for this family* — resolved in favour of survival. (It is genuinely able-to-fail: a *delocalizing* warp fails to
>    localize the zero mode → `p=3` → `RETURN_NOGO`, kept reachable.)
> 2. **But gravity does NOT come for free.** A drain that *breathes* radiates an **unavoidable, bounded monopole/dipole `c_s`-wave
>    `∝ ε0 = 1−𝒯₀(0)`**, tied to the gravity strength itself. GR forbids this (Birkhoff's theorem relies on mass conservation, which the
>    drain breaks on the brane) — so it is a clean **falsifiable prediction / departure from GR.** *You cannot have the Newtonian sink
>    without the radiation — they are the same flux.*
>
> This is the first end-to-end confirmation that **gravity-as-a-drain is self-consistent** (Newtonian `1/r²` + the GR-matched quadrupole
> + global conservation, no internal no-go) **and** yields a concrete testable difference from GR. It sharpens — does **not** close —
> `pde_ledger` open-item #9 (the full nonlinear return closure is the downstream track-3 work). Detail:
> `software/stage1_solver/{directives/pathA_29_brane_bulk_return.md, tools/pathA_29_*, reports/pathA_29_*}`.

---

## 6. Pedigree — explored before, never wired together

Individual pieces of this picture have respectable prior lives; the novelty is the *combination*. (A synthesis of pedigreed
parts is a reason to take it seriously — it is **not** evidence it is true; pieces being individually sound is exactly what would
make a *wrong* synthesis seductive.)
- **Light as elasticity of space:** MacCullagh's rotational ether (1839), Kelvin's gyrostatic ether (1890), Cosserat/micropolar
  elasticity (1909).
- **A particle as trapped field energy:** Wheeler's geon (1955).
- **Gravity/EM as emergent from a medium:** analog gravity — Unruh (1981), Volovik, Barceló–Liberati–Visser.
- **Our space as a wall in a higher-dimensional bulk:** Rubakov–Shaposhnikov, *"Do we live inside a domain wall?"* (1983); the
  superfluid **³He A–B interface** (Volovik) as a real, structured phase interface with surface physics.

### 6.1 Our OWN prior internal work on the throat — already explored, never wired together (don't rediscover it again)

The throat-as-self-bound-particle picture above has been **partially built twice in this repo** and then fragmented. Read these before
re-deriving (surveyed 2026-06-24):
- **Paper-7 — `notes/inner_throat/`** (`inner_throat_4d.md`, `4d_next_steps.md`, `inner_throat_hard_mode.md`): already builds the
  **force balance** `∂_a E = ∂_L E = 0`, `E_total = E_fluid + U_wave + E_geom`, `m_eff = E_total/c²`, in **4D with the calibrated n=5
  EOS**, plus a **working 4D GNLS solver + equilibrium-scan harness** and a drain diagnostic `Φ_w=∫ρv_w`. **STALLED:** it could never
  **self-select the aspect ratio `L/a`** (conservative-static tension+vacuum cost doesn't pin it; "requires dynamics"). The
  drain/backpressure sector was **specced but never built**, and the trapped wave was a **scalar surrogate**, not a brane-shear mode.
- **Lepton functional — `notes/lepton_mass_notes.md`**: the trapped-wave energy `F(a)=A/a + B/a² + C a³`, a virial identity, and
  **`m = (18/11)·E_trapped-wave`** with an *exact* finite-throat standing-wave spectrum. **But** the mass **tower is falsified** (1:9:25
  ≠ 1:207:3477) and the absolute scale is open — quantitative trapped-wave-as-mass, **no** lepton spectrum.
- **`software/stage1_solver/` throat (C0–C0g)** — a **different object**: a **frozen-geometry** gauged-quintic-GPE + Maxwell *profile*
  BVP solved only to feed the Path-A **gravity-sector closure** (`m̂0²·S_port` scale map / the `54/5` GR quadrupole). Throat geometry was
  a **frozen input — self-selection explicitly forbidden**; it hit a **fold at τ≈0.029** (a wall-elasticity knob, **not** the soliton's
  stability boundary) and was **deliberately re-scoped OFF the critical path** (`decisions/13` §5e–§5f). **Reusable:** the validated
  GPE+Maxwell solver, the stationary-GPE/soliton-ground-state **literature survey**
  (`reports/pathA_throat_solver_literature_synthesis.md`), and the conditioning/fold diagnostics — **not** a soliton existence result.
- **THE RECURRING WALL (both efforts):** *nothing in the built models self-selects the throat size/shape.* The 2026-06-24 candidate fix
  = the **trapped-wave pressure + drain backpressure + multi-brane return** forces (§4, §5) — **UNTESTED**. The cheapest existence test,
  **never run in either place**, is a **4D Derrick/virial scaling check** — do that *before* any big resume.

---

## 7. Open questions / conceptual audit (what is NOT yet pinned)

Honest list of what the picture still owes us. (Concept-level; the math test plan is §8.)
1. **Does the polar-arrow order actually form a stable wall (§2)?** Polar (`+w ≠ −w`), not headless — and stable, with a finite
   surface tension. The current `U(ρ)∝ρ⁵` is single-well; the wall must come from the *orientation* field. **Honest prior (GLM
   tertiary, 2026-06-23): there is likely a three-way no-win** — {light needs an orientational OP} × {a stable `±w` wall needs
   disconnected vacua, i.e. an easy-axis anisotropy} × {emergent `w` needs full rotational symmetry} can hold at most *two at once*
   for a classical OP. The anisotropy that buys stability is what *imposes* `w` — so the "**why space is 3D falls out for free**"
   bonus is the leg most likely to be lost (realistic best case: a stable wall with `w` put in by hand).
2. **Does the wall carry light, i.e. MacCullagh *rotational*-elastic shear, not mere Frank orientation waves (§2, §3)?** Not "does a
   shear modulus exist" but the sharper fork: are the in-plane modes genuine **rotational-elastic** (curl-only, traction, 2
   transverse + no longitudinal = a photon), or merely **Frank/orientation** waves (torque, not light), or **Cauchy** elasticity
   (a stray longitudinal third mode)? Make-or-break for the light sector; turns on whether the arrows couple to a genuine medium
   displacement (see the §3 deep-dive). Honest prior: a real coin-flip, hard at the curl-only + gauge-closure step.
3. **Dark-energy cycle signs & rates (§5):** is the net bulk↔brane flow *inward* (expansion)? Does influx mean *area growth* (not
   *densification*, which would drift `c_s`)? Does the matter→dark-energy crossover timing come out sane?
4. **Electron/proton ↔ charge-sign preference (§4):** does defect structure prefer a `±w` direction via bulk interaction, and how
   does that square with antiparticles + the matter-abundance asymmetry?
5. **Spin:** not yet placed in the picture. Where does it live (circulation in the trapped wave? throat geometry? the arrows?)?
6. **Charge magnitude:** the picture gives charge *universality* and *two signs* naturally; what fixes the actual *magnitude* `e`
   (the tension? the throat? the bulk)? — the calibration target.
7. **The defect's two roles at once:** a throat is simultaneously a **drain** (gravity/mass-energy inflow) and a **puncture**
   (charge). One shared energy budget or two co-located things? (Keeps the single-medium picture honest at the one place it
   matters most.)
8. **Cone lock `λγ`:** does light (brane shear, `c_γ`) end up at the same speed as the gravitational ripple (`c_s`)?
9. **Why does light come in packets (localization + quantization)?** A separate, deeper frontier (see §3): a linear wave spreads, so
   "photon as a lump" needs either nonlinear self-trapping (the optical-soliton analog — the best lead, but in tension with photon
   non-interaction) or genuine field quantization (`E=ℏω`, which the classical model doesn't produce). Deliberately out of scope for
   `pathA_24`; the standing puzzle the user has long flagged as missing a clean analog. *Candidate answer (2026-06-24, §4):* the THROAT
   provides the localization — a trapped mode "larger than the neck" is below waveguide cutoff → evanescent → genuinely bound; the lump
   is light caught in a throat, not a free wave that must self-trap in open space.
10. **Wall-energy cosmology (the σ-tension squeeze).** A stable wall with tension `σ_brane` is a classic cosmological hazard
    (Zel'dovich–Kobzarev–Okun overclosure). And the tension big enough to set the charge scale (§4) vs. tiny enough not to overclose
    differ by many orders of magnitude — a likely contradiction the dark-energy story (§5) must resolve by being the *flow*, not the
    static tension. (Tested in `pathA_24` G10.)
11. **Throat self-selection — the recurring wall (§6.1).** Neither Paper-7 nor the `stage1_solver` solve self-selects the throat
    size/aspect-ratio `L/a` (a named unsolved residual / a frozen input, respectively). The 2026-06-24 candidate force is the
    **trapped-wave pressure + drain backpressure + multi-brane return** (§4, §5) — UNTESTED. **Cheapest existence test = a 4D
    Derrick/virial scaling check** (never run); do it before any big resume — a Derrick obstruction would kill the throat-soliton cheaply.
12. **Does light's confinement (shear-free bulk) and its 2 polarizations actually fall out (§4)?** The intrinsic-(3+1)D / extrinsic-
    curvature picture *predicts* 2 polarizations + brane confinement for free, and re-frames the old "leak" as curvature, not loss — but
    this must be shown, not asserted (and it must still resolve the C5 longitudinal-mode obstruction at the throat, not in vacuum).

---

## 8. Next: the consistency-gate program for the GNLS polar-smectic superfluid (post-pivot plan)

**This replaces the old little-arrows T1–T5 derivation ladder** (T1 was run and FALSIFIED — see the §2 v3 update; the T2–T5
rungs are moot on that baseline). The plan and detail now live in **`docs/medium_requirements_and_prior_art.md`** (requirements
list A/B/C, prior-art survey, candidate structure, consistency gates). Summary:

**Goal:** test whether the **GNLS polar-smectic superfluid** (GNLS medium — KEPT — + polar orientation field + a non-local/roton
layering driver) is **internally self-consistent across all requirements**, i.e. whether one postulated structure satisfies the
whole list at once, or whether two requirements form a **no-go**. (Analog, not derivation — rule §0.6. Falsification = a no-go.)

**The make-or-break consistency gates (replace the T1–T5 ladder):**
- **Gate L — light on the smectic layer (THE CRUX):** the in-plane polar order must give MacCullagh rotational stiffness that is
  2-transverse + no-longitudinal, **bounded-below / inertially anchored** (defeat the Kelvin-gyrostat negative-energy instability),
  and **leak-free into the inter-layer bulk** (don't kill magnetism). A plain smectic is liquid in-plane — the rigidity must come
  from the arrows. Highest risk; the most likely no-go.
- **Gate S — magnetism preserved** (in-plane stiffness confined to the layer; inter-layer bulk shear-free).
- **Gate B — brane↔gravity compatibility** (the layering interaction is finite-`k`; must not disturb long-wavelength `c_s`/flow or
  the existing GR-quadrupole bundle `χ_Q`/`P0`).
- **Gate Q — two charge signs** (polar `±w` flip across the layer → two mirror, mass-independent puncture directions).
- **Gate K — cone-lock `c_γ≈c_s`** (likely a calibration gap, not a derivation — consistent with `λγ`).
- **Gate T — throat/mass** (a defect puncturing the layers → finite trapped-wave throat).

**Inherited walls (concede, don't fight):** dynamics/`G`/`α` (calibrate-predict — universal wall, not ours); emergent-axis/why-3D
(smectic helps but full isotropic 3+1D localization is structurally hard); Lorentz/preferred-frame (toy-analog scope only).

**Methodology shift:** from "freeze a MINIMAL postulate + test DERIVATION" → **"specify the FULL candidate structure (postulated
freely) + test CONSISTENCY / hunt a no-go,"** same pipeline (Codex design-review xhigh → GLM tertiary → execute, each tri-reviewed:
orchestrator arbiter re-run + transliteration-fidelity audit + adversarial review on clean agents + user gate; Mathematica leads,
SymPy cross-checks; dual-engine; units restored; Codex codes, Claude reviews). T0 (GNLS + polar OP, frozen `8fa41ac51e88`) is
preserved and extended with the layering driver.

**First action when we come back:** (1) read this doc + `docs/medium_requirements_and_prior_art.md`; (2) draft a fresh directive
encoding the GNLS polar-smectic structure + the consistency gates above (Gate L first — it's the crux/most-likely no-go);
(3) Codex design-review (xhigh) → GLM tertiary → execute gate-by-gate, tri-reviewed, user-gated.

---

## Cross-references
- EM re-founding detail + the MacCullagh/shear math history: `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md`.
- The directive to be reworked into the T1–T5 ladder: `software/stage1_solver/directives/pathA_24_brane_existence_defect_structure.md`.
- Program state / "you are here": `STATUS.md`; `software/stage1_solver/decisions/13_emergent_constants_derivation.md` §0.
