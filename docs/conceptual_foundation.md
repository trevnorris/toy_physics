# Conceptual Foundation — what each part of the model physically IS

**Read this FIRST.** This is the canonical plain-language statement of the model's physical vision — what the medium, the
brane, the four force-sectors, and the particle/defect actually *are*, in the model's **own native terms**. It exists so the
conceptual picture is never lost between work sessions and never has to be re-explained from scratch. It is a **living document**:
update it whenever the conceptual picture sharpens. It is deliberately separate from the math/derivation/verdict machinery
(those live in `software/stage1_solver/decisions/*`, the `pathA_*` directives, and `STATUS.md`).

**Status:** v2, 2026-06-23. Reconstructed from the pre-`/compact` "what is light" (drumhead/surface-tension) and "why does the
brane exist" discussions, plus the EM re-founding work (pathA_23) and the puncture/charge/mass picture. **v2 adds the agreed
working hypothesis** — the **little-arrows (polar-orientation) brane mechanism** (§2) — and the **bulk↔brane cycle / dark-energy**
extension (§5), plus the post-`/compact` test plan (§8).

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
   `pathA_24`; the standing puzzle the user has long flagged as missing a clean analog.
10. **Wall-energy cosmology (the σ-tension squeeze).** A stable wall with tension `σ_brane` is a classic cosmological hazard
    (Zel'dovich–Kobzarev–Okun overclosure). And the tension big enough to set the charge scale (§4) vs. tiny enough not to overclose
    differ by many orders of magnitude — a likely contradiction the dark-energy story (§5) must resolve by being the *flow*, not the
    static tension. (Tested in `pathA_24` G10.)

---

## 8. Next: testing the little-arrows hypothesis (the post-`/compact` plan)

**Goal:** find out, with math/code, whether the **little-arrows (polar-orientation) mechanism** actually delivers — first the
brane, then light, then charge, then (ambitiously) the dark-energy cycle. Falsification is the point: a clean FAIL at any rung is
a first-class, welcome result, never to be rescued.

**The test ladder (cheapest-concept-fatal-FIRST; one rung at a time, user gate between):**

- **T1 — Does a polar-orientation field form a stable wall? (the brane make-or-break.)** Add the *minimal* polar "arrow" order
  parameter to the medium and ask: does it produce a codimension-1 wall with **finite surface tension** and **bound zero modes**
  (confinement)? Two sub-questions that are themselves make-or-break: **(i) stability** — a vector whose vacua fill a sphere gives
  *no* topologically stable wall (`+w` can unwind to `−w`), so what minimal structure gives a genuinely stable (or adequately
  long-lived) `±w` wall — a polar double-well / anisotropy / disconnected vacua (π₀=ℤ₂)?; **(ii) is `w` emergent or pre-existing**
  — does the alignment *define* `w` (spontaneous), or did we smuggle in a preferred `w`-axis (anisotropy)? Report honestly which.
  FAIL modes: no stable wall; needs ad-hoc structure; `w` had to be imposed.
- **T2 — Does the wall carry in-plane shear → light? (the light make-or-break.)** Compute the spectrum of fluctuations localized
  on the wall. Does it contain **two transverse in-plane shear polarizations at a single speed** `c_γ=√(μ_brane/ρ_brane)`,
  cleanly separated from the (scalar) bending mode `u_w` and any trapped scalar modes? Wrong count / wrong dispersion / no shear =
  FAIL (this is the pathA_23 Stage-2 wall, now with a concrete candidate to test).
- **T3 — Bulk stays shear-free? (magnetism consistency.)** Confirm the bulk (arrows `±w`) acquires no in-plane shear rigidity, so
  Magnus/magnetism is untouched. Consistency gate.
- **T4 — Puncture into `±w` domains → two charge signs + finite throat.** Model a throat through the wall connecting to the `±w`
  bulk; confirm two mirror-image puncture orientations (the two charge signs, pure direction — no winding), and a finite throat
  radius from the tension-vs-holding-open balance.
- **T5 — The dark-energy cycle (ambitious; after T1–T4).** Set up the two opposing flows (matter-driven brane→bulk drainage;
  tension-throttled bulk→brane areal leak) and test the §5 claims: net inward leak, area-growth (not densification), the
  matter→dark-energy acceleration crossover. Signs and rates are everything here.

**Methodology (the standing process):** this supersedes/reworks the generic-domain-wall directive `pathA_24` — **rework pathA_24
to encode the little-arrows mechanism + this T1–T5 ladder**, then run it through the standing pipeline: Codex design-review
(`-c model_reasoning_effort=xhigh`) → GLM tertiary (foundational) → execute **stage-by-stage**, each **tri-reviewed** (orchestrator
re-run + transliteration-fidelity audit + adversarial review on clean agents) with a user gate. Engine: Mathematica leads (wall
profiles, fluctuation spectra, stability), SymPy cross-checks; dual-engine; units restored; Codex codes, Claude reviews.
Conditional-verdict rule on any postulated ingredient.

**First action when we come back:** (1) read this doc top-to-bottom for grounding; (2) rework `pathA_24` to test T1 (the brane
make-or-break) with the little-arrows order parameter; (3) Codex design-review (xhigh) → GLM → execute T1. Don't skip ahead —
T1 (stable wall) gates everything; if the polar field can't make a stable wall, the rest of the ladder is moot and that itself is
the result.

---

## Cross-references
- EM re-founding detail + the MacCullagh/shear math history: `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md`.
- The directive to be reworked into the T1–T5 ladder: `software/stage1_solver/directives/pathA_24_brane_existence_defect_structure.md`.
- Program state / "you are here": `STATUS.md`; `software/stage1_solver/decisions/13_emergent_constants_derivation.md` §0.
