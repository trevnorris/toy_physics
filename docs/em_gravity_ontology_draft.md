# Native EM + gravity ontology — assembled checkpoint (2026-07-11)

> [!CAUTION]
> **ARCHIVED HISTORICAL CHECKPOINT — NOT THE CURRENT ONTOLOGY.** This document is preserved for provenance only. It contains superseded statements, including older identifications of trapped standing-wave content with observed mass. The canonical model specification is [Ontology and Closure Ledger of the Toy Analog Physics Model](toy_model_ontology_summary.md), which treats the trapped mode as structural support and keeps active, passive, and inertial mass operationally distinct.

**Purpose.** This is the picture we (user + Claude) have hammered out after the gauge/Higgs program was found to be a *mis-import*. It is a CHECKPOINT, not the final doc: the plan is to (a) mine the review logs for banked verdicts, (b) have Codex + Grok check the mined verdicts against THIS ontology, and only then (c) rewrite the real doc from scratch. Everything here is tagged **[COMPUTED]** (a run result we trust), **[HYPOTHESIS]** (proposed, not yet earned), **[OPEN]** (unresolved), or **[GRAVEYARD]** (killed — do not resurrect).

---

## 0. The core correction that reorganizes everything

**Light is NOT electromagnetism in this model.** In standard physics the photon *is* the EM gauge field, so "light = EM." That identity does **not** survive translation here, for a reason internal to the model:

- charge exists only as a **puncture** (charge sign = the `±w` direction of the puncture);
- a puncture is a **throat** = a trapped/standing configuration of the medium = it has **mass**;
- therefore **anything carrying charge is massive**;
- **light is massless** ⇒ light carries no charge ⇒ light cannot be an excitation of the charge/puncture (EM) sector. **[COMPUTED-style logical result; user's proof]**

Consequence: there is **no background EM gauge field** for light to be a ripple of. **Light is a ripple of the *brane*; EM is what *throats do to each other*.** Different substrates. The entire "EM must be a first-class gauge `U(1)` whose photon is light" program (gauge structure, second-class vs first-class, Higgs/deconfinement) was **answering a question this model does not ask.** → see Graveyard.

---

## 1. What each thing IS (native substrates)

| Thing | What it is | Status |
|---|---|---|
| **The medium** | ONE compressible superfluid in 4+1D, in two phases: ordered **brane** (our 3D space) + **bulk** | foundation |
| **Light** | untrapped **brane shear** — a brane elastic (rotational/MacCullagh) wave; 2 transverse modes | **[COMPUTED]** `pathA_36` transverse modes are correct |
| **Mass** | **trapped** brane shear = a geon = a **throat** | foundation |
| **Charge** | the **`±w` puncture direction** of a throat | foundation |
| **EM force** | the interaction between throats' **4D bodies** — static = charge/Coulomb, moving = magnetism | foundation + `pathA_38` |
| **Current** | a **throat in translation** (see §3) | definition (lock into docs) |

**Light IS allowed to be a brane elastic mode** — because light is not trying to reproduce Maxwell's lock; it is just free shear. **EM is NOT allowed to be a brane elastic mode** — that is the dead end (§5).

---

## 2. Why light bends in gravity but NOT in the EM force

- **Gravity** = brane **compression/density change** → it literally alters the medium that carries light → the shear wave refracts → **light bends** (lensing). Same medium, so light feels it. **[COMPUTED-consistent]**
- **EM force** = a **throat–throat 4D-body interaction** → it does *not* change the brane's shear stiffness → a passing shear wave has nothing to refract off → **light does not bend in an E/B field.** Different substrate.
- **Honest caveat (matches reality):** light still **scatters off charges as physical defects** (Thomson/Compton) — a throat is a real obstacle in the brane — which is distinct from the EM-force deflection. A light beam crossing a magnet isn't bent, but electrons *do* scatter light. **[COMPUTED-consistent]**

---

## 3. Current, and how charge/light couple — natively (NO standard-physics smuggling)

**Current = a throat in translation**: the puncture *pattern* moving through the brane, its 4D body sweeping through the bulk. **[definition — put in docs verbatim]**
- It is **pattern-motion, not stuff-motion**: the throat is made of medium; when it "moves," the pattern translates while medium flows *through* it (a whirlpool crossing a pond). Consistent with "nothing is static, it's all flow."
- The **elementary** current is ONE moving throat. A wire is a *population* of moving `−w` throats against stationary `+w` throats.
- **NOT** "charges flowing through a neutral background" (that's the standard-physics image).

**How shaking a charge radiates light [HYPOTHESIS — mechanism to verify]:** a throat is a defect embedded *in* the brane; accelerating it **stirs the brane's shear field** → the brane radiates shear = light. Run backward: a light wave (shear) sweeping past a throat **stirs it** → light drives charges (antenna). Same one mechanical coupling (throat-embedded-in-brane), both directions. Reproduces "charges make light / light drives charges" **without** light being the EM field.

**Origin of magnetism's velocity-dependence = [OPEN].** It must be *derived* from the moving 4D throat-body (the native `pathA_39` re-do). It is **NOT** to be asserted via a "moving throat looks compressed in the bulk" story — that is a cousin of a killed idea (§5) and has not earned a place.

---

## 4. Where EM lives, and how it gets `1/r²` (the load-bearing geometry)

**EM is a BULK interaction focused around the brane — it does NOT propagate IN the brane.** The distinction is the whole escape from the wall:
- **"propagates in the brane"** = a mode built from the brane's own elastic stiffness (curl/div) → `pathA_36` second-class → **dead (§5).**
- **"bulk interaction bound near the brane"** = a mode that is fundamentally a **bulk** DOF but whose profile is *peaked at the brane* (the `h`-branon, `sech²(w/ℓ)`). **[COMPUTED]** `pathA_38`.

**`1/r²` comes from localization = a bound zero mode, NOT fractional compression.** **[COMPUTED + reviewer-endorsed]**
- A mode trapped in `w` (width `ℓ`, `sech²` profile) lives in an effectively 3D sheet → its field spreads in 3D → **`1/r²`**. Same physics as an electron bound in a quantum well / light in a fiber. Not math magic — a bound state doesn't leak into the dimension it's bound out of.
- **Falsifiable discipline:** a *true normalizable zero mode* → `1/r²` at ALL distances, **no crossover**. A *leaky resonance* → crosses over to `1/r³` at large `r` (that would be a fudge). Observation shows no Coulomb crossover ⇒ it must be a true bound mode; the `sech²` says it is.
- **[HYPOTHESIS]** Gravity may be the *opposite* case — a leaky mode that DOES crossover at large `r`, which is where a MOND-like `1/r` tail would come from. Charge bound, gravity leaky = a real checkable distinction.

---

## 5. The lock (Maxwell structure) — the central OPEN problem

Maxwell locks curl-stiffness to div-stiffness in a fixed (relativistic) ratio; that lock IS the gauge structure and the source of the dual sign. A generic elastic continuum has these as **independent** knobs → no lock → not Maxwell. This is why every brane-elastic attempt returned second-class.

- **[GRAVEYARD] Lock from the brane's ELASTIC fields** — dead. `pathA_36` (second-class longitudinal) + `pathA_35` (`FAIL_COUPLE_STRESS_NOGO`) already computed this fails.
- **[HYPOTHESIS — leading, to test] Lock as a geometric/topological identity** from the brane-in-bulk **embedding + throat topology.** In real physics the "invisible" first-class constraints ARE geometric identities (`∇·B=0`, Bianchi) — they cost no energy and carry no DOF, which is *why* they're invisible, and which is why they dodge the `pathA_36` (dynamical-stiffness) problem. A geometric identity is **derived, not placed** (answers the "we're just placing it" objection). The throat is a topological defect — exactly where such identities get non-trivial content (winding/flux).
  - Structural opening: the brane's in-plane compression speed is **not** the bulk `c_s` (they're plausibly different modes) — so there is room for a brane-associated sector distinct from gravity's compression. *(from requirements-review; verify in mining.)*
- **[HYPOTHESIS] The lock lives in the bulk via a state-change.** The constituents have an internal state: **unlocked on the brane, locked in the bulk.** The lock is present exactly where EM lives (throat 4D body, bulk) and absent where light lives (brane). Analogy class = **photoinduced phase transitions** / **nematic→smectic** (excitation/ordering locks previously-independent stiffnesses). No perfect known material; the *category* is ordinary.
  - **Possible unification:** the same brane/bulk state distinction that supplies the lock may ALSO be the potential well that localizes EM to the brane (§4). One mechanism, two effects. **[FLAG: over-determination risk]** — if one state-change is forced to do the lock AND the localization AND the drain, watch for the model's known "one thing doing too many jobs" killer.

---

## 6. Gravity plumbing, spin, gravitomagnetism

- **The throat is a passive topological puncture — it does NOT pump. [corrected]** Inflow (gravity) is a **persistent, dissipationless superflow** set by the winding around the defect (like frictionless circulation around a superfluid vortex) — no engine needed because superflow doesn't dissipate. The exact drain-vs-return-leak energetics = **[OPEN]** (do not assert a driver).
- **Spin = [OPEN].** Candidate: **quantized circulation around the throat** (vortex-like — not the throat material rigidly spinning, which it can't do). A circulating `±w` throat = a current loop = a magnetic moment (ties spin↔magnetism). **Do NOT collapse spin to "only gravitomagnetism"** — a reviewer flagged that as too hasty. *(verify in mining.)*
- **Gravitomagnetism ≠ EM magnetism, no conflict.** Gravitomagnetism = moving **mass** (gravity/flow sector, already in the PN ladder). EM magnetism = moving **charge** (throat sector). Same "static + moving-source" math shape, different couplings/sectors — parallel structures.

---

## 7. GRAVEYARD (killed — do not resurrect)

- **[GRAVEYARD] "Light = EM" / EM as a mediating gauge field with its own photon.** Wrong ontology here (§0). The whole first-class-`U(1)` / second-class / Higgs / deconfinement quest answers a question the model doesn't ask.
- **[GRAVEYARD] Higgs/deconfinement reframe** (former §9) — fatal, unanimous 3-panel: `pathA_36`'s spectrum is gapless (not Higgsed); Anderson-Higgs forbids a massless bulk photon on a Higgsing brane; the constraint-class "flip" is circular/degenerate.
- **[GRAVEYARD] geon/MacCullagh EM direction** (`em_sector_reconsideration` v2) — computed no-go; dual sign only `BY_TUNING` (`B_eff→0`, which kills gravity).
- **[GRAVEYARD] Fractional-length compression `R ∝ r^{2/3}`** as the `1/r²`-vs-`1/r³` fix — killed (unmotivated, fine-tuned, scale-breaking, leaves fingerprints). Localization (§4) is the right route.
- **[GRAVEYARD] Brane elastic fields forming a first-class lock** — `pathA_36`/`pathA_35`. (§5)
- **[GRAVEYARD] "Moving throat looks compressed in the bulk" as the mechanism for magnetism** — not verified; magnetism's velocity-origin is an OPEN derivation (§3), not a compression story.
- **[GRAVEYARD] Spin → only gravitomagnetism** — flagged too hasty; keep spin open (§6).

---

## 8. The concrete next question (once ontology is validated)

Derive **`1/r²` + the dual sign + magnetism directly from the throat–throat 4D interaction** — i.e. the native `pathA_39` re-do — with **no gauge field anywhere**, and test whether the **lock** (§5) can come from geometry/topology (the leading hypothesis) rather than being placed.
