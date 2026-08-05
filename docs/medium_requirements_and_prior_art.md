# Medium requirements + prior-art survey

**Purpose.** A clear, constituent-level list of what the one medium must BE and DO for the conceptual model to work, plus a survey
of whether existing physics already builds a system matching the list — and, critically, **what stopped each prior program from
going all the way** (the failure modes are the most useful signal).

**Status:** v1, 2026-06-23. Born from the pathA_24 **T1 falsification** (a *static* polar-vector domain wall cannot self-localize:
the O(4) vacuum manifold `S³` is connected, so `+w` unwinds to `−w` with zero barrier — `T1_FAIL_NO_STABLE_WALL`, tri-reviewed
genuine, commit `2fa91886`). Two user reframes drove this doc:
1. **Dynamic, not static.** The medium is intrinsically messy / self-interacting / time-dependent; "lowest-energy static
   configuration" was the wrong question. A dynamically-sustained brane (³He A–B interface, vortex, active wall) is the natural
   object, not a static energy minimum.
2. **Model the particles, let the rest emerge.** Stop postulating a continuum order-parameter field tuned to the answer (the
   recurring "math drifts from concept" trap, [[project-single-medium-concept-vs-math-drift]]); instead specify what the
   constituent particles are/do and let the brane/light/charge **emerge** (³He spirit). Pedigreed emergence > postulated field.

Canonical conceptual source: `docs/conceptual_foundation.md` (read first). This doc is the *requirements spec + prior-art map*.

---

## The requirements list

### A. What the constituents ARE (the "tiny particles")
- **A1.** Discrete particles with their own cohesive inter-particle forces, whose collective state is a **compressible superfluid**
  (a condensate with a coherent phase → frictionless flow + quantized circulation).
- **A2.** They fill **4 spatial + 1 time** dimensions.
- **A3.** Each particle carries a **polar internal orientation** (a genuine head≠tail vector, `+w ≠ −w`), and orientations
  **interact** (alignment lowers energy). The orientation is a property *of the particles* — carried by them, not a separate field.
- **A4.** The system is **intrinsically dynamic / self-interacting / time-dependent** — persistent flow and fluctuation are the
  normal state; static equilibrium is not the right idealization.
- **A5.** *(Lens)* The constituents may obey richer rules than written; a puzzle in the emergent physics is a hint to **enrich the
  constituent model**, not automatically a flaw.

### B. What must EMERGE (each with its hard signature)
| # | Phenomenon | Required emergent signature |
|---|---|---|
| **B1** | Gravity | inflow/**drain** toward defects (sinks); test bodies carried by the flow; an effective attractive metric |
| **B2** | Magnetism | bulk **swirl/vorticity** (quantized vortices, Magnus) — and the **bulk must stay shear-free** |
| **B3** | Gravity-ripple speed `c_s` | compression/phonon wave, density-dependent |
| **B4** | The brane | a **stable, dynamically-sustained, codim-1 3-surface** with surface tension, that confines us (bound zero modes), **bulk on both sides**, axis **emergent, not hand-imposed** |
| **B5** | Why 3D | orientation picks **one** axis (defines `w`) → perpendicular surface is 3D |
| **B6** | Light | on the brane: **2 transverse, non-dispersive, curl-only (rotational-elastic / MacCullagh) waves, no longitudinal mode**; rigidity **on the brane, not the bulk**; needs constituents that **resist local rotation** with orientation **rigidly tied to displacement** (traction, not torque) |
| **B7** | Electric charge | throat **puncture direction `±w`** — binary, two signs, quantized/universal, **mass-independent** |
| **B8** | Mass | **trapped standing wave** (geon) in a throat; **finite** throat radius (tension vs holding-open) → finite self-energy; **extended**, never a point |
| **B9** | Spin | *(open)* must live somewhere — circulation in the trapped wave? throat geometry? the orientation? |
| **B10** | Dark energy | bulk↔brane medium **cycle**: net inward → **area growth** (not densification) → decel→accel crossover |
| **B11** | Cone-lock | `c_γ ≈ c_s` (λγ→1) |
| **B12** | Light packets | *(open frontier)* localization + quantization — soliton self-trapping vs field quantization |

### C. The non-negotiable rules
- **C1.** ONE medium — everything emergent; no separate fundamental field layered on.
- **C2.** No point particles — all defects extended/finite.
- **C3.** Model the **particles**, let macro behavior emerge — don't postulate a continuum field tuned to the target.
- **C4.** Falsifiable; a clean failure is welcome.

---

## Prior-art survey (COMPLETE — 2026-06-23)

Five parallel research agents (web + scholar). Each: (i) match grade vs the list; (ii) **what stopped that program** (the key
signal); (iii) routability for us. Per-candidate detail in the session transcript; distilled below.

### Per-candidate grades
| Candidate | Best-matched reqs | Grade | One-line on what stopped it |
|---|---|---:|---|
| **Volovik / ³He** (Fermi-point emergent universe) | C1/C3 philosophy; B1-kinematics; B4 wall; B2/B9 | **6/10** | gravity is **kinematics, no Einstein dynamics, `G` ill-defined**; A–B wall **externally maintained, not self-localizing**; `l̂` axis boundary-pinned; Lorentz only IR + preferred frame; 3+1 is input |
| **Analog/emergent gravity** (Unruh, BLV, BEC) | B1, B3, (B11) | **6/10** | **kinematics-without-dynamics** — acoustic metric solves Navier–Stokes, NOT Einstein eqs; no `1/r²` attraction, no derived `G`; only crack = Sakharov/Volovik induced gravity (fermionic sector + Λ-catastrophe) |
| **Wen string-nets** (emergent photon + electron) | **B6, B7, B8, B9** | **8.5/10** | **lattice substrate** → chiral-fermion no-go (Nielsen–Ninomiya), fermion doubling (4 species), **no gravity**, Lorentz fine-tuned, `α`+`c` are free knobs not predictions |
| **Orientational/rotational matter** (polar nematics, smectics, micropolar, MacCullagh) | A3, B5; templates for B4, B6 | **7/10** (ingredients) | classical/dissipative/non-vacuum; polar **domain walls confirm T1** (boundary-stabilized, micron-wide, spread); MacCullagh light = right form but **died of negative-energy/angular-momentum instability** (Kelvin gyrostat) |
| **Braneworlds / solitons** (Rubakov–Shaposhnikov, Dvali–Shifman, Q-balls) | B8 (fermion loc.); partial B4 | **4/10** | all **postulated, not emergent** (C1/C3); **light-localization is the field's famous open problem**; Q-ball needs conserved U(1)+non-single-well (excludes our GNLS); oscillons decay |

### The unified findings (what the whole survey says)

1. **A universal wall — and we already hit it: KINEMATICS without DYNAMICS.** Every emergent-gravity program (analog gravity AND
   ³He) reproduces the effective metric + geodesics but **not the Einstein equations and not a first-principles `G`**. This is
   *exactly* our Gate-4 result (`g_G` calibrated-on-`G`, not derived). ⇒ Our calibrate-predict concession is **not a solver
   shortfall — it is the correct, honest response to a 40-year field-wide wall.** Strong vindication of the methodology.

2. **The "emergent axis / why-3D" obstruction is real and STRUCTURAL — independently confirmed twice.** (a) The leading
   domain-wall-Standard-Model program (Davies–George–Volkas / Lifshitz) found that with **isotropic** 4-space the 3+1D zero modes
   **do not exist** — they had to **break 4D rotational symmetry by hand**. (b) ³He's `l̂` axis is spontaneous only in idealized
   bulk; in practice it is **boundary/field-pinned**, and the A–B interface is **externally maintained**. ⇒ T1 was not a modeling
   slip; the **"why space is 3D, for free" leg is the one prior art says is hardest/most-likely-lost** (matches the GLM three-way
   no-win). Honestly downgrade it from goal to long-shot.

3. **The one concrete NEW escape for the brane: the SMECTIC / lamellar mechanism.** Instead of a polar-vector *domain wall* (which
   T1 killed and real ferroelectric nematics re-kill), a **smectic** spontaneously breaks continuous *translation* in one
   direction → an **emergent stack of codim-1 layers as the GROUND STATE**, layer-normal emergent, stabilized by
   compression+curvature elasticity (de Gennes `B|∂u|²+K|∇²u|²`) — **qualitatively immune to the spread-and-unwind failure**
   because it is a broken-symmetry order, not a metastable kink. Caveats: gives a *stack* (many layers, not one); the order is
   translational/density (in-plane light B6 still unsolved on it); needs hosting on a *dynamic* medium to avoid LC dissipation.
   This is the survey's **strongest structural lead** for B4+B5 with emergence intact.

4. **Light (B6) is the crux, with two RIVAL pedigreed routes, each with a documented fatal-looking catch.** (A) **Continuum
   MacCullagh** rotational-elastic shear (our route): right signature (2 transverse, no longitudinal, from curl-of-displacement),
   but historically **died of a negative-energy / angular-momentum-non-conservation instability** — any import MUST prove the
   rotational stiffness is **bounded-below and inertially anchored** (the precise make-or-break check). (B) **Lattice gauge
   constraint** (Wen): rigorously works (8.5/10) but **requires a lattice** → preferred frame + fermion doubling + no chiral
   fermions + no gravity — rival to our continuum picture; adopting it means inheriting all of that. We are deliberately the
   continuum one.

5. **Charge (B7) gets independent corroboration.** Wen's "charge = quantized flux at the **end of an extended, orientation-bearing,
   non-isolatable string**" maps onto our "charge = puncture **direction** of an extended throat, never a point" — both get
   quantization + two signs + mass-independence from charge being a property of an *extended, non-isolatable* object. Convergent
   support for the puncture-charge picture (mine the analogy; do NOT import the lattice photon it rides on).

### Implications for our program
- **Methodology vindicated:** calibrate-predict + "demonstrate the bridge, don't derive the cosmos" is the right frame; the
  dynamics/`G`/`α` wall is universal, not ours. ([[project-analog-framework-goal]], [[feedback-calibrate-predict-methodology]].)
- **Where surplus can live (our structural differences from prior art):** MacCullagh continuum light (B6), puncture-charge (B7),
  geon-mass (B8) — these are exactly the places ³He/Wen/analog-gravity *don't* go. Inherited walls: dynamics/`G` (B1), emergent
  axis/why-3D (B5).
- **Most actionable next build (if we continue the brane):** a **smectic/lamellar** brane mechanism (a fresh T0), with the
  **MacCullagh negative-energy stability** as the make-or-break light check.
- **Closest spiritual ancestor:** Volovik. We differ (and can earn surplus) precisely by being a *polar-constituent + brane-shear*
  continuum rather than a *chiral-`l̂`* fermionic liquid.

---

## ⭐⭐ ATTRIBUTION PASS — the S9→S11b light sector, named piece by piece (2026-08-04)

⭐ **Purpose: paper attribution, and a second source of checks.** An independent literature search was run
against the S9–S11b step records. ⭐ It **re-confirmed** this survey's own identifications from an entirely
separate route (MacCullagh, Volovik, Rubakov–Shaposhnikov all came back), ⭐ and added names this survey
did not have.

⛔⛔ **THE HEADLINE, AND IT IS UNCOMFORTABLE: the S9/S10/S11 light algebra is NOT new.** Curl-only stiffness
⇒ `D−1` transverse modes at `c² = μ_R/ρ_br` + a zero-restoring-force longitudinal slot, and lifting only
the longitudinal with a compression modulus, is **standard 19th-century continuum mechanics.**
⇒ ⭐ **Whittaker, _A History of the Theories of Aether and Electricity_, Ch. V** is the reference; a reader
of the step records would currently take the mode structure as earned here, and it is not.
⚠ **This does not invalidate the rebuild** — we still had to derive it, and standard ground is *solid*
ground. ⭐ **But S9, S10 and S11 must cite it.** Citing it strengthens the toy-analog framing rather than
weakening it. ⇒ [[feedback-framing-split]]

| our piece | status | attribution |
|---|---|---|
| S9/S10 curl-only ⇒ transverse count, `c² = μ_R/ρ_br` | **STANDARD** | MacCullagh 1839; Whittaker Ch. V |
| S11 longitudinal lifted by a compression modulus | **STANDARD** | Cauchy–Navier elasticity |
| **S11b-A interface law** | **STANDARD** | ⭐⭐ **structural acoustics** — fluid-loaded plate, radiation impedance, added mass |
| S11b-B passivity / Onsager–Casimir, "odd" couplings need a drive | **KNOWN** | Fruchart, Scheibner, Vitelli, *Annu. Rev. Condens. Matter Phys.* **14** 471 (2023) |
| the `h`-branon as the brane's own transverse fluctuation | **KNOWN** | Cembranos, Dobado, Maroto, PRL **90** 241301 (2003) |
| defects ⇒ long-range fields with ⛔ no gauge field | **KNOWN** | Eshelby, *Solid State Phys.* **3** (1956) |
| MacCullagh **+ topological defects ⇒ EM including charge** | **KNOWN** | ⭐⭐ **Unzicker**, `arXiv:gr-qc/0011064`, *ZAMM* **102** (2022) — ⚠ **our charge mechanism, already published** |
| `Z₂` sign = ± orientation of a puncture | **ADJACENT** | domain walls, spin-ice monopoles, signed vortices — the exact construction not found under a name |

### ⛔ NOT FOUND — ⚠ a search-failure statement, ⛔ NOT a proof of originality

- MacCullagh stiffness restricted to an **ordered thin phase with a strictly shear-free bulk**, so light is
  **confined** and cannot leave — every piece is known; this **confinement architecture** was not found.
- The **stray longitudinal as the deliberate physical anchor for charge** — the interface mathematics is
  standard; this **role** for the mode was not found.
- **Gravity as drain-flow between throats while light is MacCullagh on the same sheet** — Volovik and
  Unzicker are adjacent, ⛔ not this.

⇒ ⭐ **The novel-looking residual, if any, is the ARCHITECTURE — ⛔ not any individual derivation.**

### ⛔⛔ AND THE CHECK THIS SURVEY ALREADY NAMED, WHICH THE LIGHT SECTOR NEVER RAN

⚠⚠ **This document flagged it on 2026-06-23 (see item 4 above and "most actionable next build"):**
MacCullagh light *"died of a negative-energy / angular-momentum-non-conservation instability — any import
MUST prove the rotational stiffness is **bounded-below and inertially anchored** (the precise make-or-break
check)."* ⇒ Kelvin-gyrostat instability, `arXiv:1907.04144`.

⛔ **S9, S10 and S11 built the MacCullagh light sector and did not run it.** They correctly compute the
**consequences** of curl-only stiffness; ⛔ nothing establishes that a medium **can have** curl-only
stiffness **stably**. ⇒ ⭐ That is a **supplied premise**, and this survey already identified it as the one
most likely to break. ⇒ [[feedback-whose-law-is-it]]

⭐ **It belongs in the step records as a stated limitation, and on the plan as a step.**

## ⭐⭐ THE 4D SPATIAL BULK — separately searched (2026-08-04)

⭐ **Asked narrowly: has anyone hypothesised the 4-dimensional SPATIAL structure, ⛔ not the light sector?**

⭐ **Yes — the geometry, since 1982–83. ⛔ Not the combination, and ⛔ not the charge mechanism.**

| ingredient | label | attribution |
|---|---|---|
| 3-sheet embedded in a bulk with an open extra dimension | **ADJACENT** | ⭐⭐ **Rubakov–Shaposhnikov (1983)**, *Phys. Lett. B* **125** 136, *"Do we live inside a domain wall?"*; **Akama (1982)**; RS2; DGP (2000) |
| bulk is a **material medium** with density and flow | **KNOWN** separately | superfluid-vacuum theory — Volovik, Sbitnev, Consoli, Huang, Fedi |
| matter = **topological puncture THROUGH** the sheet | **ADJACENT** | domain-wall localization traps matter **ON** the wall (potential / zero modes) — ⛔ the opposite construction |
| **charge sign = puncture orientation ±`w`** (`Z₂`) | ⭐⭐ **NOT FOUND** | searched charge/orientation/puncture/throat/membrane + `Z₂`, R7-branes, topological `Z₂` |
| gravity = **medium flowing into** the defect | **KNOWN** separately | Cahill (2003); ⭐ Hamilton–Lisle *river model*, *Am. J. Phys.* **76** 519 (2008); Sbitnev; Unruh acoustic sinks — ⛔ all **3D**, no off-sheet bulk |
| **shear on the sheet, none in the bulk** ⇒ light confined | **ADJACENT** | ⛔ no "shear only on a 3-sheet, fluid bulk, traps light" package found |

### ⭐⭐ THE DIFFERENCES — this is the paper's positioning section

| | prior work | ⭐ ours |
|---|---|---|
| the extra dimension | one more spatial leg of **relativistic spacetime**; **metric** bulk | a **material 4-space** with density and flow |
| the sheet | a field-theoretic **domain wall** | a **shear-stiff sheet** in a **non-shear** bulk |
| matter | **trapped ON** the wall by a potential / zero modes | a **hole THROUGH** the sheet into `±w` |
| charge | a **U(1) / SM quantum number** | **puncture orientation** |
| gravity | **5D Einstein / induced Einstein–Hilbert** | the medium **draining off** the sheet |

⇒ ⭐ **Not found assembled.** Braneworlds stop at an empty metric bulk with potential-trapped matter;
superfluid-inflow models stop in 3D with no extra spatial direction and no punctures. ⚠ What stopped each:
**GR-as-geometry** on the braneworld side, **no extra spatial direction** on the superfluid side.

### ⚠ Two cautions before any of this reaches a paper

- ⛔ **Cahill's inflow work is regarded as fringe.** Citing it as precedent cuts both ways. ⭐ Prefer
  **Hamilton–Lisle**, who make the same flow point **inside standard GR**.
- ⭐⭐ **The strongest available claim is NOT "this is novel."** This survey already records that
  **light-localization is the braneworld field's famous open problem** (see the per-candidate grades,
  Braneworlds/solitons row). ⇒ ⭐ **our shear-on-sheet / no-shear-in-bulk mechanism is a DIFFERENT ANSWER
  TO A KNOWN OPEN PROBLEM** — a far stronger and more checkable claim than novelty in general.

⇒ ⭐ Add to the reading list: **Rubakov–Shaposhnikov (1983)**, **Akama (1982)**, **DGP (2000)**,
**Hamilton–Lisle (2008)**, **Maartens** braneworld *Living Reviews* as the baseline to contrast against.

## ⭐⭐⭐ THE DRUM-HEAD CHARGE MECHANISM — searched against its CURRENT formulation (2026-08-04)

⚠ **An earlier pass searched the OLD wording ("charge = `Z₂` puncture orientation") and got NOT FOUND.**
⭐ The current mechanism — **Coulomb from the elastic energy of superposed membrane deflections** — is a
**heavily studied** soft-matter problem. ⇒ [[project-puncture-deflection-charge-mechanism]]

⭐ **Its name:** *membrane-mediated* (or *curvature-mediated*) interactions between inclusions; *capillary
interactions* / the **Cheerios effect**. **STANDARD/KNOWN.** ⛔ Never cast as electric charge.

### ⭐⭐ THE SIGN IS NOT UNIVERSAL — IT IS SET BY THE BOUNDARY-CONDITION CLASS

| regime | inclusion type | like-direction | unlike |
|---|---|---|---|
| tension + gravity (Cheerios) | capillary **monopole** (net vertical force) | ⛔ **ATTRACT** | repel |
| tension-dominated, **force-free** | same force-pattern | ⭐ **REPEL** | attract |
| **bending only** (`σ=0`) | rigid **cones** (slope `α`) | **REPEL always** — ⚠ but energy `∝ α₁²+α₂²`, so **ORIENTATION-BLIND at leading order** | same |
| **bending + tension** | cones, **equal** orientation | ⭐ **REPEL at all `R`** (Weikl–Kozlov–Helfrich 1998) | — |
| bending + tension | cones, **opposite** orientation | — | ⚠ repel **near**, attract **far** |

### ⭐⭐⭐ THE BC ANSWER — this RESOLVES `R1_REQUIRED(bc_selection)` FROM OUTSIDE THE PROJECT

⚠ Our build stalled because the exterior **relaxes to zero deflection** unless something pins it, and
nothing selected the pinning class. ⭐ **Membrane physics has settled it:**

> **For a particle that PIERCES the sheet, the physically preferred BC is contact-line + prescribed contact
> angle, with the sheet free to set its height so that net force and torque vanish — i.e. a cone-like SLOPE
> BC, ⛔ not a free force monopole.**
> ⛔ *Fixed height* applies only if the puncture is **clamped to an external scaffold**; ⛔ *fixed force*
> only if something **external loads** the particle. ⚠ **Neither describes our puncture.**

⇒ ⭐⭐ **The physically-correct BC for a piercing particle is exactly the class that gives LIKE-REPELS.**
⇒ ⭐ Map onto our REPLACE/ADD variants before using this — the correspondence is **not yet checked**.

### ⛔⛔ TWO PROBLEMS THIS LEAVES, both sharp

1. ⛔ **Pure bending is ORIENTATION-BLIND at leading order** (`∝ α₁²+α₂²`) ⇒ `+w` and `−w` interact
   identically ⇒ **no charge sign at all** without tension. ⚠ And with tension, **opposite** orientation
   **repels near and attracts far** — ⛔ not Coulomb's behaviour.
2. ⛔⛔ **THE POWER LAW IS WRONG IN 2D.** ⚠ **No standard 2D-sheet-in-3D interaction gives force `∝ 1/R²`:**
   bending cones `U∼a⁴/R⁴` (`F∼1/R⁵`), anisotropic `U∼a⁴/R²` (`F∼1/R³`), tension+cones `U∼K₀(R/λ)`
   (**exponential**), Cheerios `∼K₀(R/L_c)` (exponential). ⇒ ⭐ **The sign wants the CONE BC; Coulomb's
   `1/R²` wants a MONOPOLE. In 2D you cannot have both.**

### ⭐⭐⭐ WHERE THE 4D SETTING DOES REAL WORK — AND THE CALCULATION TO RUN

> **"3-sheet in 4-space: NOT FOUND … Green functions of `∇²` on a 3-manifold are mathematically different
> (`G ∼ 1/R` for tension Poisson), so powers would change, but this case is UNSTUDIED. ⛔ Do not assume
> 2D-in-3D signs/powers carry over without re-derivation."**

⚠ A 3-dimensional sheet has a `1/R` Green's function where a 2D sheet has `log R` ⇒ **every power in that
table changes.** ⛔ **Which way is NOT established and ⛔ must not be assumed** — ⚠ the temptation to
hand-wave "so it becomes `1/R²`" is exactly the kind of step this project bans.

⇒ ⭐⭐ **THE CALCULATION:** redo the membrane-inclusion interaction for a **3-sheet in 4-space**, for **both
BC classes** (slope/cone and monopole), and obtain **the sign and the power law**.
⭐ It is bounded, well-posed, dual-engine-able, and has a **known 2D answer to calibrate against** — the
best-conditioned calculation available to the charge sector. ⭐ And it is **observational**: Coulomb's
`1/R²` is measured, so this either reproduces it or contradicts an experiment.
⇒ [[feedback-analog-find-consistent-structure]]

### ⭐ Reading list — the charge mechanism

1. **Bitbol, Constantin, Fournier**, *Membrane-mediated interactions*, `arXiv:1903.05712` — ⭐ power-law
   table + BC survey; read first.
2. **Weikl, Kozlov, Helfrich**, *Phys. Rev. E* **57** 6988 (1998) — ⭐⭐ **sign versus orientation and tension.**
3. **Goulian, Bruinsma, Pincus**, *Europhys. Lett.* **22** 145 (1993) — `1/R⁴` cones at zero tension.
4. **Vella & Mahadevan**, *Am. J. Phys.* **73** 817 (2005) — Cheerios; ⚠ like menisci **attract**.
5. **Evans, Turner, Sens**, `cond-mat/0301144` — force distributions, tension, `K₀` repulsion.
6. **Dommersnes & Fournier**; **Kim–Neu–Oster** — multipoles, anisotropy.

### ⭐ Reading list from the attribution pass, most useful first

1. **Whittaker**, *A History of the Theories of Aether and Electricity*, Ch. V — MacCullagh, Green, Kelvin.
2. **Unzicker**, `arXiv:gr-qc/0011064` and *ZAMM* **102** (2022) — closest living "MacCullagh + defects = EM";
   ⭐ read before further charge-sector work.
3. **Volovik**, *The Universe in a Helium Droplet* — closest living "one medium ⇒ gravity + EM".
4. **Barceló, Liberati, Visser**, "Analogue Gravity," *Living Rev. Relativ.*
5. **Cembranos, Dobado, Maroto**, PRL **90** 241301 (2003) — branons.
6. **Fruchart, Scheibner, Vitelli**, *Annu. Rev. CMP* **14** 471 (2023) — odd elasticity, passivity (S11b-B).
7. Any **structural acoustics** text on fluid-loaded plates — ⭐⭐ S11b-A *is* this subject, so the
   literature holds an **external cross-check** on our interface law, stronger than any review leg.
8. **Eshelby**, *Solid State Phys.* **3** (1956) — defects ⇒ long-range elastic fields, no gauge field.

---

### Sources
Volovik *cond-mat/9806010*, *0709.1258*, *The Universe in a Helium Droplet*; Barceló–Liberati–Visser *Analogue Gravity* (Living
Reviews) + *gr-qc/0106002*; Levin–Wen *cond-mat/0407140* (RMP 77, 871), Gu–Wen *hep-th/0507118*, Levin–Wen *cond-mat/0404617*;
ferroelectric nematics Chen–Clark *PNAS 2020* (arXiv 2003.03020), Mertelj–Lavrentovich *Nat.Commun. 2022* (arXiv 2206.06600);
MacCullagh ether Darrigol/O'Raifeartaigh (EPJ-H 2010), Kelvin-gyrostat instability arXiv 1907.04144; Rubakov–Shaposhnikov (1983),
Dvali–Shifman (1996–97), Coleman *Q-balls* (1985), Davies–George–Volkas arXiv 0705.1584 + arXiv 1008.2054, Zel'dovich–Kobzarev–Okun
(1975).

---

## The candidate structure — the **GNLS polar-smectic superfluid** (decided 2026-06-23)

### The reframe (what we are actually doing)
We are building a mathematical **analog**, NOT deriving the universe. We do **not** need to derive the medium/brane/constituents
from first principles — the constituents could be arbitrarily deep (e.g. 4D structures living in their own 5D space; unknowable).
The question is sharper and finite: **is there a single self-consistent superfluid structure that satisfies ALL the requirements
(A/B/C) at once?** We may **postulate the structure freely**; the only test is internal **consistency**.
- **This makes falsification STRONGER, not weaker.** The failure mode is no longer "our minimal postulate wasn't enough" (which
  invites endless structure-adding). It is a genuine **NO-GO: two or more requirements are mutually incompatible in *any*
  structure.** The three-way no-win (T1) was the first such no-go. **Hunting no-gos between requirements is the new adversarial
  target.** (Goal framing: [[project-analog-framework-goal]], [[feedback-framing-split]].)

### The structure
One compressible superfluid `ψ=√ρ e^{iθ}` (the **GNLS — KEPT as the carrier**) whose ground state spontaneously develops, at once:
(i) a **1D density layering** (smectic) → the stable, emergent-axis brane; and (ii) a **polar in-plane orientation** within the
layer (the "little arrows") → light + charge. Working name: the **GNLS polar-smectic superfluid**.

**KEPT — the GNLS medium, unchanged** (the existing program rides on this; survey-validated; T0 already froze GNLS + polar OP):
`ψ=√ρe^{iθ}`; quantum pressure `(ℏ²/8mρ)(∇ρ)²`; single-well EOS `U(ρ)=(K/4)ρ⁵` → `c_s`, `P(ρ)` (B3); flow / circulation / vortices →
gravity-drain (B1) + magnetism (B2); gauge coupling `(q_*/m)A_i`.

**ADDED — two new ingredients, same medium:**
1. The **polar orientation field** (arrows), carried by `ψ` (ρ-weighted, advected) — already in T0 → light (B6), charge (B7).
2. A **layering ("smectic") driver** — a **non-local / finite-range (roton-type) interaction** OR a **polar↔density coupling** that
   makes the dispersion develop a finite-`k` (roton) minimum that softens → a stripe/smectic ground state. The local single-well
   `U(ρ)` **cannot layer on its own**; this is genuinely new structure (the way real dipolar-BEC / soft-core supersolids stripe).

**Honest change from the original picture:** the **density NOW MODULATES** (the smectic layering). The earlier little-arrows claim
("the two states live in orientation, not density; `U(ρ)` stays single-well; we never modulate density") is **refined**: `U(ρ)` stays
single-well *locally*, but a new non-local term gives a finite-`k` density modulation. "Density doesn't modulate" is no longer true.
**Pedigree:** stripe/smectic superfluids with orientational order are real matter (dipolar-BEC supersolids; spin–orbit-coupled stripe
condensates) — we assemble known pieces, not invent physics.

### Requirement → mechanism map
| Req | Mechanism in the GNLS polar-smectic superfluid |
|---|---|
| B4 brane / B5 why-3D | **smectic layering** — codim-1 surface as the *ground state*, one emergent axis (fixes both T1 failures) |
| B6 light | **in-plane polar order** supplies MacCullagh rotational stiffness on the layer (the arrows' original job) |
| B7 charge | polar order flips `±w` across the layer → two mirror puncture directions |
| B1 gravity / B3 `c_s` | superfluid **flow → acoustic metric** + compression sound (kinematics; `G` calibrated — universal wall) |
| B2 magnetism | superfluid **vortices**; in-plane stiffness confined to the layer; inter-layer bulk stays liquid |
| B8 mass | defect **puncturing the layers** = trapped-wave throat, finite radius |
| B9 / B10 / B12 | spin from the orientational order / inter-layer flow cycle / nonlinear self-trapping (open, as before) |

### The consistency gates (the NEW make-or-break targets — hunt for a no-go), by risk
- **Gate L — light on the layer (THE CRUX, highest risk).** A plain smectic is **liquid in-plane (no shear → no light)**, so the
  rigidity MUST come from the polar order, and it must be (a) **2-transverse + NO longitudinal**, (b) **bounded-below / inertially
  anchored** (defeat the Kelvin-gyrostat negative-energy instability that killed the historical MacCullagh ether), and (c) **no
  shear leak into the inter-layer bulk** (the pathA_23 magnetism-killing leak). All three at once, or it's the no-go.
- **Gate S — magnetism preserved.** In-plane stiffness confined to the layer; inter-layer bulk shear-free (Magnus intact). Tied to Gate L(c).
- **Gate B — brane↔gravity compatibility.** The layering interaction is a finite-`k` feature; it must NOT disturb long-wavelength
  `c_s`/flow/drain or the existing GR-quadrupole bundle (`χ_Q`, `P0`). Verify, don't assume.
- **Gate Q — two charge signs.** Polar order flips `±w` across the layer → two mirror, mass-independent puncture directions. Looks consistent; verify.
- **Gate K — cone-lock `c_γ≈c_s` (B11).** In-plane rotational modulus vs bulk compressibility → two speeds, equal nowhere
  automatically (survey). Almost certainly a **calibration gap**, not a derivation — consistent with `λγ`'s current status. Flag, don't bank.
- **Gate T — throat/mass (B8).** A defect puncturing the layers → trapped-wave throat, finite radius; consistency with all the above.

### Inherited walls (NOT gates we can pass — honest scope; concede them)
- **Dynamics / `G` / `α`:** every emergent program gets kinematics, not dynamics (no Einstein eqs, no first-principles `G`). We
  **calibrate** these and predict surplus — the correct response to a universal wall, not a shortfall ([[feedback-calibrate-predict-methodology]]).
- **Emergent axis / why-3D:** the smectic makes the axis emergent (better than the domain wall), but the Davies–George–Volkas / ³He
  warning says fully-isotropic 3+1D localization is structurally hard — treat "why-3D for free" as a hoped-for bonus, not a goal.
- **Lorentz / preferred frame:** a medium has a rest frame; exact Lorentz invariance won't emerge. Acceptable under strict
  toy-analog framing; do not claim the real vacuum.

### Methodology shift (binding for the next directive)
From "freeze a MINIMAL postulate + test whether it DERIVES X" → to **"specify the FULL candidate structure (postulated freely — it's
an analog) + test internal CONSISTENCY across all requirements / hunt for a no-go,"** under the same discipline (dual-engine,
units restored, falsification-first, no-rescue, conditional-verdict; Codex codes / Claude reviews + arbiter re-run + fidelity audit
+ adversarial review + user gate). The **consistency gates above replace the T1–T5 derivation ladder.** T0 (GNLS + polar OP, frozen)
is preserved and extended with the layering driver.
