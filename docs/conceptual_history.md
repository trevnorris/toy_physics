# Conceptual History — the derivation ladder and how we got here

**The vision this records:** `docs/conceptual_foundation.md` (the plain-language, read-first statement of what the
medium, brane, four force-sectors, and defect physically ARE, in the model's own native terms).

This companion doc holds the **dated derivation history** that used to live inside `conceptual_foundation.md` and had
grown too large to keep there: the document's own version changelog, the per-gate / per-stage compute records, the dated
"⭐ vN UPDATE" blocks, and the **superseded-hypothesis detours** (the little-arrows domain wall, the density-smectic
route, the shear-surface no-go, the flow-framing retirement, etc.). Content is ordered roughly **chronologically**.
**Nothing here was deleted** — it was *moved* out of the foundation doc so that doc can stay a clean vision. Where a block
references another "§N" or a "v5 block in §2", those pointers refer to the foundation doc's live section or to the block as
it appears below.

---

## Document version changelog (the foundation doc's own status blocks)

These are the top-of-file `**Status**` paragraphs that tracked each version bump of `conceptual_foundation.md`.

**Status:** v3, 2026-06-23. (v2 added the little-arrows brane hypothesis.) **v3 records a pivot:** the little-arrows *domain-wall*
brane was **FALSIFIED** at pathA_24 T1 (a static polar-vector wall on the connected `S³` vacuum manifold spreads/unwinds —
tri-reviewed genuine, commit `2fa91886`). A prior-art survey + two user reframes led to the **GNLS polar-smectic superfluid**
candidate (the arrows now ride on a **smectic density-layering** brane, not a self-localizing wall) and a **methodology shift**:
we are building an **analog and testing whether one self-consistent structure satisfies all requirements** — NOT deriving the
universe. Full requirements list, prior-art survey, candidate structure, and the consistency-gate plan now live in
`docs/medium_requirements_and_prior_art.md` (read alongside this doc). Changes below: §0 reframe rule; §2 brane refinement (density
now modulates); §8 new plan.

**Status update:** v5, 2026-07-03. **v5 records a second brane pivot.** Two more brane routes have since died (both kept as
first-class falsifications): the v3 **smectic density-layering** route is **CLOSED** (`RC_DENSITY_SMECTIC_LIGHT_NOGO` — a density
modulation forms a codim-1 brane only by pinning the arrows along the layer normal, which *starves* in-plane light), and the
fallback **postulated shear-surface** brane hit `FAIL_COUPLE_STRESS_NOGO` (pathA_35 Gate L, provisional — the frozen light package
has no scalar potential `φ` to kill its longitudinal zero mode). The **new candidate** is a **material-state closure**: the brane
is the medium's *ordered, shear-supporting phase*; the bulk is the *same medium de-structured / shear-free*; throats are
phase-conversion sites. See the ⭐ v5 block in §2 and the updated plan in §8. This *formalizes* the "draining = de-structuring"
language already in §5. Source: `notes/brane_bulk_handoff.md`.

**Status update:** v6, 2026-07-03. **v6 = the medium DEMONSTRABLY carries light, and the "third wave / C5" is where the sectors
CONNECT.** pathA_36 (dual-engine + tri-reviewed) got **2 transverse photons at `c_γ²=μ_R/ρ_br`** out of the brane's MacCullagh shear
— the first time the one medium provably supports *light* alongside gravity. The residual longitudinal "third wave" was NOT
resolved by a postulated gauge field (that route is a provable no-go, tied to pathA_25) — instead it is **reframed** as the
gravity/electric channel, with the live front being **C5-with-a-throat** (does the electric field + Gauss's law EMERGE from the
puncture–brane coupling). Full detail in the ⭐ v6 block in §3.

**Status update:** v7, 2026-07-03. **v7 = the EM-charge MECHANISM pivot — charge is a 4D throat-body interaction, NOT a brane
flow/deformation; and gravity is a strictly one-way drain.** The v6 "third wave" reframe tried to make the electric field a
*longitudinal brane deformation/flow* around a puncture. Three sub-routes were explored and set aside (each first-class): elastic
**deformation** → wrong energy scaling (gradient energy `(∂u)²∼1/r⁶ ≠` Coulomb `∼1/r⁴`); a one-fluid **flow** → charge=mass by
construction (`ρv` *is* the mass current); a two-stream **counterflow** → needs a second nondissipative component a `T=0`
single-condensate lacks in flat 3D. The resolution (user): **charge is not a flow — it is the interaction of the throats' 4D bodies
beyond the mouth** (the throat's body has more area than its 3D mouth and extends into the bulk). This yields a clean four-way
taxonomy and **re-homes charge into the throat/bulk sector**. Correction folded in: **gravity is strictly one-way** — medium that
enters a throat de-structures into a bulk state and **cannot return to the brane through that throat**; the return is the separate
global `S_leakage`. Full detail in the ⭐ v7 block in §3 (+ §4). The v6 *achievement* (2 transverse photons) stands; only the v6
*mechanism-guess* (longitudinal deformation = the electric field) is superseded.

*(The v8 block, 2026-07-04, is the current magnetism picture and is retained in `conceptual_foundation.md` §3, not here.)*

---

## 2026-06-23 — v3 brane pivot: the little-arrows domain wall FALSIFIED; smectic refinement

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

### The original little-arrows brane hypothesis (the conceptual ancestor of §2, superseded by the v5 material-state closure)

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

---

## 2026-06-24 — our OWN prior internal work on the throat (survey)

*(This was §6.1 of the foundation doc — a map of prior internal throat work so it is not rediscovered.)*

### Our OWN prior internal work on the throat — already explored, never wired together (don't rediscover it again)

The throat-as-self-bound-particle picture has been **partially built twice in this repo** and then fragmented. Read these before
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

## 2026-06-25 — pathA_29 milestone: brane↔bulk return (RETURN_RESIDUAL_PREDICTION)

*(This was a milestone block inside §5 Cosmology.)*

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

## 2026-06-26 — the moving-throat PDE push (Gates 1–5)

*(This was a milestone block inside §5 Cosmology.)*

> **⭐ MILESTONE (2026-06-26) — the MOVING-THROAT PDE push is underway, and the first four gates hold.** Completing the actual
> moving-throat PDE (the distributed-geometry field theory whose solved branch must *return* the GR-quadrupole and the EM sector) is a
> **~6-gate ladder**, not one step — master checklist: `research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md`.
> What matters conceptually: we are checking that the **distributed geometry lift** (promote the throat from collective coordinates
> `(a,L)` to a level-set field `Σ=r−R(Ω,w,t)`) genuinely *contains the structures we already trust* as its limits.
> - **Gate 1 (`pathA_30`) ✅** — the frozen-wall reduction reproduces the exact finite-throat acoustic resonator (mouth DtN
>   `−(ω/c_s)tan(ωL/c_s)`, half-shifted mode ladder). The geometry lift's plumbing is sound. (`DN_UNITTEST_BC_DEPENDENT`: the
>   Dirichlet/Neumann boundary assignment is a labeled calibration input, not yet derived.)
> - **Gate 2 (`pathA_31`) ✅** — switching on the breathing mode, the distributed lift's Hellmann–Feynman force **reproduces the legacy
>   `(a,L)` collective dynamics** (mass + stiffness from one operator projection; the 2-mode truncation is clean for any order-unity
>   wall stiffness). So the old finite-dimensional throat dynamics are genuinely the lowest-mode limit of the field theory.
>   (`BREATHING_CALIBRATED`: the wall material constants are calibration inputs.)
> - **Gate 3 (`pathA_32`) ✅** — switching on the quadrupole channels, the **ℓ=2 quadrupole sector first appears**, and the three
>   directional lanes that the throat can wobble along **collapse to one common set of coefficients on an isotropic (round) reference**.
>   That collapse is exactly what spatial isotropy *means* — there is no preferred direction baked into the field theory; the reduction
>   is **SO(3)-covariant** (rotation-symmetric). And it is genuinely able-to-fail: deliberately squashing the reference throat splits the
>   lanes back apart. So the quadrupole — the lowest radiating gravitational shape — enters the PDE the way it should.
>   (`ISOTROPY_CALIBRATED`: the wall material constants are calibration inputs.)
> - **Gate 4 (`pathA_33`) ✅** — the quadrupole **radiation / normalization** gate. Here the model draws an honest line between what it
>   *earns* and what it *calibrates*. The **shape** of the radiating quadrupole's response (its low-frequency fingerprint) is **earned** —
>   it is **derived from the genuine outgoing-wave boundary condition** (a wave leaving to infinity), not put in by hand. The overall
>   **magnitude** — the famous `54/5` factor and the gravitational coupling `G` itself — is honest **calibration**, not a derivation of
>   Newton's constant. The split is literally `54/5 = 2·27/5`: the `27` is earned from the fingerprint, the `2/5` and `G` are the
>   calibrated bridge. So the PDE delivers the *form* of gravitational-wave emission while we still tune its strength to nature.
>   (`QUAD_CALIBRATED`.)
>
> Honest status: all four are **necessary plumbing**, all land **CALIBRATED** (the toy-model contract — calibrate freely, derive later
 if we can). **Gate 5 (`pathA_34`) is now DONE = `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (EARNED)** — the honest landing: at linear
> order the ℓ=0/1↔ℓ=2 cross-link runs through the return admittances `Z0_ret,Z1_ret`, which pathA_29 shows are a *premise*, so the
> linear theory cannot pin them and **Gate 6 must supply that selector** (the first concrete Gate-6 input). The remaining decisive
> tests are **Gate 6** (the full *nonlinear* closure — the likely wall; its 4D solve is **sim-deferred**) with the **PN-ladder
> match-back** as the falsifier. **⭐ Since the full simulation is out of reach (`docs/development_plan.md`), the program now pivots to
> completing the OTHER sectors to "simulation-ready" — next front = the BRANE sector (§2: the crux that carries light), then
> EM/light, then throat.** Process note (hard-won): Gates 2, 3, and 4 all first runs *looked* clean, but the adversarial review caught
> pass-by-construction defects fidelity missed — Gate-3's first attempt even faked one of its two independent engines (a vacuous "x minus
> x" comparison), and Gate-4's first attempt had a units check that cheated by back-solving a free knob to force the answer. The
> dimensional-rigor win is worth flagging: a units check that can **actually fail** caught a quietly dropped scale factor in the handoff —
> exactly the kind of natural-units bookkeeping slip that hides order-of-magnitude errors. So every verdict-bearing control must compute
> from inputs and carry an able-to-fail probe. The model is holding together so far, which is the only reason to keep pushing.

---

## 2026-07-03 — v5 brane pivot: the material-state closure (density-smectic CLOSED, shear-surface NOGO)

*(This was the ⭐ v5 UPDATE block at the top of §2. It records the current brane candidate — the material-state closure — as
well as the two brane-route deaths that motivated it. The one-line current statement is retained in the foundation doc's
FOUR-SECTOR CHAIN spine; the full mechanism + honest costs + the make-or-break ladder are here.)*

> **⭐ v5 UPDATE (2026-07-03) — second brane pivot: the material-state closure. Read this before the v3 smectic block above.**
> Since v3, **two more brane routes have died** (both first-class falsifications, both kept):
> - The v3 **smectic density-layering** route is **CLOSED** — `RC_DENSITY_SMECTIC_LIGHT_NOGO`: a density modulation *can* form a
>   codim-1 brane, but the same coupling that opens the codim-1 window **pins the arrows along the layer normal** (`P_∥=0`) → it
>   **starves the in-plane shear that carries light**. Density-brane and light-brane are mutually exclusive.
> - The fallback **postulated shear-surface** brane (pathA_35: concede the axis, *postulate* a light-confining shear surface, test
>   whether light lives on it) hit **`FAIL_COUPLE_STRESS_NOGO`** (Gate L, provisional): the frozen MacCullagh light package has **no
>   scalar potential `φ`** to remove its longitudinal zero mode (the "C5" obstruction — below), and the escapes either re-admit
>   hidden `P` modes or collide with the required `u_w` gap.
>
> **The new candidate — a material-state closure (from `notes/brane_bulk_handoff.md`).** Stop asking "what *shape* localizes the
> brane" and instead treat brane-vs-bulk as **two phases of the one medium**:
> - the **brane** is the medium's **ordered, shear-supporting phase** (the phase in which the substructure locks into a state that
>   can carry MacCullagh shear = light);
> - the **bulk** is the **same medium de-structured** — a shear-free phase; this is now *why* the bulk carries no light and preserves
>   magnetism (§3), a **consequence of one order parameter** instead of a separate postulate;
> - a **throat** is an open defect where **ordered brane material de-structures into the bulk phase** (the steady drain = gravity);
>   **return** is bulk material **re-ordering** into the brane phase. Brane↔bulk becomes a **conserved phase cycle**, which **kills
>   the old "bulk vacuum sucks brane material in" paradox** (an empty sink could never explain re-entry). It formalizes the exact
>   "draining = de-structuring" language already in §5.
> - Introduce a brane-order field `χ_B ∈ [0,1]` (name provisional): `χ_B=1` ordered/shear-supporting, `χ_B=0` de-structured/bulk.
>   Shear/light exists only where `χ_B≈1`.
>
> **Why this is more than a re-description — it may escape BOTH prior brane deaths (to be tested).** A `χ_B` domain wall (ordered
> slab at `w≈0`, de-structured outside) is a **genuinely different object** from both corpses:
> - It **escapes T1** (the little-arrows wall that *spread and unwound* on a **connected** `S³` vacuum): a scalar `χ_B` with a real
>   **double-well** free energy `f_B(n,χ_B)` has **disconnected** minima → a **φ⁴-kink-like, topologically stable wall** — exactly
>   the structure T1's own negative control (the stable φ⁴ kink) showed it lacked.
> - It **escapes the density no-go**: `χ_B` is a **single domain wall in an abstract order field, not a periodic density stack**, so
>   it has **no layer normal to pin the arrows against** — the in-plane shear lives in a separate field that `χ_B` merely *gates on*.
>
> **The exciting sub-hypothesis — `χ_B`'s phase may supply the missing C5 `φ`.** The shear-surface death was precise: MacCullagh's
> curl-only energy `½μ_R(∇×u)²` is gauge-invariant under `u→u+∇χ` but its kinetic term is not, forcing a **constrained physical
> longitudinal zero mode** (`∂_t²(∇·u)=0`); Maxwell escapes via a scalar potential (`φ→φ−∂_tχ`) but MacCullagh has none, and the
> only on-brane scalar (`u_w`) **must stay gapped** (a massless out-of-plane mode = an excluded fifth force), so it *cannot* be `φ`.
> But if the brane-order parameter is **complex** (amplitude `χ_B` + **phase `θ`**), that phase is a **new on-brane scalar with
> genuine mechanical provenance that is NOT `u_w`** — a natural candidate for the exact `φ`-analog the frozen light package lacked.
> If `θ` couples to the displacement like Maxwell's scalar potential, it removes the longitudinal zero mode and **could revive
> shear-surface light** as a *fresh G0*. Whether it actually does is the test — **not a claim**.
>
> **Honest costs (do not fool ourselves).** The double-well `f_B` is **still postulated** — the GNLS potential `U(ρ)∝ρ⁵` is
> single-well, so "two coexisting phases" is the *same* degenerate-vacua obstacle the brane has always had (§ "The obstacle" below),
> now conceded under the analog license (§0.6), **not derived**. It is **real drift**: `f_B` + the interface stiffness `κ_B` + the
> `χ_B`-shear gating + (if complex) the phase sector — the gauntlet must count it and decide whether `χ_B` is
> substructure-of-one-medium or a smuggled second field. And it may simply **relocate** the problems (the
> wall-exists-*and*-carries-shear rung is where I expect it to bite).
>
> **The make-or-break ladder (each rung able-to-fail), tested next (§8):** **W** — a double-well `χ_B` wall that is stable (escapes
> T1) *and* light-permitting (escapes the density no-go); **φ** — does `χ_B`'s phase supply a legitimate C5 `φ` (removes the
> longitudinal mode with no fifth force and no `u_w` collision)?; **Gate L** — re-run the six light sub-hurdles on the
> `χ_B`-founded G0. Strictly W→φ→Gate-L in order; a failure at W kills the route cheaply.

*(Program note: since v5, rung W has since **passed** — the FOUR-SECTOR CHAIN spine in the foundation doc now records `χ_B` as a
"genuine self-localizing Z₂ material-state wall, rung-W-passed: stable + light-permitting.")*

---

## 2026-07-03 — v6: the medium demonstrably carries light; the C5 "third wave" (mechanism-guess later superseded by v7)

*(This was the ⭐ v6 UPDATE block at the top of §3. The *achievement* — 2 transverse photons, pathA_36 — still stands and is
carried in the foundation doc's spine + v8 block. The v6 *mechanism-guess* (the longitudinal deformation = the electric field)
was superseded by the v7 charge-as-4D-throat-body mechanism, which is retained live in the foundation doc.)*

> **⭐ v6 UPDATE (2026-07-03) — the medium DEMONSTRABLY carries light; and the "third wave / C5" is where the sectors CONNECT.
> Read this — it reframes the light sector and defines the live front.**
>
> **Achievement (pathA_36, dual-engine + tri-reviewed).** The frozen shear-surface light package (MacCullagh curl shear on the
> brane + brane inertia) yields **exactly 2 transverse photon polarizations at `c_γ²=μ_R/ρ_br`** — the FIRST time in the program the
> one medium demonstrably supports *light* alongside gravity. That is real and banked.
>
> **The C5 "third wave" (the honest obstruction).** MacCullagh (curl-only) stiffness leaves the **longitudinal** (compression, `∇·u`)
> brane motion with no shear-restoring force — a "third mode" that ordinary elasticity would turn into a spurious 3rd photon. The
> attempt to make the **brane-order phase `θ`** act as Maxwell's scalar potential `φ` and remove it **fails in VACUUM**
> (`FAIL_CAUCHY_STRAY_LONGITUDINAL` / `BY_TUNING`, EARNED — arbiter reproduced, Codex + adversarial `EARNED` from-scratch): a *stable*
> order-parameter phase has **positive** gradient-stiffness (that positivity is *what makes the brane stable*), which is the **wrong
> sign** for Maxwell's `φ`; the only escape (negative stiffness) is a spatial-modulation instability = the **same obstruction that
> killed the density-smectic route (pathA_25)**. So the order-parameter phase **provably cannot self-supply the gauge structure**, and
> the two no-gos are *one* wall.
>
> **THE REFRAME (the live direction — do NOT bolt on an abstract gauge field).** The vacuum no-go answers the *wrong* question. In
> real EM the longitudinal/scalar sector is **pure gauge (nothing) in vacuum** and only becomes physical — the **Coulomb field** —
> in the presence of **charge**. pathA_36 analyzed the shear sector in VACUUM (no throats) — i.e. it diagnosed the one mode whose
> entire character comes from charge, *with no charge present*. Natively, the longitudinal (compression) brane mode is **NOT a
> spurious photon**: it is (a) the medium's **density/compression = the `c_s` gravity sector**, and (b) its static-around-a-charge
> part = the **ELECTRIC FIELD = the brane's tension-deformation around a throat** (§4 already: "electric field energy = brane-tension
> deformation around the puncture"). So the "third wave" is exactly **where light (shear) meets gravity (density) and charge
> (throat).**
>
> **The missing ingredient is CHARGE, not a gauge field.** Maxwell's gauge structure — Gauss's law, exactly 2 propagating photons —
> should **EMERGE** from the **throat–brane coupling**: the longitudinal brane deformation around a puncture *is* the electric
> (Coulomb) field; Gauss's law = "brane tension sourced by punctures." And the longitudinal/transverse split is naturally **two
> sectors at two speeds** (`c_s` compression/gravity vs `c_γ` shear/light), not "3 modes of one field" (pathA_36 used a single
> inertia, which blurred this). If the longitudinal mode couples to charge as *radiation* (not just as the static Coulomb field),
> that is a **falsifiable departure** from pure Maxwell (cf. the pathA_29 monopole wave) — possibly a feature.
>
> **Why this matters (the potential massive win).** If the electric field / Gauss constraint **emerges from the puncture–brane
> coupling**, the forces are not merely *coexisting* in one medium — they are **CONNECTED**: light = shear, gravity =
> density/compression, electric field = static brane-deformation around charge — all the *same* longitudinal/transverse structure of
> the *one* medium. The "third-photon problem" dissolves into "the sectors are linked."
>
> **NEXT (post-compact) = the C5-WITH-A-THROAT gate.** Put a puncture (charge) into the flat-brane light analysis and test whether
> (i) the longitudinal brane response to the throat is the **static electric (Coulomb) field** (a constraint, not a free wave);
> (ii) **Gauss's law emerges** from "brane deformation sourced by the puncture"; (iii) the longitudinal (compression) sector is the
> **`c_s` gravity/density** channel, distinct from `c_γ` light; (iv) whether any longitudinal **radiation** coupling to charge appears
> (a falsifiable departure). Build via the pipeline (directive → Codex design-review → dual-engine → tri-review). This **supersedes**
> the pathA_36 vacuum-gauge remediation (the vacuum result stands as a documented waypoint, not the final word). Detail +
> deliverables: [[project-brane-existence-defect-structure]]; directive `software/stage1_solver/directives/pathA_36_c5_phase_potential.md`.

---

## The consistency-gate program for the GNLS polar-smectic superfluid (post-pivot plan, largely superseded)

*(This was §8 of the foundation doc — the forward plan / consistency-gate ladder. It is dominated by the GNLS-polar-smectic and
material-state-closure programs and is now largely superseded by the completed four-sector chain; the current "what's next"
lives in the foundation doc's FOUR-SECTOR CHAIN spine ("THE CLOSING STEP — the CONSISTENCY KNIT"). Retained here for the record.)*

> **⭐ v5 UPDATE (2026-07-03) — this GNLS-polar-smectic gate program is itself now SUPERSEDED at the brane-existence level** (its
> density route is CLOSED and the postulated shear-surface fallback hit `FAIL_COUPLE_STRESS_NOGO`; see the §2 v5 block). The live
> plan is the **material-state closure** and its **W→φ→Gate-L ladder**: **(W)** does a double-well `χ_B` wall localize the brane
> while staying light-permitting; **(φ)** does `χ_B`'s phase supply the C5 scalar the shear-surface freeze lacked; **(Gate L)**
> re-run the light sub-hurdles on the `χ_B`-founded G0. Same pipeline as below (Codex design-review xhigh → GLM tertiary → execute,
> tri-reviewed, user-gated). **Immediate step:** a **GLM reframe pre-checks rung W** — is the double-well `χ_B` wall stable *and*
> light-permitting, or is there a fresh no-go between codim-1 localization and in-plane shear? — before the expensive directive+G0
> build. The gate list below (S/B/Q/K/T) is retained as the **downstream requirements** any surviving light-capable structure must
> still satisfy.

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

## 2026-07-05 — pathA_42: the charge-coupled scalar MAPPED (SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED)

The v8/pathA_39 results established that the assembled EM sector is **not exact Maxwell** — it carries an extra **charge-coupled
scalar** (the propagating `h`-branon). pathA_42 (dual-engine SymPy + Mathematica, tri-reviewed) gave that scalar its **physical
meaning** and mapped its character:

- **What the scalar IS:** the **electric field's flex INTO `w`** — a smooth **bulge** of the brane into the hidden dimension around
  a charge (a swelling of the embedding), **not** a throat / hole / puncture. It is the same embedding-direction (`h`) family as the
  Gate-L bending scalar, but a **distinct ledger object**.
- **It radiates:** when a charge **accelerates**, the scalar carries a dipole radiation channel alongside the ordinary EM (Poynting)
  one — derived scalar far-zone flux `∝ d²ω⁴q_h²/(M_h c_E³)` vs EM `∝ Q_E²d²ω⁴/c_γ³`. On the current ledger branch this sits at
  radiation exponent 3 (the exponent-1 branch is conditional on a future pinned-`K_h` fact).
- **Verdict `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED`:** the *character* of the departure is now mapped (a real, falsifiable
  5th-force-like scalar), but its **magnitude is SIM-GATED** — the free `h`-sector kinetic normalization, `q_h/Q_E` universality,
  `c_E=c_γ` (preferred-frame), and the `u_L`/`C_hu` mixing that could transfer mass-coupling into charge-sourced modes are all
  **not yet pinned**. So whether the scalar is a **falsifiable observable break** or is **naturally hidden** remains open (a
  `HARD_WALL`; Guard A refused all laundering fixtures — no numeric magnitude was emitted).
- **EP note:** the mass channel `h_EP` is `EARNED_SAFE` only on the *decoupled floor*; nonzero `C_hu` mixing (`C_hu·q_h·q_M`)
  reintroduces a mass coupling → `MIXED_SCALAR_EP_RISK`, so full EP safety additionally requires universality to be earned.

Detail: `software/stage1_solver/{directives/pathA_42_charge_coupled_scalar.md, reports/pathA_42_charge_coupled_scalar.md,
reports/pathA_42_charge_coupled_scalar_results.yaml, tools/pathA_42_charge_coupled_scalar.*}`.

---

## Cross-references
- The vision this history records: `docs/conceptual_foundation.md` (read first).
- EM re-founding detail + the MacCullagh/shear math history: `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md`.
- The directive reworked into the T1–T5 ladder: `software/stage1_solver/directives/pathA_24_brane_existence_defect_structure.md`.
- Requirements spec + prior-art map: `docs/medium_requirements_and_prior_art.md`.
- Program state / "you are here": `STATUS.md`; `software/stage1_solver/decisions/13_emergent_constants_derivation.md` §0.
