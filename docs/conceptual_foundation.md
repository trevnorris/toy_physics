# Conceptual Foundation — what each part of the model physically IS

**Read this FIRST.** This is the canonical plain-language statement of the model's physical vision — what the medium, the
brane, the four force-sectors, and the particle/defect actually *are*, in the model's **own native terms**. It exists so the
conceptual picture is never lost between work sessions and never has to be re-explained from scratch. It is a **living document**:
update it whenever the conceptual picture sharpens. It is deliberately separate from the math/derivation/verdict machinery
(those live in `software/stage1_solver/decisions/*`, the `pathA_*` directives, and `STATUS.md`).

**Derivation ladder + dated history: `docs/conceptual_history.md`.** This doc is the *vision*. The version-by-version changelog,
the dated per-gate / per-stage compute records, and the superseded-hypothesis detours (the little-arrows domain wall, the
density-smectic route, the shear-surface no-go, the flow-framing retirement, etc.) have been **moved** to that companion history
doc so this one stays a clean read-first statement of what the model IS. Nothing was deleted — only relocated.

**Status:** the current picture is below — read the ⭐ **FOUR-SECTOR CHAIN** spine, then §0–§5, then the ⭐ **v8** magnetism block
in §3. All four force-sectors (gravity, light, electric charge, magnetism) are now EARNED from one brane+bulk; what remains is the
consistency knit. Version history (v3 little-arrows falsification → v5 material-state closure → v6 light → v7 charge → v8 magnetism
→ pathA_42 scalar): `docs/conceptual_history.md`.

---

## ⭐ THE FOUR-SECTOR CHAIN — how ONE brane+bulk yields gravity + light + electric + magnetism (2026-07-04 status map)

**This is the spine of the whole program: one medium, organized into a brane+bulk, from which the four forces emerge as distinct
things that the SAME structure does. Read this as the derivation chain + where each link stands.**

**THE FOUNDATION (the brane+bulk) — now DEFINED (was imposed a year ago).** One compressible superfluid in 4+1D, in two phases: an
**ordered, shear-supporting phase = the brane** (our 3D space, a codim-1 wall; order field `χ_B` — a genuine self-localizing Z₂
material-state wall, rung-W-passed: stable + light-permitting) and a **de-structured, shear-free phase = the bulk** (`w≠0`). **Throats**
are open punctures where brane-ordered medium de-structures into the bulk (phase-conversion sites). Everything below hangs off this one
structure.

1. **GRAVITY = the throat drain (one-way inflow).** A throat swallows medium (de-structures it into the bulk; no return through the
   throat — return = the separate global `S_leakage`). The inflow `v_r` IS the field; a *change* in it propagates at `c_s`. Localizes to
   **`1/r²`** via a normalizable transverse zero mode through the finite slab (RS-like). Parity: the **EVEN** channel. **STATUS: EARNED**
   — PN ladder 1PN→4PN + 2.5PN GR-matched (`research/4d_*pn*`) + `pathA_29` localization (`RETURN_RESIDUAL_PREDICTION`).
2. **LIGHT = the brane's in-plane shear wave** (MacCullagh rotational-elastic). 2 transverse polarizations at `c_γ²=μ_R/ρ_br`; confined
   to the brane because the bulk is shear-free. **STATUS: EARNED** — `pathA_36` (dual-engine + tri-reviewed): the medium demonstrably
   carries 2 transverse photons (first time).
3. **ELECTRIC CHARGE = the interaction of two throats' 4D BODIES beyond the mouth.** Sign = `±w` puncture direction. Mediator = the
   **gapless transverse-embedding/orientation-lock Goldstone `h`** (wall displacement into `w`, arrows locked to the normal). Localizes
   to **`1/r²` Coulomb** by the SAME principle as gravity (normalizable `sech²` zero mode); like-repel/unlike-attract from `G₀>0`.
   Parity: the **ODD** channel — which is *exactly why charge ≠ mass* (`⟨f₀,S_grav_even⟩=0` by `w→−w`). **STATUS: EARNED** — `pathA_38`
   = `THROAT_ELECTRIC_LOCALIZED_COULOMB` (dual-engine + tri-reviewed). DERIVED{p=2, sign, parity} + CALIBRATED{`Q_E`, that a throat
   sources a nonzero monopole} + SIM-DEFERRED{operator-level parity mixing, source compactness, nonlinear interior}.
4. **MAGNETISM = the velocity-dependent (`O(V)`) part of the SAME 4D throat-body interaction** — the moving/rotating version of the
   electric charge of #3 (charge-coupled by construction; `±w` sets its sign). ⭐ **Two swirls must NOT be conflated:**
   **gravitomagnetism** = the 3D flow-swirl of the mass-inflow (EVEN, mass-coupled, `c_s`, weak — frame-dragging / Gravity Probe B;
   part of the *gravity* sector, one flow field) vs **EM magnetism** = the 4D throat-body swirl (ODD/`±w`, charge-coupled, strong,
   because the swirl concentrates toward the finite mouth radius). **STATUS: ⭐⭐ EARNED — SECTOR COMPLETE (`pathA_39`, all 4 stages;
   this is the FOURTH and last force-sector, so ALL FOUR now come from ONE brane+bulk).** The four computed results, all dual-engine
   + tri-reviewed CLEAN: (0+1) **exact isolated Maxwell is EXCLUDED** — the EM *field* carries a charge-coupled scalar admixture (a
   propagating `h`-branon, forced by pathA_38's `M_h>0`; a 5th-force-like departure); (2) a charge-coupled `R^-2` current-current
   magnetic force, like currents attracting (transverse) + an unavoidable attractive longitudinal scalar-current admixture; (3)
   `FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING` — because charge is `P_w`-ODD (`±w`), a MOVING charge is symmetry-allowed to mix the
   EM-magnetic and gravitomagnetic sectors at `O(V)`: **the two swirls above are ONE 4D-throat swirl, not operator-level separable
   under motion** (magnitude sim-deferred); (4) `FIELD_SCALAR_VECTOR_DEPARTURE` — assembling the whole linearized EM sector as one
   quadratic action gives a **transverse-vector + propagating charge-coupled-scalar (`h`-branon) multiplet, NOT exact Maxwell**. Net:
   magnetism is a real fourth force + a *characterized departure* from clean Maxwell — first-class for a toy analog. Detail: the ⭐ v8
   block below + `directives/pathA_39_magnetism_close_maxwell.md`, `..._stage3_operator_parity.md`, `..._stage4_field_classification.md`.

**THE CLOSING STEP (after magnetism) — the CONSISTENCY KNIT → "fully formed."** Four sectors each passing in isolation is NOT the same
as one self-consistent model. The knit: (a) the sharp cross-sector test **`λγ = c_γ/c_s = 1`** (the GW170817 cone-lock — do light and
gravity travel at the same speed in ONE parameter set? — still OPEN); (b) the NG cross-consistency gauntlet (all sector couplings
mutually consistent, no contradiction); (c) **assemble the end-to-end chain into the central `pde_ledger`** (the calibrated PDE
delivering GR + EM). Completing this turns *"the four forces each work"* into *"one fully-formed brane+bulk model."* Honest boundary
([[project-simulation-deferred-complete-pde-strategy]]): this completes the **spec**; the full nonlinear sim stays deferred, so
sim-dependent questions stay posed and a no-go is still possible.

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

> **Current brane picture + brane-pivot history:** the current candidate is the **material-state closure** — the brane is the one
> medium's *ordered, shear-supporting phase* (order field `χ_B`), the bulk is the *same medium de-structured / shear-free*, and a
> throat is a phase-conversion site (see the ⭐ FOUR-SECTOR CHAIN spine above). The full derivation of how we got here — the v3
> little-arrows-wall falsification, the density-smectic close, the shear-surface no-go, and the material-state closure with its
> W→φ→Gate-L ladder — is in `docs/conceptual_history.md`. The ontology below (what a brane IS, bulk on both sides, the single-well
> obstacle, and the keystone question) is timeless and stays here.

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

> **⭐ Slab caveat (2026-07-07, from the `ledger_stage006` build — do not let "wall admitted" read as "slab stable").** With the
> material-state closure's `χ_B ∈ [0,1]` double-well (minima at 0 and 1), "bulk on both sides" makes the brane a finite
> `χ_B = 1` **slab** bounded by `χ_B = 0` bulk above and below — a **kink–antikink pair**, not a single wall. The double-well
> alone provides **no mechanism selecting the slab width** (kink–antikink pairs generically attract); the width `W_slab` is an
> un-earned input — the old `Z/W/B_ℓ` profile scale in new clothes — and maps onto the known `L/a` self-selection open item
> (§7 #11, "requires dynamics", sim-deferred). `ledger_stage006` verifies single-kink *admission* (EL residual + surface tension
> `σ_wall = √(2a_Bκ_B)/6`) and registers `W_slab` as `FREE-UNREDUCED` in the parameter register.

**The obstacle (honest).** The current medium potential `U(ρ) = (K/4)ρ⁵` is a **single well** — one stable state — so as written
the medium *cannot* form a wall. That is exactly why the brane currently has to be *put in by hand* (a confinement potential
`V_conf`, the `w`-profiles `Z/W/B_ℓ`, the `k_w u_w²` restoring term — all *inputs*, not derived). To **derive** the brane, the one
medium needs structure that gives it **two coexisting states.**

*(The original **little-arrows** working hypothesis for the wall's internal structure — polar constituents defining `w`, why-space-
is-3D, the two `±w` domains, charge falling out, and the alternatives considered — is the superseded conceptual ancestor of the
material-state closure. It is preserved in full in `docs/conceptual_history.md`.)*

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

> **Light sector — banked, with the C5 "third wave" history moved out.** pathA_36 (dual-engine + tri-reviewed) got **exactly 2
> transverse photon polarizations at `c_γ²=μ_R/ρ_br`** out of the brane's MacCullagh shear — the FIRST time the one medium
> demonstrably carries *light* alongside gravity (banked). The v6 "third wave / C5" reframe that followed — the longitudinal
> (compression) mode is the `c_s` gravity/density channel, its static-around-a-charge part is the electric field, and the
> C5-with-a-throat plan to make Gauss's law EMERGE from the puncture — has since been executed and is subsumed by the v7 charge
> mechanism (below) and the ⭐ v8 magnetism block. That dated waypoint is preserved in `docs/conceptual_history.md`.

> **⭐ v7 UPDATE (2026-07-03) — the EM-charge MECHANISM: charge is a 4D THROAT-BODY interaction, not a brane flow. This supersedes
> the v6 mechanism-guess (the v6 *achievement* — 2 photons — still stands). Read this: it is the current picture.**
>
> **The four-way taxonomy (each force = a distinct thing the ONE medium does — keep them clean):**
> - **Gravity** = a throat *swallowing* medium — a **one-way inflow/drain.** Medium that passes the mouth **transforms
>   (de-structures)** into a bulk state and **cannot return to the brane through that throat**; the return is the separate **global
>   `S_leakage`** process, not a local throat loop. (This retires any "recirculation down-and-back-up the throat" picture.)
> - **Light** = shear/ripple **along** the brane (2 transverse polarizations, `c_γ`). *Banked* (pathA_36).
> - **Magnetism** = **swirl / vorticity of** the brane (a moving charge → swirl).
> - **Charge** = the throats' **4D bodies interacting beyond the mouth** — a *geometric* interaction of the extended throat structures
>   in the bulk. **NOT a brane flow, not a deformation, not a swirl.**
>
> **Why charge is NOT a flow (the honest evolution — three sub-routes, each set aside as first-class):** (a) electric field = elastic
> brane **deformation** → wrong energy scaling (generic first-gradient elasticity forces `u_L∼1/r²`, so gradient energy
> `(∂u)²∼1/r⁶`, but Coulomb needs `∼1/r⁴`; self-energy `∼1/a³` not `∼1/a`); (b) electric field = a one-fluid **flow** → charge=mass
> by construction (`ρv` *is* the mass current — a flow with energy moves mass = gravity); (c) a two-stream **counterflow**
> (mass-neutral) → needs a *second nondissipative component* the `T=0` single-condensate medium lacks in flat 3D (one U(1) →
> `ker(J_m)={0}`). The way out is **not another flow** — it is that charge lives in the throat's **4D body**, where flat-3D
> reasoning does not apply.
>
> **Charge specifics.** `±w` (throat puncture direction) = the **sign**. **Like-throats resist 4D overlap → repel; opposite (`+w` &
> `−w`) merge → attract, and can annihilate** (the two throats join and the brane heals = pair annihilation; pair *production* =
> nucleating a `+w`/`−w` pair, so **net charge is conserved** while the *count* of charges can grow). Charge magnitude is
> **universal** (a topological property of the puncture — cannot drift with energy/mass). The **classical electron radius `r_e`** is
> the **throat-body size**: setting the electric self-energy `∼ke²/a` equal to the rest energy `m_ec²` (the trapped-geon energy that
> holds the throat open) gives `a∼ke²/(m_ec²)=r_e`. "You cannot push two electrons closer than `r_e`" = "you cannot force two 4D
> throat-bodies to overlap" — the finite-size resolution of the point-charge self-energy divergence.
>
> **⭐ THE MAKE-OR-BREAK (must be nailed in the derivations): does it give `1/r²`?** An electric field is **long-range** (`1/r²` force
> between distant charges), so this cannot be only a short-range *contact* collision — the throat's 4D body must mediate a
> **long-range** interaction that comes out `1/r²` **in the 3D brane.** The catch: a source acting through the **4D bulk** naively
> gives the **wrong** falloff (4-space → `1/R²` *potential* → `1/R³` force). It becomes correct 3D Coulomb only if the throat's
> influence is **brane-localized** (does not spread freely into the bulk). Encouragingly, **that is the same localization mechanism
> that already made *gravity* come out `1/r²` through the finite slab** (pathA_29's localizing flat-slab family) — so there is a
> demonstrated route, but it **must be shown for charge, not assumed.** This is the crux of the charge derivation.
>
> **Future item (deferred — do NOT calculate yet, but do not lose):** the same `1/r²`-vs-localization mechanism may, at **large
> separations (galaxy scale, widely-separated bodies)**, turn `1/r²` into `1/r` — a candidate native handle on **galactic rotation
> curves / a MOND-like modification.** Parked for later.
>
> **Program implication.** This **re-homes the charge crux into the throat/bulk sector** (the 4D throat interior — Gate-T territory,
> where mass = the geon, throat size = `r_e`, and annihilation already live), and **out of the flat-brane far-field.** The
> flow-based flat-brane charge gate (`pathA_37`, in its flow/counterflow form) is **retired** — it tested the wrong mechanism. The
> charge crux is now: **does a throat's 4D body mediate a brane-localized interaction that comes out `1/r²` in 3D, with `±w` giving
> the sign?** Detail: [[project-brane-existence-defect-structure]], [[feedback-native-em-mechanisms]].

> **⭐ v8 UPDATE (2026-07-04) — the MAGNETISM mechanism refined + the FIRST computed result (`pathA_39`). Magnetism is the
> moving/rotating version of the v7 charge mechanism, and its first gate stage is banked. Read this for the current view of
> magnetism (it supersedes the "swirl = bulk vorticity, close Maxwell as a gauge field" framing).**
>
> **Magnetism is the velocity-dependent (`O(V)`) part of the SAME 4D throat-body interaction that gives electric charge** — the
> moving/rotating throat, not a separate mechanism.
>
> **⭐⭐ WHERE MAGNETISM LIVES — IN THE THROAT (do NOT forget this): both electric charge AND magnetism are 4D-THROAT-BODY
> phenomena.** Electric charge = the **static** interaction of two throats' 4D bodies (pathA_38); magnetism = **the SAME throat-body
> interaction under MOTION** (pathA_39). Both happen in the throat's 4D body beyond the mouth and are felt on the brane via
> localization — **NOT** a brane-surface effect like light (in-plane shear), **NOT** a bulk-flow effect like gravity (the drain).
> The `±w` puncture direction sets the sign for both. **Do NOT relapse to "magnetism = bulk vorticity of the medium" — that is
> *gravitomagnetism* (weak, mass-coupled, part of the gravity sector), a DIFFERENT thing.** The mnemonic: gravity = flow, light =
> brane shear, **charge & magnetism = the throat (static vs moving).**
>
> This applies the v7 electric-sector route under motion: charge-coupled by construction, `±w` sets the sign. It is NOT "does the
> curl of the light field close Maxwell as a gauge field" — that framing was tried and is a *predetermined* fail (the shear field's
> longitudinal stiffness `B_eff=ρ_B0²/χ_c>0` from pathA_36 blocks clean gauge closure).
>
> **⭐ TWO distinct "swirls" — do NOT conflate them (the load-bearing clean-up):**
> - **Gravitomagnetism** = the **3D flow-swirl of the mass-inflow** — the medium literally circulating as it drains into a moving or
>   spinning throat. **Mass-coupled, EVEN under `w→−w`, rides `c_s`, WEAK** (ordinary frame-dragging — Gravity Probe B saw only
>   milli-arcseconds/year). It is **part of the gravity sector** (gravity's radial inflow + this swirl = one flow field), NOT a
>   separate force. The bulk is shear-free, so the medium **cannot "push" sideways** — which is *why* frame-dragging is tiny.
> - **EM magnetism** = the **4D throat-body swirl** — the same swirl carried past the mouth into the throat's 4D body, flavored by
>   `±w`. **Charge-coupled, ODD, STRONG.** Felt on the brane by the SAME localization that gives charge its `1/r²` (a brane-trapped
>   field, not a "shadow" reached for in the bulk).
> - **One swirl, two regimes:** 3D near-mouth flow = gravitomagnetism (weak); 4D throat-body = EM magnetism (strong). **Why EM ≫
>   gravity:** the swirl **concentrates** as it funnels toward the mouth (`v_φ ~ 1/r`, saturating at the finite mouth radius `r_e`) —
>   the same localization/concentration that made charge strong. Consistency checks that match reality: magnetic moment ∥ spin ∥
>   angular momentum; the EM sense flips with `±w` while the gravitomagnetic sense follows mass (electron vs positron → same
>   frame-dragging, opposite magnetic moment).
>
> **The longitudinal mode is NOT a spurious "third photon" — it is the gravity/density (`c_s`) channel.** In a ONE-medium model, EM
> is not an isolated sector: the transverse ripple of the shear field is light (`c_γ`); its longitudinal (compression) part is a
> density wave (`c_s`) = the gravity/density channel. "Not exact isolated Maxwell" is EXPECTED and is the *point* of unification — the
> question is whether the coupling stays clean or produces an observable departure. Do not force it to zero; characterize it.
>
> **⭐ FIRST COMPUTED RESULT (`pathA_39` Stage 0+1 — dual-engine SymPy+Mathematica, `ENGINE_AGREE`, tri-reviewed CLEAN: arbiter
> reproduced + `ADVERSARIAL_SOUND` + `FIDELITY_CLEAN`): `FAIL_OBSERVABLE_SCALAR_ADMIXTURE`.** The electric charge has a **nonzero
> residue on a propagating scalar pole** — a **charge-coupled scalar** the medium carries that clean EM does not. **Exact isolated
> Maxwell is EXCLUDED**, computationally. The robust floor is the **propagating `h`-branon** (charge residue `∝ q_h²/M_h`), FORCED by
> pathA_38's own `M_h>0` (its static Coulomb mediator propagates dynamically) — so even the *electric* sector's dynamic completion is
> not clean Maxwell (the static `1/r²` Coulomb + sign + charge≠mass still stand, as the `ω→0` limit). This is the honest
> **"characterized departure"** landing we expected: not a clean Maxwell copy, but a computed, falsifiable 5th-force-like scalar.
>
> **⭐ SECOND COMPUTED RESULT (`pathA_39` Stage 2 — dual-engine SymPy+Mathematica, tri-reviewed CLEAN): `MAGNETIC_FORCE_DERIVED`.**
> Integrating out the declared moving-throat sources — via a genuine propagator inversion + exchange integral, with the sign and
> falloff FALLING OUT (not asserted) — gives a brane-localized current-current kernel with `R^-1` potential and `R^-2` point-force
> falloff. **The transverse channel gives the correct magnetic sign — like currents attract for `μ_R>0`** (matching magnetostatics).
> The longitudinal channel is ALSO attractive for `B_eff>0`, so the old cancellation/crossover story is retired: **magnetism rides
> alongside an UNAVOIDABLE attractive scalar-current admixture** (the same Stage 0+1 departure, uncancelable at the force level).
> Preferred-frame unless `c_E=c_γ`; conditional on the sim-deferred amplitudes `aT,aL`. ⚠️ Its first pass was a pass-by-construction
> RIG (hand-written kernel, cosmetic anti-readback, a WRONG asserted longitudinal sign that faked a `Ξ` crossover) — CAUGHT by the
> adversarial + fidelity legs, genuinely remediated (real integral; the Biot–Savart fixture now correctly `FAIL_TARGET_READBACK`,
> rejected on derivation-process not expression-shape), re-tri-reviewed CLEAN. ⭐ **Stage 3 (operator parity, dual-engine +
> tri-reviewed CLEAN) = `FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING`:** because charge is `P_w`-ODD (`±w`), a MOVING charge's `O(V)`
> operator is symmetry-allowed in the combined-parity-mixing class, so it generically mixes the EM-magnetic (odd) and
> gravitomagnetic (even) sectors — **the two swirls are ONE 4D-throat swirl, not operator-level separable under motion** (via a
> spurion-EFT selection-rule gate, since the naive `δÔ_i≡0` on the reduced slab; magnitude sim-deferred). ⭐ **Stage 4 (field
> classification, dual-engine + tri-reviewed CLEAN) = `FIELD_SCALAR_VECTOR_DEPARTURE` — SECTOR CLOSE:** assembling the full
> linearized EM sector as ONE quadratic action `Q(ω,k)` over `(u_T1,u_T2,u_L,h)` and COMPUTING its DOF (Dirac: 4 real / 2 under a
> Maxwell counterfactual), scalar-sector stability, and per-pole charge residues gives a **transverse-vector + propagating
> charge-coupled-scalar (`h`-branon) multiplet, NOT exact Maxwell**. **The sector landed as a fourth force WITH a characterized
> departure, not a clean Maxwell copy — ALL FOUR force-sectors now earned from ONE brane+bulk; first-class for a unified toy model,
> computed (and adversarially verified) rather than asserted.** Detail: `directives/pathA_39_magnetism_close_maxwell.md` +
> `reports/pathA_39_{scalar_admixture_screen,magnetic_force}.md`.
>
> **⭐ pathA_42 UPDATE (2026-07-05) — the charge-coupled scalar MAPPED.** The extra charge-coupled scalar the EM field carries
> (the propagating `h`-branon of the results above) now has a physical meaning: it is the **electric field's flex INTO `w`** — a
> smooth **bulge** of the brane into the hidden dimension around a charge, **not** a throat/hole/puncture. It **radiates when a
> charge accelerates** (a scalar dipole channel alongside the ordinary EM one). Verdict `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED`
> (dual-engine + tri-reviewed): the departure's *character* is mapped (a real, falsifiable 5th-force-like scalar), but its
> **magnitude is SIM-GATED** — the free `h`-sector kinetic normalization, `q_h/Q_E` universality, `c_E=c_γ` (preferred-frame), and
> the `u_L`/`C_hu` mixing that decides whether it is a falsifiable observable break vs naturally hidden are not yet pinned. Detail:
> `directives/pathA_42_charge_coupled_scalar.md`, `reports/pathA_42_charge_coupled_scalar.md`.

> **Downstream payoff (parked — do not compute yet):** the 4D throat/vortex structure gives a candidate **resistivity-free magnetic
> reconnection** — 1D vortex lines are codim-3 in 4D, so they cannot link/knot and can reconnect by a smooth excursion through the
> bulk (no dissipative cut). This is where the parked MHD work ([[project-mhd-reconnection-parked]]) reconnects once the EM sector is
> settled.

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
| **Magnetism** | Two swirls, kept distinct (⭐ v8 block above): **EM magnetism** = the velocity-dependent (`O(V)`) part of the 4D throat-body interaction (the *moving* version of electric charge — charge-coupled, `±w` sign, strong) vs **gravitomagnetism** = the 3D flow-swirl of the mass-inflow (mass-coupled, weak frame-dragging; part of the *gravity* sector). The bulk is shear-free (it cannot "push"), which is *why* frame-dragging is tiny and light stays confined to the brane. | throat-body swirl (`±w`) / flow vorticity |
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

> **⭐ v7 (2026-07-03) — what the electric FORCE is (not just the charge label).** The bullets above define charge's *sign* (`±w`) and
> *magnitude* (universal, topological). The **electric force/field** is the **interaction of two throats' 4D bodies beyond the mouth**
> (§3 ⭐ v7) — **not** a brane flow or deformation (both were explored and set aside: deformation fails on energy scaling, a flow is
> charge=mass, a counterflow needs a second component the medium lacks in flat 3D). `r_e` (classical electron radius) = the throat-body
> size, from the self-energy balance `ke²/a = m_ec²`; two *like* throats cannot be pushed closer than `r_e` (their 4D bodies resist
> overlap → repulsion), *opposite* throats merge/annihilate. **The make-or-break** — getting `1/r²` in 3D from the 4D-body interaction
> via **brane-localization** (the same mechanism that gave pathA_29 gravity its `1/r²`) — is the crux, detailed in §3 ⭐ v7. This
> **re-homes charge into the throat/bulk sector** (with mass, `r_e`, and annihilation), out of the flat-brane far-field.

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

> **Gravity-drain milestones (dated records moved to history).** `pathA_29` (2026-06-25, `RETURN_RESIDUAL_PREDICTION`,
> tri-review-verified) is the first end-to-end confirmation that **gravity-as-a-drain is self-consistent**: modelling the brane↔bulk
> return as a finite slab with a DC-sink completion, long-range **`1/r²` Newtonian gravity SURVIVES the finite slab** (normalizable
> transverse zero mode → `p=2`) **and** a *breathing* drain radiates an unavoidable, bounded monopole/dipole `c_s`-wave — a clean
> **falsifiable departure from GR** (you cannot have the Newtonian sink without the radiation; they are the same flux). The
> **moving-throat PDE push** (from 2026-06-26) then confirmed the distributed-geometry lift contains the structures we trust:
> Gates 1–4 (`pathA_30/31/32/33`) all land CALIBRATED (frozen-wall resonator; breathing `(a,L)` dynamics; the SO(3)-covariant ℓ=2
> quadrupole; the `54/5 = 2·27/5` quadrupole radiation split — shape earned, magnitude/`G` calibrated), and Gate 5 (`pathA_34`) =
> `FAIL_UNDERDETERMINED_NOT_PREDICTIVE` (EARNED) — linear theory can't pin the ℓ=0/1↔ℓ=2 return admittances, so **Gate 6's nonlinear
> closure (sim-deferred) must supply that selector**. Full dated records, per-gate detail, and the process lessons:
> `docs/conceptual_history.md`.

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

### 6.1 Our OWN prior internal work on the throat (survey moved to history — don't rediscover it again)

A survey of our **own prior internal throat work** — Paper-7's `notes/inner_throat/` 4D force-balance solver + equilibrium-scan
harness (`∂_a E = ∂_L E = 0`, `m_eff = E_total/c²`, calibrated n=5 EOS); the `notes/lepton_mass_notes.md` lepton functional
(`m = (18/11)·E_trapped-wave`) with its **falsified** mass tower (1:9:25 ≠ 1:207:3477); and the `stage1_solver` frozen-geometry
gauged-quintic-GPE+Maxwell throat (re-scoped OFF the critical path) — plus **THE RECURRING WALL** that *nothing in the built models
self-selects the throat size/aspect-ratio* (cheapest untested existence check = a **4D Derrick/virial scaling check**), has been
**moved to `docs/conceptual_history.md`** so it is not rediscovered.

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
12. **The wall's provenance fork — route (a) vs route (c) (2026-07-07, sharpened by `ledger_stage006`).** The brane is currently
    **route (a): a postulated order field** — `χ_B` with a postulated double-well, `DRIFT(6)`, honestly labeled
    (`ACTION_SPECIFIED_CLASSIFIED`; "wall made explicit as a postulated field, NOT derived from the one medium"). The deeper move
    — **route (c): χ_B as the medium's own orientational order (`χ_B = |P_∥|²`)** — would convert postulate into derivation and
    kill the drift, but its immediate neighbors are already falsified (the v3 arrows wall, the pathA_25 density-smectic
    `FAIL_LIGHT_STARVED`, the pathA_35 shear-surface no-go: the recurring killer is ONE coupling doing two jobs — creating the
    wall AND pinning P — which starves light). The exact composite has never itself been formally gated: a **named, high-risk
    future gate**, in the same deferred-nonlinear territory as Route A / the NG5 reductions. Until then the four-sector claim is
    honestly "four sectors from one medium **plus an imposed wall**." (Also the λγ tie-in: earning the wall's moduli from the
    same dynamics that sets `c_s` is what could make `λγ` a prediction instead of a calibration — register edge R10.)
    **The route-(c) physical picture (user, 2026-07-08 — banked as the intuition the gate will test):** each medium constituent
    carries an internal orientable micro-configuration (picture a Phlat-Ball / concentric gimbal rings — a structure with distinct
    folded/unfolded states; formally, the bulk-wide `Pⁱ` "little arrows" with the `(PⁱPⁱ−1)²` soft-unit constraint). **The phase IS
    the arrows' collective state**: on the brane the configurations are aligned/locked (arrows to the normal) → macroscopic shear
    rigidity → light; in the bulk the same constituents are individually intact but collectively disordered → zero macroscopic
    shear — exactly a demagnetized ferromagnet (every atom has a moment; the bar has no field). The throat is the
    aligner/de-aligner: the drain flow "throws" material hard enough to break alignment (`Γ_B<0` de-structuring — the shock that
    demagnetizes), the return re-orders (`Γ_B>0`). If route (c) lands, this picture becomes the derivation (χ_B's postulated
    double-well collapses into arrow physics and the drift shrinks); if it no-goes, the two-field ontology (P + independent χ_B,
    pin P7) survives and the picture stays a heuristic. **User steer (2026-07-08): hammer at this identification later in the
    ledger — the ledger is not complete until ONE clean equation set explains brane AND bulk** (the Part-VII unified parent action
    is the floor; route (c) is the prize above it).
13. **Does light's confinement (shear-free bulk) and its 2 polarizations actually fall out (§4)?** The intrinsic-(3+1)D / extrinsic-
    curvature picture *predicts* 2 polarizations + brane confinement for free, and re-frames the old "leak" as curvature, not loss — but
    this must be shown, not asserted (and it must still resolve the C5 longitudinal-mode obstruction at the throat, not in vacuum).

---

## 8. Next — the consistency knit (plan history moved out)

With all four force-sectors EARNED, what remains is the **consistency knit** (see the ⭐ FOUR-SECTOR CHAIN spine, "THE CLOSING
STEP"): (a) the cross-sector cone-lock `λγ = c_γ/c_s = 1` (the GW170817 test — do light and gravity travel at one speed in ONE
parameter set?); (b) the NG cross-consistency gauntlet (all sector couplings mutually consistent); (c) assembling the end-to-end
chain into the calibrated central `pde_ledger` (GR + EM). Honest boundary: this completes the *spec*; the full nonlinear sim stays
deferred, so sim-dependent questions stay posed and a no-go is still possible.

The earlier **consistency-gate program** that got the brane sector here — the GNLS polar-smectic gate list (Gate L/S/B/Q/K/T), the
material-state `W→φ→Gate-L` ladder, and the retired little-arrows T1–T5 ladder — is preserved in `docs/conceptual_history.md`.
Requirements spec + prior-art map: `docs/medium_requirements_and_prior_art.md`.

---

## Cross-references
- EM re-founding detail + the MacCullagh/shear math history: `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md`.
- The directive to be reworked into the T1–T5 ladder: `software/stage1_solver/directives/pathA_24_brane_existence_defect_structure.md`.
- Program state / "you are here": `STATUS.md`; `software/stage1_solver/decisions/13_emergent_constants_derivation.md` §0.
