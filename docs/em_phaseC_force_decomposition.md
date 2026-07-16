# EM force decomposition — the defect/antidefect mechanism and the gravity↔EM hierarchy

**Status:** conceptual hypothesis + test plan. **Nothing here is derived.** This is the leading physical
picture for how the model could produce a *signed* two-body electric/magnetic force (and possibly the
gravity↔EM strength hierarchy) from one native force package (body fields + boundary operator + flux law +
typed native roots) — and the honest, able-to-fail plan for
adjudicating it in **U1 Phase C (deliverable 7.5)** and the downstream U1+U2 two-body pass.

**Provenance:** a private conceptual discussion (user + Claude) surfaced the idea; a Codex physics audit
(`_scratch/_HANDOFF_post_throat_bulk_flow_force_decomposition.md`, 2026-07-16) supplied the rigorous
decomposition and corrections; Claude synthesized this canonical document. Persisted here because it is a
conceptual insight we do not want to lose. Companion docs: `docs/model_definition_audit.md` (§F to-define
list, §G the ~10⁴² hierarchy), `docs/conceptual_foundation.md` (native mechanisms), `docs/em_analog_next_phase_handoff.md`.

Read-first caveat (matches the whole EM sector): the recurring failure mode is IMPORTING standard EM at
under-defined gaps. This doc names a mechanism and a test; it does **not** assert the mechanism holds. Every
outcome below — including "the throat-current mechanism fails" — must remain a live landing.

---

## 1. The central correction: raw outflow does NOT explain attraction

Gravity in this model is the charge-even **drain** — medium flows radially inward on the brane toward each
throat (universal, `1/R²` far force within the postulated slab branch). Its **attractive sign is a property of
the accepted gravity construction** — the two-body far-field sign is the target-blind `pathA_21c` signature,
while `pathA_29` establishes the drain's inward *sink* character and localized falloff, **not** the two-body
attraction sign by itself — and it is **not** a consequence of the drain being unsigned: even parity in `±w`
only establishes that gravity is *insensitive to charge sign*, not that the force is attractive. The medium then continues *through the mouth into the bulk*. It is tempting to reason: "both throats
are bulk **outflows** (sources); two sources interact → attraction."

**That reasoning is wrong here, for two independent reasons:**

1. **Potential flow does not fix a universal sign.** Potential flow and Bernoulli are *inviscid* and do apply
   in a superfluid — but they do **not** give a universal attraction sign for maintained sources/sinks. The
   result depends on boundary work, on what is held fixed, and on how the force on the defect is defined.
   (Viscous jet entrainment is a *separate* mechanism, and *that* is the one genuinely unavailable in the ideal
   medium.) The native stress/action decides the sign, not a pressure slogan.
2. **The streams need not even overlap.** A `+w` throat's post-mouth stream primarily occupies one half-bulk,
   a `-w` throat's the other. If their bulk supports do not overlap and no common mode spans the brane, the
   orientation-odd cross interaction may simply be **zero**.

For a global orientation sign `s ∈ {+1,-1}` (the `±w` puncture) to produce a *force*, some **shared** field or
boundary-value problem must see it: the wall displacement, sleeve strain, orientation lock, collar stress, a
bulk return field, or a mode coupled across both sides. Only the cross term of the two complete solutions can
determine attraction or repulsion. And there is no literal time ordering: the two-body steady state is **one
global boundary-value problem** ("enters the bulk *before* interacting" is not a mechanism).

---

## 2. The leading hypothesis: defect/antidefect wall-healing

The stronger candidate is that **`+w` and `-w` throats are conjugate geometric-defect candidates** (a possible
defect/antidefect pair — note `±w` is presently a `Z₂` *orientation*, **not** a proven topological charge or
winding; see the topology gate in §7), each carrying core energy **plus a long-range region of signed
wall/sleeve strain**:

- a `+w` throat imposes a signed sleeve/wall continuation into one side;
- a `-w` throat imposes the conjugate continuation;
- an **opposite (`+w,-w`) pair** *might* have its signed far distortions **cancel** as it approaches (wall heals,
  shared deformation energy drops → attraction);
- an **equal (`+w,+w`) pair** *might* have them **add** (energy rises → repulsion).

**This sign is a hypothesis, not a theorem.** Kink–antikink attraction is plausible but is *not* universal for
elastic defects — fixed-displacement, fixed-force, capillary-like, and topological-soliton boundary conditions
can each give a different sign. Start from the sign-agnostic generic decomposition:

```
U(R; s1,s2) = 2·E_core + U_00(R) + s1·s2·U_11(R)
```

and only write the asymptotic **landing**

```
U_00(R) ~ −A_g/R ,   U_11(R) ~ +A_q/R      (A_g, A_q > 0)   [Coulomb-sign landing]
```

**after** demonstrating (a) a gapless localized mode, (b) a nonzero monopole, and (c) the native sign from the
full calculation. Local annihilation/healing at `R ≈ a` does **not** by itself imply a long-range `1/R` tail
(range gate, §7).

Two cautions before calling any of this "attraction":

- **The throat is open and driven, so a scalar pair potential may not exist.** The general two-body force is
  `F_pair = F_var + F_flux + F_𝓑 + F_rad`; the `E_pair`/`U(R)` above capture only the *conservative/static
  variational* contribution. Flux, constraint, return-path, and radiative forces must be included before naming
  the result — a drop in wall energy can be outweighed or reversed by momentum carried through the mouths.
- This is nonetheless the piece that could rescue electrostatics from the current "no-lock limit just attracts"
  problem (`docs/model_definition_audit.md` §F, I2: the electric sign presently rests on an assumed lock `G0>0`;
  a healthy unconstrained scalar with a conventional linear source attracts like-sources after being integrated
  out).

**Essential distinction (do not conflate):** a large barrier to *create or destroy* a throat is
separation-**independent** once the bodies are far apart → it contributes to **mass** (`E_core`), producing
**no force**. Attraction requires the `R`-**dependent** shared deformation energy that decreases as a conjugate
pair approaches. Annihilation and charge/topology conservation are separate downstream questions
(see the `±w` Z₂-topology + dynamical conservation note).

---

## 3. The four-channel decomposition (the honest bookkeeping)

Split a stationary one-throat configuration into orientation-even and orientation-odd parts on a genuinely
`R_w`-symmetric background:

```
Phi_s = Phi_even + s·Phi_odd ,   v_s = v_even + s·v_odd ,   s ∈ {+1,-1}
```

`v_even` = common radial drain; `v_odd` = the `R_w`-odd axial post-mouth component. The two-body interaction
then organizes into **{static / moving} × {charge-even / orientation-pair}**:

| field overlap                    | `s` parity | velocity order | physical classification      |
|----------------------------------|-----------:|---------------:|------------------------------|
| radial / common drain            | even       | `V⁰`           | **gravity**                  |
| oriented post-mouth body/flow    | `s1·s2`    | `V⁰`           | **static electric** candidate |
| moving common drain              | even       | `V1·V2`        | **gravitomagnetic** candidate |
| moving oriented post-mouth body  | `s1·s2`    | `V1·V2`        | **EM-magnetic** candidate     |

Consequences (naming discipline):

- A stationary orientation-dependent attraction/repulsion is a **static electric** candidate — **NOT** magnetism.
- A term earns the name **magnetic** only if it is generated by *moving* the oriented throat-body, is bilinear
  in the two bodies' velocities at leading order, and **vanishes when either velocity is zero**.
- The one-body coupling to the bulk is only a prerequisite; a two-body force needs a genuine **cross term**
  between the two induced configurations (shared mediator or boundary/constraint reaction).

**Full four-orientation Hadamard decomposition** (run this, do not presume the sign). Because the throat is
open/driven, the **primary object is the channel-resolved force**, not a potential — decompose each force channel
`ch ∈ {var, flux, 𝓑, rad}` separately:

```
F^ch_{A,ab}(R) = (1/4) Σ_{s1,s2}  s1^a · s2^b · F^ch_A(s1,s2) ,   a,b ∈ {0,1}
  (a,b)=(0,0)   charge-even / mass channel
  (1,0),(0,1)   should vanish on a symmetric, exchange-symmetric background
  (1,1)         orientation-pair channel (electric candidate)
```

Collapse a channel to a scalar potential decomposition `U_ab(R)` **only** after its conservative contribution
passes an integrability check — some cells (flux, constraint, radiative) may have **no** scalar potential at all.
Nonzero `(1,0)`/`(0,1)` components diagnose **ambient asymmetry** (the seam between the postulated two-sided
`R_w`-symmetric background and the one-sided slab that `pathA_29` actually computes), **not** a clean charge
interaction.

Moving throats — expand the complete interaction, do **not** insert a current by hand:

```
U_12 = U_mass(R) + s1·s2·U_orientation(R)
     + V1_i·K_mass_ij(R)·V2_j                # gravitomagnetic
     + s1·s2·V1_i·K_orientation_ij(R)·V2_j   # EM-magnetic candidate
     + other allowed one-sided/parity-mixed structures + higher orders
```

The desired `j ~ s·V` behavior is an **output classification of `K_orientation`**, not an input. It may instead
be zero, nonlocal, anisotropic, higher-order in `V`, or carry extra tensor structure.

---

## 4. The gravity↔EM hierarchy (~10⁴²) — a plausible shared origin

This does **not** solve the hierarchy, but it is the first mechanism that plausibly ties the *sign* and the
*strength* together, and it is the right *kind* of mechanism (a `(v/c)²` factor could never make 40+ orders).

1. Brane tension / orientation-lock energy presents a large barrier to material crossing the wall.
2. That barrier selects a very small steady throughput `ṁ`, **suppressing the charge-even drain** = weak gravity.
3. The signed throat orientation is an **order-one geometric defect candidate** (topological only if the topology gate passes), not the tiny transmitted flux.
4. The same stiffness that throttles throughput **stores large energy in the fixed signed defect** = strong electric.

If the native throat solve produced a tunneling/nucleation-like throughput `ṁ ~ exp(−S_wall)`, then a gravity
coefficient quadratic in throughput carries `exp(−2·S_wall)` — an exponential, which is how you manufacture many
orders of magnitude. A purely classical drain need not produce such an exponential, so this **cannot be assumed**.

**The load-bearing fork (what makes it falsifiable, not a slogan) — the mouth "ensemble":**

- **Fixed conventional source** (normal force): response scales `~ J²/K`. Stiffer wall → **weaker** coupling.
- **Fixed displacement / winding / topological defect**: stored energy scales `~ K·q_top²`. Stiffer wall →
  **stronger** coupling.

"The brane barrier makes the electric force strong" is viable **only if** the throat imposes a *fixed signed
geometric/topological boundary condition*, not an arbitrary external scalar source. U1/U2 must **derive** which
ensemble applies — `pathA_38`'s reduced source projection does not settle it. Also possible: a
localization-normalization ratio `d/ℓ` (slab scale vs wall thickness) entering the two couplings; power and
scale relation currently underived → assigning `10⁴²` to it now would be numerology.

Operationally, the hierarchy must be read from the **same solved throat**:

```
F_orientation / F_gravity = A_q / A_g
```

with `ṁ`, the signed boundary datum, wall stiffness, zero-mode normalizations, and return closure fixed by the
native solution (not fit afterward). At present `ṁ`, `Q_E`, the nonlinear mouth condition, and return closure
are all OPEN → **the hierarchy is not yet predicted.** (Cross-ref `docs/model_definition_audit.md` §G, which
frames `F_e/F_g` as the sharpest held-out both-sectors test — currently FIT, not tested.)

**Species qualifier:** `A_q/A_g` is **not** a universal constant — it carries the two bodies' effective mass and
charge source strengths. `F_e/F_g ≈ 4×10⁴²` is the *electron* pair; a *proton* pair is `~10³⁶`. The target is
the whole family of these ratios, not a single number, and any derivation must reproduce the family, not one point.
Moreover `A_q/A_g` is a well-defined *ratio* only if **both** channels land on the **same** `1/R` potential; if
they do not, compare the far radial-force coefficients directly and report the result as **distance-dependent /
non-Coulombic** rather than a single number.

---

## 5. What is already in place vs what is missing

**In place (the framework does not ban the mechanism):**

- `docs/em_u1_body_definition.md` defines the full localized throat-body (one-sided sleeve/`χ` tube, geon
  content, bound drain + motion-induced flow) and records that the drain's `v_w` component is `R_w`-odd.
- U1's boundary map spans the ensembles: E1/E2 impermeable no-slip/free-slip, E3 permeable translating texture,
  E4 nonholonomic throat↔brane-shear lock, E5 dissipative slip. The open-system law
  `F_var + F_flux + F_constraint + F_Rayleigh + F_rad = 0` is the correct container (no closed point particle).
- B2 supplies the bookkeeping layer (native mass/momentum/energy currents, intake response, closure branches,
  radiative residues) and genuinely derived the wall `χ–u` Hessian block at the frozen symbolic base.

**Missing (why current results do not answer the question):**

- `pathA_29`: conditional slab-localized charge-even drain falloff. Does **not** compute the cross interaction
  of two oriented post-mouth streams.
- `pathA_38`: conditional `h`-mediated static falloff/source projection. Does **not** derive it from the full
  drain/sleeve interaction; the nonlinear mouth boundary action (what the throat actually holds fixed) is absent;
  full native static sign OPEN.
- `pathA_39`: transverse/longitudinal moving-field kernels **given** declared sources `∝ s·V`. Does **not**
  derive those sources from throat translation.
- Stage-3: `O(V)` parity mixing is symmetry-**allowed** with nonzero witness structures. Does **not** compute
  the realized finite-mouth coefficient.
- B2: post-mouth closure is OPEN (whether intake grows the body, exits along the sleeve, or uses a return path);
  many response slots typed `UNRESOLVED`. B2 does **not** prove a stream-stream interaction sign.

No existing result yet says the post-throat bulk flow is the physical source of electricity or magnetism.

---

## 6. Phase-C adjudication plan (the decisive one-body → two-body pass)

U1 Phase C (7.5) already demands the body action + surface terms, native bulk sources `J_A`, the full matrix
operator perturbation `δO_AB` (incl. off-diagonal mixing), constraint/virtual-work channels, far multipoles of
the *total coupled* response, `s`-parity + leading `V` order, and a mass/drain-vs-orientation/charge
decomposition with a mass-only ablation. Make these diagnostics explicit so the mechanism cannot hide inside an
aggregate coupling map.

**Ownership — do NOT overassign to Phase C (it is a *one-body* deliverable):** (A) and the one-body portions of
(D) are **Phase C**. (B), (C), (E), and the two-body ablations are the **downstream U1+U2 two-body consumer**
(`docs/em_u1_body_definition.md` — the static/moving two-body signs are explicitly downstream, not Phase C).
(F) hierarchy extraction is **beyond current U1** unless an upstream amendment supplies a native
throughput-selection law: present U1 treats `ṁ` as a *fixed OPEN input* and only checks whether motion forces it
to vary, so Phase C **cannot** derive `ṁ ~ exp(−S_wall)` or select its magnitude within current scope.

- **(A) One-body post-mouth parity census.** For each field/profile and each U1 endpoint, derive the `s`-even /
  `s`-odd components of the stationary solution (axial drain, sleeve, wall, surface-flux). Tag each parity claim
  as body-only, ambient-postulate-dependent, or one-sided-asymmetry-map.
- **(B) Static two-body orientation decomposition.** One action/stress/flux prescription; compute all four
  `{s1,s2}`; extract `U_00,U_10,U_01,U_11`; report the drain, post-mouth bulk-flow, sleeve/wall, `h`,
  longitudinal-mode, and control-surface/return contributions **separately**.
- **(C) Moving two-body kernel decomposition.** Translate the complete one-body solution, solve the induced
  response, extract every `O(V1·V2)` tensor structure; project onto `K_mass` and `s1·s2·K_orientation`. **Do not
  project onto `s·V` until the full response exists.** Then classify the orientation channel: exactly
  convection-like / convection-like + departures / departure-only / null / unresolved-on-named-missing-closure.
- **(D) Required ablations** (each able to fail at its own assert): common drain/mass source; `w`-odd
  sleeve/orientation structure; post-mouth axial flow; `h` coupling; brane transverse shear; longitudinal
  compression; E4 constraint; outer control-surface flux/return. **The magnetic candidate must vanish when the
  orientation structure is removed; removing translational velocity must eliminate every magnetic term while
  leaving the static orientation term.**
- **(E) Physical comparison** (a Maxwell landing needs more than `1/R²` + attraction): vanishing at `V=0`;
  reversal under one *derived* current reversal; transverse/Darwin tensor structure; **relative coefficient** to
  the static orientation force; propagation poles/speeds; longitudinal/scalar residue; preferred-frame pieces;
  dependence on OPEN boundary/closure choices.
- **(F) Shared wall-barrier / hierarchy extraction.** From the **same** throat solution determine (not fit):
  whether the wall barrier selects `ṁ` and by what functional form; whether the signed datum is fixed-source /
  fixed-displacement / topological; `A_g` and `A_q`; which localization normalizations enter each; whether
  core/barrier energy is strictly `R`-independent; and `A_q/A_g` **with no new parameter inserted after either
  force is read off.** Run controls that exchange fixed-source ↔ fixed-defect ensembles — a sign or hierarchy
  that flips under this exchange is conditional on unresolved mouth physics and must not be promoted.

---

## 7. Evaluation rubric — predeclare the outcomes, INCLUDING the kill

The Phase-C derivation must **accept every outcome**, but model evaluation must still be willing to conclude an
outcome is *too structurally different to call electromagnetism*. The biggest risk now is not wrong algebra — it
is letting "characterized departure" quietly mean **every** outcome counts as success (the unfalsifiability /
Lorentz-ether trap). Predeclare the questions and the landings.

**Predeclared questions:**

- Does a generic charge-odd `O(V)` source exist after the mass/drain channel is ablated?
- Is it `∝ s·V`, or does it carry additional tensor structure?
- Does the *same native force package* — same body fields, boundary operator, flux law, and typed native roots
  (not necessarily one action; E4/flux/Rayleigh are not all action-derived) — give the static sign **and** the moving sign?
- Does the transverse term have the Maxwell/Darwin tensor form **and relative coefficient** — not merely `1/R²`?
- How much residue lies in the propagating `h` / longitudinal-scalar modes?
- Are `c_E`, `c_γ`, and the coupling coefficients **tied**, or independently fitted? (= the surplus-vs-knobs ledger.)
- Do preferred-frame terms remain? (= the Lorentz-violation falsifier.)
- Is there any approximate additive conserved charge in the multi-throat theory?
- Which results rely on the postulated ambient `R_w` symmetry vs the one-sided slab actually computed?

**Predeclared landings (name them honestly, with quantitative cut-lines where possible):**

- **Maxwell-like low-energy EFT** (Darwin relative coefficient within a predeclared tolerance).
- **Scalar–vector characterized departure** (right sign/structure, wrong relative coefficient beyond tolerance).
- **Signed two-body force analog, but NOT EM** (orientation sign exists but tensor structure / propagation fail).
- **Null or wrong-sign result — the throat-current mechanism FAILS.**

**Predeclared gates (each able to fail):**

- **Topology gate.** Is `±w` a genuine topological charge/winding, or only a `Z₂` orientation? Test whether the
  two sectors are disconnected, whether a finite-energy interpolation exists, and whether opposite throats can
  actually annihilate. Until shown, call them *conjugate geometric-defect candidates*, not established
  defect/antidefect topology.
- **Range gate.** Classify the potential explicitly — `1/R`, Yukawa, bulk `1/R²`, dipole, or null. Do not assume
  the desired `1/R` potential (`1/R²` force) prematurely; local healing at `R≈a` is not a long-range tail.
- **Stability gate.** A correct force sign obtained via a ghost, negative-norm mode, or unstable background is a
  **failure**, not a success (cf. the historical `P–u`-coupled long-wavelength instability on the now-*retired*
  brane-polar-field baseline — operative U1 Phase A has no `P–u` block).
- **Momentum-closure gate.** Verify total body-plus-medium momentum balance. Equal-and-opposite body forces need
  not hold on their own in an open system, but any imbalance must be carried by a **named** flux/return/radiation
  channel — not left implicit.

**Add (discipline):** make the "not EM" thresholds *numeric* wherever possible, so the departure is able to fail
per finding, not just describable after the fact.

---

## 8. Interpretation guards (failure modes)

1. Do not call a `V⁰` orientation effect magnetic — it is the static electric candidate.
2. Do not assume same-side flow ⇒ attraction and opposite-side ⇒ repulsion — it may be attraction-vs-zero or a return-closure artifact.
3. Do not import viscous jet entrainment into an ideal superfluid — the native stress/action decides the sign.
4. Do not infer `s·V` from geometry alone — a moving orientation can yield charge-even, higher-order, or no transverse source.
5. Do not let the gravity drain masquerade as charge — preserve the mass-only ablation, while allowing genuinely `w`-odd drain-derived cross terms to be classified (not banned by slogan).
6. Do not ignore the ambient seam — the executable `pathA_29` background is one-sided; clean `s1·s2` parity uses the postulated two-sided `R_w` background. Run both branches or give an explicit branch map.
7. Do not assume one body's interaction with the bulk is a two-body force — it must descend from the cross interaction of both bodies' fields or a shared constraint/flux channel.
8. Do not tune the post-mouth contribution independently after seeing the result — shared profile/coefficients must stay shared so the hierarchy test retains meaning.
9. Do not mistake opposite coordinate flux for a physical source-sink pair — both orientations are outlets; a common signed field/boundary interaction must be demonstrated.
10. Do not use a large throat-creation barrier as a force by itself — only its separation-dependent pair deformation attracts; the isolated core is mass.
11. Do not infer that greater wall stiffness always means stronger charge — true for some fixed-defect ensembles, false for conventional fixed-source response. Derive the mouth ensemble first.
12. Do not claim the hierarchy from a suggestive large ratio — power, normalization, and parameter relation must come from the same native throat solution; `d/ℓ` or an exponential is only a candidate structure until derived.

---

## 9. Bottom line

The conceptual idea is worth preserving and explicitly testing; it does **not** require abandoning U1/B2 or
inventing a new field. The refined leading hypothesis is **not** "two outflows attract." It is: **opposite throat
orientations are conjugate geometric-defect *candidates* whose shared deformation energy *may* cancel (opposite
→ heal → attract; like → add → repel), while the same wall tension *may* throttle the common drain** — a
plausible *common* physical origin for opposite-charge attraction **and** a strong-electric/weak-gravity
hierarchy. The sign is a hypothesis, not a theorem; the pair-energy picture is only the conservative/static piece
of an open-system force (`F_var+F_flux+F_𝓑+F_rad`); and it all remains conditional until the topology, mouth
boundary ensemble, throughput selection, and pair energy are derived from one native force package (body fields,
boundary operator, flux law, typed native roots — not a single action, since E4/flux/Rayleigh are not
action-derived). It must stay able to fail (null / wrong-sign is a real landing).

**Standing review question for the Phase-C contract + U1/U2 two-body consumer:**
Does the operative contract explicitly force the four-way decomposition `{static/moving} × {charge-even/
orientation-pair}`, distinguish raw post-mouth flux from defect/antidefect wall energy, determine the
fixed-source-vs-fixed-defect mouth ensemble, and extract `A_q/A_g` from one shared throat solution? If not, add
the smallest non-prejudicial amendment that inserts these diagnostics **without** assuming their sign, hierarchy,
source law, or success.

**Do not** re-open or reinterpret the completed B2 stage-0 HALT artifact (its schema `status` is
`AWAITING_ORCHESTRATOR_APPROVAL`) to force this answer. B2 supplies bookkeeping; Phase C and the two-body
calculation supply the physics.
