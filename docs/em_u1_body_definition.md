# U1 — the `±w` throat-body as a dynamical object (conceptual definition)

**Status: DRAFT v3.1 — 2026-07-12.** v1 → Codex design-review (`U1_NEEDS_FIXES`, 12 items) → v2 → Grok independent review (`U1_V2_NEEDS_FIXES`, 3 BLOCKER + 6 MAJOR + 3 MINOR; verified the v2 folds genuine) → v3 → Codex confirm-pass (`U1_V3_FOLD_ISSUES`, 2 localized mis-folds) → **v3.1 fixes both** (fold logs: §10). Next: Codex re-confirm. This is the keystone definition of the EM-analog next phase (`docs/em_analog_next_phase_handoff.md`; `docs/model_definition_audit.md` §F): what the `±w` throat-body **is** as a *moving* object in the medium — composition, kinematics, response law, couplings — **in the model's own terms**, so the IMPORTED EM pieces (`j=sV`, the electric sign, the exchange machinery) can later be *computed* instead of borrowed. Tags: **[COMMITTED]** / **[DEFINED-HERE]** / **[POSTULATE]** (new input, declared) / **[HYPOTHESIS]** (must be computed) / **[COMPUTABLE(...)]** (an output, conditional on the parenthesized inputs) / **[OPEN]** / **[NOT-U1]**.

**Anti-import rule (standing):** nothing here may borrow a point-charge, a Maxwell current, an exchange sign, closed-body classical mechanics, or `E=mc²`-style inertia. Every dynamical property is **derived from the model's own action/fields** or explicitly tagged [POSTULATE]/[OPEN]. Whatever the downstream signs/current-law come out to — Maxwell-like or a departure — that is the answer.

---

## 0. Inherited commitments (the walls U1 must build between)

1. **[COMMITTED]** One compressible superfluid in 4+1D; our 3D space is a brane (ordered-phase domain wall at `w≈0`); the brane carries MacCullagh rotational shear (postulated stiffness `μ_R`); the bulk is superfluid (supersolid lattice = a live [HYPOTHESIS] in the ontology, not assumed here).
2. **[COMMITTED]** A particle is a **throat**: a puncture in the brane whose wall-structure extends into the bulk in one `w`-direction. **The ± of that direction is the charge `s`** — a direction, *not* a winding. Not a point particle; not a bulk vortex.
3. **[COMMITTED]** The throat **drains**: one-way radial inflow = the gravity attribute; `1/r²` via the finite-slab zero mode (`pathA_29`). NOTE: `pathA_29`'s *executable* slab family is **one-sided** (brane at `w=0`, return/absorber at `w=d`) — relevant to §1's ambient-symmetry declaration. That slab zero mode (gravity localization) is a **different object** from the throat's in-brane translational zero mode (§7.1) — do not conflate. Mass = trapped configuration (geon) **[POSTULATED mechanism]**.
4. **[COMMITTED]** Static electric force = the two throats' 4D bodies-beyond-the-mouth (`pathA_38`, scalar-`h`-mediated, conditional). Magnetism = the **moving** 4D throat-body (`pathA_39` committed the object; its dynamics were never defined — that gap is U1). Distinct from gravitomagnetism. `pathA_39` Stage 3 found **symmetry-allowed unprotected parity mixing under motion**, partly as **operator perturbations on the mediators** (not only sources) — §7.5 must carry both. `pathA_39` also **separated mass-sourced (`q_M→u_L`) from charge-sourced residues** — §7.5 must preserve that split.
5. **[COMMITTED — computed]** No native first-class `U(1)`: `NATIVE_P_NO_EMERGENT_GAUSS`. Charge is natively `Z₂`. U1 does not try to fix this (U3, audit §F).
6. Throat geometry: mouth radius `a`, `w`-extent `L`, with `L/a = 37/20` (**a frozen branch ansatz, not self-selected**; `ANSATZ_LEDGER.md`, audit §A). Carried symbolically; not re-derived.

---

## 1. Composition, conjugation, and the ambient background **[DEFINED-HERE]**

The `±w` body is the **complete localized field configuration** at the puncture:
- the **sleeve** — the wall-field (`χ`) configuration wrapping into the bulk: a tube of ordered-phase material, mouth radius `~a`, extent `~L` in `w` on ONE side (`s=±1`);
- the **geon content** — the trapped (nonlinear) configuration constituting the mass **[POSTULATED mechanism; profile = declared open input, §7.0]** — part of the body *from the start*;
- the **bound flow field** — drain inflow + motion-induced flow.

Nothing else is available to *be* the body. Every dynamical property is a functional of (`χ, ρ, v`) — or a tagged new input.

**Body conjugation (a definition):** the `s=−1` body **ansatz** is *defined* as the `w`-mirror of the `s=+1` body ansatz: `Φ₋(x,w) := R_w Φ₊(x,−w)`.

**Ambient `R_w` symmetry (a separate, declared input — NOT implied by the above):** whether the *environment* is `w`-mirror-symmetric is a property of the background, and the committed executable gravity family (`pathA_29`, one-sided slab) does **not** have it. U1 declares: **[POSTULATE — U1 background] the conceptual background is `w`-symmetric about the brane (medium on both sides; a `+w` and a `−w` body see equivalent environments)**, with the one-sided `pathA_29` slab treated as a computational device whose reconciliation with this postulate is an **[OPEN]** seam (§11). Every parity statement below says which of the two (body conjugation vs ambient symmetry) it uses; gate G5 is scoped accordingly.

Parity consequences (computed per-observable, never per-mechanism):
- equal drain rates for `±s` follows from *ambient* `R_w` **[so: POSTULATE-dependent]**;
- the drain's `v_w` component is `R_w`-odd — **sign-odd drain-derived structure is allowed**; a charge-dependent cross term is *not* automatically a definition failure;
- the slogan "drain = sign-shared, orientation = sign-odd" survives only as a heuristic about specific even observables under the ambient postulate.

**Frame declaration:** the U1 lab frame is the **medium rest frame** (asymptotic bulk at rest); `V` is velocity relative to it. This is native for a superfluid analog (and matches `pathA_39`'s preferred-frame diagnostics); no Lorentz-covariant form may be smuggled in later by notation.

## 2. Kinematics + the U2 seam `𝔅`, made operational **[DEFINED-HERE]**

"The throat moves with velocity `V`" means: the **pattern** (the §1 configuration) translates in-brane. Whether and how the surrounding medium is *tied* to that translation is exactly the unknown; U1 must not pre-answer it in either direction. *(v2's "the medium does not ride along; no material plug is dragged" was itself a pre-answer — it describes only the permeable endpoint below.)*

**`𝔅` — the sleeve-surface boundary/constraint operator (the U2 object), declared as an input.** It is not assumed reducible to one scalar (the `ζ` of `directive_native_vortex_sign.md` is a one-parameter caricature). U1 computes every response coefficient as a function of `𝔅`, at minimum at these **operationally-defined endpoints**:

| Endpoint | Fields constrained | Condition (schematic) | Variational class |
|---|---|---|---|
| **E1 impermeable + no-slip (fluid)** | bulk `v` at the sleeve surface `Σ` | `v·n̂=V·n̂` and `v_∥=V_∥` on `Σ` | holonomic BC (Lagrangian) |
| **E2 impermeable + free-slip** | bulk `v` normal only | `v·n̂=V·n̂` on `Σ`; tangential stress-free | holonomic BC (Lagrangian) |
| **E3 permeable / phase-texture** | none beyond the pattern itself | medium passes through the translating `χ`-texture; drain flux unimpeded | bulk action only (Lagrangian) |
| **E4 nonholonomic shear-lock** | throat velocity ↔ **brane shear** `u_T` (curl) at the collar | velocity-level constraint `g(V, u̇_T)=0` (the ontology §3/§9 "roll") | **constraint** (multiplier / virtual work) |
| **E5 dissipative** | tangential slip with loss | slip ∝ tangential stress (Rayleigh function `R`) | **non-variational** (Rayleigh) |

**Terminology guard (Grok M1):** the ontology's decisive "roll vs slip" pivot is **E4** — a *nonholonomic throat–shear constraint* — NOT the fluid no-slip E1. v2 conflated them; they are distinct endpoints and U2 must adjudicate them separately. Where `𝔅` physically lives (mouth collar, bulk — possibly lattice-hosting — or both) is U2's question; U1 does not fix the location.

## 3. The one-body law — a four-channel force balance, not a bare Hessian **[DEFINED-HERE]**

The body is an **open system** (one-way intake; the return closure and drain-vs-charge budget are inherited OPEN), and `𝔅` includes non-variational endpoints. A single `S_eff` Hessian is therefore **not** a well-posed definition of its dynamics (Grok B1). U1 defines the one-body law as a **force balance over the collective coordinates** `q_A ∈ {X, p}`:

`F_A^var + F_A^flux + F_A^𝔅 + F_A^rad = 0`,  with `F_A^𝔅 = F_A^constraint(𝔅) + F_A^Rayleigh(𝔅)`

1. **`F^var` — the variational channel.** From a precisely scoped **`L_eff`** (a Lagrangian *density in time*; not "`S_eff` per unit time" — steady `V` makes an action diverge): the model's action on a **declared control volume** `Ω_c` **plus its explicit surface action on `∂Ω_c` and `Σ`**, with the body's fields embedded via collective coordinates (§7.1, including the zero-mode quotient). `M_AB(𝔅)` (inertia tensor) and `Ω_AB(𝔅)` (symplectic/Berry form) are defined as the `V²`-Hessian and first-order form **of `L_eff` only** — flow and geon contributions arise in the same calculation, cross terms retained; no `ρa³L` scaling, no `E/c²` term pre-asserted (whether an `E/c²`-like relation *emerges* is a result).
2. **`F^flux(closure)` — the control-surface channel.** Momentum/mass/energy carried across `∂Ω_c` by the drain and return (§5 bookkeeping). **Partition rule (anti-double-count):** any flux contribution scaling as `∝V̇` is reported **either** inside `M_AB` **or** inside the intake coefficient `C_ṁ` — never both; the report must state the partition.
3. **`F^𝔅` — the boundary/constraint channel**, decomposed:
   - **`F^constraint(𝔅)`** — **Lagrange–d'Alembert / multiplier reaction forces** for constraint-type endpoints. This is where the ontology's decisive **E4 nonholonomic shear-lock** reaction enters the one-body equation (the E4 constraint is velocity-level; its reaction is NOT obtainable from `δL_eff` and NOT dissipative — it needs its own slot). **E1's** no-slip reactions live in `F^var` when E1 is imposed as a holonomic BC on the fields inside `L_eff`, otherwise here as multiplier reactions — the report must state which. Zero at E2/E3.
   - **`F^Rayleigh(𝔅)`** — the genuinely dissipative channel (E5; Rayleigh function), zero at E1–E4.
   Every §7 deliverable that depends on `𝔅` reports its result under this same decomposition (which channel each contribution sits in).
4. **`F^rad` — the radiative channel.** Hereditary in general (memory kernels); **the local form `F = M V̇ + Ω V` is explicitly the sub-radiative, slow-acceleration truncation** — not the definition (Grok M6). §7.7 computes the per-channel residues and may promote memory kernels.

**Response-structure gating (from v2, kept):** zero winding kills exactly the **intrinsic-circulation Magnus** term; it does not forbid first-order terms in general. `Ω_AB` is computed with three separately-gated channels: (i) intrinsic-circulation Magnus (expected zero, verified — the expectation is *tested*, not assumed); (ii) translation–translation Berry curvature (wall-texture/geon-phase terms can supply it without circulation); (iii) translation–tilt gyroscopic coupling. A nonzero Berry/gyroscopic term is a first-class result, not a nuisance.

**Dissipation content:** no *viscous* drag is assumed in the ideal subcritical bulk, but the committed field content has several radiative channels — compression (`u_L`/`c_s`), brane shear (`u_T`/`c_γ`), the **charge-coupled `h` mode**, wall/`χ` modes, internal geon modes. §7.7 enumerates the native action's modes and computes which have nonzero moving-body residues. An open sink is not automatically conservative — hence §5.

## 4. Open-system bookkeeping — the control surface **[DEFINED-HERE]**

U1 declares a control surface `∂Ω_c` around the body and does explicit mass/momentum/energy flux bookkeeping:
- whether intake grows the body's mass, exits along the sleeve, or rides the (open) return path is **[OPEN — declared, carried symbolically]**; U1 *accounts* for it, it does not resolve it;
- the **intake-momentum response** `C_ṁ(𝔅, closure)` (reverse-rocket-like) is **[COMPUTABLE(𝔅, closure)]**, reported under the §3 partition rule;
- first pass holds `ṁ` fixed **[S2]**, with a computed check whether the equations force `O(V)` variation of `ṁ` (G7).

## 5. The tilt observable `p` **[HYPOTHESIS]**

- **Invariant observable:** `p` = the in-brane projection of the sleeve's tangent (equivalently the leading in-brane multipole of the sleeve `χ`-configuration). **`p` is used everywhere in §7 — the collective orientation coordinate is `p` itself** (v1's signed `θ = s·p` inserted charge-oddness by convention and is banned from the deliverables; Grok M3). Distinguish steady `p(V)` from dynamic susceptibility `p(ω,V)`.
- **Parity is computed, not assumed.** Under body conjugation + ambient `R_w`: `p(−s,V) = p(+s,V)` — charge-**even**. If the steady response is analytic and the in-brane background isotropic, the leading nonzero `p(V)` is **linear in `V`**.
- **Tilt alone cannot establish a current.** Current-likeness is decided ONLY by the §7.5 coupling package — a charge-even `p` can still yield a charge-**odd** mediator coupling once the mediators' own `R_w` transformations enter (even×odd = odd), and vice versa; the check is channel-by-channel. No commitment forces or forbids tilt; it stays distinct from the graveyard "moving throat looks compressed" *only as an orientation mode, not a compression proxy*.
- Honest outcomes (all first-class): the computed coupling's `O(V)` part ∝ `sV` (the import becomes computed); different structure (**characterized departure**); or zero (this magnetism source dies natively).

## 6. First-pass scoping **[DEFINED-HERE — declared, revisitable]**

- **S1 — Rigid collective coordinates** `{X, p}`; deformation modes deferred; validity regime reported (§7.1).
- **S2 — Fixed `ṁ`** with the forced-variation check (G7).
- **S3 — Single body.** Two-body computations (static → electric sign; moving → magnetism sign + current law) are downstream consumers of U1+U2 **[NOT-U1]**.
- **S4 — All outputs conditional** on: the §7.0 declared inputs, `𝔅`, the §4 closure declaration, and the §1 ambient postulate. Conditionality stated per-output, never silently dropped.
- **S5 — Local-EOM truncation** (§3.4) declared wherever used.

## 7. Deliverables (each **[COMPUTABLE(declared inputs, 𝔅, closure, ambient)]**)

**7.0 — Base configuration + declared inputs (FIRST).** A native one-body ansatz: sleeve profile (`χ` shape, thickness, density contrast), geon content (profile declared open where unknown), drain flow — with the **explicit list of action coefficients and background profiles needed beyond `(ρ, a, L, ṁ, μ_R)`** (χ-kinetic coefficient, wall stiffness, sleeve bending/anchoring stiffness, EOS/wave operators, IR scheme, …). Anything unavailable = **[OPEN input]**; downstream outputs conditional on it.

**7.1 — Collective-coordinate embedding** of `{X, p; s}` into the native action → `L_eff[X, Ẋ, p, ṗ; s; 𝔅]`, **including the explicit zero-mode quotient** (projection/moduli-fixing recipe for the translational Goldstone — the field Hessian is degenerate along it, so "integrate out" without the quotient is singular; distinct from `pathA_29`'s slab zero mode). Includes: canonical momenta (with first-order/Berry terms), generalized force defined as a field-functional variation, the rigid-approximation validity regime, internal-mode stability. `s` is fixed within the one-body sector (no global charge-conservation claim — non-native per §0.5).

**7.2 — Inertia tensor** `M_AB(𝔅)` from `L_eff` (§3.1), under the §3 partition rule; anisotropy reported.

**7.3 — Symplectic form** `Ω_AB(𝔅)`: the three first-order channels gated separately (§3).

**7.4 — Tilt statics + susceptibility:** existence, computed parity, leading `V`-dependence, stiffness/anchoring dependence of `p(V)` and `p(ω,V)`.

**7.5 — ⭐ THE COUPLING PACKAGE (the object that can replace `j=sV`).** Not a bare `δS/δΦ` — four declared components (Grok B3), for each mediator `Φ_A ∈ {h, u_T, u_L, wall/χ modes}`:
 a. **domain + surface terms:** `S_body` defined on `Ω_c` with explicit boundary terms on `Σ` and `∂Ω_c` (surface contributions are often the actual exterior source);
 b. **bulk source** `J_A = δS_body/δΦ_A` **and the FULL MATRIX kernel perturbation** `O_AB → O_AB + δO_AB[V, p, s; 𝔅]` — **explicitly including the off-diagonal mediator/mode-mixing entries** (pathA_39 Stage 3's committed mixing is a coupled operator, not per-mediator diagonal perturbations; computing only diagonal `δO_A` would miss it);
 c. **constraint/virtual-work channel** when `𝔅` is E4/E1 (nonholonomic or no-slip forces enter by multipliers, not `δS/δΦ`);
 d. **far-field multipoles of the TOTAL COUPLED linear response** (sources + the full mixed kernel `O_AB+δO_AB` together — the object the two-body computations actually consume), with **`s`-parity channel-by-channel** and `O(V)` structure; G8's ablations act on this total coupled response.
 **Mass/charge split (Grok M4):** every response is decomposed into **drain/mass-sourced vs orientation/charge-sourced** parts (pathA_39 separated `q_M` from charge residues); the mass-only ablation gate G11 enforces it.
 *(pathA_39 declared `j_T=q_A^T V`, `j_L=q_L V`; 7.5 replaces the declarations with computed objects.)*

**7.6 — Intake response** `C_ṁ(𝔅, closure)` (§4), partitioned per §3.

**7.7 — Radiative residues:** per-channel moving-body residues; accelerating-body forms; memory kernels where the local truncation fails.

**Gates (v3; process checks labeled as such — they document honesty, they do not gate physics):**
- **G1** existence + normalizability of the translational zero mode (else the collective frame fails);
- **G2** finite effective coefficients under the **declared IR scheme** (scheme in §7.0; regulator-dependence beyond it = fail);
- **G3** kinetic Hessian positive, or the instability explicitly classified;
- **G4** Berry-curvature result validated by a **control with a known nonzero answer** (e.g. an ansatz carrying imposed winding must yield the known Magnus form through the same pipeline) — distinguishing a derived zero from an ansatz artifact;
- **G5** `R_w` covariance end-to-end **on the declared ambient-symmetric background** (scoped by the §1 postulate; on a one-sided background this gate is replaced by the declared asymmetry map);
- **G6** *(process)* `𝔅`-endpoint sensitivity reported for every output;
- **G7** the forced-`O(V)`-`ṁ` check (S2 self-consistency);
- **G8** coupling ablation **per channel of 7.5a–d**: each claimed coupling vanishes when its originating term is switched off (computed provenance, not grep);
- **G9** flux closure across `∂Ω_c` with a **quantitative residual bound**; any imbalance must be assigned to a *named* OPEN ledger entry (the return path), not free-floating;
- **G10** *(process)* tilt-zero honestly recorded if `p(V)≡0`;
- **G11** **mass-only ablation:** with the orientation (charge) structure removed, the claimed *charge*-channel `O(V)` multipoles of 7.5d must vanish (or the contamination is quantified) — the gravity drain must not masquerade as an EM current.

## 8. What U1 is NOT **[NOT-U1]**

- **Not U2:** U2 = *deriving* `𝔅` (which endpoint — or mixture — the sleeve's actual structure enforces, and where it lives; the ontology allows collar and/or bulk-lattice hosts). U1 declares `𝔅` and computes conditional on it.
- **Not the signs:** no electric sign, no magnetism sign, no current law asserted; downstream consumers use 7.5 + U2.
- **Not a Gauss/`U(1)` rescue:** additive charge/conservation stay non-native; `Z₂` is what the body has.
- **Not a full nonlinear throat sim** (sim-deferred): collective-coordinate/effective-action level, with unavailable profiles declared open, never silently guessed.

## 9. Open seams **[OPEN]**

- The geon profile / mass mechanism (7.0 input); whether an `E/c²`-like inertia relation *emerges*.
- Bulk↔brane return closure + drain-vs-charge energy budget (§4 carries them symbolically).
- **Reconciling the §1 ambient-`R_w` postulate with the one-sided `pathA_29` slab family** (stack/dual-return geometry vs computational device).
- Whether the supersolid-bulk lattice modifies `𝔅` (U2's question).
- Which §7.0 coefficients tie to already-frozen ansatz values vs remain free.

## 10. Fold logs

**v1 → v2 (Codex design-review, 12 items — all folded):** (1 BLOCKER) U1↔U2 circularity → `𝔅` input + conditional outputs; (2 BLOCKER) `ρa³L`/`E/c²` imports → Hessian-of-native-action definition; (3 BLOCKER) sign-split axiom too strong → `R_w` covariance, per-observable parity, drain-rate equality tagged; (4 BLOCKER) no-Magnus insufficient → `Ω_AB` computed, 3 channels; (5 BLOCKER) tilt charge-parity wrong → invariant `p`, parity computed; (6 BLOCKER) nothing replaces `j=sV` → source deliverable; (7) declared-inputs deliverable 7.0; (8) full radiative-mode enumeration; (9) geon in composition + control surface; (10) gates reworked; (11) `L/a` relabeled; (12) missing body properties folded into 7.1/7.4/7.5.

**v2 → v3 (Grok independent review — all folded):**
| # | Severity | Catch | Fold |
|---|---|---|---|
| B1 | BLOCKER | `S_eff` Hessian ill-posed for an open drain + dissipative `𝔅`; flux/inertia double-count; EOM has no non-variational slot | four-channel force balance; `L_eff` scoped on `Ω_c`+surfaces; partition rule; local EOM = declared truncation (§3) |
| B2 | BLOCKER | `R_w` conjugation assumed ambient `w`-symmetry; `pathA_29` slab is one-sided; G5 tautological-or-false | body conjugation (definition) split from ambient symmetry ([POSTULATE] + open seam); G5 scoped (§1, §9, G5) |
| B3 | BLOCKER | `J_A=δS_body/δΦ_A` not closed: undefined domain, missing surface terms, `δO` operator channel, constraint forces | 7.5 = four-component coupling package (a–d); G8 ablates per channel |
| M1 | MAJOR | `𝔅` endpoints not operational; fluid no-slip conflated with the ontology's nonholonomic shear-lock; "no plug dragged" pre-answers E3 | operational endpoint table E1–E5 incl. E4 shear-lock; terminology guard; kinematics language neutralized (§2) |
| M2 | MAJOR | zero-mode quotient unspecified; slab vs translational zero mode conflatable | explicit quotient recipe in 7.1; §0.3 disambiguation |
| M3 | MAJOR | `θ` re-imports the charge-odd convention | `p` used everywhere; signed `θ` banned (§5, §7) |
| M4 | MAJOR | mass-current vs charge-current not gated; drain could masquerade as EM current | 7.5 mass/charge split + gate G11 |
| M5 | MAJOR | G4/G5/G9 vacuous or convention-tests; G6/G10 are process | G4 known-nonzero control; G5 scoped to ambient postulate; G9 quantitative residual + named ledger entry; G6/G10 labeled process (§7) |
| M6 | MAJOR | local `F=MV̇+ΩV` pre-answers radiation structure | declared sub-radiative truncation; memory kernels in 7.7 (§3.4, S5) |
| m1 | MINOR | `S_eff` diverges at steady `V` | `L_eff` wording (§3.1) |
| m2 | MINOR | channel-(i) "expected zero" must stay tested | explicit in §3 gating |
| m3 | MINOR | frame undeclared | medium-rest-frame declaration (§1) |

**v3 → v3.1 (Codex confirm-pass, 2 items — both folded):**
| # | Severity | Catch | Fold |
|---|---|---|---|
| C1 | BLOCKER | E4's nonholonomic reaction had no slot (`F^diss=0` at E4 recreated Grok-B1 for the decisive endpoint) | `F^𝔅 = F^constraint (Lagrange–d'Alembert multipliers; E4 here, E1 placement declared) + F^Rayleigh (E5)`; decomposition carried into §7 (§3.3) |
| C2 | MAJOR | §7.5b operator perturbation written diagonal; pathA_39 Stage-3 mixing is a coupled operator | full matrix kernel `O_AB+δO_AB` incl. off-diagonal mixing; 7.5d + G8 act on the total coupled response (§7.5b/d) |
