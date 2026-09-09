# S11c-d — SHARED PHYSICS (the profile-conditioned transverse↔thickness mixing/leakage for a localized-interface profile)

**S11c-d** is the fourth sub-step of the S11c curved-interface program (`directives/S11c_decisions.md` row `:52`). It
consumes S11c-c2's closed operator and closed off-diagonal kernel and produces the **profile-conditioned linear
mixing** between the uniform transverse (light) sector and the thickness/breathing (self-energy) sector, for an
**explicitly named localized-interface profile class**, together with the two distinct photon-kill channels and the
flux-normalized leakage FORM. This document is the physics authority for the two blind S11c-d engines and their
comparator. Tag prefix `S11CD_`.

The SymPy engine reads the inherited model through `ledger_fold.load_model` over the atomic frozen base
`scripts/S11c_b_exports.py` with the c1 and c2 deltas folded on top (§7), binding only its declared `IMPORT_KEYS`; the
Wolfram engine imports nothing and re-derives every consumed object from the sibling specs
(`S9_export_chain_rebuild_directive.md:16-18` is the only cross-engine control). Blindness is the control: an
agreement is independent construction, not a copy.

⭐ This is an **orchestrator-written physics spec** (the physics authority both blind engines read). Per `CLAUDE.md`
G1/G2 it is physics-bearing and gets **two legs — Codex `gpt-5.6-sol` xhigh + Grok — reviewed UNTIL CLEAR** (spec
row, ⛔ not the decision-list one-pass): leg→fold→leg until nothing outstanding changes what is computed or may be
claimed; both reports before any commit; the reviewed baseline is preserved before any repair overwrites it. The
**build directive** that follows this spec gets its **own** two decision legs before any builder (the TRIGGER, `G2`).
⚠ **Spec v1.**

⭐ **The profile-class + regime decision (§1c–§1d) was settled by a three-way physics consult** (orchestrator + `gpt-6-astra`
xhigh + `grok-4.6`, `directives/_legs/S11c_d_profile_class_consult.md`) and the user's approval: **localized interface,
Born in contrast `η` with sharpness `σ_W` and kinematics `kL` kept live** — ⛔ **not "weak-gradient."** The two
engines' load-bearing corrections (the grade-conflation that "weak-gradient" hides, and the strong-edge bridge as an
unresolved obligation rather than a reduction) are folded into §1d, §3d, §5.

---

## 0 · Scope

**In scope.**
1. **NAME the profile class** = a **localized interface** in the background thickness field `W₀(x)` (§1c): smooth,
   asymptotically constant, with a localized gradient of finite integrated jump `∫W₀′ dx = W₊ − W₋ ≠ 0` — ⛔ **not** a
   defect bump (`∫W₀′ dx = 0`, a distinct object, §1c/§5c).
2. The **profile-conditioned transverse↔thickness mixing response** (§3a): the linear mixing driven by `∇μ_R ≠ 0`
   (`∝ k·a`), organized as **one insertion of `S11CC2_CLOSED_COUPLING_KERNEL` between the diagonal sector responses
   (resolvents / Green operators) of `S11CC2_CLOSED_SLAB_OPERATOR`** (distorted-wave Born, §2).
3. The **two DISTINCT photon-kill channels** (`N13`, §3b): continuum conversion (into the thickness continuum / bulk
   escape) **and** capture into a **bound** thickness/breathing pole; and the confinement question (`N10`) as a
   computed object.
4. The **order bookkeeping** (`N12`, §3c): converted amplitude `O(εη)`, absolute converted flux `O(ε²η²)`, **and** the
   flux-normalized dimensionless conversion fraction `O(η²)` — emit **both** the absolute and fractional labels.
5. The **N6 independent shape/coordinate-route control + one-sided corruption** (§5a), the profile-**form** ablation
   (§5c), and the flux-normalized falsification **FORM** (§5d; the numeric bound withheld, `N7`).

**Out of scope (named, not solved).**
- The **order-unity slit edge** (finite-`ΔW`, non-perturbative) — a **downstream obligation named in §3d**, ⛔ not a
  mechanical reduction of the weak result: a nonzero Born coefficient gives **no** general lower bound on order-unity
  conversion (§3d counterexample). Its construction (piecewise-uniform jump / matched interior-exterior / an effective
  interface response with explicit undetermined parameters) is **new** and belongs to the S11c-e strong-edge stage,
  where its weak limit must reproduce S11c-d.
- The **falsification magnitude** — needs the throat interior `R1` (`V3_STEP_PLAN.md:1179`); only the FORM is
  computable now (`N7`). The withheld `O(1)`/grating reductio is diffed **orchestrator-side**, ⛔ never a builder
  target.
- A **global dispersion `ω(k)`** for generic `W₀(x)` (`N5`/`N10`, forbidden); the **nonlinear-light program** (`N10`);
  the **periodic→Bloch** and **slowly-varying→WKB** classes (§1d rejects them for this endgame; a periodic profile is
  a later evaluation of the **same** kernel, not an S11c-d class).
- The **full cross-engine self-energy operand residual** and c1's four giant families (≥64 GB,
  `DEFERRED_HEAVY_RUNS.md`). S11c-d must be constructible and cross-engine-testable on this box for its own mixing
  object, and **name — not silently absorb** — anything it cannot close here (§1b).

---

## 1 · Complete inherited setup — SUPPLIED and unfalsifiable in this build

Everything in §1 is an input. The mixing response, the two channels, the leakage FORM, and every control disposition
of §§3–5 are **outputs**; ⛔ none is stated here.

### 1a · Inheritance and the consumed c2 exports

The DOFs, sector split, background ansatz, `(ε,η,σ_W)` power counting (`N12`), and admissibility are exactly S11c-a
§§1–2 / S11c-b §§1–2 / S11c-c2 §1, inherited by pointer. S11c-d consumes one already-built, per-engine-reviewed model,
**S11c-c2** (`scripts/S11c_c2_exports.py`; step record `steps/S11c_c2_self_energy_fold.md`; disposition
`_measurements/S11c_c2_N6_reconcile_disposition.md`):

- **`S11CC2_CLOSED_SLAB_OPERATOR`** — the closure-modified variable-coefficient slab operator over `{u,θ,e_W}`, per
  `(anchoring α, density ρ)`, assembled two-face, with `S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS` and
  `S11CC2_CLOSED_SLAB_OPERATOR_PARITY_BLOCKS`. Its **diagonal blocks** supply the local sector spectra / resolvents
  (the transverse response `G_T`, the thickness response `G_H` and its poles).
- **`S11CC2_CLOSED_COUPLING_KERNEL`** — the re-extracted **off-diagonal** transverse↔`{θ,e_W,u_L}` block(s), with
  `S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS`. This is the **mixing vertex** S11c-d inserts.
- the self-energy increment `S11CC2_SELF_ENERGY_INCREMENT`, its operands, and the six §3d re-adjudication objects — as
  provenance for the disposition of §1b, ⛔ not re-opened here.
- inherited by pointer through c2: the S11c-a T-a..T-i face substrate, `S11CB_MU_THETA_OPERATOR`,
  `background_density_map`, the `Λ_{A,V,X}` closure channels, the two-momentum DtN kernel `dtn_kernel`, and the
  constants/profile carriers (`W_0`, `W_bg`, `w1_profile`, `L_W`, `sigma_W`, `eta_bg`, `mu_R`, `rho_m`, `rho_br`, …).

### 1b · What is per-engine-SOUND vs cross-engine-UNCLOSED in the c2 import — SUPPLIED HONESTLY (rule 6/16)

⭐ **This is the load-bearing honesty section: S11c-d's physics IS the gradient-driven off-diagonal kernel, so c2's
carried cross-engine operand DEBT is MATERIAL to this consumer — ⛔ NOT dismissible on covariance alone**
(`steps/S11c_c2_self_energy_fold.md:214-219`).

⭐ **PER-ENGINE SOUND (SymPy, 2-leg):** the self-energy fold wiring + A/C/D1–D6 and the emitted
**`S11CC2_CLOSED_SLAB_OPERATOR` + `S11CC2_CLOSED_COUPLING_KERNEL` + increment VALUES**. These are the operands
S11c-d consumes.

⭐ **CROSS-ENGINE, dual-engine confirmed (the N6 thread only, on this box):** **operator covariance (Reading B)** — the
material builder implements the declared frame change `Φ` (`R_cov` no-nonzero in both engines). ⚠ The matched
covariance-channel cross-engine zeros are **`(0)−(0)`** — a dual-engine confirmation of the **vanishing** statement,
⛔ **NOT operand agreement.** Preserved together: `R_N6 = I_E − I_{M→E}` **nonzero (18/288)** AND `R_cov` **no-nonzero**
(consistent under Reading B).

⛔ **CROSS-ENGINE UNCLOSED — S11c-d must NOT treat these as closed (a supplied, unfalsifiable-in-this-build premise it
names honestly, `M2`/rule 16):**
- **the cross-engine OPERAND DEBT** — the surfaced blind-WL-vs-imported **carrier (40)**, constitutive **source (76)**,
  and **Φ (18)** residuals are **UNADJUDICATED**, and the leftover SHAPE was **not inspected**. The v3 collapse
  instrument was CLOSED **NOT-SOUND** (the post-EL graded coefficient table is the wrong object,
  `_measurements/S11c_c2_N6_reconcile_disposition.md`); the only sound reconcile instrument is **upstream of the EL
  differentiation**. ⛔ **Do NOT let "representational-difference-UNADJUDICATED" become "known to be just thickness"**:
  the alternatives include a genuine constitutive-convention mismatch OR an implementation error. Because S11c-d's
  mixing rides on `S11CC2_CLOSED_COUPLING_KERNEL`, this DEBT is a **live premise on the object S11c-d builds** —
  S11c-d **names** it as unclosed and ⛔ does not pre-adjudicate it.
- **the two S11c-b sign conventions** that multiply the substituted `δp_s` slots and do **not** cancel from c2's
  residual (`steps/S11c_b_variable_coefficient_operator.md:112-115`): the **face-generalized-force** convention (PY
  `+diff` vs WL `−linearVirtualVariation`) and the **#90 closure-fold** sign. (The kinetic `−K/+K` convention is a bulk
  term independent of the response slots.)
- **the six §3d re-adjudications** carried per-engine SymPy, of which the c1-UNDECIDED imports (background density,
  `t_s` scalar-vs-4-vector, DtN whole-form, flat-resolvent leg-labeling, ENERGY) stay **cross-engine-UNDECIDED**.
- **the 3 N6 premise caveats** (Reading B does not close them): (1) is `Φ` itself physically correct (not merely
  reproduced); (2) does the face velocity `V` transform correctly (`V_E≡V_M` is builder agreement); (3)
  extracted-block / omitted-block leakage.
- **F and G are WITHDRAWN interpretations** — the uniform-limit decoupling (F) and directionality (G) rest on the
  retired `verify_F`/`verify_EG` instruments; numeric-probe re-grounding is a standing OWED debt, **PAUSED
  INDEFINITELY / non-blocking**. ⛔ Do not resurface it as a BLOCKER; ⛔ but pausing did not discharge it — the F/G
  conclusions do **not** stand, only the increment VALUES do.

⛔ Folding any of the above to force cross-engine closure is the exact defect this rebuild exists to catch (rule 1/6).

### 1c · The background profile class — SUPPLIED framing (`N5`/`N12`/`N14`/`N15`)

⭐ **The named class is a LOCALIZED INTERFACE.** The background thickness profile is

```text
W₀(x) = W̄₀·[ 1 + η·f(x/L_W) ] ,     f asymptotically constant,  f(−∞) ≠ f(+∞)  ⇒  ∫ W₀′(x) dx = W₊ − W₋ ≠ 0 ,
σ_W = η·W̄₀/L_W  (the first-jet / sharpness grade, kept live) ,   L_W independent of η .
```

- ⭐ **An interface is NOT a bump.** A single edge has a localized **gradient** but **different asymptotic
  backgrounds** (`∫W₀′ ≠ 0`); a defect bump returns to its original background (`∫W₀′ = 0`). Their low-momentum-transfer
  coupling content can differ (§5c). ⛔ A bump must **not** silently stand in for a single edge. A smooth finite barrier
  with two transitions is the corresponding localized model for **two** edges.
- A specific shape (e.g. `f(ξ) = (1+tanh ξ)/2`) is a **representative**, ⛔ not the class — keep `f` general within the
  localized-interface class (or name a representative as a representative). The mixing vertex depends on `∇w₁`.
- **Admissibility (`N12`).** Name which quantities vary (`W₀(x)`, `μ_R(x)`, `ρ_br⁰(x)`), the profile's **anchoring**
  (material-advected vs lab/Eulerian-held — inherited as the anchoring `α`, `N4`), and the **stationary equations or
  the named force** that holds it — an inadmissible background silently sources spurious coupling.
- **Names (`N14`).** Every spatially-varying field/kernel/observable gets a **fresh** injective standard name; ⛔ never
  reuse an imported S11b/c1/c2 key (`W_0`, `e_W`, `rho_br`, `v_0`, `slab_operator`, `coupling_kernel`, …) for a varying
  object — `F9`'s object comparison proves a false equal, and for `rho_br` it is the `∇Σ_E⁰=0` freeze that drops the
  advective channel.
- **Invariants (`N15`).** Inherit S11b's invariants with variable coefficients; **emit any new gradient-of-background
  invariants as RESULTS** (new constants if they appear) — e.g. the edge form factor / integrated-jump moment. ⛔
  Neither smuggle them in by an `W₀→W₀(x)` substitution into the uniform energy, nor forbid them.

### 1d · The regime — SUPPLIED (Born in CONTRAST, sharpness LIVE)

⭐ **Born in the contrast `η`; the sharpness `σ_W` and the kinematic `kL_W` are INDEPENDENT live grades.** The three
scales are distinct and must not be conflated:

```text
η    = contrast          (Born / scatterer strength) ,
σ_W  = η·W̄₀/L_W = first-jet / sharpness   (the adiabaticity axis; the WKB axis) ,
kL_W = kinematic         (is the edge sharp on the WAVE) .
```

⛔⛔ **"Weak-GRADIENT" is the wrong name and is FORBIDDEN in the build framing** — it identifies gradient with
contrast. A builder handed "weak-gradient" will Taylor-expand the profile in `σ_W` or send `L_W → ∞`, which **kills the
Fourier form factor of a localized edge and slides the object toward WKB** — the one class that does not match an edge.
⭐ Say instead: **localized interface, Born in contrast `η`, `σ_W` live** — ⛔ do not expand the profile in
derivatives, ⛔ do not send `L_W → ∞`. Weak **contrast** does not require WKB: provided the frozen approximation holds,
`kL_W` need not be large.

⚠ The retained rectangle is the inherited `(η^{≤1}, σ_W^{≤1})` (`N12`). ⛔ **Do NOT set `η → O(1)` inside the c2
operators**: `S11CC2_CLOSED_SLAB_OPERATOR` / `S11CC2_CLOSED_COUPLING_KERNEL` are **first-shape-order** truncations, and
`η = O(1)` uses an incomplete operator exactly where the missing `O(η²)` pieces compete — that is a bookkeeping
violation, not "going non-perturbative." Keeping `η` a formal symbol does **not** restore the missing error control
(§3d).

---

## 2 · The mixing object and the distorted-wave-Born organization — SUPPLIED framing

The object is the **linear mixing** between the uniform transverse sector and the thickness sector, supported only
where `∇μ_R ≠ 0` (the localized gradient). With a **uniform** background the two sectors decouple identically (S11b);
mixing is sourced by the gradient (`∝ k·a`). The natural organization is **distorted-wave Born**: the diagonal blocks
of `S11CC2_CLOSED_SLAB_OPERATOR` supply the incoming/outgoing sector fields, and the intersector conversion is one
insertion of `S11CC2_CLOSED_COUPLING_KERNEL`:

```text
mixing response  =  [  − G_H⁺ · K_HT · ψ_T^{in,D}  ]_{εη}          (one kernel insertion; the reverse uses the other block)
where   G_H⁺   = outgoing Green operator / resolvent of the thickness-diagonal block of S11CC2_CLOSED_SLAB_OPERATOR ,
        K_HT   = the transverse→thickness block of S11CC2_CLOSED_COUPLING_KERNEL (the off-diagonal vertex) ,
        ψ_T^{in,D} = the incoming transverse solution of the uncoupled DIAGONAL problem .
```

⭐ The **scattering labels are frequency and asymptotic incoming/outgoing channels** — these do **not** require, and
do **not** imply, a global `ω(k)` for the inhomogeneous slab (`N5`). Two framing obligations:
- ⭐ **Re-expand to the retained `η` order.** Keeping the profile-dependent diagonal solutions unexpanded is a useful
  partial resummation, but its extra powers are **not** a complete higher-order prediction — retain the specified
  `(η^{≤1}, σ_W^{≤1})` grade.
- ⭐ **Specify a REGULAR scattering domain.** Exclude unresolved threshold/resonance enhancements and interaction
  lengths that make repeated conversion appreciable — a weak pointwise gradient can accumulate a large conversion over
  a sufficiently long coherent region.
- ⚠ **Sub-choice deferred to the build directive** (they coincide at leading `O(η²)`): plain Born on the **uniform**
  sector modes vs distorted-wave Born on the **profile-dressed** diagonal Green operator `G_H⁺(ω;W₀)`. The leading
  conversion is the same; the build directive fixes which, and re-expands to the retained grade either way.

---

## 3 · The construction (OUTPUTS)

Every object below is computed for both anchorings `α` and both density representatives `ρ`; it carries its computed
`(ε,η,σ_W)` multigrade and restored `[L,T,M]` dimension, and states no component value, sign, order, parity, or grade
in this document.

### 3a · The profile-conditioned mixing response (continuum channel)

Emit the mixing response of §2 for the localized-interface profile, both intersector directions (both exported
off-diagonal blocks), per `(α,ρ)`. Its outgoing-channel projections are the **conversion / scattering amplitudes**.

```text
⇒ S11CD_MIXING_RESPONSE (per (α,ρ), both directions) , S11CD_MIXING_RESPONSE_TERM_ORIGINS ,
  S11CD_CONVERSION_AMPLITUDE  (outgoing-channel projections) .
```

### 3b · The two DISTINCT photon-kill channels (`N13`) + confinement

`N13`: "confinement of light" = **survival of the transverse polarization channel**. Conversion into a **bound**
breathing/thickness mode kills the photon **exactly as** bulk radiation does — the two are **distinct emitted
objects**, ⛔ not one "energy stays in the slab" statement.

- **(i) continuum conversion** — transverse → the thickness **continuum** / bulk escape (the §3a amplitude projected on
  the radiating/continuum channel).
- **(ii) bound-mode capture** — transverse → a **bound** thickness/breathing pole (a resolvent **pole** of the
  thickness-diagonal block of `S11CC2_CLOSED_SLAB_OPERATOR`; a weak attractive well binds a mode in 1D, so this is a
  real photon-kill channel at **small** `η`, and is **not** a Bloch band). Emit the pole/residue object and the
  transverse→bound coupling.
- The **confinement question** (`N10`): whether transverse-channel survival is **unconditional** — emit the computed
  object, ⛔ not a claim.

```text
⇒ S11CD_CONTINUUM_CONVERSION , S11CD_BOUND_MODE_CAPTURE (pole + residue + coupling) , S11CD_CONFINEMENT_CONDITION .
```

### 3c · The leakage order bookkeeping — the export representation (`N12`; the ε²-cancellation correction)

With incident field amplitude `ε` and contrast `η`, in the **linear** theory:

```text
converted amplitude       ψ_H      = O(ε η)   ,
absolute converted flux    J_conv   = O(ε² η²) ,
incident flux              J_in     = O(ε²)    ,
dimensionless conversion FRACTION   C = J_conv / J_in = O(ε⁰ η²) = O(η²)     (the incident ε² CANCELS) .
```

⭐ **Emit BOTH labels — absolute (`O(ε²η²)`) and fractional (`O(η²)`) — and keep both visible.** A normalized /
per-photon conversion rate in a linear theory **cannot** retain the incident `ε²`; the flux-**normalized** dimensionless
conversion FORM (the S11c-e observable, `N7`) is therefore `O(η²)`, while `N12`'s "`O(ε²η²)` leakage rate" is the
**absolute** converted power. ⚠ Naming both prevents an accidental intensity dependence, and prevents a mis-ordered
`O(ε²η²)`/`O(η²)` term from being read as the excluded nonlinear-light program (`N10`/`N12`). The **magnitude** is
`R1`-blocked (out of scope); only the FORM is computed here.

```text
⇒ S11CD_CONVERSION_FRACTION_FORM (flux-normalized, O(η²)) , S11CD_CONVERSION_POWER (absolute, O(ε²η²)) .
```

### 3d · The strong-edge bridge — NAMED as a downstream obligation, ⛔ NOT solved, ⛔ NOT a reduction

⭐⭐ **S11c-d establishes the WEAK matching coefficient only; it does NOT establish the order-unity edge form.** A lab
slit edge is **localized + order-unity contrast**. S11c-d honestly delivers `F′(0)` (the `O(η)` slope of the
conversion); the lab bounds `F(1)`.

⛔⛔ **A nonzero Born coefficient supplies NO general positive lower bound on strong-edge conversion.** Illustrative
counterexample (⛔ not a model of the slab): a lossless two-mode coupler with dimensionless integrated coupling `G`,

```text
A_H = −i ε sin(η G) ,     C = sin²(η G) = η²G² + O(η⁴) ,     C → 0  at finite η G .
```

The converted amplitude starts at `O(εη)` and the flux at `O(ε²η²)` **exactly as required**, yet the exact conversion
can return to **zero** at finite coupling. ⛔ **Do NOT evaluate the Born coefficient at `η=1` and compare that number
to the lab** — that identifies `F′(0)` with `F(1)`, the invalid extrapolation `N7` names.

⭐ **The downstream obligation, stated now (owned by the S11c-e strong-edge stage):** a **justified finite-contrast
response** — a piecewise-uniform (finite-`ΔW`) matching, matched interior/exterior solutions, or an effective interface
response whose undetermined parameters remain **explicit** — whose **weak limit reproduces S11c-d**. This is a **NEW
construction** (each uniform side is `η`-exact; the finite-`ΔW` matching is not a first-jet kernel). ⚠ Born can miss
repeated conversion/reconversion, diagonal reflection, resonance shifts, and altered channel availability — these
change **frequency and angular dependence**, not just magnitude. ⇒ if the finite-contrast response cannot be
established in scope, the honest S11c-e outcome is a **conditional constraint on edge-response parameters** (or a
deferred numerical exclusion), ⛔ **not** a shape-independent exclusion — and the unknown interior coupling need **not**
factor as `C_edge(ω,ϑ) = C_interior·F(ω,ϑ)` (repeated scattering can put it inside resonance denominators).

```text
⇒ S11CD_WEAK_MATCHING_COEFFICIENT (F′(0)) , S11CD_STRONG_EDGE_OBLIGATION (named premise, ⛔ not solved here) .
```

---

## 4 · Objects to compute and emit

Per anchoring `α∈{L,M}` and density representative `ρ∈{ρ_4D,ρ_br}`, multigraded and dimensioned:

- The **profile-conditioned mixing response** + its term-origins + the conversion amplitude — §3a.
- The **two photon-kill channels** (continuum conversion; bound-mode pole + residue + coupling) + the confinement
  condition — §3b.
- The **leakage bookkeeping** — the flux-normalized conversion FRACTION FORM (`O(η²)`) and the absolute converted power
  (`O(ε²η²)`) — §3c.
- The **weak matching coefficient** `F′(0)` and the named strong-edge obligation — §3d.
- The **control outputs** of §5, each emitted as the object and its literal residual.
- Any **new gradient-of-background invariant** surfaced (`N15`, e.g. the edge form factor).

Every result carries its `(ε,η,σ_W)` order (`N12`) and its restored `[L,T,M]` dimension. ⛔ No result is reported
without both.

---

## 5 · Independent routes and controls

⭐ Every control re-enters the chain **at the ACTION / the imported operands**, ⛔ never at a result. Each emits the
object and its literal residual; ⛔ none asserts a target value. A **coefficient** rescale tests arithmetic; only a
**form** change tests physics.

### 5a · The N6 control — the independent shape/coordinate route + one-sided corruption (`decisions:94-104`)

⭐ **The genuine control (rule 14), ⛔ NOT the uniform limit.** Derive the off-diagonal mixing two independent ways and
compare, then corrupt one route only:

```text
route 1 (level-set / graph):    derive the off-diagonal mixing by direct level-set / graph linearization of the
                                 localized-interface faces ;
route 2 (flattened material):   derive it AGAIN after flattening the faces into material coordinates, then transform
                                 Eulerian ↔ material EXACTLY into the common Eulerian face basis ;
S11CD_REP_INVARIANCE_RESIDUAL[α,ρ] = route1 − route2      (the representation-invariance measurement) .
```

`N6` is the physics requirement that these are the **same operator in two representations** — the uncorrupted residual
is the measurement, its **computed value is the finding**, ⛔ no target value is supplied, and the diff is adjudicated
on our side (⛔ never a builder exit condition). Then the **one-sided independence corruption** (still at fixed `α,ρ`):
mutate **one route only at its source** and require a **nonzero** residual while the **uncorrupted** route is unmoved —

- **(i) tilt probe (`N3`):** reverse **one** face's first-jet slope term in the outward normal `n̂_s` on **one** route
  (through that route's own carrier factory, ⛔ not by editing the imported operator, ⛔ not by altering only the DtN
  kernel jet).
- **(ii) N4 advection probe:** omit / flip **one** route's advective-density term (`u·∇ρ₄/ρ₄`); ⚠ this term is
  **structurally absent for `RHO4_CONSTANT`** (`∇ρ₄=0`) and present for `RHOBR_CONSTANT` — the live probe is
  `RHOBR_CONSTANT`; for `RHO4_CONSTANT` emit the **computed absence**, ⛔ never an `A−A`.

⚠ There are **≥2 same-order channels** (tilt `N3`; advection `N4`); the one-sided corruption is the **independence test
between them**, ⛔ so "the gradient channel" (singular) is wrong. ⛔ `∇W₀→0` / `η→0` is **NOT** an accepted corruption
(it is the vacuous uniform limit renamed, `N6`). ⛔ Corrupting one **anchoring** is not this test (it only shows two
distinct physical setups differ). The `M→E` map lives **inside** route 2's native builders (covector inverse-transpose)
— ⛔ no separate `T` on the differenced object; the field redefinition `Δρ` relates two descriptions of **one**
perturbation, ⛔ never the two anchorings.

```text
⇒ S11CD_REP_INVARIANCE_{ROUTE1,ROUTE2,RESIDUAL}[α,ρ] , S11CD_CONTROL_INDEPENDENCE_{BASE,CORRUPTED,RESIDUAL}[α,ρ,probe] .
```

### 5b · The uniform limit — REGRESSION smoke-test only

`W₀→W̄₀` (`η→0`): the off-diagonal mixing must vanish (S11b decoupling). ⭐ **Secondary smoke test only** — it cannot
see the coefficient, sign, or parity of the gradient coupling (S11b: coupling identically zero); it is a useful check
for a forbidden gradient-**independent** term. ⛔ It is **not** the N6 control (§5a).

### 5c · The profile-FORM ablation + the edge-vs-bump discriminant

- **Form ablation.** Perturb the **FORM** of `f(ξ)` (the localized-gradient shape) and require the mixing / leakage to
  **move** — ⛔ a coefficient rescale of `η` is insufficient; only a form change tests the coupling.
- **Edge-vs-bump discriminant.** Emit the mixing object for the **interface** (`∫W₀′ ≠ 0`) and — as the discriminant —
  for a **bump** (`∫W₀′ = 0`), and their difference: the low-momentum-transfer content differs, so a bump must **not**
  stand in for the edge. ⛔ Do not freeze a specific representative shape as "the slit."

```text
⇒ S11CD_PROFILE_FORM_ABLATION , S11CD_EDGE_MINUS_BUMP .
```

### 5d · The falsification FORM control (`N7`)

Emit the **flux-normalized dimensionless conversion FORM** (§3c, `O(η²)`; `∝ k·a`, supported where `∇μ_R ≠ 0`, carrying
the `σ_W` form factor). ⛔⛔ The `O(1)`-fraction / diffraction-grating reductio and the withheld numeric lab bound are
**orchestrator-side**, ⛔ never in the builder-facing acceptance text (the builder iterates toward any target it can
see); the magnitude is `R1`-blocked. ⚠ A slit edge is an order-unity localized gradient — the FORM here is the
**weak-contrast** object (§3d states which), ⛔ not a non-perturbative lab number.

---

## 6 · Method, dimensions, and script obligations

- **Method.** Balance laws + the binding material virtual-displacement rule + variational derivatives with held-fixed
  fields named + prescribed external virtual work (S11b), ⛔ never an irreversible response kernel in an ordinary
  action. The diagonal responses and the off-diagonal vertex are the c2 exports consumed verbatim; the mixing is one
  kernel insertion (§2), re-expanded to the retained grade.
- **Dimensions.** Restore `[L,T,M]` on every emitted object, dimensional consistency able-to-fail
  ([[feedback_dimensional_consistency_check]]); `(ε,η,σ_W)` multigrade on every object (`N12`).
- **Rest-frame limit.** Inherit `N11a` inert; S11c-d constructs **no** convective operator. Every result inherits the
  c1/S11b smallness domain (`|q_out·v_bulk_normal_0/ω|≪1` + boundary-layer/subsonic), ⛔ never aliasing
  `v_bulk_normal_0` to `v_0` (`N14`/`N11`).
- **Script obligations.** The three build-skill clauses bind the build directive (`.claude/skills/build/SKILL.md`): a
  script PRINTS computed objects and never states conclusions; PRINT the residual, do not assert it; interpretation is
  the step record. ⛔ No hand-typed CAS object standing in for a computed one; every control re-enters at the
  ACTION / imported operands. ⛔ No tautological residual (rule 2 corollary 3): the §3a/§3c export representations are
  **not** checks; the §5 residuals are emitted with both operands, and a two-route residual is emitted only where an
  **independent** second route exists.
- **Serialize CAS jobs; watch RSS.** c2's self-energy `.out` was ~499 MB and the full cross-engine residual is the
  ≥64 GB work; S11c-d's mixing is a **derived insertion** on the closed operator — measure the process that runs,
  defer heavy controls in-band→out-of-band (`DEFERRED_HEAVY_RUNS.md`), ⛔ never two memory-heavy CAS jobs concurrently.
  Detached launch (harness reaps `run_in_background`). Mathematica: 2-seat licence, `--sandbox danger-full-access`,
  serialize dual ablations.

---

## 7 · Names, F9 reservations, chain output, and export schema

**F9 / `N14` reservations.** Every new object gets a **fresh** injective `mechanical_lower_camel` name; ⛔ never reuse
an imported S11c-c2/c1/S11c-b/S11b key (`closed_slab_operator`, `closed_coupling_kernel`, `slab_operator`,
`coupling_kernel`, `dtn_kernel`, `mu_theta_operator`, `w1_profile`, `W_0`, `e_W`, `rho_br`, `v_0`, …) for a new S11c-d
object.

**Chain output (`N1`/`N8`; topology = the two-leg-gated `directives/export_ledger_bind_closure_design.md` §D1–§D3).**
The SymPy engine reads the inherited model via the **positional** `load_model` call over the frozen base
`scripts/S11c_b_exports.py` with the c1 and c2 deltas folded on top (signature `load_model(base_path, *delta_paths)`,
`scripts/ledger_fold.py:102`; ⛔ NOT keyword form), binding only its declared `IMPORT_KEYS`, and writes
`scripts/S11c_d_exports.py` as its **own-rows delta** (§D2, ⛔ not the accumulated whole-model file). ⛔ The exact
`IMPORT_KEYS` **root set** (minimal roots whose recursive closure covers the §1a consume-set) is fixed at the **build
directive** against the real export files, ⛔ not enumerated-then-frozen here; its two decision legs verify it, and that
the guard (`check_consumer`/`assert_lookups_equal_manifest`/`assert_delta_is_minimal`) passes on the fold — ⚠ noting
the guard passes on **key existence**, so it will **not** catch a wrong-provenance binding; that is the directive's +
legs' responsibility. `BUILD_INPUT_DIGESTS` pins, per §D3, `{this sub-step's SymPy audit, scripts/S11c_b_exports.py,
scripts/S11c_c1_exports.py, scripts/S11c_c2_exports.py, this spec, scripts/ledger_fold.py}`. ⛔ Never `git add -f` a big
`.out`; ⛔ never annex an `*_exports.py`.

**The comparator (`N8`, frozen `T7` contract).** The S11c-d comparator joins the two blind engines' emitted objects by
name, pairs residual operands, is three-valued, rejects a native boolean, and PRINTS/decides nothing (rule 2). ⚠ Its
load-bearing residual is on the **mixing response / conversion amplitude** (§3a) — which **rides on the carried
cross-engine operand DEBT** (§1b); the comparator **SURFACES** the DEBT and the §1b representation questions (the
staged representational bridge, [[feedback_reconcile_representational_bridge]], ⛔ never a blanket collapse), ⛔ does
not pre-adjudicate them. The full per-object symbolic residual + c1's four giants remain deferred
(`DEFERRED_HEAVY_RUNS.md`); S11c-d names, ⛔ does not pre-adjudicate, whatever it cannot close on this box.

**The blind Wolfram engine** re-derives the §§1–2 supplied inputs, the S11c-a face substrate, the S11c-b slab-operator
and c2 closed-operator/kernel rows it consumes, and the localized-interface mixing — importing nothing (the only
cross-engine control). ⛔ The denylist stays cut (`N9`/rule 12); blindness is enforced by absence.

---

## 8 · Supplied versus computed; builder report

**SUPPLIED (unfalsifiable in this build):** all of §1 (the c2 exports and their per-engine-SOUND vs
cross-engine-UNCLOSED disposition — the operand DEBT, the two S11c-b signs, the six §3d re-adjudications, the 3 N6
premise caveats, the withdrawn F/G), the §1c localized-interface class and its admissibility, the §1d regime
(Born-in-contrast, `σ_W`/`kL` live), the §2 distorted-wave-Born organization, `N11a`, `N12`, `N13`.

**COMPUTED (outputs, ⛔ none stated here):** the profile-conditioned mixing response and conversion amplitude (§3a);
the two photon-kill channels and the confinement condition (§3b); the flux-normalized conversion FRACTION FORM and the
absolute converted power (§3c); the weak matching coefficient (§3d); every control residual (§5); any new
gradient-of-background invariant (`N15`); the `(ε,η,σ_W)` orders and `[L,T,M]` dimensions.

**Builder report.** The build directive states, per emitted object, which line computed it (`.claude/skills/build`);
declares the profile class (§1c), the regime grades (§1d), the insertion organization (§2), and the leakage
bookkeeping (§3c) it implemented; and reports the literal residuals of §5 — ⛔ never a prose conclusion. The
disposition of every §5 residual and the strong-edge obligation (§3d) is read on **our** side, in the step record,
⛔ not asserted by the script (rule 5).
