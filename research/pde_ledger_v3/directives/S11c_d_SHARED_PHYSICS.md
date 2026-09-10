# S11c-d — SHARED PHYSICS (the profile-conditioned transverse↔thickness mixing/leakage for a localized-interface profile)

**S11c-d** is the S11c-d sub-step of the S11c curved-interface program — the **fifth build unit** after the c1/c2
split (`directives/S11c_decisions.md` row `:52`). It consumes S11c-c2's closed operator and closed off-diagonal
kernel and produces the **profile-conditioned linear mixing** between the uniform transverse (light) sector and the
thickness/breathing (self-energy) sector, for an **explicitly named localized-interface profile class**, together with
the two distinct photon-kill channels and the flux-normalized leakage FORM. This document is the physics authority for
the two blind S11c-d engines and their comparator. Tag prefix `S11CD_`.

The SymPy engine reads the inherited model through `ledger_fold.load_model` over the atomic frozen base
`scripts/S11c_b_exports.py` with the c1 and c2 deltas folded on top (§7), binding only its declared `IMPORT_KEYS`; the
Wolfram engine imports nothing and re-derives every consumed object from the sibling specs
(`S9_export_chain_rebuild_directive.md:16-18` is the only cross-engine control). Blindness is the control: an
agreement is independent construction, not a copy.

⭐ This is an **orchestrator-written physics spec** (the physics authority both blind engines read). Per `CLAUDE.md`
G1/G2 it is physics-bearing and gets **two legs — Codex `gpt-5.6-sol` xhigh + Grok — reviewed UNTIL CLEAR** (spec
row, ⛔ not the decision-list one-pass); both reports before any commit; the reviewed baseline is preserved before any
repair overwrites it. The **build directive** that follows this spec gets its **own** two decision legs before any
builder (the TRIGGER, `G2`). ⚠ **Spec v3** — folded from two two-leg gates (`_legs/S11c_d_shared_physics_review_*`,
both rounds NOT-SOUND, convergent). Round-1 fixed: nonexistent export rows, a reinstated withdrawn F, an `η`/`σ_W`
collapse, the reserved-key thickness name, undefined interface scattering, a 1D-well over-claim, an `F′(0)`
fraction-slope error. Round-2 fixed: the leftover export-key casing + increment-in-consume-set; an over-constrained
`λ`; an unnamed interface axis + backwards WKB direction; the thickness profile not itself constrained to an
interface; "supported only where `∇μ_R≠0`" (which drops tilt/advection and re-implies F); bound capture folded into a
continuum flux; a class-wide bound-existence claim; a non-canonical S-matrix; and — the deepest — a `route1 − Φ(route2)`
N6 that applies the constitutive map `Φ` to scattering amplitudes/modes/flux (the c2 type error relocated) and has no
constructible mutation site.

⭐ **The profile-class + regime decision (§1c–§1d) was settled by a three-way physics consult** (orchestrator + `gpt-6-astra`
xhigh + `grok-4.6`, `_legs/S11c_d_profile_class_consult{,_astra,_grok}.md`) and the user's approval: **localized
interface, Born in contrast `η` with sharpness `σ_W` and kinematics `Q_nL_W` kept live** (§1d).

---

## 0 · Scope

**In scope.**
1. **NAME the profile class** = a **localized interface** in the inherited background profiles (§1c): smooth,
   asymptotically constant, with a finite integrated jump — required on **both** the thickness `w₁` (`Δw₁ ≠ 0`, so the
   asymptotes `W₋ ≠ W₊`, §2) and the modulus `m₁` (`Δm₁ ≠ 0`, the mixer), kept **independent**; ⛔ **not** a defect
   bump (zero jump, §5c).
2. The **profile-conditioned transverse↔thickness mixing response** (§3a): the linear mixing (the **full** imported
   off-diagonal vertex — tilt `∇w₁`, modulus-gradient `∇m₁`, N4 advection), built on the **two-asymptote distorted
   basis** of `s11cc2ClosedSlabOperator` (its `W₋`- and `W₊`-asymptotic diagonal operators) with one insertion of
   `s11cc2ClosedCouplingKernel` (§2).
3. The **two DISTINCT photon-kill channels** (`N13`, §3b): continuum conversion (into the thickness continuum / bulk
   escape) **and** a **profile-functional, conditional** bound thickness/breathing pole (computed existence, ⛔ not
   assumed); and the confinement question (`N10`) as a computed object.
4. The **order bookkeeping** (`N12` + the inherited `(ε,η,σ_W)` multigrade, §3c): the full multigraded amplitude
   including its uniform component `A_0`, the **induced** amplitude `ΔA = A − A_0`, the transverse-channel flux
   pairing, and the flux-normalized conversion fraction on the named homotopy `λ ≡ η` (`L_W`, shapes fixed).
5. The **N6 insertion-construction independence control** (two scattering-coordinate constructions on the imported
   kernel + one-sided profile-atom corruption, §5a), the profile-**form** ablation (§5c), and the flux-normalized
   falsification **FORM** (§5d; the numeric bound withheld, `N7`).

**Out of scope (named, not solved).**
- The **order-unity slit edge** (finite-`ΔW`, non-perturbative) — a **downstream obligation named in §3d**, ⛔ not a
  mechanical reduction: a nonzero Born coefficient gives **no** general lower bound on order-unity conversion (§3d
  counterexample). Its construction (piecewise-uniform jump / matched interior-exterior / an effective interface
  response with explicit undetermined parameters) is **new** and belongs to the S11c-e strong-edge stage, where its
  weak limit must reproduce S11c-d.
- The **falsification magnitude** — needs the throat interior `R1` (`V3_STEP_PLAN.md:1179`); only the FORM is
  computable now (`N7`). The withheld `O(1)`/grating reductio is diffed **orchestrator-side**, ⛔ never a builder
  target.
- A **global dispersion `ω(k)`** for generic profiles (`N5`/`N10`, forbidden); the **nonlinear-light program**
  (`N10`); the **periodic→Bloch** and **slowly-varying→WKB** classes (§1d; ⛔ a periodic profile is a **distinct** N5
  object, not this interface S-matrix on a periodic profile).
- The **kernel-level Eulerian↔material representation N6** — that was **c2's** control and is carried as
  cross-engine-UNCLOSED DEBT (§1b); S11c-d ⛔ does not re-derive the material **closed operator** (it is not in the
  consume-set), it tests **insertion-construction** independence (§5a).
- The **full cross-engine self-energy operand residual** and c1's four giant families (≥64 GB,
  `DEFERRED_HEAVY_RUNS.md`). S11c-d must be constructible and cross-engine-testable on this box for its own mixing
  object, and **name — not silently absorb** — anything it cannot close here (§1b).

---

## 1 · Complete inherited setup — SUPPLIED and unfalsifiable in this build

Everything in §1 is an input. The mixing response, the two channels, the leakage FORM, and every control disposition
of §§3–5 are **outputs**; ⛔ none is stated here.

### 1a · Inheritance and the consumed c2 exports (the REAL export rows)

The DOFs, sector split, background ansatz, `(ε,η,σ_W)` power counting (`N12` + S11c-a §2), and admissibility are
exactly S11c-a §§1–2 / S11c-b §§1–2 / S11c-c2 §1, inherited by pointer. S11c-d consumes one already-built,
per-engine-reviewed model, **S11c-c2** (`scripts/S11c_c2_exports.py`; step record
`steps/S11c_c2_self_energy_fold.md`; disposition `_measurements/S11c_c2_N6_reconcile_disposition.md`). The c2 delta's
own-rows that S11c-d binds (exact `IMPORT_KEYS` root set fixed at the build directive, §7), by their **real** write-keys:

- **`s11cc2ClosedSlabOperator`** — the closure-modified variable-coefficient slab operator over `{u,θ,e_W}`, per
  `(anchoring α, density ρ)`. Its **diagonal blocks** supply the local sector operators / resolvents (the transverse
  response and the thickness response with its poles), including the two **asymptotic** diagonal operators (§2).
- **`s11cc2ClosedCouplingKernel`** — the re-extracted **off-diagonal** transverse↔`{θ,e_W,u_L}` block(s): the **full**
  mixing vertex (tilt `∇w₁`, modulus-gradient `∇m₁`, and N4 advection channels; ⛔ not `∇w₁` alone, ⛔ not `∇m₁`
  alone), multigraded `(ε,η,σ_W)`; its Fourier content is on the **momentum transfer** `Q = k_out − k_in`
  (`s11cc2FourierW1ProfileHatTransfer(−k_input + k_output, …)`), ⛔ not a single unspecified `k`.
- the **field carriers** `s11cc2Fieldtheta`, `s11cc2FieldeW`, `s11cc2Fieldu{1,2,3}`; the **profile coefficients**
  `s11cc2Coefficientw1Profile`, `s11cc2Coefficientm1Profile`; the **Fourier-of-profile** carriers
  `s11cc2FourierW1ProfileHatTransfer`, `s11cc2FourierW1ProfileJetHat*` (+ their `*Dimension` rows); and the inherited
  constants/kernels reachable through the fold (`W_0`, `mu_R`, `eta_bg`, `sigma_W`, `L_W`, `rho_m`, `rho_br`,
  `Lambda_{A,V,X}_0`, `tau_{A,V,X}`, `omega`, `c_s0`, `dtn_kernel`, `background_density_map`, …).

⛔ **There are NO c2 term-origin, parity, self-energy-increment, or six-§3d-re-adjudication export rows** — those are
**step-record provenance** (`steps/S11c_c2_self_energy_fold.md`; the increment was dropped to EMIT-only, `:56`), ⛔ not
importable objects. S11c-d neither binds nor assumes them; where it needs such provenance it derives it anew as an
S11c-d output.

### 1b · What is per-engine-SOUND vs cross-engine-UNCLOSED in the c2 import — SUPPLIED HONESTLY (rule 6/16)

⭐ **This is the load-bearing honesty section: S11c-d's physics IS the gradient-driven off-diagonal kernel, so c2's
carried cross-engine operand DEBT is MATERIAL to this consumer — ⛔ NOT dismissible on covariance alone**
(`steps/S11c_c2_self_energy_fold.md:214-219`).

⭐ **PER-ENGINE SOUND (SymPy, 2-leg):** the self-energy fold wiring + A/C/D1–D6 and the emitted
**`s11cc2ClosedSlabOperator` + `s11cc2ClosedCouplingKernel` VALUES**. These two rows are the operands S11c-d consumes.
(The self-energy **increment** value is per-engine-SOUND **parent provenance** — ⛔ it is EMIT-only, not an S11c-d
import.)

⭐ **CROSS-ENGINE, dual-engine confirmed (the N6 thread only, on this box):** **operator covariance (Reading B)** — the
material builder implements the declared frame change `Φ` (`R_cov` no-nonzero in both engines), on the specified
control/premise subset. ⚠ The matched covariance-channel cross-engine zeros are **`(0)−(0)`** — a dual-engine
confirmation of the **vanishing** statement, ⛔ **NOT operand agreement.** ⚠ **`R_N6 = I_E − I_{M→E} = 18/288` is the
per-engine SymPy RAW result and was ENTIRELY SCHEMA-UNMATCHED** — there is **no** direct cross-engine comparison of
`R_N6` itself (`_measurements/S11c_c2_N6_reconcile_disposition.md`); only the Reading-B covariance vanishing + the
control/premise subset are dual-engine. Both are preserved (`R_N6` nonzero raw, `R_cov` no-nonzero) as consistent
under Reading B.

⛔ **CROSS-ENGINE UNCLOSED — S11c-d must NOT treat these as closed (a supplied, unfalsifiable-in-this-build premise it
names honestly, `M2`/rule 16):**
- **the cross-engine OPERAND DEBT** — the surfaced blind-WL-vs-imported **carrier (40)**, constitutive **source (76)**,
  and **Φ (18)** residuals are **UNADJUDICATED**, and the leftover SHAPE was **not inspected**. The v3 collapse
  instrument was CLOSED **NOT-SOUND**; the only sound reconcile is **upstream of the EL differentiation**. ⛔ **Do NOT
  let "representational-difference-UNADJUDICATED" become "known to be just thickness"** — the alternatives include a
  genuine constitutive-convention mismatch OR an implementation error. Because S11c-d's mixing rides on
  `s11cc2ClosedCouplingKernel`, this DEBT is a **live premise on the object S11c-d builds** — S11c-d **names** it as
  unclosed and ⛔ does not pre-adjudicate it.
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
  conclusions do **not** stand, only the increment VALUES do. ⚠ **In particular, "the sectors decouple at uniform
  background" is exactly the withdrawn F** — S11c-d treats the uniform amplitude as a **computed** object `A_0` (§3c),
  ⛔ never as an assumed zero.

⛔ Folding any of the above to force cross-engine closure is the exact defect this rebuild exists to catch (rule 1/6).

### 1c · The background profile class — SUPPLIED framing (`N5`/`N12`/`N14`/`N15`)

⭐ **The named class is a LOCALIZED INTERFACE, stated on the INHERITED profiles along a NAMED in-plane normal** (S11c-a
§2a, `:171-192`). ⛔ Never write the varying thickness as `W_0`/`W₀(x)` — `W_0` and `mu_R` are **reserved constant
ledger keys** (`N14`); the varying fields are `W_bg`, `μ_R,bg`, via the inherited profiles `w₁ (=
s11cc2Coefficientw1Profile)`, `m₁ (= s11cc2Coefficientm1Profile)`. The in-plane coordinate `y` is a 3-vector; **name a
unit interface normal `n̂`** (e.g. `y¹`), with the profiles depending on `ξ ≡ n̂·y/L_W` and uniform along the edge:

```text
ξ ≡ n̂·y/L_W ,   W_bg ≡ W̄₀[1 + η w₁(ξ)] ,   μ_R,bg ≡ μ̄_R[1 + η m₁(ξ)] ,   σ_W ≡ η W̄₀/L_W ,
∂_{yᵢ}W_bg = σ_W n̂ᵢ w₁′(ξ) ,   ∂_{yᵢ}μ_R,bg = (μ̄_R/W̄₀) σ_W n̂ᵢ m₁′(ξ) ,   W̄₀ ≡ W_0 ,  μ̄_R ≡ mu_R .
```

- ⭐ **`w₁` and `m₁` are INDEPENDENT `O(1)` profiles** (S11c-a `:190`; ⛔ no engine may tie them, ⛔ no `m₁ = m₁[w₁]`).
  The **interface** condition (asymptotically constant, finite jump) is required on **both**: `Δw₁ ≡ w₁(+∞) − w₁(−∞) ≠
  0` (the thickness interface, so the asymptotes `W₋ ≠ W₊`, §2) **and** `Δm₁ ≡ m₁(+∞) − m₁(−∞) ≠ 0` (the modulus
  interface, the mixer). A specific shape (e.g. `(1+tanh ξ)/2`) is a **representative**, ⛔ not the class.
- ⭐ **An interface is NOT a bump.** The invariant discriminant is the **zero-transfer moment** `m̂₁′(Q_n=0) = Δm₁`
  (and `ŵ₁′(0) = Δw₁`): **nonzero** for an interface, **zero** for a bump. A bump returns to its background and has a
  different low-`Q` vertex; ⛔ a bump must **not** silently stand in for a single edge (§5c).
- **Momentum transfer / form factor.** Scattering conserves the tangential (edge-parallel) momentum; the **normal
  momentum transfer** is `Q_n ≡ k_{out,n} − k_{in,n}` (along `n̂`). The localized-edge form factor is `m̂₁′(Q_nL_W)`
  (and `ŵ₁′(Q_nL_W)`) — one reduced 1-D Fourier convention along `n̂`; the edge-parallel momenta stay live parameters.
- **Admissibility (`N12`).** Name which quantities vary (`W_bg`, `μ_R,bg`, `ρ_br,bg⁰` — both density representatives
  with their two asymptotic values `W_±`, `μ_±`), the anchoring (material-advected vs lab/Eulerian-held, inherited
  `α`, `N4`), and the stationary equations or the named force that holds the background — an inadmissible background
  silently sources spurious coupling.
- **Names (`N14`).** Every spatially-varying field/kernel/observable gets a **fresh** injective standard name; ⛔ never
  reuse an imported constant key (`W_0`, `mu_R`, `e_W`, `rho_br`, `v_0`, `slab_operator`, `coupling_kernel`, …) for a
  varying object.
- **Invariants (`N15`, at the RIGHT layer).** S11c-d consumes c2's operator/kernel **verbatim** — it ⛔ may **not**
  invent a new local constitutive constant (that is the variable-coefficient-operator stage's job, `N15`;
  `decisions:143`). It emits **profile moments and Fourier form factors DERIVED from the imported kernel** (e.g. the
  integrated-jump moment `Δm₁`, `m̂₁′(Q_nL_W)`) — profile/scattering **data**, ⛔ not new constitutive invariants. If
  the scattering exposes a **missing** local invariant, record it as an **upstream N15 debt**, ⛔ do not add it in d.

### 1d · The regime — SUPPLIED (Born in CONTRAST, sharpness LIVE)

⭐ **Born in the contrast `η`; the sharpness `σ_W` and the kinematic `Q_nL_W` are LIVE.** The three quantities are
distinct and must not be conflated (S11c-a `:189-198`):

```text
η      = contrast          (Born / scatterer strength;  a background bookkeeper) ,
σ_W    = η·W̄₀/L_W = first-jet / sharpness   (an INDEPENDENT background bookkeeper — ⛔ never a common order with η) ,
Q_nL_W = kinematic         (is the edge sharp on the WAVE — a PARAMETER, ⛔ not a grade) .
```

⚠ **The imported operator/kernel are ALREADY a first-`σ_W` shape expansion** — so "weak contrast" is well-defined;
what is **forbidden** is: ⛔ an **additional** `σ_W → 0` / `L_W → ∞` limit, ⛔ **expanding away the form factor**
`m̂₁′(Q_nL_W)`, and ⛔ setting **`η → O(1)` inside the first-shape-order operators** (an incomplete operator exactly
where the missing `O(η²)` pieces compete). ⭐ **The form-factor kinematics (corrected):** `Q_nL_W → 0` is the
**zero-transfer / sudden** limit, where `m̂₁′(Q_nL_W) → Δm₁` (the integrated jump, maximal conversion); `|Q_nL_W| ≫ 1`
(reached by `L_W → ∞`, i.e. `σ_W → 0`) is the **WKB / adiabatic** regime, where the form factor is suppressed
(Riemann–Lebesgue) and conversion is small — the class N5 says does not match a sharp edge. ⭐ Keep **Born in contrast
`η`, `σ_W` and the `Q_nL_W` form factor LIVE.** Weak **contrast** does not require WKB.

---

## 2 · The mixing object and the two-asymptote distorted-wave organization — SUPPLIED framing

The object is the **linear mixing** between the transverse sector and the thickness sector, carried by the **full**
imported off-diagonal vertex (tilt `∇w₁`, modulus-gradient `∇m₁`, N4 advection). ⚠ **Constitutive** mixing is *driven
by* `∇μ_R,bg ≠ 0` (`V3_STEP_PLAN.md:1179`, the modulus subchannel); ⛔ but this is **not** the whole support — a
`μ_R`-only projection **drops** the tilt and advection channels §1a requires and is an **ablation**, not the object.
The uniform-background amplitude is a **computed** object `A_0` (§3c), ⛔ not an assumed zero (the withdrawn F, §1b).

Because the interface has different asymptotes `W₋ ≠ W₊` (`Δw₁ ≠ 0`), the two ends carry **distinct** diagonal
operators, wavenumbers, and flux velocities — so name the two-asymptote distorted basis:

```text
L₋ = the W₋-asymptotic diagonal operator ,   L₊ = the W₊-asymptotic diagonal operator   (from s11cc2ClosedSlabOperator) ;
FIX one incident end and the open OUTGOING channels at each end, with an outgoing/Jost boundary prescription ;
mixing response  =  ( thickness outgoing resolvent )  ·  K_HT  ·  ( incoming transverse distorted solution )   at O(ε·[first shape order]) ,
where K_HT = the transverse→thickness block of s11cc2ClosedCouplingKernel (the FULL vertex), the reverse uses the other block .
```

⭐ The **scattering labels are frequency and asymptotic incoming/outgoing channels** — these do **not** require, and do
**not** imply, a global `ω(k)` (`N5`). Framing obligations:
- ⚠ **At the retained first-shape order and away from thresholds, uniform-mode Born and the two-asymptote
  distorted-wave construction coincide** (`G = G₀ + O(η,σ_W)`, `K = O(first shape)` ⇒ the continuum object is the
  uniform-mode matrix element of `K` at that order); ⛔ do not claim a retained-order difference. The two-asymptote
  basis is the **threshold / regular-domain organization** and fixes the `L±` kinematics — build the directive's
  regression on uniform modes only away from thresholds.
- ⭐ **Re-expand the continuum response to the retained background grade** (first order in each of `η`, `σ_W`); keeping
  the profile-dependent diagonal solutions unexpanded is a partial resummation, not a complete higher-order
  prediction. ⚠ **Exception:** a **bound-pole** spectral solve (§3b) is a **resummation** `G = (G₀⁻¹ − V)⁻¹` that the
  first-order continuum re-expansion **cannot** create — that spectral solve is **exempt** from the continuum
  re-expansion and carries its own stated error limitation. The "regular scattering domain" (exclude threshold /
  resonance enhancement and long coherent regions where repeated conversion accumulates) governs the **continuum**
  Born validity only; ⛔ it does not delete the §3b bound object.
- ⛔ **State no overall sign** in the insertion (an `i0`/Lippmann–Schwinger convention); let each engine's resolvent
  convention fix it. The comparator's load-bearing residual is the mixing amplitude, §7.

---

## 3 · The construction (OUTPUTS)

Every object below is computed for both anchorings `α` and both density representatives `ρ`; it carries its computed
`(ε,η,σ_W)` multigrade and restored `[L,T,M]` dimension, and states no component value, sign, order, parity, or grade
in this document.

### 3a · The profile-conditioned mixing response (continuum channel) — a CANONICAL S-matrix object

Emit the mixing response of §2 for the localized-interface profile, on the two-asymptote distorted basis, per
`(α,ρ)`. ⭐ **Fix the canonical object for blind comparison:** **one** named incident end (or the complete channel
matrix with both ends emitted separately), **one** amplitude normalization (flux-normalized to unit incoming flux, so
the amplitude carries the explicit `√(v_out/v_in)` factor — ⛔ do not offer "unit-flux OR `J_out/J_in`" as
interchangeable, they differ), and explicit branch labels `(incoming end, outgoing end, channel, ω, k_∥)`. Emit the
field-normalized raw amplitude separately if desired. Its continuum-channel projections are the **conversion /
scattering amplitudes**.

```text
⇒ S11CD_MIXING_RESPONSE (per (α,ρ), fixed incident end or full matrix) , S11CD_CONVERSION_AMPLITUDE (flux-normalized, labelled channels) ,
  S11CD_ASYMPTOTIC_OPERATORS (L₋, L₊, open-channel set, v_in/v_out) .
```

### 3b · The two DISTINCT photon-kill channels (`N13`) — the bound pole is a PROFILE-FUNCTIONAL CONDITIONAL

`N13`: "confinement of light" = **survival of the transverse polarization channel**. Conversion into a **bound**
breathing/thickness mode kills the photon **exactly as** bulk radiation does — the two are **distinct emitted
objects**, ⛔ not one "energy stays in the slab" statement.

- **(i) continuum conversion** — transverse → the thickness **continuum** / bulk escape (the §3a amplitude projected on
  the radiating/continuum channel).
- **(ii) bound-mode capture — a PROFILE-FUNCTIONAL, COMPUTED CONDITIONAL, ⛔ no class-wide existence claim.** A class
  with a fixed asymptotic jump contains both monotone steps **and** profiles with localized overshoots/wells, which
  have **different** pole sets — so existence/absence ⛔ cannot be asserted class-wide. ⛔ Do **not** invoke the 1D
  weak-well theorem (it needs equal asymptotes + an attractive self-adjoint well; the interface has different
  asymptotes and the thickness operator is multi-component, `ω`-dependent, nonlocal, possibly **non-Hermitian** from
  the outgoing bulk response). ⇒ emit a **profile-functional Jost / Evans determinant** and its zeros, with left/right
  residues and the nonlinear-eigenvalue (`∂_ωL_H`) normalization and spectral coupling — **PERMIT COMPUTED ABSENCE**.
  A **true bound pole** is a **physical-sheet, normalizable, zero-width** pole with **every radiation channel closed**
  (⛔ "below both continua" alone is insufficient for the non-Hermitian operator); distinguish it from a **second-sheet
  resonance**. A true bound state carries **no** asymptotic outgoing flux, so a **capture rate** requires a named
  protocol (wave-packet / switching, or a resonance width); ⛔ a pole residue + coupling alone is **not** a conversion
  probability — either define the protocol or emit only the spectral overlap / coupling. ⚠ A found pole is a localized
  **bound/resonance** object, ⛔ **not** a Bloch band (and an empty pole set is not a band either — the interface is
  nonperiodic regardless).
- The **confinement question** (`N10`): whether transverse-channel survival is **unconditional** — emit the computed
  object (a possibly-empty pole set + the continuum conversion), ⛔ not a claim.

```text
⇒ S11CD_CONTINUUM_CONVERSION , S11CD_BOUND_MODE_SPECTRAL_TEST (Jost/Evans determinant, zeros, residues, coupling; may be empty) ,
  S11CD_CAPTURE_PROTOCOL (or spectral-overlap-only) , S11CD_CONFINEMENT_CONDITION .
```

### 3c · The order bookkeeping — `(ε,η,σ_W)` multigrade + the named homotopy; the flux is the transverse continuum channel

⭐ **Every object is multigraded `(ε,η,σ_W)` from its actual data dependency (S11c-a `:195`); ⛔ no engine may assign a
common order to `η` and `σ_W`, and `Q_nL_W` is a kinematic parameter, ⛔ not a grade.** Retain first order in wave
(`ε¹`) and first shape order in **each** background bookkeeper (`η^{≤1}`, `σ_W^{≤1}`).

- **Emit the FULL multigraded amplitude including its uniform component `A_0`** (the `(η⁰,σ_W⁰)` part). ⛔ Do **not**
  assume `A_0 = 0` (that is the withdrawn F, §1b). The **leakage-relevant** object is the **induced amplitude**
  `ΔA ≡ A − A_0` (the part sourced by the gradient) — emit it as a named object.
- **The named homotopy.** `λ ≡ η` at **fixed** `L_W` and fixed shapes `w₁`, `m₁` (so `σ_W = η W̄₀/L_W` **tracks**; this
  is a physical **path**, ⛔ not a freeze of the formal `(ε,η,σ_W)` multigrade). On it the induced amplitude is
  `ΔA = O(ελ)`, the induced converted flux `O(ε²λ²)`, `J_in = O(ε²)`, and the flux-normalized fraction
  `C = J_conv/J_in = O(λ²) = O(η²)` (the incident `ε²` **cancels** in the linear theory). ⭐ Emit **both** the formal
  `(ε,η,σ_W)` amplitude multigrade **and** the homotopy `λ`-orders. The `N12` labels (`O(εη)` coupling, `O(ε²η²)`
  leakage, `O(η²)` fraction) apply to `ΔA` (⛔ do not attach them to the total `A` unless `A_0` is computed zero).
  ⚠ An `O(ε²λ²)` **observable** of a **linear** mixing is not the excluded nonlinear-light program (an `O(ε²)`
  **vertex** would be, `N10`/`N12`).
- **Flux normalization (the N13-correct observable).** `J` is the inherited **transverse-polarization-channel** energy
  flux (the S11b quadratic energy / c2 traction–slab pairing bilinear form). `J_conv` = the **continuum** (thickness /
  bulk) outgoing transverse-channel flux **ONLY** — ⛔ **do not fold bound capture into `J_conv`** (a stationary bound
  state carries no outgoing flux; a non-flux overlap in a flux is dimensionally illegal). ⛔ Do not use total
  mechanical energy or bare bulk Poynting (they hide the bound channel). ⚠ A **total photon-loss** fraction, if named,
  is a **separately-named** sum of the continuum fraction and the bound-capture probability — **same-dimension
  probabilities**, combined only **after** the §3b capture protocol exists, ⛔ never as a term inside `J`. The
  **magnitude** is `R1`-blocked; only the FORM is computed here.

```text
⇒ S11CD_AMPLITUDE_MULTIGRADE (incl. A_0) , S11CD_INDUCED_AMPLITUDE (ΔA) , S11CD_TRANSVERSE_FLUX_PAIRING ,
  S11CD_CONVERSION_FRACTION_FORM (continuum, O(λ²)) , S11CD_CONVERSION_POWER (absolute, O(ε²λ²)) .
```

### 3d · The strong-edge bridge — NAMED as a downstream obligation, ⛔ NOT solved, ⛔ NOT a reduction

⭐⭐ **S11c-d establishes the WEAK coefficients only; it does NOT establish the order-unity edge form.** A lab slit
edge is **localized + order-unity contrast**. Emit **separately** (⛔ do not conflate an amplitude slope with a
fraction slope): the leading **amplitude** coefficient `∂_λ(A_H/ε)|₀` **and** the leading **fraction** coefficient
`lim_{λ→0} C/λ² = ½C''(0)`. ⚠ **The slope of the fraction itself vanishes** (`C = O(λ²) ⇒ C'(0) = 0`) — so "the `O(η)`
slope of the conversion" is the wrong object. ⛔ **State no value or sign for the amplitude coefficient** — it is a
computed object and **can vanish** (a form-factor node `Q_nL_W = nπ` with `Δm₁ ≠ 0`, a `k·a = 0` selection rule); the
strong-edge statements below are **conditional on a nonzero computed coefficient**. The lab bounds the **strong
conversion fraction** `C_strong` at order-unity contrast, ⛔ **not** "`F(1)`" of an amplitude.

⛔⛔ **A nonzero Born coefficient supplies NO general positive lower bound on strong-edge conversion.** Illustrative
counterexample (⛔ not a model of the slab): a lossless two-mode coupler with dimensionless integrated coupling `G`,

```text
A_H = −i ε sin(η G) ,     C = sin²(η G) = η²G² + O(η⁴) ,     C = 0  at finite nonzero  η G = nπ .
```

The converted amplitude starts at `O(εη)` and the flux at `O(ε²η²)` **exactly as required**, yet the exact conversion
returns to **zero** at finite coupling. ⛔ **Do NOT evaluate the Born coefficient at `η=1` and compare that number to
the lab** — that identifies the weak coefficient with `C_strong(1)`, the invalid extrapolation `N7` names.

⭐ **The downstream obligation, stated now (owned by the S11c-e strong-edge stage):** a **justified finite-contrast
response** — a piecewise-uniform (finite-`ΔW`) matching, matched interior/exterior solutions, or an effective
interface response whose undetermined parameters remain **explicit** — whose **weak limit reproduces S11c-d**. This is
a **NEW construction** (each uniform side is `η`-exact; the finite-`ΔW` matching is not a first-jet kernel). ⚠ Born can
miss repeated conversion/reconversion, diagonal reflection, resonance shifts, and altered channel availability — these
change **frequency and angular dependence**, not just magnitude. ⇒ if the finite-contrast response cannot be
established in scope, the honest S11c-e outcome is a **conditional constraint on edge-response parameters** (or a
deferred numerical exclusion), ⛔ **not** a shape-independent exclusion — and the unknown interior coupling need **not**
factor as `C_edge(ω,ϑ) = C_interior·F(ω,ϑ)` (repeated scattering can put it inside resonance denominators).

```text
⇒ S11CD_WEAK_AMPLITUDE_COEFFICIENT (∂_λ(A_H/ε)|₀, no value/sign) , S11CD_WEAK_FRACTION_COEFFICIENT (½C''(0)) ,
  S11CD_STRONG_EDGE_OBLIGATION (named premise, ⛔ not solved here) .
```

---

## 4 · Objects to compute and emit

Per anchoring `α∈{L,M}` and density representative `ρ∈{ρ_4D,ρ_br}`, multigraded and dimensioned:

- The **mixing response** + conversion amplitude (canonical, flux-normalized, labelled) + the asymptotic operators /
  channels / velocities — §3a.
- The **two photon-kill channels** (continuum conversion; the **profile-functional** bound-pole spectral test —
  Jost/Evans zeros / residues / coupling, possibly empty; the capture protocol or spectral-overlap-only) + the
  confinement condition — §3b.
- The **leakage bookkeeping** — the full multigraded amplitude incl. `A_0`, the induced amplitude `ΔA`, the
  transverse-flux pairing, the continuum conversion FRACTION FORM (`O(λ²)`), and the absolute converted power
  (`O(ε²λ²)`) — §3c.
- The **weak amplitude coefficient** (no value/sign) and **weak fraction coefficient**, and the named strong-edge
  obligation — §3d.
- The **control outputs** of §5, each emitted as the object and its literal residual.
- **Profile moments / form factors** derived from the imported kernel (`N15` data, ⛔ no new constitutive constants).

Every result carries its `(ε,η,σ_W)` order (and, on the homotopy, its `λ`-order) and its restored `[L,T,M]` dimension.
⛔ No result is reported without both.

---

## 5 · Independent routes and controls

⭐ Every control re-enters the chain **at the ACTION / the imported operands**, ⛔ never at a result. Each emits the
object and its literal residual; ⛔ none asserts a target value. A **coefficient** rescale tests arithmetic; only a
**form** change tests physics.

### 5a · The N6 control — INSERTION-construction independence on the imported kernel + one-sided profile-atom corruption

⚠ **The kernel-level Eulerian↔material representation N6 was c2's control** (and is cross-engine-UNCLOSED DEBT, §1b);
S11c-d ⛔ does **not** re-derive the material **closed operator** (it is not in the consume-set, and no face-normal
carrier factory is imported). S11c-d's N6 tests the **INSERTION / scattering-construction** coordinate-independence on
the **imported** closed kernel:

```text
route 1 (Eulerian scattering coords):   the §2 mixing insertion (distorted waves + flux pairing) built in the Eulerian
                                        in-plane scattering coordinates, on the imported closed kernel ;
route 2 (material-flattened coords):    the SAME insertion built after flattening the interface faces to material
                                        in-plane coordinates and transforming the SCATTERING problem back into the
                                        common Eulerian basis (the coordinate map lives INSIDE the construction) ;
S11CD_INSERTION_INVARIANCE_RESIDUAL[α,ρ] = route1 − route2      (differenced DIRECTLY in the common basis) .
```

⛔⛔ **No `Φ` on the amplitude / modes / flux.** `Φ` is the constitutive field map (c2's `R_cov` acted on the source
`μ`); it is **not** defined on scattering amplitudes, incoming/outgoing modes, or the flux bilinear form — writing
`route1 − Φ(route2)` is the c2 type error relocated (a θ-shift that annihilates). Difference the two constructions
directly in the common Eulerian basis (⛔ no separate final transform on the amplitude, as c2 §5c). If a
naturality/covariance residual is emitted at all, it is the naturality of the **maps that actually act on scattering
data** (the incoming-data map, the measure) — ⛔ never a global `Φ` on the mixing amplitude.

⭐ Both residuals are **computed measurements** — their **values are the findings**, ⛔ no target is supplied, ⛔ never a
builder exit condition; the diff is adjudicated on **our** side; keep the carrier/source/Φ operand DEBT (§1b)
explicitly unclosed.

**The one-sided independence corruption** acts on the SOURCES THAT EXIST at this step — the **imported kernel's
explicit profile atoms** (⛔ not by editing the whole operator, ⛔ not a fresh face level-set rebuild), still at fixed
`α,ρ`, **PRINT** the residual (⛔ do not require a nonzero value):
- **(i) tilt probe (`N3`):** reverse the **thickness first-jet slope atom** (`w₁′`/the `∇w₁` factor in the kernel) on
  **one** route only.
- **(ii) N4 advection probe:** omit / flip the **advective-density atom** (`u·∇ρ₄/ρ₄`) on **one** route only; ⚠ this
  atom is **structurally absent for `RHO4_CONSTANT`** (`∇ρ₄=0`) and present for `RHOBR_CONSTANT` — the live probe is
  `RHOBR_CONSTANT`; for `RHO4_CONSTANT` emit the **computed absence**, ⛔ never an `A−A`.

⚠ There are **≥2 same-order channels** (tilt `N3`; advection `N4`); the one-sided corruption is the **independence
test between them**. ⛔ `∇W_bg→0` / `η→0` is **NOT** an accepted corruption (the vacuous uniform limit renamed, `N6`);
⛔ corrupting one **anchoring** is not this test; ⛔ `Δρ` never bridges `LAB_HELD ↔ MATERIAL_ADVECTED`.

```text
⇒ S11CD_INSERTION_INVARIANCE_RESIDUAL[α,ρ] , S11CD_CONTROL_INDEPENDENCE_{BASE,CORRUPTED,RESIDUAL}[α,ρ,probe] .
```

### 5b · The uniform regression — zero the JETS, keep the asymptotes live

⭐ The forbidden object the smoke test must catch is a **gradient-independent** coupling term — which a bare `η→0`
would erase automatically (`K_bad = c·η → 0` whatever `c`). ⇒ the regression sets the profile **jets** to zero
(`w₁′ = m₁′ = 0`, i.e. constant `w₁`, `m₁`) while retaining **live** `η` and the **arbitrary constant asymptotes**
`W_±`, `μ_±` — and emits the computed coupling, which **must** be a computed object (⛔ not asserted "must vanish":
whether the closure-induced coupling decouples at uniform background is the withdrawn F, §1b). Keep the full
`η = σ_W = 0` reference as an **additional** regression. The uniform regression cannot see the gradient coupling's
coefficient/sign/parity — it is **secondary**, ⛔ not the N6 control (§5a).

```text
⇒ S11CD_UNIFORM_REGRESSION (jets→0, asymptotes live; + the η=σ_W=0 reference) .
```

### 5c · The profile-FORM ablation + the edge-vs-bump discriminant

- **Form ablation.** Perturb the **FORM** of `m₁(ξ)` (and `w₁(ξ)`) and emit the **baseline operand, the altered-form
  operand, and their residual** — ⛔ do not assert that the mixing "must move"; whether it moved is adjudicated on our
  side (a coefficient rescale of `η` would be arithmetic, only a form change tests the coupling).
- **Edge-vs-bump discriminant — the zero-transfer moment.** The canonical discriminant is the zero-transfer moment
  `m̂₁′(Q_n=0) = Δm₁` (and `ŵ₁′(0) = Δw₁`): **nonzero** for an interface, **zero** for a bump. Emit the low-`Q`
  (zero-transfer) content of the mixing vertex and its dependence on that moment; a bump (zero moment) has a different
  low-`Q` vertex and ⛔ must not stand in for the edge. ⛔ Do not freeze a specific representative shape as "the slit."

```text
⇒ S11CD_PROFILE_FORM_ABLATION (baseline, altered, residual) , S11CD_ZERO_TRANSFER_MOMENT_DEPENDENCE .
```

### 5d · The falsification FORM control (`N7`; the boundary refinement is N2-permitted)

Emit the **flux-normalized dimensionless conversion FORM** (§3c, continuum, `O(λ²)`; `∝ k·a`, carrying the
`m̂₁′(Q_nL_W)` form factor). ⚠ The FORM is the **full** imported vertex projected on the conversion channel — ⛔ not a
`∇μ_R`-only projection. The decision-list table places the FORM and confinement interpretation at S11c-e; N2 permits a
spec-stage boundary refinement, so computing the FORM here is not a defect — but ⛔ d's FORM must **not** become e's
withheld `O(1)` target. ⛔⛔ The `O(1)`-fraction / diffraction-grating reductio and the withheld numeric lab bound are
**orchestrator-side**, ⛔ never in the builder-facing acceptance text; the magnitude is `R1`-blocked. ⚠ A slit edge is
an order-unity localized gradient — the FORM here is the **weak-contrast** object (§3d), ⛔ not a non-perturbative lab
number.

---

## 6 · Method, dimensions, and script obligations

- **Method.** Balance laws + the binding material virtual-displacement rule + variational derivatives with held-fixed
  fields named + prescribed external virtual work (S11b), ⛔ never an irreversible response kernel in an ordinary
  action. The diagonal responses and the off-diagonal vertex are the c2 exports consumed verbatim; the mixing is one
  kernel insertion on the two-asymptote distorted basis (§2), re-expanded to the retained background grade (the
  bound-pole spectral solve exempt, §2/§3b).
- **Dimensions.** Restore `[L,T,M]` on every emitted object, dimensional consistency able-to-fail
  ([[feedback_dimensional_consistency_check]]); `(ε,η,σ_W)` multigrade on every object (`N12`).
- **Rest-frame limit.** Inherit `N11a` inert; S11c-d constructs **no** convective operator. Every result inherits the
  c1/S11b smallness domain (`|q_out·v_bulk_normal_0/ω|≪1` + boundary-layer/subsonic; large `k c_s0/|ω|` is
  **necessary, ⛔ not sufficient**), ⛔ never aliasing `v_bulk_normal_0` to `v_0` (`N14`/`N11`).
- **Script obligations.** The three build-skill clauses bind the build directive (`.claude/skills/build/SKILL.md`): a
  script PRINTS computed objects and never states conclusions; PRINT the residual, do not assert it; interpretation is
  the step record. ⛔ No hand-typed CAS object standing in for a computed one; every control re-enters at the
  ACTION / imported operands. ⛔ No tautological residual (rule 2 corollary 3): the §3a/§3c export representations are
  **not** checks; the §5 residuals are emitted with both operands, and a two-route residual is emitted only where an
  **independent** second route exists.
- **Serialize CAS jobs; watch RSS.** c2's self-energy `.out` was ~499 MB and the full cross-engine residual is the
  ≥64 GB work; S11c-d's mixing is a **derived insertion** on the closed operator (plus a spectral solve) — measure the
  process that runs, defer heavy controls in-band→out-of-band (`DEFERRED_HEAVY_RUNS.md`), ⛔ never two memory-heavy CAS
  jobs concurrently. Detached launch (harness reaps `run_in_background`). Mathematica: 2-seat licence,
  `--sandbox danger-full-access`, serialize dual ablations.

---

## 7 · Names, F9 reservations, chain output, and export schema

**F9 / `N14` reservations.** Every new object gets a **fresh** injective `mechanical_lower_camel` name; ⛔ never reuse
an imported S11c-c2/c1/S11c-b/S11b key (`s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel`, `slab_operator`,
`coupling_kernel`, `dtn_kernel`, `mu_theta_operator`, `W_0`, `mu_R`, `e_W`, `rho_br`, `v_0`, …) for a new S11c-d
object.

**Chain output (`N1`/`N8`; topology = the two-leg-gated `directives/export_ledger_bind_closure_design.md` §D1–§D3).**
The SymPy engine reads the inherited model via the **positional** `load_model` call over the frozen base
`scripts/S11c_b_exports.py` with the c1 and c2 deltas folded on top (signature `load_model(base_path, *delta_paths)`,
`scripts/ledger_fold.py:102`; ⛔ NOT keyword form), binding only its declared `IMPORT_KEYS`, and writes
`scripts/S11c_d_exports.py` as its **own-rows delta** (§D2, ⛔ not the accumulated whole-model file). ⛔ The exact
`IMPORT_KEYS` **root set** (minimal roots whose recursive closure covers the §1a consume-set — the two c2 closed
rows + the field/coefficient/Fourier carriers + the reachable constants) is fixed at the **build directive** against
the real export files, ⛔ not enumerated-then-frozen here; its two decision legs verify it, and that the guard
(`check_consumer`/`assert_lookups_equal_manifest`/`assert_delta_is_minimal`) passes on the fold — ⚠ noting the guard
passes on **key existence**, so it will **not** catch a wrong-provenance binding (⛔ the `s11cc2Fieldtheta`-vs-
`s11cc2FieldTheta` casing and the increment-vs-operator distinction are exactly this hazard); that is the directive's +
legs' responsibility. `BUILD_INPUT_DIGESTS` pins, per §D3, `{this sub-step's SymPy audit, scripts/S11c_b_exports.py,
scripts/S11c_c1_exports.py, scripts/S11c_c2_exports.py, this spec, scripts/ledger_fold.py}`. ⛔ Never `git add -f` a
big `.out`; ⛔ never annex an `*_exports.py`.

**The comparator (`N8`, frozen `T7` contract).** The S11c-d comparator joins the two blind engines' emitted objects by
name, pairs residual operands, is three-valued, rejects a native boolean, and PRINTS/decides nothing (rule 2). ⚠ Its
load-bearing residual is on the **mixing amplitude** (§3a) — which **rides on the carried cross-engine operand DEBT**
(§1b); the comparator **SURFACES** the DEBT and the §1b representation questions (the staged representational bridge,
[[feedback_reconcile_representational_bridge]], ⛔ never a blanket collapse), ⛔ does not pre-adjudicate them. The full
per-object symbolic residual + c1's four giants remain deferred (`DEFERRED_HEAVY_RUNS.md`); S11c-d names, ⛔ does not
pre-adjudicate, whatever it cannot close on this box.

**The blind Wolfram engine** re-derives the §§1–2 supplied inputs, the S11c-a face substrate, the S11c-b slab-operator
and c2 closed-operator/kernel rows it consumes, and the localized-interface mixing — importing nothing (the only
cross-engine control). ⛔ The denylist stays cut (`N9`/rule 12); blindness is enforced by absence.

---

## 8 · Supplied versus computed; builder report

**SUPPLIED (unfalsifiable in this build):** all of §1 (the two c2 export operand rows and their per-engine-SOUND vs
cross-engine-UNCLOSED disposition — the operand DEBT, the raw/schema-unmatched `R_N6`, the two S11c-b signs, the six
§3d re-adjudications, the 3 N6 premise caveats, the withdrawn F/G, and that no term-origin/parity/increment/§3d rows
are importable), the §1c localized-interface class on the inherited `w₁`/`m₁` along the named normal `n̂` and its
admissibility, the §1d regime (Born-in-contrast, `σ_W`/`Q_nL_W` live), the §2 two-asymptote distorted-wave
organization, `N11a`, `N12`, `N13`.

**COMPUTED (outputs, ⛔ none stated here):** the mixing response and canonical conversion amplitude and asymptotic
operators (§3a); the two photon-kill channels (the **profile-functional conditional** bound-pole spectral test) and the
confinement condition (§3b); the full multigraded amplitude incl. `A_0`, the induced amplitude, the transverse-flux
pairing, the continuum conversion FRACTION FORM and the absolute converted power (§3c); the weak amplitude coefficient
(no value/sign) and weak fraction coefficient (§3d); every control residual — the insertion-invariance and one-sided
corruption residuals (§5a), the uniform regression (§5b); the profile moments / form factors (`N15` data); the
`(ε,η,σ_W)`/`λ` orders and `[L,T,M]` dimensions.

**Builder report.** The build directive states, per emitted object, which line computed it (`.claude/skills/build`);
declares the profile class (§1c), the regime grades (§1d), the distorted-basis insertion organization (§2), and the
leakage bookkeeping (§3c) it implemented; and reports the literal residuals of §5 — ⛔ never a prose conclusion. The
disposition of every §5 residual and the strong-edge obligation (§3d) is read on **our** side, in the step record,
⛔ not asserted by the script (rule 5).
