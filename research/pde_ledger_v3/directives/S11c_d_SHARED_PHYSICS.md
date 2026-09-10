# S11c-d — SHARED PHYSICS v4 (profile-conditioned transverse↔thickness scattering at a thickness interface)

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

⭐ This is a **Codex `gpt-5.6-sol`-authored physics spec**, re-authored as v4 under `CLAUDE.md` rule 15 after the
orchestrator-authored v2/v3 folds bred new defects. It is the physics authority both blind engines read. Per G1/G2 it
is physics-bearing and gets **two non-author legs — a fresh Claude agent + Grok — reviewed UNTIL CLEAR** (spec row,
⛔ not the decision-list one-pass); both reports precede any commit, and a reviewed baseline is preserved before a
repair overwrites it. The **build directive** that follows this spec gets its own two decision legs before any builder
(the G2 TRIGGER). **Spec v4** retains the v3 §0–§8 structural frame and folds every round-3 finding; no wording from a
prior version is authoritative where the round-3 derivation corrected it.

⭐ **The profile-class + regime decision (§1c–§1d) was settled by a three-way physics consult** (orchestrator + `gpt-6-astra`
xhigh + `grok-4.6`, `_legs/S11c_d_profile_class_consult{,_astra,_grok}.md`) and the user's approval: **localized
interface, Born in contrast `η` with sharpness `σ_W` and kinematics `Q_nL_W` kept live** (§1d).

---

## 0 · Scope

**In scope.**
1. **NAME the profile class** = a **localized interface** in the inherited background profiles (§1c): smooth,
   asymptotically constant, with short-range retained jets and a finite integrated **thickness** jump
   `Δw₁ ≠ 0`, hence `W₋ ≠ W₊`. This is the `N5` class gate. The independent modulus profile `m₁` may be constant,
   a localized bump, or an interface; `Δm₁` is a per-profile discriminant of the modulus subchannel, ⛔ not a class
   gate (§1c/§5c).
2. The **profile-conditioned transverse↔thickness mixing response** (§3a): the linear mixing (the **full** imported
   off-diagonal vertex — tilt `∇w₁`, modulus-gradient `∇m₁`, N4 advection), built on the **two-asymptote distorted
   basis** of the full `W₋`- and `W₊`-asymptotic block operators obtained from `s11cc2ClosedSlabOperator` and
   `s11cc2ClosedCouplingKernel`; the uniform and two end baselines are computed before any sector-diagonal
   simplification (§2).
3. The **two DISTINCT photon-kill channels** (`N13`, §3b): continuum conversion (into the thickness continuum / bulk
   escape) **and** a **profile-functional, conditional** bound thickness/breathing pole (computed existence, ⛔ not
   assumed); and the confinement question (`N10`) as a computed object.
4. The **order bookkeeping** (`N12` + the inherited `(ε,η,σ_W)` multigrade, §3c): the full multigraded amplitude,
   its uniform baseline, its separate zero-jet-contrast and first-jet pieces, the physical total conversion flux and
   its baseline/interference slots, and the separately named induced-field quadratic form on the homotopy
   `λ ≡ η` (`L_W`, shapes fixed).
5. The **scattering-coordinate covariance regression** on the imported closed kernel and its one-sided
   shape-sensitivity mutations (§5a). It does **not** discharge c2's kernel-level N6/N3/N4 independence debt. Also in
   scope are the profile-**FORM** ablation (§5c) and the flux-normalized falsification **FORM** (§5d; the numeric bound
   withheld, `N7`).

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
- The **kernel-level Eulerian↔material representation N6 and independent N3/N4 shape construction** — those were
  **c2's** control and remain carried as cross-engine-UNCLOSED DEBT (§1b). S11c-d does not import native
  pre-extraction material/Eulerian operands or channel-origin provenance and therefore cannot discharge them by a
  coordinate rewrite of the already-constructed closed kernel (§5a).
- The **full cross-engine self-energy operand residual** and c1's four giant families (≥64 GB,
  `DEFERRED_HEAVY_RUNS.md`). S11c-d must be constructible and cross-engine-testable on this box for its own mixing
  object, and **name — not silently absorb** — anything it cannot close here (§1b).

---

## 1 · Complete inherited setup — SUPPLIED and unfalsifiable in this build

Everything in §1 is an input. The mixing response, the two channels, the leakage FORM, and every control disposition
of §§3–5 are **outputs**; ⛔ none is stated here.

### 1a · Inheritance and the consumed c2 exports (the REAL export rows)

The DOFs, sector split, background ansatz, `(ε,η,σ_W)` power counting (`N12` + S11c-a §2), and admissibility are
exactly S11c-a §§1–2 / S11c-b §§1–2 / S11c-c2 §1, inherited by pointer. S11c-d consumes the already-built,
per-engine-reviewed **S11c-c2** model (reconstructed from the frozen base + the c1 delta + the c2 delta via
`load_model`, §7 — c2 is a bind-closure own-rows delta, ⛔ not an accumulated whole-model file)
(`scripts/S11c_c2_exports.py`; step record
`steps/S11c_c2_self_energy_fold.md`; disposition `_measurements/S11c_c2_N6_reconcile_disposition.md`). The c2 delta's
own-rows that S11c-d binds (exact `IMPORT_KEYS` root set fixed at the build directive, §7), by their **real** write-keys:

- **`s11cc2ClosedSlabOperator`** — the closure-modified variable-coefficient slab operator over `{u,θ,e_W}`, per
  `(anchoring α, density ρ)`. Together with the closed coupling kernel, it supplies the local sector blocks and the
  **full** two-end asymptotic block operators/resolvents, including every computed constant off-diagonal baseline
  (§2).
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
importable objects. S11c-d neither binds nor assumes them. It may emit provenance of its own downstream projections,
but cannot retroactively reconstruct c2 term origins or claim channel isolation from the closed kernel (§5a).

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
confirmation of the **vanishing** statement, ⛔ **NOT operand agreement.** ⚠ **`R_N6 = I_E − I_{M→E}` is nonzero as a
per-engine SymPy RAW census (18 of 288 columns nonzero, ⛔ a count — not an algebraic residual value) and was ENTIRELY
SCHEMA-UNMATCHED** — there is **no** direct cross-engine comparison of `R_N6` itself
(`_measurements/S11c_c2_N6_reconcile_disposition.md`); only the Reading-B covariance vanishing + the control/premise
subset are dual-engine. Both are preserved (`R_N6` nonzero raw census, `R_cov` no-nonzero) as consistent
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
  conclusions do **not** stand, only the increment VALUES do. ⚠ **The withdrawn F is specifically the c2 *increment*
  interpretation** (whether the closure-induced coupling in the closed kernel decouples at uniform background — the
  retired `verify_F` question). ⛔ It does **not** retract S11b's uniform decoupling, which stands as prior art
  (`M3` oracle). ⇒ for c2's closed kernel, S11c-d treats the uniform amplitude/baseline as a **computed** object
  `A_0`/`K_0` (§3c), ⛔ never as an assumed zero (⛔ do not type `K_0=0` from S11b) — the computation, not the F label,
  settles it.

⛔ Folding any of the above to force cross-engine closure is the exact defect this rebuild exists to catch (rule 1/6).

### 1c · The background profile class, density branches, and Fourier convention — SUPPLIED framing

⭐ **The named `N5` class is a LOCALIZED THICKNESS INTERFACE along a named in-plane normal.** The inherited profile
definition begins at `S11c_a_SHARED_PHYSICS.md:171`. The in-plane coordinate `y∈R³`; choose a unit normal `n̂` to the
edge, set `ξ ≡ n̂·y/L_W`, and take every background coefficient uniform in the two directions parallel to the edge.
`W_0` and `mu_R` remain reserved constant ledger keys; the fresh varying fields and inherited independent profiles are

```text
W̄₀ ≡ W_0 ,   μ̄_R ≡ mu_R ,
W_bg(y)   ≡ W̄₀[1 + η w₁(ξ)] ,
μ_R,bg(y) ≡ μ̄_R[1 + η m₁(ξ)] ,
σ_W       ≡ η W̄₀/L_W ,
∂_{yᵢ}W_bg   = σ_W n̂ᵢ w₁′(ξ) ,
∂_{yᵢ}μ_R,bg = (μ̄_R/W̄₀) σ_W n̂ᵢ m₁′(ξ) .
```

`w₁ (= s11cc2Coefficientw1Profile)` and `m₁ (= s11cc2Coefficientm1Profile)` are independent dimensionless `O(1)`
profiles; no engine ties them or substitutes `m₁=m₁[w₁]`. Require finite one-sided limits for both. The class gate is

```text
Δw₁ ≡ w₁(+∞) − w₁(−∞) ≠ 0 ,
W₋ ≡ W̄₀[1+ηw₁(−∞)] ,   W₊ ≡ W̄₀[1+ηw₁(+∞)] ,   hence W₋ ≠ W₊ .
```

The modulus profile remains independent and in-class when it is constant, a localized bump, or an interface. Its
`Δm₁ ≡ m₁(+∞)−m₁(−∞)` is a computed per-profile discriminant of the modulus subchannel, not membership in the
thickness-interface class. Carry `μ₋≡μ̄_R[1+ηm₁(−∞)]` and `μ₊≡μ̄_R[1+ηm₁(+∞)]` independently of `W₋,W₊`. A
representative such as `(1+tanh ξ)/2` is a profile instance, never the class.

**Short-range retained jets.** For each `f∈{w₁,m₁}` and every positive derivative order `r` actually consumed by the
retained imported kernel, require `f^(r)∈L¹(R)` and the regularity needed for the kernel operations and integrations by
parts, together with `f^(r)(ξ)→0` at both ends so the constant asymptotic operators are defined. This premise applies
to every consumed derivative, not merely `f′`; it is what licenses the
Riemann–Lebesgue statement in §1d. The unequal-asymptote profiles `f` themselves do not decay and are not treated as
ordinary `L¹` Fourier functions. A zero-jet profile factor is kept in coordinate space. If it is transformed, fix the
profile origin and use the explicit canonical subtraction

```text
f(ξ) = f(−∞) + Δf H(ξ) + f_loc(ξ) ,
f_loc(ξ) ≡ f(ξ)−f(−∞)−Δf H(ξ) ,
```

with the additional half-line tail premise that makes `f_loc∈L¹`. Retain the delta/principal-value distributional
transform of the constant-plus-Heaviside part. If that tail premise is unavailable, the zero-jet step remains in
coordinate space; it is not assigned an ordinary transform.

**Exact reduced transform and the c2 carrier map.** Put

```text
Q ≡ k_out − k_in ,   Q_n ≡ n̂·Q ,   Q_∥ ≡ Q − n̂Q_n ,   s ≡ Q_nL_W ,
f̂_red(s) ≡ ∫_{−∞}^{∞} dξ exp(−isξ) f(ξ) ,
f(ξ) = (1/2π)∫_{−∞}^{∞} ds exp(+isξ) f̂_red(s) .
```

The first line is an ordinary integral only for a localized/subtracted `f`; for a full step it is distributional.

⚠⚠ **The imported c2 carrier's normalization is NOT supplied — it must be COMPUTED, ⛔ not asserted.** The c2 engine's
transform convention is an **unnormalised forward transform with a normalised inverse** (`c1` uses `DiracDelta(k−k′)`
with no `(2π)³` coefficient; the `(2π)⁻³` sits on the application/inverse — the c2 self-energy fold's own stated
convention, `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`), so `s11cc2FourierW1ProfileHatTransfer` /
`s11cc2FourierW1ProfileJetHat{i}` are ⛔ **not** equal to the normalized `(2π)⁻³∫d³y e^{−iQ·y}(·)`. ⛔ **Do not supply
a numeric `(2π)`/`δ²(Q_∥)` map for the imported carrier** — a supplied constant map about an imported object is an
`M2` leak, and a wrong one (e.g. off by `(2π)³`) is a defect **both engines would share**. ⇒ **each engine COMPUTES
the 3-D→1-D reduction factor** by reducing the one-profile insertion **as it actually appears in
`s11cc2ClosedCouplingKernel`** against the imported kernel's own flat identity (its `DiracDelta³(k_out−k_in)` piece),
using the S11c-d reduced transform `f̂_red` above and the geometry `∂_{yᵢ}f = n̂ᵢ f′(ξ)/L_W`; it **emits that reduction
factor as an object with BOTH operands** (the imported carrier, and the derived reduced form), so the tangential
`δ²(Q_∥)`, the `2π`, and the dimensional content **fall out of that computation**, ⛔ not a supplied constant.

**Reconstruction round-trips against the DERIVED map, ⛔ not a tautology.** The 3-D↔1-D reconstruction each engine emits
compares the imported 3-D carrier to its **computed** reduced form — ⛔ **not** a defined
`A_3D ≡ [L_W/(2π)]δ²(Q_∥)A_edge` (that is `A−A`: it manufactures the 3-D object from the reduced one and checks
nothing). The comparator joins the **reduced kernels so obtained** (each engine's computed reduction), ⛔ **not** a
pre-factored `[L_W/(2π)]` coefficient. ⚠ Separately, when forming flux/rate, neither engine squares an unstripped
`δ²(Q_∥)`: the standard tangential box/continuum normalization is removed before the per-unit-edge-area limit (a
flux-normalization convention S11c-d fixes, distinct from — and ⛔ not a substitute for — the computed carrier
reduction above).

**Branchwise density admissibility (`N12`/`N4`), verbatim from S11c-a §2b.** Let
`ρ_4D,ref⁰ ≡ rho_br/W_0` and carry the full fields, their gradients, and both asymptotic values live:

```text
RHO4_CONSTANT:
    ρ_4D,bg⁰(y) ≡ ρ_4D,ref⁰                         (constant),
    ρ_br,bg⁰(y) ≡ ρ_4D,bg⁰(y) W_bg(y)               (varies) ;

RHOBR_CONSTANT:
    ρ_br,bg⁰(y) ≡ rho_br                            (constant),
    ρ_4D,bg⁰(y) ≡ rho_br/W_bg(y)                    (varies) .
```

The live asymptotic and gradient maps are

```text
RHO4_CONSTANT:
    ρ_4D,±⁰ = ρ_4D,ref⁰ ,        ∇ρ_4D,bg⁰ = 0 ,
    ρ_br,±⁰ = ρ_4D,ref⁰W_± ,     ∇ρ_br,bg⁰ = ρ_4D,ref⁰∇W_bg ;

RHOBR_CONSTANT:
    ρ_br,±⁰ = rho_br ,           ∇ρ_br,bg⁰ = 0 ,
    ρ_4D,±⁰ = rho_br/W_± ,       ∇ρ_4D,bg⁰ = −rho_br∇W_bg/W_bg² .
```

In particular, the `RHOBR_CONSTANT` branch retains `∇ρ_4D,bg⁰`; it does not make `ρ_br,bg⁰` vary. For each density
representative and each inherited anchoring `α∈{LAB_HELD,MATERIAL_ADVECTED}`, state the stationary equations or the
named external force that holds the background. An inadmissible background can source spurious coupling.

**Names and invariants (`N14`/`N15`, at the right layer).** Every new varying field, kernel, and observable gets a
fresh injective standard name; imported constants such as `W_0`, `mu_R`, `e_W`, `rho_br`, `v_0`, `slab_operator`, and
`coupling_kernel` are never reused. S11c-d consumes c2's operator/kernel verbatim and does not invent a local
constitutive constant. It emits profile moments and form factors derived from that kernel as scattering data. A
missing local invariant is recorded as upstream `N15` debt, not added here.

### 1d · The regime — SUPPLIED (contrast Born; sharpness and two kinematic axes live)

Keep the background grades and kinematics distinct:

```text
η          = zero-jet contrast and Born/scatterer-strength bookkeeper ,
σ_W        = ηW̄₀/L_W = independent first-jet/sharpness bookkeeper ,
s=Q_nL_W   = momentum-transfer argument of a profile form factor (not a grade) ,
k_aL_W     = local wavelength/edge-sharpness parameter for channel a (not a grade) .
```

The imported operator/kernel are already truncated at first order in each of `η` and `σ_W`. Do not take an additional
`σ_W→0` expansion, do not Taylor-expand away the `s` dependence, and do not set `η` to order unity inside this
first-shape-order operator. Do not take `L_W→0` at fixed `η`, because that drives `σ_W→∞` outside the retained model.

`s→0` is only the **zero-transfer limit**; `(f′)̂_red(0)=Δf`. It is not named “sudden” and supplies no statement about
where the total conversion is largest: overshoot profiles, other vertex channels, and their current weights remain
live. For `f′∈L¹`, large `|s|` suppresses that derivative form factor by
Riemann–Lebesgue; it does not by itself establish suppression of the full vertex or total conversion. A form-factor
zero is a zero of the computed form factor, if one exists; no class-generic node location is supplied.

WKB/adiabaticity is separate: it requires local wavelength and modal-gap control, for example
`|k_a(ξ)|L_W≫1` away from turning points together with `L_W|k_a(ξ)−k_b(ξ)|≫1` (or the equivalent derivative/gap
criterion) for coupled channels. That is not inferred from large `|s|` and is not performed here as an extra
`σ_W→0` expansion. Born in contrast, first-jet sharpness, transfer kinematics, and local adiabaticity remain separate
live quantities.

---

## 2 · The mixing object and the two-asymptote distorted-wave organization — SUPPLIED framing

The object is the linear scattering/mixing generated by the **full imported closed operator**, with
`s11cc2ClosedCouplingKernel` used as its canonical off-diagonal extraction. It includes the tilt, modulus-gradient,
and N4 advection content carried by that kernel. The modulus-gradient projection is only one subchannel; projecting
onto it alone drops other retained terms and is an ablation, not the S11c-d object.

Before defining channels, each engine block-decomposes the imported operator at the retained grades. In the display
below `𝓛(y)` denotes the operator/kernel, not multiplication by a local matrix. For a nonlocal normal kernel, an end
limit means simultaneous translation of both normal arguments to that end followed by the translation-invariant
asymptotic symbol at fixed `(ω,k_∥,k_n)`:

```text
𝓛(y;ω,k_∥) = [ L_TT(y)   K_TH(y) ] ,
              [ K_HT(y)   L_HH(y) ]

K₀  = the complete (η⁰,σ_W⁰) off-diagonal block at the uniform reference background ,
K₋  = lim_{ξ→−∞} K(ξ) at the left constant background ,
K₊  = lim_{ξ→+∞} K(ξ) at the right constant background ,

𝓛₋^full = lim_{ξ→−∞}𝓛(ξ) = [L_TT,−  K_TH,−; K_HT,−  L_HH,−] ,
𝓛₊^full = lim_{ξ→+∞}𝓛(ξ) = [L_TT,+  K_TH,+; K_HT,+  L_HH,+] .
```

`K₀`, `K₋`, and `K₊`, in both directions, are **computed and emitted objects**. Their values are not supplied here.
In particular, no sector-decoupling value is inherited from the withdrawn F. The incoming, outgoing, evanescent,
and threshold channels at each end are the left/right modes of the **full** pencils `𝓛₋^full` and `𝓛₊^full`, with the
outgoing/physical-sheet prescription and the modal current of §3a. Bare diagonal `L_TT,±` and `L_HH,±` define the
asymptotic channels only under a separately recorded conditional reduction based on the computed off-diagonal
baselines.

At fixed `(ω,k_∥)`, construct the complete left/right scattering problem for the retained variable-coefficient
operator between those two full asymptotic channel spaces. The continuum mixing response is its first-background-
grade expansion with one retained off-diagonal insertion, including every term induced by expansion of the
resolvent, the asymptotic modes, and the kernel. Schematically, for a homotopy parameter used only to expose the
logic,

```text
G = G₀ + λG₁ + … ,   K = K₀ + λK₁ + … ,   ψ = ψ₀ + λψ₁ + … ,
GKψ = G₀K₀ψ₀ + λ(G₁K₀ψ₀ + G₀K₁ψ₀ + G₀K₀ψ₁) + … .
```

Consequently, a uniform-mode matrix element of `K₁` is not identified with the distorted-wave result merely from
power counting. The uniform-mode/diagonal-sector simplification is permitted only when the relevant computed
`K₀`, `K₋`, and `K₊` baseline disposition and the regular-domain hypotheses make that reduction valid; the engines
emit the premises and the reduced/full residual rather than an expected value.

Scattering labels are frequency, conserved `k_∥`, end, direction, and asymptotic channel. They neither require nor
imply a global `ω(k)` for the profile (`N5`). Re-expand every **continuum** response to first order in each of `η` and
`σ_W`; keeping profile-dependent solutions unexpanded would be a partial resummation rather than a controlled
higher-order prediction. The continuum Born domain excludes thresholds, resonance enhancement, modal-gap closures,
and long coherent regions in which repeated conversion accumulates.

**Bound-pole exemption and its limitation.** The pole solve of §3b uses the resolvent of the retained
first-shape-order imported operator and is exempt from the continuum re-expansion because re-expansion cannot create a
pole. This exemption does not improve the parent-theory accuracy: omitted `O(η²,σ_W²)` operator terms are the same
order as a weak binding energy `E∼λ²`. Therefore pole existence and location are, by default, spectrum of the
**truncated operator**, not a controlled parent-theory channel. Promotion to a physical claim requires the
threshold/separation and remainder conditions stated in §3b.

No overall insertion sign is supplied: each engine derives it from its own consistent Fourier,
Lippmann–Schwinger, and outgoing-resolvent convention, while the comparator joins the resulting canonical amplitudes
and currents (§3a/§7).

---

## 3 · The construction (OUTPUTS)

Every object below is computed for both anchorings `α` and both density representatives `ρ`; it carries its computed
`(ε,η,σ_W)` multigrade and restored `[L,T,M]` dimension. Apart from the supplied power-counting contract, this
document supplies no component value, sign, parity, grade, or expected residual.

### 3a · The profile-conditioned mixing response (continuum channel) — a CANONICAL S-matrix object

Emit the mixing response of §2 for every `(α,ρ)` as the **complete left/right open-channel S-matrix** at fixed
`(ω,k_∥)`:

```text
S_{(e_out,b)←(e_in,a)}(ω,k_∥) ,
e_in,e_out ∈ {−,+} ,
a = every incoming open mode of 𝓛^full_{e_in} ,
b = every outgoing open mode of 𝓛^full_{e_out} .
```

Both incident ends and every reflected/transmitted open channel are mandatory; a one-end response is not the
canonical payload. Emit closed/evanescent modes as matching data, not flux channels. If the full asymptotic blocks
mix the inherited sectors, label transverse-like and thickness-like channels by computed spectral projectors or
continuous continuation from the reference sector basis, and emit that classification map. A classification that is
not well-defined at a degeneracy is reported as a domain limitation rather than silently replaced by bare-sector
labels.

**Modal energy current and normalization.** Derive `J` from the S11b quadratic energy current, evaluated on the
imported closed operator and including its closed/nonlocal bulk contribution; do not import or cite a c2
traction–slab-pairing EMIT tag. For each asymptotic end, compute right modes `r_a` and adjoint/left modes `l_a` of
`𝓛_e^full(ω,k_n,k_∥)`. For a simple mode use the nonlinear-pencil normalization

```text
N_{ab}^{(e)} ≡ ⟨l_a, (∂_ω𝓛_e^full) r_b⟩ ,
⟨l_a, (∂_ω𝓛_e^full) r_a⟩ = 1 ,
J_{ab}^{(e)} ≡ the polarized S11b normal energy-current bilinear 𝓙_n[l_a,r_b]
                evaluated with the imported closed operator .
```

For a degenerate mode space, emit the matrices `N^{(e)}` and `J^{(e)}` and choose a current-orthogonal channel basis.
Also emit the identity relating the derived current, `∂_{k_n}𝓛_e^full`, and `∂_ω𝓛_e^full` for the engine's Fourier
convention. A bare group velocity or a typed `√(v_out/v_in)` factor is not a substitute for this multi-component,
frequency-dependent, possibly non-Hermitian current. Flux-normalize propagating modes with their computed signed
currents, and emit the field-normalized amplitudes and normalization maps separately.

The converted continuum flux uses the **thickness-channel current**:

```text
s₋ ≡ −1 ,   s₊ ≡ +1                                      (outward end orientations) ,
J_T,in[a]    ≡ −s_{e_in} 𝓙_n,T,in[a] ,
J_H,out[a]   ≡ Σ_{e∈{−,+}} s_e 𝓙_n,H,e,out[S_{H,e←T}a] ,
C_{T→H}[a]   ≡ J_H,out[a] / J_T,in[a] .
```

This definition never calls the converted current “transverse.” Cross terms inside a non-diagonal current block are
retained before choosing the current-orthogonal basis; the end-orientation factors, rather than an absolute-value
shortcut, define incoming and outgoing flux.

```text
⇒ S11CD_ASYMPTOTIC_FULL_OPERATORS_AND_BASELINES (𝓛₋^full,𝓛₊^full,K₀,K₋,K₊,modes,classifiers) ,
  S11CD_MODAL_FLUX_BILINEAR (left/right modes,∂_ω𝓛,current matrices,normalization maps) ,
  S11CD_COMPLETE_CHANNEL_S_MATRIX (both incident ends; field- and flux-normalized forms) ,
  S11CD_CONVERSION_AMPLITUDE , S11CD_CONTINUUM_T_TO_H_FLUX_FUNCTIONAL .
```

### 3b · The two DISTINCT photon-kill channels (`N13`) — the bound pole is a PROFILE-FUNCTIONAL CONDITIONAL

`N13`: "confinement of light" = **survival of the transverse polarization channel**. Conversion into a **bound**
breathing/thickness mode kills the photon **exactly as** bulk radiation does — the two are **distinct emitted
objects**, ⛔ not one "energy stays in the slab" statement.

- **(i) continuum conversion** — transverse → the thickness **continuum** / bulk escape (the §3a amplitude projected on
  radiating/continuum channels and evaluated with each converted channel's own current).
- **(ii) bound-mode capture — a PROFILE-FUNCTIONAL, COMPUTED CONDITIONAL, ⛔ no class-wide existence claim.** A class
  with a fixed asymptotic thickness jump contains both monotone steps **and** profiles with localized
  overshoots/wells, which
  have **different** pole sets — so existence/absence ⛔ cannot be asserted class-wide. ⛔ Do **not** invoke the 1D
  weak-well theorem (it needs equal asymptotes + an attractive self-adjoint well; the interface has different
  asymptotes and the thickness operator is multi-component, `ω`-dependent, nonlocal, possibly **non-Hermitian** from
  the outgoing bulk response).

The canonical bound output is **not** a determinant value. For the retained operator pencil `𝓛(ω)`, emit its pole
set and, for each isolated pole `ω_*`, the normalized resolvent residue and Riesz projector:

```text
𝓡_* ≡ Res_{ω=ω_*} 𝓛(ω)⁻¹ ,
𝓟_* ≡ (1/2πi) ∮_{Γ_*} 𝓛(ω)⁻¹(∂_ω𝓛(ω)) dω ,
𝓛(ω_*)r_* = 0 ,   l_*†𝓛(ω_*) = 0 ,   ⟨l_*,(∂_ω𝓛)(ω_*)r_*⟩ = 1
```

with the corresponding generalized finite-rank form for multiple poles. Emit the spectral coupling/overlap obtained
by applying `𝓡_*` or `𝓟_*` to the transverse incident source. A Jost/Evans or Fredholm determinant may be emitted only
as a noncanonical diagnostic: an Evans/Jost function is defined up to a nonvanishing analytic factor, and a
determinant for the nonlocal operator additionally requires explicit analytic-Fredholm/compactness or trace-class
premises plus a fixed normalization. Its literal value is not a comparator key.

A true bound pole is on the physical sheet, normalizable, zero-width, and has every radiation channel closed;
“below both continua” alone does not decide this for the non-Hermitian operator. Emit the sheet, normalizability,
width, and all-channel-closure tests separately and distinguish second-sheet resonances. Because no concrete
preparation/switching/damping protocol is supplied in this spec, S11c-d emits **spectral overlap only**, not a capture
probability or rate. A later capture observable requires such a protocol and a same-dimension probability before it
can be combined with continuum loss.

**Truncated-model status and promotion domain.** Every pole just described is first reported as a pole of the
first-shape-order imported operator. The omitted `O(η²,σ_W²)` terms compete with weak binding `E∼λ²`, so neither a
weak pole nor an empty weak-pole set is automatically a controlled parent-theory statement. Promotion requires a
closed contour `Γ_*` separated from every asymptotic threshold, branch point/cut, other pole, and modal-gap closure,
together with a bound on the omitted operator remainder `R_{≥2}` strong enough that
`sup_{ω∈Γ_*}||𝓛_ret(ω)⁻¹R_{≥2}(ω)||<1` (or an equivalent analytic perturbation bound). That condition preserves the
enclosed Riesz rank and controls pole motion. Without it—or a higher-order/full-operator validation—the emitted pole
set and residues remain truncated-model data. Any localized pole/resonance is not a Bloch band; the profile is
nonperiodic independently of the pole-set disposition.

**Transverse survival/confinement functional.** For each transverse incident vector `a_T`, use both reflected and
transmitted transverse blocks of the complete S-matrix:

```text
J_T,out[a_T] = J_T,out,−[ S_{T,−←T}a_T ] + J_T,out,+[ S_{T,+←T}a_T ] ,
P_T,surv[a_T] ≡ J_T,out[a_T] / J_T,in[a_T] .
```

Each `J_T,out,e` is the §3a current quadratic form restricted by the computed transverse classifier. This is the
`N13` confinement object; continuum `T→H` conversion and bound spectral overlap are emitted beside it, not substituted
for it. The spec supplies no expected disposition of survival.

```text
⇒ S11CD_CONTINUUM_CONVERSION ,
  S11CD_BOUND_POLE_SET_AND_RIESZ_DATA (poles,residues,projectors,sheet/closure tests; may be empty) ,
  S11CD_BOUND_SPECTRAL_OVERLAP , S11CD_TRANSVERSE_SURVIVAL_FUNCTIONAL .
```

### 3c · The order bookkeeping — amplitude components, physical flux, and induced-field quadratic form

Every object is multigraded `(ε,η,σ_W)` from its actual data dependency (S11c-a `:195`). Retain `ε¹` and first order
in each independent background bookkeeper, including the mixed retained grade. `s=Q_nL_W` and `k_aL_W` are kinematic
parameters, not grades.

For each converted channel, apply the following decomposition to the **field-normalized** outgoing amplitude vector
before contraction with the current form; emit its flux-normalized image through the §3a normalization map as a
separate representation:

```text
A_H = A_00 + A_η + A_σ + A_ησ
    = ε[a_00 + ηa_10 + σ_Wa_01 + ησ_Wa_11] ,

A₀              ≡ A_00                         (uniform-reference amplitude) ,
A_zero-jet       ≡ A_η                          (zero-jet contrast component) ,
A_first-jet      ≡ A_σ + A_ησ                   (components containing σ_W) ,
ΔA               ≡ A_H − A₀ = A_η+A_σ+A_ησ .
```

Thus `ΔA` is the field induced relative to the uniform reference, but it is **not** synonymous with the
gradient-sourced field: it also contains the zero-jet contrast component `A_η`. `A₀`, every component, and their sum
are computed; no baseline value is supplied from the withdrawn F.

Use the named physical path `λ≡η` at fixed `L_W` and fixed shapes, so `σ_W=(W̄₀/L_W)λ` tracks. This path does not
identify the formal `η` and `σ_W` grades. Write and emit

```text
A_H(λ) = ε[a₀ + λa₁ + λ²a₂]  at the retained rectangle ,
B_H(λ) = B₀ + λB₁ + λ²B₂ + … ,
J_T,in(λ) = ε²[j₀ + λj₁ + λ²j₂ + …] ,
```

where `B_H(λ)` is the outgoing thickness-channel current form inherited from §3a, including the dependence of the
asymptotic modes and current normalization on the background. The physical continuum conversion flux is
`J_H,out(λ)=B_H(λ)[A_H(λ),A_H(λ)]`, with explicitly emitted coefficients including

```text
J_H^(0) = ε² B₀[a₀,a₀] ,
J_H^(1) = ε²{B₀[a₀,a₁]+B₀[a₁,a₀]+B₁[a₀,a₀]} ,
J_H^(2) = ε²{B₀[a₁,a₁]+B₀[a₀,a₂]+B₀[a₂,a₀]
             +B₁[a₀,a₁]+B₁[a₁,a₀]+B₂[a₀,a₀]} .
```

These are coefficients of the retained rectangular truncation. If the computed baseline `a₀` participates, omitted
pure second-order amplitude/current terms can also enter the parent-theory `λ²` flux through baseline interference;
emit that truncation status rather than promoting the displayed retained coefficient. Under a computed disposition
that removes the relevant baseline/interference slots, the leading induced-field `B₀[a₁,a₁]` term does not require an
uncomputed second-order amplitude.

Accordingly, the physical flux expansion contains a baseline `O(ε²)` slot and an interference/current-variation
`O(ε²λ)` slot before its `O(ε²λ²)` slot. The physical conversion fraction is

```text
C_{T→H,total}(λ) ≡ J_H,out(λ)/J_T,in(λ) ,
```

with the quotient expanded consistently; it is not assigned a leading `λ` order until the computed baseline and
interference disposition permits one.

Separately emit the coherently subtracted **induced-field quadratic form**

```text
Q_H,induced(λ) ≡ B_H(λ)[ΔA(λ),ΔA(λ)] ,
C_H,induced-field(λ) ≡ Q_H,induced(λ)/J_T,in(λ) .
```

This is a useful quadratic diagnostic with the homotopy slots `O(ε²λ²)` and `O(λ²)`, respectively, but it is neither
the physical total conversion fraction nor, in general, `J_H,out[A_H]−J_H,out[A₀]`; the latter contains coherent
interference and variation of the current form. The `N12` physical labels “`O(εη)` coupling,
`O(ε²η²)` leakage, `O(η²)` fraction” are attached to the physical `T→H` observable only conditionally on the
computed baseline/interference disposition. The formal component grades and the induced-field quadratic labels are
emitted regardless. An `O(ε²λ²)` quadratic observable of a linear vertex is not the excluded nonlinear-light
program.

Bound spectral overlap is not inserted into any continuum current. A stationary bound state has no asymptotic
outgoing current, and this spec supplies no capture protocol; total photon-loss probability is therefore not formed
here. Total mechanical energy and bare bulk Poynting are not substitutes for the channel-resolved currents.

```text
⇒ S11CD_AMPLITUDE_COMPONENTS_AND_MULTIGRADE (A₀,A_zero-jet,A_first-jet,ΔA) ,
  S11CD_TOTAL_CONVERSION_FLUX_EXPANSION (baseline,interference,quadratic slots) ,
  S11CD_TOTAL_T_TO_H_CONVERSION_FRACTION ,
  S11CD_INDUCED_FIELD_QUADRATIC_FORM , S11CD_N12_BASELINE_INTERFERENCE_OPERANDS .
```

### 3d · The strong-edge bridge — NAMED as a downstream obligation, ⛔ NOT solved, ⛔ NOT a reduction

S11c-d establishes weak coefficients of the retained operator; it does not establish the order-unity edge response.
A lab slit edge is localized and order-unity in contrast. Emit separately

```text
∂_λ(ΔA_H/ε)|₀ ,
C_{T→H,total}(0) ,   ∂_λC_{T→H,total}|₀ ,   ½∂²_λC_{T→H,total}|₀ ,
lim_{λ→0} C_H,induced-field(λ)/λ² .
```

The total-fraction Taylor coefficients retain the baseline/interference content of §3c. The identification
`lim C_{T→H,total}/λ² = ½C''_{T→H,total}(0)` is made only under the computed baseline/interference disposition that
makes the quotient regular with those lower slots absent. An amplitude coefficient, a total-fraction coefficient,
and the induced-field quadratic coefficient are different objects. No value or sign is supplied for any of them. A
projected coefficient is evaluated with the computed form factor and may encounter a zero of that form factor, if
one exists; no universal node is prescribed. Every statement about a hypothetical finite-contrast continuation is
conditioned on the actually computed weak coefficient. The lab bounds the strong conversion fraction
`C_strong(1)`, not an amplitude evaluated at unit contrast and not a Born coefficient substituted at `η=1`.

⛔⛔ **A nonzero Born coefficient supplies NO general positive lower bound on strong-edge conversion.** Illustrative
baseline-free counterexample (not a model of the slab): a lossless two-mode coupler with dimensionless integrated
coupling `G`,

```text
A_H = −i ε sin(η G) ,     C = sin²(η G) = η²G² + O(η⁴) ,     C = 0  at finite nonzero  η G = nπ .
```

This example has the baseline-free N12 amplitude/flux orders yet its exact conversion returns to zero at finite
coupling. Therefore even a computed nonzero weak coefficient supplies no finite-contrast lower bound. Do not evaluate
the Born coefficient at `η=1` and compare it with the lab; that identifies weak data with `C_strong(1)`.

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
⇒ S11CD_WEAK_AMPLITUDE_COEFFICIENT , S11CD_TOTAL_FRACTION_TAYLOR_COEFFICIENTS ,
  S11CD_INDUCED_FIELD_WEAK_QUADRATIC_COEFFICIENT ,
  S11CD_STRONG_EDGE_OBLIGATION (named premise, ⛔ not solved here) .
```

---

## 4 · Objects to compute and emit

Per anchoring `α∈{L,M}` and density representative `ρ∈{ρ_4D,ρ_br}`, multigraded and dimensioned:

- The **mixing response** + conversion amplitude (canonical, flux-normalized, labelled) + the asymptotic operators /
  computed `K₀,K₋,K₊`, full channels, left/right modes, and modal current bilinears — §2/§3a.
- The **two photon-kill channels** (continuum conversion; the **profile-functional** bound-pole spectral test —
  canonical pole set + normalized Riesz residues/projectors + spectral overlap, possibly empty) and the explicit
  reflected-plus-transmitted transverse survival functional — §3b. No capture probability is emitted without a
  concrete protocol.
- The **leakage bookkeeping** — the full amplitude decomposition (`A₀`, zero-jet contrast, first jet, `ΔA`), physical
  total converted-flux baseline/interference/quadratic slots, total `C_{T→H}`, the induced-field quadratic form, and
  the operands needed for the conditional `N12` disposition — §3c.
- The **weak amplitude, total-fraction Taylor, and induced-field quadratic coefficients**, and the named strong-edge
  obligation — §3d.
- The **control outputs** of §5, each emitted as the object and its literal residual.
- **Profile moments / form factors** derived from the imported kernel, their reduced one-dimensional normalization,
  and reconstruction of every c2 three-dimensional carrier (`N15` data, ⛔ no new constitutive constants).

Every result carries its `(ε,η,σ_W)` order (and, on the homotopy, its `λ`-order) and its restored `[L,T,M]` dimension.
⛔ No result is reported without both.

---

## 5 · Independent routes and controls

⭐ Every control re-enters the chain **at the ACTION / the imported operands**, ⛔ never at a result. Each emits the
object and its literal residual; ⛔ none asserts a target value. A **coefficient** rescale tests arithmetic; only a
**form** change tests physics.

### 5a · Scattering-coordinate covariance regression — downstream only; c2 N6 debt remains open

The independent kernel-level Eulerian↔material shape construction required by `N6`, including independent N3/N4
content, belonged to c2 and remains the unclosed debt of §1b. S11c-d imports only the already-constructed closed
operator/kernel: it has neither native pre-extraction Eulerian/material operands nor term-origin provenance that
isolates tilt from advection. Rewriting that same kernel in another chart is therefore **not** an independent
derivation of its shape content and does not discharge kernel-level N6, N3, or N4.

S11c-d nevertheless performs a downstream **scattering-coordinate covariance regression** at fixed anchoring `α`
and density representative `ρ`. The actual chart is the in-plane material/Eulerian map

```text
x^i = X^i + u^i(X,t)
```

with the anchoring held fixed. It is not a flattening of the slab interface faces. Route E constructs the §2–§3
scattering problem directly in Eulerian `x` coordinates. Route M rewrites the same imported closed operator in `X`,
including the chart Jacobian, covector/mode maps, conserved tangential measure, and the S11b current, and maps the
scattering data **inside the construction** back to the common Eulerian channel basis. Emit

```text
S11CD_SCATTERING_COORDINATE_COVARIANCE_ROUTE_E[α,ρ] ,
S11CD_SCATTERING_COORDINATE_COVARIANCE_ROUTE_M[α,ρ] ,
S11CD_SCATTERING_COORDINATE_COVARIANCE_RESIDUAL[α,ρ]
    ≡ ROUTE_E − ROUTE_M  in the common Eulerian channel/current basis .
```

No `Φ` acts on amplitudes, modes, the S-matrix, or flux. `Φ` is c2's constitutive source map and is not a map on
scattering data. There is no final `routeE−Φ(routeM)` operation. Any covariance map emitted here is one actually
derived from `x=X+u` for incoming data, measures, modes, and the current.

The available one-sided mutations are **shape-sensitivity probes**, not clean channel-isolation or independence
tests, because the imported kernel exports no term-origin rows and a common `w₁′` factor can move several channels.
At fixed `(α,ρ)`, emit baseline, mutated operand, and residual for:

- reversal of the explicit thickness first-jet factor `w₁′`/`∇w₁` on one coordinate route only;
- omission or reversal of an explicitly identifiable `u·∇ρ_4D,bg⁰/ρ_4D,bg⁰` factor on one route only.

For `RHO4_CONSTANT`, emit the computed structural absence of the latter factor, with its source expression and
density-gradient operand; do not manufacture an `A−A` residual. For `RHOBR_CONSTANT`, keep the varying
`ρ_4D,bg⁰=rho_br/W_bg` and its gradient live. Neither mutation is labelled “tilt isolated” or “advection isolated.”

`∇W_bg→0` and `η→0` are rejected as one-sided mutations because they leave the interface family and create the
vacuous uniform regression. Changing one physical anchoring is also not a coordinate test, and `Δρ` does not bridge
`LAB_HELD` to `MATERIAL_ADVECTED`. Each covariance and shape-sensitivity operand/residual is printed before any guard;
its value is not a builder target. The §1b carrier/source/Φ debt and kernel-level N6/N3/N4 independence remain
explicitly unclosed whatever this downstream regression prints.

```text
⇒ S11CD_SCATTERING_COORDINATE_COVARIANCE_{ROUTE_E,ROUTE_M,RESIDUAL}[α,ρ] ,
  S11CD_SHAPE_SENSITIVITY_{BASE,MUTATED,RESIDUAL}[α,ρ,probe] ,
  S11CD_RHO4_ADVECTION_FACTOR_ABSENCE_OPERANDS[α] .
```

### 5b · The jet-zero regressions — two separate constant backgrounds

Setting `w₁′=m₁′=0` everywhere makes each profile a single constant and forces its two ends equal; it cannot retain
the unequal asymptotes of one interface. Therefore run three separate uniform constructions:

1. the **left constant background everywhere**, using the interface's computed left values
   `(W₋,μ₋,density-map₋)` at both ends and retaining the corresponding live zero-jet contrast;
2. the **right constant background everywhere**, using `(W₊,μ₊,density-map₊)` at both ends and retaining its live
   zero-jet contrast;
3. the uniform reference with `η=σ_W=0`.

In each construction, recompute the full operator, `K`, its asymptotic modes, currents, and S-matrix/coupling object;
do not insert a stated coupling value. These regressions expose what the imported first-shape-order operator computes
on each constant background. They do not validate a gradient coefficient, sign, or parity and do not close the
withdrawn F or the c2 N6 debt. In the common inherited field basis, emit the reference coupling as baseline, each
left/right constant coupling as an operand, and the direct residual `K_uniform,end−K_uniform,reference`; keep the
end-specific modes/currents beside those triplets rather than subtracting unlike channel spaces.

```text
⇒ S11CD_UNIFORM_REFERENCE_REGRESSION (inputs,computed coupling,modes,current) ,
  S11CD_UNIFORM_LEFT_{BASE,OPERAND,RESIDUAL} ,
  S11CD_UNIFORM_RIGHT_{BASE,OPERAND,RESIDUAL} .
```

### 5c · The profile-FORM ablation + the edge-vs-bump discriminant

- **Form ablation.** Perturb the **FORM** of `m₁(ξ)` (and `w₁(ξ)`) and emit the **baseline operand, the altered-form
  operand, and their residual**, with no prescribed residual disposition. A coefficient rescale of `η` tests
  arithmetic; the form change tests the profile dependence.
- **Thickness edge-vs-bump discriminant.** Emit the **two operands separately** — `(w₁′)̂_red(0)` (the reduced form
  factor at zero transfer) **and** `Δw₁ ≡ w₁(+∞)−w₁(−∞)` (the jump from the asymptotic limits) — and their residual;
  ⛔ do not assert the identity `(w₁′)̂_red(0)=Δw₁` as a single payload (it is a Fourier theorem, ⛔ not a scattering
  result). `Δw₁` is the class discriminator: the in-class thickness interface has unequal thickness ends, whereas a
  thickness bump returns to the same thickness. The ablation changes the thickness profile from the interface family
  to a bump family and emits the resulting zero-transfer operands/residual without an expected value.
- **Independent modulus control.** Emit the analogous **two operands** `(m₁′)̂_red(0)` and `Δm₁ ≡ m₁(+∞)−m₁(−∞)` (and
  their residual) as the per-profile discriminator of the modulus subchannel — ⛔ again not the identity as a single
  payload. A constant or modulus bump is allowed inside the thickness-interface class, so its modulus moment does not
  reclassify the thickness object. Emit the full vertex's zero-transfer dependence; do not replace it by a
  modulus-only projection or freeze one representative shape as “the slit.”

```text
⇒ S11CD_PROFILE_FORM_ABLATION (baseline,altered,residual) ,
  S11CD_THICKNESS_EDGE_BUMP_DISCRIMINANT , S11CD_MODULUS_SUBCHANNEL_DISCRIMINANT .
```

### 5d · The falsification FORM control (`N7`; the boundary refinement is N2-permitted)

Emit the **computed flux-normalized dimensionless conversion projection** from the complete channel S-matrix,
modal-current forms, and the **full** imported vertex (§3a–§3c). Do not substitute a typed expected shape such as
`∝k·a`, and do not replace the full projection by its `∇μ_R` subchannel. Emit both the total physical conversion FORM
(including its baseline/interference slots) and the separately named induced-field quadratic FORM.

The decision-list table places the FORM and confinement interpretation at S11c-e; N2 permits this spec-stage boundary
refinement. S11c-d's weak-contrast FORM does not become S11c-e's withheld order-unity target. The diffraction-grating
reductio and numeric lab bound remain orchestrator-side; the magnitude is `R1`-blocked. A slit edge is an order-unity
localized gradient, so no non-perturbative lab number is inferred from this retained operator.

---

## 6 · Method, dimensions, and script obligations

- **Method.** Balance laws + the binding material virtual-displacement rule + variational derivatives with held-fixed
  fields named + prescribed external virtual work (S11b), ⛔ never an irreversible response kernel in an ordinary
  action. The closed operator and its off-diagonal extraction are the c2 exports consumed verbatim. Compute
  `K₀,K₋,K₊`; define modes and currents from the full asymptotic block pencils; then construct the complete two-ended
  scattering matrix and re-expand its continuum response to the retained background rectangle (§2). The bound-pole
  solve is exempt only in the limited truncated-model sense of §2/§3b.
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
  independently executed second construction exists. Section 5a's two chart implementations are explicitly a
  downstream covariance regression, not an independent derivation of the imported kernel's N3/N4 content.
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
load-bearing residuals include the complete flux-normalized S-matrix/mixing amplitudes, modal currents, total and
induced-field conversion forms, survival functional, and bound Riesz data (§3). These d objects **remain conditional
on and PROPAGATE** the material c2 operand debt (§1b). A final projected amplitude residual cannot directly surface
the carrier (40), source (76), or `Φ` (18) families: the c2 disposition identifies a schema non-join, and projection
or channel summation can cancel operand differences. Agreement on a d projection therefore compares only that d
projection and does not close any upstream operand family.

Join channel records by `(α,ρ,profile,ω,k_∥,incident end/channel,outgoing end/channel)` after the common current
normalization. Join bound records by profile, sheet, isolating contour, and algebraic-multiplicity label, comparing
pole locations and normalized Riesz data rather than any determinant value. These are schema rules, not expected
physics values.

Direct surfacing may be claimed only if d separately defines and emits common-basis projections of the carrier,
source, and `Φ` operand families, preserves their family identity before channel summation, and supplies a sound
cross-engine schema join for each. The present consume-set supplies no such c2 provenance rows, so this spec makes no
direct-surfacing claim. The full per-object symbolic residual and c1's four giant families remain deferred
(`DEFERRED_HEAVY_RUNS.md`); S11c-d names and propagates, but does not pre-adjudicate, what it cannot close on this box.

**The blind Wolfram engine** re-derives the §§1–2 supplied inputs, the S11c-a face substrate, the S11c-b slab-operator
and c2 closed-operator/kernel rows it consumes, and the localized-interface mixing — importing nothing (the only
cross-engine control). ⛔ The denylist stays cut (`N9`/rule 12); blindness is enforced by absence.

---

## 8 · Supplied versus computed; builder report

**SUPPLIED (unfalsifiable in this build):** all of §1 (the two c2 export operand rows and their per-engine-SOUND vs
cross-engine-UNCLOSED disposition — the operand DEBT, the raw/schema-unmatched `R_N6`, the two S11c-b signs, the six
§3d re-adjudications, the 3 N6 premise caveats, the withdrawn F/G, and that no term-origin/parity/increment/§3d rows
are importable); the §1c localized **thickness** interface (`Δw₁≠0`) with independent unrestricted-class `m₁`, the
branchwise density maps, short-range-jet domain, exact reduced Fourier/3-D carrier convention, and admissibility; the
§1d regime (Born in `η`, with `σ_W`, `Q_nL_W`, and local `kL_W`/gap kinematics separate and live); the §2 requirement
to compute the uniform/end baselines and define channels from full asymptotic block pencils; and `N11a`, `N12`,
`N13`. No baseline value is supplied.

**COMPUTED (outputs, with no expected result supplied):** `K₀,K₋,K₊`, the two full asymptotic operators, complete
left/right modes/classifiers, explicit S11b-derived modal currents, and the complete S-matrix (§2/§3a); the continuum
`T→H` current functional, total transverse survival functional, and the profile-functional truncated-operator pole
set + normalized Riesz residues/projectors + spectral overlaps (§3a–§3b); the separated zero-jet/first-jet amplitude
components, total physical flux baseline/interference/quadratic slots, total conversion fraction, induced-field
quadratic form, and operands for the conditional `N12` disposition (§3c); the weak Taylor coefficients and
strong-edge handoff (§3d);
the downstream coordinate-covariance and shape-sensitivity operands/residuals (§5a), three uniform regressions (§5b),
profile-form controls and thickness/modulus discriminants (§5c), and the computed falsification projection (§5d);
all derived profile moments/form factors and their 3-D reconstruction; and every output's `(ε,η,σ_W)`/`λ` order and
`[L,T,M]` dimension.

**Builder report.** The build directive states, per emitted object, which line computed it (`.claude/skills/build`);
declares the profile class (§1c), the regime grades (§1d), the distorted-basis insertion organization (§2), and the
leakage bookkeeping (§3c) it implemented; and reports every §5 baseline, altered operand, and literal residual before
guards—never a prose conclusion. The disposition of every §5 residual, every baseline/interference slot, and the
strong-edge obligation (§3d) is read on **our** side in the step record, not asserted by the script (rule 5).
