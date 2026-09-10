# S11c-d — SymPy build directive (profile-conditioned transverse↔thickness scattering; SymPy engine)

⭐ **THIN directive.** All physics is the cleared spec `directives/S11c_d_SHARED_PHYSICS.md` (v10, **CLEARED**,
committed `399a8516`; round-10 dual-engine gate both legs SOUND). This directive ⛔ does **not** restate or re-derive
the physics (a re-wording is weaker and drifts, [[feedback_point_at_the_obligation_never_restate_it]]) — it POINTS at
the spec and fixes only the **build-mechanical** layer: the import wiring, the exact `IMPORT_KEYS` root set, the
**deferred §1c 3-D→1-D Fourier-reduction element census against the real closed rows** (§3 — the spec deferred this
here on purpose, `S11c_d_SHARED_PHYSICS.md:114-116,253-255`), the emit tags, and the three script clauses. Model:
**CODE build → `gpt-6-astra` high**.

- **Deliverable script:** `scripts/S11c_d_mixing_scattering_sympy_audit.py` (tag prefix `PY_S11CD_*`; ledger
  write-keys `s11cd*` lowerCamel).
- **Deliverable export:** `scripts/S11c_d_exports.py` (own-rows delta; §5).
- **Deliverable builder report:** `_measurements/S11c_d_sympy_builder_report.md` (the frozen import/reduction map +
  per-object which-line-computed-it + literal §5 residuals — ⛔ no verdicts).
- **Physics authority (SUPPLIED, unfalsifiable in this build):** the whole S11c-d spec `399a8516`. This directive
  governs the **SymPy** engine only; the blind Wolfram engine + the T7 comparator are separate downstream artifacts
  with their own directives.

---

## 0 · The three mandatory script clauses (verbatim — `.claude/skills/build/SKILL.md`, non-negotiable)
> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit`/`Print` payload
> must be a CAS object — an expression, a solved root, a boolean from a symbolic test. ⛔ Never prose
> describing a result.
> **2. PRINT the residual; do NOT assert it.** `assert residual == 0` **is the builder writing down the
> expected output**, and it turns an informative value into a binary crash. Compute → emit → *then* assert.
> **3. Interpretation belongs to the STEP RECORD.** ⛔ The script does not editorialise.

**Structural rule (verbatim):** *The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE
ACTION and the ANSATZ. Every other expression involving them must be REACHED BY COMPUTATION. Every control re-enters
the chain at the ACTION / the imported operands, ⛔ never at a result.*

**Corollaries (all bind):**
- **A hand-typed CAS object is still hand-typed** — delete a `Solve`/`subs`/`series` and the emit must MOVE; else it
  is a typed answer wearing algebra. The **construction operands of `𝓛` are the §1c-reduced closed rows** (§2/§3), ⛔
  never a hand-typed off-diagonal payload, K₀/K₋/K₊ value, amplitude, or baseline.
- **A tag name names the OBJECT**, ⛔ never its value/sign/order/parity/grade (spec §7 F9/N14).
- **No tautological residual.** The §3a S-matrix / §3c bookkeeping export *representations* are **not** checks. A §5
  residual is emitted with **both operands**; a **two-route** residual is emitted **only where an independently
  executed second construction exists**. ⚠ **§5a's Route E / Route M are a downstream coordinate-covariance
  regression, ⛔ NOT an independent derivation of the imported kernel's N3/N4 content** (the imported kernel exports no
  term-origin rows — `S11c_d_SHARED_PHYSICS.md:711-717,844-848`); ⛔ do not dress `routeE−routeM` as a kernel-N6
  closure. For `RHO4_CONSTANT`, emit the **computed structural absence** of the advective factor with its source
  expression + density-gradient operand — ⛔ do **not** manufacture an `A−A` residual (§5a).
- **Emission is NEVER conditional on a payload's VALUE** — only on which package/quantity (a value present-and-identical
  = INVARIANT is a result; absent = indistinguishable from never-computed).
- ⛔ Every anti-example uses **placeholder** symbols. ⛔ Run `reduction/derived_or_declared.py` on the deliverable
  after the build (triage, ⛔ not a verdict); run `reduction/engine_output_checks.py` on its `.out`.

---

## 1 · What to build — POINTERS to the cleared spec (⛔ do not restate the physics)
Build the S11c-d profile-conditioned mixing/scattering object exactly per `S11c_d_SHARED_PHYSICS.md`:
- **§0** scope (in/out); **§1a** the two consumed c2 closed rows + carriers (SUPPLIED); **§1b** the per-engine-SOUND
  vs **cross-engine-UNCLOSED** disposition — ⛔ do **not** treat the carried cross-engine operand agreement on the
  coupling kernel as closed (qualitative premise, §6 here); **§1c** the localized-interface profile class + branchwise
  density maps + the `f̂_red` reduced-transform convention (the §3 census implements the reduction it defers here);
  **§1d** the regime (Born in `η`; `σ_W`, `Q_nL_W`, `k_aL_W` **separate and live** — ⛔ no extra `σ_W→0`, ⛔ no
  Taylor-away of `s`, ⛔ no `η`→O(1), ⛔ no `L_W→0`).
- **§2** the **reduced-representation rule (GOVERNING §§3–6)** + the two-asymptote distorted-wave organization —
  construct `𝓛` from the engine's **§1c-reduced** closed slab operator with its **reduced** closed coupling kernel as
  the canonical off-diagonal extraction; compute `K₀,K₋,K₊` + the two **full** asymptotic block pencils `𝓛₋^full`,
  `𝓛₊^full`; the distorted-wave insertion (`GKψ`) logic; re-expand every **continuum** response to first order in
  each of `η` and `σ_W`; the bound-pole re-expansion **exemption** with its truncated-model status.
- **§3a** the complete two-ended S-matrix (both incident ends, every reflected/transmitted open channel) + modal flux
  bilinear derived from the **S11b quadratic energy current on the reduced operator** (⛔ do **not** import/cite a c2
  traction–slab-pairing EMIT tag; ⛔ no typed `√(v_out/v_in)`); **§3b** the two **DISTINCT** photon-kill channels
  (continuum conversion + the **PROFILE-FUNCTIONAL, computed-conditional** bound pole — ⛔ do **not** invoke the 1-D
  weak-well theorem, `N13`; emit pole set + normalized Riesz residues/projectors + sheet/normalizability/width/
  all-channel-closure tests + **spectral overlap only**, ⛔ no capture probability/rate; a determinant/Jost value is a
  noncanonical diagnostic, ⛔ not a comparator key) + the reflected-plus-transmitted **transverse survival functional**;
  **§3c** the order bookkeeping (`A₀`/zero-jet/first-jet/`ΔA` multigrade; physical flux **baseline/interference/
  quadratic** slots on the `λ≡η` homotopy; the induced-field quadratic form; the **conditional** `N12` labels);
  **§3d** the strong-edge bridge NAMED as a downstream S11c-e obligation (⛔ NOT solved, ⛔ NOT a reduction) + its
  weak coefficients.
- **§4** objects to emit; **§5** the controls (5a downstream scattering-coordinate covariance regression + one-sided
  shape-sensitivity mutations; 5b the **three** jet-zero uniform regressions; 5c profile-**FORM** ablation +
  edge-vs-bump + modulus discriminants; 5d the computed falsification **FORM**); **§6** method/dimensions/script
  obligations; **§7** names/F9/export schema/comparator contract; **§8** supplied-vs-computed + builder report.

⛔ **HELD PHYSICS the build MUST carry (spec; ⛔ do not regress — each is a cleared-round finding):**
- **The reduced-representation rule (§2, governing §§3–6):** **every** downstream object — currents, modes,
  normalization, resolvent, poles, S-matrix, conversion fraction, survival functional, profile moments/form factors,
  and the §5a routes — is built from the engine's **§1c-reduced** closed operator/kernel, ⛔ **NEVER** the unreduced
  3-D rows. The unreduced rows are operands of the §1c reduction (§3) **ONLY**. ⛔ No side reduction while `𝓛` runs on
  unreduced content.
- **The withdrawn F ⇒ compute, never type a zero.** The uniform amplitude/baseline `A_0`/`K_0` and every `K₀,K₋,K₊`
  (both directions) are **COMPUTED and EMITTED** objects (§2/§3c); ⛔ **never** typed `=0` and ⛔ never inherit a
  sector-decoupling value from the withdrawn c2 F (the withdrawn F is the c2 *increment* interpretation, ⛔ it does
  not retract S11b's uniform decoupling — but the **computation**, not the F label, settles `K_0`).
- **The FULL imported off-diagonal vertex** (tilt `∇w₁`, modulus-gradient `∇m₁`, N4 advection) — ⛔ not `∇w₁` alone,
  ⛔ not a modulus-only projection; **full** = the complete retained reduced operator.
- **The reduced-slab-operator-off-diagonal-block vs reduced-kernel residual is a SURFACED FINDING** (§2), ⛔ not
  silently overridden by the canonical designation.
- **No global `ω(k)`** (`N5`); the continuum Born domain **excludes** thresholds/resonance/modal-gap closures/long
  coherent regions; any localized pole is **not** a Bloch band.
- **The c2 operand DEBT is a LIVE PREMISE on the object built** (the mixing rides on the coupling kernel): S11c-d
  **NAMES and PROPAGATES** it, ⛔ does **not** resurface the upstream operand families or pre-adjudicate — a
  d-projection residual cannot close an upstream c2 operand family (schema non-join). §6 renders it **qualitatively**;
  the comparator propagates it.

---

## 2 · Import wiring — VERIFIED against the real 3-parent fold (build-mechanical)
Read the inherited model via the **positional** `load_model` call (spec §7; signature
`load_model(base_path, *delta_paths)`, `scripts/ledger_fold.py:102` — ⛔ **NOT** the keyword form, which TypeErrors):

```python
fold, audit = load_model("scripts/S11c_b_exports.py", "scripts/S11c_c1_exports.py", "scripts/S11c_c2_exports.py")
```

**Verified on the real files (reproduce; ⛔ do not hard-code):** base **2441** rows + c1 delta **44** + c2 delta
**70** → fold **2555** (strictly additive); `audit.overwrites == []`; all pairwise exact-key intersections empty
(base∩c1 = base∩c2 = c1∩c2 = ∅). `check_consumer(fold, IMPORT_KEYS)` (`scripts/ledger_fold.py:199`) resolves cleanly
with the two c2 closed rows as roots (closure 213 keys, 273 symbol edges, 45 dimension edges, no ambiguity).
(Measurements: `_measurements/S11c_d_sympy_build_directive_census.md`.)

**`IMPORT_KEYS` — the rule (⛔ not a fixed number to copy blindly).** The fold guard `assert_lookups_equal_manifest`
requires `IMPORT_KEYS` to equal **exactly the set of keys the build looks up by `fold[key]`** — an **undeclared**
lookup **and** a **declared-but-unused** key each fail. Fix `IMPORT_KEYS` = the build's actual direct-lookup set,
subject to:
- **MUST include the two c2 CLOSED rows as roots** — `s11cc2ClosedSlabOperator` and `s11cc2ClosedCouplingKernel`
  (verified exact casing) — the §2/§3 construction operands; their `check_consumer` closure (213 keys) pulls the
  hats, momenta, integration variables, and much of the reachable consume-set.
- **MAY additionally declare** any carrier / coefficient / constant the build reads **directly by key** (for its
  symbol, its `dimension_key`, or the §5c FORM ablation, which perturbs the profile FORM): the field carriers
  `s11cc2Fieldtheta` (⚠ **lowercase `theta`** — verified), `s11cc2FieldeW`, `s11cc2Fieldu{1,2,3}`; the profile
  coefficients `s11cc2Coefficientw1Profile`, `s11cc2Coefficientm1Profile`; the Fourier carriers
  `s11cc2FourierW1ProfileHatTransfer`, `s11cc2FourierW1ProfileJetHat{1,2,3}` (+ their `*Dimension` companions); the
  momenta `s11cc1_k_output_{1,2,3}`, `s11cc1_k_input_{1,2,3}`, `s11cc2MiddleMomentum{1,2,3}`,
  `s11cc2OutgoingNormalMomentum`; and the reachable constants (`W_0`, `mu_R`, `eta_bg`, `sigma_W`, `L_W`, `rho_m`,
  `rho_br`, `Lambda_{A,V,X}_0`, `tau_{A,V,X}`, `omega`, `c_s0`, `background_density_map`, `rho_br_bg_rho4_constant`,
  …) — **iff** the build looks them up by key.
- **MUST NOT declare** any key the build does not look up (declared-but-unused fails).
⇒ The build reports its final `IMPORT_KEYS`; the two decision legs verify it against `assert_lookups_equal_manifest`
on the **real** 3-parent fold (⛔ not by reading — run it).

**⛔ Provenance the guard will NOT catch — the directive + legs must enforce:**
- ⛔⛔ **Bind the c2 CLOSED rows, ⛔ NEVER the S11c-b OPEN `slab_operator` / `coupling_kernel`.** Both the open S11c-b
  rows (`slab_operator`, `coupling_kernel` — pressure-slot / open / pre-closure) and the c2 closed rows
  (`s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel` — `class='DERIVED'`, `step='S11c-c2'`) exist in the fold;
  `check_consumer` / `assert_lookups_equal_manifest` pass on **either** (key-existence only). Binding the open rows
  silently uses **pre-closure physics** and re-invites re-closing c2's work — the exact N14/F9 false-equal hazard
  (c2's own `s11c_c1_face_response`-vs-`face_response` precedent, `S11c_c2_sympy_build_directive.md:109-115`). S11c-d
  **consumes the closed rows verbatim** as the §3 reduction operands and ⛔ does **not** re-close the open rows.
- ⚠ **Exact carrier casing:** the real key is `s11cc2Fieldtheta` (lowercase), ⛔ not `s11cc2FieldTheta`. The two
  `Closed*` rows carry **no `*Dimension` companion** (no `dimension_key`); every carrier / hat / coefficient /
  momentum key does.
- ⛔ `dtn_kernel` is reachable through the fold (c1's two-momentum DtN kernel, carrying the 3-D momentum `DiracDelta`s
  — §3) but is ⛔ **NOT** bound as the reduction-convention operand (§3).

**Case structure (per closed row).** Each `Closed*` value is a nested `Tuple`-of-pairs association list of **4
case-entries** keyed `Tuple(Str(α), Str(ρ))`, α∈{`LAB_HELD`,`MATERIAL_ADVECTED`}, ρ∈{`RHO4_CONSTANT`,`RHOBR_CONSTANT`},
payload under `Str('VALUE')`. Slab components `{U, THETA, E_W}`; coupling = outer sector layer
`{TRANSVERSE_TO_THICKNESS, THICKNESS_TO_TRANSVERSE}` then components `{THETA, E_W, DIV_U}`. Iterate all 4 `(α,ρ)`
cases (spec §3/§4).

---

## 3 · The deferred §1c 3-D→1-D Fourier reduction — PINNED against the real closed rows (build-mechanical)
The spec **defers** the reduction MECHANICS to this directive on purpose (`S11c_d_SHARED_PHYSICS.md:114-116,240-258` —
the Fourier recipe forced 5 spec rounds precisely because it was specifying the HOW; the spec supplies only `f̂_red` +
the reduce-your-own-rows requirement + the comparator-join control, [[feedback_many_review_rounds_signals_recipe_creep]]).
This section pins the **CENSUS** — *which* convention-bearing 3-D elements the SymPy engine must reduce — against the
**real imported rows**. It pins **WHAT** to reduce; ⛔ it does **NOT** supply **HOW** (the `(2π)`/`δ²` arithmetic is
the engine's **computed** reduction — spec §1c).

**The requirement (spec §1c, verbatim intent).** The SymPy engine must **COMPUTE and EMIT** the 3-D→1-D reduction of
**every** convention-bearing 3-D element of **both** its imported closed rows — the closed operator **and** the closed
coupling kernel — applying its **own realized Fourier convention** to the §1c interface geometry `f=f(n̂·y/L_W)`,
`∂_{yᵢ}f=n̂ᵢf′(ξ)/L_W`, with **no** convention-bearing element left in a 3-D convention. The **only** supplied
convention is the reduced one-dimensional template `f̂_red(s) ≡ ∫dξ e^{−isξ}f(ξ)`, `s≡Q_nL_W`, `Q≡k_out−k_in`. ⛔ Do
**NOT** supply/type a numeric `[L_W/(2π)]` or `(2π)²L_W` map or a `δ²(Q_∥)` map; ⛔ do **NOT** bind `dtn_kernel` as
the reduction-convention operand (it is c1's two-momentum DtN kernel carried inside the closed rows, ⛔ not the
3-D→1-D reduction operand).

**The census — the convention-bearing 3-D elements ACTUALLY PRESENT in the real closed rows (⛔ reduce every one;
this enumerates WHAT, ⛔ never the reduction result).** Verified per row against `s11cc2ClosedSlabOperator` and
`s11cc2ClosedCouplingKernel` (measurements file):

- **Profile hats — all four c2 hats, as 3-argument APPLIED functions** (the c1-level snake_case hats
  `s11cc1_w1_profile_hat_transfer` / `s11cc1_w1_profile_jet_hat_*` are **ABSENT** — 0 — in both closed rows):
  `s11cc2FourierW1ProfileHatTransfer` (occurrences slab/coupling ≈ 88/550), `s11cc2FourierW1ProfileJetHat{1,2,3}`
  (≈ 100/136 each). **Two argument shapes** (identical in both rows), **both to be reduced**:
  - **transfer argument** `(k_out_i − k_in_i)` (the momentum-transfer `Q_i`);
  - **middle-leg (intermediate-momentum) arguments** `(k_out_i − k_mid_i)` and `(k_mid_i − k_in_i)`, with
    `k_mid = s11cc2MiddleMomentum{1,2,3}` (≈ 24 middle-leg occurrences each).
  ⇒ reduce **every hat at every argument** (transfer AND both middle-leg forms) — ⛔ do not reduce the transfer hats
  and leave the middle-leg hats in a 3-D convention.
- **3-D integral measures — explicit `sp.Integral` nodes** (⛔ not `DiracDelta`-represented), all limits `−∞..∞`, over
  - **`d³y` real-space position integrals** — variables `s11cc2Y{1,2,3}`; **and**
  - **`d³k` momentum-convolution integrals** — over `s11cc1_k_output_{1,2,3}`, `s11cc1_k_input_{1,2,3}`, **and**
    `s11cc2MiddleMomentum{1,2,3}` (output, input, and middle momenta).
  (Literal `Integral(` counts slab/coupling ≈ 218/526; distinct `.atoms(Integral)` ≈ 26/90.) ⇒ reduce **every** `d³y`
  and `d³k` measure; ⛔ leave none in a 3-D convention.
- **No 3-D momentum `DiracDelta` in either closed row** (verified 0 by literal count and `.atoms(DiracDelta)`) — c2
  stripped them. The 3-D momentum deltas live on `dtn_kernel` (3 `DiracDelta` at `(k_output_i − k_input_i)`, one per
  component, 0 `Integral`), which is ⛔ **not** the reduction operand. So the SymPy reduction operates on **applied
  hats + explicit `Integral` measures**, ⛔ not on delta-sharpened momenta.
- The normal-momentum carrier `s11cc2OutgoingNormalMomentum` is present; carry/reduce it per the §1c edge-normal
  geometry alongside the other momenta.

⚠ This census is fixed against the **committed** c2 rows (digests in `BUILD_INPUT_DIGESTS`, §5); the two decision legs
re-run it on the real fold and confirm the build reduces **every** enumerated element with **no** convention-bearing
element left in a 3-D convention.

**Both-operand emission (⛔ not an `A−A` tautology).** For **each** convention-bearing element of **each** row, emit a
reduction record with **BOTH** operands: **(i)** the 3-D element **as it actually appears** in that engine's own
closed operator/kernel, and **(ii)** the corresponding **reduced** one-dimensional object obtained by applying that
engine's own realized convention. ⛔ Do **NOT** manufacture the 3-D operand by defining `A_3D ≡ [L_W/(2π)]δ²(Q_∥)A_edge`
and subtracting it from itself (`A−A`, spec §1c). Each engine additionally computes and emits its own 3-D↔1-D
**reconstruction round-trip** for both rows from those exposed operand pairs.

**Flux `δ²(Q_∥)` stripping (distinct step).** When forming flux/rate, ⛔ never square an unstripped tangential
`δ²(Q_∥)`: remove that factor before the per-unit-edge-area limit using the engine's **computed** reduction + the
standard tangential box/continuum normalization. This flux step is **distinct from and ⛔ not a substitute for** the
computed carrier reduction above.

**Comparator note (context, ⛔ not this build's job).** The downstream T7 comparator joins the two engines' **reduced**
operators and reduced kernels (sourced from each row's both-operand reduction record, covering every
convention-bearing element) **and** the `𝓛` / S-matrix / currents / conversion forms built from them — ⛔ not a
pre-factored coefficient, ⛔ not a side reduction while `𝓛` uses unreduced content (spec §2/§7). The SymPy build's job
is to EMIT the reduction records + the reduced construction; it does **not** run the comparator.

---

## 4 · Objects to emit (spec §3/§4) + provenance
Per anchoring `α∈{LAB_HELD,MATERIAL_ADVECTED}` and density representative `ρ∈{RHO4_CONSTANT,RHOBR_CONSTANT}`, each
object carrying its computed `(ε,η,σ_W)` (and, on the `λ≡η` homotopy, its `λ`-) order **and** its restored `[L,T,M]`
dimension (⛔ no object reported without both — [[feedback_dimensional_consistency_check]]):

- **§3a:** `S11CD_ASYMPTOTIC_FULL_OPERATORS_AND_BASELINES` (`𝓛₋^full,𝓛₊^full,K₀,K₋,K₊`, modes, classifiers),
  `S11CD_MODAL_FLUX_BILINEAR` (left/right modes, `∂_ω𝓛`, current matrices, normalization maps),
  `S11CD_COMPLETE_CHANNEL_S_MATRIX` (both incident ends; field- and flux-normalized), `S11CD_CONVERSION_AMPLITUDE`,
  `S11CD_CONTINUUM_T_TO_H_FLUX_FUNCTIONAL`.
- **§3b:** `S11CD_CONTINUUM_CONVERSION`, `S11CD_BOUND_POLE_SET_AND_RIESZ_DATA` (poles, residues, projectors,
  sheet/closure tests; **may be empty**), `S11CD_BOUND_SPECTRAL_OVERLAP`, `S11CD_TRANSVERSE_SURVIVAL_FUNCTIONAL`.
- **§3c:** `S11CD_AMPLITUDE_COMPONENTS_AND_MULTIGRADE` (`A₀,A_zero-jet,A_first-jet,ΔA`),
  `S11CD_TOTAL_CONVERSION_FLUX_EXPANSION` (baseline/interference/quadratic slots),
  `S11CD_TOTAL_T_TO_H_CONVERSION_FRACTION`, `S11CD_INDUCED_FIELD_QUADRATIC_FORM`,
  `S11CD_N12_BASELINE_INTERFERENCE_OPERANDS`.
- **§3d:** `S11CD_WEAK_AMPLITUDE_COEFFICIENT`, `S11CD_TOTAL_FRACTION_TAYLOR_COEFFICIENTS`,
  `S11CD_INDUCED_FIELD_WEAK_QUADRATIC_COEFFICIENT`, `S11CD_STRONG_EDGE_OBLIGATION` (named premise, ⛔ not solved).
- **§5 controls** (each = the object **and** its literal residual, both operands): `S11CD_SCATTERING_COORDINATE_
  COVARIANCE_{ROUTE_E,ROUTE_M,RESIDUAL}[α,ρ]`, `S11CD_SHAPE_SENSITIVITY_{BASE,MUTATED,RESIDUAL}[α,ρ,probe]`,
  `S11CD_RHO4_ADVECTION_FACTOR_ABSENCE_OPERANDS[α]` (5a); `S11CD_UNIFORM_REFERENCE_REGRESSION`,
  `S11CD_UNIFORM_LEFT_{BASE,OPERAND,RESIDUAL}`, `S11CD_UNIFORM_RIGHT_{BASE,OPERAND,RESIDUAL}` (5b);
  `S11CD_PROFILE_FORM_ABLATION` (baseline, altered, residual), `S11CD_THICKNESS_EDGE_BUMP_DISCRIMINANT`,
  `S11CD_MODULUS_SUBCHANNEL_DISCRIMINANT` (5c — ⛔ emit `(w₁′)̂_red(0)` and `Δw₁` as **two separate operands** + their
  residual, ⛔ never the identity as a single payload); the computed flux-normalized dimensionless conversion **FORM**
  projection (5d — from the complete S-matrix + modal currents + the **FULL** imported vertex; ⛔ do **not** substitute
  any typed expected shape and ⛔ do not replace the full projection by a single-gradient subchannel — spec §5d).
- **§3 reduction records:** the both-operand reduction record + reconstruction round-trip per convention-bearing
  element per row (§3 above).
- **Profile moments / form factors** derived from the imported kernel + their reduced 1-D normalization (`N15` data,
  ⛔ no new constitutive constants).

**F9 / N14:** every new object gets a **fresh injective** `mechanical_lower_camel` write-key; ⛔ never reuse an
imported key (`s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel`, `slab_operator`, `coupling_kernel`,
`dtn_kernel`, `mu_theta_operator`, `W_0`, `mu_R`, `e_W`, `rho_br`, `v_0`, …). Per emitted object the builder report
states **which line computed it**.

---

## 5 · Export schema (spec §7 / `export_ledger_bind_closure_design.md` §D1–D3)
Write `scripts/S11c_d_exports.py` as an **own-rows delta** (⛔ not the accumulated whole-model file). Membership = the
bind-closure (D1); `assert_delta_is_minimal` requires the delta's key-set = S11c-d's own bind-closure ∪ infra.
`BUILD_INPUT_DIGESTS` pins `{this SymPy audit, scripts/S11c_b_exports.py, scripts/S11c_c1_exports.py,
scripts/S11c_c2_exports.py, this spec, scripts/ledger_fold.py}` (§D3). ⛔ Never `git add -f` a big `.out`; ⛔ never
annex an `*_exports.py`.

⭐⭐ **EMIT ≠ EXPORT — export ONLY what S11c-e binds (D1); everything else is EMIT-only (→ `.out`).** §4 lists what to
**emit** (PRINT to stdout, for review + the T7 comparator). The **export delta** is far smaller: per
`directives/S11c_decisions.md:52`, S11c-d hands S11c-e **scattering amplitudes / resonances / local spectrum**, and
S11c-e's declared scope is the **flux-normalized dimensionless conversion observable + leakage + confinement** whose
weak limit must reproduce S11c-d (`S11c_d_SHARED_PHYSICS.md:60-65`). ⇒ **EXPORT the objects S11c-e's leakage
observable + strong-edge weak limit binds** — the mixing response / complete channel S-matrix / conversion amplitude,
the continuum `T→H` flux functional + transverse survival functional (the `N13` confinement object), and the §3d
weak coefficients that the finite-contrast response must reproduce — plus only their recursive new
coordinate/function/dimension bind-closure. ⚠ There is **no S11c-e manifest yet** (e unbuilt), so this membership is
grounded in e's **DECLARED scope**, ⛔ not a verified bind; S11c-e's actual `IMPORT_KEYS` confirms it when built (a
S11c-d export e never binds is the flag then; D1). The two decision legs settle the exact export membership against
that declared scope.

⛔ **EMIT-ONLY (→ `.out`, ⛔ NOT the ledger export):** the §3 reduction records/reconstruction round-trips, the bound
Riesz data, the amplitude-component/flux-slot bookkeeping, and **every §5 control operand/residual** (comparison/emit
representations — the T7 reads them from stdout; ⛔ nothing downstream binds them).

⭐ **Store the exported objects in a TRANSPARENT compact encoding** — an ordinary algebraically-equivalent factored
SymPy expression (`sp.factor`/`collect`/CSE), ⛔ **NOT `sp.expand`ed** and ⛔ **NOT an opaque `UnevaluatedExpr`/hold**
(a downstream `diff`/pole-solve must stay evaluable). PRINT the expanded form to `.out` only; ⚠ **require a casewise
semantic-equivalence check** `canonicalize(expanded_emitted_root − decode(compact_export_root)) == 0` per case (the
serialization round-trip proves storage identity, ⛔ NOT equivalence after a new compaction — D4).

---

## 6 · Supplied / withheld (leak discipline) — the c2 DEBT rendered QUALITATIVELY
Everything in spec §1 is **SUPPLIED and unfalsifiable in this build** — state so in the builder report so a passing
build does not read as if it verified it. The **outputs** are §§3–5; ⛔ this directive and the script state **no**
component value, sign, order, parity, grade, baseline, or expected residual, and carry **no acceptance criterion
referencing an expected value** (the diff happens on **OUR** side; a genuine cross-engine disagreement is a
**finding**, ⛔ not a build failure — [[feedback_no_commit_before_legs_report]]).

⚠ **The c2 cross-engine operand DEBT is rendered QUALITATIVELY (⛔ no counts to the builder).** The two consumed
closed rows (`s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel`) are per-engine-SOUND but their **cross-engine
operand agreement is UNCLOSED** — a **material, supplied-unfalsifiable** premise on the object S11c-d builds (the
mixing rides on the coupling kernel). The build **consumes them verbatim** as the §3 reduction operands, **NAMES**
the debt as unclosed, and **PROPAGATES** it; ⛔ it does **NOT** resurface the upstream carrier/source/`Φ` operand
families and ⛔ does **NOT** pre-adjudicate them (a d-projection residual cannot close an upstream c2 operand family).
⛔⛔ **Do NOT present any c2-import cross-engine status quantity as an S11c-d target** — the c2 import's operand /
covariance / naturality status quantities are the **c2 IMPORT's** status (orchestrator-side detail the builder does
not need), ⛔ **NOT** an S11c-d expected output; S11c-d's engines re-derive the mixing and never compute them. ⛔ Do
**NOT** resurface the c2 F/G re-grounding as a blocker (paused, non-blocking).

⛔ **WITHHELD (the one acceptance criterion — orchestrator-side):** the S11c-d **falsification numeric bound / the
`O(1)` grating reductio** (§0/§5d — magnitude is `R1`-blocked; only the FORM is computable now, `N7`). The directive
and script emit only the computed flux-normalized dimensionless conversion **FORM**; the numeric bound + the reductio
are diffed **on our side**, ⛔ never a builder target.

⚠ **Before launch, leak-gate this directive** ([[feedback_grep_acceptance_dodgeable]] / build-skill §195): `rg` for
co-occurring step symbols in proximity (e.g. a profile/amplitude symbol beside a typed coefficient, `sin` beside
`eta`/`G`, a residual "must vanish"/"must equal"/"=0" phrasing) and READ every hit; probe for SYMBOLS co-occurring,
⛔ never enumerate assembled-expression spellings. A leak gate returning zero hits on a packet that names the step's
own symbols is evidence the **probes** are wrong.

---

## 7 · Run discipline
Detached launch (`setsid` + a completion marker + `Monitor`; the harness reaps `run_in_background`). S11c-d's mixing
is a **derived insertion** on the closed operator (plus a spectral pole-solve) — c2's self-energy `.out` was ~499 MB
and the full cross-engine residual is the ≥64 GB work, so **measure the process that runs**; defer heavy controls
in-band→out-of-band (`DEFERRED_HEAVY_RUNS.md`); ⛔ **never run two memory-heavy CAS jobs concurrently**. Astra needs
`--sandbox danger-full-access` to run SymPy.

⛔⛔ **THE BUILDER'S JOB ENDS AT build → verify-own-deliverable → report.** ⛔ Do **NOT** launch review legs,
comparators, the Wolfram engine, or any downstream step, and ⛔ do **NOT** read the `/build` or `/review-legs` skills
and act on them — **those are the ORCHESTRATOR's, run in a separate process.** A builder that reviews (or spawns
reviewers for) its own output is the self-check trap (`CLAUDE.md` L-R8/G1: whatever writes does not review). The
builder writes only: `scripts/S11c_d_mixing_scattering_sympy_audit.py`, `scripts/S11c_d_exports.py`, its `.out`, and
one `_measurements/S11c_d_sympy_builder_report.md` (the frozen import/reduction map + per-object
which-line-computed-it + literal §5 residuals — ⛔ no verdicts). The orchestrator verifies the deliverable (exists,
non-empty, plausible token count — ⛔ not the exit status) and launches the two independent legs (Codex-written →
**fresh Claude agent + Grok**).
