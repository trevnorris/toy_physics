# ledger_stage016 — the ℓ=2 SO(3) covariance theorem (Check II-G3a)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 sector). The EARNED-FIRST leg of a 2-way split of `pathA_32`:
this stage carries the **SO(3) covariance-theorem component (1/2) of the joint `ISOTROPY_CALIBRATED`** — the EARNED angular
structure; the grouped-P2 lane isotropy + the calibration partition (which completes the joint and exports the ℓ=2 port
kernel) is **stage 017** (II-G3b).

**Verdict.** `ISOTROPY_CALIBRATED` (JOINT, 2-stage) — landed here as a **PARTIAL** component (016 EARNED, 017 PENDING).
Ledger earned-label: `L2_SO3_COVARIANCE_THEOREM_EARNED`.

**Status.** Exact closed-form / symbolic (float-free): genuine S²-integral Gram, the real `−Δ_S²` Laplace–Beltrami
operator, Rayleigh quotients, eigenfunction residuals, exact `{L,M,T}` dimension vectors, SHA-256 hashes of exact
expression strings; `expect_zero`/`expect_bool` residual asserts, no `scipy`/`numpy`/floats/tolerances. Dual-engine SymPy
**82 PASS** / Mathematica **91 PASS**, exit 0, CWD-independent.

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_32_grouped_p2_isotropy_{sympy.py,.wl}` +
> `reports/pathA_32_grouped_p2_isotropy.md` (the 016 slice = report :3, :9–11, :22/:24/:26, :29–40, :51–56, :60). The report
> is cited for provenance only; the derivation below is self-contained.

---

## 1. What this stage earns

The ℓ=2 angular sector of the frozen throat's response is one 5-dimensional real irreducible representation of SO(3). The
stage derives that fact three ways and lands it as the covariance the grouped-lane isotropy (017) rides.

### 1.1 The real ℓ=2 basis and orthonormality (`Gram = I₅`)
The five real ℓ=2 spherical harmonics
```
Y20  = √(5/16π) (3cos²θ − 1)
Y21c = −√(15/4π)  sinθ cosθ cosφ
Y21s = −√(15/4π)  sinθ cosθ sinφ
Y22c =  √(15/16π) sin²θ cos2φ
Y22s =  √(15/16π) sin²θ sin2φ
```
are orthonormal on the sphere: with the genuine measure `∫_{S²} · dΩ = ∫₀^π ∫₀^{2π} · sinθ dφ dθ`,
```
Gram_ij = ∫_{S²} Y_i Y_j dΩ = δ_ij   ⇒   Gram = I₅.
```
This is a computed integral, not asserted — a non-orthonormal basis (or a wrongly-normalized harmonic) breaks
`gram_is_identity` and drives the `covariant` gate to `FAIL_NOT_COVARIANT`.

### 1.2 The computed `−Δ_S²` eigenvalue `λ_m = 6` (the crux)
The real Laplace–Beltrami operator
```
−Δ_S² f = −[ (1/sinθ) ∂θ(sinθ ∂θ f) + (1/sin²θ) ∂φ² f ]
```
acts on each harmonic. Two independent computations certify the eigenvalue is COMPUTED (not typed):
- the **Rayleigh quotient** `λ_m = ∫ Y_m(−Δ_S²)Y_m dΩ / ∫ Y_m² dΩ = 6` for every m;
- the **eigenfunction residual** `(−Δ_S²)Y_m − λ_m Y_m = 0` for every m (`residuals_zero`) — proving each `Y_m` is a
  genuine eigenfunction, not merely a Rayleigh number.
So `λ_m = ℓ(ℓ+1) = 6` is the SAME across all five m: the ℓ=2 sector is one 5-dim SO(3)-irrep with a single eigenvalue —
the **angular degeneracy** that IS the covariance.

### 1.3 The K₂ angular stiffness uses the computed eigenvalue
The angular part of the ℓ=2 stiffness is
```
K₂ = K̃ + λ_m · T̃_Ω,      M₂ = M̃,
```
with `λ_m` the LIVE computed eigenvalue from §1.2 (assembled as `build_K2(lambdas[name])`, not a typed literal 6). Two
non-vacuous mechanisms make "the K₂ coefficient IS the computed λ" genuinely able-to-fail:
- a **residual-on-the-K₂-coefficient**: extract the `T̃_Ω` coefficient from the assembled K₂ and assert
  `(−Δ_S²)Y − (coeff_in_K2)·Y = 0` — a WRONG coefficient typed into K₂ (e.g. 2) breaks it;
- the **bare `forced_eigenvalue_probe`**: forcing the coefficient to 2 ≠ 6 → `coefficient_equals = False` →
  `FAIL_NOT_COVARIANT` (self-ablation `forced = 6` suppresses). The probe is built on the bare `K̃ + forced·T̃_Ω` (the
  angular self-overlap is the Gram diagonal `= 1` — a 016 fact), NOT via 017's lane assembly.

> ⚠ **De-count note.** The source's `k_coeff_equal = (k_coeff_used − lambdas == 0)` with `k_coeff_used := rayleigh :=
> lambdas` is a vacuous `λ − λ ≡ 0` self-compare (always True, unable to fail). It is **de-counted** here — its role is
> carried by the residual + the bare forced probe above (the fidelity leg confirmed it survives only as a documentation
> string, not a counted check).

### 1.4 Angular dimensional consistency
On pathA_32's own convention — VOLUME densities on the wall measure `dV = a²·dw·dΩ` (dimension `L³`), with `β₂`
dimensionless and `β₂' = L⁻¹`:
```
[μ_η]=M L⁻³, [T_w]=M L⁻¹ T⁻², [K_η]=M L⁻³ T⁻², [T_Ω]=M L⁻³ T⁻²
M₂ = μ_η β₂² dV  ⇒ [M₂] = M
K₂ = (T_w β₂'² + K_η β₂² + λ·T_Ω β₂²) dV  ⇒ each term M T⁻², [K₂] = M T⁻², [K₂/M₂] = T⁻².
```
The dimensions are SOURCED from the density integrals (not back-solved from K₂/M₂). Corrupting the sourced `[T_Ω]`
(and its assembled `T̃_Ω` scalar) by one power of L makes the K₂ term-sum inhomogeneous → `FAIL_DIMENSIONAL`; the same
holds for a one-power corruption of `[μ_η]`, `[T_w]`, or `[K_η]`.

### 1.5 Dimension-object enumeration (step (a), frozen before `.wl` emission)

**Membership rule and working.** The enumeration counts a dimension object once at the clean, stable binding from which
it can be read; walker-local aliases such as `measureDim` are listed by their returned top-level
`baselineDim["Dims"][…]` read path rather than counted a second time.

- Named dimension-vector literal helpers are **out**: `zeroDim`, `expectedM`, `expectedK`, and `expectedRatio` are,
  respectively, the walker's neutral element and independently asserted targets, not declarations of additional physical
  quantities. Each is nevertheless listed below.
- The computed numeric angular coefficient `lambdaRef` is **out**: it has no dimension-vector declaration and acquires
  `zeroDim` only through the walker's numeric fall-through branch. It is nevertheless listed below because it is a live
  factor in the dimension walk.
- Clean symbol→vector rule-table entries are **in**, all twelve of them. This includes the four entries whose values are
  aliases of helper bindings (`dOmegaDim`/`beta2Dim` via `zeroDim`, `Mtilde` via `expectedM`,
  `Ktilde`/`TomegaTilde` via `expectedK`): the association keys are distinct live declarations consumed by the walker.
- Clean computed per-quantity results are **in**, all nine returned under stable top-level `baselineDim["Dims"]`. The
  arity self-check's clean `evalDimensional` re-run returns the same nine-entry `Dims` association into a transient
  `Module` local and is **out** as a duplicate evaluation result; it is listed separately below.
- Expected-target and neutral-element objects are **out** individually, as recorded in the first five rows; this includes
  the unbound measure target `{3,0,0}` as well as the four named helpers.
- The unbound `{1,0,0}` corruption displacement is **out** and listed separately: it is mutation machinery, not a clean
  quantity.
- The four mutation-scoped rule maps and their four deliberately-corrupted evaluation results are **out** individually:
  they exist only to make the able-to-fail probes fail. The `mu_eta_density` result completes with a full nine-entry
  `Dims` map but `Ok=False`; the other three take the caught inhomogeneity path and have empty `Dims` maps.

| `.wl` object and definition locus | artifact status | coverage reason / read locus |
|---|---|---|
| `zeroDim` (`.wl:23`) | not emitted | Private neutral element returned for numeric/zero expressions and empty sums by `dimOf` (`.wl:140,152,159`); its two clean aliases are listed separately as rule-table entries. |
| `expectedM` (`.wl:24`) | not emitted | Independent expected-target helper read by the clean gate (`.wl:270-272`) and by assertions (`.wl:522,527`); the `Mtilde` rule entry that aliases it is a separate listed declaration. |
| `expectedK` (`.wl:25`) | not emitted | Independent expected-target helper read by the clean gate (`.wl:271-272`) and assertions (`.wl:523-528`); the two rule entries that alias it are listed separately. |
| `expectedRatio` (`.wl:26`) | not emitted | Independent expected-target helper read by the clean gate (`.wl:272`) and ratio assertion (`.wl:529`), not a walked quantity. |
| unbound measure expected target `{3,0,0}` (`.wl:270,521`) | not emitted | Independent typed target for the walked `measureDim`/`dims["measure"]`; target divergence is checked against the computed result, but the target is not itself a physical declaration or walker result. |
| unbound corruption displacement `{1,0,0}` (`.wl:308,310`) | not emitted | Mutation-only increment used to construct the four corrupted rule maps (and the coordinated `TomegaTilde` corruption); it must never reach the clean artifact. |
| `lambdaRef = lambdas["20"]` (`.wl:227`), used in `k2Ref` (`.wl:229`) and `kOmegaTermExpr` (`.wl:235,263-264`) | not emitted | Computed numeric angular coefficient, not a symbol→vector declaration. Its dimension is assigned by the **fall-through**, not by a declaration: `dimOf`'s branch `TrueQ[expr == 0] \|\| NumericQ[expr], zeroDim` (`.wl:140`) sees the live value `6` and returns the neutral dimension. |
| `dimRules[aDim]`, literal entry in `makeDimRules[]` (`.wl:238-239`), bound at `.wl:315` | emitted as `dim_rules.a` | Clean symbol declaration consumed by `dimOf` (`.wl:141`) and live at top level through `dimRules` (`.wl:315`). |
| `dimRules[dwDim]`, literal entry in `makeDimRules[]` (`.wl:238-240`), bound at `.wl:315` | emitted as `dim_rules.dw` | Clean differential-length declaration consumed by the measure walk (`.wl:231,259`) and live through `dimRules`. |
| `dimRules[dOmegaDim]`, alias entry to `zeroDim` (`.wl:238-241`), bound at `.wl:315` | emitted as `dim_rules.d_omega` | Distinct clean solid-angle declaration; although its value aliases the neutral helper, the walker reads this key through `dimOf` (`.wl:141`). |
| `dimRules[beta2Dim]`, alias entry to `zeroDim` (`.wl:238-242`), bound at `.wl:315` | emitted as `dim_rules.beta2` | Distinct clean profile declaration consumed in the mass and stiffness integrands (`.wl:232,234-235`) and read by `dimOf`. |
| `dimRules[beta2PrimeDim]`, literal entry in `makeDimRules[]` (`.wl:238-243`), bound at `.wl:315` | emitted as `dim_rules.beta2_prime` | Clean radial-derivative declaration consumed in `kTwTermExpr` (`.wl:233`) and read by `dimOf`. |
| `dimRules[muEtaDensity]`, literal entry in `makeDimRules[]` (`.wl:238-244`), bound at `.wl:315` | emitted as `dim_rules.mu_eta` | Clean pathA_32 volume-density declaration consumed in `m2IntegralExpr` (`.wl:232`) and read by `dimOf`. |
| `dimRules[TwDensity]`, literal entry in `makeDimRules[]` (`.wl:238-245`), bound at `.wl:315` | emitted as `dim_rules.T_w` | Clean pathA_32 wall-tension density declaration consumed in `kTwTermExpr` (`.wl:233`) and read by `dimOf`. |
| `dimRules[KetaDensity]`, literal entry in `makeDimRules[]` (`.wl:238-246`), bound at `.wl:315` | emitted as `dim_rules.K_eta` | Clean pathA_32 wall-stiffness volume-density declaration consumed in `kEtaTermExpr` (`.wl:234`) and read by `dimOf`. |
| `dimRules[TOmegaDensity]`, literal entry in `makeDimRules[]` (`.wl:238-247`), bound at `.wl:315` | emitted as `dim_rules.T_Omega` | Clean angular-stiffness volume-density declaration consumed in `kOmegaTermExpr` (`.wl:235`) and read by `dimOf`. |
| `dimRules[Mtilde]`, alias entry to `expectedM` (`.wl:238-248`), bound at `.wl:315` | emitted as `dim_rules.M_tilde` | Distinct clean reduced-scalar declaration used by `m2Core` (`.wl:228`) and read by `dimOf`; aliasing the target value does not erase the symbol binding. |
| `dimRules[Ktilde]`, alias entry to `expectedK` (`.wl:238-249`), bound at `.wl:315` | emitted as `dim_rules.K_tilde` | Distinct clean radial-stiffness scalar declaration used by `buildK2` (`.wl:216,229`) and read by `dimOf`. |
| `dimRules[TomegaTilde]`, alias entry to `expectedK` (`.wl:238-250`), bound at `.wl:315` | emitted as `dim_rules.T_Omega_tilde` | Distinct clean reduced angular-stiffness declaration used by `buildK2` (`.wl:216,229`) and read by `dimOf`. |
| `baselineDim["Dims"]["measure"]`, computed as `measureDim` (`.wl:259,277-278`), bound at `.wl:316` | emitted as `baseline_dims.measure` | Clean computed result of `dimOf[measureExpr, rules]`; returned from `evalDimensional` and live through top-level `baselineDim`. |
| `baselineDim["Dims"]["M2_integral"]`, computed as `m2IntegralDim` (`.wl:260,277-279`), bound at `.wl:316` | emitted as `baseline_dims.M2_integral` | Clean computed mass-integral result, live through top-level `baselineDim`. |
| `baselineDim["Dims"]["T_w_beta_prime_sq"]`, computed as `kTwDim` (`.wl:261,277-280`), bound at `.wl:316` | emitted as `baseline_dims.T_w_beta_prime_sq` | Clean computed wall-tension term result, live through top-level `baselineDim`. |
| `baselineDim["Dims"]["K_eta_beta_sq"]`, computed as `kEtaDim` (`.wl:262,277-281`), bound at `.wl:316` | emitted as `baseline_dims.K_eta_beta_sq` | Clean computed wall-stiffness term result, live through top-level `baselineDim`. |
| `baselineDim["Dims"]["lambda_T_Omega_beta_sq"]`, computed as `kOmegaDim` (`.wl:263,277-282`), bound at `.wl:316` | emitted as `baseline_dims.lambda_T_Omega_beta_sq` | Clean computed angular-stiffness term result, live through top-level `baselineDim`. |
| `baselineDim["Dims"]["K2_integral"]`, computed as `k2IntegralDim` (`.wl:264,277-283`), bound at `.wl:316` | emitted as `baseline_dims.K2_integral` | Clean computed homogeneous sum result, live through top-level `baselineDim`. |
| `baselineDim["Dims"]["actual_M2"]`, computed as `actualM2Dim` (`.wl:265,277-284`), bound at `.wl:316` | emitted as `baseline_dims.actual_M2` | Clean computed result for the live bare `m2Core=Mtilde` expression (`.wl:228,265`), live through top-level `baselineDim`. |
| `baselineDim["Dims"]["actual_K2"]`, computed as `actualK2Dim` (`.wl:266,277-285`), bound at `.wl:316` | emitted as `baseline_dims.actual_K2` | Clean computed result for `k2Ref=buildK2[lambdaRef]` (`.wl:227-229,266`), live through top-level `baselineDim`. |
| `baselineDim["Dims"]["actual_K2_over_M2"]`, computed as `actualRatioDim` (`.wl:267,277-286`), bound at `.wl:316` | emitted as `baseline_dims.actual_K2_over_M2` | Clean computed ratio result, live through top-level `baselineDim`. |
| `dimProbeResult["Dims"]`, full clean re-run result local to `runAritySelfCheck[]` (`.wl:614-616`) | not emitted | A second clean `evalDimensional[lambdaRef,m2Core,k2Ref,dimRules]` returns the full nine-entry association from `.wl:277-286`, duplicating the already listed stable `baselineDim["Dims"]` quantities. It dies at the `Module` boundary after only the `"Ok"` key is checked (`.wl:624`), so it is not a new artifact quantity. |
| `corruptRulesFor["mu_eta_density", dimRules]` (`.wl:300-313`, call at `.wl:318-320`) | not emitted | Mutation-only copy shifts `muEtaDensity` by one `L` (`.wl:302,308`); it is deliberately corrupted and exists only inside the density-probe construction. |
| `corruptRulesFor["T_w_density", dimRules]` (`.wl:300-313`, call at `.wl:318-320`) | not emitted | Mutation-only copy shifts `TwDensity` by one `L` (`.wl:303,308`); it is deliberately corrupted and exists only inside the density-probe construction. |
| `corruptRulesFor["K_eta_density", dimRules]` (`.wl:300-313`, call at `.wl:318-320`) | not emitted | Mutation-only copy shifts `KetaDensity` by one `L` (`.wl:304,308`); it is deliberately corrupted and exists only inside the density-probe construction. |
| `corruptRulesFor["T_Omega_density", dimRules]` (`.wl:300-313`, call at `.wl:318-320`) | not emitted | Mutation-only copy shifts both `TOmegaDensity` and `TomegaTilde` by one `L` (`.wl:305,308-310`); it is deliberately corrupted and exists only inside the density-probe construction. |
| `densityCorruptions["mu_eta_density"]["Dims"]` (`.wl:318-321`) | not emitted | Mutation evaluation of the prior corrupt map. AB5 measured the exact value as `<\|"measure" -> {3, 0, 0}, "M2_integral" -> {1, 1, 0}, "T_w_beta_prime_sq" -> {0, 1, -2}, "K_eta_beta_sq" -> {0, 1, -2}, "lambda_T_Omega_beta_sq" -> {0, 1, -2}, "K2_integral" -> {0, 1, -2}, "actual_M2" -> {0, 1, 0}, "actual_K2" -> {0, 1, -2}, "actual_K2_over_M2" -> {0, 0, -2}\|>`: the walk completes and returns all nine entries with `Ok=False`, rather than taking the catch. It remains out because it is deliberately corrupted and only its verdict data are consumed (`.wl:538-545`). |
| `densityCorruptions["T_w_density"]["Dims"]` (`.wl:318-321`) | not emitted | Mutation evaluation of the prior corrupt map; the inhomogeneity catch returns an empty `Dims` association (`.wl:296-297`), and only its `Ok`/error verdict data are consumed (`.wl:538-545`). |
| `densityCorruptions["K_eta_density"]["Dims"]` (`.wl:318-321`) | not emitted | Mutation evaluation of the prior corrupt map; the inhomogeneity catch returns an empty `Dims` association (`.wl:296-297`), and only its `Ok`/error verdict data are consumed (`.wl:538-545`). |
| `densityCorruptions["T_Omega_density"]["Dims"]` (`.wl:318-321`) | not emitted | Mutation evaluation of the prior corrupt map; the inhomogeneity catch returns an empty `Dims` association (`.wl:296-297`), with `Ok`/error fields feeding `dimProbe` (`.wl:323-340`). |

### 1.6 Physics-leg verdict (step (c1)) — derived from the model, independently of the conversion

Run on a fresh leg against the stage's existing declarations, deriving each dimension from
`docs/model_map.md` §2 rather than checking a claim. ⚠ Its naming determinations were produced
**without sight of the `.wl` build's proposals**, and vice versa; where the two agree below, that is two
independent parties, not one reviewing the other.

**Measure — everything hangs on this.** Stage016 integrates on `dV = a²·dw·dΩ` (`.wl:231`,
`sympy:371`), the throat wall's own **3-volume** in the 4D bulk. So every `*_density` here is a
per-wall-3-volume density, the same measure class as `ρ_br = M L⁻³`. stage013 uses `4π∫₀^{L0}dw`, a
**line** measure. The two differ by exactly the `a²` Jacobian.

**(1) Per-quantity verdict — 21 of 21 CORRECT on this stage's own convention.** Load-bearing routes:
`[μ_η]=M L⁻³` from `M₂=∫μ_ηβ₂²dV` with `[q₂]=L` ⇒ `[M₂]=M`, cross-checked as `m·ρ·δ = M·L⁻⁴·L`;
`[T_w]=M L⁻¹T⁻²` as energy per wall-3-volume (a 3-brane tension), cross-checked against the model's own
`[μ_R]=[c_γ²][ρ_br]`, with `√(T_w/μ_η)=L T⁻¹` a genuine wall speed; `[K_η]=[T_Ω]=M L⁻³T⁻²`.
⚠ **`[β₂]=1` is a CONVENTION, not a model consequence** — only the product `q₂·β₂` is fixed (`=L`), and
neither engine ever declares `q₂`. ⚠ **`T_Ω`'s independence is UNDETERMINED**: an isotropic wall gives
`T_Ω = T_w/a²`, an *identical* dimension, so no dimensional check here can support or refute the
register's decision to count `T_Ω` as a separate CALIB knob (`parameter_register.md:182`).

**(2) D4 name determinations.** ✅ `a` is the same throat radius as stage018's — a physical-radius
identity, not merely both being `L` (independently reached by both legs). ✅ `K_η`, `μ_η`, `T_w` and
`T_Ω` are **reduction levels** of the 013/016/023 families, ⛔ **not** renamed apart (§7).
⛔ **Sharpest silent-merge hazard, found by this leg and not by the build:** `T̃_Ω` (016) and stage023's
`T_Ω` both carry `M T⁻²` and are the same ℓ(ℓ+1)-shaped radial reduction one ℓ apart. They coincide
**iff** the radial profile is ℓ-independent — which the model does not assert and R42 records as
PENDING ⇒ **different quantity until R42 exists.** This is the `c_s0`/`c_S` shape from stage018: same
dimension, same word-family, attractive wrong merge, and the comparator is blind to it by construction.
The current naming keeps them apart; the hazard is recorded so a later D4 pass cannot merge them.

**(3) Coverage — 12 of 21 records are declared literals in BOTH engines** (`sympy:355-366` ↔
`.wl:239-250`, identical tuples), independently counted by the build and this leg. Of the 9 computed,
**3 are self-referential** (`actual_M2`, `actual_K2`, `actual_K2_over_M2` walk a declaration back to the
constant that defines it) and 6 are genuine algebra over the literals. **0 are computed from any
physical input.** For scale: stage013 5/15, stage018 6/10, the 037 spike 8/21 — this is the largest
non-independent fraction measured so far. ⇒ **A clean cross-engine result here certifies transcription
fidelity only.** ⚠ 7 of the 12 declarations carry no ablation coverage at all.

⭐⭐ **Sharper still, from the adversarial leg: the artifact has 21 records but only 12 FREE VALUES.**
All nine `baseline_dims.*` records are pure functions of the twelve `dim_rules.*` declarations through
five hard-coded expression lines (`.wl:231-236`), so cross-engine agreement on those nine tests only
that two transliterations walk the same hard-coded trees. `dim_rules.M_tilde` and
`baseline_dims.actual_M2` are **literally the same object printed twice**, and the comparator counts
them as 2 of 21 compared. Read the census as **12 free / 9 derived**, not 21 independent.

**(4) ⭐ Structurally uncheckable — the deliverable.** No amount of conversion fixes these:
- **H1** — `M̃`/`K̃`/`T̃_Ω` are asserted against the very constants that define them (`sympy:364` vs
  `:723`). Neither engine anywhere writes `M̃=∫μ_ηβ₂²dV`; it lives only in print strings ⇒ register
  R35's *"dual-engine dim-verified"* (`parameter_register.md:302`) is **overstated** — what is verified
  is that a symbol declared `M` has dimension `M`.
  ⚠ **Scope corrected by measurement (AB4).** An earlier draft of this item claimed that retyping
  `EXPECTED_M`/`EXPECTED_K` "leaves all three passing". **That is false for the Mathematica engine:**
  setting `expectedM = {7,3,-5}` exits **1** at the top-level `baselineDim["Ok"]` preflight
  (`.wl:317`) with `FAIL  baseline dimensional check`, before any assertion runs — because the
  *independently walked* `M2_integral` stays `{0,1,0}` and contradicts the retyped target. So the
  self-reference is real but **narrower** than stated: it makes the three bare-scalar assertions
  vacuous, it does **not** make the declaration globally unfalsifiable. The narrower claim was not
  isolated by this ablation and remains **UNMEASURED**; isolating it needs an edit that moves the
  target and the walked value together.
- **H2** — a **two-parameter family** of declarations passes every check with all four corruption probes
  still firing. Nothing enforces `[β₂']=[β₂]·L⁻¹` and `q₂` is undeclared, so with `[β₂]=L^p`,
  `[β₂']=L^q` the witness `(p,q)=(1,−1)` gives `[μ_η]=M L⁻⁵` fully green and the register's published
  dimensions silently false. Closure: assert `dim(beta2_prime)==dim(beta2)−Dim(1,0,0)` and declare `q₂`.
- **H3** — `K_η` and `T_Ω` are dimensionally identical, so their two corruption probes fire through one
  `Add`-homogeneity failure: two teeth, one fact. Setting `T_Ω := T_w/a²` passes and would remove a
  counted CALIB knob, shifting the Part-VII count.
- **H4** — the dimensional gate is invariant under this stage's own physics: `λ_m` is dimensionless, so
  the angular theorem and the dimensional block never touch each other in either direction.
- **H5** — `participates_in_verdict` restates `mutation_fires` (`sympy:474-477`, `:481`): two counted
  teeth, one fact.
- **H6** — the export path terminates in a hand-copied boolean: everything 018–024 ride on takes its
  dimension from the three H1 literals and its cross-stage validation from stage017's typed
  `CITED_016_DIMENSIONAL_OK = True`. **The three least-checked records are the only ones exported.**
- **H7 — record wiring is structurally uncheckable for equal-valued records.** The 21 records occupy
  only 10 distinct exponent triples. Every unordered pair within each of these five groups can have
  its emission bindings swapped with no observable effect on the transcript, the in-file assertions,
  or a name-keyed comparator:
  `{1,0,0}` → {`dim_rules.a`, `dim_rules.dw`};
  `{0,0,0}` → {`dim_rules.d_omega`, `dim_rules.beta2`};
  `{-3,1,-2}` → {`dim_rules.K_eta`, `dim_rules.T_Omega`};
  `{0,1,0}` → {`dim_rules.M_tilde`, `baseline_dims.M2_integral`, `baseline_dims.actual_M2`}; and
  `{0,1,-2}` → {`dim_rules.K_tilde`, `dim_rules.T_Omega_tilde`,
  `baseline_dims.T_w_beta_prime_sq`, `baseline_dims.K_eta_beta_sq`,
  `baseline_dims.lambda_T_Omega_beta_sq`, `baseline_dims.K2_integral`,
  `baseline_dims.actual_K2`}. The five singleton groups are
  `dim_rules.beta2_prime` → `{-1,0,0}`, `dim_rules.mu_eta` → `{-3,1,0}`,
  `dim_rules.T_w` → `{-1,1,-2}`, `baseline_dims.measure` → `{3,0,0}`, and
  `baseline_dims.actual_K2_over_M2` → `{0,0,-2}`. AB1 directly demonstrated the
  `K_tilde`/`T_Omega_tilde` swap. Correct name→binding wiring is established instead by
  `_scratch/stage016/liveness_probe.out`, which independently moved and restored all 21 exact
  association slots; that is a one-time probe artifact, **not a standing gate** rerun by acceptance.

⚠ Four doc-level claims did not survive this derivation and are held for independent verification
before any doc is edited (tasks #36/#38), since overturning a recorded finding needs harder evidence
than confirming one: the *"`K_η=T_wβ²` does not transfer"* catch, the *"same physical constant, two
dimensions"* framing, the *"computed"* provenance verb at `parameter_register.md:182-184`, and the
`{B̃,Z̃}` row at `:185` whose cited stage has no dimension machinery.

---

## 2. The able-to-fail battery (016-owned)

The verdict runs a SCOPED gate chain (016's gates only — `dimensional_ok`, `covariant`, `tautology_clear`,
`able_to_fail_ok`; 017's stability/denominator/lane-collapse/dynamic gates are NOT computed here). 016's aggregate battery
is its three probes:
| probe | mutation → verdict | self-ablation |
|---|---|---|
| `wrong_eigenvalue` | force K₂ coeff = 2 ≠ computed 6 → `FAIL_NOT_COVARIANT` | `forced = 6` suppresses |
| `tautology_hash_collision` | identical harmonic inputs → non-distinct SHA-256 hashes → `FAIL_TAUTOLOGICAL` | `{20,21c,22c}` distinct suppresses |
| `dimensional_corrupt_T_Omega` | corrupt sourced `[T_Ω]` (+`T̃_Ω`) by one L → `FAIL_DIMENSIONAL` | restore dims suppresses |
Each probe flag is a COMPUTED conjunction (its mutation fires its own verdict ∧ its self-ablation suppresses the fail);
**neutering any one flips the 016 aggregate `able_to_fail_ok` false** (`neutering_any_probe_flips_false`). No probe flag is
a stamped literal.

---

## 3. Honest scope

- **ANGULAR earned / RADIAL calibrated (the `ISOTROPY_CALIBRATED`-not-`PASS` reason).** 016 DERIVES the angular structure
  (Gram, `λ_m=6`, the K₂ angular-stiffness form). The RADIAL profile `β₂(w)` and the radial/support scalars remain FROZEN
  calibration inputs (not derived from the Gate-1 `R0` support equation). That frozen-ness is why the joint is CALIBRATED
  rather than PASS; the calibration PARTITION is stage 017's.
- **Deferred (Gate 4/5/6, sim-deferred).** The 54/5 quadrupole normalization, the outgoing odd-N coefficients, and the
  solved nonlinear branch data (`G = GENUINE_BLOCKED`) are downstream work, not 016's.

---

## 4. Consumed / exported

- **Consumed — PROVENANCE + self-contained dimensional integrity (NOT a cross-stage dual-site relation).** ⚠ pathA_32's
  VOLUME-density / dimensionless-`β₂` convention differs from stage 013's LINE-density / `β=L⁻¹` convention (related by an
  `∫a²dΩ ≈ L²` bridge, not equal), so the stage-013 relation `K_η = T_w β²` does NOT transfer and the density dims do NOT
  equal the 013 register. Therefore: the physical wall constants `μ_η`, `T_w` (counted at 013) + the frozen radial
  reference `R0(w)`/`β₂(w)` + the radial scalars `M̃`/`K̃`/`T̃_Ω` + the **Gate-1 D/N boundary provenance** (stage 012 R28
  `BC_DEPENDENT`, stage 011 R26) are cited as PROVENANCE; the genuine able-to-fail integrity is pathA_32's OWN
  `[M₂]=M`/`[K₂]=M T⁻²` dimensional consistency (§1.4). `c_S` is NOT consumed (the covariance theorem is speed-free).
  `B̃`/`Z̃` support scalars are NOT 016-consumed (they are 017's D-lane).
- **Register.** ZERO new counted knobs (a structural covariance edge, like 011/012/014). New tracked-not-counted symbols
  first-appearing here — `T_Ω`/`T̃_Ω` (the angular-stiffness density) and `β₂(w)` (the frozen ℓ=2 radial profile) — have
  their COUNTING **deferred to 017's calibration partition** (016 uses but does not count them). Structural edge **R34**
  (the ℓ=2 SO(3) covariance provenance; discharges nothing). Part-II counted CALIB set unchanged at
  `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) = 4.
- **Exported.** The ℓ=2 SO(3) covariance theorem (`Gram=I₅`, `λ_m=6`, the K₂ angular-stiffness form `K̃ + λ_m·T̃_Ω`) →
  stage 017 (which rides it for the grouped-lane isotropy, then exports the ℓ=2 port kernel to 018–021 + 022/023 + 024).
  The port-kernel export is 017's, not 016's.

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O. The `.wl` is a **genuinely independent native route**
(native `Integrate`/`FullSimplify` S² integrals, native `D[...]` Laplacian, native `Hash[...,"SHA256"]`, its own
`dimOf`/`scopedVerdict`/`verdictFromGates`) — no `Get`/`Import`/`Export`/YAML, no mirroring of the `.py`; the source pair's
scratch-YAML engine-agreement handoff is severed. Agreement is transcript-level: both engines emit `Gram=I₅`, `λ_m=6`, the
three probe verdicts, and the 016 aggregate. Arity self-check + unevaluated-leakage scan carried.

**Directive review** used the Codex→Grok→Codex bookend: two BLOCKING Codex folds (the aggregate-neuter meta-test wording;
narrowing the consumed-packet claim) + a Grok compute-verify pass that caught the volume-vs-line dimensional-convention
mismatch (the cross-stage dual-site was dimensionally non-transferable → reframed to provenance + self-contained
integrity), plus five hardening nits (the `k_coeff_equal` de-count strengthened to a residual-on-the-assembled-K₂-
coefficient; the forced probe made bare; `B̃`/`Z̃` trimmed; live `lambdas` in the dim re-target); a final Codex confirm →
`DIRECTIVE_CLEAN`.

**Tri-review** on fresh agents: `FIDELITY_CLEAN` (an independent read hand-re-derived the harmonics, the Laplace–Beltrami
operator, all Rayleigh `λ=6` with zero residuals, and the full L/M/T dim walk incl. the T_Ω inhomogeneity — no dropped
check, no transliteration error, the `k_coeff_equal` de-count correct) + `ADVERSARIAL_CLEAN` (per-tooth ablation: all 82
SymPy checks + 4 native `.wl` ablations fired at their own assert; the `k_coeff_equal` de-count, aggregate battery,
sourced-density dims, and mirror-`.wl` all clean) — with ONE low-severity stamped-literal (`participates_in_verdict`)
found and remediated to a computed verdict-propagation check; a fresh-agent `REVERIFY_CLEAN` confirmed the now-computed
tooth fires under its own ablation (residual = 1, exit 1 at its own assert, both engines), no regression.
