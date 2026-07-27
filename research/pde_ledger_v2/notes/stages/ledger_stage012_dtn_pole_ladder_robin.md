# ledger_stage012_dtn_pole_ladder_robin

## Status

**Part II — Gravity. II-G1b (build-order 012).** Reshape of the **DtN pole-ladder + Robin-falsifier leg** of gate
**`pathA_30`** (the gate splits into stages 011 + 012 per the finalized Part-II split). Source top-line, verbatim:
**`DN_UNITTEST_BC_DEPENDENT`** — a JOINT verdict, and **this stage COMPLETES it.** Stage 011 carried the
reduction-certificate component (`REDUCTION_CERTIFIED`: the frozen longitudinal operator `L_s` genuinely produced + its
validity window). This stage takes that cited `L_s`, solves the Dirichlet/Neumann boundary-value problem on the frozen
throat, and derives the resonator's Dirichlet-to-Neumann map, its half-shifted resonance ladder, its round-trip phase,
and the Robin counterfactual that proves the D/N determination is not hardcodable — then lands the verdict at
`BC_DEPENDENT` because the boundary condition itself is imposed (`bc_derivation_emitted=False`, a banked calibration
input). Unlike stage 011 (a bridge-strip **plus a de-rig**), this stage is a **near-pure bridge-strip**: the physics is
already genuinely computed in the source and the `.wl` transfer-matrix route already carries it, so the reshape severs
the scratch-YAML bridge and keeps the native content.

- **CONSUMED (cited, dual-site integrity, NOT re-derived):** stage 011's frozen const-coefficient Helmholtz operator
  `L_s = ψ''(s) + (ω/c_S)²ψ(s)`, the domain `[0, L0]` (cap `R0(L0)=0`), and the speed `c_S` (`k = ω/c_S`,
  `c_S² = 5Kρ*⁴/m` = Part I edge R1 at `ρ*`). The `dsolve` of the cited `L_s` OPENS this stage.
- **POSTULATED (labeled):** the D/N boundary pair — a Dirichlet condition at the brane/mouth `s=0` and a Neumann
  condition at the regular pinch-off cap `s=L0` — is the IMPOSED unit-test boundary condition (`bc_provenance=imposed`),
  a banked CALIBRATION input, NOT derived here (edge R28). `α` (the Robin cap admittance, `[α]=L⁻¹`) is a
  control-construction symbol that builds the falsifiable counterfactual, not the physics.
- **EARNED:** the outward-mouth DtN `Z₀₀ = −(ω/c_S)tan(L0ω/c_S)` (derived via the D/N coefficient-matrix `LUsolve`, then
  compared to the typed target — a genuine derived-vs-typed comparison); the half-shifted pole ladder
  `ω_n = πc_S(n+½)/L0`; the static small-ω series; the round-trip `R_rt = 1`; and the **Robin falsifier** (the six
  `counterfactual_guard` members, each a computed residual).
- **HONEST SCOPE (banked calibration — NOT a de-rig target, NOT a tooth):** `bc_derivation_emitted = False` — the
  explicit mouth/cap `V_wall` gradient derivation is not emitted, so the verdict stays `DN_UNITTEST_BC_DEPENDENT` rather
  than `DN_UNITTEST_PASS`. Earning `PASS` (the `V_wall` derivation) is a genuinely deferred upgrade (edge R28).

Ledger-local earned-label (NOT a source verdict token): `DTN_POLE_LADDER_ROBIN_FALSIFIER_EARNED`. The completed joint
verdict composes as `DN_UNITTEST_BC_DEPENDENT = (011: REDUCTION_CERTIFIED, cited from ledger_stage011) ∧ (012: DtN
ladder EARNED + bc_derivation_emitted=False → BC_DEPENDENT landing, computed here)`.

## Purpose

Stage 011 pinned the interior operator on the frozen throat. This stage poses and solves the boundary-value problem
that turns that operator into a resonator response. The physical question: given the frozen const-coefficient Helmholtz
interior, what is the wall-mouth admittance seen by the brane, where do its resonances sit, and is that admittance a
genuine consequence of the imposed boundary pair (not a hardcoded answer)? The result — the DtN `−(ω/c_S)tan(L0ω/c_S)`,
the half-shifted ladder, and the imposed-BC provenance — is banked forward to stage 013 (harmonic β lift) and stage 017
(calibration input). The source pair computed this through a scratch-YAML bridge between the two engines; the reshape's
burden is to sever that bridge, keep both engines genuinely independent (the `.py`'s coefficient-matrix `LUsolve` route
and the `.wl`'s transfer-matrix route), and keep the Robin falsifier + the honest labels un-stampable.

## The derivation (both engines, own routes)

- **The general solution (opens 012).** `dsolve` of the cited `L_s = ψ'' + (ω/c_S)²ψ = 0` on `[0,L0]` (`k = ω/c_S`)
  gives the fundamental pair `ψ(s) = C₁ sin(ωs/c_S) + C₂ cos(ωs/c_S)`. This is the opening move — stage 011 stopped just
  before it (setting up `L_s`); stage 012 cites `L_s` and calls `dsolve` on it.

- **The D/N boundary-value problem → the DtN (EARNED).** Impose the Dirichlet mouth trace `ψ(0) = ψ_M` and the Neumann
  cap `ψ'(L0) = 0`. The 2×2 coefficient matrix in `(C₁,C₂)` has determinant `dn_det = −ω cos(L0ω/c_S)/c_S`; solving it
  by `LUsolve` gives `bc_applied_solution = ψ_M·(sin(ωs/c_S)tan(L0ω/c_S) + cos(ωs/c_S))`. The outward-mouth
  Dirichlet-to-Neumann admittance is `Z₀₀ = −∂ₛψ|₀ / ψ|₀ = −(ω/c_S)·tan(L0ω/c_S)`. The check
  `dtn_matches_target := expr_equal(dtn − (−k·tan(k·L0)), 0)` compares this DERIVED `dtn` (from the `LUsolve`) to the
  typed target — a genuine derived-vs-typed comparison, not an `X≡X` (the adversarial leg confirmed even hardcoding the
  target is caught, since perturbing the Neumann RHS recomputes the derived DtN and the comparison fails).

- **The half-shifted pole ladder (EARNED).** The DtN denominator's `cos(L0ω/c_S)` factor gives the pole equation
  `cos(L0ω/c_S) = 0`, whose roots are the HALF-SHIFTED ladder `ω_n = πc_S(n+½)/L0`. `halfshift := (pole_residual == 0)`
  is computed by substituting `ω_n` back into the pole equation (residual `0`). The half-shift is the D/N fingerprint:
  a D/D cap would give the integer ladder (see the Robin counterfactual).

- **The static small-ω series (EARNED).** The small-ω expansion `series(dtn, ω, 0, 6) = −L0ω²/c_S² − L0³ω⁴/(3c_S⁴) +
  O(ω⁶)` is the leading static compliance. It is distinct from the actual DC limit `limit(dtn, ω→0⁺) = 0` (the DtN
  vanishes at zero frequency, leading order `∝ ω²`). No separate static solve — the series comes from the dynamic DtN.

- **The round trip (EARNED).** With Dirichlet reflection `r_D = −1` and Neumann reflection `r_N = +1`, the round-trip
  factor is `−exp(2i L0ω/c_S)`. Substituting the D/N ladder gives `R_rt = 1` (`φ0 = 0 mod 2π`) — `round_trip_closes`
  is computed by substitution, confirming the ladder is exactly the resonance condition.

- **The Robin counterfactual (EARNED; the FALSIFIER).** Without a counterfactual, the D/N determination would be
  hardcodable. The Robin cap `∂ₛψ − α ψ = 0` is a one-parameter family interpolating Neumann (α→0) and Dirichlet
  (α→∞); its DtN reads off the `robin_denominator_core = k cos(kL0) − α sin(kL0)`. The `counterfactual_guard` has
  EXACTLY six members, each a computed symbolic residual/comparison (NOT a stamped truthiness):
  - `robin_determinant_emitted` — tied to the computed Robin coefficient-matrix determinant `robin_det`
    (`≡ −robin_denominator_core`), via `nonzero_q(robin_det)` — NOT `bool(hstr(...))` (the source's stampable form was
    hardened here);
  - `recovers_DN_at_alpha0 := expr_equal(dn_from_robin − dtn, 0)` (α→0 recovers the D/N DtN `−(ω/c_S)tan`);
  - `recovers_DD_at_alpha_inf := expr_equal(dd_from_robin − dd_target, 0)` (α→∞ recovers the D/D DtN `k·cot(kL0)`);
  - `halfshift_destroyed_for_DD` — the D/D pole denominator is `sin(kL0)`, whose zeros are the INTEGER ladder
    `πc_S·j/L0` (NOT the half-shifted `(j+½)`), computed on the integer-vs-half-shifted samples;
  - `numeric_alpha_distinct` — a concrete `α = 2/L0` gives a Robin DtN distinct from BOTH the D/N and D/D limits
    (exact rationals, no floats);
  - `dtn_mismatch := (robin_dtn − dtn ≢ 0)` for generic `α`.
  `dd_zero_mode_removable` (`= 1/L0`) is a computed artifact, held OUT of the guard (asserted a non-member). The guard
  feeds `all(counterfactual_guard.values())` → else `FAIL_COUNTERFACTUAL`; breaking `robin_denominator_core` flips the
  recomputed `recovers_DD_at_alpha_inf` and is caught.

- **The 012 dimensional legs + corrupt-`[K]` probe (EARNED).** With `[energy]=(2,1,-2)` (order `L,M,T`), four-volume
  `(4,0,0)`, `[P]=(-2,1,-2)`, `[ρ]=(-4,0,0)`, the chain `[K]=[P]−5[ρ]=(18,1,-2)` gives `[c_S²]=(2,0,-2)`,
  `[c_S]=(1,0,-1)`, `[k]=[ω/c_S]=(-1,0,0)`, `[α]=(-1,0,0)` (so `[α c_S]=(0,0,-1)=[ω]`). Then `[tan_argument]=[k L0]=1`
  (dimensionless) and `[Z00_prefactor]=[−k]=L⁻¹`. The corrupt-`[K]` probe adds one power of `L` to `[K]`, which
  propagates `[c_S²]→[c_S]→[k]`, so `[tan_argument]` is no longer dimensionless AND `[Z00_prefactor] ≠ L⁻¹` →
  `DN_UNITTEST_FAIL_DIMENSIONAL` (`mutation_fires=True`), with a two-verdict self-ablation (with mutation →
  `FAIL_DIMENSIONAL`, without → the clean 012 verdict, `fail_suppressed=True`). This is the 012 half of the shared dim
  block; the `cs_squared` leg was stage 011 (excluded here). Note the propagation model differs from the source (which
  sourced `[c_S]` independently and fired via the `cs_squared` leg): the reshape reroutes through `5Kρ*⁴/m` to keep the
  probe alive in 012 without re-deriving the 011 leg (directive §2h; clean-walk values match the report).

The Mathematica source has exactly **27 dimension-valued objects**, and emits **19 of 27**. `dimensionAxes`, the
`baseRules`/`dimRules` maps, and the returned `buildDimensionalBlock[]` association are containers rather than
additional dimension values; aliases are counted once. The clean and corrupt `walk` results are counted separately
because `.wl:281` and `.wl:287` evaluate the walker with different `K` dimensions. This is the source-coverage record
required by `DIMENSION_REWRITE.md` §4-a.

| `.wl` object and definition locus | artifact status | coverage reason / read locus |
|---|---|---|
| `zeroDim` (`.wl:36`) | not emitted | Private neutral element returned by `dimOf` for constants, empty sums, and dimensionless trig functions (`.wl:220,232,239`); also read by `expectedTanDim` and the `Z00` walk guard (`.wl:261,270`). |
| `LengthDim` / `lengthDim` (`.wl:252`, returned at `.wl:303`) | not emitted | Primitive carrier used to define four-volume, density, `alpha`, and the expected `Z00` dimension, and bound to `L0` in each walk (`.wl:254,256,260,262,264`); no standalone artifact name is needed for those derived checks. |
| `EnergyDim` / `energyDim` (`.wl:253`, returned at `.wl:304`) | emitted as `energy_dim` | Read directly at the `DIM` print site (`.wl:573`). |
| `FourVolumeDim` / `fourVolumeDim` (`.wl:254`, returned at `.wl:305`) | emitted as `four_volume_dim` | Read directly at the `DIM` print site (`.wl:574`). |
| `PressureDim` / `pressureDim` (`.wl:255`, returned at `.wl:306`) | emitted as `pressure_dim` | Read directly at the `DIM` print site (`.wl:575`). |
| `RhoDim` / `rhoDim` (`.wl:256`, returned at `.wl:307`) | emitted as `rho_dim` | Read directly at the `DIM` print site (`.wl:576`). |
| `KDim` / `kDim` (`.wl:257`, returned at `.wl:308`) | emitted as `K_dim` | Read directly at the `DIM` print site (`.wl:577`) and supplied to the clean walk (`.wl:281`). |
| `OmegaDim` / `omegaDim` (`.wl:258`, returned at `.wl:309`) | emitted as `omega_dim` | Read directly at the `DIM` print site (`.wl:578`) and bound to `omega` in each walk (`.wl:264`). |
| `MassDim` / `massDim` (`.wl:259`, returned at `.wl:310`) | emitted as `mass_dim` | Read directly at the `DIM` print site (`.wl:579`) and bound to `m` in each walk (`.wl:264`). |
| `AlphaDim` / `alphaDim` (`.wl:260`, returned at `.wl:311`) | emitted as `alpha_dim` | Read directly at the `DIM` print site (`.wl:580`) and bound to `alpha` in each walk (`.wl:264`). |
| `ExpectedTanDim` / `expectedTanDim` (`.wl:261`, returned at `.wl:312`) | not emitted | Expected-target helper, read by the clean/corrupt dimensional predicates (`.wl:283,289,294`) and their direct assertions (`.wl:596,603,725`). |
| `ExpectedZ00Dim` / `expectedZ00Dim` (`.wl:262`, returned at `.wl:313`) | not emitted | Expected-target helper, read by the clean/corrupt dimensional predicates (`.wl:284-285,290-291,295`) and their direct assertions (`.wl:597-598,604,726`). |
| `cleanWalk["CsSquaredDim"]` (`.wl:265,272`, instantiated at `.wl:281`) | emitted as `clean_walk.cs_squared_dim` | Computed by `dimOf[5 K rhoStar^4/m, baseRules]` and printed directly (`.wl:581,594`). |
| `cleanWalk["CsDim"]` (`.wl:266,273`, instantiated at `.wl:281`) | emitted as `clean_walk.cs_dim` | Derived from the live `CsSquaredDim` and printed directly (`.wl:582,594`). |
| `cleanWalk["KDimFromOmegaOverCs"]` (`.wl:274`, instantiated at `.wl:281`) | emitted as `clean_walk.k_dim` | Computed by `dimOf[omega/cS, dimRules]` and printed directly (`.wl:583,594`). |
| `cleanWalk["TanArgumentDim"]` (`.wl:268,275`, instantiated at `.wl:281`) | emitted as `clean_walk.tan_argument_dim` | Computed by `dimOf[k L0, dimRules]`, printed (`.wl:584,595`), and asserted (`.wl:283,596`). |
| `cleanWalk["Z00PrefactorDim"]` (`.wl:269,276`, instantiated at `.wl:281`) | emitted as `clean_walk.z00_prefactor_dim` | Computed by `dimOf[-k, dimRules]`, printed (`.wl:585,595`), and asserted (`.wl:284,597`). |
| `cleanWalk["Z00Dim"]` (`.wl:270,277`, instantiated at `.wl:281`) | emitted as `clean_walk.z00_dim` | Computed by `dimOf[-k Tan[k L0], dimRules]`, printed (`.wl:586,595`), and asserted (`.wl:285,598`). |
| `cleanWalk["AlphaCSDim"]` (`.wl:278`, instantiated at `.wl:281`) | not emitted | The walker computes and stores `[alpha cS]`, but the returned field has no read locus: no assertion, prose print, `DIM` print, or sidecar references it. This is a live-but-dead dimension value. |
| `CorruptKDim` / `corruptKDim` (`.wl:286`, returned at `.wl:316`) | emitted as `corrupt_K_dim` | Supplied to the corrupt walk (`.wl:287`) and printed directly (`.wl:587,600`). |
| `corruptWalk["CsSquaredDim"]` (`.wl:265,272`, instantiated at `.wl:287`) | not emitted | The live walker result is read internally to derive corrupt `CsDim` (`.wl:266`), but the returned field has no assertion, print, or sidecar read. This is a live-but-dead dimension value. |
| `corruptWalk["CsDim"]` (`.wl:266,273`, instantiated at `.wl:287`) | emitted as `corrupt_walk.cs_dim` | Derived from the live corrupt `CsSquaredDim` and printed directly (`.wl:588,601`). |
| `corruptWalk["KDimFromOmegaOverCs"]` (`.wl:274`, instantiated at `.wl:287`) | emitted as `corrupt_walk.k_dim` | Computed by `dimOf[omega/cS, dimRules]` and printed directly (`.wl:589,601`). |
| `corruptWalk["TanArgumentDim"]` (`.wl:268,275`, instantiated at `.wl:287`) | emitted as `corrupt_walk.tan_argument_dim` | Computed by `dimOf[k L0, dimRules]`, printed (`.wl:590,602`), and asserted by the mutation gate and tooth (`.wl:289,294,603,725`). |
| `corruptWalk["Z00PrefactorDim"]` (`.wl:269,276`, instantiated at `.wl:287`) | emitted as `corrupt_walk.z00_prefactor_dim` | Computed by `dimOf[-k, dimRules]`, printed (`.wl:591,602`), and asserted by the mutation gate and tooth (`.wl:290,295,604,726`). |
| `corruptWalk["Z00Dim"]` (`.wl:270,277`, instantiated at `.wl:287`) | not emitted | Probe-only derived result, read by `corruptDimensionalOk` (`.wl:291`), which feeds the mutated verdict and asserted self-ablation (`.wl:297,606,727`). |
| `corruptWalk["AlphaCSDim"]` (`.wl:278`, instantiated at `.wl:287`) | not emitted | The walker computes and stores `[alpha cS]`, but the returned field has no read locus: no assertion, prose print, `DIM` print, or sidecar references it. This is a live-but-dead dimension value. |

- **The 012-scoped verdict.** Computed from the 012 rungs — `DN_UNITTEST_FAIL_DIMENSIONAL` → `FAIL_POLE_LADDER` (if not
  `dtn_matches_target and halfshift`) → `FAIL_COUNTERFACTUAL` (if not `all(counterfactual_guard.values())`) → else
  `bc_derivation_emitted` False → `DN_UNITTEST_BC_DEPENDENT` (the landing). The `DN_UNITTEST_PASS` branch exists but is
  never taken (the imposed BC's `V_wall` derivation is deferred). `round_trip_closes` is a standalone EARNED check + a
  tooth, not a separate verdict rung. The joint `DN_UNITTEST_BC_DEPENDENT` is printed as the composition citing stage
  011's `REDUCTION_CERTIFIED`, NOT typed as 012-earned.

## Consumed inputs

**Cited — no file reads; genuine DUAL-SITE citation-integrity (two independently-constructed sites, `siteA − siteB ≡
0`, plus explicit frozen-export anchors):**
- **From `ledger_stage011` (II-G1a):** the frozen const-coefficient Helmholtz operator `L_s = ψ'' + (ω/c_S)²ψ` on
  `[0,L0]`. Site A = the re-typed export form; **site B = reconstructed from the `dsolve`'d fundamental pair
  `{sin(ωs/c_S), cos(ωs/c_S)}`** — positing a generic monic 2nd-order operator `y'' + a y' + b y` and solving the 2×2
  system that requires both solutions be annihilated, recovering `a=0`, `b=(ω/c_S)²` (the operator recovered from its
  null space). Assert `L_s^A − L_s^B ≡ 0` + the frozen-export anchor `consumed_L_s − (ψ''+(ω/c_S)²ψ) ≡ 0`. A corruption
  of site A's export sign OR of site B's fundamental pair (`sin/cos → sinh/cosh`, which recovers `b=−(ω/c_S)²`) fails
  BOTH engines — the null-space reconstruction is a genuinely independent construction, not a `k²`-rename.
- **From `ledger_stage005` (I-2, edge R1), via stage 011:** the sound-speed law `c_S² = 5Kρ⁴/m_GNLS` (EOS exponent-5).
  Site A the literal `5Kρ⁴/m`, site B the EOS route `∂_ρ(Kρ⁵)/m`; evaluated at `ρ*`; every one-site exponent
  corruption (`5→4`, `5→6`) fails BOTH engines, a coordinated both-site drift caught by the frozen-export anchor.
  NOT re-derived. "R1 at ρ*" caveat; EOS exponent-5 stays `IMPOSED`.
- **The domain `[0, L0]`** (cap `R0(L0)=0`) is cited from stage 011; `L0` is already registered (ACTION-geometry) — not
  re-counted or re-solved here.

## Exports

- The outward-mouth DtN `Z₀₀ = −(ω/c_S)tan(L0ω/c_S)`, the resonance ladder `ω_n = πc_S(n+½)/L0`, and the `BC_DEPENDENT`
  provenance (the imposed-BC calibration slot) → stage 013 (harmonic β lift) + stage 017 (calibration input). (Distinct
  from stage009/010's *bulk* slab Helmholtz mode → 024/026 — this is the **frozen-throat** longitudinal resonator.)
- Register: **zero new counted knobs**. `α` (Robin cap admittance, `[α]=L⁻¹`) is a control-construction symbol (tracked,
  not counted, like stage010's `k_warp`); it builds the falsifiable counterfactual, not the physics. `L0` already
  registered at stage 011 (ACTION-geometry). One **structural/imposed** edge: R28 (the D/N boundary is IMPOSED → the
  `BC_DEPENDENT` landing; a `PENDING`/`IMPOSED` obligation, NOT a reduction — discharges nothing; its deferred discharge
  is the mouth/cap `V_wall` derivation that would earn `DN_UNITTEST_PASS`, analogous to R23's obligation).

## Verification

- **Reshape (blueprint §5) — bridge-severing:** stripped the `.py` scratch-YAML/`_sympy_exprs.wl` export, the
  Mathematica-YAML re-read, the `digest_mapping`/yaml helpers, the `mma_exports`/`expression_digest`, the
  report/feed writers, and the `RESULTS_YAML` write inside `main`; and the `.wl`'s `Get[sympyExprFile]`, the `sympy*`
  halves of the agreement equalities, and the `Export`. The `.wl`'s native-vs-native self-consistency checks
  (`robinAlpha0==dtnTransfer`, `robinAlphaInf==ddTransfer`) were PRESERVED. Both engines standalone, print-only, **zero
  file I/O**, float-free (the numeric-`α` check uses exact rationals), ledger idioms. **Clean CONSUME-from-011 cut:** no
  frozen-reduction operator assembly, no reduction certificate, no `operator_is_helmholtz`/`speed_is_cs`/`domain_is_L0`
  de-rig, no `cs_squared` standalone dim leg, no pinch-off domain solve (all stage 011) — 012 cites `L_s`/`[0,L0]`/`c_S`
  and opens at the `dsolve`.
- **The `.wl` independent route:** the KEPT transfer-matrix engine derives the DtN by propagating the Dirichlet mouth
  state `{ψ_M, ψ'_M}` through `transferMatrix[L0]` and imposing the Neumann cap — a genuinely DIFFERENT decomposition
  from the `.py`'s coefficient-matrix `LUsolve`; both derive `Z₀₀ = −(ω/c_S)tan(L0ω/c_S)` and agree (the agreement is
  COMPUTED — corrupting the transfer matrix breaks its own DtN assert). It derives its own pole ladder / static series /
  round trip / Robin counterfactual (native `robin*`/`ddTransfer` route) and its own `tan_argument`/`Z00` dim legs +
  corrupt-`[K]` probe; native `Solve`/`Series`/`Coefficient`; arity self-check present; `ConditionalExpression`
  non-pole wrappers normalized to their value branch (no load-bearing message silenced).
- **Dual-engine:** SymPy **84 PASS / 0 FAIL** · Mathematica **90 PASS / 0 FAIL** (the +6 = the `.wl`'s arity self-check
  block net of two SymPy-only DtN-detail checks), both exit 0, CWD-independent; runner transcripts under
  `scripts/output/` + `mathematica/output/`.
- **Tri-review (fresh agents):** arbiter re-run via the runners (both engines, repo root + foreign CWD, reproduced);
  **`FIDELITY_CLEAN`** (full 012-slice coverage diff, every value hand-re-derived; the 6-member Robin guard confirmed
  each a genuine residual; `dtn_matches_target` genuine derived-vs-typed; `bc_derivation_emitted=False` an honest scope
  flag; the `L_s` null-space dual-site genuinely independent; the consume-from-011 boundary clean); **`ADVERSARIAL_CLEAN`**
  (16-mutant matrix — the Robin guard un-stampable (SM12 flips `robin_determinant_emitted`, a fiat-flip does nothing),
  the `L_s` sinh/cosh corruption fires in both engines (SM2/MM3), the LUsolve-is-decorative attack fails (SM7: even a
  hardcoded target is caught by the Neumann-RHS tooth), `bc_derivation_emitted` load-bearing (SM5/SM5b reveal the
  deferred `DN_UNITTEST_PASS`), corrupt-`[K]` propagates (SM8), all four verdict rungs reachable, the `.wl` transfer
  route independent with computed agreement (MM1), arity clean (MM2, no silent skip), tallies reconcile 84/90, zero
  file I/O).
- **Documented nit (non-blocking, both legs CLEAN):** the `.wl`'s `robin_determinant_emitted` witnesses the typed
  `robinDenominatorCore` (`nonzeroQ`), not a derived coefficient-matrix determinant — the transfer-matrix route has no
  coefficient matrix, so this is directive-§1.5-A-permitted (the `robin_denominator_core ≢ 0` form) and sound; the
  derived Robin DtN is independently exercised by `recovers_DN`/`recovers_DD` and the native `alpha0Agree`/`ddAgree`
  self-checks, and the adversarial MM1 confirms the derived route is caught if corrupted. It is slightly asymmetric vs
  the `.py` (which additionally ties `robin_det` to the core), recorded here for transparency; no remediation cycle was
  triggered (no rig, no vacuous tooth, no escape).
- **Teeth (8, all fire):** (1) Robin-swap — BREAK `robin_denominator_core` → `recovers_DD_at_alpha_inf` (recomputed
  residual) flips → `FAIL_COUNTERFACTUAL` (a fiat forced-true is explicitly forbidden as a no-op); (2) corrupt the
  derived `dtn` → `dtn_matches_target` fails → `FAIL_POLE_LADDER`; (3) integer ladder (drop half-shift) → `halfshift`
  fails; (4) corrupt `r_D`/`r_N` → `round_trip_closes` fails; (5) break the α→∞ path → `recovers_DD_at_alpha_inf` fails;
  (6) numeric-α degeneracy (α→0) → `numeric_alpha_distinct` fails; (7) dual-site consumed-input corruption (`c_S²`
  exponent, `L_s` export sign, `L_s` site-B sinh/cosh) → integrity assert in BOTH engines; (8) corrupt-`[K]` → dim legs
  off → `FAIL_DIMENSIONAL` + self-ablation.

## Provenance

- Source gate: `software/stage1_solver/tools/pathA_30_dn_unit_test_{sympy.py,.wl}` (012 slice; sources unchanged);
  `software/stage1_solver/reports/pathA_30_dn_unit_test.md` (`## DtN Derivation` :26–33, `## Pole Ladder` :35–39,
  `## Static Limit` :41–45, `## Round Trip` :47–49, `## Robin Counterfactual` :51–58, `## BC Provenance` :60–64,
  `## Dimensional Check` `tan_argument`/`Z00_prefactor` at :79).
- Reshape directive + tri-review artifacts: `research/pde_ledger_v2/_scratch/ledger_stage012_*` + `_scratch/adv_stage012/`.
  Running-start source map: `research/pde_ledger_v2/notes/stage012_pathA30_dtn_ladder_source_map.md`.
- Split row: `research/pde_ledger_v2/notes/part2_gravity_atomic_split.md` (id 012). Completes the pathA_30 fold opened
  by stage 011 (the frozen-reduction certificate); together they carry the joint `DN_UNITTEST_BC_DEPENDENT`.
