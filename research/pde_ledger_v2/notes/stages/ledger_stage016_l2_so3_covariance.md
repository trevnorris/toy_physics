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
