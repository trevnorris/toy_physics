**CLEAR**

Packet identifier (as supplied, not independently recomputed):  
`S11 tail/Abel and inverse stability T1–T4 fidelity contract v1`  
`aggregate_sha256 = 79a8697aebca1ba1c1ac29e694083daed6c1df512259d07fc59444b27e2b8a48`

This is a statement-fidelity review of the bounded T1–T4 contract. It is not a numerical S11c error certificate, scattering-convergence proof, or clearance of unfinished S11c work. No required fidelity corrections were found. Stronger Fourier/Sobolev rates, outgoing-domain construction, and uniform folded moments remain application obligations, not blockers.

---

## Assumptions and conventions checked

- Amplitudes are `ℝ → ℂ`. Error uses complex norms, not real-part projection.
- Regulator is `exp(-a |x|)` about the fixed origin `x = 0`, with `a ≥ 0`.
- Kept domain `s` is an arbitrary measurable set; `{|x| ≤ R}` is the physical cutoff, including empty `R < 0`.
- Measures are arbitrary on `ℝ`; Lebesgue is the physical case; Dirac is a legitimate exact witness.
- Multiplier `b` is a.e. strongly measurable and pointwise bounded by `B`.
- Integrability of `g` and of `|x| ‖g x‖` are explicit premises of the top-level Abel theorem.
- Inverse theorem: supplied continuous linear equivalence `A : X ≃L[𝕜] Y` and bounded `E`, with `‖A⁻¹‖ ≤ κ`, `‖E‖ ≤ ε`, `κ ε < 1`. Completeness of `X` is required for construction. No finite-dimensional, self-adjoint, real-frequency, or unweighted `L²` hypothesis is present.
- Native identification uses selected `EdgeReduction.__init__` statements only: Gaussian Fourier mass, Abel half-lines, `s = L_W (k − k′)`, measure `L_W / (2π)`, width `h = a / L_W`, development input `L_W = 10`, origin `ξ = 0`.
- Raw `delta_mass` Meijer-G is recorded, not certified. The `2π` mass is the separate Gaussian computation.
- Assessment Fourier/Sobolev remarks are context, not theorems of this increment.

---

## T1 — actual integrals

**Object.** `tailMass` is `∫_{sᶜ} ‖g‖` and `firstMoment` is `∫ |x| ‖g x‖` (`Tail.lean`). Both use `‖·‖` on `ℂ`.

**Integrability split is correct.** `omitted_integral_bound` is the primitive comparison and does **not** certify integrability of `b g`. `truncation_error_bound` does: measurable `s`, integrable `g`, a.e. strongly measurable bounded `b`, then `bounded_product_integrable` and `setIntegral_compl`. That avoids treating Mathlib’s default `0` on a non-integrable Bochner integral as a physical tail.

**Combined Abel/truncation.** `abelFactor a x = Real.exp (-a * |x|)`. `abel_truncation_bound` concludes

```text
‖∫ b g − ∫_s (exp(-a |x|) b g)‖ ≤ B (tailMass + a * firstMoment)
```

under `MeasurableSet s`, `a ≥ 0`, `0 ≤ B`, integrable `g`, integrable first moment, measurable bounded `b`. The proof is the triangle inequality of full Abel error (`abel_pairing_bound`) plus the regulated tail (`truncation_error_bound` on `c = abelFactor • b`, using `|abelFactor| ≤ 1`). Premises for both pieces are present; `hB0 : 0 ≤ B` is redundant given a nonempty domain and `‖b x‖ ≤ B`, not a hidden restriction.

**Convention checks.**
- Negative `R`: `measurable_cutoff` is `{|x| ≤ R}` as a closed set for every real `R`; empty kept set is allowed.
- Zero regulator: `abel_zero` gives factor `1`; bound collapses to a pure tail.
- Phase factor: `real_phase_norm` shows `‖exp(-s x I)‖ = 1`; unused by the main bound, consistent with a modulus majorant.
- Physical specialization to Lebesgue is documentation, not a silent extra hypothesis.

No T1 blocker. Optional sharper form (first moment only on `s`) is not the stated contract.

---

## T2 — actual inverse existence

**Construction.** `perturbEquiv` uses `Units.oneSub (-A.symm.comp E)` under `‖A⁻¹ E‖ < 1`, then `.trans A`. Forward map is `A ∘ (I + A⁻¹ E) = A + E`, identified by `perturbEquiv_eq`. `CompleteSpace X` is on the construction only; `Y` is not assumed complete; no `FiniteDimensional`, adjoint, or `L²` structure.

**Estimates vs existence.** Intermediate lemmas (`inverse_norm_bound`, `inverse_difference_bound`, solution/observation bounds) assume a supplied equivalence `B` with `B = A + E`. `exists_controlled_inverse` constructs that `B` from `relative_error_small` and joins

- `B = A + E`,
- `‖B⁻¹‖ ≤ κ / (1 − κ ε)` with `0 < 1 − κ ε` from `κ ε < 1`,
- `‖B⁻¹ − A⁻¹‖ ≤ κ² ε / (1 − κ ε)`.

**Identities.**
- Resolvent: `B⁻¹ − A⁻¹ = − B⁻¹ ∘ E ∘ A⁻¹`.
- Source: `B⁻¹(f + df) − A⁻¹ f = B⁻¹(df − E A⁻¹ f)`.
- Observation: `‖C(diff)‖ ≤ ‖C‖` times the solution bound, for a **fixed** bounded `C : X →L Z`.

These match the intended Neumann/resolvent calculus. A graph/outgoing realization of the differential operator is an application premise, stated in `Stability.lean` and the fidelity record, not imported as a discharged fact.

No T2 blocker. Packaging solution/observation into `exists_controlled_inverse` would be optional packaging, not a missing existence proof.

---

## T3 — compact native correspondence

Selected constructor statements in `scripts/S11c_d_mixing_scattering_sympy_audit.py` (`EdgeReduction.__init__`, Fourier block ~1618–1623, Abel block ~1624–1647) are:

- Gaussian: `exp(-a u² + i q u)`, mass `∫ dq ∫ du = fourier_mass`, checked as `2π`.
- Abel: `exp(-i s ξ)` with `exp(∓ a ξ)` on the two half-lines.
- `negative = 1/(a − i s)`, `positive = 1/(a + i s)`.
- `constant = 2a/(a² + s²)`, even `a/(a² + s²)`, odd `−i s/(a² + s²)`.
- `prescribe` uses `left * constant + jump * positive` at fixed origin `ξ = 0`.
- `hat` inserts the `y₂ = L_W ξ` Jacobian and divides by `fourier_mass`.

Hand algebra: with `s = L_W q`, `h = a/L_W`, measure `L_W/(2π)`, the constant term is the Poisson kernel `h / (π (h² + q²))` and the step/positive term is `(P_h − i Q_h)/2`. Development input has `"L_W": "10"`. Lean `approved_width` is `a/10` by `rfl`.

The source-check instrument executes those AST statements on a `SimpleNamespace`; it does not import the production driver or build the closed operator. Recorded native operands match the half-line formulas. `delta_mass` is the unevaluated Meijer-G; it is **not** the `2π` check. Four nonzero controls (width `a L_W`, reversed PV sign, missing `1/2`, discarded constant) are tied to these operands.

Docs partition local differential / ordinary integrable / distributional / solve-channel contributions as **open application obligations**. They do not claim an operator-wide estimate or class-wide rate. Assessment items such as `‖P_h * − I‖_{L²→L²} = 1` and `h^α` Sobolev bounds are explicitly deferred.

No T3 blocker. This is tested selected-source translation, not a certified kernel or trial/test pairing.

---

## T4 — controls

The instrument defines 12 paired statement mutations and 16 positives (`EXAMPLES` + `EXTRA_POSITIVES`; canonical-source `MUTATIONS` is empty, as documented).

Recorded contract report: each of the 12 mutants has exit status 1, exactly one diagnostic, declaration `contract_control`, and `⊢ False`. No recorded timeout or import/instrument failure. The 16 positives are recorded PASS.

| Control | Mathematical content |
|---|---|
| Tail / complex tail | Dirac mass at 2 outside `[-1,1]`; error `1` for `1` and for `i` |
| Abel moment | Dirac at 1: `1 − e⁻¹ > 0` |
| Width / sign / `\|x\|` | `1/10` not `10`; `a = −1` grows; `x = −1` still contracts |
| Inverse family | `A = 1`, `E = −1/2`: inverse `2`, difference `1`, source `3/2`, observation `9/2`; omitted operator error `0` is false |
| Critical margin | `E = −1` has no real inverse |
| Extra positives | atomic integrability and moment; `a = 0`; zero amplitude; `1 · ½ < 1` |

Dirac integrals are actual Bochner integrals, not scalar placeholders. Inverse/solution/observation controls are exact scalar arithmetic, not edits of `exists_controlled_inverse`. That limitation is stated and is adequate for the load-bearing **formula factors** (omitted tail/moment, Abel scale/sign, lost `1/(1−κε)`, omitted `E` or `‖C‖`, singular margin). It is **not** a proof that every wrong edit of the general theorems would be detected, and the docs do not claim that.

Generic Dirac/scalar witnesses are not Lebesgue profile moments or an outgoing operator margin. The fidelity record says so.

No T4 blocker.

---

## Limits of this review

Read: `MANIFEST.json`, policy, T1–T4 coverage/fidelity/verification, assessment, four modules, audit root, contract/source/validation reports, source-check instrument, native constructor/prescribe/hat fragments, development input, lakefile/toolchain/manifest pins.

Did **not** run Lean, did **not** run the native or contract instruments, did **not** recompute file or aggregate hashes, did **not** inspect Mathlib lemma source (`Units.oneSub`, `setIntegral_compl`, `setIntegral_dirac`, Bochner comparison lemmas). Builds, 41 standard-axiom rows, 12 rejections, 16 positives, and object hashes are treated as author evidence. Mathlib is the recorded pin `db584cd6d46c92f209a44c0f1c829460d327499d` plus direct import hashes, not a rebuilt library.

Unfinished S11c items (variable coefficients, nonlinear-pole projectors, limiting absorption, full scattering certificate, D5, systematic CAS bridge) were not used as expansion or as blockers.

**Disposition:** T1–T4 statement fidelity is **CLEAR** for this bounded conditional contract.
