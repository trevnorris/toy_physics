# S1–S4 fidelity review

**Verdict: CLEAR** for the bounded S1–S4 contract in this isolated packet.

This is an independent reading of the supplied Lean, native AST, instruments, and recorded reports. It does not rest on the author’s PASS label. Compilation, mutation, and native execution were not re-run.

## What the packet actually proves

The finite object is the selected `construct` fragment in `_measurements/S11c_d_finite_scattering.py` (lines 371–393): row/column balancing, `lstsq`, unscaling, `A x - b` and the scaled residual, modal observation `P V⁻¹(T c - t_in)`, same-end current assembly, and construction of `incoming_flux`. Lean treats invertibility, inverse-norm bounds, observation maps, current entries, and denominator margins as premises. That matches the native split: `lstsq` rank, singular values, and `balancedCondition` stay outside the kernel.

### S1 — residual, scaling, perturbation

`S11ScatteringSensitivity.residual` is `A x - b`. `residual_identity` / `residual_bound` give

\[
x - A.\mathrm{symm}\, b = A.\mathrm{symm}(Ax-b), \qquad
\|x - A.\mathrm{symm}\, b\| \le \kappa\,\|Ax-b\|
\]

with `κ` an explicit bound on `‖A.symm‖` and `A : X ≃L[ℂ] Y` supplied. `residual_zero_iff` characterises exact solves. Invertibility is never taken from rank, singular values, or a condition number.

`balanced A R D z = R(A(D.symm z))` is `K = R A D⁻¹`. `balanced_residual` is the native identity `K(Dc) - Rb = R(Ac - b)`. `unscale_error` is the `D⁻¹` conversion from `z`-error to `c`-error.

`perturbed_residual_bound` reuses reviewed T2 `inverse_norm_bound` and `solution_error_bound` for supplied `E`, `B = A + E`, and `κε < 1`. Both `A` and `B` are equivalences; `B` is not reconstructed from floating-point data. The bound keeps the numerical residual of `B` at `b+db`, the RHS term `‖db‖`, and `ε‖A.symm b‖`. T2’s `CompleteSpace` existence theorem is not silently used.

`κε ≥ 1` is only an excluded estimate domain. `critical_inverse_margin` is a scalar witness that `0·x = 1` has no solution; it does not say every `κε ≥ 1` operator is singular.

### S2 — affine observation and full current

`amplitude C d x = Cx + d` matches `a = P V⁻¹(T D⁻¹ z - t_in)` once `C` is taken in the same unknown as the residual. `amplitude_difference` drops the fixed offset; `flux_*` is evaluated on the full amplitude, so `d` remains in the expansion point. `mass` is `∑ᵢ‖aᵢ‖`; `observationBound` is `∑ᵢ‖Cᵢ‖`; `pair_bound` is the entrywise bound `|conj(x)ᵀ J y| ≤ β mass(x) mass(y)`. That is an `ℓ¹` / entrywise theorem, not a Euclidean spectral-norm theorem.

`flux_change_bound` keeps both interference terms and the quadratic `mass(e)²`. `current_change_bound` / `flux_and_current_change` take a separate entrywise bound `γ` on `H`; there is no derived physical estimate of `H`. No Hermitian, positive-current, or nonzero-baseline hypothesis is imposed on the sensitivity algebra. Empty index types are allowed; `empty_amplitude_positive` is `mass = 0` on `Fin 0`.

`Pipeline` composes residual → amplitude mass → fixed-`J` flux → fraction. That is the required chain.

### S3 — signed fractions and margins

`fraction_difference` is the exact split into numerator change and denominator change. `fraction_error_bound` uses a common positive lower bound `d` on `|j|` and `|jh|`, with the `d²` factor on the denominator term. Signed numerators and denominators are admitted through absolute values; positivity and conservation are not inferred. `denominator_margin` / `denominator_nonzero` turn a reference magnitude and a strictly smaller error into a nonzero perturbed denominator. `zero_margin_witness` is the case `|0-1| = |1|`. `denominator_cases` is the disjoint split `j = 0 ∨ 0 < |j|`. Lean’s totalized division is constrained by those hypotheses.

Native incident current is the signed diagonal `-Re(J_{ii})` on incoming modes, required positive in `construct` after the selected AST. Selected execution stops at `incoming_flux`. The saved `outgoingFluxRatio` and the positivity `require` are source-only.

### S4 — identification, controls, audit

Compact identification is in the coverage/fidelity records plus the selected AST ranges (10 + 10 + 2 statements). The audit root prints 51 declarations; the recorded axiom lists are `propext`, `Classical.choice`, `Quot.sound` only.

Recorded run6, as reports in the packet: 8 isolated objects, 18 paired rejections, 24 positives (23 distinct; cross-term and quadratic-error share `full_error_witness`), plus native revalidation = 51 records. Each embedded mutant output is exit 1, one `contract_control` error, one `⊢ False`, and no warning/import/resource text. Controls are instance-level canonical arithmetic, not source-replacement mutants. Native run1’s 13 identities, 14 wrong-formula controls, and 9 bound examples are hash-revalidated, not re-executed. Run1’s formal ERROR is kept as a failure.

`lakefile.toml` has no `S11ScatteringSensitivity` target. Direct Mathlib import records are 13 pairs; the manifest lists 15 package pins. Author-side claims of 605 preserved files and 243 old objects are catalogued in `_measurements/S11_lean_sensitivity_preserved_inputs.json`; those files and objects were not re-hashed here.

## Blocking findings

None for the bounded S1–S4 statements, hypotheses, native correspondence, or recorded control contract as supplied in this packet.

## Optional / nonblocking

1. `Pipeline.normalized_residual_bound` still asks the user to supply `margin ≤ |jh|`. `denominator_margin` already proves `|j| - η ≤ |jh|`; a derived form with only `η < |j|` would be a convenience lemma, not a new obligation.

2. Unscaling and `H` stay off the residual-to-fraction pipeline. An application in coefficient space must use `C_c` and `unscale_error`; an application in balanced unknowns must use `C_z = C_c D⁻¹`. `H` is a separate additive term.

3. Formal mutants are wrong constants on canonical witnesses. They reject the intended omitted factor or sign; they do not mutate the proof terms of `residual_bound` or `flux_change_bound`.

4. `complex_current_witness` uses a Hermitian `J`. The theorems do not assume that; the off-diagonal control is still the diagonal-only value 2 versus the full value 4.

5. Native bound examples use unscaled `M`, coefficient `∞`-norm, inverse row-sum `0.625`, and `C_c` with `∑_{i,j}|C_{ij}|`. They do not numerically exercise `‖K⁻¹‖` or `C_z`. The Lean scaling identities are separate, as documented.

## Fidelity and verification limits

This review read the isolated packet only. Lean, NumPy, and scientific constructors were not executed. Network, original-repository, and subagent access were not used.

Outside this packet, and therefore not independently re-established here: run6 `.olean` bytes and raw logs (mutant diagnostics were taken from embedded `output` fields), cgroup sample files behind the resource-receipt hashes, Mathlib/Physlib caches, and the 605/243 historical files/objects named only by hash. Native floating-point comparisons (tolerance `1e-12` scaled, wrong-formula separation `> 1e-8`) remain outside Lean’s kernel.

The packet does not establish a physical inverse bound, a saved physical residual-to-flux certificate, certified numerical intervals, discretization/quadrature/tail error, continuum scattering convergence, channel completeness, conservation, parent-theory accuracy, or a native zero-denominator policy. Origin rephasing in `S11c_d_finite_scattering_domain.py:observable` (`a ↦ D a`, `J ↦ D⁻ᴴ J D⁻¹`) is provenance for the already reviewed F1–F4 pullback; it is not executed by the selected sensitivity instrument.

Policy still requires a second independent non-author review before fidelity review is complete. No commit is authorized.
