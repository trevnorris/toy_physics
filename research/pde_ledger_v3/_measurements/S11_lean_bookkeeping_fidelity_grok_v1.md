## Verdict

**CLEAR** for the bounded P1–P4 retained-order bookkeeping contract.

The Lean statements match the supplied finite rectangle, the real path, the full degree-six contraction, the field quotient, and the observable split. The recorded run7 controls and the selected native AST checks support those statements at the level the contract claims. Nothing below requires a proof or contract change.

## What the statements establish

**P1.** `rectangle` keeps `eta` and `sigma` independent: `a00 + eta·a10 + sigma·a01 + (eta*sigma)·a11`. `delta` is that expression without `a00`, so it contains the zero-jet term `eta·a10` together with the first-jet terms. `rectangle_path` specializes to `eta = t`, `sigma = r*t` and gives `a0 = a00`, `a1 = a10 + r·a01`, `a2 = r·a11`. Real scalars enter through `ofReal` casts. `quadratic_exact` and `flux_epsilon_squared` use `conj_ofReal`, so conjugation fixes those real powers and the real epsilon factor. `path_does_not_identify_rectangle` puts `a10 = -r·v`, `a01 = v`, `a11 = 0` in the kernel of that path for every real `r` and `t`. The ratio `r = W̄₀/L_W` is the §1c/§3c identification; the theorems quantify over a supplied real `r`.

**P2.** `pair` is `star x ⬝ᵥ (B *ᵥ y)`. For degree-two `current` and `amplitude`, the seven coefficients are the ordered sums `q_d = Σ_{i+j+k=d} pair B_j a_i a_k`. The counts are 1, 3, 6, 7, 6, 3, 1, totaling 27. Both cross orders, the `B1`/`B2` variations, and the `a2` slots are present. On the scalar fixture `a = (1,2,3)`, `B = (1,4,5)` the coefficients are `(1, 8, 31, 72, 107, 96, 45)`, matching `first_witness`, `second_witness`, `higher_witness`, and the native degree list. `retained_with_remainder` rewrites degrees 3–6 as the exact factor `t^3 * (q3 + t q4 + t^2 q5 + t^3 q6)`. `real_flux_expansion` takes real parts. The formal `current` is the supplied degree-two polynomial. There is no Hermitian hypothesis. `Fin 0` is allowed, and the empty-amplitude positive gives flux zero.

**P3.** Over any field, with `j0 ≠ 0`,

`c0 = n0/j0`, `c1 = (n1 - j1 c0)/j0`, `c2 = (n2 - j1 c1 - j2 c0)/j0`

satisfy the three coefficient equations, and `quotient_unique` says they are the only solution. `quotient_residual` is the exact identity

`(j0 + t j1 + t^2 j2)(c0 + t c1 + t^2 c2) - (n0 + t n1 + t^2 n2) = t^3 (j1 c2 + j2 c1) + t^4 (j2 c2)`.

`zero_leading_obstruction` says `0 * u ≠ n0` whenever `n0 ≠ 0`. `zero_model_two_solutions` gives two distinct constants, `0` and `1`, for the all-zero leading equation `0 * u = 0`. Higher-order singular quotients stay outside the contract. `negative_leading` checks `c0 2 (-2) = -1`. The complex positive instantiates the same equations in `ℂ` on real numerals. `epsilon_cancels` assumes only `eps ≠ 0`. `scaled_denominator_nonzero` is the separate claim that `eps ≠ 0` and `den ≠ 0` imply `eps^2 * den ≠ 0`.

**P4.** `subtracted_total` keeps both cross terms. `subtracted_eq_induced_iff` says total-minus-baseline equals the induced flux exactly when the real part of those two crosses vanishes. With `a0 = 0`, `induced_low_coefficients` gives `q0 = q1 = 0` and `induced_coefficient` gives `q2 = pair B0 a1 a1`, independent of `a2`, `B1`, and `B2`. The witness `B0 = 1`, `a0 = 1`, `a1 = 0`, `a2 = 3` gives `q2 = 6`. The scalar flux witness is total-minus-baseline `3` against induced flux `1`. No theorem states a physical baseline, a conversion number, conservation, channel completeness, a scattering solve, or a strong-edge bound.

## Native link and controls

`current.multiply` convolves every supplied bidegree by matrix product, with no grade cutoff. The instrument runs that function, together with `quadratic`, `subtract`, `amplitude_bookkeeping`, `quotient`, `adjoint`, and `lambda_series`, on small synthetic values. `lambda_series` sends bidegree `(p,q)` to degree `p+q` with weight `ratio^q`. For ratio `2` and components `(1,2,3,4)` the path is `(1,8,8)`; the wrong power `ratio^2 * 4 = 16` is rejected. The two-channel fixture compares all seven degrees with an independent triple loop, checks the evaluated contraction at `t = 1/2`, and rejects the diagonal-only current. The scalar path gives total-minus-baseline degree 2 equal to `26` and induced degree 2 equal to `4`.

`quotient` solves one scalar recurrence per diagonal entry of each incident column. It is not a matrix inverse. The two-column fixture gives `c0 = (1,1)`, `c1 = (1,1)`, `c2 = (1,3)` with zero recurrence residuals, and a separate fixture keeps `i/2`. There is no zero guard. The isolated input `1/0` is recorded as nonfinite output and is not tied to a saved physical denominator.

Run7’s report has 44 passing records: one native hash revalidation, seven fresh objects, sixteen rejections, and twenty positives. The audit root’s 33 `#print axioms` lines depend only on `propext`, `Classical.choice`, and `Quot.sound`. Each of the sixteen mutants exits 1 with exactly one `contract_control` error, unsolved `⊢ False`, and no warning. Three `q2 = 31` pairs share `second_witness`, leaving eighteen distinct positive statements. The four extra positives cover the complex coefficient equations, the `n0 ≠ 0` obstruction, the two distinct zero-model constants, and empty flux. Native run1’s report has 30 identities, 11 wrong-formula controls, and one invalid-domain witness. Run7 rechecks its hashes and does not rerun it. The run1 guard exit 1 is recorded as an unrelated formal failure and is not counted as a mutation. Fifteen package commits match `lake-manifest.json`. Six direct Mathlib source/object pairs are named. `verify.py` does not register this contract, and `lakefile.toml` has no bookkeeping library target; the packet does not claim either registration.

## Blocking findings

None.

## Optional observations

These do not change the bounded claims.

- `path_does_not_identify_rectangle` shows membership in the path kernel for every `v`, including `v = 0`. Coverage already declines to infer a nonzero kernel vector when the index type is empty. For a nonempty type the same identity is the non-identification, by taking `v ≠ 0`.
- `leading_denominator_cases` is excluded middle. The real case split is the `j0 ≠ 0` block, the `n0 ≠ 0` obstruction, and the all-zero leading witness.
- Under Lean’s totalized division, `epsilon_cancels` remains true when `den = 0`, because both sides become `0`. The coverage already requires `scaled_denominator_nonzero` before reading the result as a defined fraction.
- The Lean `path_ratio` mutant rejects the false rectangle value `25` by `norm_num` on `path_witness`. The direct wrong-power check is the native comparison of path degree 2 with `16`. That split matches the disclosed witness-or-arithmetic control rule.
- `parent_second_witness` changes the retained coefficient `a2`. It shows how a second-order amplitude moves `q2` when the baseline is nonzero. It is not an omitted pure `η²` or `σ²` parent term, and the contract does not promote it to one.
- The native two-channel current fixtures are Hermitian, with nonzero off-diagonal entries. The Lean contraction has no Hermitian hypothesis. The off-diagonal rejection still separates the full product from its diagonal truncation.
- The Lean `ℂ` positive uses real numerals. Retention of the non-real value `i/2` is the native fixture.
- `Quotient.lean` imports `Mathlib.Data.Complex.Basic` without using `ℂ` in that file. The complex positive lives in the control script.

## Fidelity and verification limits

This review read the packet and compared the statements with §3c–§3d and the selected native functions. It did not execute Lean or NumPy, recompute SHA-256 against trees outside the packet, or inspect the other reviewer’s report. The author PASS label was not treated as evidence.

Outside this packet, and therefore unchecked here, are the Mathlib sources and oleans named by the six direct-import hashes, the git cleanliness of the fifteen pins, the 567 historical files, the 236 preserved objects, the cgroup sample log, the pre-repair scratch logs, the engine product `J.multiply`, and the numeric `W_0`/`L_W` JSON. Within the packet, the tested `current.multiply` is unrestricted convolution, and `r` remains a supplied real.

The native report is hash-revalidated run1 evidence on binary-exact fixtures. It does not put NumPy translation or arbitrary floating-point evaluation inside the Lean kernel. The contract does not establish a parent-theory Taylor bound, a physical baseline or conversion value, conservation, channel completeness, a scattering solution, or a strong-edge extrapolation.
