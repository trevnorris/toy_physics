The F1–F4 contract is faithful to its bounded claim. The statements quantify over every finite complex amplitude space, keep the full conjugate contraction and both cross terms, and put Hermiticity, positivity, normalization, and balance only where the theorems name those hypotheses.

This is a statement-fidelity review of packet `5945c5b4…84ae`. The recorded run8 logs were read. Lean was not rebuilt, and the native modules were not executed. The packet has no compiled objects and no Mathlib cache, so this is not a fresh kernel replay.

## F1 — full finite current form

`pair` in `lean/s11/S11ScatteringFlux/Current.lean` is `star x ⬝ᵥ (J *ᵥ y)` on every `Fintype` index type, including `Fin 0`. `pair_coordinates` expands it to `∑ j, ∑ i, star (x i) * J i j * y j`, so every matrix entry enters. There is no diagonal, positivity, transverse, nondegeneracy, or basis-completeness hypothesis on `pair`, `flux`, `flux_add`, or `pair_pullback`.

`flux J x` is defined as `(pair J x x).re` for every supplied `J`. `hermitian_flux_real` proves `(pair J x x).im = 0` only from `Jᴴ = J`. `pair_hermitian` is the off-diagonal form of that hypothesis: `pair J x y = star (pair J y x)` for all `x, y`. The non-Hermitian witness is the constant matrix `Complex.I` on `realAmplitude = 1`, where the raw imaginary part is `1`. The paired mutant claiming that imaginary part is `0` is rejected. Taking the real part in `flux` does not erase that raw defect.

`flux_add` keeps both cross contractions for every matrix:

`flux J (x+y) = flux J x + flux J y + (pair J x y + pair J y x).re`.

`flux_add_iff_cross_zero` makes additivity equivalent to that real sum vanishing. Nothing in the statement deletes a closed matching mode or assumes a diagonal sector split.

`empty_flux` sends every current and every amplitude on `Fin 0` to `0`. `null_flux_nonzero_amplitude` sends `diag(1,-1)` on `both = ![1,1]` to flux `0` while proving `both ≠ 0`. Zero flux is therefore not collapsed to the zero vector.

The same expansion is what `direct_quadratic` writes in `_measurements/S11c_d_continuum_currents.py`: `conj(x[i]) * J[i,j] * y[j]` over both indices. Summation order differs and the summand does not.

## F2 — pullback and covariance

`pullback J C = Cᴴ * J * C` is stated for rectangular `C : Matrix n m`. `pair_pullback` and `flux_pullback` identify the whole form: `pair (pullback J C) x y = pair J (C*ᵥx) (C*ᵥy)`. `pullback_comp` is `(C*D)ᴴ J (C*D)`, with no invertibility assumption. A rectangular selection therefore changes coordinates on the column span only. No theorem promotes it to an invertible change of a complete physical channel basis.

`scattering_flux_covariant` is narrower and the quantifiers match the comment. `Cout` and `Dout` are square on the outgoing index, `S : Matrix n m`, `Cin : Matrix m k`, and the only inverse hypothesis is the supplied `Cout * Dout = 1`. The identity is

`flux (pullback J Cout) ((Dout * S * Cin) *ᵥ x) = flux J (S *ᵥ (Cin *ᵥ x))`.

That is outgoing-flux covariance for the represented input `Cin x`. It does not transform an incident denominator, assert that `S` solves the PDE, or invoke the §5a Eulerian/material route. `[DecidableEq n]` is the Mathlib requirement for the identity matrix.

Shared physics §3a defines the current matrix as the polarized bilinear `J_ab`, possibly non-Hermitian, and says a bare `√(v_out/v_in)` factor is not a substitute. The Lean pullback inserts no such factor. `basis_flux` shows the real scaling `C = 2` takes flux from `1` to `4`.

## F3 — signs, fraction, conditional balance

Section 2 identifies the minus end with the left background and the plus end with the right. Section 3a then sets `s₋ = -1`, `s₊ = +1`,

`J_in = -s_e J_n`, and `J_out = Σ_e s_e J_n,e`.

`Balance.lean` matches that reading. `outward left = -1`, `outward right = 1`, `incident e j = -outward e * j`, and `outgoing l r = -l + r`. So `incident left` is the identity, `incident right` negates, and the two-end total is the signed sum. `oriented_outgoing_witness` gives `outgoing 3 1 = -2`. The mutant `4` is the unsigned sum `3+1`.

`end_coverage` exhausts the two `End` constructors. `flux_sign_coverage` is `j < 0 ∨ j = 0 ∨ 0 < j`. The coverage note correctly leaves exclusivity to the real order laws. There is no separate local disjointness lemma. That does not create an overlap in the stated trichotomy.

`fraction num den` is `none` exactly when `den = 0`, and `some (num/den)` otherwise (`fraction_undefined_iff`, `fraction_defined`). The zero test is exact. `fraction 7 0 = none` rejects the Lean junk value `some 0`. `null_flux_nonzero_amplitude` keeps a nonzero vector whose flux is zero, so a vanishing denominator is not an amplitude-zero shortcut.

`fraction_nonneg` needs `0 ≤ num` and `0 < den`. `fraction_le_one_iff` needs `0 < den` and is equivalent to `num ≤ den`. Negative denominators stay in the algebra: `fraction 2 (-4) = some (-1/2)`, and the mutant `some (1/2)` is rejected. `fraction 4 2 = some 2` rejects an unconditional cap at `1`.

`conditional_balance` assumes `converted + survived + defect = incoming` and `incoming ≠ 0`, and concludes

`converted/incoming + survived/incoming = 1 - defect/incoming`.

`conservation_requires_zero_defect` is the converse direction under those same hypotheses: the normalized sum equals `1` exactly when `defect = 0`. There is no `S`, no `SᴴS = 1`, no omitted-channel construction, and no bound-capture term. Section 3c says bound spectral overlap is not inserted into the continuum current and that a total photon-loss probability is not formed. The theorem leaves that accounting identity as an application premise. `balance_positive` instantiates it at `1+1+2 = 4`. `zero_defect_positive` instantiates defect `0`.

## F4 — controls and adjudication

`_measurements/S11_lean_flux_contract_check.py` builds fresh statements. It does not patch canonical sources. `verify.adjudicate` accepts a rejection only when the exit code is `1`, the bad-instrument pattern is absent, and the single error sits in `contract_control` with the substring `⊢ False` exactly once and with no warning. Import failures, unreduced constructor goals, timeouts, and resource failures fail that test.

The run8 report records all eleven mutants in that form. Each output is one `error: unsolved goals` followed by `⊢ False`.

| Pair | Proved statement | Rejected statement |
|---|---|---|
| interference | `flux coherent both = 4` | `= 2`, the diagonal-only value |
| conjugation | `flux oneCurrent phaseAmplitude = 1` | `= -1`, the product without left conjugation |
| basis metric | pulled flux `= 4` | `= 1`, the unscaled current |
| orientation | `outgoing 3 1 = -2` | `= 4` |
| signed incident | `incident.left (-3) = -3` | `= 3` |
| zero denominator | `fraction 7 0 = none` | `= some 0` |
| normalization | `fraction 2 4 = some (1/2)` | `= some 2` |
| negative domain | `fraction 2 (-4) = some (-1/2)` | `= some (1/2)` |
| unit bound | `fraction 4 2 = some 2` | `= some 1` |
| conservation | `(1)/4 + (1)/4 = 1/2` | `= 1` |
| raw reality | imaginary part `= 1` | `= 0` |

The eleven matching positives, plus `empty_positive`, `indefinite_null_positive`, `balance_positive`, and `zero_defect_positive`, are the fifteen passes. With four builds and the native check, the report’s 31 records are the full set. The audit root prints axioms for 34 declarations. Every list is among `propext`, `Classical.choice`, and `Quot.sound`. `end_coverage` depends only on `propext`.

These controls are meaningful for this contract. The false right-hand sides are the adjacent wrong conventions: drop the off-diagonal block, drop conjugation, ignore the basis factor, add absolute end currents, flip the left incident sign, treat division by zero as zero, invert or absolutize a ratio, clamp at one, or declare the raw imaginary part absent. The zero-denominator mutant is specifically `⊢ False`, which is the repair described after earlier runs left `none = some 0`.

The conservation pair is the loosest of the eleven. Its script entry has an empty lemma list, so `norm_num` refutes `(1)/4+(1)/4 = 1` by arithmetic and does not call `conditional_balance`. The theorem is still exercised by `balance_positive` and `zero_defect_positive`, and `conservation_requires_zero_defect` is among the audited statements. That is enough for this numeric contract. It is not a source mutation of the balance proof.

## Native grade-(0,0) link

`_measurements/S11_lean_flux_source_check.py` extracts `multiply` and `quadratic` from the currents module, `adjoint` and `gram` from the response module, and `RectangularModeJets.multiply` from `scripts/S11c_d_mixing_scattering_sympy_audit.py`. It strips decorators, does not import the modules, and runs only the small fixtures. The recorded report has seven identities and five wrong-formula controls. The adjoint mutation is one of those five: `RemoveConjugation` deletes `.conj()` from the extracted adjoint, whose source is `{g: v.conj().T …}`.

`quadratic(a, J)` is `adjoint(a) * J * a`, hence `a.conj().T @ J @ a`, the complex value of `pair`, with `flux` the real part. On the Hermitian fixture `a = [1, i]`, `J = [[2, i], [-i, 3]]`, that contraction equals `3`. The diagonal part alone equals `5`. `J.T` gives `7`. The non-conjugate product `a.T @ J @ a` equals `-1`. For `x = [1, 0]` and `y = [0, i]`, each cross term equals `-1`, their sum is `-2`, and `q(x)+q(y)+cross = 3` while `q(x)+q(y) = 5`. Both interference terms are required on that fixture. The indefinite vector `[1, 1]` on `diag(1,-1)` gives `0`. The empty channel shape `(0, 1)` against a `(0, 0)` matrix gives the `1×1` zero, the matrix shape of the scalar `empty_flux` result.

`gram` is executed with `RectangularModeJets.multiply`. That function drops a product whose grade has a component greater than 1. On the recorded `(0,0)` inputs the guard keeps the single product, and the independent check is `c.conj().T @ j @ c` for the complex basis `[[2, i], [0, 2]]`. The stale-metric comparison is real: `q(c a, J) = 10` and `q(a, J) = 3`. This grade-(0,0) congruence matches `pullback`. It is not a higher-grade certificate, and the grade cutoff must not be reused as one. Production `response.J` is `boundary.J`, whose definition is not in this packet. The instrument discloses the engine multiply it actually ran.

`open_metrics` multiplies the outgoing block by `s = v['orientation']` and the incoming block by `-s`. The anchor is an `ast.unparse` substring check of that function, not an execution and not a proof that the stored orientation numbers are `-1` and `+1`. Those literals are the §3a convention encoded by `outward`. The native orientation evidence stops at the relative sign.

## Evidence limits that stay limits

Run8’s axiom log, eleven `⊢ False` transcripts, source hashes, and four object hashes agree across the contract report and `_measurements/S11_lean_flux_validation.json`. The toolchain pin is Lean 4.33.0. Seven direct Mathlib imports have source and object hashes. The rest of the external cache is the pinned baseline, not a fresh Mathlib build. Python and NumPy versions are absent from the native report. The fixtures used here are small Gaussian integers, exact in binary arithmetic, so the missing versions do not change these identities. They also do not turn the native script into a pinned numerical scattering run.

`lean/INSTALL_D5_BULK_REGISTRATION.json` is `PASS_REGISTRATION_ONLY`: catalog and tooling checks, with an explicit refusal to replay the D5B proof. The resource receipt is a separate containment record. Neither one enlarges F1–F4.

Shared physics §3a and §3c fix the sign and bilinear conventions above. They do not show that a scattering solution, a unitary `S`, a positive physical current, or a complete channel basis has been constructed.

## Optional improvements

These do not change the bounded claim.

- The conservation mutant could apply `conditional_balance` to a nonzero defect and ask for unity, instead of refuting the arithmetic with an empty lemma list.
- A Lean numeric witness for two nonzero amplitude cross terms would sit beside the native fixture. The general identity is already `flux_add`.
- `incident_right` is proved and audited. A separate numeric mutant is not present. `outgoing 3 1 = -2` already uses both end signs.
- A local lemma that the three sign cases are pairwise disjoint would match the coverage note’s appeal to the order laws.
- Keep `RectangularModeJets.multiply` inside the grade-(0,0) link. Its cutoff is a different operation from `current.multiply` once a grade component exceeds 1.

## Blocking findings

None. The hypotheses that the bounded claim requires are present on the theorems that use them: `Jᴴ = J` for raw reality, `Cout * Dout = 1` for outgoing covariance, `den = 0` for an undefined fraction, `0 ≤ num` and `0 < den` for a nonnegative ratio, `num ≤ den` for the unit bound, and `converted + survived + defect = incoming` with nonzero incident flux for normalized balance. No silent diagonal, transverse, positivity, or completeness restriction enters the finite-dimensional algebra.

CLEAR
