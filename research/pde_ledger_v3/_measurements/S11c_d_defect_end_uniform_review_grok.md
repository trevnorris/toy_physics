I'll read the method and evidence guide first, then the supplied source and evidence, and assess only the bounded weak-end comparison.The named guides are not at the workspace root. I'll locate the packet files from the work directory itself.The method and evidence guide are in the packet. Next I’ll check the source interfaces that the comparison actually joins.The bounded comparison method is sound. It compares the corrected translated end symbol with the saved selected uniform equation on the physical wave relation, and it stops when a join, a denominator, or a selected residual is unresolved. That clearance is only for this method. It does not accept a later run, a loss claim, calibration, or the draining model.

## 1. Representation joins

The objects that have to meet are already different in the saved interfaces.

The weak cells are indexed by rows `(U0, U1, U2, THETA_BALANCE, E_W_BALANCE)` and fields `(u_1, u_2, u_3, theta, e_W)`, with four grades and two spatial sides (`source/ends-worker.py` lines 19–20 and 567–578). Each nonzero cell stores the coefficient of one `epsilon_shape`. Address 8034 shows that split directly: `consumerOriginal` still contains `epsilon_shape`, while `consumerField` is the extracted coefficient (`controls/new-control-omit-height-contact.json`). The pairing factor is outside the matrix: `B = 2π ∫ hat(v)(-p)^T E hat(u) dp` (`controls/translated-weak-end-conclusion.json`).

The uniform side is the material-end restriction bundle. `P` is `CLOSED_PENCIL_LEGS[0]` after binding, `L` is `curl[:, 1:3]` over two zero `theta`/`eW` rows, and `D` and the five-row residual are built from that same `P` (`source/uniform-worker.py` lines 215–248 and 323–357). Field order in the lift is `(u1, u2, u3, theta, eW)`. The algebraic radical is `frequency['radical']`, mapped to `uniformAlgebraicRadical`; physical depth is a separate symbol. The conversion checked there is `Q = s(cs) q`, with both wave polynomials required to be degree 2 and free of a linear term. A saved point shows the split numerically in names only: normal `-sqrt(595)/10`, physical depth `3*sqrt(266)*I/200`, algebraic radical `3*sqrt(399)*I/199` (`uniform/selected-result.json`, first LEFT block). Those floats are not operands.

Spatial side joins through the tanh ends: weak `minus`/`plus` heights are `0` and `1/2`, and `W_0 = 1`, so the plus height is `W/2` (`ends/reversed-limits.json`, `controls/translated-weak-end-conclusion.json`). Both bulk faces stay inside each spatial end. The lower-face normal in the minus trace is `-I q_o` (`ends/native-constant-height-minus.json`). Memory in the weak symbol is the saved `beta = (30+9i)/109`, from `rho_m = 1/10`, `Lambda_A0 = 1/100`, `tau_A = 1/10`, and `omega = 3` (`ends/depth-domain.json`, `source/ends-worker.py` lines 349–351).

The written rules block a fitted match. A row or field change has to come from those native maps, as an exact permutation or a constant nonzero unit factor with its inverse. A missing map stops the comparison. Momentum-dependent row mixes, dropped rows, and invented scales are outside the method (`method.md`, coordinate-join section). Epsilon is joined as the already stored single-epsilon coefficient; a zero cell is left alone. The `2π` factor stays on the integral. `p` and `q` join to the saved physical normal momentum and physical outgoing depth, and the stored `Q = s(cs) q` identity is reused only with its own inputs after `omega = 3` and `cs` match. The Euclidean lift Gram stays the bookkeeping Gram of `D`. It is not turned into a physical energy pairing.

The same raw materials have to match before binding, including live memory. An unexpected free symbol is an unresolved source comparison.

## 2. Grades, full Delta, and the selected residuals

Assembly order is the right one. `E_e` is summed over `G = {00, 10, 01, 11}` with explicit zeros, and only then evaluated at `eta = 1/100` and `sigma = 1/1000` (`method.md` lines 12–13 and 84–91). That sigma is the saved relation `sigma_W = W_0 eta_bg / L_W`. Early substitution is refused, so a missing grade cannot disappear into a number.

`Delta_e = Ephys,e - P_old,e` is a new full five-by-five output at one shared physical origin, real frequency, and physical `q`. The method keeps a nonzero or unknown Delta as a result. A retained four-grade identity is not treated as equality to the untruncated finite end law, and a pencil that is already numerically bound is not regraded.

The selected tests are the five-row products

`A_e = (Ephys,e - P_old,e) L_old`, `R_e = Ephys,e L_old - L_old D_old,e`.

They use the restored `P L - L D` bundle and its operands. The new products are the `Ephys` multiplications. On the wave surface these two residuals are linked by that restored identity, so they are one selected statement recorded in two places. They are still the right statement: agreement of `Ephys` with the saved pencil on the two lift columns, and reproduction of the saved `D`. A nonzero `A` or `R` stops the applicability claim. The corrected operator is left as saved. There is no mode search and no new speed scan.

Grade-wise `E_ab L` can be kept for diagnosis. A zero translated response in grades `01` and `11` does not delete source or consumer grades. An optional old-source re-expansion is not required for `A` and `R`.

## 3. Wave-surface and grazing certificates

The modulus written in the method is the outgoing relation at `omega = 3` and tangentials `(1/5, 1/10)`: `|k_∥|^2 = 1/20`, so `q^2 = 9/cs^2 - p^2 - 1/20` (`ends/depth-domain.json`). The branch is the same one the uniform substitution uses: positive real `q` when the radicand is positive, positive imaginary `q` when it is negative (`source/uniform-worker.py` lines 594–604).

For a nongrazing rational entry, the method clears the declared denominator, divides the numerator by that quadratic, and keeps numerator, denominator, quotient, remainder, and the reconstruction residual. Coefficients that still depend on `q` stay unresolved. There is no generic nonzero-denominator claim across an unseen zero set, and no numerical tolerance.

Grazing uses the cells’ already closed `q = 0` values (`source/ends-worker.py` lines 573–574) together with the inherited two-path limits. Those limits are `t → 0+` on the radiating and evanescent paths, with leading-term orders, and they certify the selected projections (`source/uniform-worker.py` lines 651–715). Full-`P` finiteness is only diagnostic there. The method therefore extends `A` and `R` to grazing through those closed and inherited operands. A direct `q = 0` substitution into a raw singular `P_old` is outside the method. If that limit join is unavailable, grazing correspondence stays unresolved even when the nongrazing identity holds. The old limit calculation is not repeated.

The listed branch values match that wave relation and the saved schedule. `cs = sqrt(3/2) = sqrt(6)/2` is LEFT modal ratio 1, with `p = ±sqrt(595)/10`. `cs = sqrt(150/101) = 5*sqrt(606)/101` is RIGHT modal ratio 1, with `p = ±sqrt(601)/10`. They are inherited correspondence targets. The sample LEFT row also carries normal `-sqrt(595)/10` at a nonzero evanescent depth, so that normal is the selected momentum of `D` across the LEFT speeds. The operator limit at `q = 0` remains the separate grazing join.

A nonzero polynomial remainder is not, by itself, a positive-sheet witness. A reported witness has to carry its outgoing root, domain, and an exact finite nonzero certificate.

## 4. What a zero selected residual would mean

If the joins hold and `A = R = 0` on the stated domain, the corrected weak-end symbol and the saved pencil give the same source-selected transverse action at real `omega = 3`: `Ephys L = L D` in that frame. The inherited schedule points above can be cited as algebraic correspondence of that selected equation.

That result is silent about a frequency neighborhood, a differentiated group velocity, a new current, a full mode count, an outgoing basis, and plane-wave action of the nonuniform weak operator. The saved Gram and face currents remain observations about the original uniform law. The projected point table already records `completeOutgoingBasisClaim: false` and an extra full-nullity estimate at the sample point; this method does not turn that table into a census. Defect loss and protection stay out of scope even if both ends agree. Full matrix equality is a different question from the selected residual, and the uniform channels were conditional on the nonuniform and direct mixed-grade composition. Agreement on the lift leaves every other sector where `Delta` puts it.

## 5. Controls, restoration, and a mismatch

The three controls are sensitivity checks on this comparison. They do not establish the physics.

1. Omit one saved term whose product with the inherited lift is nonzero, and keep an exact residual at a declared admissible point, with its local, pressure, and grade ancestry. The transverse lift has zero `theta` and `eW` rows, as in the sample normalized basis. Address 8034 is an `e_W` pressure address, so its column does not move `E L`. If no term with a nonzero lift action is available, the coverage gap is the result.

2. Reverse the spectral normal momentum inside the inherited curl lift while the physical momentum in the equation stays fixed. The control passes only with a changed five-row residual. If symmetry makes it silent, it is inapplicable. A curl-entry ablation is then only an addressing check.

3. On a nonzero pressure-bearing `THETA` or `eW` unit column, replace the positive outgoing depth by the opposite sheet at one rational point in the `cs` interval, and test the closed row action. That unit column is not a mode. The selected transverse lift is allowed to stay sheet-blind because the saved `D` was required to be free of the acoustic radical.

Completed weak-end controls stay restored when their operands and map match. Their late failure is real and stays in the record: `controls/final-serialization.stderr` dies on a SymPy `Zero` during the final JSON write, and `controls/partial-final-checks.txt` is the truncated bytes. `controls/saved-weak-end-controls-return.json` is the earlier durable return. The method does not rerun that job to erase the failure.

Uniform `P`, `L`, Gram, `D`, residual, and the wave conversion are restored from the completed restriction and lift returns named in `uniform/opaque-interface-receipts.json`. Those payloads were not decoded here, and the method does not treat the hashes as the symbols. `uniform/selected-result.json` is a metadata projection (`packet-index.json` gives it a different hash from the source checks file) and is correspondence evidence only. A substantive nonzero selected residual, an unmatched binding, or a missing coordinate join stops the comparison and is kept.

## Nonblocking notes

The spatial end, the harmonic partner `Pminus`, and the normal-momentum sign are three different saved labels. The atomic operand for `A` and `R` is the material-end bundle whose `P`, `L`, `D`, and five-row residual were stored together. The curl-entry fallback needs an explicit movement requirement on an entry that actually enters `curl[:, 1:3]`; otherwise that fallback is inapplicable too. The receipt subset in this packet does not list `LEFT|RIGHT/exact-limits/±1`. Grazing stays unresolved until those completed returns are source-joined, and the method already forbids recomputing them. The later worker’s 4 GiB figure is a tooling parameter; the old uniform containment check required 8 GiB. Neither note changes the algebra of this comparison.

## Coverage and limits

I read `method.md`, `evidence-guide.md`, `review-prompt.md`, `packet-index.json`, `physical-input.json`, the uniform plan, input manifest, opaque receipts, and selected-result header plus one LEFT momentum block. I read the weak-end method, depth domain, translation response, first-shape join, phase arguments, reversed limits, flat closure, constant-height minus trace opening, symbol schema opening, and the control return, conclusion, interval predicate, stderr, partial checks, operation index, and the opening of the height-contact control. I read `uniform-worker.py` through its science path, the continuation through `continue_science`, `ends-worker.py` through symbol assembly, `uniform-source.py`, `pairing-check.py`, `current-source.py`, the frequency-source end binding, and the control-continuation tail.

I did not decode any opaque pickle. I did not read all 200 cells, the full selected-result table, `exact-evidence.jsonl`, `engine.py`, `c1.py`, `c2.py`, or `uniform-response.py`. No peer report was in the packet. The branch-target arithmetic above is a consistency check of the written wave relation against the stated speeds. It is not a new evaluation of the saved pencils.

CLEAR FOR THIS BOUNDED WEAK-END UNIFORM-COMPARISON METHOD