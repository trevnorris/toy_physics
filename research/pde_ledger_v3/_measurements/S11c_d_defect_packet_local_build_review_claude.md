**Verdict: NEEDS REVISION**

I found one blocking defect: the worker's coefficient-certificate check fails on the saved data, so the single authorized run would stop at cell 2, before any quadrature.

## Blocking finding

**`worker.py:158` compares quotient coefficients to their certificates by position, but the saved certificate lists are not in coefficient order.**
- The check is `all(v['finite'] is True and v['value']==c['value'] for c,v in zip(coeff,cert))`.
- Cell 2 (n=0, grade (0,1), degree 4) has coefficients in order c0..c4 (`local-cells.json:1695-1730`). Its certificates run c0, c2, c3, c1, c4 (`:1732-1812`).
- At position 1, `189/550450+…` is compared with `-112/165135-…`. That is unequal, so `require` raises.
- Cell 6 (the Leibniz control cell, `:3242-3337`) is permuted the same way: c0, c2, c3, c1.
- Cells 0 and 1 pass only by luck: cell 0 is a single constant, and both coefficients of cell 1 are equal.
- The failure is fail-closed and preserved, but it would consume the one authorized science execution. The 20 tooling tests never exercise this predicate, so they did not catch it.
- **Correction:** drop the positional pairing. Either compare the two lists as multisets, or match each certificate to its coefficient by value. Also check that each certificate's `real + i·imaginary` equals its `value`. Add a stdlib test over the saved cells for this check. This changes the worker hash, so the gate and review record must be regenerated.

## Coverage

I did not run anything.

**Read in full:**
- build.md, method.md and review-prompt.md
- worker.py and numerical-library.py
- quadrature-library.py and storage-library.py
- launcher.py, tooling-tests.py and the test record
- manifest.json, execution-authority.json and evidence-guide.md
- physical-plan, physical-input, source-units, source-text-joins, rules/extraction and the Fourier journal receipt

**Read in part:**
- **Selected cells:** all 16 were checked structurally for zero flag, degree, epsilonPower, grade order and sourceChildren. Cells 0, 1, 2, 6, 7, 8 and 9 were read in detail.
- **Child batches:** 3 of 24 returns and 1 input, only for the first child or the fields the worker reads. I did not read every child's `nativeCoefficient`.
- **Partition:** the source header, first children and tail only.
- **Context, registry and contract:** context at the top level and `physical`; merged.json keys plus `gamma_14`; the contract dimension keys the coefficients need; the LEFT origin file's head and tail.
- **Runtime helpers:** only the raw-helper segments the worker imports.

**Not read:**
- all-local-cells.json (the 400 cells). I did not verify the "selected equals filtered" claim.
- the three rule files
- the 13 MB source line, so its hash join is unverified
- the guard and supervisor
- the gate, which does not exist yet

So the saved operands were read selectively. Further runtime-predicate failures beyond finding 1 cannot be excluded.

## Checks that held up
- **Gaussian tail inequalities:** I re-derived them.
  - The majorant recursion matches `|ip − (x+5/2)/64| ≤ 325/128 + |x|/64`, valid since kappa < 5/2.
  - The `I0`, `I1` and `Ij` recurrences are correct.
  - `|uv| = exp(−x²/64 − 25/256)` holds exactly, including under the conjugation mutant.
  - The code matches the build text. The Leibniz tail adds the `a'` bound.
- **All-real constant references:** the moments `z0`, `z0²−1/128`, `z0³−3z0/128` are correct and match the code. They cover only cells 0 and 8, the two constant nonzero cells.
- **Route independence:** A builds Q by recursion and evaluates u and v separately. B uses the explicit polynomials and the reduced product. `Q2` and `Q3` are correct. Open squared half-panels have correct Jacobians. The B budget is the global `1e-13/32`, and R+8 is present.
- **Cell structure:** there are 8 zero and 8 nonzero cells. Cell 6 is degree 3 and nonconstant, so the Leibniz control is meaningful. The coefficient `a'` I checked by hand against the saved `firstDerivativePolynomial`.
- **Controls:** cell 0 is nonzero and constant. The conjugation mutant collapses to about 0, so its movement is about the baseline magnitude, far above 10x the envelope. Both controls are saved before the final guards.
- **Unit argument:** this is sound. I hand-checked THETA children 1 and 116: child 116 closes to (−3,−1,1) only with the published registry value `gamma_14 = (2,0,0)`. The shared pairing unit `M L⁻² T⁻¹` follows. The worker applies this predicate to every selected e_W child at runtime, and I could not check the rest.
- **Evidence ordering:** evidence is emitted before each `require` and before the failure guards, including the SQLite records and the control receipts.

## Limitations (not blocking)
1. **Registry origins are not joined.** The four `*-origin.json` files are copied and hashed but never read. Only agreement across the four merged-registry entries is checked, so the build's "producer/object provenance" claim is not enforced in code.
2. **Quotient-to-original join is indirect.** The worker never directly checks `cell['coefficient']` against the quotient polynomial. It relies on the saved cancelled-zero identities, whose operands it does not inspect. My hand checks of cells 0, 1, 2, 6, 7, 8 and 9 agree.
3. **Summand sum is not recomputed.** The sum of child summands is not compared to the identity's right-hand side. This is inherited by design.
4. **`a'` is not cross-checked.** The code's `profile_derivative` is not compared with the saved `firstDerivativePolynomial`.
5. **Speed independence is name-based.** The hidden-speed check looks only at symbol names. The saved cells were bound with the old `c_s0=10`.
6. **Zero-cell counts are not enforced.** Nothing checks the 8/8 split, and a cell flagged nonzero but with all-zero coefficients passes silently.
7. **build.md is not hash-pinned.** It is covered only indirectly, through the review record.
8. **Shared and empirical limits.** Both routes share the coefficient vector and `poly_eval`. The B errors are embedded empirical estimates.
9. **Output size.** I estimate a SQLite journal of a few GB. This is unquantified.
10. **Scope.** No numerical result is claimed, and no local value is a pressure, packet, current or loss value.