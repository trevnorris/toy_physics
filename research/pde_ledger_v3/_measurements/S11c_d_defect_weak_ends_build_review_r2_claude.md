**Independent assessment of the translated weak-end BUILD: no blocker found by static reading.**

This was a static read only. Nothing was executed, so I did not recalculate any saved result. Some run-time assertions cannot be settled from the packet (see limitations).

## Inspected files and regions
- `input/build-guide.md`: read in full, first.
- `input/method.md`: read in full.
- `input/worker.py`: read lines 1–637 in full.
- `input/launcher.py`: read in full.
- `input/runtime-source/inert-helpers.py`: read in full. It defines `Journal`, `decode`, `containment`, `one_symbol` and `expanded_sinh_arguments`, which `worker.py` extracts by AST.
- `input/inputs.json`:
  - Header, `sourcePins` (grepped), `resources` and alias presence.
  - The pins cover the helper, the supervisor, the shared guard, the hook, the launcher, the method, the method record and the authority.
  - I did not read the other ~420 alias entries line by line.
- `input/execution-authority.json`: read in full.
- `input/runtime-source/supervisor.py`: grepped for `--stage` handling only. `defect_weak_ends` is passed through as an arbitrary stage label.
- Saved operands, sampled:
  - `saved-inputs/weak/plus-global-height-PV-certificate.json` and the minus certificate.
  - `weak/new-weak-duality-and-order.json`, `weak/global-profile-envelope.json` and `weak/all-coefficient-certificates.json`.
  - Field `79ea957c…` (constant) and field `833d6feb…` (degree 5): operands, reconstruction input and return, and derivative class.
  - `reference/native-first-shape-input.json`, `reference/plus-new-native-trace.json`, `reference/minus-final-native-slot-routing.json` and `reference/minus-native-height-constant.json`.
  - `reference/restored-profile-arguments.json` and `saved/reference/retained-response-census.json`.
  - `full/new-pressure-wave-arguments.json` and `full/extended-binding-context.json`.
  - `inventory/address-representatives.json`: first record, plus grep of the height responses, wave multipliers and normal maps.

## Findings that support clearance
- **Phase and sign joins:**
  - The restricted `native_phase_exponent` grammar matches the saved `sourcePhase` and `profilePhase` strings, including the 2- and 3-tuple generator targets and `sum(genexp)`.
  - The source shift gives `i(l-k)a`, and the profile coefficient comes out the same.
  - The flipped-`I` mutation moves the coefficient by −2, which is nonzero.
- **PV/Dirichlet algebra:**
  - With `A0=1/(2π)` the limits are 0 and `W/2`, and both are checked against the native `H±`.
  - The reversed phase is derived from `-translated_exponent`. Both height exchanges are checked and nothing is hard-coded to zero.
  - The oscillatory diagonal term is retained, and the exponential is not differentiated.
- **Profile consistency:** by hand, `(W/4)δ + (W/2i)A/Q` is the forward transform of `(W/4)(1+tanh)` under the native `e^{-iQy}` convention. The plus end is therefore `W/2` and the minus end is 0.
- **Holder certificates:**
  - The PV certificates reduce to `Bh` and `i·f·q·Bh`. The plus normal coefficient is positive and the minus one is negative.
  - They carry `normalCoefficientBoundedDirectly: true` and `qTimesTestAssumedC1: false`.
- **Pressure fields:**
  - The certificate polynomial is the quotient in `T`. Substituting `T→tanh(x/10)` matches the old reconstruction right-hand side term by term.
  - The old literal zero return is inherited and no function is replayed.
- **Constant-height specialization:**
  - Lab-height × normal gives `i q η H±` on both faces.
  - The affine slot gives `F0 − product·F0 = F0 + η H Bh`.
  - `η²` is carried separately.
  - Grazing uses the closed `q+β` expressions, with `Re β=30/109`.
- **Address and control joins:**
  - The wave multipliers are polynomials in `composition_p`, mapped to the diagonal `p`.
  - The height responses contain exactly one `reference_height_hat(0)`, and the flat and height formulas reduce to `F0` and `H·Bh`.
  - The controls select from actual nonzero plus-end terms, and the movement enters the full end cell.
  - The control point `p=1, q=2, cs²=180/101` gives radicand 4, which is in scope.
  - The lower-normal control uses a minus-face flat normal-slot address with `-term`.
- **Runtime wiring:**
  - Order is gate, then argv, then helper pins, then `containment` (native `RLIMIT_AS` 4 GiB, swap 0, pids 32, one CPU), and only then the SymPy import.
  - There is no retry and no deadline, and the journal and evidence log persist before the guards.
  - The launcher pre-verifies via `runpy` with stdlib only, arms the hook before the coordinator proceeds, and runs the unchanged guard → supervisor → worker command.

## Blockers
None found, mathematical or tooling.

## Nonblocking limitations
1. **Structural equality at `worker.py:331`.** `actual_argument == D(proof['right'])` is a structural equality, not a cancel-to-zero check. The srepr shows matching coefficient and Add structure, so I expect it to hold. If it fails, it fails loud and early. That would be a pure schema/tooling failure, and it would cost the one authorized run. I checked 2 of the 34 fields.
2. **The "physical first-height DtN vanishes" check is tautological.** Setting `qi=qo=q` in `(−qi+qo)` already makes `first_diag` zero at `worker.py:440`, so the check at `:460` adds no further content. It is nevertheless correct for the stated claim.
3. **`omit-height-contact` subtracts a literal `W/4`** (`:571`) rather than reading the contact from the decomposition. The value equals the method's contact, and the control still enters the actual cell.
4. **Not statically established:**
   - Control movements at the point being nonzero.
   - All 13,260 address joins passing.
   - Heavy `sp.cancel` cost and memory under the 4 GiB limit. File sizes of about 18 MB per row suggest it fits.
   - The minus-face `normalJet` equaling `−i q` (I read only the minus slot and height-constant files, not the minus trace).
5. **The functional-analytic steps are assessed arguments, not machine proofs:** Riemann–Lebesgue, the uniform-in-`a` bound, and the compact-`cs` uniformity. This follows the method's own statement.
6. **Uninspected:** `input/tests.py` and the ~400 remaining saved aliases. I did not read the wave/factor identities, the cell identities, or the full address arrays beyond samples. I did not check that the shipped hashes match the files. `shared-guard.py` was not read. The supervisor was only grepped for stage handling.

No scientific READY gate exists yet, and clearance of this build does not establish any scattering, plane-wave, current or loss result.

**CLEAR FOR THIS TRANSLATED WEAK-END BUILD**