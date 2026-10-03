## Independent build assessment: translated weak-end worker

**Verdict: CLEAR FOR THIS TRANSLATED WEAK-END BUILD.** This is based on a sampled review, not a full read of every saved operand. The coverage and limits are listed below.

This is an assessment of method wiring only. It is not a physics result. I ran nothing, and I did not recalculate any saved result.

### What I checked

**Joins from saved inputs.**
- The context joins hold. `extended-binding-context.json` has `numeric.omega` = 3 (lines 415–418) and a `frequencyOverride` from 1 to 3. The `physical` and `saved.physicalInput` blocks keep the historical omega = 1, which is what `worker.py:298` compares against.
- The worker's own constant for β at Ω = 3 is (30+9i)/109. I re-derived it by hand from Λ = 1/100, ρ = 1/10, τ = 1/10 (`worker.py:349-350`).
- The 28 aliases the worker loads by name are all present in `inputs.json`. I confirmed this by counting matches.
- `all-local-cells.json` has 800 endpoint records, which is 400 cells × 2 ends. `new-pressure-wave-arguments.json` has the `sourceMomentum`/`savedMultiplier` records the address join needs.
- Every address factor proof carries a published `cancelled` zero. `inherit_zero` (`worker.py:136`) only checks that zero and calls nothing, so there is no replay. The worker never compares the two sides of an old cancelled identity.

**Quotient versus numerator, and T versus tanh.**
- The saved field polynomial files hold only the numerator. The certificates and `P0` hold the full quotient.
- The worker evaluates the certificate polynomial (`worker.py:333`) and substitutes T → tanh(x/10). It then compares that to the saved reconstruction `right` side (`worker.py:337-340`).
- The saved `left` is the factored form and `right` is the expanded form. The worker compares only against `right`, as it should.
- The endpoint stage uses the genuine saved T and x symbols from the operands file.

**Phase and translation.**
- `native_phase_exponent` (`worker.py:209`) is a restricted AST grammar with no `eval`. It accepts exactly the `sp.exp(sp.I*(sum(zip…)−sum(zip…)))` source form and the `sp.exp(-sp.I*sum(zip…))` profile form.
- The saved strings are at `new-weak-duality-and-order.json:19-20`, and they are the same strings the grammar was written for.
- With kout = (l, e1, e2) and kin = (k, e1, e2), the shift x, y → x+a, y+a gives exactly i(l−k)a. The edge coordinates cancel.
- The profile exponent comes out as −i(l−k)y. Both signs are then joined (`worker.py:395-397`).
- Flipping the profile's imaginary unit moves the relation by −2, which is nonzero, so the sign control is sensitive.

**Dirichlet limits and heights.**
- The PV limit is W/4 + (W/2i)(iπA₀·sgn·c). With A₀ = 1/(2π) this gives 0 and W/2, and the worker checks that against H.
- The oscillatory diagonal term is carried explicitly (`worker.py:410-413`).
- The reversed-phase limits are derived from the reversed exponent. Both exchanges are checked (`worker.py:415-420`), and nothing is hard-coded to zero.
- The saved native heights are η·w₁/2 on the plus face and −η·w₁/2 on the minus face. At the plus end w₁ = 1, so H₊ = 1/2 = W/2.
- The trace product is i·q·η·H on both faces, because the height and normal signs flip together.
- The address-level height transform at 0 is joined to H.

**Holder and normal coefficients.**
- The PV certificates reduce to B_h at q(k) = q(l) = q, so the pressure coefficient is −iμq/(q+β).
- The plus-face normal coefficient is +μq²/(q+β) and the minus-face one is −μq²/(q+β). The worker's check against i·sign·q·B_h matches.
- Both certificates have `normalCoefficientBoundedDirectly: true` and `qTimesTestAssumedC1: false`, so no C1(qY) assumption is made.

**First-shape grade split.**
- The saved first-shape `right` operand has an η height term and a σ slope term with no mixed term.
- The height term is i·η·ω·(qo−qi)·ĥ/(10·qo). On common depth it becomes exactly zero by cancellation inside the expression, so there is no raw 0/0.
- The slope term k·σ·ŝ/(10·qi·qo) is kept. It carries no pointwise-zero assertion and no j(0) = 0 assumption.
- Full reconstruction and the grade-free coefficient check (`worker.py:449-452`) are fatal if a grade is omitted.

**Native slot and flat joins.**
- The final affine slot is −affine_height·jet_slot + physical_target. Inserting the normal jet times F0 gives F0 + η·H·B_h at both ends of both faces.
- The flat closure μ/(q+β) is joined to the saved census flat expression with ω = 3.
- The census flat has only the symbols `reference_unrestricted_frequency` and `reference_qi`, and both are bound.

**Grades and coverage.**
- The grade-triple count of 16 is mathematically right: each component has 4 options, so 4 × 4 = 16.
- The worker checks the target-grade sum and every address against its coverage record.
- Responses (0,1) and (1,1) get explicit zeros with their analytic reasons, never an early drop.
- The 200 end-symbol cells follow the expected 2 × 5 × 5 × 4 grid.
- Source derivatives are taken at p, and the normal factor uses the output q after the diagonal support.

**Controls.**
- The point p = 1, q = 2 with cs² = 180/101 gives a radicand of exactly 4, which is in [1, 4] and outgoing.
- Controls are drawn from nonzero addresses by existing metadata, and nonzero movement is certified exactly at runtime. A control that cannot be found is fatal.
- The contact-omission control is an actual W/4·B_h·base·normal term at the plus end.
- The reversed-phase control inserts both exchanged endpoint limits.
- The lower-normal control uses a surviving FLAT normal-slot address on the minus face.
- Each control is pushed through the real full-cell symbol.

**Runtime and launcher.**
- `containment()` runs before the sympy import (`worker.py:634-636`). It requires 4 GiB, swap 0, pids 32, one CPU, one thread, and no CPU rlimit.
- The helper is exec'd from exactly the `HELPERS` function and class names out of `inert-helpers.py`. All ten are present there.
- The worker's `verify_gate` is its own function. It requires the exact status, hashes, pins, authority, a literal CLEAR verdict from both reviewers, and a fresh build record.
- The launcher is hook-first. It starts the watcher and waits for its `waiting` state before sending the "1" byte. After that the coordinator runs the unchanged-path pooled-guard command.
- There is no computation deadline and no retry. The 30-second `select` is only an arming handshake.
- Evidence is fsynced before the guards run. `main` writes `posthashes.json` and `checks.json`, and `scientificAcceptance` is False.

### Blockers

None found. I found no mathematical, implementation, provenance, readiness or scope blocker.

### Non-blocking limitations

1. **Structural-equality fragility.** `worker.py:340` requires `actual_argument == D(proof['right'])` as a structural SymPy equality. A representational mismatch (for example unreduced complex rational coefficients) would fail the run, not pass it wrongly. There is no retry, so a false failure would consume the one authorized run. The guide says no test restored a scientific object. Sampled forms look consistent (unreduced `(a+bI)/n` coefficients in both the saved `right` and the certificates), but I could not execute it.
2. **Response cache key.** The cache key `(proof, face, slot)` (`worker.py:528-529`) does not include the response grade. It is safe only if each proof's `responseCoefficient` determines the grade. The address check ties the coefficient to the address, and a mismatch would surface later as an unbound or wrong symbol.
3. **2π factor.** `weakPairingFactor` is a recorded constant, not derived from saved values. The profile normalization is pinned by A₀ = 1/(2π) and the native half-height checks. The (2π)⁻³·(2π)² edge arithmetic is inherited, not recomputed here.
4. **Memory and size risk.** The worker emits whole-row JSON for about 2652 detailed addresses per row, and builds sums of many terms per cell. The 4 GiB cap and the volume have not been exercised.
5. **Analytic steps are not machine-proved.** Schwartz convergence, uniform-in-a bounds, Riemann–Lebesgue and Dirichlet are reviewed arguments from `method.md`. The worker records them as assessed, and it states that. I checked that the method does not differentiate exp(iQa) in the uniform bound. It splits off A(0) and uses the bounded sine integral.
6. **Controls are formula-level.** The control row-identity zeros (`worker.py:602`) are algebraic tautologies. They demonstrate sensitivity only, as the method states.

### Files and regions inspected

- `input/build-guide.md` (whole), `input/method.md` (whole)
- `input/worker.py` (whole, lines 1–657)
- `input/launcher.py` (whole)
- `input/runtime-source/inert-helpers.py` (whole)
- `input/inputs.json`: key lines for worker, launcher, helper, c2 source, physical input, method path, record paths, scope, and resources; plus alias-presence counts for the 28 loaded names
- `input/full/extended-binding-context.json` (lines 1–60, 280–320, 355–449, 560–612)
- `input/full/new-pressure-factor-arguments.json` (head, factor 0, proof names)
- `input/full/new-pressure-wave-arguments.json` and `input/full/all-local-cells.json` (key counts and key names only)
- `input/pressure/new-weak-duality-and-order.json` (phase and order keys)
- `input/pressure/plus-global-height-PV-certificate.json` (whole) and `minus-global-height-PV-certificate.json` (first lines)
- `input/pressure/global-profile-envelope.json` (first lines)
- `input/reference/native-first-shape-input.json`, `plus-native-height-constant.json`, `plus-new-native-trace.json`, `minus-new-native-trace.json` (head), `minus-final-native-slot-routing.json`
- `input/reference/retained-response-census.json` (flat, height, slope, mixed entries)
- `input/address-representatives.json` (first record) and `input/address-review-index.jsonl` (only matched, long lines omitted)
- `input/saved-inputs/weak/` field 79ea957c… (operands, reconstruction-input, derivative-class), field 6e77d90a… (reconstruction-input, polynomial), and `all-coefficient-certificates.json` (sample lines)

### Not inspected (relevant material)

- `input/tests.py`, `input/execution-authority.json`, `input/saved-input-map.json`, `input/review-prompt.md` and `input/control-candidate-metadata.json`.
- `runtime-source/shared-guard.py` and `supervisor.py`. I only confirmed the worker and launcher pin them by hash and path. I could not establish that they are unchanged from any prior version.
- Most of the 422 saved aliases: the 13,260 address records beyond the sampled ones, the remaining 33 field files, the other local cells, and the slope and mixed census details.
- `native/c1.py`, `native/c2.py`, the old source workers, the non-first-field operands, and most address and control metadata.
- Any simplifier-dependent behavior of SymPy on the actual payloads. Nothing was executed.

**CLEAR FOR THIS TRANSLATED WEAK-END BUILD**