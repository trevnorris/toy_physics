**NEEDS REVISION**

I read the sources and operands only. Nothing was executed, so this is source-level review, not runtime proof.

## Blocker

**1. Wrong symbol name in the middle-sheet map, which makes the run fail.**
- **Location:** `worker.py:315`, `rootmap={...'H':m-k...}`.
- **Cause:** the saved `physical-sheet.json["heightRoute"]` has no symbol `H`. Its symbols are `increment_transfer`, `k` and `increment_effective_bulk_speed`. `H` appears only in the lower-boundary modes. The `.get(x.name,x)` lookup therefore leaves `increment_transfer` unmapped.
- **Consequence:** `native_middle` still contains `(increment_transfer+k)`. `J.zero('middle-propagating-sheet')` at `worker.py:322` then compares it with `sqrt(realrad)`, which is written in `m`, and the residual is non-zero. The single authorised run dies there. Everything after it never executes:
  - the sheet-domain checks,
  - the reused-domain and new-bounds evidence,
  - collisions,
  - dimensions,
  - all four controls, including `control-native-middle-root`.
- **Why the tests missed it:** `tooling-tests.py` never exercises `run_science`.
- **Minimum correction:** key the map on `'increment_transfer'` instead of `'H'`. Optionally add `require(oldroot.atoms(Symbol) ⊆ mapped names)` so any future unmapped saved symbol is caught.

## Required but not mathematical

**2. `height_constant` is hard-coded.** `worker.py:210` sets `'height_constant':sp.S.Zero` rather than extracting it from the saved height.
- The saved heights (`eta_bg*w1_profile/2` and its negative) are purely linear, so zero is the correct value.
- The native `trace_two` and `trace_three` assignments consume it, so it should be derived.
- Add `J.zero(label+'-height-constant', labheight.subs(profile_symbol,0), 0)` and pass that value.

**3. Gate join omits the launcher.** `verify_gate` (`worker.py:49–69`) ties the build record only to the worker, manifest, shared guard and supervisor hashes.
- `launch.py` is pinned only by the gate's own `launcherSha256`, so a launcher edited after this review would still pass.
- Add `launcherSha256` to the build-record join and to the gate.
- Pinning the tooling-test hash is optional.

## Should fix or disclose

**4. The wrong-lower-height control is not derived from the lower trace.**
- `wrongLower=(phs-tracehs)/…` (`worker.py:389`) is a hand-coded sign flip.
- `control-actual-native-trace` emits `corruptT01=-lowerT`, but nothing consumes it.
- Build the corrupt coefficient from `-lowerT` using the same `nativeTraceCoefficient` formula, so the control acts through the actual lower trace.
- At the declared point I computed by hand that all four control movements are non-zero, so a spurious zero is not a concern.

**5. Several inequalities are text only, not machine-checked.**
- The profile-product tail constant `150·e^{30π}` and `1−e^{−10π}>1/2` appear only in the reused `profileProof` string.
- The lower bound `0 ≤ x·cosh x − sinh x` is not checked.
- `domain['qMagnitudeBound']==5` is never joined. It is needed for `|μ qi|≤2`.
- `κ_δ ≥ κ_min` is checked only as a constant equality (`kappaMinimum == sqrt(879)/20`), not as an inequality over cs and δ.
- The arithmetic I re-derived by hand is correct: the `200` tail constant, the `u≥12` polynomials, the Hölder quadrant gap, the `4√7/(5 b0²)` jet variation and the `2/b0`, `2√7/b0²` height constants.
- Add `require(domain['qMagnitudeBound']==5)` and a `kappa-lower-gap` identity of the kind `closed-worker.py:237` already uses. Otherwise label these items as assessed text in the result.

**6. Minor gaps.**
- The restored direct map gives `qr` (`grazing_qs`) no sheet join. This is acceptable for a tag, but it should be stated.
- `ref[0,1]` is never compared with (5)/(6) at `k→m`.
- The final-slot check is partly tautological: `ref[0,col]` is defined as `P−T01·ref[1,col]`. Its independent content is the reference-equation residual and the closed-form (5)/(6) checks.

## Checked and found faithful

- **First shape:** the on-shell, edge and profile mapping reproduces `iμ(qo−qi)/qo·ηh + μk/(qoqi)·σj`. The lower modes are selected by grade, not replayed.
- **Closure:** `(I+aZ)F=Z` holds for unrestricted ω. I checked the saved ω=3 closed entries by hand against the worker's `F`, including `F02`.
- **Trace and mixed algebra:**
  - T01 and T12 equal `iη q h`, with face² = 1 and the native `T02=0`.
  - `phs+tracehs = C` for all `qm`.
  - The right-height contact vanishes, and `J` follows from the factorised density.
  - The odd PV piece is removed.
- **Direct term:** the restored direct density's symbol map covers every saved symbol. The tag is carried once and never integrated again.
- **Dimensions:** the keys and arithmetic reproduce `(-2,-1,1)`.
- **Evidence ordering:** operands are emitted before each guard. Failure is preserved, with posthashes and no retry.
- **Containment and launch:** containment runs before the sympy import, and there is no deadline. The hook-first launcher with one run is as described.

## Evidence if the worker passes (after fix 1)

For both faces, per unit native source, with eta and sigma independent, on cs∈[1,2], k/l∈[−3,3], and ω+i0 at the saved edges, it would establish:
- exact source and operand joins for zero, height, slope and iterated-mixed coefficients;
- the reference-trace conversion with the surviving contact;
- explicit elementary inequality certificates supporting the assessed L1, Hölder and tail arguments;
- four responsive controls.

It does not establish measure-theoretic limits (machine-proved), integral values, full source/consumer composition, a finite inverse, the 4/6 cutoffs, eta², sigma², κ=0, β=0, any defect or loss claim, or runtime restoration, memory or timing success. Source clearance is not result acceptance.