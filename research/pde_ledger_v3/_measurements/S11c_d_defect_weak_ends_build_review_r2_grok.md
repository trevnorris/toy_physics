I'll start with the build guide, then assess the worker, launcher, and inputs against the method and the saved operands, staying inside `input/` and reading only.The worker’s science path starts after the gate and phase join. I’ll read that implementation next, then check it against the saved schemas the guide names.The science routine is in one place. Next I’ll check the helper, the saved field and phase schemas, and the parts of the manifest that pin runtime and authority.The trace signs and the minus-face Holder coefficient match the method. I’m checking whether factor-cache keys, polynomial storage, and the first-shape substitution actually preserve those identities.The report still has to be written from the implementation, so I am checking the constant-height specialization and the saved first-shape operand before stating the verdict.The constant-height check in `input/worker.py` sets the independent symbol `sigma_W` to 0 and then requires the whole saved first-shape right-hand side to cancel to 0. That is a false specialization of the saved operand. The build needs revision.

## Blocker

`input/worker.py` lines 437–460 load `reference/native-first-shape-input.json` and use only its `right` side. The bindings are

```437:441:input/worker.py
    first=load('reference/native-first-shape-input.json')['right']
    first_bindings={'reference_unrestricted_frequency':sp.Integer(3),'reference_qi':q,'reference_qo':q,
                    'reference_k':p,'reference_l':p,eta.name:eta,sigma.name:sp.S.Zero}
    first_diag=first.xreplace(named_map(first,first_bindings))
    # h and j functions carry zero transfer now; j is multiplied by sigma=0.
```

`eta` and `sigma` come from the top-level pair in `input/saved-inputs/full/extended-binding-context.json` lines 597–606: `Symbol('eta_bg', real=True)` and `Symbol('sigma_W', real=True)`. The map therefore sends `sigma_W` to `S.Zero`. For each face and each endpoint the same substituted expression is passed to `zero('constant-height-physical-first-shape-'+face+'-'+side, dtn, 0)` (lines 459–460). `zero` (lines 266–269) records both sides, then requires `cancel(together(left-right)) == 0`.

The saved right-hand side (`input/saved-inputs/reference/native-first-shape-input.json` lines 4–8) is two addends:

- height: `I*eta_bg*ω*(-qi+qo)*height_hat(-k+l)/(10*qo)`
- slope: `k*ω*sigma_W*slope_hat(-k+l)/(10*qi*qo)`

On common `qi = qo = q` and `k = l = p`, the height factor `Add(-qi, qo)` becomes `Add(-q, q)`. That height addend is the object `method.md` lines 177–178 call the first-height DtN correction, and it vanishes on that diagonal without touching `sigma_W`. The slope addend has no `(l-k)` factor. After `k = l` it is still `p*ω*sigma_W*slope_hat(0)/(10*q*q)`. The code replaces `reference_height_hat` only (line 459). `reference_slope_hat` stays. The four guards pass because `sigma_W` was replaced by 0, which is what the comment on line 441 states.

`method.md` line 13 keeps `eta` and `sigma` independent in `G = {00, 10, 01, 11}`. The vanishing claim there is the height correction `Zh = i*mu*(q(l)-q(k))/q(l)*h`, not the whole first shape. The retained census slope coefficient (`input/saved-inputs/saved/reference/retained-response-census.json` lines 38–40) has the same nonvanishing shape: `reference_k*ω*slope_hat(-k+l)` over the two depth denominators, with `slope_hat` still present at zero transfer.

This is a changed scientific validation inside the new constant-height join. It is separate from the trace, affine slot, and end-symbol sums, which do not use this substituted `dtn`. Those later formulas can still be assembled while this guard certifies the first shape only on the `sigma_W = 0` slice. No SymPy process was run. The reading is the saved `srepr` plus the substitutions the worker applies before `cancel`/`together`. Whether `reference_slope_hat(0)` equals `(L/2)A(0)` was not found as a saved identity in the files opened; the slope function is still present, and the guard deletes it by zeroing `sigma_W`.

The earlier worry that `context['independentGrades']` and `context['epsilon']` exist only inside `context['saved']` is not a defect. The top-level pair is at lines 597–610, after top-level `numeric` closes at line 567, and top-level `numeric.omega` is `Integer(3)` (lines 415–418). The worker’s join at line 290 binds those top-level symbols.

## What matches the method

These checks agree with `method.md` and the saved operands that were opened. They do not clear the first-shape guard.

Endpoint storage. Pressure polynomial files store a numerator and a separate denominator. Certificates and `P0` store the quotient. The `f588…` reconstruction right-hand side is that quotient with `weak_tanh_variable` replaced by `tanh(composition_x/10)`. `inherit_zero` accepts only a cancelled `Integer(0)` and does not demand that the two sides of an old identity be structurally identical. Local cells use `full_weak_tanh_variable`. The worker takes saved `leftEndpoint` and `rightEndpoint`. The first local cell’s endpoints are both `Integer(0)`, and that cell’s degree-4 polynomial evaluates to 0 at `T = ±1` by direct rational arithmetic. One cell of 400 was checked that way.

Phase and height limits. The restricted interpreter covers the saved source and profile exponential strings. With both edges equal and a shift of only the first coordinate, the increment is `i(l-k)a`, the profile exponent is the opposite sign, and the translated diagonal exponent is `iQa`. Flipping the bound `I` moves that coefficient. The reversed exponent is that derived negation. The PV split keeps the contact `W/4` and the oscillatory diagonal term. With `A(0) = 1/(2π)` and `W = 1`, the Dirichlet limits are `0` and `W/2`. The uniform-bound text uses the sine-integral bound 4 and does not differentiate `exp(iQa)`.

Holder, trace, affine, and grazing. The plus and minus normal certificates are `± μ q²/(q+β)`, with `normalCoefficientBoundedDirectly` true and `qTimesTestAssumedC1` false in the saved certificates. Trace height and normal both flip sign, and their product is `i q η H_side` with `H_- = 0` and `H_+ = W/2`. The final slot used is `equationSolution`. The retained inverse is `1 - product`, and the excluded square is checked by `(1+product)(1-product) = 1 - product²`. The closed grazing value is the continuous value of `(1 - i q η H) F0` at `q = 0`, which is `μ/β`. The raw first shape is not evaluated at `q = 0`.

Addressed controls. Height and reversed-phase eligibility is response grade `(1,0)` and `NATIVE_HEIGHT`. Address `8034` is that route: constant nonzero source and consumer, pressure normal `1`, height response that becomes `H B_h` after diagonal substitution. Lower-normal eligibility is face `minus`, slot `normal`, response grade `(0,0)`, `NATIVE_FLAT`. Address `10062` is that route: flat response `F0`, normal factor `-I*common_outgoing_q(composition_l)`, consumer grade `(1,0)` kept. The consumer `(12-40I)*(I*tanh(composition_x/10)+I)/1744` is nonzero at the plus endpoint `T = +1` and zero at `T = -1`. The plus-end term is the one the worker requires to be symbolically nonzero, so this route is not dropped. The reverse-phase mutation inserts the exchanged derived limits into both end cells. Response grades `(0,1)` and `(1,1)` are explicit weak-limit zeros, with source and consumer cross grades retained on the flat and height responses. No old pencil or mode comparison is an acceptance test.

Runtime pins. `inputs.json` is `PREPARED_NO_SCIENCE_OR_GATE`. The worker demands a later gate with literal clearances, one run, and `durationLimits` null. `containment()` is invoked before `import sympy`. The launcher arms the hook before the scientific command, sets `automaticRetry` false, and passes the pooled guard at 4 GiB, tasks 32, inside the 16 GiB pool. Packet `shared-guard.py` records `--seconds` and does not impose a science deadline (`RuntimeMaxSec=infinity`). Packet `supervisor.py` waits on the child with no wall-clock cutoff. The 30-second `select` and the 5-second arming loop bound hook startup. The guard’s 2-second `wait` is a memory sample, and exit 124 is the host-reserve stop if `MemAvailable` stays under 4 GiB. Absence of the gate file is the declared pre-assessment state.

## Inspected files and regions

- `input/build-guide.md`, whole
- `input/method.md`, whole
- `input/worker.py`, whole, including `verify_gate`, `inherit_zero`, `native_phase_exponent`, `control_eligible`, `zero`, the field-endpoint loop, and lines 338–460 re-read for this report
- `input/launcher.py`, whole
- `input/inputs.json`, header and the resource tail
- `input/saved-input-map.json`, opening aliases only
- `input/runtime-source/inert-helpers.py`, through the start of `scientific_work` (about line 400); `containment` and `Journal`/`decode`
- `input/runtime-source/supervisor.py`, whole
- `input/runtime-source/shared-guard.py`, whole (441 lines)
- `input/saved-inputs/reference/native-first-shape-input.json`, whole
- `input/saved-inputs/full/extended-binding-context.json`, saved object, top-level numeric including `omega` and the edge momenta, and lines 567–612
- `input/saved-inputs/saved/reference/retained-response-census.json`, lines 1–80
- `input/saved-inputs/weak/analytic-conclusion.json`; `new-weak-duality-and-order.json`; `global-profile-envelope.json`; `all-coefficient-certificates.json` lines 1–109; the `f588…` polynomial, operands, derivative class, reconstruction input, and reconstruction return; plus and minus height PV certificates
- `input/saved-inputs/reference/` plus and minus new-native trace, final slot routing, and height constant; restored profile arguments, lines 1–69
- `input/saved-inputs/full/complete-weak-assembly.json`, lines 1–40; pressure factor arguments, opening and the factor-8 window; wave arguments, lines 1–60; `all-local-cells.json`, first cell through its endpoints
- `input/saved-inputs/inventory/fields.json`, lines 1–120
- `input/control-candidate-metadata.json`, through line 1011, including address `8034` and address `10062`
- `input/address-representatives.json`, lines 1–150
- `input/native/c2.py`, header and the window near line 400

## Uninspected material that limits the verdict

- The other 16 of the claimed 17 full-factor proof bodies, beyond factor 0, factor 1, factor 8, and factor 12. The response cache key omits response grade, so a shared proof id across grades was not ruled out for every id.
- Transverse and edge-derivative wave arguments after the `e_W` sample.
- The other 399 local cells’ endpoint values. Cardinality is what the worker asserts against the route grid; it was not recounted by hand.
- The other 33 field reconstruction files beyond the `f588…` sample and the certificate head.
- The `c2.py` assignment bodies for the phase strings and for `trace_two` / `trace_three`. The worker compares those strings by AST at runtime. The window near line 400 is an old second-scattering path that this worker does not call.
- `native/c1.py`, which the worker does not read.
- Byte identity of packet `shared-guard.py` and `supervisor.py` with the repository scripts named in the gate. Those repository paths were not opened.
- `address-review-index.jsonl`, `weak-address-coverage.json`, unrestricted-closure operands beyond the worker’s `inherit_zero` call, artifact-index entries beyond the receipt join, preparation tests, the middle of `inputs.json`, and the rest of `saved-input-map.json`.
- `inert-helpers.py` after the start of `scientific_work`. That body is not in the worker’s `HELPERS` exec list.
- No scientific payload, SymPy session, shell, hash recomputation, or peer report was used. The slope-hat value at zero transfer was not identified with `A(0)` from a saved equality.

## Nonblocking limitations

- The control cell identity compares a cell with the same term object that was just appended, so it does not independently recompute the movement.
- The worker requires `normalCoefficientBoundedDirectly` and does not also require the saved `qTimesTestAssumedC1` flag to be false. Both opened certificates have that flag false.
- `commonEdgesCancelInProfileAndTranslation` is stored as a literal true beside the actual edge-cancellation residual.
- `epsilon_shape` is decoded and then unused. The density tuple’s historical relation of `sigma_W` to `eta_bg` is not applied on the trace path, which matches independent grades aside from the first-shape substitution above.
- The declared control point uses `cs² ∈ (1,4)`. The value `180/101` also lies in the method interval `[1,2]`.
- The guard can still stop the job with exit 124 when host available memory stays under 4 GiB. That is the host reserve, and the scientific `subprocess` call has no timeout.
- There is no scientific job and no READY gate in this packet. That absence is the stated preparation status. It is not, by itself, a defect in this build, and a future gate is still required before any scientific run.

NEEDS REVISION