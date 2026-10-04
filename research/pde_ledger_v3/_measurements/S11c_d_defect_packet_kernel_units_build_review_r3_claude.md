**CLEAR FOR THIS BOUNDED KERNEL WAVE/MEASURE UNIT-TRANSPORT BUILD**

This is a source verdict only. It is not execution or result acceptance, and it clears nothing outside the two rest-bulk Gaussian pairings.

## What I read
- **Read in full:**
  - `build.md`, `evidence-guide.md`, `worker.py`, `library.py`, `source-contracts.json`
  - `runtime-source/…source_units_lib.py`
  - the original preflight lines 180–300, `inner_lib.py` 40–150, and `inner.py` 100–170
- **Read in part:**
  - `inputs.json` (759 KB): lines 1–62 and 12665–12723, plus alias-key greps.
  - `native_units_lib.py`: only `require`, `unit_tuple`, `add`, `scale` and the `Units` header.
  - The saved JSON files: `selected/pressure-addresses.json` (first two addresses and header), `native/effective-registry.json` (the six parameter units), `native/lower-geometry.json`, and `native/fourier-contract.json`. For the Fourier files I only grepped `fourier/transform-convention.json` for the two convention flags.
  - Saved address records 7956, 10023 and 8346, and the key fields of the others by grep.
- **Not read:** the other ~540 address records, the factor proofs, the 19 Fourier families, `launcher.py`, `supervisor.py`, `completion-hook.py`, `tooling-tests.*`, `reference-grazing`, `raw_increment`, and `fourier_lib` beyond the pinned fragment. I did not hash-check the pins. I executed nothing.
- **How far the partial coverage reaches:** the data-dependent checks (address route, depth map, source/consumer units, `pressure_dimension`) are `require`s that stop the run. A bad saved record would stop it, not be accepted silently. That is why I rely on the sampled records.

## Hand-checked results
- **Parameters:** mu=(-4,-1,1), a=(3,1,-1) and beta=(-1,0,0) follow from preflight line 235 and the registry. q+beta is homogeneous, and the bulk density keeps its four-dimensional unit.
- **Profile and H pieces:**
  - A is dimensionless, hhat is L², jhat is L, Hsub is L², the H contact is L² and the H density is L³.
  - The W/4 coefficient is L, and L plus the delta's L equals hpv.
- **Kernel densities:** J, reflected, height and quadratic are each (-1,-1,1). Each gets (-2,-1,1) after dt.
- **Templates:** the five templates give (-3,-1,1) for flat and (-2,-1,1) for the rest.
- **Addresses:**
  - Totals close to (-2,-1,1) for flat pressure 7956, mixed pressure 8346 and normal-slot height 10023 (consumer units differ by slot).
  - Slot, depth and measure fields match `address_route`, including `deltaSupport` for flat.
- **Flat control:**
  - `address_route`/`assemble_route`/`flat_delta_control` give a baseline with measure `dl`, a supported route (`dl dk` plus delta, kernel +L) with the same total, and a mutant (`dl dk`, no delta) one L-power off.
  - The mutant fails `pressure_dimension`, the same predicate used by all 544 addresses. Both assemblies are saved before the decision, and the control uses a live flat address.
- **Missing-dt J:** the unintegrated J (-1,-1,1) added to the (-2,-1,1) first term triggers the exact "inhomogeneous source addition" refusal. The refusal-node event is required.
- **Missing-dt D:** `Dvalue` is a bare name, so it is homogeneous. It reaches `assemble_route` and fails the complete predicate because the total is off by L.
- **Unitless a:** the kernel AST gives J=(-4,-2,2), propagated to (-5,-2,2) through Jvalue. This refuses at the template addition.
- **Control integrity:** `template_route_attempt` accepts only that exact refusal. Unbound names and unsupported syntax stay fatal.
- **Live addresses:** live NATIVE_MIXED_ITERATION and INHERITED_DIRECT addresses exist, so the controls' `next()` selections will not raise.

## Concrete blockers
None.

## Optional advice (none changes the verdict)
1. **Kernel call sites are pinned but not interpreted.** `worker.py:227-229` only unparses `inner-runtime-call` and `inner-proof-call`. The kernel env is hand-keyed by formal name. By eye, both calls' positional order matches the signature. At `inner_lib.py:121`, `aa,bb` fill A1 and A2. At `inner.py:139`, `aa` fills `a`, so the same name means two different things. Binding the call arguments to formals by AST, with the result saved, would remove that dependence on reading.
2. **The normal insertion is not walked.** `worker.py:242` adds `MOMENTUM` by hand. `normalSign` and `normalDepth` at `worker.py:245` come from the face label and the string `'q(l)'`, not from `a['normalMultiplier']`. The sign is covered only by the inherited adapter proof (`preflight-normal`, `IfExp`). I'd assert against `normalMultiplier` and `responseOutputDepth`.
3. **The first-shape reduction is never joined.** `worker.py:173-175` only checks `firstunit==(0,-1,1)` and emits `reduced=(-2,-1,1)`. Nothing requires the reduced value to equal anything, contrary to the build.md claim that the reductions "connect" the dimension to the 1D kernel. `remainingProfileForwardPower:-1` at `worker.py:249` is a literal, not computed from -3+2.
4. **`raw-direct-scale` is not walked.** The `raw-direct-scale` and `preflight-normal` source fragments are verified to exist but never unit-walked. The profile units come from the reference `originalPV` and `Hcontact` walks plus inherited joins.
5. **omega is not pinned.** `omega=sp.Integer(3)` is not a pinned fragment, so omega gets its unit from the registry plus the frequency-3 binding check.
6. **Some controls are partly by construction.**
   - The delta's +L lives inside `assemble_route`, so the flat control tests the route rule, not independent physics (build.md says as much).
   - The D control replaces a unit on a single `Dvalue` node. It does not reach the three-addend sum, whose equal units are verified separately at `worker.py:224`.
7. **Saved-record presence is unconfirmed.** I did not confirm that every `source/*` alias exists in `result-record['records']` under `'complete/'`, or that every `fourier/*` alias is in its artifact index. A missing one would stop the run, not corrupt it.