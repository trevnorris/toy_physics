**NEEDS REVISION**

There is one concrete blocker, in the J/D missing-`dt` controls. The dimensional interpretation itself checks out against the sources I read.

## Blocker

**The two `control-missing-dt-{J,reflected}` controls cannot fail** (`input/worker.py:294-297`).
- They compute `baseline=add(unit0,MOMENTUM)` and `mutant=add(unit0,ZERO)`, then set `responded=baseline!=mutant`. That is true for any density unit, so `require(r['responded'])` can never refuse.
- Nothing in the actual route is touched. The control never goes through the NATIVE_MIXED_ITERATION or INHERITED_DIRECT template AST, `K.assemble_route`, or `K.pressure_dimension`. It also does not use the `Jvalue`/`Dvalue` terms in `tenv` (`worker.py:233`).
- `method.md:205-209` requires that mutations "affect an applicable original expression/route, not only a displayed label; record a refusal/movement and its scope." `build.md:93-95` rejects a bare unit-tuple difference for the flat control, and the same standard should apply here.
- Fix:
  - Re-walk the actual mixed-iteration template with `Jvalue=densities[0]` (no `+MOMENTUM`). The `Add` should refuse, because `J` is `(-1,-1,1)` and the H term is `(-2,-1,1)`.
  - Assemble a live direct/D address with `Dvalue` lacking `dt`, and require that it fails `pressure_dimension`.
  - Save the walk and assembly returns before requiring the refusal.

## What checked out

I did these by hand from the saved registry and source text. I did not run anything.
- **Parameters:**
  - mu is `(-4,-1,1)` and a is `(3,1,-1)`.
  - beta is `(-1,0,0)`, so `q+beta` is homogeneous.
  - The preflight source at `S11c_d_defect_packet_preflight.py:235` matches the pinned fragments, with the bulk `rho_m` at `(-4,0,1)`.
- **Kernel:** `kernel_components` gives J and the reflected, height and quadratic terms each at `(-1,-1,1)`, so each is `(-2,-1,1)` after `dt`. The roots stay separate.
- **Call sites:** I read both of them (`inner_lib.py:121` and `inner.py:139`). Their positional bindings are correct, and the two-argument actual `A` is `5z/(2 sinh(5πz))`, i.e. L=10.
- **Profile and H:** hhat `L²`, jhat `L`, H contact `L²`, H density `L³`, and the W/4 coefficient `L`.
- **First shape:** the native first-shape expression gives `(0,-1,1)`, and `(-2,-1,1)` after the two edge deltas.
- **Template components:** the five original template components come out at `(-3,-1,1)` for flat and `(-2,-1,1)` for the others, with the normal slot adding one `L^-1`.
- **Routes:**
  - All 544 saved addresses agree with `address_route`. There are 288 flat addresses (`dl`, `k=l`, q(l)→q(l)) and 256 off-diagonal ones (`dl dk`, no flat support, q(k)→q(l)). Of the 544, 442 are exact zeros and 102 are formal.
  - Flat address 7956 totals `(-2,-1,1)` by hand.
- **Flat control:**
  - `library.py:170-186` runs the reduced baseline, the delta-in-`dl dk` supported route, and the delta-removed mutant through the same `assemble_route` and `pressure_dimension`.
  - Both returns are saved by the caller before the decision.
  - The mutant totals `(-3,-1,1)`, so it fails, and the supported route equals the baseline.
- **Normal map:** the nested `completeNormalMap` is preserved from `fullFactorProof`.

## Optional advice (non-blocking)

1. The `a`-unitless control (`worker.py:298-299`) compares only `mutated[0]!=densities[0]`. Propagating it through the J template and `pressure_dimension`, as above, would give it real gating power.
2. The worker never checks a call site's positional argument binding against the `kernel_components` parameters (`worker.py:227-229`). It only persists the call text. The sources are hash-pinned and I verified the binding by hand, but an explicit AST check would be cheap.
3. The worker does not compare `a['normalMultiplier']` with the face sign and q(l) (for example address 8619's `I*common_outgoing_q(composition_l)`). It relies on identity with the saved complete factor (`fac['mappedAddressFactor']==...`, `worker.py:262`). The same applies to `normalMap`, which is saved but not compared.
4. The X-family match at `worker.py:278` ignores `argumentDerivative` and the family's `addresses` list. The Y match checks the derivative only.
5. The height contact coefficient is taken from the mixed H record's `W/4` (`worker.py:212-214`). The only join is the arithmetic `heightcontact+L==hpv`, not a height-only source record.
6. The identification of the literal `5` in the callers' `A` with `L_W/2` is an accepted-profile-scale dependency, as are the gamma registry units (`Lambda_A_0=(-5,1,1)`). `q(l)` versus another momentum cannot be separated by dimensions alone, which `build.md` already states.

## Coverage

- **Read in full:** `build.md`, `evidence-guide.md`, `worker.py`, `library.py`, `source-contracts.json`, `method.md` sections 2 and 4, and the helper units libraries' relevant parts (`runtime-source/S11c_d_defect_packet_native_units_lib.py` lines 70-123, plus the imports of `S11c_d_defect_packet_source_units_lib.py`).
- **Read in part:**
  - `inputs.json` (lines 1-70 and 12670-end, plus targeted greps).
  - `pressure-addresses.json` (the first address and address 8619, plus counts across all 544).
  - `effective-registry.json` (lines 1-3246 of 4273; every name the worker uses is in that range).
  - `first-shape-native-transport.json` (first line).
  - `family-0.json` (first 40 lines).
  - `original-source/S11c_d_defect_packet_inner_lib.py` (lines 95-139) and `original-source/S11c_d_defect_packet_inner.py` (lines 118-152).
  - `original-source/S11c_d_defect_packet_preflight.py` (one targeted grep).
- **Not read:** `launcher.py`, `execution-authority.json`, the tooling tests and logs, the per-address and factor files, `fourier-lib` beyond the pinned product fragment, and the reference and raw-direct source bodies. I did not run any tool or script.

This is a source-level verdict only. It is not execution or result acceptance.