**Verdict: NEEDS REVISION.** One deterministic defect would burn the single authorized science run. The rest of the worker matches the bounded method on the points I checked.

## Blocking defect

**The native unit check (`dimension()` at `worker.py:498-507`, called at `worker.py:508`) cannot pass on the real native rows.**
- `dimension()` requires every `Symbol` name to be a key of `contracts['dimensions']` (`worker.py:501`).
- That table is checked at `worker.py:478` to equal the static `DIMENSION_SCHEMA` literal in `source/c2.py`. I found no `gamma_*` key in `saved-inputs/native/wave-profile.json`, and the only `gamma_` hit in `c2.py` is the code string at line 229.
- In `c2.py:228-237` the `gamma_` dimensions are derived at runtime from constraints and written into `DIMENSION_SCHEMA`. The static literal does not hold them.
- All five native rows contain `Symbol('gamma_s11cb_…')` atoms (for example 439 occurrences in `U0.json`, 313 in `THETA_BALANCE.json`, 865 in `E_W_BALANCE.json`).
- `rowunits=[dimension(native[r]) for r in ROWS]` therefore hits `require(expr.name in dimensions, 'native dimensions gamma_…')`. That `ValueError` is caught by `SourceStage('native-source-chart-units')`, which raises `SOURCE_MAP_UNRESOLVED`.

Why this blocks:
- The stop is fail-closed, so nothing false is accepted. It is also deterministic, and it comes after the heavy cyclic covariance pass.
- `scientificRunsAuthorized==1` and there is no automatic retry. The one run is spent before the end-cell, restoration, A/R, grade, grazing and control stages execute.
- `tests.py` has no `gamma` reference, and its dimension test only compares the contract to the static literal (`tests.py:67`). The 37 stdlib checks cannot see this.

Required fix: join the gamma dimensions from an actual saved source record. Either carry the derived `gamma_` dimensions in the saved dimension contract, or reproduce `c2`'s constraint derivation as an inert, checked join. Do not guess them or exempt `gamma_*`. Then add a test that runs `dimension()` over a native-row-shaped atom set containing `gamma_*`.

## What I checked that looks faithful

- **Chart and side joins (item 1):**
  - `S` is the proper cyclic permutation and `S^T E S` gives old row `i` equal to weak row `roworder[i]=(1,2,0,3,4)`. `native_lift=S*lift` (`worker.py:605`) is consistent with that.
  - `rotated_name` and `rowmap` agree (U0→U2, U1→U0, U2→U1). The classifier fails closed on unrecognised names.
  - Native rows carry `d1`, `d2` and `d3` profile jets, so the covariance test is well posed.
  - The side map comes from actual `Limit` orientation, tanh endpoint values and half-heights `W*end/2`, with no residual fitting.
- **Normalization (item 2):**
  - All 400 local cells, 200 end cells and children 122, 113 and 114 are joined to the saved operands.
  - `diff(…, eps)` is applied per unit carrier. The 2π bookkeeping is explicit and no power map is used.
- **Restoration (item 3):**
  - The 18-operation census needs the two native-reconstruction entries that `inputs.json` supplies. The worker's operation names and `sign--1` / `sign-1` blob names match the receipts and `inputs.json`.
- **Attribution (item 4):** the recursion, remainder H and A = Σ J_ab − H_phys are correct.
  - `R−A=I_old` is an exact identity given `Ephys−Pold`.
  - `attribution_usable` has the specified fatal and unavailable split.
- **Wave tests (item 5):** the q-reduction, quotient check and witness point are sound.
  - The point p=1, q=2, cs²=180/101 satisfies the wave relation.
  - `finite_constant` is conservative.
- **Controls (item 6):**
  - Each control records baseline, mutation and the changed-minus-baseline movement, and the curl-entry fallback leaves an explicit gap record.
  - Per-cell `selectedWeights` and `zeroWeight` are saved, so untested pressure/direct entries stay visible.
- **Grazing (item 5):** R0 uses the closed E values, the saved selected limits and both inherited paths. No old function is called.
- **Execution (item 7):** the posthash and encoded-summary design looks sound.

## Nonblocking notes

- `mapping` at `worker.py:566` omits `kl` and `ql`. If leg 0 contained them, `worker.py:573` would stop with `SOURCE_MAP_UNRESOLVED`. That is fail-closed, but worth confirming against the leg-0 free symbols.
- The omit-cell control removes the whole `(grade, entry)` symbol rather than one term. The saved ancestry is retained, so I read this as the "full-cell/grade omission" in item 6.
- The opposite-sheet `applicable` test at `worker.py:667` is column-level, not row-level. It is sound only because q-dependence enters through pressure sums.
- `sp.simplify` in `finite_constant` and `sp.cancel` on large expressions are risks for time and memory under the 4 GiB guard. Memory fit is unverified, as stated.

## Coverage limits

I read `method.md`, `build-guide.md` and all of `worker.py`. I sampled `inputs.json`, `opaque-interface-receipts.json`, `physical-input.json`, `c2.py` and the native row and wave-profile JSON by grep. I did not read `launcher.py`, `tests.py` in full, `engine.py`, the runtime helpers or the other source files. I did not decode any blob, so the binding of `Praw`, the lift and the old operands is unverified. I read the gamma defect from the static source and did not execute anything.

NEEDS REVISION