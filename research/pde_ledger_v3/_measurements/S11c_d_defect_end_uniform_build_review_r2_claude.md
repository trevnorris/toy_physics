**NEEDS REVISION**

I found three defects in worker.py that I would fix before this build runs. All three are narrow.

## Blocking defects

**1. A failed algebraic reconstruction is recorded as "attribution unavailable" instead of stopping.** (`worker.py:500-503`)
- The check is `wave_test(A[i,j] - (ΣηᵃσᵇJ_ab - H)_phys)`. After the exact join of P_raw to P_old (`:482-483`) and the exact A/R consistency (`:489-490`), this is an identity modulo the wave relation.
- The code treats any status other than `ZERO` the same way. It logs `ATTRIBUTION_UNAVAILABLE` and continues.
- That also catches `NONZERO_CERTIFIED` and `UNRESOLVED_NONZERO_SYMBOLIC_REMAINDER`. Both mean the reconstruction `P_raw L − retained − H` has failed, not that coefficient dependence is unsupported.
- The method (point 4) says failed algebraic reconstruction still stops. Only `UNRESOLVED_UNSUPPORTED_DEPTH_DEPENDENCE` and `newc['remainder']!=0` should be recorded statuses.
- Fix: raise on any attribution status other than `ZERO` or the unsupported-dependence status.

**2. The epsilon/Fourier "formal extraction identity" does not touch the saved symbols.** (`:129-131`, `:139`)
- `extraction` and `division` are built from freshly created symbols `normalization_trial` and `normalization_coefficient`. The identity `∂ε(ε·c·k)/c = ε·c·k/(ε·c)` is true for any inputs, so it cannot fail.
- The record is nonetheless labelled `inheritedNormalizationProofs: True`.
- The only checks that touch the 200 saved cells are `epsilonPower`, `weakPairingFactor` and "no `epsilon_shape`" (`:440-441`). Those are the saved cells' own attestations.
- The `source_fragment` checks confirm only that statements exist in the source text. Nothing joins a saved symbol to `raw_coefficient = original/(eps*wave)` or to the strong-matrix `diff(·,ε)` coefficient.
- A wrong scalar such as `i` or the wave factor would therefore surface as a nonzero A/R and be reported as `RETAINED_MISMATCH`. The task says the normalization must succeed at runtime before a residual is interpreted.
- Fix: join at least one actual saved local cell. Compare `val*(i p)^xOrder` against the stored coefficient of the native child divided by `eps*wave`, using the saved ancestry. Otherwise drop the `inheritedNormalizationProofs` label and declare the gap.

**3. Covariance and convention failures are mostly not reported as `SOURCE_MAP_UNRESOLVED`.**
- Only `source_chart_and_scale` (`:365-368`) and the row-covariance `J.zero` (`:399-400`) are wrapped.
- These failures raise plain `ValueError` and end as `FAILED_PRESERVED`:
  - an unclassified native atom (`:172`)
  - a missing mapped atom (`:393`)
  - a native source-row mismatch (`:380`)
  - a unit-schema or wave-source mismatch (`:385-387`)
  - an unsupported dimension node (`:414`)
- The method and the launcher completion message both specify `SOURCE_MAP_UNRESOLVED` for these.
- Fix: wrap `:371-415` in the same `try`/`except` as the covariance step.
- The run would still stop either way, but this mislabels the status.

## Non-blocking notes

- **Execution risk:** the covariance `J.zero` on 895 KB rows calls `cancel(together(a-b))`. If the two sides are not structurally identical it could exhaust 4 GiB, and the process would be killed with no classification. The input and return receipts are written first, so evidence survives.
- **Grazing R0:** the saved path pencils and selected join are emitted but never compared against the direct R0. That fits the method's "no A at singular P_raw", but the "join closed E to the saved limit inputs" language is satisfied only by argument-level joins (`:516-521`, `:533-538`).
- **Controls:**
  - The sheet-control applicability test (`:571`) asks only whether any pressure cell exists in the column, not whether it is the cell carrying `q`.
  - The curl-entry fallback is recorded with `directionProof:False`, so the silent normal-sign result is not lost.
  - The omit-cell control removes whole cell entries, not single terms.
- **Zero-weight cells:** `zeroWeight` and `lift[j,c]!=0` are structural tests. A symbolic zero not in canonical form would be missed.
- **Mapping checks that look correct:** these check out on direct reading.
  - The cyclic map: `S[:3,:3]·(a,b,p)=(p,a,b)`, `S^T E S` with `roworder=(1,2,0,3,4)`, and `rowmap` all agree with weak(1,2,3)=old(3,1,2).
  - The side map from tanh ends to half-heights 0 and W/2.
  - The bijection check on sides.
  - The `p→kr, q→qr` mapping for leg 0.
  - The `A=(E−P)L`, `R=EL−LD`, `R−A=I` algebra.
  - The `H` bookkeeping and the coefficient recursion.
  - The q-reduction division.
  - The 4 GiB containment, the no-CPU-deadline check, and posthashes over sources, copies and opaque blobs.
  - The `SavedCodec` and `decode` extraction, whose namespace dependencies are all supplied.

## What I read, and the limits

- **Read in full:** worker.py, method.md, build-guide.md, launcher.py, runtime-source/inert-helpers.py.
- **Read in part:** native-wave-profile-contract.json (lines 1–3036 of 4318), `SavedCodec` and `decode` in uniform-worker.py, the `fieldUnits` usage in uniform-worker.py and frequency-source.py, and the `w1_profile` atoms in saved-inputs/native/U0.json.
- **Not read:** tests.py, inputs.json, engine.py, the ends and full-weak workers, the controls evidence, and the rest of the 14 JSON inputs.
- **Consequences of what I did not read:**
  - I could not confirm that the 20 AST statement contracts match the pinned sources.
  - I could not confirm that exact-structure joins such as `:77-78` compare like types (string against string, or sympy against sympy).
  - I could not confirm that the covariance check will pass on the real rows.
- No payload was decoded and nothing was executed.

**NEEDS REVISION**