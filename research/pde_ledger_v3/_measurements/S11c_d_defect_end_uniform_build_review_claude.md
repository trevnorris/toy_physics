## Build assessment: weak-end / saved-uniform comparison

**Verdict: NEEDS REVISION.** The arithmetic and containment are largely sound. The chart join is not executed, so a wrong join would be read as a physics mismatch. The scalar normalization is only declared. Several outcomes the method requires to be recorded as statuses would instead abort the single run.

### What I read
- **Read in full:** `method.md`, `build-guide.md`, `worker.py`, `launcher.py` and `runtime-source/inert-helpers.py`.
- **Read in part:**
  - `engine.py`: the leg and pencil construction around lines 3249–3407, plus grep hits for the leg symbols.
  - `phase-arguments.json`: first 60 lines.
  - `wave-profile.json`: lines 1–3036 of 4318.
  - `physical-input.json` and `inputs.json`: grep only.
  - `tests.py`: grep only.
  - Native row JSON (13 MB lines): grep only.
  - `uniform-worker.py` and `uniform-continuation.py`: grep for the field-unit lines only.
- **Not read:** the 200 cell JSONs, `shared-guard.py`, `supervisor.py`, the rest of `uniform-worker.py`, and the 58 opaque blobs. Nothing was decoded or executed.

### Blocking defects

1. **The source-chart join is not executed (`worker.py:282–283`, `:327`).**
   - **Covariance check:** `native[row].xreplace(rotation) == native[rowmap[row]]` only shows the generic native rows are covariant under the relabelling. That holds for any coordinate permutation or reflection, so it cannot tell `weak(1,2,3)=old(3,1,2)` from the other cyclic map or from reflections.
   - **What would discriminate:** the tangent slot order (`kout=(weak_end_l, 1/5, 1/10)` against the old edge order), the sign of the carrier phase and time, profile-normal orientation, and which end is minus/plus.
   - **Where those are checked:** nowhere. They appear only as literals in the emitted dict at `worker.py:282`.
     - `ends/phase-arguments.json` is copied and never read.
     - LEFT=minus and RIGHT=plus is a hard-coded tuple at `:327`.
     - `heights` is emitted at `:360` and never compared.
   - **Consequence:** a wrong join gets through to A/R and is reported as `RETAINED_MISMATCH` or a finite mismatch. The method requires a stop with `SOURCE_MAP_UNRESOLVED` before any residual is interpreted.
   - **Fix:** add executable joins that assert the tangent order, orientation, time sign and side to the old `EdgeReduction` / uniform-source data and to the phase-arguments record.

2. **The scalar normalization is declared, not checked (`worker.py:310`).**
   - `'factor':1` is a literal. The `strong_matrix` and `construct` sources are emitted but never compared against anything.
   - `weakPairingFactor==2*pi` at `:323` compares a cell metadata field with a constant.
   - Nothing independently shows that the cell symbols exclude 2π (and any (2π)^k) and that the old pencil is the unscaled leg-0 `strong_matrix`. The origin join covers only the old side.
   - A scalar discrepancy would therefore surface as a grade mismatch, confounded with physics.
   - **Fix:** add an executable check, run before interpretation, that does not depend on the residual itself.

3. **The attribution identity uses the wrong notion of zero (`worker.py:382` against `:363`).**
   - The P_raw–P_old origin join is accepted modulo the wave relation (`wave_test`, a/b remainders).
   - `A − attrib = (Praw(origin) − Pold)·L`, but `J.zero` demands exact `cancel==0` off the wave. If the two forms differ by any multiple of `q²−radical`, this aborts.
   - It should be tested with the same wave reduction as the origin join.

4. **Required UNRESOLVED/UNAVAILABLE outcomes are implemented as aborts.**
   - The method requires status outcomes for these cases:
     - an unsupported coefficient dependence (`polynomial` raises `ValueError` at `:153`, reached from `coefficients` at `:376`);
     - the attribution zero at `:382`;
     - `require(newc['remainder']==0)` at `:383`;
     - grazing `require(regular)` at `:420`, where the method says "unresolved";
     - the control-3 `require` at `:449`;
     - the `wave_test` requirements.
   - Because the run is single-shot with no retry and the LEFT end runs first, one such failure discards RIGHT, grazing, controls and the coverage table.
   - **Fix:** wrap each in a recorded status and continue.
   - There is no `require(len(operations)==18)` either. I counted 9 restore + 1 lift + 2 restriction + 2 native + 4 limit operations, which does give 18.

5. **Memory is unmeasured under the 4 GiB address-space limit (execution risk).**
   - The five native rows are 13 MB constructors. The worker decodes them with `sympify`, rotates them with `xreplace`, and emits all five as `text` plus `srepr` at `:282`.
   - Nothing in the packet shows this fits under `RLIMIT_AS`=4 GiB.
   - Prefer an AST/atom-level covariance check and do not emit the full rows twice.

### Checked and found sound
- **Chart algebra:** `S` is proper and orthonormal. `S.T*E*S` and the row map U0→U2, U1→U0, U2→U1 are consistent with old=(w2,w3,w1). `native_lift=S*lift` is the right direction.
- **Lift and A/R:** the lift field/unit augmentation matches the old continuation. The A/R/I_old identity (`:371`) is exact, and the on-wave zero is restricted to its stated scope.
- **Coefficient recursion:** the grade recursion and H reconstruction at `:156–172` are correct. The finite-truncation versus retained-mismatch classification follows the method.
- **Wave tests:** the witness point has radicand 4, so q=2, and cs≈1.335 is inside [1,2]. UNRESOLVED is returned conservatively.
- **Engine legs:** leg 0 is the right/positive leg, with `modes.k→kright`, `qlegs[1]→qright` and the tangents mapped. This matches the worker's binding of wr, kr and qr. `PLUS_ROW_POWER_MAP` is not applied to leg 0.
- **Controls and coverage:** controls record baseline, mutated and mutated-minus-baseline operands. Silent controls are recorded as gaps. The per-cell `zeroWeight` table identifies pressure entries the lift does not test.
- **Containment and gate:** 4 GiB cgroup and `RLIMIT_AS`, zero swap, 32 pids, one CPU, one thread, no CPU deadline, the guard-manifest environment variable, strict posthashes, an encoded final summary, and an inert-only helper import. The old 8 GiB module body is never imported.

### Nonblocking
- `sp.simplify` on large expressions in `finite_constant` could be slow.
- The omit-cell control removes a whole graded cell entry, not a single term. This is acceptable under item 6 but worth stating in the record.
- `tests.py` is metadata and stand-in only, as the guide says.

NEEDS REVISION