**NEEDS REVISION**

This is a source-only review. I executed nothing and recalculated nothing. No READY gate exists. The 56 stdlib tests (`tooling-tests.log`) are tooling, not native runtime proof.

## What I read
- **Read in full:** `method.md`, `build.md`, `guide.md`, `worker.py`, `prepare.py`, `resume.py`, `geometry-tail.py`, `launcher.py`, `tooling-tests.py`, `tooling-tests.log` and `original-prepare.py`.
- **Read in part:**
  - `input-manifest.json`: the scope, source-pin head, resources, method/record paths, `attestedGuardSources/Functions` and the `priorFiles` start.
  - `original-worker.py`: only a grep of the `run` and `copy_inputs` locations.
  - `prior/source/.../S11c_d_defect_packet_preflight.py`: the grep context around `weighted_tail`, `exponential_moment` and `contributions`.
- **Not read:** `validation.py`, `constructor-reader.py`, `restore-library.py`, `numeric.py`, `geometry.py`, `request-index.py`, `evidence-store.py`, `execution-authority.json`, `static-preparation.json`, `baseline-numerical-method.md`, `review-prompt.md`, `prior-result-record.json` and every `prior/complete/**` JSON. I did not read the `prior/` supervisor and guard logs either.
- **Consequence:** the following claims are taken from the packet's own statements, not from anything I checked:
  - the numerical engine is unchanged;
  - the saved JSON values, chain hashes, the 5007-record count and the 8271-file and 158,628,313-byte totals;
  - the hash pins;
  - the "line 216" failing guard. I only saw that `original-prepare.py:216` is the per-address `require`.

## What looks sound
- **Method:** the analytic and numerical budgets are separated cleanly. The subset-only 1e-11 ceiling, the H/flat/height/PV/slope pending status, the declined re-evaluation of the enlarged formulas, and the unknown aggregate are all stated honestly.
- **Restoration:** `restore_prefix` checks the chain sequence and hashes, all 33 budget checkpoints in order (the last refused), 263 operations in the order input < decision < return, and a conditional guard inventory. Each `J.emit` comes before its `require`.
- **Unit and domain joins:** `[-2,-1,1]` is required on every entry and on `ar['totals']` before any sum. Kappa and the carrier envelope are checked below 3. The 20/10/10/5-per-component census is joined to the preflight list and the native rows.
- **Sum structure:** the K, T and carrier loops give four group records, each with D counted three times. All four records are emitted before the single final `require` (`prepare.py:206-207`).
- **Mutants:** eight mutants per group go through `check_census`, and each must refuse.
- **Gate:** it requires both literal verdicts and a review record that pins the exact worker, manifest, library, launcher, guard, supervisor and method hashes.
- **Numerical tail:** by inspection, `continue_geometry` matches the tail of `original-prepare.py` statement for statement. The `run` function line span matches the original (85 lines in each).

## Blockers
1. **The summed values are not covered by the census.**
   - `check_census` (`resume.py:34`) compares only identity fields and the encoded `outerOperand` and `middleOperand`. It never checks `total` or the sympy `outer` and `middle` objects.
   - `prepare.py:201` sums `rows[i]['total']`, which is not a census field.
   - A row whose summed `total` disagrees with its validated operands would pass every mutant and the census.
   - The `swapped-*` mutants only refuse because they also edit the operands.
   - Required fix: recompute each total from the census-validated encoded operands, or add `total` and the sympy values to the census. Add a total-only mutant.
2. **The H scaling and monotonicity premises are not bound.**
   - The method claims identical b*, CX/CY, E30 and F2 bindings (`method.md:159-165`). `tail_formula_join` joins only the four assignment statements.
   - Base `E30=E15**2` (preflight) versus enlarged `sp.Integer(3)**30`, the definition of `b`, `Cordinary`, `F` and `bstar`, and the helpers `weighted_tail` and `exponential_moment` are not AST- or value-joined between the base and the enlarged sources.
   - The H check (`prepare.py:170-172`) is a symbolic tautology (`2^-124 = 2^-122/4`). The only code-checked premise is F2 > 0.
   - Required fix: add AST/value joins for those bindings and the two helpers (`original-prepare.py` against the pinned preflight source). Otherwise withdraw the claim that those bindings are checked.
3. **Several restored attestations are emitted but never asserted.** These are `new-source-ast-join-input`, `new-profile-tail-source-joins`, `new-outgoing-source-join` and `new-Gaussian-envelope-transport` (`prepare.py:78-79`).
   - They are re-emitted with no new join to current sources. The method calls for rechecking source/AST receipts anew. Only `new-original-source-contracts` fragments are rechecked.
   - Required fix: assert the fragment and source joins for each of them, or state them as inherited-only in the method and build documents.
4. **The "unchanged" claims have no runtime or pinned evidence.** `run` equality, `continue_geometry` equality and `tail_formula_join` on real sources are checked only by `tooling-tests.py`, which is neither source-pinned nor hash-recorded in the log.
   - Required fix: pin `tooling-tests.py` and its log in the manifest. Alternatively, run the AST equality checks inside containment against the pinned originals and emit the result.

## Coverage gaps and runtime obligations
- **Non-independence:** expected and actual census rows both derive from the same saved files. The packet discloses this. The mutants test the predicate, not the values. Carrier groups are identical by construction.
- **Rational round-trip:** expected rows use the saved raw `{text,srepr}` dicts, and actual rows use `str`/`srepr` of the restored sympy value. A non-canonical saved `srepr` would cause a false refusal. That is safe, but it is untested on native data.
- **Atomic runtime checks:** the following are runtime-only. `check_tree` (8271 files and the byte total) must hold before, after and at final posthash. The 5007-record chain and the 33 budget records must verify. Per-address ancestry and the `restored-address-*` equality must hold. The strict `K==27` branch and the unique `middle` AugAssign must resolve. The whole-run order is: restoration, then the new checks, then the four group records, then the budget guard, and only then geometry and numerics.
- **Aggregates unknown:** the four strict sums (27/122 and 29/124, zero and matching carrier) and the enlarged-versus-base component regressions are unknown. A failure must stop and be preserved.
- **Gate obligations:** the actual build review record must carry both literal verdicts, and it must pin the exact `prepare`, `resume`, `tail` and test hashes. The guard and supervisor, hook-first launch, 4 GiB/16 GiB/zero-swap/one-CPU/32-task limits, no deadline and no retry are inherited and cannot be confirmed from the packet.
- **Numerical tail:** SQLite FULL, request identities, A24/A48/B50 comparisons, mutant controls, errors and ASTs are inherited unchanged. They have no runtime evidence yet. Total time, storage and RSS are unknown.

## Remaining scoped risks
- Executed-predicate attestation depends on the original pinned control flow. Branch-dependent guards are not claimed passed.
- The reused structural parser `restore_scalar` and `walk` are disclosed, but they are still re-executed.
- The analytic majorant, D triangle lemma, shift-root integral and Gaussian identity are inherited dependencies, not independent checks.
- The aggregate is subset-only. J/direct is not full packet, current or leakage. Empirical numerical indicators are not rigorous error bounds.

**NEEDS REVISION**