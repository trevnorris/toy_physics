**NEEDS REVISION**

I read the method, implementation guide, `worker.py`, the manifest and the launcher/supervisor/helper source. I sampled the saved evidence (listed under coverage limits) and read the native c2 routes. I ran nothing.

The typed-factor method is sound. The worker keeps `Rprod`, `E` and the whole closed-density tag `Dwhole` as separate objects:
- `Rprod = qi*qo*E` is checked as a relation (`worker.py:312-316`), not as a density-minus-factor test.
- The saved original-input/canonical-return pairs are inherited, not recomputed (`worker.py:318-326`).
- The saved closed-before-cancel reference equals `Rprod`. I checked by hand that `(9+30i)/(60+91i) = β` and `β² = 9i/(60+91i)`.
- Its jet is `±i*qo*Rprod`.
- The adapter signature matches native `kernel_apply`. Evaluating the native assignments gives `normal_jet = i f qo P`, and the reference pressure becomes `P(1 − i f qo H)`. That gives the retained `(1,1)` term `P` and the excluded `(2,1)` term.

The blockers below are missing or weak checks the method requires, not errors in what is already computed.

## Scientific / validation blockers

1. **The formal noncommuting check never touches the actual row.**
   - Location: `worker.py:486-496`. It expands abstract symbols `C/F/S` and compares them with `triples()`, which only checks the enumeration.
   - Method §4b requires comparison with "direct substitution of independent response placeholders in the actual row."
   - Consequence: the actual consumer grade tables (`consumers[row][slot]`) and the address list are never reassembled into `row_bound`.
   - Minimum correction: for each row and face, substitute graded placeholders for the slot atoms in `row_bound`. Then check that each retained grade equals the sum over the 16 triples of the actual `consumerOriginal` coefficients, with ordering tags kept.

2. **`taggedTotalMixed` is not used as a reconstruction target for the assembled pieces.**
   - Location: `worker.py:371-372` joins saved symbols to saved symbols, so it is a pure identity.
   - `components_for` (`worker.py:445-450`) builds the iteration with `Hwhole/Jwhole` tags plus a `Dwhole` tag, but nothing sums them with tags mapped back and compares to `rc['taggedTotalMixed']` or `jetKernels[f].mixed`.
   - Consequence: the multiplicity-one claim for the direct addend rests on a hand-built component list that is never tested.
   - Minimum correction: assemble from `components_for`, map the tags back to the bare saved symbols, and require a zero residual against both targets.

3. **Controls use hand-typed coefficients instead of saved operands.**
   - Location: `worker.py:515` (`flat3 = 3/10 / (3/2 + β)`), `:523` (`slope = (3/10)*(3/2)/(...)`), `:533` (`normflat`).
   - Consequence: the control movements are not tied to `rc['flat']`, `rc['slopeCoefficient']` or `rc['jetKernels']`. The control addresses are asserted rather than joined.
   - I verified the constants by hand: `ω/10 = 3/10`, `k = 3/2`, and the `q` values at the proposed points are consistent.
   - Minimum correction: build each coefficient by substituting `omega=3`, `qi=q(k)`, `qo=q(l)` and `k=p` into the saved objects, and `J.zero` them against the typed constants.

4. **Response operands in the addresses carry no explicit argument map.**
   - Location: `worker.py:442-450` and `:466-478`. `responseCoefficient` keeps raw `reference_k/l/qi/qo/qm` symbols. Only the strings `'q(l)'`/`'q(k)'` describe the binding.
   - Method §4b and §5 say no response-depth binding may be inferred from equal printed names. The maps exist only for the `Jwhole` and `Dwhole` tags.
   - Minimum correction: add a per-component map record, hashed and with assumptions, from `reference_k/l/qi/qo/qm/unrestricted_frequency` to `k`, `l`, `q(·)` and `3`, with explicit flat-support substitutions.

5. **Trace routing and excluded-term checks are partial.**
   - Location: `worker.py:337-339` checks only `(T·Δ)[0,2] == P`. The required statement is that there is no `[1,2]` direct contribution.
   - Location: `worker.py:363` only requires `set(pt) ⊆ {(1,1),(2,1)}`.
   - Location: `worker.py:347` hardcodes `reference = sign/2` instead of using the bound `W_0`.
   - Minimum correction: assert `T·Δ == Δ` for the whole matrix. Assert that the excluded `(2,1)` coefficient equals `−i f qo·η·ĥ·P`. Take `W_0` from the context.

6. **Required controls are absent or fail hard.**
   - "Restored constant-end/zero-profile reductions" (method §6) are not implemented, and non-applicability is not recorded.
   - Missing nonzero control entries (`worker.py:508-522`) raise via `require`. Method §6 says to record exact absence and choose an applicable entry. In a no-retry run this would discard the single authorized run.
   - Minimum correction: either implement the reductions or record their absence. Make the control selection a recorded search with the absence documented before any failure.

7. **Jet dimensions and epsilon count are hand-typed.**
   - Location: `worker.py:292` hardcodes `baseDimension` as `[1,0,0]` for `u` and zeros otherwise. `worker.py:468` hardcodes `epsilonCount: 1`, even on zero entries.
   - Minimum correction: join the dimensions to `raw/dimension-and-measure-join.json` or native metadata, or drop the labels.

## Minor validation points

- Physical sheet (`worker.py:386-388`): only the branch expressions are compared, not the Piecewise conditions that select the propagating or decaying branch.
- Control-depth check (`worker.py:499`): it checks `depth²` only, not the sign or branch condition.
- The control-sensitivity records are labelled formal tag coefficients, which is correct. Keep it that way.

## Pure tooling

- Gate, review, launcher and manifest joins are consistent. `verify_gate` pins worker, manifest, helpers, review record and method hashes. The launcher re-asserts `command`, `outputDirectory`, `completionMessageSha256` and resources.
- The launcher is hook-first. The coordinator blocks on a 30 s startup handshake, which is not a computation deadline. A watcher-arm failure makes the coordinator refuse.
- Containment: `containment()` requires a pooled guard manifest, 4 GiB `memory.max`, swap 0, `pids.max` 32, one CPU, one thread, no CPU rlimit, and `RLIMIT_AS = 4 GiB`. The guard verifies `RuntimeMaxUSec=infinity` and `Restart=no`.
- Strict exception behavior is correct. `finally` writes `posthashes.json` and `checks.json`, with integrity mismatch forcing exit 1 and `scientificAcceptance: False`.
- `Journal.emit` rewrites `artifact-index.json` on every emit. That is O(n²) I/O but not a correctness issue.
- No gate exists yet, as the packet says. The gate must be created fresh against the corrected worker and manifest hashes.

## Optional wording

None requested.

## Coverage limits

- I only sampled the saved JSON operands. I read these in detail: `binding-context`, the THETA unit and split files, the trace and slot-routing files, and the `closed-density`, `closed-before-cancel` and closure-operand files. I did not independently re-expand the full rational quotients or the full row constructors.
- I did not read `tooling-tests.py`, the reference workers, or `shared-guard.py` beyond line 200.
- I did not check whether the saved symbol assumptions (for example `Symbol('reference_qo')` across objects) will always satisfy `one_symbol` and `exact_named_map`. A mismatch fails closed.

## Runtime evidence still needed

- Every `J.zero` and `J.sinh_zero` return, with the typed-factor, bare-linear-factor and reference/jet factor joins.
- The `direct-adapter-input` and `direct-adapter-return` records for both faces.
- The full native-row partition files, with `unknownPressureAtoms` empty.
- The `*-ordered-addresses` files and `grade-coverage`, which must show 320 coverage entries (5 rows × 2 faces × 2 slots × 16 triples), including explicit zeros.
- The post-fix assembled-row reconstruction and `taggedTotalMixed` reconstruction.
- Control operands and their exact nonzero movements.
- Complete posthashes with no mismatch.

This is a source verdict only. It accepts neither runtime restoration nor any result.