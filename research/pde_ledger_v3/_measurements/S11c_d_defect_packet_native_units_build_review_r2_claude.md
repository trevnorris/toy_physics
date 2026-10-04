**CLEAR FOR THIS BOUNDED NATIVE PRESSURE-UNIT BUILD**

I found no substantive mathematics or validation blocker. I read the sources and several saved operands by hand; I could not recompute hashes, so the pin checks are not verified here. Several non-blocking gaps are listed below.

## What I checked against source

Line numbers refer to `input/worker.py` (`w:`) and `input/library.py` (`l:`).

- **Case selection.** `join_native_case` (w:116-133) with `select_case` (l:55-64) parses all labels from the full original table.
  - It requires exactly one match, then joins index, full payload, VALUE and the payload hash.
  - The decision operands are saved before the guard at w:132, and the response uses the explicit `CASES` tag.
  - The saved operands fit this: chemical `['LAB_HELD','RHO4_CONSTANT']` at index 0, responses at indices 0 and 2, velocities at 0 and 1, density at 0.
- **Provenance records.** The four `sources` keys are distinct (b:5630, c1:103, b:3649, b:545), so the journal's exclusive-create filenames cannot collide.
- **Density.** The VALUE is `Tuple(Tuple(Equalities…), rho_br*(eta_bg*w1_profile+1), Matrix)`. Leaf `[1]` equals `densityMap['rho_br_bg_rho4_constant']`, and the expected unit is [-3,0,1] (w:182-184).
- **Source and epsilon.** The saved raw is `Mul(Pow(epsilon_shape,-1), Add(…))`. It has exactly one epsilon divisor and one remainder, which is compared with the single Add in the original DELTA_P (w:227-233). I hand-checked term 1 of the Add: [-5,1,1]+[-1,-2,1]-[-3,0,1]-[-4,0,1] = [1,-1,0].
- **Flat coefficient.**
  - The native diagonal is `Mul(omega, rho_m, Pow(q_out_output,-1), 3×DiracDelta)`.
  - After delta removal it matches `flat-join.json`'s `left`, and the expected unit [-3,-1,1] is consistent (w:219-224).
  - `_LEDGER['dtn_kernel']` membership and the `self.kernel = cases(values['dtn_kernel'])` lookup are checked structurally at w:189-200.
- **Dynamic Z override.**
  - Native `kernel_bridge` (`native-c2.py:367-377`) does `raw = inputs.kernel[...]`, `diagonal = named(raw,'FLAT_DIAGONAL')`, `z0out = diagonal.xreplace(deltas→1)`, then `DIMENSION_SCHEMA[z.name] = dimension(z0out)`.
  - The worker applies the delta-removed certify return, not a hand formula, as `face_registry[zname]` (w:237).
  - The static schema has the opaque Z at [0,0,0] for both `lab_held` faces, which matches the wrong-Z control's hardcoded [0,0,0] (w:249).
  - The resolvent symbols are statically [0,0,0] and the identity operator is [0,0,0].
  - Units: feedback [3,1,-1], Z [-3,-1,1], operand [0,0,0], pressure [-2,-2,1].
- **Controls.** Both run per face on the real expressions: the Z mutation on `definition[1]`, and the brane-density mutation on the saved raw source (w:249-257). The baselines pass first, so a refusal is attributable to the mutation. The failure-evidence ordering is correct: input, walk, decision, then require.
- **Interpreter (l:87-114).**
  - It uses AST text only and executes nothing.
  - Zero, unknown, nonrational and undefined-power cases refuse as specified.
  - All factors are visited, so zero cannot hide an inhomogeneous Add or an unknown symbol.
- **Registry.** Four registry records per gamma must agree, and origin receipts are hash-checked. The worker marks `independentOfOriginalInference: False`, so the "inherited inference, not independent" limit is preserved.
- **Gate, launcher and authority.** The routes are consistent. Worker argv and output are tied to the gate command, the guard and supervisor paths are pinned, there is no deadline, and the retry flag is false. Authority scope equals the manifest scope, and the supervisor takes the `--stage` argument generically. Posthash, copy and identity checks are in `finally`, so they run on failure.

## Non-blocking gaps (local tooling, hardening only)

1. **kernel_bridge routing.** `required_statements` (w:203) omits `raw = inputs.kernel[(anchoring, face)]` and the `z = atom_named(...)` name f-string. The saved contract equals the whole native function, so the text is in the evidence. I verified both lines by hand, but the worker does not assert them.
2. **dtn lookup chain.** `values` in `Inputs.__init__` is not tied to `_LEDGER`, and the `cases()` key mapping to `(anchoring, face)` was not read. Only the `values['dtn_kernel']` statement and literal equality are checked. The face selection (w:215-218) uses labels `[0]` and `[1]` only, with no saved index. The literal is only two labels wide, and full-literal equality covers the payload.
3. **V and mu symbols.** The source's static `V` [1,-1,0] and `mu_theta` [-1,-2,1] match the separately certified velocity and chemical units only through identical hardcoded expected tuples (w:178, w:226). Substitution of the chemical and velocity amplitudes into the source belongs to the deferred transport.
4. **A-export lines.** The `a_exports` provenance lines (a:3133, 3105, 482) are never read. The worker reads only the b/c1 lines, and the saved hashes show the a and b lines are identical for the face-velocity record. Which export `Inputs` actually uses is not shown in this packet.
5. **Failure evidence.** The `require` calls on unique dtn face, definition shape and unique feedback (w:218, 239-242) fire before their operands are emitted. The operands are still recoverable from the saved copies and the traceback.
6. **Origin receipts.** They are not cross-checked against the manifest's `opaqueInputs`. The `sourcePins` posthash covers the same paths.

## Coverage and limits

- **Read in full:** `build.md`, `worker.py`, `library.py`, `evidence-guide.md`, `launcher.py`, `inputs.json`, `execution-authority.json`, `flat-join.json`, tooling record and log, `raw-helper.py` (inherited), and `LEFTNativeSource-origin.json` (first 70 lines).
- **Read partially, by targeted excerpts:** `source.json` (provenance and case fields, plus the DELTA_P and Z excerpts), `dtn-operands.json`, `binding-context.json`, `native-c2.py`, `merged.json` (first 70 lines), `consumer-unit-joins.json` (first 60 lines), the `plus-source-input.json` head, and `source-contracts.json` (one grep match, the `kernel_bridge` prefix).
- **Not read:** `method.md`, `packet-action-method.md`, `review-prompt.md`, `tooling-tests.py`, the pickles, the other three origin files, `native-chemical-amplitude.json`, `minus-source-input.json`, `native-b-exports.py`, and the shared guard, hook and supervisor beyond the `--stage` handling.
- **Not examined:** the 200 KB quoted source lines and their hashes, because the 43 tooling tests use synthetic operands only. The minus-face data are assumed symmetric to the plus-face data, not read. Gamma units are inherited, not re-derived, and the consumer-unit returns were not rechecked.

## Runtime obligations

Inspect the following when the run happens:
- **Line and registry joins:** the actual line and constructor hash joins, the full-table label census per selection, and the registry and origin records.
- **Unit walks and controls:** every unit-walk node, the dynamic-Z contract records, and both controls per face with their refusal messages.
- **Run integrity:** the evidence chain, the guard and containment logs, and the posthashes.
- **Deep-AST failures:** deep-AST parse or recursion failures on the 160 to 206 KB constructors under the 4 GiB limit.

## Deferred, not accepted

This verdict does not accept any runtime result. The following remain required:
- source, grade and field unit transport;
- the H/J/direct-profile and measure units;
- the all-544-address summand certificates and addressed controls;
- the evaluator.

`pressureSummandUnitsComplete` and `numericalEvaluatorReady` stay false.