**CLEAR FOR THIS BOUNDED SOURCE UNIT-TRANSPORT BUILD**

This verdict rests on a partial reading of the saved operands, listed under "Actual read coverage" below. I found no scientific blocker in what I read. I ran nothing, so this is a source assessment, not runtime acceptance.

## Findings on the points you asked about

**Physical parameters.**
- All 38 numeric bindings in `saved/original/binding-context.json` match `physicalInput.parameters`, except omega.
- Omega is the only exception: held 3 against the original physical input of 1. The two are kept separate (`library.py:87-96`, `worker.py:179`).
- Every bound name has a registry unit. The gamma entries sit at `effective-registry.json:4118-4263`.
- Gamma units stay marked as inference-dependent (`gammaUnitInferenceInherited`, `worker.py:387`).

**Stage2.** The V and mu_theta placeholder units are compared to the inherited velocity and chemical unit returns (`worker.py:234-237`). The chain is unbound chemical constructor, then bound LEFT operand, then amplitude (`worker.py:230`, `245-248`, `257`). I found no unit, flag or ID standing in for a missing identity.

**Fraction identities and component evidence.**
- Full = N/D uses cross multiplication with every inverse base recorded (`library.py:220-249`).
- The inverse-base certificate evaluates at the grade origin and requires an exact finite nonzero constant (`library.py:364-367`).
- `regularity()` rebuilds the components as exact rationals. Saved flags are only compared, never trusted (`library.py:55-84`).
- For plus-source, the saved fraction components agree with those rules (numerator 90900+303000I, denominator 1).

**Profile scale identity.**
- I hand-checked four saved mapped values against P_r = (1−T²)·dP_{r−1}/dT: w1_d1, m1_d1, w1_d1d1 and m1_d1d1. All four match.
- The argument form `tanh(x/10)` with `composition_x` matches the check at `library.py:295`.
- Saved profile orders reach only 2, so the recurrence covers orders 0–2.

**Controls.** Both controls follow a code path that refuses. The missing-L control divides a live mapped value by L^r and sends it through the same `profile_transport` and `ScaleVerifier`. The grade-length control makes eta a length, which creates an inhomogeneous Add in the chemical constructor. Whether they refuse at runtime is unknown.

**544-address coverage.** `worker.py:343-345` and `386` require both the pointer set and the joined address-ID set to be exactly 544 unique entries. `transport-index.json` contains 544 `selectedJsonPointer` entries.

## Nonblocking qualifications

1. **Tautological check.** The formal `L_W^r × L^−r = 0` check (`library.py:133`) is tautological. The real content is the recurrence.
2. **Strict ring normal form.**
   - Opaque negative-power atoms never cancel against a matching polynomial factor.
   - Only integer powers are supported.
   - Symbols with different assumptions are distinct atoms. The saved data mixes them (`Symbol('m1_profile_d1d1')` in map 1648 against `real=True` in map 1644).
   - Any such mismatch makes the assembly refuse. It cannot produce a false pass.
3. **Large assemblies not hand-verified.**
   - These are the ~150-term `full-bound-left`, chemical/velocity LEFT, consumer unit-left and remainder assemblies.
   - A runtime refusal is therefore possible and is the main residual risk.
   - Inputs and normal forms persist before each equality guard.
4. **Inherited from the old grade operation.**
   - The per-grade zero receipts are inherited.
   - Retained-coefficient values are inherited.
   - `fullHigherRemainder` is only joined as the definition full − Σ retained, which the build discloses.
5. **Scope of the source-jet work.** Only the e_W channel (8 jets per grade table) is transported. The theta and u jets in those tables are not needed by the selected addresses.
6. **Control 2.** It covers the chemical expression only.

## Actual read coverage

**Read fully:**
- `build.md`, `evidence-guide.md`, `worker.py`, `library.py` and the inherited `native_units_lib`.
- `binding-context.json`, `physical-input.json` and `tooling-test-record.json`.
- `plus-source-operands.json`.

**Read in part (grep, partial reads, or the lines I viewed):**
- `effective-registry.json`: lines 1–3246 of 4273, plus the gamma lines by grep.
- `plus-source-split.json`: structure and srepr heads only.
- `plus-source-regular-denominator-*`, `plus-chemical-domain-zero-grade-join-input` and `chemical-amplitude-domain.json`.
- `plus-source-jets-00.json`, `plus-native-source-join-operands.json` and `inherited-consumer-unit-returns.json` (grep only).
- `plus-unbound-source-units-return.json`.
- `source-contracts.json`: lines 1–40 only. I did not independently confirm the pinned source-line and fragment match at `library.py:105-111`.
- `native-profile-scale-join.json` and `transport-index.json` (grep, with the 544 pointer count).
- Profile maps 1221, 1226, 1264, 1273, 1644 and 1648, plus a grep survey of profile atom names.
- `S11c_d_defect_source_composition.py:296-345`, which contains `bind`.
- `inputs.json`: lines 1–150 only.

**Not read:**
- Minus-side files and all consumer group operand/split/denominator files.
- Most of the 392 profile maps, and `fields/` and `inventory/`.
- The native unit-walk files and most join input/return files.
- `launcher.py`, the execution authority, `method.md`, `tooling-tests.py` and its log. The 73 passing synthetic tests are unverified.

**File-existence check.** I checked by glob that the saved filenames the worker needs for the plus-side and consumer-group joins exist.

For the unread operands, this verdict relies on the worker's fail-closed design.

## Unknowns preserved

- Runtime acceptance and the outcome of the large assemblies.
- Wave, kernel, measure and full-summand units, and evaluator readiness.
- Gamma unit inference.