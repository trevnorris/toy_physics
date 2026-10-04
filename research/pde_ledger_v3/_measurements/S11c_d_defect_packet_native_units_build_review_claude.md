NEEDS REVISION

The mathematics and the worker's joins hold up on the parts I read. One provenance gap should be fixed before a gate. It is cheap to close and I did not find a mathematical blocker. I did not run any code or scripted validator; every check below is by reading.

## Applicable blockers

**B1. Chemical case membership is only partly checked against the original source** (`worker.py:150-153`, `119-127`)
- The worker checks the chemical case text hash, the case's `VALUE` tag against the value text, and that the value text is a substring of the pinned `S11c_b_exports.py:5630` constructor.
- It never parses the original case table to confirm which label that text sits under. The `LAB_HELD`/`RHO4_CONSTANT` label comes from the saved `chemicalSource.case` receipt. The original literals are label tables, for example `Tuple(Tuple(Tuple(Str('LAB_HELD'),Str('RHO4_CONSTANT')),Tuple(...)))` in `native-b-exports.py:466`.
- For the DTN literal the worker does parse labels (`worker.py:187-190`). The same pattern should apply to the chemical, response, velocity and density tables.
- The chemical expression carries no label of its own, so its case is unresolved from the source line alone.
- The other three are covered by content inside the expressions:
  - Response: the `…lab_held_plus_rho4_constant` resolvent and Z names.
  - Density: the `Equality(rho_br_bg_rho4_constant, W_bg*rho_br/W_0)` entry matches the selected leaf `rho_br*(eta_bg*w1_profile+1)`.
  - Velocity: the plus and minus cases are byte-identical (same hash `84dae281…`).
- The chemical certificate is therefore attached to real original text, but its stated case is not shown. Fix: parse the label tuples from the pinned original literal and require the selected case's labels to match.

## Non-blocking tooling and text items

- **T1. Z-unit override is a typed literal** (`worker.py:208`).
  - `face_registry[zname]=[-3,-1,1]` is hard-coded. The `certify` result at `worker.py:164` is discarded. The tie to the face's delta-removed expression is `U.same` at `195` plus the matching literal.
  - Values agree, so this is not wrong. It would be cleaner to take the unit from `certify` on each face's `delta_removed`, which also matches `kernel_bridge` (`native-c2.py:374-377`).
- **T2. `inputs.kernel` is not tied to the DTN literal.** The worker proves the literal occurs exactly once as a `_restore` constant in `S11c_c1_exports.py` (`worker.py:169-170`). It does not prove that the occurrence is the `dtn_kernel` key consumed at `native-c2.py:199`. I saw `'dtn_kernel'` at `native-c1-exports.py:85`.
- **T3. Gamma registry origins are hash identity plus inherited extraction only** (`worker.py:134-147`).
  - It checks the pickle hash and size, `sourceReexecuted is False`, `registryEntries>0`, four records per symbol, and agreement across them.
  - It does not tie merged records to per-registry entry counts. `LEFTNativeSource` says 259, and merged `staticEntries` is 823, but the worker compares neither.
  - It does not read the pickles. This is as designed and matches the "inferred original units" limit. It is unresolved first-occurrence provenance, not independent proof.
- **T4. Event file naming.** Events are saved under `*-complete-unit-walk` even when the walk is partial on failure (`worker.py:106`). This is cosmetic.
- **T5. Chemical symbol link.** The chemical symbol `s11cc1_mu_theta_lab_held_*` is linked to the chemical expression only through the matching literal `[-1,-2,1]` and the source Add homogeneity. There is no direct equality assertion. This is acceptable for a unit bridge.

## What I verified by reading

Unit order is (L, T, M).

- **Registry values** in `native-c2.py` line 59:
  - `Lambda_A_0 [-5,1,1]`, `Lambda_V_0 [-4,0,1]`, `omega [0,-1,0]`, `rho_m [-4,0,1]`
  - `rho_br_bg_rho4_constant [-3,0,1]`, `tau_A` and `tau_V [0,1,0]`
  - `mu_theta_lab_held_{±} [-1,-2,1]`, `V_lab_held_{±} [1,-1,0]`, `q_out_output [-1,0,0]`
  - `dtn_operator_lab_held_{±} [0,0,0]` (the static placeholder), and the resolvent and identity symbols `[0,0,0]`
- **Source** (`plus-source-input.json`): `Λ_A μ/(ρ_br ρ_m)` gives `[1,-1,0]`. The second term `V(Λ_V/ρ_m·(1−iωτ)⁻¹+1)` is homogeneous. The raw source is `Pow(eps,-1)·Add`, and its Add matches the Add inside the original `DELTA_P`, term by term and in order.
- **Flat:** the diagonal is `Mul(omega, rho_m, Pow(q_out_output,-1), 3×DiracDelta)`. After delta removal it equals the saved flat srepr and gives `[-3,-1,1]`. The original DTN `DIMENSION_L_T_M` of `[0,-1,1]` is consistent with the three deltas contributing `[3,0,0]`.
- **Inverse operand** (original text):
  - It is `Add(Λ_A ρ_m⁻² (1−iωτ_A)⁻¹ Z, identity)`.
  - The feedback coefficient is `[3,1,-1]`, so the operand is `[0,0,0]`.
  - `DELTA_P = Mul(Add(source), resolvent, Z)` has exactly one Add, giving `[-2,-2,1]`.
- **Dynamic Z override:** `kernel_bridge` (`native-c2.py:367-378`) matches the four required statements the worker checks for.
- **Controls:**
  - The bulk-density mutation refuses on the inner `Λ_V/ρ_m'+1` Add, with a dimension of `[-1,0,0]` against dimensionless.
  - The static-Z mutation refuses on `identity + feedback·Z` with `Z=[0,0,0]`.
  - Both must pass the baseline first, and both require the message to start with `inhomogeneous native Add`.
  - On failure they leave `activeOperation` set, and a failure record is written.
- **Interpreter** (`library.py`):
  - It accepts only Symbol, Integer/Rational, Add, Mul, rational Pow, the six transcendental functions, DiracDelta, and the names `I`, `pi` and `E`.
  - There is no tolerance, no simplification and no executed constructor.
  - Zero is `None`. Zero factors still evaluate every symbol, so they cannot hide unknowns. Undefined zero powers and unknown names refuse.
- **Journal ordering:** the input is saved before the walk, the walk in `finally`, the decision before `require`, and the return only on success. Posthashes, copies and identity copies are saved in `finally`.
- **Gate, argv and launcher:**
  - The gate checks the worker, manifest, source pins, guard, supervisor, launcher, library, authority, build record, both literal build verdicts, the method record and the scope.
  - The worker checks that `sys.argv` and the output directory match the gate command.
  - The launcher compares the full command, runs the hook first, and arms no deadline. `runpy` of the worker is guarded by `__name__`, so `main` does not run.
  - The authority has 1 execution, no retry and no deadline.
  - The containment helper checks `memory.max=4GiB`, swap 0, `pids.max=32`, one CPU, one thread per library, and no CPU rlimit.
- **Gate dependency:** the gate requires the literal string `CLEAR FOR THIS BOUNDED NATIVE PRESSURE-UNIT BUILD` from both reports. This verdict will not satisfy it, which is the intended outcome until B1 is fixed.

## Coverage and limits

**Read fully:**
- `build.md`, `worker.py`, `library.py`, `launcher.py`, `evidence-guide.md`
- `inputs.json`, `packet-index.json`, `execution-authority.json`
- `tooling-tests.py` and its record, the saved flat join, and the `LEFTNativeSource` origin file
- the full DTN literal
- `native-c2.py` lines 360-410 and 725-775
- the `containment` function in `raw-helper.py`

**Read partly:**
- `source.json` (904 KB): lines 26-44, 115-144, 190-229 and 392-432, plus targeted fragments of the 160-206 KB response constructors. I saw the inverse-operand definition and the Add and tail of `DELTA_P` but not the full texts.
- `plus-source-input.json`: lines 1-21 of 41.
- `native-chemical-amplitude.json`: the `raw` expression only.
- `binding-context.json`: lines 60-160.
- `merged.json`: the head and the required-symbol list. The 259-entry per-symbol registry bodies were not read.
- The C1, B and C2 exports: targeted greps only.

**Not read:** `minus-source-input.json`, `consumer-unit-joins.json`, `source-contracts.json`, `physical-input.json`, the three other origin files, `method.md`, `packet-action-method.md`, `review-prompt.md`, the tooling log, the shared guard, the supervisor, the hook, `AGENTS.md`, and the pickles (identity only, as designed).

**Not verified:** units of the roughly 200-term chemical expression and of every other static-schema symbol. I checked only the symbols listed above, and the other three origin files remain unchecked.

## Runtime obligations (not accepted by this review)

- Every constructor-text parse, all unit-walk node decisions, the registry join with no `unavailable` entries used, and every certify result must be inspected after the run.
- The run must show both controls refusing on the actual Adds, on both faces.
- The posthashes, evidence chain, guard, resource and outcome logs, and a nonzero-exit or failure record must be inspected.
- Memory and recursion behavior on the 61 MB, 206 KB and 161 KB literals must be checked at runtime.
- The result keeps `pressureSummandUnitsComplete` and `numericalEvaluatorReady` false. Source-grade/field, H/J/direct-profile/measure and all-544-address unit transport remain required.