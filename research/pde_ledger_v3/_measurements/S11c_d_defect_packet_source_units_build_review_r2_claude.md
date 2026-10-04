**Verdict: NEEDS REVISION.** I found no error in the saved science, but the new substitution provenance has gaps the worker does not close. The fixes are small and code-only. I did not run anything, so I cannot confirm that any ring equality passes on the real operands.

## Blockers

**B1. The substitution unit and value joins are incomplete.**
- `worker.py:173-176` unit-checks only `profileEqualities` and `densityMap` against the registry. The stage2 map (`worker.py:220-221`, `236-239`) replaces `s11cc1_V_lab_held_<face>` and `s11cc1_mu_theta_lab_held_<face>` with the saved amplitudes. Nothing asserts the unit of the replaced symbol equals the unit of its replacement.
- The saved data does satisfy this. The registry gives V = (1,-1,0) and μθ = (-1,-2,1). These match `native/plus-native-velocity-units-return.json` and `new-unbound-chemical-homogeneity-return.json`. The check is simply absent from the worker.
- The `numeric-binding-provenance` record (`worker.py:181-185`) takes each value from `context['numeric']` and its unit from the registry. Only `L_W` is joined to `physicalInput.parameters` (`worker.py:186`).
- The binding uses omega = 3. `physical-input.json` has `omega: "1"`. The scope frequency is 3, but no code ties `numeric['omega']` to `frequency`.
- The other values I spot-checked do match `physicalInput.parameters`: B_rho_3, C, tau_A/V/X, kappa_theta(_W), G_theta_u, rho_m, Lambda_A_0, and several gammas.
- Required fix: assert both unit equalities, assert every `numeric` value equals its `physicalInput` parameter, and assert `numeric['omega'] == frequency`.

**B2. Saved flags stand in for an identity on the six grade groups.**
- `library.py:55-56` accepts `finite` and `signedNonzero` booleans. The exact component literals sit in the same file (for example `plus-source-regular-denominator-fraction.json`: 90900, 303000, 1, 0). They are never rederived or checked against the numerator.
- For the six rational groups, `split['denominatorAtZero']` is never joined to `split['denominator']` (`worker.py:276` uses it only in `regularity`). The chemical domain does get this join (`worker.py:252`, `at_grade_zero`). The regular-domain premise for the Taylor coefficients is therefore unlinked from the actual denominator.
- `full_remainder_operands` (`library.py:211-217`) builds only `full − retained` and `N − D·retained`. `full == N·D⁻¹` is never tested, so `build.md` overstates when it says the remainders are joined to the full expression, numerator and denominator. The two assemblies match the saved remainder definitions, so they do not independently show the retained terms are Taylor coefficients. That dependency is stated honestly as inherited.
- Required fix: rederive nonzero from the component literals, add the `denominatorAtZero` join, and add the `full == N·D⁻¹` check. Alternatively, narrow the wording.

## What looks correct

- **Scale identity.** The recurrence P_r = (1−T²)·P′_{r−1} is correct (`library.py:243-249`).
  - I checked it by hand on `profile-map-1361`: w″ = T³−T and m″ = −2T⁴+8T²/3−2/3 both match the saved mapped values.
  - The producer at `source/S11c_d_defect_source_composition.py:417` computes `L**r · diff(...)`, which is what the identity certifies.
  - The base assignments sit at line 410, columns 4 and 38, as pinned.
  - L = 10, so the missing-L control (`worker.py:340-353`) does change coefficients. It runs through the same `verify` path, and its refusal text contains "actual mapped profile".
- **Native control.** The `eta_bg` → length control (`worker.py:355-362`) should raise. The chemical `Add` mixes terms with and without η·w1_profile, and the engine message contains "inhomogeneous".
- **Stored magnitude versus physical quantity.** The two are kept distinct. The coefficient unit is the group unit minus the jet unit, with the jet unit checked against the registry. This follows from the inherited unbound homogeneity plus checked substitutions. The normalized denominator gets no unit.
- **Address coverage.** There are 544 selected addresses, all with channel `e_W`. The pointer-set and joined-ID checks are both present. The address statuses 102 + 106 + 336 sum to 544.
- **Chemical LEFT joins.** The ring cancels ε⁻¹·ε structurally. The negative-power bases in the saved bound source are (η·w1+1)⁻¹, ρ_m⁻¹ and (1−iωτ)⁻¹. They reduce to a shared opaque atom or to constants, so the expected and saved forms should compare.

## Non-blocking qualifications

1. **Persistence.** `regularity` emits its evidence only after its guards pass (`library.py:64`). The tanh-argument ring forms (`library.py:230-235`) are lost if the x/L guard refuses. A ring `unsupported constructor` error leaves no decision file. Inputs still persist through `copy_inputs` and `J.start`.
2. **Inherited-walk registry.** The registry embedded in the inherited walk inputs is not asserted equal to `effective-registry.json`. They look identical for the keys I sampled.
3. **Consumer binding.** The unbound-to-bound binding for consumer coefficients is inherited. `prior(name+'-native-unit-coefficient', right=origin)` (`worker.py:273`) never checks its left operand against `item['coefficient']`, unlike the source side.
4. **Field derivation.** `oneDimensional` is not re-derived from `original` and the map. The check that fields contain only `composition_x` is also inherited.
5. **Ring limits.** Multi-term bases cannot cancel against their own inverses, and products of non-constant bases would be refused. This is a fail-closed risk of a false refusal at runtime, not a false pass.
6. **L join strength.** The x/L join identifies L by its value 10 only.
7. **Smaller points.**
   - The `formal==ZERO` check in `profile_transport` is tautological.
   - The 64/16 location counts are not asserted.
   - Gamma units are not tagged as inferred per binding in the journal. They are noted only in the final limits.

## Coverage and limits

- **Read in full:** `build.md`, `evidence-guide.md`, `worker.py`, `library.py`, `runtime-source/..._native_units_lib.py`, `binding-context.json`, the `plus-source` regular-denominator fraction and components files, profile maps 1219, 1361, 1379 and 2300, and the consumer `delta_p_plus` epsilon inputs 00, 10 and 11.
- **Read in part:**
  - `inputs.json`: head, tail and key greps only.
  - The chemical native input: lines 1-298 of 4277 (the full constructor line, but only part of the registry).
  - `plus-source-split.json` and the plus stage2 operands/bound/epsilon/velocity records: fragments only.
  - `transport-index.json` and `selected.json`: counts and samples.
  - The effective registry: greps.
  - The source-contract profile definitions and source lines 400-430.
- **Not read:** about 390 of the 392 profile maps, all minus-face records, the other five grade groups' saved operands, the 34 field records, the inventory and whole-coefficient transforms, `launcher.py`, and the tooling tests.
- **Not established:** no ring equality on a real operand was executed, so I cannot say the runtime joins pass. My static findings do not extend to the unread records.

Wave, kernel, measure and full-summand units and evaluator readiness remain deferred. Gamma units remain dependent on the inferred registries.