**Verdict: CLEAR FOR THIS BOUNDED PRESSURE-ASSEMBLY READINESS METHOD**

I found no method defect. The items below are implementation obligations and efficiency concerns. This verdict authorizes faithful build preparation only, not a worker, an evaluated convolution, current or loss.

## Unit method: checked, with no method defect

I hand-checked the unit arithmetic against the supplied native sources. No code was run.

- **Parameter units.** `rho_m` is (-4,0,1) and `rho_br` is (-3,0,1) (`c2.py:59`). The saved native resolvent definition is `I + Λ_A0/(ρ_m²(1-iωτ_A))·Z` (`consumer-census.json:87`). This gives a=(3,1,-1) and β=aμ=(-1,0,0), as the method states.
- **Source.** The raw source is `Λ_A μ_θ/(ρ_br_bg_rho4_constant ρ_m(1-iωτ_A)) + V(Λ_V/(ρ_m(1-iωτ_V))+1)` (`plus-source-input.json`). Each term has unit (1,-1,0).
- **Slot coefficients.** The pressure-slot coefficient is (-1,1,0) and the normal-jet coefficient is (0,1,0) (`consumer-unit-joins.json`).
- **Why the unbound join is needed.** The numeric templates hide the units. `unrestricted-closure-operands.json` shows a=1/(100(1/100−iω/1000)), which is 1/(1−iω/10) because Λ_A and ρ_m² cancel. The saved Λ_A/ρ_m and τ_A are both 1/10 but have different units. Λ_V0 and Λ_X0 are 0, so the λ_V/ρ_m term vanishes numerically.
- **W and L origins.** `raw-increment.py:443,448` and `inner-library.py:53-60` show the explicit a·μ²·W·L/4 and W·L/(4i) factors. Hand dimension sums give (-1,-1,1) for the J, reflected, height and quadratic densities, and (-2,-1,1) for the H, slope and height templates.
- **Profile scaling.** `native-profile-scale-join.json` and `c2.py:264-265` both put L^r on each profile derivative.
- **Pairing unit.** The result is (-2,-1,1) for both the flat and off-diagonal routes, and it matches `local-unit-conclusion.json`. That file correctly records `pressureSummandUnitsEstablished:false`.
- **Collision arrangement.** The lines, cut offsets and slab construction match `packet-action-method.md` §3–4 and `inner-library.py:31-42`. The k+l=0, ±2κ, l=k and l=±κ coincidences are complete for the root and profile structure. The k,Q height domain is carried separately.

## Implementation obligations

These are not method defects. Each is needed for a faithful build.

1. **Operator-symbol units.**
   - The literal schema gives `s11cc1_dtn_operator_*` unit (0,0,0). Taken literally, `1+aZ` has unequal addends and would trip the method's own stop rule.
   - `c2.py:377` overrides Z with the unit of the delta-removed FLAT_DIAGONAL. That is ρ_m·ω/q (`c1.py:602-603`), with unit (-3,-1,1).
   - The build must use this override and join it to `c1.py:1088-1093`. It must not use the literal 0, and it must not adjust anything to force agreement.
2. **Gamma units are inferred, not independent.**
   - `c2.py:221-237` solves each gamma unit from the first term where it appears, and `mu_theta` is processed first.
   - Gamma-bearing e_W terms are homogeneous by construction. They include e_W_d1 and the L⁻¹ terms with gamma_w_bg_13/14 and gamma_mu_r_bg_13/14.
   - The certificates should label those addresses "inherited inferred registry unit". The independent check is that each gamma is consistent across all its occurrences in the other rows. `THETA-row.json` contains gamma_w_bg_13 13 times.
3. **Inherited nan entries.** `units-and-scope.json` has nan entries in `nativeConsumerDimensions` and labels its dimensions "inherited declarations". The new proof must not inherit those entries.
4. **Exact zeros.** `plus-source-jets-00.json:44,1459` marks `coefficientDimensionIndependentlyVerified:false` and says the required units are "inferred expectations". The 442 exact-zero addresses can carry only an expected unit with their receipts.
5. **Symbolic transform record.** The saved symbolic transform record (`saved['transform']`) is not supplied. `Fourier-profile-lemma.json` is numeric at L=10. The L origin must come from `raw-increment.py:443,448`, evaluated as a formula join, not as a recomputed transform.
6. **Route independence.**
   - The old `inner-library.py:290-293` assembles every route from one shared H and the B-route J/D. A full outer Route A/Route B build must not do that.
   - Each outer route needs its own H, J/D and Fourier chain. Cache keys must include the outer route, so Route B never consumes a Route A value even when arguments coincide.
7. **Request keys.** `fourier-library.py:154-158` accepts only rationals and ±√595/10 as arguments. Arbitrary outer nodes need an exact canonical argument descriptor, including how the 30-digit and 50-digit arguments relate.
8. **Old-bank reuse.** The old-bank values are opaque SQLite receipts only. Outer nodes are generally irrational, so reuse will be essentially zero. Count the old banks as no savings unless a full-request match is verified.
9. **Exact geometry.** All line coefficients lie in ℚ(√595). The build can and should compare endpoints exactly.

## Efficiency concerns

- **Storage.** The inherited bank sizes are 12,338,020,352 bytes for 38 inner points (about 325 MB each, 4371 s) and 4,961,550,336 bytes for 190 Fourier requests (about 26 MB each, 2770 s). These are scale evidence only. They indicate that the full-evidence format at outer-node scale is far beyond 4 GiB containment, and the method's own storage gate is the correct response.
- **Inner work per outer node.** `inner-library.py:134-138` uses t-panels of width at most 1/2 over [-T,T] with T=122, so at least 488 panels per node. Inner integrals are not shared across distinct (k,l).
- **Fourier work for Y(l).** The nested l-bounds depend on k, so Y(l) requests scale roughly with the number of outer nodes. The Q-route reaches |l|≈102.

## Remaining limits

- Not established here: any runtime, any numerical agreement, any pressure value, and the correctness of the common analytic kernel convention.
- Adaptive Route B counts and unit-tag-conservative reuse limits remain unknown until execution.
- The existing no-deadline guard, supervisor and containment rules are unchanged.

## Evidence read

- **Method files:** `method.md`, `evidence-guide.md`, `review-prompt.md`, `packet-action-method.md`.
- **Source:** `c1.py` (dimension definitions and lines 1040–1140), `c2.py` (schema, `infer_dimensions`, `at_source`, `kernel_bridge`, `reference_pressure_kernels`, `build_face`), `raw-increment.py` (lines 425–555), `inner-library.py` (full), `fourier-library.py` (lines 1–273), and `preflight.py` lines 322–358.
- **Saved evidence:**
  - Native census: `consumer-census.json` lines 53–117.
  - Consumer inputs: `consumer-unit-joins.json`, `units-and-scope.json`, `binding-context.json`, `inherited-source-normalization.json`, and the key lists of `plus-source-input.json` and `native-chemical-amplitude.json`.
  - First-shape and closure operands: `first-shape-native-transport.json`, `unrestricted-closure-operands.json`, `raw-ordered-before-cancel.json`.
  - Selected addresses and fields: `pressure-addresses.json` (head), `pressure/fields.json` (head), `pressure/whole-tags.json` (head).
  - Sampled field and source-ancestry files: one `field/*-reconstruction-input.json`, and `plus-source-jets-00.json`, `plus-source-quotient-00-input.json`, `plus-source-jet-reconstruction-00-input.json`, `plus-native-source-join-status.json`.
  - Preflight results and inventory: `inherited-unit-contract.json`, `physical-plan.json`, `Fourier-profile-lemma.json`, `checks.json` (K, U, T lines), `native-profile-scale-join.json`, `fourier-and-unit-provenance.json` (partial), the first `numeric-factor-plus-pressure-NATIVE_FLAT` input and arguments, and gamma entries of `native-units/units/merged.json`.
  - Banks and local result: `local-unit-conclusion.json`, `fourier-01/` and `inner-01/` `checks.json` and receipts.

**Not read:**
- The 13 MB slab constructor and the full `THETA-row.json`.
- The remaining 33 field polynomials and returns, and the other seven source-ancestry tables.
- The consumer grade-split files, the other 19 numeric-factor adapters, and the H/J/D envelope and bound files.
- `closed-grazing.py`, `reference-grazing.py`, `source-composition.py`, `local.py` and most of `preflight.py`.

I did not check every one of the 102 formal addresses; that is the build's job. No tool failures occurred.