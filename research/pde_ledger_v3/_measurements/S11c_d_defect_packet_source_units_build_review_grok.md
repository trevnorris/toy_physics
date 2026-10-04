I'll assess this bounded source unit-transport build from the frozen packet only, starting with the build notes and the worker/library, then the saved inputs those point to.The remaining question is whether the retained grades, jet reconstructions, and address census are real identity joins or stand-ins. I am reading those saved operands and the worker checks that consume them.**NEEDS REVISION**

The homogeneous-rational coefficient argument does not hold for these saved domain and unit operands. A source-only reading is not a runtime certificate. Gamma dimensions stay inference-dependent. Wave, kernel, measure, full-summand, and evaluator units stay deferred. No pressure value is required.

## Blocker 1 — one retained grade is an inhomogeneous sum, and pure numeric jet coefficients are then given different units

The plus unbound source return is `["1","-1","0"]` in `saved/native/plus-unbound-source-units-return.json`. `worker.py` line 232 copies that unit onto every nonzero retained grade. Lines 262–269 then set each independent-jet coefficient unit to that source unit minus `jet_unit`.

The retained plus `(0, 0)` grade that receives this unit is the sum in `saved/composition/plus-source-jets-00.json` lines 6–8 (the same sum is the split’s retained `(0, 0)`). Its rational coefficients multiply these registry entries from `saved/native/effective-registry.json`:

| symbol in that sum | registry lines | unit |
|---|---|---|
| `e_W` | 353–356 | `[0, 0, 0]` |
| `e_W_d1d1` | 368–371 | `[-2, 0, 0]` |
| `e_W_t` | 713–717 | `[0, -1, 0]` |
| `theta` | 1673–1676 | `[0, 0, 0]` |
| `theta_d1d1` | 1678–1681 | `[-2, 0, 0]` |
| `u_1_d1` | 2043–2046 | `[0, 0, 0]` |

`eta_bg` is `[0, 0, 0]` (lines 728–732) and `sigma_W` is `[0, 0, 0]` (lines 1612–1616). No `L_W` factor stands in front of `e_W_d1d1` or `e_W_t`. Under `Units.dimension` in `runtime-source/S11c_d_defect_packet_native_units_lib.py` lines 100–101, an `Add` whose live units include `[0, 0, 0]`, `[-2, 0, 0]`, and `[0, -1, 0]` raises `inhomogeneous native Add`. The worker never walks this grade. The `[1, -1, 0]` value is the pre-binding constructor unit joined at `worker.py` line 200, not a unit of this bound sum.

The jet table then treats the normalized numbers as the coefficients:

- `e_W` in `plus-source-jets-00.json` lines 26–43: coefficient and field are both `5/218 + 3*I/436`, `jetDimension` `[0, 0, 0]`, `requiredCoefficientDimension` `[1, -1, 0]`.
- `e_W_d1d1` in the same file, lines 100–120: coefficient and field are both `-1/109 - 3*I/1090`, `jetDimension` `[-2, 0, 0]`, `requiredCoefficientDimension` `[3, -1, 0]`. `coefficientDimensionIndependentlyVerified` is false.

`worker.py` line 263 only checks that those stored triples match `source unit − jet unit`. That gives two different physical units to two normalized numeric coefficients of one sum. The denominator path does the separation the jet path does not: `worker.py` line 193 sets `normalizedBoundDenominatorAssignedUnit` false, and `library.py` `field_proof` lines 115–122 does not dimension a constant polynomial denominator. Address `7956` (`transport` source map lines 13813–13819) points its source field at `fields/79ea957c…-polynomial.json`, whose polynomial is that same `5/218 + 3*I/436` over denominator `Integer(1)`. The unit on that address is the jet’s required unit `[1, -1, 0]`. The field id checked at `worker.py` line 292 is an extra byte identity.

The same pattern is in the `(1, 0)` grade: `plus-source-jets-10.json` lines 27–42 give `e_W` the coefficient `-I*w1_profile*(12 - 40*I)/1744` and required unit `[1, -1, 0]`, while `e_W_d1d1` remains length⁻² in the same sum (`plus-source-jet-reconstruction-10-input.json` lines 2–8). The reconstruction return’s `cancelled` is `Integer(0)`, so that file really joins two writings of this grade. It does not make the sum homogeneous.

## Blocker 2 — the quotient receipts are `0 = 0`, and the stored remainder is never joined

`worker.py` lines 242–245 join saved-zero for retained `(0, 0)` only, then for grades `00`, `10`, `01`, and `11` require only that the quotient file’s right side and `cancelled` are `Integer(0)`. They then accept the booleans `nativeShapeCoefficientsCalled is False` and `zeroGradeReused is True`.

`plus-source-quotient-00-input.json` and `plus-source-quotient-01-input.json` are both `{left: Integer(0), right: Integer(0)}`. The same 124-byte shape is reused for other quotient and epsilon slots (for example `THETA_BALANCE-d_w_delta_p_plus-epsilon-11-input.json` and `THETA_BALANCE-delta_p_plus-epsilon-10-input.json` in `inputs.json`).

`plus-source-split.json` lines 42–51 store the actual remainder operands. `fullHigherRemainder` and `quotientRingNumeratorRemainder` are unsimplified symbolic differences, not `Integer(0)`. `excludedPure` `(2, 0)` is a nonzero `w1_profile**2` polynomial; `(0, 2)` is `Integer(0)`. The worker never reads those three fields. The producer fragment embedded in the source map, `grade_split` at lines 19273–19279, emits each quotient file from `terms.get(g)` after cancellation, so a stored `0 = 0` does not still carry the remainder expression.

`selectedCount: 544` is checked at `worker.py` line 250. The address loop at lines 284–317 certifies `index['addressLocations']` and then sets `BOUNDED_SOURCE_GRADE_PROFILE_UNIT_TRANSPORT_COMPLETE` without requiring `len(addresses) == 544` or equality of the joined address ids with the selected ids. In the frozen map, `addressLocations` begins at `/selected/0` (line 13814, address `7956`) and the last object read is `/selected/543` (line 19245, address `10576`). Pointers `1`, `3–6`, `69–71`, `319–321`, and `537–540` were also present. Contiguity of all 544 pointers was not read, so that remains unknown. The inventory count is not that identity.

## Controls and persistence

The chemical control is applicable to this constructor. `saved/native/new-unbound-chemical-homogeneity-input.json` line 2 is an `Add` whose first summand contains `eta_bg` (`B_rho_3 * epsilon_shape * eta_bg * theta * w1_profile`) and whose second summand is `Mul(B_rho_3, epsilon_shape, theta)`. With `eta_bg` changed from `[0, 0, 0]` to `[1, 0, 0]`, those two terms differ by one length, and `Units.dimension` raises `inhomogeneous native Add`. `worker.py` lines 314–316 emit the decision operands before the require. `main` lines 334–343 write the failure record, posthashes, and `checks.json` before the posthash require.

The missing-`L` control does not apply to the transport function. `worker.py` lines 303–307 store `{'LExponent': 0}` and then evaluate `add(scale(registry['L_W'], 0), (-order, 0, 0))`. `library.py` `profile_transport` lines 96–97, which is where `L_W^order` is actually required to cancel the derivative unit, is not called. For any saved order greater than 0 the side formula is `(-order, 0, 0)`, which differs from zero whether or not `profile_transport` still checks `L_W`. Registry `L_W` is `[1, 0, 0]` (lines 23–27). `profile-map-1219.json` carries `L = Integer(10)`, matching `binding-context.json` parameter `L_W` `"10"`, and its map is empty because that coefficient is the pure number above. Derivative profile symbols that were read, including `m1_profile` at line 798 and `w1_profile_d2d2d2` at lines 4003–4007, are `[0, 0, 0]`.

## What this reading did and did not cover

Read: `build.md`, `worker.py`, `library.py`, the native unit engine, the plus `(0, 0)` and `(1, 0)` jet tables and the `(1, 0)` reconstruction pair, the plus quotient `00` and `01` inputs, the plus split remainder tail, the effective registry rows cited above, the chemical constructor head, the velocity-walk head, consumer epsilon `00` for `delta_p_plus`, profile map `1219`, field polynomials `79ea…` and `9937…`, and the address-location endpoints in `S11c_d_defect_packet_pressure_transport_source_map.json` (manifest sha `51b94650…`, `inputs.json` lines 3035–3038).

The stage2, chemical-domain, and live-density joins at `worker.py` lines 199–219 are textual identity requires. Their operand bodies were not re-derived. Reciprocal-base completeness inside all six walks was not re-read; `library.py` lines 35–44 refuse a negative power that lacks its `/base` event. Whether every pointer from `/selected/0` through `/selected/543` is present is unknown. The velocity walk’s first factors are `Rational(1, 2)` at `[0, 0, 0]`, `W_0` at `[1, 0, 0]`, and `e_W_t` at `[0, -1, 0]`; that is a component walk, not a unit for the normalized denominator `Integer(1)`.

`worker.py` line 317 keeps `pressureSummandUnitsComplete`, `waveKernelMeasureUnitsComplete`, and `numericalEvaluatorReady` false. This source transport is not full readiness.