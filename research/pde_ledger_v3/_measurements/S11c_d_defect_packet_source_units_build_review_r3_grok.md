I'll assess this bounded source unit-transport build from the frozen packet only: first the build notes and source, then the saved operands named in the evidence guide. No peer reports and no scientific reruns.The verdict is still open. I am reading the remaining load-bearing saved operands—derivative profile maps, the plus-source remainder head, the other grade-group denominators, the unread consumer joins, and the 544-address set—before deciding.The six grade groups and six unit operations are named in the packet. Next I am reading the denominator certificates, the unread consumer joins, and derivative profile maps.The remaining load-bearing checks are the minus-source remainder tail, one cancelled-zero return, and whether the 544-address set is actually stored. I will read those operands before any verdict.CLEAR FOR THIS BOUNDED SOURCE UNIT-TRANSPORT BUILD

Source assessment only. The worker was not run, no producer was replayed, and this is not runtime acceptance. No scientific blocker showed up in the code or the operands that were read. Wave, kernel, measure, full-summand units, evaluator readiness, and any pressure value stay deferred.

## What the build actually checks

A stored coefficient is a magnitude in the original base units. The physical coefficient is that magnitude together with the unit of the original expression. `stage_unit_join` compares each stage-2 placeholder’s registry unit with the inherited unit return of the pre-binding expression. It does not walk the bound numeric amplitude, and the normalized denominator is recorded with no unit (`worker.py` lines 221 and 307).

Held frequency and the original physical parameter stay distinct. `binding-context.json` has `physicalInput.parameters.omega` equal to `"1"` and `frequency` / `numeric.omega` equal to `Integer(3)`. `numeric_binding` (`library.py` lines 87–96) uses the held value only for the name `omega`. The registry unit is left on the symbol. `dimensionlessValueInOriginalUnits: true` is an annotation beside that value check.

Stage-2 units match the inherited returns. Effective registry `s11cc1_V_lab_held_{plus,minus}` is `[1,-1,0]` (lines 1418–1427) and `s11cc1_mu_theta_lab_held_{plus,minus}` is `[-1,-2,1]` (lines 1523–1531). `plus-native-velocity-units-return.json` is `["1","-1","0"]` against expected `[1,-1,0]`. `new-unbound-chemical-homogeneity-return.json` is `["-1","-2","1"]` against `[-1,-2,1]`. `unit_tuple` accepts those string exponents. Epsilon is `Symbol('epsilon_shape')` with no assumption keywords; `sigma_W` in the same registry is `[0,0,0]`.

The four consumer LEFT bindings match the held-frequency arithmetic. Unbound coefficients are in `inherited-consumer-unit-returns.json`. With `Lambda_A_0=1/100`, `rho_m=1/10`, `tau_A=1/10`, `W_0=1`, and held `omega=3`:

- Both `delta_p` faces bind to `-I*(100/109)*epsilon*(3/100-I/10)`, the saved LEFT in `THETA_BALANCE-delta_p_{plus,minus}-native-unit-coefficient-input.json`.
- `d_w_delta_p_plus` binds to `-I*(25/109)*epsilon*eta*w*(3/50-I/5)`.
- `d_w_delta_p_minus` is that expression with the opposite sign.

The consumer unit used later is `coefficientDimension`: `[-1,1,0]` for `delta_p` and `[0,1,0]` for `d_w`. The slot total `[-3,-1,1]` is not substituted for it.

`binding_operand` follows the pinned `bind` in `source-contracts.json` (line 51, `executed: false`): density map, then profile equalities to a fixed point, then numeric replacement. `rho_br_bg_rho4_constant` is in both maps. Density-map right-hand side `rho_br*(eta*w+1)` is the substitution that runs; its unit is `rho_br`’s `[-3,0,1]`. The profile right-hand side `W_bg*rho_br/W_0` is the same quantity after `W_bg = W_0*(eta*w+1)`. `certify` walks both right-hand sides; only the density-map side is substituted.

## Denominators, fractions, and remainders

Plus-source `D(0,0)` is `90900+303000*I`. `plus-source-regular-denominator-fraction.json` has numerator components `90900` and `303000`, denominator components `1` and `0`, and `signedNonzero` `[[true,true],[true,false]]`. Those flags are the nonzero pattern of the components. The imaginary denominator component is exactly zero, which the regularity rule allows. `plus-source-regular-denominator-fraction-reconstruction-return.json` and `plus-source-quotient-00-return.json` are the cancelled literal `Integer(0)`. Minus-source denominator and `denominatorAtZero` are the same constant (`minus-source-split.json`). `d_w_delta_p_plus-split.json` has numerator `-I*epsilon*eta*w*(6-20*I)`, denominator and `denominatorAtZero` `436`. `(6-20*I)/436` equals `(25/109)*(3/50-I/5)`, so the quotient remainder `Integer(0)` is that cancellation. The normalized numeric denominator is not given a unit.

`full_remainder_operands` builds `full` minus the four retained grades, and `N-D*retained`. Excluded pure grades stay outside that subtraction. Minus-source `excludedPure["(2, 0)"]` is the `w1_profile**2` polynomial and `"(0, 2)"` is `Integer(0)`. Its `fullHigherRemainder` is `N/D` minus the retained rectangle, and `quotientRingNumeratorRemainder` is `N` minus `D` times that rectangle. The `d_w` higher remainder is the unsimplified two-term difference of the `25/109` LEFT and the `1/436` grade piece; the quotient remainder is `Integer(0)`. Per-grade `0=0` receipts stay inherited: `grade_split` and `quotient_recurrence` are pinned with `executed: false`, and the worker sets `quotientRecurrenceRerun` false. `Pow(symbol,0)` becomes 1 for a nonzero symbol. Literal `0**0` and a nonpositive power of a zero base refuse. A negative power of a multi-term polynomial stays an opaque atom. Inputs are emitted before equality guards; a normal form is emitted when construction succeeds and before `require`. A constructor failure keeps the input and the outer failure record.

## Scale identities and the two controls

Pinned bases in `source-contracts.json` are `w=(1+tanh(x/profile_length))/2` and `m=(1-tanh(x/profile_length)**2)/3`. On saved maps, `T=tanh(composition_x/10)` with `L` `Integer(10)`:

- `profile-map-1306.json`: `w1_profile_d1` is `1/2-T^2/2`, and `m1_profile_d1` is `-10/3*(1/5-T^2/5)*T`, which is `-2T/3+2T^3/3`. Both match `(1-T^2)` times the derivative of the pinned base.
- `profile-map-1305.json`: `w1_profile_d1d1` is `(T^2-1)*T`, and `m1_profile_d1d1` is `-2/3*(T^2-1)*(3T^2-1)`. Both match the recurrence from those order-1 maps. Transverse `d2d2` and `d3d3` maps are `Integer(0)`.

The argument `composition_x/10` is the same ring element as `composition_x` times `L_W` to the power `-1`. Formal order `L` cancels. The missing-`L` control divides a live positive-order longitudinal map by `10^r`. These polynomials are not identically zero, so the mutated coefficients miss the same certificate, whose refusal text contains `actual mapped profile`. The chemical control is also non-vacuous: the constructor in `new-unbound-chemical-homogeneity-input.json` is an `Add` whose first term contains `eta_bg` and whose second term is `B_rho_3*epsilon*theta` with no `eta_bg`. Setting `eta_bg` to `[1,0,0]` makes that `Add` inhomogeneous. The engine message is `inhomogeneous native Add at ...`.

## Address span

`transport-index.json` `addressLocations` opens at `/selected/0`, `addressId` 7956, and closes at `/selected/543`, `addressId` 10576. Pointers 1, 319, 499, and 542 land on the same 10-line stride, so the stored pointer span is 0 through 543. The first `selected.json` record is the same `addressId` 7956, channel `e_W`, row `THETA_BALANCE`. The worker requires both the pointer set and the joined address-id set; it does not use `counts.addresses` or `selectedCount` as the proof. Field ids are sha256 of the field srepr and are checked as an extra identity. One field triple, `f6ca2ddb…`, has degree 4, five constant rational-complex coefficients, denominator `160114896`, and reconstruction return cancelled `Integer(0)`. `field_proof` pins that shape. It does not recompute the polynomial.

No required unit, flag, or id was found standing in for a missing identity on the operands read. Gamma units stay inference-dependent through the `gamma_` name flag. `coefficientDimensionIndependentlyVerified` is not set true. `jet_unit` joins saved `timeOrder` / `spatialOrders` to the registry derivative rule and does not re-parse the name suffix. `jet_spec` remains `executed: false`.

## Nonblocking qualifications

`dimensionlessEpsilonRemoved` and `dimensionlessValueInOriginalUnits` are constant annotations next to the real unit and value checks. Name-keyed substitution is slightly wider than the pinned `bind`, which refuses two assumption-variants of one name. The four consumer sreprs that were expanded do not contain such a pair, and their LEFT values match. The strict ring equates normal forms. It does not rewrite a higher remainder into a cancelled Taylor tail. Old per-grade zeros remain dependent on the accepted grade operation.

## Coverage limits

The `selected.json` body was not enumerated record by record, and not every one of the 544 address ids was read for uniqueness. The six inherited operation registries were not compared byte for byte with `effective-registry.json`; the worker requires full dict equality, and the stage-2 keys above were read in the effective registry. Not every one of the 392 profile maps was opened. The scale catalog forces repeated longitudinal atoms onto one polynomial, and orders 0, 1, and 2 for `w` and `m` were checked on the maps tied to the opened grade-`(0,1)` jets. An unread longitudinal order above 2 was not shown. The full chemical epsilon-join was not expanded by hand. The live repo file named by the profile contract was not opened; the contract text in the packet matches the two base assignments, with `executed: false`. Substitution right-hand-side unit walks were not executed. `verify_gate` still requires a future READY gate and both literal CLEAR reports before any run. `native/result-record.json` remains the prior bridge status with transport and evaluator pending.

`pressureSummandUnitsComplete`, `waveKernelMeasureUnitsComplete`, and `numericalEvaluatorReady` stay false in the worker’s own result object. No pressure value or field response was assessed.