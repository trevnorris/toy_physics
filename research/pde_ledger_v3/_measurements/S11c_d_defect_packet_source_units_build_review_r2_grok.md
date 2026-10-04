I'll assess this bounded source unit-transport build from the frozen packet only, starting with the build notes, worker, library, and the saved operands those files name.The packet is under `input/`. I'll read the build, worker, library, and the operand map before judging the joins and unit provenance.The scale maps and the chemical constructor are in hand. I am reading the remaining joins that can still change the verdict: stage-2 unit provenance, the remainder operands, the grade-length control’s registry row, and the 544-address set.The remainder constructors are in the supported ring, and `eta_bg` is dimensionless in the registry. Next I am checking whether every profile candidate for one jet is the same coefficient, and whether the stage-2 symbols carry the amplitude units.CLEAR FOR THIS BOUNDED SOURCE UNIT-TRANSPORT BUILD.

This is a source reading of the frozen packet. The worker was not run, the ring was not executed, and no pressure value or field response was computed. Runtime acceptance stays outside this verdict.

The two objects stay separate. A stored rational is a magnitude in the original base units. The physical coefficient is that magnitude together with units taken from the original expression. `worker.py` sends substitution constructors and the grade-length control to the unit engine. It does not send a bare bound numerical sum, and it does not invent a unit for the normalized denominator. `normalizedBoundDenominatorUnit` is `None` in the grade finish record. Zero fields keep `transportedUnit: None` and still carry the required unit.

## What the joins actually check

Normalization-LEFT and stage 2 use saved trees. For each face, `binding_operand` rebuilds the epsilon-stripped chemical and velocity expressions in the pinned bind order: density map, then profile equalities to stability, then one numeric substitution. `assembly` compares the strict ring normal form with the saved LEFT operand. Stage 2 is one simultaneous name substitution of `s11cc1_V_lab_held_{face}` and `s11cc1_mu_theta_lab_held_{face}` only, then the same binding. Old `bind`, source, grade, profile, and quotient functions are not called.

The placeholder units match the inherited amplitude-side units on the operands that are actually saved:

- `s11cc1_V_lab_held_plus` and `_minus` are `[1,-1,0]` in `saved/native/effective-registry.json` and in `plus-unbound-source-units-input.json`. `plus-native-velocity-units-return.json` is the same triple.
- `s11cc1_mu_theta_lab_held_plus` and `_minus` are `[-1,-2,1]`. `new-unbound-chemical-homogeneity-return.json` is the same triple.
- The unbound source return is `[1,-1,0]`.

`eta_bg`, `sigma_W`, `epsilon_shape`, `w1_profile`, and `m1_profile` are `[0,0,0]`. `L_W` is `[1,0,0]`, and the substituted magnitude is `Integer(10)`. Numeric `omega` is `Integer(3)`, the same as `frequency`. Each numeric binding stores `originalUnit` from the registry beside the magnitude. The rational is not used as a unit.

The full remainder is expression assembly, not a new Taylor proof. `full_remainder_operands` builds `full - Σ retained grades` and `numerator - denominator·retained` for `(0,0),(1,0),(0,1),(1,1)`. That is the same definition as the pinned `grade_split` fragment in `saved/transport-index.json`. In `plus-source-split.json`, the saved `fullHigherRemainder` is those four grade polynomials with minus signs, plus one negative power of the multi-term denominator. The grade-`(0,0)` piece is the same polynomial as `plus-source-jets-00.json` `source`. `quotientRingNumeratorRemainder` is the numerator polynomial minus the denominator times those grades. Excluded pure grades stay in `excludedPure`; `(0,2)` is `Integer(0)`. Per-grade `cancelled == Integer(0)` receipts, including `plus-source-quotient-00-return.json`, stay inherited.

Jet units are recomputed. For `e_W` in `plus-source-jets-00.json`, orders give `[0,0,0]`, the registry row is `[0,0,0]`, and the required coefficient is the source unit `[1,-1,0]`. The saved `jetDimension` and `requiredCoefficientDimension` match. `coefficientDimensionIndependentlyVerified` is false and is not used as the proof. `e_W_d1` and `e_W_d1d1` follow the same order rule. Consumer slot units are the inherited 3-vectors: `delta_p_plus` and `delta_p_minus` are `[-1,1,0]`; both `d_w_delta_p_*` slots are `[0,1,0]`. `THETA_BALANCE-delta_p_plus-epsilon-00-input.json` has one explicit `epsilon_shape` factor on the right, which `strip_epsilon` removes. The normalized denominator fraction for the plus source is the nonzero complex `90900+303000·I` over `Integer(1)`, with the zero imaginary denominator component allowed. No unit is attached to it.

Scale identities hold on the saved maps that were opened, with `T = tanh(composition_x/10)` and `L = Integer(10)`:

- `w1_profile` in `profile-map-1264.json` and `profile-map-1527.json` is `(1+T)/2`.
- `w1_profile_d1` in `profile-map-1360.json` and `profile-map-1601.json` is `1/2 - T²/2`.
- `m1_profile_d1` there is `-10·(1/5 - T²/5)·T/3`, which collects to `-2T/3 + 2T³/3`.
- `w1_profile_d1d1` in `profile-map-1605.json` is `(T²-1)·T`.
- `m1_profile_d1d1` there is `-2/3·(T²-1)·(3T²-1)`, which collects to `-2/3 + (8/3)T² - 2T⁴`.
- Transverse `d2`, `d3`, `d2d2`, `d3d3`, and `d2d3` maps in `profile-map-1363.json`, `1602`, `1603`, `1604`, and `1605` are `Integer(0)`.

Those are the pinned bases and the recurrence `P_r = (1-T²) dP_{r-1}/dT`. The explicit `10` in the order-1 `m` map cancels inside the ring. The argument is `composition_x/10`, not an `L` label. `profile()` and `sp.diff` are not called.

Both controls are aimed at real operands. The missing-`L` control divides a saved longitudinal derivative by `Integer(10)` to the derivative order and reuses `profile_transport`. The order-1 and order-2 maps above have nontrivial rational coefficients, so that division breaks the recurrence. The refusal text is `actual mapped profile fails native L-scaled recurrence`. The grade control sets `eta_bg` to `[1,0,0]` on the chemical constructor registry, where `eta_bg` is `[0,0,0]`. That constructor contains both `B_rho_3·epsilon_shape·eta_bg·theta·w1_profile` and `B_rho_3·epsilon_shape·theta`. The unit engine’s Add rule refuses with `inhomogeneous native Add`. Gamma rows are not re-derived. They sit after `zeta_c` as string triples, and they differ by family: `gamma_s11cb_mu_r_bg_01` is `["0","0","0"]`, while `gamma_s11cb_mu_r_bg_06` is `["2","0","0"]` and `gamma_s11cb_w_bg_04` is `["-2","-2","1"]`. Their numeric magnitudes are separate rationals such as `29/101`.

## Flags and the 544 addresses

`selectedCount`, `counts.addresses`, `counts.sourceLocations`, and `counts.consumerLocations` are not the coverage test. `addressLocations` is a uniform 10-line record from `/selected/0` (`addressId` 7956) through `/selected/543` (`addressId` 10576). Sampled pointers 0, 1, 2, 69, 169, 269, 369, 527, and 543 lie on that grid. `saved/inventory/selected.json` starts at `addressId` 7956 and ends at `addressId` 10576, and it declares `selectedCount` 544. The worker requires length 544, unique ids, the pointer set `{/selected/0..543}`, and equality of the joined id set with the selected id set. Field ids are `sha256` of the field `srepr`, checked again by `field_proof`, and marked as extra identity. `originalToProfileInputCertified: false` on consumer index entries is left false; the worker still strips epsilon and matches the profile input.

## Qualifications

These do not change the tagged units or the checked identities.

The worker checks stage-2 expression identity and the unbound walks. The placeholder registry units agree with those walks on the saved rows. There is no separate `require` that restates that agreement.

`dimensionlessValueInOriginalUnits` is a constant annotation on every numeric binding. The unit that is stored is `originalUnit`.

The physical parameter sheet lists `omega` as `"1"`. The substituted magnitude and the frequency are `Integer(3)`. Substitution uses the numeric magnitude.

`profile_polynomial` can refuse an `x/L` mismatch before the scale certificate record is written. The certificate input is already emitted, and a top-level failure record keeps the traceback. Successful assembly and successful scale matches emit both normal forms, including the coefficient residual, before the equality guard.

The ring collects sums and products and expands nonnegative integer powers. A negative power of a multi-term polynomial stays an opaque atom. Algebraically equal rationals written with different denominators can fail. They do not pass by tolerance. `0**0` is 1, which is the grade-`(0,0)` monomial.

`input/library.py` loads a same-directory sibling unit engine. The packet stores that engine at `runtime-source/S11c_d_defect_packet_native_units_lib.py`, and the worker loads the pinned measurements library by absolute path. The engine text that was read is the packet `runtime-source` copy. Its Add refusal string is `inhomogeneous native Add`.

Symbol assumptions differ on some saved profile atoms (`real=True` versus a bare symbol). That can make `same()` or a ring key fail. It does not skip the recurrence, which keys the atom name.

## Coverage limits

Read in full: `build.md`, `evidence-guide.md`, `worker.py`, `library.py`, the unit-engine portion that implements `dimension`, `saved/source-contracts.json`, `binding-context.json`, the inherited unit returns named above, the consumer unit returns, the plus-source denominator fraction, the grade-00 quotient return, and the epsilon-00 consumer input.

Read in the load-bearing slices: the chemical constructor and its `eta_bg` row, registry rows for grades, epsilon, `L_W`, `e_W` and its derivatives, profile atoms, stage-2 symbols, and the sampled gamma triples, `plus-source-split.json` through both remainder keys, `plus-source-jets-00.json` for `e_W`, the profile maps named above, `transport-index.json` at the address grid, the consumer-key region, the `grade_split` fragment, and the endpoint records of `selected.json`.

Not opened term by term: the other five remainder bodies, the other jet tables, and the rest of the 392 profile maps. Duplicate maps share the zero sha `6c01d058…` or the order-0 `w` sha `c3d2b1dc…`. Interior address pointers between the sampled indices were not read one by one. Source-pin hashes in `inputs.json` were not recomputed.

Wave, kernel, measure, `h/j`, H/J/direct kernels, flat delta, paired PV, Fourier measure, and complete summand units remain deferred. `pressureSummandUnitsComplete`, `waveKernelMeasureUnitsComplete`, and `numericalEvaluatorReady` stay false. Gamma units remain inference-dependent. Old per-grade zeros remain dependent on the accepted grade operation.