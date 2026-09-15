# S11c-d finite-action momentum quadrature

The complete run finished in 4 h 2 min with exit zero and empty stderr. All 80
native factorized integrals, 70 bound sources, six nested profiles and three
Abel transfer pairs retain their source and ordered-limit joins. Both old grids
agree with the native evaluator: all 300 action components and 1,920 nonlocal
terms pass, with maximum scaled action difference 1.11e-16 and term difference
3.60e-17. Actual momentum-weight mutations respond in all evaluated layouts.

Three concentration-aware momentum refinements use panel orders 8/12/16 and
outer orders 24/40/64, with source/profile orders fixed at 128, momentum bounds
+/-2, source bounds +/-32, profile bounds +/-10 and regulator 0.2. Full action
changes fall from 1.51e-3 to 3.64e-5 for the first field and from 4.04e-4 to
1.05e-5 for the second. These are finite-grid differences, not an established
full-action limit. More momentum resolution is needed before tail/regulator work.

All 5,398 tags, 2,697 fresh write keys and 22,822 metadata paths replay. All
43 frozen sources, five main artifacts and 1,326 saved row/group/grid/partial
packets have verified hashes. Bound and numerical packets remain byte-identical
before/after emission. Peak process RSS was 235.4 MiB; estimated phase workspace
was 6.00 MiB and batch cache 89 KiB, within their separate 32 MiB budgets.

The 2,215,633-byte transcript is prepared for DataLad/git-annex publication at
`scripts/out/S11c_d_momentum_action.out`; SHA256
`d910ed2b3d47397c950f4965024850c8ac9bf4ac5f35913e75c3dc16039a9900`.
All intermediate data remain in repository scratch. Next: isolate the remaining
momentum changes and refine them, reusing accepted operands. Infinite-domain
tails/interchange, Abel weak limits, boundary matching, scattering and bound
poles remain work. Approved inputs and the solver/export contract are unchanged.
