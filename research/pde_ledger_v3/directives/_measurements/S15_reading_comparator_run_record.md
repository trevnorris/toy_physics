# S15 reading comparator for items 1–3, 5 and 6

The repaired build compared all 336 named quantities. In each engine 309 original quantities matched byte for byte and all 27 reconstructed signed sums matched exactly. For all three components of the net shell crossing, the engines agree exactly on the flux, its geometric departure and both drain-velocity contributions. The field-response departure agrees at 30 sampled points per component. [Run manifest](S15_reading_comparator_run_manifest.json) · [Summary transcript](../../scripts/out/S15_reading_comparator.out)

| Cross-engine comparison | Exact identity | Agreement at tested points |
|---|---:|---:|
| Original quantity | 303 | 33 |
| Geometric departure | 303 | 33 |
| Response departure | 279 | 57 |
| Separation derivative | 279 | 57 |
| Departure partition residual | 336 | 0 |
| Item 3 family integral | 108 | 228 |
| Item 3 emitted family result | 132 | 204 |
| Velocity contribution of drain 1 | 324 | 12 |
| Velocity contribution of drain 2 | 318 | 18 |

All 52 initial numerical flags agreed on reevaluation at 45 and 75 decimal digits. The 7,560 additional branch evaluations agreed within their empirical error estimates. Among the 144 emitted closed forms, 96 matched exactly and 48 at tested points. Among the 192 emitted general-field integrals, 36 matched exactly and 156 at tested points. The archive retains 4,812 unpaired internal diagnostic occurrences. [Detail archive](../../scripts/out/S15_reading_comparator_details.out)

All 4,704 per-object dimensional collections matched exactly. Of 406 SymPy argument and bound checks, 393 match signatures in the early printed Wolfram table, three match later Wolfram source declarations, and ten are SymPy-only bound or trigonometric checks. The numerical sample results establish agreement at the recorded points. The item-6 audit identifies the source, component and observed argument for each of the three later declarations. [Corrected audit and per-object coverage](S15_reading_comparator_coverage.json)

Of the 48 sampled emitted results with exceptional exponent planes, 24 have untested simultaneous selections. Counted separately by object, 378 feasible selections were sampled and 858 were untested. The inventory gives each object’s equations and the exact sets of equations that held at a sampled point or went untested. These counts describe coverage of the stored grid. [Coverage inventory](S15_reading_comparator_coverage.json) · [Final review disposition](S15_reading_comparator_build_review_r2/README.md)

Decision `2603dc7c`, reading engines `cc62dc9f`, original engines `3192123b`, initial comparator `e49186a6`, repaired calculation and first repair review `74d51a1f`. The scoped comparator and its reporting corrections are accepted for use in the S15 record decision list. Item 4 stays deferred. The conversion from shell momentum to throat acceleration remains open.
