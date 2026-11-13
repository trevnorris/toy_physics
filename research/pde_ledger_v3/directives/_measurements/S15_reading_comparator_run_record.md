# S15 reading comparator, items 1–3, 5 and 6

Builder result, 2026-11-12. All 336 required quantities were present in both accepted streams. For each engine, 309 reading quantities matched their original payloads byte for byte, and all 27 signed surface sums matched the independently assembled original terms exactly. Across engines, the net momentum crossing the two surfaces, its geometric departure and both drain-velocity contributions agree exactly in all three Cartesian components. The change due to the fields’ response agrees at the sampled points in all three components. [Run manifest and quoted transcript](S15_reading_comparator_run_manifest.json)

| Cross-engine comparison (336 quantities per row) | Exact identity | Agreement at tested points | Incomplete |
|---|---:|---:|---:|
| Original quantity Q | 303 | 30 | 3 |
| Geometric departure | 303 | 30 | 3 |
| Response departure | 279 | 54 | 3 |
| Separation derivative | 279 | 54 | 3 |
| Residual of the departure partition | 336 | 0 | 0 |
| Item 3 family integral before evaluation | 108 | 211 | 17 |
| Item 3 emitted family result | 132 | 204 | 0 |
| Velocity contribution of drain 1 | 324 | 6 | 6 |
| Velocity contribution of drain 2 | 318 | 8 | 10 |

**Differences and limits.** No arithmetic disagreement persisted at a checked point. Of the 354 initial numerical flags, 353 closed-form samples agreed on reevaluation at 45 and 75 decimal digits. The remaining velocity residual was 128 in binary64. Setting the target-centre velocity and both radius rates to zero exactly gave zero in both native operands, with higher derivatives retained. Incomplete rows have one or more quadrature cases with insufficient convergence. The three incomplete Q and departure rows are the three components of the drain-2-only regional momentum rate. The transcript names every other incomplete object. [Counts, cases and residuals](S15_reading_comparator_run_manifest.json)

**Item 3.** Of 144 emitted closed forms, 96 matched exactly and 48 at tested points. Of 192 emitted general-field integrals, 36 matched exactly and 156 at tested points. The elementary exponent boundaries agreed across all 336 inventories. For six mixed-drain regional quantities, Wolfram explicitly printed six simultaneous conditions alongside eight elementary boundaries. All 7,560 additional evaluations across the family integral, closed operand and emitted result agreed within their error estimates. Other intersections and the full continuous domain remain unproved for sampled results. The 1,044 Wolfram-only divisor-support tests across 288 native records returned False on their respective supports. The per-occurrence archive also identifies 4,816 internal kernel records without a cross-engine pairing. [Branch and kernel checks](S15_reading_comparator_run_manifest.json)

**Dimensions and conventions.** All 4,704 per-object dimension collections matched exactly. The three geometric and 72 response departures with different literal-zero spellings also compared exactly. The existence conditions are phrased differently. Both require regular fields and responses, regular centre and radius motion, domination for differentiation under the integral, and pressure regularity on the densities used. SymPy printed 406 observed-argument or bound checks. Of these, 393 matched Wolfram’s declared argument slots, three source-velocity response slots lacked a matching printed signature, and ten were SymPy-only bound or trigonometric checks. Two zero lower angular bounds returned an empty observed dimension set. [Typed comparisons](S15_reading_comparator_run_manifest.json)

**Reproduction and next step.** Decision list `2603dc7c`, reading engines `cc62dc9f`, original engines `3192123b`, original comparator `5a019394`. The manifest records commands, fingerprints, guard outcomes and quoted results. Use fresh guard and output directories for a reproduction. The local [summary transcript](../../scripts/out/S15_reading_comparator.out) gives the final outcomes under `final_field`. The [detail archive](../../scripts/out/S15_reading_comparator_details.out) is a tar archive containing native operands, object records, precision and branch checks. Independent build review is next, before the S15 step record. The builder changes await that review and have no commit ID. Item 4 remains deferred. The mapping from shell momentum to throat acceleration remains open. [Reviewed comparator remit](../S15_reading_comparator_decision_list.md)
