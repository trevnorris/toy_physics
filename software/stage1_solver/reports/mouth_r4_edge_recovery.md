# Radius 4 · Edge feasibility and frozen order at step 108

**A feasible edge step exists. Three order nodes blocked the paired restart.** Step 108 of `1f46c4ad` holds radius 4 at depth 0.81914036 and exit radius 3.30637727.

All eight vertices of cell (25, 11, 10) touch only melted cells. Four are held and four move freely. At depth 0.28460443, local point (0, 0, 0.399257224707046), the chart minimum is 2.39601188 × 10⁻⁶ and the physical inverse determinant is 3.46480550 × 10⁻⁶. The largest free chart derivative is −361.33790486 on the radial label at (25, 11, 10). [All 24 derivatives and signed controls](../../../research/pde_ledger_v3/_measurements/mouth_r4_edge_descent.json)

Moving the free melted labels doubled the edge margin with zero crystal work. Adding an entrance displacement lowered total energy by 3.75506093 × 10⁻⁶ at full step. The signed entrance force moved from −0.09873119 to −0.06499535. The peak magnitude moved to (24, 2, 14) and fell to 0.08680906. The whole-chart lower bound stayed positive at 1.17476133 × 10⁻⁶.

Raising one zero-order coefficient by 10⁻⁹ gave the following physical work. A 65-digit calculation traced all three reversals to a reduced crystal contribution at fixed bath fields. [All order probes](../../../research/pde_ledger_v3/_measurements/mouth_r4_edge_order_work.json) · [Independent energy components](../../../research/pde_ledger_v3/_measurements/mouth_r4_edge_order_review.json)

| Joined order node | Step 100, immediately after the field update | Step 108 |
|---:|---:|---:|
| 36353 | +9.36471 × 10⁻¹³ | −1.37832 × 10⁻¹² |
| 45532 | +4.18493 × 10⁻¹⁵ | −1.75634 × 10⁻¹⁷ |
| 62811 | +1.80694 × 10⁻¹⁴ | −3.93708 × 10⁻¹⁴ |

Peak crystal force rebounded from 0.01756028 to 0.57995197 **before** the step-94 field update. That update moved it to 0.58105854. The eight crystal steps after the step-100 field update raised the order residual from 3.78 × 10⁻⁹ to 0.00137832. [Force and field-update history](../../../research/pde_ledger_v3/_measurements/mouth_r4_edge_static.json)

The frozen rule uses 27 Lobatto points and one moving point per edge. It measured the edge dip at 2.39602878 × 10⁻⁶. The cell's unit-barrier slope for opening the edge changed from −27.4430 to −3886.8663. Faces and interiors still require the whole-cell certificate. [Coverage and restoring derivative](../../../research/pde_ledger_v3/_measurements/mouth_r4_aligned_stencil.json)

Physical energy, crystal forces and current readouts changed by zero. At ε = 10⁻⁸ the gauge energy changed by −1.44548 × 10⁻⁹. The signed shifts for depth, exit, mean width, depth asymmetry, waist and polar deformation were (−5.4404, +3.1135, +9.8348, −3.1037, +8.7525, −1.5478) × 10⁻¹⁰. All shifts were 100 times smaller at ε = 10⁻¹⁰. Among the 52,265 positive physical increments, the maximum gauge-to-physical work ratios were 6.92827 × 10⁻¹¹ and 6.92827 × 10⁻¹³. The same three total increments stayed negative at both strengths. [Frozen invariance, force shifts and signed gate](../../../research/pde_ledger_v3/_measurements/mouth_r4_edge_static.json)

Independent review concluded, “Scoped clear for the aligned-stencil order sweep and six frozen shape-force shifts at step 108. The strict Part 3 order-sign gate fails at both strengths.” The full derivative and Hessian audit, recovery burst and paired relaxation remain gated. The next discriminating test would hold the step-108 labels fixed, equilibrate the order and carrier fields, and repeat the all-node sweep before another crystal step.
