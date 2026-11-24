# Radius 4 · Projected crystal motion at fixed mouth shape

**The crystal moved past the original corner. The primary solve stopped in the melted region.** At fixed depth 0.81914036 and exit radius 3.30637727, 17 active steps lowered the energy and preserved orientation. The smallest accepted fraction was 0.000244140625. The peak force magnitude was 0.07792692 before the final field update. Afterward the signed radial force at depth 0.02587313 near the entrance was **−0.07797064**, against a magnitude target of **0.000030**. [Accepted steps](../../../research/pde_ledger_v3/_measurements/mouth_r4_projected_stop.json) · [Final forces](../../../research/pde_ledger_v3/_measurements/mouth_r4_projected_fields.json)

| Held-state measurement | Value |
|---|---:|
| Inward current / dense, subsonic share of positive incident current | 10.70776412 / 66.75885% |
| Absolute inlet-to-bulk current imbalance | 1.35 × 10⁻¹³ |
| First fully melted cross-section, start → stop | 0.31047756 → 0.33635069 |
| Physical order residual with its bounds enforced / carrier residual | 4.82 × 10⁻⁹ / 1.21 × 10⁻⁸ |
| Held total / crystal energy | 92185.41380128905 / 0.32314161609 |
| Minimum sampled physical inverse volume, ordered / fully melted | 0.03571794 / 0.09863992 |

The last field update moved two cells into the fully melted set and one out. With active labels fixed, 105 inner steps reduced the unscaled melted-label residual from 0.00033830566 to **5.93109911 × 10⁻⁸**, above the **3.0 × 10⁻⁸** target. The line search could no longer resolve sufficient energy decrease. The unsuccessful inner extension was saved separately from the retained state. [Stop record](../../../research/pde_ledger_v3/_measurements/mouth_r4_projected_stop.json)

The original corner’s physical inverse volume recovered from 0.00050144 to 0.07878394. The largest inner residual was on radial label (25, 11, 17), depth 0.28460443, in fully melted material. Both maps passed the whole-cell orientation test with chart lower bound 0.00139076. [Orientation measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_projected_orientation.json)

Raising order by 10⁻⁹ at node (24, 9, 17) raised physical energy by 1.17249 × 10⁻¹¹ but removed two cells from the numerical gauge, lowering total energy by 1.08318 × 10⁻¹⁰. This discontinuity rules out a local minimum of the gauged energy. The inner line search held order fixed, so its failure needs a separate explanation. [One-sided work](../../../research/pde_ledger_v3/_measurements/mouth_r4_projected_fields.json)

Only the primary gauge ran. The six stationary shape forces and cross-gauge differences remain unmeasured. Next, compare the melted-label gradient with high-precision signed energy changes on the saved inner iterate.

The solver and feasibility tests are recorded at `d72f1d27`, following `4c3c0a59`. Independent review concluded, “Scoped clear for the 17 fixed-shape primary descents, the changed-mask inner-solve stop, the retained-state one-sided order descent and the separate orientation measurements on both saved maps. Paired equilibrium and stationary shape forces remain unresolved.” [Input tests](mouth_r4_corner_feasibility.md)
