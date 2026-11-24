# Radius 4 · a material cell regained order during held relaxation

The held 24 × 24 run at ε = 10⁻⁸ stopped after **11 accepted internal crystal steps**. The largest free-label force was **0.17798342**, against the **0.000030** target. The last accepted fraction was 1/4096. At the 0.0001 floor, both search directions proposed turning part of the same material cell inside out. The saved state kept positive orientation. The 10⁻¹⁰ held run had not been admitted when the primary stopped, so there is no paired equilibrium, stationary six-force table or checkpoint 11. The starting-state derivative checks are recorded in `091e8c4d`. [Derivative controls](mouth_r4_slaved_work.md) · [Endpoint, step history and signed trials](../../../research/pde_ledger_v3/_measurements/mouth_r4_barrier_stop.json)

| Measurement on the last accepted state | Value |
|---|---:|
| Held depth / exit radius | 0.81914036 / 3.30637727 |
| Inward current / dense and subsonic share, including swirl | 10.70788500 / 66.75849% |
| Feed-to-bulk current discrepancy | 1.08 × 10⁻¹³ |
| First sampled fully melted cross-section, before → after the order update | 0.28460443 → 0.31047756 |
| Energy change over the 11 steps at fixed fields | −0.0289912554 |
| Total / crystal energy | 92185.44213692489 / 0.3048567027 |
| Gauge energy / crystal energy | 4.28935 × 10⁻⁸ |

The largest remaining force was **+0.17798342 on a radial reference label near the entrance**, at depth 0.02587313. The order field had equilibrated before the eleven crystal steps. With the final crystal labels its stationarity residual rose from 1.83 × 10⁻⁹ to **0.01771717**. The carrier gradient stayed exactly unchanged and the recomputed currents matched the saved post-rematch readings. [Final fields, signed label forces and independent current check](../../../research/pde_ledger_v3/_measurements/mouth_r4_barrier_held_fields.json)

The order update took 759 cells out of the fully melted set. In cell (15, 9, 1), all eight vertices had zero order before that update. One vertex regained order, reaching χ = 0.0005496956. At the diagonally opposite corner, at radius 3.30532109 and depth 0.25873130, χ stayed zero. On the saved state that corner had physical inverse determinant 0.0005014362, label-gradient norm 0.8904783 and formal elastic energy density 1.49228 × 10⁷. Its weighted crystal energy density was zero. [Order at the eight vertices](../../../research/pde_ledger_v3/_measurements/mouth_r4_held_order_cell_history.json)

At the smallest allowed trial, the corner's **chart determinant** was −0.00502231 with the accumulated search direction and −0.00504954 after restarting it. Both choices of how to extend the melted labels gave the same values. Each of the four nodal labels determining this corner belongs to at least one cell with positive order. A melt-only solve therefore cannot change this determinant at fixed active labels. Order is positive arbitrarily close to the corner inside this cell, so the proposed inversion would have reached ordered material. [Independent corner reconstruction](../../../research/pde_ledger_v3/_measurements/mouth_r4_held_corner_support.json)

The next discriminating calculation is whether an active-label displacement constrained by that corner's orientation can still lower energy at the held geometry. Physical mouth collapse and free-fall response remain open. Independent review concluded, “The retained state remains positively oriented and draining. Both tested search directions would invert positive-order material near a zero-order corner, and a melt-only extension cannot alter that corner determinant at fixed active trace.”
