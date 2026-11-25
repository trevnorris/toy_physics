# Radius 4 · Held Newton–CG stop at step 108

The primary run at ε = 10⁻⁸ accepted 74 steps beyond the reviewed step 34 of `e60e9b01`. All 74 lowered the fixed-field energy. Four order and carrier updates brought the total energy change to −0.00958300642. The smallest peak active crystal force was 0.01756028 at step 90. It rose to 0.09873119 by step 108, against a target of 0.000030. The shape stayed fixed at depth 0.81914036 and exit radius 3.30637727. [Signed work, residuals and trial fractions](../../../research/pde_ledger_v3/_measurements/mouth_r4_newton_stop_work.json)

| Held quantity | Step 34 | Step 108 |
|---|---:|---:|
| Total energy | 92185.40038481 | 92185.39080181 |
| Crystal energy | 0.31718485 | 0.31204470 |
| Numerical gauge energy | 5.36264 × 10⁻⁸ | 7.49996 × 10⁻⁸ |
| Peak active crystal force | 0.04750481 | 0.09873119 |
| Euclidean force norm over all free labels | 0.42932219 | 0.31018037 |
| Minimum sampled physical inverse determinant, whole grid | 0.01863644 | 0.00231393 |
| Minimum sampled physical inverse determinant, ordered material | 0.35957709 | 0.43358712 |

The largest remaining force was −0.09873119 on the axial label component of node (24, 2, 10), at radius 3.95280734, depth 0.05174626 and χ = 0.20550515. Eight crystal steps followed the last order and carrier update. The order optimality residual rose from 3.78 × 10⁻⁹ at step 100 to 0.00137832 at step 108. The saved mouth carried 10.70772910 inward current, 66.758962% in dense subsonic material including swirl. Feed and bulk currents agreed within 2.62 × 10⁻¹³. The first fully melted cross-section was at depth 0.33635069. [Frozen field and current measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_newton_held_stop.json)

**The limiting point lies inside an edge of a fully melted cell.** Both the projected Newton direction and its Jacobi fallback exhausted 14 step fractions each. At the smallest eligible fraction, 1/8192, both would invert cell (25, 11, 10) at local coordinates (0, 0, 0.3984375), radius 3.85106331 and depth 0.28460443. The accepted chart determinant there was +2.43628 × 10⁻⁶. The two trial values were −2.68242 × 10⁻⁶ and −2.86841 × 10⁻⁷, well beyond the 9.29 × 10⁻¹⁴ roundoff margin. All eight order coefficients were zero. The determinant barrier samples the vertices. On the accepted cell the smallest vertex value was 0.00728637, over 3,000 times the minimum of +2.39601 × 10⁻⁶ along the limiting edge. The whole-cell orientation check caught the trial inversions. [Trial determinants](../../../research/pde_ledger_v3/_measurements/mouth_r4_newton_stop_work.json) · [Independent edge calculation](../../../research/pde_ledger_v3/_measurements/mouth_r4_newton_edge.json)

All six shape coordinates remain held. The second gauge, six stationary shape forces and paired one-sided order test require an internally stationary primary. The next question is whether a direction that preserves orientation along the limiting edge can relieve the entrance load. The survival and free fall of a released mouth remain open.

Independent review concluded, “Scoped clear for the fixed-shape primary stop at step 108. The accepted chart is positively oriented and both allowed line searches encounter resolved inversions on an edge of a fully melted cell. Crystal and order remain above target. The record supports neither a paired stationary readout nor a physical collapse.”
