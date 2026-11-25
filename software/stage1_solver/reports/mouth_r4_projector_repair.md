# Radius 4 · Production corner projection through step 288

The repaired projector reproduced the feasible quarter-step from `49220309` and preserved positive orientation across the held mouth. The primary ε = 10⁻⁸ continuation uses the accepted step-288 state. Its first operation re-equilibrates the order and carrier fields at fixed labels. The ε = 10⁻¹⁰ solve starts from the primary stationary endpoint after its native residual, orientation and order-sign checks pass. The target is 0.000030. Depth is held at 0.81914036 and exit radius at 3.30637727.

| Pre-resume check | Measured result |
|---|---:|
| Row permutations matching the independent 90-digit reference | 24 / 24 |
| Largest direction component error across permutations | 1.03637 × 10⁻¹³ |
| Largest independently evaluated KKT residual norm | 5.22459 × 10⁻¹⁸ |
| Frozen energy, crystal force, field residual and current differences from the step-287 readout | 0 |
| Production step fraction and number of protected sites | ¼, five sites |
| Signed total work, production | −0.00010131370353644309 |
| Absolute work difference from diagnostic trial | 5.42 × 10⁻²⁰ |
| Peak crystal force, 287 → 288 | 0.04602957770136634 → 0.034479169454073 |
| Certified physical inverse-determinant lower bound, all 28,080 cells | +8.543232352946087 × 10⁻⁸ |

The projector normalizes each constraint and checks the returned direction, multipliers, stationarity, complementarity and descent. Eight permutations needed the extended-precision fallback. In the production replay, the three-row projection passed in binary64. The four- and five-row projections used extended precision. A new site restarts the line search along the jointly projected direction. The frozen comparison covered 93,555 physical nodal force entries, 81,125 physical free forces, 81,125 gauged free forces, 31,185 order loads and ten flow/balance readings. Recomputed total energy equaled the archived direct evaluation. Its difference from the separately accumulated trajectory energy was one binary64 spacing, 1.45519 × 10⁻¹¹. [Regression, signed work and independent readout](../../../research/pde_ledger_v3/_measurements/mouth_r4_projector_repair.json)

Independent verdict, “The 24-permutation regression, native step-288 replay and frozen step-287 invariance all pass independent review. The held primary may continue from the saved step 288. Secondary admission still requires its primary stationary endpoint.” The six stationary forces and cross-gauge comparison require paired equilibrium. [Projector](../tools/mouth_r4_corner_projection.py) · [Regression fixture and checks](../tests/test_r4_corner_projector.py)
