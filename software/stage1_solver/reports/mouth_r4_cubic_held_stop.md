# Held radius 4 on the approved cubic collar

Both step-166 fixed-label preparations passed independent review. Each gauge passed all 52,333 one-sided order tests after the fluid and order fields equilibrated. The primary ε = 10⁻⁸ run accepted steps 167–287 on the cubic reference fixed in `f99601aa`. All 121 crystal steps lowered energy, by 0.00992625528789 in total. Six field updates lowered it by a further 0.00002976632095. Depth was held at 0.81914036 and the exit radius at 3.30637727.

| Quantity | Primary, last accepted state |
|---|---:|
| Total energy | 92185.38355095353 |
| Physical crystal energy | 0.3073392682562248 |
| Gauge energy | 8.614831490634566 × 10⁻⁸ |
| Peak crystal force, step 166 → 287 | 0.31946134 → 0.04602958 |
| Internal force target | 0.000030 |
| Order optimality residual, fresh at 286 → fixed order at 287 | 5.45190 × 10⁻⁹ → 7.08540 × 10⁻⁵ |
| Certified lower bound on physical inverse determinant, all 28,080 cells | +4.10954101 × 10⁻⁷ |
| Inward current | 10.70771691 |
| Dense and subsonic current share | 66.75899% |
| Feed minus bulk current | +4.47642 × 10⁻¹³ |
| First fully melted cross-section, depth | 0.33635069 |

The last fresh order and fluid update was step 286. All 51,839 free zero-order nodes gained energy under a +10⁻⁹ order perturbation. The maximum absolute gauge-to-physical work ratio was 8.9486475242 × 10⁻¹⁰. The order and carrier arrays were byte-identical at 287. Only the crystal labels moved. The peak force at 287 was radial and outward at node (24, 1, 19), depth 0.02587313 and χ = 0.32715616.

On step 288, a backtrack identified a corner and the solver projected its trial displacement against the local orientation constraints. The projected direction failed the solver’s feasibility tolerance. The exception omitted the offending constraint and its residual. Cell (25, 10, 6) sets the smallest certified chart bound in the accepted state. Identifying the failing trial corner requires additional evidence. The held run stopped before a stationary six-force measurement or admission of the secondary gauge. All six shape coordinates remain held.

Independent verdict, “Scoped clear for a numerical projection stop with radius 4 held; paired equilibrium is unestablished.” The reviewer checked all 127 accepted and field-matched state hashes, the signed work, the exception and the last fresh field readings. [Stop measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_cubic_held_stop.json) · [Standing reference](mouth_r4_bicubic_reference.md)
