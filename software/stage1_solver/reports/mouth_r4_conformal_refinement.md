# A conforming split preserves the radius-four crystal and reaches the same positivity barrier

Starting from the feasible material state recorded at `bd44031a`, we split the polar edge of triangle 0 in prism (5, 18, 0) at an existing integration boundary, 0.12129504969343591 radians. Six incident prisms became twelve children. Eighteen neighboring transition cells matched their faces and edges. The three components of the inserted site's displacement were free. The entrance stayed held at 4 with one quantum of circulation, the same order energy, phase drop 0.4510084282857874 and carrier potential \(U(n)=(n^5-5n+4)/20\).

| Before descent, refined minus original | Maximum absolute difference |
|---|---:|
| Reference labels, placement and their first derivatives | 0 exactly |
| Total energy and whole-crystal elastic energy | 0, 0 |
| Elastic energy in the target and each of 35 surrounding prisms | 0 |
| All 10,074 original generalized forces, including each checked regional contribution | 0 |
| Crystal load on spatial order, full field gradient, eight current/balance readouts | 0, 0, 0 |

The target prism retained elastic energy 0.0000011219810152387763. The new site's three generalized forces were +0.0013280451, +0.0001196192 and +0.0000113272. Signed work on a separate smooth displacement checked all three derivatives to within 1.05 × 10⁻¹⁰. The guard certified all twelve child supports and rejected an imposed χ = 0.001 at a folded child corner with chart determinant −0.18983935. [Measurement and regional comparisons](../../../research/pde_ledger_v3/_measurements/mouth_r4_conformal_refinement.json)

**First joint material step.** With the carrier and spatial order held fixed during the crystal move, the original line search tested the coarse and new coordinates together. It rejected 34 positive-order overlaps. The first accepted fraction was 2⁻³⁴ = 5.820766091346741 × 10⁻¹¹, against the required 10⁻⁴. Crystal energy decreased by 3.7914116290949096 × 10⁻¹⁴. The two assembled fixed-field totals rounded to the same value. Doubling the accepted step gave a negative chart determinant at positive order in the upper child of the target prism. All 6,918 prisms passed the guard on the accepted state.

| Point in the upper child on the accepted state | χ | Physical J | W | χW |
|---|---:|---:|---:|---:|
| Location of the last rejected proposal's overlap | 7.00946 × 10⁻¹⁹ | 2.87323 × 10⁻¹⁴ | 30.78497 | 2.15786 × 10⁻¹⁷ |
| Tightest verified volume bound, near depth 0.20214 | 8.29976 × 10⁻²¹ | 2.11022 × 10⁻¹⁴ | 31.09361 | 2.58070 × 10⁻¹⁹ |

The refined region's positive-order volume infimum lies in [0, 2.1102191078959685 × 10⁻¹⁴]. The upper endpoint comes from a pointwise physical determinant. The order and weighted energy density each approach zero at a regular boundary point. Small χW is representable in float64. The stall diagnosis rests on the strictly rejected larger step, the measured energy decrease and the independent physical volume readback.

**Decision.** The specified stop gate fired before checkpoint 11. The fixed-field state had exit radius 3.27757267, depth 0.82008883, and inward feed current 10.57354576. The dense, subsonic share of incident inward disk current was 66.1109%. The sampled fully melted disk began at depth 0.28460443. The last six signed shape forces still belong to checkpoint 10. The physical fate of the held throat and its readiness for free release remain open.

Independent review concluded “The one-pass conforming refinement preserved the saved R4 state and reproduced the positive-order step barrier. The specified Eulerian pivot gate is met before checkpoint 11.” The next test is to represent material reference coordinates on the spatial grid and check the same ordered strain, elastic work and currents on this saved state before evolving it. Radius 3 remains blind. [Exact child support](../tools/mouth_r4_conforming_prisms.py) · [Joint displacement test](../tools/mouth_r4_refinement_descent.py) · [Independent constitutive readback](../tools/mouth_r4_refinement_review.py)
