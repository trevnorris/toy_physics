# Radius 4 · The step-383 residual and the outer field updates

The pinned solver on `ba838359` reproduced both archived residual pairs and the forcing threshold bit for bit from step 382. All 35 Krylov curvatures were positive. The audit accepted zero steps. Both residual paths used the same 81,144 free degrees of freedom and the same Euclidean normalization. The recurrence combined products on individual search directions. The direct check applied the operator to the accumulated candidate.

| Relative residual | CG iteration 30 | CG iteration 35 |
|---|---:|---:|
| Recursive | 0.016907035796648553 | 0.01783585354604564 |
| Direct, production | 0.028488689210726767 | 0.031006079753728905 |
| Direct, promoted physical gradient and displaced coordinates | 0.02231901163588519 | 0.02559630122271862 |
| Required maximum | 0.018404342361420984 | 0.018404342361420984 |

At 30 the restart kept the candidate, replaced the residual with the direct value, reapplied the diagonal preconditioner, set the search direction to the preconditioned residual, and recomputed their inner product. At 35 the second direct check worsened and the solve stopped. Summing the same five post-restart updates in extended precision changed the recursive residual by only 4.20012 × 10⁻²² in norm.

The physical crystal product uses a central difference of two binary64 gradients. Rounding the two displaced label fields accounts for a vector difference of 2.20646 × 10⁻⁶ at iteration 35. Evaluating and differencing the gradients at the same rounded fields contributes another vector of norm 8.98887 × 10⁻⁷. These errors have a nonzero cross term. Promoting both operations reduced the five-direction additivity discrepancy from 3.17289 × 10⁻⁶ to 1.71617 × 10⁻⁹. Widening the promoted probe by up to 64× kept the candidate residual near 0.025596, above tolerance. The gauge action is analytic and an independent gauge-gradient difference agreed with it to 3.75 × 10⁻¹⁰ in product norm.

The five largest mismatch components at iteration 35 lie at depths 0.02587313 and 0.05174626 on ordered entrance nodes. Each has an incident-cell physical inverse-determinant lower bound above 0.4245. Melt-only degrees of freedom carry 96.23% of the candidate's squared displacement norm and only 6.76 × 10⁻²³ of the residual gap in norm.

All ten completed crystal–field cycles lowered total energy. From step 325 onward, the largest post-update crystal force was radial at depth 0.20698504. At 377 the order/carrier update reduced the order optimality residual from 0.00175105 to 3.09655 × 10⁻¹¹ while increasing the peak crystal force from 0.0000321265 to 0.0003298342. Of the free order nodes, 151 gained positive order and 118 returned to zero. The last five crystal steps had frozen fields. Simultaneous stationarity remains open. The six shape coordinates are held pending strategy review.

[Full-precision CG trace, operator tests, nodal bounds and cycle ledger](../../../research/pde_ledger_v3/_measurements/mouth_r4_step383_forensics.json) · [Accepted step-382 state](mouth_r4_projector_held_stop.md)
