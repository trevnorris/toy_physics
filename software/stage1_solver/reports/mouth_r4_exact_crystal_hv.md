# Radius 4 at step 382: exact crystal Hessian action

The archived step-383 Newton–CG solve stopped because its finite-difference Hessian action depended on the direction and step width. The binary64 forward tangent follows the same physical gradient through the cubic reference map, the moving radial weld cut, both quadrature segments and the cut-work adjoint. It uses the existing analytic 39-point gauge action. The reference state comes from `ba838359` and the frozen Krylov vectors from the audit recorded at `578a5526`. The six shape coordinates were held. The replay accepted zero steps.

| Step-382 check | Measurement |
|---|---:|
| Relative additivity error, sum / difference | 6.649273094452398e-16 / 2.3556302400669105e-14 |
| Normalized symmetry error | 2.088056361145881e-18 |
| Largest relative error against binary128 on the three archived directions, free coordinates | 1.7560865314268252e-14 |
| Largest relative error against binary128 on four moving-weld perturbations, free coordinates | 2.7898154237107134e-15 |
| CG iterations / restarts | 29 / 0 |
| Recursive / direct relative residual | 0.01711870515041706 / 0.01711870515041708 |
| Norm of the difference between residual vectors, divided by the norm of the right-hand side | 1.4319917848364876e-15 |
| Eisenstat–Walker tolerance | 0.018404342361420984 |
| Median joint Hessian product, new / archived | 0.6022583693265915 s / 1.0687327338171992 s |

Frozen-state differences were exactly zero on all energies, all 93,555 nodal physical forces, both 81,144-entry packed force banks, 31,185 order loads, all field residuals and ten flow and balance readings. The direct total energy was 92185.38326662472.

The binary128 oracle recomputed the first variation of the archived assembly at two symmetric probe widths. Each localized perturbation moved nine quadrature-ray intersections in its witness cell. The widths agreed within 8.97e-17 on the free coordinates. The replay timings cover 30 products, excluding the separate 24-product preconditioner assembly. The assembled diagonal's historical display field retained the default finite-difference width during the replay. The exact calls had no probe width. The display now reports zero.

Independent review concluded “Phase 1 passes all four zero-step gates on saved R4 step 382.” The outstanding physical question is whether the labels and order field can reach simultaneous stationarity at this held shape. [Full measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_exact_hv382.json)
