# Held radius 4 · Conjugate-gradient stop at step 383

The ε = 10⁻⁸ primary accepted 94 crystal updates after the reviewed step 288 in `bef7526a`. Each lowered energy and preserved positive orientation. At attempted step 383, the directly recomputed Newton linear-system residual rose after a CG restart. The solver stopped before projecting a direction or trying a displacement. Elapsed time was 2 h 1 min 5 s. The held solver had no outer-step or elapsed-time cap.

| Relative CG residual at step 383 | Iteration 30 | Iteration 35 |
|---|---:|---:|
| Recursive estimate | 0.01690704 | 0.01783585 |
| Direct evaluation of the proposed direction | 0.02848869 | 0.03100608 |
| Required upper bound | 0.01840434 | 0.01840434 |

The largest signed force at 382 was radial, −0.00008798533 at node (25, 8, 14), radius 3.88592755 and depth 0.20698504. Its χ = 0 vertex affects adjoining ordered material. Across the run, 456 projected directions passed the descent and KKT checks. The sole extended-precision fallback at step 377 rejected a Newton direction. The projected Jacobi alternative lowered energy by 6.31634 × 10⁻⁹ at full step.

| Primary state at accepted step 382 | Read-only measurement |
|---|---:|
| Total / physical crystal energy | 92185.38326662472 / 0.3071402521934724 |
| Gauge energy / ratio to physical crystal energy | 9.62329302 × 10⁻⁸ / 3.13319174 × 10⁻⁷ |
| Peak crystal force / gauge-scaled melted-label residual | 0.00008798533 / 1.57114056 × 10⁻⁷ |
| Order optimality residual, fresh at 377 → frozen at 382 | 3.09654524 × 10⁻¹¹ → 0.00214183487 |
| Carrier amplitude / phase residual | 2.17239581 × 10⁻⁹ / 8.97927518 × 10⁻¹¹ |
| Inward current / dense and subsonic share | 10.707716628605539 / 66.75898947% |
| Absolute feed-to-bulk current discrepancy | 4.08562073 × 10⁻¹⁴ |
| Certified physical inverse determinant, all 28,080 cells | ≥ 4.46074987 × 10⁻⁷, limiting cell (25, 10, 9) |

The crystal residual fell below 0.000030 at steps 325, 358 and 362 with order and carrier fields frozen. Re-equilibrating those fields brought the crystal residual back to 0.00048973, 0.00024684 and 0.00031951, respectively. The last fresh field update was step 377. All 51,504 free zero-order nodes gained total energy under δχ = +10⁻⁹, with maximum absolute gauge-to-physical work ratio 1.09648 × 10⁻¹⁰. The order and carrier arrays stayed byte-identical through step 382. The first fully melted cross-section moved from depth 0.33635069 at step 288 to 0.36222382 at step 325 and held there through the stop.

The ε = 10⁻¹⁰ run required a stationary primary and was never admitted. The six stationary force components and cross-gauge differences remain unmeasured. Depth 0.81914036 and exit radius 3.30637727 were held throughout. Free release remains gated. The open numerical question is why the directly evaluated Hessian action and the CG recurrence disagree at step 383.

Independent readback verdict, “Step 382 is positively oriented and draining. Crystal force and constrained order residual both exceed target. The primary stopped inside CG at attempted step 383; the secondary gauge and stationary shape forces remain unmeasured.” [Measurements and residual history](../../../research/pde_ledger_v3/_measurements/mouth_r4_projector_held_stop.json) · [Step-288 repair](mouth_r4_projector_repair.md)
