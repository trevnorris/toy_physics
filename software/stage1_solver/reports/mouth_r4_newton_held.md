# Radius 4 · Held crystal optimizer

The guarded primary run at ε = 10⁻⁸ was paused at accepted step 29 of `393aec91`. The saved fields and crystal labels have SHA-256 `d2e764212c08bdcacfbc18121812ceb484f2eaffcd0b44f46991a9cf03120202`. The mouth is held at depth 0.81914036 and exit radius 3.30637727. At the pause the smallest sampled inverse volume determinant was 0.09535070, the intake was 10.70774948, and 66.7589% of inward current was dense and subsonic. The first fully melted section was at depth 0.33635069. These are held-state measurements.

| Original run | Steps 1–21 | Steps 22–29 |
|---|---:|---:|
| Wall time (min/step) | 1.9812 | 13.3991 |
| Inner-iteration journals on accepted trials (s, share) | 1,443.4602 · 57.82% | 6,153.7220 · 95.68% |
| Inner-iteration journals at field rematches (s) | 36.3386 | 0.0014 |
| Completed inner solves on rejected trials | 0 | 0 |
| Combined time outside those journals (s) | 1,016.4633 | 277.8572 |
| Field rematches / minimum outer backtracks | 6 / 145 | 2 / 66 |

The first window includes startup. The journal intervals give lower bounds for time inside the inner solver. In the late window, field rematches took at most 244.6351 s and projection plus early backtracking at most 22.0253 s. The original code emitted no separate timers for crystal gradients, projected directions, and order and carrier updates. The unfinished trial after step 29 is excluded. The complete per-step energy, maximum and Euclidean force norms, backtrack counts or bounds, and interval timestamps are in the [29-step measurement](../../../research/pde_ledger_v3/_measurements/mouth_r4_newton_history.json).

From step 21 to 29, energy fell by 0.0006312051 while the maximum free crystal force rose from 0.0749740 to 0.2970619 and its Euclidean norm rose from 0.5630454 to 0.7607605. The five largest nodal force magnitudes on the preserved state were:

| Label node | Radius | Depth | Polar angle (rad) | Order χ | Force magnitude |
|---|---:|---:|---:|---:|---:|
| 23,3,18 | 3.90936 | 0.07762 | 1.36253 | 0.111951 | 0.297283 |
| 23,2,18 | 3.92339 | 0.05175 | 1.36253 | 0.196127 | 0.293038 |
| 24,3,17 | 3.93875 | 0.07762 | 1.26135 | 0.112680 | 0.198383 |
| 24,2,17 | 3.95280 | 0.05175 | 1.26135 | 0.197283 | 0.194893 |
| 23,3,15 | 3.90971 | 0.07762 | 1.06525 | 0.113293 | 0.108814 |

<details>
<summary>All 29 accepted steps, including timings and residuals</summary>

Wall time runs to the acceptance marker. The two inner columns give recorded iteration time. The energy and force columns include the following field rematch where one was run. Euclidean norms use the independent active crystal coordinates.

| Step | Wall (s) | Inner on accepted trial (s) | Inner at field rematch (s) | Backtracks | Total ΔE | Max force | Euclidean force norm |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 229.19 | 21.869 | 12.640 | ≥4 | -0.00318026 | 0.089194 | 0.714910 |
| 2 | 27.69 | 24.509 | 0.000 | ≥5 | -0.00110299 | 0.068025 | 0.654809 |
| 3 | 23.94 | 20.659 | 0.000 | ≥5 | -0.00081232 | 0.148683 | 0.699538 |
| 4 | 34.67 | 31.407 | 0.000 | ≥5 | -0.00089030 | 0.079494 | 0.647657 |
| 5 | 201.57 | 26.639 | 23.695 | ≥6 | -0.00037432 | 0.150935 | 0.660061 |
| 6 | 33.28 | 30.020 | 0.000 | ≥6 | -0.00037260 | 0.069495 | 0.606932 |
| 7 | 41.56 | 38.378 | 0.000 | ≥6 | -0.00036064 | 0.063971 | 0.590132 |
| 8 | 51.00 | 47.790 | 0.000 | ≥6 | -0.00034310 | 0.117009 | 0.608019 |
| 9 | 202.60 | 40.812 | 0.001 | ≥7 | -0.00016397 | 0.066029 | 0.589830 |
| 10 | 53.02 | 49.713 | 0.000 | ≥7 | -0.00014622 | 0.463806 | 1.030973 |
| 11 | 60.31 | 57.009 | 0.000 | ≥7 | -0.00017135 | 0.819192 | 1.342607 |
| 12 | 84.59 | 81.225 | 0.000 | ≥7 | -0.00015276 | 0.077818 | 0.582524 |
| 13 | 221.73 | 53.312 | 0.002 | ≥7 | -0.00014828 | 0.170621 | 0.657499 |
| 14 | 126.65 | 123.434 | 0.000 | ≥8 | -0.00008645 | 0.061083 | 0.561494 |
| 15 | 90.25 | 87.009 | 0.000 | ≥8 | -0.00006227 | 1.245798 | 2.118056 |
| 16 | 122.42 | 119.151 | 0.000 | ≥8 | -0.00008487 | 0.219381 | 0.638767 |
| 17 | 255.74 | 113.642 | 0.000 | ≥8 | -0.00007291 | 0.226260 | 0.684300 |
| 18 | 98.36 | 94.928 | 0.000 | ≥8 | -0.00007559 | 0.086313 | 0.569930 |
| 19 | 219.00 | 215.718 | 0.000 | ≥9 | -0.00003833 | 0.612220 | 1.028475 |
| 20 | 74.72 | 71.459 | 0.000 | ≥9 | -0.00003180 | 0.549874 | 1.007729 |
| 21 | 243.97 | 94.778 | 0.001 | ≥9 | -0.00004061 | 0.074974 | 0.563045 |
| 22 | 565.85 | 562.576 | 0.000 | ≥9 | -0.00003713 | 0.059867 | 0.550989 |
| 23 | 719.79 | 716.458 | 0.000 | ≥9 | -0.00003437 | 1.195554 | 1.792252 |
| 24 | 294.80 | 291.404 | 0.000 | ≥9 | -0.00003587 | 0.861176 | 1.355622 |
| 25 | 504.57 | 373.526 | 0.000 | ≥9 | -0.00003650 | 0.083279 | 0.562729 |
| 26 | 2072.47 | 2069.213 | 0.000 | ≥9 | -0.00003509 | 0.149874 | 0.600333 |
| 27 | 152.01 | 148.804 | 0.000 | ≥8 | -0.00006335 | 0.626605 | 1.109961 |
| 28 | 293.35 | 290.025 | 0.000 | ≥7 | -0.00008652 | 0.157859 | 0.611176 |
| 29 | 1828.73 | 1701.716 | 0.001 | 6 | -0.00030237 | 0.297062 | 0.760761 |

</details>

The step-29 evaluator reproduced total energy 92185.40445819321, physical crystal energy 0.318747421489235, gauge energy 4.9558578274692115 × 10⁻⁸, active forces, all field residuals and every current reading with absolute difference zero. The largest relative errors in the Hessian-vector products were 1.20 × 10⁻⁹ in ordered material, 1.11 × 10⁻⁹ at the melting front and 7.36 × 10⁻¹³ in melted material. The largest relative symmetry error was 4.34 × 10⁻⁹. Independent review concluded, “Scoped clear for the five-step benchmark from saved step 29.” [Hessian and invariance measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_newton_audit.json)

The benchmark uses joint Newton–CG for the active and melted labels, the measured packed Hessian diagonal, and the physical energy of `936ceca4`, `1847f480` and `393aec91`. An order/carrier update is due after a fivefold drop in the peak force, at its target, or after 20 accepted steps. For this five-step comparison we specified a 20% improvement in both mean time per step and final peak force, a lower final Euclidean norm, and a last-three median force below the initial 0.2970619. Timing includes preconditioner setup inside the steps. Initial held-geometry setup is reported separately. Paired stationary forces remain open.
