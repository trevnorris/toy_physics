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

The test uses joint Newton–CG for the active and melted crystal labels, with order and carriers held through a burst. The frozen physical energy, smooth gauge and welded geometry are those of `936ceca4`, `1847f480` and `393aec91`. Advancement requires an independent Hessian check and a measured decrease in time and residual on a five-step benchmark from this exact state. The paired stationary forces remain open.
