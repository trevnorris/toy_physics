# Radius-four crystal relaxation reached its hard stop

The held mouth remained open and draining, but both crystal update orders reached the negligible-descent threshold at the boundary of weakly ordered material. Checkpoint 10 is the last completed whole-mouth state. The run saved a feasible diagnostic state before the six force probes for checkpoint 11. An independent review verified the material events and the saved state. [Measurements, signed steps and review](../../../research/pde_ledger_v3/_measurements/mouth_r4_feasible_stall.json) · [Solver and starting-state controls at commit 368ab67f](mouth_r4_feasible_crystal.md)

## What stopped the crystal

With order fixed, 29 accepted displacements lowered total energy by 0.0293338420. The last two step fractions were 3.72529 × 10⁻⁹ and 2.32831 × 10⁻¹⁰. Rematching the fields then lowered energy by another 0.00430507271. Order at node 51603 rose from 0.0000364030 to 0.000946277 and reached 0.000947306 during interleaving.

Three interleaved crystal displacements passed the exact-support and strict energy-decrease tests. The final two accepted fractions were 3.72529 × 10⁻⁹ and 2.91038 × 10⁻¹¹. They lowered energy by 2.40230 × 10⁻¹² and 1.01780 × 10⁻¹³, with zero additional resolved energy change on either field update. Both drops were below the predeclared 4.65661 × 10⁻¹⁰ threshold. This triggered the stop. The largest remaining material residual was 0.289840 against a tolerance of 0.000030.

All 256 overlap rejections occurred in prism (depth 5, radial 18, polar 0), triangle 0. All 6,912 prisms in the saved state passed the exact-support guard. Following two material points distinguishes compression at the original location from the new limit at the edge of the order field.

| Material point evaluated on the final feasible state | Order χ | Physical volume ratio J | Elastic energy density W |
|---|---:|---:|---:|
| Original fold witness | 5.08830 × 10⁻⁶ | 0.00366956 | 5.23053 |
| Location of the final blocked proposal | 9.79319 × 10⁻²¹ | 7.41281 × 10⁻¹⁴ | 29.8372 |

The second point lies at physical radius 3.55318821, depth 0.20214128 and polar angle 0.17885882 radians. It tightens the bound on the infimum of J over positive order to [0, 7.412813899151272 × 10⁻¹⁴]. The lower bound permits zero as a limit. Every point with positive order passed the strict orientation check. The measured witness is above the separate machine-zero threshold of 1.42109 × 10⁻¹⁴.

## What happened to the drain

| Reading | Checkpoint 10 | Cleared shape before crystal relaxation | Stalled diagnostic state |
|---|---:|---:|---:|
| Exit radius | 3.30638 | 3.27757 | 3.27757 |
| Depth | 0.81914 | 0.82009 | 0.82009 |
| Inward feed current | 10.70857 | 10.57368 | 10.57355 |
| Dense, subsonic share of inward incident-disk current | 66.7562% | 66.1106% | 66.1109% |
| First sampled depth fully melted across the disk through the exit | 0.28460443 | 0.28460443 | 0.28460443 |

Most of the current change occurred on the shape step. The final feed and bulk currents agree within 1.25 × 10⁻¹⁴. The sampled fully melted front stayed at the same depth while partial order changed. The last measured middle-width and exit forces are checkpoint 10's +0.70549933 and −0.64123640. Their subsequent trend is unresolved. [Readback of saved fields](../tools/mouth_r4_feasible_stall_readback.py)

The open question is whether the narrow region of positive order can carry its load with a different material mesh. Local remeshing of this prism is the smallest test of that hypothesis, provided it preserves the surrounding active crystal strain, energy and work. A reference map evolved on the spatial grid is the broader alternative. Radius 4 remains held and radius 3 remains blind.
