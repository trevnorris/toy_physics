# Radius 4 · Stage 2 stopped at fixed shape

**The held mouth still drains. Both Eulerian label relaxations stalled before the first shape step.** Checkpoint 10 remains the last whole-mouth state. Stage 1 cleared the two-metric evaluator at `fb0347a8`. The Stage 2 solver at `01622952` used one quantum of circulation, phase drop 0.4510084282857874 and the same material law. [Stage 1](mouth_r4_clearance.md) · [Stage 2 measurements and saved-state hashes](../../../research/pde_ledger_v3/_measurements/mouth_r4_eulerian_stage2.json)

| Fixed-shape reading | Gauge 10⁻⁸ | Gauge 10⁻¹⁰ |
|---|---:|---:|
| Energy-lowering label steps | 50 | 36 |
| Last accepted step fraction | 3.55 × 10⁻¹⁵ | 1.42 × 10⁻¹⁴ |
| Total energy at stop | 92185.49685485 | 92185.49240873 |
| Energy change during label descent | −0.03143979 | −0.03588591 |
| Largest free-label force (target 0.000030) | 0.12656537 | 0.11469149 |
| Limiting cell (radial, depth, polar index) | (18, 8, 19) | (18, 8, 5) |
| Polar angle at shared depth 0.23285817 | 1.465706 | 0.241735 |
| Melted-side physical inverse volume ratio | 1.78657 × 10⁻¹³ | 6.99635 × 10⁻¹³ |
| Physical label-gradient norm at corner | 0.565343 | 0.592681 |
| Formal crystal energy density at corner | 5.61549 × 10²⁵ | 6.38766 × 10²⁴ |
| Inverse volume ratio on adjacent partly ordered side | 1.34711 | 1.79085 |

All eight vertices in each limiting cell had exactly zero crystalline order, so its crystal energy density was zero. Independent one-sided reconstructions verified the table. The final rejected primary trial had a positive chart determinant below its roundoff margin. The secondary trial crossed zero. Both last accepted maps passed the orientation checks. Their crystal energies differ by 0.00444611 at different, internally unbalanced states. Gauge independence at equilibrium remains unmeasured.

| Signed shape force | Original checkpoint 10 | Gauge 10⁻⁸ | Gauge 10⁻¹⁰ |
|---|---:|---:|---:|
| Log depth | +0.150565 | +0.008780 | +0.075389 |
| Scaled log exit radius | −0.641236 | −0.674023 | −0.681712 |
| Mean width | +0.705499 | +0.573277 | +0.546254 |
| Depth asymmetry | −0.666834 | −0.576995 | −0.565523 |
| Waist profile | +0.465603 | +0.397602 | +0.373970 |
| Polar–equatorial distortion | +0.009007 | +0.020186 | +0.022610 |

The last two columns give virtual work at fixed fields after transfer from the original material map, field rematching and partial label relaxation. Five forces agreed with two-sided energy differences within 3.10 × 10⁻⁸. For log depth, only decreasing depth kept the labels oriented at the two probe sizes. [Signed force components and one-sided tests](../../../research/pde_ledger_v3/_measurements/mouth_r4_eulerian_stall_work.json)

Exit radius 3.30637727 and depth 0.81914036 were held. The saved carrier and order fields were identical across gauges. The inward current was 10.70852202 (checkpoint 10 gave 10.70857269). Feed and bulk agreed within 5.65 × 10⁻¹¹. Dense, subsonic material carried 66.75637% of inward incident disk current. The first sampled cross-section fully melted through to the exit stayed at depth 0.28460443.

Independent review concluded “Neither a coupled checkpoint 11 nor a physical collapse was established.” At both stalled corners, moving a neighboring label belonging exclusively to melted cells changes the local determinant. The next discriminating test is to hold all active labels fixed and measure whether an orientation-preserving extension through the melt exists, including its gauge cost. Radius 4 stays held for that decision and radius 3 stays blind.
