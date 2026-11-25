# Radius 4 · Step-108 order rematch and force audit

**Both order rematches passed. The log-depth force failed the wider finite-difference check at both gauges.** The source is step 108 of `0afc9e1e`. Depth was held at 0.81914036 and exit radius at 3.30637727. Both gauges started from identical crystal labels and physical fields.

| Fixed-label measurement | ε = 10⁻⁸ | ε = 10⁻¹⁰ |
|---|---:|---:|
| Max bounded-order residual | 1.94581 × 10⁻⁹ | 1.94579 × 10⁻⁹ |
| Max amplitude / phase residual | 7.39149 × 10⁻¹⁰ / 6.15743 × 10⁻¹¹ | 7.39149 × 10⁻¹⁰ / 6.97470 × 10⁻¹¹ |
| Positive / zero / negative +10⁻⁹ order perturbations | 52,061 / 0 / 0 | 52,061 / 0 / 0 |
| Largest absolute gauge work / physical work | 7.56449 × 10⁻¹¹ | 7.56449 × 10⁻¹³ |
| Peak active crystal force magnitude | 0.09881032648 | 0.09881032561 |
| Inward current | 10.707728693785 | 10.707728693789 |
| Dense, subsonic share, including swirl | 66.7589636970% | 66.7589636970% |
| Absolute feed-to-bulk current difference | 4.86722 × 10⁻¹³ | 3.89022 × 10⁻¹³ |

Two formerly failing sites, 36353 and 62811, acquired order of 0.00011816 and 0.00010912. They were the only newly ordered material-grid sites. Site 45532 stayed melted and its perturbation work became positive. The first fully melted cross-section stayed at depth 0.33635069. Maximum order difference between gauges was 2.02933 × 10⁻⁹. [Full support changes and all-node work](../../../research/pde_ledger_v3/_measurements/mouth_r4_aligned_rematch.json)

**Historical rebound.** Peak force rose from 0.01756028 to 0.57995197 before the step-94 field update, which moved it to 0.58105854. Later, the determinant on the same melted edge fell from 0.00976947460 at step 100 to 0.00000239601 at 108. Fractions at steps 106–108 fell from 1/64 to 1/2048 and 1/4096. All 19 CG searches encountered nonpositive curvature, nine on the first iteration. Exact edge minima are recoverable on the two retained maps. [Step-by-step rebound and orientation trace](../../../research/pde_ledger_v3/_measurements/mouth_r4_rebound_trace.json)

**Derivative decision.** All 81,125 free crystal-label derivatives, 143,648 free order derivatives and five other shape modes passed at each gauge. For log depth the wide (2⁻¹⁸) central-minus-analytic errors were +0.000290135077 and +0.000290135373. Physical crystal work contributed +0.000290132315 and +0.000290132314. The narrow (2⁻²⁰) errors were +8.56 × 10⁻¹⁰ and +6.87 × 10⁻¹⁰. Both widths had a 10⁻⁷ limit. The six forces in the measurement are frozen-state readings. The entrance crystal still carries a signed force of −0.09881033 at depth 0.05174626. [Signed force and derivative banks](../../../research/pde_ledger_v3/_measurements/mouth_r4_aligned_derivatives.json)

**Edge and curvature controls.** All 16 Hessian-vector comparisons passed, with maximum relative error 2.464 × 10⁻⁷. The aligned barrier exerted an opening force on the saved edge and rejected a moving-point determinant of −7.18901 × 10⁻⁶ while all 27 fixed samples stayed positive. A closer positive probe in cell (25, 11, 10) exposed a second limitation. Its whole-chart lower bound was +6.47413 × 10⁻⁸, but the physical crystal kernel could not certify the probe. It returned a lower bound of −2.73823 × 10⁻⁸. That bound does not establish inversion. [Paired Hessian, orientation and restoring-force measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_aligned_hessian.json)

Independent review concluded, “The production solid could not certify the −0.95 positive-chart trial, and the wider log-depth force audit fails at both gauges. Gate 2 does not clear for crystal descent.” The open questions are which cells supply the wide depth-probe work and why the two orientation bounds disagree. Step 108 remains the last accepted crystal step.
