# Radius 4 clears the frozen spatial-label tests

**Stage 1 passes on the held checkpoint-10 mouth.** With the same order field, reservoir phase drop of 0.4510084282857874 and two welded reference geometries, a smooth reference map gives monotonically decreasing errors in stored energy and all six signed forces. The test permits internal relaxation from checkpoint 10. Equilibrium, continued drainage and gauge independence during relaxation are separate measurements.

| Smooth-map error, spatial minus analytic, degree 8 | 1× | 2× | 4× |
|---|---:|---:|---:|
| Stored crystal energy | −0.00208441 | −0.00051999 | −0.00012961 |
| Log depth | −0.00349613 | −0.00088794 | −0.00023299 |
| ↳ Order work | −0.00070328 | −0.00016981 | −0.00004238 |
| ↳ Remaining depth work | −0.00279285 | −0.00071814 | −0.00019061 |
| Scaled log exit | −0.00215608 | −0.00053858 | −0.00013485 |
| Mean width | −0.00951194 | −0.00237232 | −0.00059373 |
| Depth asymmetry | +0.00788669 | +0.00195995 | +0.00049005 |
| Waist profile | −0.00287456 | −0.00072915 | −0.00018334 |
| Polar–equatorial distortion | −0.00001528 | −0.00000383 | −0.00000095 |

Degree 6 reproduced every sign and trend. At 4×, the interior, weld and melt-front neighborhoods contributed −0.00014799, −0.00003026 and −0.00005474 to the depth error. The interior rate on the last refinement was 2.05. The melt-front rate changed from 2.09 to 1.58, so its asymptotic order remains open. Including the independently differenced reservoirs, the largest discrepancy between analytic force and signed-energy derivative was 5.02 × 10⁻⁸ using the smaller central depth step and extrapolated central differences for the other five coordinates. [Six forces, regional splits and both quadratures](../../../research/pde_ledger_v3/_measurements/mouth_r4_clearance_smooth.json)

**The unfolded control is sensitive to where the grid samples the material strain.** Removing harmonic extension changed both signed ordered label banks by exactly zero on all three grids. Shifting interior depth knots by a quarter cell changed degree-8 depth errors from −0.00603682, +0.00212353, −0.00393637 to +0.00501805, −0.01678618, −0.00917340. Most shift sensitivity was outside the melt-front neighborhood at 2× and inside it at 4×. The original material map has gradient jumps in both places. [Splice-free energy and force breakdown](../../../research/pde_ledger_v3/_measurements/mouth_r4_clearance_shift.json)

The smooth radial label was (\xi_x=x+0.12x(1-x^2)\sin(\pi t)), where (x,t) are normalized radius and depth. Its normalized Cartesian forward determinant is at least 0.71178. This is the chart determinant used for the unfolded control. The physical elastic determinant also contains the two stored reference metrics. Its reference-volume mean is 0.45211020 for this held boundary, precluding a uniform physical bound of 0.5. [Map bounds and volume audit](../../../research/pde_ledger_v3/_measurements/mouth_r4_clearance_regions.json)

Independent review concluded “Stage 1 cleared for the held R4 checkpoint-10 two-metric Eulerian relaxation. No equilibrium or gauge-independence verdict.” The next test advances the unpinned labels, order, carrier and six shape coordinates, and compares the two gauge strengths at checkpoint 11. Radius 3 stays reserved.
