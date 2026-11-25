# Radius 4 at step 395 · Transported shape work

The cubic two-metric evaluator reproduced step 395 from `e62a30ff` with zero measured change in all three energies, eight native residuals, six shape forces, ten flow readings and the physical orientation bound. No optimization step was taken. Both transported probes rebuilt the registered cubic and aligned gauge on the moved mesh. The whole-cell certificate covered 28,080 cells on each geometry.

| Probe norm | Transported energy change | Linear prediction | Remainder | Physical lower bound on Jinv |
|---:|---:|---:|---:|---:|
| 0.00010000 | -0.000225213356316 | -0.000227770274018 | +2.556917702e-06 | 1.201593814e-05 |
| 0.00005000 | -0.000113245885586 | -0.000113885137009 | +6.392514229e-07 | 1.2015599974e-05 |

Halving the displacement reduced the remainder to 0.250009 of the full-width value. The predeclared checks permitted a full-width remainder up to 2% of the linear work and a half-width remainder up to 35% of the first, with a 16-ULP allowance for subtracting the total energies. This measures first-order consistency on a single direction.

| Coordinate | Signed force | Transported diagonal curvature at 2⁻¹⁸ | At 2⁻²⁰ |
|---|---:|---:|---:|
| Log depth | -0.006656262 | 67.76071828 | 67.74291214 |
| Scaled log exit | -0.638270855 | 9.92867734 | 9.91504605 |
| Mean width | +0.660848084 | 167.81559325 | 167.80613096 |
| Depth asymmetry | -0.527998742 | 113.25196311 | 113.23379808 |
| Waist profile | +0.618248729 | 35.78145104 | 35.77402379 |
| Polar–equatorial | +0.006866546 | 8.80610220 | 8.79657415 |

The inverse-BFGS initialization uses the reciprocal of the absolute wide-probe curvature with a lower curvature floor of one. Internal settlement on a finite displaced shape is the open test. [Measurements and source-state digest](../../../research/pde_ledger_v3/_measurements/mouth_r4_release_zero395.json) · [Zero-step evaluator](../tools/mouth_r4_release_pilot.py) · [Nested solver](../tools/mouth_r4_release_outer.py)

**Independent local review of the numerical gate.** “Scoped clear for the Phase 1 held-state equivalence and finite transported-work gate at step 395. The reviewed evidence admits the bounded primary pilot.”

**Independent local review of the outer driver.** “Scoped clear for Phase 1 and the reviewed outer-driver logic.”
