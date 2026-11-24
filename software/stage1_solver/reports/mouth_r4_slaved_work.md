# Radius 4 · shape work with the melted labels in equilibrium

At the two held geometries from `12e02d90` the mouth has the same dimensions and the two crystal traces have different residual stresses. Starting from the orientation-preserving extensions in `743c4a1b`, each signed perturbation held the physical fields and crystalline order fixed and re-equilibrated the strictly melted labels. The derivative of this reduced energy includes the change in the physical integration volume and gradients when the mouth moves. [Starting states and extension test](mouth_r4_slaved_melt.md)

| Signed shape force on the held state | Gauge 10⁻⁸ | Gauge 10⁻¹⁰ |
|---|---:|---:|
| Log depth | +0.00877970 | +0.07538911 |
| Scaled log exit radius | −0.67402272 | −0.68171220 |
| Mean width | +0.57327658 | +0.54625351 |
| Depth asymmetry | −0.57699453 | −0.56552268 |
| Waist profile | +0.39760177 | +0.37396996 |
| Polar–equatorial deformation | +0.02018597 | +0.02261000 |

All 56 signed perturbations of the six shape coordinates and the active radial label at the previously blocking corner preserved orientation. At both step sizes and both gauge strengths the largest difference between the analytic and centered shape forces was 3.10 × 10⁻⁸ against a 10⁻⁷ tolerance. On the active radial label the largest difference was 3.83 × 10⁻¹². The largest unscaled residual on the independently relaxed melted labels was below 3.0 × 10⁻⁸. [Signed measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_barrier_work.json) · [Precision audit](../../../research/pde_ledger_v3/_measurements/mouth_r4_barrier_depth_precision.json)

The original half-step log-depth discrepancies were −2.19576 × 10⁻⁷ and −2.19560 × 10⁻⁷. Re-evaluating the reservoir order-gradient and double-well energies on the same signed states in extended precision reduced them to +2.82 × 10⁻⁹ and +2.84 × 10⁻⁹. At full step the corrected discrepancies were +8.18 × 10⁻⁹ and +1.37 × 10⁻⁸. The original values remain in the measurement record.

A separate signed displacement moved all 14,433 free label groups together, including 4,887 radial, 4,898 axial and 4,648 angular groups. At both gauges and both widths the energy derivative agreed with the projected analytical force within 1.83 × 10⁻¹¹. This checks one mixed direction through the full vector of free labels. [Mixed displacement and independent readback](../../../research/pde_ledger_v3/_measurements/mouth_r4_barrier_mixed_work.json)

The exit pull is inward and the mean-width force outward on both saved traces. A physical comparison between gauge strengths requires independently settling the interior at each gauge and measuring the forces again. The force balance required for free fall remains open.

The independent force review concluded, “All 56 signed perturbations passed the orientation and inner-stationarity checks. The 28 central differences agree with their analytic counterparts at both original widths within 1e-7 when the two extensive order terms use extended precision.”
