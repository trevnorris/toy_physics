# What carries the load at the radius-four exit?

**The crystal's internal angular displacement equilibrated. The next proposed shape still failed its material-force test near the exit.** The accepted wall had depth 0.79590446 and exit radius 3.50122632. The crystal remained ordered in the interior, the core and collar shared a mobile bond, and one circulation quantum continued into the receiving superfluid.

| Fixed-wall result | 24 × 24 solid | 33 × 32 solid |
|---|---:|---:|
| Energy released | 1.31422 × 10⁻⁶ | 1.31399 × 10⁻⁶ |
| Largest extra radial displacement | 0.0000917801 | 0.0000917826 |
| Largest free meridional material force | 8.04 × 10⁻⁶ | 4.89 × 10⁻⁶ |
| Minimum sampled principal stretch before / after | 0.32919844 / 0.32919818 | 0.30429934 / 0.30429934 |
| Largest absolute global shape force before / after | 1.218387 / 1.218384 | 1.757456 / 1.757456 |
| Polar wall force before / after | −0.002029 / +0.013042 | −0.001959 / +0.013111 |

Small displacements to either side of the angular minimum raised the energy. All sampled volumes stayed positive. The finer solid resolved a smaller stretch near the exit, already present before release. The angular adjustment shifted the polar wall force by +0.01507 on both meshes while leaving the dominant global force essentially unchanged. On their common carrier grid the fine solid carried 11.61684 inward current, 70.4769% dense and subsonic including swirl. [Angular calibration](mouth_fixed_wall_angular.md) · [Paired shape forces](../../../research/pde_ledger_v3/_measurements/mouth_fixed_wall_support.json)

**The earlier stall was on a different wall.** We recovered that proposed geometry from the recorded descent rule and held it fixed. The angular mode settled there, but the final material residual was 0.00153713 against a tolerance of 0.00003000. In the last saved sweep before that rejection, the largest force was radially inward at radius 3.46628 and 97.76% of the exit depth. An independent signed displacement of that node reproduced its force of −0.00186587. The tested angular mode moves material only between 6% and 62% of the depth. [Stall and node measurements](../../../research/pde_ledger_v3/_measurements/mouth_fixed_wall_support.json) · [Independent checks](../../../research/pde_ledger_v3/_measurements/mouth_closed_rim_review.json)

**Both control volumes enclosed the entire disk and rim.** Including the stress on the extra sloping wall inside the thicker volume accounted for 0.04484 of the 0.04646 change in axial flux on the original carrier grid. The corrected signed residuals were −0.007713 and −0.006084. Doubling the carrier and conversion-layer grid at identical physical surfaces gave −0.052753 and −0.052639. The two thicknesses now agree much more closely, but the absolute force shifts by about 0.045 under refinement. The integrated mass imbalance was below 2.1 × 10⁻⁹ on both grids. A common virtual translation includes the held constraints. Assigning a separate force to the prescribed order at the shared rim still depends on which side supplies its slope. [Full momentum and constraint-work audit](mouth_closed_rim_work.md)

Whole-mouth release remains gated. A local calculation on the last rings of bonded material at the rejected wall can now target the signed radial and axial residuals, with independently refined solid and carrier stresses on the same physical cuts. Radius 3 remains blind.
