# Axial load at the radius-4 exit

**The disk-stress discrepancy is concentrated near the rim. The full-disk force still depends on the measurement thickness.** We compared drive and rest on the same 33 × 32 bonded lattice from `432eccb9`, then rematched both fields on a denser carrier grid. The wound current changed from 11.64367 to 11.58261. [Measurements and independent review](../../../research/pde_ledger_v3/_measurements/mouth_disk_pillbox.json)

| Drive-induced axial force, incident / receiving | Native fields (13 × 10 × 14) | Refined fields (19 × 14 × 22) |
|---|---:|---:|
| Direct carrier and order trace | 1.08657 / 0.44241 | 1.07270 / 0.48723 |
| Additional direct crystalline traction | −0.00906 / 0 | −0.00929 / 0 |
| Full weak force, half-thickness 0.03429 | 0.48527 / 0.47470 | 0.47411 / 0.52192 |
| Full weak force, half-thickness 0.06857 | 0.38912 / 0.47052 | 0.43662 / 0.51916 |
| Full weak force, half-thickness 0.12 | 0.34853 / 0.46791 | 0.42170 / 0.51761 |
| Interior-weighted weak force, half-thickness 0.12 | 0.14854 / 0.14921 | 0.14926 / 0.14917 |

Phase motion including order inertia supplies **1.01881 / 0.40654** of the native direct traces. Without winding, the full carrier and order traces are **0.97879 / 0.40460**. The axial virtual-work test includes advected momentum, gradient stresses, native crystal work, and lateral traction on the sloping wall and receiving cylinder. At native half-thickness 0.03429 those lateral terms are **0.43579 / 0.07303**. Tapering the test from unit weight at 75% radius to zero at 95% isolates interior transmission. Signed energy differences verify the volume work. Fixed-geometry curvature work cancels between drive and rest.

A shared cap-node variation gives phase currents **+11.64367245** from the duct, **−10.48792412** from directly above and **−1.15574833** from neighboring elements outside the rim. They sum to **9.03 × 10⁻⁹**. Momentum additionally probes the mechanical restraints. On the refined grid at half-thickness 0.03429, contracting the prescribed-order reaction with the incident or receiving axial slope gives **+0.01318** or **+0.10345**. The held solid nodes contribute **+0.00907** for an axial displacement at fixed Eulerian order. The separate variations do not select one rim force.

**After the exterior moved**, two solid meshes settled internally at identical shape coordinates, depth 0.795904 and exit radius 3.501226. Their direct total stresses were **0.85185 / 0.37637** (24 × 24) and **0.85088 / 0.37667** (33 × 32). On the refined solid the full weak pair changed from **0.52913 / 0.41784** at half-thickness 0.03429 to **0.41439 / 0.41374** at 0.12. Solid refinement reproduces the thickness dependence. The refined wound/unwound currents were **11.61684 / 13.78974** and the dense, subsonic incoming share was **70.48%**. Exterior forces remained.

The next test must enclose the bonded-to-free rim from all adjacent banks and measure material-boundary and prescribed-order reactions under one displacement as local resolution increases. Whole-mouth support and free fall remain open. Independent verdict: “scoped clear for the held R4 drive−rest force diagnostics, the two rematched field meshes and the conclusion that whole-disk axial traction remains unresolved.” [Held-mouth result](mouth_joint_collar_release.md)
