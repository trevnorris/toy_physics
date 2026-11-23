# Radius-four bonded exit motion and exterior momentum

**Review finding.** Moving the entire bonded disk and its adjoining wall in the specified axial profile gives a differential reaction of +0.50016167 on the doubled grid and +0.49895823 on the rim-refined grid. These are zero-displacement derivatives for the common motion of all material and field meshes. The companion signed-displacement test measures how closely two independent carrier meshes reproduce this result. The force convention is minus the derivative of the draining-minus-resting energy. [Measurements and source hashes](../../../research/pde_ledger_v3/_measurements/mouth_shared_rim_surface.json)

| Axial measurement, four-point quadrature | Doubled | Rim-refined |
|---|---:|---:|
| Negative outward momentum through the near caps and side | +0.50016167 | +0.49895823 |
| Negative outward momentum through the far caps and side | +0.45412343 | +0.45303005 |
| Material work to change the far displacement taper to the near one | +0.04603824 | +0.04592818 |
| Near outward momentum plus the native one-cell wall reaction | −0.02303706 | −0.04042788 |
| Near outward momentum plus the wall reaction at the doubled-grid physical width | −0.02303706 | −0.02011176 |

The caps enclose the whole terminal disk. The side integral includes both the radial flux and the axial projection of the sloping normal. I integrated all three directly on each taper family and evaluated the material work separately using the nodal stresses, crystalline work and factorized-vortex surface term. Both procedures count the common rim once. In the undisturbed part of the collar touched by either external enclosure, the two states have identical order. On each grid, all-sector work and negative outward flux agree within 2.2 × 10⁻⁸ at native quadrature and within 10⁻¹¹ at four-point quadrature. This equality checks the force signs, angular measure and physical displacement. It follows from integration with the same weight and supplies no independent test of material stationarity.

The historical one-cell wall reaction and the complete disk displacement measure different quantities. Refinement halves the width over which a unit wall-node displacement decays, including the first cell past the disk. Preserving the original physical width gives a near residual change of +0.00292530. Narrowing to the refined grid's one-cell width at fixed resolution adds −0.02031612, accounting for the earlier −0.01739082 net change. The far surface gives +0.00298216 and −0.02031637 for the same two operations. Direct slit traction plus the enclosing flux still changes from −0.05342299 to −0.07774057. Two resolutions of that local sampling rule cannot establish a material singularity. [Surface calculation](../tools/mouth_shared_rim_surface.py)
