# Radius 4 — frozen spatial refinement at checkpoint 10

**The energy deficit shrank at each refinement. The error in the transport-corrected depth force grew.** The joint criterion for starting relaxation failed. Radius 4 remains held at checkpoint 10 with inward current 10.708572690281159. The physical fate of the open mouth remains unresolved. Radius 3 remains blind.

We subdivided each of the original spatial intervals in radius, depth and polar angle, retaining the original material map, welded core/collar metrics, trilinear label interpolation, continuously matched order, phase drop 0.4510084282857874 and fluid potential. The 2× and 4× grids had 35,905 and 271,385 required inverse-label nodes. Each had a unique inverse. At zero gauge only the cells touching positive order contributed elastic energy.

| Spatial refinement | Crystal energy at identical points | Deficit | Signed error from points on the same side of the weld |
|---|---:|---:|---:|
| 1× | 0.201259290011 | 2.54732% | −0.00507218647 |
| 2× | 0.203059569482 | 1.67560% | −0.00338180155 |
| 4× | 0.204137140519 | 1.15382% | −0.00235371803 |

The saved energy was **0.206520015646** at the same 5,291,652 ordered material points. The two empirical convergence orders were 0.604 and 0.538. Three grids establish a decreasing deficit over the tested range. Neither a nonzero limiting deficit nor quadratic convergence has been established.

**Signed total-force errors**, calculated as the projected force minus the saved force. The depth comparison carries the original material motion through the moving order knots. All three rows use degree-6 spatial integration with the same partition of the two reference metrics.

| Grid | Log depth | Scaled log exit | Mean width | Depth asymmetry | Waist profile | Polar–equatorial |
|---|---:|---:|---:|---:|---:|---:|
| 1× | −0.00621924 | +0.00209041 | +0.00619616 | +0.00052372 | +0.01129138 | +0.00026026 |
| 2× | −0.00646568 | +0.00139192 | +0.00575042 | −0.00404245 | +0.00288983 | +0.00036416 |
| 4× | −0.01523021 | +0.00040976 | +0.00139524 | −0.00038965 | +0.00169728 | +0.00012232 |

The independent degree-8 depth residuals were −0.00396732, −0.00672632 and −0.01353225. Changing the integration of the original material state shifted its depth force by at most 0.000219 across the native degrees tested. Both comparisons use the archived force as their reference.

**Interior localization.** At 4×, using the projected inverse-label gradient with the *original* material labels and reference metric accounts for −0.00237742 of the −0.00238288 energy difference. The remaining −0.00000545 enters when the label values and metric are also transferred. The largest depth-band deficit, −0.00156557, lies between 2/24 and 5/24 of material depth. Both core and collar contribute. All sampled ordered points retained positive orientation.

**Depth-work control.** Keeping the label function fixed at the original normalized physical points, we varied physical depth and moved the order knots. At 4× the geometric stress-work error was +0.00073073 and the order-knot work error was −0.01442227. This localizes sensitivity to the stored strain energy where order varies with depth. It uses a different virtual motion from the production force measurement. The open question is how much of the production residual comes from resampling the label field as the depth grid moves. The next discriminating check would evaluate both label motions with one common integration rule, keeping the order-knot motion identical. [Full signed tables, region bins, integration controls and source hashes](../../../research/pde_ledger_v3/_measurements/mouth_r4_spatial_refinement.json)

Independent review concluded “The transport-corrected depth force moves farther from checkpoint 10 under spatial refinement at both degree six and degree eight.”
