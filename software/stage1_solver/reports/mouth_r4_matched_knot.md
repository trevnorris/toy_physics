# Radius 4 — depth force with identical moving order knots

**Radius 4 stays at checkpoint 10.** Matching the order-field motion and the physical integration points leaves a signed depth-force discrepancy of **−0.01414** on the finest tested grid. It grows between 2× and 4× at both quadrature degrees. This triggers the Stage 1 stop in the directive following `7b45bfec`. The test used the saved, held state throughout. Checkpoint 11, equilibrium, collapse and free fall remain open. Radius 3 remains blind.

The table gives the Eulerian force minus the original material force for the same log-depth variation, with \(F_q=-\partial E/\partial q\). At each sign, both evaluations see one set of physical points, weights and order values. The saved material map moves with the mouth. The spatial labels follow its inverse at the displaced knots. The signed quadrature splits the crystal on both sets of spatial and reference-metric boundaries. The fluid, order-field and exterior-crystal forces are common to both evaluations and cancel in the difference.

| Spatial grid | Degree 6 | Degree 8 |
|---|---:|---:|
| 1× | −0.00908246 | −0.00904741 |
| 2× | −0.00577563 | −0.00579545 |
| 4× | −0.01414182 | −0.01414111 |

At 4× and degree 8, adding the common field contribution gives total depth forces of **+0.15040682** in the original material representation and **+0.13626572** in the spatial representation. Halving the signed displacement from \(4\times10^{-7}\) to \(2\times10^{-7}\) changed their difference by \(1.39\times10^{-9}\).

**Where the difference enters.** In physical volume, the local coefficient of order \(\chi\) in the crystal energy is \(W_{\rm crystal}/J\). The work from moving the common order field therefore contains

\[
\Delta F_{\chi}=-\int (\partial_q\chi)\,\Delta(W_{\rm crystal}/J)\,d^4x .
\]

An independent degree-8 evaluation at the saved material points gives **+0.00148908** from geometric stress and **−0.01472409** from moving order at 4×. Degree 6 gives +0.00144152 and −0.01463491. These partials hold the normalized label function fixed. The table includes the additional change as spatial labels are sampled again on the moved grid. The coefficient of the moving order field accounts for most of the finest-grid discrepancy. Its magnitude alone decreases with refinement. The full depth-force error has a different trend over these three grids.

**Controls.** Both charts moved all physical boundaries identically within roundoff. Over 155,358 samples per shape, the maximum difference in signed displacement was \(7.11\times10^{-15}\). The entrance ring and reservoir roof stayed fixed. The exit disk moved by 0.81914036 per unit log depth. Direct central energy differences verified all six production shape partials to \(1.39\times10^{-8}\) with fixed spatial nodal coefficients. Independent inversion at both signed step sizes verified the transported labels to \(9.86\times10^{-15}\). The phase drop was 0.4510084282857874 throughout.

Independent review concluded “The matched-knot depth-force gate fails at checkpoint 10. I verified both signed maps use a shared integration rule, both external domains move together, and the six production partial derivatives pass.” The continuum limit remains undetermined. The remaining question is where the projected strain changes the work of the moving conversion layer. Summing that signed work by physical depth and bond side would localize it on this same checkpoint. [Signed energies, forces, source hashes and boundary checks](../../../research/pde_ledger_v3/_measurements/mouth_r4_matched_knot.json)
