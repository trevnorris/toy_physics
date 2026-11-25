# Radius 4 · Why the step-288 corner projection failed

At fixed shape and fields, **a joint projection gave a strictly energy-lowering, positively oriented trial**. Checkpoint 287 remains the last accepted state. Starting from the archived field update at step 286, the zero-step diagnostic reproduced every step-287 label bit for bit, all seven backtracks and the signed work −3.6625103882280213 × 10⁻⁶. The source state and solver stop are recorded in `f28f283f`. [Measurements, reproducing commands and independent review](../../../research/pde_ledger_v3/_measurements/mouth_r4_step288_projection.json)

At step 288 the eight-iteration Newton direction reached its inner tolerance. At trial fraction 1/16, the orientation guard detected an inversion in cell (25,10,4), local coordinate (0,0,1), at physical depth 0.25873130. Physical inverse determinant there was −0.00045075089. Adding its constraint exposed a defect in the four-row projection. SciPy 1.15.3's nonnegative least-squares call returned a residual of 0.0. Direct evaluation on the same matrix and multipliers gave 1.47205222 × 10⁻⁵. The computed direction violated the older constraint in cell (25,10,6) at (0,1,0.625), giving normalized directional slack −2.47542101 × 10⁻⁷ against a floor of −10⁻¹⁰. All four Gram eigenmodes survived. Reversing the equation order recovered the independently verified three-binding-row solution. This diagnosis applies to the captured input.

| Cell | Certified physical lower bound at 287 | Measured cell minimum at 287 | Measured cell minimum, rejected 1/16 trial |
|---|---:|---:|---:|
| (25,10,3) | +0.0022961410 | +0.0091391809 | +0.0003083131 |
| (25,10,4) | +0.0008042912 | +0.0020105818 | −0.0006093482 |
| (25,10,6) | +4.1095410 × 10⁻⁷ | +2.8537737 × 10⁻⁶ | +2.5184378 × 10⁻⁶ |

All three cells have χ = 0 throughout by trilinear interpolation. Each shares its four shallower vertices with neighboring elements carrying positive order. Cells (25,10,3) and (25,10,4) share the face with vertices (25 or 26, 10 or 11, 4). Cell (25,10,6) shares no vertex with either. The measured minima use closed-cell sampling and local minimization. The separate lower bounds certify positive orientation everywhere in these cells at step 287.

The diagnostic solved all retained corner constraints together, with a 70-digit multiplier check and a small positive margin to handle rounding. A trial exposed a fifth constraint in cell (25,10,3). With all five protected, a quarter step decreased physical crystal energy by 0.00010131418082, increased gauge energy by 4.77286377 × 10⁻¹⁰ and lowered total energy by **0.00010131370354**. The certified physical inverse determinant was at least **8.54324007 × 10⁻⁸** across 28,080 cells. The frozen-field peak crystal force fell from 0.04602958 to 0.03447917 on this test configuration.

Independent verdict, “Scoped clear for a feasible, energy-lowering zero-step trial and a numerical four-row projection failure. Convergence under subsequent steps remains open.” A production repair would need to check constraint feasibility and multiplier optimality on the captured case before a held run resumes. [Diagnostic](../tools/mouth_r4_step288_projection.py) · [Cell measurement](../tools/mouth_r4_step288_cells.py)
