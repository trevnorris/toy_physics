# Radius 4 · Smooth reference on the held step-166 collar

Following the crossing audit at `08a781cd`, a clamped tensor-product cubic interpolates the registered collar geometry. The reference determinant has a positive lower bound of 0.48942599 over all 72 patches. The maximum registered-knot position error is 1.11 × 10⁻¹⁶. Across internal faces, the bounds on jumps in position, first derivative and second derivative are 1.43 × 10⁻¹⁶, 1.71 × 10⁻¹⁵ and 2.20 × 10⁻¹³. The entire welded trace is preserved within 1.79 × 10⁻¹⁶. [Geometry and orientation](../../../research/pde_ledger_v3/_measurements/mouth_r4_smooth_collar_geometry.json)

All four boundary curves are straight. At the lateral wall, however, the stored axial tangent at the material label t = 2/3 is (0, 1.06508080) on the left and (0, 1.16861746) on the right. A differentiable map at fixed label cannot reproduce both. The cubic preserves the geometric curves, their unit tangents, the registration knots, and the full weld. Between knots the largest material position shifts along the lateral wall, entrance and exit are 0.00108954, 0.00222759 and 0.00470399, respectively. These are changes to the material spacing along the same physical boundary. [Independent trace audit](../../../research/pde_ledger_v3/_measurements/mouth_r4_collar_boundary_trace.json) · [Between-knot extrema](../../../research/pde_ledger_v3/_measurements/mouth_r4_smooth_collar_geometry.json)

With the original step-167 projected-Jacobi direction, both cells incident on the former crossing now lower physical crystal energy at α = 1/8192. The signed changes are −1.04280852 × 10⁻⁸ and −2.50477524 × 10⁻⁸. Smoother quadrature also changes the unchanged core's integrated energy by −3.09442434 × 10⁻¹¹, including the core portions of 2,009 rays that cross the weld. [Signed cell work](../../../research/pde_ledger_v3/_measurements/mouth_r4_smooth_collar_geometry.json) · [Core quadrature](../../../research/pde_ledger_v3/_measurements/mouth_r4_core_quadrature.json)

At the identical frozen fields, physical crystal energy rises from 0.31109285128 to 0.31536502625. The primary peak active crystal force rises from 0.10460517 to 0.32030122 and is now radial at (23, 1, 19), depth 0.02587313, χ = 0.32422060. Both gauges give zero shift in gauge energy and every saved flow reading. Inward current is 10.70772687 and 66.75897% is dense and subsonic. These comparisons hold the shape, labels and coherent fields fixed. [Two-gauge energies, forces and currents](../../../research/pde_ledger_v3/_measurements/mouth_r4_collar166_static.json)

| Frozen shape force | Linear, 10⁻⁸ | Cubic, 10⁻⁸ | Cubic, 10⁻¹⁰ |
|---|---:|---:|---:|
| Log depth | +0.03883339 | +0.00472231 | +0.00472233 |
| Scaled log exit | −0.55890386 | −0.58832640 | −0.58832639 |
| Mean width | +0.96217418 | +0.83084993 | +0.83084997 |
| Depth asymmetry | −0.70304637 | −0.59140129 | −0.59140132 |
| Waist profile | +0.81205505 | +0.77706111 | +0.77706114 |
| Polar–equatorial | −0.00503080 | −0.00555240 | −0.00555239 |

Independent review concluded, “The opt-in bicubic collar preserves the registered weld and boundary point sets, has continuous first and second derivatives across its 72 interiors, and removes the former positive-work reversal in both incident cells. This clears the geometry and local physical-work checks only.” Both gauges pass the unsplit log-depth check at h = 2⁻¹⁸ with an error of +2.72 × 10⁻⁹.

**Standing reference decision.** The user approved the measured between-knot shifts and adopted the clamped cubic rule from `7c32cbc9` for the two radius-4 gauges, subsequent shape release and future blind radius-3 evaluation. The reviewed frozen comparison is `3e73ce0b`. For radius 4 the SHA-256 of the registered coefficients in little-endian binary64 row-major order is `4c09ddff651e5b501e4b58014bd75b109be36878ed5e8d27b954709482ef7b12`. A future blind preparation uses the same registered-knot and clamped-tangent rule. Its geometry has yet to be evaluated.

Both fixed-label rematches must pass independent review before the primary ε = 10⁻⁸ Newton solve starts. The ε = 10⁻¹⁰ solve starts from the primary stationary endpoint, with native order, carrier, orientation and sign checks. Convergence requires simultaneous internal residuals at or below 0.000030. All six shape coordinates stay held through the paired stationary readout.
