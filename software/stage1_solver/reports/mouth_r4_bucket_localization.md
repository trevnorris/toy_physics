# Where the radius-4 depth-force discrepancy arises

**Radius 4 stays held at checkpoint 10.** At degree 8, the 2× to 4× spatial-label depth force drops from +0.14461346 to +0.13626572. The original material force moves from +0.15040891 to +0.15040682.

The entries below are signed differences, spatial minus original material. “Order work” measures the force from moving the same crystalline order field through the two stored strain-energy densities, $-\int (\partial_{\ln L}\chi)\,\Delta(W/J)\,dV$. “Full force” also includes the movement of the labels and the physical geometry. Integration used the matched motions from `eda4ffed`.

| Degree 8 | Interior A | Weld B | Fold or melt seam C | All cells |
|---|---:|---:|---:|---:|
| 4× stored energy | −0.00108187 | −0.00047848 | −0.00102343 | −0.00258377 |
| 4× order work | −0.00212737 | −0.00062493 | −0.01197178 | −0.01472409 |
| 2× to 4× change in full force error | +0.00848854 | +0.00115465 | −0.01798885 | −0.00834566 |

C carries 81.31% of the finest-grid order-work error. Its explicit order work improves by +0.01401688 on the last refinement, while its remaining terms change by −0.03200573. The reversal survives changing quadrature degree and assigning overlap cells to the weld first. The historical energy deficit of −0.00238288 used degree 3.

All ten leading order-work cells lie in C. Cells intersecting the eight folded prisms contribute only −0.00003158. Every positive-order cell uses directly inverted material labels at all corners. Material depths 2/24 through 5/24 contribute −0.00936658 in the core and −0.00422111 in the collar. Replacing only the inverse-label gradient reproduces the order-work error to 0.00000324. [Raw values, 24 depth bands and ten ranked cells](mouth_r4_bucket_tables.md) · [Full-precision measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_bucket_localization.json)

**The unfolded control also fails.** The last pre-fold state had a minimum chart determinant of 4.70 × 10⁻¹⁶. Halving checkpoint 10’s crystal displacement gave a minimum of 0.50430411 over all 6,912 prisms with identical shape, order, carrier and reference metrics. At degree 8 the energy error was −0.00127197 → −0.00049191 → −0.00074416. The depth-force error was −0.00603682 → +0.00212353 → −0.00393637. Degree 6 reproduced both reversals. All five other force errors shrank twice at both degrees. Splitting the near-static energy integral at moving cell and metric boundaries also reproduced the reversal. [Unfolded measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_bucket_unfolded.json)

The continuum limit remains open. Shifting the interior label grid on this unfolded state while holding the physical order knots and boundaries fixed would test sensitivity to where the grid samples the original strain gradients. Independent review concluded “The three-bucket localization is verified. The unfolded control fails the Stage 1 convergence gate for stored energy and log-depth force at both quadrature degrees. Radius 4 must stay held before Stage 2.”
