# Does the neutral patch recover the brane pressure law at its own resting spacing?

**The patch contracts by 0.248592% and still stiffens too rapidly.** We let the same cold, neutral, registered triangular layers settle at zero applied in-plane stress, keeping their separation fixed at one bond length. The bonds selected a pitch of **0.9975140800** in units of the isolated within-layer bond length. Both in-plane normal stresses and the shear stress vanished at this spacing. Zero applied stress is a trial preparation. The ambient stress of the physical brane is still open. [Calculation](../tools/neutral_patch_zero_load.py) · [Previous fixed-spacing comparison](neutral_patch_periodic_compression.md)

Starting from that spacing, we strained all three brane directions equally at conserved constituent number, keeping the fourth-direction gap fixed. Here \(n\) is density divided by the density of this new resting state. Only the pressure rise at \(n=1.01\) sets the energy scale. Every other row is a prediction. The pressure is the mean stress conjugate to that loading path, in solver units.

| Density change | Brane pressure change | Patch pressure change | Error in magnitude |
|---|---:|---:|---:|
| −2% | −0.01921584 | −0.01836149 | −4.446% |
| +1%, calibration | +0.01020201 | +0.01020201 | fitted |
| +2% | +0.02081616 | +0.02111620 | +1.441% |
| +5% | +0.05525631 | +0.05837168 | +5.638% |

At the former unit pitch, the +5% pressure excess was 5.902%. The new result held for 24 and 200 sites, after extending the interaction range, and when the pressure was checked against the change in energy. An independent lattice sum reproduced the predictions and found positive in-plane stiffness on a 36-by-36 grid of wavevectors at rest and at the two ends of the tested density interval, apart from the common translations. Local verdict: “Scoped clear for the cold registered, fixed-gap, equitriaxial pressure-change test and its sampled in-plane stability.” [Measurements, reviewer calculation and verdict](../../../research/pde_ledger_v3/_measurements/neutral_patch_zero_load.json)

**Material requirement.** Write \(U(u)\) for energy per two-site cell at log in-plane stretch \(u\), and \(\eta=\log n\). The work relation gives \(p=-e^\eta U'(-\eta/3)/3\) in reference-volume units. At zero load \(U'(0)=0\), so the scale-independent pressure curvature is

\[
\frac{p_{\eta\eta}(0)}{p_\eta(0)}=2-\frac{U'''(0)}{3U''(0)}.
\]

The brane pressure law \(\Delta p=(n^5-1)/5\) requires a curvature of **5**. The independent calculation gives **7.98802124**. Thus any material proposed to supply this pressure by the same fixed-gap mechanism must have \(U'''(0)/U''(0)=-9\) at equilibrium. This necessary condition can be checked on an isolated patch before using it on a mouth. It constrains the first nonlinear correction and leaves the full finite-density curve to be tested.

**0PN decision.** Carry this pressure requirement forward. The tested bonds can still contribute local attachment forces. Predicting free fall requires an independent total sideways force per actual wall area, including where the material on both sides transfers its reaction. The compression test supplies no such force. Working base `1dc45d09`. Independent local review complete. Formal review and commit remain pending. [Balanced-load test and missing input](mouth_deep_support_calibration.md)
