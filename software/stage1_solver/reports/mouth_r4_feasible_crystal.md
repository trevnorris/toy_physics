# Can a shorter crystal displacement keep the radius-four throat open?

The subsequent interleaved relaxation reached its hard stop before checkpoint 11. [Final material and flow measurements](mouth_r4_feasible_stall.md)

The starting shape has exit radius 3.2775726740, depth 0.8200888334 and inward current 10.5736799525. Its spatial fields have settled on the proposed shape. Checkpoint 10 remains the reference for the six forces until the crystal displacement, field update and signed shape probes all finish. The entrance is held at radius 4 with one quantum of circulation, the original order energy and phase drop 0.4510084282857874. [Original shape step and overlap](mouth_r4_exact_support.md)

## What folded in the previous trial

A readback followed the same material point in prism (5, 18, 0), triangle 0, across the displacement rejected by the old stop rule. At the start its order was 2.89624821 × 10⁻⁷, its signed volume ratio was +0.106350818, and its elastic energy density was 1.71879507. On the rejected displacement, its signed volume ratio became −0.001803079 while its order was still positive at 1.91818435 × 10⁻⁷. The sign change occurred as crystal coordinates moved through a fixed spatial order field.

Over the whole starting shape, the infimum of physical J where χ > 0 lies in [0.00012058322, 0.078801065], with these displayed bounds rounded outward. The diagnostic bounds this quantity because χ > 0 is an open condition and its minimum need not be attained. It reports pinning at machine precision only with an actual positive-order witness whose upper J bound is at most 64 machine epsilons. [Readout, force controls and exact values](../../../research/pde_ledger_v3/_measurements/mouth_r4_feasible_crystal.json)

## The displacement that previously stopped the run

On internal iteration 19 (counting from zero), the full proposal matched the archived rejected crystal coordinates byte for byte. Half and quarter steps still put positive order on inverted material in the same prism. The **one-eighth step** passed the whole-prism guard and lowered total energy by 0.00014446895655 at fixed shape and fields. The next direction admitted a **one-thirty-second step**, lowering energy by another 0.00008151137808. These are internal displacements at exit radius 3.2775726740. They establish a feasible descent past the earlier stop. The six shape forces are measured only after the subsequent field update.

## Acceptance and stopping rules

With shape, order, carrier amplitude and phase fixed, only the crystal energy changes. The line search halves the trial displacement until the whole ordered region is positively oriented and the energy strictly decreases. The recorded α is relative to the trust-limited direction. The accepted crystal-energy change gives the total energy change for that displacement without subtracting the large, identical reservoir background.

Two consecutive accepted steps with α < 10⁻⁴ and energy drops no greater than 32 floating-point spacings of the total energy trigger the interleaved test. At the starting energy this corresponds to 4.66 × 10⁻¹⁰. Exhaustion of the tested descent directions also triggers that test. The fallback updates order at the last feasible crystal placement, then alternates individual crystal steps with field updates. Only order coefficients whose positive interpolation support intersects J ≤ 0 must remain zero at fixed placement. All other free coefficients respond to the written material energy. Two successive unproductive interleaved updates, or a verified near-zero ordered volume on successive complete checkpoints despite interleaving, stop the run with a saved local constitutive readout. The coarse force target is 0.0004.
