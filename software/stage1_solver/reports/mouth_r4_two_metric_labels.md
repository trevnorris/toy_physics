# Radius 4 at checkpoint 10 — one welded label field and two natural geometries

**Stage 1 failed checkpoint preservation.** On the 27 × 55 × 21 physical grid, one continuous label field \(\xi(x)\) carries the two original stress-free geometries. Each side uses \(G_s=(\partial X_s/\partial\xi)\nabla_x\xi\) and crystal energy \(\int\chi\det G_s\,W(G_s^{-1})\,d^4x\), including the hoop stretch. All 5,006 required inverse labels reproduced their positions within \(1.12\times10^{-16}\).

**Welded ring.** The saved core/collar hoop stretches were 0.99904009 / 0.78284373 at χ = 0.18527171. The reconstructed ring gave 0.99914255 / 0.78292402 at χ = 0.18493922. With 32× local refinement the stretch errors fell to \(7.03\times10^{-7}\) and \(5.51\times10^{-7}\).

**Static energy.** At 5,291,652 identical ordered points, crystal energy fell from **0.206520015646 to 0.201259290011** and total energy from **92185.53614224629 to 92185.53088152065**. Of the 0.005260725634 deficit, 0.005072186467 came from points staying on the same reference side. The largest absolute matched-point element change was +0.000668557931 in prism (0, 16, 5), triangle 1. Independent physical-cell integration gave 0.200929363829 and 0.200841549833 at degrees 9 and 16. All eight held-field current and balance readings matched exactly, including the inward current of 10.708572690281159. [Element-wise energies](../../../research/pde_ledger_v3/_measurements/mouth_r4_two_metric_elements.tsv)

**Six total forces.** The analytic derivative transports label and order coefficients with the grid. Degree-6 crystal integration splits the collar’s reference triangles. Central energy differences at \(10^{-4}\) and \(5\times10^{-5}\) give the signed extrapolation. Its largest discrepancy is \(4.21\times10^{-8}\), compared with \(3.70\times10^{-7}\) for the raw half-step depth force.

| Shape coordinate | Saved | Projected analytic | Signed extrapolation | Projected − saved |
|---|---:|---:|---:|---:|
| Log depth | +0.15056488 | +0.25719448 | +0.25719451 | +0.10662959 |
| Scaled log exit | −0.64123640 | −0.63914599 | −0.63914600 | +0.00209041 |
| Mean width | +0.70549933 | +0.71169550 | +0.71169548 | +0.00619616 |
| Depth asymmetry | −0.66683446 | −0.66631074 | −0.66631074 | +0.00052372 |
| Waist profile | +0.46560293 | +0.47689431 | +0.47689435 | +0.01129138 |
| Polar–equatorial shape | +0.00900742 | +0.00926768 | +0.00926768 | +0.00026026 |

**Depth motion.** Re-inverting the saved map on the displaced depth knots changes the force by −0.11284884. The correction varied by \(1.5\times10^{-9}\) over the three smallest steps at this quadrature degree. The resulting **+0.14434564** leaves **−0.00621924** against the saved depth force. The earlier +0.41901 compared fixed and moving *order knots*, already included here. The other five shape motions preserve the normalized knots.

**Gauge and orientation.** We minimized the harmonic gauge in fully melted cells with registered wall/end labels and natural conditions on the remaining free surfaces. At ε = 0, 10⁻¹⁰ and 10⁻⁸ the frozen hoop stretches agree. The largest ε adds 6.309 × 10⁻⁹ energy and changes the largest active-node and shape-force components by 1.1243 × 10⁻⁸ and 3.9721 × 10⁻⁹. Every cell has positive chart-determinant Bernstein coefficients. Maximum sampled physical condition numbers were 4.068 in ordered material and 4.909 in melt. Ordered inversions were rejected on both sides (\(J_{\rm inv}=-1.84672,-1.36192\)). [Signed energies, controls and independent readouts](../../../research/pde_ledger_v3/_measurements/mouth_r4_two_metric_labels.json)

Independent review concluded “The matched-point elastic-energy deficit persists.” The next read-only test is to preserve one-sided strain within ordered cells and repeat the identical-point comparison. Radius 4 stays held before checkpoint 11. Radius 3 remains blind.
