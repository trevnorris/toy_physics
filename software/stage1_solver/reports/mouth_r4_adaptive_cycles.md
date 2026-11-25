# Radius 4 at step 392 · Held-shape benchmark

**The primary paused at cycle 10.** Re-equilibrating order and carrier fields at fixed step-382 labels gave a peak crystal force of 0.00035689590. Each new cycle accepted one crystal step and then re-equilibrated the fields. The matched force reached 0.00003307470 at cycle 7, rose twice, and ended at 0.00005684467. The target was 0.000030. All ten CG solves met their forcing tolerances. The largest direct-versus-recursive residual difference was 2.78 × 10⁻¹⁷. [Exact Hessian audit](mouth_r4_exact_crystal_hv.md)

### Crystal and field work by cycle

Crystal forces, rebounds and signed energy changes are in millionths of model units. Order residuals are in thousandths and matched field residuals in billionths. The two work columns give the fixed-field crystal step and subsequent field update. Each row starts at the preceding matched force.

| Cycle / step | Crystal before → after field update | Rebound | Frozen order | Matched fields | Crystal ΔE | Fields ΔE | α / CG iterations | Seconds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 / 383 | 176.56 → 313.24 | +136.68 | 1.42984 | 2.405 | -2.1046 | -1.3364 | 0.5 / 13 | 181.5 |
| 2 / 384 | 155.51 → 203.40 | +47.90 | 1.74816 | 2.565 | -2.1962 | -1.7404 | 0.5 / 14 | 181.2 |
| 3 / 385 | 152.49 → 169.39 | +16.91 | 1.28582 | 2.643 | -1.3704 | -0.6140 | 0.25 / 16 | 179.8 |
| 4 / 386 | 84.42 → 112.35 | +27.93 | 1.16750 | 2.753 | -1.8795 | -0.6055 | 0.5 / 17 | 150.6 |
| 5 / 387 | 56.08 → 76.08 | +20.01 | 2.38193 | 2.926 | -1.9575 | -1.9643 | 0.5 / 18 | 161.7 |
| 6 / 388 | 11.90 → 51.81 | +39.91 | 1.56526 | 3.077 | -1.1794 | -0.7510 | 1 / 18 | 155.6 |
| 7 / 389 | 10.76 → 33.07 | +22.31 | 1.45462 | 3.309 | -1.7962 | -0.7741 | 1 / 20 | 178.6 |
| 8 / 390 | 11.63 → 56.40 | +44.76 | 2.55630 | 3.628 | -1.4239 | -1.9002 | 1 / 54 | 210.1 |
| 9 / 391 | 42.80 → 69.61 | +26.81 | 1.54293 | 3.787 | -0.7499 | -0.6408 | 0.25 / 19 | 180.6 |
| 10 / 392 | 38.51 → 56.84 | +18.34 | 1.44216 | 3.953 | -1.9727 | -0.5628 | 0.5 / 19 | 191.0 |

### Where the adjustment occurred

Peak sites are radial-label forces at material-grid nodes. A node with χ = 0 can share a crystal label with an ordered element. Reversals count the previous cycle’s toggles that flipped back. The cosines compare successive displacements on the always-ordered labels, the full ordered labels and the free order field. The last column gives the certified lower bound on the physical inverse determinant, in millionths, over all 28,080 cells.

| Step | Peak node · depth | Local χ | Off→on / on→off | Reversed / previous | Cosines, persistent / full / χ | Bound on J_inv |
|---:|---|---:|---:|---:|---:|---:|
| 383 | (25,8,14) · 0.20699 | 0 | 57 / 44 | — | — | 12.6348 |
| 384 | (25,8,14) · 0.20699 | 0 | 60 / 37 | 63/101 (62.4%) | 0.804 / -0.022 / 0.528 | 11.9366 |
| 385 | (25,8,14) · 0.20699 | 0 | 33 / 36 | 66/97 (68.0%) | 0.366 / -0.055 / 0.443 | 1.2750 |
| 386 | (25,8,14) · 0.20699 | 0 | 57 / 32 | 65/69 (94.2%) | 0.968 / 0.777 / 0.641 | 7.6274 |
| 387 | (25,8,14) · 0.20699 | 0 | 53 / 37 | 68/89 (76.4%) | 0.314 / -0.013 / 0.589 | 2.9185 |
| 388 | (11,13,19) · 0.33635 | 0.00012582067 | 40 / 34 | 66/90 (73.3%) | 0.538 / -0.055 / 0.275 | 7.5000 |
| 389 | (25,8,13) · 0.20699 | 0.0010336256 | 58 / 36 | 68/74 (91.9%) | -0.307 / -0.468 / 0.773 | 3.4667 |
| 390 | (12,13,19) · 0.33635 | 0.00018083764 | 57 / 37 | 65/94 (69.1%) | -0.292 / -0.038 / 0.362 | 1.1016 |
| 391 | (13,14,19) · 0.36222 | 0 | 41 / 38 | 66/94 (70.2%) | -0.067 / -0.257 / 0.227 | 1.0163 |
| 392 | (13,14,19) · 0.36222 | 0 | 37 / 39 | 67/79 (84.8%) | 0.967 / 0.869 / 0.922 | 0.6423 |

The full run took 34 min 40 s, including the starting field match. Total energy fell by 0.00002751952 to 92185.38323739906. At step 392, physical crystal energy was 0.3071531764135857 and gauge energy was 9.92227192134481 × 10⁻⁸. Inward current was 10.70771658845, with 66.75898974% dense and subsonic and a feed-to-bulk difference of 2.36 × 10⁻¹³. The first fully melted cross-section stayed at depth 0.3622238246. All 51,361 free zero-order sites passed the +10⁻⁹ order perturbation at the endpoint. The largest ratio of absolute gauge work to physical work was 1.614 × 10⁻¹⁰. [Full-precision cycles, sign tests, currents, certificate and checkpoint digest](../../../research/pde_ledger_v3/_measurements/mouth_r4_adaptive_cycles.json)

The force rose after every field update, including three crystal steps that individually reached the target. Persistent-label displacement cosines changed sign in the second half. The benchmark leaves open whether accelerating those labels or solving the coupled crystal and order response would be more effective. All six shape coordinates are held. The secondary gauge and stationary shape-force comparison await primary convergence.

**Independent review.** “Scoped clear for pausing the held primary at step 392. All ten cycles lowered energy and preserved whole-cell orientation. The post-rematch crystal force rose twice in cycles 6–10 and ended above the simultaneous target. The paired gauge and shape release remain gated.”
