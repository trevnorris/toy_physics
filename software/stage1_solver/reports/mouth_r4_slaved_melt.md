# Radius 4 · the melted label extension at fixed shape

At the two held states from `12e02d90`, the Stage 2 solver already eliminated the strictly melted label unknowns. For each trial of the active labels it solved $K_{mm}u_m=-K_{ma}u_a$ with the registered boundary conditions. The scalar gauge strength cancels from this equation. The earlier shape probes transported labels to the new depth grid without re-solving that equation. The fixed-geometry label search reached its stall before it used a shape probe. [Earlier stop](mouth_r4_eulerian_stage2.md) · [Static measurements and source records](../../../research/pde_ledger_v3/_measurements/mouth_r4_slaved_melt.json)

| Held-state check | Gauge 10⁻⁸ | Gauge 10⁻¹⁰ |
|---|---:|---:|
| Physical inverse volume at the blocking corner before the repeated harmonic solve | 1.78657 × 10⁻¹³ | 6.99635 × 10⁻¹³ |
| Same corner after the repeated solve | 1.77642 × 10⁻¹³ | 6.99635 × 10⁻¹³ |
| Gauge energy change on repeating the solve | +8.27 × 10⁻²⁵ | 0 |
| Largest change in active force including the harmonic gauge | 8.47 × 10⁻²² | 1.03 × 10⁻²⁵ |
| Consecutive final active steps below 0.0001 in the original search | 20 | 15 |
| Same corner after changing the numerical gauge | 0.320984 | 0.430024 |
| Barrier gauge energy / crystal energy after the change | 4.40703 × 10⁻⁸ | 4.95383 × 10⁻¹⁰ |
| Active step fraction accepted with the barrier | 1/128 | 1/128 |
| Reduced energy change on that step | −0.00174382 | −0.00138379 |
| Remaining free-label force after that step (target 0.000030) | 0.125078 | 0.113078 |

The repeated harmonic solve moved melt-only labels by at most 2.22 × 10⁻¹⁶. The changes in physical crystal energy, every active crystal displacement force and every flow reading were exactly zero. Across all accepted steps the archived unscaled harmonic residual stayed below 6.0 × 10⁻¹⁴. Re-solving the harmonic equation at both signed widths for every shape coordinate produced a negative determinant on one side of every pair at both gauges. These measurements locate the failure in the stationary harmonic extension itself. [Signed shape and active-label measurements](../../../research/pde_ledger_v3/_measurements/mouth_r4_slaved_harmonic_work.json)

The fallback adds $d-1-\log d$ to the unit numerical gauge, with $d$ the determinant of the spatial labels relative to the identity chart. It assigns one eighth of each fully melted cell’s physical volume to each corner. Both solutions passed a separate orientation test throughout every cell. Holding the active trace fixed while applying the barrier changed physical elastic energy and active crystal forces by exactly zero. The largest change in the order derivative at a zero-order node was +0.01199905 or +0.01546004. This is the derivative for *order returning* to that material. It has to be included when the order and carrier are next rematched. At positive-order nodes its change was exactly zero.

Both initial tests used exit radius 3.30637727 and depth 0.81914036. The saved fields carried inward current 10.70852202, with 66.75637% carried by dense, subsonic material and the sampled fully melted cross-section beginning at depth 0.28460443. The internal fields must settle at both gauges before releasing the shape. Measuring free translation requires a jointly balanced shape.

Independent review of the two fixed-shape tests concluded, “The two direct solves preserve the physical held state, the original active trials already used harmonic elimination, and the measured determinant barrier restores corner orientation and permits one energy-lowering active step at both gauges.”
