# Radius 4 · Paired internal equilibrium at step 395

**Both gauges passed the simultaneous internal residual target of 0.000030 with all six shape coordinates held.** The primary took three crystal steps from step 392. After the first field update its peak was 0.00003257692. Step 394 encountered nonpositive CG curvature and accepted an energy-lowering step at 1/64. Step 395 brought the matched peak to 0.00001691319. At the secondary gauge the primary endpoint satisfied the native tests after a field update with zero change in energy at that gauge and zero additional crystal steps. The benchmark used the exact Hessian and held-shape solver from `63db8cc0`. [Starting state and force history](mouth_r4_adaptive_cycles.md)

| Held endpoint quantity | Gauge 10⁻⁸ | Gauge 10⁻¹⁰ |
|---|---:|---:|
| Maximum internal residual | 1.69131862 × 10⁻⁵ | 1.69060531 × 10⁻⁵ |
| Constrained order residual | 3.00133 × 10⁻¹¹ | 3.16068 × 10⁻¹¹ |
| Maximum carrier residual | 4.01341 × 10⁻⁹ | 4.01341 × 10⁻⁹ |
| Total energy | 92185.38323628501 | 92185.38323618693 |
| Physical crystal energy | 0.30715233478933646 | 0.30715233478933646 |
| Gauge energy | 9.90719914 × 10⁻⁸ | 9.90719914 × 10⁻¹⁰ |
| Gauge energy / crystal energy | 3.22550019 × 10⁻⁷ | 3.22550019 × 10⁻⁹ |
| Certified global lower bound on physical inverse determinant | 1.20152618 × 10⁻⁵ | 1.20152618 × 10⁻⁵ |
| First fully melted cross-section, depth | 0.3622238246 | 0.3622238246 |
| Positive work for +10⁻⁹ at free zero-order nodes | 51,359 / 51,359 | 51,359 / 51,359 |
| Largest absolute gauge work / physical work in that test | 2.25008 × 10⁻¹⁰ | 2.25008 × 10⁻¹² |
| Inward current | 10.70771658766 | 10.70771658766 |
| Dense, subsonic share including swirl | 66.75898975% | 66.75898975% |
| Absolute feed-to-bulk current discrepancy | 3.19744 × 10⁻¹⁴ | 3.19744 × 10⁻¹⁴ |

The shape, crystal labels and all three saved field arrays agree bit for bit across gauges. The maximum order difference and count of changed order-support entries are both zero. The physical orientation certificate covers all 28,080 cells. The integrated weak exit-disk current differs from the feeding current by 1.24 × 10⁻¹⁴. [Full-precision energies, both stationary states, signed work and flow](../../../research/pde_ledger_v3/_measurements/mouth_r4_held_pair395.json)

### Forces on the six held coordinates

Positive force favors increasing the listed coordinate. The comparison uses the original Checkpoint 10 measurements. [Checkpoint 10 force record](../../../research/pde_ledger_v3/_measurements/mouth_r4_checkpoint10_boundary.json)

| Shape coordinate | Checkpoint 10 | Step 395, 10⁻⁸ | Step 395, 10⁻¹⁰ |
|---|---:|---:|---:|
| Log depth | +0.150564883 | -0.006656262 | -0.006656253 |
| Scaled log exit radius | -0.641236402 | -0.638270855 | -0.638270849 |
| Mean width | +0.705499333 | +0.660848084 | +0.660848133 |
| Depth asymmetry | -0.666834458 | -0.527998742 | -0.527998782 |
| Waist profile | +0.465602934 | +0.618248729 | +0.618248764 |
| Polar–equatorial distortion | +0.009007417 | +0.006866546 | +0.006866547 |

The largest cross-gauge force difference is 4.866 × 10⁻⁸. The maximum discrepancy against central energy differences at the two probe widths is 4.439 × 10⁻⁹. Widening the middle and narrowing the exit both lower the energy locally. The held geometry has depth 0.81914036 and exit radius 3.30637727. Shape relaxation and free release remain separate, unperformed tests.

**Independent review.** “Scoped clear for paired internal stationarity at the held radius-4 shape. Both target gauges pass native force, work-sign, flow and whole-body orientation checks at step 395. All six boundary forces retain measurable signed work, so shape equilibrium and free release remain open.”
