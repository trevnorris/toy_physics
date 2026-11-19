# Can the two trapped-wave excitations be distinguished at fixed intake?

**Two peaks survive in the weakly damped example. They merge in the more strongly damped example.** The steady carrier solver gave the same mass current, 0.0373999726942, at the inlet, outlet of the throat and far boundary on two grids. The drive and release were calculated in a separate linear wall–wave model at the one-excitation reference geometry. Its coherent supporting wave stays on. Both damping examples use the same radius, depth, bend mass (12), stiffness (6) and wave coupling. Only the assumed bend damping changes. The axisymmetric current solver supplies a control. An angular pressure perturbation can have zero mean current change and still exert a bending force. That force remains unmeasured. [Carrier solver](../src/stage1_solver/open_sleeve.py) · [Trapped-mode preparation](mouth_trapped_wave_spectrum.md)

| Bend damping | Drive angular frequency | Wall displacement amplitude | Lag behind force |
|---:|---:|---:|---:|
| 0.1 | 0.456605832128 | 5.68290633648 × 10⁻⁶ | 179.857136506° |
| 0.1 | 3.93100103488 | 5.68290625433 × 10⁻⁶ | 179.857136512° |
| 6.5 | 0.456605832128 | 5.60972368175 × 10⁻⁶ | 170.793905472° |
| 6.5 | 3.93100103488 | 5.6097236005 × 10⁻⁶ | 170.793905376° |

Each row used a force amplitude of 0.001 in model units. Displacement is the component at the driven frequency where the bend profile is largest. Halving the force halved that component. At the exact sum and difference frequencies the two complex responses coincide. With both forces applied in phase, the trapped wave balances their load and the wall settles to rest. Switching both off exposes the stored wave. The dipole’s initial normalized amplitude was 0.00301370455753 in both preparations. [Drive and release](../tools/mouth_driven_bend.py) · [Two-tone cancellation](../tools/mouth_bend_pair_release.py)

| Bend damping | Trapped-wave energy decay rate from wall ringdown | Peaks in conditional spectrum |
|---:|---:|---:|
| 0.1 | 0.000177499200907 | 2 |
| 6.5 | 0.00713343083853 | 1 |

Counting peaks requires a readout law. We assumed equal integrated signals for two separately selected transitions, one cold loss channel, a common frequency shift and zero extra phase noise. This keeps the conservative angular-frequency spacing at 0.00406283876392 and gives widths proportional to the measured decay rate, by factors of one and three. The amplitudes, phases and decay rates above are numerical measurements in the stated wall model. A fresh independent calculation reproduced the wall responses, ringdown rates and peak counts. The draining material’s damping and its two-transition spectrum remain open. [Calculation, ringdown samples and independent check](../../../research/pde_ledger_v3/_measurements/mouth_driven_bend.json)

**Next physical input.** Impose a measured bend displacement on the draining carrier and integrate the opposing angular traction over the wall. Subtract the steady background and account separately for work by the supporting wave and the feed. The cycle-averaged work lost at matched intake would constrain damping independently of the two observed lines. [Distributed wall action and force](../../../research/pde/paper/pde.tex#L996)

Working base `df13acaa`. Formal review and commit are pending.
