# Fixed-wall angular release of the radius-four crystal

At the previously verified depth 0.79590446 and exit radius 3.50122632, the ordered interior lowered its energy by a small polar–equatorial displacement. Both solid meshes resolved a local minimum along that coordinate. The wall, the bonded interface between core and collar, the receiving reservoirs and one quantum of coherent circulation were preserved. The material can move through the core–collar join. [Preparation and frozen-shape controls](mouth_joint_polar_shape.md)

| Fixed-wall measurement | 24 × 24 | 33 × 32 |
|---|---:|---:|
| Equilibrium internal amplitude | −0.00087021140 | −0.00087023495 |
| Peak extra radial displacement | 0.0000917801 | 0.0000917826 |
| Change in total energy | −0.0000013142194 | −0.0000013139852 |
| Signed angular energy derivative at the minimum | −1.34 × 10⁻¹¹ | −7.94 × 10⁻¹⁵ |
| Max. free meridional material force | 8.04 × 10⁻⁶ | 4.89 × 10⁻⁶ |
| Smallest sampled principal stretch before / after | 0.32919844 / 0.32919818 | 0.30429934 / 0.30429934 |
| Smallest volume ratio after release | 0.11845998 | 0.12705384 |
| Inward current / dense and total-subsonic share | 11.616869 / 70.4773% | 11.616844 / 70.4769% |

The material displacement is the amplitude times \(x^3(1-x)h(z)P_2(\cos\theta)\), where \(x\) is fractional radius. The depth factor vanishes beyond 62% of the depth. Its peak occurs at 34% of the depth, at three quarters of the local radius. The longitude-independent motion therefore vanishes at the exit rim. The smallest stretch on each mesh lies near that rim, at 97.2% and 98.3% of the depth respectively. The mesh-dependent difference in these samples predates this relaxation. Two meshes do not establish a continuum lower bound on stretch. [Material law and signed checks](../tools/mouth_fixed_wall_angular.py)

Changing the equilibrium amplitude by ±0.0001 or ±0.0002 raised the energy on each mesh. On coarse, fixing the meridional displacements and terminal fields gave an angular energy drop within 1.4 × 10⁻¹² of the re-equilibrated result. On fine the rechecked meridional field already satisfied the force criterion and its nodal displacements stayed fixed. The force from a further meridional perturbation of the angular correction agreed with a signed difference to 5.1 × 10⁻¹⁰ on the finer control steps. All mobile terminal order samples lay between 0.03792 and 0.93725. The stored shaft reference ranged from 0.5 to 1, with order 0.998997 on its mid-depth axis. These results apply at the verified wall geometry. [Guarded energies, force residuals and stretches](../../../research/pde_ledger_v3/_measurements/mouth_fixed_wall_angular.json)
