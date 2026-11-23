# Internal polar work on the radius-four drain

**The bonded solid can lower its energy at fixed wall shape.** At the same nonspherical wall, the generalized force on an independent internal quadrupole was **−0.0030204568842106005** on the 24 × 24 material mesh and **−0.003019839931944075** on the 33 × 32 mesh. A displacement of −0.0001 lowered the energy on both by about 2.85 × 10⁻⁷. The finer mesh used the transferred geometry and a fresh material and field relaxation. The front overlap was zero in both probes. [Guarded measurements](../../../research/pde_ledger_v3/_measurements/mouth_polar_material_work.json)

The crystal moves radially by `s x³(1−x) h(z) P₂(cos θ)`, where `x = r/b(z)` on the meridian before its polar deformation. The depth window spans 0.06–0.62 of the duct length. Both bonded pieces move continuously across their join. The displacement vanishes on the axis, wall and axial cuts. Its force is conjugate to the coefficient `s`. The calculation evaluates each moved lattice site's occupancy in the fixed physical order field and checks that the probe has zero overlap with the mobile terminal order. [Virtual-work calculation](../tools/mouth_polar_internal_work.py)

| Accepted state | Wall quadrupole amplitude | Analytic `dE/ds` | Centered `dE/ds` at step 0.0001 |
|---|---:|---:|---:|
| Coarse spherical control | 0 | 3.541871657075646 × 10⁻¹⁵ | 6.147500036329001 × 10⁻¹² |
| Coarse first nonspherical state | 0.00016963940257820912 | 0.0029318856002603835 | 0.0029318898103246776 |
| Coarse last accepted state | 0.0001743942780611947 | 0.0030204568842106005 | 0.003020456891081834 |
| Fine mesh at the same exterior | 0.0001743942780611947 | 0.003019839931944075 | 0.0030198399367779146 |

The accepted coarse and fine material residuals were **1.1218360803789962 × 10⁻⁵** and **4.709288860968228 × 10⁻⁶**. On the subsequent coarse trial, 70 alternating sweeps ended at residual **0.0006830708042412943**. The subsequent fine trial stopped during material polishing at **0.0008128832206617774**, about 27.1 times the joint tolerance of 0.00003. That iterate exceeded both the joint tolerance and the tighter polishing target of 0.00001. Intermediate polishing residuals were unrecorded. Whether the wall and independent internal polar displacements can equilibrate together remains open. [Accepted states and stopping residual](../../../research/pde_ledger_v3/_measurements/mouth_polar_material_work.json)

The receiving-side order was also checked separately. Re-solving its native order and exterior crystal terms at wall amplitudes −0.025 and +0.025 gave signed slopes within 1.7 × 10⁻⁷ of the corresponding fixed-field derivative. At the historical held shape, a separate feasible receiver-order variation gave driven carrier work **+2.752109651282808 × 10⁻⁶** per unit amplitude at step 0.0001. The native receiving equilibrium therefore specifies the reference order used by the carrier. Full stationarity of that field would also require its carrier contribution. [Receiver variations](../../../research/pde_ledger_v3/_measurements/mouth_polar_material_work.json)
