# Order and transverse-wave motion in the held, draining waist

We let the crystalline fraction vary across and along each held throat. The entrance remained ordered, the side and terminal half-order surfaces stayed in their registered positions, and the carrier passed only through the terminal disk. Both reservoirs, material parameters and the initial degree-four action \(J=1\) were common to both radii. The static comparison re-solved the density, phase, interior order and occupied outgoing shear mode for each profile. The separate exterior-order work is needed for total radial and longitudinal force. [Numerical records](../../../research/pde_ledger_v3/_measurements/mouth_live_waist_order.json) · [Source of the held geometries](mouth_sealed_shape_relax.md)

The interior order energy is \(a_B\int[\delta_B^2|\nabla\chi|^2+\chi^2(1-\chi)^2]dV\), with \(a_B=75\) and \(\delta_B=0.04\). For a stationary order, eliminating the Route A momentum gives the reciprocal carrier term \(-nI(\nabla\theta\cdot\nabla\chi)^2/[2(n+I|\nabla\chi|^2)]\) in addition to the ordinary phase energy. Both wave kinetic and strain coefficients depend on \(\chi\). The derivative of \(J\,\mathrm{Re}\,\omega\) supplies the cycle-averaged wave reaction, with an equal angular occupation over the degree-four multiplet to keep an axisymmetric order field. The loss of crystalline order at a value of 0.5 leaves finite transverse wave coefficients in this law. [Field and order solver](../tools/mouth_live_waist_order.py)

| Entrance and waist state | Fine-grid center order, flow only | Center order with \(J=1\) | Interior plus drive energy derivative for widening | Same derivative for increasing log depth |
|---|---:|---:|---:|---:|
| 1, broad | 0.9996614 | 0.9996618 | +4.59555 | +2.62652 |
| 1, intermediate | 0.8978337 | 0.5177808 | +1.18481 | +0.19441 |
| 1, tight | 0.5017517 | 0.5003962 | −0.22459 | −0.51938 |
| 2, broad | 0.9996977 | 0.9996982 | +13.09977 | +10.86736 |
| 2, intermediate | 0.7221515 | 0.7221205 | +6.86235 | +4.76415 |
| 2, tight | 0.5000018 | 0.5000013 | +4.35743 | +3.08833 |

The columns use signed changes in the middle width, \(b(t)\mapsto b(t)+4t(1-t)q\), or in \(\log L\) at unchanged radius in material coordinates. We solved at \(q=\pm0.001,\pm0.002\) for both variations. The middle parameter has length units. It displaces the interior waist while both opening radii stay fixed. To obtain the total mouth load, add the outside order work with the same sign convention and take the negative energy derivative. The field was initialized from the ordered preparation and the occupation raised to one. A radius-2 mode stays concentrated on the downstream side even as the waist reaches half-order.

The narrow radius-1 profiles develop a connected patch where the ordered fraction sits at its lower bound. The wave pencil loses both mass and stiffness on completely uncovered wave degrees of freedom. We removed those degrees and continued the occupied mode on its positive support. The fine mesh has 388 and 269 free order nodes at zero at the intermediate and tight configurations. The counts of completely inactive wave nodes are 4 and 7. Increasing the order over each entire zero patch at two small amplitudes cost energy. Supporting checks also looked at individual sites and halves of the patch. These constrained stationary candidates permit conditional shape work. None of these perturbations establishes stability in every material direction. [One-sided field checks](mouth_live_waist_fields.md) · [Independent bounded checks](../../../research/pde_ledger_v3/_measurements/mouth_live_waist_independent_bounded.json)

The two order meshes used 48 and 72 radial intervals, with 14 and 20 intervals in the corresponding carrier and wave assembly. Fine-grid interior derivatives are shown above. The tight radius-1 case changes sign when exterior wall work is included. One original launcher disappeared after admitting the coarse tight job. Its child returned exit zero with all 18 records and the native resource-limit evidence. The archive records the missing parent exit separately. Both fine radius-1 candidates have native zero guard and worker exits. [Combined force and mesh comparisons](mouth_coupled_waist_balance.md)

### Can the order and wave move together in real time?

The time control retains one waist-centered trial displacement \(\chi=\chi_*+a(t)\eta(r,w)\) and one standing-wave coordinate \(Q(t)\), initialized with either all strain or all kinetic energy. The basis \(\eta\) vanishes on the pins and at a saturated site. The carrier's prepared density and phase are held for this short-time projection. The fixed background carries a steady current and still exerts reciprocal order work. At fixed phase the projected inertial density is \(I n/(n+I|\nabla\chi|^2)\). The velocity-linear coefficient is a total derivative for this one order coordinate on a stationary background. The normalized wave mass and stiffness vary as \(1+ma\) and \(\omega_*^2(1+ka)\), respectively, from the measured kinetic and strain overlap with \(\eta\). In the gauge removing the total derivative the conserved Hamiltonian is
\[
H=\frac{p_a^2}{2M_a(a)}+E_{\rm order}(a)+E_{\rm carrier,held}(a)
 + \frac{P^2}{2(1+ma)}+\frac{\omega_*^2(1+ka)}{2}Q^2 .
\]
We integrated both canonical pairs for two shear-wave periods. Energy exchanges between the wave and order. The initial action sets the preparation only. Open-domain leakage was measured separately on the actual mode. The real-time control uses a conservative projection of that mode with frozen spatial eigenfunction, azimuthally averaged stresses, and held geometry. It does not evolve the full carrier, mixing between wave modes, all order displacements or the moving wall. [Hamiltonian and timestep controls](../tools/mouth_live_waist_transient.py)

| State | Smallest to largest center order, both initial shear phases |
|---|---:|
| Radius 1, broad, refined | 0.99963512–0.99968832 |
| Radius 1, intermediate | 0.46935505–0.56596435 |
| Radius 1, tight | 0.48221669–0.51846528 |
| Radius 2, broad | 0.99969602–0.99970080 |
| Radius 2, intermediate | 0.71809350–0.72543745 |
| Radius 2, tight, refined | 0.49455870–0.50544149 |

The coarse runs gave projected order frequencies of approximately 123–296, compared with shear frequencies of 5.86–9.15. Holding wave action instantaneously fixed during a fast order oscillation would change the experiment. At the tighter time tolerance the largest energy drift was less than \(7.7\times10^{-14}\) of the prepared energy and the two tolerances agreed on the order trace to within \(5.2\times10^{-6}\) of its amplitude. At both tight radii the center remained near half-order through two periods. On the two narrower radius-1 states this time projection held the saturated sites at zero throughout. Reactivating those sites would require a larger space of order motions and a dynamical rule at the endpoints.
