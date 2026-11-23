# Radius-four conversion layer coupled to both reservoirs

A dense, subsonic annulus survives when the conversion layer and the two remote reservoirs respond to each other. The earlier five-coordinate balance shifts even with the meridional drain drive turned off.

A single conforming mesh connects the reservoirs through the disk. Density amplitude and meridional phase vary throughout. The crystalline order relaxes within 0.12 on either side of the disk, with its rim and the outer edges of that band matched to the material preparation. We retain the ordered interior, a sealed lateral wall, remote phase difference 0.4510084283 and one externally prepared quantum of circulation through both fluids. The comparison removes the entire sharp carrier energy, including the sharp terminal contact, and adds the joined carrier Hamiltonian plus the change in order and elastic energy within the band. The skeleton is held at the indicated mechanical state for each test. [Carrier/layer Hamiltonian](../tools/mouth_supercore_joined_front.py) · [Native receiving order](../tools/mouth_supercore_terminal_order.py)

| Geometry | Depth / exit | Wound / unwound intake | Intake retained | Dense, subsonic share¹ | Inside / receiving axial stress² |
|---|---:|---:|---:|---:|---:|
| Registered | 0.81583 / 3.57310 | 11.95131 / 14.09377 | 84.80% | 72.55% | +0.97391 / +0.39908 |
| Released iterate 203 | 0.79368 / 3.49430 | 11.64120 / 13.81086 | 84.29% | 68.96% | +1.08835 / +0.44067 |
| Released iterate 309 | 0.79332 / 3.49473 | 11.64411 / 13.81350 | 84.30% | 68.96% | +1.08765 / +0.44054 |

¹ Positive incoming current in the last inside cell layer, with density above half bath density and total speed including swirl below local sound speed, on the common 13 × 10 × 14 mesh. At iterate 203 the integrated feed and receiving currents agreed to 2.51 × 10⁻¹¹. The iterate had sufficient material order for a separate energy-lowering displacement below the conversion band. Iteration 309 is the last accepted state before a subsequent mechanical trial failed.

² Integrated on opposite sides of the same exact disk and differenced against zero meridional bias at the same winding. These stresses include the coherent carrier and crystalline order. The separate skeletal elastic stress belongs to the material calculation. On the registered state the raw disk currents from the two sides changed from 12.06569 / 10.81532 to 11.98879 / 11.08548 upon mesh refinement to 18 × 13 × 24. The corresponding stresses changed from +0.97391 / +0.39908 to +0.97442 / +0.44248. A unique pointwise contact traction remains unresolved despite weak current conservation. [Disk-trace calculation](../tools/mouth_joined_front_trace.py)

| Geometry | Five shape-force changes from replacing the sharp interface | Largest part induced by drain bias |
|---|---|---:|
| Registered | +2.9219, +1.4401, -5.7377, -4.4239, -2.2684 | 0.00760 |
| Released iterate 203 | +3.7117, +2.0711, -7.9922, -6.0998, -3.2788 | 0.00707 |
| Released iterate 309 | +3.7219, +2.0774, -8.0192, -6.1179, -3.2936 | 0.00711 |

The signed force differences largely come from allowing the terminal order profile to adjust at rest. At iterate 203 the original wound line search lost precision before returning a field. Integrating the directional gradient along small trial steps reached the original stationarity tolerance using a previously converged wound field as the initial guess. The successful unwound control was reused. The independent review checked the disk measure, both field gradients, pressure signs and conserved currents. [Matched force audit](../tools/mouth_joined_front_audit.py) · [Stationarity check](../tools/mouth_joined_front_stationarity.py) · [Independent review](mouth_collar_continuum_review.md) · [Verified measurements and source hashes](../../../research/pde_ledger_v3/_measurements/mouth_supercore_joined_front.json)
