# How much of a held bend the draining material follows at 0PN

**Fresh material partly follows the bend and exerts a restoring force on the wall.** We tested one curved mouth of radius 1 at inward intake 0.25. The entrance stayed fixed while we displaced the end of its first bend sideways by ±0.000005 mouth radii. The strip used a measured axial material path and evolved the sideways motion of its sites. An individual axial tag took 6.92032503 model time units from depth 0.5 to mid-depth 1.175, compared with the earlier cohort mean of 6.92036821. [Axial path](../tools/mouth_entrance_tag_path.py) · [Material and contact equations](../tools/mouth_entrance_tag_transport.cpp)

| Where fresh wall-contacting material crossed | Mean sideways motion as a fraction of the imposed end displacement |
|---|---:|
| 16 local spacings upstream of mid-depth | 12.6697% |
| Upstream edge of the force gauge | 22.2259% |
| Mid-depth | 27.7018% |

Each entry averages 16 fresh crossings at that location over arrivals 313–320, using the difference between opposite steps. Relative to an unbent run at the same intake, the material exerted **−0.01255338958** on the mid-depth wall for the positive step and **+0.01255338960** for the negative step, in patch force units. The percentages above measure displacement. The force is the reaction integrated across the local gauge.

**Controls.** All seven guarded runs completed, and both consecutive eight-arrival windows gave a restoring force in each. In the last window, the positive-step force per unit displacement differed from baseline by +0.648 parts per million at half amplitude, +1.741 with half the outer time step, and +0.015 with twice the transverse period at fixed coverage. Both alternative near-centered inlet preparations changed it by less than 0.002 parts per million. Extending the outlet and measuring later (arrivals 361–368) changed it by −3.243 parts per million. The two baseline windows differed by 1.022 parts per million. The seven earlier fixed-step attempts left the calibrated contact range and supplied no force measurement. The diagnostic located the baseline failure in unforced initialization. [Guarded outcomes, literal readings and independent reviews](../../../research/pde_ledger_v3/_measurements/mouth_entrance_tag_release.json)

**What this establishes for free fall.** We separately released the entrance and both bends in a weak oscillating brane stream. Both model walls had the same axial geometry, intake, zero trapped-wave occupation and first-bend effective tension 0.02. Only the deeper effective tension changed. Each gave the same first-bend area holding force at the prescribed displacement. The free response includes the carrier and the assumed wall support. The local strip reaction measured above was evaluated separately. [Release equations](../tools/mouth_entrance_release_gate.py)

| Whole-mouth result on mesh 34 | Deep tension 0.02 | Deep tension 0.04 |
|---|---:|---:|
| First-bend area holding force at +0.000005 radii | +1.99657 × 10⁻⁷ | +1.99657 × 10⁻⁷ |
| Steady pull per incident speed | 0.2500039295 | 0.2500039317 |
| Slow entrance inertia with the finite-bath spring removed | 15.57260 | 15.31017 |
| Free entrance acceleration per incident speed, in phase at angular frequency 0.005 | 0.01610712 | 0.01640658 |

Doubling deeper support increased free acceleration by **1.85922%** on mesh 34 and **1.85940%** on mesh 26. These two examples establish that the first-bend measurement alone leaves the acceleration undetermined. The next independent input is the **total wall force and the motion of the entrance, middle and flare under a weak sideways load on the deeper bend** at the same intake. A material law calibrated independently across both bends could also supply that input. Once the model reproduces the independent response, use it to predict the flexible-mouth release, then compare accelerations across radii. The earlier rigid-mouth radius comparison is a separate result. [Rigid mouth and radius dependence](mouth_carrier_0pn.md) · [Earlier flexible-mouth tests](mouth_free_bend_response.md)

Working base `49dca23e`. Independent local reviewers checked the completed strip measurement and the conditional two-depth release. This working checkpoint preserves the local review scopes. Formal review is pending.
