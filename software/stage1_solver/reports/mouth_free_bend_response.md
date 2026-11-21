# Does a draining mouth fall together when its bends can move?

An isolated draining mouth can accelerate along a weak brane stream with its entrance and deeper flare free to respond. Raising trial wall tension from 0.02 to 0.04 gave essentially the same held force and an **18.10% higher in-phase entrance acceleration**. The written material model has yet to select its physical lateral support. [Verified measurements and review scopes](../../../research/pde_ledger_v3/_measurements/mouth_free_bend_response.json)

**Preparation.** One radius-1 flare drains at intake 0.25 from the feeding slab down the fourth coordinate into the bulk. Translation and two sideways bends are free. The axial and radial profile is prescribed, with zero bare wall mass and zero trapped-wave occupation. Both meshes have bath radius 9 and roof 10. Mesh 34 has a wider displacement interpolation and a smaller pulse timestep. All readings use model units.

**Release without wall tension.** The full linear field and wall equations admit an unforced displacement proportional to the real part of $v e^{(0.12162116+0.03582566i)t}$. The amplitude grows by a factor of $e$ every 8.2223 time units. Mesh 26 gave 8.2224. Opposite weak pushes, including half strength at the entrance, confirmed growing relative displacement after the push ended at time 8. At time 32 the two meshes passed the work check to within $2.63\times10^{-9}$ of cumulative absolute input work. Two earlier runs to time 120 failed the energy check. Cancellation in the indefinite energy is a possible numerical cause. Late readings from those attempts are excluded. [Full field and wall equations](../tools/mouth_free_bend_release.py) · [Signed-pulse controls](../tools/mouth_free_bend_certificate.py)

**Supported mouth.** Adding tension times flared wall area makes the full quadratic energy positive for the represented carrier field and three wall coordinates at either trial tension. Separate experiments measured the holding force in a steady brane stream, the push response in resting material and the free acceleration in an oscillating brane stream. The bulk reservoir was undriven. [Held force](../tools/mouth_free_bend_held_flow.py) · [Entrance inertia](../tools/mouth_free_bend_inertia.py) · [Free incident drive](../tools/mouth_free_bend_incident.py)

| Mesh 34 reading | Tension 0.02 | Tension 0.04 |
|---|---:|---:|
| Force on held entrance per incident speed | 0.25000393 | 0.25000399 |
| Slow entrance inertia, after subtracting the finite-bath spring | 15.57260 | 13.37587 |
| Held force coefficient divided by slow inertia | 0.01605409 | 0.01869067 |
| Direct entrance acceleration per incident speed, in phase at angular frequency 0.005 | 0.01610712 | 0.01902300 |
| Direct flare acceleration per incident speed, in phase at the same frequency | 0.01639548 | 0.01913612 |

The holder applies the opposite force. Independent holder calculations agreed within $10^{-13}$ on each mesh. Constraining rigid motion gives inertia 12.26045. A small finite-bath centering stiffness limits extrapolation to zero frequency. Residual towing drag has unresolved sign.

**Stopping point.** The missing input is the lateral force transmitted between entrance and flare during drainage. Apply two independent small bend loads on the isolated mouth. Read displacement and phase at the entrance, mid-flare and deeper flare as material is replaced. A calibrated material stress, attachment coverage and renewal law could also supply this response. Fix this response independently, then test how the physical geometry, intake and total inertia vary with radius before predicting an unequal pair. [Material evidence](neutral_patch_neighbor_turn.md) · [Two-mouth criterion](../../../research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md#what-two-inward-drains-must-reproduce-in-the-newtonian-limit)

Working base `49dca23e`. Two independent local reviewers examined the scopes listed in the measurements. Formal review and commit are pending.
