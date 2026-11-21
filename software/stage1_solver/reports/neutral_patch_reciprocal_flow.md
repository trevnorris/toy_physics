# Momentum received by wall-facing sites from a companion flow at 0PN

The held rim gained force toward the companion in all twelve fine-grid comparisons at eight mouth radii. The wall-facing sites could receive momentum in either direction. We measured that force as incoming shear minus the change in fluid momentum around the sites. Both layers returned their drag reaction to the same local fluid. [Calculation](../tools/neutral_patch_reciprocal_flow.py) · [Guarded measurements](../../../research/pde_ledger_v3/_measurements/neutral_patch_reciprocal_flow.json)

| Fluid at patch entry | Material behind wall-facing sites | Force on those sites, central / layered bonds | Force reaching held rim, central / layered |
|---|---|---:|---:|
| Flow already established | Zero shear traction | +0.001125 / +0.001156 | +0.001652 / +0.001646 |
| Flow already established | Held at rest | −0.000218 / −0.000146 | +0.000743 / +0.000729 |
| Initially at rest | Zero shear traction | −0.000268 / −0.000213 | +0.000324 / +0.000314 |
| Initially at rest | Held at rest | −0.000326 / −0.000263 | +0.000298 / +0.000288 |

Forces are paired minus isolated for one equatorial patch, in microscopic binding-energy units per preferred site spacing. Positive points toward the companion. Mouth radius 1, intake 0.25, grip and incident flow history were fixed. Each result averages four thermal histories over one transit from 95% to 1% shear order. Doubling separation left 0.249348–0.249367 of the rim increment. [Incident field and preparation](neutral_patch_companion_rim.md)

**What the sites did.** With central bonds and established flow, tagged site 5 moved toward the companion in both cases. Against a shear-free backing its incremental speed was 0.0008781 versus 0.0009676 in adjacent fluid. Against a resting backing the readings were 0.0003684 versus 0.0003481, in patch lengths per model time. Shear alone gave the wrong force sign in the first case. The surface contribution was −0.000898 and the fluid momentum rate was −0.002023, giving +0.001125 to the wall-facing layer. Two fluid regions and the site reactions gave the same load.

**Material assumptions.** The archived curved-mouth field supplies the drive at the outer edge. The added viscous fluid resolves the layer normal and one transverse coordinate, averaging along the force direction. Density, kinematic viscosity and site drag are each one in patch units. A patch length is 0.01 mouth radii. The bonds and rim connection remain elastic. This local fluid omits elastic memory, individual three-dimensional wakes and fluctuating solvent stress. The established state starts from the empty-channel steady flow for each backing. The second preparation starts all fluids at zero sideways velocity. [Transport and controls](../tools/neutral_patch_shear_channel.py) · [Related fluid–structure methods](https://arxiv.org/abs/1009.5648)

Finite backing resistances of 0.5 and 4 also gave opposite wall-layer force signs with established flow. An independent calculation predicted the signs at the limiting backings for both bond laws. Grid refinement shifted those rim readings by less than 0.8% and wall-layer readings by up to 7.8%. A weak shear pulse with traction and velocity measured behind the sites could determine the material response independently. Carrier mass per volume and single-site drag also need calibration. Free acceleration further requires active coverage, the remaining momentum channels and translating-mouth inertia.

Working base `49dca23e`. Independent local review covered the reciprocal force and its scope. Formal review and commit are pending.
