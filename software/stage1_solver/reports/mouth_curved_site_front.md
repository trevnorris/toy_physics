# How far does one bound site move the conversion layer?

**One weakly bound site shifts the conversion layer all the way to the axis and changes the sideways pressure load on the held mouth.** At radius 1 and intake 0.25 we used the same archived carrier and deep conversion front with and without the site. The physical site strength, reciprocal change in carrier momentum and free-mouth acceleration remain open. [Empty mouth and order equation](mouth_free_order_depth.md) · [Earlier local-site test](mouth_order_slip.md)

**Extent.** The site favors order over a normalized Gaussian footprint of width 0.01 mouth radii on the material side of the wall. The geometry has three brane directions and the depth coordinate w. We differentiated the stationary equation with the carrier held fixed, using angular symmetry of the empty mouth to resolve the mean, dipole and quadrupole of a single localized site. A positive coupling of magnitude h enters the energy as −h times the order sampled by the site. The table gives the change in depth along the fourth coordinate of the half-order surface, per unit h, on the finest mesh. The wall column is an angular average. At linear order the axis and wall mean depend only on the angularly averaged part of the single-site response.

| Location of site | Distance to observed axis point | Shift at axis / h | Mean shift at wall / h |
|---|---:|---:|---:|
| Upstream, depth 1.175 | 1.00547 radii | +1.24873 | +1.26083 |
| At half order, depth 1.25195 | 1.00160 radii | +2.71019 | +2.82456 |

Giving the same empty shape the positive local restoring stiffness found upstream reduced the axis response for the half-order site to 2.20 × 10⁻⁷ per h. This control changes the linear restoring force and holds the original background by construction. Three meshes, a doubled footprint width and an independent calculation on the middle mesh support a response across a mouth radius. The absolute displacement scales with the weak coupling. These measurements leave the pointwise shift on the opposite wall undetermined. [Front calculation](../tools/mouth_curved_site_front.py) · [Displacement across the mouth](figures/mouth_curved_site_front.png)

**Force.** The neutral wall has zero normal derivative of order, so the order field exerts normal pressure there. Its angular imbalance on a curved wall gives a sideways force. On mesh 48 the external force required to hold the wall and its attached site together, per h, was **+0.276621** upstream and **+0.605166** at half order. Positive points toward the side occupied by the site. The force the material exerts on the assembly points the other way. The bare wall terms were +0.207901 and +0.324356. A volume balance including advection and the remote boundaries reproduced the combined loads to within 0.000021 per h. Opposing contributions from the collar and flare make that check important. These are forces from the order field and its site interaction on a held assembly with prescribed flow. [Wall and volume stress calculations](../tools/mouth_curved_site_stress.py)

**Which side wins?** We placed two equal-mass sites on opposite sides of this *one mouth*, with the upstream site on +x and the half-order site on −x. For site energy m μ(⟨χ⟩), the earlier adsorption calculation gives h = −m μ′ at the empty order averaged over each footprint. If the ratio of the absolute half-order slope to its upstream value exceeds **0.45710** on mesh 48, the net material force points toward the upstream site's side. The two proposed contact laws fall on opposite sides of this threshold.

| Contact law | Half-order / upstream slope magnitude | Net material force per one site's mass |
|---|---:|---:|
| All contacts weaken | 0.99995 | +0.01798 |
| Mean attraction persists | 0.26021 | −0.00715 |

Exchanging sides reverses the result. Opposed sites at equal depth cancel, as does a population distributed uniformly around each depth. These comparisons use first-order superposition and the common energy reference of the two adsorption models. An additional intrinsic site energy that varies with order could change the slopes and the direction. The signs survived the mesh, footprint and host-size controls. [Opposite-site calculation](../tools/mouth_curved_site_pair.py)

**0PN implication.** A collective change in order can transmit a lateral load even when a single site cannot be carried through conversion. To test gravity at fixed intake, compare the held mouth alone and in a weak companion flow. At each depth measure how attachment strength varies with angle, weighting it by the cosine toward the companion, and how the companion changes the pressure exerted per existing site. Equivalently, subtract the empty-wall load from the occupied-wall load in each flow, then compare those differences. The full momentum balance also needs the carrier and downstream flux. The present tests impose the site positions on an isolated held mouth. [Measurements and independent reviews](../../../research/pde_ledger_v3/_measurements/mouth_curved_site_front.json)

Working checkpoint `a91a633f`. Formal repository review is pending.
