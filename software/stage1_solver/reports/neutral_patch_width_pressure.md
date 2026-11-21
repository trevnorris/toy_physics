# When material leaves a brane patch, does the slab thin or does its density fall?

**Both outcomes are possible.** In two candidate materials, losing particles without squeezing the faces made the slab thinner and lowered its depth-mean number density. Squeezing the faces enough to produce the *same loss* made it thinner and denser. A combination of the two loads kept the density almost constant. These are measurements on neutral materials outside a throat.

**Pressure over a growing region.** We used two settled 1,152-site bilayers with different bond laws, free to move in two lateral directions and in depth. We applied equal and opposite pressures to the two faces over a stripe, allowing sites to move between the loaded region and its surroundings. The pressure did work on the volume between the actual faces. In a fixed central window we measured their separation, their midpoint, and the particle number. Density was the particle number divided by the enclosed volume. Under a weak squeeze over each finite stripe both materials thinned, lost some local particles, and became denser. Reversing and halving the loads, shifting the stripe, and an independent force calculation checked the response. [Pressure experiment](../tools/neutral_patch_width_pressure.py) · [Material and wrinkle preparation](neutral_patch_wrinkle.md)

| Width of the pressure stripe (Gaussian standard deviation in bond-length units) | 1 | 2 | 4 | Whole cell |
|---|---:|---:|---:|---:|
| Central thinning relative to the whole-cell squeeze, first bond law | 67.7% | 87.9% | 96.7% | 100% |
| Second bond law | 56.0% | 78.6% | 94.0% | 100% |
| Expected from averaging a purely local linear response with the same sensor | 70.7% | 89.4% | 97.0% | 100% |

These finite cells resisted even a uniform change in width. Much of the observed size trend comes from how much of the central window was under full pressure. The globally conserved site count allows lateral redistribution. This strip test supplies no exchange rate with a bulk bath or radial falloff exponent. [Size and amplitude checks](../tools/neutral_patch_width_readout.py)

**Equal column losses with a reservoir.** A separate uniform-patch test let both the width and the number of sites per projected area relax at a common far chemical potential. We matched a column loss of 50 parts per million in logarithmic units. The values below are changes in natural log times one million. The two preparations are stretched sheets kept continuously covered. Vacuum has lower grand potential, so these results apply to the local material branch.

| How the same material loss was produced | Thickness change, first / second law | Mean density change, first / second law |
|---|---:|---:|
| Face compression | −99.99 / −91.66 | +49.99 / +41.66 |
| Raising the local energy per particle, with no face pressure | −14.57 / −13.70 | −35.43 / −36.30 |

Separately measuring the two responses predicted a combination of face compression and higher local energy per particle that would keep density constant to first order. Across the two laws, halving both loads without refitting left a density change of only **0.0206% or 0.0220% of the column change**. Choosing this mixture demonstrates a possibility. Its occurrence in a drain still requires an independent mechanism. [Measurements, independent scripts and review scopes](../../../research/pde_ledger_v3/_measurements/neutral_patch_width_pressure.json)

**What this gives 0PN.** In the existing slowly varying, fixed-intake carrier model, prescribing a 1% narrower annulus at eight radii lowers its density by a further 9.81 × 10⁻¹⁰ in bath-density units. That calculation provides only the carrier response. The physical slab’s thinning profile and its force on a second mouth remain open. The missing inputs are the net normal load from both sides of each face and the energy change that drives exchange of constituents at that same intake. Their measured ratio would determine whether a drain made of either candidate material follows the nearly constant-density path. A freely draining calculation must measure width and depth-resolved density with the imposed squeeze removed. [Active-drain force measurement and width equation](active_drain_width.md)

Working base `1dc45d09`. Independent local reviews are recorded with the measurements. Formal project review and commit pending.
