# Can flow hold a constituent inside a converting mouth?

**One candidate surface confined a cold constituent to the conversion region while brane flowed past. A second surface let it through at the same intake, launch speed and drag rate.** This supplies a possible way to retain material near a mouth. The calculation held the guide fixed and measured motion along its surface. The sideways force on a free mouth remains open.

**Material measurement.** We placed one constituent above a cold triangular sheet with the pitch and energy scale of the independent neutral compression test. We held the area and host plane fixed, allowed the host sites to move in the plane, and relaxed the added site in all three local coordinates. We subtracted the empty host energy at each order. Both trial laws use the same bonds, reduce the lateral variation in their potential by the same factor at fixed height and host positions, and set the site energy to zero far from the sheet. They differ in what happens to the average face attraction when the material loses shear order. [Surface calculation](../tools/mouth_surface_site_transfer.py) · [Neutral calibration](neutral_patch_zero_load.md)

| Loss of shear order from 0.875299 to 0.5 | Reversible work per unit site mass | Half-order / upstream mean normal holding force |
|---|---:|---:|
| All surface contacts weaken | +0.02053739 | 0.571233 |
| The unregistered fraction retains the spatially averaged attraction | +0.01146821 | 1.000000 |

The work includes host relaxation on the sampled stable adsorption branch. Its order integral agreed with independent endpoint energies. Enlarging the host from 25 to 64 sites gave +0.02053755 and +0.01147122. The force ratio refers to a separate weak probe scanned uniformly at height 1.1 microscopic lengths above a clamped plane at fixed site density. Scrambling the lateral positions redistributed the local force while preserving its area average for unchanged pair potentials. The earlier single-well values of +0.00665233 and zero describe a different binding environment. An intrinsic energy per constituent that depends on order would contribute to the transfer work and leave fixed-order force scans unchanged. [Single-site energy reference](mouth_unpinned_feature_passage.md)

**Release into the same flow.** We used the radius-one, intake-0.25 carrier from depth 1.175 to 1.252010. One longitudinal coordinate follows wall arclength $s$. The site remains on the equilibrated adsorption branch with specific energy $μ(s)$. In this trial the stream exerts the linear drag $γ(u-v)$ per unit mass. The rates below are examples in inverse model time. The constituent starts at the upstream carrier speed for every launch. [Carrier readout](../tools/mouth_bound_site_carrier_work.py) · [Cold release](../tools/mouth_bound_site_retention.py)

| Trial rate $γ$ | All contacts weaken | Mean attraction persists |
|---:|---|---|
| 0 or 3 | Returned to entry | Returned to entry |
| 5.5 | Confined between the two cuts | Crossed half order in 3.69874 model time units |
| 12 | Crossed half order | Crossed half order |

The wall phase satisfies $u=dθ/ds$. The effective mechanical energy $H=v²/2+μ−γθ$ obeys $dH/dt=−γv²$. For the weakening-contact law the independently relaxed midpoint lies below both endpoint effective potentials whenever $4.53345 < γ < 6.33741$ on mesh 34. This guarantees an interior minimum in that interval of rates. At 5.5 the actual launch fell below both escape energies by 0.00008300. Its energy-balance error was below 2.0 × 10⁻¹¹ at that event. Both carrier meshes, both path resolutions and the larger material cell gave the same outcomes. Confinement here assumes a cold, stationary, one-dimensional guide. The calculation did not release the host or test escape by detachment.

**Energy and next measurement.** In the crossing at 5.5, the prescribed stream delivered 0.0111011 specific work. A reciprocal drag would also produce at least 0.0129854 of slip heat, so the flow would have to supply at least **0.0240865 per unit participating mass**. This is a conditional energy demand. The unloaded carrier’s equal-frequency reservoirs supply zero net steady power, and actual local transfer to a site has yet to be measured. The next tests are a fixed-height normal-force scan at the two orders and an independent longitudinal slip-decay measurement in a uniform-order sample with a known time scale. A loaded-versus-empty energy balance at measured participating mass flux would then test whether the flow can sustain retention or passage before using either to predict free fall. [Work and heat check](../tools/mouth_bound_site_work_readout.py) · [Reservoir balance](mouth_bend_flow_reaction.md) · [Measurements and review](../../../research/pde_ledger_v3/_measurements/mouth_bound_site_transfer.json)

Working checkpoint `a91a633f`. Independent local scientific reviews cover the surface work, force comparison and conditional release. Formal repository review is pending.
