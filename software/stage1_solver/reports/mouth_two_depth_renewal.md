# Force transmitted by renewing attachments on one draining mouth

At 0PN the brane and attachment equations transmit a measurable sideways force between the entrance and deeper wall of the cubic mouth. The force remaining after a push depends on whether the contacts carry their relaxed configuration downstream. Both tested laws still admit an unforced growing bend. [Measurements, guard outcomes and review](../../../research/pde_ledger_v3/_measurements/mouth_two_depth_renewal.json)

**Preparation.** The radius was 1 and the inward intake was 0.25. Separate sideways pushes acted at depths 1.175 and 2.25. The entrance and both bends moved freely. We used the saved deep shear-order profile, unit ordered-brane shear modulus and unit tangential grip from the earlier material tests. Each contact relaxed on the time for material to travel from depth 0.5 to 2.0, measured here as 28.431 model units. In one trial its unloaded configuration traveled with the carrier from a straight inlet. In the other it renewed at a fixed site. Both use quasistatic brane shear, zero bare wall mass and zero occupied wave. The steady flow, order profile and axial geometry are maintained. [Renewal law](mouth_bend_renewal.md) · [Matched carrier preparation](mouth_two_depth_load.md)

**Motion.** For the traveling contacts, the rising oscillatory drive gave the following mesh-34 displacement amplitudes per unit force and signed phase lags. The angular frequency was 0.40 and the envelope growth rate was 0.30. Negative lag denotes a phase lead. These phases describe this drive and do not measure propagation time.

| Gauge | Middle push | Flare push |
|---|---:|---:|
| Entrance | 1.690, +131.1° | 0.285, −40.2° |
| Middle | 1.738, +101.2° | 0.294, −59.0° |
| Flare | 0.184, −99.4° | 0.531, +104.5° |

**Force during and after the pulse.** A separate smooth pulse peaked at time 4 and ended at 8. The table gives the actual integrated force exerted by the contacts on the wall between the waist and middle gauge, divided by peak applied force. Positive points along the push. The downstream wall takes the opposite resultant in this internal shear model.

| Pulse and time | Renew in place | Travel with material |
|---|---:|---:|
| Middle push at 4 | +0.17838 | +0.17727 |
| Middle push at 16 | +0.05348 | +0.02401 |
| Middle push at 28.375, near one transit | +0.07444 | +0.00456 |
| Flare push at 4 | −0.01823 | −0.01810 |
| Flare push at 16 | −0.02203 | −0.00929 |
| Flare push at 28.375, near one transit | −0.05144 | +0.00034 |

We replayed the measured positions through the carrier alone, subtracted the applied load from the required carrier load, and compared the inferred material-on-wall force with direct integration of attachment stress. Equivalent station forces agreed within $4.83\times10^{-9}$ of the peak push throughout the fine-mesh pulse. The table uses the spatial traction distribution from the material equations. The three positions alone give forces conjugate to the three measured motions. The order-weighted stress trace ended at depth 1.924 on the finer mesh, above the flare gauge. A carrier calibration from the coarser mesh could not resolve the smallest flare-station signal. [Force subtraction](../tools/mouth_two_depth_renewal.py) · [Spatial traction and work](../tools/mouth_two_depth_transfer.py)

**Gravity implication.** Each law transmitted load yet retained a growing unforced bend, with an amplitude e-folding time of 17.77 or 17.89 model units for the identified mode. The pulse budget includes signed work from prescribed reference transport. Its physical energy supply, contact lifetime, wall grip and shear inertia still need independent selection. Measuring traction after the push ends would distinguish the two renewal laws. A physical free-fall coefficient remains open. Two supplemental fixed-site spectra rejected an additional near-zero candidate. An independent stress formulation verified the positive growth rate and pulse forces on both meshes. [Fixed-site verification](../tools/mouth_two_depth_local_stress.py) · [Cross-mesh comparison](../tools/mouth_two_depth_renewal_compare.py)

Working base `49dca23e`. Two independent local reviews are recorded. Formal review and commit are pending.
