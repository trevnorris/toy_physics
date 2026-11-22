# Can outgoing waves account for the resistance of a moving drain?

**The outgoing wave condition passes independent checks on a sphere. Applying it to the moving, draining mouth requires a matched exterior for the carrier and its steady flow.** In the preceding finite-box study, the corrected push and steady tow agreed to five significant figures. The two independent interior work readings were about 2.4% and 0.8% below the push in the largest box. [Finite-box audit](mouth_finite_wall_force_audit.md)

### What the outgoing test should measure

Using the solver’s $e^{+i\omega t}$ convention, a homogeneous acoustic dipole in $d$ spatial dimensions has radial potential $g(s)=s^{-(d-2)/2} H^{(2)}_{d/2}(ks)$. Its exact outgoing logarithmic derivative approaches $s g'/g=-(d-1)$ as $ks\to0$. At $\omega=0.000225$, $c=1$ and matching radius 18, the real parts are −1.999984 in 3D and −2.999992 in 4D. A plane-wave boundary misses this leading static gradient. [Calibration, inputs and literal outputs](../../../research/pde_ledger_v3/_measurements/mouth_outgoing_check.json) · [NIST Hankel limits](https://dlmf.nist.gov/10.7)

| Calibration in resting material | 3D radiated power per squared speed | 4D radiated power per squared speed |
|---|---:|---:|
| Rigid sphere driven at prescribed velocity | $\omega^4$ | $\omega^5$ |
| Point source or sink of fixed flux displaced at that velocity | $\omega^2$ | $\omega^3$ |

The difference follows from translating the original drain flow. That produces a field proportional to displacement $X=V/(i\omega)$. At fixed velocity the field has one extra inverse power of frequency, and the power has two. Direct acoustic energy flux at three observation radii, extending into the wave zone, reproduced both limits over six frequencies. An independent annulus calculation recovered the rigid-sphere resistance at matching radii 9, 13.5 and 18. These are calibrations of acoustic power. For an open drain, equating that power to total actuator work additionally requires the change in energy supplied by the two reservoirs. [Calculation](../tools/mouth_outgoing_dipole_calibration.py)

### What worked on the actual mouth

* **Moving boundary.** The archived solve used a chart moving with the mouth. The outgoing field in laboratory coordinates is $h_{\rm lab}=h_{\rm chart}-X\partial_x f_0$. Releasing the old boundary nodes also exposes a surface drive in the weak equation. Removing this artificial drive for an outer test function recovered direct volume advection to within $1.8\times10^{-11}$ in all three boxes. At zero intake the same correction removes a large drive in completely uniform material. This tests the equation used to excite the field. The corrected momentum definition used to measure force stays fixed. [Moving-boundary test and readings](../../../research/pde_ledger_v3/_measurements/mouth_translating_open_cut.json)
* **Stationary throughflow.** We matched the leading exterior phase of a 3D sink in a unit-thickness slab and a 4D source in a reflecting half-space, using the same prescribed intake of 0.25. The two successful preparations at box radii 9 and 18 required far phase differences of 0.4525058154 and 0.4525054766. An interior surface six radii from the outlet measured currents of 0.24995620 and 0.24995623. Their local flux varied with angle by 2.89% and 2.67% RMS relative to a uniform monopole. These stationary fields provide useful starting points for harmonic matching.
* **Middle box.** At radius 13.5 the phase adjustment did not reach the prescribed intake tolerance within ten updates. Its final intake discrepancy was not recorded. This supplies no stationary field or force at that radius. The earlier batch stopped at an interior residual tolerance stricter than the native carrier solver. [All outcomes and archived field hashes](../../../research/pde_ledger_v3/_measurements/mouth_outgoing_check.json)

### Where Step 1 stands

The receiving mesh ends on a cylindrical side and a flat roof. Selecting degree one in the three brane coordinates still allows variation along the fourth direction. On a spherical matching surface in the reflecting half-space this includes odd 4D angular degrees 1, 3, 5 and higher. The coherent field also has a propagating sound solution and a spatially decaying density solution. Applying one acoustic dipole coefficient separately to every face would discard this structure. The production outgoing push was held at this point.

The next implementation is a surface operator that matches both field components and their normal fluxes to the coherent exterior. Use the laboratory-frame correction and the two retained stationary preparations. Change the matching radius at fixed frequency, and compare actuator work, interior flux and both reservoir powers before lowering the frequency. For asymptotic 0PN free fall the roadmap requires the resistance divided by $M\omega$ to tend to zero on the relevant orbital scaling. A finite-band bound can instead support the separations for which it was measured. [0PN roadmap](../../../docs/zero_pn_free_fall_roadmap.md)

**Independent review.** The external-wave reviewer reported “Scoped clear for the outgoing calibration, moving-cut identity, two stationary preparations, failed middle preparation and stated limits in mouth_outgoing_check.md.” A second reviewer verified both drives and distinguished the three-radius annulus test, which used the rigid-sphere drive, from the flux measurements on both drives.
