# What tension does the ordered-to-superfluid interface transmit?

**One tension measured on a neutral patch predicted the normal holding force when the curvature was halved and when the patch was shaped like a straight throat section.** This determines the response of the model’s shear-order interface. Applying it to a draining flare also requires identifying the material surface that carries the wall load and the normal stress on its outer face.

### Material and preparation

Both sides contain the same coherent material at rest, at equal density and phase. The brane side supports transverse shear and the disordered bulk side loses that stiffness. No shear displacement or trapped wave was excited. We used the written order energy $W(\chi)=a_B\chi^2(1-\chi)^2$ and gradient energy $\kappa_B|\nabla\chi|^2/2$, with the unspecified mixing and throat terms inactive in this neutral preparation. The spherical calculation relaxed the order field and used the exact uniform coherent solution. The straight-collar control called the native Stage 1 coherent carrier solver on the same elements as the order field. Neither test used the retired classical particle material. [Material action](../../../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md) · [Spherical solver](../tools/mouth_curved_order_patch.py) · [Native carrier control](../tools/mouth_curved_collar_patch.py)

We fixed the order at its two distant pure-phase values, held the half-order surface at the stated curvature, and let the profile relax on either side. The signed pressure is the external force per unit three-dimensional interface area needed to hold it there, positive along the normal from ordered material to disordered material. The uniform carrier supplies equal far pressures. Denote the interface width by $\delta_B=\sqrt{\kappa_B/(2a_B)}$ and use the sum of the three principal curvatures for $K$.

| Patch and $K\delta_B$ | Prediction from spherical calibration, $p_{\rm hold}/a_B$ | Measured $p_{\rm hold}/a_B$ |
|---|---:|---:|
| Flat, 0 | 0 | 0 |
| Spherical, +0.06 | Calibration | +0.02000032021 |
| Spherical, +0.03 | +0.01000016010 | +0.01000003991 |
| Straight collar, +0.03 | +0.01000016010 | +0.01000000110 |

The negative-curvature runs reversed the force. The two full-curvature spherical readings give $\sigma_{\rm fit}/(a_B\delta_B)=0.33333867$. The independent flat-interface energy and integrated stress anisotropy each give $1/3$ to the reported precision. The calibrated slope predicted the spherical half-curvature force within 0.00120%. Recomputing the energy after moving the spherical pin gave the same force independently. The collar bends in two directions and is straight in the third. Its finite-element energy check also reproduces the holding force. The small geometry-dependent corrections arise within the square-gradient law itself. These data do not require a separate bending modulus at the two sampled curvatures. The energy scale $a_B$ and length $\delta_B$ still need independent calibration in mouth units. [Spherical measurements](../../../research/pde_ledger_v3/_measurements/mouth_curved_order_patch.json) · [Collar measurements and independent sphere calculation](../../../research/pde_ledger_v3/_measurements/mouth_curved_collar_patch.json)

### Backing pressure in the draining mouth

At equal density $n$, both pure phases of the neutral law have homogeneous pressure $nU'(n)-U(n)$. The equal-depth order well supplies no additional difference. For a draining wall, use the normal momentum flux at the outer face under a specified coherent reservoir state and an interface law for any crossing material. The flow through the conduit already enters the inner force measurement. The holding pressure in the table measures interfacial work and enters the balance once.

All eight radial-force bands in the present mouth mesh touch fluid on the conduit side only. The receiver chamber begins downstream of them. The next test is a two-sided wall patch at the same supply. Locate the half-order surface, apply a known weak normal load on the receiving side, and read the wall reaction and the normal material current together. This determines whether the calibrated layer transmits that load and tests the reflecting boundary assumed by the mouth solver. A subsequent free-shape calculation needs one exterior geometry and pressure law carried across depths and radii. The affine pressure used to exclude the five trial shapes was only a diagnostic. [Geometry and pressure audit](../../../research/pde_ledger_v3/_measurements/mouth_curved_patch_backing.json) · [0PN roadmap](../../../docs/zero_pn_free_fall_roadmap.md)

**Independent local reviews.** The interface reviewer reported “Scoped clear for the reported neutral-interface calibration, native straight-collar transfer, and draining-mouth backing gap.” The geometry reviewer reported “Scoped clear for the final report’s identification of the neutral constitutive law, coherent solver, equal homogeneous pressures, sign of the applied holding traction, and lack of lateral backing data in all eight mouth bands.”
