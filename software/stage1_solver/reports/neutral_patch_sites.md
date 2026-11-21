# Can a neutral patch form anchors and predict their holding force?

**A neutral substrate with short-range binding formed persistent interlayer pairs. Its independent tensile calibration predicted the force to hold a surviving pair one eighth of a pitch off center.** Two material laws with the same normal pull curve gave substantially different sideways forces. This calculation supplies a candidate microscopic mechanism for a 0PN attachment. The physical mouth’s force law remains open. [Previous native anchor search](mouth_anchor_order_spectrum.md)

**What was assumed.** Two parallel planes, one model length apart, each carried 12 neutral mobile constituents with free edges. Both planes started from independent random in-plane positions. All pairs interacted, with no assigned neighbors. Same-plane binding used a Morse well of depth 1, equilibrium length 1 and inverse range 7. Cross-plane binding used depth 0.6, equilibrium length 0.94 and inverse range 10. For tangential separation \(q\) and gap \(d\), the central law used the distance \(R=\sqrt{d^2+q^2}\). The layered alternative used \(R=\sqrt{d^2+q^4/a^2}\), with \(a=1\) assumed equal to the same-plane equilibrium length. Both reduce to the same normal interaction at alignment. A real calibration would measure normal gap, within-plane and cross-plane stretch curves and an off-axis response for the angular dependence. [Material and formation calculation](../tools/neutral_patch_sites.py)

**What persisted.** After energy relaxation we chose the eligible interlayer pair closest to the patch center. Eligibility required mutual nearest neighbors across the planes and binding deeper than 60% of the cross-plane well. We followed the same identities at all 4,000 Brownian steps over two model time units at \(k_BT=0.008\), then relaxed the patch again. The following counts give the initial strict pairs and the subset that remained eligible throughout, including the final relaxation.

| Preparation | Central law | Layered law |
|---|---:|---:|
| Calibration patch (191) | 11 of 11 persisted | 0 of 5 persisted |
| First test patch (817) | 2 of 3 persisted | 0 eligible initially |
| Second test patch (1321) | 9 of 9 persisted | 9 of 11 persisted |

In patch 1321, 8 of the 9 central pairs and all 11 layered pairs eligible after relaxation were new matches relative to the random initial positions. Both preselected pairs survived the strict test, including an independent check at half the integration step. In patch 817 the central tag survived at the 50% threshold and crossed the 60% threshold during observation. No layered tag qualified there. The force test below uses patch 1321. These observations establish persistence over the sampled interval in a finite patch.

**Prediction and held-out reading.** The reference pitch was the median same-plane nearest-neighbor distance on patch 191. Separately, an isolated neutral pair supplied its normal tension curve \(T(R)\) through energy changes under normal extension. The largest sampled tension was 2.99986007 model force units. We saved \(H=T(R)\,\partial R/\partial q\) before applying the clamp to either test patch. At fixed midpoint and normal gap, the clamp imposed the actual relative position \((q,0)=(p/8,0)\). Other particles relaxed. We measured the total energy change under small opposite translations of the two tagged sites. Repeating that reading at identical coordinates with just their direct interaction removed gave the surrounding-material contribution. [Prediction and clamp](../tools/neutral_patch_hold.py)

| Interlayer law | Reference pitch | Predicted pair force | Measured pair force | Total clamp force, including surroundings |
|---|---:|---:|---:|---:|
| Central | 0.997524 | +0.371114 | +0.371104 | +0.884152 |
| Layered | 0.993908 | +0.01139775 | +0.01139778 | −1.389453 |

The pair forces were 12.37% and 0.380% of the sampled tensile force scale. At one sixteenth pitch, the measured forces fell to 0.185580 and 0.00142454, respectively. In the layered example the surroundings reversed the sign of the total clamp reading. The subtraction isolates one direct interaction in an additive model. [Readings, guarded outcomes and review scopes](../../../research/pde_ledger_v3/_measurements/neutral_patch_sites.json)

**What the next material test must decide.** At fixed normal preload, a weak sideways displacement gives a linear pair reaction for the central law and a cubic one for this layered law. An independent sideways measurement can distinguish them even when the normal pull curves agree. Transferring the result to the mouth requires evidence that the same sites and their force scale persist as flowing brane material loses shear order. This neutral test includes no intake, conversion or free-mouth motion. The earlier narrow and wide crossed-strand fixtures have a different microscopic construction and force normalization. [Earlier mouth clamp](mouth_two_site_clamp.md)

Working base `49dca23e`. Formal review and commit are pending.
