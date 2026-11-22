# How much shear does the receiving material retain over a mouth transit?

**The dense receiving material relaxed nearly all the sideways stress that reversed with the applied strain. Both crystalline controls retained substantial stress.** All 38 trajectories passed the prespecified late-window structural checks at bath temperature 1.088. The hold lasted 840.17332 native units, one nominal transit if the site mass, pitch and bond energy from the earlier impact calibration apply to this material. [Short preparation controls](mouth_independent_receiver.md) · [Clock calibration](mouth_pressure_matched_receiver.md)

| Preparation | Undriven late density | Original neighbors absent | RMS displacement¹ | Final density-wave peak | Retained shear, ±0.02 | Retained shear, ±0.01 |
|---|---:|---:|---:|---:|---:|---:|
| Sample 87, receiver A | 1.29618 | 89.54% | 18.58 | 0.09163 | −0.61% | −0.14% |
| Sample 87, receiver B | 1.28940 | 88.93% | 18.26 | 0.09258 | +0.70% | +1.60% |
| Sample 87, receiver C | 1.30496 | 90.54% | 19.21 | 0.08368 | −0.62% | −1.55% |
| Independent sample 41, receiver | 1.30252 | 87.78% | 17.90 | 0.10002 | −0.34% | −2.84% |
| Sample 87, crystal A | 1.67431 | 70.96% | 1.29 | 0.74546 | +51.07% | +54.82% |
| Sample 87, crystal C | 1.67533 | 68.22% | 1.30 | 0.85148 | +41.68% | +33.82% |

¹ Undriven material displacement in bond lengths, corrected for box dilation and center motion. All 128 sites remained connected at the sampled times. One early block of receiver A under tensile pressure averaged density 1.19783, below the late-window threshold of 1.2. Its eight late blocks met the criterion. [State criteria](../tools/mouth_receiver_transit_gate.py)

**Signed shear.** For each amplitude we subtracted the negative-strain stress from the positive-strain stress, then divided the last-quarter difference by the prompt difference. For the three histories from sample 87, the ratio of the mean differences was **−0.10% [−3.27%, +1.44%]** at full strain and **+0.06% [−5.87%, +3.70%]** at half strain. Both approximate 95% Fieller intervals are conditional on one founding geometry and fit inside the prespecified ±25% margin. Sample 41 was a single independent prediction. The half-strength signal common to both strain signs, measured against the undriven baseline and normalized to the prompt response, ranged from −13.52% to −3.23% for sample 87 and reached +12.96% for sample 41. A small signed difference leaves room for that offset.

**Pressure.** With external pressure steps of ±0.25, sample 87 had positive volume compliance **0.09117 [0.04693, 0.13542]** and a realized mean-pressure difference **0.49223 [0.37015, 0.61432]** (approximate 95% intervals). Sample 41 gave 0.07970 and 0.42114. These are finite-time measurements including a tensile branch. The independent audit reconstructed 608 block readings and 652 stress evaluations across six undriven and 32 driven trajectories. [Measurements and independent reviews](../../../research/pde_ledger_v3/_measurements/mouth_receiver_transit.json)

**Interface preparation.** With zero imposed mean pressure, the normal pressure along the laboratory fourth direction ranged from −0.05175 to +0.31546 in the receivers, versus +1.02328 and +1.23619 in the crystals. The crystal's undriven sideways stress in the tested direction was about −1.96, versus −0.015 to +0.013 in the receivers. On the *actual periodic face* chosen to join each cell, the magnitude of the total tangential traction was 4.94 and 4.93 for the crystals, versus 0.01–0.06 for the receivers. Units are bond energy per bond length to the fourth power. The common geometry therefore needs its own zero-load measurement. An instantaneous fit to the crystal's in-plane periods added 1.75 and 9.80 bond energies per receiver site on the two selected starting geometries. In the subsequent separate-material test, 11 of 12 histories met the criteria after the fit. Neither pair resolved compatible traction within the required margin. A force during conversion and the physical mouth acceleration remain open. [Common-plane preparation and traction measurements](mouth_receiver_common_face.md) · [Free-fall roadmap](../../../docs/zero_pn_free_fall_roadmap.md)

Working baseline `a91a633f`. Formal repository review pending.
