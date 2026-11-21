# Calibrating a held 0PN mouth with a shear pulse

A weak sideways pulse distinguished four hidden conditions behind the wall-facing sites. We held one radius-one draining mouth at intake 0.25 with the same rim grip. Opposite pulses and a half-strength pair drove the brane-facing edge for 2.3835 model time units. Two internal gauges measured lateral velocity and shear force per area, half a patch length and one patch length behind the sites. The lateral fluid increment started at zero. Each drive had a matched unforced control and four thermal histories. [Pulse experiment](../tools/neutral_patch_shear_pulse.py) · [Measured traces](figures/neutral_patch_shear_pulse.png)

| Hidden backing | Resistance inferred from the pulse, in patch units |
|---|---:|
| Stress-free | 0 to 0.012 |
| Resistance 0.5 | 0.475 to 0.526 |
| Resistance 4 | 3.72 to 4.31 |
| Held at rest | At least 78 |

These bounds allow a 1% RMS discrepancy in the early traction reading, with known layer thickness, density, viscosity and input velocity. They describe resolution in this numerical experiment. Using the nearer gauge’s velocity as input, we fitted its traction in the first 65% of the transit. The deeper gauge and remaining time were held out. Velocity peaked at the deeper gauge 1.62–2.62 model time units after the drive peaked. Accounting for stored momentum predicted late traction within 0.041% relative RMS for the four finite-layer preparations. An instantaneous force proportional to speed missed it by 36–609%. [Blind inversion](../tools/neutral_patch_pulse_inverse.py)

**A different material behind the same gauges.** We extended the initially resting viscous medium beyond the nominal backing. With the far boundary 10.5 patch lengths below the sites, a homogeneous half-space with the supplied density and viscosity predicted late traction within 0.24% relative RMS without fitting a resistance. A one-patch-length layer with its best early-fit terminal resistance missed late traction by 98.94%. Moving the far boundary between 6.5 and 10.5 lengths changed the recorded histories by less than 0.003%. This is a measurable signature of momentum diffusing into depth. [Independent checks and guarded measurements](../../../research/pde_ledger_v3/_measurements/neutral_patch_shear_pulse.json)

**Connection to free fall.** Four additional calibrated drive shapes determined the response during conversion. A separate release measured the contribution from pre-existing lateral flow. Combining these responses with the calculated companion velocity predicted the signs of both wall-layer drag and held-rim force in all 32 comparisons at 8 and 16 mouth radii. Eight cases were run for the first time after the predictions were frozen. Maximum relative errors were 0.07070% for drag and 0.000877% for rim force. This checks transfer between drives within the same material laws and four thermal histories. [Forecast](../tools/neutral_patch_pulse_forecast.py)

The physical receiving material still needs an independent transverse stress law as order fades. That would determine how much momentum it stores or returns. Predicting the acceleration of a free mouth also needs whole-mouth force and inertia. Working base `49dca23e`. Independently reviewed locally; formal review and commit pending.
