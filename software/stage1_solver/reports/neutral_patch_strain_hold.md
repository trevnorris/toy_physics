# How long a converting mouth patch bears sideways force

At inward intake 0.24999999999821823 we sampled a 24-site patch at four depths of the same draining mouth. At each depth we fixed shear order and normal gap, displaced the two layer centers oppositely by a total strain of 0.001 and held for 5.29667 carrier time units. Sites could rearrange. Both candidate bond laws have the same normal force. The layered law is less sensitive to small sideways offsets. Cross-layer binding scales with order in both. [Archived material preparation](neutral_patch_conversion.md) · [Experiment](../tools/neutral_patch_order_hold.py)

| Fixed shear order | Initial force / strain, central | Initial force / strain, layered | Force retained late, central | Force retained late, layered | Pulse error, central / layered |
|---:|---:|---:|---:|---:|---:|
| 0.95 | 18.0714 | 16.1112 | 30.2% | 73.0% | 1.36% / 5.27% |
| 0.50 | 9.5632 | 7.4144 | 53.1% | 79.6% | 1.14% / 7.85% |
| 0.10 | 1.8856 | 1.1967 | 86.5% | 80.0% | 1.98% / 14.49% |
| 0.01 | 0.1839 | 0.1074 | 91.0% | 73.2% | 5.74% / 22.35% |

Initial force per unit strain is in patch units. Retention is the last tenth of the hold averaged over eight thermal histories, relative to the initial mean. Only the 95%-ordered central material fell below half force for at least 5% of a transit. It reached that threshold at 2.9251 time units.

**Independent pulse.** The measured hold response $K$ predicts $F(t)=\int_0^t K(t-s)\,\dot\gamma(s)\,ds$ if it transfers between drives. We froze the prediction before applying a smooth strain pulse lasting 1.32417 time units. Eight new thermal histories used the same four starting frames, which coincide at 95% order. The table scores the internal force. The actuator also overcomes bath drag $6\dot\gamma$ during motion. Opposite and half-strength drives, zero-strain controls and a halved integration step checked the force reading. [Forecast and scoring](../tools/neutral_patch_strain_forecast.py)

Six of eight full-record errors met the predeclared 10% relative RMS criterion. Using only the initial stiffness did better in three cases. The integrated force after release was negative in all eight measurements and predictions. At half order with central bonds the predicted and measured impulses were −2.251 and −2.202 per unit peak strain. At 1% order they were −0.00934 and −0.00112. [Force histories](figures/neutral_patch_strain_hold.png) · [Recovery](../tools/neutral_patch_strain_recovery.py)

**Low-order replication.** For the two layered cases at 10% and 1% order, 64 new holds on the same four frames gave late retention of **80.1% and 56.2%**. On 64 independent pulses the original forecasts scored 12.14% and 20.98%. The enlarged hold forecasts scored **6.61% and 16.01%**, and the initial-stiffness-only forecasts scored 5.53% and 11.57%. The enlarged forecasts gave the correct sign of recovery impulse and overestimated its magnitude by factors of 2.05 and 9.62. The 1%-order error stayed above 15% when any one starting frame was omitted. [Replicated forecasts and sampling checks](../tools/neutral_patch_layered_hold_ensemble.py)

**Motion of the material patch.** In a noise-free control, we relaxed one preparation at each order. The maximum relative RMS error of its predicted odd pulse response was $6.8\times10^{-6}$. Holding the patch’s orientation fixed removed 97.1–98.0% of the central material’s force loss over a transit and 6.4–12.8% of the layered material’s. The extra grip also changes initial stiffness. This control concerns rotation of material in the finite patch while the mouth stays still. [Stationary control](../tools/neutral_patch_resting_order_hold.py) · [Angle grip](../tools/neutral_patch_resting_angle_hold.py)

**Open 0PN requirement.** Measure the same low-order patch after different waiting times before strain is applied, and measure how strongly the surrounding brane constrains its orientation. These would test whether one force-decay history can transfer to a patch driven by a companion. The current forces belong to the modeled patch. Its coupling to a free mouth remains open. [Measurements, guard validation and review](../../../research/pde_ledger_v3/_measurements/neutral_patch_strain_hold.json)

Working base `49dca23e`. The independent local review covers the fixed-order forces, frozen forecasts, recovery impulses and stationary orientation control. Formal review and commit are pending.
