# What a tagged attachment predicts after the 0PN middle push

**Fresh material can carry a residual sideways force after the originally loaded material has left the waist-to-middle section.** Which force the mouth transmits still depends on how attachments remember their unloaded position. We tested the same radius-one cubic mouth at inward intake 0.25. [Guarded readings and source hashes](../../../research/pde_ledger_v3/_measurements/mouth_middle_tag_grip.json)

**Following the same pair.** At middle depth 1.175 we initialized the earlier neutral 24-site patch from its saved relaxed coordinates. Its nine central-bond pairs included the preselected sites 5 and 12. We carried the patch with the computed flow, holding its lower layer to the local station by the existing elastic connections. The cross-layer well strength followed the local shear order divided by its starting value, 0.87530. We displaced the two tagged sites slightly in opposite directions, released them, and subtracted the reversed tap with the same thermal noise. Separate frozen copies tested their stiffness to a fresh displacement. Both read the direct sideways force of site 12 on site 5. The tap also strained site 5’s support. [Tagged-pair calculation](../tools/mouth_middle_tag_grip.py) · [Neutral material and support](neutral_patch_conversion.md)

| Fine-grid observation across four thermal histories | At half order | At 1% order |
|---|---:|---:|
| Time and downstream depth from the middle start | 1.592 at depth 1.2520 | 4.720 at depth 1.4275 |
| Original pairs continuously registered at sampled steps, each history | 9 of 9 | 9 of 9 |
| Tagged-pair stiffness under a fresh displacement, relative to entry | 57.09–57.25% | 1.14–1.15% |
| Signed force remaining from the released tap, relative to its initial signal | +0.0836–0.0854% | about +0.0003% |

The released signal first stayed below half for at least 0.05 model time units starting at about 0.003. It subsequently changed sign. Pair lifetimes are censored at the final reading. We observed rapid unloading of a connected pair. The independent speed of its unloaded attachment to the wall is still unknown.

**Force during replacement.** A separate forward solve applied the original middle push, peaking at time 4 and ending at 8. The carrier took 10.40093 time units from the waist to the middle. By time 18.40093 it had replaced all the material in that wall segment. The table gives the integrated attachment force on that segment in percent of the peak push, positive along the push. Both trial laws assume unit ordered grip and a relaxation time of one waist-to-depth-2 transit, 28.43103. [Force calculation](../tools/mouth_middle_replacement_force.py) · [Original renewal laws](mouth_two_depth_renewal.md)

| Time | Rest position renews in place | Rest position travels with carrier |
|---|---:|---:|
| 16, during replacement | +5.35% | +2.40% |
| 18.5, just after full replacement | +7.62% | +4.04% |
| 28.375 | +7.44% | +0.46% |

At time 16, integration on the two sides of the incoming material front gave **+16.00%** on the fresh side and **−13.60%** on the older side for the traveling law. After full passage, fresh material occupied the entire segment. Elastic stress and motion can pass a load between cohorts. For in-place renewal, the front tracks only the carrier. Both original whole-mouth preparations admit a growing bend. The reported forces are internal loads at finite times.

**What would select a prediction.** Track the wall-side and brane-side ends in depth relative to an independently marked wall station. On matched copies at several depths, hold a weak relative sideways displacement and measure the total force at the wall as fresh material passes. Reverse and halve the displacement and include a zero-drive control. The transported unloaded position and the force per unit area over time would calibrate the wall law independently of the middle push. The tagged-pair clock and lifetime alone cannot fix that law.

Working base `49dca23e`. Independent local review is recorded in the measurements. Formal review and commit are pending.
