# Does the formation wave straighten the first bend?

**The ringing wave pushed the held flare in the direction of the imposed bend when the material’s shear modulus was independent of bend strain.** We used the same radius-one, intake-0.25 mouth that rang after a small depthward stroke of the collar. We launched the wave while the mouth was straight, then shifted the flare sideways by 0.01 mouth radii over times 12–20. The offset grew smoothly with depth across the first bend, from \(w=0.55\) to \(1.15\). The comparison applied the same bend to a quiet mouth. [Formation and ringing](mouth_formation_wall_work.md) · [Bend after formation](../tools/mouth_ringing_bend_after.py)

| Wave-induced force, Hann-weighted mean over times 24–48 | Entrance | Flare |
|---|---:|---:|
| Quiet, offset +0.01 | 0 | 0 |
| Ringing, offset +0.01, finer mesh | \(-7.9513\times10^{-11}\) | \(+7.9513\times10^{-11}\) |
| Ringing, offset +0.01, reference mesh | \(-7.9362\times10^{-11}\) | \(+7.9362\times10^{-11}\) |
| Ringing, reversed offset, reference mesh | \(+7.9362\times10^{-11}\) | \(-7.9362\times10^{-11}\) |

The signs use the same sideways axis at both ends. A positive flare force at positive offset tends to increase the bend. Holding the flare against that contribution takes an opposite force. The entries are generalized wave forces on the translated sections and their attached material. Equal and opposite follows because the wave action depends on their relative displacement. Separate physical traction readings at the walls would also require the carrier, material attachments and reservoir momentum.

**Checks on the sign.** We let the previously round wave acquire its first angular distortion as the sections moved. The force includes the change in momentum of the moving wave displacement. Short averages can have either sign during ringing. A smooth average over many oscillations gave the outward force above. Halving the time step, moving the receiving cut deeper and shortening the bend ramp preserved its sign. An independent triangular calculation gave \(+7.926\times10^{-11}\) on the flare at the same offset to first order. Order and intake were fixed and the wave followed the attached material. [Four-dimensional wave action](../tools/mouth_ringing_bend.py) · [Guarded comparisons and independent reviews](../../../research/pde_ledger_v3/_measurements/mouth_ringing_bend.json)

**A material measurement could change the outcome.** As a separate hypothesis, let the shear modulus increase or decrease under slow local section shear \(s\), with \(\mu(s)/\mu(0)=1+\lambda s^2+O(s^4)\) and fixed transverse inertia. Every \(\lambda\) gives the same wave in the straight mouth. On the finer mesh the measured leading-order force at offset 0.01 is \(F_{\rm flare}\simeq 7.9513\times10^{-11}-2.4870\times10^{-9}\lambda\) in model force units. The wave contribution changes to a restoring sign above **\(\lambda\simeq0.0320\)** for this preparation and averaging window. The threshold applies to this one bend and wave. The total stiffness, inertia and free-fall acceleration remain open. [Material sensitivity calculation](../tools/mouth_ringing_bend_modulus.py)

The next independent test is to hold a neutral material sample at zero and opposite weak shears, at the order levels present in the collar. Measuring its incremental shear modulus with the same small probe at each strain, then repeating at half strain, would determine the sign and size of the quadratic response without using the mouth force to set it.

Working base **1dc45d09**. Three local reviewers checked the force, time evolution and material condition. Formal review and commit pending.
