# Can converting material power a wave at the mouth?

**A converting patch can do positive work on an oscillation if its attachments change their relaxed lengths strongly enough.** Ordinary loss of shear stiffness absorbed energy at the tested frequency. The physical attachment response and energy supply remain unmeasured.

Material crossed the same planar patch at intake 0.25 as its shear order fell from 0.85 to 0.45. We prescribed a small, spatially uniform sinusoidal slip at the radius-one mouth’s scalar-mode frequency, 2.14006206. The order field responded throughout the cycle. Fresh contacts entered without stored extension. We measured the contact force times the slip velocity, heat from conversion and contact renewal, and energy crossing both ends. This is a local work measurement. Applying it to the trapped mode requires its actual slip at the converting material. [Patch calculation](../tools/mouth_conversion_wave_patch.py) · [Earlier contact law](mouth_conversion_contact_transfer.md)

| Material response | Work on probe per cycle, amplitude 0.0001 | At 0.0002 |
|---|---:|---:|
| Ordinary order-dependent shear stiffness | −6.38 × 10⁻²² | −1.02 × 10⁻²⁰ |
| Contact with conversion changing its relaxed length, β = 1 | −1.46 × 10⁻⁸ | −5.84 × 10⁻⁸ |
| Stronger change, β = 10 | +4.41 × 10⁻⁸ | +1.77 × 10⁻⁷ |

**Why the sign changes.** A contact carries its earlier strain downstream. Under the strong trial law, conversion amplifies that strain enough for the contact to pull in the direction the probe is now moving. Halving the frequency made the same material absorb energy at both amplitudes. An independent calculation found a weak-probe sign crossing at β = 7.78380 for this transit and frequency. This predicts a measurable material property. Hold the ends of a tagged, slightly prestrained contact at fixed slip. At the crossing, its elastic extension relative to its evolving force-free position grows by a factor of **51.99** while it passes through the patch. The β = 10 test corresponds to a factor of 162.43. [Two independent checks, controls and raw measurements](../../../research/pde_ledger_v3/_measurements/mouth_conversion_wave_patch.json)

**Who pays.** The bare order potential has equal energies at its homogeneous endpoints. For the table we gave the incoming material an internal energy supply, of which **0.0258199 per entering constituent** was released between these two sections. This is a specified material preparation. At amplitude 0.0002 and β = 10, the probe received **2.40565 × 10⁻⁷ per passing constituent**. Most supplied energy became conversion heat. Most energy transferred into the contacts left as stored strain. Net mechanical work also left through the patch boundaries. Removing the supplied energy still produced positive probe work, but then mechanical work entered through the boundaries and paid for it.

**0PN decision.** Both weak probes gained work in the strong-contact example. Increasing the amplitude changed the order-boundary work from outward at 0.00022 to inward at 0.00023 while probe work stayed positive at both readings. There is no measured saturation or selected throat-wave energy. The next independent test is to tag a contact at fixed slip with opposite small initial prestrains, locate its zero-force position during conversion, and measure the inlet-to-outlet extension ratio and available energy per constituent. Those measurements would test the proposed source before using it in a free-fall calculation. [Previous whole-mouth energy test](mouth_wave_drain_power.md)

Working base `1dc45d09`. Local review found “Scoped clear for work on the prescribed oscillatory probe and signed cycle-energy accounting on this flowing planar patch.” Formal review and commit pending.
