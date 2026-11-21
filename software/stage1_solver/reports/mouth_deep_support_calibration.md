# Can the draining flow determine the deeper wall support?

**Phase 1 of the 0PN free-fall roadmap stopped at independent calibration.** At radius 1 and intake 0.25 we tested a wave-empty cubic mouth with an isotropic membrane of constant tension and uniform pressure on its other side. These two unknown constants could in principle be inferred from the stationary shape before observing any sideways motion. The model failed to predict the carrier pressure in two unused portions of the same deep bend. [Calculation](../tools/mouth_deep_support_calibration.py) · [Field, guard and force records](../../../research/pde_ledger_v3/_measurements/mouth_deep_support_calibration.json)

For a compact radial displacement of each spherical section, let δA and δV be the changes in three-dimensional wall area and four-dimensional enclosed volume. Steady local balance requires the carrier’s outward work per δV to equal the backing pressure plus the tension times δA/δV. The pressure zero is relative to resting carrier. We computed the carrier work from both the reflecting wall’s stress and the interior variation of the stationary fixed-intake potential. A direct area difference checked δA. We fitted the outer two of four disjoint bands, all inside the deeper bend, and predicted the inner two.

| Depth band in mouth radii | Wall pressure from carrier | Predicted from the two fitted bands |
|---|---:|---:|
| 1.22–1.40, calibration | −0.00157264 | −0.00157264 |
| 1.42–1.59, prediction | −0.00153316 | −0.00111260 |
| 1.62–1.79, prediction | −0.00110886 | −0.00086475 |
| 1.81–1.99, calibration | −0.00072697 | −0.00072697 |

On mesh 48 the two unused bands required additional signed pressure of −0.00042056 and −0.00024411 relative to that fit. The two force calculations gave the same signs and magnitudes on meshes 26, 34 and 48. The fitted tension of 0.00387236 therefore has no independent standing as the total lateral support. An imposed shape can be maintained by additional stress, a depth-dependent backing, or an external holder. These measurements alone cannot separate them. Even a successful radial fit could leave a sideways bending energy proportional to the squared offset undetected on the straight preparation.

**How to distinguish the two earlier support models.** We calculated two weak sinusoidal loads at angular frequency 0.005, pushing neighboring stations sideways with equal and opposite forces. For the flare–middle pair, a force amplitude of 10⁻⁹ on each station gave an in-phase flare displacement relative to the middle of **6.80694 × 10⁻⁹** for deep tension 0.02 and **3.34853 × 10⁻⁹** for 0.04 on mesh 34. Mesh 26 gave 6.80691 × 10⁻⁹ and 3.34852 × 10⁻⁹. The first-bend tension stayed at 0.02 throughout. The two-mesh forecasts include the complex displacement of the entrance, middle and flare for both loads. Recovering an assumed tension from its own predicted motion checks the readout only. These two loads probe internal bends. If the attachment also carries momentum in common translation, an independent common-motion load is needed. [Earlier local strip and single-mouth release](mouth_entrance_tag_release.md)

**Decision.** To advance phase 2, supply the total depth-dependent material reaction independently. A candidate next test is to take the same material used in the attachment calculation, extend it across the deeper flare, calibrate its force per represented wall area, and check the flare–middle displacement against a separate applied load. Coverage, force units and additional wall stresses must be accounted for before treating that as total support. Until then the physical single-mouth release, common response across radii and unequal-drain prediction remain open. Baseline checkpoint `1dc45d09`.

Independent local review: “Scoped verdict: the phase‑1 radial falsification and conditional balanced-push forecasts are supported. Stop before phase 2.”
