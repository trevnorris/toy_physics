# Where does the passing brane acquire its sideways motion?

**In this candidate material, the passing sites acquired almost all their sideways speed before reaching the deep flare.** We tracked them from depth 1.175 to 1.65 on the radius-one cubic mouth at intake 0.25. Its meridional arc was represented by 56 material rows, free to slide sideways. We prescribed their axial locations, transverse spacing, normal gap and surviving within-layer bonds. One transverse pitch is 0.00997514 mouth radii. [Material and carrier inputs](../tools/mouth_mobile_lining.py)

Relaxing small offsets restored the prepared staggered pattern. With fresh sites entering uniformly across a pitch at zero sideways speed, every site has a reflected partner. Their forces cancel on each free row even as individual sites scatter. The imposed axial flow can supply energy to that scattering. [Symmetry calculation](../tools/mouth_lining_symmetric_passage.py)

| Flux-weighted reading on the finest symmetric run | Result |
|---|---:|
| Sites within one eighth of a transverse pitch of a site on the reference row at depth 1.55 | 24.18% (uniform reference 25%) |
| RMS sideways speed arriving there | 3.87405 pitches per model time |
| RMS change of individual sideways velocities from there to 1.65 | 0.002764 pitches per model time |
| RMS outlet position difference when sites coast freely after 1.55 | 0.004939 pitch |

Across four calculations the eighth-pitch fraction ranged from 23.44% to 24.48% at the gauge. It reached 27.27% one hundredth of a depth unit farther on in the finest run. The distribution was broad and varied with depth.

**Bound with moving rows.** At the tested fine-mesh gauge, the tabulated contact law limits the change in any site's absolute sideways velocity through the last stretch to 0.034900 pitches per model time. Bounding the underlying continuous Morse forces gives a conservative ceiling of 0.160 over the recorded transit. Both allow arbitrary sideways row motions and require the specified axial flow, gap, interaction range and weakening of contact with shear order. [Force bound](../tools/mouth_lining_impulse_bound.py)

**Roadmap decision.** Shifting the entire fine strip and feed together changed the calculated wall momentum from +3.07543 to −3.41958 in units of one row mass times one pitch per model time. Physics requires the same reading. The momentum ledger balanced, but the signed recoil failed this check. Step 2 still needs a measured force transmitted to the curved mouth. Under opposite weak displacements of one held mouth, measure traction where the brane is still ordered alongside the front-to-back difference in attached material, its exchange rate and its flow between two depths. These readings distinguish sites that carry the bend from passing sites that reproduce its shape. [Unrestricted crossing](../tools/mouth_mobile_lining.cpp) · [Surface counting relation](mouth_orientation_and_motion.md)

Working checkpoint `a91a633f`. Independent local scientific review completed within the recorded scopes. Formal repository review is pending. [Guarded results and literal verdicts](../../../research/pde_ledger_v3/_measurements/mouth_mobile_lining.json)
