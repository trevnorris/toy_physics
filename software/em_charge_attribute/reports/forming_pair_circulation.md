# Can a forming pair acquire opposite wave circulation?

**Yes, as a mechanical response in the scalar-wave model.** A shared standing wave with zero initial momentum can develop opposite orbital angular momenta above and below the brane midplane. Reversing a small deformation reverses which side turns which way. The symmetric opening considered here gives both assignments equal status.

**Formation test.** Use the existing neutral support equation \(u_{tt}=\nabla^2 u-Vu\), in dimensionless units and per unit area in the two spectator directions. Resolve one brane coordinate \(x\) and the depth \(w\). Choose \(V=x^2+(w^2-d^2)^2/4\). Its single minimum splits into minima at \(w=\pm d\), with \(d=2(3\tau^2-2\tau^3)\), \(\tau=t/4\), on \(0\le t\le4\). A control keeps \(d=0\). This is prescribed confinement. Material drainage, wall backreaction, particle production and charge selection are outside this calculation. [Parent wave Hamiltonian](../../../notes/inner_throat/inner_throat_hard_mode.md#L459) · [Mixed-plane circulation](../../../docs/native_light_em_and_vortex_throat_interpretation.md#L1036) · [Orbital angular momentum of light](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.45.8185)

Let \(G=\exp[-(x^2+w^2)/2]\). Start with \(u=xG\) and \(\dot u=\epsilon(w^2-1/2)G\), for \(\epsilon=0,\pm0.01,+0.005\). In each half measure the angular momentum of wave momentum density \(-\dot u\nabla u\) about that half’s center of energy, in common laboratory axes. The table gives upper-half angular momentum divided by its initial value for \(\epsilon=+0.01\). The wave is even in \(w\), which enforces opposite lower-half readings. [Two-dimensional solver](../forming_pair_circulation.py) · [Independent mode reduction](../forming_pair_circulation.wl)

| Confinement | Halfway, finite differences / spectral | End, finite differences / spectral |
|---|---:|---:|
| Splitting | 1.1334 / 1.1346 | 0.3107 / 0.3059 |
| Held fixed | 1.0811 / 1.0822 | −0.4568 / −0.4645 |

There is transient amplification even in the control. Reversing the seed reverses both readings, and zero seed gives zero circulation to numerical precision. At the end of the split, about 40% of the energy lies within one depth unit of the minima. These are half-space measurements of a shared wave. Long-lived circulation confined to two separate mouths remains unverified. [Guarded outputs and refinement checks](../../../docs/_measurements/forming_pair_circulation.json)

**Starting with zero angular momentum everywhere.** At \(t=0\), instead set \(\dot u=0\) and \(u=[x+\eta(w^2-1/2)]G\). At \(\eta=0.01\), the initial angular momentum rates about the two energy centroids are \(+0.005595272867\) and \(-0.005595272867\). About the common origin the upper-half torque is \(3\sqrt{\pi}\eta/8\). Stress across the midplane contributes \(+\sqrt{\pi}\eta/2\) and confinement contributes \(-\sqrt{\pi}\eta/8\). Both contributions reverse below. This initial torque also exists in the unsplit control. [Direct acceleration balance](../forming_pair_initial_torque.py) · [Independent stress integral](../forming_pair_initial_torque.wl)

**Implication for charge.** Reflection in \(x\), together with the irrelevant overall sign of the real wave, exchanges the two seeds while preserving the entire opening history. Neither assignment is selected. A physical formation law would have to produce and retain a circulation correlated with mouth orientation and give it a conserved electric response. The written neutral-wave Hamiltonian has zero direct response to a static electric probe. [Mouth matching and open physical requirement](mirror_internal_matching.md)

Working calculation on base e0a94594, uncommitted pending independent review.
