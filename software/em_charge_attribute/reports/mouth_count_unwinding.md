# Can an open draining mouth lose its signed geometric count?

**Yes, for the count carried by the layer where the brane loses shear order.** A smooth deformation of a receiving-side flare takes that count from \(+1\) to zero with the neck open, the mouth still pointing into \(+w\), and the same positive number intake \(D\). The reflected flare goes from \(-1\) to zero. The exterior lip stays fixed. This is a finite-energy path through trial configurations. Whether a freely supported throat takes that path depends on its wave and rim dynamics. [Original definition of the layer count](stationary_mouth_charge.md) · [Wave support](trapped_wave_charge.md)

**Geometry.** Past a neck of radius \(R>0\), use the closed material sections \(r=|\mathbf x|\) of an expanding flare. In the conversion region \(2R\le r\le3R\), take its depth to be

\[
 w=H_*+\epsilon\left(\frac{r^2}{2R}-\beta z\right),\qquad z=x_3,\qquad 0\le\beta\le4,\quad\epsilon>0 .
\]

Choose \(H_*\) deep enough to keep the whole patch on the same side. Taper the tilt smoothly to zero in the fully ordered and fully disordered regions before joining the fixed neck and receiver. Continue the upper-face normal through the waist. Its in-brane projection points along \(\mathbf x/R-\beta\hat z\). On each section these vectors form a sphere of radius \(r/R\) centered at \(-\beta\hat z\). Their directions cover the unit sphere once while the origin is inside, and have degree zero once it is outside. At equality the projection vanishes at one point while the full surface normal remains regular. [Degree by signed preimages](https://www.ime.usp.br/~gorodski/teaching/mat5799-2015/gorodski-smooth-manifolds-2013.pdf)

**Count through the transition.** Weight the section degrees by the fraction of order lost, \(-d\chi=30u^2(1-u)^2du\), with \(u=r/R-2\). SymPy integrated the signed solid angle. Mathematica counted the preimages of a fixed direction. Their upper-mouth results agree. [Scripts and guarded outputs](../../../docs/_measurements/mouth_count_unwinding.json)

| Tilt \(\beta\) | 0 | \(9/4\) | \(5/2\) | \(11/4\) | 4 |
|---|---:|---:|---:|---:|---:|
| Count in the conversion layer | 1 | \(459/512\) | \(1/2\) | \(53/512\) | 0 |

The midpoint is a weighted sum across a layer of finite width. Each regular section still has integer degree. During the passage, the section through the zero has no defined direction at that point. Omitting that single section defines the displayed integral over the remaining sections.

**Energy and flow at the zero.** The physical area factor \(W=\sqrt{1+|\nabla w|^2}\) reaches 1 there. All three principal curvatures have magnitude \(\epsilon/R\). Tension, curvature-squared energy and smooth order-gradient energy stay finite on the compact patch and its smooth joins. A smooth nonzero wave trial field transported by the deformation also has finite energy in the support Hamiltonian. For a material-flow witness take surface number density \(n_0/W\) and downstream coordinate speed \(v^r=D/(4\pi n_0r^2)\). The number rate is \(D\) through every cut, including during a smooth tilt. For a stationary shape the patch kinetic energy in units \(mD^2/(8\pi n_0R)\) is \(1/6+\epsilon^2(1+\beta^2/18)\). This trial current bounds the minimum kinetic energy at the same boundary flux. Minimizing over divergence-free currents with those boundary data gives an irrotational flow for the fixed density and metric. The intake is prescribed here. Neither construction solves the coupled pressure, conversion and wave-support equations. [Wall energy](../../../research/pde_ledger_v2/notes/stage044_v2_unfreeze_prep.md#L19) · [Order balance](../../../docs/moving_interface_conversion_work_proposal.md#L9) · [Wave Hamiltonian](../../../notes/inner_throat/inner_throat_hard_mode.md#L479)

**Condition for protection.** For the initially once-covered mouth, closed conversion sections, the full order drop, no boundary escape and a nonzero in-brane projection of the normal throughout conversion suffice to preserve the unit. Equivalently the material face must have a nonzero tangential depth gradient there, since both magnitudes equal \(\sqrt{1-(N\cdot\hat w)^2}\). Strict local progress into the bulk *along the face relative to its motion* would suffice. Positive total intake allows a patch to level off. The model still needs to establish whether the physical rim enforces this local condition or imposes a substantial energetic cost for violating it. The electric coupling is a separate open requirement. [Current conversion-law scope](../../../research/pde_ledger_v3/steps/S12_order_source.md)

Working calculation on base e0a94594, uncommitted pending independent review.
