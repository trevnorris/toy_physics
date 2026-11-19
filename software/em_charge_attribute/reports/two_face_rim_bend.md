# How far does a rim load bend the opposite face?

**Result.** In the project’s minimal two-interface model, a net normal rim load produces a midplane displacement proportional to \(1/r\) in the three brane directions. A positive restoring stiffness for thickness makes the difference between the two face displacements decay over a separate length. At large distance both faces then undergo the same signed bend. The actual draining slab’s stiffnesses and its net load on this mode remain open. [Existing two-interface action](../../../docs/opposite_orientation_throat_coincidence.md#5-quadratic-two-interface-action) · [Width selection](../../../research/pde_ledger_v3/V3_STEP_PLAN.md#L315)

**Test.** Subtract the original one-mouth state at fixed radius, orientation and inward preparation. Let \(a\) be the rim radius in the benchmark, or a matching radius outside a physical core. On \(r\ge a\) use weak spherical perturbations about a homogeneous, reflection-symmetric exterior at rest at infinity. The imposed loads \(F_+,F_-\) are integrated over the rim, both positive toward \(+w\). They include any force transmitted from the excluded region. Upper-face loading means \((F_+,F_-)=(F,0)\); zero lower load allows the lower face to move. This prescribes a test force; a drain's intake alone gives neither load.

Write \(m=(u_++u_-)/2\), \(t=(u_+-u_-)/2\); the change in thickness is \(2t\). The static energy of the template is

\[
E=\tfrac12\int_{r\ge a}[\kappa_m|\nabla m|^2+\kappa_t|\nabla t|^2+B_t t^2]\,d^3x
-(F_++F_-)m(a)-(F_+-F_-)t(a).
\]

Here the gradient stiffnesses are positive and the midplane has no pinning to a fixed \(w\). Variation gives \(\nabla^2 m=0\) and \((\nabla^2-1/\lambda^2)t=0\), with \(\lambda=\sqrt{\kappa_t/B_t}\) when \(B_t>0\). The inner flux in each mode is \(-4\pi a^2\) times its stiffness times its radial slope. Applying the two loads gives

\[
m(r)=\frac{F_++F_-}{4\pi\kappa_m r},\qquad
t(r)=\frac{F_+-F_-}{4\pi\kappa_t r}\frac{e^{-(r-a)/\lambda}}{1+a/\lambda},\qquad
u_\pm=m\pm t .
\]

For a single upper load, \(dm/dr=-F/(4\pi\kappa_m r^2)\). Mirroring the entire preparation to \((0,-F)\) reverses the midplane bend and preserves the thickness change. Equal and opposite face loads have zero midplane response. A self-supported throat could therefore have a smaller or vanishing net bend even with substantial local stress.

**Opposite-face motion.** For two identical independent area tensions \(\sigma\), expanding their areas at small slope gives \(\kappa_m=\kappa_t=2\sigma\) for the gradient terms. As an illustration take \(\lambda=a\). The lower face rises by \(0.333333333\), \(0.689275193\), and \(0.873242123\) times the upper-face rise at \(r/a=1,2,3\). Without thickness restoration, the same equal-tension benchmark leaves the opposite face still everywhere, despite its \(1/r\) midplane signal. With unequal mode stiffnesses even the near-rim displacement ratio depends on both. [Displacement figure](two_face_rim_bend.png) · [SymPy and Mathematica checks](../../../docs/_measurements/two_face_rim_bend.json)

**Physical discriminator.** On a free-faced slab under a calibrated one-face load, measure \((u_+-u_-)/(u_++u_-)\) at several radii. This model predicts attenuation by \(e^{-(r_2-r_1)/\lambda}\) between any two readings with nonzero sum. The fixed-depth carrier solver cannot measure this. The order energy gives an isolated interface tension, while the effective width curvature and normal bulk response of the active slab remain open. A physical electric interpretation also requires the net oriented throat load on an unpinned mode with long-wavelength stiffness proportional to wave number squared. The sign of the force on a second mouth depends on its own coupling. This normal response supplies no in-brane translational inertia for 0PN. [Single-interface tension](../../../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md#L133) · [Active slab limitations](../../stage1_solver/reports/active_drain_width.md) · [Mouth force and bend](one_mouth_electric_bend.md)

Working result on base e0a94594, pending independent review.
