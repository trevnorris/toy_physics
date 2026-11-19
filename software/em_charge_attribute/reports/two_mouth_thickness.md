# Can preferred slab thickness make like mouths repel?

**Two-mouth test.** In the weak-deformation, reflection-symmetric slab, restoring the preferred thickness lets a mouth on one face affect the other face. For two identical, passively supported mouths, like orientations attract at every disjoint separation with nonzero load. At large separation opposite mouths repel. Bringing the mouths together gave the following external work, calculated from the difference of the printed energies in common dimensionless units. [SymPy and Mathematica calculations](../../../docs/_measurements/two_mouth_thickness.json)

| Orientations | Work from \(8a\) to \(4a\), divided by \(P_0^2\) |
|---|---:|
| \(++\) and \(--\) | −0.1157707837 |
| \(+-\) and \(-+\) | +0.1352499088 |

**Setup.** The midplane displacement is \(m=(u_++u_-)/2\). The half-thickness change is \(t=(u_+-u_-)/2\). Reflection flips \(m\) and preserves \(t\). For orientation \(s=\pm1\), the outward mouth depth is \(z=\langle sm+t\rangle\), averaged on a spherical surface of radius \(a\) in the three brane directions. This is a finite-radius load distribution on a continuous slab. The two identical cores have energy \(V(L)=TL+A/L\) at fixed wave action with reflecting length \(L=Y-z>0\). Each reflecting end has the same spring to a mirrored, bulk-fixed reference. The rim depths, reflecting ends and loads relax at each separation. Discharge is absent from this elastic test. [Two-face equations](two_face_rim_bend.md) · [Wave and rear support](supported_mouth_bend.md)

**Sign.** The two stiffnesses are \(\kappa_m,\kappa_t>0\). For thickness relaxation length \(\lambda>0\) and separation \(d>2a\), the mutual displacement per unit load is \(g=s_1s_2G_m+G_t\), where

\[
G_m=\frac{1}{4\pi\kappa_m d},\qquad
G_t=\left[\frac{\sinh(a/\lambda)}{a/\lambda}\right]^2
       \frac{e^{-d/\lambda}}{4\pi\kappa_t d}.
\]

Let \(S\) be the isolated compliance. At equilibrium \(P=V'(L)=z/(S+g)\). The pair field energy at equal outward depths is \(z^2/(S+g)\). Differentiating the *total* energy at the relaxed rim and reflector leaves \(dU/dg=-P^2\), hence the outward force is \(P^2\,dg/dd\). Both kernels decrease with distance. For like mouths both therefore pull inward, for any positive stiffnesses in this model. Thickness gives the same direct contribution for either orientation at a given load. Its force decays exponentially up to powers of distance. The loads themselves do readjust. At long range the leading orientation-dependent outward force is \(-s_1s_2P_0^2/(4\pi\kappa_m d^2)\). Removing thickness restoration with equal gradient stiffnesses makes the opposite-face response vanish. Measuring thickness perpendicular to locally parallel faces gives the same quadratic thickness energy at weak slope.

**Checks and limits.** The example used \(a=\lambda=4\pi\kappa_m=4\pi\kappa_t=k_b=A=1\), \(P_0=1/20\), \(T=1+P_0\), and a common reference giving isolated \(L=1\). SymPy compared relaxed energies. Mathematica integrated the drive work and the force over the receiving shell. All eight rows agreed to ten decimal places. A physical mouth might constrain material strain differently. Its attachment, nonlinear rim geometry and the source of an intrinsic electric sign remain open. Working calculation on base e0a94594, uncommitted pending independent review.
