# Can the mouth geometry distinguish charge signs independently of the drain?

Two mouths can approach along the three shared brane directions while material leaves one face toward \(+w\) and the other toward \(-w\). The force would have to be transmitted along the slab. In the earlier guide calculation moving the upstream steering reaction onto the mouth reversed its lateral response at the same isolated bend and intake. Such recoil cannot establish the required intrinsic sign. [Flow-guide result](bend_flow_guide.md)

**A static one-mouth check.** Approximate the punctured face by a spherical graph \(w=H(r)\) for \(r\ge a\). Use a common reference normal toward \(+w\) for the two reflected graphs. Assume the face is flat far away and turns vertical at the lip, so \(H'(a)\to-s\infty\), with \(s=+1\) for an upward mouth. Define the signed turning by the integral of the product of the three principal curvatures. With \(\alpha=\arctan H'\), its value is

\[
{\cal I}=\int k_r k_\theta^2\,dA
=4\pi\int_{-s\pi/2}^{0}\sin^2\alpha\,d\alpha=s\pi^2.
\]

Thus \({\cal I}/\pi^2=\pm1\) regardless of radius or radial profile for this class. The intact face contributes zero if it starts and ends flat and is a regular radial graph. The unit values depend on the vertical tangent. A \(45^\circ\) lip has only \(1/2-1/\pi\) of the vertical value. SymPy used the shape operator and Mathematica used the area swept out by the unit normal. [Guarded checks](../../../docs/_measurements/oriented_rim_geometry.json)

**Does this determine what a distant mouth feels?** The profiles \(H=s(C/a)\arcsin(a/r)\) all have the same signed turning at a fixed orientation, yet their far bend is \(sC/r\) for any \(C>0\). As an equilibrium control, minimize the area of one face at uniform tension with a vertical waist. Its first integral gives \(r^2H'/\sqrt{1+H'^2}=-s a^2\). The far bend is then \(s a^2/r\), so doubling the radius quadruples its coefficient. The geometric label survives both checks. Its conversion into a common electric strength remains open. This control omits the second face and the slab’s thickness energy. [Two-face model](two_face_rim_bend.md)

**Decisive force test.** Derive which slab stress or constraint, if any, a turned lip fixes on a surrounding sphere independent of discharge. For well-separated mouths, compute the work of moving the centers along the brane, including rim work, with one common attachment law. For the orientation-dependent part, the energy measured from infinite separation must rise as like mouths approach and fall by the same amount for the mirrored pair. The target is \(+K/d\) and \(-K/d\), respectively, with the same \(K>0\) across radii and discharge preparations. This is a test of a proposed stress law. Simply adding up far-field bending energy omits the rim work, which already reversed the sign in the passive scalar example. [Passive attachment calculation](wall_load_to_bend.md)

Working calculation on base e0a94594, pending independent review.
