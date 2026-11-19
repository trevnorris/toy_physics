# Does one inward drain leave a gravity-range bend?

**In the open-receiver benchmark the slab bends locally. Its \(1/r\) displacement coefficient is exactly zero when the mouth recoil and all face pressures act on the same slab midpoint.** The surviving displacement falls as \(r^{-4}\), giving a slope of order \(r^{-5}\). The physical converting mouth still needs a measurement or a law identifying where its normal reaction goes. [Guarded SymPy and Mathematica results](../../../docs/_measurements/isolated_mouth_bend.json)

**One-mouth test.** Take \(+w\) toward the receiving bulk. Work to leading order in weak steady flow about flat faces. At distances \(r\ge a\) the feed occupies a reflecting slab of thickness \(H\), and the receiver occupies a four-dimensional half-space. Both fluids are incompressible and at rest at infinity. Match the total inward and outward mass rate \(D>0\). Keep all turning, conversion and compact mechanical contacts inside the cuts at \(a\). Their local fields have not been matched. Use the unpinned linear midpoint of stiffness \(\kappa_m>0\) and set its displacement \(m(\infty)=0\). [Two-face action](../../em_charge_attribute/reports/two_face_rim_bend.md) · [Receiver geometry and discharge](../../em_charge_attribute/reports/rounded_end_bulk_flux.md)

On the feed side take the potential \((A+B e^{-kr}\cos kz)/r\), where \(0\le z\le H\), \(k=\pi/H\), \(A=D/(4\pi\rho_f H)\), and \(\rho_f\) is the feed density. The adjustable mode \(B\) carries normal momentum near the mouth at fixed intake. The receiver has density \(\rho_4\) and the rounded outlet's uniform or quadrupole-shaped discharge, parameterized by \(\epsilon\). The signed outward \(w\)-momentum rates on the feeding cylinder and receiving hemisphere, including pressure, are

\[
\Phi_f(r)=8\pi\rho_f AB\,e^{-kr}/r,\qquad
\Phi_b(r)=\frac{2D^2}{\pi^3\rho_4 r^3}
\left(\frac13-\frac{2\epsilon a^2}{5r^2}+\frac{\epsilon^2a^4}{7r^4}\right).
\]

**Where the momentum goes.** Write \(\Phi=\Phi_f+\Phi_b\). In a stationary, unsupported assembly with complete transfer to the common midpoint, compact recoil is \(-\Phi(a)\). The intact faces between \(a\) and \(r\) receive \(\Phi(a)-\Phi(r)\). Thus the total midpoint load inside \(r\) is \(Q_m(r)=-\Phi(r)=-4\pi\kappa_m r^2m'(r)\). The internal mouth cut cancels when both sides are included. SymPy integrated momentum crossing the outer cuts. Mathematica integrated pressure on both feed faces and the receiving plane and reconstructed the displacement from spherical shells. [Cut calculation](../tools/isolated_mouth_bend.py) · [Face calculation](../tools/isolated_mouth_bend.wl) · [Three-dimensional Green function](https://www.cfm.brown.edu/people/dobrush/am34/Mathematica/ch6/laplace.html)

\[
\lim_{r\to\infty} r\,m(r)=0,\qquad
\lim_{r\to\infty} r^4m(r)=-\frac{D^2}{24\pi^4\kappa_m\rho_4}.
\]

**Physical requirement.** In the general budget, \(Q_m=F_{w,\mathrm{applied}}-\dot P_w-\Phi_f-\Phi_b-T\), where \(T\) is any additional signed outward momentum rate through a separate channel, excluding midpoint traction itself. A steady isolated source with \(T\to T_\infty\ne0\) would have \(m(r)\sim -T_\infty/(4\pi\kappa_m r)\). This control specifies a transfer. The existing fixed-face solver supplies neither that transfer nor a free midpoint reading. The decisive datum is a nonzero plateau in \(-4\pi\kappa_m r^2m'(r)\) over expanding cuts, together with its independently accounted opposite reaction. Until then this branch supplies no additional Newtonian force or free-fall coefficient. [Fixed-face control](../../em_charge_attribute/reports/joined_mouth_normal.md) · [0PN force and inertia](../../../docs/zero_pn_gravity.md)

Working calculation on checkpoint 3666e56f, pending independent review.
