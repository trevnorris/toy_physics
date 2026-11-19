# How fast does the receiving bulk’s recoil fall with distance?

**The existing quiescent half-space receiver gives a \(d^{-5}\) reaction away from a positive bend source for an upper mouth.** This is a one-drain test with an imposed distant bend. The receiver carries a prescribed mass rate \(D\) through a compact inner cut. Its reflecting screen lies at \(w=\epsilon f(\mathbf x)\) outside target radius \(a\). Use steady, constant-density potential flow with slip, subtract ambient pressure on the same geometry, and take zero net excess momentum flux at infinity in each bent state. A lower receiver reverses the sign. [Flat reference flow and pressure](rounded_end_bulk_flux.md)

**Force balance.** Use the mouth-center tangent \((\mathbf e_x,\epsilon f_x(0))\) throughout the balance. The outward area vector at the lower wall of the upper receiver is \((\epsilon\nabla f,-1)d^3x\). The impermeable wall transmits pressure. At first order the projection of its momentum flux is
\[
J'(0)=\int_{r>a}p_0(r)[f_x(\mathbf x)-f_x(0)]\,d^3x .
\]
Here \(J\) is the momentum rate **supplied through the inner cut into the receiver**. The screen receives that momentum and a separate inlet supplier takes the opposite recoil. The unknown first pressure variation multiplies a zero flat-screen tangent projection. A regular shift of the projected rim likewise gives zero at this order. If inlet and screen belong to the same translating assembly their two contributions cancel. An actual second discharge, far incident carrier flow or separate support changes the balance.

**Why the inverse-square term vanishes.** Let the applied height have exterior potential \(f=A/|\mathbf x-d\mathbf e_x|\) for a distant positive spherical source of radius \(b\). The mean gradient on a target-centered shell equals the target gradient for \(r<d-b\), varies through the source region, and vanishes for \(r>d+b\). Only the distant pressure tail contributes. For the point-source limit, taken after linearization of a smooth source,
\[
J'(0)=\frac{A}{d^2}C(d),\qquad
C(d)=-4\pi\int_d^\infty r^2p_0(r)\,dr
=\frac{2D^2}{\pi^3\rho d^3}
\left[\frac13-\frac{2e}{5}\frac{a^2}{d^2}+\frac{e^2}{7}\frac{a^4}{d^4}\right].
\]
The dimensionless \(e=0\) or \(1/4\) specifies the previous discharge profiles. For a positive source within radius \(b\) about \(d\mathbf e_x\), \(d>a+b\), superposition bounds \(J'(0)/f_x(0)\) between \(C(d+b)\) and \(C(d-b)\).

**Numerical check.** A uniform source ball with \(a=b=1\) was integrated two ways, over the screen and over the source. From \(d=32\) to \(64\), the force exponents were 5.0019032186 for \(e=0\) and 5.0015853439 for \(e=1/4\). The Python angular quadrature reported roundoff at its requested tolerance. An independent 35-digit Mathematica integral reproduced all eight scaled forces to the twelve displayed decimal places after rounding. Working result on base e0a94594. [Calculations and evidence](../../../docs/_measurements/drain_bend_response.json)
