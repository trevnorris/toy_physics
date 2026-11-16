# Slow force on a companion from a translating drain

Let \(C=\delta F_{s,x}/A_{t,x}\), with the target at \(x=0\), the stationary companion at \(x=-d\), and positive \(x\) pointing toward the target from the companion. The held force is subtracted at the same mouth positions. Both mouths drain inward. Work to first order in the imposed acceleration and leading order in the intake rates \(J_t,J_s\). The calculation takes the slow limit, removes the lateral reservoir, then takes large separation at fixed radii and slab depth \(H=1\). The solver's brane coordinates have three spatial dimensions. The coordinate normal to the brane spans the slab. [Force definition and field checks](drain_driven_initial_acceleration.md#predicting-the-companion-load-at-separation-6)

## Response on the brane side of the mouths

Let \(W_i\) be a mouth's intake per unit total rate, with \(\int W_i\,d^3x=1\). At uniform density the continuum profile is proportional to \(g_i=[1-r^2/a_i^2]_+^3\) inside its radius \(a_i\), and zero outside. Let \(D_N=-\Delta\) on the slab with reflecting faces, \(\gamma\) denote evaluation at the mouth face, and \(\gamma^*\) insert a face source in the weak equation. A unit translation changes the face source by \(S_i=-\partial_x W_i\). Its phase response per unit current is \(p_i=D_N^{-1}\gamma^*S_i\).

The carrier term \((|\psi|^8-1)\psi/4\) has derivative \(2\operatorname{Re}\delta\psi\) at \(\psi=1\). The translation equations give the real velocity tangent \(u_R=-2T_g p_t\) and imaginary acceleration tangent \(h_I=2D_N^{-1}u_R\), per unit target rate. Here \(T_g=(-\Delta+4)^{-1}\) has \(\partial_n u_R=-(g_t+g_s)u_R\) on the mouth face and reflection below. This condition comes from the solver's restoring law for the density amplitude at the openings. To bilinear order the companion samples \(J_sJ_t\int W_s\partial_x\gamma h_I\,d^3x\). The explicit mixed-stress variation in the mouth law is higher order at uniform density. Integrating by parts gives the reciprocal formula

\[
\frac{C}{J_sJ_t}=-4\langle p_s,T_g p_t\rangle_{\rm slab}.
\]

This uses the phase, density amplitude and mouth momentum flux. No force between the two centers enters these equations. [Carrier operator](../software/stage1_solver/tools/lateral_reservoir_field.py) · [Intake and restoring laws](../software/stage1_solver/tools/lateral_mouth_radius.py) · [Motion tangents and force components](../software/stage1_solver/tools/lateral_mouth_acceleration.py)

## Density response and finite apertures

First set \(g=0\) in the amplitude boundary condition, keeping the two intake profiles. At distances large compared with the depth, the depth-uniform mode gives a factor \(1/H\). With \(k\) the brane wave number, the radial kernel whose axial second derivative gives \(C H/(J_sJ_t)\) has transform

\[
\frac{4\widehat W_s(k)\widehat W_t(k)}{k^4(k^2+4)}
=\frac{1}{k^4}
 -\left(\frac14+\frac{\langle r^2\rangle_s+\langle r^2\rangle_t}{6}\right)\frac{1}{k^2}
 +O(1).
\]

The inverse transforms of \(k^{-4}\) and \(k^{-2}\) are \(-d/(8\pi)\) and \(1/(4\pi d)\), up to an irrelevant constant in the former. At distinct points on the line joining the mouths their axial second derivatives are zero and \(1/(2\pi d^3)\), respectively. Thus the inverse cube comes from a spatial derivative of the carrier response. The \(1/4\) comes from the \(k^2+4\) density operator. Each mouth contributes its second moment because flow is emitted and sampled over finite footprints. Integrating the stated radial profile with the three-dimensional volume element gives \(\langle r^2\rangle_i=3a_i^2/11\).

## Restoring condition at each opening

The same-order amplitude correction can be obtained from two isolated-mouth problems. Write \(T_0=(D_N+4)^{-1}\), \(T_{g_i}\) for the operator with just mouth \(i\)'s restoring condition, and \(w_i=\gamma T_{g_i}p_i\). Define

\[
\mu_i=\int (x-X_i)\,g_i\, w_i\,d^3x .
\]

To see how this enters, use \(T_0-T_g=T_0\gamma^*g\gamma T_g\) in the reciprocal formula. The screened transfer between the two localized restoring conditions decays exponentially with separation, so to leading algebraic order each contributes separately. The boundary operator \(g_i(I+\gamma T_0\gamma^*g_i)^{-1}\) is symmetric. Consequently the term at mouth \(i\) is \(4\int (\gamma T_0 p_j) g_i w_i\,d^3x\) for the other mouth \(j\). Across mouth \(i\), \(p_j\) is a distant phase dipole and \(\gamma T_0 p_j\simeq p_j/4\). Its constant part integrates to zero against the odd function \(g_i w_i\). Its axial slope is \(-1/(2\pi H d^3)\). Each mouth therefore adds \(-\mu_i/(2\pi H d^3)\) to the response per product of rates. For a radially decreasing intake, \(S_i\) is positive on the positive-\(x\) half of the mouth. The odd dipole solution there has positive sign under the Poisson and screened Robin inverses, giving \(\mu_i>0\).

For numerical evaluation the Neumann surface multipliers are \(B_m(k)=\coth(H\sqrt{k^2+m^2})/\sqrt{k^2+m^2}\). The isolated face problem is

\[
 w_i+B_2(g_i w_i)=\tfrac14(B_0-B_2)S_i .
\]

The calculation projects this equation onto radial dipole basis functions for each mouth, using the spherical Bessel transform of order one. It varies the number of basis functions, quadrature nodes and wave-number cutoff. An additional calculation solves the linearized equations directly in each archived finite cylinder with the calibrated mouth shapes. [Coefficient calculation and checks](../software/stage1_solver/tools/lateral_companion_inverse_cube.py) · [Hankel transform convention](https://dlmf.nist.gov/10.22.v)

## Result

In solver units the leading weak-flow, large-separation coefficient is

\[
\boxed{\displaystyle
 C(d)=-\frac{J_sJ_t}{2\pi H\,d^3}
 \left[\frac14+\frac{a_s^2+a_t^2}{22}+\mu_s+\mu_t\right]
 + o(d^{-3})\quad\text{at bilinear order}.}
\]

The negative sign refers to the extra load when the target accelerates away from its companion. The steady attraction is a separate observable. The continuum formula uses the uniform-density intake profile. The finite-grid check retains the spatial calibration from the nonlinear runs and compares at their measured steady rates. Higher powers of the rates can change the coefficient in the finite-flow model. [Measured comparison](drain_driven_initial_acceleration.md#deriving-the-size-of-the-inverse-cube-response)
