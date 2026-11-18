# How fast must a trapped shear wave follow the rim at 0PN?

**What the existing equations bound.** The wall action gives the traction transmitted to the rim by a specified wall deformation. The baseline carrier calculation measures the feeding-side reaction to prescribed motion. Neither gives the relative lateral restoring stiffness \(k\) of a trapped shear envelope or its effective inertia \(M_w\) on the same draining state. To see the freedom, add \(\lambda b(a)|X-Y|^2/2\), \(b(a)>0\), to a centered configuration with radius \(a\), rim center \(X\) and envelope center \(Y\). Every \(\lambda>0\) preserves the centered energy, radial equilibrium and radial curvature while changing the relative stiffness to \(\lambda b(a)\) in this example. Thus radial support and stability allow an arbitrarily slow response within the open effective laws. The rapid oscillation sustaining the throat, its leakage lifetime and a shear crossing time do not supply a lower bound on the lateral restoring rate. [Mouth traction and open inputs](one_mouth_translation_equation.md#work-exchanged-when-the-opening-shifts) · [Radius and length stability calculation](pathA_26_derrick.md) · [Measured carrier reaction](../../../docs/drain_driven_initial_acceleration.md#L229)

**Slowly varying gravity.** In the conservative two-coordinate example the wave obeys \(M_w\ddot Y=k(X-Y)\), with the exterior drive on the rim. Define \(\omega_f^2=k/M_w\). For a nonzero harmonic component of the rim motion at angular frequency \(\nu<\omega_f\), with no free relative oscillation, the fractional acceleration mismatch is

\[
\frac{|\ddot Y_\nu-\ddot X_\nu|}{|\ddot X_\nu|}
=\frac{\nu^2}{\omega_f^2-\nu^2}.
\]

For frequencies through \(\nu_g\) and allowed fractional error \(\epsilon\) per component, the exact condition is \(\omega_f\ge\nu_g\sqrt{1+1/\epsilon}\). At an illustrative 1% tolerance the multiplier is \(\sqrt{101}\). For a prescribed circular comparison of separation \(d\) and relative speed \(v\), the leading companion flow turns at \(\nu_g=v/d\) by kinematics. A noncircular comparison needs its drive spectrum. The tracking criterion assumes stable linear response on the same inward branch and that the retained drive and envelope rates lie below the fast standing-wave scale. [Two independent algebraic checks](../../../docs/_measurements/wave_follow_bound.json)

**Release from rest.** A sudden constant rim force excites the free relative oscillation at \(\Omega^2=k(1/M_r+1/M_w)\). For positive rim-associated field inertia \(M_r\), zero initial slip and velocity, and nonzero force, \(X_{\rm cm}\) denotes the displacement of the center of inertia. For \(t>0\),

\[
\frac{|X-Y|}{|X_{\rm cm}|}=\frac{2(1-\cos\Omega t)}{(\omega_f t)^2}
\le\frac{4}{(\omega_f t)^2}.
\]

This bounds separation *relative to the distance fallen*. The accelerations keep oscillating in the lossless example. If the envelope must also stay within an absolute length \(\ell\) of the rim, the same solution gives \(|X-Y|\le 2|A_{\rm cm}|/\omega_f^2\), where \(A_{\rm cm}\) is the constant center-of-inertia acceleration. The physical release would need its own initial state, damping and leakage response. [Push histories and assumptions](trapped_wave_push.md)

**A bound using only data at the mouth.** For this two-coordinate law, a weak push from rest gives the ratio of Laplace-transformed total net force to rim acceleration, \(m(s)=\widetilde F(s)/(s^2\widetilde X(s))=M_r+M_w\omega_f^2/(\omega_f^2+s^2)\). Here \(s>0\) is a real weighting rate, and \(m(0)\) denotes its slow limit. On the same isolated drain, measure the calibrated force and geometric center. Subtract or independently account for other external loads and dynamic channels. Form

\[
R=\frac{m(0)-m(s)}{m(0)-m(2s)},\qquad
\omega_f^2=4s^2\,\frac{1-R}{4R-1},\qquad \tfrac14<R<1.
\]

The inferred rate decreases with \(R\), so a resolved upper confidence bound for \(R\) inside this interval yields a lower bound on \(\omega_f\), conditional on this law. Test the prediction at an unused rate and at half drive. Unresolved differences, additional modes or losses preclude this three-reading inference without further information. Other mechanisms could reproduce the same rim response, so attribution to the trapped wave still needs an independent identification. The archived radius-one run shows less than 0.1% variation in the *feeding-side* reaction between the slow limit and harmonic frequency 0.01 at reservoir extent 48. It supplies no measured value of \(R\). The physical comparison therefore remains open. Working result on base cda9a8eb, pending independent review. [Archived frequency comparison](../../../docs/drain_driven_initial_acceleration.md#L237) · [Run manifest](../../../docs/_measurements/wave_follow_bound.json)
