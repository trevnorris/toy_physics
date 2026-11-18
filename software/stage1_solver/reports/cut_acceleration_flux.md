# Sideways momentum crossing a cut just beyond a slowly accelerating mouth

Orient the cut \(C\) downstream and take \(+x\) along the imposed rim acceleration \(A\). Material exports sideways momentum per unit time
\[
\Phi_x^C=\int_C [j_m v_x-\sigma_{xN}]\,dS,\qquad
j_m=\rho(\mathbf v-\mathbf V_C)\cdot\mathbf N .
\]
Here \(v_x\) is the transported lab velocity, \(\mathbf V_C\) the cut velocity and \(\sigma\) the tensile-positive stress. In an instantaneous rest frame the common sideways boost is zero. In the neutral carrier on a flat section normal to \(w\), ordinary pressure has zero \(xw\) projection. The density-gradient stress contributes \(-\hbar^2 n\,\partial_x\partial_w\ln n/(4m)\) to the flux density. A tilted cut admits a pressure projection. If the cut also intersects the wall or bound shear field, include its transmitted momentum and traction once for the complete load. Momentum stored between this cut and the rim also enters a rim-force inference. [Moving-cut balance](one_mouth_moving_momentum.md) · [Wall reaction](rim_end_load.md)

**Nearest archived interior cut.** The finite neutral-carrier receiver has a straight entrance of length \(L\). At \(w=L/4\), 128 cells per mouth radius, the following are linear coefficients per unit imposed acceleration in solver mass units. Both start from the same inward number rate \(0.125\). The whole rig accelerates. The second outlet condition additionally imposes relative sideways bath velocity \(A\tau\), with \(\tau=1\) solver time unit.

| Far bath | Carried | Density-gradient stress | Total | Total minus matching zero-drain control |
|---|---:|---:|---:|---:|
| Uniform phase | −0.0000786533 | +0.167814125 | +0.167735471 | +0.0000356997 |
| Controlled sideways drift | +0.001102305 | +0.168415396 | +0.169517700 | +0.0018179285 |

The 64-cell totals were \(+0.167653903\) and \(+0.169432737\). The table samples the carrier on one interior plane, at finite resolution. Neither receiver has been matched to the physical mouth. [Original guarded fields and loads](../../../docs/_measurements/open_sleeve_lateral_balance.txt) · [Sums and control subtraction](../../../docs/_measurements/cut_acceleration_flux.json)

**Does the whole-cut integral cancel the missing response?** Use a straight receiver as a local uniqueness test. In dimensionless units \(m=n_0=v_{w,0}=1,\ c_s^2=35/24,\ \hbar^2/(4m^2)=1/24\). Let \(Y=j_1(r)\cos\theta/j_1(a)\) on the round three-dimensional section of radius \(a\), with \(j_1\) the spherical Bessel function and \(j_1'(a)=0\) at its first positive stationary point. This dipole has transverse Laplacian \(-Y\) and zero radial derivative at the side. Set \(w=0\) at the measurement cut. Add to a common linear acceleration response the density and velocity-potential perturbations
\[
(\delta n,\delta\Phi)=bA Y\left(-3e^{-2w}+8e^{-3w},\;-2e^{-2w}+3e^{-3w}\right).
\]
The addition solves the homogeneous linear carrier equations because the common acceleration forcing cancels in the difference. It decays downstream, preserves reflecting side conditions and changes the normal current pointwise by zero at the cut. Per \(bA\), it adds \(\partial_xY\) by carriage and \(3\partial_xY/4\) by stress. Integrating over the entire ball gives \(7\pi a^2/3\), since \(Y=\cos\theta\) on the rim. Thus even the normal-current profile and a common far response leave a finite force difference if mouth matching leaves \(b\) free. The mouth-level law would have to select the amplitude \(b\). [Parent receiver modes](../../../docs/_measurements/sleeve_outlet_boundary_data.txt) · [Independent hydrodynamic and complex-field checks](../../../docs/_measurements/cut_acceleration_flux.json)

**Physical input.** On the actual downstream section, determine the acceleration-dependent flux-weighted sideways velocity and the integrated sideways stress, including any wall traction across its perimeter. A direct total-current reading suffices. Use opposite small accelerations on the same inward branch, accounting for static displacement, velocity and any history dependence before extracting the slow coefficient. The normal intake and far backpressure alone leave its magnitude and sign open. Working result on base cda9a8eb, pending independent review.
