"""Axial-potential shear model in a homogeneous exterior with assigned flux.

One stationary source is enclosed by spheres oriented radially outwards.
Stress is the Maxwell field stress of the supplied exterior fields.
"""

import sympy as sp


def emit(label, expr):
    print(f"{label} = {expr}", flush=True)
    return expr


def check(label, left, right, failures):
    diff = left - right
    if isinstance(diff, sp.MatrixBase):
        diff = diff.applyfunc(sp.simplify)
        passed = diff.is_zero_matrix
    else:
        diff = sp.simplify(diff)
        passed = diff == 0
    print(f"{label}: left={left}, right={right}, difference={diff}", flush=True)
    if not passed:
        failures.append(label)


def free_waves(b, eps, failures):
    k, y = sp.symbols("k y", positive=True)
    wavevector = sp.Matrix([0, 0, k])
    W = sp.Matrix(sp.symbols("Wx Wy Wz", real=True))
    H = sp.Matrix(sp.symbols("Hx Hy Hz", real=True))
    # Real instantaneous standing-wave amplitudes. Volume averaging gives
    # the same spatial normalization to the kinetic and gradient forms.
    curl_amplitude = wavevector.cross(W)
    kinetic = b*H.dot(H)/2
    gradient = curl_amplitude.dot(curl_amplitude)/(2*eps)
    inertia = sp.hessian(kinetic, H)
    stiffness = sp.hessian(gradient, W)
    operator = stiffness-y*inertia
    emit("plane-wave frame k along z; operator K-y*M (y=omega^2)", operator)
    emit("coefficient eigenvalues omega^2 with multiplicities",
         (inertia.inv()*stiffness).eigenvals())
    longitudinal = (stiffness*W)[2]
    emit("longitudinal equation left b*ddot Wz; right", -longitudinal)
    divB = sp.I*wavevector.dot(b*H)
    emit("longitudinal magnetic divergence to constrain at initial time", divB)
    divBdot = -sp.I*wavevector.dot(stiffness*W)
    emit("time derivative of magnetic divergence from Euler-Lagrange", divBdot)
    emit("time-independent gradient shift; delta D, for amplitude alpha",
         sp.I*wavevector.cross(sp.I*wavevector*sp.Symbol("alpha")))
    check("curl W has zero bulk divergence", wavevector.dot(curl_amplitude), 0, failures)
    return inertia, stiffness


def power_balance(failures):
    E = sp.Matrix(sp.symbols("Ex Ey Ez", real=True))
    H = sp.Matrix(sp.symbols("Hx Hy Hz", real=True))
    # Rows of each derivative array are the differentiation directions.
    def derivatives(name):
        return sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"d{'xyz'[i]}_{name}{'xyz'[j]}"))

    def curl(deriv):
        return sp.Matrix([deriv[1, 2]-deriv[2, 1], deriv[2, 0]-deriv[0, 2],
                          deriv[0, 1]-deriv[1, 0]])

    dE, dH = derivatives("E"), derivatives("H")
    # Ddot = curl H kinematically. Bdot = -curl E by varying W.
    energy_rate = E.dot(curl(dH))-H.dot(curl(dE))
    emit("field energy flux", E.cross(H))
    flux_div = 0
    for i in range(3):
        differentiated = dE.row(i).T.cross(H) + E.cross(dH.row(i).T)
        flux_div += differentiated[i]
    check("local energy rate versus negative divergence of field energy flux",
          energy_rate, -flux_div, failures)


def sphere_charts(eps, failures):
    theta, t, seam = sp.symbols("theta t seam", real=True)
    r, inner, outer = sp.symbols("r r_in r_out", positive=True)
    Q, rotation_scale = sp.symbols("Q lambda", real=True)
    # An axial one-form on the sphere is a(theta)*dphi. Its circulation
    # density a is r*sin(theta)*W_phi. Its theta derivative is the radial
    # flux density per dtheta*dphi. Regular pole caps fix the two constants.
    radial_D = Q/(4*sp.pi*r**2)
    flux_density = r**2*sp.sin(t)*radial_D
    aN, aS = [sp.integrate(flux_density, (t, pole, theta)) for pole in (0, sp.pi)]
    WN, WS = [sp.factor(a/(r*sp.sin(theta))) for a in (aN, aS)]
    emit("north and south azimuthal physical W components", (WN, WS))
    emit("W_phi at each regular pole", (sp.limit(WN, theta, 0, dir="+"),
                                        sp.limit(WS, theta, sp.pi, dir="-")))
    radial_curls = [sp.simplify(sp.diff(a, theta)/(r**2*sp.sin(theta))) for a in (aN, aS)]
    polar_curls = [sp.simplify(-sp.diff(r*azimuthal, r)/r) for azimuthal in (WN, WS)]
    emit("north and south (radial, polar) curls; azimuthal curls vanish by symmetry",
         list(zip(radial_curls, polar_curls)))
    cap_fluxes = [2*sp.pi*sp.integrate(d*r**2*sp.sin(t), limits)
                  for d, limits in zip(radial_curls, [(t, 0, seam), (t, seam, sp.pi)])]
    transition = sp.simplify(WN-WS)
    emit("north minus south W_phi on overlap", transition)
    circulation = sp.factor(2*sp.pi*r*sp.sin(seam)*transition.subs(theta, seam))
    emit("north and south cap fluxes at arbitrary seam", cap_fluxes)
    check("total sphere flux and circulation of overlap difference",
          sp.simplify(sum(cap_fluxes)), circulation, failures)
    common = sp.symbols("C_shared", real=True)
    emit("globally smooth W; north and south Stokes boundary circulations", (common, -common))
    emit("globally smooth W; sum of oriented boundary circulations", common-common)
    shell_integrand = 4*sp.pi*r**2*radial_curls[0]**2/(2*eps)
    emit("static radial field shell energy for r_out>r_in>0",
         sp.factor(sp.integrate(shell_integrand, (r, inner, outer))))

    angle = sp.factor(rotation_scale*transition)
    emit("same-axis ordinary rotation angle difference (north minus south)", angle)
    emit("angle difference in turns on equator", sp.factor(angle.subs(theta, sp.pi/2)/(2*sp.pi)))
    emit("radial derivative of same-axis angle difference at fixed polar angle", angle.diff(r))
    emit("zero assigned flux (radial D, chart circulation, material mismatch)",
         [obj.subs(Q, 0) for obj in (radial_curls[0], circulation, angle)])
    return r, Q, radial_curls[0], circulation


def surface_stress(eps, r, Q, radial_D, failures):
    # Use c=cos(theta). Azimuthal symmetry permits aligning z with the
    # uniform applied field without restricting its direction relative to a sphere.
    c, phi, E0 = sp.symbols("c phi E0", real=True)
    radial = sp.Matrix([sp.sqrt(1-c**2)*sp.cos(phi), sp.sqrt(1-c**2)*sp.sin(phi), c])
    field = radial_D*radial/eps + sp.Matrix([0, 0, E0])
    normal_E = sp.trigsimp(radial.dot(field).expand())
    field_sq = sp.trigsimp(field.dot(field).expand())
    traction = eps*(field*normal_E-radial*field_sq/2)
    emit("electric Maxwell stress traction z-component on the outward sphere",
         sp.factor(traction[2]))
    ring = traction.applyfunc(lambda term: sp.integrate(sp.expand(term), (phi, 0, 2*sp.pi)))
    force = ring.applyfunc(lambda term: sp.simplify(r**2*sp.integrate(term, (c, -1, 1))))
    emit("integrated electric stress (Fx,Fy,Fz), applied field E0 along z", force.T)
    emit("zero applied field control", force.subs(E0, 0).T)
    emit("zero source control", force.subs(Q, 0).T)
    check("mirror Q -> -Q versus applied-field reversal", force.subs(Q, -Q).T,
          force.subs(E0, -E0).T, failures)


def compact_phase(inertia, stiffness, Q, circulation, failures):
    # Extra constitutive hypotheses: a dimensionless period-2*pi phase,
    # fixed nonzero calibration a, and optional spatial phase stiffness.
    a = sp.symbols("a", nonzero=True, real=True)
    phase_stiffness = sp.symbols("Kphase", positive=True)
    turns = sp.symbols("n", integer=True)
    phase_change = a*circulation
    emit("transition phase increment in radians", phase_change)
    emit("single-valued transition exp(i*alpha); assigned flux versus integer n",
         sp.solve(phase_change-2*sp.pi*turns, Q))
    # A regular phase of nonvanishing modulus on each entire disk has a real
    # lift on that disk. Its tangential derivative integrates to zero.
    phi = sp.symbols("phi", real=True)
    north, south = [sp.Function(name)(phi) for name in ("alpha_N", "alpha_S")]
    endpoints = (north-south).subs(phi, 2*sp.pi)-(north-south).subs(phi, 0)
    emit("difference of cap phase circulations from regular lifts", endpoints)
    periodic = {lift.subs(phi, 2*sp.pi): lift.subs(phi, 0) for lift in (north, south)}
    emit("cap phase circulation difference when both lifts extend smoothly", endpoints.subs(periodic))

    # On a flat bulk patch grad(alpha) is longitudinal for k parallel to z.
    # No temporal phase action is supplied. Only the transverse W block is tested.
    W = sp.Matrix(sp.symbols("Wx Wy Wz", real=True))
    phase_gradient = sp.Matrix([0, 0, sp.Symbol("grad_alpha_longitudinal", real=True)])
    covariant = phase_gradient-a*W
    energy = phase_stiffness*covariant.dot(covariant)/2
    add_stiffness = sp.hessian(energy, list(W[:2]))
    frequency_matrix = inertia[:2, :2].inv()*(stiffness[:2, :2]+add_stiffness)
    roots = list(frequency_matrix.eigenvals())
    emit("phase-stiff flat region; transverse omega^2 with multiplicities", frequency_matrix.eigenvals())
    emit("phase stiffness set to zero; transverse omega^2", [root.subs(phase_stiffness, 0) for root in roots])
    k = sp.symbols("k", positive=True)
    emit("phase-stiff transverse limit omega as k approaches zero",
         [sp.sqrt(sp.limit(root, k, 0)) for root in roots])
    check("zero phase stiffness recovers shear transverse operator",
          frequency_matrix.subs(phase_stiffness, 0), inertia[:2, :2].inv()*stiffness[:2, :2], failures)

    # W=0 is prescribed for this kinematic ring comparison. The energy
    # derivative records its load if a positive phase stiffness is present.
    radius = sp.symbols("R", positive=True)
    Wphi = sp.symbols("Wphi", real=True)
    ring_alpha = turns*phi
    ring_density = phase_stiffness*(ring_alpha.diff(phi)/radius-a*Wphi)**2/2
    emit("prescribed W=0 ring; winding number",
         sp.integrate(ring_alpha.diff(phi), (phi, 0, 2*sp.pi))/(2*sp.pi))
    theta = sp.symbols("theta", real=True)
    Wzero = sp.zeros(3, 1)  # physical spherical components in both charts
    flux_density = (sp.diff(radius*sp.sin(theta)*Wzero[2], theta)
                    - sp.diff(radius*Wzero[1], phi))
    emit("prescribed globally W=0; enclosing flux",
         sp.integrate(flux_density, (theta, 0, sp.pi), (phi, 0, 2*sp.pi)))
    emit("prescribed W=0 ring; phase energy per unit section area",
         sp.integrate(radius*ring_density.subs(Wphi, 0), (phi, 0, 2*sp.pi)))
    emit("prescribed W=0 ring; material energy-density derivative with respect to Wphi",
         ring_density.diff(Wphi).subs(Wphi, 0))


def main():
    b, eps = sp.symbols("b eps", positive=True)
    failures = []
    inertia, stiffness = free_waves(b, eps, failures)
    power_balance(failures)
    r, Q, radial_D, circulation = sphere_charts(eps, failures)
    surface_stress(eps, r, Q, radial_D, failures)
    compact_phase(inertia, stiffness, Q, circulation, failures)
    if failures:
        raise AssertionError(", ".join(failures))


if __name__ == "__main__":
    main()
