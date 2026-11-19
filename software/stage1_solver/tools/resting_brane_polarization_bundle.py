"""Berry curvature of conditional circular polarization lines on a saved mouth section.

The added hypothesis locks the transverse plane to the projected normal.
Use A=i*e^dagger de, outward orientation, and helicity h relative to that
normal. No electric coupling or helicity selection enters the construction.
"""

import sympy as sp


def emit(label, value):
    print(f"{label} = {value}", flush=True)
    return value


def check(label, left, right, failures):
    difference = sp.simplify(left-right)
    print(f"{label}: left={left}, right={right}, difference={difference}", flush=True)
    if difference != 0:
        failures.append(label)


def berry_one_form(vector, coordinates):
    return sp.Matrix([
        sp.trigsimp(sp.expand_complex(sp.I*vector.conjugate().dot(vector.diff(coord))))
        for coord in coordinates
    ])


def circular_lines(theta, phi, failures):
    direction = sp.Matrix([sp.sin(theta)*sp.cos(phi), sp.sin(theta)*sp.sin(phi), sp.cos(theta)])
    tangent = direction.diff(theta)
    azimuthal = (direction.diff(phi)/sp.sin(theta)).applyfunc(sp.simplify)
    area = sp.trigsimp(direction.dot(tangent.cross(direction.diff(phi))))
    curvatures = {}
    for h in (1, -1):
        vector = (tangent+sp.I*h*azimuthal)/sp.sqrt(2)
        emit(f"h={h} polarization norm", sp.trigsimp(vector.conjugate().dot(vector).expand()))
        connection = berry_one_form(vector, (theta, phi))
        emit(f"h={h} frame A(theta,phi)", list(connection))
        # At the north and south poles the frame's phi phases are opposite.
        # e -> exp(i chi)*e shifts this connection by -d chi.
        potentials = []
        for pole, pole_theta, winding in (("north", 0, h), ("south", sp.pi, -h)):
            phase = sp.exp(sp.I*winding*phi)
            pole_vector = (phase*vector.subs(theta, pole_theta)).applyfunc(
                lambda item: sp.trigsimp(sp.expand_complex(item)))
            emit(f"h={h} {pole} chart phase multiplier and regular pole vector",
                 (phase, list(pole_vector)))
            potential = connection-sp.Matrix([0, winding])
            emit(f"h={h} {pole} chart (A_theta,A_phi)", list(potential))
            potentials.append(potential)
        curvature = potentials[0][1].diff(theta)-potentials[0][0].diff(phi)
        emit(f"h={h} chart curvature F_theta_phi", curvature)
        sphere_integral = sp.integrate(curvature, (theta, 0, sp.pi), (phi, 0, 2*sp.pi))
        seam_integral = sp.integrate((potentials[0]-potentials[1])[1], (phi, 0, 2*sp.pi))
        check(f"h={h} integrated curvature versus chart circulation", sphere_integral, seam_integral, failures)
        emit(f"h={h} direction sphere first Chern number", sphere_integral/(2*sp.pi))
        curvatures[h] = sp.simplify(curvature/area)
    # The trace of a matrix-connection commutator vanishes, also in the degenerate doublet.
    emit("complete doublet trace curvature F_theta_phi", sum(curvatures.values())*area)
    return curvatures


def jacobian(raw, c, phi):
    # c=cos(theta) decreases on the positively oriented theta path.
    # Integrate columns (phi,c) over c=-1..1 for the signed sphere area.
    volume = sp.factor(sp.trigsimp(raw.dot(raw.diff(phi).cross(raw.diff(c))).expand()))
    length_sq = sp.factor(sp.trigsimp(raw.dot(raw).expand()))
    return volume/length_sq**sp.Rational(3, 2), sp.sqrt(length_sq)


def mouth_section(phi, area_multiples, failures):
    c = sp.symbols("c", real=True)
    t, delta = sp.symbols("t delta", positive=True)
    beta = sp.symbols("beta", nonnegative=True)
    a = sp.symbols("a", real=True, nonzero=True)
    radial = sp.Matrix([sp.sqrt(1-c**2)*sp.cos(phi), sp.sqrt(1-c**2)*sp.sin(phi), c])
    raw = t*radial-sp.Matrix([0, 0, beta])
    density, norm = jacobian(raw, c, phi)
    mirror_density, _ = jacobian(-raw, c, phi)
    emit("mouth signed area density, variables phi and c=cos(theta)", density)
    check("antipodal signed area density versus negative original", mirror_density, -density, failures)
    # The azimuth in the original map is phi. The signed area primitive is n_z.
    primitive = raw[2]/norm
    check("mapped n_z derivative versus direct triple product", primitive.diff(c), density, failures)
    endpoints = [sp.factor(primitive.subs(c, value), deep=True) for value in (-1, 1)]
    emit("mouth n_z at south and north input poles, regular sections", endpoints)
    tilt_cases = [("beta<t", {t: beta+delta}), ("beta>t", {beta: t+delta})]
    for label, substitution in tilt_cases:
        lower, upper = [sp.simplify(endpoint.subs(substitution)) for endpoint in endpoints]
        signed_area = 2*sp.pi*(upper-lower)
        emit(f"{label} integrated signed solid angle", signed_area)
        for h, multiple in area_multiples.items():
            F = multiple*signed_area
            C = sp.simplify(F/(2*sp.pi))
            fixed = sp.simplify(-multiple*signed_area/(2*sp.pi))
            flipped = sp.simplify(-area_multiples[-h]*signed_area/(2*sp.pi))
            emit(f"{label} h={h}: Chern (original,mirror fixed h,mirror flipped h)",
                 (C, fixed, flipped))
            emit(f"{label} h={h}: assigned original flux only if a*W_tangent=A", sp.simplify(F/a))
    singular = raw.subs(beta, t)
    emit("beta=t unnormalized projected normal at north input pole", list(singular.subs(c, 1)))
    unit = singular/norm.subs(beta, t)
    emit("beta=t limiting unit normal approaching north along fixed phi",
         [sp.simplify(sp.limit(component, c, 1, dir="-")) for component in unit])


def relative_phase(failures):
    # Additional surface constitutive hypothesis. W and A are locally smooth
    # at vortex cores in regular charts. Eta shifts under both gauge changes,
    # so v is a globally defined one-form on the ordered punctured surface.
    a = sp.symbols("a", real=True, nonzero=True)
    kappa, amplitude_sq, v_sq = sp.symbols("kappa absPsi_sq norm_v_sq", nonnegative=True)
    energy_density = kappa*amplitude_sq*v_sq/2
    emit("added local phase energy per unit surface area", energy_density)
    chi_W, chi_h = sp.symbols("dchi_W dchi_h", real=True)
    # Transform aW -> aW+dchi_W, A -> A+dchi_h.
    delta_eta = sp.Symbol("d_delta_eta", real=True)
    emit("d(delta eta) required by v invariance", sp.solve(delta_eta-chi_W+chi_h, delta_eta))
    q = sp.Symbol("q", real=True)
    C, winding = sp.symbols("C_h sum_m", integer=True)
    x1, x2 = sp.symbols("x1 x2", real=True)
    eta, W1, W2, A1, A2 = [sp.Function(name)(x1, x2) for name in ("eta", "W1", "W2", "A1", "A2")]
    v1 = eta.diff(x1)-a*W1+A1
    v2 = eta.diff(x2)-a*W2+A2
    dv = (v2.diff(x1)-v1.diff(x2)).expand()
    curl_W = W2.diff(x1)-W1.diff(x2)
    curl_A = A2.diff(x1)-A1.diff(x2)
    coef_W, coef_A = dv.coeff(W2.diff(x1)), dv.coeff(A2.diff(x1))
    check("smooth dv from local v versus chart curvatures", dv, coef_W*curl_W+coef_A*curl_A, failures)
    # Positive small loops go around vortex cores as boundaries of the
    # excised disks. The boundary orientation of the punctured surface reverses them.
    stokes_boundary = -2*sp.pi*winding
    smooth_dv = coef_W*q+coef_A*2*sp.pi*C
    emit("integral of dv on punctured S, shrinking cores; boundary and bulk",
         (stokes_boundary, smooth_dv))
    flux, = sp.solve(stokes_boundary-smooth_dv, q)
    emit("closed-surface a*q/(2*pi) including positively oriented vortex windings",
         sp.simplify(a*flux/(2*sp.pi)))
    emit("nowhere-zero coherence on closed S; no vortices, a*q/(2*pi)",
         sp.simplify((a*flux/(2*sp.pi)).subs(winding, 0)))
    emit("phase energy at vanishing coherence for arbitrary finite norm_v_sq",
         energy_density.subs(amplitude_sq, 0))


def main():
    theta, phi = sp.symbols("theta phi", real=True)
    failures = []
    multiples = circular_lines(theta, phi, failures)
    mouth_section(phi, multiples, failures)
    relative_phase(failures)
    if failures:
        raise AssertionError(", ".join(failures))


if __name__ == "__main__":
    main()
