"""Fourier response of an isotropic polarizable elastic slab at rest.

Use exp(i*k*z-i*w*t), q=k^2, y=w^2 and k>0. Each transverse block
represents one polarization. A and u are driven by electric current J and
force density F. Spectral weights are coefficients of 1/(y_mode-y).
"""

import argparse

import sympy as sp


q, k, w, y, x = sp.symbols("q k w y x", positive=True)
eps, nu, rho, K, mu, s = sp.symbols("eps nu rho K mu s", positive=True)
g0, f = sp.symbols("g0 f", real=True)
g, T, S = sp.symbols("g T S", real=True)
z = sp.symbols("z", real=True)
failures = []


def show(label, expr):
    print(f"{label} = {expr}", flush=True)
    return expr


def compare(label, left, right):
    difference = sp.factor(left - right)
    print(f"{label}: left={left}, right={right}, difference={difference}", flush=True)
    if difference != 0:
        failures.append(label)


def even_frequency(expr):
    return sp.expand(expr).subs(w**2, y)


def coefficient(expr, variable, degree):
    return sp.factor(sp.expand(expr).coeff(variable, degree))


def next_coefficient(poly, trial, degree, variable=q):
    eq = coefficient(poly.subs(y, trial), variable, degree)
    return [sp.factor(value) for value in sp.solve(eq, z)]


def weight(numerator, determinant, branch, numerator_order, denominator_order):
    # The input branch contains powers through the first nonzero numerator.
    top = coefficient(numerator.subs(y, branch), q, numerator_order)
    bottom = coefficient(determinant.diff(y).subs(y, branch), q, denominator_order)
    return sp.factor(-top / bottom)


def transverse_operator():
    # Euler-Lagrange derivatives of one transverse Fourier oscillator.
    A, u, Ad, ud, Add, udd, J, F = sp.symbols("A u Ad ud Add udd J F")
    lag = (eps*Ad**2 - nu*q*A**2 + rho*ud**2 - T*u**2)/2 - g*u*Ad + J*A + F*u
    coordinates = sp.Matrix([A, u])
    velocities = [Ad, ud]
    accelerations = [Add, udd]
    time_vars = list(coordinates) + velocities
    rates = velocities + accelerations
    fourier = {Ad: -sp.I*w*A, ud: -sp.I*w*u, Add: -w**2*A, udd: -w**2*u}
    equations = []
    for coord, vel in zip(coordinates, velocities):
        momentum = lag.diff(vel)
        dt = sum(momentum.diff(var)*rate for var, rate in zip(time_vars, rates))
        equations.append((dt - lag.diff(coord)).subs(fourier))
    return sp.Matrix(equations).jacobian(coordinates)


def pinned_flux(matrix, low_slope):
    # With F=0 the matter row fixes u/A. Form u/E while both numerator
    # and denominator still use the same w, then set w=sqrt(y) everywhere.
    u_per_E = sp.factor((-matrix[1, 0]/(sp.I*w*matrix[1, 1])).subs(w, sp.sqrt(y)))
    show("pinned u/E for F=0", u_per_E)
    show("pinned low u/E", u_per_E.subs(y, low_slope*q).subs(q, 0))
    # Real plane-wave fluxes are nu*w*k*|A|^2/2 and mu*w*k*|u|^2/2.
    # Both gradients act on fields with the same k. The dipole term is local.
    flux_ratio = sp.factor(mu*y*u_per_E**2/nu)
    show("pinned shear/EM time-averaged flux", flux_ratio)
    low_flux = show("pinned low (shear/EM flux)/q",
                    sp.factor(sp.limit(flux_ratio.subs(y, low_slope*q)/q, q, 0)))

    # On the free branch the source-free Maxwell row gives a second u/A.
    # It can disagree with the driven matter-row ratio away from the branch.
    u_per_A = -matrix[0, 0]/matrix[0, 1]
    maxwell_flux = (mu*u_per_A*sp.conjugate(u_per_A)/nu).subs(w, sp.sqrt(y))
    compare("low flux coefficient from matter and Maxwell rows", low_flux,
            sp.factor(sp.limit(maxwell_flux.subs(y, low_slope*q)/q, q, 0)))
    sample = {q: 1, s: 1, mu: 1, eps: 1, rho: 1, g0: 1, nu: sp.Rational(5, 6)}
    print(f"pinned finite-q parameters = {sample}", flush=True)
    sample_modes = sp.solve(even_frequency(matrix.det()).subs(sample), y)
    sample_y = show("pinned finite-q low y", min(sample_modes))
    matter_value = sp.factor(flux_ratio.subs(sample).subs(y, sample_y))
    maxwell_value = sp.factor(maxwell_flux.subs(sample).subs(y, sample_y))
    compare("finite-q on-shell flux ratio from matter and Maxwell rows", matter_value, maxwell_value)
    if matter_value.is_positive is not True:
        failures.append("finite-q on-shell shear flux")


def main(pinned_flux_only=False):
    M = transverse_operator()
    pinned = M.subs({g: g0, T: s+mu*q})
    vp, = next_coefficient(even_frequency(pinned.det()), z*q, 1)
    if pinned_flux_only:
        pinned_flux(pinned, vp)
        return

    show("transverse Fourier operator, rows (J,F), columns (A,u)", M)
    D = show("D, transverse determinant", sp.collect(even_frequency(M.det()), y))
    adj = M.adjugate().subs(w**2, y)
    show("response numerators (A,u)/(J,F), common denominator D", adj)
    nA, nuu = adj[0, 0], adj[1, 1]

    # Gauss's law and matter force balance keep phi as a constraint.
    phi, charge, u, F = sp.symbols("phi charge u F")
    EL = -sp.I*k*phi
    long_matrix, long_source = sp.linear_eq_to_matrix(
        [sp.I*k*(eps*EL + g*u) - charge, (S-rho*y)*u - g*EL-F], [phi, u])
    responseL = (long_matrix.inv()*long_source).applyfunc(sp.factor)
    phi_static = responseL[0].diff(charge).subs(y, 0)
    long_gap, = sp.solve(long_matrix.det(), y)
    show("longitudinal (phi,u) for (charge,F)", responseL)
    show("longitudinal free y(k)", sp.factor(long_gap))
    specify = {S: s+K*q, T: s+mu*q, g: g0+f*q}
    kernel = sp.factor(phi_static.subs(k**2, q).subs(specify))
    show("static phi/charge with g(k)", kernel)
    show("static A/J", sp.factor((nA / D).subs(y, 0)))

    acoustic = {s: 0, f: 0}
    Da = D.subs(specify).subs(acoustic)
    Aa = nA.subs(specify).subs(acoustic)
    Ua = nuu.subs(specify).subs(acoustic)
    show("unpinned dipole static phi/charge", sp.factor(kernel.subs(acoustic)))
    low, = next_coefficient(Da, z*q**2, 2)
    gap, = [r for r in sp.solve(Da.subs(q, 0), y) if r != 0]
    high, = next_coefficient(Da, gap+z*q, 1)
    show("unpinned dipole low y/q^2", low)
    show("unpinned dipole upper (intercept, slope in q)", (gap, high))
    show("unpinned dipole low (A/J,u/F) weights / q",
         (weight(Aa, Da, low*q**2, 1, 0), weight(Ua, Da, low*q**2, 1, 0)))
    show("unpinned dipole upper (A/J,u/F) weights at q=0",
         (weight(Aa, Da, gap, 0, 0), weight(Ua, Da, gap, 0, 0)))

    # a and b are squared bare EM and shear speeds, d is f^2/(eps*rho).
    a, b, d, v = sp.symbols("a b d v", positive=True)
    grad = {s: 0, g0: 0}
    Dg = sp.factor((D.subs(specify).subs(grad)/(eps*rho))
                   .subs({nu: a*eps, mu: b*rho, f**2: d*eps*rho}))
    Ag = sp.factor((nA/rho).subs(specify).subs(grad).subs(mu, b*rho))
    Ug = sp.factor((nuu/eps).subs(nu, a*eps))
    show("gradient static phi/charge", sp.factor(kernel.subs(grad)))
    poly_x = sp.factor(Dg.subs(y, x*q)/q**2)
    show("gradient exact x=y/q roots (a=nu/eps,b=mu/rho,d=f^2/(eps*rho))",
         sp.solve(poly_x, x))
    a1, = next_coefficient(Dg, a*q+z*q**2, 3)
    b1, = next_coefficient(Dg, b*q+z*q**2, 3)
    ya, yb = a*q+a1*q**2, b*q+b1*q**2
    show("gradient distinct q^2 coefficient in y, (EM-like,shear-like)", (a1, b1))
    show("gradient distinct (eps*A/J shear weight,rho*u/F EM weight) / q",
         (weight(Ag, Dg, yb, 2, 1), weight(Ug, Dg, ya, 2, 1)))
    show("gradient distinct (eps*A/J EM weight,rho*u/F shear weight) at q=0",
         (weight(Ag, Dg, ya, 1, 1), weight(Ug, Dg, yb, 1, 1)))
    detuning = sp.symbols("detuning", real=True)
    near = {a: v+detuning*k/2, b: v-detuning*k/2, q: k**2}
    split = next_coefficient(Dg.subs(near), v*k**2+z*k**3, 6, k)
    show("gradient near-degenerate (a-b=detuning*k) k^3 coefficients in y", split)
    for splitting in split:
        curve = v*k**2+splitting*k**3
        slope = coefficient(Dg.diff(y).subs(near).subs(y, curve), k, 3)
        normalized_weights = [sp.factor(-coefficient(num.subs(near).subs(y, curve), k, 3)/slope)
                              for num in (Ag, Ug)]
        show("gradient (eps*A/J,rho*u/F) leading weights for k^3 coefficient "
             + str(splitting), normalized_weights)
        show("gradient coincident (k^3 coefficient,eps*A/J weight,rho*u/F weight)",
             [value.subs(detuning, 0) for value in (splitting, *normalized_weights)])

    pin = {f: 0}
    Dp = D.subs(specify).subs(pin)
    Ap = nA.subs(specify).subs(pin)
    Up = nuu.subs(specify).subs(pin)
    yp, = [r for r in sp.solve(Dp.subs(q, 0), y) if r != 0]
    slope_p, = next_coefficient(Dp, yp+z*q, 1)
    show("pinned static phi/charge", sp.factor(kernel.subs(pin)))
    show("pinned lim q*phi/charge", sp.factor((q*kernel).subs(pin).cancel().subs(q, 0)))
    show("pinned (low slope in q,upper intercept,upper slope in q)", (vp, yp, slope_p))
    pinned_flux(pinned, vp)
    show("pinned low (A/J weight, u/F weight / q)",
         (weight(Ap, Dp, vp*q, 0, 0), weight(Up, Dp, vp*q, 1, 0)))
    show("pinned high (A/J,u/F) weights at q=0",
         (weight(Ap, Dp, yp, 0, 0), weight(Up, Dp, yp, 0, 0)))
    show("pinned longitudinal y", sp.factor(long_gap.subs(specify).subs(pin)))
    show("pinned longitudinal phi/charge weight times q",
         sp.factor(-sp.residue(q*responseL[0].diff(charge).subs(k**2, q).subs(specify).subs(pin),
                               y, long_gap.subs(specify).subs(pin))))
    show("decoupled (A/J,u/F,u/J), any T", [sp.factor(entry.subs(g, 0)/D.subs(g, 0))
                                         for entry in (nA, nuu, adj[1, 0])])
    show("decoupled static phi/charge", kernel.subs({g0: 0, f: 0}).cancel())

    # Replace the pinned matter coordinate by internal p and add neutral U.
    # The ordering is A,U,p. The supplied matter block is mass-orthogonal.
    mc, mp, muc, mup = sp.symbols("mC mP muC muP", positive=True)
    h = sp.symbols("h", real=True)
    two = sp.zeros(3)
    two[0, 0] = M[0, 0]
    two[0, 2] = M[0, 1].subs(g, g0)
    two[2, 0] = M[1, 0].subs(g, g0)
    two[2, 2] = M[1, 1].subs({rho: mp, T: s+mup*q})
    two[1, 1] = muc*q-mc*w**2
    two[1, 2] = two[2, 1] = h*q
    D2 = even_frequency(two.det())
    A2 = even_frequency(two.minor_submatrix(0, 0).det())
    leading = coefficient(D2.subs(y, x*q), q, 2)
    show("mobile reference leading D/q^2 at y=x*q", leading)
    show("mobile reference q=0 determinant", sp.factor(D2.subs(q, 0)))
    show("mobile reference h=0 A/J after cancellation", sp.factor((A2/D2).subs(h, 0)))
    common = sp.symbols("cC2", positive=True)
    mobile = {muc: mc*common}
    correction, = next_coefficient(D2.subs(mobile), common*q+z*q**2, 3)
    show("mobile reference distinct common branch q^2 shift, cC2=muC/mC", correction)
    show("mobile reference distinct common branch A/J weight / q",
         weight(A2.subs(mobile), D2.subs(mobile), common*q+correction*q**2, 2, 1))
    crossing = {muc: mc*vp, q: k**2}
    D2_cross = D2.subs(crossing)
    mixed = next_coefficient(D2_cross, vp*k**2+z*k**3, 6, k)
    show("mobile reference coincident h!=0 k^3 coefficients", mixed)
    for splitting in mixed:
        curve = vp*k**2+splitting*k**3
        num = coefficient(A2.subs(crossing).subs(y, curve), k, 3)
        den = coefficient(D2_cross.diff(y).subs(y, curve), k, 3)
        show("mobile reference coincident A/J weight for " + str(splitting), sp.factor(-num/den))

    # Compare independent representations only after emitting the responses.
    compare("operator adjugate check (0,0)", even_frequency((M*M.adjugate())[0, 0]), D)
    compare("pinned optical gap from longitudinal and transverse systems",
            long_gap.subs(specify).subs(pin).subs(q, 0), yp)
    compare("pinned static and propagating susceptibilities", sp.factor(vp/nu),
            (q*kernel).subs(pin).cancel().subs(q, 0))
    compare("two-component leading product against pinned block",
            leading, (muc-mc*x)*coefficient(Dp.subs(y, x*q), q, 1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pinned-flux", action="store_true",
                        help="only measure the pinned displacement and shear/EM flux ratio")
    main(pinned_flux_only=parser.parse_args().pinned_flux)
    if failures:
        raise AssertionError(", ".join(failures))
