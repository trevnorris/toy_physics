"""Carrier trace, conjugation and weak electric readout at a stationary rim."""

import sympy as sp


def record(label, result, expected):
    value = sp.factor(result)
    print(f"{label} = {value}", flush=True)
    if sp.simplify(value - expected) != 0:
        raise ArithmeticError(label)


def carrier_trace():
    x, y, u, v, h = sp.symbols("x y u v h", real=True)
    p = sp.Matrix([x, y])
    d = sp.Matrix([u, v])  # normal covariant derivative in units with A_i=0
    form = sp.Matrix([[0, 1], [-1, 0]])

    def flux(trace, derivative):
        return (trace.T * form * derivative)[0]

    j = flux(p, d)
    # h is real. C [D_N psi] = h psi, with the same inward N on both traces.
    c = sp.symbols("C", positive=True)
    record("current jump in single-trace contact (hbar/m removed)",
           flux(p, d + h*p/c) - j, 0)
    a, b, f, g = sp.symbols("a b f g", real=True)
    internal = sp.Matrix([[a, b], [f, g]])
    record("normal current under constant internal map, relative to J",
           flux(internal*p, internal*d)/j, internal.det())
    for s in (1, -1):
        # The lab derivative AND the inward normal transform under reflection.
        for label, op in [("pullback", sp.eye(2)), ("conjugation", sp.diag(1, -1))]:
            inward = s * flux(op*p, s*op*d)
            record(f"mouth {s:+} {label} intake / original", inward/j, op.det())

    # Write the variation of h(n,chi) psi in independent conjugate coordinates.
    z, zb, dz, dzb, dchi, hn, hc, t = sp.symbols("z zb dz dzb dchi h_n h_chi t")
    density = (z+t*dz)*(zb+t*dzb)
    force = (h + hn*(density-z*zb) + hc*t*dchi)*(z+t*dz)
    tangent = sp.expand(force.diff(t).subs(t, 0))
    record("contact coefficients, direct", tangent.coeff(dz), h+z*zb*hn)
    record("contact coefficients, conjugate", tangent.coeff(dzb), hn*z**2)
    record("contact coefficients, order", tangent.coeff(dchi), hc*z)
    record("contact response to a unit phase variation",
           tangent.subs({dz:sp.I*z, dzb:-sp.I*zb, dchi:0}), sp.I*h*z)


def probe():
    n, chi, radius, q, phi, phase_dot, hb = sp.symbols("n chi R q phi theta_dot hbar", real=True)
    H, stiffness, drive, s = sp.symbols("H K J s", real=True)
    z, t, alpha = sp.symbols("z t alpha", real=True)
    speed, length, amp, inertia = sp.symbols("c L amp rho_w", positive=True)
    k = sp.pi/(2*length)
    wave = amp*sp.sin(k*z)*sp.cos(speed*k*t+alpha)
    for side in (1, -1):
        lab_wave = wave.subs(z, side*z)
        record(f"mouth {side:+} support node", lab_wave.subs(z, 0), 0)
        record(f"mouth {side:+} free-end slope", lab_wave.diff(z).subs(z, side*length), 0)
    # Integrate in inward arclength for the upper representative. The same
    # squares and positive volume element result for the mirror.
    wave_energy = sp.integrate(inertia*(wave.diff(t)**2 + speed**2*wave.diff(z)**2)/2,
                               (z, 0, length))
    record("support energy at arbitrary temporal phase", sp.trigsimp(wave_energy),
           inertia*amp**2*speed**2*sp.pi**2/(16*length))
    # Per carrier volume, -n(hbar theta_dot + q phi) is the temporal action.
    # Keep surface and stored-mode contributions separate measures.
    carrier = -n*(hb*phase_dot + q*phi)
    rim = stiffness*H**2/2 - drive*s*H
    contact = sp.Function("g")(n, chi, radius)
    record("carrier electric source per volume", -carrier.diff(phi), q*n)
    record("scalar contact direct electric source per area", contact.diff(phi), 0)
    record("annular rim direct electric source per normalized weight", rim.diff(phi), 0)
    record("wave direct electric source", wave_energy.diff(phi), 0)
    record("rim load conjugate to H at H=0", -rim.diff(H).subs(H, 0), drive*s)
    rate = sp.symbols("Gamma", real=True)
    record("charge source from order exchange, both fractions", q*n*rate + q*(-n*rate), 0)

    n1, n2, t1, t2, v1, v2 = sp.symbols("n1 n2 t1 t2 V1 V2", real=True)
    w1, w2, gradient = sp.symbols("w1 w2 gradient", real=True)
    def excess(side):
        return sum(vol*q*(density-ref)*(1+gradient*side*depth)
                   for vol, density, ref, depth in [(v1,n1,t1,w1), (v2,n2,t2,w2)])

    up, down = excess(1), excess(-1)
    print("each mirror residual for uniform probe =", sp.expand(up.subs(gradient, 0)), flush=True)
    contrast = sp.expand((up-down)/2)
    record("uniform-probe half difference", contrast.subs(gradient, 0), 0)
    record("depth-gradient half difference per gradient", contrast.diff(gradient),
           q*(v1*w1*(n1-t1)+v2*w2*(n2-t2)))


if __name__ == "__main__":
    carrier_trace()
    probe()
