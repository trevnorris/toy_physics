"""Translation current on a short wall segment and on a polarized shear patch."""

import sympy as sp


def check(label, value, reference):
    difference = sp.simplify(sp.trigsimp(value - reference))
    print(f"{label}: computed={value}; comparison={reference}; difference={difference}", flush=True)
    if difference != 0:
        raise ArithmeticError(label)


def wall():
    w, b, mu, tension, s = sp.symbols("w b mu T s", positive=True)
    x, y = sp.symbols("X Y", real=True)
    k = s * sp.sqrt(mu/tension)
    q = (x*sp.sinh(k*(b-w)) + y*sp.sinh(k*w))/sp.sinh(k*b)
    current = -tension*sp.diff(q,w)
    inventory = sp.integrate(mu*s*q, (w,0,b))
    check("wall equation", tension*q.diff(w,2), mu*s**2*q)
    check("wall incoming minus outgoing", current.subs(w,0)-current.subs(w,b), s*inventory)
    check("wall cut current", current.subs(w,b), tension*k*(x-y*sp.cosh(k*b))/sp.sinh(k*b))

    # Characterize the termination only by the force it requires per displacement.
    for label, load in [("free", 0), ("reflecting extension", tension*k*sp.tanh(k*b))]:
        end = sp.solve(current.subs(w,b)-load*y, y)[0]
        flux_per_a = sp.simplify(load*end/(s**2*x))
        store_per_a = s*inventory.subs(y,end)/(s**2*x)
        check(label+" slow stored rate per acceleration", sp.limit(store_per_a,s,0), mu*b)
        check(label+" slow cut flux per acceleration", sp.limit(flux_per_a,s,0),
              0 if load == 0 else mu*b)


def wave():
    # One intrinsic transverse displacement. Rotate the transverse integral to
    # its even profile about the packet center so integral h*h_x vanishes.
    rho, modulus, f, fz, h, hx, velocity, slope, omega = sp.symbols(
        "rho B f fz h hx V slope omega", real=True)
    phase = sp.symbols("phase", real=True)
    co, si = sp.cos(phase), sp.sin(phase)
    ut = -f*(velocity*hx*co + omega*h*si)
    ux = f*hx*co
    uz = (fz*h - f*hx*slope)*co
    kinetic = rho*ut**2/2
    strain = modulus*(ux**2+uz**2)/2
    p = sp.diff(kinetic,velocity)
    flow = -sp.diff(strain,slope)
    def average(expr):
        return sp.expand(sp.integrate(expr, (phase,0,2*sp.pi))/(2*sp.pi))
    check("cycle-averaged packet momentum density", average(p), rho*f**2*hx**2*velocity/2)
    check("cycle-averaged flux before transverse integration", average(flow),
          modulus*f*fz*h*hx/2 - modulus*f**2*hx**2*slope/2)
    z, length = sp.symbols("z L", positive=True)
    for n in (0,1):
        mode = sp.sqrt(2/length)*sp.sin((2*n+1)*sp.pi*z/(2*length))
        weight = sp.integrate(mode**2, (z,0,length/4))/sp.integrate(mode**2, (z,0,length))
        check(f"first-quarter translation fraction, mode {n}", weight,
              sp.Rational(1,4)-sp.sin((2*n+1)*sp.pi/4)/((2*n+1)*sp.pi))
        print(f"mode {n} fraction = {float(weight):.12f}", flush=True)


if __name__ == "__main__":
    wall()
    wave()
