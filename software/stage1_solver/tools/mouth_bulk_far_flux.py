"""First-order axial momentum on a 4D sphere in the steady neutral carrier."""

import sympy as sp


def check(label, actual, expected):
    error = sp.factor(actual-expected)
    print(f"{label}: {sp.factor(actual)} | {sp.factor(expected)} | residual {error}", flush=True)
    if error != 0:
        raise ArithmeticError(label)


def main():
    r, m, k, n, s = sp.symbols("r m k n s", positive=True)
    u, t, mu, e = sp.symbols("u t mu e", real=True)
    logn, b = (sp.Function(name)(r) for name in ("L", "b"))
    area = 2*sp.pi**2*r**3
    mass_rate = area*m*n*s

    # mu is the axial direction cosine on S^3. The first-order flow is
    # radial u*mu plus tangent t*(e_x-mu*e_r). s is the base radial speed.
    def average(value):
        return sp.factor(2/sp.pi*sp.integrate(
            sp.expand(value)*sp.sqrt(1-mu**2), (mu, -1, 1)))

    density = n*(1+e*b*mu)
    radial = s+e*u*mu
    axial = s*mu+e*(u*mu**2+t*(1-mu**2))
    carried = average(sp.diff(m*density*radial*axial, e).subs(e, 0))

    # Q = -k [Laplacian(log n) + |grad log n|^2/2], k=hbar^2/(4m).
    # On S^3 the l=1 angular Laplacian of mu is -3 mu.
    delta_q = -k*(b.diff(r, 2)+3*b.diff(r)/r-3*b/r**2+logn.diff(r)*b.diff(r))*mu
    # Steady phase equation, angle-independent chemical potential, dp=n dh.
    pressure = average(n*(-m*s*u*mu-delta_q)*mu)
    # Hessian of the logarithmic density in the axial, radial directions.
    hx_r = b.diff(r, 2)*mu**2+(b.diff(r)/r-b/r**2)*(1-mu**2)
    quantum = average(-k*n*(hx_r + b*logn.diff(r, 2)*mu**2))
    for label, term in (("carriage", carried), ("pressure", pressure), ("density stress", quantum)):
        print(f"{label} per area = {term}", flush=True)
    result = sp.factor(area*(carried+pressure+quantum))
    material = mass_rate*(u+3*t+s*b)/4
    gradient = area*k*n*(logn.diff(r)*b.diff(r)-logn.diff(r, 2)*b)/4
    check("integrated export", result, material+gradient)
    check("constant base density, density-gradient contribution",
          (result-material).subs(logn, 0).doit(), 0)
    # Algebraic controls; not global solutions for the radial background.
    check("unit stream at the cut", result.subs({b:0, u:1, t:1}).doit(), mass_rate)
    check("decaying potential dipole at the cut", result.subs({b:0, u:-3, t:1}).doit(), 0)


if __name__ == "__main__":
    main()
