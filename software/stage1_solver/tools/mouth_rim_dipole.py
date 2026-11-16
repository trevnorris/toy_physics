"""Project the parent shape equation onto the geometric mouth translation.

Use the polar axis along the small applied force. cos(theta) multiplies the
actual center displacement, so no spherical-harmonic normalization is hidden
in the endpoint force. All expressions are linearized about the spherical rim.
"""

import sympy as sp


def main():
    u, eps = sp.symbols("u epsilon", real=True)
    radius, length = sp.symbols("a L", positive=True)
    force, dilation, dipole, deformation = sp.symbols("F b q k", real=True)

    def solid_angle(f):
        return 2*sp.pi*sp.integrate(f, (u, -1, 1))

    def project(f):
        return sp.simplify(3*solid_angle(u*f)/(4*sp.pi))

    quadrupole = (3*u**2-1)/2
    eta = dilation + dipole*u + deformation*quadrupole
    R = radius + eps*eta
    # The volume centroid numerator at the reference sphere vanishes.
    vol = solid_angle(radius**3/3)
    center_change = sp.simplify(solid_angle(sp.diff(u*R**4/4,eps).subs(eps,0))/vol)
    traction = 3*force*u/(4*sp.pi*radius**2)
    work = sp.simplify(solid_angle(traction*radius**2*eta))
    normal_force = sp.simplify(solid_angle(radius**2*traction*u))

    w, t = sp.symbols("w t", real=True)
    mu, Tw, To, K = [sp.Function(n)(w) for n in ("mu","T_w","T_Omega","K")]
    q = sp.Function("q")(w,t)
    displacement = q*u
    lap_sphere = sp.diff((1-u**2)*sp.diff(displacement,u),u)
    projected = project(mu*sp.diff(displacement,t,2)
                        -sp.diff(Tw*sp.diff(displacement,w),w)
                        -To*lap_sphere + K*displacement)
    expected = mu*sp.diff(q,t,2)-sp.diff(Tw*sp.diff(q,w),w)+(K+2*To)*q

    # An algebraic trace test, not two solutions with identical full initial
    # conditions. Both profiles have identical displacement and axial traction
    # at the mouth and the same value and slope at the opposite end.
    b = sp.Symbol("delta_q", real=True)
    delta = b*(w/length)**2*(1-w/length)**2
    traces = [sp.simplify(g.subs(w,at)) for at in (0,length)
              for g in (delta,sp.diff(delta,w))]
    delta_acc = sp.simplify(((sp.diff(Tw*sp.diff(delta,w),w)-(K+2*To)*delta)/mu).subs(w,0))
    print(f"first-order center shift = {center_change}")
    print(f"integrated applied axial force = {normal_force}; first-order work = {work}")
    print(f"projected angular eigenvalue = {project(-lap_sphere)/q}")
    print(f"projected equation residual against dipole operator = {sp.simplify(projected-expected)}")
    print(f"endpoint displacement and slope differences at 0,L = {traces}")
    print(f"endpoint acceleration difference from the axial operator = {delta_acc}")
    assert (center_change==dipole and normal_force==force and work==force*dipole
            and sp.simplify(projected-expected)==0 and traces==[0,0,0,0]
            and sp.simplify(delta_acc-2*b*Tw.subs(w,0)/(length**2*mu.subs(w,0)))==0)


if __name__ == "__main__":
    main()
