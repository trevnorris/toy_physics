"""Interaction energy and moving-cut power of two prescribed ideal monopoles.

The pair and both isolated controls use the same domain, excluding two
disjoint translating spheres. Rates and radii are fixed; no translating
dipole or physical mouth-response condition is supplied. Terms linear in
the prescribed center velocities measure the work of the static force.
"""
import sympy as sp


def emit(label, expr):
    print(f"{label}: {sp.factor(expr)}", flush=True)


def check(label, lhs, rhs):
    residual = sp.factor(lhs-rhs)
    print(f"{label}: {sp.factor(lhs)} | {sp.factor(rhs)} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(label)


def channel(n):
    d, rho, r = sp.symbols("d rho r", positive=True)
    q1, q2, v1, v2 = sp.symbols("q1 q2 V1 V2", real=True)
    omega = 2*sp.pi**sp.Rational(n, 2)/sp.gamma(sp.Rational(n, 2))
    # Green identity: the companion potential is harmonic on the other ball.
    # On hole 1, phi1 is constant and the companion normal flux is zero.
    # On hole 2, the mean phi1 equals its value at center 2, and the domain
    # normal points into the hole. The boundary at infinity vanishes.
    phi1 = -sp.integrate(q1/(rho*omega*r**(n-1)), (r, d, sp.oo))
    energy = sp.factor(-q2*phi1)
    energy_rate = sp.diff(energy, d)*(v2-v1)
    u21 = -q2/(rho*omega*d**(n-1))
    u12 = q1/(rho*omega*d**(n-1))

    # Cross part of the momentum tensor on sphere 1. g is its own radial
    # speed, c_x and c_r are the full companion components on the surface.
    # The cancellation is pointwise, before the harmonic mean is used.
    g, cx, cr, mu, va, vb, area = sp.symbols("g c_x c_r mu V_i V_j area", real=True)
    cross_momentum = rho*(g*mu*cr+cx*g-g*cr*mu)
    check("cross momentum integrand", cross_momentum, rho*g*cx)
    # From unsteady Bernoulli, dot(phi_i)=-V_i . grad(phi_i).
    # J_cross=-rho(dot(phi1)*v2+dot(phi2)*v1). Reynolds subtracts E_cross*V_i.
    lab_flux = rho*(va*g*mu*cr+vb*cx*g)
    moving_term = -rho*g*cr*va*mu
    cross_flux = sp.factor(lab_flux+moving_term)
    check("moving-cut energy integrand", cross_flux, rho*g*vb*cx)
    # A regular harmonic velocity has its center value as spherical average.
    flux1 = (area*cross_flux).subs({g:q1/(rho*area), cx:u21, vb:v2})
    flux2 = (area*cross_flux).subs({g:q2/(rho*area), cx:u12, vb:v1})
    force1 = -(area*cross_momentum).subs({g:q1/(rho*area), cx:u21})
    force2 = -(area*cross_momentum).subs({g:q2/(rho*area), cx:u12})
    mechanical = sp.factor(force1*v1+force2*v2)
    supply1 = sp.factor(flux1+force1*v1)
    supply2 = sp.factor(flux2+force2*v2)
    print(f"dimension = {n}", flush=True)
    emit("common-domain interaction energy", energy)
    emit("interaction energy rate", energy_rate)
    emit("force on center 1", force1)
    emit("force on center 2", force2)
    emit("mechanical work rate", mechanical)
    emit("energy from mouth 1 into exterior", flux1)
    emit("energy from mouth 2 into exterior", flux2)
    check("moving exterior energy balance", energy_rate, flux1+flux2)
    emit("inner energy required at mouth 1", supply1)
    emit("inner energy required at mouth 2", supply2)
    check("whole interaction energy requirement", supply1+supply2, energy_rate+mechanical)

    # The l=1 normal projection of a regular harmonic companion has angular
    # mean <mu*c_r>=u_center/n. Higher harmonics integrate to zero.
    th = sp.Symbol("theta", real=True)
    norm = sp.integrate(sp.sin(th)**(n-2), (th,0,sp.pi))
    projection = sp.integrate(sp.cos(th)**2*sp.sin(th)**(n-2), (th,0,sp.pi))/norm
    wrong_no_unsteady = -projection*(q1*v1*u21+q2*v2*u12)
    wrong_no_reynolds = flux1+flux2+projection*(q1*v1*u21+q2*v2*u12)
    for label, wrong in [("omit unsteady pressure", wrong_no_unsteady),
                         ("omit moving-cut transport", wrong_no_reynolds)]:
        defect = sp.factor(energy_rate-wrong)
        emit(label+", energy defect", defect)
        if defect == 0:
            raise ArithmeticError(label)
    # Common translation changes local power allocation, not separation energy.
    check("common translation interaction storage", energy_rate.subs(v2,v1),0)
    check("one drain absent", energy.subs(q2,0),0)
    return energy, mechanical, supply1+supply2, (q1,q2,rho,d,v1,v2)


def join(brane, bulk):
    d1, d2, rho3, rho4, sep = sp.symbols("D1 D2 rho3 rho4 d", positive=True)
    v1, v2 = sp.symbols("V1 V2", real=True)
    terms=[]
    for data, sign, density in [(brane,-1,rho3),(bulk,1,rho4)]:
        energy, work, supplied, (q1,q2,rho,d,va,vb)=data
        replace={q1:sign*d1,q2:sign*d2,rho:density,d:sep,va:v1,vb:v2}
        terms.append(tuple(sp.factor(x.subs(replace)) for x in (energy,work,supplied)))
    for label,index in [("brane-plus-bulk interaction energy",0),
                        ("brane-plus-bulk mechanical power",1),
                        ("brane-plus-bulk inner energy requirement",2)]:
        emit(label, sum(row[index] for row in terms))


if __name__ == "__main__":
    print("SymPy version =", sp.__version__, flush=True)
    join(channel(3),channel(4))
