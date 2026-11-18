"""Energy on translating cuts of the ideal 3D intake and 4D receiver.

Constant density, inviscid irrotational exterior, fixed throughput and far
lab stream. Retain quadratic weak lateral motion. B(t) is an unspecified
decaying dipole; no boundary condition or inertia is assigned to the mouth.
"""
import sympy as sp


def emit(label, value):
    print(f"{label}: {sp.factor(value)}", flush=True)


def check(label, lhs, rhs):
    residual = sp.factor(lhs-rhs)
    print(f"{label}: {sp.factor(lhs)} | {sp.factor(rhs)} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(label)


def exterior(n):
    r, a, outer, rho = sp.symbols("r a R rho", positive=True)
    eps, theta = sp.symbols("eps theta", real=True)
    q, u, v, b, bd, pinf, internal = sp.symbols("q U V B Bdot p_inf e", real=True)
    omega = 2*sp.pi**sp.Rational(n, 2)/sp.gamma(sp.Rational(n, 2))
    area = omega*r**(n-1)
    angular_norm = sp.integrate(sp.sin(theta)**(n-2), (theta, 0, sp.pi))

    def average(expr):
        return sp.factor(sp.integrate(sp.expand(expr)*sp.sin(theta)**(n-2),
                                      (theta, 0, sp.pi))/angular_norm)

    def coefficient(expr, order):
        return sp.expand(expr).coeff(eps, order)

    mu, sn = sp.cos(theta), sp.sin(theta)
    g = q/(rho*area)
    # Lab velocity, and relative crossing speed for a cut centered on X(t).
    vr = g+eps*mu*(u-(n-1)*b/r**n)
    vt = -eps*sn*(u+b/r**n)
    vx = mu*vr-sn*vt
    crossing = vr-eps*v*mu
    relative2 = crossing**2+(vt+eps*v*sn)**2
    speed2 = vr**2+vt**2
    p_steady = pinf+rho*(eps**2*(u-v)**2-relative2)/2
    p_unsteady = -eps*rho*bd*mu/r**(n-1)
    pressure = p_steady+p_unsteady
    transport = rho*(internal+speed2/2)*crossing
    work = pressure*vr  # lab material velocity, not relative crossing speed
    energy = area*average(coefficient(transport+work, 2))
    energy0 = area*average(coefficient(transport+work, 0))
    momentum = area*average(coefficient(rho*vx*crossing+pressure*mu, 1))
    print(f"dimension = {n}", flush=True)
    emit("mass flux", area*rho*average(crossing))
    emit("energy base", energy0)
    emit("quadratic transported energy", area*average(coefficient(transport, 2)))
    emit("quadratic pressure work", area*average(coefficient(work, 2)))
    emit("quadratic total energy flux", energy)
    emit("linear momentum flux", momentum)

    # Stored kinetic energy relative to the radial flow and the same-volume
    # uniform lab stream, on a finite spherical shell with fixed radii.
    localized_density = average(coefficient(rho*speed2/2, 2))-rho*u**2/2
    stored = sp.integrate(area*localized_density, (r, a, outer))
    stored_rate = sp.diff(stored, b)*bd
    emit("localized kinetic energy in shell", stored)
    emit("localized kinetic energy rate", stored_rate)
    check("shell energy conservation", energy.subs(r, outer)-energy.subs(r, a), -stored_rate)

    # FORM control: omit the time derivative from Bernoulli pressure while
    # allowing the dipole to evolve; the resulting energy defect must survive.
    wrong_flux = area*average(coefficient(transport+p_steady*vr, 2))
    defect = sp.factor(wrong_flux.subs(r, outer)-wrong_flux.subs(r, a)+stored_rate)
    emit("omit unsteady pressure, shell defect", defect)
    if defect == 0:
        raise ArithmeticError("unsteady-pressure form control was invisible")
    check("steady total energy", energy.subs(bd, 0), q*u**2/2)
    return energy0+energy.subs(bd, 0), momentum.subs(bd, 0), (q, u, pinf, rho, internal)


def joined(feed, bulk):
    rate = sp.Symbol("D", positive=True)
    u, w, v, hf, hb, boost = sp.symbols("U W V H_feed H_bulk boost", real=True)
    def substitute(data, signed, speed, head):
        energy, momentum, (q, local_u, pressure, rho, internal) = data
        values = {q:signed, local_u:speed, internal:head-pressure/rho}
        return sp.factor(energy.subs(values)), sp.factor(momentum.subs(values))
    ef, pf = substitute(feed, -rate, u, hf)
    eb, pb = substitute(bulk, rate, w, hb)
    arrival = sp.factor(-ef-eb)
    force = sp.factor(-pf-pb)
    power = sp.factor(force*v)
    remainder = sp.factor(arrival-power)
    emit("net energy arrival", arrival)
    emit("ideal steady mouth force", force)
    emit("mechanical power delivered to center or holder", power)
    emit("remaining energy rate", remainder)
    check("mouth-frame head difference", remainder,
          rate*(hf-hb+((u-v)**2-(w-v)**2)/2))
    check("remaining energy under observer change",
          remainder.subs({u:u-boost, w:w-boost, v:v-boost}), remainder)
    check("energy arrival under observer change",
          arrival.subs({u:u-boost, w:w-boost}), arrival-boost*force)
    emit("held brane stream, remaining rate", remainder.subs({w:0, v:0}))
    emit("tow through resting media, remaining rate", remainder.subs({u:0, w:0}))
    emit("equal head, W=0, V=U diagnostic", remainder.subs({hb:hf, w:0, v:u}))


if __name__ == "__main__":
    print("SymPy version =", sp.__version__, flush=True)
    joined(exterior(3), exterior(4))
