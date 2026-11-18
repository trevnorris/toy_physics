"""Rim work from the structural dipole equation, with densitized coefficients."""

import sympy as sp


def check(name, left, right):
    residual = sp.simplify(left - right)
    print(f"{name}: {left} | {right} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(name)


def work_and_reciprocity():
    w, t, u, lam = sp.symbols("w t u lam", real=True)
    angular = 2 * sp.pi * sp.integrate(u**2, (u, -1, 1))
    check("angular norm", angular, 4 * sp.pi / 3)
    tension, mu, stiffness = [sp.Function(key)(w) for key in ("T", "mu", "K1")]
    q = sp.Function("q")(w, t)
    vel, slope = q.diff(t), q.diff(w)
    energy = angular * (mu * vel**2 + tension * slope**2 + stiffness * q**2) / 2
    source = mu * vel.diff(t) - (tension * slope).diff(w) + stiffness * q
    power = sp.expand(energy.diff(t) - angular * source * vel)
    check("energy rate less source work", power, (angular * tension * slope * vel).diff(w))

    def operator(profile, parameter):
        return -(tension * profile.diff(w)).diff(w) + (stiffness + parameter * mu) * profile

    p, r = [sp.Function(key)(w) for key in ("p", "r")]
    concomitant = tension * (r * p.diff(w) - p * r.diff(w))
    check("same-frequency source identity", sp.expand(p * operator(r, lam) - r * operator(p, lam)),
          concomitant.diff(w))
    check("two-frequency overlap identity", sp.expand(p * operator(r, lam) - r * operator(p, 0)),
          concomitant.diff(w) + lam * mu * p * r)


def free_end_control():
    w, s, length, tension, mu = sp.symbols("w s L T mu", positive=True)
    k = s * sp.sqrt(mu / tension)
    profile = sp.cosh(k * (length - w)) / sp.cosh(k * length)
    check("unit mouth trace", profile.subs(w, 0), 1)
    check("free remote traction", tension * profile.diff(w).subs(w, length), 0)
    check("elastic equation", tension * profile.diff(w, 2), mu * s**2 * profile)
    angular = 4 * sp.pi / 3
    response = sp.simplify(-angular * tension * profile.diff(w).subs(w, 0))
    overlap = sp.integrate(profile, (w, 0, length))
    check("drive force per displacement", response, sp.simplify(angular * s**2 * mu * overlap))
    check("slow force per acceleration", sp.limit(response / s**2, s, 0), angular * mu * length)
    check("fast force per acceleration", sp.limit(response / s**2, s, sp.oo), 0)


if __name__ == "__main__":
    work_and_reciprocity()
    free_end_control()
