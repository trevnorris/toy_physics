"""Vary the conservative couplings and test the actual normalized intake on the resting mouth."""

import numpy as np
import sympy as sp

from held_crossflow_budget import archived_packet


def confinement_and_inertia():
    q, radius_offset, polar, eps, n0, n1 = sp.symbols("q d mu eps n0 n1", real=True)
    v = sp.Function("V")
    energy = (n0 + eps*n1*polar)*v(radius_offset-q*polar)
    mixed = sp.diff(energy, eps, q).subs(q, 0)
    wall_drive = sp.simplify(-2*sp.pi*sp.integrate(mixed, (polar, -1, 1)))
    print("wall force per r^2 dr dw and unit dipolar density drive =", wall_drive)

    # Unit-kinetic-normalized conservative mode, with zero direct X kinetic
    # term. Test the consequence *if* this reduction exists for translation.
    t, s, stiffness, coupling, frequency = sp.symbols("t s K c Om", positive=True)
    x, y = (sp.Function(name)(t) for name in ("X", "Y"))
    lagrangian = (sp.diff(y, t)**2 - frequency**2*y**2 - stiffness*x**2)/2 + coupling*x*y
    def required_force(coord):
        return sp.diff(sp.diff(lagrangian, coord.diff(t)), t) - sp.diff(lagrangian, coord)
    mode = required_force(y).subs(y.diff(t, 2), s**2*y)
    y_driven, = sp.solve(mode, y)
    response = sp.factor(required_force(x).subs(y, y_driven)/x)
    slow_mass = sp.simplify(response.diff(s, 2).subs(s, 0)/2)
    fast = sp.limit(response, s, sp.oo)
    print("required force per displacement in Laplace space =", response)
    print("slow induced coefficient, fast force/displacement limit =", slow_mass, fast)
    if (sp.simplify(wall_drive-4*sp.pi*n1*v(radius_offset).diff(radius_offset)/3) != 0
            or slow_mass != coupling**2/frequency**4 or fast != stiffness):
        raise ArithmeticError("Confinement variation or field elimination failed")


def intake_variations():
    # Moments of the baseline current against h^k on the 3D opening.
    eps, j, m1, m2, m3 = sp.symbols("eps J m1 m2 m3", real=True, nonzero=True)
    paired = (m1+2*eps*m2+eps**2*m3)/(1+2*eps*m1/j+eps**2*m2/j)
    curl = sp.factor(-paired.diff(eps).subs(eps, 0))
    print("amplitude-then-phase mixed work minus reverse work =", curl, flush=True)
    if sp.simplify(curl + 2*(m2-m1**2/j)) != 0:
        raise ArithmeticError("Normalized intake derivative failed")

    path = archived_packet("mouth_feed_weak_run.json", "n32_restoring_slope_u0.npz")
    with np.load(path, allow_pickle=False) as data:
        u = abs(data["field"][-1])
        current = data["current"]
        dx = float(data["dx"])
        if float(data["speed"]) != 0 or u.shape != (32, 32, 32) or np.any(u <= 0):
            raise ValueError("Expected the positive-density, resting periodic preparation")
    axis = (np.arange(32)+0.5-16)*dx
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij", sparse=True)
    weight = np.maximum(1-x*x-y*y-z*z, 0)**3
    h = x*weight
    vol = dx**3
    rate = vol * np.sum(current)
    np.testing.assert_allclose(current, 0.25*weight/(vol*weight.sum()), atol=1e-14, rtol=1e-12)
    normal_shape = weight/u**2
    def trial_current(test, delta):
        raw = normal_shape*(u*(1+delta*test))**2
        return rate*raw/(vol*raw.sum())
    def mixed(test, delta):
        dj = (trial_current(test, delta) - trial_current(test, -delta))/(2*delta)
        return -vol*np.sum(test*dj), vol*dj.sum()
    variance_result = -2*vol*np.sum(current*(h-vol*np.sum(h*current)/rate)**2)
    print(f"rate={rate:.16g} mixed_variance={variance_result:.16g}", flush=True)
    for delta in (1e-4, 5e-5):
        result, rate_change = mixed(h, delta)
        uniform, _ = mixed(np.ones_like(h), delta)
        print(f"eps={delta:g} mixed_fd={result:.16g} rate_derivative={rate_change:.3g} "
              f"uniform_control={uniform:.3g}", flush=True)
        if (variance_result >= 0 or not np.isclose(result, variance_result, rtol=2e-8, atol=0)
                or max(abs(rate_change), abs(uniform)) > 1e-10):
            raise ArithmeticError("Saved-state variation or fixed-rate control failed")


if __name__ == "__main__":
    confinement_and_inertia()
    intake_variations()
