"""Local, dimensionless 0PN example with an open receiving half-space.

Use the smooth carrier benchmark and a single illustrative order kinetics in
both preparations. Order is passive in this example. The full throat coupling
and upstream match remain unspecified.
"""
import sympy as sp

from sleeve_outlet_boundary_data import beta, c2, k, lam, linear_equations, v, w, x


def main():
    # Use one of the previous benchmark's decaying modes, without its
    # prescribed zero normal-current perturbation at the cut.
    eps, U, alpha, offset = sp.symbols("eps U alpha b", real=True)
    params = {v: 1, c2: sp.Rational(35, 24), beta: sp.Rational(1, 24)}
    mode = sp.sin(k*x)*sp.exp(-lam*w)
    potential_amplitude = sp.Symbol("P", real=True)
    eqs = linear_equations(mode, potential_amplitude*mode)
    reduced = [sp.simplify((eq/mode).subs(params).subs({k: 1, lam: 2})) for eq in eqs]
    coefficient, = sp.solve(reduced[0], potential_amplitude)
    print("mode [continuity, phase equation] operands =", reduced)
    print("potential amplitude for unit density =", coefficient)
    density_change = alpha*U*mode.subs({k: 1, lam: 2})
    potential_change = U*x + coefficient*density_change
    bulk = linear_equations(density_change, potential_change).subs(params).applyfunc(sp.simplify)
    print("driven bulk [continuity, phase equation] residuals =", list(bulk))

    n = 1 + eps*density_change
    ux = eps*sp.diff(potential_change, x)
    uw = 1 + eps*sp.diff(potential_change, w)
    current = sp.diff(n*uw, eps).subs(eps, 0).expand()
    print("at mouth [density, axial velocity, normal current] derivatives =",
          [sp.simplify(sp.diff(field, eps).subs({eps: 0, w: 0})) for field in (n, uw, n*uw)])
    print("deep [density, lateral velocity, axial velocity] derivatives =",
          [sp.limit(sp.diff(field, eps), w, sp.oo) for field in (n, ux, uw)])
    print("current derivative integrated on symmetric strip =",
          sp.integrate(current.subs(w, 0), (x, -sp.pi, sp.pi)))

    # Net de-ordering per constituent chi*(1-chi) gives material derivative
    # D_t Y=1 for Y=log[(1-chi)/chi]. Here the steady local profiles realize
    # that dynamical law, with a common decaying correction at the far end.
    y1 = sp.integrate(-sp.diff(potential_change, w), w)
    y1 -= sp.limit(y1, w, sp.oo)
    Y = w-offset + eps*y1
    advected = ux*sp.diff(Y, x) + uw*sp.diff(Y, w)
    print("order coordinate perturbation =", y1)
    print("advected order coordinate [base, derivative] =",
          [sp.simplify(advected.subs(eps, 0)), sp.simplify(sp.diff(advected, eps).subs(eps, 0))])

    # The level's position follows Y=0. The excess net loss per volume
    # peaks where the normal derivative of n*chi*(1-chi) vanishes.
    chi = 1/(1 + sp.exp(Y))
    loss = n*chi*(1-chi)
    at_peak = {eps: 0, w: offset}
    def at_base(expression):
        return sp.simplify(expression.subs(at_peak))
    level_shift = at_base(-sp.diff(Y, eps)/sp.diff(Y, w))
    peak_slope = sp.diff(loss, w)
    peak_shift = at_base(-sp.diff(peak_slope, eps)/sp.diff(peak_slope, w))
    print("at isolated w=b [fraction, loss, loss slope, loss curvature] =",
          [at_base(field) for field in (chi, loss, peak_slope, sp.diff(peak_slope, w))])
    print("[half-order shift, net-loss peak shift] =", [level_shift, peak_shift])
    if any(bulk) or at_base(sp.diff(advected, eps)) != 0:
        raise ValueError("Nonzero local equation residual")


if __name__ == "__main__":
    main()
