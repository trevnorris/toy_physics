"""Constituent and order budgets per unit three-area in the existing flat example.

The isolated component of conversion_region_weak_flow has n=m=u_w=1,
no lateral velocity and no relative order current. The cuts are fixed.
"""
import sympy as sp


def main():
    w, t, offset = sp.symbols("w t b", real=True)
    length = sp.Symbol("L", positive=True)
    n = sp.S.One
    # In this test the velocity potential is in units with m=1.
    potential = w - (sp.Rational(1, 2) + sp.Rational(35, 96))*t
    uw = sp.diff(potential, w)
    current = n*uw
    chi = 1/(1 + sp.exp(w-offset))
    ordered = chi*current
    bulk = (1-chi)*current
    loss_density = n*chi*(1-chi)
    print("constituent [intake, discharge, storage, side current] =",
          [current.subs(w, 0), current.subs(w, length), sp.diff(sp.integrate(n, (w, 0, length)), t), 0])
    print("local [ordered divergence, net loss density] =",
          [sp.diff(ordered, w), loss_density])
    order_residual = sp.simplify(sp.diff(ordered, w) + loss_density)
    # Integrate the local loss in z=exp(w-b) between the two cuts.
    z = sp.Symbol("z", positive=True)
    loss_per_z = sp.cancel(loss_density.subs(sp.exp(w-offset), z)/z)
    loss = sp.integrate(loss_per_z, (z, sp.exp(-offset), sp.exp(length-offset)))
    order_drop = ordered.subs(w, 0) - ordered.subs(w, length)
    bulk_gain = bulk.subs(w, length) - bulk.subs(w, 0)
    residuals = [sp.simplify(order_drop-loss), sp.simplify(bulk_gain-loss)]
    print("integrated [order drop, bulk gain, loss] =", [order_drop, bulk_gain, loss])
    print("local and integrated residuals =", [order_residual, *residuals])
    # A finite segment with its conversion midpoint halfway between the cuts.
    sample = {offset: sp.log(3), length: 2*sp.log(3)}
    readings = [ordered.subs(w, 0), bulk.subs(w, 0), ordered.subs(w, length),
                bulk.subs(w, length), loss]
    print("sample [in ordered, in disordered, out ordered, out disordered, net converted] =",
          [sp.simplify(val.subs(sample)) for val in readings])
    if order_residual != 0 or any(residuals):
        raise ValueError("Constituent fractions do not satisfy the order budget.")


if __name__ == "__main__":
    main()
