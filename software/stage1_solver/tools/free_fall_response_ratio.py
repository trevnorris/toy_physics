"""Compare slow free translation for two admissible lateral field normalizations."""

import sympy as sp


def check(label, actual, expected):
    print(f"{label} = {actual}", flush=True)
    if sp.simplify(actual - expected) != 0:
        raise ArithmeticError(label)


def main():
    z, length, rate, speed, density, scale = sp.symbols("z L s c mu lam", positive=True)
    base_inertia, pull = sp.symbols("M C", positive=True)
    # q(z,s) per unit imposed rim displacement. The continuation reflects only
    # the sideways perturbation. The separate material drain remains open.
    shape = sp.cosh(rate * (length - z) / speed) / sp.cosh(rate * length / speed)
    line_mass = scale * density
    tension = line_mass * speed**2
    check("field equation residual", sp.simplify(tension * shape.diff(z, 2)
          - line_mass * rate**2 * shape), 0)
    check("unit rim motion", shape.subs(z, 0), 1)
    check("far lateral traction", tension * shape.diff(z).subs(z, length), 0)

    rim_force = sp.simplify(-tension * shape.diff(z).subs(z, 0))
    momentum = sp.integrate(line_mass * rate * shape, (z, 0, length))
    check("rim force minus momentum rate per unit displacement",
          sp.simplify(rim_force - rate * momentum), 0)
    check("static stiffness", sp.limit(rim_force, rate, 0), 0)
    line_inertia = sp.limit(rim_force / rate**2, rate, 0)
    check("slow lateral field inertia", line_inertia, scale * density * length)

    # No direct stream drive on this incremental lateral sector. The
    # common exterior delivers pull * U; the holder balances it at X=0.
    accel_per_flow = sp.factor(sp.limit(
        pull * rate**2 / (base_inertia * rate**2 + rim_force), rate, 0))
    for normalization in (1, 2):
        check(f"slow acceleration per speed at scale {normalization}",
              accel_per_flow.subs(scale, normalization),
              pull / (base_inertia + normalization * density * length))
    change = sp.factor(accel_per_flow.subs(scale, 2) - accel_per_flow.subs(scale, 1))
    print(f"change at fixed held pull = {change}", flush=True)
    if change.is_negative is not True:
        raise ArithmeticError("expected smaller slow acceleration for positive pull")

    # Repeat for a sector that also receives a separately specified steady
    # force lambda * B * U. This is an affine family of effective laws.
    added_pull, mode_inertia = sp.symbols("B m", real=True)
    ratio = (pull + scale * added_pull) / (base_inertia + scale * mode_inertia)
    check("scale sensitivity times squared total inertia",
          sp.factor(sp.diff(ratio, scale) * (base_inertia + scale * mode_inertia)**2),
          added_pull * base_inertia - pull * mode_inertia)


if __name__ == "__main__":
    main()
