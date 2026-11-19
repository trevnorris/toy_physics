"""Static radial exterior of the two-interface slab, driven by signed rim loads."""

import argparse
from pathlib import Path

import sympy as sp


def check(label, value, expected):
    print(f"{label} = {value}", flush=True)
    if sp.simplify(value - expected) != 0:
        raise ArithmeticError(label)


def plot_faces(path, radius, profiles):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    distances = np.linspace(1, 5, 300)
    fig, ax = plt.subplots(figsize=(6.6, 3.8), layout="constrained")
    for (label, color, style), profile in zip([
        ("Loaded face", "#9c4317", "-"), ("Midplane", "#505968", "--"),
        ("Opposite face", "#106d77", "-"),
    ], profiles):
        values = sp.lambdify(radius, profile, "numpy")(distances)
        ax.plot(distances, values, style, color=color, linewidth=2.2, label=label)
    ax.set(xlabel="Distance from mouth center (rim radii)", ylabel="Displacement / midplane rise at rim",
           title="Equal face tensions, thickness decay length = rim radius", xlim=(1, 5), ylim=(0, 1.6))
    ax.grid(axis="y", alpha=0.22)
    ax.legend(frameon=False)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()
    r, a, length, km, kt = sp.symbols("r a ell K_m K_t", positive=True)
    fu, fl = sp.symbols("F_up F_low", real=True)
    m, t = (sp.Function(i)(r) for i in ("m", "t"))
    jac = 4 * sp.pi * r**2
    density = jac * (km*m.diff(r)**2 + kt*t.diff(r)**2 + kt*t**2/length**2)/2
    work = fu*(m+t) + fl*(m-t)
    fields = []
    for field, basis, label in [(m, 1/r, "middle"),
                                (t, sp.exp(-(r-a)/length)/r, "half width")]:
        # Vary the radial exterior energy and include the inner boundary sign.
        flux = sp.diff(density, field.diff(r))
        operator = sp.diff(flux, r) - sp.diff(density, field)
        check(label + " homogeneous equation", sp.simplify(operator.subs(field, basis).doit()), 0)
        unit_load = -flux.subs(field, basis).doit().subs(r, a)
        mode = sp.factor(sp.diff(work, field) / unit_load) * basis
        fields.append(mode)
        print(f"{label} profile = {mode}", flush=True)
        check(label + " load residual", sp.simplify(-flux.subs(field, mode).doit().subs(r, a)
              - sp.diff(work, field)), 0)
    middle, half = fields
    f = sp.Symbol("F", positive=True)
    top = {fu:f, fl:0}
    mirrored = {fu:0, fl:-f}
    check("one-face middle slope times r^2", sp.factor(middle.subs(top).diff(r)*r**2),
          -f/(4*sp.pi*km))
    check("mirrored middle sum", middle.subs(top)+middle.subs(mirrored), 0)
    check("mirrored full width difference", 2*(half.subs(top)-half.subs(mirrored)), 0)
    check("balanced opposite loads middle", middle.subs({fu:f, fl:-f}), 0)
    contrast = sp.factor((half/middle).subs(top))
    print(f"face difference over face sum = {contrast}", flush=True)
    lower = (middle-half).subs(top)
    upper = (middle+half).subs(top)
    check("unlocked equal-tension lower face", sp.limit(lower.subs(kt,km),length,sp.oo), 0)
    check("finite restoring curvature distant face ratio", sp.limit(lower/upper,r,sp.oo), 1)
    check("rim ratio for equal tensions", sp.factor((lower/upper).subs({kt:km,r:a})), a/(a+2*length))
    rp = sp.Symbol("dr", positive=True)
    check("contrast attenuation over dr", sp.simplify(contrast.subs(r,r+rp)/contrast),
          sp.exp(-rp/length))
    for factor in (1, 2, 3):
        example = sp.simplify((lower/upper).subs({kt:km,length:a,r:factor*a}, simultaneous=True))
        print(f"equal tensions, ell=a, r/a={factor}: lower/upper = {float(example):.9f}", flush=True)
    if args.plot:
        norm = middle.subs(top).subs(r,a)
        profiles = [sp.simplify((field/norm).subs({kt:km, length:1, a:1}))
                    for field in (upper, middle.subs(top), lower)]
        plot_faces(args.plot, r, profiles)


if __name__ == "__main__":
    main()
