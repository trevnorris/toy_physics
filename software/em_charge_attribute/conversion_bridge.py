"""Weak electric gauge test along the saved relaxed-mouth paths, per signed unit."""
import json
from pathlib import Path

import sympy as sp


repo = Path(__file__).resolve().parents[2]
paths = json.loads((repo / "docs/_measurements/relaxed_neck_paths.json").read_text())
x, y = sp.symbols("x y", real=True)

# x is the distance travelled by the zero from the inner face, in units R.
# y is the independent radial coordinate measured from the same face.
def order(y):
    return 10*(1-y)**3 - 15*(1-y)**4 + 6*(1-y)**5


chi = order(y)
test_gauge = x*(1-x)
results = {}
for label, primitive, include_point in (
    ("layer alone", 1-chi, False),
    ("completed", 1-chi, True),
    ("localized primitive", chi*(1-chi), True),
):
    # The solid angle on radial sections has degree +1 only outside the zero.
    regular = sp.integrate(primitive.diff(y), (y, x, 1))
    point = primitive.subs(y, x) if include_point else sp.Integer(0)
    total = sp.simplify(regular + point)
    # A spatially uniform gauge test on the collar, zero at both time ends.
    # Its derivative pairs directly with the charge in A wedge J.
    variation = sp.integrate(test_gauge.diff(x)*total, (x, 0, 1))
    results[label] = (regular, point, total, variation)
    print(f"{label}: gauge pairing = {variation}, charge at entry = {total.subs(x,0)}, "
          f"charge just inside exit = {total.subs(x,1)}", flush=True)
    for name, path in paths.items():
        tilts = [sp.Rational(str(state["beta"])).limit_denominator(10000)
                 for state in path["states"] if 2 <= state["beta"] <= 3]
        if not tilts or tilts[0] != 2 or tilts[-1] != 3 or tilts != sorted(set(tilts)):
            raise ValueError(f"Invalid layer crossing in {name}")
        measured = [tuple(sp.factor(v.subs(x, t-2)) for v in (regular, point)) for t in tilts]
        print(f"{name} regular/point at tilts {tilts}: {measured}", flush=True)

# Independently use the weak divergence and the two end sections.
for label, (regular, point, total, variation) in results.items():
    weak_divergence = sp.integrate(test_gauge * total.diff(x), (x, 0, 1))
    assert sp.simplify(variation + weak_divergence) == 0, label
    if label != "layer alone":
        assert sp.simplify(total.diff(x)) == 0, label
