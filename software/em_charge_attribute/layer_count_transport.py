"""Balance the solid-angle current at the moving zero on the saved relaxed flares."""
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp


def check(name, value, target):
    value = sp.factor(value)
    print(f"{name} = {value}", flush=True)
    if sp.simplify(value - target) != 0:
        raise AssertionError(name)


# Use x/R as material chart coordinates. Smooth the angular area form
# in slope space before taking its distributional limit at the zero.
x, y, z, beta, speed = sp.symbols("x y z beta v", real=True)
a, length = sp.symbols("a L", positive=True)
coords = (x, y, z)
offset = sp.Matrix([x, y, z-beta])
vel = sp.Matrix([0, 0, speed])  # v = d beta/dt, freely chosen path clock
field = offset / (4*sp.pi*(offset.dot(offset) + a**2)**sp.Rational(3, 2))


def div(vector):
    return sum(vector[i].diff(c) for i, c in enumerate(coords))


kernel = 3*a**2 / (4*sp.pi*(offset.dot(offset) + a**2)**sp.Rational(5, 2))
check("divergence of smoothed solid angle / kernel", div(field)/kernel, 1)
enclosed = 4*sp.pi*length**2*field[2].subs({x: 0, y: 0, z: beta+length})
check("defect index from a surrounding sphere", sp.limit(enclosed, a, 0, dir="+"), 1)
check("radial integral derivative", enclosed.diff(length) -
      4*sp.pi*length**2*kernel.subs({x: 0, y: 0, z: beta+length}), 0)

# Components of -d chi ^ omega for a stationary order profile. Pulling
# omega back by the translated slope produces both density and current.
chi = sp.Function("chi")(*coords)
grad = sp.Matrix([chi.diff(c) for c in coords])
layer_density = -grad.dot(field)
layer_flux = grad.dot(vel)*field + layer_density*vel
check("flux normal to an order surface", grad.dot(layer_flux), 0)
residual = speed*layer_density.diff(beta) + div(layer_flux)
check("layer continuity residual / kernel", sp.cancel(residual/kernel), grad.dot(vel))
whole_density = layer_density + (1-chi)*kernel
whole_flux = layer_flux + (1-chi)*kernel*vel
check("layer plus weighted moving zero residual",
      sp.cancel(speed*whole_density.diff(beta) + div(whole_flux)), 0)

# Same radial loss as the saved paths. The attachment modes are supported
# on r/R < 2. Thus every stored state uses the same layer geometry at its beta.
u = sp.symbols("u", real=True)
loss_weight = 30*u**2*(1-u)**2
remaining = sp.integrate(loss_weight, (u, beta-2, 1))
check("residual per unit tilt within conversion", remaining.diff(beta),
      -loss_weight.subs(u, beta-2))
check("integrated local change for upper orientation",
      sp.integrate(remaining.diff(beta), (beta, 2, 3)), -1)
check("integrated local change for mirror",
      sp.integrate(-remaining.diff(beta), (beta, 2, 3)), 1)
check("residual per unit tilt at halfway", remaining.diff(beta).subs(beta, sp.Rational(5, 2)),
      -sp.Rational(15, 8))

repo = Path(__file__).resolve().parents[2]
saved = json.loads((repo / "docs/_measurements/relaxed_neck_paths.json").read_text())
for name, path in saved.items():
    states = path["states"]
    pieces = []
    previous = 1
    for state in states:
        value = sp.Rational(Fraction(state["beta"]).limit_denominator(10**4))
        if abs(float(value)-state["beta"]) > 1.e-12:
            raise ValueError("Unresolved tilt in saved path")
        layer = sp.Integer(1) if value < 2 else remaining.subs(beta, value)
        if value > 3 or layer > previous or layer < 0:
            raise ValueError("Invalid path range or ordering")
        previous = layer
        if value >= 2:
            pieces.append(f"{value}:{layer}")
    print(f"{name} ({len(states)} shapes), tilt:layer = {', '.join(pieces)}")
    assert states[0]["beta"] == 0 and states[-1]["beta"] == 3 and previous == 0
