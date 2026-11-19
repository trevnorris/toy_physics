"""Work of an attached wall and its orientation-odd far-field interaction."""

import json
import re
from pathlib import Path

import sympy as sp


def verify(label, expr, expected):
    result = sp.factor(expr)
    print(f"{label} = {result}", flush=True)
    if sp.simplify(result - expected) != 0:
        raise AssertionError(f"{label}: expected {expected}")


# Differentiate with isolated mouth preparations fixed. The local attachment
# energy can depend on additional internal coordinates. All free coordinates
# are stationary at the point of evaluation.
x, za, zb, Sa, Sb = sp.symbols("x z_a z_b S_a S_b", positive=True)
ua, ub = sp.symbols("u_a u_b", real=True)
ha, hb = sp.symbols("h_a h_b", real=True)
h = sp.Matrix([ha, hb])
G = sp.Matrix([[Sa, x*za*zb], [x*za*zb, Sb]])
field_energy = h.dot(G.inv()*h)/2
attachment = sp.Function("V_a")(ha, ua) + sp.Function("V_b")(hb, ub)
energy = field_energy + attachment
variables = (ha, hb, ua, ub)
rates = sp.symbols("h_a_prime h_b_prime u_a_prime u_b_prime")
along_branch = energy.diff(x) + sum(energy.diff(v)*vp for v, vp in zip(variables, rates))
stationary = {attachment.diff(v): -field_energy.diff(v).subs(x, 0) for v in variables}
ya, yb = sp.symbols("y_a y_b", real=True)
leading = along_branch.subs(x, 0).subs(stationary).subs({ha: Sa*ya, hb: Sb*yb})
verify("stationary attachment pair energy slope in x", leading, -za*zb*ya*yb)

# A finite mechanical contact independently illustrates the distinction
# between fixing a spring reference and supplying conjugate reservoir work.
ka, kb, pa, pb = sp.symbols("k_a k_b H_a H_b", positive=True)
sa, sb = sp.symbols("s_a s_b", real=True)
springs = ka*(ha-sa*pa)**2/2 + kb*(hb-sb*pb)**2/2
isolated = sp.solve([sp.diff(field_energy+springs, v).subs(x, 0) for v in h], list(h))
spring_coeff = field_energy.diff(x).subs(x, 0).subs(isolated)/(sa*sb)
verify("two spring references", spring_coeff, -za*zb*ka*pa*kb*pb/((1+Sa*ka)*(1+Sb*kb)))
verify("infinite stiffness at fixed references", sp.limit(sp.limit(spring_coeff, ka, sp.oo), kb, sp.oo),
       -za*zb*pa*pb/(Sa*Sb))
verify("source-work reservoir form at held field", (-field_energy).diff(x).subs(x, 0),
       za*zb*ha*hb/(Sa*Sb))

# Endpoint virtual work at fixed radius and shape during an infinitesimal
# local bend. Keep both displacements live, even when they later coincide.
radius, ell, sigma, bend = sp.symbols("a ell sigma k_bend", positive=True)
near, far = sp.Function("near")(ha), sp.Function("far")(ha)
pull = 4*sp.pi*sigma*radius**2 + 16*sp.pi*bend/9
wall = sa*pull*(far-near)
verify("wall load per unit bend", -wall.diff(ha), sa*pull*(near.diff(ha)-far.diff(ha)))
verify("comoving wall ends", (-wall.diff(ha)).subs(far, near).doit(), 0)

manifest = json.loads(Path("docs/_measurements/unequal_mouth_electric.json").read_text())
readings = manifest["runs"]["annular_mouth_radius_sympy"]["stdout_claim_lines"]
escape = {}
for line in readings:
    match = re.fullmatch(r"radius (\d+), escape = ([\d.]+), self response = [\d.]+", line)
    if match:
        escape[int(match[1])] = sp.Rational(match[2])

# f is the tension fraction of the small-mouth pull. Nonnegative tension
# and bending give 0 <= f <= 1. Print limits at both single-term endpoints.
f = sp.symbols("f", nonnegative=True)
ratio = sp.factor(pull.subs(radius, 2*radius)/pull)
tension = 4*sp.pi*sigma*radius**2
verify("wall load ratio at doubled radius", ratio, 1+3*tension/pull)
gain = (escape[1]/escape[2])/(1+3*f)
print("equal transfer, large/small exterior strength, curvature and tension limits =",
      [float((1/gain).subs(f, end)) for end in (0, 1)])
print("required relative mechanical transfer, curvature and tension limits =",
      [float(gain.subs(f, end)) for end in (0, 1)])
verify("relative transfer times relative emitted load", gain*(1+3*f)*escape[2]/escape[1], 1)
