"""Far-field force at fixed mouth data in the signed bending channel."""

import sympy as sp


def check(label, calculated, expected):
    print(f"{label}: {sp.factor(calculated)}", flush=True)
    if sp.factor(calculated - expected) != 0:
        raise AssertionError(f"{label}: expected {expected}, got {calculated}")


# x = 1/(4 pi kappa d). The derivative of the interaction energy in x
# at zero equals 4 pi kappa times the limiting outward d^2 F.
x, s1, s2 = sp.symbols("x s1 s2", real=True)
S1, S2, z1, z2, phi1, phi2, q1, q2, j1, j2, g = sp.symbols(
    "S1 S2 z1 z2 phi1 phi2 q1 q2 j1 j2 g", positive=True
)
lam = sp.symbols("lambda", nonnegative=True)
response = sp.Matrix([[S1, x * z1 * z2], [x * z1 * z2, S2]])


def odd_force(energy):
    derivative = sp.diff(energy, x).subs(x, 0)
    return sp.factor(sum(
        i * j * derivative.subs({s1: i, s2: j})
        for i in (-1, 1) for j in (-1, 1)
    ) / 4)


def force_with_fixed_source(total1, total2, work1, work2):
    y = sp.Matrix([s1 * total1, s2 * total2])
    work = sp.Matrix([s1 * work1, s2 * work2])
    h = response * y
    return odd_force(y.dot(h) / 2 - work.dot(h))


h = sp.Matrix([s1 * phi1, s2 * phi2])
y = response.inv() * h
fixed_bend = odd_force(-y.dot(h) / 2)
check("held bend", fixed_bend, z1 * z2 * phi1 * phi2 / (S1 * S2))
check("held bend if reservoir work is omitted", odd_force(y.dot(h) / 2), -fixed_bend)
fixed_load = force_with_fixed_source(q1 + g, q2 + g, g, g)
fixed_source = force_with_fixed_source(j1 + g, j2 + g, j1 + g, j2 + g)
mixed = force_with_fixed_source(q1 + g, q2 + g, g + lam * q1, g + lam * q2)
check("fixed integrated load", fixed_load, z1 * z2 * (q1 * q2 - g**2))
check("fixed energy source", fixed_source, -z1 * z2 * (j1 + g) * (j2 + g))
check("fractional source work", mixed, z1 * z2 * (
    (1 - 2 * lam) * q1 * q2 - lam * g * (q1 + q2) - g**2
))
check("zero work fraction recovers fixed load", mixed.subs(lam, 0), fixed_load)
check("full work fraction recovers fixed source", mixed.subs(lam, 1),
      fixed_source.subs({j1: q1, j2: q2}))
small = mixed.subs({q2: q1, z2: z1})
large = mixed.subs({q1: q2, z1: z2})
check("two-radius determinant", small * large - mixed**2,
      -(1 - lam)**2 * (z1 * z2 * g * (q1 - q2))**2)
check("equal-data sign threshold", small.subs(lam, (q1 - g) / (2 * q1)), 0)


# A single exterior law and preparation for every pair. Units here use
# a = 1, 4 pi kappa = 1 and the unit datum at the small mouth.
def sharp_sphere_triplet(expression, data_at_radius):
    results = []
    for r1, r2 in ((1, 1), (1, 2), (2, 2)):
        p1, load1 = data_at_radius(sp.Integer(r1))
        p2, load2 = data_at_radius(sp.Integer(r2))
        results.append(expression.subs({
            S1: sp.Rational(1, r1), S2: sp.Rational(1, r2), z1: 1, z2: 1,
            phi1: p1, phi2: p2, q1: load1, q2: load2,
            j1: load1, j2: load2, g: 0,
        }))
    return tuple(results)


for name, expression, law, expected in [
    ("equal bend", fixed_bend, lambda r: (1, 1), (1, 2, 4)),
    ("bend proportional to radius", fixed_bend, lambda r: (r, 1), (1, 4, 16)),
    ("bend inversely proportional to radius", fixed_bend, lambda r: (1/r, 1), (1, 1, 1)),
    ("equal integrated load", fixed_load, lambda r: (1, 1), (1, 1, 1)),
    ("equal load per area", fixed_load, lambda r: (1, r**2), (1, 4, 16)),
    ("equal integrated energy source", fixed_source, lambda r: (1, 1), (-1, -1, -1)),
]:
    result = sharp_sphere_triplet(expression, law)
    print(f"{name}, small-small / small-large / large-large: {result}", flush=True)
    if result != expected:
        raise AssertionError((name, result, expected))
