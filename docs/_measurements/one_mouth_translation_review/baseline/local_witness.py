"""Exact two-coordinate release check. Toy parameters use fixed reference units."""
from sympy import Matrix, diff, simplify, symbols

f, U, c, k, b = symbols("f U c k b", real=True)
mobility = Matrix([[1, k], [k, 1]])  # positive definite for -1 < k < 1
flow = c * U * Matrix([1, b])
H = -c * U
support_channel = Matrix([1, b])
probe_channel = Matrix([1, 0])
held = mobility * (flow + H * support_channel)
same = (mobility * (flow + f * support_channel))[0]
changed = (mobility * (flow + f * probe_channel))[0]
print("held acceleration =", held.T)
for name, acceleration in (("matched actuator", same), ("different actuator", changed)):
    print(name, "X acceleration at f=H =", simplify(acceleration.subs(f, H)))
    measured = diff(acceleration, U)
    predicted = -diff(acceleration, f) * diff(H, U)
    print("  flow-only slope =", measured, "support and probe prediction =", predicted,
          "difference =", simplify(measured - predicted))
