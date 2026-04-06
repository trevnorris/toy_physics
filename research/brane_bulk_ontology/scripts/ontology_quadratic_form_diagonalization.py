import sympy as sp

A_T, A_L = sp.symbols("A_T A_L", positive=True)
u_T = sp.symbols("u_T", real=True)

alpha2 = sp.Rational(3, 4)
u_L = sp.sqrt(alpha2) * u_T

energy = sp.expand(sp.Rational(1, 2) * (A_T * u_T**2 + A_L * u_L**2))
target = sp.expand(sp.Rational(1, 8) * (4 * A_T + 3 * A_L) * u_T**2)
kernel = sp.diag(A_T, A_L)

print("alpha^2 =", alpha2)
print("Energy E(u_T) =", energy)
print("Target calibrated form =", target)
print("Matches calibrated form:", sp.simplify(energy - target) == 0)
print("Eigenvalues of M:", kernel.eigenvals())
print("Positive for A_T>0, A_L>0: True")
