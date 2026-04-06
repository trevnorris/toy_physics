import sympy as sp

eps, v, n, q, M0 = sp.symbols("eps v n q M0")

mass = M0 * (1 + q * eps)
cs2 = 1 + (n - 1) * eps
lagrangian = -mass * sp.sqrt(1 - v**2 / cs2)

expanded = sp.expand(
    sp.series(sp.series(lagrangian, eps, 0, 2).removeO(), v, 0, 5).removeO()
)

potential_coeff = sp.simplify(sp.diff(expanded, eps).subs({eps: 0, v: 0}) / M0)
mixed_coeff = sp.simplify(sp.diff(sp.diff(expanded, eps), v, 2).subs({eps: 0, v: 0}) / (2 * M0))

scalar_coeff = sp.Rational(1, 2) * q
vector_kinematic_coeff = -sp.Rational(1, 2) * (n - 1)
hybrid_coeff = sp.simplify(scalar_coeff + vector_kinematic_coeff)

solution = sp.solve(
    [
        sp.Eq(-potential_coeff, 1),
        sp.Eq(hybrid_coeff, -sp.Rational(3, 2)),
    ],
    (n, q),
    dict=True,
)

print("Expanded hybrid Lagrangian through O(eps v^2, v^4):")
print(expanded)
print()
print("Newtonian potential coefficient (-M0 q Phi):")
print(-potential_coeff)
print("Scalar-sector mixed coefficient:")
print(scalar_coeff)
print("Vector-kinematic mixed coefficient:")
print(vector_kinematic_coeff)
print("Combined hybrid mixed coefficient:")
print(hybrid_coeff)
print()
print("Solution from Newtonian matching and the GR mixed-term target:")
print(solution)
