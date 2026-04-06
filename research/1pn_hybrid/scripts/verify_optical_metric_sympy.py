import sympy as sp

r, beta = sp.symbols("r beta", positive=True)
N = sp.Function("N")

turning_point_relation = sp.Eq(sp.Symbol("b", positive=True) ** 2, N(r) ** 4 * r**2)
circular_condition = sp.diff(N(r) ** 4 * r**2, r)
log_condition = sp.simplify(circular_condition / (2 * N(r) ** 3 * r))

N_simple = 1 + beta / r
rph_solution = sp.solve(sp.Eq(sp.diff(N_simple**4 * r**2, r), 0), r)[0]
bph_solution = sp.simplify((N_simple**2 * r).subs(r, rph_solution))

print("Optical-metric turning-point relation:")
print(turning_point_relation)
print()
print("Circular-orbit condition d/dr[N(r)^4 r^2] = 0 becomes:")
print(sp.Eq(log_condition, 0))
print("Equivalent logarithmic form:")
print(sp.Eq(sp.diff(sp.log(N(r)), r), -sp.Rational(1, 2) / r))
print()
print("For N(r) = 1 + beta/r:")
print("r_ph =", sp.simplify(rph_solution))
print("b_ph =", sp.simplify(bph_solution))
