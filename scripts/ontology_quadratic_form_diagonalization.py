import sympy as sp

# Parameters
A_T, A_L = sp.symbols('A_T A_L', positive=True)
alpha2 = -sp.Rational(2, 5)  # alpha^2 = -2/5
alpha = sp.sqrt(alpha2)      # this will be imaginary

# Basis vector: [u_T, u_L] with u_L = alpha * u_T
u_T = sp.symbols('u_T')
u = sp.Matrix([u_T, alpha * u_T])

# Euclidean kernel on (T,L)
M = sp.diag(A_T, A_L)

E = sp.simplify((u.T * M * u)[0])
print("Energy E(u_T) =", E)

# Diagonalize the quadratic form
evals = sp.simplify(M.eigenvals())
print("Eigenvalues of M:", evals)

