import sympy as sp

# ----------------------------------------------------------------------
# 1. Symbols and constants
# ----------------------------------------------------------------------
r, w, a, L, t = sp.symbols('r w a L t', positive=True, real=True)
cs, A = sp.symbols('c_s A', positive=True, real=True)

# Bessel root x_01: first zero of J0
x = sp.symbols('x', real=True)
J0 = sp.besselj(0, x)
x01 = sp.nsolve(J0, 2.4)  # initial guess near first zero
x01 = sp.N(x01)
print("First J0 root x_01 ≈", x01)

# ----------------------------------------------------------------------
# 2. Wavenumbers and mode frequency
# ----------------------------------------------------------------------
k_r = x01 / a
k_w = sp.pi / L

omega = cs * sp.sqrt(k_r**2 + k_w**2)
print("Mode frequency omega =", omega)

# ----------------------------------------------------------------------
# 3. Fundamental 4D mode h(r, w, t)
# ----------------------------------------------------------------------
h = A * sp.besselj(0, k_r * r) * sp.sin(k_w * w) * sp.exp(-sp.I * omega * t)
print("Fundamental 4D mode h(r, w, t) =")
sp.pprint(h)

# ----------------------------------------------------------------------
# 4. Aspect ratio L/a from Paper IV result: L/a = sqrt(2)*pi/x01
# ----------------------------------------------------------------------
L_over_a = sp.sqrt(2) * sp.pi / x01
print("Predicted aspect ratio L/a ≈", sp.N(L_over_a))

