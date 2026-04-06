import sympy as sp

r, w, a, L, t = sp.symbols("r w a L t", positive=True, real=True)
c_s, A = sp.symbols("c_s A", positive=True, real=True)
omega_11 = sp.symbols("omega_11", real=True)

x = sp.symbols("x", real=True)
x01 = sp.N(sp.nsolve(sp.besselj(0, x), 2.4))

k_r = x01 / a
k_w = sp.pi / L
omega_11_sq = sp.simplify(c_s**2 * (k_r**2 + k_w**2))
h_11 = A * sp.besselj(0, k_r * r) * sp.sin(k_w * w) * sp.cos(omega_11 * t)
L_over_a = sp.sqrt(2) * sp.pi / x01

print("First J0 root x_01 ≈", x01)
print("Radial wavenumber k_r =", k_r)
print("Axial wavenumber k_w =", k_w)
print("Dispersion relation omega_11^2 =", omega_11_sq)
print("Fundamental throat mode h_11(t, r, w) =")
sp.pprint(h_11)
print("Predicted aspect ratio L/a ≈", sp.N(L_over_a))
