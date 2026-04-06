import sympy as sp

G, M, r, b, c, n, z = sp.symbols("G M r b c n z", positive=True)
Z_E, Z_R = sp.symbols("Z_E Z_R", positive=True)

alpha_n = sp.simplify((n - 1) / 2)
delta_rho_over_rho = -G * M / (c**2 * r)
delta_cs_over_c = sp.simplify(alpha_n * delta_rho_over_rho)
N_of_r = sp.simplify(1 - delta_cs_over_c)

lensing_integrand = alpha_n * G * M * b / (c**2 * (b**2 + z**2) ** sp.Rational(3, 2))
delta_theta = sp.simplify(sp.integrate(lensing_integrand, (z, -sp.oo, sp.oo)))

shapiro_exact = sp.simplify(
    alpha_n
    * G
    * M
    / c**3
    * sp.integrate(1 / sp.sqrt(b**2 + z**2), (z, -Z_E, Z_R))
)
shapiro_asym_n5 = 2 * G * M / c**3 * sp.log(4 * Z_E * Z_R / b**2)

print("SymPy optics cross-check for the preferred 1PN branch")
print(f"alpha_n = (n-1)/2 = {sp.sstr(alpha_n)}")
print(f"delta_rho/rho = {sp.sstr(delta_rho_over_rho)}")
print(f"N(r) = {sp.sstr(N_of_r)}")
print()

print("General-n lensing coefficient:")
print(f"Delta theta_n = {sp.sstr(delta_theta)}")
print(f"At n=5: {sp.sstr(sp.simplify(delta_theta.subs(n, 5)))}")
print(
    "Checks 4GM/(b c^2):",
    sp.simplify(delta_theta.subs(n, 5) - 4 * G * M / (b * c**2)) == 0,
)
print()

print("General-n Shapiro coefficient:")
print(f"Exact one-way delay = {sp.sstr(shapiro_exact)}")
print(f"At n=5: {sp.sstr(sp.simplify(shapiro_exact.subs(n, 5)))}")
print(f"Weak-field asymptotic n=5 formula = {sp.sstr(shapiro_asym_n5)}")
print()

print("Mass-based redshift (q=1):")
print(f"Delta m / m = {sp.sstr(delta_rho_over_rho)}")
print("Magnitude matches the weak-field GR redshift potential.")
