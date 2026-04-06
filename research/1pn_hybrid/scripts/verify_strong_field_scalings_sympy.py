import sympy as sp

a, a0, c, sigma = sp.symbols("a a0 c sigma", positive=True)

Zth = sp.simplify(c / a**2)
Mmodel = a * (1 + sigma / (a**2 + a0**2))
k_eff = sp.simplify(a * sp.diff(sp.log(Mmodel), a))

print("Phenomenological throat-impedance scaling:")
print("Z_th(a) =", Zth)
print()
print("Representative mass model:")
print("M(a) =", Mmodel)
print("k_eff(a) =", k_eff)
print("Limit[k_eff, a->oo] =", sp.limit(k_eff, a, sp.oo))
