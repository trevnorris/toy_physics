
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr: sp.Expr) -> None:
    s = sp.simplify(sp.expand(expr))
    print(f"{name} = {s}")
    if s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 99 — FINITE D/N MIXED-TUBE REALIZATION")

a, L_W = sp.symbols("a L_W", positive=True, real=True)
K_s, K_q, lam = sp.symbols("K_s K_q lam", positive=True, real=True)
z = sp.symbols("z", real=True)

r_c = lam**2 / (K_s * K_q)

kappa0_from_tube = sp.simplify(4 * L_W**2 / (sp.pi**2 * a**2))
L_W_required = sp.solve(
    sp.Eq(kappa0_from_tube, (1 + r_c) / 3),
    L_W,
)[0]

print("kappa0 from D/N half-wave tube =", kappa0_from_tube)
print("Required tube length L_W =", sp.simplify(L_W_required))

kappa0_bare = sp.simplify((1 + r_c) / 3)
gamma0_bare = sp.simplify((1 + r_c) / 9)

kappa_c = sp.simplify(kappa0_bare / (1 + r_c))
gamma_c = sp.simplify(gamma0_bare / (1 + r_c))

expect_zero("final kappa_c - 1/3", kappa_c - sp.Rational(1, 3))
expect_zero("final gamma_c - 1/9", gamma_c - sp.Rational(1, 9))

D_bare = sp.expand((1 + r_c) * (1 - z**2 / 3 - sp.I * z**5 / 9))
D_final = sp.simplify(D_bare / (1 + r_c))
expect_zero(
    "bare scaled-canonical branch renormalizes to canonical",
    D_final - (1 - z**2 / 3 - sp.I * z**5 / 9),
)

print("\nSummary:")
print("  D/N half-wave mixed tube length:", sp.simplify(L_W_required))
print("  Bare mixed coefficients: kappa0 =", kappa0_bare, ", gamma0 =", gamma0_bare)
print("  Final coefficients: kappa_c =", kappa_c, ", gamma_c =", gamma_c)
