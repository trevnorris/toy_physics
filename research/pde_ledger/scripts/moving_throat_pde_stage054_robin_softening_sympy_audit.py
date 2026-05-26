#!/usr/bin/env python3
"""
Stage 37 SymPy audit: Robin-compliance softening of the lowest support lane.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 37 — ROBIN-COMPLIANCE SOFTENING")

s, L, k, h, y, eta = sp.symbols("s L k h y eta", positive=True, real=True)
A, B = sp.symbols("A B", real=True)
pi = sp.pi

psi = A * sp.cos(k * s) + B * sp.sin(k * s)
# Bottom Neumann condition.
B_expr = sp.solve(sp.Eq(sp.diff(psi, s).subs(s, L), 0), B)[0]
psi_bn = sp.simplify(psi.subs(B, B_expr))
char_eq = sp.simplify(sp.diff(psi_bn, s).subs(s, 0) - h * psi_bn.subs(s, 0))

print("B from Neumann bottom =", B_expr)
print("Robin characteristic equation factor =", sp.factor(char_eq / A))
expect_zero("Robin equation -> k tan(kL) - h", char_eq / A + h - k * sp.tan(k * L))

# Dimensionless form.
expect_zero("dimensionless form", (k * sp.tan(k * L) - h).subs({k: y / L, h: eta / L}) - (y * sp.tan(y) - eta) / L)

banner("EXACT SOFTENING FACTOR")

KX, TX, KW, x = sp.symbols("K_X T_X K_W x", positive=True, real=True)
Kphi = KX + TX * y**2 / L**2
KW_expr = KX + pi**2 * TX / (4 * L**2)
AK = sp.simplify(KW_expr / Kphi)

print("K_W^(eff) =", KW_expr)
print("K_phi,0^(eff) =", Kphi)
print("A_K =", AK)

x_def = sp.Eq(x, pi**2 * TX / (L**2 * KW))
print("x definition:", x_def)
KX_from_x = sp.simplify(KW * (1 - x / 4))
TX_from_x = sp.simplify(x * L**2 * KW / pi**2)
expect_zero("K_W identity", KW - (KX_from_x + pi**2 * TX_from_x / (4 * L**2)))
AK_x = sp.simplify((KW / (KX + TX * y**2 / L**2)).subs({KX: KX_from_x, TX: TX_from_x}))
print("A_K in x,y form =", AK_x)
expect_zero("A_K x-form", AK_x - 1 / (1 - x / 4 + x * y**2 / pi**2))

AK_sym = sp.simplify(1 / (1 - x / 4 + x * y**2 / pi**2))
AK_DN = sp.simplify(AK_sym.subs(y, pi / 2))
AK_soft = sp.simplify(sp.limit(AK_sym, y, 0, dir="+"))
print("A_K(y=pi/2) =", AK_DN)
print("A_K(y->0+) =", AK_soft)
expect_zero("DN limit", AK_DN - 1)
expect_zero("soft-mouth limit", AK_soft - 4 / (4 - x))

# Monotonicity certificate: A_K is strictly decreasing in y on (0, pi/2)
# with 0 < x < 4. Verify the closed form of dA_K/dy; positivity of the
# prefactor 2*x*y/pi^2 on the assumed domain then implies dA_K/dy < 0,
# which together with the endpoint values brackets A_K in [1, 4/(4-x)].
dAK_dy = sp.diff(AK_sym, y)
dAK_dy_expected = -2 * x * y / (pi**2 * (1 - x / 4 + x * y**2 / pi**2) ** 2)
expect_zero("dA_K/dy closed form", dAK_dy - dAK_dy_expected)
print("Prefactor 2*x*y/pi^2 > 0 on 0<x<4, 0<y<pi/2 => dA_K/dy < 0 (monotone decreasing).")

banner("PURE SOFTENING THRESHOLD")

zeta_req = sp.symbols("zeta_req", positive=True, real=True)
ineq_rhs = sp.simplify(1 / zeta_req - 1 + x / 4)
y_req_sq = sp.simplify(pi**2 * ineq_rhs / x)
print("Condition A_K >= zeta_req implies y^2 <=", y_req_sq)

# Maximum reachable pure-softening ratio.
AK_max = sp.simplify(4 / (4 - x))
print("A_K,max =", AK_max)

x_floor = sp.simplify(sp.solve(sp.Eq(AK_max, zeta_req), x)[0])
print("x floor at saturation =", x_floor)
expect_zero("x floor = 4 - 4/zeta_req", x_floor - (4 - 4 / zeta_req))

print("\nFINAL LEDGER")
print("  Robin compliance replaces the D/N half-wave by y tan y = eta.")
print("  The support softening factor is A_K = 1 / [1 - x/4 + x y^2/pi^2].")
print("  It ranges from 1 up to 4/(4-x).")
print("  Pure softening rescue is possible only if zeta_req <= 4/(4-x).")
