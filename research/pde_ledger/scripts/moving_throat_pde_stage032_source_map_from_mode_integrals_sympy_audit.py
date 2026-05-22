#!/usr/bin/env python3
"""
Moving-throat PDE Stage 15 SymPy audit.

Checks:
1. Exact finite-throat axial integrals for the N/N wall basis and D/N half-wave.
2. Local-kernel reduction of wall/internal/source couplings.
3. Exact Schur-complement decomposition of the reduced wall operator.
4. Exact selected-source map on the natural D/N source branch.
5. Elimination of the abstract selected-branch source-map factor.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr):
    if isinstance(expr, sp.MatrixBase):
        simp = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(simp)
        if any(entry != 0 for entry in simp):
            raise AssertionError(f"{name} is not zero")
    else:
        simp = sp.simplify(sp.expand(expr))
        print(f"{name} = {simp}")
        if simp != 0:
            raise AssertionError(f"{name} is not zero")


banner("STAGE 15.1 — EXACT FINITE-THROAT AXIAL INTEGRALS")

L, s = sp.symbols("L s", positive=True, real=True)
u0 = 1 / sp.sqrt(L)
u1 = sp.sqrt(sp.Integer(2) / L) * sp.cos(sp.pi * s / L)
f0 = sp.sqrt(sp.Integer(2) / L) * sp.sin(sp.pi * s / (2 * L))

kappa0 = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))
kappa1 = sp.simplify(sp.integrate(u1 * f0, (s, 0, L)))

print("u0 . u0 =", sp.simplify(sp.integrate(u0 * u0, (s, 0, L))))
print("u1 . u1 =", sp.simplify(sp.integrate(u1 * u1, (s, 0, L))))
print("u0 . u1 =", sp.simplify(sp.integrate(u0 * u1, (s, 0, L))))
print("f0 . f0 =", sp.simplify(sp.integrate(f0 * f0, (s, 0, L))))
print("kappa0 =", kappa0)
print("kappa1 =", kappa1)
print("sigma =", sp.simplify(kappa0**2 + kappa1**2))
print("sigma/kappa0^2 =", sp.simplify((kappa0**2 + kappa1**2) / kappa0**2))

expect_zero("u0.u0 - 1", sp.integrate(u0 * u0, (s, 0, L)) - 1)
expect_zero("u1.u1 - 1", sp.integrate(u1 * u1, (s, 0, L)) - 1)
expect_zero("u0.u1", sp.integrate(u0 * u1, (s, 0, L)))
expect_zero("f0.f0 - 1", sp.integrate(f0 * f0, (s, 0, L)) - 1)
expect_zero("kappa0 - 2 sqrt(2)/pi", kappa0 - 2 * sp.sqrt(2) / sp.pi)
expect_zero("kappa1 + 4/(3 pi)", kappa1 + sp.Rational(4, 3) / sp.pi)
expect_zero("sigma - 88/(9 pi^2)", kappa0**2 + kappa1**2 - sp.Rational(88, 9) / sp.pi**2)
expect_zero("sigma/kappa0^2 - 11/9", (kappa0**2 + kappa1**2) / kappa0**2 - sp.Rational(11, 9))

banner("STAGE 15.2 — LOCAL-KERNEL MODE REDUCTION")

q0, q1 = sp.symbols("q0 q1", real=True)
U0, U1 = sp.symbols("U0 U1", real=True)
phi, W, Qstf = sp.symbols("phi W Qstf", real=True)
gU, gB, gW, gR, gQ = sp.symbols("g_U g_B g_W g_R g_Q", real=True)

eta = q0 * u0 + q1 * u1
Ufield = U0 * u0 + U1 * u1
phifield = phi * f0
Wfield = W * f0

L_etaU = sp.simplify(gU * sp.integrate(eta * Ufield, (s, 0, L)))
L_etaphi = sp.simplify(gB * sp.integrate(eta * phifield, (s, 0, L)))
L_etaW = sp.simplify(gW * sp.integrate(eta * Wfield, (s, 0, L)))
L_UW = sp.simplify(-gR * sp.integrate(Ufield * Wfield, (s, 0, L)))
L_src = sp.simplify(gQ * Qstf * sp.integrate(eta * f0, (s, 0, L)))

print("L_etaU =", L_etaU)
print("L_etaphi =", L_etaphi)
print("L_etaW =", L_etaW)
print("L_UW =", L_UW)
print("L_src =", L_src)

expect_zero(
    "L_etaU - gU (q.U)",
    L_etaU - gU * (q0 * U0 + q1 * U1),
)
expect_zero(
    "L_etaphi - gB (v.q) phi",
    L_etaphi - gB * (kappa0 * q0 + kappa1 * q1) * phi,
)
expect_zero(
    "L_etaW - gW (v.q) W",
    L_etaW - gW * (kappa0 * q0 + kappa1 * q1) * W,
)
expect_zero(
    "L_UW + gR (v.U) W",
    L_UW + gR * (kappa0 * U0 + kappa1 * U1) * W,
)
expect_zero(
    "L_src - gQ Q (v.q)",
    L_src - gQ * Qstf * (kappa0 * q0 + kappa1 * q1),
)

banner("STAGE 15.3 — EXACT SCHUR-COMPLEMENT DECOMPOSITION")

D0, D1 = sp.symbols("D0 D1", real=True)
Aphi, AU, AW = sp.symbols("A_phi A_U A_W", real=True)

v = sp.Matrix([kappa0, kappa1])
I2 = sp.eye(2)
Deta = sp.diag(D0, D1)

# Internal block x = (U0,U1,phi,W)
Kint = sp.Matrix([
    [AU, 0, 0, -gR * kappa0],
    [0, AU, 0, -gR * kappa1],
    [0, 0, Aphi, 0],
    [-gR * kappa0, -gR * kappa1, 0, AW],
])
Bmat = sp.Matrix([
    [gU, 0, gB * kappa0, gW * kappa0],
    [0, gU, gB * kappa1, gW * kappa1],
])
Sigma = sp.simplify(Bmat * Kint.inv() * Bmat.T)

sigma = sp.simplify((v.T * v)[0])
Xi = sp.simplify(gU**2 / AU)
alpha = sp.simplify(gB**2 / Aphi + (AU * gW + gR * gU)**2 / (AU * (AU * AW - gR**2 * sigma)))
Sigma_target = sp.simplify(Xi * I2 + alpha * (v * v.T))

print("Sigma =")
sp.pprint(Sigma)
print("Xi =", Xi)
print("alpha =", alpha)
expect_zero("Sigma - [Xi I + alpha vv^T]", Sigma - Sigma_target)

banner("STAGE 15.4 — NATURAL D/N SOURCE MAP")

A, DK, alpha0 = sp.symbols("A DeltaK alpha0", positive=True, real=True)
B = A + DK
sigma_sym = sp.symbols("sigma", positive=True, real=True)
delta_kappa = sp.symbols("delta_kappa", real=True)
Kprod = sp.symbols("Kprod", positive=True, real=True)

R = sp.sqrt((DK + alpha0 * delta_kappa)**2 + 4 * alpha0**2 * Kprod)
lambda_minus = sp.simplify((A + B - alpha0 * sigma_sym - R) / 2)
s_minus = sp.simplify(sp.Rational(1, 2) * (sigma_sym + ((DK + alpha0 * delta_kappa) * delta_kappa + 4 * alpha0 * Kprod) / R))

subs_nat = {
    sigma_sym: kappa0**2 + kappa1**2,
    delta_kappa: kappa0**2 - kappa1**2,
    Kprod: kappa0**2 * kappa1**2,
}
s_minus_nat = sp.simplify(s_minus.subs(subs_nat))
mhat_sq = sp.simplify(s_minus_nat / kappa0**2)

print("mhat_-^2 =", mhat_sq)
expect_zero("mhat_-^2(alpha=0) - 1", sp.simplify(mhat_sq.subs(alpha0, 0) - 1))

# Non-trivial identity on the natural-D/N kappa products.
expect_zero(
    "delta_kappa^2 + 4*Kprod - sigma^2 (natural)",
    (delta_kappa**2 + 4 * Kprod - sigma_sym**2).subs(subs_nat),
)

# Interior consistency: s_minus_nat at alpha0 = 1, DK = 1 equals the
# closed form obtained using the natural-D/N identity above.
R_nat = sp.sqrt(DK**2 + 2 * alpha0 * DK * delta_kappa + alpha0**2 * sigma_sym**2)
s_minus_nat_simplified = sp.simplify(
    (sigma_sym + (DK * delta_kappa + alpha0 * sigma_sym**2) / R_nat) / 2
)
expect_zero(
    "s_minus_nat - s_minus_nat_simplified (interior identity)",
    sp.simplify((s_minus_nat - s_minus_nat_simplified.subs(subs_nat))),
)
expect_zero(
    "s_minus_nat at (alpha0=1, DK=1) interior point",
    sp.simplify(
        s_minus_nat.subs({alpha0: sp.Integer(1), DK: sp.Integer(1)})
        - s_minus_nat_simplified.subs({alpha0: sp.Integer(1), DK: sp.Integer(1)}).subs(subs_nat)
    ),
)
expect_zero("limit_{alpha->oo} mhat_-^2 - 11/9", sp.simplify(sp.limit(mhat_sq, alpha0, sp.oo) - sp.Rational(11, 9)))

banner("STAGE 15.5 — ELIMINATION OF THE ABSTRACT SOURCE-MAP FACTOR")

beta0 = sp.symbols("beta0", positive=True, real=True)

# Closed-form product under the natural D/N substitution.
P0_minus_nat = sp.simplify((beta0 * s_minus / lambda_minus).subs(subs_nat))
Nprod_nat = sp.simplify((s_minus / kappa0**2).subs(subs_nat) * P0_minus_nat)

# (i) At alpha0 = 0 the elimination must yield beta0 * kappa0^2 / A.
Nprod_alpha0 = sp.simplify(Nprod_nat.subs(alpha0, 0))
expect_zero(
    "Nprod(alpha=0) - beta0 * kappa0^2 / A",
    Nprod_alpha0 - beta0 * kappa0**2 / A,
)

# (ii) As alpha0 -> oo the product vanishes (lambda_minus diverges, s_minus is finite).
Nprod_inf = sp.limit(Nprod_nat, alpha0, sp.oo)
expect_zero("limit_{alpha->oo} Nprod_nat", Nprod_inf)

print("All Stage 15 checks passed.")
