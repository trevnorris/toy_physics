#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("STAGE 192 — EXACT PAIRWISE RATIO OPTIMIZER AND MIXED-RAY WINNER THEOREM")

# ---------------------------------------------------------------------------
# I. Exact algebraic form of the pairwise certified objective
# ---------------------------------------------------------------------------
subbanner("I. Exact algebraic form of the pairwise certified objective")

ki, kj, H0, r = sp.symbols("k_i k_j H0 r", positive=True, real=True)
u, v, w = sp.symbols("u v w", real=True)

kij = sp.simplify((ki + kj * r) / sp.sqrt(1 + r**2))
kappa = sp.simplify((u + 2 * v * r + w * r**2) / (1 + r**2))

tau = sp.simplify(2 * H0 / (kij + sp.sqrt(kij**2 - 2 * H0 * kappa)))

A = sp.simplify(ki**2 - 2 * H0 * u)
B = sp.simplify(2 * ki * kj - 4 * H0 * v)
C = sp.simplify(kj**2 - 2 * H0 * w)
S = sp.sqrt(A + B * r + C * r**2)

tau_expected = sp.simplify(2 * H0 * sp.sqrt(1 + r**2) / (ki + kj * r + S))

print("k_ij(r) =")
sp.pprint(kij)
print("kappa_ij(r) =")
sp.pprint(kappa)
print("A,B,C =")
sp.pprint((A, B, C))
print("tau_ij(r) =")
sp.pprint(tau)
expect_zero("explicit algebraic tau form", tau - tau_expected)
expect_zero(
    "discriminant numerator reduction",
    sp.simplify((1 + r**2) * (kij**2 - 2 * H0 * kappa) - (A + B * r + C * r**2)),
)

# ---------------------------------------------------------------------------
# II. Exact stationary numerator theorem via the denominator functional
# ---------------------------------------------------------------------------
subbanner("II. Exact stationary numerator theorem")

Phi = sp.simplify((ki + kj * r + S) / sp.sqrt(1 + r**2))
N = sp.simplify(2 * (kj - ki * r) * S + B + 2 * (C - A) * r - B * r**2)
Phi_prime_expected = sp.simplify(N / (2 * (1 + r**2) ** sp.Rational(3, 2) * S))

tau_from_Phi = sp.simplify(2 * H0 / Phi)
print("Phi_ij(r) =")
sp.pprint(Phi)
print("N_ij(r) =")
sp.pprint(N)
expect_zero("tau = 2H0 / Phi", tau - tau_from_Phi)
expect_zero("Phi derivative law", sp.diff(Phi, r) - Phi_prime_expected)

# ---------------------------------------------------------------------------
# III. Exact quartic elimination theorem
# ---------------------------------------------------------------------------
subbanner("III. Exact quartic elimination theorem")

J = sp.simplify(B + 2 * (C - A) * r - B * r**2)
Q = sp.expand(J**2 - 4 * (kj - ki * r) ** 2 * (A + B * r + C * r**2))
Qpoly = sp.Poly(Q, r)

print("J_ij(r) =")
sp.pprint(J)
print("Q_ij(r) =")
sp.pprint(Qpoly.as_expr())
print(f"degree(Q_ij) = {Qpoly.degree()}")
expect_zero("quartic degree minus 4", Qpoly.degree() - 4)
expect_zero("quartic factorization identity", Q - (J - 2 * (kj - ki * r) * S) * (J + 2 * (kj - ki * r) * S))
expect_zero("N - (J + 2(kj-kir)S)", N - (J + 2 * (kj - ki * r) * S))

# ---------------------------------------------------------------------------
# IV. Diagonal-neutral curvature reduction -> gradient-optimal ratio
# ---------------------------------------------------------------------------
subbanner("IV. Diagonal-neutral curvature reduction")

kappa_star = sp.symbols("kappa_star", real=True)
r_grad = sp.simplify(kj / ki)

tau_diag = sp.simplify(tau.subs({u: kappa_star, v: 0, w: kappa_star}))
kappa_diag = sp.simplify(kappa.subs({u: kappa_star, v: 0, w: kappa_star}))
print("tau on diagonal-neutral curvature branch =")
sp.pprint(tau_diag)
expect_zero("diagonal-neutral kappa(r) - kappa_*", kappa_diag - kappa_star)
expect_zero("gradient-optimal stationarity on diagonal-neutral branch", sp.simplify(sp.diff(tau_diag, r).subs(r, r_grad)))

# ---------------------------------------------------------------------------
# V. Pair-symmetry reduction -> equal-mix critical ray
# ---------------------------------------------------------------------------
subbanner("V. Pair-symmetry reduction")

tau_sym = sp.simplify(tau.subs({kj: ki, w: u}))
print("tau on pair-symmetric branch =")
sp.pprint(tau_sym)
expect_zero("pair-symmetry under r -> 1/r", sp.simplify(tau_sym - tau_sym.subs(r, 1 / r)))
expect_zero("equal-mix stationarity on pair-symmetric branch", sp.simplify(sp.diff(tau_sym, r).subs(r, 1)))

# ---------------------------------------------------------------------------
# VI. Explicit quartic coefficients in the original variables
# ---------------------------------------------------------------------------
subbanner("VI. Explicit quartic coefficients in the original variables")

Qorig = sp.expand(Q.subs({A: ki**2 - 2 * H0 * u, B: 2 * ki * kj - 4 * H0 * v, C: kj**2 - 2 * H0 * w}))
Qorig_poly = sp.Poly(Qorig, r)
print("Q_ij(r) in original variables =")
sp.pprint(Qorig_poly.as_expr())
print("coefficients highest -> constant =")
for coeff in Qorig_poly.all_coeffs():
    sp.pprint(sp.factor(coeff))

banner("STAGE 192 SYMPY AUDIT COMPLETED SUCCESSFULLY")
