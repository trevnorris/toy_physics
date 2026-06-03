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


banner("STAGE 208 — PAIRWISE MIXED RAYS AND OFF-DIAGONAL HESSIAN SYNERGY")

# ---------------------------------------------------------------------------
# I. Exact pairwise mixed-ray family and oriented slope law
# ---------------------------------------------------------------------------
subbanner("I. Exact pairwise mixed-ray family and oriented slope law")

ki, kj, r = sp.symbols("k_i k_j r", positive=True, real=True)

s = sp.Matrix([1, r]) / sp.sqrt(1 + r**2)
g = sp.Matrix([-ki, -kj])
Kij = sp.simplify((g.T * s)[0])
kij = sp.simplify(-Kij)

print("mixed-ray direction s_ij(r) =")
sp.pprint(s)
print("oriented slope K_ij(r) =")
sp.pprint(Kij)
print("positive slope magnitude k_ij(r) =")
sp.pprint(kij)
expect_zero(
    "mixed slope law",
    Kij + (ki + r * kj) / sp.sqrt(1 + r**2),
)

kij_prime = sp.simplify(sp.diff(kij, r))
print("d k_ij / dr =")
sp.pprint(kij_prime)
expect_zero(
    "mixed slope derivative law",
    kij_prime - (kj - ki * r) / (1 + r**2) ** sp.Rational(3, 2),
)

r_grad = sp.simplify(kj / ki)
kij_grad = sp.simplify(kij.subs(r, r_grad))
print("gradient-optimal ratio r_grad =")
sp.pprint(r_grad)
print("k_ij(r_grad) =")
sp.pprint(kij_grad)
expect_zero("gradient-optimal ratio stationarity", kij_prime.subs(r, r_grad))
expect_zero("gradient-optimal slope", kij_grad - sp.sqrt(ki**2 + kj**2))
expect_zero("gradient gain over primitive i", kij_grad**2 - ki**2 - kj**2)

# ---------------------------------------------------------------------------
# II. Exact mixed curvature decomposition and cross-Hessian leverage
# ---------------------------------------------------------------------------
subbanner("II. Exact mixed curvature decomposition and cross-Hessian leverage")

hii, hij, hjj = sp.symbols("h_ii h_ij h_jj", real=True)
Hpair = sp.Matrix([[hii, hij], [hij, hjj]])
H1ij = sp.simplify((s.T * Hpair * s)[0])

w_i = sp.simplify(1 / (1 + r**2))
w_j = sp.simplify(r**2 / (1 + r**2))
w_x = sp.simplify(2 * r / (1 + r**2))
H1_expected = sp.simplify(w_i * hii + w_x * hij + w_j * hjj)

print("pairwise curvature H_1,ij(r) =")
sp.pprint(H1ij)
print("weights (w_i, w_x, w_j) =")
sp.pprint((w_i, w_x, w_j))
expect_zero("mixed curvature decomposition", H1ij - H1_expected)
expect_zero("diagonal neutrality when h_ij = 0", H1ij.subs(hij, 0) - (w_i * hii + w_j * hjj))

w_x_prime = sp.simplify(sp.diff(w_x, r))
print("d w_x / dr =")
sp.pprint(w_x_prime)
expect_zero("cross-weight derivative", w_x_prime - 2 * (1 - r**2) / (1 + r**2) ** 2)
expect_zero("cross-weight stationarity at r=1", w_x_prime.subs(r, 1))
expect_zero("cross-weight maximum value", w_x.subs(r, 1) - 1)

# ---------------------------------------------------------------------------
# III. Canonical gradient-optimal and equal-mix mixed rays
# ---------------------------------------------------------------------------
subbanner("III. Canonical gradient-optimal and equal-mix mixed rays")

# Gradient-optimal ray
s_grad = sp.simplify(sp.Matrix([ki, kj]) / sp.sqrt(ki**2 + kj**2))
K_grad = sp.simplify((g.T * s_grad)[0])
H_grad = sp.simplify((s_grad.T * Hpair * s_grad)[0])

print("gradient-optimal mixed direction =")
sp.pprint(s_grad)
print("gradient-optimal slope =")
sp.pprint(K_grad)
print("gradient-optimal curvature =")
sp.pprint(H_grad)
expect_zero("gradient-optimal direction from r_grad", s.subs(r, r_grad) - s_grad)
expect_zero("gradient-optimal K", K_grad + sp.sqrt(ki**2 + kj**2))
expect_zero(
    "gradient-optimal curvature law",
    H_grad - (ki**2 * hii + 2 * ki * kj * hij + kj**2 * hjj) / (ki**2 + kj**2),
)

# Equal-mix synergy ray
s_eq = sp.Matrix([1, 1]) / sp.sqrt(2)
K_eq = sp.simplify((g.T * s_eq)[0])
H_eq = sp.simplify((s_eq.T * Hpair * s_eq)[0])

print("equal-mix direction =")
sp.pprint(s_eq)
print("equal-mix slope =")
sp.pprint(K_eq)
print("equal-mix curvature =")
sp.pprint(H_eq)
expect_zero("equal-mix direction from r=1", s.subs(r, 1) - s_eq)
expect_zero("equal-mix slope law", K_eq + (ki + kj) / sp.sqrt(2))
expect_zero("equal-mix curvature law", H_eq - (hii + 2 * hij + hjj) / 2)

# ---------------------------------------------------------------------------
# IV. Exact mixed curvature envelopes from entrywise envelopes
# ---------------------------------------------------------------------------
subbanner("IV. Exact mixed curvature envelopes from entrywise envelopes")

hii_lo, hij_lo, hjj_lo = sp.symbols("hii_lo hij_lo hjj_lo", real=True)
hii_hi, hij_hi, hjj_hi = sp.symbols("hii_hi hij_hi hjj_hi", real=True)

kappa_lo = sp.simplify((hii_lo + 2 * r * hij_lo + r**2 * hjj_lo) / (1 + r**2))
kappa_hi = sp.simplify((hii_hi + 2 * r * hij_hi + r**2 * hjj_hi) / (1 + r**2))

print("lower mixed curvature envelope =")
sp.pprint(kappa_lo)
print("upper mixed curvature envelope =")
sp.pprint(kappa_hi)
expect_zero(
    "lower envelope weighted form",
    kappa_lo - (w_i * hii_lo + w_x * hij_lo + w_j * hjj_lo),
)
expect_zero(
    "upper envelope weighted form",
    kappa_hi - (w_i * hii_hi + w_x * hij_hi + w_j * hjj_hi),
)

# ---------------------------------------------------------------------------
# V. Fixed-ratio mixed-ray certified bracket data
# ---------------------------------------------------------------------------
subbanner("V. Fixed-ratio mixed-ray certified bracket data")

H0 = sp.symbols("H0", positive=True, real=True)
Delta_lo = sp.simplify(kij**2 - 2 * kappa_lo * H0)
Delta_hi = sp.simplify(kij**2 - 2 * kappa_hi * H0)

T_lo = sp.simplify(2 * H0 / (kij + sp.sqrt(Delta_lo)))
T_hi = sp.simplify(2 * H0 / (kij + sp.sqrt(Delta_hi)))

print("Delta_lo(r) =")
sp.pprint(Delta_lo)
print("Delta_hi(r) =")
sp.pprint(Delta_hi)
print("tau_lo(r) =")
sp.pprint(T_lo)
print("tau_hi(r) =")
sp.pprint(T_hi)
expect_zero(
    "quadratic root relation for tau_lo",
    sp.simplify(H0 - kij * T_lo + sp.Rational(1, 2) * kappa_lo * T_lo**2),
)
expect_zero(
    "quadratic root relation for tau_hi",
    sp.simplify(H0 - kij * T_hi + sp.Rational(1, 2) * kappa_hi * T_hi**2),
)

# ---------------------------------------------------------------------------
# VI. Canonical screen rays separate gradient gain from cross-Hessian leverage
# ---------------------------------------------------------------------------
subbanner("VI. Canonical screen rays separate gradient gain from cross-Hessian leverage")

expect_zero(
    "gradient and equal-mix coincide iff k_i-k_j = 0 via ratio difference",
    sp.simplify(r_grad - 1) - (kj - ki) / ki,
)

w_x_grad = sp.simplify(w_x.subs(r, r_grad))
print("cross-Hessian weight on gradient-optimal ray =")
sp.pprint(w_x_grad)
print("cross-Hessian weight on equal-mix ray =")
sp.pprint(w_x.subs(r, 1))

# A useful exact comparison identity
expect_zero(
    "gradient-optimal cross weight formula",
    w_x_grad - 2 * ki * kj / (ki**2 + kj**2),
)

banner("STAGE 208 SYMPY AUDIT COMPLETED SUCCESSFULLY")
