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


banner("STAGE 194 — FULL INTERIOR SIMPLEX OPTIMIZER AND FINITE ALGEBRAIC CANDIDATE REDUCTION")

# ---------------------------------------------------------------------------
# I. Exact interior ratio objective and stationary numerators
# ---------------------------------------------------------------------------
subbanner("I. Exact interior ratio objective and stationary numerators")

r, s = sp.symbols("r s", nonnegative=True, real=True)
ki, kj, kk = sp.symbols("k_i k_j k_k", positive=True, real=True)
H0 = sp.symbols("H0", positive=True, real=True)

A, B, C, D, E, F = sp.symbols("A B C D E F", real=True)
Delta = A + B * r + C * s + D * r**2 + E * r * s + F * s**2
sqrtDelta = sp.sqrt(Delta)

Phi = (ki + kj * r + kk * s + sqrtDelta) / sp.sqrt(1 + r**2 + s**2)
tau = sp.simplify(2 * H0 / Phi)

Mr = sp.simplify((1 + r**2 + s**2) * kj - r * (ki + kj * r + kk * s))
Ms = sp.simplify((1 + r**2 + s**2) * kk - s * (ki + kj * r + kk * s))
Lr = sp.simplify((1 + r**2 + s**2) * sp.diff(Delta, r) - 2 * r * Delta)
Ls = sp.simplify((1 + r**2 + s**2) * sp.diff(Delta, s) - 2 * s * Delta)
Nr = sp.simplify(2 * Mr * sqrtDelta + Lr)
Ns = sp.simplify(2 * Ms * sqrtDelta + Ls)

print("Phi(r,s) =")
sp.pprint(Phi)
print("tau(r,s) =")
sp.pprint(tau)
print("M_r =")
sp.pprint(Mr)
print("M_s =")
sp.pprint(Ms)
print("L_r =")
sp.pprint(Lr)
print("L_s =")
sp.pprint(Ls)

expected_dPhi_dr = Nr / (2 * (1 + r**2 + s**2) ** sp.Rational(3, 2) * sqrtDelta)
expected_dPhi_ds = Ns / (2 * (1 + r**2 + s**2) ** sp.Rational(3, 2) * sqrtDelta)
expect_zero("exact dPhi/dr compiler", sp.diff(Phi, r) - expected_dPhi_dr)
expect_zero("exact dPhi/ds compiler", sp.diff(Phi, s) - expected_dPhi_ds)

# ---------------------------------------------------------------------------
# II. Exact elimination system and degree counts
# ---------------------------------------------------------------------------
subbanner("II. Exact quartic-sextic elimination system and degree counts")

Ccross = sp.expand(Ms * Lr - Mr * Ls)
Sr = sp.expand(Lr**2 - 4 * Mr**2 * Delta)
Ss = sp.expand(Ls**2 - 4 * Ms**2 * Delta)

print("C_cross(r,s) =")
sp.pprint(Ccross)
print("S_r(r,s) =")
sp.pprint(Sr)
print("S_s(r,s) =")
sp.pprint(Ss)

expect_zero("cross-elimination identity", Ms * Nr - Mr * Ns - Ccross)
expect_zero("square elimination identity (r)", Nr * (Nr - 4 * Mr * sqrtDelta) - Sr)
expect_zero("square elimination identity (s)", Ns * (Ns - 4 * Ms * sqrtDelta) - Ss)

poly_C = sp.Poly(Ccross, r, s)
poly_Sr = sp.Poly(Sr, r, s)
poly_Ss = sp.Poly(Ss, r, s)

print(f"deg C_cross = {poly_C.total_degree()}")
print(f"deg S_r = {poly_Sr.total_degree()}")
print(f"deg S_s = {poly_Ss.total_degree()}")

if poly_C.total_degree() != 4:
    raise AssertionError("C_cross is not quartic")
if poly_Sr.total_degree() != 6:
    raise AssertionError("S_r is not sextic")
if poly_Ss.total_degree() != 6:
    raise AssertionError("S_s is not sextic")

# ---------------------------------------------------------------------------
# III. Diagonal-isotropic curvature reduction
# ---------------------------------------------------------------------------
subbanner("III. Diagonal-isotropic curvature reduction")

u = sp.symbols("u", real=True)
A_iso = ki**2 - 2 * H0 * u
B_iso = 2 * ki * kj
C_iso = 2 * ki * kk
D_iso = kj**2 - 2 * H0 * u
E_iso = 2 * kj * kk
F_iso = kk**2 - 2 * H0 * u
Delta_iso = sp.simplify(Delta.subs({A: A_iso, B: B_iso, C: C_iso, D: D_iso, E: E_iso, F: F_iso}))

k_rs = sp.simplify((ki + kj * r + kk * s) / sp.sqrt(1 + r**2 + s**2))
tau_iso = sp.simplify(tau.subs({A: A_iso, B: B_iso, C: C_iso, D: D_iso, E: E_iso, F: F_iso}))

a_iso_expr = (ki + kj * r + kk * s) ** 2 - 2 * H0 * u * (1 + r**2 + s**2)
expect_zero("Delta_iso reduction", Delta_iso - a_iso_expr)
expect_zero(
    "tau_iso depends only on k_rs",
    tau_iso - 2 * H0 / (k_rs + sp.sqrt(k_rs**2 - 2 * H0 * u)),
)

# ---------------------------------------------------------------------------
# IV. Full triple-symmetry reduction and the equal-mix barycenter
# ---------------------------------------------------------------------------
subbanner("IV. Full triple-symmetry reduction and the equal-mix barycenter")

k = sp.symbols("k", positive=True, real=True)
ud, ux = sp.symbols("u_d u_x", real=True)

sym_subs = {
    ki: k,
    kj: k,
    kk: k,
    A: k**2 - 2 * H0 * ud,
    B: 2 * k**2 - 4 * H0 * ux,
    C: 2 * k**2 - 4 * H0 * ux,
    D: k**2 - 2 * H0 * ud,
    E: 2 * k**2 - 4 * H0 * ux,
    F: k**2 - 2 * H0 * ud,
}

Nr_sym = sp.simplify(Nr.subs(sym_subs).subs({r: 1, s: 1}))
Ns_sym = sp.simplify(Ns.subs(sym_subs).subs({r: 1, s: 1}))
expect_zero("equal-mix stationary numerator Nr(1,1)", Nr_sym)
expect_zero("equal-mix stationary numerator Ns(1,1)", Ns_sym)

# ---------------------------------------------------------------------------
# V. Candidate-count bound statement
# ---------------------------------------------------------------------------
subbanner("V. Algebraic candidate-count bound")

bezout_bound = poly_C.total_degree() * poly_Sr.total_degree()
print(f"Bezout bound from (quartic, sextic) system = {bezout_bound}")
if bezout_bound != 24:
    raise AssertionError("Unexpected Bezout bound")

print("\nAll Stage 211 identities verified.")
