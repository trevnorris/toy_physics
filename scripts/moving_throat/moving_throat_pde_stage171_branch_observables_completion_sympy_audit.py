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


banner("STAGE 171 — BRANCH-OBSERVABLE COMPLETION AND THE EXACT FIRST-ORDER OBSERVABLE COMPILER")

chi0s, deltaUs, epsetas = sp.symbols(
    "chi0_star deltaU_star epsiloneta_star", positive=True, real=True
)

# Coherent-branch coefficients
Ctr = sp.simplify(
    chi0s * deltaUs
    / ((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs))
)
Cstar = sp.simplify(1 / Ctr)
Bstar = sp.simplify(2 * (1 + chi0s + deltaUs) / deltaUs)
Astar = sp.simplify(2 * chi0s / ((1 + chi0s) * (1 + deltaUs)))

subbanner("I. Exact coherent-branch coefficient identities")
print("C_tr,* =")
sp.pprint(Ctr)
print("C_* =")
sp.pprint(Cstar)
print("B_* =")
sp.pprint(Bstar)
print("A_tr,* =")
sp.pprint(Astar)
expect_zero("C_* C_tr,* - 1", Cstar * Ctr - 1)
expect_zero("A_tr,* - B_* C_tr,*", Astar - Bstar * Ctr)

# First-order branch-observable packet
# dr = delta ln R_tr
# dn = delta ln N_*
# de = delta ln epsilon_eta

dr, dn, de = sp.symbols("dln_Rtr dln_Nstar dln_epseta", real=True)
obs = sp.Matrix([dr, dn, de])

subbanner("II. Observable packet -> tangent quotient packet")
C_obs_to_quot = sp.diag(-Cstar, 1, 1)
quot = sp.simplify(C_obs_to_quot * obs)
print("Delta_obs^(1) =")
sp.pprint(obs)
print("C_obs->quot =")
sp.pprint(C_obs_to_quot)
print("Delta_quot^(1) = C_obs->quot * Delta_obs^(1) =")
sp.pprint(quot)
print("det(C_obs->quot) =")
sp.pprint(sp.factor(sp.simplify(C_obs_to_quot.det())))
C_quot_to_obs = sp.simplify(C_obs_to_quot.inv())
print("C_quot->obs =")
sp.pprint(C_quot_to_obs)
expect_zero("C_quot->obs * C_obs->quot - I", C_quot_to_obs * C_obs_to_quot - sp.eye(3))
expect_zero("C_obs->quot * C_quot->obs - I", C_obs_to_quot * C_quot_to_obs - sp.eye(3))

Sigma_tr, Sigma_nt, Sigma_eta = quot

subbanner("III. Tangent quotient packet -> defect packet")
C_quot_to_def = sp.Matrix(
    [
        [-Ctr, 0, 0],
        [Astar, 1, 0],
        [-Astar, -1, -epsetas / (1 - epsetas)],
    ]
)
def_from_quot = sp.simplify(C_quot_to_def * quot)
print("C_quot->def =")
sp.pprint(C_quot_to_def)
print("Delta_def^(1) from tangent quotient packet =")
sp.pprint(def_from_quot)

Theta_from_quot = sp.simplify(def_from_quot[0])
Xi_from_quot = sp.simplify(def_from_quot[1])
R_from_quot = sp.simplify(def_from_quot[2])

expect_zero("Theta_from_quot - dln(R_tr)", Theta_from_quot - dr)
expect_zero("Xi_from_quot - (dln(N_*) - B_* dln(R_tr))", Xi_from_quot - (dn - Bstar * dr))
expect_zero(
    "R_from_quot - (-(epseta_*/(1-epseta_*)) dln(epseta) - Xi)",
    R_from_quot - (-epsetas / (1 - epsetas) * de - Xi_from_quot),
)

subbanner("IV. Direct observable packet -> defect packet")
C_obs_to_def_expected = sp.Matrix(
    [
        [1, 0, 0],
        [-Bstar, 1, 0],
        [Bstar, -1, -epsetas / (1 - epsetas)],
    ]
)
C_obs_to_def = sp.simplify(C_quot_to_def * C_obs_to_quot)
print("C_obs->def (from factorization) =")
sp.pprint(C_obs_to_def)
print("C_obs->def (expected direct compiler) =")
sp.pprint(C_obs_to_def_expected)
expect_zero("factorized compiler - expected compiler", C_obs_to_def - C_obs_to_def_expected)

Delta_def = sp.simplify(C_obs_to_def * obs)
print("Delta_def^(1) = C_obs->def * Delta_obs^(1) =")
sp.pprint(Delta_def)
print("det(C_obs->def) =")
sp.pprint(sp.factor(sp.simplify(C_obs_to_def.det())))

Theta, Xi, Rcal = Delta_def
expect_zero("Theta - dln(R_tr)", Theta - dr)
expect_zero("Xi - (dln(N_*) - B_* dln(R_tr))", Xi - (dn - Bstar * dr))
expect_zero(
    "R - (B_* dln(R_tr) - dln(N_*) - epseta_*/(1-epseta_*) dln(epseta))",
    Rcal - (Bstar * dr - dn - epsetas / (1 - epsetas) * de),
)

subbanner("V. Exact inverse observable compiler")
C_def_to_obs_expected = sp.Matrix(
    [
        [1, 0, 0],
        [Bstar, 1, 0],
        [0, -(1 - epsetas) / epsetas, -(1 - epsetas) / epsetas],
    ]
)
C_def_to_obs = sp.simplify(C_obs_to_def.inv())
print("C_def->obs =")
sp.pprint(C_def_to_obs)
print("C_def->obs (expected) =")
sp.pprint(C_def_to_obs_expected)
expect_zero("inverse compiler - expected inverse", C_def_to_obs - C_def_to_obs_expected)
expect_zero("C_def->obs * C_obs->def - I", C_def_to_obs * C_obs_to_def - sp.eye(3))
expect_zero("C_obs->def * C_def->obs - I", C_obs_to_def * C_def_to_obs - sp.eye(3))

obs_roundtrip = sp.simplify(C_def_to_obs * Delta_def)
expect_zero("observable roundtrip - Delta_obs^(1)", obs_roundtrip - obs)

subbanner("VI. Complementary selected-branch observable")
dln_E = sp.simplify(-epsetas / (1 - epsetas) * de)
print("delta ln(1 - epseta) =")
sp.pprint(dln_E)
expect_zero("(R + Xi) - delta ln(1-epseta)", (Rcal + Xi) - dln_E)

subbanner("VII. Zero-set equivalence")
# Because both compilers are invertible, zero defect <-> zero observables <-> zero quotient packet.
# Verify by substituting zero defects and recovering zero observables exactly.
zero_def = sp.Matrix([0, 0, 0])
zero_obs_from_zero_def = sp.simplify(C_def_to_obs * zero_def)
zero_quot_from_zero_obs = sp.simplify(C_obs_to_quot * zero_obs_from_zero_def)
expect_zero("zero observables from zero defect", zero_obs_from_zero_def)
expect_zero("zero quotient packet from zero observables", zero_quot_from_zero_obs)

banner("STAGE 171 LEDGER")
print("1. The Stage-167 branch-observable packet")
print("      (delta ln R_tr, delta ln N_*, delta ln epseta)")
print("   is an exact first-order compiler image of the Stage-170 tangent quotient packet.")
print("2. The same observable packet compiles exactly to")
print("      (Theta_1, Xi_1, R_1).")
print("3. The coefficient identity A_tr,* = B_* C_tr,* is the algebraic hinge that makes")
print("   the nontracking composite subtract the universal tracking feed-through exactly.")
print("4. Therefore the first-order zero-defect test is equivalently the zero-set of")
print("      delta ln R_tr, delta ln N_*, delta ln epseta.")
