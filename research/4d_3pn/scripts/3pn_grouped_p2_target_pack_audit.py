#!/usr/bin/env python3
"""
3pn_grouped_p2_target_pack_audit.py

Exact audit for the 3PN grouped-P2 target pack.

What this script does
---------------------
1. Reconstructs the exact 3PN COM residual beyond the frozen one-body/self-static seed.
2. Verifies the grouped-P2 axisymmetric inverse map
       (u2^(20),u2^(21),u2^(22)) <-> (ubar2,a2,b2).
3. Verifies the exact co-rotating pair-frame source coefficients for the grouped real P2 bundle.
4. Builds the first exact time-local O(omega^2) grouped front-end ansatz and prints the
   isotropic collapse.

This script is intentionally a target-pack / front-end audit rather than a full throat-side
3PN derivation.  Its job is to freeze the exact data vector and the first explicit grouped-P2
kinematic scaffold that the next constitutive matching step has to use.
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


def expect_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(sp.expand(expr))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Part I — Exact GR 3PN COM residual data pack
# ---------------------------------------------------------------------------

banner("PART I — EXACT 3PN COM RESIDUAL DATA PACK")

nu = sp.symbols("nu", real=True)

Delta_l = {
    1: 3 * nu * (3 * nu - 1) * (4 * nu - 1) / 16,
    2: sp.Integer(0),
    3: sp.Integer(0),
    4: sp.Integer(0),
    5: sp.Integer(0),
    6: nu * (38 - 116 * nu - 57 * nu**2) / 16,
    7: nu**2 * (20 - 69 * nu) / 16,
    8: 3 * nu**2 * (3 - 11 * nu) / 16,
    9: 5 * nu**3 / 16,
    10: nu * (129 - 98 * nu + 52 * nu**2) / 16,
    11: nu * (-3 + 52 * nu + 124 * nu**2) / 16,
    12: nu * (-5 + 11 * nu + 48 * nu**2) / 12,
    13: -nu * (244 + 3 * sp.pi**2 + 1272 * nu + 96 * nu**2) / 192,
    14: nu * (452 + 3 * sp.pi**2 - 384 * nu - 224 * nu**2) / 64,
    15: nu * (-908 + 63 * sp.pi**2) / 96,
}

nonzero_slots = [i for i, expr in Delta_l.items() if sp.simplify(expr) != 0]
print("Nonzero residual slots:", nonzero_slots)
print("Count =", len(nonzero_slots))
if len(nonzero_slots) != 11:
    raise AssertionError("Unexpected number of nonzero COM residual slots.")

nu_eq = sp.Rational(1, 4)
print("\nEqual-mass specialization nu = 1/4:")
for i in nonzero_slots:
    print(f"  Delta l_{i} =", sp.simplify(Delta_l[i].subs(nu, nu_eq)))


# ---------------------------------------------------------------------------
# Part II — Grouped-P2 axisymmetric inverse map
# ---------------------------------------------------------------------------

banner("PART II — GROUPED-P2 AXISYMMETRIC INVERSE MAP")

u2_20, u2_21, u2_22 = sp.symbols("u2_20 u2_21 u2_22", real=True)
ubar2, a2, b2 = sp.symbols("ubar2 a2 b2", real=True)

ubar2_expr = (u2_20 + 2 * u2_21 + 2 * u2_22) / 5
a2_expr = (2 * u2_20 - u2_21 - u2_22) / 10
b2_expr = (u2_21 - u2_22) / 2

u2_20_back = ubar2 + 4 * a2
u2_21_back = ubar2 - a2 + b2
u2_22_back = ubar2 - a2 - b2

expect_zero("u2_20 recovered", u2_20_back.subs({ubar2: ubar2_expr, a2: a2_expr, b2: b2_expr}) - u2_20)
expect_zero("u2_21 recovered", u2_21_back.subs({ubar2: ubar2_expr, a2: a2_expr, b2: b2_expr}) - u2_21)
expect_zero("u2_22 recovered", u2_22_back.subs({ubar2: ubar2_expr, a2: a2_expr, b2: b2_expr}) - u2_22)

print("\nParameter count:")
print("  full grouped O(w^2) payload   : 3 real numbers")
print("  exact COM residual data vector : 11 nonzero slots")
print("  overdetermination              : 8 equations")
print("  isotropic branch (a2=b2=0)     : 1 datum, 10 checks")


# ---------------------------------------------------------------------------
# Part III — Co-rotating grouped source coefficients
# ---------------------------------------------------------------------------

banner("PART III — CO-ROTATING GROUPED SOURCE COEFFICIENTS")

ux, uy, d, r, U = sp.symbols("ux uy d r U", real=True)
u2 = ux**2 + uy**2

C20 = sp.sqrt(sp.Rational(2, 3)) * (d**2 - u2 / 2 - U)
C21c = sp.sqrt(2) * d * ux
C21s = sp.sqrt(2) * d * uy
C22c = (ux**2 - uy**2) / sp.sqrt(2)
C22s = sp.sqrt(2) * ux * uy

# Co-rotating pair-frame kinematics.
ddot = (u2 - U) / r
uxdot = -d * ux / r
uydot = -d * uy / r
Udot = -d * U / r


def total_derivative(expr: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(expr, d) * ddot
        + sp.diff(expr, ux) * uxdot
        + sp.diff(expr, uy) * uydot
        + sp.diff(expr, U) * Udot
    )

C20_dot = sp.simplify(total_derivative(C20))
C21c_dot = sp.simplify(total_derivative(C21c))
C21s_dot = sp.simplify(total_derivative(C21s))
C22c_dot = sp.simplify(total_derivative(C22c))
C22s_dot = sp.simplify(total_derivative(C22s))

expect_zero(
    "C20_dot formula",
    C20_dot - sp.sqrt(sp.Rational(2, 3)) * d * (3 * u2 - U) / r,
)
expect_zero(
    "C21c_dot formula",
    C21c_dot - sp.sqrt(2) * ux * (u2 - d**2 - U) / r,
)
expect_zero(
    "C21s_dot formula",
    C21s_dot - sp.sqrt(2) * uy * (u2 - d**2 - U) / r,
)
expect_zero(
    "C22c_dot formula",
    C22c_dot + sp.sqrt(2) * d * (ux**2 - uy**2) / r,
)
expect_zero(
    "C22s_dot formula",
    C22s_dot + 2 * sp.sqrt(2) * d * ux * uy / r,
)

A20 = sp.simplify(C20_dot**2)
A21 = sp.simplify(C21c_dot**2 + C21s_dot**2)
A22 = sp.simplify(C22c_dot**2 + C22s_dot**2)

expect_zero("A20 grouped norm", A20 - 2 * d**2 * (3 * u2 - U)**2 / (3 * r**2))
expect_zero("A21 grouped norm", A21 - 2 * u2 * (u2 - d**2 - U)**2 / r**2)
expect_zero("A22 grouped norm", A22 - 2 * d**2 * u2**2 / r**2)

print("\nExact grouped source norms:")
print("  A20 =", A20)
print("  A21 =", sp.factor(A21))
print("  A22 =", sp.factor(A22))


# ---------------------------------------------------------------------------
# Part IV — First exact time-local O(w^2) front-end scaffold
# ---------------------------------------------------------------------------

banner("PART IV — FIRST EXACT TIME-LOCAL O(w^2) FRONT-END SCAFFOLD")

c2_20, c2_21, c2_22 = sp.symbols("c2_20 c2_21 c2_22", real=True)

L_front = sp.simplify(sp.expand(sp.Rational(1, 2) * (c2_20 * A20 + c2_21 * A21 + c2_22 * A22)))
print("L_front =")
sp.pprint(L_front)

expect_zero(
    "compact grouped front-end form",
    L_front
    - (
        c2_20 * d**2 * (3 * u2 - U)**2 / (3 * r**2)
        + c2_21 * u2 * (u2 - d**2 - U)**2 / r**2
        + c2_22 * d**2 * u2**2 / r**2
    ),
)

# Axisymmetric variables.
ubar_c2, a_c2, b_c2 = sp.symbols("ubar_c2 a_c2 b_c2", real=True)
L_axis = sp.expand(
    L_front.subs(
        {
            c2_20: ubar_c2 + 4 * a_c2,
            c2_21: ubar_c2 - a_c2 + b_c2,
            c2_22: ubar_c2 - a_c2 - b_c2,
        }
    )
)
print("\nL_front in (ubar,a,b) variables =")
sp.pprint(L_axis)

u2_sym, v2, d2 = sp.symbols("u2_sym v2 d2", real=True)
ciso = sp.Symbol("ciso")
L_iso_u = sp.simplify(
    ciso * (
        d2 * (3 * u2_sym - U)**2 / (3 * r**2)
        + u2_sym * (u2_sym - d2 - U)**2 / r**2
        + d2 * u2_sym**2 / r**2
    )
)
L_iso = sp.simplify(sp.expand(L_iso_u.subs(u2_sym, v2 - d2)))
L_iso_expected = ciso * (
    3 * v2**3
    - 3 * d2 * v2**2
    - 6 * U * v2**2
    + 12 * U * d2 * v2
    - 6 * U * d2**2
    + 3 * U**2 * v2
    - 2 * U**2 * d2
) / (3 * r**2)
expect_zero("isotropic collapse", L_iso - L_iso_expected)

print("\nIsotropic front-end scaffold:")
print("  L_iso =", L_iso_expected)


# ---------------------------------------------------------------------------
# Part V — Final ledger
# ---------------------------------------------------------------------------

banner("PART V — FINAL LEDGER")
print("1. The exact COM residual beyond the frozen seed has 11 nonzero slots.")
print("2. The grouped P2 O(w^2) payload has exactly 3 real data: (u2^(20),u2^(21),u2^(22)).")
print("3. Therefore the grouped-P2 3PN inverse problem is 8-fold overdetermined before any")
print("   isotropy assumption, and 10-fold overdetermined on the minimal isotropic branch.")
print("4. The co-rotating grouped source coefficients and their exact derivative norms are now")
print("   explicit, giving the first exact time-local grouped-P2 front-end scaffold.")
print("5. The remaining live task is to compute the actual throat-side dictionary from this")
print("   grouped front end into the solved 3PN ordinary/Hamiltonian target chart.")
