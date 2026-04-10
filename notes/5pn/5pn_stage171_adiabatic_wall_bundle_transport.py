#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = '=' * 88
    print('\n' + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = '-' * 88
    print('\n' + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        simp = expr.applyfunc(lambda e: sp.simplify(sp.expand(e)))
        print(f"{name} =")
        sp.pprint(simp)
        if any(e != 0 for e in simp):
            raise AssertionError(f"{name} is not zero")
    else:
        simp = sp.simplify(sp.expand(expr))
        print(f"{name} = {simp}")
        if simp != 0:
            raise AssertionError(f"{name} is not zero")


"""
5pn_stage171_adiabatic_wall_bundle_transport.py

Stage 171 — impose the adiabatic/frozen-wall condition
    delta ln Theta_w = 0
on the Stage-149/150 isotropic bundle transport.

What this script does
---------------------
1. Replays the exact Stage-149 inversion formulas.
2. Imposes the adiabatic-wall constraint dlnTheta = 0.
3. Derives the simplified first-order transport laws for
      rho_w, c_{s,w}, ell, a, L_W, c_s, Z_q, T_m, v_{w0}, g_s, g_q, lambda.
4. Proves that the parent compensation invariants
      r_c, mathfrak r, mathfrak g
   still have zero first-order drift on the adiabatic branch.

Interpretation
--------------
The adiabatic wall removes wall-depth / thermal-fraying drift as an isotropic source,
but it does not by itself solve the weak-axisymmetric quotient-orbit problem.
"""

banner("STAGE 171 — ADIABATIC WALL BUNDLE TRANSPORT")

# Bundle observables

dlnTheta, dlnKs, dlnKq, dlnP0 = sp.symbols('dlnTheta dlnKs dlnKq dlnP0', real=True)

# Remaining microscopic drifts

dlnrhow, dlncsw, dlnell = sp.symbols('dlnrhow dlncsw dlnell', real=True)
dlna, dlnLW, dlncs, dlnZq = sp.symbols('dlna dlnLW dlncs dlnZq', real=True)
dlnTm, dlnvw0 = sp.symbols('dlnTm dlnvw0', real=True)
dlngs, dlngq, dlnlambda = sp.symbols('dlngs dlngq dlnlambda', real=True)

dlnrc, dlnr, dlng = sp.symbols('dlnrc dlnr dlng', real=True)

subbanner("1. Exact Stage-149 inversion formulas")

subs149 = {
    dlnrhow: sp.Rational(1, 2) * dlnTheta,
    dlna: sp.Rational(1, 2) * dlnKs - sp.Rational(1, 4) * dlnTheta,
    dlncs: sp.Rational(1, 2) * dlnKs - sp.Rational(1, 4) * dlnTheta + sp.Rational(1, 5) * dlnP0,
    dlnZq: dlnKq - sp.Rational(2, 5) * dlnP0,
}

dlncsw_expr = dlnTheta
# frozen n=5 wall-EOS branch

dlnell_expr = -dlnTheta
# healing-locked shell / explicit Family-1 branch

dlnLW_expr = subs149[dlna]

print("delta ln rho_w =", subs149[dlnrhow])
print("delta ln c_{s,w} =", dlncsw_expr)
print("delta ln ell =", dlnell_expr)
print("delta ln a =", subs149[dlna])
print("delta ln L_W =", dlnLW_expr)
print("delta ln c_s =", subs149[dlncs])
print("delta ln Z_q =", subs149[dlnZq])

subbanner("2. Adiabatic wall condition: delta ln Theta_w = 0")

adiabatic = {dlnTheta: 0}

expect_zero("delta ln rho_w on adiabatic wall", subs149[dlnrhow].subs(adiabatic))
expect_zero("delta ln c_{s,w} on adiabatic wall", dlncsw_expr.subs(adiabatic))
expect_zero("delta ln ell on adiabatic wall", dlnell_expr.subs(adiabatic))

print("delta ln a |ad =", sp.simplify(subs149[dlna].subs(adiabatic)))
print("delta ln L_W |ad =", sp.simplify(dlnLW_expr.subs(adiabatic)))
print("delta ln c_s |ad =", sp.simplify(subs149[dlncs].subs(adiabatic)))
print("delta ln Z_q |ad =", sp.simplify(subs149[dlnZq].subs(adiabatic)))

subbanner("3. Exact adiabatic transport of v_{w0} and T_m")

# Stage-148 transport laws
vw0_expr = sp.simplify(
    (sp.Rational(1,2) * (dlnZq - dlnrhow) + sp.Rational(3,2) * dlncsw + dlncs - sp.Rational(5,2) * dlna)
    .subs(subs149)
    .subs({dlncsw: dlncsw_expr})
    .subs(adiabatic)
)
Tm_expr = sp.simplify(
    (sp.Rational(1,2) * (dlnZq - dlnrhow) + sp.Rational(3,2) * dlncsw - dlncs - sp.Rational(3,2) * dlna)
    .subs(subs149)
    .subs({dlncsw: dlncsw_expr})
    .subs(adiabatic)
)

expect_zero(
    "delta ln v_{w0} + 3/4 dlnK_s - 1/2 dlnK_q",
    vw0_expr - (-sp.Rational(3,4) * dlnKs + sp.Rational(1,2) * dlnKq),
)
expect_zero(
    "delta ln T_m + 5/4 dlnK_s - 1/2 dlnK_q + 2/5 dlnP0",
    Tm_expr - (-sp.Rational(5,4) * dlnKs + sp.Rational(1,2) * dlnKq - sp.Rational(2,5) * dlnP0),
)
print("delta ln v_{w0} |ad =", vw0_expr)
print("delta ln T_m |ad =", Tm_expr)

ratio_expr = sp.simplify(vw0_expr - Tm_expr)
product_expr = sp.simplify(vw0_expr + Tm_expr)
expect_zero("delta ln(v_{w0}/T_m) - (1/2 dlnK_s + 2/5 dlnP0)", ratio_expr - (sp.Rational(1,2) * dlnKs + sp.Rational(2,5) * dlnP0))
expect_zero("delta ln(v_{w0} T_m) + 2 dlnK_s - dlnK_q + 2/5 dlnP0", product_expr - (-2 * dlnKs + dlnKq - sp.Rational(2,5) * dlnP0))
print("delta ln(v_{w0}/T_m) |ad =", ratio_expr)
print("delta ln(v_{w0} T_m) |ad =", product_expr)

subbanner("4. Exact adiabatic transport of g_s, g_q, lambda")

gs_expr = sp.simplify(Tm_expr + 2 * subs149[dlna].subs(adiabatic) + dlnell_expr.subs(adiabatic))
gq_expr = sp.simplify(subs149[dlnZq].subs(adiabatic) - sp.Rational(3,2) * dlnLW_expr.subs(adiabatic))
lambda_expr = sp.simplify(vw0_expr + 2 * subs149[dlna].subs(adiabatic) + dlnell_expr.subs(adiabatic) + sp.Rational(1,2) * dlnLW_expr.subs(adiabatic))

expect_zero(
    "delta ln g_s + 1/4 dlnK_s - 1/2 dlnK_q + 2/5 dlnP0",
    gs_expr - (-sp.Rational(1,4) * dlnKs + sp.Rational(1,2) * dlnKq - sp.Rational(2,5) * dlnP0),
)
expect_zero(
    "delta ln g_q + 3/4 dlnK_s - dlnK_q + 2/5 dlnP0",
    gq_expr - (-sp.Rational(3,4) * dlnKs + dlnKq - sp.Rational(2,5) * dlnP0),
)
expect_zero(
    "delta ln lambda - 1/2(dlnK_s + dlnK_q)",
    lambda_expr - sp.Rational(1,2) * (dlnKs + dlnKq),
)
print("delta ln g_s |ad =", gs_expr)
print("delta ln g_q |ad =", gq_expr)
print("delta ln lambda |ad =", lambda_expr)

subbanner("5. Parent compensation invariants remain tangent")

# Parent compensation invariants

dlnrc_expr = sp.simplify(2 * lambda_expr - dlnKs - dlnKq)
dlnr_expr = sp.simplify(lambda_expr - sp.Rational(1,2) * (dlnKs + dlnKq))
dlng_expr = sp.simplify(gq_expr + sp.Rational(1,2) * dlnKs - gs_expr - sp.Rational(1,2) * dlnKq)

expect_zero("delta ln r_c |ad", dlnrc_expr)
expect_zero("delta ln mathfrak r |ad", dlnr_expr)
expect_zero("delta ln mathfrak g |ad", dlng_expr)

banner("FINAL STAGE-171 LEDGER")
print("1. The adiabatic wall condition sets delta ln Theta_w = 0 and therefore")
print("   freezes rho_w, c_{s,w}, and ell at first order.")
print("2. The remaining isotropic bundle drift family is three-dimensional:")
print("      (delta ln K_s, delta ln K_q, delta ln P_0).")
print("3. On this adiabatic family, the exact transport laws are")
print("      delta ln a = delta ln L_W = 1/2 delta ln K_s,")
print("      delta ln c_s = 1/2 delta ln K_s + 1/5 delta ln P_0,")
print("      delta ln Z_q = delta ln K_q - 2/5 delta ln P_0,")
print("      delta ln v_{w0} = -3/4 delta ln K_s + 1/2 delta ln K_q,")
print("      delta ln T_m  = -5/4 delta ln K_s + 1/2 delta ln K_q - 2/5 delta ln P_0.")
print("4. The parent compensation invariants r_c, mathfrak r, mathfrak g still have zero first-order drift.")
print("5. So the adiabatic wall removes isotropic wall-depth / thermal-fraying motion, but it does not by itself decide the weak-axisymmetric quotient-orbit question.")
