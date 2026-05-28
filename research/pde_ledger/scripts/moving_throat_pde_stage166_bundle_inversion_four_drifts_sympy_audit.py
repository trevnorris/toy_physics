#!/usr/bin/env python3
"""
moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py

SymPy audit for Stage 166.

Checks:
1. Solve the exact logarithmic inversion system
      Theta_w, K_s, K_q, P_0  ->  rho_w, a, c_s, Z_q.
2. Verify the forward transport identities after substitution.
3. Verify the equivalent bundle form using P_0 = N_0 / D_0.
4. Report the explicit Family-1 frozen-wall simplification and rho_w^(chi).
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


banner("STAGE 166 — EXACT BUNDLE INVERSION OF THE LAST FOUR DRIFTS")

dTheta, dKs, dKq, dP = sp.symbols("dTheta dKs dKq dP", real=True)
drho, da, dcs, dZ = sp.symbols("drho da dcs dZ", real=True)
dN0, dD0 = sp.symbols("dN0 dD0", real=True)

# Exact logarithmic branch laws used in the note.
eq1 = sp.Eq(dTheta, 2 * drho)
eq2 = sp.Eq(dKs, 2 * da + drho)
eq3 = sp.Eq(dKq, dZ + 2 * dcs - 2 * da)
eq4 = sp.Eq(dP, 5 * (dcs - da))

sol = sp.solve((eq1, eq2, eq3, eq4), (drho, da, dcs, dZ), dict=True)[0]

print("drho =", sp.simplify(sol[drho]))
print("da   =", sp.simplify(sol[da]))
print("dcs  =", sp.simplify(sol[dcs]))
print("dZ   =", sp.simplify(sol[dZ]))

banner("General inversion forms (paper Sec. 2)")
expect_zero("drho general", sol[drho] - sp.Rational(1, 2) * dTheta)
expect_zero("da general", sol[da] - (sp.Rational(1, 2) * dKs - sp.Rational(1, 4) * dTheta))
expect_zero("dcs general", sol[dcs] - (sp.Rational(1, 2) * dKs - sp.Rational(1, 4) * dTheta + sp.Rational(1, 5) * dP))
expect_zero("dZ general", sol[dZ] - (dKq - sp.Rational(2, 5) * dP))

banner("Forward verification")
expect_zero("Theta law", eq1.lhs.subs(sol) - eq1.rhs.subs(sol))
expect_zero("Ks law", eq2.lhs.subs(sol) - eq2.rhs.subs(sol))
expect_zero("Kq law", eq3.lhs.subs(sol) - eq3.rhs.subs(sol))
expect_zero("P0 law", eq4.lhs.subs(sol) - eq4.rhs.subs(sol))

banner("Equivalent full-bundle form with P_0 = N_0 / D_0")
sol_bundle = {
    drho: sol[drho],
    da: sol[da],
    dcs: sp.simplify(sol[dcs].subs(dP, dN0 - dD0)),
    dZ: sp.simplify(sol[dZ].subs(dP, dN0 - dD0)),
}
print("dcs(bundle) =", sol_bundle[dcs])
print("dZ(bundle)  =", sol_bundle[dZ])
expect_zero(
    "bundle identity for dcs",
    sol_bundle[dcs] - (sp.Rational(1, 2) * dKs - sp.Rational(1, 4) * dTheta + sp.Rational(1, 5) * (dN0 - dD0)),
)
expect_zero(
    "bundle identity for dZ",
    sol_bundle[dZ] - (dKq - sp.Rational(2, 5) * (dN0 - dD0)),
)

banner("Frozen-wall corollary")
sol_frozen = {k: sp.simplify(v.subs(dTheta, 0)) for k, v in sol.items()}
print("drho|frozen =", sol_frozen[drho])
print("da|frozen   =", sol_frozen[da])
print("dcs|frozen  =", sol_frozen[dcs])
print("dZ|frozen   =", sol_frozen[dZ])
expect_zero("frozen drho", sol_frozen[drho])
expect_zero("frozen da", sol_frozen[da] - sp.Rational(1, 2) * dKs)
expect_zero("frozen dcs", sol_frozen[dcs] - (sp.Rational(1, 2) * dKs + sp.Rational(1, 5) * dP))
expect_zero("frozen dZ", sol_frozen[dZ] - (dKq - sp.Rational(2, 5) * dP))

banner("Explicit Family-1 wall density")
Theta_chi = sp.Float("4.06863235008162", 30)
rho_chi = sp.sqrt(Theta_chi / 25)
print("rho_w^(chi) =", sp.N(rho_chi, 18), "* lambda_mu^(-1)")

print("\nCarry-forward formulas:")
print("  delta ln rho_w = 1/2 delta ln Theta_w")
print("  delta ln a     = 1/2 delta ln K_s - 1/4 delta ln Theta_w")
print("  delta ln c_s   = 1/2 delta ln K_s - 1/4 delta ln Theta_w + 1/5 delta ln P_0")
print("  delta ln Z_q   = delta ln K_q - 2/5 delta ln P_0")
print("  with  delta ln P_0 = delta ln N_0 - delta ln D_0  on the isotropic bundle.")
