#!/usr/bin/env python3
"""
moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py

SymPy-backed audit for the finite orbit–quotient closure theorem.

Checks:
1. The exact finite log-ratio equalities of the three monomial invariants are
   exactly M_* Delta x = 0.
2. Solving the finite fibre equations for (Delta_eta, Delta_T, Delta_mu)
   reproduces the Stage 186 similarity-orbit laws.
3. Substituting those exact finite laws annihilates all invariant log-ratio
   equations.
4. The same three monomials are therefore complete finite orbit invariants, so
   the weak-axisymmetric defect lives in the exact three-dimensional quotient.
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


banner("STAGE 187 — EXACT ORBIT–QUOTIENT CLOSURE")

# Finite log-ratio variables between two positive microscopic states.
DL, DC, DG, DU, DEta, DW, DM, DT = sp.symbols(
    "Delta_lambda Delta_c Delta_gamma Delta_U Delta_eta Delta_W Delta_mu Delta_T",
    real=True,
)
chi, deltaU, E, F = sp.symbols("chi0_star deltaU_star E_star F_star", positive=True, real=True)

# Exact monomial-ratio equations in logarithmic form.
row_tr = (1 + deltaU) * (DG + DC - DU) + (1 + chi) * (DT - DU)
row_nt = 2 * (1 + E) * DL + 2 * E * DG + (F - E) * DU - DEta - (2 + E) * DW + DM - F * DT
row_eta = 2 * DC - DU - DEta

print("Exact finite log-ratio equations:")
print("row_tr  =", row_tr)
print("row_nt  =", row_nt)
print("row_eta =", row_eta)

# Derive the rows from the actual Stage-187 monomials at finite level.
# Positive primitive ratios (xtilde/x) for the eight microscopic variables.
rL, rC, rG, rU, rEta, rW, rM, rT = sp.symbols(
    "r_lambda r_c r_gamma r_U r_eta r_W r_mu r_T", positive=True
)
log_subs = {
    DL: sp.log(rL),
    DC: sp.log(rC),
    DG: sp.log(rG),
    DU: sp.log(rU),
    DEta: sp.log(rEta),
    DW: sp.log(rW),
    DM: sp.log(rM),
    DT: sp.log(rT),
}
Ctr_ratio = (rG * rC / rU) ** (1 + deltaU) * (rT / rU) ** (1 + chi)
Cnt_ratio = (rL ** 2 * rM / (rEta * rW ** 2)) * (
    rG ** 2 * rL ** 2 / (rU * rW)
) ** E * (rT / rU) ** (-F)
eps_ratio = rC ** 2 / (rU * rEta)
expect_zero(
    "log C_tr ratio - row_tr",
    sp.expand_log(sp.log(Ctr_ratio), force=True) - row_tr.subs(log_subs),
)
expect_zero(
    "log C_nt ratio - row_nt",
    sp.expand_log(sp.log(Cnt_ratio), force=True) - row_nt.subs(log_subs),
)
expect_zero(
    "log epsilon_eta ratio - row_eta",
    sp.expand_log(sp.log(eps_ratio), force=True) - row_eta.subs(log_subs),
)

# Same matrix as Stage 186, but now acting on finite log-ratios.
M = sp.Matrix([
    [0, 1 + deltaU, 1 + deltaU, -(2 + chi + deltaU), 0, 0, 0, 1 + chi],
    [2 * (1 + E), 0, 2 * E, F - E, -1, -(2 + E), 1, -F],
    [0, 2, 0, -1, -1, 0, 0, 0],
])
Dx = sp.Matrix([DL, DC, DG, DU, DEta, DW, DM, DT])
Mx = sp.expand(M * Dx)

expect_zero("matrix row 1 - exact row_tr", Mx[0] - row_tr)
expect_zero("matrix row 2 - exact row_nt", Mx[1] - row_nt)
expect_zero("matrix row 3 - exact row_eta", Mx[2] - row_eta)

minor = sp.Matrix([
    [0, 0, 1 + chi],
    [-1, 1, -F],
    [-1, 0, 0],
])
print("det selected minor (Delta_eta, Delta_mu, Delta_T) =", sp.simplify(minor.det()))
expect_zero("selected minor determinant", minor.det() - (1 + chi))

# Exact finite solve.
sol = sp.solve([row_tr, row_nt, row_eta], [DEta, DT, DM], dict=True)[0]
print("\nExact finite fibre solution:")
for k in (DEta, DT, DM):
    print(f"{k} = {sp.simplify(sol[k])}")

DEta_expected = 2 * DC - DU
DT_expected = DU - (1 + deltaU) * (DG + DC - DU) / (1 + chi)
DM_expected = (
    2 * DC - DU + 2 * DW - 2 * DL
    - E * (2 * DG + 2 * DL - DU - DW)
    - F * (1 + deltaU) * (DG + DC - DU) / (1 + chi)
)

expect_zero("Delta_eta finite law", sol[DEta] - DEta_expected)
expect_zero("Delta_T finite law", sol[DT] - DT_expected)
expect_zero("Delta_mu finite law", sol[DM] - DM_expected)

# Substitute back and verify the exact invariant-fibre equations vanish.
expect_zero("row_tr after solve", row_tr.subs(sol))
expect_zero("row_nt after solve", row_nt.subs(sol))
expect_zero("row_eta after solve", row_eta.subs(sol))

banner("Finite orbit interpretation")
print("The three monomial equalities reduce exactly to the same rank-3 matrix condition")
print("M_* Delta x = 0, but now Delta x is a finite log-ratio vector rather than an")
print("infinitesimal drift. Therefore the Stage 186 similarity orbit integrates exactly:")
print("its fibres are the full level sets of (C_tr, C_nt, epsilon_eta).")

banner("Carry-forward formulas")
print("  Delta_eta = 2 Delta_c - Delta_U")
print("  Delta_T   = Delta_U - ((1+deltaU_*)/(1+chi_*))(Delta_gamma + Delta_c - Delta_U)")
print("  Delta_mu  = 2 Delta_c - Delta_U + 2 Delta_W - 2 Delta_lambda")
print("              - E_*(2 Delta_gamma + 2 Delta_lambda - Delta_U - Delta_W)")
print("              - F_*((1+deltaU_*)/(1+chi_*))(Delta_gamma + Delta_c - Delta_U)")
print("\nConclusion:")
print("  Equal values of (C_tr, C_nt, epsilon_eta) are exactly equivalent to lying on")
print("  the same finite similarity orbit G_*. The weak-axisymmetric defect therefore")
print("  lives in the exact three-dimensional orbit quotient.")
