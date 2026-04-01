#!/usr/bin/env python3
"""
moving_throat_pde_stage169_similarity_orbit_closure_sympy_audit.py

SymPy audit for Stage 169: exact microscopic similarity orbit and weak-axisymmetric
closure theorem.

Checks:
1. Build the exact 3x8 monomial-drift matrix M_*.
2. Verify the convenient 3x3 minor det = 1 + chi_*.
3. Solve the linear compatibility ledger for (tau_1, kappa_eta, mu_1).
4. Construct the finite five-parameter similarity orbit and verify exact
   preservation of the three direct microscopic monomials.
5. Verify that the linearization of the finite orbit reproduces the Stage-168
   compatibility formulas exactly.
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

banner("STAGE 169 — EXACT MICROSCOPIC SIMILARITY ORBIT")

chi, delta = sp.symbols("chi delta", positive=True, real=True)
E, F = sp.symbols("E F", real=True)

# Microscopic grouped weak-axisymmetric drifts
lam1, c1, gam1, kU, kEta, kW, mu1, tau1 = sp.symbols(
    "lam1 c1 gam1 kU kEta kW mu1 tau1", real=True
)

# Direct microscopic monomial drifts from Stage 168
dlog_Ctr = (1 + delta) * (gam1 + c1 - kU) + (1 + chi) * (tau1 - kU)
dlog_Cnt = (2 + 2*E) * lam1 + 2*E * gam1 + mu1 + (F - E) * kU - kEta - (2 + E) * kW - F * tau1
dlog_eta = 2*c1 - kU - kEta

print("delta log C_tr =", sp.expand(dlog_Ctr))
print("delta log C_nt =", sp.expand(dlog_Cnt))
print("delta log eps_eta =", sp.expand(dlog_eta))

M = sp.Matrix([
    [0, 1 + delta, 1 + delta, -(2 + chi + delta), 0, 0, 0, 1 + chi],
    [2 + 2*E, 0, 2*E, F - E, -1, -(2 + E), 1, -F],
    [0, 2, 0, -1, -1, 0, 0, 0],
])
print("\nM_* =")
sp.pprint(M)

minor = M[:, [7, 4, 6]]  # (tau1, kEta, mu1)
det_minor = sp.simplify(sp.factor(minor.det()))
print("\nminor det(tau1,kEta,mu1) =", det_minor)
if det_minor != 1 + chi:
    raise AssertionError("Unexpected rank-3 minor determinant.")

banner("Linear compatibility solve")
sol = sp.solve(
    [sp.Eq(dlog_Ctr, 0), sp.Eq(dlog_eta, 0), sp.Eq(dlog_Cnt, 0)],
    [tau1, kEta, mu1],
    dict=True,
)[0]

tau_formula = sp.simplify(sol[tau1])
kEta_formula = sp.simplify(sol[kEta])
mu_formula = sp.simplify(sol[mu1])

print("tau1 =", tau_formula)
print("kEta =", kEta_formula)
print("mu1  =", mu_formula)

expect_zero(
    "tau1 - [kU - ((1+delta)/(1+chi)) (gam1+c1-kU)]",
    tau_formula - (kU - (1 + delta) * (gam1 + c1 - kU) / (1 + chi)),
)
expect_zero("kEta - (2 c1 - kU)", kEta_formula - (2*c1 - kU))
expect_zero(
    "mu1 - Stage-168 form",
    mu_formula
    - (
        (2*c1 - kU)
        + 2*kW
        - 2*lam1
        - E * (2*gam1 + 2*lam1 - kU - kW)
        - F * (1 + delta) * (gam1 + c1 - kU) / (1 + chi)
    ),
)

banner("Finite five-parameter similarity orbit")
Lam, C, Gam, U, W = sp.symbols("Lam C Gam U W", real=True)

Eta_exp = 2*C - U
Tau_exp = U - (1 + delta) * (Gam + C - U) / (1 + chi)
Mu_exp = Eta_exp + 2*W - 2*Lam - E * (2*Gam + 2*Lam - U - W) + F * (Tau_exp - U)

print("K_eta exponent =", Eta_exp)
print("T_U exponent   =", Tau_exp)
print("mu_W exponent  =", sp.simplify(Mu_exp))

# Exact monomial preservation
Ctr_orbit = (1 + delta) * (Gam + C - U) + (1 + chi) * (Tau_exp - U)
Cnt_orbit = (2 + 2*E) * Lam + 2*E * Gam + Mu_exp + (F - E) * U - Eta_exp - (2 + E) * W - F * Tau_exp
Eta_orbit = 2*C - U - Eta_exp

expect_zero("finite orbit preserves C_tr", Ctr_orbit)
expect_zero("finite orbit preserves C_nt", Cnt_orbit)
expect_zero("finite orbit preserves eps_eta", Eta_orbit)

banner("Linearization reproduces compatibility ledger")
subs_basis = {Lam: lam1, C: c1, Gam: gam1, U: kU, W: kW}
expect_zero("linearized tau formula", Tau_exp.subs(subs_basis) - tau_formula)
expect_zero("linearized kEta formula", Eta_exp.subs(subs_basis) - kEta_formula)
expect_zero("linearized mu formula", Mu_exp.subs(subs_basis) - mu_formula)

print("\nConclusion:")
print("  The three Stage-168 compatibility equations are exactly the tangent-space")
print("  equations of a finite five-parameter multiplicative similarity orbit.")
print("  The coherent weak-axisymmetric zero-defect branch is therefore codimension 3,")
print("  not an isolated fine-tuning point.")
