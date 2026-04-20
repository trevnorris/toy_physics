#!/usr/bin/env python3
"""
moving_throat_pde_stage143_bare_mixed_port_slippage_sympy_audit.py

SymPy-backed audit for Stage 143.

Checks:
1. Exact compensated-branch identity
     d gamma_W - (1/3) d kappa_W
     = (d gamma_0 - (1/3) d kappa_0)/(1+r_c)
   for the concrete core model.
2. Collapse under the Stage-142 canonical-even gate d kappa_W = 0.
3. Pure-scale harmlessness: if d gamma_0 = d kappa_0 / 3 then d gamma_W = 0.
4. Final defect laws after inserting the tangential-mouth transport.
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

banner("STAGE 143 — BARE MIXED-PORT SLIPPAGE THEOREM")

eps = sp.symbols("eps", real=True)
rc, drc = sp.symbols("r_c dr_c", real=True)
dk0, dg0 = sp.symbols("dκ0 dγ0", real=True)

k0_star = (1 + rc) / 3
g0_star = (1 + rc) / 9

kW = ((k0_star + eps * dk0) / (1 + rc + eps * drc)).series(eps, 0, 2).removeO()
gW = ((g0_star + eps * dg0) / (1 + rc + eps * drc)).series(eps, 0, 2).removeO()

dkW = sp.expand(kW.coeff(eps, 1))
dgW = sp.expand(gW.coeff(eps, 1))

print("dκ_W =", sp.simplify(dkW))
print("dγ_W =", sp.simplify(dgW))

identity = dgW - sp.Rational(1, 3) * dkW - (dg0 - sp.Rational(1, 3) * dk0) / (1 + rc)
expect_zero("exact compensated-branch slippage identity", identity)

# Canonical-even gate dκ_W = 0
dgW_gate = sp.simplify(((dg0 - sp.Rational(1, 3) * dk0) / (1 + rc)))
print("dγ_W under dκ_W = 0 =", dgW_gate)

# Pure-scale harmlessness
dgW_purescale = sp.simplify(dgW_gate.subs(dg0, dk0 / 3))
expect_zero("pure-scale harmlessness", dgW_purescale)

banner("Tangential DtN susceptibility and final defect law")

UpsilonPi = sp.symbols("Upsilon_Pi", real=True)
dSigma0, dS, sigma_star = sp.symbols("dSigma0 dS sigma_*", real=True)

dPi_tan = sp.Float("0.832409471081635") * dSigma0 - sp.Float("1.16275838754222") * dS
dgW_tan = sp.simplify(UpsilonPi * dPi_tan / (1 + rc))
DeltaQ = sp.simplify(-9 * sigma_star * dgW_tan / (1 - sigma_star))
NQm1 = sp.simplify(9 * sigma_star * dgW_tan / (1 - sigma_star))

print("dPi_tan =", dPi_tan)
print("dγ_W =", dgW_tan)
print("Δ_Q =", DeltaQ)
print("N_Q - 1 =", NQm1)

banner("Carry-forward formulas")
print("1) dγ_W - (1/3)dκ_W = (dγ_0 - (1/3)dκ_0)/(1+r_c)")
print("2) with dκ_W = 0: dγ_W = (dγ_0 - (1/3)dκ_0)/(1+r_c)")
print("3) if dγ_0 = dκ_0/3, then dγ_W = 0")
print("4) if dγ_0 - dκ_0/3 = Upsilon_Pi dPi_tan, then")
print("   Δ_Q = -9 σ_* Upsilon_Pi dPi_tan / [(1-σ_*)(1+r_c)]")
print("   N_Q-1 = +9 σ_* Upsilon_Pi dPi_tan / [(1-σ_*)(1+r_c)]")
