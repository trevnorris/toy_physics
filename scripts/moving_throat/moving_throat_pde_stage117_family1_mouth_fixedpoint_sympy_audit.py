
#!/usr/bin/env python3
"""
moving_throat_pde_stage117_family1_mouth_fixedpoint_sympy_audit.py

SymPy audit for the Family-1 shell + first mixed D/N tube reduction.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

Pi = sp.symbols("Pi", positive=True, real=True)
M_s, M_q = sp.symbols("M_s M_q", real=True)

def S(Pi, kappa):
    return sp.simplify(
        Pi * (kappa * sp.tanh(kappa) + Pi * (sp.exp(-Pi) / sp.cosh(kappa) - 1))
        / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2))
    )

banner("STAGE 117 — FAMILY-1 FIXED-POINT REDUCTION")
kk = sp.symbols("kk", positive=True, real=True)
S_shell = sp.simplify(sp.limit(S(Pi, kk), kk, 0))
S_q = sp.simplify(S(Pi, sp.pi/2))
Pi_eq = sp.simplify(M_s + M_q * S_q)

print("S_shell =", S_shell)
print("S_q(Pi) =")
sp.pprint(S_q)
print("Fixed-point law Pi =")
sp.pprint(Pi_eq)

Pi_star = sp.Float("1.50882951349316")
S_star = sp.N(S_q.subs(Pi, Pi_star), 30)
print("S_q(Pi_star) =", S_star)
print("Canonical gain line: M_s = Pi_star - S_q(Pi_star) M_q")
sp.pprint(sp.N(Pi_star - S_star*M_q, 30))
