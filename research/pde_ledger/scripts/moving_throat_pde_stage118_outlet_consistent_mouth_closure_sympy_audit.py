
#!/usr/bin/env python3
"""
moving_throat_pde_stage118_outlet_consistent_mouth_closure_sympy_audit.py

SymPy audit for the outlet-consistent mouth closure:
    M_s = 4 Sigma_m,   M_q = -Sigma_m.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

Pi, Sigma_m = sp.symbols("Pi Sigma_m", positive=True, real=True)

def S(Pi, kappa):
    return sp.simplify(
        Pi * (kappa * sp.tanh(kappa) + Pi * (sp.exp(-Pi) / sp.cosh(kappa) - 1))
        / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2))
    )

banner("STAGE 118 — OUTLET-CONSISTENT ONE-PARAMETER CLOSURE")
S_q = sp.simplify(S(Pi, sp.pi/2))
Pi_eq = sp.simplify(Sigma_m * (4 - S_q))
print("Pi =")
sp.pprint(Pi_eq)

Pi_star = sp.Float("1.50882951349316")
S_star = sp.N(S_q.subs(Pi, Pi_star), 30)
Sigma_var = sp.symbols("Sigma_var", positive=True, real=True)
Sigma_star = sp.N(sp.solve(sp.Eq(Pi_star, Sigma_var * (4 - S_star)), Sigma_var)[0], 30)
M_s_star = sp.N(4 * Sigma_star, 30)
M_q_star = sp.N(-Sigma_star, 30)

print("S_q(Pi_star) =", S_star)
print("0 < S_q(Pi_star) < 1 ->", bool(0 < S_star < 1))
print("Sigma_m^* =", Sigma_star)
print("M_s^* =", M_s_star)
print("M_q^* =", M_q_star)

residual = sp.N(Pi_star - Sigma_star * (4 - S_star), 30)
print("Pi_star - Sigma_star*(4 - S_star) =", residual)
if abs(residual) > 1e-12:
    raise AssertionError("Outlet-consistent threshold did not close.")
