
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

def expect_zero(name: str, expr):
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 341 — ORBIT PACKET TO THREE DIRECT OBSERVABLES")

# Direct coherent branch observables
dln_Rtr, dln_Rtarget, dln_epseta = sp.symbols("dln_Rtr dln_Rtarget dln_eps_eta", real=True)
epseta = sp.symbols("epsilon_eta", positive=True, real=True)

# Physical defect triplet used later in the notes
Theta1 = dln_Rtr
Xi1 = -dln_Rtarget - epseta * dln_epseta / (1 - epseta)
R1 = dln_Rtarget

subbanner("I. Exact direct-observable defect map")
print("Theta_1 =")
sp.pprint(Theta1)
print("Xi_1 =")
sp.pprint(Xi1)
print("R_1 =")
sp.pprint(R1)

# Exact inverse map
inv_dln_Rtr = Theta1
inv_dln_Rtarget = R1
inv_dln_epseta = - (1 - epseta) * (R1 + Xi1) / epseta

subbanner("II. Exact inverse map")
print("d ln R_tr =")
sp.pprint(inv_dln_Rtr)
print("d ln R_target =")
sp.pprint(inv_dln_Rtarget)
print("d ln epsilon_eta =")
sp.pprint(inv_dln_epseta)

expect_zero("roundtrip Theta_1", inv_dln_Rtr - dln_Rtr)
expect_zero(
    "roundtrip Xi_1",
    (-inv_dln_Rtarget - epseta * inv_dln_epseta / (1 - epseta)) - Xi1,
)
expect_zero("roundtrip R_1", inv_dln_Rtarget - dln_Rtarget)

subbanner("III. Zero-set equivalence")
print("Exact zero-set theorem:")
print("  Theta_1 = Xi_1 = R_1 = 0")
print("    iff")
print("  d ln R_tr = d ln R_target = d ln epsilon_eta = 0")

banner("STAGE 341 LEDGER")
print("1. The full weak-axisymmetric orbit-lock problem is equivalent to invariance of three direct observables.")
print("2. Those observables are (R_tr, R_target, epsilon_eta).")
print("3. Vanishing of (Theta_1, Xi_1, R_1) is exactly equivalent to vanishing drift of those three observables.")
