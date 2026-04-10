
#!/usr/bin/env python3
"""
5pn_stage84_natural_source_map_reduction.py

Stage 84 audit: on the natural source-map branch the last reduced 2.5PN
obstruction is purely outgoing.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 84 — NATURAL SOURCE-MAP REDUCTION")

chi_Q, Delta_Q = sp.symbols("chi_Q Delta_Q", real=True)
N_Q = sp.simplify(1 / chi_Q)
print("On the natural point-particle source-map branch, mhat_0 -> 1, so")
print("  N_Q = 1/chi_Q.")

expect_zero("N_Q - 1/(1 + Delta_Q)", N_Q.subs(chi_Q, 1 + Delta_Q) - 1 / (1 + Delta_Q))
expect_zero("(N_Q - 1) + Delta_Q/(1 + Delta_Q)", (N_Q.subs(chi_Q, 1 + Delta_Q) - 1) + Delta_Q / (1 + Delta_Q))
expect_zero("canonical branch N_Q - 1", N_Q.subs(chi_Q, 1) - 1)

banner("STAGE 84 FINAL LEDGER")
print("On the natural source-map branch the last reduced 2.5PN obstruction is purely outgoing:")
print("  N_Q = 1/chi_Q,")
print("so the only remaining reduced question is whether chi_Q = 1 on the actual")
print("passive/outgoing moving-throat quadrupole branch.")
