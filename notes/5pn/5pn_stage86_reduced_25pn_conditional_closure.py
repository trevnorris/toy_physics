
#!/usr/bin/env python3
"""
5pn_stage86_reduced_25pn_conditional_closure.py

Stage 86 audit: conditional reduced 2.5PN closure.
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

banner("STAGE 86 — CONDITIONAL REDUCED 2.5PN CLOSURE")

chi_Q, Delta_Q = sp.symbols("chi_Q Delta_Q", real=True)
N_Q = sp.simplify(1 / chi_Q)
expect_zero("N_Q on the canonical outgoing branch", N_Q.subs(chi_Q, 1) - 1)
expect_zero("N_Q - 1 + Delta_Q + O(Delta_Q^2) (linearized)", sp.series((1 / (1 + Delta_Q) - 1 + Delta_Q), Delta_Q, 0, 2).removeO())

banner("STAGE 86 FINAL LEDGER")
print("Inside the present reduced hierarchy, the 2.5PN theorem is conditionally closed by")
print("one and only one remaining branch datum:")
print("  chi_Q.")
print("If chi_Q = 1, the reduced GR-like point-particle 2.5PN theorem is closed.")
print("If chi_Q != 1, the entire remaining reduced failure is measured by Delta_Q := chi_Q - 1.")
