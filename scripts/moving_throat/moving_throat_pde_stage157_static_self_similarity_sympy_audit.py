#!/usr/bin/env python3
"""
moving_throat_pde_stage157_static_self_similarity_sympy_audit.py

SymPy-backed audit for the static self-similarity decomposition of the remaining
linear grouped loading defect.

Checks:
1. Exact weighted decomposition of delta_D = D01/D0.
2. Exact wall-referenced decomposition of Xi_load.
3. Exact weighted-average formulas for BdG, conservative Maxwell/mixed, and
   outgoing-transfer logarithmic slopes.
4. Static self-similarity theorem Xi_load = 0.
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


banner("STAGE 157 — STATIC SELF-SIMILARITY DECOMPOSITION")

# ---------------------------------------------------------------------------
# 1. Exact weighted decomposition of the conservative static slope
# ---------------------------------------------------------------------------
K, B0, Z0 = sp.symbols("K B0 Z0", real=True, nonzero=True)
K1, B01, Z01 = sp.symbols("K1 B01 Z01", real=True)
N0, N01 = sp.symbols("N0 N01", real=True, nonzero=True)

D0 = K - B0 - Z0
D01 = K1 - B01 - Z01

deltaK = K1 / K
deltaB = B01 / B0
deltaZ = Z01 / Z0
deltaN = N01 / N0
deltaD = D01 / D0

omegaK = K / D0
omegaB = B0 / D0
omegaZ = Z0 / D0

expect_zero("weight identity", omegaK - omegaB - omegaZ - 1)
expect_zero(
    "delta_D weighted decomposition",
    deltaD - (omegaK * deltaK - omegaB * deltaB - omegaZ * deltaZ),
)

Xi_load = sp.simplify(deltaN - deltaD)
Xi_wall_ref = sp.simplify((deltaN - deltaK) + omegaB * (deltaB - deltaK) + omegaZ * (deltaZ - deltaK))
expect_zero("Xi_load wall-referenced form", Xi_load - Xi_wall_ref)
print("Xi_load =", Xi_load)

# ---------------------------------------------------------------------------
# 2. BdG weighted-average logarithmic slope
# ---------------------------------------------------------------------------
banner("BdG bundle weighted logarithmic slope")

c1, c2, w1, w2 = sp.symbols("c1 c2 w1 w2", real=True, nonzero=True)
dc1, dc2, dw1, dw2 = sp.symbols("dc1 dc2 dw1 dw2", real=True)

B1 = c1**2 / w1**2
B2 = c2**2 / w2**2
B0_two = sp.simplify(B1 + B2)
B01_two = sp.simplify(2 * c1 * dc1 / w1**2 - 2 * c1**2 * dw1 / w1**3 + 2 * c2 * dc2 / w2**2 - 2 * c2**2 * dw2 / w2**3)

rhoB1 = sp.simplify(B1 / B0_two)
rhoB2 = sp.simplify(B2 / B0_two)

deltaB_weighted = sp.simplify(rhoB1 * (2 * dc1 / c1 - 2 * dw1 / w1) + rhoB2 * (2 * dc2 / c2 - 2 * dw2 / w2))
expect_zero("BdG weighted-average formula", B01_two / B0_two - deltaB_weighted)

# ---------------------------------------------------------------------------
# 3. Conservative Maxwell/mixed static weighted slope
# ---------------------------------------------------------------------------
banner("Conservative Maxwell/mixed weighted logarithmic slope")

Q1, Q2, Delta1, Delta2 = sp.symbols("Q1 Q2 Delta1 Delta2", real=True, nonzero=True)
Q1p, Q2p, Delta1p, Delta2p = sp.symbols("Q1p Q2p Delta1p Delta2p", real=True)

Z1 = Q1 / Delta1
Z2 = Q2 / Delta2
Z0_two = sp.simplify(Z1 + Z2)
Z01_two = sp.simplify((Delta1 * Q1p - Q1 * Delta1p) / Delta1**2 + (Delta2 * Q2p - Q2 * Delta2p) / Delta2**2)

rhoZ1 = sp.simplify(Z1 / Z0_two)
rhoZ2 = sp.simplify(Z2 / Z0_two)

deltaZ_weighted = sp.simplify(rhoZ1 * (Q1p / Q1 - Delta1p / Delta1) + rhoZ2 * (Q2p / Q2 - Delta2p / Delta2))
expect_zero("Z weighted-average formula", Z01_two / Z0_two - deltaZ_weighted)

# ---------------------------------------------------------------------------
# 4. Outgoing-transfer weighted slope
# ---------------------------------------------------------------------------
banner("Outgoing-transfer weighted logarithmic slope")

P1, P2 = sp.symbols("P1 P2", real=True, nonzero=True)
P1p, P2p = sp.symbols("P1p P2p", real=True)

N1s = P1**2 / Delta1**2
N2s = P2**2 / Delta2**2
N0_two = sp.simplify(N1s + N2s)
N01_two = sp.simplify(2 * P1 * P1p / Delta1**2 - 2 * P1**2 * Delta1p / Delta1**3 + 2 * P2 * P2p / Delta2**2 - 2 * P2**2 * Delta2p / Delta2**3)

rhoN1 = sp.simplify(N1s / N0_two)
rhoN2 = sp.simplify(N2s / N0_two)

deltaN_weighted = sp.simplify(rhoN1 * (2 * P1p / P1 - 2 * Delta1p / Delta1) + rhoN2 * (2 * P2p / P2 - 2 * Delta2p / Delta2))
expect_zero("N weighted-average formula", N01_two / N0_two - deltaN_weighted)

# ---------------------------------------------------------------------------
# 5. Static self-similarity theorem
# ---------------------------------------------------------------------------
banner("Static self-similarity theorem")

delta_star = sp.symbols("delta_star", real=True)
Xi_self_similar = sp.simplify(Xi_load.subs({K1: delta_star*K, B01: delta_star*B0, Z01: delta_star*Z0, N01: delta_star*N0}))
expect_zero("Xi_load under static self-similarity", Xi_self_similar)

# Wall-referenced defect-field form
SigmaB, SigmaZ, SigmaN = sp.symbols("SigmaB SigmaZ SigmaN", real=True)
Xi_sigma = sp.simplify(Xi_load.subs({B01: (deltaK + SigmaB)*B0, Z01: (deltaK + SigmaZ)*Z0, N01: (deltaK + SigmaN)*N0}))
expect_zero("Xi_load mismatch-field form", Xi_sigma - (SigmaN + omegaB * SigmaB + omegaZ * SigmaZ))
print("Xi_load mismatch-field prototype =", sp.simplify(SigmaN + omegaB * SigmaB + omegaZ * SigmaZ))

print("\nCarry-forward formulas:")
print("  delta_D = omegaK deltaK - omegaB deltaB - omegaZ deltaZ")
print("  omegaK - omegaB - omegaZ = 1")
print("  Xi_load = (deltaN-deltaK) + omegaB(deltaB-deltaK) + omegaZ(deltaZ-deltaK)")
print("  delta_B = weighted average of 2 dln(c_alpha/varpi_alpha)")
print("  delta_Z = weighted average of dln(Q_r/Delta_r)")
print("  delta_N = weighted average of 2 dln(P_r/Delta_r)")
print("  static self-similarity => Xi_load = 0")
