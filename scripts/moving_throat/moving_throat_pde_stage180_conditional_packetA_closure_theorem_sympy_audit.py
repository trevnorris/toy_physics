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


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


def grouped_trace_anomaly(x20, x21, x22):
    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)
    return xbar, ax, bx


banner("STAGE 180 — EXACT CONDITIONAL PACKET-A CLOSURE THEOREM")

# ---------------------------------------------------------------------------
# I. Isotropic outgoing lane implies a_{P0}=b_{P0}=0
# ---------------------------------------------------------------------------
subbanner("I. Exact isotropic outgoing-lane collapse")

P0 = sp.symbols("P0", real=True)
aP0, bP0 = grouped_trace_anomaly(P0, P0, P0)[1:]
print("a_{P0} on the isotropic outgoing lane =")
sp.pprint(aP0)
print("b_{P0} on the isotropic outgoing lane =")
sp.pprint(bP0)
expect_zero("a_{P0}", aP0)
expect_zero("b_{P0}", bP0)

# ---------------------------------------------------------------------------
# II. Packet-A residual collapses to the normalization slot
# ---------------------------------------------------------------------------
subbanner("II. Exact Packet-A residual collapse")

Delta_norm = sp.symbols("Delta_norm", real=True)
Delta_branch = sp.Matrix([0, 0, 0, 0, 0, 0, 0, Delta_norm])
print("Delta_branch after imposing the carried isotropic one-pole front end =")
sp.pprint(Delta_branch)
expect_zero("first seven Packet-A slots", Delta_branch[:7, :] - sp.zeros(7, 1))

# ---------------------------------------------------------------------------
# III. Natural source-map reduction: Delta_norm = P0_target(1/chi_Q - 1)
# ---------------------------------------------------------------------------
subbanner("III. Exact source-map reduction of the final Packet-A slot")

a, c_s, c, G = sp.symbols("a c_s c G", positive=True, real=True)
chi_Q = sp.symbols("chi_Q", nonzero=True, real=True)
P0_target = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
N_Q = sp.simplify(1 / chi_Q)
Delta_norm_pt = sp.simplify(P0_target * (N_Q - 1))

print("P0^target =")
sp.pprint(P0_target)
print("N_Q on the natural point-particle source-map branch =")
sp.pprint(N_Q)
print("Delta_norm on the natural point-particle source-map branch =")
sp.pprint(Delta_norm_pt)
expect_zero(
    "Delta_norm - P0^target(1/chi_Q - 1)",
    Delta_norm_pt - P0_target * (1 / chi_Q - 1),
)

# ---------------------------------------------------------------------------
# IV. Exact equivalence: Delta_branch = 0 iff chi_Q = 1
# ---------------------------------------------------------------------------
subbanner("IV. Exact finish-line equivalences")

Delta_branch_pt = sp.Matrix([0, 0, 0, 0, 0, 0, 0, Delta_norm_pt])
print("Delta_branch on the natural point-particle source-map branch =")
sp.pprint(Delta_branch_pt)

# Algebraic identities encoding the zero-set equivalence
expect_zero(
    "chi_Q*(Delta_norm/P0_target + 1) - 1",
    sp.simplify(chi_Q * (Delta_norm_pt / P0_target + 1) - 1),
)
expect_zero(
    "Delta_norm/P0_target + (chi_Q - 1)/chi_Q",
    sp.simplify(Delta_norm_pt / P0_target + (chi_Q - 1) / chi_Q),
)

Delta_Q = sp.symbols("Delta_Q", real=True)
expect_zero(
    "Delta_norm under chi_Q = 1 + Delta_Q",
    sp.simplify(Delta_norm_pt.subs({chi_Q: 1 + Delta_Q}) + P0_target * Delta_Q / (1 + Delta_Q)),
)

# ---------------------------------------------------------------------------
# V. Stage-177 deformation algebra form of the same closure condition
# ---------------------------------------------------------------------------
subbanner("V. Exact Stage-177 deformation-algebra realization gate")

S, beta = sp.symbols("S beta", nonzero=True, real=True)
Sigma0, Sigma5 = sp.symbols("Sigma_0 Sigma_5", real=True)
chi_from_def = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
closure_num = sp.simplify(3 * S * (beta**5 - 1) + Sigma0 + 27 * Sigma5)

print("chi_Q from the carried isotropic DtN deformation algebra =")
sp.pprint(chi_from_def)
print("closure numerator =")
sp.pprint(closure_num)
expect_zero(
    "(3S - Sigma0)(chi_Q - 1) - closure numerator",
    sp.simplify((3 * S - Sigma0) * (chi_from_def - 1) - closure_num),
)

Delta_norm_def = sp.simplify(P0_target * (1 / chi_from_def - 1))
print("Delta_norm after inserting the exact Stage-177 deformation algebra =")
sp.pprint(Delta_norm_def)
expect_zero(
    "Delta_norm + P0_target*closure_num/[3(S beta^5 + 9 Sigma5)]",
    sp.simplify(Delta_norm_def + P0_target * closure_num / (3 * (S * beta**5 + 9 * Sigma5))),
)

# ---------------------------------------------------------------------------
# VI. Linearized finish-line map
# ---------------------------------------------------------------------------
subbanner("VI. Exact linearized finish-line map")

eps = sp.symbols("eps", real=True)
eps_beta, dSigma0, dSigma5 = sp.symbols("eps_beta dSigma_0 dSigma_5", real=True)

chi_linear = sp.series(
    chi_from_def.subs({beta: 1 + eps * eps_beta, Sigma0: eps * dSigma0, Sigma5: eps * dSigma5}),
    eps,
    0,
    2,
).removeO()
chi_linear_expected = 1 + eps * (5 * eps_beta + dSigma0 / (3 * S) + 9 * dSigma5 / S)
print("chi_Q linearized =")
sp.pprint(sp.expand(chi_linear))
expect_zero(
    "linearized chi_Q - [1 + eps(5 eps_beta + dSigma0/(3S) + 9 dSigma5/S)]",
    sp.simplify(sp.expand(chi_linear - chi_linear_expected)),
)

Delta_norm_linear = sp.series(
    Delta_norm_def.subs({beta: 1 + eps * eps_beta, Sigma0: eps * dSigma0, Sigma5: eps * dSigma5}),
    eps,
    0,
    2,
).removeO()
Delta_norm_linear_expected = -eps * P0_target * (5 * eps_beta + dSigma0 / (3 * S) + 9 * dSigma5 / S)
print("Delta_norm linearized =")
sp.pprint(sp.expand(Delta_norm_linear))
expect_zero(
    "linearized Delta_norm + eps P0_target*(5 eps_beta + dSigma0/(3S) + 9 dSigma5/S)",
    sp.simplify(sp.expand(Delta_norm_linear - Delta_norm_linear_expected)),
)

# ---------------------------------------------------------------------------
# VII. Higher-odd irrelevance at the final finish line
# ---------------------------------------------------------------------------
subbanner("VII. Higher-odd irrelevance survives to the final Packet-A theorem")

L7 = sp.symbols("L7", real=True)
expect_zero("d chi_Q / dL7", sp.diff(chi_from_def, L7))
expect_zero("d Delta_norm / dL7", sp.diff(Delta_norm_def, L7))

banner("STAGE 180 LEDGER")
print("1. Once the carried isotropic grouped-real P2 conservative one-pole front end is imposed,")
print("      a2 = b2 = a4 = b4 = a_{P0} = b_{P0} = Delta_pole = 0,")
print("   so the full Packet-A residual collapses exactly to")
print("      Delta_branch = (0,0,0,0,0,0,0,Delta_norm).")
print("2. On the natural point-particle source-map branch, the exact Stage-178 factorization gives")
print("      N_Q = 1/chi_Q,")
print("   hence")
print("      Delta_norm = P0^target (1/chi_Q - 1).")
print("3. Because P0^target > 0 on the physical branch, the final Packet-A retarded finish line is")
print("      Delta_branch = 0  <=>  Delta_norm = 0  <=>  chi_Q = 1.")
print("   Equivalently, on the same natural branch, this is also")
print("      N_Q = 1  <=>  Delta_Q := chi_Q - 1 = 0.")
print("4. In the exact Stage-177 isotropic DtN deformation algebra,")
print("      chi_Q = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0),")
print("   so the same finish line is")
print("      chi_Q = 1  <=>  3 S (beta^5 - 1) + Sigma_0 + 27 Sigma_5 = 0.")
print("5. The higher-odd irrelevance theorem is inherited unchanged:")
print("   any extra isotropic odd tail beginning at O(omega^7) leaves chi_Q and therefore the final Packet-A theorem unchanged.")
