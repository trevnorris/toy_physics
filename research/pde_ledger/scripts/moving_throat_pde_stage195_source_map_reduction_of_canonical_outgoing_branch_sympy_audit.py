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


banner("STAGE 178 — EXACT SOURCE-MAP REDUCTION OF THE CANONICAL OUTGOING BRANCH")

# ---------------------------------------------------------------------------
# I. Carry-forward isotropic retarded invariant tuple
# ---------------------------------------------------------------------------
subbanner("I. Carry-forward isotropic retarded invariant tuple")

a, c_s, c, G = sp.symbols("a c_s c G", positive=True, real=True)
P0, chi_Q, mhat0 = sp.symbols("P0 chi_Q mhat_0", positive=True, real=True)

Omega_Q = sp.simplify(sp.Integer(3) * c_s / (2 * a))
P0_target = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
Gamma5_target = sp.simplify(2 * G / (5 * c**5))
Gamma5 = sp.simplify(chi_Q * a**5 * P0 / (27 * c_s**5))
Gamma5_alt = sp.simplify(chi_Q * 9 * P0 / (32 * Omega_Q**5))

N_Q = sp.symbols("N_Q", real=True)
N_Q_def = sp.simplify(P0 / P0_target)

print("Omega_Q =")
sp.pprint(Omega_Q)
print("P0_target =")
sp.pprint(P0_target)
print("Gamma5_target =")
sp.pprint(Gamma5_target)
print("Gamma5 =")
sp.pprint(Gamma5)
expect_zero("Gamma5 - chi_Q*9*P0/(32 Omega_Q^5)", Gamma5 - Gamma5_alt)
expect_zero(
    "Gamma5/Gamma5_target - chi_Q*N_Q",
    sp.simplify(Gamma5 / Gamma5_target).subs(P0, N_Q * P0_target) - chi_Q * N_Q,
)

# ---------------------------------------------------------------------------
# II. Exact observable odd-closure factorization and Packet-A normalization collapse
# ---------------------------------------------------------------------------
subbanner("II. Exact odd-closure factorization and Packet-A normalization collapse")

Delta_norm = sp.simplify(mhat0**2 * P0 - P0_target)
Delta_norm_NQ = sp.simplify(Delta_norm.subs(P0, N_Q * P0_target))
odd_closure = sp.simplify(mhat0**2 * chi_Q * N_Q - 1)
NQ_from_odd = sp.simplify(1 / (mhat0**2 * chi_Q))
Delta_norm_from_odd = sp.simplify(Delta_norm_NQ.subs(N_Q, NQ_from_odd))

print("Delta_norm =")
sp.pprint(Delta_norm_NQ)
print("odd closure =")
sp.pprint(odd_closure)
print("N_Q from odd closure =")
sp.pprint(NQ_from_odd)
print("Delta_norm after imposing odd closure =")
sp.pprint(Delta_norm_from_odd)
odd_condition_residual = sp.simplify(
    (mhat0**2 * Gamma5 - Gamma5_target).subs(P0, N_Q * P0_target)
)
expect_zero(
    "observable odd condition factorizes as Gamma5_target*(mhat0^2 chi_Q N_Q - 1)",
    odd_condition_residual - Gamma5_target * (mhat0**2 * chi_Q * N_Q - 1),
)
expect_zero("Delta_norm from odd closure - P0_target*(1/chi_Q - 1)", Delta_norm_from_odd - P0_target * (1 / chi_Q - 1))

Delta_Q = sp.symbols("Delta_Q", real=True)
Delta_norm_DeltaQ = sp.simplify(Delta_norm_from_odd.subs(chi_Q, 1 + Delta_Q))
expect_zero(
    "Delta_norm in terms of Delta_Q",
    Delta_norm_DeltaQ + P0_target * Delta_Q / (1 + Delta_Q),
)

# ---------------------------------------------------------------------------
# III. Natural source-map reduction
# ---------------------------------------------------------------------------
subbanner("III. Natural source-map reduction")

NQ_pt = sp.simplify(NQ_from_odd.subs(mhat0, 1))
Delta_norm_pt = sp.simplify(Delta_norm_NQ.subs({mhat0: 1, N_Q: NQ_pt}))

print("N_Q on the natural point-particle source-map branch =")
sp.pprint(NQ_pt)
print("Delta_norm on the natural point-particle source-map branch =")
sp.pprint(Delta_norm_pt)
expect_zero("N_Q(point-particle) - 1/chi_Q", NQ_pt - 1 / chi_Q)
expect_zero("Delta_norm(point-particle) - P0_target*(1/chi_Q - 1)", Delta_norm_pt - P0_target * (1 / chi_Q - 1))
expect_zero(
    "N_Q(point-particle)-1 in terms of Delta_Q",
    sp.simplify((NQ_pt - 1).subs(chi_Q, 1 + Delta_Q)) + Delta_Q / (1 + Delta_Q),
)

# ---------------------------------------------------------------------------
# IV. Stage 194 deformation algebra after source-map reduction
# ---------------------------------------------------------------------------
subbanner("IV. Stage 194 deformation algebra after source-map reduction")

S, beta = sp.symbols("S beta", nonzero=True, real=True)
Sigma0, Sigma5 = sp.symbols("Sigma_0 Sigma_5", real=True)

chi_from_def = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
NQ_from_def = sp.simplify(1 / chi_from_def)
Delta_norm_from_def = sp.simplify(P0_target * (NQ_from_def - 1))
Delta_norm_from_def_expected = sp.simplify(
    -P0_target * (3 * S * (beta**5 - 1) + Sigma0 + 27 * Sigma5) / (3 * (S * beta**5 + 9 * Sigma5))
)

print("chi_Q from Stage 194 DtN deformation algebra =")
sp.pprint(chi_from_def)
print("N_Q after natural source-map reduction =")
sp.pprint(NQ_from_def)
print("Delta_norm(point-particle) after inserting DtN deformation algebra =")
sp.pprint(Delta_norm_from_def)
expect_zero(
    "N_Q deformation law",
    NQ_from_def - (3 * S - Sigma0) / (3 * (S * beta**5 + 9 * Sigma5)),
)
expect_zero("Delta_norm deformation law", Delta_norm_from_def - Delta_norm_from_def_expected)

# First-order linearization about the canonical compact branch beta=1, Sigma0=0, Sigma5=0.
eps_beta, dSigma0, dSigma5, t = sp.symbols("eps_beta dSigma_0 dSigma_5 t", real=True)
NQ_minus_1_linear = sp.simplify(
    (NQ_from_def - 1)
    .subs({beta: 1 + t * eps_beta, Sigma0: t * dSigma0, Sigma5: t * dSigma5})
    .series(t, 0, 2)
    .removeO()
    .coeff(t, 1)
)
Delta_norm_linear = sp.simplify(P0_target * NQ_minus_1_linear)

print("linearized N_Q - 1 =")
sp.pprint(sp.expand(NQ_minus_1_linear))
print("linearized Delta_norm(point-particle) =")
sp.pprint(sp.expand(Delta_norm_linear))
expect_zero(
    "linearized N_Q - 1",
    NQ_minus_1_linear + 5 * eps_beta + dSigma0 / (3 * S) + 9 * dSigma5 / S,
)
expect_zero(
    "linearized Delta_norm(point-particle)",
    Delta_norm_linear + P0_target * (5 * eps_beta + dSigma0 / (3 * S) + 9 * dSigma5 / S),
)

# ---------------------------------------------------------------------------
# V. Canonical compact outgoing branch
# ---------------------------------------------------------------------------
subbanner("V. Canonical compact outgoing branch")

expect_zero("chi_Q(canonical) - 1", chi_from_def.subs({beta: 1, Sigma0: 0, Sigma5: 0}) - 1)
expect_zero("N_Q(canonical) - 1", NQ_from_def.subs({beta: 1, Sigma0: 0, Sigma5: 0}) - 1)
expect_zero(
    "Delta_norm(canonical)",
    Delta_norm_from_def.subs({beta: 1, Sigma0: 0, Sigma5: 0}),
)
expect_zero(
    "P0(canonical source-map-reduced) - P0_target",
    sp.simplify((NQ_from_def * P0_target).subs({beta: 1, Sigma0: 0, Sigma5: 0}) - P0_target),
)

banner("STAGE 178 LEDGER")
print("1. The isotropic retarded grouped-P2 one-pole branch carries the exact odd ratio")
print("      Gammabar_5 / Gammabar_5^target = chi_Q N_Q,")
print("   where N_Q := Pbar_0 / P0^target and P0^target = 54 G c_s^5/(5 a^5 c^5).")
print("2. The observable odd closure therefore factorizes exactly as")
print("      mhat_0^2 chi_Q N_Q = 1.")
print("3. The Packet-A normalization residual satisfies")
print("      Delta_norm = P0^target (mhat_0^2 N_Q - 1),")
print("   and after imposing the observable odd closure it collapses exactly to")
print("      Delta_norm = P0^target (1/chi_Q - 1).")
print("4. On the natural point-particle source-map branch mhat_0 -> 1, the conservative")
print("   carrier ratio itself reduces to")
print("      N_Q = 1/chi_Q,")
print("   so the last reduced isotropic defect is purely outgoing.")
print("5. Inserting the Stage 194 isotropic DtN deformation algebra gives")
print("      chi_Q = 3 (S beta^5 + 9 Sigma_5)/(3 S - Sigma_0),")
print("      N_Q = (3 S - Sigma_0)/(3 (S beta^5 + 9 Sigma_5)),")
print("      Delta_norm^(pt) = -P0^target [3 S(beta^5-1)+Sigma_0+27 Sigma_5]/[3(S beta^5+9 Sigma_5)].")
print("6. Linearized about the canonical compact passive/outgoing branch,")
print("      N_Q - 1 = -5 eps_beta - dSigma_0/(3 S) - 9 dSigma_5/S + O(2).")
print("7. On the canonical compact outgoing branch itself (beta=1, Sigma_0=Sigma_5=0, hence chi_Q=1),")
print("   one gets exactly")
print("      N_Q = 1,   Delta_norm = 0.")
print("8. Stage 195 therefore removes the last source-map ambiguity from the isotropic retarded")
print("   finish line: once the odd observable closure is imposed, the remaining Packet-A normalization")
print("   defect is already the purely outgoing scalar chi_Q-1, and on the natural source-map branch")
print("   the conservative carrier ratio itself is N_Q = 1/chi_Q.")
