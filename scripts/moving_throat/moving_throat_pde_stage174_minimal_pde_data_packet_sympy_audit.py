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


def grouped_inverse(xbar, ax, bx):
    x20 = sp.simplify(xbar + 4 * ax)
    x21 = sp.simplify(xbar - ax + bx)
    x22 = sp.simplify(xbar - ax - bx)
    return x20, x21, x22


def response_moments_from_D(D0, D2, D4):
    nu2 = sp.simplify(-D2 / D0)
    nu4 = sp.simplify((D2**2 - D0 * D4) / D0**2)
    return nu2, nu4


def prefactor_moments(D0, D2, D4, N0, N2, N4):
    P0 = sp.simplify(N0 / D0)
    P2 = sp.simplify((D0 * N2 - 2 * D2 * N0) / D0**2)
    P4 = sp.simplify(
        (D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0)
        / D0**3
    )
    return P0, P2, P4


banner("STAGE 174 — MINIMAL PDE DATA PACKET AND THE EXACT HOME-STRETCH THEOREM")

# ---------------------------------------------------------------------------
# I. Exact single-lane compiler from (D0,D2,D4,N0,N2,N4)
# ---------------------------------------------------------------------------
subbanner("I. Exact response / prefactor compiler from one grouped lane")

omega = sp.symbols("omega", real=True)
D0, D2, D4, N0, N2, N4 = sp.symbols("D0 D2 D4 N0 N2 N4", nonzero=True, real=True)

nu2, nu4 = response_moments_from_D(D0, D2, D4)
P0, P2, P4 = prefactor_moments(D0, D2, D4, N0, N2, N4)

Y = sp.simplify(D0 / (D0 + D2 * omega**2 + D4 * omega**4))
Y_series = sp.series(Y, omega, 0, 6).removeO()
Y_expected = sp.expand(1 + nu2 * omega**2 + nu4 * omega**4)
print("Y(omega) = D0 / D(omega) =")
sp.pprint(Y)
print("nu2 =")
sp.pprint(nu2)
print("nu4 =")
sp.pprint(nu4)
expect_zero("response series compiler", Y_series - Y_expected)

Pref = sp.simplify(D0 * (N0 + N2 * omega**2 + N4 * omega**4) / (D0 + D2 * omega**2 + D4 * omega**4) ** 2)
Pref_series = sp.series(Pref, omega, 0, 6).removeO()
Pref_expected = sp.expand(P0 + P2 * omega**2 + P4 * omega**4)
print("Pref(omega) = D0 N(omega) / D(omega)^2 =")
sp.pprint(Pref)
print("P0 =")
sp.pprint(P0)
print("P2 =")
sp.pprint(P2)
print("P4 =")
sp.pprint(P4)
expect_zero("prefactor series compiler", Pref_series - Pref_expected)

Delta_pole_one = sp.simplify(nu4 - 4 * nu2**2)
print("Delta_pole^(one-lane) =")
sp.pprint(Delta_pole_one)
expect_zero(
    "Delta_pole^(one-lane) + (3 D2^2 + D0 D4)/D0^2",
    Delta_pole_one + (3 * D2**2 + D0 * D4) / D0**2,
)

# ---------------------------------------------------------------------------
# II. Exact grouped weighted trace/anomaly calculus
# ---------------------------------------------------------------------------
subbanner("II. Exact grouped weighted trace/anomaly calculus")

x20, x21, x22 = sp.symbols("x20 x21 x22", real=True)
xbar, ax, bx = grouped_trace_anomaly(x20, x21, x22)
x20_rec, x21_rec, x22_rec = grouped_inverse(xbar, ax, bx)

print("xbar =")
sp.pprint(xbar)
print("a_x =")
sp.pprint(ax)
print("b_x =")
sp.pprint(bx)
expect_zero("x20 inverse", x20_rec - x20)
expect_zero("x21 inverse", x21_rec - x21)
expect_zero("x22 inverse", x22_rec - x22)

Ggrp = sp.diag(1, 2, 2)
ebar = sp.Matrix([1, 1, 1])
ea = sp.Matrix([4, -1, -1])
eb = sp.Matrix([0, 1, -1])
Pbar = sp.Matrix([[1, 2, 2], [1, 2, 2], [1, 2, 2]]) / 5
Pa = sp.Matrix([[16, -8, -8], [-4, 2, 2], [-4, 2, 2]]) / 20
Pb = sp.Matrix([[0, 0, 0], [0, 2, -2], [0, -2, 2]]) / 4

expect_zero("ebar^T G ea", (ebar.T * Ggrp * ea)[0])
expect_zero("ebar^T G eb", (ebar.T * Ggrp * eb)[0])
expect_zero("ea^T G eb", (ea.T * Ggrp * eb)[0])
expect_zero("Pbar + Pa + Pb - I", Pbar + Pa + Pb - sp.eye(3))
expect_zero("Pbar^2 - Pbar", Pbar * Pbar - Pbar)
expect_zero("Pa^2 - Pa", Pa * Pa - Pa)
expect_zero("Pb^2 - Pb", Pb * Pb - Pb)
expect_zero("Pbar Pa", Pbar * Pa)
expect_zero("Pbar Pb", Pbar * Pb)
expect_zero("Pa Pb", Pa * Pb)

# ---------------------------------------------------------------------------
# III. Packet A -> Delta_branch
# ---------------------------------------------------------------------------
subbanner("III. Exact compilation of Delta_branch from Packet A")

D20_0, D20_2, D20_4, N20_0, N20_2, N20_4 = sp.symbols(
    "D20_0 D20_2 D20_4 N20_0 N20_2 N20_4", nonzero=True, real=True
)
D21_0, D21_2, D21_4, N21_0, N21_2, N21_4 = sp.symbols(
    "D21_0 D21_2 D21_4 N21_0 N21_2 N21_4", nonzero=True, real=True
)
D22_0, D22_2, D22_4, N22_0, N22_2, N22_4 = sp.symbols(
    "D22_0 D22_2 D22_4 N22_0 N22_2 N22_4", nonzero=True, real=True
)
G, c, cs, a, mhat0 = sp.symbols("G c c_s a mhat_0", positive=True, real=True)

nu2_20, nu4_20 = response_moments_from_D(D20_0, D20_2, D20_4)
nu2_21, nu4_21 = response_moments_from_D(D21_0, D21_2, D21_4)
nu2_22, nu4_22 = response_moments_from_D(D22_0, D22_2, D22_4)

P0_20, P2_20, P4_20 = prefactor_moments(D20_0, D20_2, D20_4, N20_0, N20_2, N20_4)
P0_21, P2_21, P4_21 = prefactor_moments(D21_0, D21_2, D21_4, N21_0, N21_2, N21_4)
P0_22, P2_22, P4_22 = prefactor_moments(D22_0, D22_2, D22_4, N22_0, N22_2, N22_4)

ubar2, a2, b2 = grouped_trace_anomaly(nu2_20, nu2_21, nu2_22)
ubar4, a4, b4 = grouped_trace_anomaly(nu4_20, nu4_21, nu4_22)
Pbar0, aP0, bP0 = grouped_trace_anomaly(P0_20, P0_21, P0_22)

Delta_pole = sp.simplify(ubar4 - 4 * ubar2**2)
Delta_norm = sp.simplify(mhat0**2 * Pbar0 - 54 * G * cs**5 / (5 * a**5 * c**5))
Delta_branch = sp.Matrix([a2, b2, a4, b4, aP0, bP0, Delta_pole, Delta_norm])

print("Delta_branch =")
sp.pprint(Delta_branch)

# Isotropic collapse check
D0c, D2c, D4c, N0c, N2c, N4c = sp.symbols("D0c D2c D4c N0c N2c N4c", nonzero=True, real=True)
iso_subs = {
    D20_0: D0c, D20_2: D2c, D20_4: D4c, N20_0: N0c, N20_2: N2c, N20_4: N4c,
    D21_0: D0c, D21_2: D2c, D21_4: D4c, N21_0: N0c, N21_2: N2c, N21_4: N4c,
    D22_0: D0c, D22_2: D2c, D22_4: D4c, N22_0: N0c, N22_2: N2c, N22_4: N4c,
}
expect_zero("a2 on isotropic branch", a2.subs(iso_subs))
expect_zero("b2 on isotropic branch", b2.subs(iso_subs))
expect_zero("a4 on isotropic branch", a4.subs(iso_subs))
expect_zero("b4 on isotropic branch", b4.subs(iso_subs))
expect_zero("aP0 on isotropic branch", aP0.subs(iso_subs))
expect_zero("bP0 on isotropic branch", bP0.subs(iso_subs))

nu2_iso, nu4_iso = response_moments_from_D(D0c, D2c, D4c)
P0_iso, _, _ = prefactor_moments(D0c, D2c, D4c, N0c, N2c, N4c)
expect_zero("ubar2 on isotropic branch - nu2_iso", ubar2.subs(iso_subs) - nu2_iso)
expect_zero("ubar4 on isotropic branch - nu4_iso", ubar4.subs(iso_subs) - nu4_iso)
expect_zero("Pbar0 on isotropic branch - P0_iso", Pbar0.subs(iso_subs) - P0_iso)

# Isotropic one-pole + normalization target gives Delta_branch = 0
P0_target = sp.simplify(54 * G * cs**5 / (5 * a**5 * c**5) / mhat0**2)
onepole_norm_subs = {
    D4c: -3 * D2c**2 / D0c,
    N0c: D0c * P0_target,
}
Delta_branch_iso_test = sp.simplify(Delta_branch.subs(iso_subs).subs(onepole_norm_subs))
expect_zero("Delta_branch on isotropic one-pole normalized branch", Delta_branch_iso_test)

# ---------------------------------------------------------------------------
# IV. Exact orbit-packet interconversion
# ---------------------------------------------------------------------------
subbanner("IV. Exact orbit-packet interconversion and zero-set equivalence")

chi0s, Fstar = sp.symbols("chi0_star F_star", positive=True, real=True)
mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)
Rtr_orb, Rnt_orb, Reta_orb = sp.symbols("Rtr_orb Rnt_orb Reta_orb", positive=True, real=True)
qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)

R_from_m = sp.Matrix([
    sp.simplify(mT ** (1 + chi0s)),
    sp.simplify(mMu / (mK * mT**Fstar)),
    sp.simplify(1 / mK),
])
q_from_R = sp.Matrix([
    sp.simplify(sp.log(Rtr_orb)),
    sp.simplify(sp.log(Rnt_orb)),
    sp.simplify(sp.log(Reta_orb)),
])
m_from_q = sp.Matrix([
    sp.simplify(sp.exp(qtr / (1 + chi0s))),
    sp.simplify(sp.exp(-qeta)),
    sp.simplify(sp.exp(qnt - qeta + Fstar * qtr / (1 + chi0s))),
])
R_from_q = sp.Matrix([
    sp.simplify(sp.exp(qtr)),
    sp.simplify(sp.exp(qnt)),
    sp.simplify(sp.exp(qeta)),
])
q_from_m = sp.Matrix([
    sp.simplify((1 + chi0s) * sp.log(mT)),
    sp.simplify(sp.log(mMu) - sp.log(mK) - Fstar * sp.log(mT)),
    sp.simplify(-sp.log(mK)),
])

print("R from mismatch packet =")
sp.pprint(R_from_m)
print("q from mismatch packet =")
sp.pprint(q_from_m)
print("m from quotient packet =")
sp.pprint(m_from_q)
print("R from quotient packet =")
sp.pprint(R_from_q)

expect_zero("R_from_m after m_from_q", R_from_m.subs({mT: m_from_q[0], mK: m_from_q[1], mMu: m_from_q[2]}) - R_from_q)
expect_zero("q_from_m after m_from_q", q_from_m.subs({mT: m_from_q[0], mK: m_from_q[1], mMu: m_from_q[2]}) - sp.Matrix([qtr, qnt, qeta]))
expect_zero("q_from_R after R_from_q", q_from_R.subs({Rtr_orb: R_from_q[0], Rnt_orb: R_from_q[1], Reta_orb: R_from_q[2]}) - sp.Matrix([qtr, qnt, qeta]))
expect_zero("m_from_q at orbit lock", m_from_q.subs({qtr: 0, qnt: 0, qeta: 0}) - sp.Matrix([1, 1, 1]))
expect_zero("R_from_q at orbit lock", R_from_q.subs({qtr: 0, qnt: 0, qeta: 0}) - sp.Matrix([1, 1, 1]))
expect_zero("q_from_m at orbit lock", q_from_m.subs({mT: 1, mK: 1, mMu: 1}))

Delta_orbit = sp.Matrix([qtr, qnt, qeta])
print("Delta_orbit =")
sp.pprint(Delta_orbit)

# ---------------------------------------------------------------------------
# V. Final packet theorem ledger
# ---------------------------------------------------------------------------
subbanner("V. Final packet theorem ledger")

print("Packet A — grouped bundle data:")
print("  (D_A0,D_A2,D_A4,N_A0,N_A2,N_A4) for A in {20,21,22}, plus mhat_0")
print()
print("Packet B — any one of:")
print("  (m_T,m_K,m_mu),  (R_tr,R_nt,R_eta),  or  (q_tr,q_nt,q_eta)")
print()
print("Everything else in the reduced endgame is an exact compiler output of these packets.")
print("The reduced verdict is exactly the pair of zero-set tests:")
print("  Delta_branch = 0")
print("  Delta_orbit  = 0")

banner("STAGE 174 LEDGER")
print("1. Packet A compiles exactly to the grouped response moments, the outgoing")
print("   prefactor moments, and the final branch residual packet Delta_branch.")
print("2. Packet B can be provided in any of three exact forms and compiles exactly")
print("   to the orbit residual packet Delta_orbit = (q_tr,q_nt,q_eta).")
print("3. On the isotropic one-pole normalized test branch, Delta_branch vanishes exactly.")
print("4. On orbit lock, Delta_orbit vanishes exactly.")
print("5. So the completed moving-throat PDE only has to supply the finite data needed")
print("   to evaluate Delta_branch and Delta_orbit on the actual branch.")
