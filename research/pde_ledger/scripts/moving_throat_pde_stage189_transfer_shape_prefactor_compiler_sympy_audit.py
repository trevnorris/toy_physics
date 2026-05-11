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


banner("STAGE 172 — TRANSFER-SHAPE / OUTGOING-PREFACTOR COMPILER")

# ---------------------------------------------------------------------------
# I. Stage 188 observable packet -> transfer-shape packet
# ---------------------------------------------------------------------------
Bstar, epsetas = sp.symbols("B_star epsiloneta_star", positive=True, real=True)
dr, dn, de = sp.symbols("dln_Rtr dln_Nstar dln_epseta", real=True)
obs = sp.Matrix([dr, dn, de])

subbanner("I. Observable packet -> transfer-shape / selected-branch packet")
dlnT2 = sp.simplify(dn - Bstar * dr)
dlnOneMinus = sp.simplify(-epsetas / (1 - epsetas) * de)
dlnRtarget = sp.simplify(dlnOneMinus - dlnT2)
trf = sp.Matrix([dlnT2, dlnOneMinus, dlnRtarget])
C_obs_to_trf = sp.Matrix(
    [
        [-Bstar, 1, 0],
        [0, 0, -epsetas / (1 - epsetas)],
        [Bstar, -1, -epsetas / (1 - epsetas)],
    ]
)
print("Delta_obs^(1) =")
sp.pprint(obs)
print("C_obs->trf =")
sp.pprint(C_obs_to_trf)
print("Delta_trf^(1) =")
sp.pprint(trf)
expect_zero("C_obs->trf * Delta_obs^(1) - Delta_trf^(1)", C_obs_to_trf * obs - trf)
expect_zero("selected-branch compatibility row", trf[2] + trf[0] - trf[1])
print("rank(C_obs->trf) =", C_obs_to_trf.rank())

# Match to defect notation from Stage 188.
Theta = dr
Xi = dlnT2
Rcal = dlnRtarget
expect_zero("Xi - (dln N_* - B_* dln R_tr)", Xi - (dn - Bstar * dr))
expect_zero("R_1 - dln R_target", Rcal - dlnRtarget)
expect_zero("(R_1 + Xi_1) - dln(1-epseta)", (Rcal + Xi) - dlnOneMinus)

# ---------------------------------------------------------------------------
# II. One-port continuum transfer shape and selected-branch identity
# ---------------------------------------------------------------------------
subbanner("II. One-port continuum transfer-shape identity")
Zw, rho, OmW2, epsW, epseta, Lambda0 = sp.symbols(
    "Z_W rho OmegaW2 epsilonW epsiloneta Lambda_0", positive=True, real=True
)
Rtarget = sp.symbols("R_target", positive=True, real=True)
T2_direct = sp.simplify(Zw * (1 + rho) ** 2 / (OmW2 * (1 - epsW) ** 2))
Rtarget_formula = sp.simplify(Lambda0 * OmW2 * (1 - epseta) * (1 - epsW) ** 2 / (Zw * (1 + rho) ** 2))
T2_selected = sp.simplify(Lambda0 * (1 - epseta) / Rtarget_formula)
print("T_A^2 (direct continuum form) =")
sp.pprint(T2_direct)
print("R_target (selected-branch form) =")
sp.pprint(Rtarget_formula)
expect_zero("Lambda_0 (1-epseta) / R_target - T_A^2", T2_selected - T2_direct)

chi0, deltaU = sp.symbols("chi0 deltaU", positive=True, real=True)
epssplit = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
T2_coh = sp.simplify(Zw * (1 + chi0) ** 2 / (OmW2 * (1 - epssplit) ** 2))
print("T^2 (coherent local D/N specialization) =")
sp.pprint(T2_coh)

# ---------------------------------------------------------------------------
# III. Exact response/prefactor compiler
# ---------------------------------------------------------------------------
subbanner("III. Exact isotropic grouped response / prefactor compiler")
omega = sp.symbols("omega", real=True)
D0, D2, D4, N0, N2, N4 = sp.symbols("D_0 D_2 D_4 N_0 N_2 N_4", nonzero=True, real=True)
Dcons = D0 + D2 * omega**2 + D4 * omega**4
Nseries = N0 + N2 * omega**2 + N4 * omega**4

Y_series = sp.series(D0 / Dcons, omega, 0, 6).removeO()
u2 = sp.simplify(-D2 / D0)
u4 = sp.simplify((D2**2 - D0 * D4) / D0**2)
print("Y(omega) = D0 / Dcons =")
sp.pprint(Y_series)
expect_zero("Y(omega) - (1 + u2 w^2 + u4 w^4)", Y_series - (1 + u2 * omega**2 + u4 * omega**4))

Pref_series = sp.series(D0 * Nseries / Dcons**2, omega, 0, 6).removeO()
P0 = sp.simplify(N0 / D0)
P2 = sp.simplify((D0 * N2 - 2 * D2 * N0) / D0**2)
P4 = sp.simplify((D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3)
print("Pref(omega) = D0 N(omega) / Dcons(omega)^2 =")
sp.pprint(Pref_series)
expect_zero("Pref(omega) - (P0 + P2 w^2 + P4 w^4)", Pref_series - (P0 + P2 * omega**2 + P4 * omega**4))

Kbl = sp.symbols("K_bl", nonzero=True, real=True)
Teff2 = sp.simplify(N0 / Kbl)
expect_zero("P0 - (K_bl/D0) T_eff^2", P0 - (Kbl / D0) * Teff2)

# ---------------------------------------------------------------------------
# IV. Weak-axisymmetric prefactor slope
# ---------------------------------------------------------------------------
subbanner("IV. Weak-axisymmetric prefactor slope")
eps, lam, D1, N1 = sp.symbols("eps lambda D_1 N_1", real=True)
P_A_series = sp.series((N0 + eps * lam * N1) / (D0 + eps * lam * D1), eps, 0, 2).removeO()
P1 = sp.simplify((N1 * D0 - N0 * D1) / D0**2)
expect_zero("P_A - (P0 + eps lambda P1)", P_A_series - (P0 + eps * lam * P1))
ln_ratio = sp.series(sp.log(P_A_series / P0), eps, 0, 2).removeO()
expect_zero("log(P_A/P0) - eps lambda (P1/P0)", ln_ratio - eps * lam * (P1 / P0))
expect_zero("P1/P0 - (N1/N0 - D1/D0)", sp.simplify(P1 / P0 - (N1 / N0 - D1 / D0)))

# ---------------------------------------------------------------------------
# V. Compact outgoing l=2 fingerprint compiler
# ---------------------------------------------------------------------------
subbanner("V. Compact outgoing l=2 fingerprint compiler")
a, cs = sp.symbols("a c_s", positive=True, real=True)
Aout = sp.simplify(a**2 / (9 * cs**2))
Bout = sp.simplify(4 * a**4 / (81 * cs**4))
G5out = sp.simplify(a**5 / (27 * cs**5))
Yhat = 1 + Aout * omega**2 + Bout * omega**4 + sp.I * G5out * omega**5
Pref_trunc = P0 + P2 * omega**2 + P4 * omega**4
out_series = sp.expand(Pref_trunc * Yhat)
out_series = sp.Poly(out_series, omega).as_dict()
# reconstruct up to O(omega^5)
out_trunc = sum(coef * omega**power[0] for power, coef in out_series.items() if power[0] <= 5)
K0 = P0
K2 = sp.simplify(P2 + Aout * P0)
K4 = sp.simplify(P4 + Aout * P2 + Bout * P0)
Gamma5 = sp.simplify(G5out * P0)
expect_zero(
    "outgoing branch expansion - (K0 + K2 w^2 + K4 w^4 + i Gamma5 w^5)",
    out_trunc - (K0 + K2 * omega**2 + K4 * omega**4 + sp.I * Gamma5 * omega**5),
)

# ---------------------------------------------------------------------------
# VI. Normalization equivalence and constant-prefactor branch
# ---------------------------------------------------------------------------
subbanner("VI. Normalization equivalence and constant-prefactor branch")
G, c, mhat0 = sp.symbols("G c mhat_0", positive=True, real=True)
P0_target = sp.simplify(54 * G * cs**5 / (5 * a**5 * c**5 * mhat0**2))
expect_zero(
    "mhat0^2 Gamma5 - 2G/(5c^5) at P0_target",
    (mhat0**2 * Gamma5 - 2 * G / (5 * c**5)).subs(P0, P0_target),
)
expect_zero(
    "Gamma5 - a^5 P0/(27 c_s^5)",
    Gamma5 - a**5 * P0 / (27 * cs**5),
)

N2_const = sp.simplify(2 * D2 * N0 / D0)
N4_const = sp.simplify((D2**2 + 2 * D0 * D4) * N0 / D0**2)
expect_zero("P2 on constant-prefactor branch", P2.subs(N2, N2_const))
expect_zero("P4 on constant-prefactor branch", P4.subs({N2: N2_const, N4: N4_const}))
expect_zero("K2 on constant-prefactor branch - A P0", K2.subs({N2: N2_const}) - Aout * P0)
expect_zero("K4 on constant-prefactor branch - B P0", K4.subs({N2: N2_const, N4: N4_const}) - Bout * P0)

banner("STAGE 172 LEDGER")
print("1. The Stage 188 branch-observable packet compiles exactly to the direct transfer packet")
print("      (d ln T^2, d ln(1-epseta), d ln R_target),")
print("   with one built-in selected-branch compatibility relation.")
print("2. The isotropic grouped operator/transfer data")
print("      (D0,D2,D4,N0,N2,N4)")
print("   compile exactly to")
print("      (u2,u4,P0,P2,P4).")
print("3. Multiplying by the compact outgoing l=2 fingerprint gives")
print("      (K0,K2,K4,Gamma5),")
print("   and the first odd coefficient depends only on the static prefactor P0.")
print("4. The universal normalization target is exactly equivalent to")
print("      mhat0^2 P0 = 54 G c_s^5 / (5 a^5 c^5).")
print("5. On the constant-prefactor branch, the even outgoing coefficients collapse to the")
print("   pure compact-fingerprint values K2 = A P0 and K4 = B P0.")
