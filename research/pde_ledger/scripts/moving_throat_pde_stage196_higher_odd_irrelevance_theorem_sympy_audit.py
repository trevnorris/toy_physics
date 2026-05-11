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


banner("STAGE 179 — EXACT HIGHER-ODD IRRELEVANCE THEOREM")

# ---------------------------------------------------------------------------
# I. Response-side higher-odd difference identity
# ---------------------------------------------------------------------------
subbanner("I. Response-side higher-odd difference identity")

omega = sp.symbols("omega", real=True)
Omega_Q, chi_Q, tau_Q = sp.symbols("Omega_Q chi_Q tau_Q", nonzero=True, real=True)

sigma_can = sp.simplify(sp.Rational(9, 8) / Omega_Q**5)
X_Q = sp.simplify(omega**2 / Omega_Q**2 + sp.I * chi_Q * sigma_can * omega**5)
H_Q = sp.I * tau_Q * omega**7

Y5 = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - X_Q))
Y7 = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - X_Q - H_Q))
Ydiff_expected = sp.simplify(H_Q / (4 * (1 - X_Q) * (1 - X_Q - H_Q)))

print("sigma_Q^can =")
sp.pprint(sigma_can)
print("Yhat_Q^(ret,5) =")
sp.pprint(Y5)
print("Yhat_Q^(ret,>=7) =")
sp.pprint(Y7)
expect_zero("exact response-side difference identity", sp.simplify(Y7 - Y5) - Ydiff_expected)

Y5_series = sp.series(Y5, omega, 0, 8).removeO()
Y7_series = sp.series(Y7, omega, 0, 8).removeO()
Ydiff_series = sp.series(sp.simplify(Y7 - Y5), omega, 0, 8).removeO()

print("Yhat_Q^(ret,5) series through O(omega^7) =")
sp.pprint(sp.expand(Y5_series))
print("Yhat_Q^(ret,>=7) series through O(omega^7) =")
sp.pprint(sp.expand(Y7_series))
print("difference series through O(omega^7) =")
sp.pprint(sp.expand(Ydiff_series))

expect_zero("higher-odd difference through O(omega^5)", sp.series(sp.simplify(Y7 - Y5), omega, 0, 7).removeO())
expected_Y5_through5 = sp.expand(1 + omega**2 / (4 * Omega_Q**2) + omega**4 / (4 * Omega_Q**4) + sp.I * chi_Q * 9 * omega**5 / (32 * Omega_Q**5))
expect_zero(
    "Yhat_Q^(ret,>=7) through O(omega^5) - canonical one-pole form",
    sp.series(Y7, omega, 0, 6).removeO() - expected_Y5_through5,
)
expect_zero(
    "first higher-odd correction really starts at omega^7",
    Ydiff_series - sp.I * tau_Q * omega**7 / 4,
)

# ---------------------------------------------------------------------------
# II. DtN-side higher-odd difference identity
# ---------------------------------------------------------------------------
subbanner("II. DtN-side higher-odd difference identity")

z = sp.symbols("z", real=True)
L0, L2, L4, L5, L7 = sp.symbols("L0 L2 L4 L5 L7", nonzero=True, real=True)

D5 = sp.simplify(L0 + L2 * z**2 + L4 * z**4 + sp.I * L5 * z**5)
L_ge7 = sp.I * L7 * z**7
D7 = sp.simplify(D5 + L_ge7)

Ydef5 = sp.simplify(L0 / D5)
Ydef7 = sp.simplify(L0 / D7)
Ydef_diff_expected = sp.simplify(-L0 * L_ge7 / (D5 * D7))

print("Yhat_2^(def,5) =")
sp.pprint(Ydef5)
print("Yhat_2^(def,>=7) =")
sp.pprint(Ydef7)
expect_zero("exact DtN-side difference identity", sp.simplify(Ydef7 - Ydef5) - Ydef_diff_expected)
expect_zero("DtN higher-odd difference through O(z^5)", sp.series(sp.simplify(Ydef7 - Ydef5), z, 0, 7).removeO())

Ydef7_series = sp.series(Ydef7, z, 0, 8).removeO()
print("Yhat_2^(def,>=7) series through O(z^7) =")
sp.pprint(sp.expand(Ydef7_series))
expect_zero(
    "DtN normalized response through O(z^5)",
    sp.series(Ydef7, z, 0, 6).removeO()
    - (1 - L2 * z**2 / L0 + (L2**2 / L0**2 - L4 / L0) * z**4 - sp.I * L5 * z**5 / L0),
)

# ---------------------------------------------------------------------------
# III. Stage 194 canonical-even matching and chi_Q stability
# ---------------------------------------------------------------------------
subbanner("III. Stage 194 canonical-even matching and chi_Q stability")

S, beta = sp.symbols("S beta", nonzero=True, real=True)
Sigma0, Sigma2, Sigma4, Sigma5 = sp.symbols("Sigma_0 Sigma_2 Sigma_4 Sigma_5", real=True)

L0_stage194 = -3 * S + Sigma0
L2_stage194 = S * beta**2 / 3 + Sigma2
L4_stage194 = S * beta**4 / 9 + Sigma4
L5_stage194 = S * beta**5 / 9 + Sigma5

Y_stage194_hi = sp.simplify(
    L0_stage194 / (L0_stage194 + L2_stage194 * z**2 + L4_stage194 * z**4 + sp.I * L5_stage194 * z**5 + sp.I * L7 * z**7)
)
Y_stage194_hi_series = sp.series(Y_stage194_hi, z, 0, 8).removeO()

Sigma2_match = sp.simplify(-(3 * S * beta**2 - 3 * S + Sigma0) / 9)
Sigma4_match = sp.simplify(-(3 * S * beta**4 - 3 * S + Sigma0) / 27)

Y_stage194_matched = sp.simplify(Y_stage194_hi_series.subs({Sigma2: Sigma2_match, Sigma4: Sigma4_match}))
Y_stage194_through5 = sp.series(Y_stage194_matched, z, 0, 6).removeO()
z5_coeff = sp.simplify(sp.expand(Y_stage194_matched).coeff(z, 5))
chi_from_series = sp.simplify(-27 * sp.I * z5_coeff)

chi_from_def = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
expected_through5 = sp.expand(1 + z**2 / 9 + 4 * z**4 / 81 + sp.I * chi_from_def * z**5 / 27)

print("Sigma_2 matching law =")
sp.pprint(Sigma2_match)
print("Sigma_4 matching law =")
sp.pprint(Sigma4_match)
print("matched normalized DtN response through O(z^7) =")
sp.pprint(sp.expand(Y_stage194_matched))
print("z^5 coefficient in the matched normalized DtN response =")
sp.pprint(z5_coeff)
print("chi_Q extracted from the matched z^5 coefficient =")
sp.pprint(chi_from_series)
print("chi_Q from the matched Stage 194 deformation algebra =")
sp.pprint(chi_from_def)

expect_zero("canonical-even plus chi_Q matching through O(z^5)", sp.simplify(Y_stage194_through5 - expected_through5))
expect_zero("chi_Q extractor - deformation algebra formula", chi_from_series - chi_from_def)
expect_zero("d chi_Q^(series extractor) / dL7", sp.diff(chi_from_series, L7))

# ---------------------------------------------------------------------------
# IV. Stage 195 source-map reduction is unchanged
# ---------------------------------------------------------------------------
subbanner("IV. Stage 195 source-map reduction is unchanged")

a, c_s, c, G = sp.symbols("a c_s c G", positive=True, real=True)
P0_target = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
N_Q_pt = sp.simplify(1 / chi_from_series)
Delta_norm_pt = sp.simplify(P0_target * (N_Q_pt - 1))

print("N_Q on the natural point-particle source-map branch =")
sp.pprint(N_Q_pt)
print("Delta_norm on the natural point-particle source-map branch =")
sp.pprint(Delta_norm_pt)
expect_zero("d N_Q / dL7", sp.diff(N_Q_pt, L7))
expect_zero("d Delta_norm / dL7", sp.diff(Delta_norm_pt, L7))
expect_zero("Delta_norm - P0_target*(1/chi_Q - 1)", Delta_norm_pt - P0_target * (1 / chi_from_series - 1))

# ---------------------------------------------------------------------------
# V. Representative higher-odd correction enters only at z^7 / omega^7
# ---------------------------------------------------------------------------
subbanner("V. Representative first higher-odd coefficient")

L7_coeff_in_Y = sp.simplify(sp.expand(Ydef7_series).coeff(z, 7).coeff(L7, 1))
print("coefficient multiplying L7 at z^7 in the normalized DtN response =")
sp.pprint(L7_coeff_in_Y)
expect_zero("L7 enters normalized response first at z^7", L7_coeff_in_Y + sp.I / L0)

banner("STAGE 179 LEDGER")
print("1. For the generalized isotropic grouped-P2 retarded module")
print("      Yhat_Q^(ret,>=7) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i chi_Q sigma_Q^can omega^5 - i T_{>=7}(omega))")
print("   with T_{>=7}(omega)=O(omega^7), one has the exact difference identity")
print("      Yhat_Q^(ret,>=7) - Yhat_Q^(ret,5) = H_Q/[4(1-X_Q)(1-X_Q-H_Q)] = O(omega^7).")
print("   So every coefficient through O(omega^5) is unchanged.")
print("2. At the DtN level, if")
print("      D_{>=7}(z) = L0 + L2 z^2 + L4 z^4 + i L5 z^5 + L_{>=7}(z),   L_{>=7}(z)=O(z^7),")
print("   then")
print("      Yhat_2^(def,>=7) - Yhat_2^(def,5) = -L0 L_{>=7}/(D5 D_{>=7}) = O(z^7).")
print("   So the normalized DtN response through O(z^5) depends only on (L0,L2,L4,L5).")
print("3. After imposing the Stage 194 canonical-even matching, the same exact formulas survive:")
print("      Sigma_2 = -(3 S beta^2 - 3 S + Sigma_0)/9,   Sigma_4 = -(3 S beta^4 - 3 S + Sigma_0)/27,")
print("      chi_Q = 3(S beta^5 + 9 Sigma_5)/(3 S - Sigma_0),")
print("   and chi_Q is completely independent of all higher odd DtN data beginning at O(z^7).")
print("4. Therefore the Stage 195 source-map-reduced Packet-A finish line is unchanged:")
print("      mhat_0^2 chi_Q N_Q = 1,   N_Q = 1/chi_Q on the natural point-particle branch,   Delta_norm = P0^target(1/chi_Q - 1).")
print("   No extra higher odd datum re-enters Delta_norm at 2.5PN order.")
print("5. The only live retarded obstruction at reduced point-particle 2.5PN order is therefore")
print("      Delta_Q = chi_Q - 1.")
print("   Any extra isotropic retarded structure beginning at O(omega^7) or higher lives beyond the theorem order.")
