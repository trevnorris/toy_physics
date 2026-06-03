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


banner("STAGE 194 — EXACT OUTGOING l=2 DTN FINGERPRINT, FIXING OF chi_Q, AND THE ISOTROPIC DEFORMATION ALGEBRA")

# ---------------------------------------------------------------------------
# I. Exact outgoing spherical l=2 DtN fingerprint
# ---------------------------------------------------------------------------
subbanner("I. Exact outgoing spherical l=2 DtN fingerprint")

z = sp.symbols("z", real=True)

j2 = ((sp.Integer(3) / z**3 - sp.Integer(1) / z) * sp.sin(z) - sp.Integer(3) * sp.cos(z) / z**2)
y2 = -((sp.Integer(3) / z**3 - sp.Integer(1) / z) * sp.cos(z) + sp.Integer(3) * sp.sin(z) / z**2)
h2 = sp.simplify(j2 + sp.I * y2)
Lambda_out = sp.simplify(z * sp.diff(h2, z) / h2)
Yhat_out = sp.simplify(-sp.Integer(3) / Lambda_out)

Lambda_out_series = sp.series(Lambda_out, z, 0, 8).removeO()
Yhat_out_series = sp.series(Yhat_out, z, 0, 8).removeO()

Lambda_out_expected = sp.expand(
    -3
    + z**2 / 3
    + z**4 / 9
    + sp.I * z**5 / 9
    - 2 * z**6 / 27
    - sp.I * z**7 / 27
)
Yhat_out_expected = sp.expand(
    1
    + z**2 / 9
    + 4 * z**4 / 81
    + sp.I * z**5 / 27
    - 11 * z**6 / 729
    - sp.I * z**7 / 243
)

print("h_2^(1)(z) =")
sp.pprint(h2)
print("Lambda_2^out(z) =")
sp.pprint(Lambda_out)
print("Yhat_2^out(z) =")
sp.pprint(Yhat_out)
expect_zero("Lambda_out series - expected", Lambda_out_series - Lambda_out_expected)
expect_zero("Yhat_out series - expected", Yhat_out_series - Yhat_out_expected)
expect_zero("static DtN slot + 3", sp.simplify(Lambda_out_series.subs(z, 0) + 3))

# ---------------------------------------------------------------------------
# II. Exact matching to the retarded grouped-P2 one-pole module
# ---------------------------------------------------------------------------
subbanner("II. Exact matching to the retarded grouped-P2 one-pole module")

a, c_s, omega = sp.symbols("a c_s omega", positive=True, real=True)
chi_Q = sp.symbols("chi_Q", real=True)

Omega_Q = sp.simplify(sp.Integer(3) * c_s / (2 * a))
sigma_Q_can = sp.simplify(sp.Integer(9) / (8 * Omega_Q**5))

Yret = sp.simplify(
    sp.Rational(3, 4)
    + sp.Rational(1, 4)
    / (1 - omega**2 / Omega_Q**2 - sp.I * chi_Q * sigma_Q_can * omega**5)
)
Yret_series = sp.series(Yret, omega, 0, 6).removeO()
Yret_expected = sp.expand(
    1
    + a**2 * omega**2 / (9 * c_s**2)
    + 4 * a**4 * omega**4 / (81 * c_s**4)
    + sp.I * chi_Q * a**5 * omega**5 / (27 * c_s**5)
)

Yout_omega_expected = sp.expand(Yhat_out_expected.subs(z, a * omega / c_s).series(omega, 0, 6).removeO())

print("Omega_Q =")
sp.pprint(Omega_Q)
print("sigma_Q^can =")
sp.pprint(sigma_Q_can)
print("Y_Q^ret(omega) =")
sp.pprint(Yret)
expect_zero("sigma_Q^can - 4 a^5/(27 c_s^5)", sigma_Q_can - 4 * a**5 / (27 * c_s**5))
expect_zero("Yret series - expected low-frequency form", Yret_series - Yret_expected)
expect_zero(
    "retarded module matches outgoing DtN branch when chi_Q = 1",
    sp.expand(Yret_series.subs(chi_Q, 1) - Yout_omega_expected),
)

odd_coeff_ret = sp.expand((sp.diff(Yret_series, omega, 5).subs(omega, 0)) / sp.factorial(5))
odd_coeff_out = sp.expand((sp.diff(Yout_omega_expected, omega, 5).subs(omega, 0)) / sp.factorial(5))
print("odd coefficient of retarded module =")
sp.pprint(odd_coeff_ret)
print("odd coefficient of outgoing DtN branch =")
sp.pprint(odd_coeff_out)
expect_zero("odd-coefficient matching fixes chi_Q - 1", sp.simplify((odd_coeff_ret - odd_coeff_out) * 27 * c_s**5 / (sp.I * a**5) - (chi_Q - 1)))

# ---------------------------------------------------------------------------
# III. Exact isotropic DtN deformation algebra
# ---------------------------------------------------------------------------
subbanner("III. Exact isotropic DtN deformation algebra")

S, beta = sp.symbols("S beta", nonzero=True, real=True)
Sigma0, Sigma2, Sigma4, Sigma5 = sp.symbols("Sigma_0 Sigma_2 Sigma_4 Sigma_5", real=True)

L0 = sp.simplify(-3 * S + Sigma0)
L2 = sp.simplify(S * beta**2 / 3 + Sigma2)
L4 = sp.simplify(S * beta**4 / 9 + Sigma4)
L5 = sp.simplify(S * beta**5 / 9 + Sigma5)

Ydef = sp.simplify(L0 / (L0 + L2 * z**2 + L4 * z**4 + sp.I * L5 * z**5))
Ydef_series = sp.series(Ydef, z, 0, 6).removeO()
Ydef_expected = sp.expand(1 - L2 * z**2 / L0 + (L2**2 / L0**2 - L4 / L0) * z**4 - sp.I * L5 * z**5 / L0)

print("L0 =")
sp.pprint(L0)
print("L2 =")
sp.pprint(L2)
print("L4 =")
sp.pprint(L4)
print("L5 =")
sp.pprint(L5)
print("Yhat_2^def(z) =")
sp.pprint(Ydef)
expect_zero("Ydef series - exact compiler", Ydef_series - Ydef_expected)

sol_even = sp.solve(
    [
        sp.Eq(-L2 / L0, sp.Rational(1, 9)),
        sp.Eq(L2**2 / L0**2 - L4 / L0, sp.Rational(4, 81)),
    ],
    [Sigma2, Sigma4],
    dict=True,
)[0]

Sigma2_expected = sp.simplify(-(3 * S * beta**2 - 3 * S + Sigma0) / 9)
Sigma4_expected = sp.simplify(-(3 * S * beta**4 - 3 * S + Sigma0) / 27)

print("Sigma_2 (canonical-even match) =")
sp.pprint(sol_even[Sigma2])
print("Sigma_4 (canonical-even match) =")
sp.pprint(sol_even[Sigma4])
expect_zero("Sigma_2 formula", sol_even[Sigma2] - Sigma2_expected)
expect_zero("Sigma_4 formula", sol_even[Sigma4] - Sigma4_expected)

chi_from_def = sp.simplify((-L5 / L0) / (sp.Rational(1, 27)))
chi_from_def_even = sp.simplify(chi_from_def.subs(sol_even))
chi_expected = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
chi_minus_one_expected = sp.simplify((3 * S * (beta**5 - 1) + Sigma0 + 27 * Sigma5) / (3 * S - Sigma0))

print("chi_Q from deformed DtN branch =")
sp.pprint(chi_from_def_even)
expect_zero("chi_Q deformation law", chi_from_def_even - chi_expected)
expect_zero("chi_Q - 1 deformation law", sp.simplify(chi_from_def_even - 1) - chi_minus_one_expected)

# ---------------------------------------------------------------------------
# IV. Carry-forward corollary: canonical invariant tuple on the canonical branch
# ---------------------------------------------------------------------------
subbanner("IV. Carry-forward corollary: canonical invariant tuple")

G, c = sp.symbols("G c", positive=True, real=True)
Kbar0 = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
Kbar2 = sp.simplify(Kbar0 / (4 * Omega_Q**2))
Kbar4 = sp.simplify(Kbar0 / (4 * Omega_Q**4))
Gammabar5 = sp.simplify(9 * Kbar0 / (32 * Omega_Q**5))

print("Kbar_0 =")
sp.pprint(Kbar0)
print("Kbar_2 =")
sp.pprint(Kbar2)
print("Kbar_4 =")
sp.pprint(Kbar4)
print("Gammabar_5 =")
sp.pprint(Gammabar5)
expect_zero("Kbar_2 - 6 G c_s^3/(5 a^3 c^5)", Kbar2 - 6 * G * c_s**3 / (5 * a**3 * c**5))
expect_zero("Kbar_4 - 8 G c_s/(15 a c^5)", Kbar4 - 8 * G * c_s / (15 * a * c**5))
expect_zero("Gammabar_5 - 2 G/(5 c^5)", Gammabar5 - 2 * G / (5 * c**5))

banner("STAGE 194 LEDGER")
print("1. The exact outgoing spherical l=2 DtN operator has the small-z fingerprint")
print("      Lambda_2^out = -3 + z^2/3 + z^4/9 + i z^5/9 - 2 z^6/27 - i z^7/27 + ...")
print("   and the normalized outgoing branch therefore has")
print("      Yhat_2^out = 1 + z^2/9 + 4 z^4/81 + i z^5/27 - 11 z^6/729 - i z^7/243 + ...")
print("2. Matching the carried retarded grouped-P2 one-pole module to that exact outgoing DtN")
print("   fingerprint fixes the last canonical outgoing scalar exactly:")
print("      chi_Q = 1.")
print("3. The first explicit isotropic DtN deformation family")
print("      Lambda_2^def = S Lambda_2^out(beta z) + Sigma_0 + Sigma_2 z^2 + Sigma_4 z^4 + i Sigma_5 z^5 + ...")
print("   preserves the canonical even moments iff")
print("      Sigma_2 = -(3 S beta^2 - 3 S + Sigma_0)/9,")
print("      Sigma_4 = -(3 S beta^4 - 3 S + Sigma_0)/27,")
print("   and then the exact outgoing-normalization factor is")
print("      chi_Q = 3 (S beta^5 + 9 Sigma_5) / (3 S - Sigma_0).")
print("4. Therefore the only isotropic DtN-side data that can move chi_Q while keeping the")
print("   canonical even fingerprint fixed are beta, Sigma_0, and Sigma_5.")
print("5. If one also carries the natural point-particle source-map branch, the canonical outgoing")
print("   branch returns the invariant tuple")
print("      Kbar_0 = 54 G c_s^5/(5 a^5 c^5),")
print("      Kbar_2 = 6 G c_s^3/(5 a^3 c^5),")
print("      Kbar_4 = 8 G c_s/(15 a c^5),")
print("      Gammabar_5 = 2 G/(5 c^5).")
print("6. Stage 194 therefore freezes the retarded successor to the Stage 193 conservative surface:")
print("   the canonical compact passive/outgoing branch has chi_Q = 1, and any remaining isotropic")
print("   PDE-facing defect is an explicit DtN deformation problem in (beta, Sigma_0, Sigma_5).")
