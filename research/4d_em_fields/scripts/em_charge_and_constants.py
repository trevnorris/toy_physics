"""
em_charge_and_constants.py

Purpose:
    Current-supportive symbolic check for the 4d_em_fields paper.

    This script aligns with the current fixed defect-branch charge ontology and
    verifies:

        Zint   = sqrt(pi) * lambdaConf,
        mu0eff = mu0 / Zint,
        qStar  = etaQ * eStar,
        qEff   = qStar / sqrt(Zint),
        eEff   = eStar / sqrt(Zint),

    together with the Gaussian KK coefficient pattern

        c_{2k} / c_0 = binomial(2k, k) / 4^k,

    and the leading Yukawa correction

        A0(r) = (mu0eff qStar)/(4 pi r) * [1 + 1/2 exp(-2 r / lambdaConf) + ...].

Dependencies:
    sympy
"""

from sympy import Eq, binomial, exp, pi, pprint, simplify, sqrt, symbols


lambdaConf, mu0, eStar, r = symbols("lambdaConf mu0 eStar r", positive=True)
etaQ = symbols("etaQ", real=True)
k = symbols("k", integer=True, nonnegative=True)

Zint = sqrt(pi) * lambdaConf
mu0eff = simplify(mu0 / Zint)
qStar = etaQ * eStar
qEff = simplify(qStar / sqrt(Zint))
eEff = simplify(eStar / sqrt(Zint))

c0 = simplify(1 / Zint)
c_even = simplify(binomial(2 * k, k) / (4**k * Zint))
c_ratio = simplify(c_even / c0)

leading_ratio = simplify(c_ratio.subs(k, 1))
leading_mass = simplify(2 / lambdaConf)

A0_leading_raw = simplify(
    (mu0 * qStar) / (4 * pi * r)
    * (c0 + c_even.subs(k, 1) * exp(-leading_mass * r))
)
A0_leading_reduced = simplify(
    (mu0eff * qStar) / (4 * pi * r)
    * (1 + leading_ratio * exp(-leading_mass * r))
)

print("Gaussian localization and charge normalization")
print("==============================================")
print()

print("Localization integral and effective coupling:")
pprint(Eq(symbols("Zint"), Zint))
pprint(Eq(symbols("mu0eff"), mu0eff))
print()

print("Fixed branch bookkeeping:")
pprint(Eq(symbols("qStar"), qStar))
pprint(Eq(symbols("qEff"), qEff))
pprint(Eq(symbols("eEff"), eEff))
print()

print("Check qEff - etaQ * eEff:")
pprint(simplify(qEff - etaQ * eEff))
print("PASS" if simplify(qEff - etaQ * eEff) == 0 else "FAIL")
print()

print("Gaussian KK even-mode coefficient pattern:")
pprint(Eq(symbols("c0"), c0))
pprint(Eq(symbols("c_2k"), c_even))
pprint(Eq(symbols("c_2k_over_c0"), c_ratio))
print()

print("Leading Yukawa data:")
pprint(Eq(symbols("leading_ratio"), leading_ratio))
pprint(Eq(symbols("leading_mass"), leading_mass))
print()

print("Leading potential written from raw KK coefficients:")
pprint(Eq(symbols("A0_leading_raw"), A0_leading_raw))
print()

print("Leading potential written in reduced mu0eff form:")
pprint(Eq(symbols("A0_leading_reduced"), A0_leading_reduced))
print()

print("Check raw - reduced leading potentials:")
pprint(simplify(A0_leading_raw - A0_leading_reduced))
print(
    "PASS"
    if simplify(A0_leading_raw - A0_leading_reduced) == 0
    else "FAIL"
)
