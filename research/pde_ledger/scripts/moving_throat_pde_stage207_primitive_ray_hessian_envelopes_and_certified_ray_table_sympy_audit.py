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


banner("STAGE 190 — PRIMITIVE-RAY HESSIAN ENVELOPES AND CERTIFIED RAY TABLE")

# ---------------------------------------------------------------------------
# I. Exact primitive directional Hessian reduction
# ---------------------------------------------------------------------------
subbanner("I. Exact primitive directional Hessian reduction")

Hll, Hcc, Hgg, HUU, HWW = sp.symbols("Hll Hcc Hgg HUU HWW", real=True)
Hlc, Hlg, HlU, HlW = sp.symbols("Hlc Hlg HlU HlW", real=True)
Hcg, HcU, HcW = sp.symbols("Hcg HcU HcW", real=True)
HgU, HgW, HUW = sp.symbols("HgU HgW HUW", real=True)

H = sp.Matrix(
    [
        [Hll, Hlc, Hlg, HlU, HlW],
        [Hlc, Hcc, Hcg, HcU, HcW],
        [Hlg, Hcg, Hgg, HgU, HgW],
        [HlU, HcU, HgU, HUU, HUW],
        [HlW, HcW, HgW, HUW, HWW],
    ]
)

labels = ["lambda", "c", "gamma", "U", "W"]
basis = [sp.Matrix([1, 0, 0, 0, 0]), sp.Matrix([0, 1, 0, 0, 0]), sp.Matrix([0, 0, 1, 0, 0]), sp.Matrix([0, 0, 0, 1, 0]), sp.Matrix([0, 0, 0, 0, 1])]
diag_entries = [Hll, Hcc, Hgg, HUU, HWW]

for label, e, d in zip(labels, basis, diag_entries):
    expr_plus = sp.simplify((e.T * H * e)[0] - d)
    expr_minus = sp.simplify(((-e).T * H * (-e))[0] - d)
    print(f"primitive curvature entry for {label} =")
    sp.pprint((e.T * H * e)[0])
    expect_zero(f"+e_{label} diagonal reduction", expr_plus)
    expect_zero(f"-e_{label} diagonal reduction", expr_minus)

# ---------------------------------------------------------------------------
# II. Off-diagonal Hessian entry first appears only on mixed rays
# ---------------------------------------------------------------------------
subbanner("II. Off-diagonal Hessian entry first appears only on mixed rays")

a, b = sp.symbols("a b", real=True)
ei = basis[0]  # lambda
ej = basis[1]  # c
mixed_expr = sp.simplify(((a * ei + b * ej).T * H * (a * ei + b * ej))[0])
expected_mixed = sp.simplify(a**2 * Hll + 2 * a * b * Hlc + b**2 * Hcc)
print("mixed two-coordinate curvature (lambda,c) =")
sp.pprint(mixed_expr)
expect_zero("mixed-ray quadratic form", mixed_expr - expected_mixed)
expect_zero(
    "off-diagonal first appearance",
    mixed_expr - (a**2 * Hll + b**2 * Hcc) - 2 * a * b * Hlc,
)

# ---------------------------------------------------------------------------
# III. Canonical primitive orientation gives negative slope magnitude
# ---------------------------------------------------------------------------
subbanner("III. Canonical primitive orientation gives negative slope magnitude")

Gam = sp.symbols("Gam", real=True)
eps = -sp.sign(Gam)
K_oriented = sp.simplify(eps * Gam)
print("canonical primitive sign eps =")
sp.pprint(eps)
print("oriented primitive slope =")
sp.pprint(K_oriented)
expect_zero("canonical orientation law", K_oriented + sp.Abs(Gam))

# ---------------------------------------------------------------------------
# IV. Sign-adapted primitive dependent exponents from Stage 204
# ---------------------------------------------------------------------------
subbanner("IV. Sign-adapted primitive dependent exponents from Stage 204")

chi0s, deltaUs, Estar, Fstar = sp.symbols(
    "chi0_star deltaU_star E_star F_star", positive=True, real=True
)
a_star = sp.simplify((1 + deltaUs) / (1 + chi0s))

sl, sc, sg, sU, sW = sp.symbols("s_lambda s_c s_gamma s_U s_W", real=True)

sigma_delta = sp.simplify(-a_star * (sg + sc - sU))
sigma_T = sp.simplify(sU + sigma_delta)
sigma_Keta = sp.simplify(2 * sc - sU)
sigma_mu = sp.simplify(
    2 * sc
    - sU
    + 2 * sW
    - 2 * sl
    - Estar * (2 * sg + 2 * sl - sU - sW)
    + Fstar * sigma_delta
)

primitive_specs = {
    "lambda": {
        "vec": sp.Matrix([1, 0, 0, 0, 0]),
        "expected": sp.Matrix([0, 0, 0, -2 - 2 * Estar]),
    },
    "c": {
        "vec": sp.Matrix([0, 1, 0, 0, 0]),
        "expected": sp.Matrix([-a_star, -a_star, 2, 2 - Fstar * a_star]),
    },
    "gamma": {
        "vec": sp.Matrix([0, 0, 1, 0, 0]),
        "expected": sp.Matrix([-a_star, -a_star, 0, -2 * Estar - Fstar * a_star]),
    },
    "U": {
        "vec": sp.Matrix([0, 0, 0, 1, 0]),
        "expected": sp.Matrix([a_star, 1 + a_star, -1, -1 + Estar + Fstar * a_star]),
    },
    "W": {
        "vec": sp.Matrix([0, 0, 0, 0, 1]),
        "expected": sp.Matrix([0, 0, 0, 2 + Estar]),
    },
}

svec = sp.Matrix([sl, sc, sg, sU, sW])
raw_rows = sp.Matrix([sigma_delta, sigma_T, sigma_Keta, sigma_mu])

eps_lambda, eps_c, eps_gamma, eps_U, eps_W = sp.symbols(
    "eps_lambda eps_c eps_gamma eps_U eps_W", real=True
)
eps_map = {
    "lambda": eps_lambda,
    "c": eps_c,
    "gamma": eps_gamma,
    "U": eps_U,
    "W": eps_W,
}

for label, spec in primitive_specs.items():
    eps_i = eps_map[label]
    subs_dict = dict(zip(svec, eps_i * spec["vec"]))
    row = raw_rows.subs(subs_dict)
    expected = eps_i * spec["expected"]
    print(f"sign-adapted primitive exponents for {label} =")
    sp.pprint(row)
    expect_zero(f"primitive exponent row ({label})", row - expected)

# ---------------------------------------------------------------------------
# V. Certified monotone primitive bracket compiler
# ---------------------------------------------------------------------------
subbanner("V. Certified monotone primitive bracket compiler")

H0, k, cL, cU = sp.symbols("H0 k cL cU", positive=True, real=True)
TauL = sp.simplify(2 * H0 / (k + sp.sqrt(k**2 - 2 * cL * H0)))
TauU = sp.simplify(2 * H0 / (k + sp.sqrt(k**2 - 2 * cU * H0)))
print("tau_lo(H0,k;cL) =")
sp.pprint(TauL)
print("tau_hi(H0,k;cU) =")
sp.pprint(TauU)
expect_zero(
    "lower comparison quadratic",
    sp.simplify(H0 - k * TauL + sp.Rational(1, 2) * cL * TauL**2),
)
expect_zero(
    "upper comparison quadratic",
    sp.simplify(H0 - k * TauU + sp.Rational(1, 2) * cU * TauU**2),
)
Width = sp.simplify(TauU - TauL)
print("primitive monotone bracket width =")
sp.pprint(Width)

# ---------------------------------------------------------------------------
# VI. Certified turning primitive bracket compiler
# ---------------------------------------------------------------------------
subbanner("VI. Certified turning primitive bracket compiler")

a_turn, b_turn = sp.symbols("a_turn b_turn", positive=True, real=True)
TauTurnLo = sp.sqrt(2 * H0 / a_turn)
TauTurnHi = sp.sqrt(2 * H0 / b_turn)
print("tau_lo^(tp) =")
sp.pprint(TauTurnLo)
print("tau_hi^(tp) =")
sp.pprint(TauTurnHi)
expect_zero(
    "turning lower quadratic",
    sp.simplify(H0 - sp.Rational(1, 2) * a_turn * TauTurnLo**2),
)
expect_zero(
    "turning upper quadratic",
    sp.simplify(H0 - sp.Rational(1, 2) * b_turn * TauTurnHi**2),
)

# ---------------------------------------------------------------------------
# VII. Primitive certified row packet summary
# ---------------------------------------------------------------------------
subbanner("VII. Primitive certified row packet summary")

k_lambda, k_c, k_gamma, k_U, k_W = sp.symbols(
    "k_lambda k_c k_gamma k_U k_W", positive=True, real=True
)
kappa_lo_lambda, kappa_hi_lambda = sp.symbols("kappa_lo_lambda kappa_hi_lambda", real=True)
kappa_lo_c, kappa_hi_c = sp.symbols("kappa_lo_c kappa_hi_c", real=True)
kappa_lo_gamma, kappa_hi_gamma = sp.symbols("kappa_lo_gamma kappa_hi_gamma", real=True)
kappa_lo_U, kappa_hi_U = sp.symbols("kappa_lo_U kappa_hi_U", real=True)
kappa_lo_W, kappa_hi_W = sp.symbols("kappa_lo_W kappa_hi_W", real=True)

primitive_row_data = {
    "lambda": (k_lambda, kappa_lo_lambda, kappa_hi_lambda),
    "c": (k_c, kappa_lo_c, kappa_hi_c),
    "gamma": (k_gamma, kappa_lo_gamma, kappa_hi_gamma),
    "U": (k_U, kappa_lo_U, kappa_hi_U),
    "W": (k_W, kappa_lo_W, kappa_hi_W),
}

for label, (k_i, c_lo_i, c_hi_i) in primitive_row_data.items():
    tau_lo_i = sp.simplify(2 * H0 / (k_i + sp.sqrt(k_i**2 - 2 * c_lo_i * H0)))
    tau_hi_i = sp.simplify(2 * H0 / (k_i + sp.sqrt(k_i**2 - 2 * c_hi_i * H0)))
    print(f"primitive certified row ({label}) -> tau_lo, tau_hi =")
    sp.pprint((tau_lo_i, tau_hi_i))

banner("STAGE 190 SYMPY AUDIT PASSED")
