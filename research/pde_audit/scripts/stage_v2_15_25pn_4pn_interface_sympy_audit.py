#!/usr/bin/env python3
"""
Stage V2-15: 2.5PN / 4PN quadrupole-normalization interface audit.

This script verifies the algebraic bridge between:
  1. the 2.5PN Burke-Thorne / outgoing STF quadrupole coefficient, and
  2. the conservative 4PN hereditary/tail coefficient.

It is intentionally symbolic. It does not assume a solved moving-throat PDE.
It checks that the 4PN tail does not introduce a new quadrupole normalization
once the canonical STF quadrupole branch is fixed; any remaining non-P2 datum is
isolated as a scalar tail-transport factor Theta_tail.
"""

from __future__ import annotations

import sympy as sp


def check_zero(name: str, expr: sp.Expr, rows: list[str]) -> bool:
    simplified = sp.factor(sp.cancel(sp.simplify(expr)))
    ok = simplified == 0
    rows.append(f"{name}: {'PASS' if ok else 'FAIL'}")
    if not ok:
        rows.append(f"    residual = {simplified}")
    return ok


def check_equal(name: str, lhs: sp.Expr, rhs: sp.Expr, rows: list[str]) -> bool:
    return check_zero(name, lhs - rhs, rows)


def dim_mul(*dims: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(sum(d[i] for d in dims) for i in range(3))  # type: ignore[return-value]


def dim_pow(dim: tuple[int, int, int], n: int) -> tuple[int, int, int]:
    return tuple(n * x for x in dim)  # type: ignore[return-value]


def dim_div(a: tuple[int, int, int], b: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(a[i] - b[i] for i in range(3))  # type: ignore[return-value]


def fmt_dim(d: tuple[int, int, int]) -> str:
    labels = ["M", "L", "T"]
    return " ".join(f"{lab}^{pow_}" for lab, pow_ in zip(labels, d) if pow_ != 0) or "1"


def main() -> None:
    rows: list[str] = []
    rows.append("Stage V2-15: 2.5PN / 4PN interface audit")
    rows.append("=" * 62)

    G, c, c_s, M, a = sp.symbols("G c c_s M a", positive=True, nonzero=True)
    theta_tail = sp.symbols("Theta_tail", positive=True, nonzero=True)
    P_eff = sp.symbols("P_eff", nonzero=True)  # P_eff = mhat0^2 * S_port * N0/D0
    omega = sp.symbols("omega")
    P0, P2, P4 = sp.symbols("P0 P2 P4")

    gamma_GR = sp.Rational(2, 5) * G / c**5
    P_eff_target = sp.Rational(54, 5) * G * c_s**5 / (a**5 * c**5)
    Gamma5_port = a**5 / (27 * c_s**5)
    gamma_eff = P_eff * Gamma5_port

    rows.append("\nDefinitions")
    rows.append(f"gamma_GR = {gamma_GR}")
    rows.append(f"P_eff_target = {P_eff_target}")
    rows.append(f"Gamma5_port = {Gamma5_port}")
    rows.append(f"gamma_eff = {gamma_eff}")

    checks = []
    checks.append(check_equal(
        "2.5PN target converts P_eff into gamma_GR",
        gamma_eff.subs(P_eff, P_eff_target),
        gamma_GR,
        rows,
    ))

    C_tail_GR = G**2 * M / (5 * c**8)
    C_tail_from_gamma = G * M * gamma_GR / (2 * c**3)
    checks.append(check_equal(
        "GR tail coefficient equals (GM/2c^3)*gamma_GR",
        C_tail_from_gamma,
        C_tail_GR,
        rows,
    ))

    C_tail_toy = theta_tail * G * M * gamma_eff / (2 * c_s**3)
    ratio_toy = sp.cancel(C_tail_toy / C_tail_GR)
    ratio_expected = sp.cancel(theta_tail * (c / c_s) ** 3 * (gamma_eff / gamma_GR))
    checks.append(check_equal(
        "toy/GR tail ratio factorizes into tail transport times quadrupole ratio",
        ratio_toy,
        ratio_expected,
        rows,
    ))

    checks.append(check_equal(
        "on c_s=c, Theta_tail=1, and 2.5PN target, toy tail equals GR tail",
        C_tail_toy.subs({P_eff: P_eff_target, c_s: c, theta_tail: 1}),
        C_tail_GR,
        rows,
    ))

    theta_required = (c_s / c) ** 3
    target_ratio = sp.cancel(C_tail_toy.subs(P_eff, P_eff_target) / C_tail_GR)
    checks.append(check_equal(
        "required Theta_tail condition after 2.5PN target",
        target_ratio.subs(theta_tail, theta_required),
        1,
        rows,
    ))

    A2 = a**2 / (9 * c_s**2)
    A4 = 4 * a**4 / (81 * c_s**4)
    Yout = 1 + A2 * omega**2 + A4 * omega**4 + sp.I * Gamma5_port * omega**5
    Pref = P0 + P2 * omega**2 + P4 * omega**4
    product = sp.expand(Pref * Yout)
    odd_w5_coeff = sp.expand(product).coeff(omega, 5).coeff(sp.I)
    checks.append(check_equal(
        "leading i*omega^5 coefficient ignores P2 and P4",
        odd_w5_coeff,
        P0 * Gamma5_port,
        rows,
    ))
    checks.append(check_equal(
        "tail interface ignores P2",
        sp.diff(C_tail_toy, P2),
        0,
        rows,
    ))
    checks.append(check_equal(
        "tail interface ignores P4",
        sp.diff(C_tail_toy, P4),
        0,
        rows,
    ))

    delta_Q, delta_tail = sp.symbols("delta_Q delta_tail")
    residual_ratio = sp.expand((1 + delta_Q) * (1 + delta_tail) - 1)
    checks.append(check_equal(
        "fractional 4PN residual = delta_Q + delta_tail + delta_Q*delta_tail",
        residual_ratio,
        delta_Q + delta_tail + delta_Q * delta_tail,
        rows,
    ))

    dim_G = (-1, 3, -2)
    dim_c = (0, 1, -1)
    dim_c_s = dim_c
    dim_M = (1, 0, 0)
    dim_a = (0, 1, 0)
    dim_gamma = dim_div(dim_G, dim_pow(dim_c, 5))
    dim_Ctail = dim_div(dim_mul(dim_pow(dim_G, 2), dim_M), dim_pow(dim_c, 8))
    dim_Peff_target = dim_div(dim_mul(dim_G, dim_pow(dim_c_s, 5)), dim_mul(dim_pow(dim_a, 5), dim_pow(dim_c, 5)))
    dim_gamma_eff = dim_mul(dim_Peff_target, dim_div(dim_pow(dim_a, 5), dim_pow(dim_c_s, 5)))
    dim_bridge = dim_mul(dim_div(dim_mul(dim_G, dim_M), dim_pow(dim_c, 3)), dim_gamma)

    rows.append("\nDimension ledger, exponent order (M,L,T)")
    dims_expected = {
        "gamma_GR_dim": dim_gamma,
        "C_tail_GR_dim": dim_Ctail,
        "P_eff_target_dim": dim_Peff_target,
        "gamma_eff_from_target_dim": dim_gamma_eff,
        "bridge_dim_(GM/2c^3)*gamma": dim_bridge,
    }
    for key, val in dims_expected.items():
        rows.append(f"{key}: {fmt_dim(val)}")

    dim_checks = []
    dim_checks.append(dim_gamma_eff == dim_gamma)
    rows.append(f"dimension gamma_eff == gamma_GR: {'PASS' if dim_checks[-1] else 'FAIL'}")
    dim_checks.append(dim_bridge == dim_Ctail)
    rows.append(f"dimension bridge == C_tail: {'PASS' if dim_checks[-1] else 'FAIL'}")
    dim_checks.append(dim_div(dim_mul(dim_G, dim_M), dim_pow(dim_c, 3)) == (0, 0, 1))
    rows.append(f"dimension GM/c^3 is time: {'PASS' if dim_checks[-1] else 'FAIL'}")

    rows.append("\nKey symbolic consequences")
    rows.append(f"C_tail_GR = {C_tail_GR}")
    rows.append(f"C_tail_toy = {C_tail_toy}")
    rows.append(f"C_tail_toy / C_tail_GR = {sp.factor(ratio_toy)}")
    rows.append(f"target_ratio after P_eff target = {sp.factor(target_ratio)}")
    rows.append(f"Theta_tail required after P_eff target = {theta_required}")
    rows.append(f"omega^5 odd coefficient of Pref*Yout = {odd_w5_coeff}")

    rows.append("\nSummary")
    rows.append(f"symbolic_checks_total: {len(checks)}")
    rows.append(f"symbolic_checks_passed: {sum(bool(x) for x in checks)}")
    rows.append(f"dimension_checks_total: {len(dim_checks)}")
    rows.append(f"dimension_checks_passed: {sum(bool(x) for x in dim_checks)}")

    if all(checks) and all(dim_checks):
        rows.append("FINAL_STATUS: PASS with one explicit tail-transport gate Theta_tail*(c/c_s)^3 = 1")
    else:
        rows.append("FINAL_STATUS: FAIL; inspect residuals above")

    print("\n".join(rows))


if __name__ == "__main__":
    main()
