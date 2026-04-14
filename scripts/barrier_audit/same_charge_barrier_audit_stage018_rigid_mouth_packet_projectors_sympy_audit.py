#!/usr/bin/env python3

from __future__ import annotations

import sympy as sp


def expect_zero(expr, name: str) -> None:
    expr = sp.simplify(sp.expand(expr))
    if isinstance(expr, sp.MatrixBase):
        if any(sp.simplify(sp.expand(v)) != 0 for v in expr):
            raise AssertionError(f"{name} is not zero:\n{expr}")
    else:
        if expr != 0:
            raise AssertionError(f"{name} is not zero: {expr}")


def main() -> None:
    c = sp.symbols("c_eta", positive=True, finite=True)
    R1, E1, qeta, L = sp.symbols("R1 E1 q_eta L", real=True)

    x = sp.Matrix([R1, E1])

    # Exact rigid-mouth direct packet map.
    Mrm = sp.Matrix([[-1, -c], [0, 1]])
    q = sp.simplify(Mrm * x)
    Xi1 = sp.simplify(q[0])
    qeta_expr = sp.simplify(q[1])

    # The map is an involution.
    expect_zero(Mrm * Mrm - sp.eye(2), "Mrm^2 - I")
    Srm = sp.simplify(Mrm.inv())
    expect_zero(Srm - Mrm, "Srm - Mrm")
    expect_zero(Srm * q - x, "inverse map")

    # Packet projectors.
    Qnt = sp.diag(1, 0)
    Qeta = sp.diag(0, 1)
    expect_zero(Qnt + Qeta - sp.eye(2), "Qnt + Qeta - I")
    expect_zero(Qnt * Qnt - Qnt, "Qnt^2 - Qnt")
    expect_zero(Qeta * Qeta - Qeta, "Qeta^2 - Qeta")
    expect_zero(Qnt * Qeta, "Qnt Qeta")
    expect_zero(Qeta * Qnt, "Qeta Qnt")

    # Direct-space projectors.
    Pnt = sp.simplify(Srm * Qnt * Mrm)
    Peta = sp.simplify(Srm * Qeta * Mrm)
    Pnt_expected = sp.Matrix([[1, c], [0, 0]])
    Peta_expected = sp.Matrix([[0, -c], [0, 1]])
    expect_zero(Pnt - Pnt_expected, "Pnt - expected")
    expect_zero(Peta - Peta_expected, "Peta - expected")
    expect_zero(Pnt + Peta - sp.eye(2), "Pnt + Peta - I")
    expect_zero(Pnt * Pnt - Pnt, "Pnt^2 - Pnt")
    expect_zero(Peta * Peta - Peta, "Peta^2 - Peta")
    expect_zero(Pnt * Peta, "Pnt Peta")
    expect_zero(Peta * Pnt, "Peta Pnt")

    x_nt = sp.simplify(Pnt * x)
    x_eta = sp.simplify(Peta * x)
    x_nt_expected = sp.Matrix([-Xi1, 0])
    x_eta_expected = sp.Matrix([-c * qeta_expr, qeta_expr])
    expect_zero(x_nt - x_nt_expected, "x_nt - expected")
    expect_zero(x_eta - x_eta_expected, "x_eta - expected")
    expect_zero(x_nt + x_eta - x, "x_nt + x_eta - x")

    # Compensated strip: q_nt = 0, q_eta free.
    x_strip = sp.simplify(Srm * sp.Matrix([0, qeta]))
    expect_zero(x_strip - sp.Matrix([-c * qeta, qeta]), "x_strip - expected")
    expect_zero(Mrm * x_strip - sp.Matrix([0, qeta]), "packet on compensated strip")
    strip_norm_sq = sp.simplify((x_strip.T * x_strip)[0])
    expect_zero(strip_norm_sq - (1 + c**2) * qeta**2, "strip norm law")

    # Codimension-two orbit-lock theorem.
    expect_zero(qeta_expr - E1, "q_eta = E1")
    expect_zero(Xi1 + R1 + c * E1, "Xi1 = -R1 - c E1")
    # Xi1 = 0 and R1 = 0 implies E1 = 0.
    E_from_Xi_R = sp.solve([sp.Eq(Xi1, 0), sp.Eq(R1, 0)], [R1, E1], dict=True)
    if E_from_Xi_R != [{R1: 0, E1: 0}]:
        raise AssertionError(f"orbit-lock solve unexpected: {E_from_Xi_R}")

    # Static-only and full orbit-lock corrections.
    dx_static = sp.simplify(-x_nt)
    dx_orbit = sp.simplify(-x)
    dx_blind = sp.simplify(dx_orbit - dx_static)
    expect_zero(dx_static - sp.Matrix([Xi1, 0]), "dx_static - expected")
    expect_zero(dx_blind - sp.Matrix([c * qeta_expr, -qeta_expr]), "dx_blind - expected")

    print("Stage 018 audit checks passed.")
    print()
    print("Rigid-mouth packet map:")
    print(f"M_rm =\n{sp.pretty(Mrm)}")
    print(f"q_rm = M_rm x_rm = {q}")
    print()
    print("Canonical direct-space projectors:")
    print(f"P_nt =\n{sp.pretty(Pnt)}")
    print(f"P_eta =\n{sp.pretty(Peta)}")
    print()
    print("Canonical decomposition:")
    print(f"x_nt  = {x_nt}")
    print(f"x_eta = {x_eta}")
    print()
    print("Compensated strip:")
    print(f"x_strip(q_eta) = {x_strip}")
    print(f"||x_strip||^2  = {strip_norm_sq}")
    print()
    print("Corrections:")
    print(f"delta x_static = {dx_static}")
    print(f"delta x_blind  = {dx_blind}")
    print(f"delta x_orbit  = {dx_orbit}")


if __name__ == "__main__":
    main()
