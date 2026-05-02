from __future__ import annotations

import sympy as sp
from sympy.physics.wigner import gaunt


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m == 0:
        return sp.Integer(1)
    same_sign = sp.simplify(gaunt(2, 2, 2, 0, m, m))
    if same_sign != 0:
        raise AssertionError(f"Real-harmonic same-sign cross term should vanish for m={m}: {same_sign}")
    return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)


def grouped_trace_anomaly(x20: sp.Expr, x21: sp.Expr, x22: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    xbar = sp.simplify((x20 + 2*x21 + 2*x22) / 5)
    ax = sp.simplify((2*x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)
    return xbar, ax, bx


def main() -> None:
    eps = sp.symbols('eps', real=True)
    M1, K1w = sp.symbols('M1 K1w', real=True)
    D0, N0 = sp.symbols('D0 N0', positive=True)

    lam20 = real_y20_square_ratio(0)
    lam21 = real_y20_square_ratio(1)
    lam22 = real_y20_square_ratio(2)
    assert_zero('Y20 overlap lane 20', lam20 - 1)
    assert_zero('Y20 overlap lane 21', lam21 - sp.Rational(1, 2))
    assert_zero('Y20 overlap lane 22', lam22 + 1)

    dM20 = eps * lam20 * M1
    dM21 = eps * lam21 * M1
    dM22 = eps * lam22 * M1

    dK20 = eps * lam20 * K1w
    dK21 = eps * lam21 * K1w
    dK22 = eps * lam22 * K1w

    Mbar, aM, bM = grouped_trace_anomaly(dM20, dM21, dM22)
    Kbar, aK, bK = grouped_trace_anomaly(dK20, dK21, dK22)
    assert_zero('wall inertia trace', Mbar)
    assert_zero('wall inertia b=3a', bM - 3*aM)
    assert_zero('wall stiffness trace', Kbar)
    assert_zero('wall stiffness b=3a', bK - 3*aK)

    # Wall-only contributions to the Stage-5 weak-axisymmetric gates.
    D01_20, D01_21, D01_22 = dK20, dK21, dK22
    D21_20, D21_21, D21_22 = -dM20, -dM21, -dM22
    D41_20, D41_21, D41_22 = sp.Integer(0), sp.Integer(0), sp.Integer(0)

    dKsym, dMsym = sp.symbols('dKsym dMsym', real=True)
    B01, B21, B41 = sp.symbols('B01 B21 B41')
    Z01, Z21, Z41 = sp.symbols('Z01 Z21 Z41')
    D01_full = dKsym - B01 - Z01
    D21_full = -(dMsym + B21 + Z21)
    D41_full = -(B41 + Z41)
    K1_full = sp.expand(D21_full + D01_full / 9)
    Hev_full = sp.expand(D41_full - sp.Rational(2, 3) * D21_full - D01_full / 27)
    wall_only = {B01: 0, B21: 0, B41: 0, Z01: 0, Z21: 0, Z41: 0}
    K1_wall = sp.expand(K1_full.subs(wall_only))
    Hev_wall = sp.expand(Hev_full.subs(wall_only))
    assert_zero("wall-only K1 specialization", K1_wall - (-dMsym + dKsym / 9))
    assert_zero("wall-only H_even specialization", Hev_wall - (sp.Rational(2, 3) * dMsym - dKsym / 27))

    K1_gate_20 = sp.simplify(D21_20 + D01_20 / 9)
    K1_gate_21 = sp.simplify(D21_21 + D01_21 / 9)
    K1_gate_22 = sp.simplify(D21_22 + D01_22 / 9)

    Hev_20 = sp.simplify(D41_20 - sp.Rational(2, 3) * D21_20 - D01_20 / 27)
    Hev_21 = sp.simplify(D41_21 - sp.Rational(2, 3) * D21_21 - D01_21 / 27)
    Hev_22 = sp.simplify(D41_22 - sp.Rational(2, 3) * D21_22 - D01_22 / 27)

    # Solve the wall-only even gates.
    wall_matrix = sp.Matrix([
        [sp.diff(K1_wall, dKsym), sp.diff(K1_wall, dMsym)],
        [sp.diff(Hev_wall, dKsym), sp.diff(Hev_wall, dMsym)],
    ])
    sol_even = sp.solve(
        [sp.Eq(K1_wall, 0), sp.Eq(Hev_wall, 0)],
        [dKsym, dMsym],
        dict=True,
    )
    assert_zero("wall-only even-gate determinant", wall_matrix.det() - sp.Rational(1, 27))
    assert_nonzero("mutated wall-only determinant should fail", wall_matrix.det() + sp.Rational(1, 27))

    # If the outgoing transfer remains isotropic at first order, Xi_load is driven only by D01.
    Xi20 = sp.simplify(-D01_20 / D0)
    Xi21 = sp.simplify(-D01_21 / D0)
    Xi22 = sp.simplify(-D01_22 / D0)
    Xibar, aXi, bXi = grouped_trace_anomaly(Xi20, Xi21, Xi22)

    # Prefactor shifts for P0 = N0 / D0 with N01 = 0.
    dP20 = sp.simplify(-N0 * D01_20 / D0**2)
    dP21 = sp.simplify(-N0 * D01_21 / D0**2)
    dP22 = sp.simplify(-N0 * D01_22 / D0**2)
    Pbar, aP, bP = grouped_trace_anomaly(dP20, dP21, dP22)
    if sol_even != [{dKsym: 0, dMsym: 0}]:
        raise AssertionError(f'Wall-only even gates should solve trivially: {sol_even}')
    assert_zero('Xi trace', Xibar)
    assert_zero('Xi b=3a', bXi - 3*aXi)
    assert_zero('prefactor trace', Pbar)
    assert_zero('prefactor b=3a', bP - 3*aP)

    lines: list[str] = []
    lines.append("Weak-axisymmetric throat-action bridge\n")
    lines.append("Assume a weak axisymmetric Y20 perturbation of the promoted wall action:\n")
    lines.append("  mu_eta -> mu_eta + eps * delta_mu(w) Y20\n")
    lines.append("  T_w    -> T_w    + eps * delta_Tw(w) Y20\n")
    lines.append("  T_Omega-> T_Omega+ eps * delta_TO(w) Y20\n")
    lines.append("  K_eta  -> K_eta  + eps * delta_Keta(w) Y20.\n")
    lines.append("The exact grouped-lane signature from the Y20 triple overlap is\n")
    lines.append("  lambda_(20) = 1, lambda_(21) = 1/2, lambda_(22) = -1.\n")
    lines.append("\n")
    lines.append("So any wall-sector shift X_A has the form X_A = X^(0) + eps * lambda_A * X1.\n")
    lines.append(f"Grouped trace/anomaly of the wall inertia shift: Mbar={sp.sstr(Mbar)}, a_M={sp.sstr(aM)}, b_M={sp.sstr(bM)}\n")
    lines.append(f"Grouped trace/anomaly of the wall stiffness shift: Kbar={sp.sstr(Kbar)}, a_K={sp.sstr(aK)}, b_K={sp.sstr(bK)}\n")
    lines.append("Hence every first-order axisymmetric wall anisotropy obeys b = 3 a exactly.\n")
    lines.append("\n")
    lines.append("Wall-only contributions to the live weak-axisymmetric gates are\n")
    lines.append(f"  K1_(20,21,22)   = ({sp.sstr(K1_gate_20)}, {sp.sstr(K1_gate_21)}, {sp.sstr(K1_gate_22)})\n")
    lines.append(f"  H_even_(20,21,22) = ({sp.sstr(Hev_20)}, {sp.sstr(Hev_21)}, {sp.sstr(Hev_22)})\n")
    lines.append("with generic lane formulas\n")
    lines.append("  K1_wall = -delta_M + delta_K/9\n")
    lines.append("  H_even,wall = 2 delta_M / 3 - delta_K / 27.\n")
    lines.append(f"Solving K1_wall = 0 and H_even,wall = 0 gives: {sol_even}\n")
    lines.append("Therefore pure wall anisotropy closes the even gates only on the trivial branch delta_K = delta_M = 0.\n")
    lines.append("\n")
    lines.append("If the outgoing-transfer bundle stays isotropic at first order (N01 = 0), then\n")
    lines.append("  Xi_load^(A) = - D01^(A) / D0 = - eps * lambda_A * K1w / D0.\n")
    lines.append(f"Its grouped defects are Xibar={sp.sstr(Xibar)}, a_Xi={sp.sstr(aXi)}, b_Xi={sp.sstr(bXi)}\n")
    lines.append("so Xi_load again obeys b = 3 a.\n")
    lines.append("\n")
    lines.append("The induced prefactor shift dP0^(A) = -N0 D01^(A) / D0^2 has grouped defects\n")
    lines.append(f"  Pbar={sp.sstr(Pbar)}, a_P={sp.sstr(aP)}, b_P={sp.sstr(bP)}\n")
    lines.append("and again lies on the same weak-axisymmetric line b = 3 a.\n")
    lines.append("\n")
    lines.append("Conclusion: promoting S_Sigma gives a parent-level origin for the grouped signature (20,21,22) ~ (1,1/2,-1), but by itself it does not solve the full weak-axisymmetric 5PN gates. The Maxwell/mixed and support sectors still have to participate.\n")
    lines.append("\nSTATUS: PASS\n")

    print("".join(lines))


if __name__ == "__main__":
    main()
