#!/usr/bin/env python3
"""
Stage V2-11: grouped real P2 projector calculus audit.

Symbolically verifies the exact weighted projector calculus for the grouped real
P2 lanes (20,21,22), the weak-axisymmetric splitting fingerprint, and the
first-order transport formulas used by later normalization stages.
"""
from __future__ import annotations

import sympy as sp

PASS = "PASS"
FAIL = "FAIL"


def smat(M: sp.Matrix) -> sp.Matrix:
    return M.applyfunc(lambda z: sp.factor(sp.cancel(z)))


def sexpr(z: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(z))


def is_zero(obj) -> bool:
    if isinstance(obj, sp.MatrixBase):
        return all(sexpr(e) == 0 for e in obj)
    return sexpr(obj) == 0


def check(name: str, expr, expected) -> tuple[str, str, object]:
    residual = expr - expected
    if isinstance(residual, sp.MatrixBase):
        residual = smat(residual)
    else:
        residual = sexpr(residual)
    return (name, PASS if is_zero(residual) else FAIL, residual)


def status(name: str, ok: bool, detail="") -> tuple[str, str, object]:
    return (name, PASS if ok else FAIL, detail)


def main() -> None:
    checks: list[tuple[str, str, object]] = []

    # Five-lane to grouped-three embedding. Grouped lanes identify
    # 21c=21s=x21 and 22c=22s=x22.
    E = sp.Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [0, 0, 1],
    ])
    G = E.T * E
    checks.append(check("group_metric_from_five_lane_embedding", G, sp.diag(1, 2, 2)))

    e_bar = sp.Matrix([1, 1, 1])
    e_a = sp.Matrix([4, -1, -1])
    e_b = sp.Matrix([0, 1, -1])
    basis = [("bar", e_bar), ("a", e_a), ("b", e_b)]

    norm = {}
    for ni, vi in basis:
        norm[ni] = sexpr((vi.T * G * vi)[0])
        for nj, vj in basis:
            if ni != nj:
                checks.append(check(f"G_orthogonality_{ni}_{nj}", (vi.T * G * vj)[0], 0))
    checks.append(check("norm_e_bar", norm["bar"], 5))
    checks.append(check("norm_e_a", norm["a"], 20))
    checks.append(check("norm_e_b", norm["b"], 4))

    def projector(v: sp.Matrix) -> sp.Matrix:
        den = (v.T * G * v)[0]
        return smat(v * (v.T * G) / den)

    P_bar = projector(e_bar)
    P_a = projector(e_a)
    P_b = projector(e_b)
    projectors = [("bar", P_bar), ("a", P_a), ("b", P_b)]

    P_bar_expected = sp.Matrix([[1, 2, 2], [1, 2, 2], [1, 2, 2]]) / 5
    P_a_expected = sp.Matrix([[16, -8, -8], [-4, 2, 2], [-4, 2, 2]]) / 20
    P_b_expected = sp.Matrix([[0, 0, 0], [0, 2, -2], [0, -2, 2]]) / 4
    checks.append(check("P_bar_formula", P_bar, P_bar_expected))
    checks.append(check("P_a_formula", P_a, P_a_expected))
    checks.append(check("P_b_formula", P_b, P_b_expected))
    checks.append(check("projector_sum_identity", P_bar + P_a + P_b, sp.eye(3)))
    for name, P in projectors:
        checks.append(check(f"P_{name}_idempotent", P * P, P))
        checks.append(check(f"P_{name}_G_self_adjoint", P.T * G, G * P))
    for i, (ni, Pi) in enumerate(projectors):
        for j, (nj, Pj) in enumerate(projectors):
            if i != j:
                checks.append(check(f"P_{ni}_P_{nj}_orthogonal", Pi * Pj, sp.zeros(3, 3)))

    # Coordinate extraction and inverse map.
    x20, x21, x22 = sp.symbols("x20 x21 x22")
    x = sp.Matrix([x20, x21, x22])
    xbar = sexpr((e_bar.T * G * x)[0] / norm["bar"])
    a_x = sexpr((e_a.T * G * x)[0] / norm["a"])
    b_x = sexpr((e_b.T * G * x)[0] / norm["b"])
    checks.append(check("xbar_formula", xbar, (x20 + 2*x21 + 2*x22)/5))
    checks.append(check("a_x_formula", a_x, (2*x20 - x21 - x22)/10))
    checks.append(check("b_x_formula", b_x, (x21 - x22)/2))
    checks.append(check("inverse_reconstruction", smat(xbar*e_bar + a_x*e_a + b_x*e_b), x))

    weighted_norm = sexpr((x.T * G * x)[0])
    decomposed_norm = sexpr(5*xbar**2 + 20*a_x**2 + 4*b_x**2)
    anis_vec = smat(P_a*x + P_b*x)
    anis_norm = sexpr((anis_vec.T * G * anis_vec)[0] / 5)
    anis_norm_expected = sexpr(4*a_x**2 + sp.Rational(4,5)*b_x**2)
    checks.append(check("weighted_norm_decomposition", weighted_norm, decomposed_norm))
    checks.append(check("normalized_anisotropy_norm", anis_norm, anis_norm_expected))

    # Direct isotropy implication: b=0 => x21=x22 and a=0 => x20=x22.
    iso_residuals = [sexpr((x21 - x22) - 2*b_x), sexpr((x20 - x22) - (5*a_x + b_x))]
    checks.append(status("isotropy_relations_from_a_b", all(r == 0 for r in iso_residuals), iso_residuals))

    # Weak axisymmetric splitting fingerprint.
    lam = sp.Matrix([sp.Integer(1), sp.Rational(1, 2), sp.Integer(-1)])
    lam_bar = sexpr((e_bar.T * G * lam)[0] / norm["bar"])
    lam_a = sexpr((e_a.T * G * lam)[0] / norm["a"])
    lam_b = sexpr((e_b.T * G * lam)[0] / norm["b"])
    checks.append(check("axisymmetric_lambda_trace_zero", lam_bar, 0))
    checks.append(check("axisymmetric_lambda_a", lam_a, sp.Rational(1, 4)))
    checks.append(check("axisymmetric_lambda_b", lam_b, sp.Rational(3, 4)))
    checks.append(check("axisymmetric_b_equals_3a", lam_b, 3*lam_a))

    # First-order transport formulas.
    eps, lams = sp.symbols("eps lambda")
    D0, D2, dD0, dD2, N0, dN0 = sp.symbols("D0 D2 dD0 dD2 N0 dN0")
    u_iso = -D2/D0
    D0A = D0 + eps*lams*dD0
    D2A = D2 + eps*lams*dD2
    uA_series = sp.series(-D2A/D0A, eps, 0, 2).removeO()
    u1_expected = sexpr(-(dD2 + u_iso*dD0)/D0)
    checks.append(check("first_order_u2_transport", sexpr(uA_series), sexpr(u_iso + eps*lams*u1_expected)))

    P_iso = N0/D0
    N0A = N0 + eps*lams*dN0
    PA_series = sp.series(N0A/D0A, eps, 0, 2).removeO()
    P1_expected = sexpr((dN0 - P_iso*dD0)/D0)
    checks.append(check("first_order_P0_transport", sexpr(PA_series), sexpr(P_iso + eps*lams*P1_expected)))

    u1, P1 = sp.symbols("u1 P1")
    u_vec = sp.Matrix([u_iso + eps*lam[i]*u1 for i in range(3)])
    P_vec = sp.Matrix([P_iso + eps*lam[i]*P1 for i in range(3)])
    u_a = sexpr((e_a.T * G * u_vec)[0] / norm["a"])
    u_b = sexpr((e_b.T * G * u_vec)[0] / norm["b"])
    p_a = sexpr((e_a.T * G * P_vec)[0] / norm["a"])
    p_b = sexpr((e_b.T * G * P_vec)[0] / norm["b"])
    checks.append(check("axisymmetric_u_b_equals_3_u_a", u_b, 3*u_a))
    checks.append(check("axisymmetric_P_b_equals_3_P_a", p_b, 3*p_a))

    n_fail = sum(1 for _, st, _ in checks if st != PASS)
    print("Stage V2-11 grouped real P2 projector calculus audit")
    print("=" * 72)
    print(f"checks_total: {len(checks)}")
    print(f"checks_passed: {len(checks)-n_fail}")
    print(f"checks_failed: {n_fail}")
    print()
    for name, st, detail in checks:
        print(f"{st:4s}  {name}")
        if st != PASS:
            print(f"      detail: {detail}")
    print()
    print("Key exact outputs")
    print("-" * 72)
    print(f"G = {G}")
    print(f"P_bar = {P_bar}")
    print(f"P_a = {P_a}")
    print(f"P_b = {P_b}")
    print(f"xbar = {xbar}")
    print(f"a_x = {a_x}")
    print(f"b_x = {b_x}")
    print("x = xbar*(1,1,1) + a_x*(4,-1,-1) + b_x*(0,1,-1)")
    print(f"x^T G x = {decomposed_norm}")
    print(f"A_x^2 = {anis_norm_expected}")
    print(f"axisymmetric lambda = {tuple(lam)}")
    print(f"lambda trace/a/b = ({lam_bar}, {lam_a}, {lam_b})")
    print(f"u2 first-order slope = {u1_expected}")
    print(f"P0 first-order slope = {P1_expected}")
    if n_fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
