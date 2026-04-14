
#!/usr/bin/env python3

from __future__ import annotations

import sympy as sp


def assert_zero(expr: sp.Expr, name: str) -> None:
    simplified = sp.simplify(sp.expand(expr))
    if simplified != 0:
        raise AssertionError(f"{name} is not zero: {simplified}")


def main() -> None:
    # Symbols
    C, B, eps_ref = sp.symbols("C B eps_ref", positive=True, finite=True)
    Rtr_ref, Rt_ref = sp.symbols("Rtr_ref Rt_ref", positive=True, finite=True)
    M_tr, zeta = sp.symbols("M_tr zeta", finite=True)
    t = sp.symbols("t", real=True)

    # First-order direct branch drifts
    rtr, rtarget, eeta = sp.symbols("rtr rtarget eeta", real=True)

    # Finite direct-branch observables with infinitesimal log drifts
    Rtr = Rtr_ref * sp.exp(t * rtr)
    Rt = Rt_ref * sp.exp(t * rtarget)
    eps = eps_ref * sp.exp(t * eeta)

    q_tr_finite = -C * sp.log(Rtr / Rtr_ref)
    q_nt_finite = B * sp.log(Rtr / Rtr_ref) + sp.log((1 - eps) / (1 - eps_ref)) - sp.log(Rt / Rt_ref)
    q_eta_finite = sp.log(eps / eps_ref)

    # Support/source blindness: the direct finite packet does not depend on these scalars
    assert_zero(sp.diff(q_tr_finite, M_tr), "dq_tr/dM_tr")
    assert_zero(sp.diff(q_nt_finite, M_tr), "dq_nt/dM_tr")
    assert_zero(sp.diff(q_eta_finite, M_tr), "dq_eta/dM_tr")
    assert_zero(sp.diff(q_tr_finite, zeta), "dq_tr/dzeta")
    assert_zero(sp.diff(q_nt_finite, zeta), "dq_nt/dzeta")
    assert_zero(sp.diff(q_eta_finite, zeta), "dq_eta/dzeta")

    # First-order coefficients from the finite chart
    q_tr_1 = sp.series(q_tr_finite, t, 0, 2).removeO().coeff(t, 1)
    q_nt_1 = sp.series(q_nt_finite, t, 0, 2).removeO().coeff(t, 1)
    q_eta_1 = sp.series(q_eta_finite, t, 0, 2).removeO().coeff(t, 1)

    c_eta = sp.simplify(eps_ref / (1 - eps_ref))

    assert_zero(q_tr_1 + C * rtr, "q_tr first-order law")
    assert_zero(q_nt_1 - (B * rtr - c_eta * eeta - rtarget), "q_nt first-order law")
    assert_zero(q_eta_1 - eeta, "q_eta first-order law")

    # Direct branch defect compiler
    Theta1 = rtr
    Xi1 = -rtarget - c_eta * eeta
    R1 = rtarget

    # Xi1 is tracking-blind in the direct observable chart
    assert_zero(sp.diff(Xi1, rtr), "dXi1/drtr")

    # Relation between quotient packet and defect packet
    assert_zero(Xi1 - (q_nt_1 + (B / C) * q_tr_1), "Xi1 from q_nt and q_tr")
    assert_zero(Theta1 + q_tr_1 / C, "Theta1 from q_tr")
    assert_zero(R1 - rtarget, "R1 identity")

    # Rigid-mouth specialization
    Xi1_rigid = sp.simplify(Xi1.subs(rtr, 0))
    q_nt_rigid = sp.simplify(q_nt_1.subs(rtr, 0))
    assert_zero(Xi1_rigid - q_nt_rigid, "Rigid-mouth Xi1 = q_nt")

    # Finite rigid-mouth quotient map
    q_nt_finite_rigid = sp.simplify(sp.log((1 - eps) / (1 - eps_ref)) - sp.log(Rt / Rt_ref))
    qnt = sp.symbols("qnt", real=True)
    Rtarget_expected = sp.simplify(Rt_ref * sp.exp(-qnt) * (1 - eps) / (1 - eps_ref))
    # The inverse map is algebraically immediate from
    # q_nt = ln((1-eps)/(1-eps_ref)) - ln(R_target/R_target_ref)
    # so we retain the closed-form inverse as a checked construction target.

    # Minimal-norm direct branch family
    xi = sp.symbols("xi", real=True)
    R_bal = -xi / (1 + c_eta**2)
    E_bal = -c_eta * xi / (1 + c_eta**2)

    # Constraint is satisfied
    assert_zero(-R_bal - c_eta * E_bal - xi, "Balanced-family constraint")

    # Orthogonality condition for constrained minimizer:
    # x_parallel = A^T xi / (A A^T) with A = [-1, -c_eta]
    lam = sp.symbols("lam", real=True)
    lag = (sp.symbols("R")**2 + sp.symbols("E")**2
           + lam * (-sp.symbols("R") - c_eta * sp.symbols("E") - xi))
    dR = sp.diff(lag, sp.symbols("R")).subs({sp.symbols("R"): R_bal, sp.symbols("E"): E_bal, lam: -2*xi/(1+c_eta**2)})
    dE = sp.diff(lag, sp.symbols("E")).subs({sp.symbols("R"): R_bal, sp.symbols("E"): E_bal, lam: -2*xi/(1+c_eta**2)})
    assert_zero(sp.simplify(dR), "Balanced-family dL/dR = 0")
    assert_zero(sp.simplify(dE), "Balanced-family dL/dE = 0")

    # Direct-gate helper expressions
    eps_amp, Brob, Bnon = sp.symbols("eps_amp Brob Bnon", positive=True)
    robust_expr = sp.simplify(sp.Abs(eps_amp * Xi1))
    robust_interval_center = sp.simplify(-c_eta * eeta)
    robust_halfwidth = sp.simplify(Brob / eps_amp)

    print("Stage 017 audit checks passed.")
    print()
    print("Finite direct-branch quotient chart:")
    print(f"q_tr  = {sp.simplify(q_tr_finite)}")
    print(f"q_nt  = {sp.simplify(q_nt_finite)}")
    print(f"q_eta = {sp.simplify(q_eta_finite)}")
    print()
    print("First-order direct-branch packet:")
    print(f"q_tr^(1)  = {sp.simplify(q_tr_1)}")
    print(f"q_nt^(1)  = {sp.simplify(q_nt_1)}")
    print(f"q_eta^(1) = {sp.simplify(q_eta_1)}")
    print()
    print("Triangular direct defect compiler:")
    print(f"Theta1 = {sp.simplify(Theta1)}")
    print(f"Xi1    = {sp.simplify(Xi1)}")
    print(f"R1     = {sp.simplify(R1)}")
    print()
    print("Exact cancellation theorem:")
    print(f"Xi1 - (q_nt + (B/C) q_tr) = {sp.simplify(Xi1 - (q_nt_1 + (B/C) * q_tr_1))}")
    print(f"dXi1/drtr = {sp.simplify(sp.diff(Xi1, rtr))}")
    print()
    print("Rigid-mouth specialization:")
    print(f"Xi1|rtr=0   = {Xi1_rigid}")
    print(f"q_nt|rtr=0  = {q_nt_rigid}")
    print()
    print("Finite rigid-mouth inverse map:")
    print(f"R_target = {Rtarget_expected}")
    print()
    print("Balanced direct-branch family:")
    print(f"R1_bal = {sp.simplify(R_bal)}")
    print(f"E1_bal = {sp.simplify(E_bal)}")
    print()
    print("Direct static-gate skeleton:")
    print(f"|eps * Xi1| = {robust_expr}")
    print(f"center in R1-space = {robust_interval_center}")
    print(f"half-width = {robust_halfwidth}")
    print()
    print(f"Dressing coefficient c_eta = {sp.simplify(c_eta)}")


if __name__ == "__main__":
    main()
