import sympy as sp


def assert_zero(expr: sp.Expr) -> None:
    simplified = sp.simplify(expr)
    if simplified != 0:
        raise AssertionError(f"Expression is not zero: {simplified}")


def expand_log_simplify(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.expand_log(expr, force=True))


def main() -> None:
    sp.init_printing()
    print(
        "=== Stage 217 SymPy audit: direct branch-observable static gate and the two-observable kill test ==="
    )

    # ------------------------------------------------------------------
    # 1. Exact finite quotient chart and inverse
    # ------------------------------------------------------------------
    C_star, B_star = sp.symbols("C_star B_star", positive=True, real=True)
    Rtr_ref, Rtarget_ref = sp.symbols("Rtr_ref Rtarget_ref", positive=True, real=True)
    eps_eta_ref = sp.symbols("eps_eta_ref", positive=True, real=True)
    q_tr, q_nt, q_eta = sp.symbols("q_tr q_nt q_eta", real=True)

    Rtr = sp.symbols("Rtr", positive=True, real=True)
    Rtarget = sp.symbols("Rtarget", positive=True, real=True)
    eps_eta = sp.symbols("eps_eta", positive=True, real=True)

    qtr_fwd = -C_star * sp.log(Rtr / Rtr_ref)
    qnt_fwd = (
        B_star * sp.log(Rtr / Rtr_ref)
        + sp.log((1 - eps_eta) / (1 - eps_eta_ref))
        - sp.log(Rtarget / Rtarget_ref)
    )
    qeta_fwd = sp.log(eps_eta / eps_eta_ref)

    Rtr_inv = Rtr_ref * sp.exp(-q_tr / C_star)
    eps_eta_inv = eps_eta_ref * sp.exp(q_eta)
    Rtarget_inv = (
        Rtarget_ref
        * sp.exp(-q_nt - (B_star / C_star) * q_tr)
        * (1 - eps_eta_inv)
        / (1 - eps_eta_ref)
    )

    qtr_roundtrip = expand_log_simplify(qtr_fwd.subs({Rtr: Rtr_inv}))
    qnt_roundtrip = expand_log_simplify(
        qnt_fwd.subs({Rtr: Rtr_inv, eps_eta: eps_eta_inv, Rtarget: Rtarget_inv})
    )
    qeta_roundtrip = expand_log_simplify(qeta_fwd.subs({eps_eta: eps_eta_inv}))

    assert_zero(qtr_roundtrip - q_tr)
    assert_zero(qnt_roundtrip - q_nt)
    assert_zero(qeta_roundtrip - q_eta)

    print("\nExact finite quotient chart:")
    print("q_tr =", qtr_fwd)
    print("q_nt =", qnt_fwd)
    print("q_eta =", qeta_fwd)
    print("Inverse map verified exactly.")

    # ------------------------------------------------------------------
    # 2. First-order linearization of the finite chart
    # ------------------------------------------------------------------
    eps = sp.symbols("eps", real=True)
    r1, R1, E1 = sp.symbols("r1 R1 E1", real=True)
    eps_eta_star = sp.symbols("eps_eta_star", positive=True, real=True)

    Rtr_eps = Rtr_ref * sp.exp(eps * r1)
    Rtarget_eps = Rtarget_ref * sp.exp(eps * R1)
    eps_eta_eps = eps_eta_star * sp.exp(eps * E1)

    qtr_eps = qtr_fwd.subs({Rtr: Rtr_eps})
    qnt_eps = qnt_fwd.subs(
        {Rtr: Rtr_eps, Rtarget: Rtarget_eps, eps_eta: eps_eta_eps, eps_eta_ref: eps_eta_star}
    )
    qeta_eps = qeta_fwd.subs({eps_eta: eps_eta_eps, eps_eta_ref: eps_eta_star})

    qtr1 = sp.simplify(sp.diff(qtr_eps, eps).subs(eps, 0))
    qnt1 = sp.simplify(sp.diff(qnt_eps, eps).subs(eps, 0))
    qeta1 = sp.simplify(sp.diff(qeta_eps, eps).subs(eps, 0))

    c_eta = eps_eta_star / (1 - eps_eta_star)

    assert_zero(qtr1 + C_star * r1)
    assert_zero(qnt1 - (B_star * r1 - c_eta * E1 - R1))
    assert_zero(qeta1 - E1)

    qnt1_pretty = B_star * r1 - c_eta * E1 - R1
    print("\nFirst-order direct observable drifts:")
    print("q_tr^(1) =", qtr1)
    print("q_nt^(1) =", qnt1_pretty)
    print("q_eta^(1) =", qeta1)

    # ------------------------------------------------------------------
    # 3. Exact direct-branch defect compiler
    # ------------------------------------------------------------------
    Theta1 = -qtr1 / C_star
    Xi1 = sp.simplify(qnt1 + (B_star / C_star) * qtr1)
    Rcal1 = sp.simplify(-Xi1 - c_eta * qeta1)

    assert_zero(Theta1 - r1)
    assert_zero(Xi1 + R1 + c_eta * E1)
    assert_zero(Rcal1 - R1)

    E1_inv = sp.simplify(-(1 - eps_eta_star) / eps_eta_star * (Rcal1 + Xi1))
    assert_zero(E1_inv - E1)

    Xi1_pretty = -R1 - c_eta * E1
    print("\nTriangular direct-branch map:")
    print("Theta1 =", Theta1)
    print("Xi1 =", Xi1_pretty)
    print("R1 =", Rcal1)
    print("Inverse drift map verified.")

    # ------------------------------------------------------------------
    # 4. Exact cancellations and rigid-mouth specialization
    # ------------------------------------------------------------------
    assert_zero(sp.diff(Xi1, r1))

    qnt_rigid_finite = qnt_fwd.subs({Rtr: Rtr_ref})
    desired_qnt_rigid = sp.log((1 - eps_eta) / (1 - eps_eta_ref)) - sp.log(Rtarget / Rtarget_ref)
    if qnt_rigid_finite != desired_qnt_rigid:
        raise AssertionError(f"Unexpected rigid-mouth finite q_nt form: {qnt_rigid_finite}")

    rigid = {r1: 0}
    Xi1_rigid = sp.simplify(Xi1.subs(rigid))
    qnt1_rigid = sp.simplify(qnt1.subs(rigid))
    Theta1_rigid = sp.simplify(Theta1.subs(rigid))

    assert_zero(Theta1_rigid)
    assert_zero(Xi1_rigid - qnt1_rigid)

    print("\nExact cancellations:")
    print("d Xi1 / d(ln R_tr) =", sp.diff(Xi1, r1))
    print("Rigid-mouth finite q_nt =", qnt_rigid_finite)
    print("Rigid-mouth first-order Xi1 = q_nt^(1) =", sp.simplify(Xi1_pretty.subs(rigid)))

    # ------------------------------------------------------------------
    # 5. Two-observable static strip
    # ------------------------------------------------------------------
    vareps = sp.symbols("vareps", real=True, nonzero=True)
    robust = sp.Float("0.367930328492646")
    nonempty = sp.Float("0.737619063660757")

    Xi1_twoobs = -R1 - c_eta * E1
    robust_lower = -c_eta * E1 - robust / sp.Abs(vareps)
    robust_upper = -c_eta * E1 + robust / sp.Abs(vareps)
    nonempty_lower = -c_eta * E1 - nonempty / sp.Abs(vareps)
    nonempty_upper = -c_eta * E1 + nonempty / sp.Abs(vareps)

    robust_low_val = sp.simplify(Xi1_twoobs.subs(R1, robust_lower))
    robust_up_val = sp.simplify(Xi1_twoobs.subs(R1, robust_upper))
    nonempty_low_val = sp.simplify(Xi1_twoobs.subs(R1, nonempty_lower))
    nonempty_up_val = sp.simplify(Xi1_twoobs.subs(R1, nonempty_upper))

    assert_zero(robust_low_val - robust / sp.Abs(vareps))
    assert_zero(robust_up_val + robust / sp.Abs(vareps))
    assert_zero(nonempty_low_val - nonempty / sp.Abs(vareps))
    assert_zero(nonempty_up_val + nonempty / sp.Abs(vareps))

    print("\nTwo-observable static strip:")
    print("Xi1 =", Xi1_pretty)
    print("Robust strip:   R1 in [", robust_lower, ",", robust_upper, "]")
    print("Nonempty strip: R1 in [", nonempty_lower, ",", nonempty_upper, "]")

    # ------------------------------------------------------------------
    # 6. Canonical direct-branch families
    # ------------------------------------------------------------------
    xi = sp.symbols("xi", real=True)

    # Pure target drift
    Xi_target = sp.simplify(Xi1_twoobs.subs({E1: 0, R1: -xi}))
    assert_zero(Xi_target - xi)

    # Pure dressing drift
    Xi_dressing = sp.simplify(Xi1_twoobs.subs({R1: 0, E1: -xi / c_eta}))
    assert_zero(Xi_dressing - xi)

    # Balanced minimal-norm family via Lagrange multiplier solve
    lam = sp.symbols("lam", real=True)
    eqs = [
        sp.Eq(sp.diff(R1**2 + E1**2 + lam * (-R1 - c_eta * E1 - xi), R1), 0),
        sp.Eq(sp.diff(R1**2 + E1**2 + lam * (-R1 - c_eta * E1 - xi), E1), 0),
        sp.Eq(-R1 - c_eta * E1 - xi, 0),
    ]
    sol = sp.solve(eqs, [R1, E1, lam], dict=True)
    if len(sol) != 1:
        raise AssertionError(f"Expected one balanced-family solution, got {sol}")
    sol = sol[0]

    R1_bal = sp.simplify(sol[R1])
    E1_bal = sp.simplify(sol[E1])

    assert_zero(R1_bal + xi / (1 + c_eta**2))
    assert_zero(E1_bal + c_eta * xi / (1 + c_eta**2))
    assert_zero(sp.simplify(Xi1_twoobs.subs({R1: R1_bal, E1: E1_bal}) - xi))

    R1_bal_expected = -xi / (1 + c_eta**2)
    E1_bal_expected = -c_eta * xi / (1 + c_eta**2)

    print("\nCanonical direct-branch families:")
    print("Pure target drift:    (R1, E1) =", (-xi, 0))
    print("Pure dressing drift:  (R1, E1) =", (0, -xi / c_eta))
    print("Balanced minimal norm: (R1, E1) =", (R1_bal_expected, E1_bal_expected))

    print("\nAll Stage 217 symbolic checks passed.")


if __name__ == "__main__":
    main()
