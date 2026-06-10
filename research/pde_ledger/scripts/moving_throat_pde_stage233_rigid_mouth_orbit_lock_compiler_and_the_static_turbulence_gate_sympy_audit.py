import sympy as sp


def assert_zero(expr: sp.Expr) -> None:
    simplified = sp.simplify(expr)
    if simplified != 0:
        raise AssertionError(f"Expression is not zero: {simplified}")


def assert_close(actual: float, expected: float, tol: float = 1e-15) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{actual} !~= {expected} (tol={tol})")


def main() -> None:
    sp.init_printing()
    print("=== Stage 233 SymPy audit: rigid-mouth orbit-lock compiler and the static turbulence gate ===")

    # ------------------------------------------------------------------
    # 1. Exact Stage 188 branch-observable compiler
    # ------------------------------------------------------------------
    dln_Rtr, dln_Nstar, dln_eps_eta = sp.symbols(
        "dln_Rtr dln_Nstar dln_eps_eta", real=True
    )
    B_star, eps_eta_star = sp.symbols("B_star eps_eta_star", real=True)

    Theta1 = dln_Rtr
    Xi1 = dln_Nstar - B_star * dln_Rtr
    c_eta = eps_eta_star / (1 - eps_eta_star)
    R1 = -c_eta * dln_eps_eta - Xi1

    print("\nStage 188 observable compiler:")
    print("Theta1 =", Theta1)
    print("Xi1 =", Xi1)
    print("R1 =", R1)

    # ------------------------------------------------------------------
    # 2. Rigid-mouth / track-locked specialization
    # ------------------------------------------------------------------
    rigid_mouth = {dln_Rtr: 0}
    Theta1_rm = sp.simplify(Theta1.subs(rigid_mouth))
    Xi1_rm = sp.simplify(Xi1.subs(rigid_mouth))
    R1_rm = sp.simplify(R1.subs(rigid_mouth))

    assert_zero(Theta1_rm)
    assert_zero(Xi1_rm - dln_Nstar)
    assert_zero(R1_rm + c_eta * dln_eps_eta + dln_Nstar)

    R1_rm_expected = -c_eta * dln_eps_eta - dln_Nstar
    print("\nRigid-mouth specialization dln R_tr = 0:")
    print("Theta1 ->", Theta1_rm)
    print("Xi1 ->", Xi1_rm)
    print("R1 ->", R1_rm_expected)

    # ------------------------------------------------------------------
    # 3. Exact load-defect / prefactor identity
    # ------------------------------------------------------------------
    N01, N0, D01, D0 = sp.symbols("N01 N0 D01 D0", nonzero=True, real=True)

    P0 = N0 / D0
    P1 = (N01 * D0 - N0 * D01) / D0**2
    Xi_load = N01 / N0 - D01 / D0

    assert_zero(sp.simplify(P1 / P0 - Xi_load))

    print("\nOutgoing-prefactor slope identity:")
    print("P0 =", P0)
    print("P1 =", P1)
    print("Xi_load =", Xi_load)
    print("Verified: Xi_load = P1 / P0")

    # Stronger operator-rigidity specialization
    operator_rigid = {D01: 0}
    Xi_load_or = sp.simplify(Xi_load.subs(operator_rigid))
    assert_zero(Xi_load_or - N01 / N0)

    print("\nOperator-rigid specialization D01 = 0:")
    print("Xi_load ->", Xi_load_or)

    # ------------------------------------------------------------------
    # 4. Transported static ceiling and equivalent forms
    # ------------------------------------------------------------------
    epsilon, Pcrit, mhat0, Delta_norm, Tquad = sp.symbols(
        "epsilon Pcrit mhat0 Delta_norm Tquad", positive=True, real=True
    )

    gate_rhs = Pcrit * mhat0**2 / (Delta_norm + Tquad) - 1
    Pbar_expr = (Delta_norm + Tquad) / mhat0**2
    Pbar = sp.Symbol("Pbar", positive=True, real=True)
    gate_rhs_pbar = Pcrit / Pbar - 1

    assert_zero(sp.simplify(gate_rhs - (Pcrit / Pbar_expr - 1)))

    gate_rhs_cal = sp.simplify(gate_rhs.subs({Delta_norm: 0}))
    assert_zero(gate_rhs_cal - (Pcrit * mhat0**2 / Tquad - 1))

    print("\nTransported static ceiling:")
    print("|epsilon * Xi1| <=", gate_rhs)
    print("Equivalent Pbar form:", gate_rhs_pbar)
    print("Calibrated branch Delta_norm = 0 ->", gate_rhs_cal)

    # Track-locked gate
    gate_track_locked = sp.Abs(epsilon * Xi1_rm)
    print("\nTrack-locked internal-transfer gate:")
    print("|epsilon * dln N_*| <=", gate_rhs)

    # Operator-rigid internal-transfer reading
    gate_operator_rigid = sp.Abs(epsilon * Xi_load_or)
    print("If Xi1 is identified with the static load scalar on the operator-rigid branch:")
    print("|epsilon * (N01 / N0)| <=", gate_rhs)

    # ------------------------------------------------------------------
    # 5. Numerical recovery of the carried Stage 224 budgets
    # ------------------------------------------------------------------
    Pbar_num = sp.Float("0.002069792318062885")
    robust_budget = sp.Float("0.367930328492646")
    nonempty_budget = sp.Float("0.737619063660757")

    # Source file: scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py, ceilings dict.
    Pcrit_robust = sp.Float("0.0028313316855593175")
    Pcrit_nonempty = sp.Float("0.0035965105896846573")

    recovered_robust = sp.N(Pcrit_robust / Pbar_num - 1, 30)
    recovered_nonempty = sp.N(Pcrit_nonempty / Pbar_num - 1, 30)

    assert_close(float(recovered_robust), float(robust_budget), tol=1e-15)
    assert_close(float(recovered_nonempty), float(nonempty_budget), tol=1e-15)

    print("\nStage 223 / Stage 224 carried numbers:")
    print("Pbar compatibility point =", Pbar_num)
    print("Robust budget =", robust_budget)
    print("Nonempty budget =", nonempty_budget)
    print("Implied Pcrit (robust) =", Pcrit_robust)
    print("Implied Pcrit (nonempty) =", Pcrit_nonempty)
    print("Recovered robust budget =", recovered_robust)
    print("Recovered nonempty budget =", recovered_nonempty)

    print("\nAll Stage 233 symbolic and numerical checks passed.")


if __name__ == "__main__":
    main()
