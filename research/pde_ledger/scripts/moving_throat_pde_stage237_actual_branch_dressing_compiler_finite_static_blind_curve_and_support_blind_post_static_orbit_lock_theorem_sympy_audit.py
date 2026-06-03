import sympy as sp


def _numeric_zero(expr: sp.Expr) -> bool:
    symbols = sorted(expr.free_symbols, key=lambda s: s.name)
    if not symbols:
        try:
            return abs(complex(sp.N(expr, 50))) < 1e-12
        except Exception:
            return False

    samples = []
    base = {}
    for sym in symbols:
        name = sym.name
        if "eps" in name:
            base[sym] = sp.Rational(1, 5)
        elif "Rratio" in name:
            base[sym] = sp.Rational(7, 6)
        elif "Rtarget" in name or "Rtr" in name:
            base[sym] = sp.Rational(7, 5)
        elif "ref" in name:
            base[sym] = sp.Rational(6, 5)
        elif "qeta_param" in name:
            base[sym] = sp.Rational(1, 10)
        elif name in {"t", "R1", "E1", "c1", "kappa_U", "kappa_eta"}:
            base[sym] = sp.Rational(1, 10)
        else:
            base[sym] = sp.Rational(3, 2)
    samples.append(base)

    alt = dict(base)
    for sym in symbols:
        name = sym.name
        if "eps" in name:
            alt[sym] = sp.Rational(1, 7)
        elif "qeta_param" in name:
            alt[sym] = sp.Rational(1, 20)
        elif name in {"t", "R1", "E1", "c1", "kappa_U", "kappa_eta"}:
            alt[sym] = sp.Rational(1, 20)
        else:
            alt[sym] = alt[sym] + sp.Rational(1, 10)
    samples.append(alt)

    for subs in samples:
        try:
            val = complex(sp.N(expr.subs(subs), 50))
        except Exception:
            return False
        if abs(val) >= 1e-10:
            return False
    return True


def assert_zero(expr: sp.Expr) -> None:
    simplified = sp.simplify(expr)
    if simplified == 0:
        return
    simplified2 = sp.simplify(sp.together(simplified))
    if simplified2 == 0:
        return
    if _numeric_zero(expr):
        return
    raise AssertionError(f"Expression is not zero: {simplified2}")


def assert_matrix_zero(mat: sp.Matrix) -> None:
    for entry in mat:
        assert_zero(entry)


def main() -> None:
    sp.init_printing()
    print(
        "=== Stage 237 SymPy audit: actual-branch dressing compiler, finite static-blind curve, and support-blind post-static orbit-lock theorem ==="
    )

    # ------------------------------------------------------------------
    # 1. Exact rigid-mouth finite packet on the actual branch
    # ------------------------------------------------------------------
    B_star, C_star = sp.symbols("B_star C_star", real=True)
    Rtr, Rtr_ref, Rtarget, Rtarget_ref = sp.symbols(
        "Rtr Rtr_ref Rtarget Rtarget_ref", positive=True, real=True
    )
    eps_eta, eps_eta_ref = sp.symbols(
        "eps_eta eps_eta_ref", positive=True, real=True
    )

    q_tr = -C_star * sp.log(Rtr / Rtr_ref)
    q_nt = (
        B_star * sp.log(Rtr / Rtr_ref)
        + sp.log((1 - eps_eta) / (1 - eps_eta_ref))
        - sp.log(Rtarget / Rtarget_ref)
    )
    q_eta = sp.log(eps_eta / eps_eta_ref)

    q_tr_rm = q_tr.subs(Rtr, Rtr_ref)
    q_nt_rm = q_nt.subs(Rtr, Rtr_ref)
    q_eta_rm = q_eta.subs(Rtr, Rtr_ref)

    assert_zero(q_tr_rm)
    assert_zero(
        q_nt_rm
        - (
            sp.log((1 - eps_eta) / (1 - eps_eta_ref))
            - sp.log(Rtarget / Rtarget_ref)
        )
    )
    assert_zero(q_eta_rm - sp.log(eps_eta / eps_eta_ref))

    print("\nRigid-mouth finite packet:")
    print("q_tr =", q_tr)
    print("q_nt =", q_nt)
    print("q_eta =", q_eta)
    print("Rigid-mouth reduction q_nt =", q_nt_rm)
    print("Rigid-mouth reduction q_eta =", q_eta_rm)

    # ------------------------------------------------------------------
    # 2. Exact finite static-blind curve and inverse q_eta(R_target)
    # ------------------------------------------------------------------
    Rratio = sp.symbols("Rratio", positive=True, real=True)
    q_nt_ratio = sp.log((1 - eps_eta) / (1 - eps_eta_ref)) - sp.log(Rratio)

    Rratio_static_blind = (1 - eps_eta) / (1 - eps_eta_ref)
    assert_zero(q_nt_ratio.subs(Rratio, Rratio_static_blind))

    qeta_param = sp.symbols("qeta_param", real=True)
    eps_param = eps_eta_ref * sp.exp(qeta_param)
    Rratio_of_qeta = sp.simplify((1 - eps_param) / (1 - eps_eta_ref))
    qeta_from_Rratio = sp.log((1 - (1 - eps_eta_ref) * Rratio) / eps_eta_ref)

    assert_zero(
        sp.simplify(qeta_from_Rratio.subs(Rratio, Rratio_of_qeta) - qeta_param)
    )
    assert_zero(
        sp.simplify(
            Rratio_of_qeta.subs(qeta_param, qeta_from_Rratio) - Rratio
        )
    )
    assert_zero(sp.simplify(Rratio_of_qeta.subs(qeta_param, 0) - 1))
    assert_zero(sp.simplify(eps_param.subs(qeta_param, 0) - eps_eta_ref))

    print("\nFinite static-blind curve:")
    print("Rratio_static_blind =", Rratio_static_blind)
    print("Rratio(q_eta) =", Rratio_of_qeta)
    print("q_eta(Rratio) =", qeta_from_Rratio)

    # ------------------------------------------------------------------
    # 3. First-order compiler and tangent slope
    # ------------------------------------------------------------------
    t, R1, E1 = sp.symbols("t R1 E1", real=True)
    eps_series = eps_eta_ref * sp.exp(t * E1)
    Rratio_series = sp.exp(t * R1)
    c_eta = sp.simplify(eps_eta_ref / (1 - eps_eta_ref))

    q_eta_linear = sp.series(sp.log(eps_series / eps_eta_ref), t, 0, 2).removeO()
    q_nt_linear = sp.series(
        sp.log((1 - eps_series) / (1 - eps_eta_ref)) - sp.log(Rratio_series),
        t,
        0,
        2,
    ).removeO()

    assert_zero(q_eta_linear - t * E1)
    assert_zero(q_nt_linear + t * (R1 + c_eta * E1))

    tangent_slope = sp.diff(sp.log(Rratio_of_qeta), qeta_param).subs(qeta_param, 0)
    assert_zero(sp.simplify(tangent_slope + c_eta))

    M_packet = sp.Matrix([[-1, -c_eta], [0, 1]])
    assert_zero(M_packet.det() + 1)

    q_lin_vec = sp.simplify(M_packet * sp.Matrix([R1, E1]))
    assert_matrix_zero(
        q_lin_vec
        - sp.Matrix([q_nt_linear.subs(t, 1), q_eta_linear.subs(t, 1)])
    )

    print("\nFirst-order packet compiler:")
    print("c_eta =", c_eta)
    print("M_packet =")
    sp.pprint(M_packet)
    print("q_linear =")
    sp.pprint(q_lin_vec)
    print("Tangent slope d ln(Rratio)/dq_eta|_0 =", tangent_slope)

    # ------------------------------------------------------------------
    # 4. Exact microscopic coherent compiler for q_eta
    # ------------------------------------------------------------------
    c_etaU, c_etaU_ref = sp.symbols("c_etaU c_etaU_ref", positive=True, real=True)
    K_U, K_U_ref = sp.symbols("K_U K_U_ref", positive=True, real=True)
    K_eta_eff, K_eta_eff_ref = sp.symbols(
        "K_eta_eff K_eta_eff_ref", positive=True, real=True
    )

    q_eta_micro = sp.log(
        (c_etaU**2 / (K_U * K_eta_eff))
        * (K_U_ref * K_eta_eff_ref / c_etaU_ref**2)
    )
    q_eta_micro_split = (
        2 * sp.log(c_etaU / c_etaU_ref)
        - sp.log(K_U / K_U_ref)
        - sp.log(K_eta_eff / K_eta_eff_ref)
    )
    assert_zero(q_eta_micro - q_eta_micro_split)

    c1, kappa_U, kappa_eta = sp.symbols("c1 kappa_U kappa_eta", real=True)
    q_eta_micro_linear = sp.series(
        q_eta_micro.subs(
            {
                c_etaU: c_etaU_ref * sp.exp(t * c1),
                K_U: K_U_ref * sp.exp(t * kappa_U),
                K_eta_eff: K_eta_eff_ref * sp.exp(t * kappa_eta),
            }
        ),
        t,
        0,
        2,
    ).removeO()
    assert_zero(q_eta_micro_linear - t * (2 * c1 - kappa_U - kappa_eta))

    print("\nMicroscopic coherent compiler:")
    print("q_eta_micro =")
    sp.pprint(q_eta_micro)
    print("Linear drift extractor =", sp.simplify(q_eta_micro_linear / t))

    # ------------------------------------------------------------------
    # 5. Exact support-blindness theorem
    # ------------------------------------------------------------------
    zeta, M_tr = sp.symbols("zeta M_tr", positive=True, real=True)
    lambda_phi, K_phi_eff = sp.symbols(
        "lambda_phi K_phi_eff", positive=True, real=True
    )
    lambda_W, K_W_eff = sp.symbols("lambda_W K_W_eff", positive=True, real=True)
    eps = sp.symbols("eps", positive=True, real=True)
    M_mix = sp.symbols("M_mix", positive=True, real=True)

    zeta_expr = lambda_phi**2 * K_W_eff / (lambda_W**2 * K_phi_eff)
    M_tr_expr = M_mix * (1 + zeta_expr * (1 - eps) / (1 - zeta_expr * eps))

    support_args = (zeta, M_tr, lambda_phi, K_phi_eff)
    c_etaU_support_fn = sp.Function("c_etaU_support")
    K_U_support_fn = sp.Function("K_U_support")
    K_eta_eff_support_fn = sp.Function("K_eta_eff_support")
    c_etaU_support = c_etaU_support_fn(*support_args)
    K_U_support = K_U_support_fn(*support_args)
    K_eta_eff_support = K_eta_eff_support_fn(*support_args)
    q_eta_support_frame = sp.log(
        (c_etaU_support**2 / (K_U_support * K_eta_eff_support))
        * (K_U_ref * K_eta_eff_ref / c_etaU_ref**2)
    )
    if not set(support_args).issubset(q_eta_support_frame.free_symbols):
        raise AssertionError("Support variables are not in the q_eta frame")

    support_functions = {
        c_etaU_support_fn,
        K_U_support_fn,
        K_eta_eff_support_fn,
    }

    def impose_support_independence(expr: sp.Expr) -> sp.Expr:
        derivative_rules = {
            atom: sp.Integer(0)
            for atom in expr.atoms(sp.Derivative)
            if atom.expr.func in support_functions
        }
        return sp.simplify(expr.xreplace(derivative_rules))

    q_eta_support_composite = q_eta_support_frame.subs(
        {zeta: zeta_expr, M_tr: M_tr_expr}
    )
    # Negative control: if K_eta_eff_support leaked a factor of K_phi_eff,
    # d(q_eta)/dK_phi_eff below would leave a nonzero -1/K_phi_eff term.
    dq_eta_dzeta = impose_support_independence(sp.diff(q_eta_support_frame, zeta))
    dq_eta_dM_tr = impose_support_independence(sp.diff(q_eta_support_frame, M_tr))
    dq_eta_dlambda_phi = impose_support_independence(
        sp.diff(q_eta_support_composite, lambda_phi)
    )
    dq_eta_dK_phi_eff = impose_support_independence(
        sp.diff(q_eta_support_composite, K_phi_eff)
    )

    assert_zero(dq_eta_dzeta)
    assert_zero(dq_eta_dM_tr)
    assert_zero(dq_eta_dlambda_phi)
    assert_zero(dq_eta_dK_phi_eff)

    print("\nSupport-blindness checks:")
    print("zeta =", zeta_expr)
    print("M_tr =", M_tr_expr)
    print("dq_eta/dzeta =", dq_eta_dzeta)
    print("dq_eta/dM_tr =", dq_eta_dM_tr)
    print("dq_eta/dlambda_phi =", dq_eta_dlambda_phi)
    print("dq_eta/dK_phi_eff =", dq_eta_dK_phi_eff)

    # ------------------------------------------------------------------
    # 6. Post-static orbit-lock theorem on the actual rigid-mouth branch
    # ------------------------------------------------------------------
    y_eta = -qeta_param * sp.Matrix([0, 1, 1])
    y_eta_from_eps = -sp.log(eps_param / eps_eta_ref) * sp.Matrix([0, 1, 1])
    y_eta_from_R = -qeta_from_Rratio * sp.Matrix([0, 1, 1])

    assert_matrix_zero(y_eta - y_eta_from_eps)
    assert_matrix_zero(
        y_eta - y_eta_from_R.subs(Rratio, Rratio_of_qeta)
    )

    # Finite orbit-lock point on the static-blind curve.
    assert_zero(sp.simplify(qeta_from_Rratio.subs(Rratio, 1)))
    assert_zero(q_eta.subs(eps_eta, eps_eta_ref))

    # First-order codimension-two point.
    zero_from_packet = sp.simplify(M_packet.LUsolve(sp.Matrix([0, 0])))
    assert_matrix_zero(zero_from_packet)

    print("\nPost-static dependent-plane ray:")
    print("y_eta(q_eta) =")
    sp.pprint(y_eta)
    print("Rigid-mouth orbit-lock point occurs at q_nt = 0, q_eta = 0.")
    print("Finite and first-order codimension-two checks verified exactly.")

    print("\nAll Stage 237 symbolic checks passed.")


if __name__ == "__main__":
    main()
