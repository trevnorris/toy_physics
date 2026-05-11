import math
import sympy as sp


def assert_close(actual: float, expected: float, tol: float = 1e-12) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{actual} !~= {expected} (tol={tol})")


def vector_to_float_list(v: sp.Matrix, digits: int = 15) -> list[float]:
    return [float(sp.N(v[i], digits)) for i in range(v.rows)]


def mat_entry_list(M: sp.Matrix, digits: int = 15) -> list[list[float]]:
    return [[float(sp.N(M[i, j], digits)) for j in range(M.cols)] for i in range(M.rows)]


def coeff_vector(expr: sp.Expr, vars_: list[sp.Symbol]) -> sp.Matrix:
    expr = sp.expand(expr)
    return sp.Matrix([sp.N(expr.coeff(v), 30) for v in vars_])


def main() -> None:
    sp.init_printing()

    print("=== Stage 226 SymPy audit: strict 5PN even-gate package and pure-transfer subcorridor ===")

    # ------------------------------------------------------------------
    # 1. Exact bridge Xi_load = Xi1 and strict-even formulas
    # ------------------------------------------------------------------
    eps, lam = sp.symbols("eps lam", real=True)
    D0, D2, D4, N0 = sp.symbols("D0 D2 D4 N0", nonzero=True)
    D01, D21, D41, N01 = sp.symbols("D01 D21 D41 N01", real=True)

    D0A = D0 + eps * lam * D01
    N0A = N0 + eps * lam * N01

    u2 = -D2 / D0
    u4 = (D2**2 - D0 * D4) / D0**2
    P0 = N0 / D0
    P0A = N0A / D0A
    P1 = sp.simplify(sp.diff(P0A, eps).subs(eps, 0) / lam)
    Xi1 = sp.simplify(P1 / P0)
    Xi_load = sp.simplify(N01 / N0 - D01 / D0)

    assert sp.simplify(Xi1 - Xi_load) == 0

    K1 = D21 + D01 / 9
    H_even = D41 - sp.Rational(2, 3) * D21 - D01 / 27

    K1_comp = sp.simplify(K1.subs(D21, -u2 * D01))
    H_even_comp = sp.simplify(H_even.subs({D21: -u2 * D01, D41: D4 * D01 / D0}))

    U2, U4 = sp.symbols("U2 U4", real=True)
    H_even_one_pole = sp.simplify(
        H_even_comp.subs({D2: -U2 * D0, D4: D0 * (U2**2 - U4)}).subs(U4, 4 * U2**2)
    )

    assert sp.simplify(K1_comp - (sp.Rational(1, 9) - u2) * D01) == 0
    assert sp.simplify(H_even_comp - (D4 / D0 + sp.Rational(2, 3) * u2 - sp.Rational(1, 27)) * D01) == 0
    assert sp.simplify(H_even_one_pole - (-3 * U2**2 + sp.Rational(2, 3) * U2 - sp.Rational(1, 27)) * D01) == 0

    print("\nVerified exact bridge and strict-even formulas:")
    print("Xi_load = Xi1 =", Xi1)
    print("K1 on compensation surface =", K1_comp)
    print("H_even on compensation surface =", H_even_comp)
    print("H_even on one-pole branch =", H_even_one_pole)

    # ------------------------------------------------------------------
    # 2. Primitive compiler on the explicit Stage 223 compatibility branch
    # ------------------------------------------------------------------
    kappa = sp.symbols("kappa", positive=True)
    K, M = sp.symbols("K M", positive=True)
    lamB, varpi = sp.symbols("lamB varpi", positive=True)
    lamU, lamW, lamR = sp.symbols("lamU lamW lamR", positive=True)
    OmU, OmW = sp.symbols("OmU OmW", positive=True)

    xK, xM = sp.symbols("xK xM", real=True)
    xLB, xV = sp.symbols("xLB xV", real=True)
    xLU, xLW, xLR = sp.symbols("xLU xLW xLR", real=True)
    xOU, xOW = sp.symbols("xOU xOW", real=True)

    all_vars = [xK, xM, xLB, xV, xLU, xLW, xLR, xOU, xOW]
    mixed_vars = [xLU, xLW, xLR, xOU, xOW]

    C = kappa * lamB
    GU = lamU
    GW = kappa * lamW
    R = kappa * lamR

    Delta = OmU**2 * OmW**2 - R**2
    S2 = OmU**2 + OmW**2
    H = GU**2 + GW**2
    Q = GU**2 * OmW**2 + 2 * GU * GW * R + GW**2 * OmU**2
    P = OmU**2 * GW + R * GU

    B0 = C**2 / varpi**2
    B2 = C**2 / varpi**4
    B4 = C**2 / varpi**6
    Z0 = Q / Delta
    Z2 = (Q * S2 - H * Delta) / Delta**2
    Z4 = (Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3
    N0_bundle = P**2 / Delta**2

    # Dressed primitive parameters
    Kd = K * sp.exp(eps * xK)
    Md = M * sp.exp(eps * xM)
    lamBd = lamB * sp.exp(eps * xLB)
    vard = varpi * sp.exp(eps * xV)
    lamUd = lamU * sp.exp(eps * xLU)
    lamWd = lamW * sp.exp(eps * xLW)
    lamRd = lamR * sp.exp(eps * xLR)
    OmUd = OmU * sp.exp(eps * xOU)
    OmWd = OmW * sp.exp(eps * xOW)

    Cd = kappa * lamBd
    GUd = lamUd
    GWd = kappa * lamWd
    Rd = kappa * lamRd

    Deltad = OmUd**2 * OmWd**2 - Rd**2
    S2d = OmUd**2 + OmWd**2
    Hd = GUd**2 + GWd**2
    Qd = GUd**2 * OmWd**2 + 2 * GUd * GWd * Rd + GWd**2 * OmUd**2
    Pd = OmUd**2 * GWd + Rd * GUd

    B0d = Cd**2 / vard**2
    B2d = Cd**2 / vard**4
    B4d = Cd**2 / vard**6
    Z0d = Qd / Deltad
    Z2d = (Qd * S2d - Hd * Deltad) / Deltad**2
    Z4d = (Qd * (S2d**2 - Deltad) - S2d * Hd * Deltad) / Deltad**3
    N0d = Pd**2 / Deltad**2

    D01_compiler = sp.simplify(sp.diff(Kd - B0d - Z0d, eps).subs(eps, 0))
    D21_compiler = sp.simplify(sp.diff(-(Md + B2d + Z2d), eps).subs(eps, 0))
    D41_compiler = sp.simplify(sp.diff(-(B4d + Z4d), eps).subs(eps, 0))
    N01_compiler = sp.simplify(sp.diff(N0d, eps).subs(eps, 0))

    sample = {
        kappa: 2 * sp.sqrt(2) / sp.pi,
        lamB: sp.Rational(1, 2),
        lamU: sp.Rational(3, 10),
        lamW: sp.Rational(2, 5),
        lamR: sp.Rational(1, 4),
        OmU: sp.Integer(1),
        OmW: sp.Rational(7, 5),
        varpi: sp.Integer(2),
        M: sp.Integer(1),
    }

    B0_s = sp.simplify(B0.subs(sample))
    B2_s = sp.simplify(B2.subs(sample))
    B4_s = sp.simplify(B4.subs(sample))
    Z0_s = sp.simplify(Z0.subs(sample))
    Z2_s = sp.simplify(Z2.subs(sample))
    Z4_s = sp.simplify(Z4.subs(sample))
    N0_s = sp.simplify(N0_bundle.subs(sample))

    D0_compat = sp.simplify(3 * (sp.Integer(1) + B2_s + Z2_s) ** 2 / (B4_s + Z4_s))
    K_compat = sp.simplify(B0_s + Z0_s + D0_compat)
    sample_full = dict(sample)
    sample_full[K] = K_compat

    D0_s = sp.simplify((K - B0 - Z0).subs(sample_full))
    D2_s = sp.simplify((-(M + B2 + Z2)).subs(sample_full))
    D4_s = sp.simplify((-(B4 + Z4)).subs(sample_full))
    u2_s = sp.simplify(-D2_s / D0_s)
    u4_s = sp.simplify((D2_s**2 - D0_s * D4_s) / D0_s**2)

    assert_close(float(sp.N(K_compat, 16)), 24.4737548792910)
    assert_close(float(sp.N(D0_s, 16)), 24.2373099886223)
    assert_close(float(sp.N(D2_s, 16)), -1.18562046858190)
    assert_close(float(sp.N(D4_s, 16)), -0.173991572849491)
    assert_close(float(sp.N(u2_s, 16)), 0.0489171640391802)
    assert_close(float(sp.N(u4_s, 16)), 0.00957155575054425)
    assert_close(float(sp.N(D4_s / D0_s, 16)), -0.00717866681290820)
    assert_close(float(sp.N(N0_s / D0_s, 18)), 0.00206979231806289)
    assert sp.simplify(u4_s - 4 * u2_s**2) == 0

    print("\nVerified Stage 223 compatibility branch:")
    print("K_compat =", sp.N(K_compat, 16))
    print("D0       =", sp.N(D0_s, 16))
    print("D2       =", sp.N(D2_s, 16))
    print("D4       =", sp.N(D4_s, 16))
    print("u2       =", sp.N(u2_s, 16))
    print("u4       =", sp.N(u4_s, 16))
    print("P0       =", sp.N(N0_s / D0_s, 18))

    # Mixed-sector-only restriction of the primitive compiler
    zero_nonmixed = {xK: 0, xM: 0, xLB: 0, xV: 0}
    D01_mixed = sp.expand(sp.simplify(D01_compiler.subs(sample_full))).subs(zero_nonmixed)
    D21_mixed = sp.expand(sp.simplify(D21_compiler.subs(sample_full))).subs(zero_nonmixed)
    D41_mixed = sp.expand(sp.simplify(D41_compiler.subs(sample_full))).subs(zero_nonmixed)
    N01_mixed = sp.expand(sp.simplify(N01_compiler.subs(sample_full))).subs(zero_nonmixed)
    Xi1_mixed = sp.expand(sp.simplify(N01_mixed / N0_s - D01_mixed / D0_s))

    K1_mixed = sp.expand(sp.simplify(D21_mixed + D01_mixed / 9))
    H_even_mixed = sp.expand(sp.simplify(D41_mixed - sp.Rational(2, 3) * D21_mixed - D01_mixed / 27))

    # ------------------------------------------------------------------
    # 3. Stage 225 compensation versus strict 5PN even gates
    # ------------------------------------------------------------------
    K1_coeff_D01 = sp.simplify((sp.Rational(1, 9) - u2_s))
    H_even_coeff_D01 = sp.simplify((D4_s / D0_s + sp.Rational(2, 3) * u2_s - sp.Rational(1, 27)))

    assert_close(float(sp.N(K1_coeff_D01, 18)), 0.0621939470719309)
    assert_close(float(sp.N(H_even_coeff_D01, 18)), -0.0116042611571584)

    print("\nVerified compensation-vs-strict-even comparison on the sample branch:")
    print("K1 =", sp.N(K1_coeff_D01, 16), "* D01")
    print("H_even =", sp.N(H_even_coeff_D01, 16), "* D01")

    # ------------------------------------------------------------------
    # 4. Mixed-sector-only strict even-gate corridor
    # ------------------------------------------------------------------
    M_even = sp.Matrix([
        [sp.N(K1_mixed.coeff(v), 30) for v in mixed_vars],
        [sp.N(H_even_mixed.coeff(v), 30) for v in mixed_vars],
    ])

    expected_even = [
        [-0.255028994532, -0.132167046465, -0.067875763349, 0.568483003085, 0.300375205864],
        [-0.086801409924, -0.010480267714, -0.045759251298, 0.510038362482, 0.131455867026],
    ]
    for i in range(2):
        for j in range(5):
            assert_close(float(sp.N(M_even[i, j], 15)), expected_even[i][j], tol=1e-12)

    assert M_even.rank() == 2
    even_null = M_even.nullspace()
    assert len(even_null) == 3

    expected_w = [
        [-0.606454972136, 0.656652628212, 1.0, 0.0, 0.0],
        [6.983614208603, -9.174307357027, 0.0, 1.0, 0.0],
        [1.616693986742, -0.846872492318, 0.0, 0.0, 1.0],
    ]
    for b, exp_b in zip(even_null, expected_w):
        vals = vector_to_float_list(b)
        for actual, expected in zip(vals, exp_b):
            assert_close(actual, expected, tol=1e-11)

    Xi_coeff = coeff_vector(Xi1_mixed, mixed_vars)
    Xi_w = [float(sp.N((Xi_coeff.T * b)[0], 16)) for b in even_null]
    expected_Xi_w = [1.33691841376792, -13.9944400566810, -5.02163500066813]
    for actual, expected in zip(Xi_w, expected_Xi_w):
        assert_close(actual, expected, tol=1e-11)

    NS_even = sp.Matrix.hstack(*even_null)
    Q_even, _ = NS_even.QRdecomposition()
    proj_even = sp.simplify(Q_even * Q_even.T)
    sigma_even = sp.sqrt((Xi_coeff.T * proj_even * Xi_coeff)[0])
    sigma_even_f = float(sp.N(sigma_even, 20))
    assert_close(sigma_even_f, 2.67386816837173, tol=1e-11)

    print("\nVerified strict mixed-sector even-gate corridor:")
    print("Matrix =")
    sp.pprint(sp.N(M_even, 12))
    print("rank   =", M_even.rank(), " nullity =", len(even_null))
    print("w1     =", vector_to_float_list(even_null[0]))
    print("w2     =", vector_to_float_list(even_null[1]))
    print("w3     =", vector_to_float_list(even_null[2]))
    print("Xi1(w) =", Xi_w)
    print("sigma_even =", sp.N(sigma_even, 18))

    # ------------------------------------------------------------------
    # 5. Pure-transfer subcorridor
    # ------------------------------------------------------------------
    M_transfer = sp.Matrix([
        [sp.N(D01_mixed.coeff(v), 30) for v in mixed_vars],
        [sp.N(D21_mixed.coeff(v), 30) for v in mixed_vars],
        [sp.N(D41_mixed.coeff(v), 30) for v in mixed_vars],
    ])

    assert M_transfer.rank() == 3
    transfer_null = M_transfer.nullspace()
    assert len(transfer_null) == 2

    expected_t = [
        [-4.359222794718, 3.107402039105, 18.703510605854, 1.0, 0.0],
        [1.909256655687, -1.163651238154, -0.482414494705, 0.0, 1.0],
    ]
    for b, exp_b in zip(transfer_null, expected_t):
        vals = vector_to_float_list(b)
        for actual, expected in zip(vals, exp_b):
            assert_close(actual, expected, tol=1e-11)

    N01_coeff = coeff_vector(N01_mixed, mixed_vars)
    Xi_t = [float(sp.N((Xi_coeff.T * b)[0], 16)) for b in transfer_null]
    N01_t = [float(sp.N((N01_coeff.T * b)[0], 16)) for b in transfer_null]
    expected_Xi_t = [11.0106276743889, -5.66658382170817]
    expected_N01_t = [0.552361328292489, -0.284270966124842]
    for actual, expected in zip(Xi_t, expected_Xi_t):
        assert_close(actual, expected, tol=1e-11)
    for actual, expected in zip(N01_t, expected_N01_t):
        assert_close(actual, expected, tol=1e-11)

    # Verify that D01 = D21 = D41 = 0 on the pure-transfer basis
    for b in transfer_null:
        assert_close(float(sp.N((coeff_vector(D01_mixed, mixed_vars).T * b)[0], 18)), 0.0, tol=1e-12)
        assert_close(float(sp.N((coeff_vector(D21_mixed, mixed_vars).T * b)[0], 18)), 0.0, tol=1e-12)
        assert_close(float(sp.N((coeff_vector(D41_mixed, mixed_vars).T * b)[0], 18)), 0.0, tol=1e-12)

    NS_transfer = sp.Matrix.hstack(*transfer_null)
    Q_transfer, _ = NS_transfer.QRdecomposition()
    proj_transfer = sp.simplify(Q_transfer * Q_transfer.T)
    sigma_transfer = sp.sqrt((Xi_coeff.T * proj_transfer * Xi_coeff)[0])
    sigma_transfer_f = float(sp.N(sigma_transfer, 20))
    assert_close(sigma_transfer_f, 2.31561904386057, tol=1e-11)

    print("\nVerified pure-transfer subcorridor:")
    print("rank   =", M_transfer.rank(), " nullity =", len(transfer_null))
    print("t1     =", vector_to_float_list(transfer_null[0]))
    print("t2     =", vector_to_float_list(transfer_null[1]))
    print("Xi1(t) =", Xi_t)
    print("N01(t) =", N01_t)
    print("sigma_transfer =", sp.N(sigma_transfer, 18))

    # ------------------------------------------------------------------
    # 6. Transported Stage 224 budgets on the strict corridors
    # ------------------------------------------------------------------
    budget_10_both = sp.Float("0.367930328492646")
    budget_10_nonempty = sp.Float("0.737619063660757")
    budget_30_both = sp.Float("2.94889585703134")
    budget_30_nonempty = sp.Float("4.63505472371892")

    even_budgets = [
        float(sp.N(budget_10_both / sigma_even, 15)),
        float(sp.N(budget_10_nonempty / sigma_even, 15)),
        float(sp.N(budget_30_both / sigma_even, 15)),
        float(sp.N(budget_30_nonempty / sigma_even, 15)),
    ]
    transfer_budgets = [
        float(sp.N(budget_10_both / sigma_transfer, 15)),
        float(sp.N(budget_10_nonempty / sigma_transfer, 15)),
        float(sp.N(budget_30_both / sigma_transfer, 15)),
        float(sp.N(budget_30_nonempty / sigma_transfer, 15)),
    ]

    expected_even_budgets = [
        0.137602269567650,
        0.275862165676603,
        1.10285760977778,
        1.73346419189450,
    ]
    expected_transfer_budgets = [
        0.158890698998242,
        0.318540765855427,
        1.27348056877049,
        2.00164821411704,
    ]
    for actual, expected in zip(even_budgets, expected_even_budgets):
        assert_close(actual, expected, tol=1e-11)
    for actual, expected in zip(transfer_budgets, expected_transfer_budgets):
        assert_close(actual, expected, tol=1e-11)

    print("\nVerified transported ceiling budgets:")
    print("strict even-gate corridor  =", even_budgets)
    print("pure-transfer subcorridor =", transfer_budgets)

    print("\nAll Stage 226 checks passed.")


if __name__ == "__main__":
    main()
