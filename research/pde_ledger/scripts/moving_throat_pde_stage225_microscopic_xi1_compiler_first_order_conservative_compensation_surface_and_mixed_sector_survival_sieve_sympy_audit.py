
import math
import sympy as sp


def assert_close(actual: float, expected: float, tol: float = 1e-12) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{actual} !~= {expected} (tol={tol})")


def main() -> None:
    sp.init_printing()

    print("=== Stage 225 SymPy audit: microscopic Xi1 compiler and conservative compensation surface ===")

    # ------------------------------------------------------------------
    # Exact arbitrary-base first-order formulas
    # ------------------------------------------------------------------
    eps, lam = sp.symbols("eps lam", real=True)
    D0, D2, D4, N0 = sp.symbols("D0 D2 D4 N0", nonzero=True)
    D01, D21, D41, N01 = sp.symbols("D01 D21 D41 N01", real=True)

    D0A = D0 + eps * lam * D01
    D2A = D2 + eps * lam * D21
    D4A = D4 + eps * lam * D41
    N0A = N0 + eps * lam * N01

    u2A = -D2A / D0A
    u4A = (D2A**2 - D0A * D4A) / D0A**2
    P0A = N0A / D0A

    u2 = -D2 / D0
    u4 = (D2**2 - D0 * D4) / D0**2
    P0 = N0 / D0

    u2_1 = sp.simplify(sp.diff(u2A, eps).subs(eps, 0) / lam)
    u4_1 = sp.simplify(sp.diff(u4A, eps).subs(eps, 0) / lam)
    P1 = sp.simplify(sp.diff(P0A, eps).subs(eps, 0) / lam)
    Xi1 = sp.simplify(P1 / P0)

    assert sp.simplify(u2_1 - (-D0 * D21 + D2 * D01) / D0**2) == 0
    assert sp.simplify(u2_1 + (D21 + u2 * D01) / D0) == 0
    assert sp.simplify(
        u4_1
        - (
            -D0 * (D0 * D41 + D01 * D4 - 2 * D2 * D21)
            + 2 * D01 * (D0 * D4 - D2**2)
        )
        / D0**3
    ) == 0
    assert sp.simplify(Xi1 - (N01 / N0 - D01 / D0)) == 0

    print("\nVerified arbitrary-base first-order formulas:")
    print("u2^(1) =", u2_1)
    print("u4^(1) =", u4_1)
    print("Xi1    =", Xi1)

    # ------------------------------------------------------------------
    # Exact conservative compensation surface
    # ------------------------------------------------------------------
    D21_comp = sp.solve(sp.Eq(u2_1, 0), D21)[0]
    D41_comp = sp.solve(sp.Eq(sp.simplify(u4_1.subs(D21, D21_comp)), 0), D41)[0]

    assert sp.simplify(D21_comp + u2 * D01) == 0
    assert sp.simplify(D41_comp - D4 * D01 / D0) == 0
    assert sp.simplify(D41_comp - (u2**2 - u4) * D01) == 0

    # One-pole reduction: on a one-pole branch D4/D0 = u2^2 - u4 with u4 = 4 u2^2,
    # so D4 = -3 D0 u2^2. Substitute into the general surface D41 = (D4/D0) D01.
    D4_one_pole = -3 * D0 * u2**2
    one_pole_D41 = sp.simplify((D41_comp).subs(D4, D4_one_pole))
    assert sp.simplify(one_pole_D41 - (-3 * u2**2) * D01) == 0

    print("\nVerified exact first-order conservative compensation surface:")
    print("D21 =", D21_comp)
    print("D41 =", D41_comp)
    print("On a one-pole base branch: D41 =", one_pole_D41)

    # ------------------------------------------------------------------
    # Primitive logarithmic slope compiler
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

    # Base bundle data
    C = kappa * lamB
    GU = lamU
    GW = kappa * lamW
    R = kappa * lamR

    Delta = OmU**2 * OmW**2 - R**2
    S2 = OmU**2 + OmW**2
    H = GU**2 + GW**2
    Q = GU**2 * OmW**2 + 2 * GU * GW * R + GW**2 * OmU**2
    P = OmU**2 * GW + R * GU

    B0_base = C**2 / varpi**2
    B2_base = C**2 / varpi**4
    B4_base = C**2 / varpi**6
    Z0_base = Q / Delta
    Z2_base = (Q * S2 - H * Delta) / Delta**2
    Z4_base = (Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3
    N0_base = P**2 / Delta**2

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

    B0_1 = sp.simplify(sp.diff(B0d, eps).subs(eps, 0))
    B2_1 = sp.simplify(sp.diff(B2d, eps).subs(eps, 0))
    B4_1 = sp.simplify(sp.diff(B4d, eps).subs(eps, 0))
    Z0_1 = sp.simplify(sp.diff(Z0d, eps).subs(eps, 0))
    Z2_1 = sp.simplify(sp.diff(Z2d, eps).subs(eps, 0))
    Z4_1 = sp.simplify(sp.diff(Z4d, eps).subs(eps, 0))
    N0_1 = sp.simplify(sp.diff(N0d, eps).subs(eps, 0))

    Delta1 = 2 * OmU**2 * OmW**2 * (xOU + xOW) - 2 * R**2 * xLR
    S2_1 = 2 * OmU**2 * xOU + 2 * OmW**2 * xOW
    H1 = 2 * GU**2 * xLU + 2 * GW**2 * xLW
    Q1 = (
        2 * GU**2 * OmW**2 * (xLU + xOW)
        + 2 * GU * GW * R * (xLU + xLW + xLR)
        + 2 * GW**2 * OmU**2 * (xLW + xOU)
    )
    P1_raw = OmU**2 * GW * (2 * xOU + xLW) + R * GU * (xLR + xLU)

    assert sp.simplify(B0_1 - B0_base * (2 * xLB - 2 * xV)) == 0
    assert sp.simplify(B2_1 - B2_base * (2 * xLB - 4 * xV)) == 0
    assert sp.simplify(B4_1 - B4_base * (2 * xLB - 6 * xV)) == 0

    assert sp.simplify(Z0_1 - (Q1 * Delta - Q * Delta1) / Delta**2) == 0
    assert sp.simplify(
        Z2_1
        - (Delta * (-Delta * H1 - H * Delta1 + Q * S2_1 + S2 * Q1) + 2 * Delta1 * (Delta * H - Q * S2))
        / Delta**3
    ) == 0
    assert sp.simplify(
        Z4_1
        + (
            Delta**2 * H * S2_1
            + Delta**2 * S2 * H1
            + Delta**2 * Q1
            - 2 * Delta * H * S2 * Delta1
            - 2 * Delta * Q * S2 * S2_1
            - 2 * Delta * Q * Delta1
            - Delta * S2**2 * Q1
            + 3 * Q * S2**2 * Delta1
        )
        / Delta**4
    ) == 0
    assert sp.simplify(N0_1 - (2 * P * P1_raw) / Delta**2 + 2 * P**2 * Delta1 / Delta**3) == 0

    D01_compiler = sp.simplify((sp.diff(Kd - B0d - Z0d, eps).subs(eps, 0)))
    D21_compiler = sp.simplify((sp.diff(-(Md + B2d + Z2d), eps).subs(eps, 0)))
    D41_compiler = sp.simplify((sp.diff(-(B4d + Z4d), eps).subs(eps, 0)))

    assert sp.simplify(D01_compiler - (K * xK - B0_1 - Z0_1)) == 0
    assert sp.simplify(D21_compiler + (M * xM + B2_1 + Z2_1)) == 0
    assert sp.simplify(D41_compiler + (B4_1 + Z4_1)) == 0

    print("\nVerified primitive logarithmic-slope compiler.")

    # ------------------------------------------------------------------
    # Concrete Stage 223 compatibility point
    # ------------------------------------------------------------------
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

    B0_s = sp.simplify(B0_base.subs(sample))
    B2_s = sp.simplify(B2_base.subs(sample))
    B4_s = sp.simplify(B4_base.subs(sample))
    Z0_s = sp.simplify(Z0_base.subs(sample))
    Z2_s = sp.simplify(Z2_base.subs(sample))
    Z4_s = sp.simplify(Z4_base.subs(sample))
    N0_s = sp.simplify(N0_base.subs(sample))

    D0_compat = sp.simplify((3 * (M + B2_s + Z2_s) ** 2 / (B4_s + Z4_s)).subs(sample))
    K_compat = sp.simplify(B0_s + Z0_s + D0_compat)
    P0_target_compat = sp.simplify(N0_s / D0_compat)

    sample_full = dict(sample)
    sample_full[K] = K_compat

    D0_s = sp.simplify((K - B0_base - Z0_base).subs(sample_full))
    D2_s = sp.simplify((-(M + B2_base + Z2_base)).subs(sample_full))
    D4_s = sp.simplify((-(B4_base + Z4_base)).subs(sample_full))
    u2_s = sp.simplify((-D2_s / D0_s))
    u4_s = sp.simplify((D2_s**2 - D0_s * D4_s) / D0_s**2)

    assert_close(float(D0_s.evalf()), 24.2373099886223)
    assert_close(float(D2_s.evalf()), -1.18562046858190)
    assert_close(float(D4_s.evalf()), -0.173991572849491)
    assert_close(float(u2_s.evalf()), 0.0489171640391802)
    assert_close(float(u4_s.evalf()), 0.00957155575054425)
    assert_close(float((D4_s / D0_s).evalf()), -0.00717866681290820)
    assert_close(float(P0_target_compat.evalf()), 0.00206979231806289)

    assert sp.simplify(u4_s - 4 * u2_s**2) == 0

    print("\nVerified concrete Stage 223 compatibility point:")
    print("K_compat               =", sp.N(K_compat, 16))
    print("D0                     =", sp.N(D0_s, 16))
    print("D2                     =", sp.N(D2_s, 16))
    print("D4                     =", sp.N(D4_s, 16))
    print("u2                     =", sp.N(u2_s, 16))
    print("u4                     =", sp.N(u4_s, 16))
    print("P0_target,compat       =", sp.N(P0_target_compat, 18))

    # ------------------------------------------------------------------
    # Sample-point Xi1 compiler
    # ------------------------------------------------------------------
    D01_s = sp.simplify(D01_compiler.subs(sample_full))
    D21_s = sp.simplify(D21_compiler.subs(sample_full))
    D41_s = sp.simplify(D41_compiler.subs(sample_full))
    N01_s = sp.simplify(N0_1.subs(sample_full))

    Xi1_s = sp.simplify(N01_s / N0_s - D01_s / D0_s)
    Xi1_coeffs = {
        xK: float(sp.N(sp.expand(Xi1_s).coeff(xK), 16)),
        xM: float(sp.N(sp.expand(Xi1_s).coeff(xM), 16)),
        xLB: float(sp.N(sp.expand(Xi1_s).coeff(xLB), 16)),
        xV: float(sp.N(sp.expand(Xi1_s).coeff(xV), 16)),
        xLU: float(sp.N(sp.expand(Xi1_s).coeff(xLU), 16)),
        xLW: float(sp.N(sp.expand(Xi1_s).coeff(xLW), 16)),
        xLR: float(sp.N(sp.expand(Xi1_s).coeff(xLR), 16)),
        xOU: float(sp.N(sp.expand(Xi1_s).coeff(xOU), 16)),
        xOW: float(sp.N(sp.expand(Xi1_s).coeff(xOW), 16)),
    }

    assert_close(Xi1_coeffs[xK], -1.00975540977030)
    assert_close(Xi1_coeffs[xM], 0.0)
    assert_close(Xi1_coeffs[xLB], 0.00418038073077834)
    assert_close(Xi1_coeffs[xV], -0.00418038073077834)
    assert_close(Xi1_coeffs[xLU], 0.324464020216766)
    assert_close(Xi1_coeffs[xLW], 1.69086641859305)
    assert_close(Xi1_coeffs[xLR], 0.423379354382463)
    assert_close(Xi1_coeffs[xOU], -0.747843374599229)
    assert_close(Xi1_coeffs[xOW], -4.11424577297551)

    print("\nConcrete Xi1 compiler at the compatibility point:")
    for var, coeff in Xi1_coeffs.items():
        print(f"  coeff[{var}] = {coeff}")

    # ------------------------------------------------------------------
    # Wall-only family: exact generic no-go on the concrete branch
    # ------------------------------------------------------------------
    eq_wall_1 = sp.simplify((D21_s + u2_s * D01_s).subs({
        xM: xM, xK: xK, xLB: 0, xV: 0, xLU: 0, xLW: 0, xLR: 0, xOU: 0, xOW: 0
    }))
    eq_wall_2 = sp.simplify((D41_s - (D4_s / D0_s) * D01_s).subs({
        xM: xM, xK: xK, xLB: 0, xV: 0, xLU: 0, xLW: 0, xLR: 0, xOU: 0, xOW: 0
    }))

    wall_solutions = sp.solve([sp.Eq(eq_wall_1, 0), sp.Eq(eq_wall_2, 0)], [xK, xM], dict=True)
    assert len(wall_solutions) == 1
    assert sp.simplify(wall_solutions[0][xK]) == 0
    assert sp.simplify(wall_solutions[0][xM]) == 0

    print("\nWall-only family:")
    print("  compensation solutions =", wall_solutions)

    # ------------------------------------------------------------------
    # Pure BdG family: exact determinant and sample-point no-go
    # ------------------------------------------------------------------
    B0g, B2g, B4g, D0g, D4g, u2g = sp.symbols("B0g B2g B4g D0g D4g u2g")
    M_BdG = sp.Matrix([
        [-(B2g + u2g * B0g), 2 * B2g + u2g * B0g],
        [-(B4g - D4g * B0g / D0g), 3 * B4g - D4g * B0g / D0g],
    ])
    det_BdG = sp.expand(M_BdG.det())
    expected_det = -B0g * B2g * D4g / D0g - 2 * B0g * B4g * u2g - B2g * B4g
    assert sp.simplify(det_BdG - expected_det) == 0

    det_BdG_num = float(sp.N(expected_det.subs({
        B0g: B0_s, B2g: B2_s, B4g: B4_s, D0g: D0_s, D4g: D4_s, u2g: u2_s
    }), 20))
    assert_close(det_BdG_num, -5.11886996120011e-5, 1e-15)

    print("\nPure BdG family determinant:")
    print("  Delta_BdG =", det_BdG)
    print("  Sample-point Delta_BdG =", det_BdG_num)

    # ------------------------------------------------------------------
    # Mixed/U family: surviving compensated nullspace
    # ------------------------------------------------------------------
    eq_mixed_1 = sp.expand((D21_s + u2_s * D01_s).subs({
        xK: 0, xM: 0, xLB: 0, xV: 0
    }))
    eq_mixed_2 = sp.expand((D41_s - (D4_s / D0_s) * D01_s).subs({
        xK: 0, xM: 0, xLB: 0, xV: 0
    }))

    mixed_vars = [xLU, xLW, xLR, xOU, xOW]
    row1 = [float(sp.N(eq_mixed_1.coeff(v), 18)) for v in mixed_vars]
    row2 = [float(sp.N(eq_mixed_2.coeff(v), 18)) for v in mixed_vars]

    for actual, expected in zip(row1, [
        -0.241952861865934, -0.122133861432532, -0.0656784156312263,
         0.553209522700447,  0.288144673113677,
    ]):
        assert_close(actual, expected, 1e-12)
    for actual, expected in zip(row2, [
        -0.250543086743604, -0.0937748521387244, -0.0899548469020231,
         0.881694465041011,  0.325834311088034,
    ]):
        assert_close(actual, expected, 1e-12)

    A = sp.Matrix([row1, row2])
    assert A.rank() == 2

    A11 = A[:, :2]
    A12 = A[:, 2:]
    A11_inv = A11.inv()

    basis = []
    for j in range(3):
        free = sp.Matrix([1 if i == j else 0 for i in range(3)])
        solved = -A11_inv * A12 * free
        vec = [float(sp.N(solved[0], 15)), float(sp.N(solved[1], 15))] + list(free)
        basis.append(vec)

    expected_basis = [
        [-0.610255553634424, 0.671187016268095, 1.0, 0.0, 0.0],
        [7.05469842496522, -9.44615143817664, 0.0, 1.0, 0.0],
        [1.61486053113911, -0.839860892848583, 0.0, 0.0, 1.0],
    ]
    for vec, exp_vec in zip(basis, expected_basis):
        for actual, expected in zip(vec, exp_vec):
            assert_close(actual, expected, 1e-12)

    Xi1_mixed = sp.expand(Xi1_s.subs({xK: 0, xM: 0, xLB: 0, xV: 0}))
    xi_basis = []
    for vec in basis:
        subs = {xLU: vec[0], xLW: vec[1], xLR: vec[2], xOU: vec[3], xOW: vec[4]}
        xi_val = float(sp.N(Xi1_mixed.subs(subs), 20))
        xi_basis.append(xi_val)

    assert_close(xi_basis[0], 1.36026097049402)
    assert_close(xi_basis[1], -14.4310278139755)
    assert_close(xi_basis[2], -5.01037421295998)

    print("\nMixed/U family compensation matrix:")
    print(A)
    print("Null-basis vectors:")
    for idx, vec in enumerate(basis, start=1):
        print(f"  v{idx} =", vec)
    print("Xi1 on basis:")
    for idx, val in enumerate(xi_basis, start=1):
        print(f"  Xi1(v{idx}) =", val)

    # ------------------------------------------------------------------
    # Transported Stage 224 headroom on the first surviving family
    # ------------------------------------------------------------------
    sigma1 = xi_basis[0]
    xi_budgets = {
        "both_10": 0.367930328492646,
        "one_10": 0.737619063660757,
        "both_30": 2.94889585703134,
        "one_30": 4.63505472371892,
    }
    t_budgets = {key: value / sigma1 for key, value in xi_budgets.items()}

    assert_close(t_budgets["both_10"], 0.270485102839510)
    assert_close(t_budgets["one_10"], 0.542262903708006)
    assert_close(t_budgets["both_30"], 2.16788978070904)
    assert_close(t_budgets["one_30"], 3.40747461278373)

    print("\nTransported headroom on the first surviving mixed family:")
    for key, value in t_budgets.items():
        print(f"  {key}: |eps t| <= {value}")

    print("\nAll Stage 225 symbolic and numerical checks passed.")


if __name__ == "__main__":
    main()
