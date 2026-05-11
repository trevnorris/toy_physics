import math
import sympy as sp


def assert_close(actual: float, expected: float, tol: float = 1e-12) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{actual} !~= {expected} (tol={tol})")


def vector_to_float_list(v: sp.Matrix, digits: int = 15) -> list[float]:
    return [float(sp.N(v[i], digits)) for i in range(v.rows)]


def coeff_vector(expr: sp.Expr, vars_: list[sp.Symbol]) -> sp.Matrix:
    expr = sp.expand(expr)
    return sp.Matrix([sp.N(expr.coeff(v), 40) for v in vars_])


def unit_vector(v: sp.Matrix) -> sp.Matrix:
    return sp.simplify(v / sp.sqrt((v.T * v)[0]))


def main() -> None:
    sp.init_printing()
    print("=== Stage 227 SymPy audit: pure-transfer load factor and outgoing-rigidity sieve ===")

    # ------------------------------------------------------------------
    # 1. Primitive one-port compiler on the explicit Stage 223 branch
    # ------------------------------------------------------------------
    eps = sp.symbols("eps", real=True)

    kappa = sp.symbols("kappa", positive=True)
    K, M = sp.symbols("K M", positive=True)
    lamB, varpi = sp.symbols("lamB varpi", positive=True)
    lamU, lamW, lamR = sp.symbols("lamU lamW lamR", positive=True)
    OmU, OmW = sp.symbols("OmU OmW", positive=True)

    xK, xM = sp.symbols("xK xM", real=True)
    xLB, xV = sp.symbols("xLB xV", real=True)
    xLU, xLW, xLR = sp.symbols("xLU xLW xLR", real=True)
    xOU, xOW = sp.symbols("xOU xOW", real=True)

    mixed_vars = [xLU, xLW, xLR, xOU, xOW]

    C = kappa * lamB
    GU = lamU
    GW = kappa * lamW
    R = kappa * lamR

    Delta = OmU**2 * OmW**2 - R**2
    S2 = OmU**2 + OmW**2
    H_mix = GU**2 + GW**2
    Q = GU**2 * OmW**2 + 2 * GU * GW * R + GW**2 * OmU**2
    P = OmU**2 * GW + R * GU

    B0 = C**2 / varpi**2
    B2 = C**2 / varpi**4
    B4 = C**2 / varpi**6
    Z0 = Q / Delta
    Z2 = (Q * S2 - H_mix * Delta) / Delta**2
    Z4 = (Q * (S2**2 - Delta) - S2 * H_mix * Delta) / Delta**3
    N0 = P**2 / Delta**2

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
    H_mix_d = GUd**2 + GWd**2
    Qd = GUd**2 * OmWd**2 + 2 * GUd * GWd * Rd + GWd**2 * OmUd**2
    Pd = OmUd**2 * GWd + Rd * GUd

    B0d = Cd**2 / vard**2
    B2d = Cd**2 / vard**4
    B4d = Cd**2 / vard**6
    Z0d = Qd / Deltad
    Z2d = (Qd * S2d - H_mix_d * Deltad) / Deltad**2
    Z4d = (Qd * (S2d**2 - Deltad) - S2d * H_mix_d * Deltad) / Deltad**3
    N0d = Pd**2 / Deltad**2

    D01 = sp.simplify(sp.diff(Kd - B0d - Z0d, eps).subs(eps, 0))
    D21 = sp.simplify(sp.diff(-(Md + B2d + Z2d), eps).subs(eps, 0))
    D41 = sp.simplify(sp.diff(-(B4d + Z4d), eps).subs(eps, 0))
    N01 = sp.simplify(sp.diff(N0d, eps).subs(eps, 0))
    P01 = sp.simplify(sp.diff(Pd, eps).subs(eps, 0))
    Delta01 = sp.simplify(sp.diff(Deltad, eps).subs(eps, 0))

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

    D0_compat = sp.simplify(3 * (sp.Integer(1) + B2_s + Z2_s) ** 2 / (B4_s + Z4_s))
    K_compat = sp.simplify(B0_s + Z0_s + D0_compat)

    sample_full = dict(sample)
    sample_full[K] = K_compat

    zero_nonmixed = {xK: 0, xM: 0, xLB: 0, xV: 0}

    D01_mixed = sp.expand(sp.simplify(D01.subs(sample_full))).subs(zero_nonmixed)
    D21_mixed = sp.expand(sp.simplify(D21.subs(sample_full))).subs(zero_nonmixed)
    D41_mixed = sp.expand(sp.simplify(D41.subs(sample_full))).subs(zero_nonmixed)
    N01_mixed = sp.expand(sp.simplify(N01.subs(sample_full))).subs(zero_nonmixed)
    P01_mixed = sp.expand(sp.simplify(P01.subs(sample_full))).subs(zero_nonmixed)
    Delta01_mixed = sp.expand(sp.simplify(Delta01.subs(sample_full))).subs(zero_nonmixed)

    P_s = sp.simplify(P.subs(sample_full))
    Delta_s = sp.simplify(Delta.subs(sample_full))
    N0_s = sp.simplify(N0.subs(sample_full))

    # ------------------------------------------------------------------
    # 2. Pure-transfer corridor from Stage 226
    # ------------------------------------------------------------------
    M_transfer = sp.Matrix([
        [sp.simplify(D01_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(D21_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(D41_mixed.coeff(v)) for v in mixed_vars],
    ])
    assert M_transfer.rank() == 3

    transfer_null = M_transfer.nullspace()
    assert len(transfer_null) == 2
    T = sp.Matrix.hstack(*transfer_null)

    expected_t = [
        [-4.359222794718, 3.107402039105, 18.703510605854, 1.0, 0.0],
        [1.909256655687, -1.163651238154, -0.482414494705, 0.0, 1.0],
    ]
    for basis_vec, exp_basis in zip(transfer_null, expected_t):
        vals = vector_to_float_list(basis_vec)
        for actual, expected in zip(vals, exp_basis):
            assert_close(actual, expected, tol=1e-11)

    print("\nVerified Stage 226 pure-transfer corridor basis:")
    for idx, basis_vec in enumerate(transfer_null, start=1):
        print(f"t_{idx} =", vector_to_float_list(basis_vec, 15))

    # ------------------------------------------------------------------
    # 3. Exact pure-transfer theorem
    # ------------------------------------------------------------------
    Xi_transfer = sp.simplify(N01_mixed / N0_s)

    assert sp.simplify(Xi_transfer - 2 * (P01_mixed / P_s - Delta01_mixed / Delta_s)) == 0

    Lambda = sp.simplify(P / Delta)
    I = sp.simplify(R * GU / (OmU**2 * GW))
    H = sp.simplify(R**2 / (OmU**2 * OmW**2))
    assert sp.simplify(Lambda - GW / OmW**2 * (1 + I) / (1 - H)) == 0

    # Microscopic drift coordinates
    m_expr = xLW - 2 * xOW
    i_expr = xLR + xLU - xLW - 2 * xOU
    h_expr = 2 * xLR - 2 * xOU - 2 * xOW

    I_s = sp.simplify(I.subs(sample_full))
    H_s = sp.simplify(H.subs(sample_full))

    Xi_mih = sp.simplify(2 * (m_expr + I_s / (1 + I_s) * i_expr + H_s / (1 - H_s) * h_expr))
    assert sp.simplify(Xi_transfer - Xi_mih) == 0

    Xi_sample_specialized = sp.simplify(Xi_transfer - (2 * m_expr + sp.Rational(6, 19) * i_expr + 50 * h_expr / (98 * sp.pi**2 - 25)))
    assert Xi_sample_specialized == 0

    assert I_s == sp.Rational(3, 16)
    assert H_s == sp.Rational(25, 98) / sp.pi**2

    print("\nVerified pure-transfer theorem:")
    print("Xi_1 =", Xi_transfer)
    print("Lambda =", sp.simplify(Lambda.subs(sample_full)))
    print("I =", I_s)
    print("H =", H_s)
    print("Xi_1 = 2[m + I/(1+I) i + H/(1-H) h]")
    print("Xi_1 sample law =", sp.simplify(2 * m_expr + sp.Rational(6, 19) * i_expr + 50 * h_expr / (98 * sp.pi**2 - 25)))

    # ------------------------------------------------------------------
    # 4. Rigidity sieve on the pure-transfer corridor
    # ------------------------------------------------------------------
    i_row = sp.Matrix([[sp.simplify(i_expr.coeff(v)) for v in mixed_vars]])
    h_row = sp.Matrix([[sp.simplify(h_expr.coeff(v)) for v in mixed_vars]])
    m_row = sp.Matrix([[sp.simplify(m_expr.coeff(v)) for v in mixed_vars]])

    red_i = sp.simplify(i_row * T)
    red_h = sp.simplify(h_row * T)
    red_m = sp.simplify(m_row * T)

    det_ih = sp.factor(sp.simplify(sp.Matrix.vstack(red_i, red_h).det()))
    expected_det = -sp.Integer(19) * (-25 + 98 * sp.pi**2) * (200 + 147 * sp.pi**2) * (441 * sp.pi**2 + 4400) / (6 * (8670000 + 14894275 * sp.pi**2 + 2117682 * sp.pi**4))
    assert sp.simplify(det_ih - expected_det) == 0
    assert det_ih != 0

    print("\nVerified combined i=h rigidity determinant:")
    print("det[(i,h)|_pure transfer] =", det_ih)

    # ------------------------------------------------------------------
    # 5. One-dimensional rigid survivors
    # ------------------------------------------------------------------
    ui = red_i.nullspace()[0]
    uh = red_h.nullspace()[0]
    um = red_m.nullspace()[0]

    vi = unit_vector(T * ui)
    vh = unit_vector(T * uh)
    vm = -unit_vector(T * um)  # sign chosen to match the note convention

    expected_vi = [0.45280825, -0.29424612, -0.82815170, -0.04054866, 0.14458380]
    expected_vh = [0.66561963, -0.38941932, 0.46712837, 0.03609301, 0.43103536]
    expected_vm = [0.13386239, -0.10586713, -0.98242900, -0.05389175, -0.05293356]

    for actual, expected in zip(vector_to_float_list(vi), expected_vi):
        assert_close(actual, expected, tol=1e-8)
    for actual, expected in zip(vector_to_float_list(vh), expected_vh):
        assert_close(actual, expected, tol=1e-8)
    for actual, expected in zip(vector_to_float_list(vm), expected_vm):
        assert_close(actual, expected, tol=1e-8)

    Xi_coeff = coeff_vector(Xi_transfer, mixed_vars)

    Xi_vi = float(sp.N((Xi_coeff.T * vi)[0], 18))
    Xi_vh = float(sp.N((Xi_coeff.T * vh)[0], 18))
    Xi_vm = float(sp.N((Xi_coeff.T * vm)[0], 18))

    assert_close(abs(Xi_vi), 1.26576248, tol=1e-8)
    assert_close(abs(Xi_vh), 2.04509123, tol=1e-8)
    assert_close(abs(Xi_vm), 0.29342952, tol=1e-8)

    print("\nVerified rigid one-dimensional survivors:")
    print("v_i =", vector_to_float_list(vi, 15), " |Xi_1| =", abs(Xi_vi))
    print("v_h =", vector_to_float_list(vh, 15), " |Xi_1| =", abs(Xi_vh))
    print("v_m =", vector_to_float_list(vm, 15), " |Xi_1| =", abs(Xi_vm))

    # ------------------------------------------------------------------
    # 6. Corridor norms and transported 10%-loss ceilings
    # ------------------------------------------------------------------
    G = sp.simplify(T.T * T)
    proj_transfer = sp.simplify(T * G.inv() * T.T)
    sigma_transfer = float(sp.N(sp.sqrt((Xi_coeff.T * proj_transfer * Xi_coeff)[0]), 20))
    sigma_i = abs(Xi_vi)
    sigma_h = abs(Xi_vh)
    sigma_m = abs(Xi_vm)

    assert_close(sigma_transfer, 2.31561904386057, tol=1e-11)

    budget_both = 0.367930328492646
    budget_nonempty = 0.737619063660757

    ceilings = {
        "transfer": (budget_both / sigma_transfer, budget_nonempty / sigma_transfer),
        "i": (budget_both / sigma_i, budget_nonempty / sigma_i),
        "h": (budget_both / sigma_h, budget_nonempty / sigma_h),
        "m": (budget_both / sigma_m, budget_nonempty / sigma_m),
    }

    expected_ceilings = {
        "transfer": (0.15889070, 0.31854077),
        "i": (0.29067881, 0.58274682),
        "h": (0.17990900, 0.36067783),
        "m": (1.25389678, 2.51378617),
    }

    for key, (c1, c2) in ceilings.items():
        e1, e2 = expected_ceilings[key]
        assert_close(c1, e1, tol=1e-8)
        assert_close(c2, e2, tol=1e-8)

    print("\nVerified corridor norms and transported stricter 10%-loss ceilings:")
    print("sigma_transfer =", sigma_transfer, " ceilings =", ceilings["transfer"])
    print("sigma_i        =", sigma_i, " ceilings =", ceilings["i"])
    print("sigma_h        =", sigma_h, " ceilings =", ceilings["h"])
    print("sigma_m        =", sigma_m, " ceilings =", ceilings["m"])

    print("\nStage 227 audit completed successfully.")


if __name__ == "__main__":
    main()
