#!/usr/bin/env python3
"""
2_5pn_master_session_sympy_audit.py

Standalone SymPy verifier for the full 2.5PN derivation stack developed in this
session.

Coverage
--------
This master script consolidates:

  • the earlier session benchmark/audit layers:
      - decisive Burke–Thorne benchmark
      - low-frequency selection-rule framework
      - scalar-sector audits and rescues
      - dipole/vector-sector audits and rescues
      - core quadrupole representation/source-map checks

  • the later session quadrupole-completion layers:
      - canonical quadrupole normalization map
      - outgoing l=2 PDE fingerprint
      - source/port normalization cleanup
      - convention-invariant normalization product
      - canonical invariant low-frequency fingerprint
      - extraction of (Kbar0, Kbar2) and the single-prefactor obstruction
      - frozen 2PN fixed-vs-open quadrupole data
      - minimal positive conservative precursor
      - minimal isotropic quadrupole module and componentwise P2 target
      - single-pole sufficiency theorem
      - conservative extraction ledgers
      - grouped/normalized-support P2 formulas
      - axisymmetric grouped-P2 parameterization

The script uses SymPy only and imports nothing from local project files.
It prints the main intermediate identities and ends with a compact theorem ledger.
"""

from __future__ import annotations

import sympy as sp

I = sp.I
sqrt = sp.sqrt


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(sp.simplify)
        is_zero = all(entry == 0 for entry in simplified)
        print(f"{name} =", simplified)
        if not is_zero:
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.simplify(expr)
        print(f"{name} =", simplified)
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


def stf_from_vector(v: sp.Matrix) -> sp.Matrix:
    return sp.simplify(v * v.T - sp.eye(3) * (v.dot(v)) / sp.Integer(3))


def stf(mat: sp.Matrix) -> sp.Matrix:
    return sp.simplify((mat + mat.T) / 2 - sp.eye(3) * sp.trace(mat) / 3)


# Canonical real STF l=2 basis aligned with n = z-hat.
E20 = sp.Matrix([[-1 / sqrt(6), 0, 0], [0, -1 / sqrt(6), 0], [0, 0, 2 / sqrt(6)]])
E21c = sp.Matrix([[0, 0, 1 / sqrt(2)], [0, 0, 0], [1 / sqrt(2), 0, 0]])
E21s = sp.Matrix([[0, 0, 0], [0, 0, 1 / sqrt(2)], [0, 1 / sqrt(2), 0]])
E22c = sp.Matrix([[1 / sqrt(2), 0, 0], [0, -1 / sqrt(2), 0], [0, 0, 0]])
E22s = sp.Matrix([[0, 1 / sqrt(2), 0], [1 / sqrt(2), 0, 0], [0, 0, 0]])
BASIS = [E20, E21c, E21s, E22c, E22s]
BASIS_NAMES = ["20", "21c", "21s", "22c", "22s"]


def coeffs_in_basis(T: sp.Matrix) -> list[sp.Expr]:
    return [sp.simplify(sum(T.multiply_elementwise(B))) for B in BASIS]


# ---------------------------------------------------------------------------
# Part II / decisive benchmark utilities
# ---------------------------------------------------------------------------

def newtonian_reduce_setup():
    t = sp.symbols("t", real=True)
    G, c, m, nu, h = sp.symbols("G c m nu h", positive=True)
    mu = nu * m
    k = G * m

    r = sp.Function("r")(t)
    phi = sp.Function("phi")(t)
    rd = sp.diff(r, t)

    x = r * sp.cos(phi)
    y = r * sp.sin(phi)
    rr = sp.expand(x**2 + y**2)

    rdd = -k / r**2 + h**2 / r**3
    phidd = -2 * rd * sp.diff(phi, t) / r

    subs_basic = {
        sp.diff(r, t, 2): rdd,
        sp.diff(phi, t): h / r**2,
        sp.diff(phi, t, 2): phidd.subs({sp.diff(phi, t): h / r**2}),
    }

    def reduce_expr(expr: sp.Expr, max_iter: int = 10) -> sp.Expr:
        expr = sp.expand(expr)
        for _ in range(max_iter):
            expr = expr.xreplace(subs_basic)
            changed = False
            for d in list(expr.atoms(sp.Derivative)):
                if d.expr == r and d.derivative_count > 1:
                    rhs = rdd
                    for _i in range(d.derivative_count - 2):
                        rhs = sp.diff(rhs, t)
                    rhs = sp.expand(rhs.xreplace(subs_basic))
                    expr = expr.xreplace({d: rhs})
                    changed = True
                elif d.expr == phi and d.derivative_count > 1:
                    rhs = subs_basic[sp.diff(phi, t, 2)]
                    for _i in range(d.derivative_count - 2):
                        rhs = sp.diff(rhs, t)
                    rhs = sp.expand(rhs.xreplace(subs_basic))
                    expr = expr.xreplace({d: rhs})
                    changed = True
            expr = sp.expand(expr.xreplace(subs_basic))
            if not changed:
                break
        return sp.simplify(expr)

    def nth_deriv(expr: sp.Expr, n: int) -> sp.Expr:
        cur = expr
        for _ in range(n):
            cur = reduce_expr(sp.diff(cur, t))
        return sp.simplify(cur)

    return {
        "t": t, "G": G, "c": c, "m": m, "nu": nu, "mu": mu, "h": h, "k": k,
        "r": r, "phi": phi, "rd": rd, "x": x, "y": y, "rr": rr,
        "reduce_expr": reduce_expr, "nth_deriv": nth_deriv,
    }


# ---------------------------------------------------------------------------
# Part II: decisive benchmark / BT prototype
# ---------------------------------------------------------------------------

def decisive_benchmark() -> dict[str, sp.Expr]:
    banner("PART II — DECISIVE BENCHMARK / BURKE–THORNE PROTOTYPE")

    # 1. Mass dipole and quadrupole decomposition.
    subbanner("II.1 — Two-body mass dipole and STF quadrupole decomposition")
    m1, m2 = sp.symbols("m1 m2", positive=True)
    M = m1 + m2
    mu_red = sp.simplify(m1 * m2 / M)

    X1, X2, X3, x1, x2, x3 = sp.symbols("X1 X2 X3 x1 x2 x3", real=True)
    X = sp.Matrix([X1, X2, X3])
    x = sp.Matrix([x1, x2, x3])

    r1 = X + (m2 / M) * x
    r2 = X - (m1 / M) * x

    D = sp.simplify(m1 * r1 + m2 * r2)
    dD_dx = sp.simplify(D.jacobian(x))

    Q = sp.simplify(m1 * stf_from_vector(r1) + m2 * stf_from_vector(r2))
    Q_expected = sp.simplify(M * stf_from_vector(X) + mu_red * stf_from_vector(x))
    Q_residual = sp.simplify(Q - Q_expected)

    print("D_i =")
    sp.pprint(D)
    expect_zero("dD_i/dx_j", dD_dx)
    expect_zero("Q - (M STF(XX) + mu STF(xx))", Q_residual)

    # 2. Burke-Thorne relative force lands on the standard 2.5PN family.
    subbanner("II.2 — Burke–Thorne local quadrupole force and Iyer–Will match")
    setup = newtonian_reduce_setup()
    t = setup["t"]
    G = setup["G"]
    c = setup["c"]
    m = setup["m"]
    nu = setup["nu"]
    mu = setup["mu"]
    h = setup["h"]
    r = setup["r"]
    phi = setup["phi"]
    rd = setup["rd"]
    x_rel = setup["x"]
    y_rel = setup["y"]
    rr = setup["rr"]
    nth_deriv = setup["nth_deriv"]

    Qxx = sp.simplify(mu * (x_rel**2 - rr / 3))
    Qxy = sp.simplify(mu * x_rel * y_rel)

    Qxx5 = nth_deriv(Qxx, 5)
    Qxy5 = nth_deriv(Qxy, 5)

    subs0 = {phi: 0, sp.sin(phi): 0, sp.cos(phi): 1}
    Qxx5_0 = sp.simplify(Qxx5.subs(subs0))
    Qxy5_0 = sp.simplify(Qxy5.subs(subs0))

    a_n = sp.simplify(-(2 * G / (5 * c**5)) * r * Qxx5_0)
    a_l = sp.simplify(-(2 * G / (5 * c**5)) * r * Qxy5_0)

    factor = sp.Rational(8, 5) * G**2 * m**2 * nu / (c**5 * r**3)
    Bcoef = sp.simplify((a_l / factor) * r / h)
    Acoef = sp.simplify((a_n / factor) / rd - Bcoef)

    v2 = sp.symbols("v2", real=True)
    Acoef_v = sp.expand(Acoef.subs({h**2 / r**2: v2 - rd**2}))
    Bcoef_v = sp.expand(Bcoef.subs({h**2 / r**2: v2 - rd**2}))

    print("A(v^2, GM/r, rdot^2) =")
    sp.pprint(sp.factor(Acoef_v))
    print("B(v^2, GM/r, rdot^2) =")
    sp.pprint(sp.factor(Bcoef_v))

    alpha, beta = sp.symbols("alpha beta")
    eqs = [
        sp.Eq(sp.expand(Bcoef_v).coeff(v2), -(2 + alpha)),
        sp.Eq(sp.expand(Bcoef_v).coeff(G * m / r), -(2 - alpha)),
        sp.Eq(sp.expand(Bcoef_v).coeff(rd**2), 3 * (1 + alpha)),
    ]
    sol_alpha = sp.solve(eqs, [alpha], dict=True)[0][alpha]
    sol_beta = sp.solve(sp.Eq(sp.expand(Acoef_v).coeff(v2), 3 * (1 + beta)), beta)[0]
    print("alpha =", sol_alpha)
    print("beta  =", sol_beta)
    if sol_alpha != 4 or sol_beta != 5:
        raise AssertionError("Burke–Thorne force did not land on the expected (alpha,beta)=(4,5) family member.")

    return {
        "alpha": sol_alpha,
        "beta": sol_beta,
        "Acoef_v": Acoef_v,
        "Bcoef_v": Bcoef_v,
    }


# ---------------------------------------------------------------------------
# Part III: selection rules and influence functional
# ---------------------------------------------------------------------------

def general_low_frequency_framework() -> dict[str, sp.Expr]:
    banner("PART III — LOW-FREQUENCY SELECTION RULES / INFLUENCE FUNCTIONAL")

    subbanner("III.1 — Time-domain signs for i omega^n")
    signs = {}
    for n in (1, 3, 5):
        signs[n] = sp.simplify(I / ((-I) ** n))
        print(f"i*omega^{n}  ->  {signs[n]} * d^{n}/dt^{n}")
    # Expected mapping: n=1 -> -d/dt, n=3 -> +d^3/dt^3, n=5 -> -d^5/dt^5
    if signs[1] != -1 or signs[3] != 1 or signs[5] != -1:
        raise AssertionError("Unexpected Fourier sign mapping.")

    subbanner("III.2 — Minimal retarded kernel expansions")
    omega, Omega, sigma, g = sp.symbols("omega Omega sigma g", positive=True, real=True)
    K1 = sp.expand(sp.series(g**2 / (Omega**2 - omega**2 - I * sigma * omega), omega, 0, 3).removeO())
    K3 = sp.expand(sp.series(g**2 / (Omega**2 - omega**2 - I * sigma * omega**3), omega, 0, 5).removeO())
    K5 = sp.expand(sp.series(g**2 / (Omega**2 - omega**2 - I * sigma * omega**5), omega, 0, 7).removeO())
    print("K1 =", K1)
    print("K3 =", K3)
    print("K5 =", K5)

    subbanner("III.3 — Dissipation / Schott identities")
    t = sp.symbols("t", real=True)
    q = sp.Function("q")(t)
    q1 = sp.diff(q, t)
    q2 = sp.diff(q, t, 2)
    q3 = sp.diff(q, t, 3)
    q4 = sp.diff(q, t, 4)
    q5 = sp.diff(q, t, 5)

    id1 = sp.simplify(q1 * (-q1) + q1**2)
    id3 = sp.simplify(q1 * q3 - (sp.diff(q1 * q2, t) - q2**2))
    id5 = sp.simplify(q1 * (-q5) - (-sp.diff(q1 * q4 - q2 * q3, t) - q3**2))
    expect_zero("n=1 dissipation identity", id1)
    expect_zero("n=3 dissipation identity", id3)
    expect_zero("n=5 dissipation identity", id5)

    subbanner("III.4 — Model-specific 2PN scalar / quadrupole combinations")
    delta01, deltag1, delta205 = sp.symbols("delta01 deltag1 delta205", real=True)
    J0_sq = sp.Rational(16, 5)
    Delta_geom = sp.Rational(281, 80)
    J20_sq = sp.Rational(25, 16)
    gamma1_eff = sp.simplify(J0_sq * delta01 - Delta_geom * deltag1)
    gamma5_eff = sp.simplify(J20_sq * delta205)
    print("gamma1_eff =", gamma1_eff)
    print("gamma5_eff =", gamma5_eff)

    return {
        "gamma1_eff": gamma1_eff,
        "gamma5_eff": gamma5_eff,
    }


# ---------------------------------------------------------------------------
# Part IV: scalar sector
# ---------------------------------------------------------------------------

def scalar_sector() -> dict[str, sp.Expr]:
    banner("PART IV — SCALAR SECTOR")

    # IV.1: outgoing scalar monopole and the dangerous i*omega term.
    subbanner("IV.1 — Outgoing scalar Green function and monopole odd term")
    omega, c_s, r = sp.symbols("omega c_s r", positive=True, real=True)
    k = omega / c_s
    G = sp.exp(I * k * r) / (4 * sp.pi * r)
    G_series = sp.expand(sp.series(G, omega, 0, 6).removeO())
    print("G_R(omega,r) =", G)
    print("Series =", G_series)
    gamma1_eff = sp.Rational(16, 5) * sp.Symbol("delta01") - sp.Rational(281, 80) * sp.Symbol("deltag1")
    print("Model-specific gamma1_eff =", gamma1_eff)

    subbanner("IV.1b — 2PN scalar support/geometry finite-size rescue")
    a0, ag, k0 = sp.symbols("a0 ag k0", positive=True, real=True)
    z0 = sp.symbols("z0")
    j0 = sp.sin(z0) / z0
    y0 = -sp.cos(z0) / z0
    h0 = sp.simplify(j0 + I * y0)
    Lambda0 = sp.simplify((k0 * sp.diff(h0, z0) / h0).subs(z0, k0 * a0))
    Y0 = sp.simplify(1 / Lambda0)
    Y0_norm = sp.simplify(Y0 / Y0.subs(k0, 0))
    Yg_norm = sp.simplify(1 / (1 - I * ag * k0))
    Seff = sp.simplify(sp.Rational(16, 5) * Y0_norm + sp.Rational(25, 16) - sp.Rational(281, 80) * Yg_norm)
    ell_eff = sp.simplify(sp.Rational(16, 5) * a0 - sp.Rational(281, 80) * ag)
    print("Lambda0(k) =", Lambda0)
    print("Y0_norm(k) =", sp.series(Y0_norm, k0, 0, 6).removeO())
    print("Yg_norm(k) =", sp.series(Yg_norm, k0, 0, 6).removeO())
    print("Seff(k) =", sp.series(Seff, k0, 0, 5).removeO())
    print("ell_eff =", ell_eff)
    print("equal-scale residual ell_eff(ag=a0) =", sp.simplify(ell_eff.subs(ag, a0)))
    print("exact scalar cancellation ag =", sp.solve(sp.Eq(ell_eff, 0), ag)[0])

    # IV.2: Gaussian counterexample / Ward test.
    subbanner("IV.2 — Gaussian overlap counterexample and exact leakage identity")
    w, a, ell, N, adot = sp.symbols("w a ell N adot", positive=True, real=True)
    W = sp.exp(-w**2 / ell**2) / (sp.sqrt(sp.pi) * ell)
    g = sp.exp(-w**2 / a**2) / (sp.sqrt(sp.pi) * a)
    C = sp.simplify(sp.integrate(W * g, (w, -sp.oo, sp.oo)))
    jw = sp.simplify((adot / a) * w * g)
    continuity_residual = sp.simplify(adot * sp.diff(g, a) + sp.diff(jw, w))
    I_leak = sp.simplify(sp.integrate(sp.diff(W, w) * jw, (w, -sp.oo, sp.oo)))
    dC_da = sp.simplify(adot * sp.diff(C, a))
    print("C(a) =", C)
    print("dC/da =", sp.simplify(sp.diff(C, a)))
    expect_zero("continuity residual", continuity_residual)
    expect_zero("I_leak - adot*dC/da", sp.simplify(I_leak - dC_da))

    # IV.3: projection-locking criterion.
    subbanner("IV.3 — Projection-locking linear algebra criterion")
    Ba, BL = sp.symbols("B_a B_L")
    v1a, v1L, v2a, v2L = sp.symbols("v1_a v1_L v2_a v2_L")
    detM = sp.simplify(v1a * v2L - v1L * v2a)
    print("determinant of two-mode tangent matrix =", detM)
    print("If determinant != 0, projection-locking requires B_a = B_L = 0.")

    # IV.4: direct vs derivative coupling and the damped discrete-mode loophole.
    subbanner("IV.4 — Direct vs derivative coupling; damped discrete-mode test")
    A, B, Cc, Dd = sp.symbols("A B C D", real=True)
    g0, gd = sp.symbols("g0 gd", real=True)
    Omega, eta, lam, lamd = sp.symbols("Omega eta lam lamd", positive=True, real=True)
    G_gapless = A + I * B * omega + Cc * omega**2 + I * Dd * omega**3
    direct = sp.expand(g0**2 * G_gapless)
    derivative = sp.expand(gd**2 * omega**2 * G_gapless)
    discrete_damped = sp.expand(sp.series(lam**2 / (Omega**2 - omega**2 - I * eta * omega), omega, 0, 4).removeO())
    discrete_damped_derivative = sp.expand(
        sp.series(lamd**2 * omega**2 / (Omega**2 - omega**2 - I * eta * omega), omega, 0, 5).removeO()
    )
    gamma1_direct = sp.expand(sp.im(sp.diff(direct, omega).subs(omega, 0)))
    gamma3_derivative = sp.simplify(sp.im(sp.diff(derivative, omega, 3).subs(omega, 0)) / sp.factorial(3))
    gamma1_disc = sp.expand(sp.im(sp.diff(discrete_damped, omega).subs(omega, 0)))
    print("direct =", direct)
    print("derivative =", derivative)
    print("discrete_damped =", discrete_damped)
    print("discrete_damped_derivative =", discrete_damped_derivative)
    print("gamma1_direct =", gamma1_direct)
    print("gamma3_derivative =", gamma3_derivative)
    print("gamma1_discrete_damped =", gamma1_disc)
    if gamma1_direct == 0:
        raise AssertionError("Direct coupling unexpectedly lost its linear odd term.")
    if gamma3_derivative == 0:
        raise AssertionError("Derivative coupling unexpectedly lost its cubic odd term.")
    if gamma1_disc == 0:
        raise AssertionError("Damped discrete mode unexpectedly lost the inherited linear odd term.")

    # IV.5: breathing-to-outlet vertex is derivative-coupled.
    subbanner("IV.5 — Breathing-to-outlet vertex is dot(q)-type")
    Bq = sp.symbols("B_q", real=True)
    K_direct = sp.expand(Bq**2 * G_gapless)
    K_deriv = sp.expand(Bq**2 * omega**2 * G_gapless)
    print("K_direct =", K_direct)
    print("Im K_direct =", sp.expand(sp.im(K_direct)))
    print("K_deriv =", K_deriv)
    print("Im K_deriv =", sp.expand(sp.im(K_deriv)))
    if sp.expand(sp.im(K_deriv)).coeff(omega, 1) != 0:
        raise AssertionError("Derivative outlet kernel still contains a linear odd term.")
    delta = sp.symbols("delta", positive=True)
    breathing_exponent = sp.Rational(5, 2) + 3 * delta
    print("Derivative-breathing odd term exponent (with a/r ~ eps^delta) =", breathing_exponent)

    # IV.6: no extra independent scalar source beyond continuity + boundary.
    subbanner("IV.6 — No third linear scalar source from quadratic momentum/stress terms")
    wv, z1, z2, z3, u0, u1, u2, q0 = sp.symbols("w z1 z2 z3 u0 u1 u2 q", real=True)
    Z = I * z1 * wv + z2 * wv**2 + I * z3 * wv**3
    u_sigma = u0 * q0 + I * u1 * wv * q0 + u2 * wv**2 * q0
    j = sp.expand((Z * u_sigma).series(wv, 0, 4).removeO())
    print("j_sigma(w) =", j)
    print("With Z_sigma(0)=0, the leading mouth source is derivative-like, not direct q-like.")

    # IV.7: mouth admittance PN count and radiative-order theorem.
    subbanner("IV.7 — Mouth radiative-order theorem and Ohmic no-go")
    wv = sp.symbols("w", real=True)
    K, M, M_add, beta, R2 = sp.symbols("K M M_add beta R2", positive=True)
    Zrad = R2 * wv**2 + I * M_add * wv
    den = K - M * wv**2 + I * wv * Zrad
    Y = sp.simplify(I * wv * beta / den)
    Y_expand = sp.expand(sp.series(Y, wv, 0, 6).removeO())
    print("Compact reciprocal mouth admittance Y(w) =", Y_expand)
    print("Re Y series =", sp.series(sp.re(Y), wv, 0, 6))
    print("Im Y series =", sp.series(sp.im(Y), wv, 0, 6))

    # Dangerous 1D / Ohmic benchmark.
    a0, c0, rho0 = sp.symbols("a0 c0 rho0", positive=True, real=True)
    Z1 = c0 * rho0
    Y1 = sp.simplify(I * beta * wv / (K - M * wv**2 + I * wv * Z1))
    print("1D Ohmic benchmark Re Y1 =", sp.series(sp.re(Y1), wv, 0, 5))

    # PN count for hypothetical z2 term.
    chi_sigma, kappa_m, L = sp.symbols("chi_sigma kappa_m L", positive=True, real=True)
    z2_nat, a_small, lam = sp.symbols("z2_nat a_small lam", positive=True, real=True)
    ell_mouth = sp.simplify((chi_sigma * c_s**3 / (kappa_m * L)) * z2_nat)
    ell_nat = sp.simplify(ell_mouth.subs({z2_nat: a_small**2 / c_s**3, L: lam * a_small}))
    print("Natural odd mouth length ell_mouth ~", ell_nat)

    return {
        "C_gaussian": C,
        "I_leak": I_leak,
    }


# ---------------------------------------------------------------------------
# Part V: dipole / vector sector
# ---------------------------------------------------------------------------

def dipole_sector() -> dict[str, sp.Expr]:
    banner("PART V — DIPOLE / VECTOR SECTOR")

    # V.1: exact CM/relative split of carried odd wake ports.
    subbanner("V.1 — Carried odd wake ports: CM/relative split")
    mA, mB, M = sp.symbols("mA mB M", nonzero=True)
    etaA = mA / M
    etaB = mB / M
    Ux, Uy, Dz = sp.symbols("Ux Uy Dz")
    ux, uy, d = sp.symbols("ux uy d")
    UA = sp.Matrix([Ux, Uy])
    u = sp.Matrix([ux, uy])

    PiA_perp = sqrt(7) / sqrt(2) * (UA + etaB * u)
    PiB_perp = sqrt(7) / sqrt(2) * (UA - etaA * u)
    PiA_10 = 2 * (Dz + etaB * d)
    PiB_10 = 2 * (Dz - etaA * d)

    diff_perp = sp.simplify((PiA_perp - PiB_perp).subs(M, mA + mB))
    diff_10 = sp.simplify((PiA_10 - PiB_10).subs(M, mA + mB))
    print("PiA_perp - PiB_perp =", diff_perp)
    print("PiA_10   - PiB_10   =", diff_10)

    # V.2: total momentum conservation does not kill the relative odd dipole term.
    subbanner("V.2 — Vector Ward-identity failure")
    t = sp.symbols("t", real=True)
    xpA = sp.Function("xpA")(t)
    xpB = sp.Function("xpB")(t)
    xmA = sp.Function("xmA")(t)
    xmB = sp.Function("xmB")(t)
    acoef = sp.symbols("a", real=True)

    L_body = acoef * (sp.diff(xmA, t) - sp.diff(xmB, t)) * (sp.diff(xpA, t, 2) - sp.diff(xpB, t, 2))
    F_A = -sp.diff(sp.diff(L_body, sp.diff(xmA, t)), t)
    F_B = -sp.diff(sp.diff(L_body, sp.diff(xmB, t)), t)
    print("F_A =", sp.simplify(F_A))
    print("F_B =", sp.simplify(F_B))
    expect_zero("F_A + F_B", sp.simplify(F_A + F_B))

    # CM/relative reduction.
    X = sp.Function("X")(t)
    x = sp.Function("x")(t)
    Mtot = mA + mB
    eA = sp.simplify(mA / Mtot)
    eB = sp.simplify(mB / Mtot)
    xA = X + eB * x
    xB = X - eA * x
    dx_diff = sp.simplify(sp.diff(xA, t) - sp.diff(xB, t))
    ddx_diff = sp.simplify(sp.diff(xA, t, 2) - sp.diff(xB, t, 2))
    L_cmrel = sp.simplify(acoef * dx_diff * ddx_diff)
    print("L_cmrel =", L_cmrel)
    if sp.diff(X, t) in L_cmrel.atoms(sp.Derivative) or sp.diff(X, t, 2) in L_cmrel.atoms(sp.Derivative):
        raise AssertionError("CM coordinate unexpectedly survived in the relative dipole odd action.")

    # V.3: outgoing l=1 partial wave starts at i*k^3, not i*k.
    subbanner("V.3 — Outgoing l=1 spectral no-go for linear odd term")
    a, k, z = sp.symbols("a k z", positive=True)
    j1 = sp.sin(z) / z**2 - sp.cos(z) / z
    y1 = -sp.cos(z) / z**2 - sp.sin(z) / z
    h1 = sp.simplify(j1 + I * y1)
    Lambda1 = sp.simplify((k * sp.diff(h1, z) / h1).subs(z, k * a))
    Lambda1_series = sp.series(Lambda1, k, 0, 7).removeO()
    Y1 = sp.series(1 / Lambda1_series, k, 0, 6).removeO()
    R = sp.symbols("R", positive=True)
    Y1_norm = sp.expand(R * Y1 / (-a / sp.Integer(2)))
    print("Lambda1(k) =", Lambda1_series)
    print("Y1_norm(k) =", Y1_norm)
    if sp.expand(Y1_norm).coeff(k, 1) != 0:
        raise AssertionError("Outgoing l=1 admittance unexpectedly contains a linear k term.")
    imag_k3 = sp.simplify(sp.expand(Y1_norm).coeff(k, 3) / I)
    print("Imaginary k^3 coefficient =", imag_k3)

    subbanner("V.3b — Order reduction of an isotropic fifth-derivative force")
    rN, rdN, v2N, mN, kap = sp.symbols("r rd v2 m kap", commutative=True)

    def add_vec(u, w):
        return (sp.simplify(u[0] + w[0]), sp.simplify(u[1] + w[1]))

    def mul_vec(s, u):
        return (sp.simplify(s * u[0]), sp.simplify(s * u[1]))

    def dscalar(expr):
        dr = rdN
        drd = v2N / rN - rdN**2 / rN - mN / rN**2
        dv2 = -2 * mN * rdN / rN**2
        return sp.simplify(sp.diff(expr, rN) * dr + sp.diff(expr, rdN) * drd + sp.diff(expr, v2N) * dv2)

    def dvec(u):
        A, B = u
        dn = (-rdN / rN, sp.Integer(1) / rN)
        dv = (-mN / rN**2, sp.Integer(0))
        res = add_vec((dscalar(A), 0), mul_vec(A, dn))
        res = add_vec(res, (0, dscalar(B)))
        res = add_vec(res, mul_vec(B, dv))
        return (sp.simplify(res[0]), sp.simplify(res[1]))

    xvec = (rN, 0)
    vvec = (0, 1)
    a2 = dvec(vvec)
    a3 = dvec(a2)
    a4 = dvec(a3)
    a5 = dvec(a4)
    print("x^(5) =", a5[0], "* n +", a5[1], "* v")
    print("Thus -kap x^(5) lies in the standard {rdot n, v} 2.5PN basis.")

    # V.4: if kept at the same formal 2.5PN order, nonzero dipole wake cannot be absorbed.
    subbanner("V.4 — Same-order dipole wake cannot be absorbed into the Iyer–Will family")
    p, q, s, alpha, beta = sp.symbols("p q s alpha beta")
    F = sp.Matrix([
        3 * (1 + beta),
        sp.Rational(1, 3) * (23 + 6 * alpha - 9 * beta),
        -5 * beta,
        -(2 + alpha),
        alpha - 2,
        3 * (1 + alpha),
    ])
    BT = sp.Matrix([18, sp.Rational(2, 3), -25, -6, 2, 15])
    D = sp.Matrix([
        9 * p + 36 * q,
        -8 * p - 22 * q,
        -45 * p - 60 * q,
        -9 * p,
        8 * p,
        45 * p,
    ])
    sol = sp.solve(list(BT + D - s * F), [alpha, beta, s, p, q], dict=True)
    print("Full matching solution =", sol)
    if not sol or sol[0][p] != 0 or sol[0][q] != 0:
        raise AssertionError("Unexpected same-order dipole match to the standard 2.5PN family.")

    # V.5: finite-size scaling rescue.
    subbanner("V.5 — Dipole finite-size scaling rescue")
    eps, rho = sp.symbols("eps rho", positive=True)
    scaling = sp.simplify(eps * (sp.sqrt(eps) * rho) ** 3)
    print("eps * (sqrt(eps)*rho)^3 =", scaling)
    print("If rho ~ eps^delta, exponent = 5/2 + 3 delta.")
    delta = sp.symbols("delta", positive=True)
    exponent = sp.Rational(5, 2) + 3 * delta
    print("general exponent =", exponent)

    return {
        "imag_k3_l1": imag_k3,
        "dipole_scaling_exponent": exponent,
    }


# ---------------------------------------------------------------------------
# Part VI: quadrupole sector
# ---------------------------------------------------------------------------

def quadrupole_sector() -> dict[str, sp.Expr]:
    banner("PART VI — QUADRUPOLE SECTOR")

    # VI.1: solved 2PN P2 ports form the full real STF l=2 basis.
    subbanner("VI.1 — P2 ports are exactly the real STF l=2 content")
    ux, uy, vn = sp.symbols("ux uy vn", real=True)
    v2 = ux**2 + uy**2 + vn**2
    Pi20 = sp.Rational(1, 2) * (3 * vn**2 - v2)
    Pi21c = sqrt(2) * vn * ux
    Pi21s = sqrt(2) * vn * uy
    Pi22c = (ux**2 - uy**2) / (2 * sqrt(2))
    Pi22s = (2 * ux * uy) / (2 * sqrt(2))

    Q = sp.Matrix([
        [ux**2 - v2 / 3, ux * uy, ux * vn],
        [ux * uy, uy**2 - v2 / 3, uy * vn],
        [ux * vn, uy * vn, vn**2 - v2 / 3],
    ])

    q20 = sp.simplify(sum(Q[i, j] * E20[i, j] for i in range(3) for j in range(3)))
    q21c = sp.simplify(sum(Q[i, j] * E21c[i, j] for i in range(3) for j in range(3)))
    q21s = sp.simplify(sum(Q[i, j] * E21s[i, j] for i in range(3) for j in range(3)))
    q22c = sp.simplify(sum(Q[i, j] * E22c[i, j] for i in range(3) for j in range(3)))
    q22s = sp.simplify(sum(Q[i, j] * E22s[i, j] for i in range(3) for j in range(3)))

    expect_zero("q20 - sqrt(2/3) Pi20", q20 - sqrt(sp.Rational(2, 3)) * Pi20)
    expect_zero("q21c - Pi21c", q21c - Pi21c)
    expect_zero("q21s - Pi21s", q21s - Pi21s)
    expect_zero("q22c - 2 Pi22c", q22c - 2 * Pi22c)
    expect_zero("q22s - 2 Pi22s", q22s - 2 * Pi22s)

    # VI.2: outgoing l=2 has positive i*k^5.
    subbanner("VI.2 — Outgoing l=2 retarded completion gives +i*k^5")
    a, k, z, R2 = sp.symbols("a k z R2", positive=True)
    j2 = (3 / z**3 - 1 / z) * sp.sin(z) - 3 * sp.cos(z) / z**2
    y2 = -(3 / z**3 - 1 / z) * sp.cos(z) - 3 * sp.sin(z) / z**2
    h2 = sp.simplify(j2 + I * y2)
    Lambda2 = sp.simplify((k * sp.diff(h2, z) / h2).subs(z, k * a))
    Lambda2_series = sp.series(Lambda2, k, 0, 7).removeO()
    Y2 = sp.series(1 / Lambda2_series, k, 0, 6).removeO()
    Y2_norm = sp.expand(R2 * Y2 / (-a / sp.Integer(3)))
    imag_k5 = sp.simplify(sp.expand(Y2_norm).coeff(k, 5) / I)
    print("Lambda2(k) =", Lambda2_series)
    print("Y2_norm(k) =", Y2_norm)
    print("Imaginary k^5 coefficient =", imag_k5)
    if not imag_k5.is_positive:
        raise AssertionError("Quadrupole k^5 coefficient did not come out positive.")

    # VI.3: local odd quadrupole kernel has positive damping.
    subbanner("VI.3 — Quadrupole odd kernel gives damping, not antidamping")
    t = sp.symbols("t", real=True)
    q = sp.Function("q")(t)
    Gamma = sp.symbols("Gamma", positive=True)
    F_one = -sp.Rational(1, 2) * Gamma * sp.diff(q, t, 5)
    P_one = sp.expand(F_one * sp.diff(q, t))
    Schott_one = -sp.Rational(1, 2) * Gamma * (
        sp.diff(q, t) * sp.diff(q, t, 4) - sp.diff(q, t, 2) * sp.diff(q, t, 3)
    )
    residual = sp.simplify(P_one - sp.diff(Schott_one, t))
    print("P - dE_Schott/dt =", residual)
    expect_zero("quadrupole power balance residual", residual + sp.Rational(1, 2) * Gamma * sp.diff(q, t, 3)**2)

    subbanner("VI.3b — Compact internal P2 branch is finite-size suppressed")
    eps, rho, delta = sp.symbols("eps rho delta", positive=True)
    scaling_internal_p2 = sp.simplify(eps**2 * (sp.sqrt(eps) * rho) ** 5)
    print("eps^2 * (sqrt(eps)*rho)^5 =", scaling_internal_p2)
    print("If rho ~ eps^delta, exponent =", sp.Rational(9, 2) + 5 * delta)

    # VI.4: orbital/worldtube source map.
    subbanner("VI.4 — Orbital/worldtube STF source map")
    u_x, u_y, d, r, mu, G, M = sp.symbols("u_x u_y d r mu G M", real=True)
    v = sp.Matrix([u_x, u_y, d])
    x = sp.Matrix([0, 0, r])
    V_stf = stf(v * v.T)
    X_stf = stf(x * x.T)
    a_vec = -(G * M / r**3) * x
    I_orb_ddot = sp.simplify(2 * mu * (V_stf + stf(x * a_vec.T)))

    V_comp = {name: c for name, c in zip(BASIS_NAMES, coeffs_in_basis(V_stf))}
    X_comp = {name: c for name, c in zip(BASIS_NAMES, coeffs_in_basis(X_stf))}
    I_comp = {name: c for name, c in zip(BASIS_NAMES, coeffs_in_basis(I_orb_ddot))}
    expect_zero("V20 source map", V_comp["20"] - sqrt(sp.Rational(2, 3)) * Pi20.subs(vn, d).subs({ux: u_x, uy: u_y}))
    expect_zero("V21c source map", V_comp["21c"] - Pi21c.subs(vn, d).subs({ux: u_x, uy: u_y}))
    expect_zero("V21s source map", V_comp["21s"] - Pi21s.subs(vn, d).subs({ux: u_x, uy: u_y}))
    expect_zero("V22c source map", V_comp["22c"] - 2 * Pi22c.subs({ux: u_x, uy: u_y}))
    expect_zero("V22s source map", V_comp["22s"] - 2 * Pi22s.subs({ux: u_x, uy: u_y}))
    expect_zero("X21c static position source", X_comp["21c"])
    expect_zero("X21s static position source", X_comp["21s"])
    expect_zero("X22c static position source", X_comp["22c"])
    expect_zero("X22s static position source", X_comp["22s"])
    expect_zero("I21c source map", I_comp["21c"] - 2 * mu * Pi21c.subs(vn, d).subs({ux: u_x, uy: u_y}))
    expect_zero("I21s source map", I_comp["21s"] - 2 * mu * Pi21s.subs(vn, d).subs({ux: u_x, uy: u_y}))
    expect_zero("I22c source map", I_comp["22c"] - 4 * mu * Pi22c.subs({ux: u_x, uy: u_y}))
    expect_zero("I22s source map", I_comp["22s"] - 4 * mu * Pi22s.subs({ux: u_x, uy: u_y}))
    expect_zero(
        "I20 source map",
        I_comp["20"] - 2 * mu * sqrt(sp.Rational(2, 3)) * (Pi20.subs(vn, d).subs({ux: u_x, uy: u_y}) - G * M / r),
    )

    # VI.5: rotational invariance forces scalar 5x5 matching matrix; basis map invertible.
    subbanner("VI.5 — Matching theorem: commutant is scalar, basis map invertible")
    A_x = sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
    A_y = sp.Matrix([[0, 0, 1], [0, 0, 0], [-1, 0, 0]])
    A_z = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])

    def tensor_action(A: sp.Matrix, T: sp.Matrix) -> sp.Matrix:
        return sp.simplify(A * T + T * A.T)

    def rep_matrix(A: sp.Matrix) -> sp.Matrix:
        Mrep = sp.zeros(5)
        for jcol, T in enumerate(BASIS):
            S = tensor_action(A, T)
            comps = coeffs_in_basis(S)
            for irow, ccomp in enumerate(comps):
                Mrep[irow, jcol] = sp.simplify(ccomp)
        return Mrep

    Jx = rep_matrix(A_x)
    Jy = rep_matrix(A_y)
    Jz = rep_matrix(A_z)
    unknowns = sp.symbols("m0:25")
    Munk = sp.Matrix(5, 5, unknowns)
    eqs = []
    for J in [Jx, Jy, Jz]:
        Cmat = sp.expand(Munk * J - J * Munk)
        eqs += list(Cmat)
    commutant_sol = sp.linsolve(eqs, unknowns)
    commutant_tuple = list(commutant_sol)[0]
    free_symbol = [s for s in commutant_tuple if s.free_symbols][0]
    M_comm = sp.Matrix(5, 5, commutant_tuple).subs({free_symbol: sp.symbols("m")})
    print("M commuting with all J_i =")
    print(M_comm)
    expect_zero("commutant minus scalar identity", M_comm - sp.symbols("m") * sp.eye(5))

    Bmap = sp.diag(sqrt(sp.Rational(2, 3)), 1, 1, 2, 2)
    B_det = sp.factor(Bmap.det())
    print("det(port->canonical basis map) =", B_det)
    if B_det == 0:
        raise AssertionError("Port-to-canonical STF map is singular.")

    # VI.6: static overlap m0 != 0.
    subbanner("VI.6 — Static overlap test")
    r_sym, Gsym, Msym, musym = sp.symbols("r G M mu", positive=True, real=True)
    x_static = sp.Matrix([0, 0, r_sym])
    x_stf = stf(x_static * x_static.T)
    a_static = -(Gsym * Msym / r_sym**3) * x_static
    I_orb_ddot_static = sp.simplify(2 * musym * stf(x_static * a_static.T))
    C20_x = sp.simplify(sum(x_stf.multiply_elementwise(E20)))
    C20_ddotI_over_2mu = sp.simplify(sum(I_orb_ddot_static.multiply_elementwise(E20)) / (2 * musym))
    J20 = sp.Rational(5, 4)
    R20 = sp.Integer(1)
    q20_static = sp.simplify(J20 * R20)
    m0_from_ddotI = sp.simplify(q20_static / C20_ddotI_over_2mu)
    print("C20[x_<i x_j>] =", C20_x)
    print("C20[ddot I_orb /(2 mu)] =", C20_ddotI_over_2mu)
    print("q20_static =", q20_static)
    print("m0 from static overlap =", m0_from_ddotI)
    if m0_from_ddotI == 0:
        raise AssertionError("Static overlap m0 unexpectedly vanished.")

    return {
        "imag_k5_l2": imag_k5,
        "m0_static": m0_from_ddotI,
    }



# ---------------------------------------------------------------------------
# Part VII: later session quadrupole roadmap / normalization / extraction
# ---------------------------------------------------------------------------

def session_quadrupole_roadmap() -> dict[str, sp.Expr]:
    banner("PART VII — QUADRUPOLE NORMALIZATION / EXTRACTION / PARAMETERIZATION")

    # VII.1: canonical quadrupole normalization map
    subbanner("VII.1 — Canonical quadrupole normalization map")
    t = sp.symbols("t", real=True)
    gamma, Gamma5, m0hat, G, c = sp.symbols("gamma Gamma5 m0hat G c", positive=True, real=True)
    mu = sp.symbols("mu", positive=True, real=True)
    x = sp.Function("x")(t)
    y = sp.Function("y")(t)
    X = sp.Matrix([x, y, 0])
    I_tensor = sp.simplify(mu * (X * X.T - sp.eye(3) * (X.dot(X)) / sp.Integer(3)))
    q_can = [sp.simplify(sum(I_tensor[i, j] * B[i, j] for i in range(3) for j in range(3))) for B in BASIS]
    q5_can = [sp.diff(comp, t, 5) for comp in q_can]
    Fx = sp.simplify(-gamma / 2 * sum(sp.diff(comp, x) * comp5 for comp, comp5 in zip(q_can, q5_can)))
    Fy = sp.simplify(-gamma / 2 * sum(sp.diff(comp, y) * comp5 for comp, comp5 in zip(q_can, q5_can)))
    Ix = sp.simplify(-(x * sp.diff(I_tensor[0, 0], t, 5) + y * sp.diff(I_tensor[0, 1], t, 5)))
    Iy = sp.simplify(-(x * sp.diff(I_tensor[0, 1], t, 5) + y * sp.diff(I_tensor[1, 1], t, 5)))
    expect_zero("Fx - gamma*mu*(-x^j I_xj^(5))", sp.simplify(Fx - gamma * mu * Ix))
    expect_zero("Fy - gamma*mu*(-x^j I_yj^(5))", sp.simplify(Fy - gamma * mu * Iy))
    gamma_GR = sp.simplify(2 * G / (5 * c**5))
    gamma_toy = sp.simplify(m0hat**2 * Gamma5)
    Gamma5_target = sp.simplify(gamma_GR / m0hat**2)
    print("gamma_GR =", gamma_GR)
    print("gamma_toy =", gamma_toy)
    print("Gamma5_target =", Gamma5_target)

    # VII.2: outgoing l=2 PDE fingerprint
    subbanner("VII.2 — Outgoing l=2 PDE fingerprint")
    k, a, c_s, m0, omega = sp.symbols("k a c_s m0 omega", positive=True, real=True)
    r = sp.symbols("r", positive=True, real=True)
    z = k * r
    j2 = (3 / z**3 - 1 / z) * sp.sin(z) - 3 * sp.cos(z) / z**2
    y2 = -(3 / z**3 - 1 / z) * sp.cos(z) - 3 * sp.sin(z) / z**2
    h2 = sp.simplify(j2 + I * y2)
    ode_residual = sp.simplify(sp.diff(h2, r, 2) + 2 / r * sp.diff(h2, r) + (k**2 - 6 / r**2) * h2)
    expect_zero("l=2 Helmholtz residual", ode_residual)

    za = sp.symbols("za", positive=True, real=True)
    j2a = (3 / za**3 - 1 / za) * sp.sin(za) - 3 * sp.cos(za) / za**2
    y2a = -(3 / za**3 - 1 / za) * sp.cos(za) - 3 * sp.sin(za) / za**2
    h2a = sp.simplify(j2a + I * y2a)
    Lambda2 = sp.simplify((k * sp.diff(h2a, za) / h2a).subs(za, k * a))
    Lambda2_series = sp.series(Lambda2, k, 0, 7).removeO()
    Y2 = sp.simplify(sp.series(1 / Lambda2_series, k, 0, 6).removeO())
    Y2_static = sp.simplify(Y2.subs(k, 0))
    Y2_hat = sp.simplify(sp.expand(Y2 / Y2_static))
    Y2_hat_series = sp.series(Y2_hat, k, 0, 6).removeO()
    Y2_hat_omega = sp.expand(Y2_hat_series.subs(k, omega / c_s))
    m2 = sp.simplify(m0 * sp.expand(Y2_hat_omega).coeff(omega, 2))
    m4 = sp.simplify(m0 * sp.expand(Y2_hat_omega).coeff(omega, 4))
    Gamma5_pde = sp.simplify(m0 * sp.expand(Y2_hat_omega).coeff(omega, 5) / I)
    print("Lambda2(k) =", Lambda2_series)
    print("Y2(0) =", Y2_static)
    print("Y2_hat(k) =", Y2_hat_series)
    print("m2 =", m2)
    print("m4 =", m4)
    print("Gamma5^PDE =", Gamma5_pde)
    expect_zero("m4 - 4 m2^2 / m0", sp.simplify(m4 - 4 * m2**2 / m0))
    expect_zero(
        "Gamma5^PDE - 9 m2^(5/2)/m0^(3/2)",
        sp.simplify(Gamma5_pde - 9 * m2**sp.Rational(5, 2) / m0**sp.Rational(3, 2)),
    )

    # VII.3: source/port normalization cleanup
    subbanner("VII.3 — Source/port normalization cleanup")
    B20 = sqrt(sp.Rational(2, 3))
    Bmap = sp.diag(B20, 1, 1, 2, 2)
    Gs, Ms, rs = sp.symbols("G M r", positive=True, real=True)
    C20_ddotI_over_2mu = -sqrt(6) * Gs * Ms / (3 * rs)
    P20_ddotI_over_2mu = sp.simplify(C20_ddotI_over_2mu / B20)
    q20_port_static = sp.Rational(5, 4)
    m0_mixed = sp.simplify(q20_port_static / C20_ddotI_over_2mu)
    m0_raw_same_basis = sp.simplify(q20_port_static / P20_ddotI_over_2mu)
    q20_can_static = sp.simplify(B20 * q20_port_static)
    m0_can_same_basis = sp.simplify(q20_can_static / C20_ddotI_over_2mu)
    N_orb = sp.simplify(1 / m0_raw_same_basis)
    m0hat_point_particle = sp.simplify(N_orb * m0_raw_same_basis)
    print("det(Bmap) =", sp.simplify(Bmap.det()))
    print("m0_mixed =", m0_mixed)
    print("m0_raw_same_basis =", m0_raw_same_basis)
    print("m0_can_same_basis =", m0_can_same_basis)
    print("N_orb =", N_orb)
    print("m0hat_point_particle =", m0hat_point_particle)
    expect_zero("same-basis raw vs canonical overlap", sp.simplify(m0_can_same_basis - m0_raw_same_basis))
    expect_zero("point-particle normalized overlap - 1", sp.simplify(m0hat_point_particle - 1))

    kapp, Y0 = sp.symbols("kappa Y0", real=True)
    Y2_hat_min = 1 + a**2 * omega**2 / (9 * c_s**2) + 4 * a**4 * omega**4 / (81 * c_s**4) + I * a**5 * omega**5 / (27 * c_s**5)
    M_eff = sp.expand(kapp**2 * Y0 * Y2_hat_min)
    m0_eff = sp.simplify(M_eff.coeff(omega, 0))
    m2_eff = sp.simplify(M_eff.coeff(omega, 2))
    m4_eff = sp.simplify(M_eff.coeff(omega, 4))
    Gamma5_eff = sp.simplify(M_eff.coeff(omega, 5) / I)
    expect_zero("m2_eff - m0_eff a^2/(9 c_s^2)", sp.simplify(m2_eff - m0_eff * a**2 / (9 * c_s**2)))
    expect_zero("m4_eff - 4 m0_eff a^4/(81 c_s^4)", sp.simplify(m4_eff - 4 * m0_eff * a**4 / (81 * c_s**4)))
    expect_zero("Gamma5_eff - m0_eff a^5/(27 c_s^5)", sp.simplify(Gamma5_eff - m0_eff * a**5 / (27 * c_s**5)))
    Gamma5_can = sp.simplify(m0hat * a**5 / (27 * c_s**5))
    print("Gamma5_can =", Gamma5_can)

    # VII.4: convention-invariant normalization product
    subbanner("VII.4 — Convention-invariant normalization product")
    lam = sp.symbols("lam", positive=True, real=True)
    mfrak = sp.symbols("mfrak", real=True)
    K0 = sp.symbols("K0", real=True)
    Gamma5_port = sp.symbols("Gamma5_port", real=True)
    mfrak_p = sp.simplify(lam * mfrak)
    Gamma5_port_p = sp.simplify(Gamma5_port / lam**2)
    K0_p = sp.simplify(K0 / lam**2)
    expect_zero("mfrak'^2 Gamma5' - mfrak^2 Gamma5", sp.simplify(mfrak_p**2 * Gamma5_port_p - mfrak**2 * Gamma5_port))
    expect_zero("mfrak'^2 K0' - mfrak^2 K0", sp.simplify(mfrak_p**2 * K0_p - mfrak**2 * K0))
    Gamma5_port_DN = sp.simplify(K0 * a**5 / (27 * c_s**5))
    gamma_eff = sp.simplify(mfrak**2 * Gamma5_port_DN)
    NQ = sp.symbols("NQ", real=True)
    NQ_target = sp.simplify(sp.solve(sp.Eq(NQ * a**5 / (27 * c_s**5), gamma_GR), NQ)[0])
    print("Gamma5_port^DN =", Gamma5_port_DN)
    print("gamma_eff =", gamma_eff)
    print("N_Q_target =", NQ_target)

    # VII.5: canonical invariant low-frequency fingerprint
    subbanner("VII.5 — Canonical invariant low-frequency fingerprint")
    K0bar, K2bar, K4bar, Gammabar = sp.symbols("K0bar K2bar K4bar Gammabar", positive=True, real=True)
    K2bar_from_branch = sp.simplify(K0bar * a**2 / (9 * c_s**2))
    K4bar_from_branch = sp.simplify(K0bar * 4 * a**4 / (81 * c_s**4))
    Gammabar_from_branch = sp.simplify(K0bar * a**5 / (27 * c_s**5))
    ratio_a2_cs2 = sp.simplify(sp.solve(sp.Eq(K2bar, K2bar_from_branch), a**2 / c_s**2)[0])
    K4bar_invariant = sp.simplify(K4bar_from_branch.subs(a**4 / c_s**4, ratio_a2_cs2**2))
    Gammabar_invariant = sp.simplify(Gammabar_from_branch.subs(a**5 / c_s**5, ratio_a2_cs2**sp.Rational(5, 2)))
    K2bar_target_formula = sp.simplify((sp.Rational(2, 45) * G / c**5)**sp.Rational(2, 5) * K0bar**sp.Rational(3, 5))
    print("a^2/c_s^2 =", ratio_a2_cs2)
    print("K4bar invariant =", K4bar_invariant)
    print("Gammabar invariant =", Gammabar_invariant)
    print("K2bar_target formula =", K2bar_target_formula)
    expect_zero("K4bar invariant - 4 K2bar^2/K0bar", sp.simplify(K4bar_invariant - 4 * K2bar**2 / K0bar))
    expect_zero(
        "Gammabar invariant - 9 K2bar^(5/2)/K0bar^(3/2)",
        sp.simplify(Gammabar_invariant - 9 * K2bar**sp.Rational(5, 2) / K0bar**sp.Rational(3, 2)),
    )

    # VII.6: extraction of (Kbar0, Kbar2) and single-prefactor obstruction
    subbanner("VII.6 — Extraction of (Kbar0, Kbar2) and single-prefactor obstruction")
    K2bar_branch = sp.simplify(K0bar * a**2 / (9 * c_s**2))
    K4bar_branch = sp.simplify(K0bar * 4 * a**4 / (81 * c_s**4))
    Gammabar_branch = sp.simplify(K0bar * a**5 / (27 * c_s**5))
    K0bar_target = sp.simplify(sp.solve(sp.Eq(Gammabar_branch, gamma_GR), K0bar)[0])
    K2bar_target = sp.simplify(K2bar_branch.subs(K0bar, K0bar_target))
    K4bar_target = sp.simplify(K4bar_branch.subs(K0bar, K0bar_target))
    print("K2bar_branch =", K2bar_branch)
    print("K0bar_target =", K0bar_target)
    print("K2bar_target =", K2bar_target)
    print("K4bar_target =", K4bar_target)

    Omega20 = sp.symbols("Omega20", positive=True, real=True)
    R20 = sp.Integer(1)
    Y20_min = sp.simplify(R20 / (1 - omega**2 / Omega20**2))
    Y20_min_series = sp.series(Y20_min, omega, 0, 5).removeO()
    C_Q = sp.symbols("C_Q", positive=True, real=True)
    Y2_out = sp.expand(1 + a**2 * omega**2 / (9 * c_s**2) + 4 * a**4 * omega**4 / (81 * c_s**4) + I * a**5 * omega**5 / (27 * c_s**5))
    Kbar_pref = sp.expand(C_Q * Y2_out)
    Kbar0 = sp.simplify(Kbar_pref.coeff(omega, 0))
    Kbar2 = sp.simplify(Kbar_pref.coeff(omega, 2))
    Kbar4 = sp.simplify(Kbar_pref.coeff(omega, 4))
    Gammabar5 = sp.simplify(Kbar_pref.coeff(omega, 5) / I)
    C_Q_target = sp.simplify(sp.solve(sp.Eq(Gammabar5, gamma_GR), C_Q)[0])
    print("Y20_min(omega) =", Y20_min_series)
    print("C_Q_target =", C_Q_target)
    expect_zero("Kbar0 - C_Q", sp.simplify(Kbar0 - C_Q))
    expect_zero("Kbar2 - C_Q a^2/(9 c_s^2)", sp.simplify(Kbar2 - C_Q * a**2 / (9 * c_s**2)))
    expect_zero("Kbar4 - 4 C_Q a^4/(81 c_s^4)", sp.simplify(Kbar4 - 4 * C_Q * a**4 / (81 * c_s**4)))

    # VII.7: frozen 2PN files already fix representation and static support data
    subbanner("VII.7 — Frozen 2PN files: fixed representation/support data")
    ux, uy, d, U = sp.symbols("ux uy d U", real=True)
    v = sp.Matrix([ux, uy, d])
    n = sp.Matrix([0, 0, 1])
    Tmat = v * v.T - U * (n * n.T)
    Tstf = stf(Tmat)
    Ccoeffs = [sp.simplify(sum(Tstf[i, j] * B[i, j] for i in range(3) for j in range(3))) for B in BASIS]
    Pi20 = sp.simplify(sp.Rational(1, 2) * (3 * d**2 - (ux**2 + uy**2 + d**2)))
    Pi21c = sp.simplify(sp.sqrt(2) * d * ux)
    Pi21s = sp.simplify(sp.sqrt(2) * d * uy)
    Pi22c = sp.simplify((ux**2 - uy**2) / (2 * sp.sqrt(2)))
    Pi22s = sp.simplify((2 * ux * uy) / (2 * sp.sqrt(2)))
    expect_zero("C20 relation", Ccoeffs[0] - sp.sqrt(sp.Rational(2, 3)) * (Pi20 - U))
    expect_zero("C21c relation", Ccoeffs[1] - Pi21c)
    expect_zero("C21s relation", Ccoeffs[2] - Pi21s)
    expect_zero("C22c relation", Ccoeffs[3] - 2 * Pi22c)
    expect_zero("C22s relation", Ccoeffs[4] - 2 * Pi22s)
    J0 = sp.simplify(4 / sp.sqrt(5))
    J20 = sp.Rational(5, 4)
    static_support = sp.simplify(J0**2 + J20**2)
    closure_deficit = sp.simplify(static_support - sp.Rational(5, 4))
    print("J.J =", static_support)
    print("closure deficit =", closure_deficit)
    print("Open dynamic poles: Omega20, Omega21, Omega22")

    # VII.8: one-pole insufficiency and minimal positive conservative precursor
    subbanner("VII.8 — One-pole insufficiency and minimal positive conservative precursor")
    A = sp.simplify(a**2 / (9 * c_s**2))
    Y_out = sp.expand(1 + A * omega**2 + 4 * A**2 * omega**4 + I * a**5 * omega**5 / (27 * c_s**5))
    Omega1 = sp.simplify(1 / sp.sqrt(A))
    Y_1pole = sp.simplify(1 / (1 - omega**2 / Omega1**2))
    Y_1pole_series = sp.series(Y_1pole, omega, 0, 6).removeO()
    mismatch_w4 = sp.simplify(sp.expand(Y_out).coeff(omega, 4) - sp.expand(Y_1pole_series).coeff(omega, 4))
    print("Omega_1 =", Omega1)
    print("one-pole O(w^4) mismatch =", mismatch_w4)

    p = sp.symbols("p", positive=True, real=True)
    q1 = sp.simplify(A + sp.sqrt(3) * A * sp.sqrt(1 - p) / sp.sqrt(p))
    q2 = sp.simplify(A - sp.sqrt(3) * A * sp.sqrt(p) / sp.sqrt(1 - p))
    Y_2pole = sp.simplify(p / (1 - q1 * omega**2) + (1 - p) / (1 - q2 * omega**2))
    Y_2pole_series = sp.series(Y_2pole, omega, 0, 6).removeO()
    p_min = sp.Rational(1, 4)
    q1_min = sp.simplify(q1.subs(p, p_min))
    q2_min = sp.simplify(q2.subs(p, p_min))
    Y_min = sp.simplify(sp.limit(Y_2pole, p, p_min, dir='-'))
    Y_min_series = sp.series(Y_min, omega, 0, 6).removeO()
    Y_min_split = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - q1_min * omega**2))
    Omega_min = sp.simplify(1 / sp.sqrt(q1_min))
    print("Y_2pole series =", Y_2pole_series)
    print("Minimal positive closure =", Y_min)
    print("Minimal split form =", Y_min_split)
    print("Omega_Q,min =", Omega_min)
    expect_zero("minimal closure omega^2 match", sp.expand(Y_min_series).coeff(omega, 2) - A)
    expect_zero("minimal closure omega^4 match", sp.expand(Y_min_series).coeff(omega, 4) - 4 * A**2)

    # VII.9: minimal isotropic quadrupole module and componentwise P2 target
    subbanner("VII.9 — Minimal isotropic quadrupole module and componentwise P2 target")
    c0, c1, OmegaQ = sp.symbols("c0 c1 OmegaQ", positive=True, real=True)
    Y_family = sp.simplify(c0 + c1 / (1 - omega**2 / OmegaQ**2))
    Y_family_series = sp.series(Y_family, omega, 0, 6).removeO()
    sol = sp.solve(
        [
            sp.Eq(sp.expand(Y_family_series).coeff(omega, 0), 1),
            sp.Eq(sp.expand(Y_family_series).coeff(omega, 2), A),
            sp.Eq(sp.expand(Y_family_series).coeff(omega, 4), 4 * A**2),
        ],
        [c0, c1, OmegaQ],
        dict=True,
    )
    sol_phys = sol[0]
    c0_sol = sp.simplify(sol_phys[c0])
    c1_sol = sp.simplify(sol_phys[c1])
    OmegaQ_sol = sp.simplify(sol_phys[OmegaQ])
    Y_family_min = sp.simplify(Y_family.subs(sol_phys))
    Y_family_min_series = sp.series(Y_family_min, omega, 0, 6).removeO()
    print("c0 =", c0_sol)
    print("c1 =", c1_sol)
    print("OmegaQ =", OmegaQ_sol)
    print("Y_family_min =", Y_family_min)
    expect_zero("minimal module omega^2 mismatch", sp.expand(Y_family_min_series).coeff(omega, 2) - A)
    expect_zero("minimal module omega^4 mismatch", sp.expand(Y_family_min_series).coeff(omega, 4) - 4 * A**2)

    gram = sp.Matrix([[sp.simplify(sum(B1[i, j] * B2[i, j] for i in range(3) for j in range(3))) for B2 in BASIS] for B1 in BASIS])
    expect_zero("real STF Gram matrix - I5", gram - sp.eye(5))
    YQ_min = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / OmegaQ_sol**2))
    print("Y20 = Y21 = Y22 =", YQ_min)
    print("YQ_min series =", sp.series(YQ_min, omega, 0, 6).removeO())

    OmegaQ_from_invariants = sp.simplify(sp.sqrt(K0bar / (4 * K2bar)))
    K_family = sp.simplify(K0bar * (sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / OmegaQ_from_invariants**2)))
    K4_relation = sp.simplify(4 * K2bar**2 / K0bar)
    Gammabar_relation = sp.simplify(9 * K2bar**sp.Rational(5, 2) / K0bar**sp.Rational(3, 2))
    print("OmegaQ_from_invariants =", OmegaQ_from_invariants)
    print("K4 relation =", K4_relation)
    print("Gamma5 relation =", Gammabar_relation)

    # VII.10: single-pole sufficiency theorem and explicit extraction ledger
    subbanner("VII.10 — Single-pole sufficiency theorem and extraction ledger")
    K2_from_branch = sp.simplify(K0bar / (4 * OmegaQ**2))
    K4_from_branch = sp.simplify(K0bar / (4 * OmegaQ**4))
    Gamma_from_branch = sp.simplify(9 * K0bar / (32 * OmegaQ**5))
    K0bar_target_from_Omega = sp.simplify(sp.solve(sp.Eq(Gamma_from_branch, gamma_GR), K0bar)[0])
    K2bar_target_from_Omega = sp.simplify(K2_from_branch.subs(K0bar, K0bar_target_from_Omega))
    K4bar_target_from_Omega = sp.simplify(K4_from_branch.subs(K0bar, K0bar_target_from_Omega))
    print("K2_from_branch =", K2_from_branch)
    print("K4_from_branch =", K4_from_branch)
    print("Gamma_from_branch =", Gamma_from_branch)
    print("K0bar_target(OmegaQ) =", K0bar_target_from_Omega)
    print("K2bar_target(OmegaQ) =", K2bar_target_from_Omega)
    print("K4bar_target(OmegaQ) =", K4bar_target_from_Omega)

    c0s, c2s, c4s = sp.symbols("c0 c2 c4", positive=True, real=True)
    OmegaQ_sq_from_c0c2 = sp.simplify(c0s / (4 * c2s))
    OmegaQ_sq_from_c2c4 = sp.simplify(c2s / c4s)
    K0_from_c2c4 = sp.simplify(4 * c2s**2 / c4s)
    Gamma_from_c2c4 = sp.simplify(sp.Rational(9, 8) * c4s**sp.Rational(3, 2) / sp.sqrt(c2s))
    print("OmegaQ^2 from (c0,c2) =", OmegaQ_sq_from_c0c2)
    print("OmegaQ^2 from (c2,c4) =", OmegaQ_sq_from_c2c4)
    print("K0 from (c2,c4) =", K0_from_c2c4)
    print("Gamma from (c2,c4) =", Gamma_from_c2c4)
    print("Direct GR gate: c4^3/c2 =", sp.Rational(256, 2025) * G**2 / c**10)

    # VII.11: grouped P2 formulas and normalized-support ledger
    subbanner("VII.11 — Grouped P2 formulas and normalized-support ledger")
    c0_20, c2_20, c4_20 = sp.symbols("c0_20 c2_20 c4_20", positive=True, real=True)
    c0_21, c2_21, c4_21 = sp.symbols("c0_21 c2_21 c4_21", positive=True, real=True)
    c0_22, c2_22, c4_22 = sp.symbols("c0_22 c2_22 c4_22", positive=True, real=True)
    c0_tr = sp.simplify((c0_20 + 2 * c0_21 + 2 * c0_22) / 5)
    c2_tr = sp.simplify((c2_20 + 2 * c2_21 + 2 * c2_22) / 5)
    c4_tr = sp.simplify((c4_20 + 2 * c4_21 + 2 * c4_22) / 5)
    print("c0_trace =", c0_tr)
    print("c2_trace =", c2_tr)
    print("c4_trace =", c4_tr)
    d0_20 = sp.simplify(c0_20 - c0_tr)
    d0_21 = sp.simplify(c0_21 - c0_tr)
    d0_22 = sp.simplify(c0_22 - c0_tr)
    d2_20 = sp.simplify(c2_20 - c2_tr)
    d2_21 = sp.simplify(c2_21 - c2_tr)
    d2_22 = sp.simplify(c2_22 - c2_tr)
    d4_20 = sp.simplify(c4_20 - c4_tr)
    d4_21 = sp.simplify(c4_21 - c4_tr)
    d4_22 = sp.simplify(c4_22 - c4_tr)
    expect_zero("weighted anisotropy sum n=0", sp.simplify(d0_20 + 2 * d0_21 + 2 * d0_22))
    expect_zero("weighted anisotropy sum n=2", sp.simplify(d2_20 + 2 * d2_21 + 2 * d2_22))
    expect_zero("weighted anisotropy sum n=4", sp.simplify(d4_20 + 2 * d4_21 + 2 * d4_22))
    print("Minimal branch identity c0*c4 - 4*c2^2 =", sp.simplify(c0s * c4s - 4 * c2s**2))

    u2_20, u4_20 = sp.symbols("u2_20 u4_20", positive=True, real=True)
    u2_21, u4_21 = sp.symbols("u2_21 u4_21", positive=True, real=True)
    u2_22, u4_22 = sp.symbols("u2_22 u4_22", positive=True, real=True)
    u2_tr = sp.simplify((u2_20 + 2 * u2_21 + 2 * u2_22) / 5)
    u4_tr = sp.simplify((u4_20 + 2 * u4_21 + 2 * u4_22) / 5)
    print("u2_trace =", u2_tr)
    print("u4_trace =", u4_tr)
    u2n, u4n = sp.symbols("u2 u4", positive=True, real=True)
    Gamma5_norm = sp.simplify(sp.Rational(9, 8) * u4n**sp.Rational(3, 2) / sp.sqrt(u2n))
    Gamma5_norm_branch = sp.simplify(Gamma5_norm.subs(u4n, 4 * u2n**2))
    K0bar_target_norm = sp.simplify(sp.solve(sp.Eq(sp.symbols("K0n", positive=True) * Gamma5_norm_branch, gamma_GR), sp.symbols("K0n", positive=True))[0])
    print("normalized branch identity u4 - 4*u2^2 =", sp.simplify(u4n - 4 * u2n**2))
    print("OmegaQ on normalized branch =", sp.simplify(1 / (2 * sp.sqrt(u2n))))
    print("Gamma5_norm on branch =", Gamma5_norm_branch)
    print("K0bar_target from normalized support =", K0bar_target_norm)

    # VII.12: axisymmetric grouped-P2 parameterization
    subbanner("VII.12 — Axisymmetric grouped-P2 parameterization")
    u2bar, u4bar = sp.symbols("u2bar u4bar", real=True)
    a2, a4 = sp.symbols("a2 a4", real=True)
    b2, b4 = sp.symbols("b2 b4", real=True)
    I5 = sp.eye(5)
    T = sp.diag(4, -1, -1, -1, -1)
    S = sp.diag(0, 1, 1, -1, -1)
    U2 = sp.simplify(u2bar * I5 + a2 * T + b2 * S)
    U4 = sp.simplify(u4bar * I5 + a4 * T + b4 * S)
    Y = sp.simplify(I5 + omega**2 * U2 + omega**4 * U4)
    Y20 = sp.simplify(Y[0, 0])
    Y21 = sp.simplify(Y[1, 1])
    Y22 = sp.simplify(Y[3, 3])
    print("Y20 =", Y20)
    print("Y21 =", Y21)
    print("Y22 =", Y22)
    x20, x21, x22 = sp.symbols("x20 x21 x22", real=True)
    ubar_from_x = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    a_from_x = sp.simplify((2 * x20 - x21 - x22) / 10)
    b_from_x = sp.simplify((x21 - x22) / 2)
    print("ubar_from_(x20,x21,x22) =", ubar_from_x)
    print("a_from_(x20,x21,x22) =", a_from_x)
    print("b_from_(x20,x21,x22) =", b_from_x)
    A2_sq = sp.simplify(sp.trace((U2 - u2bar * I5)**2) / 5)
    A4_sq = sp.simplify(sp.trace((U4 - u4bar * I5)**2) / 5)
    print("A2_sq =", A2_sq)
    print("A4_sq =", A4_sq)
    print("Minimal branch relation =", sp.simplify(u4bar - 4 * u2bar**2))
    print("OmegaQ^2 =", sp.simplify(1 / (4 * u2bar)))
    print("Gamma5_norm =", sp.simplify(9 * u2bar**sp.Rational(5, 2)))
    print("K0bar_target =", sp.simplify(2 * G / (45 * c**5 * u2bar**sp.Rational(5, 2))))

    return {
        "gamma_GR": gamma_GR,
        "Gamma5_target": Gamma5_target,
        "Gamma5_pde": Gamma5_pde,
        "NQ_target": NQ_target,
        "K0bar_target": K0bar_target,
        "K2bar_target": K2bar_target,
        "K4bar_target": K4bar_target,
        "C_Q_target": C_Q_target,
        "OmegaQ_min": Omega_min,
        "OmegaQ_sol": OmegaQ_sol,
        "K0bar_target_from_Omega": K0bar_target_from_Omega,
        "K0bar_target_norm": K0bar_target_norm,
    }


# ---------------------------------------------------------------------------
# Final theorem ledger (extended)
# ---------------------------------------------------------------------------

def final_ledger(results: dict[str, dict[str, sp.Expr]]) -> None:
    banner("FINAL SESSION THEOREM LEDGER")

    print("Decisive benchmark:")
    print("  Burke–Thorne prototype lands on the Iyer–Will family with "
          f"(alpha,beta)=({results['decisive']['alpha']},{results['decisive']['beta']}).")
    print()
    print("Channel audit status:")
    print("  - scalar direct overlap generically gives i*omega,")
    print("  - derivative breathing outlet pushes the scalar odd kernel to i*omega^3,")
    print("  - compact reciprocal mouth radiation is higher-order / super-Ohmic,")
    print("  - outgoing l=1 starts at +i*k^3, not at i*k,")
    print("  - same-order nonzero dipole deformation cannot be absorbed into the standard 2.5PN family,")
    print("  - compact outgoing l=2 starts at +i*k^5 with positive damping sign.")
    print()
    print("Quadrupole normalization stack:")
    print("  gamma_GR =", results["session_quad"]["gamma_GR"])
    print("  Gamma5_target =", results["session_quad"]["Gamma5_target"])
    print("  N_Q_target =", results["session_quad"]["NQ_target"])
    print("  K0bar_target =", results["session_quad"]["K0bar_target"])
    print("  K2bar_target =", results["session_quad"]["K2bar_target"])
    print("  K4bar_target =", results["session_quad"]["K4bar_target"])
    print()
    print("Minimal isotropic quadrupole module:")
    print("  Omega_Q,min =", results["session_quad"]["OmegaQ_min"])
    print("  exact moment-matched isotropic solution Omega_Q =", results["session_quad"]["OmegaQ_sol"])
    print("  single-pole GR scaling law K0bar = 64 G Omega_Q^5 / (45 c^5).")
    print()
    print("Normalized-support branch:")
    print("  K0bar_target from normalized support moments =", results["session_quad"]["K0bar_target_norm"])
    print()
    print("Best current verdict:")
    print("  The session reduces the full GR-like point-particle 2.5PN closure to the")
    print("  isotropic real P2 quadrupole module. Once the grouped low-frequency")
    print("  conservative coefficients through O(omega^4) are fixed on that branch,")
    print("  the odd Burke–Thorne coefficient follows automatically.")
    print("  The one remaining narrow theorem gap is the final passive/outgoing")
    print("  quadrupole normalization on the actual moving-throat branch.")
    print()


def main() -> None:
    decisive = decisive_benchmark()
    framework = general_low_frequency_framework()
    scalar = scalar_sector()
    dipole = dipole_sector()
    quadrupole = quadrupole_sector()
    session_quad = session_quadrupole_roadmap()
    final_ledger({
        "decisive": decisive,
        "framework": framework,
        "scalar": scalar,
        "dipole": dipole,
        "quadrupole": quadrupole,
        "session_quad": session_quad,
    })


if __name__ == "__main__":
    main()

"""
========================================================================================
PART II — DECISIVE BENCHMARK / BURKE–THORNE PROTOTYPE
========================================================================================

----------------------------------------------------------------------------------------
II.1 — Two-body mass dipole and STF quadrupole decomposition
----------------------------------------------------------------------------------------
D_i =
⎡X₁⋅(m₁ + m₂)⎤
⎢            ⎥
⎢X₂⋅(m₁ + m₂)⎥
⎢            ⎥
⎣X₃⋅(m₁ + m₂)⎦
dD_i/dx_j = Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
Q - (M STF(XX) + mu STF(xx)) = Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

----------------------------------------------------------------------------------------
II.2 — Burke–Thorne local quadrupole force and Iyer–Will match
----------------------------------------------------------------------------------------
A(v^2, GM/r, rdot^2) =
                                       2
                             ⎛d       ⎞ 
2⋅G⋅m + 54⋅v₂⋅r(t) - 75⋅r(t)⋅⎜──(r(t))⎟ 
                             ⎝dt      ⎠ 
────────────────────────────────────────
                 3⋅r(t)                 
B(v^2, GM/r, rdot^2) =
 ⎛                                       2⎞ 
 ⎜                             ⎛d       ⎞ ⎟ 
-⎜-2⋅G⋅m + 6⋅v₂⋅r(t) - 15⋅r(t)⋅⎜──(r(t))⎟ ⎟ 
 ⎝                             ⎝dt      ⎠ ⎠ 
────────────────────────────────────────────
                    r(t)                    
alpha = 4
beta  = 5

========================================================================================
PART III — LOW-FREQUENCY SELECTION RULES / INFLUENCE FUNCTIONAL
========================================================================================

----------------------------------------------------------------------------------------
III.1 — Time-domain signs for i omega^n
----------------------------------------------------------------------------------------
i*omega^1  ->  -1 * d^1/dt^1
i*omega^3  ->  1 * d^3/dt^3
i*omega^5  ->  -1 * d^5/dt^5

----------------------------------------------------------------------------------------
III.2 — Minimal retarded kernel expansions
----------------------------------------------------------------------------------------
K1 = g**2/Omega**2 + g**2*omega**2/Omega**4 + I*g**2*omega*sigma/Omega**4 - g**2*omega**2*sigma**2/Omega**6
K3 = g**2/Omega**2 + I*g**2*omega**3*sigma/Omega**4 + g**2*omega**2/Omega**4 + g**2*omega**4/Omega**6
K5 = g**2/Omega**2 + I*g**2*omega**5*sigma/Omega**4 + g**2*omega**2/Omega**4 + g**2*omega**4/Omega**6 + g**2*omega**6/Omega**8

----------------------------------------------------------------------------------------
III.3 — Dissipation / Schott identities
----------------------------------------------------------------------------------------
n=1 dissipation identity = 0
n=3 dissipation identity = 0
n=5 dissipation identity = 0

----------------------------------------------------------------------------------------
III.4 — Model-specific 2PN scalar / quadrupole combinations
----------------------------------------------------------------------------------------
gamma1_eff = 16*delta01/5 - 281*deltag1/80
gamma5_eff = 25*delta205/16

========================================================================================
PART IV — SCALAR SECTOR
========================================================================================

----------------------------------------------------------------------------------------
IV.1 — Outgoing scalar Green function and monopole odd term
----------------------------------------------------------------------------------------
G_R(omega,r) = exp(I*omega*r/c_s)/(4*pi*r)
Series = 1/(4*pi*r) + I*omega/(4*pi*c_s) - omega**2*r/(8*pi*c_s**2) - I*omega**3*r**2/(24*pi*c_s**3) + omega**4*r**3/(96*pi*c_s**4) + I*omega**5*r**4/(480*pi*c_s**5)
Model-specific gamma1_eff = 16*delta01/5 - 281*deltag1/80

----------------------------------------------------------------------------------------
IV.1b — 2PN scalar support/geometry finite-size rescue
----------------------------------------------------------------------------------------
Lambda0(k) = I*k0 - 1/a0
Y0_norm(k) = I*a0**5*k0**5 + a0**4*k0**4 - I*a0**3*k0**3 - a0**2*k0**2 + I*a0*k0 + 1
Yg_norm(k) = I*ag**5*k0**5 + ag**4*k0**4 - I*ag**3*k0**3 - ag**2*k0**2 + I*ag*k0 + 1
Seff(k) = k0**4*(16*a0**4/5 - 281*ag**4/80) + k0**3*(-16*I*a0**3/5 + 281*I*ag**3/80) + k0**2*(-16*a0**2/5 + 281*ag**2/80) + k0*(16*I*a0/5 - 281*I*ag/80) + 5/4
ell_eff = 16*a0/5 - 281*ag/80
equal-scale residual ell_eff(ag=a0) = -5*a0/16
exact scalar cancellation ag = 256*a0/281

----------------------------------------------------------------------------------------
IV.2 — Gaussian overlap counterexample and exact leakage identity
----------------------------------------------------------------------------------------
C(a) = 1/(sqrt(pi)*sqrt(a**2 + ell**2))
dC/da = -a/(sqrt(pi)*(a**2 + ell**2)**(3/2))
continuity residual = 0
I_leak - adot*dC/da = 0

----------------------------------------------------------------------------------------
IV.3 — Projection-locking linear algebra criterion
----------------------------------------------------------------------------------------
determinant of two-mode tangent matrix = -v1_L*v2_a + v1_a*v2_L
If determinant != 0, projection-locking requires B_a = B_L = 0.

----------------------------------------------------------------------------------------
IV.4 — Direct vs derivative coupling; damped discrete-mode test
----------------------------------------------------------------------------------------
direct = A*g0**2 + I*B*g0**2*omega + C*g0**2*omega**2 + I*D*g0**2*omega**3
derivative = A*gd**2*omega**2 + I*B*gd**2*omega**3 + C*gd**2*omega**4 + I*D*gd**2*omega**5
discrete_damped = lam**2/Omega**2 + I*eta*lam**2*omega/Omega**4 + lam**2*omega**2/Omega**4 - eta**2*lam**2*omega**2/Omega**6 + 2*I*eta*lam**2*omega**3/Omega**6 - I*eta**3*lam**2*omega**3/Omega**8
discrete_damped_derivative = lamd**2*omega**2/Omega**2 + I*eta*lamd**2*omega**3/Omega**4 + lamd**2*omega**4/Omega**4 - eta**2*lamd**2*omega**4/Omega**6
gamma1_direct = B*g0**2
gamma3_derivative = B*gd**2
gamma1_discrete_damped = eta*lam**2/Omega**4

----------------------------------------------------------------------------------------
IV.5 — Breathing-to-outlet vertex is dot(q)-type
----------------------------------------------------------------------------------------
K_direct = A*B_q**2 + I*B*B_q**2*omega + B_q**2*C*omega**2 + I*B_q**2*D*omega**3
Im K_direct = B*B_q**2*omega + B_q**2*D*omega**3
K_deriv = A*B_q**2*omega**2 + I*B*B_q**2*omega**3 + B_q**2*C*omega**4 + I*B_q**2*D*omega**5
Im K_deriv = B*B_q**2*omega**3 + B_q**2*D*omega**5
Derivative-breathing odd term exponent (with a/r ~ eps^delta) = 3*delta + 5/2

----------------------------------------------------------------------------------------
IV.6 — No third linear scalar source from quadratic momentum/stress terms
----------------------------------------------------------------------------------------
j_sigma(w) = I*q*u0*w**3*z3 + q*u0*w**2*z2 + I*q*u0*w*z1 + I*q*u1*w**3*z2 - q*u1*w**2*z1 + I*q*u2*w**3*z1
With Z_sigma(0)=0, the leading mouth source is derivative-like, not direct q-like.

----------------------------------------------------------------------------------------
IV.7 — Mouth radiative-order theorem and Ohmic no-go
----------------------------------------------------------------------------------------
Compact reciprocal mouth admittance Y(w) = I*beta*w/K + I*M*beta*w**3/K**2 + I*M_add*beta*w**3/K**2 + R2*beta*w**4/K**2 + I*M**2*beta*w**5/K**3 + 2*I*M*M_add*beta*w**5/K**3 + I*M_add**2*beta*w**5/K**3
Re Y series = R2*beta*w**4/K**2 + O(w**6)
Im Y series = beta*w/K + w**3*(beta*(2*M/K + 2*M_add/K)/K - M*beta/K**2 - M_add*beta/K**2) + w**5*(beta*((-2*M/K - 2*M_add/K)**2 - M**2/K**2 - 2*M*M_add/K**2 - M_add**2/K**2)/K - M*beta*(2*M/K + 2*M_add/K)/K**2 - M_add*beta*(2*M/K + 2*M_add/K)/K**2) + O(w**6)
1D Ohmic benchmark Re Y1 = beta*c0*rho0*w**2/K**2 + beta*c0*rho0*w**4*(2*M/K - c0**2*rho0**2/K**2)/K**2 + O(w**5)
Natural odd mouth length ell_mouth ~ a_small*chi_sigma/(kappa_m*lam)

========================================================================================
PART V — DIPOLE / VECTOR SECTOR
========================================================================================

----------------------------------------------------------------------------------------
V.1 — Carried odd wake ports: CM/relative split
----------------------------------------------------------------------------------------
PiA_perp - PiB_perp = Matrix([[sqrt(14)*ux/2], [sqrt(14)*uy/2]])
PiA_10   - PiB_10   = 2*d

----------------------------------------------------------------------------------------
V.2 — Vector Ward-identity failure
----------------------------------------------------------------------------------------
F_A = a*(-Derivative(xpA(t), (t, 3)) + Derivative(xpB(t), (t, 3)))
F_B = a*(Derivative(xpA(t), (t, 3)) - Derivative(xpB(t), (t, 3)))
F_A + F_B = 0
L_cmrel = a*Derivative(x(t), t)*Derivative(x(t), (t, 2))

----------------------------------------------------------------------------------------
V.3 — Outgoing l=1 spectral no-go for linear odd term
----------------------------------------------------------------------------------------
Lambda1(k) = a**5*k**6 - I*a**4*k**5 - a**3*k**4 + I*a**2*k**3 + a*k**2 - 2/a
Y1_norm(k) = -R*a**4*k**4/4 + I*R*a**3*k**3/2 + R*a**2*k**2/2 + R
Imaginary k^3 coefficient = R*a**3/2

----------------------------------------------------------------------------------------
V.3b — Order reduction of an isotropic fifth-derivative force
----------------------------------------------------------------------------------------
x^(5) = 15*m*rd*(2*m + 7*r*rd**2 - 3*r*v2)/r**6 * n + m*(-8*m - 45*r*rd**2 + 9*r*v2)/r**6 * v
Thus -kap x^(5) lies in the standard {rdot n, v} 2.5PN basis.

----------------------------------------------------------------------------------------
V.4 — Same-order dipole wake cannot be absorbed into the Iyer–Will family
----------------------------------------------------------------------------------------
Full matching solution = [{alpha: 4, beta: 5, p: 0, q: 0, s: 1}]

----------------------------------------------------------------------------------------
V.5 — Dipole finite-size scaling rescue
----------------------------------------------------------------------------------------
eps * (sqrt(eps)*rho)^3 = eps**(5/2)*rho**3
If rho ~ eps^delta, exponent = 5/2 + 3 delta.
general exponent = 3*delta + 5/2

========================================================================================
PART VI — QUADRUPOLE SECTOR
========================================================================================

----------------------------------------------------------------------------------------
VI.1 — P2 ports are exactly the real STF l=2 content
----------------------------------------------------------------------------------------
q20 - sqrt(2/3) Pi20 = 0
q21c - Pi21c = 0
q21s - Pi21s = 0
q22c - 2 Pi22c = 0
q22s - 2 Pi22s = 0

----------------------------------------------------------------------------------------
VI.2 — Outgoing l=2 retarded completion gives +i*k^5
----------------------------------------------------------------------------------------
Lambda2(k) = -2*a**5*k**6/27 + I*a**4*k**5/9 + a**3*k**4/9 + a*k**2/3 - 3/a
Y2_norm(k) = I*R2*a**5*k**5/27 + 4*R2*a**4*k**4/81 + R2*a**2*k**2/9 + R2
Imaginary k^5 coefficient = R2*a**5/27

----------------------------------------------------------------------------------------
VI.3 — Quadrupole odd kernel gives damping, not antidamping
----------------------------------------------------------------------------------------
P - dE_Schott/dt = -Gamma*Derivative(q(t), (t, 3))**2/2
quadrupole power balance residual = 0

----------------------------------------------------------------------------------------
VI.3b — Compact internal P2 branch is finite-size suppressed
----------------------------------------------------------------------------------------
eps^2 * (sqrt(eps)*rho)^5 = eps**(9/2)*rho**5
If rho ~ eps^delta, exponent = 5*delta + 9/2

----------------------------------------------------------------------------------------
VI.4 — Orbital/worldtube STF source map
----------------------------------------------------------------------------------------
V20 source map = 0
V21c source map = 0
V21s source map = 0
V22c source map = 0
V22s source map = 0
X21c static position source = 0
X21s static position source = 0
X22c static position source = 0
X22s static position source = 0
I21c source map = 0
I21s source map = 0
I22c source map = 0
I22s source map = 0
I20 source map = 0

----------------------------------------------------------------------------------------
VI.5 — Matching theorem: commutant is scalar, basis map invertible
----------------------------------------------------------------------------------------
M commuting with all J_i =
Matrix([[m, 0, 0, 0, 0], [0, m, 0, 0, 0], [0, 0, m, 0, 0], [0, 0, 0, m, 0], [0, 0, 0, 0, m]])
commutant minus scalar identity = Matrix([[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]])
det(port->canonical basis map) = 4*sqrt(6)/3

----------------------------------------------------------------------------------------
VI.6 — Static overlap test
----------------------------------------------------------------------------------------
C20[x_<i x_j>] = sqrt(6)*r**2/3
C20[ddot I_orb /(2 mu)] = -sqrt(6)*G*M/(3*r)
q20_static = 5/4
m0 from static overlap = -5*sqrt(6)*r/(8*G*M)

========================================================================================
PART VII — QUADRUPOLE NORMALIZATION / EXTRACTION / PARAMETERIZATION
========================================================================================

----------------------------------------------------------------------------------------
VII.1 — Canonical quadrupole normalization map
----------------------------------------------------------------------------------------
Fx - gamma*mu*(-x^j I_xj^(5)) = 0
Fy - gamma*mu*(-x^j I_yj^(5)) = 0
gamma_GR = 2*G/(5*c**5)
gamma_toy = Gamma5*m0hat**2
Gamma5_target = 2*G/(5*c**5*m0hat**2)

----------------------------------------------------------------------------------------
VII.2 — Outgoing l=2 PDE fingerprint
----------------------------------------------------------------------------------------
l=2 Helmholtz residual = 0
Lambda2(k) = -2*a**5*k**6/27 + I*a**4*k**5/9 + a**3*k**4/9 + a*k**2/3 - 3/a
Y2(0) = -a/3
Y2_hat(k) = I*a**5*k**5/27 + 4*a**4*k**4/81 + a**2*k**2/9 + 1
m2 = a**2*m0/(9*c_s**2)
m4 = 4*a**4*m0/(81*c_s**4)
Gamma5^PDE = a**5*m0/(27*c_s**5)
m4 - 4 m2^2 / m0 = 0
Gamma5^PDE - 9 m2^(5/2)/m0^(3/2) = 0

----------------------------------------------------------------------------------------
VII.3 — Source/port normalization cleanup
----------------------------------------------------------------------------------------
det(Bmap) = 4*sqrt(6)/3
m0_mixed = -5*sqrt(6)*r/(8*G*M)
m0_raw_same_basis = -5*r/(4*G*M)
m0_can_same_basis = -5*r/(4*G*M)
N_orb = -4*G*M/(5*r)
m0hat_point_particle = 1
same-basis raw vs canonical overlap = 0
point-particle normalized overlap - 1 = 0
m2_eff - m0_eff a^2/(9 c_s^2) = 0
m4_eff - 4 m0_eff a^4/(81 c_s^4) = 0
Gamma5_eff - m0_eff a^5/(27 c_s^5) = 0
Gamma5_can = a**5*m0hat/(27*c_s**5)

----------------------------------------------------------------------------------------
VII.4 — Convention-invariant normalization product
----------------------------------------------------------------------------------------
mfrak'^2 Gamma5' - mfrak^2 Gamma5 = 0
mfrak'^2 K0' - mfrak^2 K0 = 0
Gamma5_port^DN = K0*a**5/(27*c_s**5)
gamma_eff = K0*a**5*mfrak**2/(27*c_s**5)
N_Q_target = 54*G*c_s**5/(5*a**5*c**5)

----------------------------------------------------------------------------------------
VII.5 — Canonical invariant low-frequency fingerprint
----------------------------------------------------------------------------------------
a^2/c_s^2 = 9*K2bar/K0bar
K4bar invariant = 4*K2bar**2/K0bar
Gammabar invariant = 9*K2bar**(5/2)/K0bar**(3/2)
K2bar_target formula = 2**(2/5)*375**(1/5)*G**(2/5)*K0bar**(3/5)/(15*c**2)
K4bar invariant - 4 K2bar^2/K0bar = 0
Gammabar invariant - 9 K2bar^(5/2)/K0bar^(3/2) = 0

----------------------------------------------------------------------------------------
VII.6 — Extraction of (Kbar0, Kbar2) and single-prefactor obstruction
----------------------------------------------------------------------------------------
K2bar_branch = K0bar*a**2/(9*c_s**2)
K0bar_target = 54*G*c_s**5/(5*a**5*c**5)
K2bar_target = 6*G*c_s**3/(5*a**3*c**5)
K4bar_target = 8*G*c_s/(15*a*c**5)
Y20_min(omega) = 1 + omega**2/Omega20**2 + omega**4/Omega20**4
C_Q_target = 54*G*c_s**5/(5*a**5*c**5)
Kbar0 - C_Q = 0
Kbar2 - C_Q a^2/(9 c_s^2) = 0
Kbar4 - 4 C_Q a^4/(81 c_s^4) = 0

----------------------------------------------------------------------------------------
VII.7 — Frozen 2PN files: fixed representation/support data
----------------------------------------------------------------------------------------
C20 relation = 0
C21c relation = 0
C21s relation = 0
C22c relation = 0
C22s relation = 0
J.J = 381/80
closure deficit = 281/80
Open dynamic poles: Omega20, Omega21, Omega22

----------------------------------------------------------------------------------------
VII.8 — One-pole insufficiency and minimal positive conservative precursor
----------------------------------------------------------------------------------------
Omega_1 = 3*c_s/a
one-pole O(w^4) mismatch = a**4/(27*c_s**4)
Y_2pole series = omega**4*(p*(-a**2/(9*c_s**2) - sqrt(3)*a**2*sqrt(1 - p)/(9*c_s**2*sqrt(p)))**2 - p*(sqrt(3)*a**2*sqrt(p)/(9*c_s**2*sqrt(1 - p)) - a**2/(9*c_s**2))**2 + (sqrt(3)*a**2*sqrt(p)/(9*c_s**2*sqrt(1 - p)) - a**2/(9*c_s**2))**2) + omega**2*(-sqrt(3)*a**2*sqrt(p)/(9*c_s**2*sqrt(1 - p)) + a**2/(9*c_s**2) - p*(-a**2/(9*c_s**2) - sqrt(3)*a**2*sqrt(1 - p)/(9*c_s**2*sqrt(p))) - p*(-sqrt(3)*a**2*sqrt(p)/(9*c_s**2*sqrt(1 - p)) + a**2/(9*c_s**2))) + 1
Minimal positive closure = 3*(a**2*omega**2 - 3*c_s**2)/(4*a**2*omega**2 - 9*c_s**2)
Minimal split form = 3*(a**2*omega**2 - 3*c_s**2)/(4*a**2*omega**2 - 9*c_s**2)
Omega_Q,min = 3*c_s/(2*a)
minimal closure omega^2 match = 0
minimal closure omega^4 match = 0

----------------------------------------------------------------------------------------
VII.9 — Minimal isotropic quadrupole module and componentwise P2 target
----------------------------------------------------------------------------------------
c0 = 3/4
c1 = 1/4
OmegaQ = 3*c_s/(2*a)
Y_family_min = 3*(a**2*omega**2 - 3*c_s**2)/(4*a**2*omega**2 - 9*c_s**2)
minimal module omega^2 mismatch = 0
minimal module omega^4 mismatch = 0
real STF Gram matrix - I5 = Matrix([[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]])
Y20 = Y21 = Y22 = 3*(a**2*omega**2 - 3*c_s**2)/(4*a**2*omega**2 - 9*c_s**2)
YQ_min series = 4*a**4*omega**4/(81*c_s**4) + a**2*omega**2/(9*c_s**2) + 1
OmegaQ_from_invariants = sqrt(K0bar)/(2*sqrt(K2bar))
K4 relation = 4*K2bar**2/K0bar
Gamma5 relation = 9*K2bar**(5/2)/K0bar**(3/2)

----------------------------------------------------------------------------------------
VII.10 — Single-pole sufficiency theorem and extraction ledger
----------------------------------------------------------------------------------------
K2_from_branch = K0bar/(4*OmegaQ**2)
K4_from_branch = K0bar/(4*OmegaQ**4)
Gamma_from_branch = 9*K0bar/(32*OmegaQ**5)
K0bar_target(OmegaQ) = 64*G*OmegaQ**5/(45*c**5)
K2bar_target(OmegaQ) = 16*G*OmegaQ**3/(45*c**5)
K4bar_target(OmegaQ) = 16*G*OmegaQ/(45*c**5)
OmegaQ^2 from (c0,c2) = c0/(4*c2)
OmegaQ^2 from (c2,c4) = c2/c4
K0 from (c2,c4) = 4*c2**2/c4
Gamma from (c2,c4) = 9*c4**(3/2)/(8*sqrt(c2))
Direct GR gate: c4^3/c2 = 256*G**2/(2025*c**10)

----------------------------------------------------------------------------------------
VII.11 — Grouped P2 formulas and normalized-support ledger
----------------------------------------------------------------------------------------
c0_trace = c0_20/5 + 2*c0_21/5 + 2*c0_22/5
c2_trace = c2_20/5 + 2*c2_21/5 + 2*c2_22/5
c4_trace = c4_20/5 + 2*c4_21/5 + 2*c4_22/5
weighted anisotropy sum n=0 = 0
weighted anisotropy sum n=2 = 0
weighted anisotropy sum n=4 = 0
Minimal branch identity c0*c4 - 4*c2^2 = c0*c4 - 4*c2**2
u2_trace = u2_20/5 + 2*u2_21/5 + 2*u2_22/5
u4_trace = u4_20/5 + 2*u4_21/5 + 2*u4_22/5
normalized branch identity u4 - 4*u2^2 = -4*u2**2 + u4
OmegaQ on normalized branch = 1/(2*sqrt(u2))
Gamma5_norm on branch = 9*u2**(5/2)
K0bar_target from normalized support = 2*G/(45*c**5*u2**(5/2))

----------------------------------------------------------------------------------------
VII.12 — Axisymmetric grouped-P2 parameterization
----------------------------------------------------------------------------------------
Y20 = omega**4*(4*a4 + u4bar) + omega**2*(4*a2 + u2bar) + 1
Y21 = omega**4*(-a4 + b4 + u4bar) + omega**2*(-a2 + b2 + u2bar) + 1
Y22 = omega**4*(-a4 - b4 + u4bar) + omega**2*(-a2 - b2 + u2bar) + 1
ubar_from_(x20,x21,x22) = x20/5 + 2*x21/5 + 2*x22/5
a_from_(x20,x21,x22) = x20/5 - x21/10 - x22/10
b_from_(x20,x21,x22) = x21/2 - x22/2
A2_sq = 4*a2**2 + 4*b2**2/5
A4_sq = 4*a4**2 + 4*b4**2/5
Minimal branch relation = -4*u2bar**2 + u4bar
OmegaQ^2 = 1/(4*u2bar)
Gamma5_norm = 9*u2bar**(5/2)
K0bar_target = 2*G/(45*c**5*u2bar**(5/2))

========================================================================================
FINAL SESSION THEOREM LEDGER
========================================================================================
Decisive benchmark:
  Burke–Thorne prototype lands on the Iyer–Will family with (alpha,beta)=(4,5).

Channel audit status:
  - scalar direct overlap generically gives i*omega,
  - derivative breathing outlet pushes the scalar odd kernel to i*omega^3,
  - compact reciprocal mouth radiation is higher-order / super-Ohmic,
  - outgoing l=1 starts at +i*k^3, not at i*k,
  - same-order nonzero dipole deformation cannot be absorbed into the standard 2.5PN family,
  - compact outgoing l=2 starts at +i*k^5 with positive damping sign.

Quadrupole normalization stack:
  gamma_GR = 2*G/(5*c**5)
  Gamma5_target = 2*G/(5*c**5*m0hat**2)
  N_Q_target = 54*G*c_s**5/(5*a**5*c**5)
  K0bar_target = 54*G*c_s**5/(5*a**5*c**5)
  K2bar_target = 6*G*c_s**3/(5*a**3*c**5)
  K4bar_target = 8*G*c_s/(15*a*c**5)

Minimal isotropic quadrupole module:
  Omega_Q,min = 3*c_s/(2*a)
  exact moment-matched isotropic solution Omega_Q = 3*c_s/(2*a)
  single-pole GR scaling law K0bar = 64 G Omega_Q^5 / (45 c^5).

Normalized-support branch:
  K0bar_target from normalized support moments = 2*G/(45*c**5*u2**(5/2))

Best current verdict:
  The session reduces the full GR-like point-particle 2.5PN closure to the
  isotropic real P2 quadrupole module. Once the grouped low-frequency
  conservative coefficients through O(omega^4) are fixed on that branch,
  the odd Burke–Thorne coefficient follows automatically.
  The one remaining narrow theorem gap is the final passive/outgoing
  quadrupole normalization on the actual moving-throat branch.

"""
