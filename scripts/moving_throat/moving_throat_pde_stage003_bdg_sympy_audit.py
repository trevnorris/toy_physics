
#!/usr/bin/env python3
"""
moving_throat_pde_stage3_bdg_sympy_audit.py

SymPy audit for the first matter-coupled moving-throat PDE reduction.

Scope
-----
This script backs the first genuinely coupled step after the geometry-only and
breathing-reduction notes:

  • axisymmetric (a,L)-type wall coordinates coupled to stable scalar BdG modes,
  • grouped real P2 wall coordinates coupled to stable quadrupole BdG modes,
  • exact elimination of the matter modes to obtain the effective wall kernels,
  • low-frequency coefficient shifts,
  • exact two-pole formula and perturbative pole shift in the one-mode case,
  • grouped-P2 isotropy/anomaly bookkeeping under channelwise coupling data.

This is still a controlled reduced-sector calculation, not a full solution of the
coupled GNLS + localized Maxwell + moving-wall PDE.
"""

from __future__ import annotations

import sympy as sp
from sympy.calculus.euler import euler_equations


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
        simplified = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.simplify(sp.expand(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Section I. Axisymmetric wall + scalar BdG normal-mode reduction
# ---------------------------------------------------------------------------

def axisymmetric_matrix_kernel_audit() -> None:
    banner("SECTION I — AXISYMMETRIC WALL + SCALAR BdG REDUCTION")

    # Coordinates and time
    t, omega = sp.symbols("t omega", real=True)

    # Wall coordinates Q^A = (qa, qL)
    qa = sp.Function("qa")(t)
    qL = sp.Function("qL")(t)

    # Two stable scalar matter modes X_alpha
    xa = sp.Function("xa")(t)
    xb = sp.Function("xb")(t)

    Maa, MaL, MLL = sp.symbols("Maa MaL MLL", real=True)
    Kaa, KaL, KLL = sp.symbols("Kaa KaL KLL", real=True)
    wa, wb = sp.symbols("wa wb", positive=True, real=True)
    c1a, c1b, c2a, c2b = sp.symbols("c1a c1b c2a c2b", real=True)

    Lred = (
        sp.Rational(1, 2) * Maa * sp.diff(qa, t) ** 2
        + MaL * sp.diff(qa, t) * sp.diff(qL, t)
        + sp.Rational(1, 2) * MLL * sp.diff(qL, t) ** 2
        - sp.Rational(1, 2) * Kaa * qa ** 2
        - KaL * qa * qL
        - sp.Rational(1, 2) * KLL * qL ** 2
        + sp.Rational(1, 2) * sp.diff(xa, t) ** 2
        - sp.Rational(1, 2) * wa**2 * xa**2
        + sp.Rational(1, 2) * sp.diff(xb, t) ** 2
        - sp.Rational(1, 2) * wb**2 * xb**2
        + c1a * qa * xa
        + c1b * qa * xb
        + c2a * qL * xa
        + c2b * qL * xb
    )

    EL_qa = euler_equations(Lred, qa, [t])[0]
    EL_qL = euler_equations(Lred, qL, [t])[0]
    EL_xa = euler_equations(Lred, xa, [t])[0]
    EL_xb = euler_equations(Lred, xb, [t])[0]

    subbanner("I.1 — Euler–Lagrange equations")
    expect_zero(
        "qa equation",
        EL_qa.lhs + Maa * sp.diff(qa, t, 2) + MaL * sp.diff(qL, t, 2) + Kaa * qa + KaL * qL - c1a * xa - c1b * xb,
    )
    expect_zero(
        "qL equation",
        EL_qL.lhs + MaL * sp.diff(qa, t, 2) + MLL * sp.diff(qL, t, 2) + KaL * qa + KLL * qL - c2a * xa - c2b * xb,
    )
    expect_zero(
        "xa equation",
        EL_xa.lhs + sp.diff(xa, t, 2) + wa**2 * xa - c1a * qa - c2a * qL,
    )
    expect_zero(
        "xb equation",
        EL_xb.lhs + sp.diff(xb, t, 2) + wb**2 * xb - c1b * qa - c2b * qL,
    )

    subbanner("I.2 — Frequency-space elimination")
    Mmat = sp.Matrix([[Maa, MaL], [MaL, MLL]])
    Kmat = sp.Matrix([[Kaa, KaL], [KaL, KLL]])
    Cmat = sp.Matrix([[c1a, c1b], [c2a, c2b]])
    Omat = sp.diag(wa**2, wb**2)
    Deff = sp.simplify(Kmat - omega**2 * Mmat - Cmat * (Omat - omega**2 * sp.eye(2)).inv() * Cmat.T)
    print("Effective axisymmetric kernel D0_eff(omega) =")
    sp.pprint(Deff)

    manual = sp.Matrix(
        [
            [
                Kaa - Maa * omega**2 - c1a**2 / (wa**2 - omega**2) - c1b**2 / (wb**2 - omega**2),
                KaL - MaL * omega**2 - c1a * c2a / (wa**2 - omega**2) - c1b * c2b / (wb**2 - omega**2),
            ],
            [
                KaL - MaL * omega**2 - c1a * c2a / (wa**2 - omega**2) - c1b * c2b / (wb**2 - omega**2),
                KLL - MLL * omega**2 - c2a**2 / (wa**2 - omega**2) - c2b**2 / (wb**2 - omega**2),
            ],
        ]
    )
    expect_zero("D0_eff - manual form", sp.simplify(Deff - manual))

    subbanner("I.3 — Low-frequency expansion")
    Deff_series = Deff.applyfunc(lambda z: sp.expand(sp.series(z, omega, 0, 5).removeO()))
    print("D0_eff(omega) through O(omega^4) =")
    sp.pprint(Deff_series)

    Keff = sp.Matrix(
        [
            [Kaa - c1a**2 / wa**2 - c1b**2 / wb**2, KaL - c1a * c2a / wa**2 - c1b * c2b / wb**2],
            [KaL - c1a * c2a / wa**2 - c1b * c2b / wb**2, KLL - c2a**2 / wa**2 - c2b**2 / wb**2],
        ]
    )
    Meff = sp.Matrix(
        [
            [Maa + c1a**2 / wa**4 + c1b**2 / wb**4, MaL + c1a * c2a / wa**4 + c1b * c2b / wb**4],
            [MaL + c1a * c2a / wa**4 + c1b * c2b / wb**4, MLL + c2a**2 / wa**4 + c2b**2 / wb**4],
        ]
    )
    Neff = sp.Matrix(
        [
            [c1a**2 / wa**6 + c1b**2 / wb**6, c1a * c2a / wa**6 + c1b * c2b / wb**6],
            [c1a * c2a / wa**6 + c1b * c2b / wb**6, c2a**2 / wa**6 + c2b**2 / wb**6],
        ]
    )
    target_series = Keff - omega**2 * Meff - omega**4 * Neff
    expect_zero("series match", Deff_series - target_series)

    print("So the scalar BdG sector renormalizes:")
    print("  K_AB -> K_AB - sum_alpha C_Aa C_Ba / varpi_a^2")
    print("  M_AB -> M_AB + sum_alpha C_Aa C_Ba / varpi_a^4")
    print("  N_AB ->          sum_alpha C_Aa C_Ba / varpi_a^6")


# ---------------------------------------------------------------------------
# Section II. One wall mode + one stable BdG mode: exact pole formula
# ---------------------------------------------------------------------------

def one_mode_pole_shift_audit() -> None:
    banner("SECTION II — ONE WALL MODE + ONE STABLE BdG MODE")

    M, K, varpi2, g, eps = sp.symbols("M K varpi2 g eps", positive=True, real=True)
    Om2 = sp.symbols("Omega_eta2", positive=True, real=True)
    w2 = sp.symbols("w2", real=True)

    # Exact reduced dispersion relation
    dispersion = sp.expand((K - M * w2) * (varpi2 - w2) - g**2)
    roots = sp.solve(sp.Eq(dispersion, 0), w2)
    roots = [sp.simplify(r.subs(K, M * Om2)) for r in roots]

    subbanner("II.1 — Exact roots")
    print("omega_-^2 =")
    sp.pprint(roots[0])
    print("omega_+^2 =")
    sp.pprint(roots[1])

    closed_form_minus = sp.simplify((Om2 + varpi2 - sp.sqrt((Om2 - varpi2) ** 2 + 4 * g**2 / M)) / 2)
    closed_form_plus = sp.simplify((Om2 + varpi2 + sp.sqrt((Om2 - varpi2) ** 2 + 4 * g**2 / M)) / 2)
    expect_zero("minus-root closed form", roots[0] - closed_form_minus)
    expect_zero("plus-root closed form", roots[1] - closed_form_plus)

    subbanner("II.2 — Perturbative wall-like shift")
    delta = sp.symbols("delta", positive=True, real=True)
    # Assume the matter mode is above the wall mode: varpi2 = Om2 + delta
    wall_like = sp.simplify(closed_form_minus.subs(varpi2, Om2 + delta).subs(g, eps * g))
    matter_like = sp.simplify(closed_form_plus.subs(varpi2, Om2 + delta).subs(g, eps * g))
    wall_series = sp.expand(sp.series(wall_like, eps, 0, 3).removeO())
    matter_series = sp.expand(sp.series(matter_like, eps, 0, 3).removeO())
    print("wall-like pole through O(g^2) =")
    sp.pprint(wall_series)
    print("matter-like pole through O(g^2) =")
    sp.pprint(matter_series)

    expect_zero("wall-like shift", wall_series - (Om2 - eps**2 * g**2 / (M * delta)))
    expect_zero("matter-like shift", matter_series - (Om2 + delta + eps**2 * g**2 / (M * delta)))

    print("Therefore, for varpi^2 > Omega_eta^2,")
    print("  delta Omega_eta^2 = - g^2 / [ M (varpi^2 - Omega_eta^2) ] + O(g^4)")
    print("  delta varpi^2     = + g^2 / [ M (varpi^2 - Omega_eta^2) ] + O(g^4)")


# ---------------------------------------------------------------------------
# Section III. Grouped real P2 channels with channelwise BdG couplings
# ---------------------------------------------------------------------------

def grouped_p2_coupling_audit() -> None:
    banner("SECTION III — GROUPED REAL P2 + BdG COUPLINGS")

    omega = sp.symbols("omega", real=True)
    K20, K21, K22 = sp.symbols("K20 K21 K22", real=True)
    M20, M21, M22 = sp.symbols("M20 M21 M22", real=True)
    g20, g21, g22 = sp.symbols("g20 g21 g22", real=True)
    w20, w21, w22 = sp.symbols("w20 w21 w22", positive=True, real=True)

    D20 = sp.simplify(K20 - M20 * omega**2 - g20**2 / (w20**2 - omega**2))
    D21 = sp.simplify(K21 - M21 * omega**2 - g21**2 / (w21**2 - omega**2))
    D22 = sp.simplify(K22 - M22 * omega**2 - g22**2 / (w22**2 - omega**2))

    subbanner("III.1 — Channelwise low-frequency coefficients")
    D20s = sp.expand(sp.series(D20, omega, 0, 5).removeO())
    D21s = sp.expand(sp.series(D21, omega, 0, 5).removeO())
    D22s = sp.expand(sp.series(D22, omega, 0, 5).removeO())
    print("D20 =", D20s)
    print("D21 =", D21s)
    print("D22 =", D22s)

    d020 = sp.simplify(D20s.coeff(omega, 0))
    d021 = sp.simplify(D21s.coeff(omega, 0))
    d022 = sp.simplify(D22s.coeff(omega, 0))
    d220 = sp.simplify(D20s.coeff(omega, 2))
    d221 = sp.simplify(D21s.coeff(omega, 2))
    d222 = sp.simplify(D22s.coeff(omega, 2))
    d420 = sp.simplify(D20s.coeff(omega, 4))
    d421 = sp.simplify(D21s.coeff(omega, 4))
    d422 = sp.simplify(D22s.coeff(omega, 4))

    print("d0-shifts:")
    print("  d0^(20) =", d020)
    print("  d0^(21) =", d021)
    print("  d0^(22) =", d022)
    print("d2-shifts:")
    print("  d2^(20) =", d220)
    print("  d2^(21) =", d221)
    print("  d2^(22) =", d222)
    print("d4-shifts:")
    print("  d4^(20) =", d420)
    print("  d4^(21) =", d421)
    print("  d4^(22) =", d422)

    subbanner("III.2 — Grouped trace / anisotropy diagnostics")
    d2bar = sp.simplify((d220 + 2 * d221 + 2 * d222) / 5)
    a2 = sp.simplify((2 * d220 - d221 - d222) / 10)
    b2 = sp.simplify((d221 - d222) / 2)

    print("d2bar =")
    sp.pprint(d2bar)
    print("a2 =")
    sp.pprint(a2)
    print("b2 =")
    sp.pprint(b2)

    subbanner("III.3 — Isotropic matter-coupling branch")
    iso_subs = {
        K20: sp.Symbol("K2"), K21: sp.Symbol("K2"), K22: sp.Symbol("K2"),
        M20: sp.Symbol("M2"), M21: sp.Symbol("M2"), M22: sp.Symbol("M2"),
        g20: sp.Symbol("g2"), g21: sp.Symbol("g2"), g22: sp.Symbol("g2"),
        w20: sp.Symbol("w2"), w21: sp.Symbol("w2"), w22: sp.Symbol("w2"),
    }
    expect_zero("isotropic a2", a2.subs(iso_subs))
    expect_zero("isotropic b2", b2.subs(iso_subs))
    expect_zero("isotropic D20-D21", sp.simplify((D20 - D21).subs(iso_subs)))
    expect_zero("isotropic D21-D22", sp.simplify((D21 - D22).subs(iso_subs)))

    print("So identical channelwise wall data and identical BdG couplings/frequencies preserve grouped-P2 degeneracy.")
    print("Any nonzero a2 or b2 must come from anisotropy in K_A, M_A, g_A, or varpi_A.")


# ---------------------------------------------------------------------------
# Section IV. Coupling selection rule under isotropic background
# ---------------------------------------------------------------------------

def selection_rule_audit() -> None:
    banner("SECTION IV — HARMONIC SELECTION RULE FOR THE CONFINEMENT COUPLING")

    th, ph = sp.symbols("theta phi", real=True)
    dOmega = sp.sin(th)

    Y00 = sp.Rational(1, 2) / sp.sqrt(sp.pi)
    Y20 = sp.sqrt(5) / (4 * sp.sqrt(sp.pi)) * (3 * sp.cos(th) ** 2 - 1)
    Y21c = sp.sqrt(15) / (2 * sp.sqrt(sp.pi)) * sp.sin(th) * sp.cos(th) * sp.cos(ph)
    Y22c = sp.sqrt(15) / (4 * sp.sqrt(sp.pi)) * sp.sin(th) ** 2 * sp.cos(2 * ph)

    I00_20 = sp.simplify(sp.integrate(sp.integrate(Y00 * Y20 * dOmega, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
    I20_21c = sp.simplify(sp.integrate(sp.integrate(Y20 * Y21c * dOmega, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
    I20_22c = sp.simplify(sp.integrate(sp.integrate(Y20 * Y22c * dOmega, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
    N20 = sp.simplify(sp.integrate(sp.integrate(Y20 * Y20 * dOmega, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))

    expect_zero("Y00-Y20 cross integral", I00_20)
    expect_zero("Y20-Y21c cross integral", I20_21c)
    expect_zero("Y20-Y22c cross integral", I20_22c)
    expect_zero("Y20 norm - 1", N20 - 1)

    print("With an isotropic background and scalar confinement coupling, the angular overlap is diagonal in (l,m).")
    print("So l=0 wall motion couples only to l=0 matter modes, and grouped real l=2 wall motion couples channelwise inside l=2.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    axisymmetric_matrix_kernel_audit()
    one_mode_pole_shift_audit()
    grouped_p2_coupling_audit()
    selection_rule_audit()

    banner("FINAL STAGE-3 LEDGER")
    print("Verified with SymPy:")
    print("  • the reduced (a,L) wall variables coupled to stable scalar BdG modes yield")
    print("    D0_eff(omega) = K - omega^2 M - C (Omega_m^2 - omega^2 I)^(-1) C^T;")
    print("  • low-frequency matter coupling renormalizes the scalar/geometry stiffness and inertia matrices as")
    print("    K_eff = K - C Omega_m^{-2} C^T and M_eff = M + C Omega_m^{-4} C^T;")
    print("  • a single wall mode coupled to a single BdG mode has the exact two-pole spectrum")
    print("    omega_±^2 = [Omega_eta^2 + varpi^2 ± sqrt((Omega_eta^2-varpi^2)^2 + 4 g^2/M)]/2;")
    print("  • for varpi^2 > Omega_eta^2, the wall-like pole shifts downward by")
    print("    delta Omega_eta^2 = - g^2 / [ M (varpi^2 - Omega_eta^2) ] + O(g^4);")
    print("  • grouped real P2 channels acquire channelwise self-energies")
    print("    -g_A^2/(varpi_A^2-omega^2), whose anisotropy is measured by the usual (trace, a2, b2) combinations;")
    print("  • identical grouped-P2 couplings/frequencies preserve a2=b2=0;")
    print("  • harmonic orthogonality enforces the l=0 / l=2 selection rule for the scalar confinement coupling on an isotropic background.")


if __name__ == "__main__":
    main()
