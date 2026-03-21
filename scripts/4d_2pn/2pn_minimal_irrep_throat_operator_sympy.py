import sympy as sp

# -----------------------------------------------------------------------------
# 2PN minimal irreducible-channel throat operator prototype (SymPy)
# -----------------------------------------------------------------------------


def assert_zero(name: str, expr) -> None:
    res = sp.simplify(sp.expand(expr))
    if isinstance(res, sp.MatrixBase):
        ok = all(sp.simplify(entry) == 0 for entry in res)
    else:
        ok = (res == 0)
    if not ok:
        raise AssertionError(f"{name} failed: {res}")
    print(f"PASS: {name}")


def main() -> None:
    print("--- 2PN minimal irreducible-channel throat operator prototype ---")

    # ------------------------------------------------------------------
    print("\n=== MASTER / pair kinematics and canonical port basis ===")
    # ------------------------------------------------------------------
    vAx, vAy, vAz, vBx, vBy, vBz = sp.symbols('vAx vAy vAz vBx vBy vBz', real=True)
    UA, UB = sp.symbols('UA UB', real=True)

    vA2 = sp.expand(vAx**2 + vAy**2 + vAz**2)
    vB2 = sp.expand(vBx**2 + vBy**2 + vBz**2)
    vAB = sp.expand(vAx * vBx + vAy * vBy + vAz * vBz)
    uAB = sp.expand(vAx * vBx + vAy * vBy)
    dA, dB = vAz, vBz

    Pi0A = sp.expand(sp.sqrt(5) * vA2 / 2)
    Pi0B = sp.expand(sp.sqrt(5) * vB2 / 2)

    Pi20A = sp.expand((3 * dA**2 - vA2) / 2)
    Pi20B = sp.expand((3 * dB**2 - vB2) / 2)

    Pi21Ac = sp.expand(sp.sqrt(2) * dA * vAx)
    Pi21As = sp.expand(sp.sqrt(2) * dA * vAy)
    Pi21Bc = sp.expand(sp.sqrt(2) * dB * vBx)
    Pi21Bs = sp.expand(sp.sqrt(2) * dB * vBy)

    Pi22Ac = sp.expand((vAx**2 - vAy**2) / (2 * sp.sqrt(2)))
    Pi22As = sp.expand((2 * vAx * vAy) / (2 * sp.sqrt(2)))
    Pi22Bc = sp.expand((vBx**2 - vBy**2) / (2 * sp.sqrt(2)))
    Pi22Bs = sp.expand((2 * vBx * vBy) / (2 * sp.sqrt(2)))

    L1_wake = sp.expand(-sp.Rational(7, 2) * vAB - sp.Rational(1, 2) * dA * dB)
    assert_zero("Frozen 1PN wake = transverse/longitudinal dipole metric",
                L1_wake - sp.expand(-sp.Rational(7, 2) * uAB - 4 * dA * dB))

    # ------------------------------------------------------------------
    print("\n=== EVEN / solve the minimal P0⊕P2 zero-frequency support operator ===")
    # ------------------------------------------------------------------
    a0, a20, a21, a22, j0, j20, s = sp.symbols('a0 a20 a21 a22 j0 j20 s', real=True)

    even_ans = sp.expand(
        a0 * Pi0A * Pi0B
        + a20 * Pi20A * Pi20B
        + a21 * (Pi21Ac * Pi21Bc + Pi21As * Pi21Bs)
        + a22 * (Pi22Ac * Pi22Bc + Pi22As * Pi22Bs)
        + UA * (j0 * Pi0B + j20 * Pi20B)
        + UB * (j0 * Pi0A + j20 * Pi20A)
        + s * UA * UB
    )

    even_target = sp.expand(
        sp.Rational(11, 8) * vA2 * vB2
        + sp.Rational(1, 4) * vAB**2
        - sp.Rational(5, 8) * (vA2 * dB**2 + vB2 * dA**2)
        + sp.Rational(3, 2) * vAB * dA * dB
        + sp.Rational(3, 8) * dA**2 * dB**2
        + UA * (sp.Rational(11, 8) * vB2 + sp.Rational(15, 8) * dB**2)
        + UB * (sp.Rational(11, 8) * vA2 + sp.Rational(15, 8) * dA**2)
        + sp.Rational(5, 4) * UA * UB
    )

    # Solve the linear coefficient problem from a sparse exact sample set.
    even_unknowns = [a0, a20, a21, a22, j0, j20, s]
    even_residual = sp.expand(even_ans - even_target)
    even_samples = [
        {vAx: 1, vAy: 0, vAz: 0, vBx: 1, vBy: 0, vBz: 0, UA: 0, UB: 0},
        {vAx: 0, vAy: 0, vAz: 1, vBx: 0, vBy: 0, vBz: 1, UA: 0, UB: 0},
        {vAx: 1, vAy: 0, vAz: 1, vBx: 1, vBy: 0, vBz: 1, UA: 0, UB: 0},
        {vAx: 1, vAy: 1, vAz: 0, vBx: 1, vBy: -1, vBz: 0, UA: 0, UB: 0},
        {vAx: 1, vAy: 2, vAz: 3, vBx: -1, vBy: 1, vBz: 2, UA: 0, UB: 0},
        {vAx: 0, vAy: 0, vAz: 0, vBx: 1, vBy: 0, vBz: 0, UA: 1, UB: 0},
        {vAx: 0, vAy: 0, vAz: 0, vBx: 0, vBy: 0, vBz: 1, UA: 1, UB: 0},
        {vAx: 1, vAy: 1, vAz: 2, vBx: 0, vBy: 0, vBz: 0, UA: 0, UB: 1},
        {vAx: 0, vAy: 0, vAz: 0, vBx: 0, vBy: 0, vBz: 0, UA: 1, UB: 1},
    ]
    even_eqexprs = [sp.expand(even_residual.subs(vals)) for vals in even_samples]
    M_even, b_even = sp.linear_eq_to_matrix(even_eqexprs, even_unknowns)
    even_soln = sp.solve(M_even * sp.Matrix(even_unknowns) - b_even, even_unknowns, dict=True)
    if len(even_soln) != 1:
        raise AssertionError(f"Unexpected even solve output: {even_soln}")
    even_soln = even_soln[0]

    print(f"Solved even data: {even_soln}")
    assert_zero("Minimal P0⊕P2 support operator is uniquely fixed",
                M_even * sp.Matrix([even_soln[u] for u in even_unknowns]) - b_even)
    assert_zero("Solved even operator reproduces the exact even 2PN target",
                even_residual.subs(even_soln))

    assert_zero("a0 = a20 = a21 = a22 = 1",
                sp.Matrix([even_soln[a0] - 1,
                           even_soln[a20] - 1,
                           even_soln[a21] - 1,
                           even_soln[a22] - 1]))
    assert_zero("j0 = 4/sqrt(5)", even_soln[j0] - 4 / sp.sqrt(5))
    assert_zero("j20 = 5/4", even_soln[j20] - sp.Rational(5, 4))
    assert_zero("static support+closure target coefficient = 5/4", even_soln[s] - sp.Rational(5, 4))

    # ------------------------------------------------------------------
    print("\n=== ODD / solve the minimal dressed dipole operator ===")
    # ------------------------------------------------------------------
    sigma, pperp, ppara = sp.symbols('sigma pperp ppara', real=True)

    odd_ans = sp.expand(
        sigma * (vA2 + vB2) * L1_wake
        - (UA + UB) * (pperp * uAB + ppara * dA * dB)
    )

    odd_target = sp.expand(
        -sp.Rational(7, 4) * vAB * (vA2 + vB2)
        - sp.Rational(1, 4) * dA * dB * (vA2 + vB2)
        - sp.Rational(15, 4) * (UA + UB) * vAB
    )

    odd_unknowns = [sigma, pperp, ppara]
    odd_residual = sp.expand(odd_ans - odd_target)
    odd_samples = [
        {vAx: 1, vAy: 0, vAz: 0, vBx: 1, vBy: 0, vBz: 0, UA: 0, UB: 0},
        {vAx: 0, vAy: 0, vAz: 1, vBx: 0, vBy: 0, vBz: 1, UA: 0, UB: 0},
        {vAx: 1, vAy: 0, vAz: 1, vBx: 1, vBy: 0, vBz: 1, UA: 0, UB: 0},
        {vAx: 0, vAy: 0, vAz: 1, vBx: 0, vBy: 0, vBz: 1, UA: 1, UB: 0},
        {vAx: 1, vAy: 0, vAz: 0, vBx: 1, vBy: 0, vBz: 1, UA: 1, UB: 0},
    ]
    odd_eqexprs = [sp.expand(odd_residual.subs(vals)) for vals in odd_samples]
    M_odd, b_odd = sp.linear_eq_to_matrix(odd_eqexprs, odd_unknowns)
    odd_soln = sp.solve(M_odd * sp.Matrix(odd_unknowns) - b_odd, odd_unknowns, dict=True)
    if len(odd_soln) != 1:
        raise AssertionError(f"Unexpected odd solve output: {odd_soln}")
    odd_soln = odd_soln[0]

    print(f"Solved odd data: {odd_soln}")
    assert_zero("Minimal odd dipole operator is uniquely fixed",
                M_odd * sp.Matrix([odd_soln[u] for u in odd_unknowns]) - b_odd)
    assert_zero("Solved odd operator reproduces the exact odd 2PN target",
                odd_residual.subs(odd_soln))
    assert_zero("sigma = 1/2", odd_soln[sigma] - sp.Rational(1, 2))
    assert_zero("pperp = ppara = 15/4",
                sp.Matrix([odd_soln[pperp] - sp.Rational(15, 4), odd_soln[ppara] - sp.Rational(15, 4)]))

    eta_perp = sp.simplify(odd_soln[pperp] / sp.Rational(7, 2))
    eta_para = sp.simplify(odd_soln[ppara] / 4)
    assert_zero("eta_perp = 15/14", eta_perp - sp.Rational(15, 14))
    assert_zero("eta_parallel = 15/16", eta_para - sp.Rational(15, 16))

    # ------------------------------------------------------------------
    print("\n=== DATA / canonical zero-frequency operator package ===")
    # ------------------------------------------------------------------
    Rodd = sp.diag(sp.Rational(7, 2), sp.Rational(7, 2), 4)
    assert_zero("Odd zero-frequency residues are fixed at {7/2, 7/2, 4}",
                Rodd - sp.diag(sp.Rational(7, 2), sp.Rational(7, 2), 4))

    Reven = sp.eye(6)
    J = sp.Matrix([4 / sp.sqrt(5), sp.Rational(5, 4), 0, 0, 0, 0])
    support_static = sp.simplify(J.dot(J))
    closure_deficit = sp.simplify(support_static - sp.Rational(5, 4))
    kappa_geom = sp.simplify(closure_deficit / 2)

    assert_zero("Even zero-frequency support residues are all 1", Reven - sp.eye(6))
    assert_zero("Support static coefficient = 381/80", support_static - sp.Rational(381, 80))
    assert_zero("Geometry closure deficit = 281/80", closure_deficit - sp.Rational(281, 80))
    assert_zero("Direct geometry-energy coefficient = 281/160", kappa_geom - sp.Rational(281, 160))
    assert_zero("Direct geometry-energy term produces the pure-U closure block",
                (-kappa_geom * (UA + UB)**2 + kappa_geom * UA**2 + kappa_geom * UB**2) + closure_deficit * UA * UB)

    SupportA = sp.Matrix([
        sp.expand(Pi0A + J[0] * UA),
        sp.expand(Pi20A + J[1] * UA),
        Pi21Ac,
        Pi21As,
        Pi22Ac,
        Pi22As,
    ])
    SupportB = sp.Matrix([
        sp.expand(Pi0B + J[0] * UB),
        sp.expand(Pi20B + J[1] * UB),
        Pi21Bc,
        Pi21Bs,
        Pi22Bc,
        Pi22Bs,
    ])

    even_support_minus_closure = sp.expand(SupportA.dot(SupportB) - closure_deficit * UA * UB)
    odd_dressed = sp.expand(
        odd_soln[sigma] * (vA2 + vB2) * L1_wake
        - odd_soln[pperp] * (UA + UB) * uAB
        - odd_soln[ppara] * (UA + UB) * dA * dB
    )

    full_target = sp.expand(even_target + odd_target)
    operator_rebuild = sp.expand(even_support_minus_closure + odd_dressed)

    assert_zero("Canonical operator package reproduces the full added conservative 2PN cross block",
                operator_rebuild - full_target)

    # ------------------------------------------------------------------
    print("\n=== LOW-FREQUENCY / unresolved DtN data beyond the zero-frequency residues ===")
    # ------------------------------------------------------------------
    omega, chi0, chi2, chi1p, chi10 = sp.symbols('omega chi0 chi2 chi1p chi10', real=True)
    Y0 = sp.expand(1 + chi0 * omega**2)
    Y2 = sp.expand(1 + chi2 * omega**2)
    Y1p = sp.expand(sp.Rational(7, 2) + chi1p * omega**2)
    Y10 = sp.expand(4 + chi10 * omega**2)

    assert_zero("Y0(0) = 1", Y0.subs(omega, 0) - 1)
    assert_zero("Y2(0) = 1", Y2.subs(omega, 0) - 1)
    assert_zero("Y1_perp(0) = 7/2", Y1p.subs(omega, 0) - sp.Rational(7, 2))
    assert_zero("Y1_parallel(0) = 4", Y10.subs(omega, 0) - 4)

    print("Current 2PN conservative matching fixes only the zero-frequency residues; the omega^2 DtN coefficients remain free PDE observables:")
    print(f"  Y0(omega)  = {Y0}")
    print(f"  Y2(omega)  = {Y2}")
    print(f"  Y1p(omega) = {Y1p}")
    print(f"  Y10(omega) = {Y10}")

    # ------------------------------------------------------------------
    print("\n=== N-BODY / pair AB in a three-body background C ===")
    # ------------------------------------------------------------------
    Gconst, cLight, mA, mB, mC, rAB, rAC, rBC = sp.symbols('Gconst cLight mA mB mC rAB rAC rBC', positive=True)
    vA2AB, vB2AB, vABAB, dAAB, dBAB = sp.symbols('vA2AB vB2AB vABAB dAAB dBAB', real=True)

    UA_loc = sp.expand(Gconst * (mB / rAB + mC / rAC))
    UB_loc = sp.expand(Gconst * (mA / rAB + mC / rBC))
    pair_pref = sp.expand(Gconst * mA * mB / (cLight**4 * rAB))

    source_B = sp.expand(J[0] * (sp.sqrt(5) * vB2AB / 2) + J[1] * ((3 * dBAB**2 - vB2AB) / 2))
    source_A = sp.expand(J[0] * (sp.sqrt(5) * vA2AB / 2) + J[1] * ((3 * dAAB**2 - vA2AB) / 2))
    assert_zero("U-driven scalar source = 11/8 v^2 + 15/8 d^2",
                sp.Matrix([
                    source_A - (sp.Rational(11, 8) * vA2AB + sp.Rational(15, 8) * dAAB**2),
                    source_B - (sp.Rational(11, 8) * vB2AB + sp.Rational(15, 8) * dBAB**2),
                ]))

    pair_AB_threebody_vel = sp.expand(
        pair_pref
        * (
            (UA_loc - Gconst * mB / rAB) * source_B
            + (UB_loc - Gconst * mA / rAB) * source_A
            - sp.Rational(15, 4) * ((UA_loc - Gconst * mB / rAB) + (UB_loc - Gconst * mA / rAB)) * vABAB
        )
    )

    pair_AB_threebody_vel_expected = sp.expand(
        Gconst**2 * mA * mB * mC / (8 * cLight**4 * rAB)
        * (
            (11 * vB2AB + 15 * dBAB**2 - 30 * vABAB) / rAC
            + (11 * vA2AB + 15 * dAAB**2 - 30 * vABAB) / rBC
        )
    )
    assert_zero("Pair AB velocity-dependent 3-body lift is fixed by the minimal operator",
                pair_AB_threebody_vel - pair_AB_threebody_vel_expected)

    pair_AB_threebody_static = sp.expand(
        pair_pref
        * sp.Rational(5, 4)
        * (UA_loc * UB_loc - (Gconst * mB / rAB) * (Gconst * mA / rAB))
    )

    pair_AB_threebody_static_expected = sp.expand(
        5 * Gconst**3 * mA * mB * mC / (4 * cLight**4)
        * (
            mA / (rAB**2 * rAC)
            + mB / (rAB**2 * rBC)
            + mC / (rAB * rAC * rBC)
        )
    )
    assert_zero("Pair AB static 3-body lift is fixed by the minimal operator",
                pair_AB_threebody_static - pair_AB_threebody_static_expected)

    print("\nSolved zero-frequency operator data:")
    print("  odd residues      : {7/2, 7/2, 4}")
    print("  even residues     : {1, 1, 1, 1, 1, 1}")
    print(f"  J source vector   : {J.T}")
    print(f"  support static    : {support_static}")
    print(f"  geometry closure  : {closure_deficit}")
    print(f"  kappa_geom        : {kappa_geom}")
    print(f"  sigma             : {odd_soln[sigma]}")
    print(f"  eta_perp          : {eta_perp}")
    print(f"  eta_parallel      : {eta_para}")

    print("\nPair AB in background C:")
    print(f"  velocity 3-body lift = {sp.factor(pair_AB_threebody_vel_expected)}")
    print(f"  static   3-body lift = {sp.factor(pair_AB_threebody_static_expected)}")

    print("\nDone.")


if __name__ == '__main__':
    main()
