import sympy as sp

# -----------------------------------------------------------------------------
# 2PN support-minus-closure throat EFT prototype (SymPy)
# -----------------------------------------------------------------------------

# Symbols
vAx, vAy, vAz, vBx, vBy, vBz = sp.symbols('vAx vAy vAz vBx vBy vBz', real=True)
Gconst, cLight, mA, mB, rAB = sp.symbols('Gconst cLight mA mB rAB', positive=True)


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
    print("--- 2PN support-minus-closure throat EFT prototype ---")

    # ------------------------------------------------------------------
    print("\n=== MASTER / component kinematics and solved targets ===")
    # ------------------------------------------------------------------
    vA2 = sp.expand(vAx**2 + vAy**2 + vAz**2)
    vB2 = sp.expand(vBx**2 + vBy**2 + vBz**2)
    vAB = sp.expand(vAx * vBx + vAy * vBy + vAz * vBz)

    dA = vAz
    dB = vBz

    UA = sp.expand(Gconst * mB / rAB)
    UB = sp.expand(Gconst * mA / rAB)

    L1_wake_target = sp.expand(-sp.Rational(7, 2) * vAB - sp.Rational(1, 2) * dA * dB)

    quartic_target = sp.expand(
        -sp.Rational(7, 4) * vAB * (vA2 + vB2)
        - sp.Rational(1, 4) * dA * dB * (vA2 + vB2)
        + sp.Rational(11, 8) * vA2 * vB2
        + sp.Rational(1, 4) * vAB**2
        - sp.Rational(5, 8) * (vA2 * dB**2 + vB2 * dA**2)
        + sp.Rational(3, 2) * vAB * dA * dB
        + sp.Rational(3, 8) * dA**2 * dB**2
    )

    quadratic_target = sp.expand(
        sp.Rational(11, 8) * (mA * vA2 + mB * vB2)
        - sp.Rational(15, 4) * (mA + mB) * vAB
        + sp.Rational(15, 8) * (mA * dA**2 + mB * dB**2)
    )

    added_cross_target = sp.expand(
        (Gconst * mA * mB) / (cLight**4 * rAB) * quartic_target
        + (Gconst**2 * mA * mB) / (cLight**4 * rAB**2) * quadratic_target
        + sp.Rational(5, 4) * Gconst**3 * mA**2 * mB**2 / (cLight**4 * rAB**3)
    )

    # ------------------------------------------------------------------
    print("\n=== 1PN / frozen wake as odd dipole channels ===")
    # ------------------------------------------------------------------
    DipA = sp.Matrix([sp.sqrt(sp.Rational(7, 2)) * vAx,
                      sp.sqrt(sp.Rational(7, 2)) * vAy,
                      2 * dA])
    DipB = sp.Matrix([sp.sqrt(sp.Rational(7, 2)) * vBx,
                      sp.sqrt(sp.Rational(7, 2)) * vBy,
                      2 * dB])
    L1_wake_dipole = sp.expand(-(DipA.dot(DipB)))
    assert_zero("Frozen 1PN wake is exactly an odd dipole overlap", L1_wake_dipole - L1_wake_target)

    leg_dressed_odd = sp.expand(sp.Rational(1, 2) * (vA2 + vB2) * L1_wake_dipole)

    # ------------------------------------------------------------------
    print("\n=== 2PN / even P0⊕P2 block in mouth-port language ===")
    # ------------------------------------------------------------------
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

    quartic_residual = sp.expand(quartic_target - leg_dressed_odd)
    port_quartic = sp.expand(
        Pi0A * Pi0B
        + Pi20A * Pi20B
        + Pi21Ac * Pi21Bc
        + Pi21As * Pi21Bs
        + Pi22Ac * Pi22Bc
        + Pi22As * Pi22Bs
    )
    assert_zero("New quartic residual is exactly the positive P0⊕P2 port overlap", port_quartic - quartic_residual)

    J0 = sp.simplify(4 / sp.sqrt(5))
    J20 = sp.Rational(5, 4)
    static_cross = sp.Rational(5, 4)

    eta_perp = sp.Rational(15, 14)
    eta_para = sp.Rational(15, 16)
    odd_added_target = sp.expand(
        leg_dressed_odd - sp.Rational(15, 4) * (UA + UB) * vAB
    )

    even_block_target = sp.expand(
        port_quartic
        + UA * (J0 * Pi0B + J20 * Pi20B)
        + UB * (J0 * Pi0A + J20 * Pi20A)
        + static_cross * UA * UB
    )

    # ------------------------------------------------------------------
    print("\n=== ODD / kinetic and potential dressing of the dipole wake ===")
    # ------------------------------------------------------------------
    lam = sp.symbols('lam', real=True)
    DipA_lam = sp.Matrix([
        (1 + lam * (sp.Rational(1, 2) * vA2 + eta_perp * UA)) * sp.sqrt(sp.Rational(7, 2)) * vAx,
        (1 + lam * (sp.Rational(1, 2) * vA2 + eta_perp * UA)) * sp.sqrt(sp.Rational(7, 2)) * vAy,
        (1 + lam * (sp.Rational(1, 2) * vA2 + eta_para * UA)) * 2 * dA,
    ])
    DipB_lam = sp.Matrix([
        (1 + lam * (sp.Rational(1, 2) * vB2 + eta_perp * UB)) * sp.sqrt(sp.Rational(7, 2)) * vBx,
        (1 + lam * (sp.Rational(1, 2) * vB2 + eta_perp * UB)) * sp.sqrt(sp.Rational(7, 2)) * vBy,
        (1 + lam * (sp.Rational(1, 2) * vB2 + eta_para * UB)) * 2 * dB,
    ])
    odd_cross_lam = sp.expand(-(DipA_lam.dot(DipB_lam)))
    odd_added_from_dressing = sp.expand(sp.series(odd_cross_lam, lam, 0, 2).removeO().coeff(lam, 1))
    assert_zero("Odd dipole kinetic+potential dressing reproduces the full added odd 2PN block",
                odd_added_from_dressing - odd_added_target)

    # ------------------------------------------------------------------
    print("\n=== MATRIX / exact support-minus-closure decomposition ===")
    # ------------------------------------------------------------------
    Xsymbols = sp.Matrix(sp.symbols('Pi0 Pi20 Pi21c Pi21s Pi22c Pi22s U', real=True))
    Meven = sp.Matrix([
        [1, 0, 0, 0, 0, 0, J0],
        [0, 1, 0, 0, 0, 0, J20],
        [0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 0],
        [J0, J20, 0, 0, 0, 0, static_cross],
    ])

    closure_deficit = sp.simplify(J0**2 + J20**2 - static_cross)
    assert_zero("Closure deficit is exactly 281/80", closure_deficit - sp.Rational(281, 80))

    Rsupport = sp.Matrix([
        [1, 0, 0, 0, 0, 0, J0],
        [0, 1, 0, 0, 0, 0, J20],
        [0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 0],
    ])
    Rgeom = sp.Matrix([[0, 0, 0, 0, 0, 0, sp.sqrt(closure_deficit)]])

    support_matrix = sp.expand(Rsupport.T * Rsupport)
    geom_matrix = sp.expand(Rgeom.T * Rgeom)
    assert_zero("Exact even response matrix = support PSD block minus pure-U geometry block",
                (support_matrix - geom_matrix) - Meven)

    det_scalar = sp.factor(Meven.extract([0, 1, 6], [0, 1, 6]).det())
    assert_zero("The scalar {P0,P20,U} block has determinant -281/80", det_scalar + sp.Rational(281, 80))

    # ------------------------------------------------------------------
    print("\n=== CHANNELS / support and closure sources ===")
    # ------------------------------------------------------------------
    SupportA = sp.Matrix([
        sp.expand(Pi0A + J0 * UA),
        sp.expand(Pi20A + J20 * UA),
        Pi21Ac,
        Pi21As,
        Pi22Ac,
        Pi22As,
    ])
    SupportB = sp.Matrix([
        sp.expand(Pi0B + J0 * UB),
        sp.expand(Pi20B + J20 * UB),
        Pi21Bc,
        Pi21Bs,
        Pi22Bc,
        Pi22Bs,
    ])

    ClosureA = sp.expand(sp.sqrt(closure_deficit) * UA)
    ClosureB = sp.expand(sp.sqrt(closure_deficit) * UB)

    even_support_minus_closure = sp.expand(SupportA.dot(SupportB) - ClosureA * ClosureB)
    assert_zero("Even block = six positive support channels minus one pure-U closure channel",
                even_support_minus_closure - even_block_target)

    # ------------------------------------------------------------------
    print("\n=== AUXILIARY / generic low-frequency elimination identities ===")
    # ------------------------------------------------------------------
    q, SA, SB = sp.symbols('q SA SB', real=True)
    Lsup = -sp.Rational(1, 2) * q**2 + q * (SA + SB)
    qstar = sp.solve(sp.Eq(sp.diff(Lsup, q), 0), q)[0]
    Lsup_on = sp.expand(Lsup.subs(q, qstar))
    assert_zero("Positive support auxiliary produces +SA*SB cross term",
                Lsup_on - sp.expand(sp.Rational(1, 2) * SA**2 + SA * SB + sp.Rational(1, 2) * SB**2))

    g, CA, CB = sp.symbols('g CA CB', real=True)
    Lgeom = sp.Rational(1, 2) * g**2 + g * (CA + CB)
    gstar = sp.solve(sp.Eq(sp.diff(Lgeom, g), 0), g)[0]
    Lgeom_on = sp.expand(Lgeom.subs(g, gstar))
    assert_zero("Wrong-sign geometry auxiliary produces -CA*CB cross term",
                Lgeom_on - sp.expand(-sp.Rational(1, 2) * CA**2 - CA * CB - sp.Rational(1, 2) * CB**2))

    d, DA, DB = sp.symbols('d DA DB', real=True)
    Lodd = -sp.Rational(1, 2) * d**2 + d * (DA - DB)
    dstar = sp.solve(sp.Eq(sp.diff(Lodd, d), 0), d)[0]
    Lodd_on = sp.expand(Lodd.subs(d, dstar))
    assert_zero("Odd dipole auxiliary produces -DA*DB cross term",
                Lodd_on - sp.expand(sp.Rational(1, 2) * DA**2 - DA * DB + sp.Rational(1, 2) * DB**2))

    # ------------------------------------------------------------------
    print("\n=== ASSEMBLY / full added conservative 2PN cross block ===")
    # ------------------------------------------------------------------
    added_cross_support_closure = sp.expand(
        (Gconst * mA * mB) / (cLight**4 * rAB)
        * (odd_added_from_dressing + even_support_minus_closure)
    )

    assert_zero("Full added conservative 2PN cross block equals dressed odd dipole wake plus support-minus-closure even sector",
                added_cross_support_closure - added_cross_target)

    print("\nSupport channels:")
    print("  S0  = Pi0 + (4/sqrt(5)) U")
    print("  S20 = Pi20 + (5/4) U")
    print("  S21c, S21s, S22c, S22s = pure P2 components")
    print(f"Closure channel: Cgeom = sqrt({sp.Rational(281,80)}) U")

    print("\nDone.")


if __name__ == '__main__':
    main()
