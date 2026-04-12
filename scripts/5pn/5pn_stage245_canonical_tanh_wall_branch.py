#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 245 — canonical tanh-wall branch and natural local mouth closure.
"""


def main() -> None:
    banner("STAGE 245 — CANONICAL TANH-WALL BRANCH")

    xi = sp.symbols("xi", real=True)
    hbar, m = sp.symbols("hbar m", positive=True, real=True)
    rho_w, csw = sp.symbols("rho_w csw", positive=True, real=True)
    a, ell, L, V0 = sp.symbols("a ell L V0", positive=True, real=True)

    subbanner("I. Exact tanh-wall moments")
    f = sp.simplify((1 + sp.tanh(xi)) / 2)
    fp = sp.simplify(sp.diff(f, xi))
    fpp = sp.simplify(sp.diff(fp, xi))

    If = sp.integrate(sp.simplify(fp**2), (xi, -sp.oo, sp.oo))
    t = sp.symbols("t", real=True)
    Ig = sp.integrate(t**2 * (1 - t**2), (t, -1, 1))  # with t = tanh(xi)

    print("f'(xi) =")
    sp.pprint(fp)
    print()
    print("f''(xi) =")
    sp.pprint(fpp)
    print()
    print("I_f =")
    sp.pprint(If)
    print("I_g =")
    sp.pprint(Ig)
    expect_zero("I_f = 1/3", sp.simplify(If - sp.Rational(1, 3)))
    expect_zero("I_g = 4/15", sp.simplify(Ig - sp.Rational(4, 15)))
    expect_zero("I_g/I_f = 4/5", sp.simplify(Ig / If - sp.Rational(4, 5)))

    subbanner("II. Explicit branch coefficients")
    Hw = sp.simplify(m * csw**2 / rho_w)
    TX = sp.simplify(sp.pi * a**2 * ell * If * hbar**2 / (m * rho_w))
    KX = sp.simplify(4 * sp.pi * a**2 * ell * If * Hw + sp.pi * a**2 * Ig * hbar**2 / (m * rho_w * ell))
    J1 = sp.simplify(If / Hw)
    kappa = sp.simplify(KX * L**2 / TX)
    W_wall = sp.simplify(4 * rho_w**2 * V0**2 * L**2 / (hbar**2 * csw**2 * ell**2))

    print("T_X =")
    sp.pprint(TX)
    print()
    print("K_X =")
    sp.pprint(KX)
    print()
    print("J_1 =")
    sp.pprint(J1)
    print()
    print("kappa =")
    sp.pprint(kappa)
    expect_zero(
        "canonical tanh kappa identity",
        sp.simplify(kappa - (4 * (m * csw * L / hbar)**2 + sp.Rational(4, 5) * (L / ell)**2)),
    )

    subbanner("III. Natural local mouth closure")
    Km = sp.simplify(TX / ell)
    eta = sp.simplify(Km * L / TX)
    expect_zero("eta = L/ell", sp.simplify(eta - L / ell))
    print("eta =")
    sp.pprint(eta)

    subbanner("IV. Three branch control parameters")
    chi_s = sp.symbols("chi_s", positive=True, real=True)
    Lambda_ell = sp.symbols("Lambda_ell", positive=True, real=True)
    Upsilon_w = sp.symbols("Upsilon_w", positive=True, real=True)

    kappa_ctrl = sp.simplify(4 * chi_s**2 + sp.Rational(4, 5) * Lambda_ell**2)
    eta_ctrl = sp.simplify(Lambda_ell)
    W_ctrl = sp.simplify(Upsilon_w * Lambda_ell**2)

    expect_zero(
        "control-map kappa identity",
        sp.simplify(
            kappa.subs({m * csw * L / hbar: chi_s, L / ell: Lambda_ell}) - kappa_ctrl
        ),
    )
    expect_zero("control-map eta identity", sp.simplify(eta.subs({L / ell: Lambda_ell}) - eta_ctrl))
    expect_zero(
        "control-map W_wall identity",
        sp.simplify(
            W_wall.subs({4 * rho_w**2 * V0**2 / (hbar**2 * csw**2): Upsilon_w, L / ell: Lambda_ell}) - W_ctrl
        ),
    )

    banner("STAGE 245 LEDGER")
    print("1. The canonical tanh wall has exact moments")
    print("      I_f = 1/3,   I_g = 4/15,   I_g/I_f = 4/5.")
    print("2. The corresponding thin-shell support branch obeys")
    print("      kappa = 4 (m c_(s,w) L / hbar)^2 + (4/5) (L/ell)^2,")
    print("      eta   = L/ell,")
    print("      W_wall = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2).")
    print("3. So the first explicit parent branch collapses to three dimensionless controls:")
    print("      chi_s = m c_(s,w) L / hbar,")
    print("      Lambda_ell = L/ell,")
    print("      Upsilon_w = 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2).")
    print("4. In those variables,")
    print("      kappa = 4 chi_s^2 + (4/5) Lambda_ell^2,")
    print("      eta = Lambda_ell,")
    print("      W_wall = Upsilon_w Lambda_ell^2.")


if __name__ == "__main__":
    main()
