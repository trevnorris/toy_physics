#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 244 — explicit GNLS wall-shell reduction for the first support branch.
This script derives T_X, K_X, kappa, and verifies W_wall = Xi on the matched thin-wall branch.
"""


def main() -> None:
    banner("STAGE 244 — GNLS WALL-SHELL REDUCTION")

    hbar, m = sp.symbols("hbar m", positive=True, real=True)
    rho_w, csw = sp.symbols("rho_w csw", positive=True, real=True)
    a, ell, L, V0 = sp.symbols("a ell L V0", positive=True, real=True)
    If, Ig = sp.symbols("I_f I_g", positive=True, real=True)

    Hw = sp.simplify(m * csw**2 / rho_w)
    Nphiphi = sp.simplify(4 * sp.pi * a**2 * ell * If)
    Gphiphi = sp.simplify(4 * sp.pi * a**2 * Ig / ell)

    subbanner("I. Parent quadratic shell coefficients")
    TX = sp.simplify(hbar**2 * Nphiphi / (4 * m * rho_w))
    KX = sp.simplify(Hw * Nphiphi + hbar**2 * Gphiphi / (4 * m * rho_w))

    print("T_X =")
    sp.pprint(TX)
    print()
    print("K_X =")
    sp.pprint(KX)

    expect_zero(
        "T_X explicit shell identity",
        sp.simplify(TX - sp.pi * a**2 * ell * If * hbar**2 / (m * rho_w)),
    )
    expect_zero(
        "K_X explicit shell identity",
        sp.simplify(KX - (4 * sp.pi * a**2 * ell * If * Hw + sp.pi * a**2 * Ig * hbar**2 / (m * rho_w * ell))),
    )

    subbanner("II. Exact support geometry parameter kappa")
    kappa = sp.simplify(KX * L**2 / TX)
    kappa_expected = sp.simplify(4 * (m * csw * L / hbar)**2 + (Ig / If) * (L / ell)**2)
    print("kappa =")
    sp.pprint(kappa)
    expect_zero("kappa identity", sp.simplify(kappa - kappa_expected))

    subbanner("III. Exact matched-layer gain and wall figure")
    J1 = sp.simplify(If / Hw)
    W_wall = sp.simplify(4 * sp.pi * a**2 * L**2 * J1 * V0**2 / (TX * ell))
    W_expected = sp.simplify(4 * rho_w**2 * V0**2 * L**2 / (hbar**2 * csw**2 * ell**2))
    print("W_wall =")
    sp.pprint(W_wall)
    expect_zero("W_wall collapse", sp.simplify(W_wall - W_expected))

    gphi = sp.simplify(V0 / ell)
    I1 = sp.simplify(Nphiphi / Hw)
    Xi = sp.simplify(gphi**2 * I1 * L**2 / TX)
    print()
    print("Xi =")
    sp.pprint(Xi)
    expect_zero("W_wall = Xi", sp.simplify(W_wall - Xi))

    banner("STAGE 244 LEDGER")
    print("1. The first explicit shell reduction fixes")
    print("      T_X = (hbar^2 / (4 m rho_w)) N_(phi phi),")
    print("      K_X = H_w N_(phi phi) + (hbar^2 / (4 m rho_w)) G_(phi phi).")
    print("2. On the thin wall shell with chi_phi = f'(xi), these become")
    print("      T_X = pi a^2 ell I_f hbar^2 / (m rho_w),")
    print("      K_X = 4 pi a^2 ell I_f H_w + pi a^2 I_g hbar^2 / (m rho_w ell).")
    print("3. The geometry parameter is")
    print("      kappa = 4 (m c_(s,w) L / hbar)^2 + (I_g/I_f) (L/ell)^2.")
    print("4. The wall figure of merit collapses exactly to")
    print("      W_wall = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2).")
    print("5. On this matched thin-wall branch, the Stage-41/42 support/source coupling is the same object:")
    print("      Xi = W_wall.")


if __name__ == "__main__":
    main()
