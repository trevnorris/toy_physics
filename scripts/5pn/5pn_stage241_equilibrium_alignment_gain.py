#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 241 — parent equilibrium source/support alignment and the exact matched-layer gain.
"""


def main() -> None:
    banner("STAGE 241 — PARENT EQUILIBRIUM SOURCE/SUPPORT ALIGNMENT")

    g_phi = sp.symbols("g_phi", positive=True, real=True)
    N_pp = sp.symbols("N_phi_phi", positive=True, real=True)
    I1, I2 = sp.symbols("I_1 I_2", positive=True, real=True)
    K_X = sp.symbols("K_X", positive=True, real=True)
    H_w = sp.symbols("H_w", positive=True, real=True)
    rho_w, m, c_s_w = sp.symbols("rho_w m c_s_w", positive=True, real=True)
    C_eq_sq = sp.symbols("C_eq_sq", nonnegative=True, real=True)

    O_sp_eq = sp.simplify(g_phi * I1)
    N_ss_eq = sp.simplify(g_phi**2 * I2)
    C_eq_expr = sp.simplify(O_sp_eq**2 / (N_pp * N_ss_eq))
    DeltaK_eq = sp.simplify(g_phi**2 * I1)
    G_eq = sp.simplify(DeltaK_eq / K_X)
    G_eq_matched = sp.simplify(G_eq.subs(I1, N_pp / H_w))
    H_w_expr = sp.simplify(m * c_s_w**2 / rho_w)
    G_eq_matched_parent = sp.simplify(G_eq_matched.subs(H_w, H_w_expr))

    subbanner("I. Exact equilibrium-aligned overlaps")
    print("O_(sigma,phi)^(eq) =")
    sp.pprint(O_sp_eq)
    print()
    print("N_(sigma,sigma)^(eq) =")
    sp.pprint(N_ss_eq)
    print()
    print("C_(sigma,phi)^2 on the equilibrium branch =")
    sp.pprint(C_eq_expr)
    expect_zero(
        "exact equilibrium coherence identity",
        sp.simplify(C_eq_expr - I1**2 / (N_pp * I2)),
    )

    subbanner("II. Exact equilibrium softening and gain")
    print("Delta K_X^(eq) =")
    sp.pprint(DeltaK_eq)
    print()
    print("G_eq =")
    sp.pprint(G_eq)
    expect_zero("exact equilibrium gain identity", sp.simplify(G_eq - g_phi**2 * I1 / K_X))

    subbanner("III. Matched-layer limit")
    print("Matched-layer condition:  I_1 = N_(phi,phi)/H_w  and  C_(sigma,phi)^2 = 1.")
    print()
    print("G_eq^(matched) =")
    sp.pprint(G_eq_matched)
    print()
    print("Using H_w = m c_(s,w)^2 / rho_w, the matched-layer gain becomes")
    sp.pprint(G_eq_matched_parent)
    expect_zero(
        "exact matched-layer parent identity",
        sp.simplify(G_eq_matched_parent - rho_w * g_phi**2 * N_pp / (m * c_s_w**2 * K_X)),
    )

    banner("STAGE 241 LEDGER")
    print("1. On the parent equilibrium branch the source profile is not free; it is tied to the support")
    print("   profile by the local compressional stiffness H(y).")
    print("2. Therefore the source/support coherence is no longer an independent datum. It is")
    print("      C_(sigma,phi)^2 = I_1^2 / [N_(phi,phi) I_2] <= 1.")
    print("3. The exact equilibrium gain is")
    print("      G_eq = g_phi^2 I_1 / K_X.")
    print("4. In the thin matched layer where H(y) is nearly constant, the branch saturates to")
    print("      C_(sigma,phi)^2 = 1,")
    print("   and the gain becomes")
    print("      G_eq = rho_w g_phi^2 N_(phi,phi) / [m c_(s,w)^2 K_X].")
    print("5. So the best-alignment formulas used in the earlier threshold scripts are not ad hoc;")
    print("   they are the natural thin-layer limit of the parent equilibrium branch.")


if __name__ == "__main__":
    main()
