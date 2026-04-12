#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 240 — exact parent-threshold decision theorem after folding the Stage-236 twin
criterion into the Stage-237/238/239 physical non-twin placement and microscopic-gain
phase diagram.
"""


def main() -> None:
    banner("STAGE 240 — PARENT THRESHOLD DECISION THEOREM")

    Pe_req = sp.symbols("Pe_req", positive=True, real=True)
    kappa = sp.symbols("kappa", positive=True, real=True)
    eta = sp.symbols("eta", positive=True, real=True)
    zeta_req = sp.symbols("zeta_req", positive=True, real=True)
    rho_star, g_phi = sp.symbols("rho_star g_phi", positive=True, real=True)
    N_ss, N_pp = sp.symbols("N_sigma_sigma N_phi_phi", positive=True, real=True)
    O_sp = sp.symbols("O_sigma_phi", real=True)
    C_sp_sq = sp.symbols("C_sigma_phi_sq", nonnegative=True, real=True)
    m, c_s_star = sp.symbols("m c_s_star", positive=True, real=True)
    K_X, T_X, L = sp.symbols("K_X T_X L", positive=True, real=True)
    y = sp.symbols("y", positive=True, real=True)
    pi = sp.pi

    alpha = sp.sqrt(kappa)
    Delta0 = sp.simplify(eta * (sp.cosh(alpha) - 1) / (kappa * (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))))
    Delta_inf = sp.simplify((sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) / (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha)))

    zeta_max = sp.simplify((pi**2 / 4) * (kappa + pi**2 / 4) / kappa)
    G_fail = sp.simplify(Pe_req / (kappa * Delta_inf))
    G_suff = sp.simplify(Pe_req / (kappa * Delta0))

    gphi_fail_sq = sp.simplify(m * c_s_star**2 * K_X * N_ss * G_fail / (rho_star * O_sp**2))
    gphi_suff_sq = sp.simplify(m * c_s_star**2 * K_X * N_ss * G_suff / (rho_star * O_sp**2))
    C_fail_sq = sp.simplify(m * c_s_star**2 * K_X * G_fail / (rho_star * g_phi**2 * N_pp))
    C_suff_sq = sp.simplify(m * c_s_star**2 * K_X * G_suff / (rho_star * g_phi**2 * N_pp))

    gphi_fail_sq_TX = sp.simplify(gphi_fail_sq.subs(K_X, kappa * T_X / L**2))
    gphi_suff_sq_TX = sp.simplify(gphi_suff_sq.subs(K_X, kappa * T_X / L**2))
    G_max = sp.simplify(rho_star * g_phi**2 * N_pp / (m * c_s_star**2 * K_X))

    subbanner("I. Exact microscopic phase diagram")
    print("G_fail(kappa,eta) =")
    sp.pprint(G_fail)
    print()
    print("G_suff(kappa,eta) =")
    sp.pprint(G_suff)
    print()
    print("G_suff - G_fail =")
    sp.pprint(sp.simplify(G_suff - G_fail))

    subbanner("II. Parent thresholds on g_phi and C_(sigma,phi)^2")
    print("g_(phi,fail)^2 =")
    sp.pprint(gphi_fail_sq)
    print()
    print("g_(phi,suff)^2 =")
    sp.pprint(gphi_suff_sq)
    print()
    print("C_fail^2 =")
    sp.pprint(C_fail_sq)
    print()
    print("C_suff^2 =")
    sp.pprint(C_suff_sq)
    print()
    print("Best-case gain at fixed g_phi:  G_max =")
    sp.pprint(G_max)
    print("(attained only at perfect alignment C_(sigma,phi)^2 = 1)")

    expect_zero(
        "exact amplitude-threshold identity (fail)",
        sp.simplify(gphi_fail_sq_TX - m * c_s_star**2 * T_X * N_ss * Pe_req / (rho_star * L**2 * O_sp**2 * Delta_inf)),
    )
    expect_zero(
        "exact amplitude-threshold identity (suff)",
        sp.simplify(gphi_suff_sq_TX - m * c_s_star**2 * T_X * N_ss * Pe_req / (rho_star * L**2 * O_sp**2 * Delta0)),
    )

    subbanner("III. Twin vs non-twin exact decision split")
    print("zeta_max(kappa) =")
    sp.pprint(zeta_max)
    expect_zero(
        "exact zeta_max identity",
        sp.simplify(zeta_max - (pi**2 / 4) * (kappa + pi**2 / 4) / kappa),
    )
    print()
    print("Decision tree:")
    print("  zeta_req <= 1                : symmetric lowest twin already sufficient")
    print("  1 < zeta_req <= zeta_max     : first explicit non-twin family can in principle rescue")
    print("  zeta_req > zeta_max          : first explicit non-twin family impossible")

    kappa_max = sp.simplify(sp.solve(sp.Eq(zeta_req, zeta_max), kappa)[0])
    print()
    print("Equivalent stiffness ceiling when zeta_req > pi^2/4:")
    print("kappa <=")
    sp.pprint(kappa_max)
    expect_zero(
        "exact kappa_max identity",
        sp.simplify(kappa_max - pi**4 / (4 * (4 * zeta_req - pi**2))),
    )

    banner("STAGE 240 LEDGER")
    print("1. Stage 236 still gives the sharp symmetric-twin theorem:")
    print("      zeta_req <= 1  iff  the lowest symmetric twin lane is enough.")
    print("2. Once zeta_req > 1, the first explicit physical non-twin family is controlled by")
    print("      zeta_max(kappa) = (pi^2/4)(kappa + pi^2/4)/kappa.")
    print("3. Inside that reachability window, the operator-selected branch has the exact microscopic")
    print("   phase diagram")
    print("      G_micro <= G_fail  : fail,   G_micro >= G_suff : success,")
    print("   with only the narrow band G_fail < G_micro < G_suff requiring the full root solve.")
    print("4. In parent variables this becomes exact threshold surfaces on the confinement-loading")
    print("   amplitude g_phi and the source/support coherence factor C_(sigma,phi)^2.")
    print("5. So the remaining support theorem gap is now fully ordered:")
    print("      first ask whether the twin branch is enough; if not, ask whether the first physical")
    print("      non-twin family is even reachable; only then check whether the parent branch crosses")
    print("      the explicit microscopic gain thresholds.")


if __name__ == "__main__":
    main()
