#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 234 — exact support-compensation theorem on the coherent tracking branch.
"""


def main() -> None:
    banner("STAGE 234 — TRACKING SUPPORT-COMPENSATION THEOREM")

    xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)
    eps = sp.symbols("eps", positive=True, real=True)
    zeta = sp.symbols("zeta", positive=True, real=True)
    M_mix = sp.symbols("M_mix", positive=True, real=True)
    xi_req = sp.symbols("xi_req", positive=True, real=True)
    S_req = sp.symbols("S_req", positive=True, real=True)

    G_tr = sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R**2) * xi))
    F_tr = sp.simplify(
        ((9 * delta + (9 + 2 * R**2) * xi) ** 2 * (9 * delta + (9 + 2 * R) * xi) ** 2)
        / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + (9 + 2 * R**2) * xi**2) ** 2)
    )
    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    dS_dzeta = sp.simplify(sp.diff(S, zeta))

    subbanner("I. Exact tracking-branch functions")
    print("G_tr(xi,delta;R) =")
    sp.pprint(G_tr)
    print()
    print("F_tr(xi,delta;R) =")
    sp.pprint(F_tr)

    subbanner("II. Exact support-enhancement factor")
    print("S(zeta;eps) =")
    sp.pprint(S)
    print()
    print("dS/dzeta =")
    sp.pprint(dS_dzeta)
    expect_zero(
        "exact monotonicity identity",
        sp.simplify(dS_dzeta - (1 - eps) / (1 - zeta * eps) ** 2),
    )

    zeta_req_expr = sp.simplify(sp.solve(sp.Eq(S_req, S), zeta)[0])
    print()
    print("zeta_req(S_req;eps) =")
    sp.pprint(zeta_req_expr)
    expect_zero("roundtrip S(zeta_req)-S_req", sp.simplify(S.subs(zeta, zeta_req_expr) - S_req))

    subbanner("III. Stable-side branch point and required support ratio")
    M_req = sp.simplify(G_tr.subs(xi, xi_req))
    M_crit = sp.simplify(G_tr.subs(xi, 1))
    dG_dxi = sp.factor(sp.simplify(sp.diff(G_tr, xi)))
    S_req_branch = sp.simplify(M_req / M_mix)
    zeta_req_branch = sp.simplify(zeta_req_expr.subs(S_req, S_req_branch))

    print("M_req = G_tr(xi_req,delta;R) =")
    sp.pprint(M_req)
    print()
    print("M_crit = G_tr(1,delta;R) =")
    sp.pprint(M_crit)
    print()
    print("dG_tr/dxi =")
    sp.pprint(dG_dxi)
    expect_zero(
        "exact positivity numerator for dG_tr/dxi",
        sp.simplify(
            dG_dxi
            - 9 * (2 * R**2 * xi**2 + 9 * delta**2 + 18 * delta * xi + 9 * xi**2)
            / (9 * delta + (9 + 2 * R**2) * xi) ** 2
        ),
    )
    print()
    print("S_req = M_req/M_mix =")
    sp.pprint(S_req_branch)
    print()
    print("zeta_req(branch) =")
    sp.pprint(zeta_req_branch)

    subbanner("IV. Exact no-go elimination on the reduced branch")
    zeta_gap = sp.simplify((1 / eps) - zeta_req_expr)
    print("1/eps - zeta_req =")
    sp.pprint(zeta_gap)
    expect_zero(
        "exact softening-gap identity",
        sp.simplify(zeta_gap - (1 - eps) / (eps * (1 + eps * (S_req - 2)))),
    )
    print()
    print("F(0,delta;R) =", sp.simplify(F_tr.subs(xi, 0)))
    print("lim_{xi->1^-} F_tr =", sp.limit(F_tr, xi, 1, dir="-"))

    banner("STAGE 234 LEDGER")
    print("1. The coherent tracking branch carries support only through")
    print("      S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps).")
    print("2. S is strictly increasing and exactly invertible:")
    print("      zeta_req = (S_req - 1)/(1 + eps(S_req - 2)).")
    print("3. Once a stable-side point xi_req solves F_tr(xi_req,delta;R) = R_target, the")
    print("   required load is M_req = G_tr(xi_req,delta;R), so the support threshold is")
    print("      zeta_req = (M_req/M_mix - 1)/(1 + eps(M_req/M_mix - 2)).")
    print("4. Because 0 < zeta_req < 1/eps for every finite S_req > 1 on 0 < eps < 1, there")
    print("   is no reduced-level support no-go on the coherent tracking branch.")


if __name__ == "__main__":
    main()
