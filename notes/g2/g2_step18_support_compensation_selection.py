#!/usr/bin/env python3
"""
Step 18 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Imports the exact coherent support-enhancement law from the moving-throat
   notes:
       M_tr = M_mix S(zeta;eps),
       S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps).
2. Imports the exact physical tracking-branch load law
       M_tr = G_tr(xi,delta;R_tr)
   and its critical load M_crit.
3. Verifies the exact inverse of the support-enhancement factor and the strict
   stability margin below the blocking pole zeta = 1/eps.
4. Verifies the exact Stage-31 support-compensation theorem:
       if M_mix < M_req,
       zeta_req = (S_req - 1) / [1 + eps (S_req - 2)]
   hits the target before softening, with zeta_req < zeta_crit < 1/eps.
5. Differentiates the implicit selected-branch law and verifies
       dxi_phys/dzeta > 0.
6. Interprets the result back in Step-16 language:
   coherent local support continuation changes the available baseline at fixed
   demand ratio, so the physical branch is structurally organized by the
   coherent (fixed-R_target) side rather than by direct retargeting.

Interpretation
--------------
This does NOT prove that the quartic anomaly correction is literally nothing but
zeta-motion. It proves something narrower and more useful: once the actual
coherent local D/N kernel is inserted, the support lane works by increasing the
baseline M_tr while leaving R_target fixed. So the natural PDE-side continuation
of the branch is a fixed-target/load-compensation problem, i.e. the coherent
side of Step 16. Any direct retargeting would have to come from the mixed/
outgoing microscopic variables, not from the support lane itself.
"""

from __future__ import annotations

import sympy as sp


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


def main() -> None:
    banner("STEP 18 — SUPPORT-COMPENSATION SELECTION LAW")

    zeta, eps = sp.symbols("zeta epsilon", positive=True, real=True)
    M_mix = sp.symbols("M_mix", positive=True, real=True)
    xi, delta, R_tr = sp.symbols("xi delta R_tr", positive=True, real=True)

    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    G_tr = sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R_tr**2) * xi))
    M_crit = sp.simplify(G_tr.subs(xi, 1))

    subbanner("XVIII.1 — Exact support-enhancement factor and tracking load law")

    print("S(zeta;eps) =")
    sp.pprint(S)
    print("G_tr(xi,delta;R_tr) =")
    sp.pprint(G_tr)
    print("M_crit(delta,R_tr) =")
    sp.pprint(M_crit)

    expect_zero(
        "dG_tr/dxi - positive Stage-31 form",
        sp.simplify(
            sp.diff(G_tr, xi)
            - 9 * (2 * R_tr**2 * xi**2 + 9 * delta**2 + 18 * delta * xi + 9 * xi**2)
            / (2 * R_tr**2 * xi + 9 * delta + 9 * xi) ** 2
        ),
    )
    expect_zero(
        "M_crit - G_tr - exact positive remainder",
        sp.simplify(
            M_crit
            - G_tr
            - 9
            * (1 - xi)
            * (2 * R_tr**2 * xi + 9 * delta**2 + 9 * delta * xi + 9 * delta + 9 * xi)
            / ((2 * R_tr**2 + 9 * delta + 9) * (2 * R_tr**2 * xi + 9 * delta + 9 * xi))
        ),
    )

    print("So G_tr is strictly increasing on 0 < xi < 1, and every stable point lies")
    print("strictly below the critical load M_crit.")

    subbanner("XVIII.2 — Exact inverse support law and stability margin")

    S_req = sp.symbols("S_req", positive=True, real=True)
    zeta_req = sp.simplify(sp.solve(sp.Eq(S, S_req), zeta)[0])
    print("zeta_req(S_req) =")
    sp.pprint(zeta_req)

    expect_zero(
        "1/eps - zeta_req - exact positive margin",
        sp.simplify(1 / eps - zeta_req - (1 - eps) / (eps * (1 + eps * (S_req - 2)))),
    )
    expect_zero(
        "dS/dzeta - positive Stage-31 form",
        sp.simplify(sp.diff(S, zeta) - (1 - eps) / (1 - zeta * eps) ** 2),
    )

    print("So every finite required enhancement S_req > 1 sits strictly below the")
    print("blocking pole zeta = 1/eps.")

    subbanner("XVIII.3 — Exact support-compensation theorem")

    M_req = sp.symbols("M_req", positive=True, real=True)
    S_req_load = sp.simplify(M_req / M_mix)
    zeta_req_load = sp.simplify(zeta_req.subs(S_req, S_req_load))
    print("If M_mix < M_req, the required support ratio is")
    print("zeta_req =")
    sp.pprint(zeta_req_load)

    S_crit = sp.simplify(M_crit / M_mix)
    zeta_crit = sp.simplify(zeta_req.subs(S_req, S_crit))
    expect_zero(
        "zeta_crit - zeta_req - exact positive remainder",
        sp.simplify(
            zeta_crit
            - zeta_req_load
            - (S_crit - S_req_load) * (1 - eps)
            / ((1 + eps * (S_crit - 2)) * (1 + eps * (S_req_load - 2)))
        ),
    )

    print("So whenever M_mix < M_req, the unique coherent support ratio zeta_req reaches")
    print("the target before the selected branch softens out, with zeta_req < zeta_crit < 1/eps.")

    subbanner("XVIII.4 — Exact monotone motion deeper into the tracking family")

    dxi_dzeta = sp.simplify(M_mix * sp.diff(S, zeta) / sp.diff(G_tr, xi))
    print("From M_mix S(zeta;eps) = G_tr(xi_phys,delta;R_tr),")
    print("dxi_phys/dzeta =")
    sp.pprint(dxi_dzeta)

    # certify positivity through the already verified factorization.
    dG_pos = 9 * (2 * R_tr**2 * xi**2 + 9 * delta**2 + 18 * delta * xi + 9 * xi**2) / (2 * R_tr**2 * xi + 9 * delta + 9 * xi) ** 2
    expect_zero("dxi/dzeta - M_mix (dS/dzeta)/(dG_tr/dxi)", sp.simplify(dxi_dzeta - M_mix * sp.diff(S, zeta) / dG_pos))

    print("Because dS/dzeta > 0 and dG_tr/dxi > 0 on the stable branch, we have")
    print("  dxi_phys/dzeta > 0.")
    print("So coherent support enhancement always drives the physical branch to larger")
    print("softening depth at fixed demand ratio.")

    subbanner("XVIII.5 — Branch-selection implication for the g-2 quartic layer")

    print("Stage 17 already showed that zeta does not enter R_target at all.")
    print("Stage 18 now shows that zeta only increases M_tr and moves the physical branch")
    print("monotonically deeper into the same tracking family.")
    print()
    print("So the coherent local D/N support lane is structurally a fixed-target /")
    print("load-compensation mechanism. In Step-16 language, that favors the coherent side")
    print("  delta ln R_target = 0")
    print("over the direct retargeting side")
    print("  delta ln R_target = -Lambda_1.")
    print()
    print("This is an inference about the natural PDE-side continuation of the physical")
    print("branch, not yet a proof that the full quartic anomaly correction is exhausted")
    print("by support enhancement alone. Any genuine retargeting would still have to come")
    print("from the mixed/outgoing microscopic variables isolated in Step 17.")


if __name__ == "__main__":
    main()
