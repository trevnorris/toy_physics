#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage265_267_common import banner, subbanner, expect_zero

"""
Stage 266 — exact twin-safe / Family-1 corridor theorem for the first concrete
l=0 <-> l=2 induced geometry-mixing mechanism.

This script takes the Stage-265 one-pole scalar-lane coefficients

    a_2 = mu_mix * r * c_g,
    a_4 = mu_mix * r^2 * c_g * (1-c_g)

and substitutes them into the exact Stage-264 weak-anisotropy corridor tests.

The main exact outputs are:

1. the support-demand drift threshold
       u := r(1-c_g) > 2,
2. the twin-margin shrink threshold
       u > 4,
3. the exact twin-failure threshold
       u >= 4 + 2*sqrt(2),
4. the exact Family-1 margin shrink threshold
       u > 8 c_*,
5. the exact Family-1 failure threshold
       u >= 8 c_* + 4 sqrt(c_*(4 c_* - 1)).

So the first actual geometry-mixing mechanism can threaten the isotropic branch
only when the scalar lane is both sufficiently pole-light (small c_g) and much
faster than the quadrupole carrier (large r).
"""


def main() -> None:
    banner("STAGE 266 — ONE-POLE GEOMETRY-MIXING CORRIDOR THEOREM")

    mu_mix, r, c_g, cstar, y = sp.symbols("mu_mix r c_g c_star y", positive=True, real=True)
    u = sp.symbols("u", real=True)

    a2 = sp.simplify(mu_mix * r * c_g)
    a4 = sp.simplify(mu_mix * r**2 * c_g * (1 - c_g))
    u_def = sp.simplify(r * (1 - c_g))
    alpha = sp.simplify(mu_mix * r * c_g)

    M_twin = sp.expand(1 + (4 * a2 - a4) * y + 2 * a2**2 * y**2)
    M_F1 = sp.expand((4 * cstar - 1) + (8 * cstar * a2 - a4) * y + 4 * cstar * a2**2 * y**2)

    M_twin_u = sp.expand(1 + alpha * (4 - u) * y + 2 * alpha**2 * y**2)
    M_F1_u = sp.expand((4 * cstar - 1) + alpha * (8 * cstar - u) * y + 4 * cstar * alpha**2 * y**2)

    subbanner("I. Exact weak-anisotropy corridor quadratics in u := r(1-c_g)")
    print("u := r(1-c_g) =")
    sp.pprint(u_def)
    print("M_twin(y) =")
    sp.pprint(sp.factor(M_twin_u))
    print("M_F1(y) =")
    sp.pprint(sp.factor(M_F1_u))

    expect_zero(
        "4 a_2 - a_4 - alpha(4-u_def)",
        sp.simplify((4 * a2 - a4) - alpha * (4 - u_def)),
    )
    expect_zero(
        "8 c_* a_2 - a_4 - alpha(8 c_* - u_def)",
        sp.simplify((8 * cstar * a2 - a4) - alpha * (8 * cstar - u_def)),
    )

    subbanner("II. Exact initial-drift thresholds")
    demand_drift = sp.factor(a4 - 2 * a2)
    twin_linear = sp.factor(sp.diff(M_twin, y).subs(y, 0))
    F1_linear = sp.factor(sp.diff(M_F1, y).subs(y, 0))

    print("a_4 - 2 a_2 =")
    sp.pprint(demand_drift)
    print("dM_twin/dy at y=0 =")
    sp.pprint(twin_linear)
    print("dM_F1/dy at y=0 =")
    sp.pprint(F1_linear)

    expect_zero(
        "a_4 - 2 a_2 - alpha(u_def-2)",
        sp.simplify(demand_drift - alpha * (u_def - 2)),
    )
    expect_zero(
        "M_twin'(0) - alpha(4-u_def)",
        sp.simplify(twin_linear - alpha * (4 - u_def)),
    )
    expect_zero(
        "M_F1'(0) - alpha(8 c*-u_def)",
        sp.simplify(F1_linear - alpha * (8 * cstar - u_def)),
    )

    print()
    print("Therefore:")
    print("  - support demand grows initially iff           u > 2,")
    print("  - the twin-safe margin shrinks initially iff   u > 4,")
    print("  - the Family-1 margin shrinks initially iff    u > 8 c_*.")

    subbanner("III. Exact twin-failure condition and first positive root")
    Delta_twin = sp.factor(sp.discriminant(M_twin_u, y))
    u_twin_fail = sp.simplify(4 + 2 * sp.sqrt(2))
    y_twin_crit = sp.simplify((u - 4 - sp.sqrt((u - 4) ** 2 - 8)) / (4 * alpha))

    print("Delta_twin =")
    sp.pprint(Delta_twin)
    print("u_twin^fail =")
    sp.pprint(u_twin_fail)
    print("First positive twin-boundary root y_twin,crit =")
    sp.pprint(y_twin_crit)

    expect_zero(
        "Delta_twin - alpha^2((u-4)^2-8)",
        sp.simplify(Delta_twin - alpha**2 * ((u - 4) ** 2 - 8)),
    )
    expect_zero(
        "M_twin(y_twin,crit)",
        sp.simplify(M_twin_u.subs(y, y_twin_crit)),
    )

    print()
    print("Twin-safety can fail only if u >= 4 + 2*sqrt(2).")
    print("If 0 <= u < 4 + 2*sqrt(2), the one-pole scalar-lane mixing mechanism never drives")
    print("the actual grouped-P2 branch out of the universal twin-safe strip for any y >= 0.")

    subbanner("IV. Exact Family-1 failure condition and first positive root")
    Delta_F1 = sp.factor(sp.discriminant(M_F1_u, y))
    u_F1_fail = sp.simplify(8 * cstar + 4 * sp.sqrt(cstar * (4 * cstar - 1)))
    y_F1_crit = sp.simplify((u - 8 * cstar - sp.sqrt((u - 8 * cstar) ** 2 - 16 * cstar * (4 * cstar - 1))) / (8 * cstar * alpha))

    print("Delta_F1 =")
    sp.pprint(Delta_F1)
    print("u_F1^fail =")
    sp.pprint(u_F1_fail)
    print("First positive Family-1-boundary root y_F1,crit =")
    sp.pprint(y_F1_crit)

    expect_zero(
        "Delta_F1 - alpha^2((u-8 c*)^2 - 16 c*(4 c*-1))",
        sp.simplify(Delta_F1 - alpha**2 * ((u - 8 * cstar) ** 2 - 16 * cstar * (4 * cstar - 1))),
    )
    expect_zero(
        "M_F1(y_F1,crit)",
        sp.simplify(M_F1_u.subs(y, y_F1_crit)),
    )

    print()
    print("The exact Family-1 corridor can fail only if")
    print("      u >= 8 c_* + 4 sqrt[c_*(4 c_* - 1)].")
    print("Below that threshold the first actual l=0<->l=2 mixing mechanism never drives the")
    print("actual grouped-P2 branch out of the exact Family-1 admissible strip for any y >= 0.")

    subbanner("V. Pure-pole scalar lane is automatically safe")
    expect_zero("a_4(c_g=1)", sp.simplify(a4.subs(c_g, 1)))
    expect_zero("u(c_g=1)", sp.simplify(u_def.subs(c_g, 1)))
    print("When c_g = 1 (pure scalar pole, no scalar contact piece):")
    print("  a_4 = 0,  u = 0, so the branch lies far below every drift/failure threshold.")
    print("Therefore pure one-pole scalar-lane mixing can only soften the demand; it can never trigger")
    print("loss of twin-safety or Family-1 admissibility by itself.")

    banner("STAGE 266 LEDGER")
    print("1. The first actual induced l=0<->l=2 geometry-mixing mechanism is controlled by")
    print("      u := r(1-c_g) = (Omega_Q^2/Omega_g^2)(1-c_g).")
    print("2. Exact thresholds are:")
    print("      u > 2                    -> support demand grows initially,")
    print("      u > 4                    -> twin margin shrinks initially,")
    print("      u >= 4 + 2 sqrt(2)       -> actual twin failure can occur,")
    print("      u > 8 c_*                -> Family-1 margin shrinks initially,")
    print("      u >= 8 c_* + 4 sqrt(c_*(4 c_* - 1)) -> actual Family-1 failure can occur.")
    print("3. So the first induced geometry-mixing mechanism threatens the isotropic branch only when")
    print("   the scalar lane is both sufficiently pole-light and sufficiently faster than the quadrupole carrier.")
    print("4. In the pure scalar-pole limit c_g = 1, the branch is automatically safe.")


if __name__ == "__main__":
    main()
