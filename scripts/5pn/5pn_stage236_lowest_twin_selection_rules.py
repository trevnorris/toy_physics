#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 236 — exact lowest-twin sufficiency and higher-harmonic selection rules.
"""


def main() -> None:
    banner("STAGE 236 — LOWEST TWIN SELECTION RULES")

    n = sp.symbols("n", integer=True, positive=True)
    eps = sp.symbols("eps", positive=True, real=True)
    zeta_req = sp.symbols("zeta_req", positive=True, real=True)
    S_req = sp.symbols("S_req", positive=True, real=True)
    x = sp.symbols("x", nonnegative=True, real=True)
    zeta = sp.symbols("zeta", positive=True, real=True)

    S = lambda zz: sp.simplify(1 + zz * (1 - eps) / (1 - zz * eps))
    zeta_req_from_S = sp.simplify(sp.solve(sp.Eq(S_req, S(zeta)), zeta)[0])

    zeta_0_twin = sp.Integer(1)
    S_0 = sp.simplify(S(zeta_0_twin))

    zeta_n_twin = sp.simplify(1 / ((2 * n + 1) ** 2 * (1 + x * n * (n + 1))))
    zeta_n_max = sp.simplify(zeta_n_twin.subs(x, 0))
    x_max = sp.simplify(sp.solve(sp.Eq(zeta_n_twin, zeta_req), x)[0])
    S_n_max = sp.simplify(S(zeta_n_max))

    subbanner("I. Exact lowest-twin enhancement")
    print("S_0 =")
    sp.pprint(S_0)
    expect_zero("lowest twin enhancement", sp.simplify(S_0 - 2))

    subbanner("II. Exact lowest-twin sufficiency criterion")
    print("zeta_req(S_req;eps) =")
    sp.pprint(zeta_req_from_S)
    print()
    print("Lowest twin succeeds iff zeta_req <= 1, equivalently S_req <= 2.")
    expect_zero("S(zeta=1)-2", sp.simplify(S(sp.Integer(1)) - 2))
    expect_zero("zeta_req(S_req=2)-1", sp.simplify(zeta_req_from_S.subs(S_req, 2) - 1))

    subbanner("III. Higher-harmonic immediate impossibility bound")
    print("zeta_n^(max) =")
    sp.pprint(zeta_n_max)
    print()
    print("So harmonic n is ruled out immediately if")
    print("  zeta_req > 1/(2n+1)^2.")

    subbanner("IV. Exact twin softness threshold")
    print("x_max(n;zeta_req) =")
    sp.pprint(x_max)
    expect_zero(
        "x_max identity",
        sp.simplify(x_max - (1 / ((2 * n + 1) ** 2 * zeta_req) - 1) / (n * (n + 1))),
    )

    subbanner("V. Exact higher-harmonic enhancement ceiling")
    print("S_n^(max) =")
    sp.pprint(S_n_max)
    expect_zero(
        "higher-harmonic enhancement ceiling",
        sp.simplify(S_n_max - (1 + (1 - eps) / ((2 * n + 1) ** 2 - eps))),
    )

    banner("STAGE 236 LEDGER")
    print("1. The symmetric lowest twin lane has zeta_0^(twin) = 1, so its enhancement is")
    print("      S_0 = 2")
    print("   exactly, independent of eps.")
    print("2. Therefore the lowest twin lane succeeds exactly when")
    print("      zeta_req <= 1,   equivalently   S_req <= 2.")
    print("3. Every higher D/N harmonic obeys the immediate impossibility bound")
    print("      zeta_req > 1/(2n+1)^2  =>  harmonic n impossible.")
    print("4. When that bound is not violated, the exact twin softness window is")
    print("      x <= x_max(n;zeta_req) = [1/((2n+1)^2 zeta_req) - 1]/[n(n+1)].")
    print("5. The corresponding enhancement ceiling is")
    print("      S_n^(max) = 1 + (1-eps)/((2n+1)^2 - eps),")
    print("   which is only modest for n>=1. So the explicit D/N support tower is strongly biased")
    print("   toward the lowest symmetric twin lane.")


if __name__ == "__main__":
    main()
