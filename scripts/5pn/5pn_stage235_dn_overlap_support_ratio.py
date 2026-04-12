#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 235 — explicit D/N overlap extraction of the coherent support ratio.
"""


def main() -> None:
    banner("STAGE 235 — D/N OVERLAP EXTRACTION OF ZETA")

    n = sp.symbols("n", integer=True, nonnegative=True)
    s, L = sp.symbols("s L", positive=True, real=True)
    K_W_eff = sp.symbols("K_W_eff", positive=True, real=True)
    K_phi_n_eff = sp.symbols("K_phi_n_eff", positive=True, real=True)
    x = sp.symbols("x", nonnegative=True, real=True)

    chi_n = sp.sqrt(2 / L) * sp.sin((n + sp.Rational(1, 2)) * sp.pi * s / L)
    I_n = sp.simplify(sp.integrate(chi_n, (s, 0, L)))
    I_0 = sp.simplify(I_n.subs(n, 0))
    I_ratio = sp.simplify(I_n / I_0)

    zeta_phys = sp.simplify((K_W_eff / K_phi_n_eff) * I_ratio**2)
    zeta_twin = sp.simplify(zeta_phys.subs(K_phi_n_eff, K_W_eff * (1 + x * n * (n + 1))))

    subbanner("I. Exact D/N overlap law")
    print("chi_n(s) =")
    sp.pprint(chi_n)
    print()
    print("I_n =")
    sp.pprint(I_n)
    print()
    print("I_0 =")
    sp.pprint(I_0)
    print()
    print("I_n/I_0 =")
    sp.pprint(I_ratio)
    expect_zero("exact overlap ratio law", sp.simplify(I_ratio - 1 / (2 * n + 1)))

    subbanner("II. Exact physical coherent support ratio")
    print("zeta_n^(phys) =")
    sp.pprint(zeta_phys)
    expect_zero(
        "physical zeta law",
        sp.simplify(zeta_phys - (K_W_eff / K_phi_n_eff) / (2 * n + 1) ** 2),
    )

    subbanner("III. Same-operator twin family")
    print("zeta_n^(twin) =")
    sp.pprint(zeta_twin)
    expect_zero(
        "twin zeta law",
        sp.simplify(zeta_twin - 1 / ((2 * n + 1) ** 2 * (1 + x * n * (n + 1)))),
    )
    print()
    print("zeta_0^(twin) =")
    sp.pprint(sp.simplify(zeta_twin.subs(n, 0)))
    expect_zero("lowest twin value", sp.simplify(zeta_twin.subs(n, 0) - 1))

    banner("STAGE 235 LEDGER")
    print("1. On the exact finite-throat D/N family with a uniform local source, the overlap law is")
    print("      I_n / I_0 = 1 / (2n+1).")
    print("2. Therefore the physical coherent support ratio is")
    print("      zeta_n^(phys) = (K_W^(eff) / K_{phi,n}^(eff)) / (2n+1)^2.")
    print("3. On the same-operator twin family this collapses to")
    print("      zeta_n^(twin) = 1 / [ (2n+1)^2 (1 + x n(n+1)) ].")
    print("4. The lowest symmetric twin lane has the exact universal value")
    print("      zeta_0^(twin) = 1,")
    print("   so it is the strongest coherent support branch in this explicit D/N tower.")


if __name__ == "__main__":
    main()
