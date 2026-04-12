#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage329_332_common import *

"""
Stage 331 — actual coherent local D/N support-threshold comparison.

What this stage does
--------------------
1. Evaluates the exact Stage-31 support threshold on the actual coherent local D/N branch.
2. Eliminates zeta_req in favor of the tracking-branch product Pi_tr and the mixed-only
   product scale C_mix = 8 Lambda (1-epsilon)/pi^2.
3. Proves the exact three-regime split:
      - mixed-only enough,
      - symmetric lowest twin enough,
      - non-twin asymmetry required.
4. Records the exact higher-harmonic impossibility and softness bounds.
"""

if __name__ == "__main__":
    banner("STAGE 331 — ACTUAL COHERENT LOCAL D/N SUPPORT THRESHOLD")

    chi0, deltaU, ZW, epsW, epsEta, Lambda = sp.symbols(
        "chi0 deltaU Z_W epsilon_W epsilon_eta Lambda", positive=True, real=True
    )
    xi, delta = sp.symbols("xi delta", positive=True, real=True)
    R = sp.symbols("R", positive=True, real=True)
    n = sp.symbols("n", integer=True, nonnegative=True)
    x = sp.symbols("x", positive=True, real=True)

    branch = coherent_tracking_state(chi0, deltaU, ZW, epsW, epsEta, Lambda, sp.Symbol("zeta", nonnegative=True, real=True))
    Pi = Pi_tr(xi, delta, R)
    Cmix = sp.simplify(branch["Ctr"])
    Sreq = S_req_from_branch(Pi, Cmix)
    zetareq = zeta_req_from_product(Pi, Cmix, branch["eps"])
    zeta0 = zeta_twin_same_operator(sp.Integer(0), x)
    S0 = sp.simplify(S_of_zeta(zeta0, branch["eps"]))

    subbanner("I. Exact tracking-branch support threshold")
    print("Pi_tr(xi,delta;R) =")
    sp.pprint(Pi)
    print("C_mix = 8 Lambda (1-epsilon)/pi^2 =")
    sp.pprint(Cmix)
    print("S_req = Pi_tr / C_mix =")
    sp.pprint(Sreq)
    print("zeta_req =")
    sp.pprint(zetareq)

    subbanner("II. Exact lowest symmetric twin criterion")
    print("zeta_0^(twin) =")
    sp.pprint(zeta0)
    print("S_0 =")
    sp.pprint(S0)
    expect_zero("zeta_0^(twin) - 1", zeta0 - 1)
    expect_zero("S_0 - 2", S0 - 2)
    twin_boundary = lowest_twin_condition(Pi, Cmix)
    mixed_boundary = mixed_only_condition(Pi, Cmix)
    print("mixed-only boundary Pi_tr - C_mix =")
    sp.pprint(mixed_boundary)
    print("lowest-twin boundary Pi_tr - 2 C_mix =")
    sp.pprint(twin_boundary)

    subbanner("III. Exact regime split")
    print("Mixed-only enough iff")
    print("  Pi_tr <= C_mix.")
    print("Symmetric lowest twin enough iff")
    print("  C_mix < Pi_tr <= 2 C_mix.")
    print("Non-twin asymmetry required iff")
    print("  Pi_tr > 2 C_mix.")
    zeta_req_mixed_only = sp.simplify(zetareq.subs({Pi: Cmix}))
    zeta_req_lowest_twin = sp.simplify(zetareq.subs({Pi: 2 * Cmix}))
    print("zeta_req at Pi_tr = C_mix ->")
    sp.pprint(zeta_req_mixed_only)
    print("zeta_req at Pi_tr = 2 C_mix ->")
    sp.pprint(zeta_req_lowest_twin)
    expect_zero("zeta_req(Pi=C_mix)", zeta_req_mixed_only)
    expect_zero("zeta_req(Pi=2C_mix) - 1", zeta_req_lowest_twin - 1)

    subbanner("IV. Exact higher-harmonic bounds")
    zeta_n = zeta_twin_same_operator(n, x)
    print("zeta_n^(twin) =")
    sp.pprint(zeta_n)
    print("Necessary impossibility bound for nth twin harmonic:")
    print("  zeta_req <= 1/(2n+1)^2")
    print("Exact softness threshold when not already impossible:")
    xmax = x_max_higher_harmonic(n, zetareq)
    sp.pprint(xmax)
    print("Enhancement ceiling S_n^(max) =")
    sp.pprint(S_n_max(n, branch["eps"]))

    subbanner("V. Interpretation")
    print("On the actual coherent local D/N branch, the support theorem is now completely")
    print("reduced to the tracking-branch product Pi_tr and the mixed-only product scale C_mix.")
    print()
    print("The first concrete physical support lane — the symmetric lowest twin branch — is")
    print("exactly the universal doubling branch.  So the coherent support problem on the")
    print("actual branch is no longer vague: either mixed loading already suffices, or the")
    print("lowest twin suffices, or genuinely non-twin asymmetry is required.")
