#!/usr/bin/env python3
"""
Step 22 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Imports the exact coherent support-compensation structure from the later
   moving-throat notes:
       S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps),
       C_mix       = 8 Lambda (1-eps)/pi^2,
       S_req       = Pi_tr / C_mix,
   together with the regime split
       Pi_tr <= C_mix          : mixed-only enough,
       C_mix < Pi_tr <= 2 C_mix: symmetric lowest twin enough,
       Pi_tr > 2 C_mix         : non-twin asymmetry required.
2. Repackages the support demand into the single dimensionless selector
       varrho := pi^2 Pi_tr / (16 Lambda),
   and proves the exact epsilon_* windows for the three support regimes.
3. Converts those epsilon_* windows into exact sigma windows, where
       sigma = 2 epsilon_*/(1-epsilon_*)
   is the coefficient that actually appears in the coherent quartic microledger.
4. Intersects the support-selector windows with the Step-21 primitive ranking
   crossover q_W = q_Lambda at
       epsilon_* = 1/(2+beta^2),
   and derives the exact demand thresholds telling us when support-compatible
   branches can or cannot enter the strong-blocking q_Lambda > q_W phase.

Interpretation
--------------
Step 21 left one live ranking ambiguity: whether the quartic layer is carried
more strongly by q_W or q_Lambda. Step 22 shows that the moving-throat support
selection side already constrains that strongly. For moderate or large normalized
support demand varrho, every mixed-only or symmetric-twin-compatible branch stays
on the q_W-dominant side. The q_Lambda > q_W regime is then forced either into a
low-demand corner or into a genuinely non-twin support branch.
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
    banner("STEP 22 — SUPPORT-DEMAND SELECTOR AND THE EPSILON_* WINDOW")

    eps = sp.symbols("epsilon_star", positive=True, real=True)
    beta = sp.symbols("beta", positive=True, real=True)
    varrho = sp.symbols("varrho", positive=True, real=True)
    Lambda, Pi_tr = sp.symbols("Lambda Pi_tr", positive=True, real=True)

    sigma = sp.simplify(2 * eps / (1 - eps))
    eps_cross = sp.simplify(1 / (2 + beta**2))
    sigma_cross = sp.simplify(2 * eps_cross / (1 - eps_cross))

    subbanner("XXII.1 — Exact support selector from the moving-throat branch product")

    C_mix = sp.simplify(8 * Lambda * (1 - eps) / sp.pi**2)
    varrho_def = sp.simplify(sp.pi**2 * Pi_tr / (16 * Lambda))
    Pi_from_varrho = sp.simplify(16 * Lambda * varrho / sp.pi**2)
    S_req = sp.simplify(Pi_from_varrho / C_mix)

    print("C_mix =")
    sp.pprint(C_mix)
    print("varrho := pi^2 Pi_tr / (16 Lambda)")
    print("Pi_tr(varrho) =")
    sp.pprint(Pi_from_varrho)
    print("S_req(varrho, epsilon_*) = Pi_tr/C_mix =")
    sp.pprint(S_req)
    expect_zero("S_req - 2 varrho/(1-epsilon_*)", sp.simplify(S_req - 2 * varrho / (1 - eps)))

    subbanner("XXII.2 — Exact epsilon_* windows for the three support regimes")

    eps_mix_max = sp.simplify(1 - 2 * varrho)
    eps_twin_max = sp.simplify(1 - varrho)

    print("Mixed-only enough  (Pi_tr <= C_mix)  <=>  epsilon_* <=")
    sp.pprint(eps_mix_max)
    print("Symmetric lowest twin enough  (Pi_tr <= 2 C_mix)  <=>  epsilon_* <=")
    sp.pprint(eps_twin_max)
    print("Therefore the exact regime split is:")
    print("  mixed-only enough:         0 < epsilon_* <= 1 - 2 varrho")
    print("  symmetric lowest twin:     1 - 2 varrho < epsilon_* <= 1 - varrho")
    print("  non-twin asymmetry needed: epsilon_* > 1 - varrho")

    expect_zero("Pi_tr - C_mix at epsilon_* = 1 - 2 varrho", sp.simplify(Pi_from_varrho - C_mix.subs(eps, eps_mix_max)))
    expect_zero("Pi_tr - 2 C_mix at epsilon_* = 1 - varrho", sp.simplify(Pi_from_varrho - 2 * C_mix.subs(eps, eps_twin_max)))

    print("Existence conditions:")
    print("  mixed-only window is nonempty iff varrho < 1/2")
    print("  twin-sufficient window is nonempty iff varrho < 1")
    print("  if varrho >= 1, the branch is automatically non-twin-required for every epsilon_* in (0,1)")

    subbanner("XXII.3 — Exact sigma windows")

    sigma_mix_max = sp.simplify(2 * eps_mix_max / (1 - eps_mix_max))
    sigma_twin_max = sp.simplify(2 * eps_twin_max / (1 - eps_twin_max))

    print("sigma = 2 epsilon_*/(1-epsilon_*)")
    print("Mixed-only upper sigma ceiling:")
    sp.pprint(sigma_mix_max)
    print("Twin-sufficient upper sigma ceiling:")
    sp.pprint(sigma_twin_max)
    expect_zero("sigma_mix_max - (1/varrho - 2)", sp.simplify(sigma_mix_max - (1 / varrho - 2)))
    expect_zero("sigma_twin_max - (2/varrho - 2)", sp.simplify(sigma_twin_max - (2 / varrho - 2)))

    print("So the exact sigma windows are:")
    print("  mixed-only enough:         0 < sigma <= 1/varrho - 2")
    print("  symmetric lowest twin:     1/varrho - 2 < sigma <= 2/varrho - 2")
    print("  non-twin asymmetry needed: sigma > 2/varrho - 2")

    subbanner("XXII.4 — Intersecting the support selector with the Step-21 q_W/q_Lambda crossover")

    print("Step-21 crossover:")
    print("  q_W = q_Lambda at epsilon_cross =")
    sp.pprint(eps_cross)
    print("Equivalent sigma crossover sigma_cross =")
    sp.pprint(sigma_cross)
    expect_zero("sigma_cross - 2/(1+beta^2)", sp.simplify(sigma_cross - 2 / (1 + beta**2)))

    varrho_mix_cross = sp.simplify((1 - eps_cross) / 2)
    varrho_twin_cross = sp.simplify(1 - eps_cross)

    print("Demand thresholds controlling whether support-compatible branches can enter q_Lambda > q_W:")
    print("  varrho_mix^x  = (1 - epsilon_cross)/2 =")
    sp.pprint(varrho_mix_cross)
    print("  varrho_twin^x = 1 - epsilon_cross =")
    sp.pprint(varrho_twin_cross)
    expect_zero("varrho_mix^x - (1+beta^2)/(2(2+beta^2))", sp.simplify(varrho_mix_cross - (1 + beta**2) / (2 * (2 + beta**2))))
    expect_zero("varrho_twin^x - (1+beta^2)/(2+beta^2)", sp.simplify(varrho_twin_cross - (1 + beta**2) / (2 + beta**2)))

    print("Consequences:")
    print("  - Mixed-only branch can reach q_Lambda > q_W only if varrho < varrho_mix^x.")
    print("  - Symmetric-twin-compatible branch can reach q_Lambda > q_W only if varrho < varrho_twin^x.")
    print("  - If varrho >= varrho_twin^x, then every mixed-only or twin-sufficient branch point stays on q_W > q_Lambda.")

    subbanner("XXII.5 — Exact phase classification in the (varrho, beta) plane")

    print("Phase A: varrho >= varrho_twin^x")
    print("  -> even the upper twin ceiling epsilon_* = 1 - varrho lies below the Step-21 crossover.")
    print("  -> all mixed-only and twin-sufficient branches are forced into q_W > q_Lambda.")
    print()
    print("Phase B: varrho_mix^x <= varrho < varrho_twin^x")
    print("  -> mixed-only branch is still forced into q_W > q_Lambda,")
    print("     but the symmetric-twin window straddles the crossover.")
    print("  -> q_Lambda > q_W is then possible only with support help and only in the upper part of the twin window, or beyond it on a non-twin branch.")
    print()
    print("Phase C: varrho < varrho_mix^x")
    print("  -> even the mixed-only branch can reach the strong-blocking q_Lambda > q_W phase.")

    beta_max = sp.Rational(2, 11)
    varrho_mix_min = sp.simplify(sp.limit(varrho_mix_cross, beta, 0, dir="+"))
    varrho_mix_max = sp.simplify(varrho_mix_cross.subs(beta, beta_max))
    varrho_twin_min = sp.simplify(sp.limit(varrho_twin_cross, beta, 0, dir="+"))
    varrho_twin_max = sp.simplify(varrho_twin_cross.subs(beta, beta_max))

    print("Using the constructive coherent bound 0 < beta < 2/11, these thresholds live in the narrow windows:")
    print("  varrho_mix^x  in (")
    sp.pprint(varrho_mix_min)
    print(",")
    sp.pprint(varrho_mix_max)
    print(")")
    print("  varrho_twin^x in (")
    sp.pprint(varrho_twin_min)
    print(",")
    sp.pprint(varrho_twin_max)
    print(")")

    subbanner("XXII.6 — Representative samples")

    samples = [
        (sp.Rational(1, 20), sp.Rational(1, 10)),
        (sp.Rational(1, 20), sp.Rational(3, 10)),
        (sp.Rational(1, 6), sp.Rational(3, 10)),
        (sp.Rational(1, 6), sp.Rational(3, 5)),
    ]
    for beta_val, varrho_val in samples:
        eps_cross_val = sp.N(eps_cross.subs(beta, beta_val))
        eps_mix_val = sp.N(eps_mix_max.subs(varrho, varrho_val))
        eps_twin_val = sp.N(eps_twin_max.subs(varrho, varrho_val))
        print(f"beta = {beta_val}, varrho = {varrho_val}:")
        print(f"  epsilon_cross = {eps_cross_val}")
        print(f"  mixed-only ceiling = {eps_mix_val}")
        print(f"  twin ceiling       = {eps_twin_val}")
        if eps_twin_val <= eps_cross_val:
            verdict = "all mixed/twin-compatible branches force q_W > q_Lambda"
        elif eps_mix_val <= eps_cross_val:
            verdict = "mixed-only forces q_W > q_Lambda, but the twin window can cross into q_Lambda > q_W"
        else:
            verdict = "even the mixed-only window can reach q_Lambda > q_W"
        print(f"  verdict: {verdict}\n")

    banner("STEP 22 LEDGER")
    print("Exact support selector:")
    print("  varrho := pi^2 Pi_tr / (16 Lambda)")
    print("  mixed-only enough:         0 < epsilon_* <= 1 - 2 varrho")
    print("  symmetric lowest twin:     1 - 2 varrho < epsilon_* <= 1 - varrho")
    print("  non-twin asymmetry needed: epsilon_* > 1 - varrho")
    print()
    print("Exact sigma windows:")
    print("  mixed-only enough:         0 < sigma <= 1/varrho - 2")
    print("  symmetric lowest twin:     1/varrho - 2 < sigma <= 2/varrho - 2")
    print("  non-twin asymmetry needed: sigma > 2/varrho - 2")
    print()
    print("Step-21 crossover intersection:")
    print("  q_W = q_Lambda at epsilon_* = 1/(2+beta^2)")
    print("  mixed-only can reach q_Lambda > q_W only if varrho < (1+beta^2)/(2(2+beta^2))")
    print("  twin-compatible branch can reach q_Lambda > q_W only if varrho < (1+beta^2)/(2+beta^2)")
    print()
    print("Interpretation:")
    print("  support selection already constrains the last live quartic dominance ambiguity.")
    print("  For moderate or large normalized support demand varrho, every mixed-only or")
    print("  symmetric-twin-compatible branch remains on the q_W-dominant side, and")
    print("  q_Lambda > q_W becomes a low-demand or genuinely non-twin effect.")


if __name__ == "__main__":
    main()
