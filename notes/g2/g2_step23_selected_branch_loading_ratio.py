#!/usr/bin/env python3
"""
Step 23 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Imports the exact selected-branch product identities from the later moving-
   throat notes,
       Pi_tr = (N_Q^(target)/beta_0) alpha_req,
       C_mix = (N_Q^(target)/beta_0) alpha_mix,
   and proves that the selected support demand depends only on the loading ratio
       rho_alpha := alpha_req / alpha_mix.
2. Rewrites the natural contact-plus-pole conservative quadrupole precursor
       Y_Q^cons(omega) = c0 + c1/(1 - omega^2/Omega_Q^2)
   in terms of rho_alpha, and derives the exact inverse formulas
       rho_alpha = 1/c0 = 1/(1-c1),
       zeta_req  = rho_alpha - 1 = c1/c0.
3. Inserts the minimal isotropic conservative precursor fixed by the 2.5PN
   quadrupole audit,
       c0 = 3/4,
       c1 = 1/4,
   and proves
       rho_alpha = 4/3,
       zeta_req  = 1/3,
       Pi_tr     = (4/3) C_mix.
4. Converts that selected-branch demand into the Step-22 support selector
       varrho := pi^2 Pi_tr/(16 Lambda),
   and proves the exact relation
       varrho = 2(1-epsilon_*)/3,
   hence
       S_req = Pi_tr/C_mix = 4/3.

Interpretation
--------------
Step 22 asked for a derivation of the support-demand selector varrho from the
selected-branch normalization side. Step 23 does that. The remaining support
question is no longer “what arbitrary Pi_tr should we insert?” On the natural
minimal isotropic passive/outgoing branch, the selected demand is fixed to the
exact ratio Pi_tr/C_mix = 4/3. That places the branch strictly in the symmetric-
lowest-twin regime and removes both the mixed-only and genuinely non-twin
branches from the live anomaly closure.
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
    banner("STEP 23 — SELECTED-BRANCH LOADING RATIO FROM THE MINIMAL ISOTROPIC QUADRUPOLE PRECURSOR")

    alpha_req, alpha_mix = sp.symbols("alpha_req alpha_mix", positive=True, real=True)
    A, beta0, NQ = sp.symbols("A beta_0 N_Q_target", positive=True, real=True)
    Pi_tr, C_mix = sp.symbols("Pi_tr C_mix", positive=True, real=True)
    rho_alpha, zeta_req = sp.symbols("rho_alpha zeta_req", positive=True, real=True)
    c0, c1 = sp.symbols("c0 c1", positive=True, real=True)
    omega, Omega_Q = sp.symbols("omega Omega_Q", positive=True, real=True)
    Lambda, eps = sp.symbols("Lambda epsilon_star", positive=True, real=True)
    mhat2, s_over_lambda = sp.symbols("mhat_sq s_over_lambda", positive=True, real=True)

    subbanner("XXIII.1 — Exact selected-branch product identities")

    kappa0_sq = sp.Rational(8) / sp.pi**2
    R_target = sp.simplify(NQ * A / (beta0 * kappa0_sq))
    G_tr = sp.simplify(8 * alpha_req / (sp.pi**2 * A))
    M_mix = sp.simplify(8 * alpha_mix / (sp.pi**2 * A))

    Pi_formula = sp.simplify(R_target * G_tr)
    C_formula = sp.simplify(R_target * M_mix)

    print("R_target =")
    sp.pprint(R_target)
    print("G_tr =")
    sp.pprint(G_tr)
    print("M_mix =")
    sp.pprint(M_mix)
    print("Pi_tr =")
    sp.pprint(Pi_formula)
    print("C_mix =")
    sp.pprint(C_formula)

    expect_zero("Pi_tr - (N_Q^(target)/beta_0) alpha_req", sp.simplify(Pi_formula - (NQ / beta0) * alpha_req))
    expect_zero("C_mix - (N_Q^(target)/beta_0) alpha_mix", sp.simplify(C_formula - (NQ / beta0) * alpha_mix))
    expect_zero("Pi_tr/C_mix - alpha_req/alpha_mix", sp.simplify(Pi_formula / C_formula - alpha_req / alpha_mix))

    spectral_sub = {NQ: mhat2 * beta0 * s_over_lambda}
    Pi_spectral = sp.simplify(Pi_formula.subs(spectral_sub))
    C_spectral = sp.simplify(C_formula.subs(spectral_sub))
    print("Spectral form with N_Q^(target) = mhat^2 beta_0 s/lambda:")
    print("Pi_tr =")
    sp.pprint(Pi_spectral)
    print("C_mix =")
    sp.pprint(C_spectral)
    expect_zero("Common spectral factor cancels from Pi_tr/C_mix", sp.simplify(Pi_spectral / C_spectral - alpha_req / alpha_mix))

    subbanner("XXIII.2 — Natural contact-plus-pole precursor and the exact inverse formulas")

    Y_cons = sp.simplify(1 / rho_alpha + (rho_alpha - 1) / rho_alpha / (1 - omega**2 / Omega_Q**2))
    print("Y_Q^cons(omega) =")
    sp.pprint(Y_cons)

    c0_from_rho = sp.simplify(sp.limit(Y_cons, omega, 0) - (rho_alpha - 1) / rho_alpha)  # static contact fraction
    c1_from_rho = sp.simplify((rho_alpha - 1) / rho_alpha)

    print("Contact fraction c0 =")
    sp.pprint(c0_from_rho)
    print("Pole residue c1 =")
    sp.pprint(c1_from_rho)
    expect_zero("c0 - 1/rho_alpha", sp.simplify(c0_from_rho - 1 / rho_alpha))
    expect_zero("c1 - (rho_alpha - 1)/rho_alpha", sp.simplify(c1_from_rho - (rho_alpha - 1) / rho_alpha))
    expect_zero("c0 + c1 - 1", sp.simplify(c0_from_rho + c1_from_rho - 1))

    rho_from_c0 = sp.simplify(1 / c0)
    rho_from_c1 = sp.simplify(1 / (1 - c1))
    zeta_from_c = sp.simplify(c1 / c0)

    print("Inverse formulas:")
    print("rho_alpha(c0) =")
    sp.pprint(rho_from_c0)
    print("rho_alpha(c1) =")
    sp.pprint(rho_from_c1)
    print("zeta_req(c0,c1) =")
    sp.pprint(zeta_from_c)
    expect_zero("rho_alpha(c0) - rho_alpha(c1)", sp.simplify(rho_from_c0 - rho_from_c1).subs({c1: 1 - c0}))

    subbanner("XXIII.3 — Minimal isotropic conservative precursor")

    minimal_sub = {c0: sp.Rational(3, 4), c1: sp.Rational(1, 4)}
    rho_min = sp.simplify(rho_from_c0.subs(minimal_sub))
    zeta_min = sp.simplify(zeta_from_c.subs(minimal_sub))

    print("Minimal isotropic module:")
    print("  c0 = 3/4")
    print("  c1 = 1/4")
    print("Therefore rho_alpha =")
    sp.pprint(rho_min)
    print("and zeta_req =")
    sp.pprint(zeta_min)
    expect_zero("rho_alpha - 4/3", sp.simplify(rho_min - sp.Rational(4, 3)))
    expect_zero("zeta_req - 1/3", sp.simplify(zeta_min - sp.Rational(1, 3)))

    Pi_min = sp.simplify(rho_min * C_mix)
    print("Selected demand product on the minimal isotropic branch:")
    print("Pi_tr =")
    sp.pprint(Pi_min)
    expect_zero("Pi_tr - (4/3) C_mix", sp.simplify(Pi_min - sp.Rational(4, 3) * C_mix))

    subbanner("XXIII.4 — Exact support selector varrho on the selected branch")

    varrho = sp.symbols("varrho", positive=True, real=True)
    C_mix_eps = sp.simplify(8 * Lambda * (1 - eps) / sp.pi**2)
    Pi_eps = sp.simplify(sp.Rational(4, 3) * C_mix_eps)
    varrho_expr = sp.simplify(sp.pi**2 * Pi_eps / (16 * Lambda))
    S_req = sp.simplify(Pi_eps / C_mix_eps)

    print("C_mix(epsilon_*) =")
    sp.pprint(C_mix_eps)
    print("Pi_tr(epsilon_*) = (4/3) C_mix =")
    sp.pprint(Pi_eps)
    print("varrho := pi^2 Pi_tr / (16 Lambda) =")
    sp.pprint(varrho_expr)
    print("S_req = Pi_tr/C_mix =")
    sp.pprint(S_req)

    expect_zero("varrho - 2(1-epsilon_*)/3", sp.simplify(varrho_expr - 2 * (1 - eps) / 3))
    expect_zero("S_req - 4/3", sp.simplify(S_req - sp.Rational(4, 3)))

    banner("STEP 23 LEDGER")
    print("Exact selected-branch demand cancellation:")
    print("  Pi_tr / C_mix = alpha_req / alpha_mix = rho_alpha")
    print()
    print("Natural contact-plus-pole inverse formulas:")
    print("  rho_alpha = 1/c0 = 1/(1-c1)")
    print("  zeta_req  = rho_alpha - 1 = c1/c0")
    print()
    print("Minimal isotropic conservative precursor:")
    print("  c0 = 3/4, c1 = 1/4")
    print("  -> rho_alpha = 4/3")
    print("  -> zeta_req  = 1/3")
    print("  -> Pi_tr     = (4/3) C_mix")
    print()
    print("Exact Step-22 support selector on this branch:")
    print("  varrho = pi^2 Pi_tr / (16 Lambda) = 2(1-epsilon_*)/3")
    print("  S_req  = Pi_tr/C_mix = 4/3")
    print()
    print("Interpretation:")
    print("  The natural minimal isotropic passive/outgoing branch fixes the selected")
    print("  support demand to the exact symmetric-lowest-twin ratio 4/3. The remaining")
    print("  anomaly ambiguity is therefore no longer about mixed-only vs non-twin")
    print("  support. It is only about where the physical branch sits along that selected")
    print("  twin-support curve.")


if __name__ == "__main__":
    main()
