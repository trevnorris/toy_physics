#!/usr/bin/env python3
"""
Step 25 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the conservative 3PN split suggested by the moving-throat notes:
   one isotropic grouped-P2 pole plus a static geometry completion.
2. Expands the resulting conservative quadrupole module through O(omega^4) and
   imposes the already-frozen minimal isotropic quadrupole branch identity
       K0*K4 = 4*K2^2.
3. Proves that this forces
       K_geom = 3 K_pole,
   so the normalized conservative module is exactly
       3/4 + (1/4)/(1 - omega^2/Omega_Q^2).
4. Uses the exact support/source contact-plus-pole map to convert that forced
   3/4 + 1/4 split into the selected loading verdict
       rho_alpha = 4/3,
       zeta_req  = 1/3,
       Pi_tr/C_mix = 4/3.

Interpretation
--------------
This is the first step where the old support-side phase diagram stops being the
right organizing object. Once the higher conservative payload is read in the
actual grouped-P2 + static-geometry language, the minimal isotropic conservative
module is not guessed; it is forced. The support/source ratio 4/3 then follows
as a corollary rather than as a separate free datum.
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
    banner("STEP 25 — GROUPED-P2 + STATIC-GEOMETRY SPLIT FORCES THE 3/4 + 1/4 MODULE")

    omega, Omega_Q = sp.symbols("omega Omega_Q", positive=True, real=True)
    K_geom, K_pole = sp.symbols("K_geom K_pole", positive=True, real=True)
    K0, K2, K4 = sp.symbols("K0 K2 K4", real=True)

    subbanner("XXV.1 — Minimal isotropic grouped-P2 + geometry realization")

    K_cons = sp.simplify(K_geom + K_pole / (1 - omega**2 / Omega_Q**2))
    K_series = sp.expand(sp.series(K_cons, omega, 0, 6).removeO())
    print("K_Q^cons(omega) =")
    sp.pprint(K_cons)
    print("Series through O(omega^4) =")
    sp.pprint(K_series)

    K0_expr = sp.simplify(K_series.coeff(omega, 0))
    K2_expr = sp.simplify(K_series.coeff(omega, 2))
    K4_expr = sp.simplify(K_series.coeff(omega, 4))

    print("K0 =")
    sp.pprint(K0_expr)
    print("K2 =")
    sp.pprint(K2_expr)
    print("K4 =")
    sp.pprint(K4_expr)

    expect_zero("K0 - (K_geom + K_pole)", sp.simplify(K0_expr - (K_geom + K_pole)))
    expect_zero("K2 - K_pole/Omega_Q^2", sp.simplify(K2_expr - K_pole / Omega_Q**2))
    expect_zero("K4 - K_pole/Omega_Q^4", sp.simplify(K4_expr - K_pole / Omega_Q**4))

    subbanner("XXV.2 — Minimal isotropic branch identity forces the split")

    branch_identity = sp.simplify(K0_expr * K4_expr - 4 * K2_expr**2)
    print("K0*K4 - 4*K2^2 =")
    sp.pprint(branch_identity)

    Kgeom_sol = sp.solve(sp.Eq(branch_identity, 0), K_geom)[0]
    print("Solved static geometry completion:")
    print("K_geom =")
    sp.pprint(Kgeom_sol)
    expect_zero("K_geom - 3*K_pole", sp.simplify(Kgeom_sol - 3 * K_pole))

    K0_sol = sp.simplify(K0_expr.subs(K_geom, Kgeom_sol))
    print("Then K0 =")
    sp.pprint(K0_sol)
    expect_zero("K0 - 4*K_pole", sp.simplify(K0_sol - 4 * K_pole))

    Yhat = sp.simplify((K_cons / K0_sol).subs(K_geom, Kgeom_sol))
    print("Normalized conservative module =")
    sp.pprint(Yhat)
    target_Yhat = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / Omega_Q**2))
    expect_zero("Yhat - (3/4 + (1/4)/(1 - omega^2/Omega_Q^2))", sp.simplify(Yhat - target_Yhat))

    subbanner("XXV.3 — Exact contact/pole map to the selected loading ratio")

    alpha_req, alpha_mix = sp.symbols("alpha_req alpha_mix", positive=True, real=True)
    rho_alpha, zeta_req = sp.symbols("rho_alpha zeta_req", positive=True, real=True)
    c0, c1 = sp.symbols("c0 c1", positive=True, real=True)
    Pi_tr, C_mix = sp.symbols("Pi_tr C_mix", positive=True, real=True)

    c0_expr = sp.simplify(alpha_mix / alpha_req)
    c1_expr = sp.simplify((alpha_req - alpha_mix) / alpha_req)
    rho_expr = sp.simplify(alpha_req / alpha_mix)
    zeta_expr = sp.simplify((alpha_req - alpha_mix) / alpha_mix)

    print("Contact fraction c0 =")
    sp.pprint(c0_expr)
    print("Pole fraction c1 =")
    sp.pprint(c1_expr)
    print("Loading ratio rho_alpha =")
    sp.pprint(rho_expr)
    print("Support ratio zeta_req =")
    sp.pprint(zeta_expr)

    expect_zero("c0 + c1 - 1", sp.simplify(c0_expr + c1_expr - 1))
    expect_zero("rho_alpha - 1/c0", sp.simplify(rho_expr - 1 / c0_expr))
    expect_zero("zeta_req - c1/c0", sp.simplify(zeta_expr - c1_expr / c0_expr))

    rho_min = sp.simplify(1 / sp.Rational(3, 4))
    zeta_min = sp.simplify(sp.Rational(1, 4) / sp.Rational(3, 4))
    print("Forced values from c0 = 3/4, c1 = 1/4:")
    print("rho_alpha =")
    sp.pprint(rho_min)
    print("zeta_req =")
    sp.pprint(zeta_min)
    expect_zero("rho_alpha - 4/3", sp.simplify(rho_min - sp.Rational(4, 3)))
    expect_zero("zeta_req - 1/3", sp.simplify(zeta_min - sp.Rational(1, 3)))

    Pi_ratio = sp.simplify(rho_min * C_mix)
    print("Therefore Pi_tr =")
    sp.pprint(Pi_ratio)
    expect_zero("Pi_tr/C_mix - 4/3", sp.simplify(Pi_ratio / C_mix - sp.Rational(4, 3)))

    subbanner("XXV.4 — What this means for the selected branch")

    print("If the actual conservative higher-order branch is")
    print("  one isotropic grouped-P2 pole + static geometry completion,")
    print("then the contact fraction is forced to 3/4 and the pole fraction to 1/4.")
    print("So the support/source loading ratio is not free anymore:")
    print("  rho_alpha = 4/3,  zeta_req = 1/3,  Pi_tr/C_mix = 4/3.")
    print("The live question therefore shifts away from support-side phase selection")
    print("and toward whether the actual moving-throat branch really realizes this")
    print("minimal grouped-P2 + static-geometry conservative module.")

    banner("STEP 25 LEDGER")
    print("Minimal isotropic conservative module:")
    print("  K_Q^cons(omega) = K_geom + K_pole/(1 - omega^2/Omega_Q^2)")
    print("  K0 = K_geom + K_pole")
    print("  K2 = K_pole/Omega_Q^2")
    print("  K4 = K_pole/Omega_Q^4")
    print()
    print("Minimal isotropic branch identity:")
    print("  K0*K4 = 4*K2^2")
    print("  -> K_geom = 3*K_pole")
    print("  -> K_pole = K0/4, K_geom = 3*K0/4")
    print()
    print("Forced normalized conservative response:")
    print("  Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)")
    print()
    print("Exact support/source corollary:")
    print("  rho_alpha = 4/3")
    print("  zeta_req  = 1/3")
    print("  Pi_tr/C_mix = 4/3")
    print()
    print("Interpretation:")
    print("  The minimal isotropic conservative quadrupole module is not guessed once")
    print("  the higher conservative branch is read as grouped-P2 + static geometry.")
    print("  It is forced. The old support-side loading ratio then follows as a direct")
    print("  corollary rather than an extra free parameter.")


if __name__ == "__main__":
    main()
