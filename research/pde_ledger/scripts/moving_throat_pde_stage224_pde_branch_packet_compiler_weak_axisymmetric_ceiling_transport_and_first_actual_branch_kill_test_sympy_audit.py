import math
from decimal import Decimal, getcontext

import sympy as sp


def assert_close(actual: float, expected: float, tol: float = 1e-12) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{actual} !~= {expected} (tol={tol})")


def main() -> None:
    sp.init_printing()
    getcontext().prec = 50

    print("=== Stage 224 SymPy audit: actual-branch packet compiler and weak-axisymmetric ceiling transport ===")

    # ------------------------------------------------------------------
    # Exact packet-to-lane compiler
    # ------------------------------------------------------------------
    Delta_norm, T_quad, mhat0 = sp.symbols("Delta_norm T_quad mhat0", positive=True)
    aP, bP = sp.symbols("aP bP", real=True)
    Pcrit = sp.symbols("Pcrit", positive=True)

    Pbar = (Delta_norm + T_quad) / mhat0**2
    P20 = Pbar + 4 * aP
    P21 = Pbar - aP + bP
    P22 = Pbar - aP - bP

    inv_bar = sp.simplify((P20 + 2 * P21 + 2 * P22) / 5)
    inv_a = sp.simplify((2 * P20 - P21 - P22) / 10)
    inv_b = sp.simplify((P21 - P22) / 2)

    assert sp.simplify(inv_bar - Pbar) == 0
    assert sp.simplify(inv_a - aP) == 0
    assert sp.simplify(inv_b - bP) == 0

    print("\nVerified exact grouped inverse map:")
    print("P20 =", P20)
    print("P21 =", P21)
    print("P22 =", P22)
    print("Recovered Pbar =", inv_bar)
    print("Recovered aP   =", inv_a)
    print("Recovered bP   =", inv_b)

    # ------------------------------------------------------------------
    # Exact isotropic ceiling compiler
    # ------------------------------------------------------------------
    isotropic_ceiling = sp.simplify(Pbar - Pcrit)
    expected_rhs = sp.simplify(mhat0**2 * Pcrit - T_quad)
    # Direct algebraic equivalence because mhat0**2 > 0.
    assert sp.expand((Pbar - Pcrit) * mhat0**2) == sp.expand(Delta_norm - expected_rhs)

    calibrated_lower_bound = sp.simplify(T_quad / Pcrit)

    print("\nExact isotropic compiler:")
    print("Pbar =", Pbar)
    print("Pbar <= Pcrit  <=>  Delta_norm <=", expected_rhs)
    print("Calibrated branch (Delta_norm = 0):  mhat0^2 >=", calibrated_lower_bound)

    # ------------------------------------------------------------------
    # Exact weak-axisymmetric signature and Xi1 compiler
    # ------------------------------------------------------------------
    eps, Xi1 = sp.symbols("eps Xi1", real=True)
    lam20 = sp.Integer(1)
    lam21 = sp.Rational(1, 2)
    lam22 = -sp.Integer(1)

    P20_wa = sp.expand(Pbar * (1 + eps * lam20 * Xi1))
    P21_wa = sp.expand(Pbar * (1 + eps * lam21 * Xi1))
    P22_wa = sp.expand(Pbar * (1 + eps * lam22 * Xi1))

    aP_wa = sp.simplify((2 * P20_wa - P21_wa - P22_wa) / 10)
    bP_wa = sp.simplify((P21_wa - P22_wa) / 2)

    assert sp.simplify(aP_wa - eps * Pbar * Xi1 / 4) == 0
    assert sp.simplify(bP_wa - 3 * eps * Pbar * Xi1 / 4) == 0
    assert sp.simplify(bP_wa - 3 * aP_wa) == 0

    P20_from_ab = sp.expand(Pbar + 4 * aP_wa)
    P21_from_ab = sp.expand(Pbar - aP_wa + bP_wa)
    P22_from_ab = sp.expand(Pbar - aP_wa - bP_wa)

    assert sp.simplify(P20_from_ab - P20_wa) == 0
    assert sp.simplify(P21_from_ab - P21_wa) == 0
    assert sp.simplify(P22_from_ab - P22_wa) == 0

    print("\nVerified weak-axisymmetric lane law:")
    print("P20 =", P20_wa)
    print("P21 =", P21_wa)
    print("P22 =", P22_wa)
    print("aP  =", aP_wa)
    print("bP  =", bP_wa)
    print("Exact relation: bP = 3 aP")

    # ------------------------------------------------------------------
    # Exact robust ceiling collapse
    # ------------------------------------------------------------------
    z = sp.symbols("z", real=True)
    # Positive branch: z = |eps Xi1| >= 0
    zabs = sp.symbols("zabs", positive=True)

    P20_pos = sp.expand(Pbar * (1 + zabs))
    P21_pos = sp.expand(Pbar * (1 + zabs / 2))
    P22_pos = sp.expand(Pbar * (1 - zabs))
    assert sp.simplify(P20_pos - P21_pos - Pbar * zabs / 2) == 0
    assert sp.simplify(P20_pos - P22_pos - 2 * Pbar * zabs) == 0

    P20_neg = sp.expand(Pbar * (1 - zabs))
    P21_neg = sp.expand(Pbar * (1 - zabs / 2))
    P22_neg = sp.expand(Pbar * (1 + zabs))
    assert sp.simplify(P22_neg - P21_neg - 3 * Pbar * zabs / 2) == 0
    assert sp.simplify(P22_neg - P20_neg - 2 * Pbar * zabs) == 0

    robust_a_form = sp.simplify(Pbar + 4 * sp.Abs(aP_wa))
    robust_xi_form = sp.simplify(Pbar * (1 + sp.Abs(eps * Xi1)))
    # SymPy keeps Abs symbolic; verify by substituting zabs for Abs(eps*Xi1)
    subs_abs = {sp.Abs(eps * Xi1): zabs, sp.Abs(aP_wa): Pbar * zabs / 4}
    assert sp.simplify(robust_a_form.subs(subs_abs) - Pbar * (1 + zabs)) == 0
    assert sp.simplify(robust_xi_form.subs(subs_abs) - Pbar * (1 + zabs)) == 0

    calibrated_weak_bound = sp.simplify(T_quad * (1 + sp.Abs(eps * Xi1)) / Pcrit)

    print("\nVerified robust all-lane ceiling collapse:")
    print("For eps*Xi1 > 0, max lane prefactor =", P20_pos)
    print("For eps*Xi1 < 0, max lane prefactor =", P22_neg)
    print("Common robust form = Pbar * (1 + |eps Xi1|) = Pbar + 4 |aP|")
    print("Calibrated weak-axisymmetric lower bound on mhat0^2:", calibrated_weak_bound)

    # ------------------------------------------------------------------
    # Numerical headroom at the Stage 223 compatibility point
    # ------------------------------------------------------------------
    barP0_compat = Decimal("0.002069792318062885")

    ceilings = {
        "both_10": Decimal("0.0028313316855593175"),
        "one_10": Decimal("0.0035965105896846573"),
        "both_30": Decimal("0.00817339430971383"),
        "one_30": Decimal("0.0116633929790174"),
    }

    budgets = {}
    for key, Pcrit_val in ceilings.items():
        eps_xi_budget = Pcrit_val / barP0_compat - Decimal(1)
        a_budget = (Pcrit_val - barP0_compat) / Decimal(4)
        budgets[key] = (eps_xi_budget, a_budget)

    print("\nExplicit headroom at the Stage 223 compatibility point:")
    for key, (eps_xi_budget, a_budget) in budgets.items():
        print(f"{key}: |eps Xi1| <= {eps_xi_budget} ; |aP| <= {a_budget}")

    assert_close(float(budgets["both_10"][0]), 0.367930328492646)
    assert_close(float(budgets["both_10"][1]), 1.90384841874108e-4)
    assert_close(float(budgets["one_10"][0]), 0.737619063660757)
    assert_close(float(budgets["one_10"][1]), 3.81679567905443e-4)
    assert_close(float(budgets["both_30"][0]), 2.94889585703134)
    assert_close(float(budgets["both_30"][1]), 1.52590049791274e-3)
    assert_close(float(budgets["one_30"][0]), 4.63505472371892)
    assert_close(float(budgets["one_30"][1]), 2.39840016523863e-3)

    print("\nAll Stage 224 symbolic and numerical checks passed.")


if __name__ == "__main__":
    main()
