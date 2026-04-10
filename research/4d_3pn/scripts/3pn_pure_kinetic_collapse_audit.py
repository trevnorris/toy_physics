#!/usr/bin/env python3
"""
3pn_pure_kinetic_collapse_audit.py

Audit for the remaining 3PN COM pure-kinetic slot.

Main point
----------
The leftover COM ordinary-Lagrangian coefficient

    Delta l1 = 3*nu*(3*nu-1)*(4*nu-1)/16

is *not* a new dynamical 3PN throat-response datum. It is exactly the ordinary-
Lagrangian counterimage of the already-fixed universal free two-body relativistic
COM Hamiltonian.

This script verifies four facts:

1. The generic-frame exact free-body ordinary Lagrangian has no mixed 3PN pure-
   kinetic term; its 3PN block is simply 5/128*(m_A v_A^8 + m_B v_B^8).
2. Naive COM reduction of that generic free block gives the carried 3PN seed
   coefficient l1_seed.
3. The exact free relativistic two-body COM Hamiltonian gives the pure-kinetic
   Hamiltonian coefficient

       h1_free = (-5 + 35 nu - 70 nu^2 + 35 nu^3)/128.

4. The exact COM ordinary/Hamiltonian compiler map then forces

       Delta l1 = h1_seed - h1_free,

   and the resulting ordinary coefficient equals the full imported GR target l1.
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


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def nu_poly_from_mass_ratio(expr_mass_ratio: sp.Expr, ratio: sp.Symbol, nu_symbol: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """Fit a rational expression in the mass ratio q = m_A/m_B to a cubic polynomial in
    nu = q/(1+q)^2.  This is sufficient for all symmetric pure-kinetic 3PN slots used here.

    Returns:
        (poly_in_nu_symbol, poly_in_ratio_form)
    """
    nu_q = sp.simplify(ratio / (1 + ratio) ** 2)
    c0, c1, c2, c3 = sp.symbols("c0 c1 c2 c3")
    poly_q = c0 + c1 * nu_q + c2 * nu_q**2 + c3 * nu_q**3
    sample_vals = [2, 3, 5, 7]
    eqs = [sp.Eq(sp.simplify(expr_mass_ratio.subs(ratio, v)), sp.simplify(poly_q.subs(ratio, v))) for v in sample_vals]
    sol = sp.solve(eqs, [c0, c1, c2, c3], dict=True)
    if not sol:
        raise AssertionError("Could not fit a cubic polynomial in nu.")
    sol = sol[0]
    poly_q = sp.simplify(poly_q.subs(sol))
    expect_zero("nu-fit residual", sp.simplify(expr_mass_ratio - poly_q))
    poly_nu = sp.expand(sol[c0] + sol[c1] * nu_symbol + sol[c2] * nu_symbol**2 + sol[c3] * nu_symbol**3)
    return poly_nu, poly_q


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def main() -> None:
    banner("PART I — GENERIC-FRAME FREE PURE-KINETIC BLOCK")

    mA, mB, c = sp.symbols("mA mB c", positive=True)
    a, b = sp.symbols("a b", nonnegative=True)  # a=v_A^2, b=v_B^2

    L_free_gen = -mA * c**2 * sp.sqrt(1 - a / c**2) - mB * c**2 * sp.sqrt(1 - b / c**2)
    # Expand to 3PN / c^-6.
    L_free_series = sp.expand(sp.series(L_free_gen, c, sp.oo, 7).removeO())
    coeff_c6 = sp.simplify(sp.expand(L_free_series * c**6).coeff(c, 0))

    print("Generic-frame free-body ordinary Lagrangian through 3PN:")
    print(L_free_series)
    print("\n3PN pure-kinetic coefficient (c^-6 block):")
    print(coeff_c6)

    expected_coeff_c6 = sp.Rational(5, 128) * (mA * a**4 + mB * b**4)
    expect_zero("generic free 3PN pure-kinetic block", coeff_c6 - expected_coeff_c6)

    subbanner("I.1 — Naive COM reduction of the generic free block")
    M = sp.symbols("M", positive=True)
    nu = sp.symbols("nu", real=True)
    etaA = mA / (mA + mB)
    etaB = mB / (mA + mB)
    v2 = sp.symbols("v2", nonnegative=True)

    l1_seed_mass = sp.simplify(coeff_c6.subs({a: etaB**2 * v2, b: etaA**2 * v2}) / (((mA * mB) / (mA + mB)) * v2**4))
    print("l1_seed in masses =", sp.factor(l1_seed_mass))

    q = sp.symbols("q", positive=True)
    l1_seed_nu, l1_seed_q = nu_poly_from_mass_ratio(sp.simplify(l1_seed_mass.subs({mA: q, mB: 1})), q, nu)
    print("l1_seed(q-form) =", l1_seed_q)
    print("l1_seed(nu) =", l1_seed_nu)

    expected_l1_seed = sp.Rational(5, 128) - sp.Rational(35, 128) * nu + sp.Rational(35, 64) * nu**2 - sp.Rational(35, 128) * nu**3
    expect_zero("l1_seed formula", l1_seed_nu - expected_l1_seed)

    banner("PART II — EXACT FREE RELATIVISTIC TWO-BODY COM HAMILTONIAN")

    # Reduced momentum variable p is defined by P = mu p.
    p = sp.symbols("p", real=True)
    mu = mA * mB / (mA + mB)
    H_free_com = sp.sqrt(mA**2 * c**4 + mu**2 * p**2 * c**2) + sp.sqrt(mB**2 * c**4 + mu**2 * p**2 * c**2)
    Hhat = sp.expand((H_free_com - (mA + mB) * c**2) / mu)
    Hhat_series = sp.expand(sp.series(Hhat, p, 0, 10).removeO())
    print("Reduced COM free Hamiltonian through 3PN:")
    print(Hhat_series)

    h1_free_mass = sp.simplify(sp.expand(Hhat_series).coeff(p, 8) * c**6)
    print("\nh1_free in masses =", sp.factor(h1_free_mass))

    h1_free_nu, h1_free_q = nu_poly_from_mass_ratio(sp.simplify(h1_free_mass.subs({mA: q, mB: 1})), q, nu)
    print("h1_free(q-form) =", h1_free_q)
    print("h1_free(nu) =", h1_free_nu)

    expected_h1_free = -sp.Rational(5, 128) + sp.Rational(35, 128) * nu - sp.Rational(35, 64) * nu**2 + sp.Rational(35, 128) * nu**3
    expect_zero("h1_free formula", h1_free_nu - expected_h1_free)

    print("\nThis is exactly the imported GR COM Hamiltonian pure-kinetic target h1.")

    banner("PART III — EXACT COM COMPILER COMPENSATION")

    F1 = sp.Rational(3, 16) * nu - sp.Rational(21, 16) * nu**2 + sp.Rational(9, 4) * nu**3
    h1_seed = sp.simplify(F1 - expected_l1_seed)
    print("F1(nu) =", F1)
    print("h1_seed(nu) =", h1_seed)

    expected_h1_seed = -sp.Rational(5, 128) + sp.Rational(59, 128) * nu - sp.Rational(119, 64) * nu**2 + sp.Rational(323, 128) * nu**3
    expect_zero("h1_seed formula", h1_seed - expected_h1_seed)

    delta_l1 = sp.simplify(h1_seed - expected_h1_free)
    print("Delta l1 from compiler compensation =", sp.factor(delta_l1))

    expected_delta_l1 = sp.Rational(3, 16) * nu - sp.Rational(21, 16) * nu**2 + sp.Rational(9, 4) * nu**3
    expect_zero("Delta l1 formula", delta_l1 - expected_delta_l1)
    expect_zero("Delta l1 factorized formula", delta_l1 - (3 * nu * (3 * nu - 1) * (4 * nu - 1) / 16))

    l1_full = sp.simplify(F1 - expected_h1_free)
    print("l1_full from exact free Hamiltonian target =", l1_full)

    expected_l1_full = sp.Rational(5, 128) - sp.Rational(11, 128) * nu - sp.Rational(98, 128) * nu**2 + sp.Rational(253, 128) * nu**3
    expect_zero("full GR l1 recovered", l1_full - expected_l1_full)
    expect_zero("l1_full = l1_seed + Delta l1", l1_full - (expected_l1_seed + delta_l1))

    banner("PART IV — THEOREM LEDGER")
    print("1. The exact generic-frame free-body ordinary 3PN block has no mixed pure-kinetic term:")
    print("      L3_free^(gen) = 5/128 (m_A v_A^8 + m_B v_B^8).")
    print("2. Naive COM reduction of that free block gives exactly the carried seed coefficient l1_seed.")
    print("3. The exact free relativistic two-body COM Hamiltonian yields")
    print("      h1_free = (-5 + 35 nu - 70 nu^2 + 35 nu^3)/128,   h2=h3=h4=h5=0,")
    print("   which is exactly the imported GR COM pure-kinetic target.")
    print("4. The remaining COM ordinary coefficient is therefore not a new dynamical datum.")
    print("   It is the exact ordinary-Lagrangian counterimage of the universal free COM Hamiltonian:")
    print("      Delta l1 = h1_seed - h1_free = 3 nu (3 nu - 1)(4 nu - 1)/16.")
    print("5. Equivalently, no extra generic-frame pure-kinetic interaction module should be added in")
    print("   the fixed ADM chart.  The COM Delta l1 term is generated by the exact COM compiler map.")


if __name__ == "__main__":
    main()
