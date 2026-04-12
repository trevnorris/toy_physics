#!/usr/bin/env python3
"""
Step 2 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Splits the current exact local anomaly law into:
      - the frozen lower-order sharp-boundary ledger,
      - the charge-side local transport residue,
      - the inertia-side local transport residue.
2. Proves that a naive common multiplicative dressing of the whole law
   contaminates orders below O(f^4).
3. Introduces the minimal transport-residue dressing ansatz driven by the
   moving-throat quotient path q(f) = (q_tr, q_nt, q_eta).
4. Shows that the first common correction then has the exact quartic form

       Δ(g/2) = (c3_q * λ_Q + c3_i * λ_I) f^4 + O(f^5),

   where λ_Q and λ_I are the first tangent projections of q(f) into the
   charge-side and inertia-side common-dressing functionals.
5. Solves the quartic matching condition using the Step-1 benchmark residual,
   and extracts the minimal "single common scalar" specialization.

Interpretation
--------------
This step isolates the first genuinely coupled correction that the moving-throat
quotient variables can carry without reopening the already-frozen lower-order
anomaly ledger.
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Step 2 — Common quotient-path ansatz
# ---------------------------------------------------------------------------

def step2_common_path() -> None:
    banner("STEP 2 — MINIMAL COMMON QUOTIENT PATH AND QUARTIC MATCHING")

    f, kappa = sp.symbols("f kappa", positive=True, real=True)
    pi = sp.pi

    # Frozen lower-order sharp-boundary ledger.
    eta0 = sp.Rational(11, 36)
    g_sharp_over_2 = 1 + f - (1 + eta0) * f**2

    # Step-1 transport residue coefficients.
    c3_q = sp.simplify((4 - pi) / (pi**2 * kappa))
    q4 = sp.simplify(4 * (pi - 3) / (pi**3 * kappa))

    c3_i = sp.simplify(sp.Rational(11, 6) * kappa)
    i4 = sp.simplify(-sp.Rational(55, 6) * kappa**2)

    T_q = sp.expand(c3_q * f**3 + q4 * f**4)
    T_i = sp.expand(c3_i * f**3 + i4 * f**4)
    T_loc = sp.expand(T_q + T_i)

    subbanner("II.1 — Exact transport-residue split")
    print("g_sharp(f)/2 =")
    sp.pprint(g_sharp_over_2)
    print("T_q(f)  = charge-side local transport residue =")
    sp.pprint(T_q)
    print("T_i(f)  = inertia-side local transport residue =")
    sp.pprint(T_i)
    print("T_loc(f)= T_q + T_i =")
    sp.pprint(T_loc)

    print("\nCubic coefficients:")
    print("c3_q =")
    sp.pprint(c3_q)
    print("c3_i =")
    sp.pprint(c3_i)
    print("c3_total =")
    sp.pprint(sp.simplify(c3_q + c3_i))

    # -----------------------------------------------------------------------
    # Naive whole-law dressing fails
    # -----------------------------------------------------------------------
    L1 = sp.symbols("L1", real=True)
    naive_multiplier = sp.expand(sp.series(sp.exp(L1 * f), f, 0, 3).removeO())
    g_loc_over_2_series = sp.expand(g_sharp_over_2 + T_loc)

    bad_whole = sp.expand(g_loc_over_2_series * naive_multiplier - g_loc_over_2_series)
    bad_series = sp.expand(sp.series(bad_whole, f, 0, 4).removeO())

    subbanner("II.2 — Why whole-law dressing is not acceptable")
    print("Naive whole-law correction [exp(L1 f) - 1] * g_loc(f)/2 =")
    sp.pprint(bad_series)
    print("\nThis starts below O(f^4), so it reopens already-frozen orders.")

    # Also show that dressing the whole charge factor would do the same.
    Q_whole = sp.expand((1 + f - f**2 + T_q) * naive_multiplier - (1 + f - f**2 + T_q))
    Q_whole_series = sp.expand(sp.series(Q_whole, f, 0, 4).removeO())
    print("\nWhole-charge-factor dressing contamination =")
    sp.pprint(Q_whole_series)

    # -----------------------------------------------------------------------
    # General transport-residue dressing
    # -----------------------------------------------------------------------
    subbanner("II.3 — General common-path dressing of the transport residue")

    # Quotient-path tangent data:
    s_tr, s_nt, s_eta = sp.symbols("s_tr s_nt s_eta", real=True)
    a_tr, a_nt, a_eta = sp.symbols("a_tr a_nt a_eta", real=True)
    b_tr, b_nt, b_eta = sp.symbols("b_tr b_nt b_eta", real=True)

    # Minimal analytic quotient path:
    #     q_i(f) = s_i f + O(f^2)
    qvec = sp.Matrix([s_tr * f, s_nt * f, s_eta * f])

    alpha = sp.Matrix([a_tr, a_nt, a_eta])
    beta = sp.Matrix([b_tr, b_nt, b_eta])

    lambda_Q = sp.simplify((alpha.T * sp.Matrix([s_tr, s_nt, s_eta]))[0])
    lambda_I = sp.simplify((beta.T * sp.Matrix([s_tr, s_nt, s_eta]))[0])

    Lambda_Q = sp.simplify(lambda_Q * f)
    Lambda_I = sp.simplify(lambda_I * f)

    expQ = sp.expand(sp.series(sp.exp(Lambda_Q), f, 0, 3).removeO())
    expI = sp.expand(sp.series(sp.exp(Lambda_I), f, 0, 3).removeO())

    g_common_over_2 = sp.expand(g_sharp_over_2 + expQ * T_q + expI * T_i)
    delta_common = sp.expand(g_common_over_2 - (g_sharp_over_2 + T_loc))
    delta_common_series = sp.expand(sp.series(delta_common, f, 0, 5).removeO())
    a4_common = sp.simplify(delta_common_series.coeff(f, 4))

    print("lambda_Q =")
    sp.pprint(lambda_Q)
    print("lambda_I =")
    sp.pprint(lambda_I)
    print("\nΔ(g/2)_common through O(f^4) =")
    sp.pprint(delta_common_series)
    print("\nQuartic common coefficient a4_common =")
    sp.pprint(a4_common)

    expect_zero(
        "a4_common - (c3_q*lambda_Q + c3_i*lambda_I)",
        a4_common - (c3_q * lambda_Q + c3_i * lambda_I),
    )

    # Order-counting theorem:
    n = sp.symbols("n", integer=True, positive=True)
    print("\nOrder-counting result:")
    print("  Since T_loc(f) = O(f^3), a quotient path q(f) = O(f^n) produces")
    print("  the first common correction at O(f^(3+n)).")
    print("  Therefore the quartic residual can be hit only if the path starts")
    print("  linearly: q_i(f) = s_i f + O(f^2).")

    # -----------------------------------------------------------------------
    # Single-scalar common specialization
    # -----------------------------------------------------------------------
    subbanner("II.4 — Minimal single-scalar common specialization")

    w_tr, w_nt, w_eta = sp.symbols("w_tr w_nt w_eta", real=True)
    common_weight = sp.Matrix([w_tr, w_nt, w_eta])

    Lambda_common = sp.simplify((common_weight.T * qvec)[0])
    lambda_common = sp.simplify((common_weight.T * sp.Matrix([s_tr, s_nt, s_eta]))[0])

    exp_common = sp.expand(sp.series(sp.exp(Lambda_common), f, 0, 3).removeO())
    g_single_over_2 = sp.expand(g_sharp_over_2 + exp_common * T_loc)
    delta_single = sp.expand(g_single_over_2 - (g_sharp_over_2 + T_loc))
    delta_single_series = sp.expand(sp.series(delta_single, f, 0, 5).removeO())
    a4_single = sp.simplify(delta_single_series.coeff(f, 4))

    print("Lambda_common(f) =")
    sp.pprint(Lambda_common)
    print("Δ(g/2)_single through O(f^4) =")
    sp.pprint(delta_single_series)
    print("a4_single =")
    sp.pprint(a4_single)

    expect_zero(
        "a4_single - (c3_total * lambda_common)",
        a4_single - sp.simplify((c3_q + c3_i) * lambda_common),
    )

    # -----------------------------------------------------------------------
    # Numerical benchmark using the Step-1 carry-forward baseline
    # -----------------------------------------------------------------------
    subbanner("II.5 — Numerical quartic matching against the Step-1 residual")

    kappa_atom = sp.Float("1.177746578880", 40)
    # Step-1 frozen carry-forward data:
    #   f = alpha_fs/(2 pi)
    #   g_loc exact = 2.002319304358647956...
    #   target      = 2.00231930436092
    f_atom = sp.Float("0.001161409732093", 40)
    residual_g = sp.Float("2.27204390584705e-12", 40)

    c3_q_num = sp.N(c3_q.subs({kappa: kappa_atom}), 30)
    c3_i_num = sp.N(c3_i.subs({kappa: kappa_atom}), 30)
    c3_total_num = sp.N((c3_q + c3_i).subs({kappa: kappa_atom}), 30)

    a4_resid = sp.N(residual_g / (2 * f_atom**4), 30)

    print(f"c3_q      = {c3_q_num}")
    print(f"c3_i      = {c3_i_num}")
    print(f"c3_total  = {c3_total_num}")
    print(f"a4_resid  = {a4_resid}")

    # General matching line
    lamQ_num, lamI_num = sp.symbols("lamQ_num lamI_num", real=True)
    lamI_line = sp.simplify((a4_resid - c3_q_num * lamQ_num) / c3_i_num)

    print("\nGeneral quartic matching line:")
    print("  c3_q * lambda_Q + c3_i * lambda_I = a4_resid")
    print("  lambda_I(lambda_Q) =")
    sp.pprint(lamI_line)

    lambdaI_only = sp.N(a4_resid / c3_i_num, 30)
    lambdaQ_only = sp.N(a4_resid / c3_q_num, 30)
    lambda_common_num = sp.N(a4_resid / c3_total_num, 30)

    print(f"\nInertia-only match  (lambda_Q = 0): lambda_I = {lambdaI_only}")
    print(f"Charge-only match   (lambda_I = 0): lambda_Q = {lambdaQ_only}")
    print(f"Single-common match (lambda_Q = lambda_I = Lambda_1): Lambda_1 = {lambda_common_num}")

    print("\nInterpretation:")
    print("  1. The quartic residual constrains only one linear combination of the")
    print("     quotient-path tangent data at this order.")
    print("  2. A natural single-common dressing of the entire transport residue")
    print("     fixes the required first tangent scalar to Lambda_1 ≈")
    print(f"     {lambda_common_num}.")
    print("  3. Because c3_i >> c3_q numerically, the quartic match is naturally")
    print("     carried mostly by an inertia-linked common dressing unless the")
    print("     charge-side projection is much larger.")


if __name__ == "__main__":
    step2_common_path()
