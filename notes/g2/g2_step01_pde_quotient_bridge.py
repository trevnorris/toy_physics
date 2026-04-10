#!/usr/bin/env python3
"""
Step 1 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Reconstructs the exact moving-throat quotient variables from the PDE handoff:
      - tracking monomial C_tr,*
      - nontracking monomial C_nt,*
      - dressing invariant epsilon_eta
2. Builds the exact monomial-drift matrix M_* and verifies:
      - rank(M_*) = 3
      - nullity(M_*) = 5
3. Derives an explicit right inverse R so that q = M_* δx and M_* R = I_3.
4. Derives the exact tangent similarity-orbit decomposition:
      - five pure similarity directions
      - three genuine quotient directions q = (q_tr, q_nt, q_eta)
5. Rebuilds the current staggered anomaly closure from atom_work.md and
   extracts the exact quartic series coefficient already present in that law.
6. Computes the small benchmark coefficient that would be needed only if one
   compared the measured residual against a cubic truncation.

Interpretation
--------------
This is the algebraic starting point for replacing the old staggered transport
story by a genuinely coupled common charge–inertia layer.  The key result is
that the first common correction should be built from the three quotient
coordinates q = (δln C_tr, δln C_nt, δln ε_eta), not from arbitrary sequential
updates of charge and inertia separately.
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
        simplified = expr.applyfunc(sp.simplify)
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
# Part I — Exact moving-throat quotient bridge
# ---------------------------------------------------------------------------

def moving_throat_bridge() -> None:
    banner("PART I — EXACT MOVING-THROAT QUOTIENT BRIDGE")

    # Reference-branch microscopic parameters.
    chi_s, delta_s = sp.symbols("chi_* delta_*", positive=True, real=True)
    epsW_s, eps_s = sp.symbols("epsilonW_* epsilon_*", real=True)

    # Exact coefficients from the nontracking monomial.
    E_s = sp.simplify(
        (2 * epsW_s / (1 - eps_s))
        * (11 + 9 * delta_s) / (11 * (1 + delta_s))
    )
    F_s = sp.simplify(
        2 * chi_s / (1 + delta_s)
        + 4 * epsW_s * delta_s / (11 * (1 - eps_s) * (1 + delta_s) ** 2)
    )

    print("E_* =")
    sp.pprint(E_s)
    print("F_* =")
    sp.pprint(F_s)

    # Microscopic grouped weak-axisymmetric drift vector.
    lam, c1, gam, kU, keta, kW, mu1, tau1 = sp.symbols(
        "lambda_1 c_1 gamma_1 kappa_U kappa_eta kappa_W mu_1 tau_1",
        real=True,
    )
    dx = sp.Matrix([lam, c1, gam, kU, keta, kW, mu1, tau1])

    # Exact quotient drift vector q = M_* dx.
    qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)
    q = sp.Matrix([qtr, qnt, qeta])

    M = sp.Matrix(
        [
            [0, 1 + delta_s, 1 + delta_s, -(2 + chi_s + delta_s), 0, 0, 0, 1 + chi_s],
            [2 * (1 + E_s), 0, 2 * E_s, F_s - E_s, -1, -(2 + E_s), 1, -F_s],
            [0, 2, 0, -1, -1, 0, 0, 0],
        ]
    )

    subbanner("I.1 — Monomial-drift matrix M_*")
    sp.pprint(M)
    print(f"rank(M_*) = {M.rank()}")
    print(f"nullity(M_*) = {M.shape[1] - M.rank()}")

    # A convenient right inverse: choose the quotient lift with free co-scalings
    # (lambda_1, c_1, gamma_1, kappa_U, kappa_W) set to zero.
    R = sp.Matrix(
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, -1],
            [0, 0, 0],
            [F_s / (1 + chi_s), 1, -1],
            [1 / (1 + chi_s), 0, 0],
        ]
    )

    subbanner("I.2 — Exact right inverse for the quotient directions")
    print("R =")
    sp.pprint(R)
    expect_zero("M_* R - I_3", M * R - sp.eye(3))

    # Tangent similarity-orbit formulas (q=0) using free co-scalings
    # (lambda_1, c_1, gamma_1, kappa_U, kappa_W).
    tau_orbit = sp.simplify(
        kU - (1 + delta_s) / (1 + chi_s) * (gam + c1 - kU)
    )
    keta_orbit = sp.simplify(2 * c1 - kU)
    mu_orbit = sp.simplify(
        2 * c1 - kU + 2 * kW - 2 * lam
        - E_s * (2 * gam + 2 * lam - kU - kW)
        - F_s * (1 + delta_s) / (1 + chi_s) * (gam + c1 - kU)
    )

    # Exact full decomposition with quotient motion turned back on.
    tau_full = sp.simplify(tau_orbit + qtr / (1 + chi_s))
    keta_full = sp.simplify(keta_orbit - qeta)
    mu_full = sp.simplify(mu_orbit + qnt - qeta + F_s * qtr / (1 + chi_s))

    subbanner("I.3 — Exact decomposition into orbit motion + quotient motion")
    print("kappa_eta =")
    sp.pprint(keta_full)
    print("tau_1 =")
    sp.pprint(tau_full)
    print("mu_1 =")
    sp.pprint(mu_full)

    # Verify that these solve q = M dx.
    dx_full = sp.Matrix([lam, c1, gam, kU, keta_full, kW, mu_full, tau_full])
    expect_zero("M_* dx_full - q", M * dx_full - q)

    # Build an explicit kernel basis from the five free co-scalings.
    free_syms = [lam, c1, gam, kU, kW]
    Gcols = []
    for s in free_syms:
        subs = {fs: 0 for fs in free_syms}
        subs[s] = 1
        Gcols.append(
            sp.Matrix(
                [
                    subs[lam],
                    subs[c1],
                    subs[gam],
                    subs[kU],
                    sp.simplify(keta_orbit.subs(subs)),
                    subs[kW],
                    sp.simplify(mu_orbit.subs(subs)),
                    sp.simplify(tau_orbit.subs(subs)),
                ]
            )
        )
    G = sp.Matrix.hstack(*Gcols)

    subbanner("I.4 — Five exact tangent similarity directions")
    print("G =")
    sp.pprint(G)
    expect_zero("M_* G", M * G)

    # Chosen projector onto quotient motion along the similarity directions.
    Pq = sp.simplify(R * M)
    Psim = sp.simplify(sp.eye(8) - Pq)
    subbanner("I.5 — Chosen quotient/similarity projectors")
    expect_zero("P_q^2 - P_q", sp.simplify(Pq * Pq - Pq))
    expect_zero("M_* P_sim", sp.simplify(M * Psim))

    # The physically important one-line conclusion.
    print("\nKey conclusion:")
    print("  Any first common correction should be organized by")
    print("      q = (δln C_tr, δln C_nt, δln ε_eta)")
    print("  because the five free co-scalings sit in ker(M_*).")


# ---------------------------------------------------------------------------
# Part II — Rebuild the current staggered anomaly benchmark
# ---------------------------------------------------------------------------

def staggered_anomaly_benchmark() -> None:
    banner("PART II — CURRENT STAGGERED ANOMALY BENCHMARK")

    f, kappa = sp.symbols("f kappa", positive=True, real=True)
    tau = sp.symbols("tau", positive=True, real=True)
    s = sp.symbols("s", real=True)

    # Collar-width variable used in atom_work:
    #     tau(f) = 1 - sqrt(1-f),   f = 2 tau - tau^2.
    f_of_tau = sp.simplify(2 * tau - tau ** 2)
    tau_of_f = sp.simplify(1 - sp.sqrt(1 - f))

    # Exact charge-side collar mode in the rescaled collar coordinate s ∈ [0,1].
    cbar_tau = sp.simplify(
        (2 / (2 - tau))
        * sp.integrate((1 - tau * s) * sp.cos(sp.pi * s / 2), (s, 0, 1))
    )
    Xi_tau = sp.simplify(
        (2 / (2 - tau))
        * sp.integrate(
            (1 - tau * s) * (1 - tau * s) ** 2
            * (sp.cos(sp.pi * s / 2) - cbar_tau),
            (s, 0, 1),
        )
    )
    Aq_tau = sp.simplify(tau / kappa)

    subbanner("II.1 — Exact charge-side collar integrals")
    print("cbar(tau) =")
    sp.pprint(cbar_tau)
    print("Xi(tau) =")
    sp.pprint(Xi_tau)

    # Exact Q_loc(f) reconstructed from the Appendix F definitions.
    Q_tau = sp.simplify(1 + f_of_tau - f_of_tau ** 2 + 2 * f_of_tau * Aq_tau * Xi_tau)
    Q_exact = sp.simplify(Q_tau.subs(tau, tau_of_f))
    Q_series = sp.expand(sp.series(Q_exact, f, 0, 5).removeO())

    subbanner("II.2 — Exact current staggered charge moment Q_loc(f)")
    print("Q_loc exact =")
    sp.pprint(Q_exact)
    print("Q_loc series through O(f^4) =")
    sp.pprint(Q_series)

    # Exponential-blur inertial factor through the quartic series order.
    Bexp_series = 1 - 6 * kappa * f + 30 * kappa ** 2 * f ** 2 - 120 * kappa ** 3 * f ** 3 + 360 * kappa ** 4 * f ** 4
    eta1_series = sp.Rational(11, 36) * Bexp_series

    g_over_2_series = sp.expand(Q_series - eta1_series * f ** 2)
    a3_staggered = sp.simplify(g_over_2_series.coeff(f, 3))
    a4_staggered = sp.simplify(g_over_2_series.coeff(f, 4))

    subbanner("II.3 — Exact staggered g_loc(f) series through O(f^4)")
    print("g_loc(f)/2 =")
    sp.pprint(g_over_2_series)
    print("\nRaw cubic coefficient a3 (in +a3 f^3 convention) =")
    sp.pprint(a3_staggered)
    print("Raw quartic coefficient a4_staggered (in +a4 f^4 convention) =")
    sp.pprint(a4_staggered)

    # Numerical benchmark values frozen in atom_work.md.
    f_atom = sp.Float("0.001161409732093", 30)
    kappa_atom = sp.Float("1.177746578880", 30)
    g_target = sp.Float("2.00231930436092", 30)

    # Use the exact reconstructed current staggered closure.
    tau_num = sp.simplify(tau_of_f.subs({f: f_atom}))
    cbar_num = sp.N(cbar_tau.subs({tau: tau_num}), 30)
    Xi_num = sp.N(Xi_tau.subs({tau: tau_num}), 30)
    Aq_num = sp.N(Aq_tau.subs({tau: tau_num, kappa: kappa_atom}), 30)
    Q_num = sp.N(Q_exact.subs({f: f_atom, kappa: kappa_atom}), 30)

    d_num = sp.N(kappa_atom * f_atom, 30)
    Bexp_exact_num = sp.N(
        1 - 6 * d_num + 30 * d_num ** 2 - 120 * d_num ** 3 + 360 * d_num ** 4
        - 720 * d_num ** 5 + 720 * d_num ** 6 * (1 - sp.exp(-1 / d_num)),
        30,
    )
    eta_num = sp.N(sp.Rational(11, 36) * Bexp_exact_num, 30)
    g_loc_exact = sp.N(2 * (Q_num - eta_num * f_atom ** 2), 30)
    residual = sp.N(g_target - g_loc_exact, 30)

    # Two different quartic benchmarks:
    #  (i) raw quartic coefficient already present in the exact current staggered law;
    # (ii) extra coefficient that would be needed only if we corrected a cubic truncation.
    a4_staggered_num = sp.N(a4_staggered.subs({kappa: kappa_atom}), 30)
    a4_cubic_benchmark = sp.N(residual / (2 * f_atom ** 4), 30)
    c4_part33_convention = sp.N(-a4_cubic_benchmark, 30)

    subbanner("II.4 — Numerical benchmark")
    print(f"tau(f_atom) = {tau_num}")
    print(f"cbar(f_atom) = {cbar_num}")
    print(f"Xi(f_atom) = {Xi_num}")
    print(f"A_q(f_atom) = {Aq_num}")
    print(f"Q_loc(f_atom) = {Q_num}")
    print(f"B_exp(kappa f) = {Bexp_exact_num}")
    print(f"eta_1(f_atom) = {eta_num}")
    print(f"g_loc exact = {g_loc_exact}")
    print(f"target - g_loc = {residual}")
    print(f"a4_staggered exact-series coefficient = {a4_staggered_num}")
    print(f"a4 needed against cubic truncation only = {a4_cubic_benchmark}")
    print(f"c4 in Part VIII sign convention ( ... - c4 f^4 ) = {c4_part33_convention}")

    print("\nInterpretation:")
    print("  1. The current staggered closure already contains a definite quartic series coefficient.")
    print("  2. The tiny measured residual should therefore be treated as an incremental")
    print("     common-layer correction Δg_common(f) starting at O(f^4), not as the")
    print("     quartic coefficient of the full exact staggered law itself.")
    print("  3. If one works only against the cubic truncation, the raw added quartic")
    print("     coefficient is +0.6238..., while the Part VIII notation would call this")
    print("     c4 = -0.6238... because the term is written as -c4 f^4.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    moving_throat_bridge()
    staggered_anomaly_benchmark()
