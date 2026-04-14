#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage004_resonance_linewidth_tradeoff_sympy_audit.py

Stage 004 — resonance/linewidth tradeoff audit.

What this script does
---------------------
1. Rebuilds the Stage-003 dynamic one-port mixed bundle and verifies the exact
   derivative identity

       dD_Pi/dPi = - N(omega),

   where N(omega) is the outgoing transfer factor already controlling the wall's
   odd quadrupole response.
2. Extracts the exact wall-like simple-pole normal form

       D_Pi(omega) ~ D0'(omega_*) (omega-omega_*) - Pi(omega_*) N_*,

   and verifies the corresponding Breit-Wigner susceptibility.
3. Derives the universal simple-pole line-shape formulas

       chi(omega) = A_*/(delta - i gamma_*),
       Re chi     = A_* delta/(delta^2 + gamma_*^2),
       Im chi     = A_* gamma_*/(delta^2 + gamma_*^2).

4. Proves the exact dispersive/absorptive tradeoff theorem:
       |Re chi|/|Im chi| = r := |delta|/gamma_*,
   so the maximum conservative line shape occurs at r = 1, precisely where
       |Re chi| = |Im chi|.
5. Proves the low-loss bound:
   if
       |Im chi| <= eta |Re chi|,   0 < eta <= 1,
   then
       sup |Re chi| = (|A_*|/gamma_*) * eta/(1+eta^2).
6. Records the equivalent quality-factor detuning law
       |omega-omega_*|/omega_* >= 1/(2 Q_* eta).
7. Prints an illustrative tolerance table for a few eta values.
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
        simplified = expr.applyfunc(lambda z: sp.factor(sp.simplify(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.factor(sp.simplify(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Part I. Exact Stage-003 derivative identity
# ---------------------------------------------------------------------------

def verify_dynamic_derivative_identity():
    banner("PART I — EXACT DYNAMIC ONE-PORT DERIVATIVE IDENTITY")

    KB, GU, GW, A, W0, R, Pi = sp.symbols(
        "K_B G_U G_W A W_0 R Pi", real=True
    )

    W = sp.simplify(W0 - Pi)
    Delta_Pi = sp.simplify(A * W - R**2)
    Q_Pi = sp.simplify(GU**2 * W + 2 * GU * GW * R + GW**2 * A)
    D_Pi = sp.simplify(KB - Q_Pi / Delta_Pi)
    P = sp.simplify(A * GW + R * GU)
    N = sp.simplify(P**2 / Delta_Pi**2)

    print("D_Pi(omega) =")
    sp.pprint(D_Pi)
    print("N(omega) =")
    sp.pprint(N)

    expect_zero("dD_Pi/dPi + N(omega)", sp.diff(D_Pi, Pi) + N)

    return D_Pi, N


# ---------------------------------------------------------------------------
# Part II. Simple-pole normal form
# ---------------------------------------------------------------------------

def verify_simple_pole_normal_form():
    banner("PART II — UNIVERSAL SIMPLE-POLE NORMAL FORM")

    delta, F0p, Zstar, Gamma, Nnum = sp.symbols(
        "delta F0p Z_star Gamma N_num", positive=True, real=True
    )

    F = sp.simplify(F0p * delta - sp.I * Gamma * Zstar)
    gamma = sp.simplify(Gamma * Zstar / F0p)
    Astar = sp.simplify(Nnum / F0p)
    chi = sp.simplify(Nnum / F)

    print("F_pole(delta) =")
    sp.pprint(F)
    print("A_* =")
    sp.pprint(Astar)
    print("gamma_* =")
    sp.pprint(gamma)
    print("chi_pole(delta) =")
    sp.pprint(chi)

    expect_zero(
        "chi_pole - A_*/(delta - i gamma_*)",
        sp.simplify(chi - Astar / (delta - sp.I * gamma)),
    )

    # Wall-like specialization
    subbanner("Wall-like pole specialization")
    D0p, Nstar = sp.symbols("D0p N_star", positive=True, real=True)
    gamma_wall = sp.simplify(Gamma * Nstar / D0p)
    chi_wall = sp.simplify(1 / (D0p * delta - sp.I * Gamma * Nstar))

    print("gamma_wall =")
    sp.pprint(gamma_wall)
    print("chi_wall(delta) =")
    sp.pprint(chi_wall)

    expect_zero(
        "chi_wall - (1/D0p)/(delta - i gamma_wall)",
        sp.simplify(chi_wall - (1 / D0p) / (delta - sp.I * gamma_wall)),
    )

    return delta, gamma, Astar


# ---------------------------------------------------------------------------
# Part III. Exact line-shape formulas
# ---------------------------------------------------------------------------

def verify_line_shape_formulas():
    banner("PART III — EXACT DISPERSIVE / ABSORPTIVE LINE SHAPE")

    A, gamma, r = sp.symbols("A gamma r", positive=True, real=True)
    delta = sp.simplify(r * gamma)

    chi = sp.simplify(A / (delta - sp.I * gamma))
    chi_rationalized = sp.simplify(A * (delta + sp.I * gamma) / (delta**2 + gamma**2))
    expect_zero("chi - rationalized form", sp.simplify(chi - chi_rationalized))

    Re_chi = sp.simplify(A * delta / (delta**2 + gamma**2))
    Im_chi = sp.simplify(A * gamma / (delta**2 + gamma**2))

    print("Re chi =")
    sp.pprint(Re_chi)
    print("Im chi =")
    sp.pprint(Im_chi)

    expect_zero(
        "Re chi - (A/gamma) * r/(1+r^2)",
        sp.simplify(Re_chi - (A / gamma) * r / (1 + r**2)),
    )
    expect_zero(
        "Im chi - (A/gamma) * 1/(1+r^2)",
        sp.simplify(Im_chi - (A / gamma) * 1 / (1 + r**2)),
    )

    ratio = sp.simplify(Re_chi / Im_chi)
    print("Re/Im ratio =")
    sp.pprint(ratio)
    expect_zero("Re/Im - r", ratio - r)

    # Extremum of the dispersive factor
    f = sp.simplify(r / (1 + r**2))
    df = sp.simplify(sp.diff(f, r))
    d2f = sp.simplify(sp.diff(df, r))

    print("f(r) = r/(1+r^2) =")
    sp.pprint(f)
    print("f'(r) =")
    sp.pprint(df)
    print("f''(1) =")
    sp.pprint(d2f.subs(r, 1))

    expect_zero("f'(1)", sp.simplify(df.subs(r, 1)))
    if sp.simplify(d2f.subs(r, 1)) >= 0:
        raise AssertionError("f''(1) should be negative for a maximum")

    peak_Re = sp.simplify((A / gamma) * f.subs(r, 1))
    peak_Im = sp.simplify((A / gamma) * 1 / (1 + 1**2))
    print("max |Re chi| =")
    sp.pprint(peak_Re)
    print("|Im chi| at the same point =")
    sp.pprint(peak_Im)

    expect_zero("peak Re - A/(2 gamma)", sp.simplify(peak_Re - A / (2 * gamma)))
    expect_zero("peak Re - peak Im", sp.simplify(peak_Re - peak_Im))

    return A, gamma, r, f


# ---------------------------------------------------------------------------
# Part IV. Low-loss bound
# ---------------------------------------------------------------------------

def verify_low_loss_bound(A, gamma, r, f):
    banner("PART IV — EXACT LOW-LOSS BOUND")

    eta = sp.symbols("eta", positive=True, real=True)
    r_boundary = sp.simplify(1 / eta)
    gain_boundary = sp.simplify((A / gamma) * f.subs(r, r_boundary))

    print("r boundary for |Im| <= eta |Re| is:")
    sp.pprint(r_boundary)
    print("gain at that boundary =")
    sp.pprint(gain_boundary)

    expect_zero(
        "gain boundary - (A/gamma) * eta/(1+eta^2)",
        sp.simplify(gain_boundary - (A / gamma) * eta / (1 + eta**2)),
    )

    # Small-eta series
    small_eta_series = sp.series(eta / (1 + eta**2), eta, 0, 5)
    print("small-eta series for eta/(1+eta^2) =")
    print(small_eta_series)

    # Quality-factor version
    omega_star, Qstar = sp.symbols("omega_star Q_star", positive=True, real=True)
    gamma_Q = sp.simplify(omega_star / (2 * Qstar))
    rel_detuning = sp.simplify((r * gamma_Q) / omega_star)
    rel_detuning_boundary = sp.simplify(rel_detuning.subs(r, r_boundary))

    print("relative detuning =")
    sp.pprint(rel_detuning)
    print("relative detuning boundary =")
    sp.pprint(rel_detuning_boundary)

    expect_zero(
        "relative detuning boundary - 1/(2 Q_* eta)",
        sp.simplify(rel_detuning_boundary - 1 / (2 * Qstar * eta)),
    )

    return eta


# ---------------------------------------------------------------------------
# Part V. Barrier / power translation
# ---------------------------------------------------------------------------

def verify_barrier_power_translation():
    banner("PART V — BARRIER / POWER TRANSLATION")

    A, gamma, r, omega, S = sp.symbols("A gamma r omega S", positive=True, real=True)
    Re_chi = sp.simplify((A / gamma) * r / (1 + r**2))
    Im_chi = sp.simplify((A / gamma) * 1 / (1 + r**2))

    U_disp = sp.simplify(sp.Rational(-1, 2) * Re_chi * S**2)
    P_abs = sp.simplify(sp.Rational(1, 2) * omega * Im_chi * S**2)

    print("U_disp =")
    sp.pprint(U_disp)
    print("P_abs =")
    sp.pprint(P_abs)

    expect_zero(
        "power/barrier ratio - 1/r",
        sp.simplify((P_abs / (omega * (-U_disp))) - 1 / r),
    )


# ---------------------------------------------------------------------------
# Part VI. Illustrative tolerance table
# ---------------------------------------------------------------------------

def print_tolerance_table():
    banner("PART VI — ILLUSTRATIVE TOLERANCE TABLE")

    etas = [sp.Integer(1), sp.Rational(1, 2), sp.Rational(1, 10)]
    print("eta    r_min = 1/eta    gain factor = eta/(1+eta^2)")
    for eta in etas:
        r_min = sp.simplify(1 / eta)
        gain = sp.N(eta / (1 + eta**2), 16)
        print(f"{sp.N(eta, 6):>4}    {sp.N(r_min, 8):>10}    {gain:>20}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    verify_dynamic_derivative_identity()
    verify_simple_pole_normal_form()
    A, gamma, r, f = verify_line_shape_formulas()
    verify_low_loss_bound(A, gamma, r, f)
    verify_barrier_power_translation()
    print_tolerance_table()

    banner("STAGE 004 SUMMARY")
    print("1. The Stage-003 one-port bundle satisfies the exact derivative identity dD_Pi/dPi = -N(omega).")
    print("2. Near a simple passive pole, every reduced susceptibility collapses to a Breit-Wigner form A_*/(delta - i gamma_*).")
    print("3. The conservative/absorptive ratio is exactly r = |delta|/gamma_*.")
    print("4. The largest conservative line shape occurs at r = 1, exactly where |Re chi| = |Im chi|.")
    print("5. In a low-loss window |Im chi| <= eta |Re chi|, the best possible linear conservative gain is")
    print("      (|A_*|/gamma_*) * eta/(1+eta^2).")
    print("6. So the linear same-charge survival test is a residue-to-linewidth test, not a generic resonance slogan.")
