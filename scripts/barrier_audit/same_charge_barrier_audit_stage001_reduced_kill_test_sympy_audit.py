
#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage001_reduced_kill_test.py

Stage 001 — reduced same-charge barrier audit.

What this script does
---------------------
1. Verifies the exact 5PN isotropic one-pole and prefactor identities used as
   the admissibility filter.
2. Proves the elementary no-go for Coulomb + inverse-sixth attraction:
      V(x) = 1/x - beta/x^6
   can lower the barrier but cannot erase it by itself.
3. Derives the core-resolved quadrupole small-x coefficient implied by the
   Stage-B.6 tidal map asymptotics.
4. Records the exact threshold formula for a generic mixed/localization
   attraction on the admissible interval x >= 1.
5. Runs one illustrative numerical scan for the exploratory family

      G_mix(x; ell, p) = exp(-x/ell) / x^p

   and reports:
      - the barrier height of the known potential,
      - the finite threshold gamma_crit(epsilon),
      - and the monotone reduction of a reduced WKB action as gamma increases.

Important scope note
--------------------
The numerical family is only an exploratory placeholder.  It is not claimed to be
the actual PDE-derived mixed-sector kernel.  The real mixed kernel must eventually
be computed from the moving-throat bundle.
"""

from __future__ import annotations

import math
from typing import Callable, Optional, Sequence, Tuple

import mpmath as mp
import numpy as np
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
# Part I. Exact 5PN isotropic algebra used as admissibility filter
# ---------------------------------------------------------------------------

def verify_isotropic_5pn_filter() -> None:
    banner("PART I — EXACT 5PN ISOTROPIC FILTER")

    K, M, B0, B2, B4, Z0, Z2, Z4, N0, N2, N4 = sp.symbols(
        "K M B0 B2 B4 Z0 Z2 Z4 N0 N2 N4", real=True
    )

    D0 = K - B0 - Z0
    D2 = -(M + B2 + Z2)
    D4 = -(B4 + Z4)

    u2 = sp.simplify(-D2 / D0)
    u4 = sp.simplify((D2**2 - D0 * D4) / D0**2)
    P0 = sp.simplify(N0 / D0)
    P2 = sp.simplify((D0 * N2 - 2 * D2 * N0) / D0**2)
    P4 = sp.simplify((D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3)

    print("u2 =")
    sp.pprint(u2)
    print("u4 =")
    sp.pprint(u4)
    print("P0 =")
    sp.pprint(P0)
    print("P2 =")
    sp.pprint(P2)
    print("P4 =")
    sp.pprint(P4)

    defect = sp.simplify(u4 - 4 * u2**2)
    target_defect = sp.simplify((D0 * (B4 + Z4) - 3 * (M + B2 + Z2) ** 2) / D0**2)
    expect_zero("one-pole identity", defect - target_defect)

    const_prefactor_N2 = sp.simplify(2 * D2 * N0 / D0)
    const_prefactor_N4 = sp.simplify((2 * D0 * (D2 * N2 + D4 * N0) - 3 * D2**2 * N0) / D0**2)

    print("Constant-prefactor conditions:")
    print("N2 =")
    sp.pprint(const_prefactor_N2)
    print("N4 =")
    sp.pprint(const_prefactor_N4)


# ---------------------------------------------------------------------------
# Part II. Exact inverse-sixth no-go
# ---------------------------------------------------------------------------

def verify_inverse_sixth_nogo() -> None:
    banner("PART II — COULOMB + INVERSE-SIXTH CANNOT ERASE THE BARRIER")

    x, beta = sp.symbols("x beta", positive=True, real=True)

    V = 1 / x - beta / x**6
    dV = sp.diff(V, x)
    x_star = sp.simplify((6 * beta) ** sp.Rational(1, 5))
    V_star = sp.simplify(V.subs(x, x_star))

    print("V(x) =")
    sp.pprint(V)
    print("dV/dx =")
    sp.pprint(dV)
    print("x_* =")
    sp.pprint(x_star)
    print("V(x_*) =")
    sp.pprint(V_star)

    expect_zero("stationary equation at x_*", sp.simplify(dV.subs(x, x_star)))
    expect_zero("positive-hump identity", sp.simplify(V_star - sp.Rational(5, 6) / x_star))


# ---------------------------------------------------------------------------
# Part III. Core-resolved quadrupole asymptotics
# ---------------------------------------------------------------------------

def verify_core_resolved_asymptotics() -> None:
    banner("PART III — CORE-RESOLVED QUADRUPOLE ASYMPTOTICS")

    x = sp.symbols("x", positive=True, real=True)
    c0 = sp.simplify(sp.Rational(2, 225) / sp.pi)

    # G6(x) chosen so that G6 ~ x^{-6} for x >> 1 and G6 ~ c0 x^4 for x << 1
    G6 = x**4 / (x**10 + 1 / c0)

    G6_small = sp.simplify(sp.series(G6, x, 0, 5).removeO())
    G6_large = sp.simplify(sp.limit(x**6 * G6, x, sp.oo))

    print("c0 =")
    sp.pprint(c0)
    print("G6(x) =")
    sp.pprint(G6)
    print("small-x leading term =")
    sp.pprint(G6_small)
    print("large-x coefficient lim x^6 G6(x) =")
    sp.pprint(G6_large)

    expect_zero("small-x coefficient check", sp.simplify(G6_small - c0 * x**4))
    expect_zero("large-x coefficient check", sp.simplify(G6_large - 1))


# ---------------------------------------------------------------------------
# Part IV. Generic threshold formula and WKB monotonicity
# ---------------------------------------------------------------------------

def verify_generic_threshold_and_wkb() -> None:
    banner("PART IV — GENERIC MIXED-SECTOR THRESHOLD AND WKB MONOTONICITY")

    V0, eps, gamma, G = sp.symbols("V0 eps gamma G", positive=True, real=True)
    integrand = sp.sqrt(V0 - gamma * G - eps)
    dgamma = sp.simplify(sp.diff(integrand, gamma))

    print("sqrt integrand =")
    sp.pprint(integrand)
    print("d/dgamma of integrand =")
    sp.pprint(dgamma)
    expect_zero(
        "integrand derivative identity",
        sp.simplify(dgamma + G / (2 * sp.sqrt(V0 - gamma * G - eps))),
    )

    print("\nExact threshold formula on x >= 1 at fixed reduced energy eps:")
    print("gamma_crit(eps) = sup_{x >= 1, V_known(x) > eps} (V_known(x) - eps) / G_mix(x)")


# ---------------------------------------------------------------------------
# Part V. Illustrative numerical family
# ---------------------------------------------------------------------------

def G6_core_numeric(x: float) -> float:
    c0 = 2.0 / (225.0 * math.pi)
    return x**4 / (x**10 + 1.0 / c0)

def known_potential(
    x: float,
    *,
    kappa: float,
    alpha6: float,
    alpha2: float,
    use_core_resolved_quadrupole: bool = True,
) -> float:
    rep = (1.0 + 0.5 * math.exp(-2.0 * kappa * x)) / x
    quad = alpha6 * (G6_core_numeric(x) if use_core_resolved_quadrupole else x**-6)
    scalar = alpha2 * math.exp(-4.0 * kappa * x) / x**2
    return rep - 3.0 * quad - scalar

def mixed_kernel(x: float, *, ell_mix: float, p_mix: float) -> float:
    return math.exp(-x / ell_mix) / (x**p_mix)

def total_potential(
    x: float,
    *,
    gamma_mix: float,
    kappa: float,
    alpha6: float,
    alpha2: float,
    ell_mix: float,
    p_mix: float,
    use_core_resolved_quadrupole: bool = True,
) -> float:
    return known_potential(
        x,
        kappa=kappa,
        alpha6=alpha6,
        alpha2=alpha2,
        use_core_resolved_quadrupole=use_core_resolved_quadrupole,
    ) - gamma_mix * mixed_kernel(x, ell_mix=ell_mix, p_mix=p_mix)

def barrier_profile(
    xs: np.ndarray,
    *,
    gamma_mix: float,
    kappa: float,
    alpha6: float,
    alpha2: float,
    ell_mix: float,
    p_mix: float,
    use_core_resolved_quadrupole: bool = True,
) -> np.ndarray:
    return np.array(
        [
            total_potential(
                float(x),
                gamma_mix=gamma_mix,
                kappa=kappa,
                alpha6=alpha6,
                alpha2=alpha2,
                ell_mix=ell_mix,
                p_mix=p_mix,
                use_core_resolved_quadrupole=use_core_resolved_quadrupole,
            )
            for x in xs
        ]
    )

def gamma_crit_numeric(
    *,
    epsilon: float,
    kappa: float,
    alpha6: float,
    alpha2: float,
    ell_mix: float,
    p_mix: float,
    x_min: float = 1.0,
    x_max: float = 50.0,
    n_grid: int = 50000,
) -> Tuple[float, float]:
    xs = np.linspace(x_min, x_max, n_grid)
    ratios = np.zeros_like(xs)
    for i, x in enumerate(xs):
        vk = known_potential(float(x), kappa=kappa, alpha6=alpha6, alpha2=alpha2)
        gm = mixed_kernel(float(x), ell_mix=ell_mix, p_mix=p_mix)
        if vk > epsilon and gm > 0.0:
            ratios[i] = (vk - epsilon) / gm
    idx = int(np.argmax(ratios))
    return float(xs[idx]), float(ratios[idx])

def outer_turning_point(
    *,
    gamma_mix: float,
    epsilon: float,
    kappa: float,
    alpha6: float,
    alpha2: float,
    ell_mix: float,
    p_mix: float,
    x_min: float = 1.0,
    x_max: float = 50.0,
    n_grid: int = 20000,
) -> Optional[float]:
    xs = np.linspace(x_min, x_max, n_grid)
    vals = barrier_profile(
        xs,
        gamma_mix=gamma_mix,
        kappa=kappa,
        alpha6=alpha6,
        alpha2=alpha2,
        ell_mix=ell_mix,
        p_mix=p_mix,
    ) - epsilon

    positive = vals > 0.0
    if not np.any(positive):
        return None

    last_pos = np.where(positive)[0][-1]
    if last_pos == len(xs) - 1:
        return x_max

    a = float(xs[last_pos])
    b = float(xs[last_pos + 1])
    fa = float(vals[last_pos])
    fb = float(vals[last_pos + 1])

    for _ in range(80):
        m = 0.5 * (a + b)
        fm = total_potential(
            m,
            gamma_mix=gamma_mix,
            kappa=kappa,
            alpha6=alpha6,
            alpha2=alpha2,
            ell_mix=ell_mix,
            p_mix=p_mix,
        ) - epsilon
        if fa * fm > 0.0:
            a, fa = m, fm
        else:
            b, fb = m, fm

    return 0.5 * (a + b)

def reduced_wkb_action(
    *,
    gamma_mix: float,
    epsilon: float,
    kappa: float,
    alpha6: float,
    alpha2: float,
    ell_mix: float,
    p_mix: float,
) -> float:
    x_out = outer_turning_point(
        gamma_mix=gamma_mix,
        epsilon=epsilon,
        kappa=kappa,
        alpha6=alpha6,
        alpha2=alpha2,
        ell_mix=ell_mix,
        p_mix=p_mix,
    )
    if x_out is None:
        return 0.0

    def integrand(x: float) -> float:
        value = total_potential(
            x,
            gamma_mix=gamma_mix,
            kappa=kappa,
            alpha6=alpha6,
            alpha2=alpha2,
            ell_mix=ell_mix,
            p_mix=p_mix,
        ) - epsilon
        return math.sqrt(max(value, 0.0))

    return float(mp.quad(lambda z: integrand(float(z)), [1.0, x_out]))

def run_illustrative_scan() -> None:
    banner("PART V — ILLUSTRATIVE NUMERICAL SCAN (EXPLORATORY ONLY)")

    # Illustrative, user-editable dimensionless parameters
    epsilon = 0.10
    kappa = 0.50
    alpha6 = 0.08
    alpha2 = 0.01
    ell_mix = 2.0
    p_mix = 1.0

    print("Illustrative reduced parameters:")
    print(f"epsilon = {epsilon}")
    print(f"kappa   = {kappa}")
    print(f"alpha6  = {alpha6}")
    print(f"alpha2  = {alpha2}")
    print(f"ell_mix = {ell_mix}")
    print(f"p_mix   = {p_mix}")

    xs = np.linspace(1.0, 20.0, 5000)
    V_known_vals = barrier_profile(
        xs,
        gamma_mix=0.0,
        kappa=kappa,
        alpha6=alpha6,
        alpha2=alpha2,
        ell_mix=ell_mix,
        p_mix=p_mix,
    )
    idx_max = int(np.argmax(V_known_vals))
    x_max = float(xs[idx_max])
    V_max = float(V_known_vals[idx_max])

    print("\nKnown barrier profile on x >= 1:")
    print(f"max_x = {x_max:.6f}")
    print(f"max_V = {V_max:.12f}")
    print(f"V(1)  = {V_known_vals[0]:.12f}")

    xcrit, gcrit = gamma_crit_numeric(
        epsilon=epsilon,
        kappa=kappa,
        alpha6=alpha6,
        alpha2=alpha2,
        ell_mix=ell_mix,
        p_mix=p_mix,
    )
    print("\nBarrier-removal threshold for the exploratory mixed kernel:")
    print(f"x_crit       = {xcrit:.12f}")
    print(f"gamma_crit   = {gcrit:.12f}")

    print("\nReduced WKB action decreases monotonically as gamma_mix grows:")
    test_gammas: Sequence[float] = (0.0, 2.0, 4.0, 6.0, 8.0, 10.0)
    prev_action: Optional[float] = None
    for gamma_mix in test_gammas:
        x_out = outer_turning_point(
            gamma_mix=gamma_mix,
            epsilon=epsilon,
            kappa=kappa,
            alpha6=alpha6,
            alpha2=alpha2,
            ell_mix=ell_mix,
            p_mix=p_mix,
        )
        action = reduced_wkb_action(
            gamma_mix=gamma_mix,
            epsilon=epsilon,
            kappa=kappa,
            alpha6=alpha6,
            alpha2=alpha2,
            ell_mix=ell_mix,
            p_mix=p_mix,
        )
        print(f"gamma_mix = {gamma_mix:>4.1f} | x_out = {x_out!s:>18} | S_red = {action:.12f}")
        if prev_action is not None and action > prev_action + 1e-10:
            raise AssertionError("Reduced WKB action should not increase for this attractive family.")
        prev_action = action

    print("\nInterpretation:")
    print("  - The known repulsive baseline survives on the admissible interval x >= 1.")
    print("  - The exploratory mixed attraction can reduce the barrier strongly.")
    print("  - But barrier removal requires a finite threshold gamma_mix.")
    print("  - The real theorem question is whether the PDE-derived mixed kernel can")
    print("    produce that strength while staying on the 5PN admissible branch.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    verify_isotropic_5pn_filter()
    verify_inverse_sixth_nogo()
    verify_core_resolved_asymptotics()
    verify_generic_threshold_and_wkb()
    run_illustrative_scan()

if __name__ == "__main__":
    main()
