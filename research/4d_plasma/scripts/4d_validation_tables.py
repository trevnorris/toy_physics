#!/usr/bin/env python3
"""
Numerical validation tables for the 4d_plasma paper.

This script replaces the earlier animation-heavy toy scripts with a deterministic
stdout-only harness tied directly to the manuscript:

1. Gaussian soft-wall / Hermite mode-tower identities
2. Coulomb + Yukawa tower convergence on the brane
3. Source-width mode-content tables for finite localization
4. Open-system leakage identity check
5. Projection-covariance EMF example

The goal is not to solve the full PDE. The goal is to provide readable numerical
tables that directly exercise the formulas the paper actually uses.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np
from numpy.polynomial.hermite import hermgauss, hermval


LAMBDA = 1.0
NORM_POINTS = 160
N_MAX = 18
GH_NODES, GH_WEIGHTS = hermgauss(NORM_POINTS)


def hermite_phys(n: int, y: np.ndarray | float) -> np.ndarray | float:
    coeffs = np.zeros(n + 1)
    coeffs[n] = 1.0
    return hermval(y, coeffs)


def gaussian_z(w: np.ndarray | float, lam: float = LAMBDA) -> np.ndarray | float:
    return np.exp(-(w / lam) ** 2)


def zint_exact(lam: float = LAMBDA) -> float:
    return lam * math.sqrt(math.pi)


def mode_mass_sq(n: int, lam: float = LAMBDA) -> float:
    return 2.0 * n / (lam * lam)


def raw_mode_norm_exact(n: int, lam: float = LAMBDA) -> float:
    return lam * math.sqrt(math.pi) * (2.0**n) * math.factorial(n)


def phi_mode(n: int, w: np.ndarray | float, lam: float = LAMBDA) -> np.ndarray | float:
    return hermite_phys(n, np.asarray(w) / lam) / math.sqrt(raw_mode_norm_exact(n, lam))


def coupling_weight_exact(n: int, lam: float = LAMBDA) -> float:
    if n % 2 == 1:
        return 0.0
    m = n // 2
    return (math.comb(2 * m, m) / (4.0**m)) / (lam * math.sqrt(math.pi))


def gh_integral_with_z(func, lam: float = LAMBDA) -> float:
    values = func(lam * GH_NODES)
    return float(lam * np.dot(GH_WEIGHTS, values))


def raw_mode_norm_numeric(n: int, lam: float = LAMBDA) -> float:
    return gh_integral_with_z(lambda w: hermite_phys(n, w / lam) ** 2, lam)


def modal_source_coeffs(sigma: float, n_max: int = N_MAX, lam: float = LAMBDA) -> np.ndarray:
    def source_profile(w: np.ndarray) -> np.ndarray:
        return np.exp(-0.5 * (w / sigma) ** 2)

    coeffs = []
    for n in range(n_max):
        coeffs.append(
            gh_integral_with_z(lambda w, n=n: phi_mode(n, w, lam) * source_profile(w), lam)
        )
    return np.array(coeffs)


def potential_factor(r: float, n_modes: int, lam: float = LAMBDA) -> float:
    c0 = coupling_weight_exact(0, lam)
    total = 0.0
    for n in range(n_modes):
        cn = coupling_weight_exact(n, lam)
        if cn == 0.0:
            continue
        total += cn * math.exp(-math.sqrt(mode_mass_sq(n, lam)) * r)
    return total / c0


def leakage_identity() -> tuple[float, float, float]:
    grid = np.linspace(-8.0, 8.0, 20001)
    w = grid
    W = np.exp(-(w**2))
    S = (1.0 + w) * np.exp(-0.5 * w**2)
    dS = np.gradient(S, w)
    dW = np.gradient(W, w)
    lhs = -float(np.trapz(W * dS, w))
    rhs = -float(W[-1] * S[-1] - W[0] * S[0]) + float(np.trapz(dW * S, w))
    return lhs, rhs, abs(lhs - rhs)


@dataclass
class CovarianceCase:
    label: str
    a: float
    b: float


def covariance_emf(case: CovarianceCase) -> tuple[float, float, float, float]:
    grid = np.linspace(-8.0, 8.0, 20001)
    w = grid
    W = np.exp(-(w**2)) / math.sqrt(math.pi)
    vx = 1.0 + case.a * w**2
    bz = 1.0 + case.b * w**2
    vbar_x = float(np.trapz(W * vx, w))
    bbar_z = float(np.trapz(W * bz, w))
    avg_vxb_y = float(np.trapz(W * (-(vx * bz)), w))
    ecov_y = -(vbar_x * bbar_z) - avg_vxb_y
    return vbar_x, bbar_z, avg_vxb_y, ecov_y


def print_header(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def print_mode_identity_table(lam: float = LAMBDA) -> bool:
    print_header("1. Gaussian Soft-Wall Identities")
    zint_num = gh_integral_with_z(lambda w: np.ones_like(w), lam)
    zint_ref = zint_exact(lam)
    print(f"Zint numeric = {zint_num:.12f}")
    print(f"Zint exact   = {zint_ref:.12f}")
    print(f"abs error    = {abs(zint_num - zint_ref):.3e}")
    print()
    print(f"{'n':>2} {'m_n^2':>12} {'norm_num':>18} {'norm_exact':>18} {'rel_err':>11} {'c_n':>14}")
    print("-" * 83)

    max_rel_err = 0.0
    odd_ok = True
    for n in range(7):
        norm_num = raw_mode_norm_numeric(n, lam)
        norm_ref = raw_mode_norm_exact(n, lam)
        rel_err = abs(norm_num - norm_ref) / norm_ref
        max_rel_err = max(max_rel_err, rel_err)
        c_n = coupling_weight_exact(n, lam)
        odd_ok = odd_ok and ((n % 2 == 0) or abs(c_n) < 1.0e-15)
        print(f"{n:>2d} {mode_mass_sq(n, lam):>12.6f} {norm_num:>18.10f} {norm_ref:>18.10f} {rel_err:>11.3e} {c_n:>14.8f}")

    pass_flag = abs(zint_num - zint_ref) < 1.0e-12 and max_rel_err < 1.0e-12 and odd_ok
    print("-" * 83)
    print(f"PASS: {'yes' if pass_flag else 'no '}  max_norm_rel_err={max_rel_err:.3e}  odd_mode_couplings_zero={odd_ok}")
    return pass_flag


def print_potential_table(lam: float = LAMBDA) -> bool:
    print_header("2. Static Coulomb + Yukawa Tower")
    radii = [1.0, 2.0, 4.0, 8.0]
    truncations = [1, 5, 9, 17]
    ref_modes = 81
    print(f"{'r/lambda':>8} {'N=1':>12} {'N=5':>12} {'N=9':>12} {'N=17':>12} {'N=81 ref':>12} {'rel_err_17':>12}")
    print("-" * 86)

    good = True
    for r in radii:
        row = [potential_factor(r * lam, n, lam) for n in truncations]
        ref = potential_factor(r * lam, ref_modes, lam)
        rel_err = abs(row[-1] - ref) / abs(ref)
        threshold = 2.0e-3 if r <= 1.0 else 1.0e-5
        good = good and rel_err < threshold
        print(f"{r:>8.2f} {row[0]:>12.8f} {row[1]:>12.8f} {row[2]:>12.8f} {row[3]:>12.8f} {ref:>12.8f} {rel_err:>12.3e}")

    asymptotic_ok = abs(potential_factor(8.0 * lam, 17, lam) - 1.0) < 1.0e-3
    pass_flag = good and asymptotic_ok
    print("-" * 86)
    print(f"PASS: {'yes' if pass_flag else 'no '}  truncation_converged={good}  large_r_zero_mode_limit={asymptotic_ok}")
    return pass_flag


def print_mode_content_table(lam: float = LAMBDA) -> bool:
    print_header("3. Source-Width Mode Content")
    sigmas = [0.25, 0.5, 1.0, 2.0, 4.0]
    print(f"{'sigma/lambda':>12} {'frac_n0':>12} {'frac_n2':>12} {'frac_n4':>12} {'frac_tail':>12} {'odd_leak':>12}")
    print("-" * 74)

    frac_n0 = []
    odd_small = True
    for sigma in sigmas:
        coeffs = modal_source_coeffs(sigma * lam, N_MAX, lam)
        power = coeffs**2
        power /= power.sum()
        frac0 = float(power[0])
        frac2 = float(power[2])
        frac4 = float(power[4])
        odd = float(power[1::2].sum())
        tail = float(power[6:].sum())
        frac_n0.append(frac0)
        odd_small = odd_small and odd < 1.0e-28
        print(f"{sigma:>12.2f} {frac0:>12.8f} {frac2:>12.8f} {frac4:>12.8f} {tail:>12.8f} {odd:>12.3e}")

    monotone = all(a <= b for a, b in zip(frac_n0, frac_n0[1:]))
    pass_flag = monotone and odd_small and frac_n0[-1] > 0.999 and frac_n0[0] < 0.7
    print("-" * 74)
    print(f"PASS: {'yes' if pass_flag else 'no '}  n0_fraction_monotone={monotone}  odd_modes_suppressed={odd_small}")
    return pass_flag


def print_leakage_table() -> bool:
    print_header("4. Open-System Leakage Identity")
    lhs, rhs, err = leakage_identity()
    lhs_label = "lhs = -∫W ∂w S"
    rhs_label = "rhs = -[WS] + ∫W' S"
    err_label = "abs error"
    print(f"{lhs_label:<24} {lhs:>18.12f}")
    print(f"{rhs_label:<24} {rhs:>18.12f}")
    print(f"{err_label:<24} {err:>18.3e}")
    pass_flag = err < 1.0e-8
    print(f"PASS: {'yes' if pass_flag else 'no '}")
    return pass_flag


def print_covariance_table() -> bool:
    print_header("5. Projection-Covariance EMF Example")
    cases = [
        CovarianceCase("zero-mode", 0.0, 0.0),
        CovarianceCase("structured", 0.4, 0.6),
    ]
    print(f"{'case':<12} {'vbar_x':>12} {'Bbar_z':>12} {'avg(vxB)_y':>14} {'Ecov_y':>14}")
    print("-" * 68)

    zero_ok = False
    structured_ok = False
    for case in cases:
        vbar_x, bbar_z, avg_vxb_y, ecov_y = covariance_emf(case)
        if case.label == "zero-mode":
            zero_ok = abs(ecov_y) < 1.0e-12
        else:
            structured_ok = abs(ecov_y) > 1.0e-2
        print(f"{case.label:<12} {vbar_x:>12.8f} {bbar_z:>12.8f} {avg_vxb_y:>14.8f} {ecov_y:>14.8f}")

    pass_flag = zero_ok and structured_ok
    print("-" * 68)
    print(f"PASS: {'yes' if pass_flag else 'no '}  zero_mode_vanishes={zero_ok}  structured_case_nonzero={structured_ok}")
    return pass_flag


def main() -> None:
    print("4d_plasma validation tables")
    print("==========================")
    print(f"lambda = {LAMBDA}")
    print(f"quadrature points = {NORM_POINTS}")
    print(f"retained modes = {N_MAX}")

    checks = [
        print_mode_identity_table(),
        print_potential_table(),
        print_mode_content_table(),
        print_leakage_table(),
        print_covariance_table(),
    ]

    print()
    print("Summary")
    print("-------")
    print(f"PASS: {'yes' if all(checks) else 'no'}")


if __name__ == "__main__":
    main()
