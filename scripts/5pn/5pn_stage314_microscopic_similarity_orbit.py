#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage312_314_common import *

"""
Stage 314 — exact microscopic similarity orbit.

This stage builds the finite five-parameter multiplicative similarity family that
preserves the three direct microscopic monomials
    (C_tr, C_nt, epsilon_eta)
exactly. It then proves that:
    1. the finite orbit preserves the monomials exactly,
    2. its tangent space at the identity is five-dimensional,
    3. and that tangent space is precisely the kernel of the monomial-drift matrix.

So the coherent weak-axisymmetric zero-defect theorem becomes a geometric orbit
statement rather than a loose collection of linear drift identities.
"""

if __name__ == "__main__":
    banner("STAGE 314 — EXACT MICROSCOPIC SIMILARITY ORBIT")

    # Base microscopic variables.
    lambda_W, c_etaU, gamma, K_U, K_etaeff, K_Weff, mu_W, T_U = sp.symbols(
        "lambda_W c_etaU gamma K_U K_etaeff K_Weff mu_W T_U",
        positive=True,
        real=True,
    )
    L, sigma = sp.symbols("L sigma", positive=True, real=True)
    chi0_star, deltaU_star = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
    epsW_star, eps_star = sp.symbols("epsilonW_star epsilon_star", positive=True, real=True)

    # Free orbit parameters in logarithmic form.
    u_lambda, u_c, u_gamma, u_KU, u_KW = sp.symbols(
        "u_lambda u_c u_gamma u_KU u_KW", real=True
    )

    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)

    u_Keta = sp.simplify(2 * u_c - u_KU)
    u_TU = sp.simplify(
        u_KU - (1 + deltaU_star) * (u_gamma + u_c - u_KU) / (1 + chi0_star)
    )
    u_muW = sp.simplify(
        2 * u_c
        - u_KU
        + 2 * u_KW
        - 2 * u_lambda
        - E_star * (2 * u_gamma + 2 * u_lambda - u_KU - u_KW)
        - F_star * (1 + deltaU_star) * (u_gamma + u_c - u_KU) / (1 + chi0_star)
    )

    subbanner("I. Exact finite multiplicative orbit map")
    orbit_map = {
        lambda_W: sp.simplify(sp.exp(u_lambda) * lambda_W),
        c_etaU: sp.simplify(sp.exp(u_c) * c_etaU),
        gamma: sp.simplify(sp.exp(u_gamma) * gamma),
        K_U: sp.simplify(sp.exp(u_KU) * K_U),
        K_etaeff: sp.simplify(sp.exp(u_Keta) * K_etaeff),
        K_Weff: sp.simplify(sp.exp(u_KW) * K_Weff),
        mu_W: sp.simplify(sp.exp(u_muW) * mu_W),
        T_U: sp.simplify(sp.exp(u_TU) * T_U),
    }

    print("u_Keta =")
    sp.pprint(u_Keta)
    print("u_TU =")
    sp.pprint(u_TU)
    print("u_muW =")
    sp.pprint(u_muW)

    subbanner("II. Exact monomial preservation on the finite orbit")
    mono = direct_monomials(
        lambda_W,
        c_etaU,
        gamma,
        K_U,
        K_etaeff,
        K_Weff,
        mu_W,
        T_U,
        L,
        sigma,
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
    )
    mono_orbit = direct_monomials(
        orbit_map[lambda_W],
        orbit_map[c_etaU],
        orbit_map[gamma],
        orbit_map[K_U],
        orbit_map[K_etaeff],
        orbit_map[K_Weff],
        orbit_map[mu_W],
        orbit_map[T_U],
        L,
        sigma,
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
    )

    log_ratio_Ctr = sp.simplify(sp.expand_log(sp.log(mono_orbit["Ctr"] / mono["Ctr"]), force=True))
    log_ratio_Cnt = sp.simplify(sp.expand_log(sp.log(mono_orbit["Cnt"] / mono["Cnt"]), force=True))
    log_ratio_epsEta = sp.simplify(sp.expand_log(sp.log(mono_orbit["epsEta"] / mono["epsEta"]), force=True))
    expect_zero("ln(C_tr orbit/original)", log_ratio_Ctr)
    expect_zero("ln(C_nt orbit/original)", log_ratio_Cnt)
    expect_zero("ln(epsilon_eta orbit/original)", log_ratio_epsEta)

    subbanner("III. Tangent space at the identity")
    free_params = [u_lambda, u_c, u_gamma, u_KU, u_KW]
    log_orbit = sp.Matrix([
        sp.log(orbit_map[lambda_W] / lambda_W),
        sp.log(orbit_map[c_etaU] / c_etaU),
        sp.log(orbit_map[gamma] / gamma),
        sp.log(orbit_map[K_U] / K_U),
        sp.log(orbit_map[K_etaeff] / K_etaeff),
        sp.log(orbit_map[K_Weff] / K_Weff),
        sp.log(orbit_map[mu_W] / mu_W),
        sp.log(orbit_map[T_U] / T_U),
    ])
    Torbit = sp.Matrix.hstack(*[sp.diff(log_orbit, p) for p in free_params])
    print("Tangent-basis matrix T_orbit =")
    sp.pprint(Torbit)
    print("rank(T_orbit) =", Torbit.rank())

    Mstar = sp.Matrix([
        [0, 1 + deltaU_star, 1 + deltaU_star, -(2 + chi0_star + deltaU_star), 0, 0, 0, 1 + chi0_star],
        [2 * (1 + E_star), 0, 2 * E_star, F_star - E_star, -1, -(2 + E_star), 1, -F_star],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ])
    expect_zero("M_* T_orbit", sp.simplify(Mstar * Torbit))
    print("rank(M_*) =", Mstar.rank())
    print("kernel dimension =", Torbit.rows - Mstar.rank())

    subbanner("IV. Exact finite orbit-fibre theorem inside the closure")
    print("The three monomials")
    print("  C_tr,*,   C_nt,*,   epsilon_eta")
    print("are preserved exactly by the five-parameter multiplicative family above.")
    print()
    print("Because the invariant equations solve uniquely for")
    print("  (K_etaeff, mu_W, T_U)")
    print("in terms of the five free multiplicative directions")
    print("  (lambda_W, c_etaU, gamma, K_U, K_Weff),")
    print("the finite zero-defect fibre is exactly the similarity orbit inside this closure.")
    print()
    print("Equivalently: two branch states lie on the same microscopic invariant fibre iff")
    print("their relative multiplicative drift sits in this orbit family, and the tangent")
    print("space of that family is precisely ker(M_*).")
