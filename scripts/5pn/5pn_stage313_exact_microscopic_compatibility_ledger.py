#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage312_314_common import *

"""
Stage 313 — exact microscopic compatibility ledger.

This stage converts the monomial-rigidity statement of Stage 312 into the direct
three-equation microscopic compatibility system for the grouped weak-axisymmetric
log drifts
    (lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1),
and solves it for the three dependent drifts
    (tau_1, kappa_eta, mu_1).

That is the smallest exact linear zero-defect ledger currently available inside
 the coherent reduced hierarchy.
"""

if __name__ == "__main__":
    banner("STAGE 313 — EXACT MICROSCOPIC COMPATIBILITY LEDGER")

    chi0_star, deltaU_star = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
    epsW_star, eps_star = sp.symbols("epsilonW_star epsilon_star", positive=True, real=True)

    lambda_1, c_1, gamma_1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
    kappa_U, kappa_eta, kappa_W = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
    mu_1, tau_1 = sp.symbols("mu_1 tau_1", real=True)

    slips = microscopic_slippages(lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1)
    branch = branch_adapted_slippages(
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
        slips["Sigma_chi"],
        slips["Sigma_delta"],
        slips["Sigma_eta"],
        slips["Sigma_Z"],
        slips["Sigma_eps"],
    )

    subbanner("I. Exact three-equation microscopic compatibility ledger")
    eq_tr = sp.Eq(branch["Sigma_tr"], 0)
    eq_nt = sp.Eq(branch["Sigma_nt"], 0)
    eq_eta = sp.Eq(branch["Sigma_eta"], 0)

    print("Tracking compatibility:")
    sp.pprint(eq_tr)
    print("Nontracking compatibility:")
    sp.pprint(eq_nt)
    print("Dressing compatibility:")
    sp.pprint(eq_eta)

    subbanner("II. Exact solve for the dependent drifts")
    sol = compatibility_ledger(
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
        lambda_1,
        c_1,
        gamma_1,
        kappa_U,
        kappa_W,
    )

    print("tau_1 =")
    sp.pprint(sol["tau_1"])
    print("kappa_eta =")
    sp.pprint(sol["kappa_eta"])
    print("mu_1 =")
    sp.pprint(sol["mu_1"])

    solved_slips = microscopic_slippages(
        lambda_1,
        c_1,
        gamma_1,
        kappa_U,
        sol["kappa_eta"],
        kappa_W,
        sol["mu_1"],
        sol["tau_1"],
    )
    solved_branch = branch_adapted_slippages(
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
        solved_slips["Sigma_chi"],
        solved_slips["Sigma_delta"],
        solved_slips["Sigma_eta"],
        solved_slips["Sigma_Z"],
        solved_slips["Sigma_eps"],
    )

    expect_zero("Sigma_tr on the compatibility branch", solved_branch["Sigma_tr"])
    expect_zero("Sigma_nt on the compatibility branch", solved_branch["Sigma_nt"])
    expect_zero("Sigma_eta on the compatibility branch", solved_branch["Sigma_eta"])

    subbanner("III. Explicit nontracking compatibility after eliminating tau_1 and kappa_eta")
    mu1_expanded = sp.simplify(sp.expand(sol["mu_1"]))
    sp.pprint(mu1_expanded)

    subbanner("IV. Rank statement")
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    Mstar = sp.Matrix([
        [0, 1 + deltaU_star, 1 + deltaU_star, -(2 + chi0_star + deltaU_star), 0, 0, 0, 1 + chi0_star],
        [2 * (1 + E_star), 0, 2 * E_star, F_star - E_star, -1, -(2 + E_star), 1, -F_star],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ])
    print("rank(M_*) =", Mstar.rank())
    print("So the zero-defect branch is codimension 3 inside the microscopic 8-drift space.")

    subbanner("V. Interpretation")
    print("The zero-defect branch is now a direct microscopic rigidity ledger.")
    print("Given the five free grouped drifts")
    print("  (lambda_1, c_1, gamma_1, kappa_U, kappa_W),")
    print("the three dependent drifts")
    print("  (tau_1, kappa_eta, mu_1)")
    print("are fixed exactly by the compatibility equations above.")
