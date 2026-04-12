#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 231 — exact bridge from the reduced two-selector overlap test to the
microscopic coherent-branch slippages.

What this script does
---------------------
1. Imports the Stage-230 reduced selector formulas.
2. Inserts the actual minimal coherent-continuum branch observables:
      dlnK         = kappa_eta,
      dlnM         = 0,
      dlnOmegaU    = kappa_U/2,
      dlnchi0      = Sigma_chi,
      dlnZW        = zeta_Z,
      dlnepsilon_W = Sigma_epsilon.
3. Proves that the reduced shape selector collapses exactly to
      S_shape = - Sigma_eta,
   where
      Sigma_eta = 2 c_1 - kappa_U - kappa_eta.
4. Shows that the selected-branch dressing mismatch
      R_1 + Xi_1
   is therefore proportional to S_shape.

Interpretation
--------------
This is the first exact bridge between the Stage-228/230 overlap-image selector
language and the Stage-165/166 microscopic coherent-branch slippage language.
"""


def _load_stage230_module():
    stage230_path = Path(__file__).with_name("5pn_stage230_reduced_branch_realizability_tester.py")
    spec = importlib.util.spec_from_file_location(
        "stage230_reduced_branch_realizability_tester", stage230_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load Stage 230 module.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if __name__ == "__main__":
    banner("STAGE 231 — COHERENT SELECTOR BRIDGE")

    stage230 = _load_stage230_module()
    data = stage230.branch_law_tester_data()
    obs = data["observed"]

    # Microscopic coherent-branch drifts from Stage 165.
    lam1, c1, gam1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
    kapU, kapEta, kapW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
    eps_eta = sp.symbols("epsilon_eta", positive=True, real=True)

    Sigma_chi = sp.simplify(gam1 + c1 - kapU)
    zeta_Z = sp.simplify(2 * lam1 - kapEta - kapW)
    Sigma_epsilon = sp.simplify(2 * gam1 + 2 * lam1 - kapU - kapW)
    Sigma_eta = sp.simplify(2 * c1 - kapU - kapEta)

    subs_actual = {
        obs["dlnK"]: kapEta,
        obs["dlnM"]: sp.Integer(0),
        obs["dlnOmegaU"]: kapU / 2,
        obs["dlnchi0"]: Sigma_chi,
        obs["dlnZW"]: zeta_Z,
        obs["dlnepsilon_W"]: Sigma_epsilon,
    }

    S_shape = sp.simplify(data["S_shape"].subs(subs_actual))

    subbanner("I. Exact minimal coherent-branch observable substitution")
    print("dlnK =")
    sp.pprint(kapEta)
    print("dlnM =")
    sp.pprint(sp.Integer(0))
    print("dlnOmegaU =")
    sp.pprint(kapU / 2)
    print("dlnchi0 = Sigma_chi =")
    sp.pprint(Sigma_chi)
    print("dlnZW = zeta_Z =")
    sp.pprint(zeta_Z)
    print("dlnepsilon_W = Sigma_epsilon =")
    sp.pprint(Sigma_epsilon)

    subbanner("II. Exact shape-selector collapse")
    print("S_shape(actual coherent branch) =")
    sp.pprint(S_shape)
    print("Sigma_eta =")
    sp.pprint(Sigma_eta)
    expect_zero("S_shape + Sigma_eta", sp.simplify(S_shape + Sigma_eta))

    subbanner("III. Exact selected-branch dressing mismatch bridge")
    RplusXi = sp.simplify(-eps_eta / (1 - eps_eta) * Sigma_eta)
    print("R_1 + Xi_1 =")
    sp.pprint(RplusXi)
    print("(epsilon_eta/(1-epsilon_eta)) * S_shape =")
    sp.pprint(sp.simplify(eps_eta / (1 - eps_eta) * S_shape))
    expect_zero(
        "R_1 + Xi_1 - epsilon_eta S_shape/(1-epsilon_eta)",
        sp.simplify(RplusXi - eps_eta * S_shape / (1 - eps_eta)),
    )

    subbanner("IV. Structural consequences")
    print("Support-lane variables are absent from S_shape.")
    print("The reduced overlap-image shape selector is exactly the microscopic dressing slippage.")
    print("Therefore")
    print("  S_shape = 0  <=>  Sigma_eta = 0  <=>  R_1 + Xi_1 = 0")
    print("on the physical coherent branch 0 < epsilon_eta < 1.")

    banner("STAGE 231 LEDGER")
    print("1. On the actual minimal coherent continuum branch, the reduced overlap-image")
    print("   shape selector is not abstract. It is exactly")
    print("      S_shape = -Sigma_eta.")
    print("2. The selected-branch dressing mismatch is therefore carried by the same scalar:")
    print("      R_1 + Xi_1 = (epsilon_eta/(1-epsilon_eta)) S_shape.")
    print("3. So the first reduced selector is now fully identified with one microscopic")
    print("   coherent-branch defect coordinate.")
