#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 232 — exact support-plane selector on the minimal coherent continuum branch.

What this script does
---------------------
1. Imports the Stage-229/230 support-plane compiler.
2. Evaluates it on the actual minimal coherent continuum observables:
      dlnK         = kappa_eta,
      dlnM         = 0,
      dlnOmegaU    = kappa_U/2,
      dlnchi0      = Sigma_chi,
      dlnZW        = zeta_Z.
3. Derives the exact microscopic support-plane scalar
      S_support = Sigma_sup.
4. Prints exact closed coefficients and their numerical values.

Interpretation
--------------
After Stage 231 identified the shape selector with Sigma_eta, this stage shows
that the remaining overlap-image selector is one explicit support-plane scalar
built from the wall-stiffness drift, the internal-U frequency drift, and the
coherent orbit/shape drifts (Sigma_chi, zeta_Z).
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
    banner("STAGE 232 — COHERENT SUPPORT-PLANE SELECTOR")

    stage230 = _load_stage230_module()
    data = stage230.branch_law_tester_data()
    obs = data["observed"]

    # Minimal coherent-continuum branch observables.
    lam1, c1, gam1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
    kapU, kapEta, kapW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)

    Sigma_chi = sp.simplify(gam1 + c1 - kapU)
    zeta_Z = sp.simplify(2 * lam1 - kapEta - kapW)

    coeff_K = sp.simplify(sp.diff(data["OU_pred_from_support"], obs["dlnK"]))
    coeff_M = sp.simplify(sp.diff(data["OU_pred_from_support"], obs["dlnM"]))
    coeff_chi = sp.simplify(sp.diff(data["OU_pred_from_support"], obs["dlnchi0"]))
    coeff_Z = sp.simplify(sp.diff(data["OU_pred_from_support"], obs["dlnZW"]))
    coeff_epsW = sp.simplify(sp.diff(data["OU_pred_from_support"], obs["dlnepsilon_W"]))

    expect_zero("support-plane has no dlnepsilon_W dependence", coeff_epsW)

    dlnK = kapEta
    dlnM = sp.Integer(0)
    dlnOmegaU = kapU / 2

    OU_pred_actual = sp.simplify(
        data["OU_pred_from_support"].subs(
            {
                obs["dlnK"]: dlnK,
                obs["dlnM"]: dlnM,
                obs["dlnchi0"]: Sigma_chi,
                obs["dlnZW"]: zeta_Z,
            }
        )
    )
    S_support = sp.simplify(dlnOmegaU - OU_pred_actual)

    # Positive coefficient convention for the support-plane scalar.
    cK = coeff_K
    cchi = sp.simplify(-coeff_chi)
    cZ = sp.simplify(-coeff_Z)
    Sigma_sup = sp.simplify(kapU / 2 - cK * kapEta + cchi * Sigma_chi + cZ * zeta_Z)

    subbanner("I. Exact support-plane prediction on the minimal coherent branch")
    print("Predicted dlnOmegaU =")
    sp.pprint(OU_pred_actual)
    print()
    print("S_support =")
    sp.pprint(S_support)

    subbanner("II. Exact microscopic support-plane scalar")
    print("Sigma_sup :=")
    sp.pprint(Sigma_sup)
    expect_zero("S_support - Sigma_sup", sp.simplify(S_support - Sigma_sup))

    subbanner("III. Exact coefficient data")
    print("c_K =")
    sp.pprint(cK)
    print("c_chi =")
    sp.pprint(cchi)
    print("c_Z =")
    sp.pprint(cZ)
    print()
    print("Numerically:")
    print("  c_K   ≈", sp.N(cK, 16))
    print("  c_chi ≈", sp.N(cchi, 16))
    print("  c_Z   ≈", sp.N(cZ, 16))

    subbanner("IV. Structural consequences")
    print("The support-plane selector depends on:")
    print("  - the wall-stiffness drift kappa_eta,")
    print("  - the internal-U frequency drift kappa_U/2,")
    print("  - the coherent interference-ratio drift Sigma_chi,")
    print("  - and the wall-to-mixed overlap drift zeta_Z.")
    print("It is independent of dlnepsilon_W at this stage.")

    banner("STAGE 232 LEDGER")
    print("1. On the minimal coherent continuum branch, the reduced overlap-image support selector")
    print("   is one explicit microscopic support-plane scalar:")
    print("      S_support = Sigma_sup")
    print("      = kappa_U/2 - c_K kappa_eta + c_chi Sigma_chi + c_Z zeta_Z.")
    print("2. The exact coefficients are frozen by the Stage-229 support compiler and the")
    print("   constructive packet-null branch; numerically they are")
    print("      c_K   ≈ 1.95241347905,")
    print("      c_chi ≈ 14.3144730533,")
    print("      c_Z   ≈ 1.99787596774.")
    print("3. So the second reduced selector is now a concrete microscopic support-plane test,")
    print("   not an abstract residual of the overlap image.")
