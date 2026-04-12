#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

_stage224_path = Path(__file__).with_name("5pn_stage224_overlap_image_parameterization.py")
_spec = importlib.util.spec_from_file_location("stage224_overlap_image_parameterization", _stage224_path)
_mod = importlib.util.module_from_spec(_spec)
assert _spec is not None and _spec.loader is not None
_spec.loader.exec_module(_mod)
overlap_image_parameterization = _mod.overlap_image_parameterization
constructive_packet_null_data = _mod.constructive_packet_null_data

"""
Stage 225 — exact similarity-orbit bridge in overlap variables.

What this script does
---------------------
1. Takes the Stage-224 overlap-image parameterization of the constructive packet-null family.
2. Rewrites the direct coherent weak-axisymmetric zero-defect relations entirely in overlap variables:
      Sigma_eta,
      Sigma_tr,
      Sigma_nt.
3. Verifies that the whole support-restored packet-null overlap image lies tangent to the exact
   monomial-preserving similarity orbit on the constructive branch.
4. Shows that the orbit verdict is support-blind in the explicit support pair
      (dlnC, dlnvarpi).
5. Isolates the exact orbit-layer formulas that remain after the packet-null solve.

Interpretation
--------------
The overlap image now splits cleanly into
  - a 2D orbit/shape layer controlled by (dlnchi0, dlnZW),
  - a 2D support layer controlled by (dlnC, dlnvarpi).
The packet-null family is automatically orbit-clean. The remaining PDE question is whether
any actual moving-throat overlap tangent lands inside the full 4D image, not whether it obeys
orbit lock once it is there.
"""


if __name__ == "__main__":
    banner("STAGE 225 — EXACT SIMILARITY-ORBIT BRIDGE IN OVERLAP VARIABLES")

    data = constructive_packet_null_data()
    param = overlap_image_parameterization()

    E_star = sp.simplify(data["E_star"])
    F_star = sp.simplify(data["F_star"])
    chi0 = sp.simplify(data["chi0"])
    deltaU = sp.simplify(data["deltaU"])

    dlnchi0 = param["dlnchi0"]
    dlnZW = param["dlnZW"]
    dlnC = param["dlnC"]
    dlnvarpi = param["dlnvarpi"]
    dlnOmegaW = param["dlnOmegaW"]
    dlnepsilon_eta = param["dlnepsilon_eta"]
    dlnepsilon_W = param["dlnepsilon_W"]
    dlndeltaU = param["dlndeltaU"]

    Sigma_eta = sp.simplify(dlnepsilon_eta)
    Sigma_tr = sp.simplify((1 + chi0) * dlndeltaU + (1 + deltaU) * dlnchi0)
    Sigma_nt = sp.simplify((dlnZW - 2 * dlnOmegaW) + E_star * dlnepsilon_W - F_star * dlndeltaU)

    subbanner("I. Direct monomial zero-defect relations in overlap variables")
    print("Sigma_eta =")
    sp.pprint(Sigma_eta)
    print("Sigma_tr =")
    sp.pprint(Sigma_tr)
    print("Sigma_nt =")
    sp.pprint(Sigma_nt)

    expect_zero("Sigma_eta on the overlap image", Sigma_eta)
    expect_zero("Sigma_tr on the overlap image", Sigma_tr)
    expect_zero("Sigma_nt on the overlap image", Sigma_nt)

    subbanner("II. Constructive branch constants")
    print("chi_0 =")
    sp.pprint(chi0)
    print("delta_U =")
    sp.pprint(deltaU)
    print("E_* =")
    sp.pprint(E_star)
    print("F_* =")
    sp.pprint(F_star)
    print("So the exact orbit-clean overlap laws are")
    print("  dlnepsilon_eta = 0,")
    print("  (1+chi_0) dlndeltaU + (1+delta_U) dlnchi0 = 0,")
    print("  (dlnZW - 2 dlnOmegaW) + E_* dlnepsilon_W - F_* dlndeltaU = 0.")

    subbanner("III. Support-blindness of the orbit verdict")
    expect_zero("d Sigma_eta / dlnC", sp.diff(Sigma_eta, dlnC))
    expect_zero("d Sigma_eta / dlnvarpi", sp.diff(Sigma_eta, dlnvarpi))
    expect_zero("d Sigma_tr / dlnC", sp.diff(Sigma_tr, dlnC))
    expect_zero("d Sigma_tr / dlnvarpi", sp.diff(Sigma_tr, dlnvarpi))
    expect_zero("d Sigma_nt / dlnC", sp.diff(Sigma_nt, dlnC))
    expect_zero("d Sigma_nt / dlnvarpi", sp.diff(Sigma_nt, dlnvarpi))

    print("The orbit-clean conditions depend only on")
    print("  (dlnchi0, dlnZW, dlnepsilon_eta, dlnepsilon_W, dlndeltaU, dlnOmegaW),")
    print("and are blind to the explicit support pair")
    print("  (dlnC, dlnvarpi).")

    subbanner("IV. Exact orbit-layer formulas inside the overlap image")
    expect_zero(
        "dlndeltaU + (18/11) dlnchi0",
        sp.simplify(dlndeltaU + sp.Rational(18, 11) * dlnchi0),
    )
    expect_zero(
        "dlnOmegaW - (41/44 dlnchi0 + 5/8 dlnZW)",
        sp.simplify(dlnOmegaW - (sp.Rational(41, 44) * dlnchi0 + sp.Rational(5, 8) * dlnZW)),
    )
    expect_zero(
        "dlnepsilon_W - (2 dlnchi0 + dlnZW)",
        sp.simplify(dlnepsilon_W - (2 * dlnchi0 + dlnZW)),
    )

    banner("STAGE 225 LEDGER")
    print("1. The Stage-224 overlap image lies exactly tangent to the coherent monomial-preserving")
    print("   similarity orbit on the constructive branch: Sigma_eta = Sigma_tr = Sigma_nt = 0.")
    print("2. The orbit verdict is support-blind: the explicit support drifts (dlnC, dlnvarpi) do not")
    print("   enter the orbit-clean conditions at all.")
    print("3. The overlap image therefore splits naturally into")
    print("      a 2D orbit/shape layer  (dlnchi0, dlnZW),")
    print("      a 2D support layer      (dlnC, dlnvarpi).")
    print("4. The remaining realizability problem is not orbit lock. It is whether the actual moving-throat")
    print("   overlap tangent lands inside the full 4D overlap image from Stage 224.")
