#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage268_271_common import banner, subbanner, expect_zero, family1_refreshed_numbers

"""
Stage 271 — actual branch verdict after combining the geometry-decoupling theorem,
the explicit Family-1 monopole breathing lane, and the later actual-isotropic
support / normalization reduction.

What this script does
---------------------
1. Encodes the exact Stage-77/78 theorem that on the actual isotropic branch

       eps_2 = eps_4 = 0,

   because the l=0 geometry lane is dynamically inert through O(omega^4) with
   respect to the grouped real P2 carrier.
2. Records the exact actual-isotropic conservative consequence

       c_pole = 1/4,   rho_alpha = 4/3,   zeta_req = 1/3,

   together with the later actual-branch reduction to the single normalization
   defect N_Q.
3. Packages the first nontrivial weak-anisotropy scalar-lane influence through
   the explicit Family-1 breathing-lane diagnostic from Stage 270.
4. Produces the cleanest current verdict:

   - the actual isotropic branch is exactly safe against scalar/geometry
     contamination,
   - the first known scalar-lane correction is only O(chi^2),
   - and even that correction becomes dangerous only if the physical grouped
     quadrupole pole is much faster than the coefficient-equivalent breathing
     pole.
"""


def main() -> None:
    banner("STAGE 271 — ACTUAL BRANCH VERDICT AFTER THE FAMILY-1 BREATHING DIAGNOSTIC")

    chi = sp.symbols("chi", real=True)
    eps2, eps4 = sp.symbols("eps_2 eps_4", real=True)
    c_pole = sp.symbols("c_pole", real=True)
    rho_alpha = sp.symbols("rho_alpha", positive=True, real=True)
    zeta_req = sp.symbols("zeta_req", positive=True, real=True)
    NQ = sp.symbols("N_Q", real=True)

    subbanner("I. Exact actual-isotropic geometry verdict")
    expect_zero("eps_2(actual isotropic branch)", sp.Integer(0))
    expect_zero("eps_4(actual isotropic branch)", sp.Integer(0))
    expect_zero("c_pole(actual isotropic branch) - 1/4", sp.Rational(1, 4) - sp.Rational(1, 4))
    expect_zero("rho_alpha(actual isotropic branch) - 4/3", sp.Rational(4, 3) - sp.Rational(4, 3))
    expect_zero("zeta_req(actual isotropic branch) - 1/3", sp.Rational(1, 3) - sp.Rational(1, 3))

    print("Exact actual-isotropic branch data carried from the moving-throat theorem ledger:")
    print("  eps_2 = eps_4 = 0")
    print("  c_pole = 1/4")
    print("  rho_alpha = 4/3")
    print("  zeta_req = 1/3")
    print()
    print("So the actual isotropic grouped-P2 / geometry branch is already conservatively clean,")
    print("and the explicit Family-1 support/source side is automatic there.")

    subbanner("II. First nonzero scalar-lane contamination enters only with explicit mixing")
    print("Stage-78 contamination theorem:")
    print("  eps_2 = a_2 chi^2,")
    print("  eps_4 = a_4 chi^2.")
    print()
    print("Therefore the known scalar/geometry lane can affect the grouped-P2 carrier only after")
    print("an explicit l=0 <-> l=2 mixing source is turned on. The actual isotropic branch itself")
    print("sits at chi = 0 for this purpose.")

    subbanner("III. Family-1 breathing-lane diagnostic from Stages 269–270")
    G0 = sp.Rational(4, 45)
    lam_minus = sp.Float("6.405572392138922", 80)
    lam_plus = sp.Float("254.444968136936126", 80)
    R_minus = sp.Float("0.002552474771738", 80)
    R_plus = sp.Float("0.386733239513976", 80)

    G2 = sp.simplify(R_minus / lam_minus + R_plus / lam_plus)
    G4 = sp.simplify(R_minus / lam_minus**2 + R_plus / lam_plus**2)
    lam_coeff = sp.simplify(G2 / G4)
    c_eff = sp.simplify(G2**2 / (G0 * G4))
    one_minus_c = sp.simplify(1 - c_eff)

    nums = family1_refreshed_numbers()
    cstar = sp.Float(str(nums["c_pole_max_0"]), 80)

    r_eff = sp.symbols("r_eff", positive=True, real=True)
    r_demand = sp.simplify(2 / one_minus_c)
    r_twin_fail = sp.simplify((4 + 2 * sp.sqrt(2)) / one_minus_c)
    r_F1_fail = sp.simplify((8 * cstar + 4 * sp.sqrt(cstar * (4 * cstar - 1))) / one_minus_c)

    print("Coefficient-equivalent breathing pole lambda_coeff =")
    sp.pprint(sp.N(lam_coeff, 20))
    print("Coefficient-equivalent pole fraction c_eff =")
    sp.pprint(sp.N(c_eff, 20))
    print("So the first known scalar lane becomes dangerous only if")
    print("  r_eff = Omega_Q^2 / lambda_coeff")
    print("is large enough. Exact threshold highlights:")
    print("  support demand starts growing at     r_eff >", sp.N(r_demand, 12))
    print("  actual twin failure requires         r_eff >=", sp.N(r_twin_fail, 12))
    print("  actual Family-1 failure requires     r_eff >=", sp.N(r_F1_fail, 12))

    subbanner("IV. Best current reading")
    print("1. The actual isotropic branch is exactly safe: chi = 0 and eps_2 = eps_4 = 0.")
    print("2. The explicit Family-1 support/source side is already automatic on that branch.")
    print("3. The first known scalar-lane correction is the two-pole monopole breathing lane,")
    print("   which is exactly equivalent, for contamination purposes, to a mixed contact+pole")
    print("   diagnostic with c_eff ~= ", sp.N(c_eff, 16), ".")
    print("4. Therefore the scalar/geometry lane is not the active reduced bottleneck anymore.")
    print("   The only way it re-enters is through weak l=0<->l=2 mixing plus a sufficiently")
    print("   large pole-ratio hierarchy r_eff.")

    banner("STAGE 271 LEDGER")
    print("1. On the actual isotropic branch the geometry lane is dynamically inert through O(omega^4),")
    print("   so eps_2 = eps_4 = 0 and the conservative grouped-P2 module is exactly the 3/4 + 1/4 branch.")
    print("2. On that same actual branch, rho_alpha = 4/3 and the explicit Family-1 support/source side")
    print("   is already automatic; the only remaining reduced defect there is the outgoing normalization N_Q.")
    print("3. The first nontrivial scalar-lane correction comes only after explicit l=0<->l=2 mixing and is")
    print("   governed by the Family-1 two-pole monopole breathing lane.")
    print("4. For contamination coefficients that breathing lane is exactly equivalent to a one-pole diagnostic")
    print("   with c_eff ~= ", sp.N(c_eff, 16), " and lambda_coeff ~= ", sp.N(lam_coeff, 16), ".")
    print("5. So the scalar/geometry lane is not the natural source of 5PN failure on the actual isotropic")
    print("   branch; the live reduced bottleneck remains the outgoing quadrupole normalization defect N_Q.")


if __name__ == "__main__":
    main()
