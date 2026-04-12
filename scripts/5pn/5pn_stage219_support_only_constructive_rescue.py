#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero
import importlib.util
from pathlib import Path

_stage218_path = Path(__file__).with_name('5pn_stage218_support_restored_master_matrix.py')
_spec = importlib.util.spec_from_file_location('stage218_support_restored_master_matrix', _stage218_path)
_mod = importlib.util.module_from_spec(_spec)
assert _spec is not None and _spec.loader is not None
_spec.loader.exec_module(_mod)
build_master_system = _mod.build_master_system

"""
Stage 219 — constructive support-only packet-null rescue.

What this script does
---------------------
1. Freezes the constructive orbit point (E_*,F_*) = (1/4,5/6) on the same explicit
   isotropic prototype used in Stages 215–218.
2. Restricts to the support-only slice
      alpha_GW = 0, alpha_R = 0,
   while keeping the support carriers
      (alpha_K, alpha_GU, alpha_OU, beta_B, beta_varpi).
3. Solves the *true* first-order Packet-A/B tangency equations
      u2^(1) = 0, u4^(1) = 0, Xi_load = 0
   exactly for (alpha_K, alpha_GU, alpha_OU) in terms of (beta_B, beta_varpi).
4. Exhibits the resulting two independent support-only packet-null directions.

Interpretation
--------------
The current isotropic baseline already survives once the BdG support variables are
restored. The support-blind failure was real, but it was only a slice failure.
"""


if __name__ == "__main__":
    banner("STAGE 219 — CONSTRUCTIVE SUPPORT-ONLY PACKET-NULL RESCUE")

    E_star = sp.Rational(1, 4)
    F_star = sp.Rational(5, 6)
    data = build_master_system(E_star, F_star)
    free = data["free"]
    alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi = free
    u2_1 = data["u2_1"]
    u4_1 = data["u4_1"]
    Xi_load = data["Xi_load"]

    subbanner("I. Exact support-only solve")
    eqs = [
        sp.Eq(sp.simplify(u2_1.subs({alpha_GW: 0, alpha_R: 0})), 0),
        sp.Eq(sp.simplify(u4_1.subs({alpha_GW: 0, alpha_R: 0})), 0),
        sp.Eq(sp.simplify(Xi_load.subs({alpha_GW: 0, alpha_R: 0})), 0),
    ]
    sol_list = sp.solve(eqs, [alpha_K, alpha_GU, alpha_OU], dict=True)
    if len(sol_list) != 1:
        raise AssertionError("Expected a unique support-only solve for (alpha_K, alpha_GU, alpha_OmegaU).")
    sol = {k: sp.simplify(v) for k, v in sol_list[0].items()}

    print("alpha_K =")
    sp.pprint(sol[alpha_K])
    print("alpha_GU =")
    sp.pprint(sol[alpha_GU])
    print("alpha_OmegaU =")
    sp.pprint(sol[alpha_OU])

    expect_zero("u2^(1) after support-only solve", sp.simplify(u2_1.subs({alpha_GW: 0, alpha_R: 0, **sol})))
    expect_zero("u4^(1) after support-only solve", sp.simplify(u4_1.subs({alpha_GW: 0, alpha_R: 0, **sol})))
    expect_zero("Xi_load after support-only solve", sp.simplify(Xi_load.subs({alpha_GW: 0, alpha_R: 0, **sol})))

    subbanner("II. Canonical support-only packet-null basis")
    expr_vec = sp.Matrix([sol[alpha_K], sp.Integer(0), sol[alpha_GU], sp.Integer(0), sol[alpha_OU], beta_B, beta_varpi])

    basis_betaB = sp.simplify(expr_vec.subs({beta_B: 1, beta_varpi: 0}))
    basis_betav = sp.simplify(expr_vec.subs({beta_B: 0, beta_varpi: 1}))

    print("Basis vector with beta_B = 1, beta_varpi = 0 =")
    sp.pprint(basis_betaB)
    print("Numerically =")
    sp.pprint(sp.Matrix([sp.N(x) for x in basis_betaB]))

    print("Basis vector with beta_B = 0, beta_varpi = 1 =")
    sp.pprint(basis_betav)
    print("Numerically =")
    sp.pprint(sp.Matrix([sp.N(x) for x in basis_betav]))

    # Direct verification on both basis vectors.
    for name, vec in [("beta_B basis", basis_betaB), ("beta_varpi basis", basis_betav)]:
        subs = {
            alpha_K: vec[0],
            alpha_GW: vec[1],
            alpha_GU: vec[2],
            alpha_R: vec[3],
            alpha_OU: vec[4],
            beta_B: vec[5],
            beta_varpi: vec[6],
        }
        expect_zero(f"u2^(1) on {name}", sp.simplify(u2_1.subs(subs)))
        expect_zero(f"u4^(1) on {name}", sp.simplify(u4_1.subs(subs)))
        expect_zero(f"Xi_load on {name}", sp.simplify(Xi_load.subs(subs)))

    banner("STAGE 219 LEDGER")
    print("1. On the constructive orbit point (E_*,F_*) = (1/4,5/6), the support-only slice")
    print("   already carries a 2-parameter exact packet-null family.")
    print("2. The exact solve expresses")
    print("      (alpha_K, alpha_GU, alpha_OmegaU)")
    print("   linearly in")
    print("      (beta_B, beta_varpi).")
    print("3. So the same isotropic prototype that killed the support-blind corridor is rescued")
    print("   immediately once the BdG support carriers are restored.")
    print("4. This means the current continuation should stay on the present isotropic baseline")
    print("   and solve the support-restored master matrix before considering any pivot.")
