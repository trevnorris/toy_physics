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
Stage 220 — full constructive support-restored master solve.

What this script does
---------------------
1. Freezes the constructive orbit point (E_*,F_*) = (1/4,5/6).
2. Solves the full support-restored 3x7 packet master system
      u2^(1) = 0, u4^(1) = 0, Xi_load = 0
   for the wall/support directions
      (alpha_K, alpha_GU, alpha_OU)
   in terms of the four free carriers
      (alpha_GW, alpha_R, beta_B, beta_varpi).
3. Exhibits the corresponding 4-dimensional packet-null basis.
4. Shows that the support-only rescue from Stage 219 is a genuine subfamily of the
   full master solve.

Interpretation
--------------
The current isotropic baseline now carries a full 4-parameter first-order Packet-A/B
null family. So there is no reason to pivot baselines yet; the correct next move is
still to understand which of these support-restored corridors can be realized by the
actual moving-throat branch.
"""


if __name__ == "__main__":
    banner("STAGE 220 — FULL CONSTRUCTIVE SUPPORT-RESTORED MASTER SOLVE")

    E_star = sp.Rational(1, 4)
    F_star = sp.Rational(5, 6)
    data = build_master_system(E_star, F_star)
    free = data["free"]
    alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi = free
    u2_1 = data["u2_1"]
    u4_1 = data["u4_1"]
    Xi_load = data["Xi_load"]

    subbanner("I. Exact full constructive solve")
    eqs = [sp.Eq(u2_1, 0), sp.Eq(u4_1, 0), sp.Eq(Xi_load, 0)]
    sol_list = sp.solve(eqs, [alpha_K, alpha_GU, alpha_OU], dict=True)
    if len(sol_list) != 1:
        raise AssertionError("Expected a unique full constructive solve for (alpha_K, alpha_GU, alpha_OmegaU).")
    sol = {k: sp.simplify(v) for k, v in sol_list[0].items()}

    print("alpha_K =")
    sp.pprint(sol[alpha_K])
    print("alpha_GU =")
    sp.pprint(sol[alpha_GU])
    print("alpha_OmegaU =")
    sp.pprint(sol[alpha_OU])

    expect_zero("u2^(1) after full solve", sp.simplify(u2_1.subs(sol)))
    expect_zero("u4^(1) after full solve", sp.simplify(u4_1.subs(sol)))
    expect_zero("Xi_load after full solve", sp.simplify(Xi_load.subs(sol)))

    subbanner("II. Canonical null basis in free carriers (alpha_GW, alpha_R, beta_B, beta_varpi)")
    expr_vec = sp.Matrix([sol[alpha_K], alpha_GW, sol[alpha_GU], alpha_R, sol[alpha_OU], beta_B, beta_varpi])

    basis = {
        "alpha_GW": sp.simplify(expr_vec.subs({alpha_GW: 1, alpha_R: 0, beta_B: 0, beta_varpi: 0})),
        "alpha_R": sp.simplify(expr_vec.subs({alpha_GW: 0, alpha_R: 1, beta_B: 0, beta_varpi: 0})),
        "beta_B": sp.simplify(expr_vec.subs({alpha_GW: 0, alpha_R: 0, beta_B: 1, beta_varpi: 0})),
        "beta_varpi": sp.simplify(expr_vec.subs({alpha_GW: 0, alpha_R: 0, beta_B: 0, beta_varpi: 1})),
    }

    for name, vec in basis.items():
        print(f"Basis vector for {name} =")
        sp.pprint(vec)
        print("Numerically =")
        sp.pprint(sp.Matrix([sp.N(x) for x in vec]))
        subs = {
            alpha_K: vec[0],
            alpha_GW: vec[1],
            alpha_GU: vec[2],
            alpha_R: vec[3],
            alpha_OU: vec[4],
            beta_B: vec[5],
            beta_varpi: vec[6],
        }
        expect_zero(f"u2^(1) on {name} basis", sp.simplify(u2_1.subs(subs)))
        expect_zero(f"u4^(1) on {name} basis", sp.simplify(u4_1.subs(subs)))
        expect_zero(f"Xi_load on {name} basis", sp.simplify(Xi_load.subs(subs)))

    subbanner("III. Support-only rescue is a subfamily")
    support_only_subs = {alpha_GW: 0, alpha_R: 0}
    print("Full solve restricted to alpha_GW = alpha_R = 0 =")
    sp.pprint({k: sp.simplify(v.subs(support_only_subs)) for k, v in sol.items()})

    banner("STAGE 220 LEDGER")
    print("1. On the constructive orbit point, the full support-restored packet master system")
    print("   has a 4-parameter exact null family.")
    print("2. The free carriers can be taken as")
    print("      (alpha_GW, alpha_R, beta_B, beta_varpi),")
    print("   with the wall/support drifts")
    print("      (alpha_K, alpha_GU, alpha_OmegaU)")
    print("   fixed linearly by the true packet-null conditions.")
    print("3. The Stage-219 support-only rescue is recovered by setting alpha_GW = alpha_R = 0.")
    print("4. Therefore the current isotropic baseline remains live after support restoration.")
    print("5. The next theorem gate is not a baseline pivot. It is deciding which of these")
    print("   support-restored null directions can be realized by the actual moving-throat branch.")
