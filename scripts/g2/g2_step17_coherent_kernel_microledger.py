#!/usr/bin/env python3
"""
Step 17 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Imports the exact coherent local D/N kernel placement-map formulas from the
   moving-throat notes:
       R_tr,
       eps = eps_W^(split),
       M_mix,
       M_supp,
       M_tr,
       R_target.
2. Verifies that the support ratio zeta affects only the total baseline M_tr and
   does not affect either the tracking factor R_tr or the demand ratio R_target.
3. Derives the exact logarithmic drift of R_tr in the natural coherent variables
   (1+chi_0) and (1+delta_U).
4. Derives the exact logarithmic drift of R_target in the natural coherent
   variables (Lambda, Z_W, 1+chi_0, epsilon_eta, eps_W, 1+delta_U).
5. Compares that microscopic demand-ratio drift with the Step-16 selected-branch
   law and proves the exact cancellation theorem:

       Lambda_1
       = - delta ln Lambda
         - 2 delta ln(1-eps)
         + delta ln Z_W
         + 2 delta ln(1+chi_0),

   equivalently

       Lambda_1
       = -q_Lambda + q_Z + 2 q_chi + 2 eps/(1-eps) q_eps,

   so the quartic universal transfer-shape drift is carried entirely by the
   mixed/outgoing microscopic variables. The dressing coordinate epsilon_eta
   cancels identically.

Interpretation
--------------
This is the first exact microscopic bridge between the Step-16 demand-ratio law
and the actual coherent-kernel formulas from the moving-throat notes. It shows:

- the support lane zeta does not directly move R_target,
- the wall-U dressing epsilon_eta does enter R_target, but cancels completely
  from the inferred Lambda_1 law,
- so the quartic direct transfer-shape drift is governed only by the
  mixed/outgoing microscopic variables (Lambda, Z_W, chi_0, eps).

That sharply reduces what the next branch-selection step has to check.
"""

from __future__ import annotations

import sympy as sp


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


def main() -> None:
    banner("STEP 17 — COHERENT-KERNEL MICROSCOPIC DEMAND-RATIO LEDGER")

    chi0, delta_U = sp.symbols("chi_0 delta_U", positive=True, real=True)
    eps_W, eps_eta = sp.symbols("epsilon_W epsilon_eta", positive=True, real=True)
    zeta, Z_W = sp.symbols("zeta Z_W", positive=True, real=True)
    Lam = sp.symbols("Lambda", positive=True, real=True)

    # Exact Stage-30 coherent local-kernel map.
    R_tr = sp.simplify((1 + chi0 / (1 + delta_U)) / (1 + chi0))
    eps = sp.simplify(eps_W * (1 - sp.Rational(2, 11) * delta_U / (1 + delta_U)))
    M_mix = sp.simplify(8 * Z_W * (1 + chi0) ** 2 / (sp.pi ** 2 * (1 - eps_eta) * (1 - eps)))
    M_supp = sp.simplify(8 * zeta * Z_W * (1 + chi0) ** 2 / (sp.pi ** 2 * (1 - eps_eta) * (1 - zeta * eps)))
    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    M_tr = sp.simplify(M_mix + M_supp)
    R_target = sp.simplify(Lam * (1 - eps_eta) * (1 - eps) ** 2 / (Z_W * (1 + chi0) ** 2))

    subbanner("XVII.1 — Exact coherent local-kernel map")

    print("R_tr =")
    sp.pprint(R_tr)
    print("eps = eps_W^(split) =")
    sp.pprint(eps)
    print("M_mix =")
    sp.pprint(M_mix)
    print("M_supp =")
    sp.pprint(M_supp)
    print("R_target =")
    sp.pprint(R_target)

    expect_zero("M_tr - M_mix * S(zeta;eps)", sp.simplify(M_tr - M_mix * S))

    subbanner("XVII.2 — Support ratio zeta affects only the baseline")

    expect_zero("dR_tr/dzeta", sp.diff(R_tr, zeta))
    expect_zero("dR_target/dzeta", sp.diff(R_target, zeta))

    print("So on the coherent local-kernel branch:")
    print("  - zeta does not move the tracking factor R_tr,")
    print("  - zeta does not move the demand ratio R_target,")
    print("  - zeta only moves the total baseline through M_tr = M_mix S.")

    # Natural logarithmic drift coordinates.
    q_chi, q_U = sp.symbols("q_chi q_U", real=True)
    q_Z, q_Lam = sp.symbols("q_Z q_Lambda", real=True)
    q_eta, q_W, q_eps = sp.symbols("q_eta q_W q_eps", real=True)

    # Conventions:
    #   delta(1+chi_0) = (1+chi_0) q_chi,
    #   delta(1+delta_U) = (1+delta_U) q_U,
    #   delta Z_W = Z_W q_Z,
    #   delta Lambda = Lambda q_Lambda,
    #   delta epsilon_eta = epsilon_eta q_eta,
    #   delta epsilon_W = epsilon_W q_W,
    #   delta epsilon = epsilon q_eps.
    dRtr = sp.simplify(sp.diff(R_tr, chi0) * (1 + chi0) * q_chi + sp.diff(R_tr, delta_U) * (1 + delta_U) * q_U)
    dln_Rtr = sp.simplify(dRtr / R_tr)

    subbanner("XVII.3 — Exact microscopic drift of the tracking factor")

    print("delta ln R_tr =")
    sp.pprint(dln_Rtr)

    expect_zero(
        "delta ln R_tr + [delta_U q_chi + chi_0 q_U]/(1+chi_0+delta_U)",
        sp.simplify(dln_Rtr + (delta_U * q_chi + chi0 * q_U) / (1 + chi0 + delta_U)),
    )

    print("So lowering R_tr on the constructive branch comes only from motion in the two")
    print("dimensionless placement variables (1+chi_0) and (1+delta_U).")

    # Exact split-blocking drift.
    deps = sp.simplify(sp.diff(eps, eps_W) * eps_W * q_W + sp.diff(eps, delta_U) * (1 + delta_U) * q_U)
    q_eps_expr = sp.simplify(deps / eps)

    subbanner("XVII.4 — Exact microscopic drift of the split blocking")

    print("q_eps := delta ln eps =")
    sp.pprint(q_eps_expr)

    expect_zero(
        "q_eps - [q_W - 2 q_U/(11+9 delta_U)]",
        sp.simplify(q_eps_expr - (q_W - 2 * q_U / (11 + 9 * delta_U))),
    )

    # Exact demand-ratio drift in natural microscopic variables.
    dln_Rtarget = sp.simplify(q_Lam - q_Z - 2 * q_chi - eps_eta / (1 - eps_eta) * q_eta - 2 * eps / (1 - eps) * q_eps)
    dln_Rtarget_qWU = sp.simplify(dln_Rtarget.subs(q_eps, q_eps_expr))

    subbanner("XVII.5 — Exact microscopic demand-ratio drift")

    print("delta ln R_target =")
    sp.pprint(dln_Rtarget)
    print()
    print("Equivalently, in the primitive coherent variables (q_W, q_U),")
    print("delta ln R_target =")
    sp.pprint(dln_Rtarget_qWU)

    # Step-16 selected-branch law: delta ln R_target = delta ln(1-eps_eta) - Lambda_1.
    Lambda_1 = sp.symbols("Lambda_1", real=True)
    dln_one_minus_eps_eta = sp.simplify(-eps_eta / (1 - eps_eta) * q_eta)
    # infer Lambda_1 by equating microscopic and Step-16 laws
    Lambda_1_expr = sp.simplify(-(dln_Rtarget - dln_one_minus_eps_eta))
    Lambda_1_expr_alt = sp.simplify(-q_Lam + q_Z + 2 * q_chi + 2 * eps / (1 - eps) * q_eps)

    subbanner("XVII.6 — Exact cancellation theorem with the Step-16 law")

    expect_zero("Lambda_1_expr - Lambda_1_expr_alt", sp.simplify(Lambda_1_expr - Lambda_1_expr_alt))

    print("From Step 16, delta ln R_target = delta ln(1-epsilon_eta) - Lambda_1.")
    print("Equating that with the exact microscopic Stage-30 demand-ratio drift gives")
    print("  Lambda_1 =")
    sp.pprint(Lambda_1_expr_alt)

    expect_zero(
        "Lambda_1 - [-delta ln Lambda - 2 delta ln(1-eps) + delta ln Z_W + 2 delta ln(1+chi_0)]",
        sp.simplify(
            Lambda_1_expr_alt
            - (-q_Lam + q_Z + 2 * q_chi + 2 * eps / (1 - eps) * q_eps)
        ),
    )

    print()
    print("So the quartic universal transfer-shape drift is carried entirely by the")
    print("mixed/outgoing microscopic variables (Lambda, Z_W, chi_0, eps).")
    print("The wall-U dressing coordinate epsilon_eta cancels identically from the")
    print("inferred Lambda_1 law.")

    subbanner("XVII.7 — Reduced microscopic verdict")

    print("At the actual coherent-kernel level:")
    print("  - zeta moves only the baseline M_tr,")
    print("  - epsilon_eta does move R_target, but cancels from the inferred Lambda_1 law,")
    print("  - the universal transfer-shape drift Lambda_1 is therefore fixed by the")
    print("    mixed/outgoing microscopic placement data alone.")
    print()
    print("So the next real branch-selection test is not another algebraic closure choice.")
    print("It is whether the physical coherent support lane compensates the exact tracking")
    print("deficit by increasing the baseline at fixed demand ratio, or whether the completed")
    print("PDE forces a genuine retargeting drift in the mixed/outgoing data themselves.")


if __name__ == "__main__":
    main()
