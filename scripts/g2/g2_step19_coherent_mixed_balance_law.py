#!/usr/bin/env python3
"""
Step 19 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the exact coherent local-kernel formulas for M_mix and R_target.
2. Verifies the exact mixed-only product law
       R_target * M_mix = 8 * Lambda * (1 - eps) / pi^2.
3. Computes the logarithmic drift laws for
       delta ln R_target,
       delta ln M_mix,
       delta ln(R_target M_mix).
4. Imposes the Step-18 coherent branch selection law
       delta ln R_target = 0,
   together with the Step-16/14 slaved dressing law
       q_eta = -((1-eps_eta)/eps_eta) * Lambda_1.
5. Derives the exact coherent mixed/outgoing balance pair
       Lambda_1 = -q_Lambda + q_Z + 2 q_chi + sigma q_eps,
       Delta_mix = q_Lambda - (sigma/2) q_eps,
   where sigma = 2 eps / (1-eps).
6. Solves that pair exactly for (q_Lambda, q_eps) and computes the
   minimum-Euclidean-norm representative in the reduced variable set
       x = (q_Lambda, q_Z, q_chi, q_eps).

Interpretation
--------------
This is the first place where the quartic g-2 repair splits cleanly into two
separate microscopic jobs on the coherent branch:

  * one linear combination of mixed/outgoing drifts supplies the universal
    transfer-shape correction Lambda_1,
  * another linear combination sets the mixed-only baseline drift Delta_mix,
    with the remaining baseline compensation then carried by the support factor
    S(zeta;eps).

So the anomaly problem is no longer a vague branch-choice ambiguity. It is an
exact two-equation balance law in the mixed/outgoing microledger.
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
    banner("STEP 19 — COHERENT MIXED/OUTGOING BALANCE LAW")

    # Microscopic coherent-kernel variables.
    Lam, eps, eps_eta, ZW, chi0 = sp.symbols(
        "Lambda epsilon epsilon_eta Z_W chi_0", positive=True, real=True
    )
    pi = sp.pi

    M_mix = sp.simplify(8 * ZW * (1 + chi0) ** 2 / (pi**2 * (1 - eps_eta) * (1 - eps)))
    R_target = sp.simplify(Lam * (1 - eps_eta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))

    subbanner("XIX.1 — Exact coherent mixed-only product law")

    print("M_mix =")
    sp.pprint(M_mix)
    print("R_target =")
    sp.pprint(R_target)
    expect_zero(
        "R_target * M_mix - 8 Lambda (1-eps)/pi^2",
        sp.simplify(R_target * M_mix - 8 * Lam * (1 - eps) / pi**2),
    )

    # Logarithmic drift variables.
    qLam, qZ, qchi, qeta, qeps = sp.symbols(
        "q_Lambda q_Z q_chi q_eta q_eps", real=True
    )

    dln_R = sp.simplify(
        qLam - qZ - 2 * qchi - eps_eta / (1 - eps_eta) * qeta - 2 * eps / (1 - eps) * qeps
    )
    dln_M = sp.simplify(
        qZ + 2 * qchi + eps_eta / (1 - eps_eta) * qeta + eps / (1 - eps) * qeps
    )
    dln_prod = sp.simplify(dln_R + dln_M)

    subbanner("XIX.2 — Exact drift ledger for R_target, M_mix, and the mixed product")

    print("delta ln R_target =")
    sp.pprint(dln_R)
    print("delta ln M_mix =")
    sp.pprint(dln_M)
    expect_zero(
        "delta ln(R_target M_mix) - [q_Lambda - eps/(1-eps) q_eps]",
        sp.simplify(dln_prod - (qLam - eps / (1 - eps) * qeps)),
    )

    # Coherent branch closure and slaved dressing law.
    Lambda1, Delta_mix = sp.symbols("Lambda_1 Delta_mix", real=True)
    qeta_coh = sp.simplify(-(1 - eps_eta) * Lambda1 / eps_eta)
    sigma = sp.simplify(2 * eps / (1 - eps))

    subbanner("XIX.3 — Coherent fixed-target balance pair")

    # Coherent side: delta ln R_target = 0.
    dln_R_coh = sp.simplify(dln_R.subs(qeta, qeta_coh))
    expect_zero(
        "delta ln R_target|coh - [Lambda_1 + q_Lambda - q_Z - 2 q_chi - sigma q_eps]",
        sp.simplify(dln_R_coh - (Lambda1 + qLam - qZ - 2 * qchi - sigma * qeps)),
    )

    Lambda1_law = sp.simplify(-qLam + qZ + 2 * qchi + sigma * qeps)
    Delta_mix_raw = sp.simplify(dln_M.subs(qeta, qeta_coh))
    Delta_mix_law = sp.simplify(Delta_mix_raw.subs(Lambda1, Lambda1_law))

    print("Lambda_1 law =")
    sp.pprint(Lambda1_law)
    print("Raw coherent delta ln M_mix =")
    sp.pprint(Delta_mix_raw)
    expect_zero(
        "Delta_mix - [q_Lambda - (sigma/2) q_eps]",
        sp.simplify(Delta_mix_law - (qLam - sp.Rational(1, 2) * sigma * qeps)),
    )
    print("Delta_mix law =")
    sp.pprint(sp.simplify(qLam - sp.Rational(1, 2) * sigma * qeps))

    print("So the coherent branch splits the mixed/outgoing microledger into")
    print("  Lambda_1  = -q_Lambda + q_Z + 2 q_chi + sigma q_eps")
    print("  Delta_mix =  q_Lambda - (sigma/2) q_eps")
    print("with sigma = 2 eps/(1-eps).")

    subbanner("XIX.4 — Exact solution for q_Lambda and q_eps")

    sol = sp.solve(
        [sp.Eq(Lambda1, Lambda1_law), sp.Eq(Delta_mix, qLam - sp.Rational(1, 2) * sigma * qeps)],
        [qLam, qeps],
        dict=True,
    )[0]

    qLam_sol = sp.simplify(sol[qLam])
    qeps_sol = sp.simplify(sol[qeps])

    print("q_Lambda =")
    sp.pprint(qLam_sol)
    print("q_eps =")
    sp.pprint(qeps_sol)

    expect_zero(
        "Recovered Lambda_1 constraint",
        sp.simplify(Lambda1_law.subs({qLam: qLam_sol, qeps: qeps_sol}) - Lambda1),
    )
    expect_zero(
        "Recovered Delta_mix constraint",
        sp.simplify((qLam - sp.Rational(1, 2) * sigma * qeps).subs({qLam: qLam_sol, qeps: qeps_sol}) - Delta_mix),
    )

    subbanner("XIX.5 — Special coherent sub-branches")

    # (i) Support carries the baseline compensation: Delta_mix = 0.
    qLam_support = sp.simplify(qLam_sol.subs(Delta_mix, 0))
    qeps_support = sp.simplify(qeps_sol.subs(Delta_mix, 0))
    print("Support-carried baseline closure (Delta_mix = 0):")
    print("q_Lambda =")
    sp.pprint(qLam_support)
    print("q_eps =")
    sp.pprint(qeps_support)

    # (ii) Frozen overlap/interference geometry drifts qZ=qchi=0.
    qLam_frozen = sp.simplify(qLam_sol.subs({qZ: 0, qchi: 0}))
    qeps_frozen = sp.simplify(qeps_sol.subs({qZ: 0, qchi: 0}))
    print("Frozen overlap/interference branch (q_Z = q_chi = 0):")
    print("q_Lambda =")
    sp.pprint(qLam_frozen)
    print("q_eps =")
    sp.pprint(qeps_frozen)

    subbanner("XIX.6 — Minimum-norm reduced representative")

    x = sp.Matrix([qLam, qZ, qchi, qeps])
    A = sp.Matrix([
        [-1, 1, 2, sigma],
        [1, 0, 0, -sigma / 2],
    ])
    b = sp.Matrix([Lambda1, Delta_mix])
    x_min = sp.simplify(A.T * (A * A.T).inv() * b)

    print("x_min = (q_Lambda, q_Z, q_chi, q_eps)^T =")
    sp.pprint(x_min)
    expect_zero("A x_min - b", sp.simplify(A * x_min - b))

    x_min_support = sp.simplify(x_min.subs(Delta_mix, 0))
    print("Minimum-norm support-carried-baseline representative (Delta_mix = 0) =")
    sp.pprint(x_min_support)

    print()
    print("Interpretation:")
    print("  * Lambda_1 and Delta_mix are the two independent mixed/outgoing targets.")
    print("  * The coherent support lane then carries any remaining baseline correction")
    print("    through M_tr = M_mix S(zeta;eps).")
    print("  * The quartic repair is therefore a coupled but sharply organized balance law,")
    print("    not a broad branch-choice ambiguity.")


if __name__ == "__main__":
    main()
