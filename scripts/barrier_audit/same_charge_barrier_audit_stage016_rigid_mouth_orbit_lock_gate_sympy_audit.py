from __future__ import annotations

from pathlib import Path
from sympy import (
    Abs,
    Eq,
    Matrix,
    Rational,
    Symbol,
    symbols,
    simplify,
)


def main() -> str:
    out: list[str] = []

    # Stage-171 observable compiler symbols
    dln_Rtr, dln_Nstar, dln_epseta = symbols("dln_Rtr dln_Nstar dln_epseta")
    Bstar = Symbol("Bstar", finite=True)
    epseta = Symbol("epseta", positive=True)

    # Exact observable-to-defect compiler
    C_obs_to_def = Matrix(
        [
            [1, 0, 0],
            [-Bstar, 1, 0],
            [Bstar, -1, -epseta / (1 - epseta)],
        ]
    )
    Delta_obs = Matrix([dln_Rtr, dln_Nstar, dln_epseta])
    Theta1, Xi1, R1 = list(C_obs_to_def * Delta_obs)

    out.append("=== Stage 016 symbolic audit ===")
    out.append("Observable -> defect compiler:")
    out.append(str(C_obs_to_def))
    out.append("")

    out.append("Theta1:")
    out.append(str(simplify(Theta1)))
    out.append("Xi1:")
    out.append(str(simplify(Xi1)))
    out.append("R1:")
    out.append(str(simplify(R1)))
    out.append("")

    # Track-locked specialization
    Theta_track = simplify(Theta1.subs(dln_Rtr, 0))
    Xi_track = simplify(Xi1.subs(dln_Rtr, 0))
    R_track = simplify(R1.subs(dln_Rtr, 0))

    out.append("Track-locked specialization dln R_tr = 0:")
    out.append(f"Theta1 -> {Theta_track}")
    out.append(f"Xi1    -> {Xi_track}")
    out.append(f"R1     -> {R_track}")
    out.append("")

    # Static load defect identity
    N01, N0, D01, D0 = symbols("N01 N0 D01 D0", nonzero=True)
    Xi_load = simplify(N01 / N0 - D01 / D0)
    Xi_load_rigid = simplify(Xi_load.subs(D01, 0))

    out.append("Static load defect:")
    out.append(f"Xi_load = {Xi_load}")
    out.append(f"Under D01 = 0: Xi_load -> {Xi_load_rigid}")
    out.append("")

    # Transported ceiling formula
    Delta_norm, T_quad, mhat0_sq, Pcrit, eps = symbols(
        "Delta_norm T_quad mhat0_sq Pcrit eps", positive=True
    )
    ceiling = simplify(Pcrit * mhat0_sq / (Delta_norm + T_quad) - 1)
    out.append("Transported static ceiling:")
    out.append(
        "|eps * Xi1| <= " + str(ceiling)
    )
    out.append("Under track lock Xi1 = dln N_*:")
    out.append(
        "|eps * dln N_*| <= " + str(ceiling)
    )
    out.append("Under additional D01 = 0:")
    out.append(
        "|eps * N01/N0| <= " + str(ceiling)
    )
    out.append("")

    # Numerical Stage-007 budgets at the compatibility point
    Pbar = 0.002069792318062885
    Pcrit_both_10 = 0.0028313316855593175
    Pcrit_one_10 = 0.0035965105896846573
    budget_both_10 = Pcrit_both_10 / Pbar - 1.0
    budget_one_10 = Pcrit_one_10 / Pbar - 1.0

    out.append("Stage-007 compatibility-point budgets:")
    out.append(f"robust 10% both-pole budget   = {budget_both_10:.15f}")
    out.append(f"nonempty 10% corridor budget  = {budget_one_10:.15f}")
    out.append("")

    # Consistency checks
    assert simplify(Xi_track - dln_Nstar) == 0
    assert simplify(Xi_load_rigid - N01 / N0) == 0

    out.append("Checks passed:")
    out.append("- Xi1 collapses exactly to dln N_* under dln R_tr = 0.")
    out.append("- Xi_load collapses exactly to N01/N0 under D01 = 0.")
    out.append("- The transported ceiling compiles to a direct internal-transfer bound.")
    return "\n".join(out)


if __name__ == "__main__":
    text = main()
    print(text)
