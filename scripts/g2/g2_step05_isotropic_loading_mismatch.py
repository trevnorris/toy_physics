#!/usr/bin/env python3
"""
Step 5 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the weak-axisymmetric grouped-lane expansion on the isotropic
   one-port branch and derives the exact first-order physical slopes
       u_2^(1), u_4^(1), P_1, P_1/P_0.
2. Uses the canonical compensated/even-preserving branch to prove that the
   conservative grouped response can be kept fixed to first order iff
       D_21 = -D_01/9,
       D_41 = -D_01/27,
   so the only remaining linear outlet defect is the loading scalar
       Xi_load = P_1/P_0 = N_01/N_0 - D_01/D_0.
3. Matches that scalar to the anomaly quartic target Lambda_1 and solves the
   minimal normalized loading families:
       n - d = Lambda_1,
   where n := N_01/N_0 and d := D_01/D_0.
4. Resolves the conservative static slope further into the microscopic reduced
   bundle lanes
       d = k - b - z,
   with k := K_01/D_0, b := B_01/D_0, z := Z_01/D_0,
   and computes the minimum-norm conservative split.
5. Uses the one-port transfer formula
       N_1/N_0 = 2 pi_1 - 2 delta_1
   (equivalently 2 d ln P_r - 2 d ln Delta_r)
   to obtain the minimum-norm port-level deformation that realizes a chosen
   outgoing-transfer slope n.

Interpretation
--------------
This step is the first point where the quartic anomaly correction simplifies to
one scalar static loading mismatch between outgoing-transfer strengthening and
conservative static softening. That is a genuine simplification of the old
staggered derivation.
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


def step5_isotropic_loading_mismatch() -> None:
    banner("STEP 5 — ISOTROPIC ONE-PORT LOADING MISMATCH")

    # ------------------------------------------------------------------
    # I. Exact physical slopes on the weak-axisymmetric grouped branch
    # ------------------------------------------------------------------
    subbanner("V.1 — Exact weak-axisymmetric physical slopes")

    eps, lam = sp.symbols("epsilon lambda_A", real=True)
    D0, D2, D4, N0 = sp.symbols("D_0 D_2 D_4 N_0", real=True, nonzero=True)
    D01, D21, D41, N01 = sp.symbols("D_01 D_21 D_41 N_01", real=True)

    D_A0 = D0 + eps * lam * D01
    D_A2 = D2 + eps * lam * D21
    D_A4 = D4 + eps * lam * D41
    N_A0 = N0 + eps * lam * N01

    u2_A = sp.simplify(-D_A2 / D_A0)
    u4_A = sp.simplify((D_A2**2 - D_A0 * D_A4) / D_A0**2)
    P_A0 = sp.simplify(N_A0 / D_A0)

    u2_1 = sp.simplify(sp.diff(u2_A, eps).subs(eps, 0) / lam)
    u4_1 = sp.simplify(sp.diff(u4_A, eps).subs(eps, 0) / lam)
    P1 = sp.simplify(sp.diff(P_A0, eps).subs(eps, 0) / lam)
    P1_over_P0 = sp.simplify(P1 / (N0 / D0))

    print("u_2^(1) =")
    sp.pprint(u2_1)
    print("u_4^(1) =")
    sp.pprint(u4_1)
    print("P_1 =")
    sp.pprint(P1)
    print("P_1 / P_0 =")
    sp.pprint(P1_over_P0)

    expect_zero("u_2^(1) + (D_21 + u_2 D_01)/D_0",
                u2_1 + (D21 + (-D2 / D0) * D01) / D0)
    expect_zero("P_1/P_0 - (N_01/N_0 - D_01/D_0)",
                P1_over_P0 - (N01 / N0 - D01 / D0))

    # ------------------------------------------------------------------
    # II. Canonical compensated/even-preserving branch
    # ------------------------------------------------------------------
    subbanner("V.2 — Canonical compensated/even-preserving branch")

    u2_can = sp.Rational(1, 9)
    u4_can = sp.Rational(4, 81)
    D2_can = sp.simplify(-u2_can * D0)
    D4_can = sp.simplify((u2_can**2 - u4_can) * D0)

    print("Canonical compensated even data:")
    print(f"  u_2 = {u2_can}")
    print(f"  u_4 = {u4_can}")
    print("  D_2 =")
    sp.pprint(D2_can)
    print("  D_4 =")
    sp.pprint(D4_can)

    u2_1_can = sp.simplify(u2_1.subs({D2: D2_can}))
    u4_1_can = sp.simplify(u4_1.subs({D2: D2_can, D4: D4_can}))

    print("u_2^(1) on canonical branch =")
    sp.pprint(u2_1_can)
    print("u_4^(1) on canonical branch =")
    sp.pprint(u4_1_can)

    even_preserving = sp.solve(
        [sp.Eq(u2_1_can, 0), sp.Eq(u4_1_can, 0)],
        [D21, D41],
        dict=True,
    )[0]

    print("Even-preserving conditions (u_2^(1) = u_4^(1) = 0):")
    for k, v in even_preserving.items():
        print(f"  {k} = {sp.simplify(v)}")

    expect_zero("D_21 + D_01/9", even_preserving[D21] + D01 / 9)
    expect_zero("D_41 + D_01/27", even_preserving[D41] + D01 / 27)

    # Under that branch the only remaining linear outlet defect is P1/P0.
    Xi_load = sp.simplify(P1_over_P0)
    print("Xi_load := P_1/P_0 =")
    sp.pprint(Xi_load)

    # ------------------------------------------------------------------
    # III. Quartic anomaly target and minimal normalized loading families
    # ------------------------------------------------------------------
    subbanner("V.3 — Quartic target as a loading mismatch")

    Lam1 = sp.symbols("Lambda_1", real=True)
    n, d, mu = sp.symbols("n d mu", real=True)

    print("Normalized loading slopes:")
    print("  n := N_01 / N_0")
    print("  d := D_01 / D_0")
    print("Constraint:  n - d = Lambda_1")

    pure_transfer = {n: Lam1, d: 0}
    pure_softening = {n: 0, d: -Lam1}

    min_nd = sp.solve(
        [sp.Eq(2 * n + mu, 0), sp.Eq(2 * d - mu, 0), sp.Eq(n - d, Lam1)],
        [n, d, mu],
        dict=True,
    )[0]

    print("Pure transfer family:")
    sp.pprint(pure_transfer)
    print("Pure conservative-softening family:")
    sp.pprint(pure_softening)
    print("Minimum-norm (n,d) family:")
    sp.pprint({n: sp.simplify(min_nd[n]), d: sp.simplify(min_nd[d])})

    expect_zero("n_min - Lambda_1/2", min_nd[n] - Lam1 / 2)
    expect_zero("d_min + Lambda_1/2", min_nd[d] + Lam1 / 2)

    # ------------------------------------------------------------------
    # IV. Conservative static denominator split d = k - b - z
    # ------------------------------------------------------------------
    subbanner("V.4 — Minimal microscopic split of the conservative static slope")

    k, b, z, nu = sp.symbols("k b z nu", real=True)
    print("Microscopic conservative reduced-bundle slopes:")
    print("  k := K_01 / D_0")
    print("  b := B_01 / D_0")
    print("  z := Z_01 / D_0")
    print("Constraint:  d = k - b - z")

    conservative_min = sp.solve(
        [
            sp.Eq(2 * k + nu, 0),
            sp.Eq(2 * b - nu, 0),
            sp.Eq(2 * z - nu, 0),
            sp.Eq(k - b - z, d),
        ],
        [k, b, z, nu],
        dict=True,
    )[0]

    print("Minimum-norm conservative split:")
    sp.pprint({
        k: sp.simplify(conservative_min[k]),
        b: sp.simplify(conservative_min[b]),
        z: sp.simplify(conservative_min[z]),
    })

    expect_zero("k_min - d/3", conservative_min[k] - d / 3)
    expect_zero("b_min + d/3", conservative_min[b] + d / 3)
    expect_zero("z_min + d/3", conservative_min[z] + d / 3)

    balanced_conservative = {
        k: sp.simplify(conservative_min[k].subs(d, min_nd[d])),
        b: sp.simplify(conservative_min[b].subs(d, min_nd[d])),
        z: sp.simplify(conservative_min[z].subs(d, min_nd[d])),
    }

    print("Balanced-branch conservative minimum-norm realization:")
    sp.pprint(balanced_conservative)

    # ------------------------------------------------------------------
    # V. One-port outgoing-transfer slope n = 2*pi_1 - 2*delta_1
    # ------------------------------------------------------------------
    subbanner("V.5 — Minimal one-port outgoing-transfer deformation")

    pi1, del1, xi = sp.symbols("pi_1 delta_1 xi", real=True)
    print("On a single outgoing-transfer port:")
    print("  n = 2*pi_1 - 2*delta_1")
    print("where pi_1 is the port-amplitude log slope and delta_1 the")
    print("static-denominator log slope of the port block.")

    port_min = sp.solve(
        [
            sp.Eq(2 * pi1 + xi, 0),
            sp.Eq(2 * del1 - xi, 0),
            sp.Eq(2 * pi1 - 2 * del1, n),
        ],
        [pi1, del1, xi],
        dict=True,
    )[0]

    print("Minimum-norm one-port deformation at fixed n:")
    sp.pprint({pi1: sp.simplify(port_min[pi1]), del1: sp.simplify(port_min[del1])})

    expect_zero("pi_1,min - n/4", port_min[pi1] - n / 4)
    expect_zero("delta_1,min + n/4", port_min[del1] + n / 4)

    balanced_port = {
        pi1: sp.simplify(port_min[pi1].subs(n, min_nd[n])),
        del1: sp.simplify(port_min[del1].subs(n, min_nd[n])),
    }

    print("Balanced-branch minimum-norm one-port deformation:")
    sp.pprint(balanced_port)

    # ------------------------------------------------------------------
    # VI. Numerical benchmark values for the current anomaly target
    # ------------------------------------------------------------------
    subbanner("V.6 — Numerical benchmark values")

    Lam1_num = sp.Float("0.279605891931464", 30)

    n_bal = sp.N(min_nd[n].subs(Lam1, Lam1_num), 18)
    d_bal = sp.N(min_nd[d].subs(Lam1, Lam1_num), 18)
    k_bal = sp.N(balanced_conservative[k].subs(Lam1, Lam1_num), 18)
    b_bal = sp.N(balanced_conservative[b].subs(Lam1, Lam1_num), 18)
    z_bal = sp.N(balanced_conservative[z].subs(Lam1, Lam1_num), 18)
    pi_bal = sp.N(balanced_port[pi1].subs(Lam1, Lam1_num), 18)
    del_bal = sp.N(balanced_port[del1].subs(Lam1, Lam1_num), 18)

    print(f"Lambda_1 = {Lam1_num}")
    print("Balanced minimum-norm normalized loading split:")
    print(f"  n = N_01/N_0   = {n_bal}")
    print(f"  d = D_01/D_0   = {d_bal}")
    print("Balanced minimum-norm conservative split:")
    print(f"  k = K_01/D_0   = {k_bal}")
    print(f"  b = B_01/D_0   = {b_bal}")
    print(f"  z = Z_01/D_0   = {z_bal}")
    print("Balanced minimum-norm one-port outgoing deformation:")
    print(f"  pi_1           = {pi_bal}")
    print(f"  delta_1        = {del_bal}")

    # ------------------------------------------------------------------
    # VII. Final interpretation
    # ------------------------------------------------------------------
    subbanner("V.7 — Interpretation")
    print("1. On the simplest isotropic one-port branch, the anomaly quartic scalar")
    print("   from Step 4 is exactly the weak-axisymmetric loading mismatch")
    print("      Xi_load = P_1/P_0 = N_01/N_0 - D_01/D_0.")
    print("2. On the canonical even-preserving branch, the conservative grouped")
    print("   response is transported by a single static slope D_01; the even")
    print("   response coefficients themselves stay fixed to first order.")
    print("3. The remaining quartic g-2 problem therefore collapses to a static")
    print("   competition between outgoing-transfer strengthening (n) and")
    print("   conservative static softening (d).")
    print("4. The minimum-norm balanced realization splits the needed loading")
    print("   equally between those two normalized channels.")


if __name__ == "__main__":
    step5_isotropic_loading_mismatch()
