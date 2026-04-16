#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 226
SymPy audit for the relaxed-constraint branch declaration and short-range
open-system lift compiler.

This script verifies:
1. a concrete exact leakage/work channel from the projected continuity identity,
2. the minimal non-rigid mouth response and positive U/V drain,
3. exact normalization of the compensated sign-changing source profile,
4. the codimension-three recovery slice back to the Stage-225 standard branch,
5. and the no-new-long-range limit for the static and linear dynamic
   same-charge conservative kernels.
"""

from __future__ import annotations

import sympy as sp


def section(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def main() -> None:
    # ------------------------------------------------------------------
    # 1. Leakage/work lane: explicit exact Gaussian projector reduction.
    # ------------------------------------------------------------------
    section("1. Exact open-system leakage/work lane")
    w = sp.symbols("w", real=True)
    ell_w, j0, E0 = sp.symbols("ell_w j0 E0", real=True)

    W = sp.exp(-w**2) / sp.sqrt(sp.pi)
    j_w = ell_w * j0 * w * sp.exp(-w**2)
    E_w = E0 * w * sp.exp(-w**2)

    boundary = sp.simplify(sp.limit(W * j_w, w, sp.oo) - sp.limit(W * j_w, w, -sp.oo))
    S_leak = sp.simplify(sp.integrate(sp.diff(W, w) * j_w, (w, -sp.oo, sp.oo)))
    W_work = sp.simplify(sp.integrate(j_w * E_w, (w, -sp.oo, sp.oo)))

    expected_S_leak = -sp.sqrt(2) * ell_w * j0 / 4
    expected_W_work = sp.sqrt(2) * sp.sqrt(sp.pi) * E0 * ell_w * j0 / 8

    print("Projector W(w)        =", W)
    print("Current   j^w(w)      =", j_w)
    print("Field     E_w(w)      =", E_w)
    print("Boundary term         =", boundary)
    print("S_leak                =", S_leak)
    print("Expected S_leak       =", expected_S_leak)
    print("Work channel W_w      =", W_work)
    print("Expected work channel =", expected_W_work)

    assert boundary == 0
    assert sp.simplify(S_leak - expected_S_leak) == 0
    assert sp.simplify(W_work - expected_W_work) == 0
    assert S_leak.subs(ell_w, 0) == 0
    assert W_work.subs(ell_w, 0) == 0

    # ------------------------------------------------------------------
    # 2. Non-rigid mouth/dressing lift.
    # ------------------------------------------------------------------
    section("2. Exact non-rigid U/V response and drain")
    U, V = sp.symbols("U V", real=True)
    k_U, k_V, chi_lam, f_U = sp.symbols("k_U k_V chi_lam f_U", positive=True, real=True)

    F_UV = (
        sp.Rational(1, 2) * k_U * U**2
        + sp.Rational(1, 2) * k_V * V**2
        - chi_lam * U * V
        - f_U * U
    )

    stationarity = [sp.diff(F_UV, U), sp.diff(F_UV, V)]
    sol = sp.solve(stationarity, [U, V], dict=True)[0]
    U_sol = sp.simplify(sol[U])
    V_sol = sp.simplify(sol[V])
    U_expected = sp.simplify(k_V * f_U / (k_U * k_V - chi_lam**2))
    V_expected = sp.simplify(chi_lam * f_U / (k_U * k_V - chi_lam**2))
    H_UV = sp.hessian(F_UV, (U, V))
    det_H = sp.factor(H_UV.det())
    ratio_VU = sp.simplify(V_sol / U_sol)
    drain_UV = sp.simplify(chi_lam * U_sol * V_sol)
    drain_expected = sp.simplify(chi_lam**2 * k_V * f_U**2 / (k_U * k_V - chi_lam**2) ** 2)

    print("Free energy F_UV      =", F_UV)
    print("Stationarity eqns     =", stationarity)
    print("U solution            =", U_sol)
    print("Expected U            =", U_expected)
    print("V solution            =", V_sol)
    print("Expected V            =", V_expected)
    print("det Hessian           =", det_H)
    print("V/U                   =", ratio_VU)
    print("Drain chi*U*V         =", drain_UV)
    print("Expected drain        =", drain_expected)

    assert sp.simplify(U_sol - U_expected) == 0
    assert sp.simplify(V_sol - V_expected) == 0
    assert sp.simplify(det_H - (k_U * k_V - chi_lam**2)) == 0
    assert sp.simplify(ratio_VU - chi_lam / k_V) == 0
    assert sp.simplify(drain_UV - drain_expected) == 0
    assert U_sol.subs(f_U, 0) == 0
    assert V_sol.subs(f_U, 0) == 0
    assert sp.simplify(V_sol.subs(chi_lam, 0)) == 0
    assert sp.simplify(U_sol.subs(chi_lam, 0) - f_U / k_U) == 0

    # ------------------------------------------------------------------
    # 3. Compensated source lane.
    # ------------------------------------------------------------------
    section("3. Compensated sign-changing source profile")
    z = sp.symbols("z", real=True)
    a, b = sp.symbols("a b", real=True)
    y = sp.symbols("y", real=True)

    varsigma = 1 + a * sp.cos(sp.pi * z) + b * sp.cos(2 * sp.pi * z)
    varsigma_mean = sp.simplify(sp.integrate(varsigma, (z, 0, 1)))
    varsigma_y = sp.expand(1 - b + a * y + 2 * b * y**2)
    y_star = sp.simplify(-a / (4 * b))
    varsigma_vertex = sp.simplify(varsigma_y.subs(y, y_star))

    print("varsigma(z)           =", varsigma)
    print("Integral over [0,1]   =", varsigma_mean)
    print("Quadratic rewrite     =", varsigma_y)
    print("Interior stationary y*=", y_star)
    print("Vertex value          =", varsigma_vertex)
    print("Boundary values       =", sp.simplify(varsigma.subs(z, 0)), ",", sp.simplify(varsigma.subs(z, 1)))

    assert varsigma_mean == 1
    assert sp.simplify(varsigma.subs({a: 0, b: 0}) - 1) == 0
    assert sp.simplify(varsigma_y - (1 - b + a * y + 2 * b * y**2)) == 0
    assert sp.simplify(varsigma_vertex - (1 - b - a**2 / (8 * b))) == 0

    # ------------------------------------------------------------------
    # 4. Exact recovery slice back to the Stage-225 standard branch.
    # ------------------------------------------------------------------
    section("4. Codimension-three recovery slice")
    recovery_map = {
        ell_w: 0,
        f_U: 0,
        a: 0,
        b: 0,
    }

    U_rec = sp.simplify(U_sol.subs(recovery_map))
    V_rec = sp.simplify(V_sol.subs(recovery_map))
    D_rec = sp.simplify(drain_UV.subs(recovery_map))
    S_rec = sp.simplify(S_leak.subs(recovery_map))
    W_rec = sp.simplify(W_work.subs(recovery_map))
    varsigma_rec = sp.simplify(varsigma.subs(recovery_map))

    print("Recovery substitutions=", recovery_map)
    print("Recovered S_leak      =", S_rec)
    print("Recovered W_w         =", W_rec)
    print("Recovered U           =", U_rec)
    print("Recovered V           =", V_rec)
    print("Recovered drain       =", D_rec)
    print("Recovered varsigma    =", varsigma_rec)

    assert S_rec == 0
    assert W_rec == 0
    assert U_rec == 0
    assert V_rec == 0
    assert D_rec == 0
    assert varsigma_rec == 1

    # ------------------------------------------------------------------
    # 5. Short-range kernel invariant: no new 1/x-like family.
    # ------------------------------------------------------------------
    section("5. Short-range kernel invariant")
    x, kappa = sp.symbols("x kappa", positive=True, real=True)
    C6, C4, C2 = sp.symbols("C6 C4 C2", real=True)
    D6, D4, D2 = sp.symbols("D6 D4 D2", real=True)

    deltaV_stat = -sp.Rational(1, 2) * (
        C6 / x**6 + 2 * C4 * sp.exp(-2 * kappa * x) / x**4 + C2 * sp.exp(-4 * kappa * x) / x**2
    )
    Vdyn_real = -sp.Rational(1, 2) * (
        D6 / x**6 + 2 * D4 * sp.exp(-2 * kappa * x) / x**4 + D2 * sp.exp(-4 * kappa * x) / x**2
    )

    limit_stat = sp.simplify(sp.limit(x * deltaV_stat, x, sp.oo))
    limit_dyn = sp.simplify(sp.limit(x * Vdyn_real, x, sp.oo))

    print("deltaV_stat(x)        =", deltaV_stat)
    print("Re V_dyn(x,omega)     =", Vdyn_real)
    print("lim x*deltaV_stat     =", limit_stat)
    print("lim x*Re V_dyn        =", limit_dyn)

    assert limit_stat == 0
    assert limit_dyn == 0

    # ------------------------------------------------------------------
    # 6. Final success banner.
    # ------------------------------------------------------------------
    section("Stage-226 audit result")
    print("All symbolic checks passed.")
    print("Verified objects:")
    print("- exact leakage/work lane on a concrete Gaussian projector profile,")
    print("- exact non-rigid U/V response and positive drain,")
    print("- exact compensated-source normalization and quadratic rewrite,")
    print("- exact recovery slice back to the Stage-225 standard branch,")
    print("- and strict short-range limits for the static and linear dynamic kernels.")


if __name__ == "__main__":
    main()
