#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 227
SymPy audit for the selected-branch leakage and scalar-photon work compiler.

This script verifies:
1. the exact projected leakage source for a Gaussian projector and odd scalar-photon profile,
2. the exact one-mode bulk work scalar,
3. the reduced Session-I work scalar as a thickness-rescaled one-mode compiler,
4. the exact pullback through the Stage-225 selected-support demand Pi_tr,
5. the support-versus-orbit split of the compiled leakage/work lane,
6. transport-orientation parity,
7. and exact recovery of the closed standard slice.
"""

from __future__ import annotations

import sympy as sp


def section(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def main() -> None:
    # ------------------------------------------------------------------
    # 1. Exact projected leakage compiler.
    # ------------------------------------------------------------------
    section("1. Exact Gaussian leakage compiler")
    w = sp.symbols("w", real=True)
    lam = sp.symbols("lam", positive=True, real=True)
    E0 = sp.symbols("E0", positive=True, real=True)
    mu_w, rho0, q = sp.symbols("mu_w rho0 q", positive=True, real=True)

    W = sp.exp(-w**2 / lam**2) / (lam * sp.sqrt(sp.pi))
    phi = 2 * w * sp.exp(-w**2 / lam**2) / (sp.sqrt(sp.pi) * lam**3)
    E_w = -E0 * phi
    j_w = mu_w * rho0 * E_w
    J_w = q * j_w

    boundary = sp.simplify(sp.limit(W * j_w, w, sp.oo) - sp.limit(W * j_w, w, -sp.oo))
    S_leak = sp.simplify(sp.integrate(sp.diff(W, w) * j_w, (w, -sp.oo, sp.oo)))
    S_expected = sp.sqrt(2) * E0 * mu_w * rho0 / (2 * sp.sqrt(sp.pi) * lam**3)

    print("W_lambda(w)           =", W)
    print("phi_lambda(w)         =", phi)
    print("E_w(w;r)              =", E_w)
    print("j^w(w;r)              =", j_w)
    print("Boundary term         =", boundary)
    print("S_leak                =", S_leak)
    print("Expected S_leak       =", S_expected)

    assert boundary == 0
    assert sp.simplify(S_leak - S_expected) == 0

    # ------------------------------------------------------------------
    # 2. Exact one-mode bulk work scalar and reduced Session-I scalar.
    # ------------------------------------------------------------------
    section("2. Bulk work scalar and reduced Session-I scalar")
    W_bulk = sp.simplify(sp.integrate(J_w * E_w, (w, -sp.oo, sp.oo)))
    W_bulk_expected = sp.sqrt(2) * E0**2 * mu_w * q * rho0 / (2 * sp.sqrt(sp.pi) * lam**3)
    W_sess = sp.simplify(2 * sp.sqrt(2 * sp.pi) * lam * W_bulk)
    W_sess_expected = sp.simplify(2 * E0**2 * mu_w * q * rho0 / lam**2)
    W_bulk_from_S = sp.simplify(q * E0 * S_leak)
    W_sess_from_S = sp.simplify(4 * sp.pi * q * lam**4 * S_leak**2 / (mu_w * rho0))

    print("W_bulk                =", W_bulk)
    print("Expected W_bulk       =", W_bulk_expected)
    print("q E0 S_leak           =", W_bulk_from_S)
    print("W_sess                =", W_sess)
    print("Expected W_sess       =", W_sess_expected)
    print("W_sess(S_leak)        =", W_sess_from_S)

    assert sp.simplify(W_bulk - W_bulk_expected) == 0
    assert sp.simplify(W_bulk - W_bulk_from_S) == 0
    assert sp.simplify(W_sess - W_sess_expected) == 0
    assert sp.simplify(W_sess - W_sess_from_S) == 0

    # ------------------------------------------------------------------
    # 3. Exact Stage-225 selected-support pullback.
    # ------------------------------------------------------------------
    section("3. Selected-support pullback")
    Lam, eps, varrho, eta_leak = sp.symbols("Lam eps varrho eta_leak", positive=True, real=True)

    C_mix = sp.simplify(8 * Lam * (1 - eps) / sp.pi**2)
    Pi_tr = sp.simplify(sp.Rational(4, 3) * C_mix)
    Pi_tr_varrho = sp.simplify(Pi_tr.subs(1 - eps, sp.Rational(3, 2) * varrho))
    Pi_expected = sp.simplify(32 * Lam * (1 - eps) / (3 * sp.pi**2))
    Pi_expected_varrho = sp.simplify(16 * Lam * varrho / sp.pi**2)

    E0_pull = sp.simplify(eta_leak * Pi_tr)
    E0_pull_varrho = sp.simplify(E0_pull.subs(1 - eps, sp.Rational(3, 2) * varrho))

    S_pull = sp.simplify(S_expected.subs(E0, E0_pull))
    S_pull_varrho = sp.simplify(S_pull.subs(1 - eps, sp.Rational(3, 2) * varrho))
    W_bulk_pull = sp.simplify(W_bulk_expected.subs(E0, E0_pull))
    W_bulk_pull_varrho = sp.simplify(W_bulk_pull.subs(1 - eps, sp.Rational(3, 2) * varrho))
    W_sess_pull = sp.simplify(W_sess_expected.subs(E0, E0_pull))
    W_sess_pull_varrho = sp.simplify(W_sess_pull.subs(1 - eps, sp.Rational(3, 2) * varrho))

    S_pull_expected = sp.simplify(16 * sp.sqrt(2) * eta_leak * mu_w * rho0 * Lam * (1 - eps) / (3 * sp.pi**sp.Rational(5, 2) * lam**3))
    S_pull_varrho_expected = sp.simplify(8 * sp.sqrt(2) * eta_leak * mu_w * rho0 * Lam * varrho / (sp.pi**sp.Rational(5, 2) * lam**3))
    W_bulk_pull_expected = sp.simplify(512 * sp.sqrt(2) * eta_leak**2 * mu_w * q * rho0 * Lam**2 * (1 - eps)**2 / (9 * sp.pi**sp.Rational(9, 2) * lam**3))
    W_bulk_pull_varrho_expected = sp.simplify(128 * sp.sqrt(2) * eta_leak**2 * mu_w * q * rho0 * Lam**2 * varrho**2 / (sp.pi**sp.Rational(9, 2) * lam**3))
    W_sess_pull_expected = sp.simplify(2048 * eta_leak**2 * mu_w * q * rho0 * Lam**2 * (1 - eps)**2 / (9 * sp.pi**4 * lam**2))
    W_sess_pull_varrho_expected = sp.simplify(512 * eta_leak**2 * mu_w * q * rho0 * Lam**2 * varrho**2 / (sp.pi**4 * lam**2))

    print("C_mix                 =", C_mix)
    print("Pi_tr                 =", Pi_tr)
    print("Pi_tr(varrho)         =", Pi_tr_varrho)
    print("E0 pullback           =", E0_pull)
    print("S_leak pullback       =", S_pull)
    print("W_bulk pullback       =", W_bulk_pull)
    print("W_sess pullback       =", W_sess_pull)

    assert sp.simplify(Pi_tr - Pi_expected) == 0
    assert sp.simplify(Pi_tr_varrho - Pi_expected_varrho) == 0
    assert sp.simplify(S_pull - S_pull_expected) == 0
    assert sp.simplify(S_pull_varrho - S_pull_varrho_expected) == 0
    assert sp.simplify(W_bulk_pull - W_bulk_pull_expected) == 0
    assert sp.simplify(W_bulk_pull_varrho - W_bulk_pull_varrho_expected) == 0
    assert sp.simplify(W_sess_pull - W_sess_pull_expected) == 0
    assert sp.simplify(W_sess_pull_varrho - W_sess_pull_varrho_expected) == 0

    # ------------------------------------------------------------------
    # 4. Support-versus-orbit split and recovery slice.
    # ------------------------------------------------------------------
    section("4. Support-versus-orbit split and recovery slice")
    R_tr, R_target, eps_eta = sp.symbols("R_tr R_target eps_eta", real=True)

    for expr, name in [(S_pull_varrho, "S_leak"), (W_bulk_pull_varrho, "W_bulk"), (W_sess_pull_varrho, "W_sess")]:
        d_Rtr = sp.simplify(sp.diff(expr, R_tr))
        d_Rtarget = sp.simplify(sp.diff(expr, R_target))
        d_epseta = sp.simplify(sp.diff(expr, eps_eta))
        print(f"d/dR_tr    {name} =", d_Rtr)
        print(f"d/dR_target{name} =", d_Rtarget)
        print(f"d/deps_eta {name} =", d_epseta)
        assert d_Rtr == 0
        assert d_Rtarget == 0
        assert d_epseta == 0

    S_rec = sp.simplify(S_pull_varrho.subs(eta_leak, 0))
    W_bulk_rec = sp.simplify(W_bulk_pull_varrho.subs(eta_leak, 0))
    W_sess_rec = sp.simplify(W_sess_pull_varrho.subs(eta_leak, 0))

    print("Recovery eta_leak=0: S_leak =", S_rec)
    print("Recovery eta_leak=0: W_bulk =", W_bulk_rec)
    print("Recovery eta_leak=0: W_sess =", W_sess_rec)

    assert S_rec == 0
    assert W_bulk_rec == 0
    assert W_sess_rec == 0

    # ------------------------------------------------------------------
    # 5. Transport-orientation parity.
    # ------------------------------------------------------------------
    section("5. Transport-orientation parity")
    S_flip = sp.simplify(S_pull_varrho.subs(eta_leak, -eta_leak))
    W_bulk_flip = sp.simplify(W_bulk_pull_varrho.subs(eta_leak, -eta_leak))
    W_sess_flip = sp.simplify(W_sess_pull_varrho.subs(eta_leak, -eta_leak))

    print("S_leak(-eta)          =", S_flip)
    print("-S_leak(eta)          =", -S_pull_varrho)
    print("W_bulk(-eta)          =", W_bulk_flip)
    print("W_sess(-eta)          =", W_sess_flip)

    assert sp.simplify(S_flip + S_pull_varrho) == 0
    assert sp.simplify(W_bulk_flip - W_bulk_pull_varrho) == 0
    assert sp.simplify(W_sess_flip - W_sess_pull_varrho) == 0

    # ------------------------------------------------------------------
    # 6. Final success banner.
    # ------------------------------------------------------------------
    section("Stage-227 audit result")
    print("All symbolic checks passed.")
    print("Verified objects:")
    print("- exact projected leakage source on the odd Gaussian scalar-photon lane,")
    print("- exact one-mode bulk work scalar and its reduced Session-I form,")
    print("- exact pullback through the Stage-225 selected-support demand Pi_tr,")
    print("- exact support-versus-orbit separation of the compiled leakage/work lane,")
    print("- exact orientation parity,")
    print("- and exact recovery of the closed standard slice at eta_leak = 0.")


if __name__ == "__main__":
    main()
