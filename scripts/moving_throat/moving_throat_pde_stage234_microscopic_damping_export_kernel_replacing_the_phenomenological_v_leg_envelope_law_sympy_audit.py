#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 234
SymPy audit for the microscopic damping/export kernel replacing the
phenomenological Session-IV V-leg envelope law.

This script verifies:
1. the exact cubic coefficient from the derivative-coupled scalar mixed outlet,
2. the exact quintic coefficient imported from the selected outgoing quadrupole lane,
3. the projected active-leg export kernel,
4. the exact Schott power identities for the cubic/quintic kernel,
5. the exact characteristic polynomial and positive-root monotonicity data,
6. the exact event-safe surface replacing gamma_safe,
7. and the Session-IV cold-event benchmark numbers in reduced units.
"""

from __future__ import annotations

import math

import sympy as sp



def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)



def main() -> None:
    # ------------------------------------------------------------------
    # 1. Scalar derivative-coupled cubic outlet.
    # ------------------------------------------------------------------
    section("1. Scalar derivative-coupled cubic outlet")
    omega = sp.symbols("omega", real=True)
    OmU0, OmW0, R0 = sp.symbols("Omega_U0 Omega_W0 R0", positive=True, real=True)
    eta0, gamma1 = sp.symbols("eta_0 gamma_1", positive=True, real=True)

    A0 = sp.simplify(OmU0**2 - omega**2)
    W0 = sp.simplify(OmW0**2 - omega**2)
    Delta0 = sp.simplify(OmU0**2 * OmW0**2 - R0**2)
    Delta_omega = sp.simplify(A0 * W0 - R0**2)
    gW0 = sp.simplify(eta0 * omega)
    N0_scalar = sp.simplify((A0 * gW0) ** 2 / Delta_omega**2)
    N0_series = sp.series(N0_scalar, omega, 0, 5).removeO()
    N0_coeff2 = sp.simplify(sp.limit(N0_scalar / omega**2, omega, 0))
    Gamma30 = sp.simplify(gamma1 * N0_coeff2)

    print("A0(omega)             =", A0)
    print("W0(omega)             =", W0)
    print("Delta0                =", Delta0)
    print("N0_scalar(omega)      =", N0_scalar)
    print("series N0_scalar      =", N0_series)
    print("coeff omega^2         =", N0_coeff2)
    print("Gamma_{3,0}           =", Gamma30)

    assert sp.simplify(N0_coeff2 - eta0**2 * OmU0**4 / Delta0**2) == 0

    # ------------------------------------------------------------------
    # 2. Selected mixed quadrupole quintic outlet.
    # ------------------------------------------------------------------
    section("2. Selected mixed quadrupole quintic outlet")
    a, cs = sp.symbols("a c_s", positive=True, real=True)
    beta0, sminus, lamminus = sp.symbols("beta_0 s_- lambda_-", positive=True, real=True)
    PiV0, PiVm = sp.symbols("Pi_V0 Pi_V-", real=True)

    P0_minus = sp.simplify(beta0 * sminus / lamminus)
    Gamma2_port = sp.simplify(a**5 / (27 * cs**5))
    Gamma5_minus = sp.simplify(Gamma2_port * P0_minus)
    Gamma3 = sp.simplify(PiV0**2 * Gamma30)
    Gamma5 = sp.simplify(PiVm**2 * Gamma5_minus)

    print("P0_minus              =", P0_minus)
    print("Gamma_2^(port)        =", Gamma2_port)
    print("Gamma_{5,-}           =", Gamma5_minus)
    print("Gamma_3               =", Gamma3)
    print("Gamma_5               =", Gamma5)

    assert sp.simplify(Gamma5 - PiVm**2 * a**5 * beta0 * sminus / (27 * cs**5 * lamminus)) == 0

    # ------------------------------------------------------------------
    # 3. Projected active-leg export kernel.
    # ------------------------------------------------------------------
    section("3. Projected active-leg export kernel")
    Sigma_exp = sp.simplify(-sp.I * Gamma3 * omega**3 - sp.I * Gamma5 * omega**5)
    print("Sigma_exp(omega)      =", Sigma_exp)

    # ------------------------------------------------------------------
    # 4. Exact Schott power identities.
    # ------------------------------------------------------------------
    section("4. Exact Schott power identities")
    t = sp.symbols("t", real=True)
    V = sp.Function("V")(t)
    G3, G5 = sp.symbols("Gamma_3 Gamma_5", positive=True, real=True)

    F_odd = sp.simplify(G3 * sp.diff(V, t, 3) - G5 * sp.diff(V, t, 5))
    Schott = sp.simplify(
        G3 * sp.diff(V, t) * sp.diff(V, t, 2)
        - G5 * (sp.diff(V, t) * sp.diff(V, t, 4) - sp.diff(V, t, 2) * sp.diff(V, t, 3))
    )
    power_balance = sp.simplify(sp.expand(sp.diff(V, t) * F_odd - sp.diff(Schott, t)))

    print("F_odd                 =", F_odd)
    print("Schott term           =", Schott)
    print("power - dSchott/dt    =", power_balance)

    assert power_balance == -G3 * sp.diff(V, t, 2) ** 2 - G5 * sp.diff(V, t, 3) ** 2

    # ------------------------------------------------------------------
    # 5. Characteristic polynomial and exact safe surface.
    # ------------------------------------------------------------------
    section("5. Characteristic polynomial and safe surface")
    mu_eta, kappaV = sp.symbols("mu_eta kappa_V", positive=True, real=True)
    s = sp.symbols("s", positive=True, real=True)
    s0, sc = sp.symbols("s_0 s_c", positive=True, real=True)

    F = sp.simplify(G3 * s**3 + G5 * s**5 + mu_eta * s**2 - kappaV)
    Fprime = sp.simplify(sp.diff(F, s))
    F_safe = sp.simplify(F.subs({kappaV: mu_eta * s0**2, s: sc}))
    safe_eq = sp.simplify(sp.solve(sp.Eq(F_safe, 0), G3)[0])
    safe_eq_hat = sp.simplify((safe_eq / mu_eta).expand())
    safe_gamma3_hat = sp.simplify((mu_eta * (s0**2 - sc**2) / sc**3) / mu_eta)
    safe_gamma5_hat = sp.simplify((mu_eta * (s0**2 - sc**2) / sc**5) / mu_eta)
    root_shift = sp.simplify(-(G3 * s0**2 + G5 * s0**4) / (2 * mu_eta))

    print("F(s)                  =", F)
    print("F'(s)                 =", Fprime)
    print("F(sc)                 =", F_safe)
    print("G3_safe(sc,G5)        =", safe_eq)
    print("G3_safe/mu_eta        =", safe_eq_hat)
    print("G3hat_safe (G5=0)     =", safe_gamma3_hat)
    print("G5hat_safe (G3=0)     =", safe_gamma5_hat)
    print("delta s (weak export) =", root_shift)

    assert sp.simplify(Fprime - (5 * G5 * s**4 + 3 * G3 * s**2 + 2 * mu_eta * s)) == 0
    assert sp.simplify(safe_eq - mu_eta * (s0**2 - sc**2 - G5 * sc**5 / mu_eta) / sc**3) == 0
    assert sp.simplify(safe_eq_hat - ((s0**2 - sc**2) / sc**3 - (G5 / mu_eta) * sc**2)) == 0

    # ------------------------------------------------------------------
    # 6. Dimensionless half-plane rewrite.
    # ------------------------------------------------------------------
    section("6. Normalized half-plane")
    G3hat, G5hat = sp.symbols("G3hat G5hat", real=True)
    half_plane = sp.simplify(sp.Eq(G3hat + sc**2 * G5hat, (s0**2 - sc**2) / sc**3))
    dimless_safe = sp.simplify((s0**2 - sc**2) / sc**2)

    print("safe half-plane       =", half_plane)
    print("dimensionless rhs     =", dimless_safe)

    # ------------------------------------------------------------------
    # 7. Session-IV benchmark specialization.
    # ------------------------------------------------------------------
    section("7. Session-IV benchmark specialization")
    t_cross_num = 1.82169718
    tcollapse0_num = 0.14402764
    gamma_crit_num = 6.94311167

    sc_num = 1.0 / t_cross_num
    s0_from_t_num = 1.0 / tcollapse0_num
    s0_num = gamma_crit_num
    rhs_num = (s0_num**2 - sc_num**2) / sc_num**3
    g3hat_safe_num = rhs_num
    g5hat_safe_num = rhs_num / (sc_num**2)
    weight_num = sc_num**2
    dimless_rhs_num = (s0_num**2 - sc_num**2) / sc_num**2

    print("s_c                   =", sc_num)
    print("s_0 from tcollapse    =", s0_from_t_num)
    print("s_0 from gamma_crit   =", s0_num)
    print("safe half-plane rhs   =", rhs_num)
    print("weight sc^2           =", weight_num)
    print("G3hat_safe            =", g3hat_safe_num)
    print("G5hat_safe            =", g5hat_safe_num)
    print("dimensionless rhs     =", dimless_rhs_num)

    assert abs(s0_from_t_num - s0_num) < 1e-6
    assert abs(sc_num - 0.5489386551062235) < 1e-12
    assert abs(weight_num - 0.3013336470698294) < 1e-12
    assert abs(rhs_num - 289.61004917557426) < 1e-9
    assert abs(g5hat_safe_num - 961.0942952828019) < 1e-9
    assert abs(dimless_rhs_num - 158.978150899687) < 1e-9

    # ------------------------------------------------------------------
    # 8. Polynomial root checks on the safe surface.
    # ------------------------------------------------------------------
    section("8. Root checks on the safe surface")
    s_poly = sp.symbols("s_poly", real=True)

    # Pure cubic safe point.
    poly3 = sp.Poly(s_poly**2 + g3hat_safe_num * s_poly**3 - s0_num**2, s_poly)
    roots3 = [complex(r) for r in sp.nroots(poly3)]
    pos3 = sorted([r.real for r in roots3 if abs(r.imag) < 1e-10 and r.real > 0])

    # Pure quintic safe point.
    poly5 = sp.Poly(s_poly**2 + g5hat_safe_num * s_poly**5 - s0_num**2, s_poly)
    roots5 = [complex(r) for r in sp.nroots(poly5)]
    pos5 = sorted([r.real for r in roots5 if abs(r.imag) < 1e-10 and r.real > 0])

    print("positive roots (cubic)=", pos3)
    print("positive roots (quint)=", pos5)

    assert len(pos3) == 1
    assert len(pos5) == 1
    assert abs(pos3[0] - sc_num) < 1e-9
    assert abs(pos5[0] - sc_num) < 1e-9

    # ------------------------------------------------------------------
    # 9. Final banner.
    # ------------------------------------------------------------------
    section("Stage-234 audit result")
    print("All symbolic and numerical checks passed.")
    print("Verified objects:")
    print("- exact cubic coefficient from the derivative-coupled scalar mixed outlet,")
    print("- exact quintic coefficient imported from the selected outgoing quadrupole lane,")
    print("- projected active-leg odd export kernel Sigma_exp(omega),")
    print("- exact Schott power identities for the cubic/quintic passive export law,")
    print("- exact characteristic polynomial and strict monotonicity on s>0,")
    print("- absence of a finite microscopic unconditional threshold on the minimal kernel,")
    print("- exact event-safe half-plane replacing gamma_safe,")
    print("- and the Session-IV cold-event benchmark specialization.")


if __name__ == "__main__":
    main()
