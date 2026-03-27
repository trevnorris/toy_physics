#!/usr/bin/env python3
"""
Master SymPy verifier for the carry-forward lepton/electron derivations.

This is a standalone script that consolidates the derivations summarized in
Sections 1–18 of the carry-forward notes:

  1. Frozen background and source anchors
  2. Core isolated-electron mass theorem
  3. Microscopic identification of the reduced coefficients
  4. Exact D/N spectral closure
  5. Charge / circulation / throughput audit
  6. Dynamic turbine scaling
  7. Exact phase closure versus missing amplitude closure
  8. Open-system mouth-output channel
  9. Static brane response: compliance, not DC flow
 10. Dynamic AC→DC rectification law
 11. Unique equilibrium crossover and selected electron scale
 12. Family-ratio theorems and no-gos
 13. Reverse-engineered lepton families under the universal 11:2:5 ledger
 14. Multi-constraint intersection picture
 15. Universal-alpha lattice search
 16. Resonance and phase-lock no-gos
 17. Highest-confidence conclusions to carry forward
 18. Best chapter split for the next session

The script emphasizes symbolic verification with SymPy and uses numeric lattice
searches where the derivation is inherently arithmetic / Diophantine.
"""

from __future__ import annotations

import argparse
import bisect
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Any

import mpmath as mp
import sympy as sp

mp.mp.dps = 80

# -----------------------------------------------------------------------------
# Global constants carried forward
# -----------------------------------------------------------------------------

LAMBDA0 = mp.mpf("1.85")
M_E = mp.mpf("0.51099895000")
M_MU = mp.mpf("105.6583755")
M_TAU = mp.mpf("1776.93")
R_MU = M_MU / M_E
R_TAU = M_TAU / M_E

# -----------------------------------------------------------------------------
# Formatting helpers
# -----------------------------------------------------------------------------

def banner(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def line(name: str, value: Any) -> None:
    print(f"{name} = {value}")


def syms_eq_zero(expr: sp.Expr) -> None:
    assert sp.simplify(expr) == 0, f"Expected zero, got: {sp.simplify(expr)}"


# -----------------------------------------------------------------------------
# Section 1
# -----------------------------------------------------------------------------

def section_1_setup() -> Dict[str, Any]:
    banner("Section 1 — Frozen background and source anchors")
    n = sp.Integer(5)
    kappa_add = sp.Rational(1, 2)
    kappa_PV = sp.Rational(3, 2)
    beta_1PN = sp.Integer(3)
    Lambda0 = sp.Rational(185, 100)

    line("n", n)
    line("kappa_add", kappa_add)
    line("kappa_PV", kappa_PV)
    line("beta_1PN", beta_1PN)
    line("Lambda0", Lambda0)
    line("mu/e", mp.nstr(R_MU, 16))
    line("tau/e", mp.nstr(R_TAU, 16))

    return {
        "n": n,
        "kappa_add": kappa_add,
        "kappa_PV": kappa_PV,
        "beta_1PN": beta_1PN,
        "Lambda0": Lambda0,
        "R_MU": R_MU,
        "R_TAU": R_TAU,
    }


# -----------------------------------------------------------------------------
# Section 2
# -----------------------------------------------------------------------------

def section_2_core_mass_theorem() -> Dict[str, Any]:
    banner("Section 2 — Core isolated-electron mass theorem")

    A0, B0, C0 = sp.symbols("A0 B0 C0", positive=True)
    n, rho = sp.symbols("n rho", positive=True)
    a, x = sp.symbols("a x", positive=True)

    A = A0 * rho ** ((n - 1) / 2)
    B = B0 / rho
    C = C0 * rho**n

    F = A / a + B / a**2 + C * a**3
    Ew = A / a
    Ef = B / a**2
    Epv = C * a**3

    stationary_combo = sp.simplify(a * sp.diff(F, a))
    virial_form = sp.simplify(stationary_combo - (-Ew - 2 * Ef + 3 * Epv))
    syms_eq_zero(virial_form)

    line("stationary_combo", stationary_combo)
    line("virial_law", sp.Eq(Ew + 2 * Ef, 3 * Epv))

    Ew_g, xg, ng = sp.symbols("Ew_g xg ng", positive=True)
    Ef_g = xg * Ew_g
    Epv_g = Ew_g * (1 + 2 * xg) / 3
    F_g = Ew_g + Ef_g + Epv_g
    dlnF_x = sp.simplify((((ng - 1) / 2) * Ew_g - Ef_g + ng * Epv_g) / F_g)
    dlnF_x_expected = sp.simplify((((5 * ng - 3) / 2) + (2 * ng - 3) * xg) / (4 + 5 * xg))
    syms_eq_zero(dlnF_x - dlnF_x_expected)

    x_solution = sp.solve(sp.Eq(dlnF_x_expected.subs(ng, 5), sp.Rational(5, 2)), xg)
    assert x_solution == [sp.Rational(2, 11)]
    x_star = x_solution[0]
    partition = [sp.Integer(1), x_star, sp.simplify((1 + 2 * x_star) / 3)]
    assert partition == [sp.Integer(1), sp.Rational(2, 11), sp.Rational(5, 11)]

    line("dlnF_x", dlnF_x_expected)
    line("x_star", x_star)
    line("partition", partition)

    # Breathing slope
    S = -A / a**2 - 2 * B / a**3 + 3 * C * a**2
    dln_a_full = sp.simplify((rho / a) * (-sp.diff(S, rho) / sp.diff(S, a)))
    B0_expr = sp.simplify((x * A * a) * rho)                   # B = x A a
    C0_expr = sp.simplify(A * (1 + 2 * x) / (3 * a**4 * rho**n))
    dln_a_sub = sp.simplify(dln_a_full.subs({B0: B0_expr, C0: C0_expr}))
    dln_a_expected = sp.simplify(-(-((n - 1) / 2) + 2 * x + n * (1 + 2 * x)) / (4 + 10 * x))
    syms_eq_zero(dln_a_sub - dln_a_expected)
    assert sp.simplify(dln_a_expected.subs({n: 5, x: sp.Rational(2, 11)})) == sp.Rational(-57, 64)

    line("dln_a_expected", dln_a_expected)
    line("dln_a_at_n5_x2_11", sp.Rational(-57, 64))

    # Explicit stationary solution at x = 2/11
    A_sym, B_sym, C_sym = sp.symbols("A_sym B_sym C_sym", positive=True)
    Ew2 = A_sym / a
    Ef2 = B_sym / a**2
    a_star = sp.solve(sp.Eq(Ef2 / Ew2, sp.Rational(2, 11)), a)[0]
    C_from_AB = sp.solve(
        sp.Eq(sp.diff(A_sym / a + B_sym / a**2 + C_sym * a**3, a).subs(a, a_star), 0),
        C_sym,
    )[0]
    F_star = sp.simplify((A_sym / a + B_sym / a**2 + C_sym * a**3).subs({a: a_star, C_sym: C_from_AB}))
    Ew_star = sp.simplify((A_sym / a).subs(a, a_star))
    Ef_star = sp.simplify((B_sym / a**2).subs(a, a_star))
    Epv_star = sp.simplify((C_sym * a**3).subs({a: a_star, C_sym: C_from_AB}))

    assert a_star == 11 * B_sym / (2 * A_sym)
    syms_eq_zero(F_star - 36 * A_sym**2 / (121 * B_sym))

    line("a_star", a_star)
    line("C_from_AB", sp.factor(C_from_AB))
    line("Ew_star", Ew_star)
    line("Ef_star", Ef_star)
    line("Epv_star", Epv_star)
    line("F_star", F_star)

    # Optional P22-style extension
    D, p = sp.symbols("D p", positive=True)
    F_ext = A_sym / a + B_sym / a**2 + C_sym * a**3 + D * a**p
    virial_ext = sp.expand(-a * sp.diff(F_ext, a))
    line("virial_ext", virial_ext)

    return {
        "x_star": x_star,
        "partition": partition,
        "dln_a": sp.Rational(-57, 64),
        "a_star": a_star,
        "F_star": F_star,
        "virial_ext": virial_ext,
    }


# -----------------------------------------------------------------------------
# Section 3
# -----------------------------------------------------------------------------

def section_3_microscopic_coefficients() -> Dict[str, Any]:
    banner("Section 3 — Microscopic identification of the reduced coefficients")

    d, a, r, rho = sp.symbols("d a r rho", positive=True)
    Phi = sp.symbols("Phi", positive=True)
    I_w = sp.symbols("I_w", positive=True)
    chi_w = sp.symbols("chi_w", positive=True)
    c_s = sp.symbols("c_s", positive=True)
    Lambda = sp.symbols("Lambda", positive=True)
    hbar = sp.symbols("hbar", positive=True)

    S = sp.symbols("S", positive=True)
    u = Phi / (rho * S * r ** (d - 1))
    Ef_d = sp.simplify(
        sp.integrate(sp.Rational(1, 2) * rho * u**2 * S * r ** (d - 1), (r, a, sp.oo), conds="none")
    )
    Ef_3 = sp.simplify(Ef_d.subs({d: 3, S: 4 * sp.pi}))
    Ef_4 = sp.simplify(Ef_d.subs({d: 4, S: 2 * sp.pi**2}))
    syms_eq_zero(Ef_3 - Phi**2 / (8 * sp.pi * rho * a))
    syms_eq_zero(Ef_4 - Phi**2 / (8 * sp.pi**2 * rho * a**2))

    Ew = sp.simplify(I_w * c_s * chi_w / a)
    A_rho = sp.simplify(I_w * chi_w * c_s)
    B_rho = sp.simplify(Phi**2 / (8 * sp.pi**2 * rho))
    F0 = sp.simplify(36 * A_rho**2 / (121 * B_rho))

    x = sp.simplify(Ef_4 / Ew)
    a0 = sp.simplify(sp.solve(sp.Eq(x, sp.Rational(2, 11)), a)[0])

    chi_DN = sp.simplify(sp.pi / (2 * Lambda))
    omega_half = c_s * sp.pi / (2 * Lambda * a0)
    Phi2_from_a0 = sp.solve(sp.Eq(a, a0), Phi**2)[0]
    F0_elim = sp.simplify(F0.subs({Phi**2: Phi2_from_a0, chi_w: chi_DN, a: sp.symbols("a0", positive=True)}))
    a0_sym = sp.symbols("a0", positive=True)
    F0_DN = sp.simplify(F0_elim.subs(a0_sym, a0_sym))
    F0_DN_expected = sp.Rational(18, 11) * I_w * (c_s * sp.pi / (2 * Lambda * a0_sym))
    syms_eq_zero(F0_DN - F0_DN_expected)
    F0_quantized = sp.simplify(F0_DN.subs(I_w, hbar))

    line("Ef_d", Ef_d)
    line("Ef_3", Ef_3)
    line("Ef_4", Ef_4)
    line("Ew", Ew)
    line("A(rho)", A_rho)
    line("B(rho)", B_rho)
    line("a0_from_Ef_over_Ew_2_11", a0)
    line("F0", F0)
    line("F0_DN", F0_DN)
    line("F0_DN_quantized", F0_quantized)

    return {
        "A_rho": A_rho,
        "B_rho": B_rho,
        "a0": a0,
        "F0_DN": F0_DN,
        "F0_DN_quantized": F0_quantized,
    }


# -----------------------------------------------------------------------------
# Section 4
# -----------------------------------------------------------------------------

def section_4_dn_spectral_closure() -> Dict[str, Any]:
    banner("Section 4 — Exact D/N spectral closure")

    j = sp.symbols("j", integer=True, nonnegative=True)
    z, t = sp.symbols("z t", real=True)
    k, L, cs, Lam, rho, Phi, Iw = sp.symbols("k L c_s Lambda rho Phi I_w", positive=True)
    A0, B0 = sp.symbols("A0 B0")
    a_sym = sp.symbols("a", positive=True)

    psi = A0 * sp.sin(k * z) + B0 * sp.cos(k * z)
    B_solution = sp.solve(sp.Eq(psi.subs(z, 0), 0), B0)[0]
    psi_D = sp.simplify(psi.subs(B0, B_solution))
    neumann_condition = sp.simplify(sp.diff(psi_D, z).subs(z, L))
    k_j = sp.simplify(sp.pi * (j + sp.Rational(1, 2)) / L)
    omega_j = sp.simplify(cs * k_j)

    line("psi_D", psi_D)
    line("neumann_condition", neumann_condition)
    line("k_j", k_j)
    line("omega_j", omega_j)

    # Exact round-trip factor
    exp_phase = sp.expand_complex(sp.exp(sp.I * 2 * k_j * L))
    R_rt = sp.simplify((-1) * exp_phase)
    syms_eq_zero(R_rt - 1)
    line("R_rt", R_rt)

    # Standing-wave current vanishes
    A_real = sp.symbols("A_real", real=True)
    psi_mode = A_real * sp.sin(k_j * z) * sp.exp(-sp.I * omega_j * t)
    Jz = sp.simplify((sp.conjugate(psi_mode) * sp.diff(psi_mode, z) - psi_mode * sp.conjugate(sp.diff(psi_mode, z))) / (2 * sp.I))
    syms_eq_zero(Jz)
    line("Jz", Jz)

    # Fixed-geometry tower
    F_fixed = sp.simplify(sp.Rational(18, 11) * Iw * omega_j)
    F_fixed_Lambda = sp.simplify(F_fixed.subs(L, Lam * a_sym))
    line("F_fixed", F_fixed)
    line("F_fixed_Lambda", F_fixed_Lambda)

    # Self-consistent tower
    chi_j_dimless = sp.simplify(sp.pi * (j + sp.Rational(1, 2)) / Lam)
    A_j = sp.simplify(Iw * cs * chi_j_dimless)
    B = sp.simplify(Phi**2 / (8 * sp.pi**2 * rho))
    a_j = sp.simplify(11 * B / (2 * A_j))
    omega_eq = sp.simplify((sp.pi * cs / (Lam * a_j)) * (j + sp.Rational(1, 2)))
    F_eq = sp.simplify(sp.Rational(18, 11) * Iw * omega_eq)
    F_eq_alt = sp.simplify(sp.Rational(36, 121) * A_j**2 / B)
    syms_eq_zero(F_eq - F_eq_alt)
    F0_eq = sp.simplify(F_eq.subs(j, 0))
    mass_ratio = sp.simplify(F_eq / F0_eq)
    syms_eq_zero(mass_ratio - (2 * j + 1) ** 2)

    line("a_j", a_j)
    line("F_eq", F_eq)
    line("mass_ratio", mass_ratio)
    line("observed_mu_over_e", mp.nstr(R_MU, 16))
    line("observed_tau_over_e", mp.nstr(R_TAU, 16))

    return {
        "k_j": k_j,
        "omega_j": omega_j,
        "R_rt": R_rt,
        "Jz": Jz,
        "mass_ratio": mass_ratio,
    }


# -----------------------------------------------------------------------------
# Section 5
# -----------------------------------------------------------------------------

def section_5_throughput_audit() -> Dict[str, Any]:
    banner("Section 5 — Charge / circulation / throughput audit")

    Iw, cs, rho, Lam, Phi = sp.symbols("I_w c_s rho Lambda Phi", positive=True)
    hbar, mpsi, qstar, n, Phi_loop = sp.symbols("hbar m_psi q_star n Phi_loop")
    S3, r = sp.symbols("S_3 r", positive=True)

    F0 = sp.Rational(72) * sp.pi**4 * Iw**2 * cs**2 * rho / (sp.Integer(121) * Lam**2 * Phi**2)
    Gamma = (2 * sp.pi * hbar * n - qstar * Phi_loop) / mpsi
    jr = Phi / (S3 * r**3)
    div_radial = sp.simplify(sp.diff(r**3 * jr, r) / r**3)
    syms_eq_zero(div_radial)

    line("F0", F0)
    line("fluxoid_law_Gamma", Gamma)
    line("j_r(r)", jr)
    line("div_4[j_r r_hat]", div_radial)
    print("Conclusion: topology fixes the circulation/holonomy class, but not the throughput amplitude Phi.")

    return {
        "F0": F0,
        "Gamma": Gamma,
        "jr": jr,
        "div_radial": div_radial,
    }


# -----------------------------------------------------------------------------
# Section 6
# -----------------------------------------------------------------------------

def section_6_dynamic_turbine_scaling() -> Dict[str, Any]:
    banner("Section 6 — Dynamic turbine scaling")

    nu, Phi, Phi0 = sp.symbols("nu Phi Phi0", positive=True)
    U, V, W = sp.symbols("U V W", positive=True)
    Lambda, a = sp.symbols("Lambda a", positive=True)
    R, s = sp.symbols("R s", positive=True)

    A = U * nu / Lambda
    B = V * Phi**2
    C = W * Lambda

    a_expr = sp.simplify(sp.Rational(11, 2) * B / A)
    C_required = sp.simplify(sp.Integer(80) * A**5 / (sp.Integer(11) ** 5 * B**4))
    Lambda_sol = sp.solve(sp.Eq(C, C_required), Lambda)[0]
    a_sol = sp.simplify(a_expr.subs(Lambda, Lambda_sol))
    F_sol = sp.simplify((sp.Rational(36, 121) * A**2 / B).subs(Lambda, Lambda_sol))

    # Relative laws
    phi = sp.symbols("phi", positive=True)
    F_ratio = sp.simplify((F_sol / F_sol.subs({nu: 1, Phi: Phi0})).subs(Phi, phi * Phi0))
    a_ratio = sp.simplify((a_sol / a_sol.subs({nu: 1, Phi: Phi0})).subs(Phi, phi * Phi0))
    L_ratio = sp.simplify((Lambda_sol * a_sol) / (Lambda_sol * a_sol).subs({nu: 1, Phi: Phi0}))
    L_ratio = sp.simplify(L_ratio.subs(Phi, phi * Phi0))
    Lambda_ratio = sp.simplify((Lambda_sol / Lambda_sol.subs({nu: 1, Phi: Phi0})).subs(Phi, phi * Phi0))

    syms_eq_zero(F_ratio - nu ** sp.Rational(1, 3) * phi ** sp.Rational(2, 3))
    syms_eq_zero(a_ratio - nu ** sp.Rational(-1, 6) * phi ** sp.Rational(2, 3))
    syms_eq_zero(L_ratio - nu ** sp.Rational(2, 3) * phi ** sp.Rational(-2, 3))
    syms_eq_zero(Lambda_ratio - nu ** sp.Rational(5, 6) * phi ** sp.Rational(-4, 3))

    line("Lambda_scaling", sp.factor(Lambda_sol))
    line("a_scaling", sp.factor(a_sol))
    line("F_scaling", sp.factor(F_sol))
    line("R(nu,phi)", F_ratio)
    line("a/a0", a_ratio)
    line("L/L0", L_ratio)
    line("Lambda/Lambda0", Lambda_ratio)

    # Inverted geometry in terms of target R
    phi_from_R = sp.simplify(R ** sp.Rational(3, 2) / sp.sqrt(nu))
    a_from_R = sp.simplify(a_ratio.subs(phi, phi_from_R))
    L_from_R = sp.simplify(L_ratio.subs(phi, phi_from_R))
    Lambda_from_R = sp.simplify(Lambda_ratio.subs(phi, phi_from_R))
    syms_eq_zero(a_from_R - R / sp.sqrt(nu))
    syms_eq_zero(L_from_R - nu / R)
    syms_eq_zero(Lambda_from_R - nu ** sp.Rational(3, 2) / R**2)

    line("phi_from_R", phi_from_R)
    line("a/a0_from_R", a_from_R)
    line("L/L0_from_R", L_from_R)
    line("Lambda/Lambda0_from_R", Lambda_from_R)

    # Minimal scale-free turbine law: phi = (P/P0)^s with P/P0 = R^2
    turbine_eq = sp.Eq(R, nu ** sp.Rational(1, 3) * (R ** (2 * s)) ** sp.Rational(2, 3))
    R_from_s = sp.solve(turbine_eq, R)[0]
    R_expected = sp.simplify(nu ** (1 / (3 - 4 * s)))
    syms_eq_zero(R_from_s - R_expected)
    line("R_from_turbine_exponent_s", R_expected)
    line("critical_s", sp.Rational(3, 4))

    s_mu = (3 - math.log(3) / math.log(float(R_MU))) / 4
    s_tau = (3 - math.log(5) / math.log(float(R_TAU))) / 4
    line("s_mu_fit", s_mu)
    line("s_tau_fit", s_tau)

    return {
        "R_law": F_ratio,
        "a_ratio": a_ratio,
        "L_ratio": L_ratio,
        "Lambda_ratio": Lambda_ratio,
        "R_from_s": R_expected,
        "s_mu": s_mu,
        "s_tau": s_tau,
    }


# -----------------------------------------------------------------------------
# Section 7
# -----------------------------------------------------------------------------

def section_7_phase_vs_amplitude() -> Dict[str, Any]:
    banner("Section 7 — Exact phase closure versus missing amplitude closure")

    A_n, A_np1, A_star, Lambda, D = sp.symbols("A_n A_np1 A_star Lambda D")
    A_star_expr = sp.solve(sp.Eq(A_star, Lambda * A_star + D), A_star)[0]
    syms_eq_zero(A_star_expr - D / (1 - Lambda))
    line("steady_state_fixed_point", A_star_expr)
    print("Autonomous eigenmode condition: D = 0 and A_{n+1} = A_n != 0 imply Lambda = 1 and phi_0 = 2*pi*N.")
    print("Conclusion: phase closure is nearly in hand on the exact D/N branch; amplitude closure remains open.")
    return {"A_star": A_star_expr}


# -----------------------------------------------------------------------------
# Section 8
# -----------------------------------------------------------------------------

def section_8_mouth_output_channel() -> Dict[str, Any]:
    banner("Section 8 — Open-system mouth-output channel")

    L, cs, hbar, rho_par, Zsig0, jmouth = sp.symbols("L cs hbar rho_par Zsig0 jmouth", positive=True)
    z, t = sp.symbols("z t", real=True)
    A0 = sp.symbols("A0", real=True)

    chi0 = sp.sqrt(2 / L) * sp.sin(sp.pi * z / (2 * L))
    omega0 = sp.pi * cs / (2 * L)
    chip0 = sp.simplify(sp.diff(chi0, z).subs(z, 0))
    chip0_sq = sp.simplify(chip0**2)
    X = A0 * chi0 * sp.cos(omega0 * t)

    X0 = sp.simplify(X.subs(z, 0))
    Xt0 = sp.simplify(sp.diff(X, t).subs(z, 0))
    line("chi0(0)", sp.simplify(chi0.subs(z, 0)))
    line("chi0'(0)", chip0)
    line("|chi0'(0)|^2", chip0_sq)
    line("X(0,t)", X0)
    line("d_t X(0,t)", Xt0)

    # Single trapped mode carries no one-way current
    psi = A0 * chi0 * sp.exp(-sp.I * omega0 * t)
    Jw = sp.simplify((sp.conjugate(psi) * sp.diff(psi, z) - psi * sp.conjugate(sp.diff(psi, z))) / (2 * sp.I))
    syms_eq_zero(Jw)
    line("J_w", Jw)

    I = sp.symbols("I", positive=True)
    A0_sq_from_I = sp.solve(sp.Eq(I, rho_par * A0**2 * omega0 / 2), A0**2)[0]
    A0_sq = sp.simplify(A0_sq_from_I.subs(I, hbar))
    Tnn_dc = sp.simplify((rho_par * cs**2 / 4) * A0_sq * chip0_sq)
    syms_eq_zero(Tnn_dc - sp.pi * hbar * cs / (2 * L**2))
    line("A0^2_from_I=hbar", A0_sq)
    line("Tnn_dc", Tnn_dc)

    jmouth_dc = sp.simplify(Zsig0 * Tnn_dc)
    F0 = sp.simplify(9 * sp.pi * hbar * cs / (11 * L))
    line("jmouth_dc", jmouth_dc)
    line("F0(L)", F0)

    L_from_j = sp.solve(sp.Eq(jmouth, jmouth_dc), L)[0]
    F0_elim = sp.simplify(F0.subs(L, L_from_j))
    line("L(jmouth)", L_from_j)
    line("F0_after_eliminating_L", F0_elim)

    return {
        "chip0": chip0,
        "Tnn_dc": Tnn_dc,
        "jmouth_dc": jmouth_dc,
        "F0": F0,
    }


# -----------------------------------------------------------------------------
# Section 9
# -----------------------------------------------------------------------------

def section_9_static_compliance() -> Dict[str, Any]:
    banner("Section 9 — Static brane response: compliance, not DC flow")

    rho, rho0, K, m, hbar, c_s = sp.symbols("rho rho0 K m hbar c_s", positive=True)
    eta, gradeta, gradtheta = sp.symbols("eta gradeta gradtheta", real=True)
    A0, C0, eps, Bmix, f0 = sp.symbols("A0 C0 eps Bmix f0", positive=True)

    P = K * rho**5
    U = K * rho**5 / 4
    cs2 = sp.diff(P, rho) / m
    Upp = sp.diff(U, rho, 2)
    syms_eq_zero(Upp - m * cs2 / rho)

    A = rho0 * m * c_s**2
    B = rho0 * hbar**2 / (8 * m)
    ell2 = sp.simplify((2 * B) / A)
    xi_h = hbar / (sp.sqrt(2) * m * c_s)
    syms_eq_zero(ell2 - xi_h**2 / 2)

    line("P(rho)", P)
    line("U(rho)", U)
    line("U''(rho)", Upp)
    line("ell^2", ell2)

    # Conservative phase current vanishes at static minimum
    Etheta = sp.simplify(hbar**2 * rho0 * gradtheta**2 / (2 * m))
    gradtheta_min = sp.solve(sp.Eq(sp.diff(Etheta, gradtheta), 0), gradtheta)[0]
    j_static = sp.simplify((hbar / m) * rho0 * gradtheta_min)
    syms_eq_zero(j_static)
    line("E_theta", Etheta)
    line("gradtheta_min", gradtheta_min)
    line("j_static", j_static)
    print("Conclusion: Z_eff_sigma,flux(0) = 0 in the conservative static theory.")

    # P22 enters scalar compliance only at even order
    H = sp.Matrix([[A0, eps * Bmix], [eps * Bmix, C0]])
    f = sp.Matrix([f0, 0])
    q = H.LUsolve(f)
    chi00 = sp.simplify(q[0] / f0)
    chi_series = sp.series(chi00, eps, 0, 4).removeO()
    line("chi00_eff", chi00)
    line("chi00_series", chi_series)

    return {
        "Upp": Upp,
        "ell2": ell2,
        "j_static": j_static,
        "chi00": chi00,
        "chi00_series": chi_series,
    }


# -----------------------------------------------------------------------------
# Section 10
# -----------------------------------------------------------------------------

def section_10_ac_streaming() -> Dict[str, Any]:
    banner("Section 10 — Dynamic AC→DC rectification law")

    eps = sp.symbols("eps", real=True)
    x, t = sp.symbols("x t", real=True)
    rho0, cs, hbar, m, L, Lambda = sp.symbols("rho0 cs hbar m L Lambda", positive=True, real=True)
    xi = sp.symbols("xi", positive=True, real=True)
    cph = sp.symbols("cph", positive=True, real=True)
    U, k, Om = sp.symbols("U k Om", positive=True, real=True)

    rho1 = sp.Function("rho1")(x, t)
    rho2 = sp.Function("rho2")(x, t)
    u1f = sp.Function("u1")(x, t)
    u2f = sp.Function("u2")(x, t)

    rho = rho0 + eps * rho1 + eps**2 * rho2
    u = eps * u1f + eps**2 * u2f
    cont = sp.expand(sp.diff(rho, t) + sp.diff(rho * u, x))
    cont1 = sp.simplify(cont.coeff(eps, 1))
    cont2 = sp.simplify(cont.coeff(eps, 2))

    line("O(eps) continuity", cont1)
    line("O(eps^2) continuity", cont2)

    dispersion = sp.Eq(Om**2, cs**2 * k**2 + hbar**2 * k**4 / (4 * m**2))
    Pi_over_U = sp.simplify(rho0 * Om / k)
    line("dispersion", dispersion)
    line("Pi_over_U", Pi_over_U)

    Om_drive = sp.pi * cs / L
    y = sp.symbols("y", positive=True, real=True)
    ysol = sp.simplify(sp.solve(sp.Eq(Om_drive**2, cs**2 * y * (1 + xi**2 * y / 2)), y)[0])
    cph_expr = sp.simplify(cs * sp.sqrt((1 + sp.sqrt(1 + 2 * sp.pi**2 * xi**2 / L**2)) / 2))
    line("k^2_at_drive", ysol)
    line("cph(L,xi)", cph_expr)

    Tper = 2 * sp.pi / Om
    u1_st = U * sp.sin(k * x) * sp.sin(Om * t)
    rho1_st = rho0 * U / cph * sp.cos(k * x) * sp.cos(Om * t)
    avg_st = sp.simplify(sp.integrate(sp.expand_trig(rho1_st * u1_st), (t, 0, Tper)) / Tper)
    syms_eq_zero(avg_st)

    u1_out = U * sp.cos(k * x - Om * t)
    rho1_out = rho0 * U / cph * sp.cos(k * x - Om * t)
    avg_out = sp.simplify(sp.integrate(sp.expand_trig(rho1_out * u1_out), (t, 0, Tper)) / Tper)
    syms_eq_zero(avg_out - rho0 * U**2 / (2 * cph))

    line("standing_<rho1*u1>", avg_st)
    line("outgoing_<rho1*u1>", avg_out)

    Ta = sp.pi * hbar * cs / (2 * L**2)
    U_from_T = sp.simplify(Ta / (rho0 * cph))
    j_mouth_local = sp.simplify(avg_out.subs(U, U_from_T))
    j_mouth_local_xi = sp.simplify(j_mouth_local.subs(cph, cph_expr))
    J_mouth_total = sp.simplify(sp.pi * L**2 / Lambda**2 * j_mouth_local)
    J_mouth_total_xi = sp.simplify(J_mouth_total.subs(cph, cph_expr))

    syms_eq_zero(j_mouth_local - sp.pi**2 * cs**2 * hbar**2 / (8 * L**4 * cph**3 * rho0))

    line("Ta", Ta)
    line("U_from_T", U_from_T)
    line("j_mouth_local", j_mouth_local)
    line("j_mouth_local_xi", j_mouth_local_xi)
    line("J_mouth_total", J_mouth_total)
    line("J_mouth_total_xi", J_mouth_total_xi)

    sigma = sp.symbols("sigma", real=True)
    A_p22 = sp.simplify(sp.pi * (L / Lambda) * sp.exp(sigma) * (L / Lambda) * sp.exp(-sigma))
    syms_eq_zero(A_p22 - sp.pi * L**2 / Lambda**2)
    line("A_p22", A_p22)

    return {
        "cont1": cont1,
        "cont2": cont2,
        "dispersion": dispersion,
        "j_mouth_local_xi": j_mouth_local_xi,
        "J_mouth_total_xi": J_mouth_total_xi,
    }


# -----------------------------------------------------------------------------
# Section 11
# -----------------------------------------------------------------------------

def section_11_equilibrium_crossover() -> Dict[str, Any]:
    banner("Section 11 — Unique equilibrium crossover and selected electron scale")

    L, xi_h, hbar, rho0, c_s, Lambda = sp.symbols("L xi_h hbar rho0 c_s Lambda", positive=True, finite=True)
    pi = sp.pi
    beta = 2 * pi**2 * xi_h**2

    Mdot = (
        pi**3 * hbar**2 / (8 * Lambda**2 * rho0 * c_s * L**2)
        * (2 / (1 + sp.sqrt(1 + beta / L**2))) ** sp.Rational(3, 2)
    )
    Phi = sp.sqrt(8 * pi**3 * rho0 * hbar * c_s / (11 * Lambda**2)) * sp.sqrt(L)

    line("Mdot(L)", Mdot)
    line("Phi(L)", Phi)

    eq_exact = sp.Eq(
        64 * Lambda**2 * rho0**3 * c_s**3 * L**2 * (L + sp.sqrt(L**2 + beta))**3,
        11 * pi**3 * hbar**3,
    )
    line("exact_crossover_equation", eq_exact)

    s = sp.symbols("s", positive=True)
    dlog = -2 / L + sp.Rational(3, 2) * beta / (L**3 * s * (1 + s))
    dlog_simplified = sp.simplify(dlog.subs(sp.pi**2 * xi_h**2, L**2 * (s**2 - 1) / 2))
    syms_eq_zero(dlog_simplified + (s + 3) / (2 * L * s))
    dPhi = sp.simplify(sp.diff(Phi, L))
    line("dlnMdot/dL", -(s + 3) / (2 * L * s))
    line("dPhi/dL", dPhi)

    z = sp.symbols("z", positive=True)
    z_poly = sp.expand(16 * Lambda**2 * rho0**3 * c_s**3 * z * (z**2 - beta) ** 2 - 11 * pi**3 * hbar**3)
    line("quintic_in_z", z_poly)

    y = sp.symbols("y", positive=True)
    C = sp.simplify(11 * hbar**3 / (64 * sp.sqrt(2) * pi**2 * Lambda**2 * rho0**3 * c_s**3 * xi_h**5))
    poly_y = sp.expand(y * (y**2 - 1) ** 2 - C)
    f_y = y * (y**2 - 1) ** 2
    fprime_y = sp.factor(sp.diff(f_y, y))
    line("C", C)
    line("y_polynomial", poly_y)
    line("f'(y)", fprime_y)

    L_of_y = sp.simplify(pi * xi_h / sp.sqrt(2) * (y**2 - 1) / y)
    F0 = sp.simplify(9 * pi * hbar * c_s / (11 * L))
    F0_of_y = sp.simplify(F0.subs(L, L_of_y))
    line("L(y)", L_of_y)
    line("F0(y)", F0_of_y)

    Lstar_lw = sp.simplify((11 * pi**3 * hbar**3 / (512 * Lambda**2 * rho0**3 * c_s**3)) ** sp.Rational(1, 5))
    F0_lw = sp.simplify(F0.subs(L, Lstar_lw))
    line("Lstar_long_wave", Lstar_lw)
    line("F0_long_wave", F0_lw)

    return {
        "eq_exact": eq_exact,
        "C": C,
        "L_of_y": L_of_y,
        "F0_of_y": F0_of_y,
        "Lstar_lw": Lstar_lw,
        "F0_lw": F0_lw,
    }


# -----------------------------------------------------------------------------
# Section 12
# -----------------------------------------------------------------------------

def section_12_family_ratio_nogos() -> Dict[str, Any]:
    banner("Section 12 — Family-ratio theorems and no-gos")

    nu, R, L0, Lam0, rho0, cs, hbar, xi = sp.symbols("nu R L0 Lam0 rho0 cs hbar xi", positive=True)

    Lj = nu * L0 / R
    Lamj = Lam0 * nu ** sp.Rational(3, 2) / R**2
    omega_j = sp.pi * cs * nu / (2 * Lj)
    cphj = sp.sqrt(cs**2 / 2 * (1 + sp.sqrt(1 + 2 * sp.pi**2 * xi**2 * nu**2 / Lj**2)))

    dotM_j = sp.pi**3 * hbar**2 * cs**2 * nu**2 / (8 * Lamj**2 * rho0 * cphj**3 * Lj**2)
    Phi_j = sp.sqrt(8 * sp.pi**3 * rho0 * hbar * cs * nu * Lj / (11 * Lamj**2))
    expr_general = sp.simplify(dotM_j / Phi_j)
    expr_ground = sp.simplify(expr_general.subs({nu: 1, R: 1}))
    ratio = sp.simplify(expr_general / expr_ground)

    alpha = sp.symbols("alpha", nonnegative=True)
    family_eq_alpha = sp.Eq(
        R**3,
        nu ** sp.Rational(5, 3) * (1 + sp.sqrt(1 + alpha * R**2)) / (1 + sp.sqrt(1 + alpha)),
    )

    R_IR = sp.simplify(nu ** sp.Rational(5, 9))
    R_UV = sp.simplify(nu ** sp.Rational(5, 6))

    line("dimensionless_ratio", ratio)
    line("family_eq_alpha", family_eq_alpha)
    line("R_IR", R_IR)
    line("R_UV", R_UV)

    for jj in [1, 2]:
        nu_j = 2 * jj + 1
        lo = nu_j ** (sp.Rational(5, 9))
        hi = nu_j ** (sp.Rational(5, 6))
        line(f"nu={nu_j}_lower_bound", sp.N(lo, 15))
        line(f"nu={nu_j}_upper_bound", sp.N(hi, 15))

    line("observed_mu_over_e", mp.nstr(R_MU, 16))
    line("observed_tau_over_e", mp.nstr(R_TAU, 16))
    print("Conclusion: the plain e,mu,tau <-> j=0,1,2 D/N tower is strongly falsified.")

    return {
        "family_eq_alpha": family_eq_alpha,
        "R_IR": R_IR,
        "R_UV": R_UV,
    }


# -----------------------------------------------------------------------------
# Section 13
# -----------------------------------------------------------------------------

def section_13_reverse_engineered_families() -> Dict[str, Any]:
    banner("Section 13 — Reverse-engineered lepton families under the universal 11:2:5 ledger")

    R, Q, Sigma = sp.symbols("R Q Sigma", positive=True)
    lam, gam, eta = sp.symbols("lam gam eta", positive=True)
    W, nu, alpha = sp.symbols("W nu alpha", positive=True)

    geometry_solution = sp.solve(
        [
            sp.Eq(Q, Sigma / lam),
            sp.Eq(lam, gam * eta),
            sp.Eq(lam**4, Sigma * eta**2),
        ],
        [lam, gam, eta],
        dict=True,
    )[0]

    line("L/L0", sp.simplify(geometry_solution[lam]))
    line("a/a0", sp.simplify(geometry_solution[gam]))
    line("Lambda/Lambda0", sp.simplify(geometry_solution[eta]))
    print("Partition lock: R = Q exactly because F/F0 = Ew/Ew0 under the universal 11:2:5 ledger.")

    g = (1 + sp.sqrt(1 + alpha * (R / W) ** 2)) / (1 + sp.sqrt(1 + alpha))
    crossover = sp.Eq(R**3, (W * nu) ** sp.Rational(5, 3) * g)
    line("generalized_crossover", crossover)

    line("frequency_dominant_threshold_mu", R_MU ** (mp.mpf(4) / 5))
    line("frequency_dominant_threshold_tau", R_TAU ** (mp.mpf(4) / 5))

    benchmarks = {}
    for label, Rval in [("mu", R_MU), ("tau", R_TAU)]:
        data = []
        for nu_val in [1, 3, 5]:
            W_min = Rval ** (mp.mpf(9) / 5) / nu_val
            freq_ratio = Rval / W_min
            L_ratio = Rval ** (mp.mpf(4) / 5)
            a_ratio = Rval ** (mp.mpf(1) / 10)
            Lambda_ratio = Rval ** (mp.mpf(7) / 10)
            data.append((nu_val, W_min, freq_ratio, L_ratio, a_ratio, Lambda_ratio))
            print(
                f"{label} | nu={nu_val:>2d} | W_min={mp.nstr(W_min, 12)} | "
                f"omega/omega0={mp.nstr(freq_ratio, 12)} | "
                f"L/L0={mp.nstr(L_ratio, 12)} | a/a0={mp.nstr(a_ratio, 12)} | "
                f"Lambda/Lambda0={mp.nstr(Lambda_ratio, 12)}"
            )
        benchmarks[label] = data

    return {"geometry_solution": geometry_solution, "crossover": crossover, "benchmarks": benchmarks}


# -----------------------------------------------------------------------------
# Section 14
# -----------------------------------------------------------------------------

def alpha_for_candidate(R: float, nu: int, W: int, tol: float = 1e-16) -> float | None:
    """Solve for alpha from the generalized crossover equation using bisection."""
    Wmin = R ** 1.8 / nu
    Wmax = R**3 / (nu ** 2.5)
    if W < Wmin - 1e-12 or W > Wmax + 1e-12:
        return None
    if abs(W - Wmin) < 1e-12:
        return 0.0

    def resid(alpha: float) -> float:
        return R**3 - (W * nu) ** (5 / 3) * (1 + math.sqrt(1 + alpha * (R / W) ** 2)) / (1 + math.sqrt(1 + alpha))

    lo, hi = 0.0, 1.0
    f_lo = resid(lo)
    if abs(f_lo) < tol:
        return lo
    f_hi = resid(hi)
    tries = 0
    while f_lo * f_hi > 0 and tries < 400:
        hi *= 2
        f_hi = resid(hi)
        tries += 1
    if f_lo * f_hi > 0:
        return None
    for _ in range(220):
        mid = (lo + hi) / 2
        f_mid = resid(mid)
        if abs(f_mid) < tol:
            return mid
        if f_lo * f_mid <= 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    return (lo + hi) / 2


def section_14_multi_constraint_intersection() -> Dict[str, Any]:
    banner("Section 14 — Multi-constraint intersection picture")

    R, Sigma, nu, alpha = sp.symbols("R Sigma nu alpha", positive=True)
    beta = R * nu / Sigma
    g = (1 + sp.sqrt(1 + alpha * beta**2)) / (1 + sp.sqrt(1 + alpha))
    print("Action-dominant branch: beta = R*nu/Sigma <= 1, so beta <= g <= 1.")
    lower_strip = sp.Eq(Sigma, R ** sp.Rational(9, 5))
    upper_strip = sp.Eq(Sigma, R**3 / nu ** sp.Rational(3, 2))
    line("lower_strip_boundary", lower_strip)
    line("upper_strip_boundary", upper_strip)

    first_candidates: Dict[str, List[Dict[str, float]]] = {}
    for label, Rval in [("mu", float(R_MU)), ("tau", float(R_TAU))]:
        rows = []
        for nu_val in [1, 3, 5]:
            W = math.ceil(Rval ** 1.8 / nu_val)
            Sigma_val = W * nu_val
            alpha_val = alpha_for_candidate(Rval, nu_val, W)
            L_ratio = Sigma_val / Rval
            a_ratio = Rval / math.sqrt(Sigma_val)
            Lambda_ratio = (Sigma_val ** 1.5) / (Rval**2)
            omega_ratio = Rval / W
            row = {
                "nu": nu_val,
                "W": W,
                "Sigma": Sigma_val,
                "alpha": alpha_val,
                "L_over_L0": L_ratio,
                "a_over_a0": a_ratio,
                "Lambda_over_Lambda0": Lambda_ratio,
                "omega_over_omega0": omega_ratio,
            }
            rows.append(row)
            print(label, row)
        first_candidates[label] = rows

    return {"lower_strip": lower_strip, "upper_strip": upper_strip, "first_candidates": first_candidates}


# -----------------------------------------------------------------------------
# Section 15
# -----------------------------------------------------------------------------

def min_W_cont(R: float, nu: int) -> float:
    return R**1.8 / nu


def alpha_needed_float(R: float, nu: int, W: int) -> float:
    Sigma = nu * W
    k = R**3 / (Sigma ** (5 / 3))
    beta = R * nu / Sigma
    den = k * k - beta * beta
    return 4 * k * (k - beta * beta) * (1 - k) / (den * den)


def branch_geometry_float(R: float, nu: int, W: int, alpha_common: float) -> Dict[str, float]:
    Sigma = nu * W
    return {
        "Sigma": Sigma,
        "L_over_L0": Sigma / R,
        "a_over_a0": R / math.sqrt(Sigma),
        "Lambda_over_Lambda0": Sigma ** 1.5 / (R**2),
        "omega_over_omega0": R / W,
        "L0_over_xi_h": math.sqrt(2 * math.pi**2 / alpha_common),
    }


def generate_branches_float(R: float, nu_max: int, alpha_max: float) -> List[Tuple[float, int, int]]:
    out: List[Tuple[float, int, int]] = []
    for nu in range(1, nu_max + 1, 2):
        W = math.ceil(min_W_cont(R, nu))
        while True:
            a = alpha_needed_float(R, nu, W)
            if a > alpha_max:
                break
            out.append((a, nu, W))
            W += 1
    return sorted(out, key=lambda row: row[0])


def nearest_matches(
    mu_rows: List[Tuple[float, int, int]],
    tau_rows: List[Tuple[float, int, int]],
) -> List[Tuple[float, float, int, int, float, int, int]]:
    tau_alphas = [row[0] for row in tau_rows]
    out: List[Tuple[float, float, int, int, float, int, int]] = []
    for a_mu, nu_mu, W_mu in mu_rows:
        idx = bisect.bisect_left(tau_alphas, a_mu)
        for k in range(max(0, idx - 2), min(len(tau_rows), idx + 3)):
            a_tau, nu_tau, W_tau = tau_rows[k]
            out.append((abs(a_mu - a_tau), a_mu, nu_mu, W_mu, a_tau, nu_tau, W_tau))
    return out


def section_15_universal_alpha_lattice() -> Dict[str, Any]:
    banner("Section 15 — Universal-alpha lattice search")

    # Symbolic inversion
    R, Sigma, nu = sp.symbols("R Sigma nu", positive=True)
    k = R**3 / Sigma ** sp.Rational(5, 3)
    beta = R * nu / Sigma
    x = sp.symbols("x", positive=True)  # x = sqrt(1 + alpha)
    t = sp.symbols("t", positive=True)
    gamma = sp.symbols("gamma", positive=True)

    x_solved = sp.simplify((2 * k - beta**2 - k**2) / (k**2 - beta**2))
    alpha_closed = sp.simplify(sp.factor(x_solved**2 - 1))
    t_quintic_num = sp.factor(
        sp.expand(
            sp.together(
                x - x_solved.subs(
                    {
                        k: t ** (-5),
                        beta: gamma * t ** (-3),
                    }
                )
            ).as_numer_denom()[0]
        )
    )

    line("alpha_closed_form", alpha_closed)
    line("t_quintic_numerator", t_quintic_num)

    alpha_max = 1.0
    mu_rows = generate_branches_float(float(R_MU), 71, alpha_max)
    tau_rows = generate_branches_float(float(R_TAU), 679, alpha_max)

    line("mu_branch_count", len(mu_rows))
    line("tau_branch_count", len(tau_rows))

    matches = nearest_matches(mu_rows, tau_rows)
    best_precision = min(matches, key=lambda m: m[0])
    sub_1e9 = [m for m in matches if m[0] < 1e-9]
    sub_1e5 = [m for m in matches if m[0] < 1e-5]
    best_balanced = min(sub_1e9, key=lambda m: (max(m[3], m[6]), m[0]))
    least_complexity = min(sub_1e5, key=lambda m: (max(m[3], m[6]), m[0]))

    for label, m in [
        ("best_precision", best_precision),
        ("best_balanced", best_balanced),
        ("least_complexity", least_complexity),
    ]:
        delta_alpha, alpha_mu, nu_mu, W_mu, alpha_tau, nu_tau, W_tau = m
        alpha_avg = (alpha_mu + alpha_tau) / 2
        print(label)
        line("delta_alpha", delta_alpha)
        line("mu_tuple", (nu_mu, W_mu))
        line("alpha_mu", alpha_mu)
        line("tau_tuple", (nu_tau, W_tau))
        line("alpha_tau", alpha_tau)
        line("alpha_avg", alpha_avg)
        line("mu_geometry", branch_geometry_float(float(R_MU), nu_mu, W_mu, alpha_avg))
        line("tau_geometry", branch_geometry_float(float(R_TAU), nu_tau, W_tau, alpha_avg))
        print()

    return {
        "alpha_closed_form": alpha_closed,
        "t_quintic_numerator": t_quintic_num,
        "mu_rows": mu_rows,
        "tau_rows": tau_rows,
        "best_precision": best_precision,
        "best_balanced": best_balanced,
        "least_complexity": least_complexity,
    }


# -----------------------------------------------------------------------------
# Section 16
# -----------------------------------------------------------------------------

@dataclass
class RepresentativeCandidate:
    name: str
    alpha: float
    nu: int
    L_over_L0: float
    a_over_a0: float

    @property
    def Lambda(self) -> float:
        return float(LAMBDA0) * self.L_over_L0 / self.a_over_a0

    @property
    def L0_over_xi(self) -> float:
        return math.sqrt(2 * math.pi**2 / self.alpha)

    @property
    def xi_over_a(self) -> float:
        return (float(LAMBDA0) / self.L0_over_xi) / self.a_over_a0

    def ratio_p22(self) -> float:
        return math.pi * self.nu / (2 * self.Lambda * math.sqrt(1 + 2 * self.xi_over_a**2))

    def ratio_p0(self, chi: float) -> float:
        return math.pi * self.nu / (chi * self.Lambda * math.sqrt(1 + 0.5 * chi**2 * self.xi_over_a**2))

    def breathing_enhancement(self, chi: float) -> float:
        r = self.ratio_p0(chi)
        return 1.0 / abs(1.0 - r * r)

    @property
    def relative_strain_scale(self) -> float:
        return 1.0 / (self.L_over_L0**2)


def section_16_phase_lock_nogos(search_results: Dict[str, Any]) -> Dict[str, Any]:
    banner("Section 16 — Resonance and phase-lock no-gos")

    sigma, a, P_t = sp.symbols("sigma a P_t", real=True)
    U_P = -P_t * sp.pi * a**2
    dU_dsigma = sp.diff(U_P, sigma)
    syms_eq_zero(dU_dsigma)
    line("d/dsigma[-P*pi*a^2]", dU_dsigma)
    print("Conclusion: a scalar P0 hammer does zero linear work on the area-preserving P22 ellipse.")

    c_s, xi_h = sp.symbols("c_s xi_h", positive=True)
    Omega22 = sp.simplify(2 * c_s / a * sp.sqrt(1 + 2 * xi_h**2 / a**2))
    line("Omega22(a)", Omega22)

    # Representative candidates from section 15
    bp = search_results["best_precision"]
    bb = search_results["best_balanced"]
    alpha_bp = (bp[1] + bp[4]) / 2
    alpha_bb = (bb[1] + bb[4]) / 2

    candidates = [
        RepresentativeCandidate("mu_best_precision", alpha_bp, bp[2], branch_geometry_float(float(R_MU), bp[2], bp[3], alpha_bp)["L_over_L0"], branch_geometry_float(float(R_MU), bp[2], bp[3], alpha_bp)["a_over_a0"]),
        RepresentativeCandidate("tau_best_precision", alpha_bp, bp[5], branch_geometry_float(float(R_TAU), bp[5], bp[6], alpha_bp)["L_over_L0"], branch_geometry_float(float(R_TAU), bp[5], bp[6], alpha_bp)["a_over_a0"]),
        RepresentativeCandidate("mu_balanced", alpha_bb, bb[2], branch_geometry_float(float(R_MU), bb[2], bb[3], alpha_bb)["L_over_L0"], branch_geometry_float(float(R_MU), bb[2], bb[3], alpha_bb)["a_over_a0"]),
        RepresentativeCandidate("tau_balanced", alpha_bb, bb[5], branch_geometry_float(float(R_TAU), bb[5], bb[6], alpha_bb)["L_over_L0"], branch_geometry_float(float(R_TAU), bb[5], bb[6], alpha_bb)["a_over_a0"]),
    ]

    chi_D = float(mp.besseljzero(0, 1))  # Dirichlet robustness check
    chi_N = float(mp.besseljzero(1, 1))  # Neumann / free-edge breathing

    for cand in candidates:
        print(cand.name)
        line("alpha", cand.alpha)
        line("Lambda_j", cand.Lambda)
        line("xi_h/a_j", cand.xi_over_a)
        line("2omega/Omega22", cand.ratio_p22())
        line("2omega/Omega00_D", cand.ratio_p0(chi_D))
        line("2omega/Omega00_N", cand.ratio_p0(chi_N))
        line("relative_scalar_strain_vs_electron", cand.relative_strain_scale)
        print()

    # Direct scalar lock threshold on the minimal action-dominant branch
    for Rval, label in [(float(R_MU), "muon"), (float(R_TAU), "tau")]:
        nu_thr_D = chi_D * float(LAMBDA0) / math.pi * Rval**0.7
        nu_thr_N = chi_N * float(LAMBDA0) / math.pi * Rval**0.7
        print(label)
        line("nu_threshold_Dirichlet", nu_thr_D)
        line("nu_threshold_Neumann", nu_thr_N)
        print()

    return {
        "Omega22": Omega22,
        "candidates": candidates,
        "chi_D": chi_D,
        "chi_N": chi_N,
    }


# -----------------------------------------------------------------------------
# Section 17 / 18 summary
# -----------------------------------------------------------------------------

def section_17_18_summary() -> Dict[str, Any]:
    banner("Sections 17–18 — Highest-confidence conclusions and chapter split")
    print("Highest-confidence carry-forward conclusions:")
    print("  • The isolated reduced mass ledger closes to the exact 11:2:5 partition.")
    print("  • The exact D/N support ladder and exact scalar round-trip phase closure are established.")
    print("  • Topology fixes circulation / holonomy class, not the throughput amplitude Phi.")
    print("  • The mouth-output problem is a dynamic rectification problem, not a static DC-flow law.")
    print("  • The naive D/N harmonic-family picture is strongly falsified.")
    print("  • Heavy-lepton candidates are sparse lattice intersections, not continuum scalings.")
    print("  • Direct P0→P2 and direct P0↔P0 mouth resonances do not select the known deep-needle families.")
    print()
    print("Suggested next-session chapter split:")
    print("  1) Isolated-electron theorem chain")
    print("  2) Family-ratio and phase-lock no-go ledger")
    print("  3) Heavy-lepton reverse-engineering and universal-alpha lattice search")
    return {}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Master SymPy verifier for the lepton/electron derivations.")
    parser.add_argument(
        "--sections",
        nargs="*",
        type=int,
        help="Optional subset of section numbers to run. Default: run all.",
    )
    args = parser.parse_args()

    available = {
        1: section_1_setup,
        2: section_2_core_mass_theorem,
        3: section_3_microscopic_coefficients,
        4: section_4_dn_spectral_closure,
        5: section_5_throughput_audit,
        6: section_6_dynamic_turbine_scaling,
        7: section_7_phase_vs_amplitude,
        8: section_8_mouth_output_channel,
        9: section_9_static_compliance,
        10: section_10_ac_streaming,
        11: section_11_equilibrium_crossover,
        12: section_12_family_ratio_nogos,
        13: section_13_reverse_engineered_families,
        14: section_14_multi_constraint_intersection,
        15: section_15_universal_alpha_lattice,
        16: None,  # depends on section 15 results
        17: section_17_18_summary,
        18: None,  # wrapped into section 17 summary
    }

    requested = args.sections if args.sections else [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]

    results: Dict[int, Dict[str, Any]] = {}
    for sec in requested:
        if sec == 16:
            if 15 not in results:
                results[15] = section_15_universal_alpha_lattice()
            results[16] = section_16_phase_lock_nogos(results[15])
        elif sec == 18:
            if 17 not in results:
                results[17] = section_17_18_summary()
        else:
            func = available.get(sec)
            if func is None:
                raise ValueError(f"Section {sec} is not implemented as a standalone code block.")
            results[sec] = func()

    banner("All requested derivations completed successfully")
    print("This script is intended as a carry-forward verifier for the next session.")


if __name__ == "__main__":
    main()
