
#!/usr/bin/env python3
"""Standalone SymPy driver for the full reduced-sector 4D write-up.

This script consolidates the symbolic derivations and benchmark numerics that
appear across Parts I–VIII and the appendices of the write-up. It is designed
to be runnable as a single file with SymPy installed.

Coverage:
- Hydrogenic reduction and Bohr-scale minimization from the reduced 4D action
- Two-body reduced mass upgrade
- Finite-size atomic response, constant-area P22 ellipse, and core regulator
- Lepton doublet machinery: intrinsic P22 bracing and isolated ellipticity
- Reduced Dirac g=2 bridge and common-carrier inertia matching
- Electron anomaly hierarchy:
  support transfer, eta1 = 11/36, leakage softening, local PDE blur,
  charge-side transport closure, and final benchmark residual

The script prints both symbolic identities and the numerical values used in the
paper's benchmark tables.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Any, OrderedDict
import sympy as sp


# ---------------------------------------------------------------------------
# Global helpers
# ---------------------------------------------------------------------------

sp.init_printing()


def banner(title: str) -> None:
    print()
    print("=" * 88)
    print(title)
    print("=" * 88)


def subbanner(title: str) -> None:
    print()
    print("-" * 88)
    print(title)
    print("-" * 88)


def p(label: str, expr: Any) -> None:
    print(f"{label} = {expr}")


def pn(label: str, expr: Any, digits: int = 20) -> None:
    print(f"{label} = {sp.N(expr, digits)}")


def to_float(expr: Any, digits: int = 30) -> float:
    return float(sp.N(expr, digits))


@dataclass(frozen=True)
class Benchmarks:
    # Benchmark values reused in the write-up
    alpha_fs: sp.Float = sp.Float("0.0072973525643", 40)
    g_e_mag: sp.Float = sp.Float("2.00231930436092", 40)

    # Constants for hydrogen numerics
    e_charge: sp.Float = sp.Float("1.602176634e-19", 40)
    epsilon0: sp.Float = sp.Float("8.8541878128e-12", 40)
    hbar: sp.Float = sp.Float("1.054571817e-34", 40)
    m_e: sp.Float = sp.Float("9.1093837139e-31", 40)
    m_p: sp.Float = sp.Float("1.67262192595e-27", 40)

    # Frozen reduced closure input
    L_over_a: sp.Float = sp.Float("1.85", 40)
    K22: sp.Rational = sp.Rational(8, 3)


BENCH = Benchmarks()




# ---------------------------------------------------------------------------
# Part I: Foundations
# ---------------------------------------------------------------------------

def foundations_snapshot() -> Dict[str, Any]:
    banner("PART I — FOUNDATIONS SNAPSHOT")

    eta_Q = sp.symbols("eta_Q", real=True)
    e_star, q_star, Z_int, q_eff, e_eff = sp.symbols(
        "e_star q_star Z_int q_eff e_eff", positive=True, real=True
    )
    rho, K, mu0_eff, qsrc, r, lam = sp.symbols(
        "rho K mu0_eff q_src r lambda", positive=True, real=True
    )

    q_star_law = sp.Eq(q_star, eta_Q * e_star)
    q_eff_law = sp.Eq(q_eff, q_star / sp.sqrt(Z_int))
    e_eff_law = sp.Eq(e_eff, e_star / sp.sqrt(Z_int))
    h_rho = sp.simplify(sp.Rational(5, 4) * K * rho**4)
    A0_static = sp.Eq(
        sp.Function("A_0")(r),
        mu0_eff * qsrc / (4 * sp.pi * r) * (1 + sp.Rational(1, 2) * sp.exp(-2 * r / lam)),
    )

    subbanner("Core reduction identities used everywhere downstream")
    p("q_star = eta_Q e_star", q_star_law)
    p("q_eff", q_eff_law)
    p("e_eff", e_eff_law)
    p("h(rho) = dU/drho for U = K rho^5 / 4", h_rho)
    p("Static zero-mode + leading KK potential", A0_static)

    return {
        "q_star_law": q_star_law,
        "q_eff_law": q_eff_law,
        "e_eff_law": e_eff_law,
        "h_rho": h_rho,
        "A0_static": A0_static,
    }


# ---------------------------------------------------------------------------
# Part II: Hydrogen from the existing 4D action
# ---------------------------------------------------------------------------

def hydrogenic_reduction() -> Dict[str, Any]:
    banner("PART II — HYDROGEN FROM THE EXISTING 4D ACTION")

    r, a, lam = sp.symbols("r a lam", positive=True, real=True)
    hbar, m, mu, g_C, E_perp, K, Gamma10 = sp.symbols(
        "hbar m mu g_C E_perp K Gamma10", positive=True, real=True
    )
    eps0, e_eff, e_star, Zint = sp.symbols(
        "epsilon_0 e_eff e_star Z_int", positive=True, real=True
    )
    lam_w = sp.symbols("lambda", positive=True, real=True)

    phi = sp.exp(-r / a) / sp.sqrt(sp.pi * a**3)

    norm = sp.simplify(4 * sp.pi * sp.integrate(r**2 * phi**2, (r, 0, sp.oo)))
    kinetic = sp.simplify(
        (hbar**2 / (2 * m)) * 4 * sp.pi * sp.integrate(r**2 * sp.diff(phi, r) ** 2, (r, 0, sp.oo))
    )
    avg_inv_r = sp.simplify(4 * sp.pi * sp.integrate(r**2 * (1 / r) * phi**2, (r, 0, sp.oo)))
    avg_yuk = sp.simplify(
        4 * sp.pi * sp.integrate(r**2 * sp.exp(-2 * r / lam) * phi**2 / r, (r, 0, sp.oo))
    )
    avg_phi10 = sp.simplify(4 * sp.pi * sp.integrate(r**2 * phi**10, (r, 0, sp.oo)))

    E_a = sp.simplify(
        E_perp
        + kinetic
        - g_C * avg_inv_r
        - g_C * sp.Rational(1, 2) * avg_yuk
        + K * Gamma10 * avg_phi10 / 4
    )
    E_clean = sp.simplify(E_perp + hbar**2 / (2 * m * a**2) - g_C / a)
    a_star = sp.solve(sp.Eq(sp.diff(E_clean, a), 0), a)[0]

    # Reduced-mass derivation
    m_e_sym, m_p_sym = sp.symbols("m_e m_p", positive=True, real=True)
    M = sp.symbols("M", positive=True, real=True)
    Xed, Xpd, Rd, rd = sp.symbols("Xed Xpd Rd rd", real=True)
    M_expr = sp.simplify(m_e_sym + m_p_sym)
    Rdot, rdot = sp.symbols("Rdot rdot", real=True)
    xe_dot = Rdot + m_p_sym / M_expr * rdot
    xp_dot = Rdot - m_e_sym / M_expr * rdot
    T_two = sp.expand(sp.Rational(1, 2) * m_e_sym * xe_dot**2 + sp.Rational(1, 2) * m_p_sym * xp_dot**2)
    mu_expr = sp.simplify(m_e_sym * m_p_sym / M_expr)
    T_split = sp.expand(sp.Rational(1, 2) * M_expr * Rdot**2 + sp.Rational(1, 2) * mu_expr * rdot**2)

    # Standard Coulomb matching and thickness law
    a_star_em = sp.simplify(a_star.subs(g_C, e_eff**2 / (4 * sp.pi * eps0)))
    a_star_thickness = sp.simplify(a_star_em.subs(e_eff, e_star / sp.sqrt(Zint)))
    a_star_gaussian = sp.simplify(a_star_thickness.subs(Zint, sp.sqrt(sp.pi) * lam_w))
    a_star_two_body = sp.simplify((hbar**2 / (mu * g_C)))

    # Benchmark numerics
    a0_fixed = sp.simplify(
        4 * sp.pi * BENCH.epsilon0 * BENCH.hbar**2 / (BENCH.m_e * BENCH.e_charge**2)
    )
    mu_num = sp.simplify(BENCH.m_e * BENCH.m_p / (BENCH.m_e + BENCH.m_p))
    a0_reduced = sp.simplify(
        4 * sp.pi * BENCH.epsilon0 * BENCH.hbar**2 / (mu_num * BENCH.e_charge**2)
    )
    eV = sp.Float("1.602176634e-19", 40)
    E1_fixed_J = sp.simplify(-BENCH.m_e * BENCH.e_charge**4 / (2 * (4 * sp.pi * BENCH.epsilon0) ** 2 * BENCH.hbar**2))
    E1_reduced_J = sp.simplify(-mu_num * BENCH.e_charge**4 / (2 * (4 * sp.pi * BENCH.epsilon0) ** 2 * BENCH.hbar**2))
    E1_fixed_eV = sp.simplify(E1_fixed_J / eV)
    E1_reduced_eV = sp.simplify(E1_reduced_J / eV)

    subbanner("Hydrogenic expectation values for the 1s trial state")
    p("Normalization", norm)
    p("<T>", kinetic)
    p("<1/r>", avg_inv_r)
    p("<e^{-2r/lambda}/r>", avg_yuk)
    p("∫|phi_a|^10 d^3x", avg_phi10)

    subbanner("Reduced energy and Bohr-scale minimization")
    p("E(a)", E_a)
    p("Clean limit E_clean(a)", E_clean)
    p("Stationary radius a_*", a_star)
    p("Coulomb-matched a_*", a_star_em)
    p("Thickness law a_*(Z_int)", a_star_thickness)
    p("Gaussian thickness law a_*(lambda)", a_star_gaussian)

    subbanner("Two-body upgrade")
    p("Two-body kinetic energy T_two", T_two)
    p("Center/relative split T_split", T_split)
    print("T_two - T_split =", sp.simplify(T_two - T_split))
    p("Reduced mass mu", mu_expr)
    p("Two-body hydrogen radius", a_star_two_body)

    subbanner("Hydrogen benchmark numerics")
    pn("Fixed-source Bohr radius a0 [m]", a0_fixed, 25)
    pn("Reduced-mass hydrogen radius a_H [m]", a0_reduced, 25)
    pn("Reduced mass mu [kg]", mu_num, 25)
    pn("Fixed-source ground-state energy [eV]", E1_fixed_eV, 20)
    pn("Reduced-mass ground-state energy [eV]", E1_reduced_eV, 20)

    return {
        "norm": norm,
        "kinetic": kinetic,
        "avg_inv_r": avg_inv_r,
        "avg_yuk": avg_yuk,
        "avg_phi10": avg_phi10,
        "E_a": E_a,
        "a_star": a_star,
        "a_star_em": a_star_em,
        "a_star_thickness": a_star_thickness,
        "a_star_gaussian": a_star_gaussian,
        "mu_expr": mu_expr,
        "a0_fixed": a0_fixed,
        "a0_reduced": a0_reduced,
        "E1_fixed_eV": E1_fixed_eV,
        "E1_reduced_eV": E1_reduced_eV,
    }


# ---------------------------------------------------------------------------
# Part III: Finite-size atomic response and core regulation
# ---------------------------------------------------------------------------

def finite_size_atomic_response() -> Dict[str, Any]:
    banner("PART III — FINITE-SIZE ATOMIC RESPONSE AND CORE REGULATION")

    x, y, r, g = sp.symbols("x y r g", positive=True, real=True)
    phi = sp.symbols("phi", real=True)
    a, sigma = sp.symbols("a sigma", positive=True, real=True)
    u = sp.symbols("u", real=True)

    rho = sp.sqrt(x**2 + y**2)
    Phi_point = -g / rho
    Txx = sp.simplify(sp.diff(Phi_point, x, 2))
    Tyy = sp.simplify(sp.diff(Phi_point, y, 2))
    Txy = sp.simplify(sp.diff(Phi_point, x, y))
    T2_cart = sp.simplify(Txx - Tyy)
    T2_mix = sp.simplify(2 * Txy)

    T2_radial = sp.simplify(sp.diff(-g / r, r, 2) - sp.diff(-g / r, r) / r)

    # Constant-area ellipse
    r_boundary = sp.simplify(a / sp.sqrt(sp.exp(-2 * sigma) * sp.cos(u) ** 2 + sp.exp(2 * sigma) * sp.sin(u) ** 2))
    boundary_series = sp.series(sp.simplify(r_boundary / a), sigma, 0, 3)
    theta = sp.symbols("theta", real=True)
    xell = a * sp.exp(sigma) * sp.cos(theta)
    yell = a * sp.exp(-sigma) * sp.sin(theta)
    area_integrand = sp.simplify(xell * sp.diff(yell, theta) - yell * sp.diff(xell, theta))
    area_loop = sp.simplify(sp.integrate(area_integrand, (theta, 0, 2 * sp.pi)))

    # Quadrupole moments
    phi_m = sp.symbols("phi_m", real=True)
    Mxx_p = a**2 * sp.exp(2 * sigma) / 4
    Myy_p = a**2 * sp.exp(-2 * sigma) / 4
    c, s = sp.cos(phi_m), sp.sin(phi_m)
    Mxx = sp.simplify(Mxx_p * c**2 + Myy_p * s**2)
    Myy = sp.simplify(Mxx_p * s**2 + Myy_p * c**2)
    Mxy = sp.simplify((Mxx_p - Myy_p) * s * c)
    Qc = sp.simplify(Mxx - Myy)
    Qs = sp.simplify(2 * Mxy)

    # Gaussian-smoothed core regulator
    xdim = sp.symbols("xdim", positive=True, real=True)
    Phi_sm = -g * sp.erf(r / (sp.sqrt(2) * a)) / r
    T2_eff = sp.simplify(sp.diff(Phi_sm, r, 2) - sp.diff(Phi_sm, r) / r)
    F_dimless = sp.simplify((a**3 / g) * T2_eff.subs(r, a * xdim))
    F_small = sp.series(F_dimless, xdim, 0, 5)
    F_large = sp.simplify(sp.limit(xdim**3 * F_dimless, xdim, sp.oo))

    subbanner("Point-Coulomb Hessian and quadrupole tide")
    p("T_xx - T_yy", T2_cart)
    p("2 T_xy", T2_mix)
    p("T2(r) = Phi'' - Phi'/r", T2_radial)

    subbanner("Constant-area P22 ellipse")
    p("Boundary r(theta)/a", r_boundary / a)
    p("Boundary series through O(sigma^2)", boundary_series)
    p("Area integrand x dy - y dx", area_integrand)
    p("Loop area integral", area_loop)
    p("Q_c", Qc)
    p("Q_s", Qs)

    subbanner("Finite-throat Gaussian core regulator")
    p("Phi_sm(r)", Phi_sm)
    p("T2_eff(r;a)", T2_eff)
    p("Dimensionless F(x) with T2_eff = (g/a^3) F(x)", F_dimless)
    p("Small-x series of F(x)", F_small)
    p("Large-x limit of x^3 F(x)", F_large)

    return {
        "T2_radial": T2_radial,
        "boundary_series": boundary_series,
        "area_loop": area_loop,
        "Qc": Qc,
        "Qs": Qs,
        "T2_eff": T2_eff,
        "F_dimless": F_dimless,
        "F_small": F_small,
        "F_large": F_large,
    }


# ---------------------------------------------------------------------------
# Part IV: From atomic P22 forcing to the lepton doublet
# ---------------------------------------------------------------------------

def lepton_doublet_bridge() -> Dict[str, Any]:
    banner("PART IV — FROM ATOMIC P22 FORCING TO THE LEPTON DOUBLET")

    r, theta = sp.symbols("r theta", positive=True, real=True)
    bx, by = sp.symbols("b_x b_y", real=True)
    chi, chir = sp.symbols("chi chi_r", real=True)
    nu0, hbar, a, Cmix, LambdaQ, k22 = sp.symbols(
        "nu_0 hbar a C_mix Lambda_Q k_22", positive=True, real=True
    )
    phi_alpha = sp.symbols("phi_alpha", real=True)

    x = r * sp.cos(theta)
    y = r * sp.sin(theta)
    bdotx = bx * x + by * y
    dAx = bx * chi + bdotx * chir * x / r
    dAy = by * chi + bdotx * chir * y / r

    stress_xx = sp.expand(dAx * dAx - sp.Rational(1, 2) * (dAx * dAx + dAy * dAy))
    stress_xy = sp.expand(dAx * dAy)

    I_xx = sp.simplify(sp.integrate(stress_xx, (theta, 0, 2 * sp.pi)))
    I_xy = sp.simplify(sp.integrate(stress_xy, (theta, 0, 2 * sp.pi)))
    L = sp.expand((2 * chi + r * chir) ** 2)
    target_xx = sp.simplify(sp.pi * L * (bx**2 - by**2) / 4)
    target_xy = sp.simplify(sp.pi * L * bx * by / 2)

    s_alpha = sp.simplify(hbar / nu0)
    h_alpha = sp.simplify(LambdaQ * a**2 * Cmix * hbar / (4 * nu0))
    sigma_inf = sp.simplify(2 * h_alpha / k22)

    subbanner("Mixed-sector traceless stress theorem")
    p("Integrated I_xx", I_xx)
    p("Target traceless form I_xx", target_xx)
    print("I_xx - target =", sp.simplify(I_xx - target_xx))
    p("Integrated I_xy", I_xy)
    p("Target traceless form I_xy", target_xy)
    print("I_xy - target =", sp.simplify(I_xy - target_xy))

    subbanner("Intrinsic half-flux lock and isolated P22 bracing")
    p("alpha = nu_0 s / (2 hbar)", sp.Symbol("alpha"))
    p("s_alpha from alpha = 1/2", s_alpha)
    p("h_alpha", h_alpha)
    p("sigma_infinity ≈ 2 h_alpha / k22", sigma_inf)

    return {
        "I_xx": I_xx,
        "I_xy": I_xy,
        "target_xx": target_xx,
        "target_xy": target_xy,
        "s_alpha": s_alpha,
        "h_alpha": h_alpha,
        "sigma_inf": sigma_inf,
    }


# ---------------------------------------------------------------------------
# Part V: Dirac g=2 from geometric gear reduction
# ---------------------------------------------------------------------------

def dirac_g2_bridge() -> Dict[str, Any]:
    banner("PART V — DIRAC g=2 FROM GEOMETRIC GEAR REDUCTION")

    theta = sp.symbols("theta", real=True)
    a, sigma = sp.symbols("a sigma", positive=True, real=True)
    q_eff, T = sp.symbols("q_eff T", positive=True, real=True)
    M_mu, M_part = sp.symbols("M_mu M_part", positive=True, real=True)
    m_mode, alpha = sp.symbols("m_mode alpha", positive=True, real=True)
    C2q, C2M, Omega = sp.symbols("C2q C2M Omega", positive=True, real=True)

    x = a * sp.exp(sigma) * sp.cos(theta)
    y = a * sp.exp(-sigma) * sp.sin(theta)
    loop_integrand = sp.simplify(x * sp.diff(y, theta) - y * sp.diff(x, theta))
    loop_integral = sp.simplify(sp.integrate(loop_integrand, (theta, 0, 2 * sp.pi)))

    mu_z = sp.simplify(q_eff / (2 * T) * loop_integral)
    L_ext = sp.simplify(M_mu / T * loop_integral)
    g_red = sp.simplify((m_mode / alpha) * (M_part / M_mu))

    mu_general = sp.simplify(q_eff * a**2 * Omega * C2q / 2)
    L_general = sp.simplify(M_part * a**2 * Omega * C2M)
    zeta = sp.simplify(C2q / C2M)
    common_carrier_ratio = sp.simplify(mu_general / (q_eff * L_general / (2 * M_part)))
    Mmu_expr = sp.simplify(M_part * C2M / C2q)

    subbanner("Charged ellipse loop identities")
    p("Loop integrand x dy - y dx", loop_integrand)
    p("Loop integral", loop_integral)
    p("mu_z", mu_z)
    p("L_ext", L_ext)
    print("mu_z - q_eff/(2 M_mu) L_ext =", sp.simplify(mu_z - q_eff * L_ext / (2 * M_mu)))

    subbanner("Reduced gyromagnetic ratio")
    p("g_red", g_red)
    p("g_red at m=1, alpha=1/2", sp.simplify(g_red.subs({m_mode: 1, alpha: sp.Rational(1, 2)})))
    p(
        "g_red at m=1, alpha=1/2, M_mu=M_part",
        sp.simplify(g_red.subs({m_mode: 1, alpha: sp.Rational(1, 2), M_mu: M_part})),
    )

    subbanner("Common-carrier inertia matching")
    p("mu from support moment", mu_general)
    p("L from support moment", L_general)
    p("zeta = C2q/C2M", zeta)
    p("mu / [q/(2M) L]", common_carrier_ratio)
    p("M_mu", Mmu_expr)

    return {
        "loop_integral": loop_integral,
        "g_red": g_red,
        "zeta": zeta,
        "Mmu_expr": Mmu_expr,
    }


# ---------------------------------------------------------------------------
# Part VI: Electron anomaly hierarchy
# ---------------------------------------------------------------------------

def electron_anomaly_chain() -> Dict[str, Any]:
    banner("PART VI — THE ELECTRON ANOMALY AS A CONTINUUM SELF-DRESSING PROBLEM")

    f, delta, kappa = sp.symbols("f delta kappa", positive=True, real=True)
    rho, theta = sp.symbols("rho theta", positive=True, real=True)
    A = sp.symbols("A", positive=True, real=True)

    # Experimental target and tree-level mismatch
    zeta_e = sp.simplify(BENCH.g_e_mag / 2)
    f_num = sp.simplify(BENCH.alpha_fs / (2 * sp.pi))
    a_e = sp.simplify(zeta_e - 1)
    delta_zeta_e = sp.simplify(zeta_e - 1)

    # Normalized real 22 mode and eta1 = 11/36
    dmu = rho / sp.pi
    phi22 = A * rho**2 * sp.cos(2 * theta)
    norm22 = sp.simplify(sp.integrate(sp.integrate(phi22**2 * dmu, (theta, 0, 2 * sp.pi)), (rho, 0, 1)))
    A22 = sp.solve(sp.Eq(norm22, 1), A)[0]
    phi22_n = sp.simplify(phi22.subs(A, A22))
    self_loop_factor = sp.simplify(
        sp.Rational(1, 2)
        * sp.integrate(sp.integrate((f * phi22_n) ** 2 * dmu, (theta, 0, 2 * sp.pi)), (rho, 0, 1))
    )
    eta1_min = sp.Rational(11, 18) * sp.Rational(1, 2)

    # Sharp anomaly law
    zeta_sharp = sp.expand((1 + f - f**2) / (1 + eta1_min * f**2)).series(f, 0, 4).removeO()
    g_sharp_exact = sp.simplify(2 * (1 + f - (1 + eta1_min) * f**2))

    # Boundary-leakage overlap factors
    chi_lin = sp.Piecewise((1, rho <= 1 - delta), ((1 - rho) / delta, True))
    B_lin = sp.simplify(
        6 * (
            sp.integrate(rho**5, (rho, 0, 1 - delta))
            + sp.integrate(rho**5 * (1 - rho) / delta, (rho, 1 - delta, 1))
        )
    )

    chi_exp = 1 - sp.exp(-(1 - rho) / delta)
    B_exp = sp.simplify(6 * sp.integrate(rho**5 * chi_exp, (rho, 0, 1)))
    eta1_of_delta = sp.simplify(eta1_min * B_exp)

    eta1_e_target = sp.simplify((1 + f_num - zeta_e) / (f_num**2) - 1)
    B_e = sp.simplify(eta1_e_target / eta1_min)

    delta_lin_target = sp.nsolve(sp.Eq(B_lin.subs(delta, sp.Symbol("x")), B_e), 0.0028)
    delta_exp_target = sp.nsolve(sp.Eq(B_exp.subs(delta, sp.Symbol("x")), B_e), 0.0014)

    # Local self-transport PDE
    x = sp.symbols("x", nonnegative=True, real=True)
    chi = sp.Function("chi")
    v_mix, omega_half = sp.symbols("v_mix omega_half", positive=True, real=True)
    ode = sp.Eq(v_mix * sp.diff(chi(x), x), omega_half * (1 - chi(x)))
    ode_sol = sp.dsolve(ode, ics={chi(0): 0})
    L, a = sp.symbols("L a", positive=True, real=True)
    c_s = sp.symbols("c_s", positive=True, real=True)
    k_half = sp.pi / (2 * L)
    ell_blur = sp.simplify((f * c_s) / (c_s * k_half))
    delta_loc = sp.simplify(ell_blur / a)
    kappa_expr = sp.simplify(2 * L / (sp.pi * a))

    # Local inertia blur and cubic coefficient
    eta1_loc_f = sp.simplify(eta1_min * B_exp.subs(delta, kappa * f))
    g_local_inertia = sp.simplify(2 * (1 + f - (1 + eta1_loc_f) * f**2))
    g_local_series = sp.expand(
        sp.series(
            2
            * (
                1
                + f
                - (
                    1
                    + eta1_min
                    * (1 - 6 * kappa * f + 30 * kappa**2 * f**2 - 120 * kappa**3 * f**3)
                )
                * f**2
            ),
            f,
            0,
            5,
        ).removeO()
    )
    c3_inertia = sp.simplify(sp.Rational(11, 6) * kappa)

    # Charge-side transport closure
    tau = sp.simplify(1 - sp.sqrt(1 - f))
    xx = sp.symbols("xx", nonnegative=True, real=True)
    weight = 2 * (1 - xx)
    Ncol = sp.integrate(weight, (xx, 0, tau))
    cbar = sp.simplify(sp.integrate(weight * sp.cos(sp.pi * xx / (2 * tau)), (xx, 0, tau)) / Ncol)
    phi_q = sp.simplify(sp.cos(sp.pi * xx / (2 * tau)) - cbar)
    Xi = sp.simplify(sp.integrate(weight * (1 - xx) ** 2 * phi_q, (xx, 0, tau)) / Ncol)
    Aq = sp.simplify(tau / kappa)
    Q_loc = sp.simplify(1 + f - f**2 + 2 * f * Aq * Xi)
    Xi_series = sp.series(Xi, f, 0, 4).removeO()
    Q_loc_series = sp.series(Q_loc, f, 0, 5).removeO()
    c3_q = sp.simplify((4 - sp.pi) / (sp.pi**2 * kappa))
    c3_total = sp.simplify(c3_inertia + c3_q)

    g_loc = sp.simplify(2 * (Q_loc - eta1_loc_f * f**2))

    # Benchmark numerics
    kappa_num = sp.simplify(2 * BENCH.L_over_a / sp.pi)
    delta_loc_num = sp.simplify(kappa_num * f_num)
    B_exp_loc = sp.simplify(B_exp.subs(delta, delta_loc_num))
    eta1_loc_num = sp.simplify(eta1_min * B_exp_loc)

    g_tree = sp.Integer(2)
    g_sharp_num = sp.simplify(g_sharp_exact.subs(f, f_num))
    g_local_inertia_num = sp.simplify(g_local_inertia.subs({f: f_num, kappa: kappa_num}))
    g_loc_num = sp.simplify(g_loc.subs({f: f_num, kappa: kappa_num}))
    residual = sp.simplify(BENCH.g_e_mag - g_loc_num)

    # Optional quartic coefficient needed to exactly match the benchmark target
    c4_req = sp.simplify(
        (BENCH.g_e_mag / 2 - (1 + f_num - sp.Rational(47, 36) * f_num**2 + c3_total.subs(kappa, kappa_num) * f_num**3))
        / f_num**4
    )

    subbanner("Experimental target and tree-level mismatch")
    pn("f = alpha/(2 pi)", f_num, 20)
    pn("zeta_e = |g_e|/2", zeta_e, 20)
    pn("a_e = zeta_e - 1", a_e, 20)
    pn("delta zeta_e", delta_zeta_e, 20)

    subbanner("Normalized 22 mode and eta1 = 11/36")
    p("phi_22 normalization integral", norm22)
    p("A_22", A22)
    p("Normalized phi_22", phi22_n)
    p("Canonical self-loop factor", self_loop_factor)
    p("eta1_min", eta1_min)
    p("Sharp zeta series", zeta_sharp)
    p("Sharp g(f)", g_sharp_exact)

    subbanner("Boundary leakage softening")
    p("B_lin(delta)", sp.expand(B_lin))
    p("B_exp(delta)", sp.expand(B_exp))
    p("eta1(delta)", eta1_of_delta)
    pn("Electron-target eta1", eta1_e_target, 20)
    pn("Electron-target overlap B_e", B_e, 20)
    pn("Linear collar delta for B_e", delta_lin_target, 20)
    pn("Exponential collar delta for B_e", delta_exp_target, 20)

    subbanner("Local self-transport PDE")
    p("ODE", ode)
    p("ODE solution with chi(0)=0", ode_sol.rhs)
    p("k_{1/2}", k_half)
    p("ell_blur", ell_blur)
    p("delta_loc", delta_loc)
    p("kappa = 2L/(pi a)", kappa_expr)

    subbanner("Local inertia blur and cubic coefficient")
    p("eta1(f) = (11/36) B_exp(kappa f)", eta1_loc_f)
    p("g_inertia_local(f)", g_local_inertia)
    p("g_inertia_local series", g_local_series)
    p("c3_inertia", c3_inertia)

    subbanner("Charge-side transport closure")
    p("tau(f)", tau)
    p("cbar(f)", cbar)
    p("Xi(f)", Xi)
    p("Xi(f) series", Xi_series)
    p("A_q(f)", Aq)
    p("Q_loc(f)", Q_loc)
    p("Q_loc(f) series", Q_loc_series)
    p("c3_q", c3_q)
    p("c3_total", c3_total)
    print("Final g_loc(f) := 2*(Q_loc(f) - eta1(f)*f**2)")

    subbanner("Benchmark numerics")
    pn("alpha_fs", BENCH.alpha_fs, 20)
    pn("|g_e| target", BENCH.g_e_mag, 20)
    pn("kappa", kappa_num, 20)
    pn("eta1_min = 11/36", eta1_min, 20)
    pn("delta_loc = kappa f", delta_loc_num, 20)
    pn("B_exp(delta_loc)", B_exp_loc, 20)
    pn("eta1_loc", eta1_loc_num, 20)
    pn("c3_q", c3_q.subs(kappa, kappa_num), 20)
    pn("c3_total", c3_total.subs(kappa, kappa_num), 20)
    pn("g_tree", g_tree, 20)
    pn("g_sharp", g_sharp_num, 20)
    pn("g_local_inertia", g_local_inertia_num, 20)
    pn("g_loc", g_loc_num, 20)
    pn("Residual |g_e| - g_loc", residual, 20)
    pn("Quartic coefficient needed to hit benchmark exactly", c4_req, 20)

    return {
        "f_num": f_num,
        "zeta_e": zeta_e,
        "a_e": a_e,
        "eta1_min": eta1_min,
        "eta1_e_target": eta1_e_target,
        "B_e": B_e,
        "delta_lin_target": delta_lin_target,
        "delta_exp_target": delta_exp_target,
        "delta_loc_num": delta_loc_num,
        "B_exp_loc": B_exp_loc,
        "eta1_loc_num": eta1_loc_num,
        "c3_q": c3_q.subs(kappa, kappa_num),
        "c3_total": c3_total.subs(kappa, kappa_num),
        "g_tree": g_tree,
        "g_sharp_num": g_sharp_num,
        "g_local_inertia_num": g_local_inertia_num,
        "g_loc_num": g_loc_num,
        "residual": residual,
        "c4_req": c4_req,
    }


# ---------------------------------------------------------------------------
# Final benchmark table
# ---------------------------------------------------------------------------

def final_summary(hyd: Dict[str, Any], anom: Dict[str, Any]) -> None:
    banner("APPENDIX-H STYLE BENCHMARK SUMMARY")

    summary = [
        ("alpha_fs", BENCH.alpha_fs),
        ("f = alpha_fs/(2 pi)", anom["f_num"]),
        ("|g_e| target", BENCH.g_e_mag),
        ("a_e = |g_e|/2 - 1", anom["a_e"]),
        ("L/a", BENCH.L_over_a),
        ("kappa = 2L/(pi a)", 2 * BENCH.L_over_a / sp.pi),
        ("eta1_min = 11/36", anom["eta1_min"]),
        ("eta1_e_target", anom["eta1_e_target"]),
        ("delta_loc", anom["delta_loc_num"]),
        ("B_exp(delta_loc)", anom["B_exp_loc"]),
        ("eta1_loc", anom["eta1_loc_num"]),
        ("c3_q", anom["c3_q"]),
        ("c3_total", anom["c3_total"]),
        ("g_tree", anom["g_tree"]),
        ("g_sharp", anom["g_sharp_num"]),
        ("g_local_inertia", anom["g_local_inertia_num"]),
        ("g_loc", anom["g_loc_num"]),
        ("Residual |g_e| - g_loc", anom["residual"]),
        ("Fixed-source Bohr radius a0 [m]", hyd["a0_fixed"]),
        ("Reduced-mass hydrogen radius a_H [m]", hyd["a0_reduced"]),
        ("Fixed-source ground-state energy [eV]", hyd["E1_fixed_eV"]),
        ("Reduced-mass ground-state energy [eV]", hyd["E1_reduced_eV"]),
    ]

    for key, val in summary:
        pn(key, val, 20)


def main() -> None:
    banner("FULL REDUCED-SECTOR 4D PAPER DRIVER")
    print(
        "This script symbolically reproduces the derivation chain used across the\n"
        "write-up and prints the benchmark values appearing in the appendices."
    )

    foundations = foundations_snapshot()
    hyd = hydrogenic_reduction()
    finite = finite_size_atomic_response()
    lep = lepton_doublet_bridge()
    g2 = dirac_g2_bridge()
    anom = electron_anomaly_chain()
    final_summary(hyd, anom)


if __name__ == "__main__":
    main()
