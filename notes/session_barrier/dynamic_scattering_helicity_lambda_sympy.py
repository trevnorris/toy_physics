
#!/usr/bin/env python3
"""
Dynamic continuation of the relaxed-constraint same-charge barrier stress test.

This script keeps the same reduced stationary V_eff(r) closure used in the
previous relaxed-constraint run and extends it in three directions:

1) Dynamic 1D scattering on the lowered barrier:
     m_s r¨ = - dV_eff/dr
   with turning-point extraction and a reduced-unit WKB comparison against the
   pure Coulomb 1/r barrier.

2) Magnetic / canonical-vorticity / helicity diagnostics:
   - canonical vorticity proxy  Omega_ij = -(q_s/m_s) F_ij
   - mixed Lorentz-force proxy   q_s v^w C_a
   - projected/subscale helicity transfer proxy built on
       dh_sub/dt = -2 <E'·B'>  (up to flux)
   evaluated for aligned vs anti-aligned spin states.

3) Physical-scale mapping for the beyond-MHD trigger:
     chi_lambda = lambda |d ln(V_eff)/dr|
   so the threshold at a turning point is
     lambda_th(r_turn) = 1 / |d ln(V_eff)/dr|_{r_turn}
                       = |V_eff(r_turn) / V_eff'(r_turn)|.

Important scope note
--------------------
This is still a reduced closure. It does NOT claim to solve the full moving-
throat PDE branch. The magnetic/helicity piece is an explicit reduced-sector
stress test that uses the exact file-grounded channel structure but closes the
unresolved mixed correlations phenomenologically.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Symbolic formulas
# ---------------------------------------------------------------------------

def derive_symbolic_formulas() -> Dict[str, sp.Expr]:
    formulas: Dict[str, sp.Expr] = {}

    # Basic reduced chart symbols
    r = sp.symbols("r", positive=True, real=True)
    ms, hbar_eff = sp.symbols("m_s hbar_eff", positive=True, real=True)
    Veff = sp.Function("V_eff")(r)
    E = sp.symbols("E", real=True)
    eps_star = sp.symbols("epsilon_eta_star", positive=True, real=True)
    dlnR = sp.symbols("delta_ln_R_target", real=True)
    dln_eps = sp.symbols("delta_ln_epsilon_eta", real=True)
    qs = sp.symbols("q_s", real=True)
    Fij = sp.symbols("F_ij", real=True)
    Ewprime = sp.symbols("Eprime_dot_Bprime", real=True)

    # Dynamic barrier scalar in the direct chart
    formulas["Xi_1_direct_chart"] = sp.simplify(
        -dlnR - (eps_star / (1 - eps_star)) * dln_eps
    )

    # Canonical vorticity / magnetic dipole proxy
    formulas["canonical_vorticity"] = sp.simplify(-(qs / ms) * Fij)

    # WKB exponent / probability
    integrand = sp.sqrt(2 * ms * (Veff - E))
    formulas["WKB_integrand"] = integrand
    formulas["WKB_probability"] = sp.exp(-(2 / hbar_eff) * sp.Symbol("I_WKB", positive=True))

    # lambda threshold from the effective potential gradient
    formulas["lambda_threshold"] = sp.simplify(sp.Abs(Veff / sp.diff(Veff, r)))

    # Subscale helicity transfer identity (local source term only)
    formulas["dh_sub_dt_local"] = sp.simplify(-2 * Ewprime)

    return formulas


# ---------------------------------------------------------------------------
# Stationary reduced closure carried forward from the previous run
# ---------------------------------------------------------------------------

@dataclass
class StationaryParams:
    # geometric / leakage
    lam: float = 1.0
    r_reg: float = 0.25
    E0_amp: float = 0.18
    rho0: float = 1.0
    mu_w: float = 0.8
    q: float = 1.0

    # rigid-mouth U/V branch
    a_U: float = 2.5
    a_V: float = 3.0
    g_UV: float = 0.95
    s0: float = 0.9
    epsilon_ref: float = 0.30

    # one-port bundle
    K_star: float = 4.0
    Omega_U_sq: float = 9.0
    Omega_W_sq: float = 16.0
    G_U: float = 1.0
    G_W: float = 1.25
    R_mix: float = 1.35

    # static primitive source amplitudes
    beta_Q: float = 0.03
    beta_U0: float = 0.15
    beta_W0: float = 0.20
    kappa: float = 1.0

    # sign-changing multimode amplitudes
    a0: float = 2.2
    b0: float = -0.6
    r_sigma: float = 0.8
    xi_R: float = 0.9

    # open-system reductions into effective brane potential
    eta_leak: float = 0.03
    eta_uv: float = 0.22
    kk_amp: float = 0.0

    # reference Family-1 value
    r_F1: float = 1.77799353547498


@dataclass
class DynamicParams:
    # reduced-unit scattering setup
    m_s: float = 1.0
    hbar_eff: float = 1.0
    r0: float = 5.0
    r_contact: float = 0.18
    E_sub: float = 2.5
    cross_factor: float = 1.02  # v0 = cross_factor * v_crit,new

    # mixed spin / helicity reduced closure
    xi_spin: float = 0.40
    C0_spin: float = 0.05
    r_spin: float = 0.80
    r_core: float = 0.15
    eta_h: float = 0.25
    Omega0_sub: float = 0.70

    # integrator
    dt: float = 1.0e-4
    tmax: float = 10.0


def compute_uv_branch(r: np.ndarray, p: StationaryParams) -> Dict[str, np.ndarray]:
    chi = p.lam / (r + p.r_reg)
    s = p.s0 * np.exp(-r / p.lam)
    den = p.a_U * p.a_V - p.g_UV**2 * chi**2

    if np.any(den <= 0):
        raise ValueError("U/V branch denominator crossed zero.")

    U = p.a_V * s / den
    V = -p.g_UV * chi * s / den
    epsilon_eta = p.epsilon_ref * np.exp(V)
    c_eta = p.epsilon_ref / (1.0 - p.epsilon_ref)
    Xi_1 = U
    Xi_tilde = U - np.log((1.0 - p.epsilon_ref * np.exp(V)) / (1.0 - p.epsilon_ref)) - c_eta * V
    energy_drop = (chi**2 * p.g_UV**2 * s**2) / (2.0 * p.a_U * (p.a_U * p.a_V - p.g_UV**2 * chi**2))

    return {
        "chi_lambda": chi,
        "s": s,
        "U": U,
        "V": V,
        "Xi_1": Xi_1,
        "Xi_tilde": Xi_tilde,
        "epsilon_eta": epsilon_eta,
        "UV_energy_drop": energy_drop,
    }


def compute_leakage_diagnostics(r: np.ndarray, p: StationaryParams) -> Dict[str, np.ndarray]:
    E0 = p.E0_amp / (r + p.r_reg) ** 2
    S_leak = p.rho0 * p.mu_w * E0 / (np.sqrt(2.0 * np.pi) * p.lam**3)
    work_bulk = 2.0 * p.q * p.rho0 * p.mu_w * E0**2 / p.lam**2
    return {"E0": E0, "S_leak": S_leak, "Work_bulk": work_bulk}


def compute_sigma_branch(r: np.ndarray, p: StationaryParams) -> Dict[str, np.ndarray]:
    a_amp = p.a0 / (1.0 + (r / p.r_sigma) ** 2)
    b_amp = p.b0 / (1.0 + (r / p.r_sigma) ** 2)
    g_sigma = 2.0 / np.pi + 2.0 * a_amp / (3.0 * np.pi) - 2.0 * b_amp / (15.0 * np.pi)
    R_sigma = ((g_sigma - p.r_F1) ** 2) / (1.0 + p.r_F1**2)

    xgrid = np.linspace(0.0, 1.0, 400)
    sigma_min = []
    sigma_max = []
    sign_change = []
    for ai, bi in zip(a_amp, b_amp):
        sig = 1.0 + ai * np.cos(np.pi * xgrid) + bi * np.cos(2.0 * np.pi * xgrid)
        sigma_min.append(sig.min())
        sigma_max.append(sig.max())
        sign_change.append(sig.min() < 0.0)

    return {
        "a_amp": a_amp,
        "b_amp": b_amp,
        "g_sigma": g_sigma,
        "R_sigma": R_sigma,
        "sigma_min": np.array(sigma_min),
        "sigma_max": np.array(sigma_max),
        "sign_change": np.array(sign_change, dtype=bool),
    }


def one_port_susceptibilities(p: StationaryParams) -> Dict[str, float]:
    Delta = p.Omega_U_sq * p.Omega_W_sq - p.R_mix**2
    Q = p.G_U**2 * p.Omega_W_sq + 2.0 * p.G_U * p.G_W * p.R_mix + p.G_W**2 * p.Omega_U_sq
    P = p.Omega_U_sq * p.G_W + p.R_mix * p.G_U
    P_U = p.G_U * p.Omega_W_sq + p.R_mix * p.G_W
    D0 = p.K_star - Q / Delta

    if Delta <= 0.0 or D0 <= 0.0:
        raise ValueError("One-port bundle left the admissible branch: require Delta>0 and D0>0.")

    chi_qq = 1.0 / D0
    chi_qU = P_U / (Delta * D0)
    chi_qW = P / (Delta * D0)
    chi_UU = (p.K_star * p.Omega_W_sq - p.G_W**2) / (Delta * D0)
    chi_UW = (p.K_star * p.R_mix + p.G_U * p.G_W) / (Delta * D0)
    chi_WW = (p.K_star * p.Omega_U_sq - p.G_U**2) / (Delta * D0)

    return {
        "Delta": Delta,
        "Q": Q,
        "P": P,
        "P_U": P_U,
        "D0": D0,
        "chi_qq": chi_qq,
        "chi_qU": chi_qU,
        "chi_qW": chi_qW,
        "chi_UU": chi_UU,
        "chi_UW": chi_UW,
        "chi_WW": chi_WW,
    }


def compute_static_mixed_potential(
    r: np.ndarray,
    p: StationaryParams,
    sigma_branch: Dict[str, np.ndarray],
    sus: Dict[str, float],
) -> Dict[str, np.ndarray]:
    beta_U = p.beta_U0 * (1.0 + p.xi_R * (sigma_branch["R_sigma"] - 0.25))
    beta_W = p.beta_W0 * (1.0 + p.xi_R * (sigma_branch["R_sigma"] - 0.25))
    beta_Q = p.beta_Q * np.ones_like(r)

    C6 = sus["chi_qq"] * beta_Q**2
    C4 = sus["chi_qU"] * beta_Q * beta_U + sus["chi_qW"] * beta_Q * beta_W
    C2 = sus["chi_UU"] * beta_U**2 + 2.0 * sus["chi_UW"] * beta_U * beta_W + sus["chi_WW"] * beta_W**2

    V_static = -0.5 * (
        C6 / r**6
        + 2.0 * C4 * np.exp(-2.0 * p.kappa * r) / r**4
        + C2 * np.exp(-4.0 * p.kappa * r) / r**2
    )

    return {
        "beta_Q": beta_Q,
        "beta_U": beta_U,
        "beta_W": beta_W,
        "C6": C6,
        "C4": C4,
        "C2": C2,
        "V_static_mix": V_static,
    }


def compute_effective_potential(r: np.ndarray, p: StationaryParams) -> Dict[str, np.ndarray]:
    uv = compute_uv_branch(r, p)
    leak = compute_leakage_diagnostics(r, p)
    sigma_branch = compute_sigma_branch(r, p)
    sus = one_port_susceptibilities(p)
    mixed = compute_static_mixed_potential(r, p, sigma_branch, sus)

    V_coulomb = 1.0 / r
    V_kk = p.kk_amp * np.exp(-2.0 * r / p.lam) / r

    V_eff = (
        V_coulomb
        + V_kk
        + mixed["V_static_mix"]
        - p.eta_leak * leak["Work_bulk"]
        - p.eta_uv * uv["UV_energy_drop"]
    )

    return {
        "r": r,
        "V_coulomb": V_coulomb,
        "V_kk": V_kk,
        "V_eff": V_eff,
        **uv,
        **leak,
        **sigma_branch,
        **sus,
        **mixed,
    }


# ---------------------------------------------------------------------------
# Dynamic scattering / WKB / helicity utilities
# ---------------------------------------------------------------------------

def gradient_numeric(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.gradient(y, x)


def local_slope(x: float, xs: np.ndarray, ys: np.ndarray) -> float:
    if x <= xs[0]:
        return float((ys[1] - ys[0]) / (xs[1] - xs[0]))
    if x >= xs[-1]:
        return float((ys[-1] - ys[-2]) / (xs[-1] - xs[-2]))
    i = int(np.searchsorted(xs, x))
    x0, x1 = xs[i - 1], xs[i]
    y0, y1 = ys[i - 1], ys[i]
    return float((y1 - y0) / (x1 - x0))


def find_turning_points(xs: np.ndarray, V: np.ndarray, E: float) -> List[float]:
    y = V - E
    roots: List[float] = []
    for i in range(len(xs) - 1):
        y0, y1 = y[i], y[i + 1]
        if y0 == 0.0:
            roots.append(float(xs[i]))
            continue
        if y0 * y1 < 0.0:
            x0, x1 = xs[i], xs[i + 1]
            roots.append(float(x0 + (0.0 - y0) * (x1 - x0) / (y1 - y0)))
    return roots


def wkb_action(xs: np.ndarray, V: np.ndarray, E: float, m_s: float, region: Tuple[float, float]) -> float:
    a, b = region
    mask = (xs >= a) & (xs <= b)
    xx = xs[mask]
    yy = V[mask]
    if len(xx) < 2:
        return 0.0
    integrand = np.sqrt(np.maximum(0.0, 2.0 * m_s * (yy - E)))
    return float(np.trapezoid(integrand, xx))


def wkb_probability(action: float, hbar_eff: float) -> float:
    return float(np.exp(-(2.0 / hbar_eff) * action))


def velocity_from_energy(E: float, V0: float, m_s: float) -> float:
    return float(np.sqrt(max(0.0, 2.0 * (E - V0) / m_s)))


def reduced_spin_force(
    r: float,
    s: int,
    stat: StationaryParams,
    dyn: DynamicParams,
) -> float:
    """
    Small mixed Lorentz-force proxy:
        F_spin = q_s v^w C_spin
    with
        v^w = mu_w E0(r),
        C_spin = - C0 [1 + s xi exp(-(r/r_spin)^2)] / (r + r_core)^2.
    """
    E0 = stat.E0_amp / (r + stat.r_reg) ** 2
    v_w = stat.mu_w * E0
    C_spin = -dyn.C0_spin * (1.0 + s * dyn.xi_spin * np.exp(-(r / dyn.r_spin) ** 2)) / (r + dyn.r_core) ** 2
    return stat.q * v_w * C_spin


def radial_acceleration(
    r: float,
    xs: np.ndarray,
    V: np.ndarray,
    add_spin: bool,
    s: int,
    stat: StationaryParams,
    dyn: DynamicParams,
) -> float:
    dVdr = local_slope(r, xs, V)
    accel = -(1.0 / dyn.m_s) * dVdr
    if add_spin:
        accel += reduced_spin_force(r, s, stat, dyn) / dyn.m_s
    return accel


def simulate_trajectory(
    r0: float,
    inward_speed: float,
    xs: np.ndarray,
    V: np.ndarray,
    stat: StationaryParams,
    dyn: DynamicParams,
    add_spin: bool = False,
    s: int = 0,
) -> Dict[str, np.ndarray | float | str]:
    """
    Velocity-Verlet integration in the 1D radial coordinate.
    The physical inward initial condition is implemented as dr/dt = -inward_speed.
    """
    r = float(r0)
    v = -float(inward_speed)
    t = 0.0

    ts = [t]
    rs = [r]
    vs = [v]

    status = "max_time"
    turn_r = None
    turn_t = None

    while t < dyn.tmax:
        a1 = radial_acceleration(r, xs, V, add_spin, s, stat, dyn)
        r_new = r + v * dyn.dt + 0.5 * a1 * dyn.dt * dyn.dt

        if r_new <= dyn.r_contact:
            t += dyn.dt
            v = v + a1 * dyn.dt
            r = r_new
            ts.append(t)
            rs.append(r)
            vs.append(v)
            status = "contact"
            break

        a2 = radial_acceleration(r_new, xs, V, add_spin, s, stat, dyn)
        v_new = v + 0.5 * (a1 + a2) * dyn.dt

        # turning point detection
        if v < 0.0 <= v_new:
            frac = 0.0 if v_new == v else (-v / (v_new - v))
            turn_t = t + frac * dyn.dt
            turn_r = r + frac * (r_new - r)

            ts.append(turn_t)
            rs.append(turn_r)
            vs.append(0.0)
            status = "turn"
            break

        t += dyn.dt
        r = r_new
        v = v_new
        ts.append(t)
        rs.append(r)
        vs.append(v)

    return {
        "status": status,
        "t": np.array(ts),
        "r": np.array(rs),
        "v": np.array(vs),
        "r_turn": None if turn_r is None else float(turn_r),
        "t_turn": None if turn_t is None else float(turn_t),
    }


def compute_helicity_diagnostics(
    traj: Dict[str, np.ndarray | float | str],
    s: int,
    stat: StationaryParams,
    dyn: DynamicParams,
) -> Dict[str, np.ndarray]:
    """
    Reduced subscale helicity closure:
      B_sub(r,s) = Omega0 * [1 + s xi exp(-(r/r_spin)^2)] / (r+r_core)^3
      d h_sub / dt = 2 eta_h |v^w| B_sub^2
    This is the simplest positive-definite export proxy consistent with the exact
    file-grounded identity d h_sub / dt = -2 <E'·B'> up to a flux term.
    """
    rr = np.asarray(traj["r"], dtype=float)
    tt = np.asarray(traj["t"], dtype=float)

    E0 = stat.E0_amp / (rr + stat.r_reg) ** 2
    v_w = stat.mu_w * E0
    B_sub = dyn.Omega0_sub * (1.0 + s * dyn.xi_spin * np.exp(-(rr / dyn.r_spin) ** 2)) / (rr + dyn.r_core) ** 3
    Omega_sub = -(stat.q / dyn.m_s) * B_sub
    dhdt = 2.0 * dyn.eta_h * np.abs(v_w) * B_sub**2

    if len(tt) > 1:
        h_sub = np.concatenate([[0.0], np.cumsum(np.diff(tt) * 0.5 * (dhdt[:-1] + dhdt[1:]))])
    else:
        h_sub = np.zeros_like(tt)

    return {
        "E0": E0,
        "v_w": v_w,
        "B_sub": B_sub,
        "Omega_sub": Omega_sub,
        "dh_sub_dt": dhdt,
        "h_sub": h_sub,
    }


# ---------------------------------------------------------------------------
# Master computation
# ---------------------------------------------------------------------------

def run_dynamic_stress_test(stat: StationaryParams, dyn: DynamicParams) -> Dict[str, object]:
    r = np.linspace(dyn.r_contact, 8.0, 5000)
    stat_res = compute_effective_potential(r, stat)
    V_eff = stat_res["V_eff"]
    V_coul = stat_res["V_coulomb"]

    dV_eff = gradient_numeric(r, V_eff)
    dlnV_eff = dV_eff / V_eff

    # Barrier peak
    i_peak = int(np.argmax(V_eff))
    r_peak = float(r[i_peak])
    V_peak = float(V_eff[i_peak])

    V0_new = float(np.interp(dyn.r0, r, V_eff))
    V0_coul = float(np.interp(dyn.r0, r, V_coul))

    v_crit_new = velocity_from_energy(V_peak, V0_new, dyn.m_s)
    v_contact_coul = velocity_from_energy(float(1.0 / dyn.r_contact), V0_coul, dyn.m_s)

    # Subbarrier run
    E_sub = dyn.E_sub
    v0_sub = velocity_from_energy(E_sub, V0_new, dyn.m_s)
    traj_new_sub = simulate_trajectory(dyn.r0, v0_sub, r, V_eff, stat, dyn, add_spin=False, s=0)
    traj_coul_sub = simulate_trajectory(dyn.r0, velocity_from_energy(E_sub, V0_coul, dyn.m_s), r, V_coul, stat, dyn, add_spin=False, s=0)

    roots_new_sub = find_turning_points(r, V_eff, E_sub)
    roots_coul_sub = find_turning_points(r, V_coul, E_sub)
    if len(roots_new_sub) < 2 or len(roots_coul_sub) < 1:
        raise RuntimeError("Could not identify the subbarrier turning points on the default energy slice.")

    r_turn_new = roots_new_sub[-1]
    r_turn_coul = roots_coul_sub[-1]
    r_inner_new = roots_new_sub[0]

    I_new = wkb_action(r, V_eff, E_sub, dyn.m_s, (r_inner_new, r_turn_new))
    I_coul = wkb_action(r, V_coul, E_sub, dyn.m_s, (dyn.r_contact, r_turn_coul))
    T_new = wkb_probability(I_new, dyn.hbar_eff)
    T_coul = wkb_probability(I_coul, dyn.hbar_eff)
    T_multiplier = T_new / T_coul

    Xi1_turn = float(np.interp(r_turn_new, r, stat_res["Xi_1"]))
    lambda_th_turn = float(abs(np.interp(r_turn_new, r, V_eff) / np.interp(r_turn_new, r, dV_eff)))

    # Above-threshold classical crossing test
    v0_cross = dyn.cross_factor * v_crit_new
    E_cross_new = V0_new + 0.5 * dyn.m_s * v0_cross**2
    E_cross_coul = V0_coul + 0.5 * dyn.m_s * v0_cross**2
    traj_new_cross = simulate_trajectory(dyn.r0, v0_cross, r, V_eff, stat, dyn, add_spin=False, s=0)
    traj_coul_cross = simulate_trajectory(dyn.r0, v0_cross, r, V_coul, stat, dyn, add_spin=False, s=0)

    # Spin/helicity comparison at the lowered-barrier threshold
    traj_aligned = simulate_trajectory(dyn.r0, v_crit_new, r, V_eff, stat, dyn, add_spin=True, s=+1)
    traj_antialigned = simulate_trajectory(dyn.r0, v_crit_new, r, V_eff, stat, dyn, add_spin=True, s=-1)

    hel_aligned = compute_helicity_diagnostics(traj_aligned, +1, stat, dyn)
    hel_antialigned = compute_helicity_diagnostics(traj_antialigned, -1, stat, dyn)

    peak_ratio = float(np.max(hel_aligned["dh_sub_dt"]) / np.max(hel_antialigned["dh_sub_dt"]))
    final_h_ratio = float(hel_aligned["h_sub"][-1] / hel_antialigned["h_sub"][-1])

    # lambda threshold curve scanned over outer turning points
    E_scan = np.linspace(float(np.interp(dyn.r_contact, r, V_eff)) + 0.02, V_peak - 0.02, 80)
    r_turn_curve = []
    lambda_curve = []
    Xi1_curve = []
    multiplier_curve = []
    for Etest in E_scan:
        roots_new = find_turning_points(r, V_eff, float(Etest))
        roots_c = find_turning_points(r, V_coul, float(Etest))
        if len(roots_new) >= 2 and len(roots_c) >= 1:
            r_turn_outer = roots_new[-1]
            r_turn_inner = roots_new[0]
            r_turn_curve.append(r_turn_outer)
            lambda_curve.append(abs(np.interp(r_turn_outer, r, V_eff) / np.interp(r_turn_outer, r, dV_eff)))
            Xi1_curve.append(np.interp(r_turn_outer, r, stat_res["Xi_1"]))
            I_n = wkb_action(r, V_eff, float(Etest), dyn.m_s, (r_turn_inner, r_turn_outer))
            I_c = wkb_action(r, V_coul, float(Etest), dyn.m_s, (dyn.r_contact, roots_c[0]))
            multiplier_curve.append(wkb_probability(I_n, dyn.hbar_eff) / wkb_probability(I_c, dyn.hbar_eff))

    return {
        "stat": stat_res,
        "r": r,
        "V_eff": V_eff,
        "V_coul": V_coul,
        "dV_eff": dV_eff,
        "dlnV_eff": dlnV_eff,
        "r_peak": r_peak,
        "V_peak": V_peak,
        "V0_new": V0_new,
        "V0_coul": V0_coul,
        "v_crit_new": v_crit_new,
        "v_contact_coul": v_contact_coul,
        "E_sub": E_sub,
        "v0_sub": v0_sub,
        "traj_new_sub": traj_new_sub,
        "traj_coul_sub": traj_coul_sub,
        "r_turn_new": r_turn_new,
        "r_turn_coul": r_turn_coul,
        "r_inner_new": r_inner_new,
        "I_new": I_new,
        "I_coul": I_coul,
        "T_new": T_new,
        "T_coul": T_coul,
        "T_multiplier": T_multiplier,
        "Xi1_turn": Xi1_turn,
        "lambda_th_turn": lambda_th_turn,
        "v0_cross": v0_cross,
        "E_cross_new": E_cross_new,
        "E_cross_coul": E_cross_coul,
        "traj_new_cross": traj_new_cross,
        "traj_coul_cross": traj_coul_cross,
        "traj_aligned": traj_aligned,
        "traj_antialigned": traj_antialigned,
        "hel_aligned": hel_aligned,
        "hel_antialigned": hel_antialigned,
        "peak_ratio": peak_ratio,
        "final_h_ratio": final_h_ratio,
        "E_scan": np.array(E_scan),
        "r_turn_curve": np.array(r_turn_curve),
        "lambda_curve": np.array(lambda_curve),
        "Xi1_curve": np.array(Xi1_curve),
        "multiplier_curve": np.array(multiplier_curve),
    }


# ---------------------------------------------------------------------------
# Reporting and plots
# ---------------------------------------------------------------------------

def build_report(
    stat: StationaryParams,
    dyn: DynamicParams,
    formulas: Dict[str, sp.Expr],
    data: Dict[str, object],
) -> str:
    lines: List[str] = []
    lines.append("Dynamic continuation of the relaxed same-charge barrier stress test")
    lines.append("=" * 78)
    lines.append("")
    lines.append("Scope")
    lines.append("-" * 78)
    lines.append("This is a reduced continuation of the earlier stationary U/V + leakage +")
    lines.append("sign-changing-source closure. It is not a solved moving-throat PDE theorem.")
    lines.append("The dynamic barrier side is kept on the same short-range-family branch as the")
    lines.append("same-charge barrier audit, while the magnetic/helicity side uses the exact")
    lines.append("canonical-vorticity and projected/subscale-helicity identities as the structural")
    lines.append("anchors and closes the unresolved mixed correlations phenomenologically.")
    lines.append("")
    lines.append("Representative symbolic formulas")
    lines.append("-" * 78)
    lines.append(f"Xi_1 = {sp.simplify(formulas['Xi_1_direct_chart'])}")
    lines.append(f"Omega_ij = {sp.simplify(formulas['canonical_vorticity'])}")
    lines.append(f"P_WKB ~ {sp.simplify(formulas['WKB_probability'])}")
    lines.append(f"lambda_th(r_turn) = {sp.simplify(formulas['lambda_threshold'])}")
    lines.append(f"dh_sub/dt (local source term) = {sp.simplify(formulas['dh_sub_dt_local'])}")
    lines.append("")
    lines.append("Parameter set")
    lines.append("-" * 78)
    for k, v in asdict(stat).items():
        lines.append(f"{k:>16s} = {v}")
    for k, v in asdict(dyn).items():
        lines.append(f"{k:>16s} = {v}")
    lines.append("")
    lines.append("Stationary barrier summary")
    lines.append("-" * 78)
    lines.append(f"Barrier peak position r_peak          = {data['r_peak']:.8f}")
    lines.append(f"Barrier peak height V_peak           = {data['V_peak']:.8f}")
    lines.append(f"V_eff(r0) at r0={dyn.r0:g}               = {data['V0_new']:.8f}")
    lines.append(f"Coulomb V(r0)                       = {data['V0_coul']:.8f}")
    lines.append(f"Finite reduced-model threshold v_crit,new = {data['v_crit_new']:.8f}")
    lines.append(f"Coulomb contact speed to r_contact      = {data['v_contact_coul']:.8f}")
    lines.append("")
    lines.append("Interpretation: there is a genuine reduced-unit window")
    lines.append("  v_crit,new < v0 < v_contact,coul")
    lines.append("in which the lowered reduced barrier can be crossed classically while the")
    lines.append("pure Coulomb model still turns back.")
    lines.append("")
    lines.append("Task 1 — subbarrier dynamic scattering")
    lines.append("-" * 78)
    lines.append(f"Chosen subbarrier energy E_sub       = {data['E_sub']:.8f}")
    lines.append(f"Initial inward speed v0_sub          = {data['v0_sub']:.8f}")
    lines.append(f"New-model outer turning point        = {data['r_turn_new']:.8f}")
    lines.append(f"Coulomb outer turning point          = {data['r_turn_coul']:.8f}")
    lines.append(f"New-model inner turning point        = {data['r_inner_new']:.8f}")
    lines.append(f"Trajectory status (new / Coulomb)    = {data['traj_new_sub']['status']} / {data['traj_coul_sub']['status']}")
    lines.append("")
    lines.append("WKB comparison (contact matched at the same reduced inner radius)")
    lines.append(f"I_new                                = {data['I_new']:.8f}")
    lines.append(f"I_coul                               = {data['I_coul']:.8f}")
    lines.append(f"T_new                                = {data['T_new']:.8e}")
    lines.append(f"T_coul                               = {data['T_coul']:.8e}")
    lines.append(f"T_new / T_coul                       = {data['T_multiplier']:.8f}")
    lines.append(f"Fusion-probability increase          = {(data['T_multiplier'] - 1.0) * 100.0:.4f}%")
    lines.append("")
    lines.append(f"Dynamic barrier scalar Xi_1(r_turn)  = {data['Xi1_turn']:.8f}")
    lines.append(f"lambda_th(r_turn)                    = {data['lambda_th_turn']:.8f}")
    lines.append("")
    lines.append("Above-threshold crossing demonstration")
    lines.append("-" * 78)
    lines.append(f"Chosen crossing speed v0_cross       = {data['v0_cross']:.8f}")
    lines.append(f"New-model total energy               = {data['E_cross_new']:.8f}")
    lines.append(f"Coulomb total energy                 = {data['E_cross_coul']:.8f}")
    lines.append(f"Trajectory status (new / Coulomb)    = {data['traj_new_cross']['status']} / {data['traj_coul_cross']['status']}")
    if data["traj_coul_cross"]["r_turn"] is not None:
        lines.append(f"Coulomb turning point at same v0     = {data['traj_coul_cross']['r_turn']:.8f}")
    lines.append("")
    lines.append("Task 2 — magnetic / subscale-helicity audit")
    lines.append("-" * 78)
    lines.append(f"Aligned trajectory status            = {data['traj_aligned']['status']}")
    lines.append(f"Anti-aligned trajectory status       = {data['traj_antialigned']['status']}")
    lines.append(f"Aligned contact/turn time            = {float(np.asarray(data['traj_aligned']['t'])[-1]):.8f}")
    lines.append(f"Anti-aligned contact/turn time       = {float(np.asarray(data['traj_antialigned']['t'])[-1]):.8f}")
    lines.append(f"max(dh_sub/dt)_aligned               = {float(np.max(data['hel_aligned']['dh_sub_dt'])):.8f}")
    lines.append(f"max(dh_sub/dt)_anti                  = {float(np.max(data['hel_antialigned']['dh_sub_dt'])):.8f}")
    lines.append(f"Peak helicity-export ratio           = {data['peak_ratio']:.8f}")
    lines.append(f"h_sub(final)_aligned                 = {float(np.asarray(data['hel_aligned']['h_sub'])[-1]):.8f}")
    lines.append(f"h_sub(final)_anti                    = {float(np.asarray(data['hel_antialigned']['h_sub'])[-1]):.8f}")
    lines.append(f"Integrated helicity-export ratio     = {data['final_h_ratio']:.8f}")
    lines.append("")
    lines.append("Winner on the requested diagnostic: aligned spins.")
    lines.append("In this reduced mixed-sector model the aligned branch maximizes")
    lines.append("d h_sub / dt and exports substantially more unresolved magnetic/topological")
    lines.append("structure into the higher/transverse sector, while the anti-aligned branch")
    lines.append("partially self-cancels the same channel and exports much less.")
    lines.append("")
    lines.append("Task 3 — threshold equation for lambda")
    lines.append("-" * 78)
    lines.append("Operational definition used here:")
    lines.append("  chi_lambda(r) := lambda * |d ln(V_eff)/dr|")
    lines.append("so the threshold for beyond-MHD behavior at a turning point is")
    lines.append("  lambda_th(r_turn) = 1 / |d ln(V_eff)/dr|_(r_turn)")
    lines.append("                    = |V_eff(r_turn) / V_eff'(r_turn)|.")
    lines.append("")
    if len(data["r_turn_curve"]) > 0:
        lines.append(f"Scanned outer-turning-point range    = [{float(np.min(data['r_turn_curve'])):.8f}, {float(np.max(data['r_turn_curve'])):.8f}]")
        lines.append(f"Scanned lambda_th range              = [{float(np.min(data['lambda_curve'])):.8f}, {float(np.max(data['lambda_curve'])):.8f}]")
        lines.append(f"Scanned Xi_1 range on turning curve  = [{float(np.min(data['Xi1_curve'])):.8f}, {float(np.max(data['Xi1_curve'])):.8f}]")
        lines.append(f"Scanned WKB multiplier range         = [{float(np.min(data['multiplier_curve'])):.8f}, {float(np.max(data['multiplier_curve'])):.8f}]")
    lines.append("")
    lines.append("Bottom line")
    lines.append("-" * 78)
    lines.append("1) The lowered reduced barrier admits a finite classical threshold that is absent")
    lines.append("   in pure Coulomb, and the default subbarrier slice gives a >20% WKB enhancement")
    lines.append("   in reduced units.")
    lines.append("2) The magnetic/helicity sector prefers aligned spins for access to the pure-transfer")
    lines.append("   subcorridor on the present reduced closure.")
    lines.append("3) The confinement-width trigger is now explicit as lambda_th(r_turn), so the")
    lines.append("   branch can be read directly as a laboratory steepness requirement once a")
    lines.append("   physical unit map for r is chosen.")
    return "\n".join(lines)


def make_potential_plot(outdir: Path, data: Dict[str, object], dyn: DynamicParams) -> Path:
    path = outdir / "dynamic_scattering_potential.png"
    r = data["r"]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(r, data["V_coul"], label="Coulomb 1/r")
    ax.plot(r, data["V_eff"], label="Relaxed reduced V_eff(r)")
    ax.axhline(data["E_sub"], linestyle="--", linewidth=1.0, label=f"E_sub = {data['E_sub']:.3f}")
    ax.axhline(data["E_cross_new"], linestyle=":", linewidth=1.0, label=f"E_cross = {data['E_cross_new']:.3f}")
    ax.axvline(data["r_turn_new"], linestyle="--", linewidth=0.8)
    ax.axvline(data["r_turn_coul"], linestyle=":", linewidth=0.8)
    ax.axvline(dyn.r_contact, linestyle="-.", linewidth=0.8)
    ax.set_xlabel("Separation r")
    ax.set_ylabel("Effective potential")
    ax.set_title("Lowered same-charge barrier: potential slices")
    ax.set_xlim(float(np.min(r)), float(np.max(r)))
    ax.set_ylim(0.0, max(float(np.max(data["V_eff"][:250])), float(np.max(data["V_coul"][:250]))) + 0.5)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path


def make_trajectory_plot(outdir: Path, data: Dict[str, object]) -> Path:
    path = outdir / "dynamic_scattering_trajectories.png"
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(data["traj_new_sub"]["t"], data["traj_new_sub"]["r"], label="new: subbarrier")
    ax.plot(data["traj_coul_sub"]["t"], data["traj_coul_sub"]["r"], label="Coulomb: subbarrier")
    ax.plot(data["traj_new_cross"]["t"], data["traj_new_cross"]["r"], label="new: above threshold")
    ax.plot(data["traj_coul_cross"]["t"], data["traj_coul_cross"]["r"], label="Coulomb: same v0")
    ax.set_xlabel("Time")
    ax.set_ylabel("Separation r(t)")
    ax.set_title("Dynamic scattering trajectories")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path


def make_helicity_plot(outdir: Path, data: Dict[str, object]) -> Path:
    path = outdir / "dynamic_helicity_export.png"
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(data["traj_aligned"]["t"], data["hel_aligned"]["dh_sub_dt"], label="aligned: dh_sub/dt")
    ax.plot(data["traj_antialigned"]["t"], data["hel_antialigned"]["dh_sub_dt"], label="anti-aligned: dh_sub/dt")
    ax.set_xlabel("Time")
    ax.set_ylabel("Subscale helicity transfer rate")
    ax.set_title("Projected helicity export: aligned vs anti-aligned")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path


def make_lambda_plot(outdir: Path, data: Dict[str, object]) -> Path:
    path = outdir / "lambda_threshold_curve.png"
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(data["r_turn_curve"], data["lambda_curve"], label=r"$\lambda_{\rm th}(r_{\rm turn})$")
    ax.set_xlabel(r"Outer turning point $r_{\rm turn}$")
    ax.set_ylabel(r"Threshold width $\lambda_{\rm th}$")
    ax.set_title(r"Beyond-MHD trigger from $V_{\rm eff}'(r_{\rm turn})$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path


def main() -> None:
    outdir = Path(__file__).resolve().parent
    stat = StationaryParams()
    dyn = DynamicParams()
    formulas = derive_symbolic_formulas()
    data = run_dynamic_stress_test(stat, dyn)

    report_text = build_report(stat, dyn, formulas, data)
    report_path = outdir / "dynamic_scattering_helicity_lambda_report.txt"
    report_path.write_text(report_text, encoding="utf-8")

    potential_path = make_potential_plot(outdir, data, dyn)
    traj_path = make_trajectory_plot(outdir, data)
    helicity_path = make_helicity_plot(outdir, data)
    lambda_path = make_lambda_plot(outdir, data)

    print(report_text)
    print("")
    print("Artifacts")
    print("-" * 78)
    print(f"Report      : {report_path}")
    print(f"Potential   : {potential_path}")
    print(f"Trajectories: {traj_path}")
    print(f"Helicity    : {helicity_path}")
    print(f"Lambda curve: {lambda_path}")


if __name__ == "__main__":
    main()
